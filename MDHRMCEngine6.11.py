"""
MDHRMCEngine implements MD-HRMC model using LAMMPS ReaxFF and PDF constraints.
"""
from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from fullrmc.Engine import Engine 
from fullrmc.Globals import LOGGER
from fullrmc.Core.boundary_conditions_collection import transform_coordinates
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.EnergyConstraints import EnergyConstraint

try:
    from lammps import lammps
    LAMMPS_IMPORT = True
except ImportError:
    LAMMPS_IMPORT = False

# Element masses (in g/mol)
ELEMENT_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.380
}    

class MDHRMCEngine(Engine):
    def __init__(self, path, freshStart=False):
        super(MDHRMCEngine, self).__init__(path=path, freshStart=freshStart)
        
        # Check LAMMPS availability
        if not LAMMPS_IMPORT:
            raise ImportError(LOGGER.error("LAMMPS Python interface not found"))
            
        # Constants - Use proper units for LAMMPS real units
        self.KB = 0.001987204258  # Boltzmann constant in kcal/mol/K
        
        # MD-HRMC parameters
        self.md_temperature = None
        self.md_steps = None
        self.hrmc_initial_temp = None
        self.hrmc_final_temp = None 
        self.hrmc_steps = None
        self.cycles_md_hrmc = None
        self.current_cycle = 1

        # Counters for statistics
        self._tried = 1 # Start at 1
        self._accepted = 0

        # Store constraints
        self.pdf_constraint = None
        self.energy_constraint = None
                
        # LAMMPS setup
        self.__lmp = None
        self.__initialized = False
        self.__element_types = {}
        self.__cached_box = None

        # Initialize counters (will be set in set_pdb)
        self.N_total = 0
        self.N_nonH = 0
        self.omega = None

        # Initialize output directory and files
        self._output_dir = os.path.dirname(path) if path else None 
        self._files = {}

        # Setup output files if path is provided
        if self._output_dir:
            self.setup_output_files(self._output_dir)

        # Store last accepted energy
        self.__last_energy = None        

    def setup_output_files(self, output_dir):
        """Setup output directory and files."""
        try:
            self._output_dir = output_dir
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # Define files and headers
            self._files = {
                'HRMC_acceptance.txt': "# Step Probability DeltaChi2 DeltaU Effective-Temperature",
                'HRMC_energy.txt': "# Step Energy DeltaU",
                'MD_energy.txt': "# Cycle PE Temperature"
            }

            # Create files with headers
            for filename, header in self._files.items():
                filepath = os.path.join(output_dir, filename)
                with open(filepath, 'w') as f:
                    f.write(header + '\n')
                    
        except Exception as e:
            LOGGER.error(f"Failed to setup output files: {str(e)}")
            raise

    def _write_data(self, filename, data):
        """Write data to file using context manager."""
        try:
            if self._output_dir is None:
                LOGGER.error("Output directory not set")
                return
                
            filepath = os.path.join(self._output_dir, filename)
            LOGGER.debug(f"Attempting to write to {filepath}")

            # Verify file exists
            if not os.path.exists(filepath):
                with open(filepath, 'w') as f:
                    f.write(self._files[filename] + '\n')
                    
            # Append data
            with open(filepath, 'a') as f:
                f.write(data + '\n')

            LOGGER.debug(f"Successfully wrote data to {filepath}")    
                
        except Exception as e:
            LOGGER.error(f"Error writing to {filename}: {str(e)}")

    def __getstate__(self):
        """Custom serialization."""
        state = super(MDHRMCEngine, self).__getstate__()
        # Don't serialize LAMMPS instance
        state['_MDHRMCEngine__lmp'] = None
        state['_MDHRMCEngine__initialized'] = False
        state['_MDHRMCEngine__last_energy'] = self.__last_energy
        return state
        
    def __setstate__(self, state):
        """Handle deserialization."""
        self.__dict__.update(state)
        # Reinitialize LAMMPS if needed
        self.__lmp = None
        self.__initialized = False

    def cleanup(self):
        """Clean up resources."""
        if hasattr(self, '_MDHRMCEngine__lmp') and self._MDHRMCEngine__lmp:
            self._MDHRMCEngine__lmp.close()
            self._MDHRMCEngine__lmp = None

    def increment_tried(self):
        """Increment number of tried moves starting from 1."""
        if self._tried == 0:  # If 0, start from 1
            self._tried = 1   # This is the first move
        else:
            self._tried += 1  # Increment tried moves

    def increment_accepted(self):
        """Increment number of accepted moves."""  
        self._accepted += 1

    def set_pdb(self, pdb, *args, **kwargs):
        """Override to calculate omega after structure is loaded."""
        # Call parent implementation first
        super(MDHRMCEngine, self).set_pdb(pdb, *args, **kwargs)
        
        # Skip omega calculation if no structure loaded yet
        if pdb is None:
            return
            
        # Calculate atom counts and omega
        elements = self.allElements
        self.N_total = len(elements)
        self.N_nonH = sum(1 for elem in elements if elem != 'H')
        
        if self.N_nonH == 0:
            raise ValueError("Structure must contain non-hydrogen atoms")
            
        # Calculate omega
        self.omega = (self.N_total**2 / self.N_nonH) * 7 * 23.0609
        
        LOGGER.info(f"Structure loaded with {self.N_total} total atoms ({self.N_nonH} non-hydrogen)")
        LOGGER.info(f"Calculated omega: {self.omega}")

    def add_constraints(self, constraints):
        """Add constraints with proper initialization."""
        try:
            for constraint in constraints:
                if isinstance(constraint, PairDistributionConstraint):
                    self.pdf_constraint = constraint
                elif isinstance(constraint, EnergyConstraint):
                    self.energy_constraint = constraint
            
            # Call parent's add_constraints AFTER storing references
            super(MDHRMCEngine, self).add_constraints(constraints)
            
        except Exception as e:
            LOGGER.error(f"Error adding constraints: {str(e)}")
            raise
                
    def initialize_used_constraints(self):
        """Initialize constraints safely."""
        try:
            # 先初始化所有约束数据
            super(MDHRMCEngine, self).initialize_used_constraints()
        
            # 确保能量约束被正确初始化
            if hasattr(self, 'energy_constraint') and self.energy_constraint is not None:
                # 检查last_accepted_energy是否为None
                if self.energy_constraint._EnergyConstraint__last_accepted_energy is None:
                    LOGGER.info("手动初始化能量约束...")
                    # 显式调用energy_constraint的initialize方法
                    self.energy_constraint.initialize()
                
            # 验证初始化是否成功
            if hasattr(self, 'energy_constraint') and self.energy_constraint._EnergyConstraint__last_accepted_energy is None:
                LOGGER.error("能量约束初始化失败，last_accepted_energy仍为None")
                raise ValueError("能量约束初始化失败")
            
        except Exception as e:
            LOGGER.error(f"约束初始化出错: {str(e)}")
            raise
            
    def setup_md(self, temperature, steps, reaxff_params, control_params):
        """Setup MD parameters."""
        self.md_temperature = temperature
        self.md_steps = steps
        
        # Initialize LAMMPS if needed
        if self.__lmp is None:
            self.__lmp = lammps(cmdargs=["-screen", "none", "-log", "none"])
            
            # Get unique elements and create type mapping
            elements = sorted(set(self.allElements))
            self.__element_types = {elem: idx+1 for idx, elem in enumerate(elements)}
        
            # Basic LAMMPS setup
            lmp = self.__lmp
            lmp.command("units real")
            lmp.command("atom_style charge")
            lmp.command("boundary p p p")

            # Create simulation box
            box = self.boundaryConditions.get_vectors()
            box_lengths = np.array([np.linalg.norm(v) for v in box])
            
            LOGGER.info(f"Creating simulation box with dimensions: {box_lengths}")
            lmp.command(f"region box block 0 {box_lengths[0]} 0 {box_lengths[1]} 0 {box_lengths[2]}")
            lmp.command(f"create_box {len(self.__element_types)} box")
            
            # Set element masses
            for elem, type_id in self.__element_types.items():
                if elem not in ELEMENT_MASSES:
                    raise ValueError(f"Mass not defined for element: {elem}")
                mass = ELEMENT_MASSES[elem]
                lmp.command(f"mass {type_id} {mass}")
            
            # Setup ReaxFF
            LOGGER.info("Setting up ReaxFF...")
            lmp.command(f"pair_style reaxff {control_params} safezone 6 mincap 100")
            lmp.command("neighbor 2.0 bin")
            lmp.command("neigh_modify delay 0 every 1 check yes")
            
            # Setup pair coefficients
            elem_string = " ".join(elements)
            lmp.command(f"pair_coeff * * {reaxff_params} {elem_string}")

    def setup_hrmc(self, initial_temp, final_temp, steps, cycles):
        """Setup HRMC parameters."""
        self.hrmc_initial_temp = initial_temp
        self.hrmc_final_temp = final_temp
        self.hrmc_steps = steps
        self.cycles_md_hrmc = cycles
        
        # Count atoms for omega calculation
        elements = self.allElements
        self.N_total = len(elements)
        self.N_nonH = sum(1 for elem in elements if elem != 'H')
        
        # Calculate T_chi_0 and cooling rate
        self.T_chi_0 = self.KB * initial_temp / self.omega
        self.cooling_rate = (final_temp/initial_temp)**(1.0/cycles)

    def run(self, numberOfSteps):
        """
        运行HRMC步骤。
        遵循Engine标准的移动接受/拒绝流程，并使用PDF和能量约束计算误差。
        
        :Parameters:
            #. numberOfSteps (int): 要运行的HRMC步骤数量
        """
        try:
            for step in range(numberOfSteps):
                # 获取要移动的组
                group_index = self.groupSelector.select_index()
                if group_index is None:
                    continue
                    
                group = self.groups[group_index]
                
                # 转换索引为numpy数组
                real_indexes = np.array(group.indexes, dtype=np.int32)
                relative_indexes = np.array([
                    self._atomsCollector.get_relative_index(idx) 
                    for idx in group.indexes
                ], dtype=np.int32)
                
                # 存储旧坐标
                old_positions = self.realCoordinates[relative_indexes].copy()
                
                # 应用移动（遵循Engine流程）
                moved_positions = group.moveGenerator.move(old_positions)
                
                # 应用周期性边界条件
                if self.isPBC:
                    box_coords = transform_coordinates(
                        transMatrix=self.reciprocalBasisVectors,
                        coords=moved_positions
                    )
                    box_coords = box_coords % 1.0
                    moved_positions = transform_coordinates(
                        transMatrix=self.basisVectors,
                        coords=box_coords
                    )
                
                # 移动前计算约束
                for c in self.constraints:
                    c.compute_before_move(
                        realIndexes=real_indexes,
                        relativeIndexes=relative_indexes
                    )
                
                # 更新坐标
                self.realCoordinates[relative_indexes] = moved_positions
                
                # 移动后计算约束
                box_coords = transform_coordinates(
                    transMatrix=self.reciprocalBasisVectors,
                    coords=moved_positions
                )
                for c in self.constraints:
                    c.compute_after_move(
                        realIndexes=real_indexes,
                        relativeIndexes=relative_indexes,
                        movedBoxCoordinates=box_coords
                    )
                
                # 检查是否应该接受移动（使用PDF和能量约束的误差）
                if self.should_step_get_rejected(real_indexes, relative_indexes):
                    # 拒绝移动
                    self.realCoordinates[relative_indexes] = old_positions
                    self.groupSelector.move_rejected(group_index)
                    for c in self.constraints:
                        c.reject_move(
                            realIndexes=real_indexes,
                            relativeIndexes=relative_indexes
                        )
                else:
                    # 接受移动
                    self.increment_accepted()
                    self.groupSelector.move_accepted(group_index)
                    for c in self.constraints:
                        c.accept_move(
                            realIndexes=real_indexes,
                            relativeIndexes=relative_indexes
                        )
                
                # 递增尝试次数
                self.increment_tried()
                
        except Exception as e:
            LOGGER.error(f"Error in HRMC step: {str(e)}")
            raise
        
    def should_step_get_rejected(self, realIndexes, relativeIndexes):
        """
        实现MD-HRMC接受准则，直接使用PDF约束和能量约束的误差值。
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
            
        :Returns:
            #. reject (bool): 是否应该拒绝移动
        """
        try:
            print(f"[HRMC DEBUG] pdf_standardError_before={self.pdf_constraint.standardError}, after={self.pdf_constraint.afterMoveStandardError}")
            print(f"[HRMC DEBUG] energy_current={self.energy_constraint.data}, energy_last={self.energy_constraint._EnergyConstraint__last_accepted_energy}")
            # 获取移动的原子信息
            moved_elem = self.allElements[realIndexes[0]]
            is_hydrogen = (moved_elem == 'H')
            
            # 直接从约束中获取误差值
            pdf_old_error = self.pdf_constraint.standardError
            pdf_new_error = self.pdf_constraint.afterMoveStandardError
            delta_chi_sq = pdf_new_error - pdf_old_error
            
            # 获取能量差异
            current_energy = self.energy_constraint.data
            last_energy = self.energy_constraint._EnergyConstraint__last_accepted_energy
            delta_U = current_energy - last_energy

            # 计算当前T_chi（温度退火）
            T_chi = self.T_chi_0 * (self.cooling_rate**self.current_cycle)
            
            # 根据原子类型计算接受概率
            if is_hydrogen:
                # 对H原子: P = min[1, exp(-(ΔU/ω)/T_χ)]
                delta = (delta_U / self.omega) / T_chi
            else:
                # 对非H原子: P = min[1, exp(-(Δχ² + ΔU/ω)/T_χ)]
                delta = (delta_chi_sq + delta_U / self.omega) / T_chi

            # 处理数值溢出
            if delta > 700:
                self._write_data('HRMC_acceptance.txt',
                    f"{self._tried} 0.0 {delta_chi_sq:.6f} {delta_U:.6f} {T_chi:.6f}")
                self._write_data('HRMC_energy.txt',
                    f"{self._tried} {current_energy:.6f} {delta_U:.6f}")
                return True
                
            if delta < -700:
                self._write_data('HRMC_acceptance.txt', 
                    f"{self._tried} 1.0 {delta_chi_sq:.6f} {delta_U:.6f} {T_chi:.6f}")
                self._write_data('HRMC_energy.txt',
                    f"{self._tried} {current_energy:.6f} {delta_U:.6f}")
                return False
                            
            # 计算接受概率
            prob = min(1.0, np.exp(-delta))

            # 记录数据
            self._write_data('HRMC_acceptance.txt', 
                f"{self._tried} {prob:.6f} {delta_chi_sq:.6f} {delta_U:.6f} {T_chi:.6f}")
            self._write_data('HRMC_energy.txt',
                f"{self._tried} {current_energy:.6f} {delta_U:.6f}")

            # 做出接受决定
            accepted = np.random.random() <= prob
            if accepted:
                self.__last_energy = current_energy
                
            return not accepted

        except Exception as e:
            LOGGER.error(f"Error in HRMC acceptance: {str(e)}")
            return True

    def _run_md_step(self):
        """Run MD simulation step."""
        try:
            lmp = self.__lmp
            positions = self.realCoordinates
            elements = self.allElements

            # Suppress output
            lmp.command("echo none")  # Suppress LAMMPS output
            lmp.command("log none")  # Suppress LAMMPS log output

            # Important: Set atom creation style to silent before adding atoms
            lmp.command("thermo_style one")
            lmp.command("thermo_modify lost ignore flush no")

            # Reset system
            lmp.command("delete_atoms group all")
            
            # Update box
            box = self.boundaryConditions.get_vectors()
            box_lengths = np.array([np.linalg.norm(v) for v in box])
            lmp.command(f"change_box all x final 0 {box_lengths[0]} "
                       f"y final 0 {box_lengths[1]} "
                       f"z final 0 {box_lengths[2]} units box")

            # Create atoms with minimal output - create all atoms at once
            atom_commands = []
            for idx, (pos, elem) in enumerate(zip(positions, elements), start=1):
                atom_type = self.__element_types[elem]
                atom_commands.append(f"create_atoms {atom_type} single {pos[0]} {pos[1]} {pos[2]} units box")
            
            # Execute all create_atoms commands at once
            for cmd in atom_commands:
                lmp.command(cmd)
                
            # Set all charges at once using group all
            lmp.command("set group all charge 0.0")

            # Reset computes
            try:
                lmp.command("uncompute pe")
                lmp.command("uncompute temp")
            except:
                pass # Ignore if computes don't exist yet

            # Define computes before using them
            lmp.command("compute pe all pe")  # Define potential energy compute
            lmp.command("compute temp all temp")  # Define temperature compute

            # Setup QEq
            lmp.command("fix QEQ all qeq/reaxff 1 0.0 10.0 1e-6 reaxff maxiter 1000")

            # Apply minimize
            lmp.command("min_style cg")  # Use conjugate gradient
            lmp.command("min_modify dmax 0.1 line quadratic")            
            lmp.command("minimize 0.0e-4 1.0e-4 10 10") #0 0 is equivalent to not performing minimization.            
            
            # Setup timesetp
            lmp.command("timestep 0.1")

            # Only output final thermo values
            lmp.command("thermo_style custom step pe temp")
            lmp.command(f"thermo {self.md_steps}")
            lmp.command("log lammps.log append") 

            # Run NVT MD
            lmp.command(f"velocity all create {self.md_temperature} 4928459 rot yes dist gaussian")
            lmp.command(f"fix 1 all nve")
            lmp.command(f"fix NVT all temp/berendsen {self.md_temperature} {self.md_temperature} $(100*dt)")
            lmp.command(f"run {self.md_steps}")

            # Get final positions
            x = lmp.extract_atom("x", 3)
            final_positions = np.array([[x[i][0], x[i][1], x[i][2]] 
                                      for i in range(len(positions))])

            # Validate positions
            if not np.all(np.isfinite(final_positions)):
                raise ValueError("MD produced invalid positions (NaN or inf values)")
                
            # Apply PBC to ensure positions are within box
            box_lengths = np.array([np.linalg.norm(v) for v in self.boundaryConditions.get_vectors()])
            for i in range(3):
                final_positions[:, i] = np.mod(final_positions[:, i], box_lengths[i])

            # Update engine coordinates
            self.realCoordinates[:] = final_positions

            # Force full recomputation of constraints
            LOGGER.info("Recomputing constraints after MD...")
            self.pdf_constraint._runtime_initialize()
            self.energy_constraint._runtime_initialize()

            # These values are automatically available through extract_compute
            pe = lmp.extract_compute("pe", 0, 0)  # Potential energy
            temp = lmp.extract_compute("temp", 0, 0)  # Temperature

            # Write to MD_energy.txt
            self._write_data('MD_energy.txt', 
                f"{self.current_cycle} {pe:.6f} {temp:.6f}")

        except Exception as e:
            LOGGER.error(f"MD step failed: {str(e)}")
            raise

        finally:
            # Reset LAMMPS output
            lmp.command("log none")
            
    def run_md_hrmc(self):
        """Run complete MD-HRMC simulation."""
        try:
            for cycle in range(1, self.cycles_md_hrmc + 1):
                self.current_cycle = cycle
                LOGGER.info(f"\nStarting MD-HRMC cycle {cycle}/{self.cycles_md_hrmc}")

                # 运行HRMC
                LOGGER.info("运行HRMC模拟...")
                self.run(numberOfSteps=self.hrmc_steps)    

                # 每个循环后，计算并写入G(r)数据
                exp_data = self.pdf_constraint.experimentalData
                exp_r = exp_data[:, 0]
                
                # 计算当前G(r)
                hist = self.pdf_constraint.data
                rho0 = len(self.allElements) / self.volume
                sim_gr = self.pdf_constraint._PairDistributionConstraint__get_total_Gr(
                    data=hist,
                    rho0=rho0
                )

                # 运行MD
                LOGGER.info("运行MD模拟...")

                # 确保LAMMPS已初始化
                if self.__lmp is None:
                    self.__initialize_lammps()
                    
                self._run_md_step() 

                # Write G(r) data in same format as experimental data
                with open(os.path.join(self._output_dir, f'cycle_{cycle+1}_gr.txt'), 'w') as f:
                    for r, gr in zip(exp_r, sim_gr):
                        f.write(f"{r:.2f} {gr:.8e}\n")

                # Plot PDF comparison
                self.plot_pdf_comparison(
                    save_path=os.path.join(self._output_dir, f"cycle_{cycle}_pdf_comparison.png"),
                    show=False
                )

                # Export PDB after each cycle
                self.energy_constraint._export_pdb(
                    positions=self.realCoordinates,
                    elements=self.allElements,
                    box_lengths=np.array([np.linalg.norm(v) for v in self.boundaryConditions.get_vectors()]),
                    filename=os.path.join(self._output_dir, f"cycle_{cycle}.pdb")
                )

                # Save state
                tmp_lmp = self.__lmp  # Store LAMMPS object temporarily
                self.__lmp = None  # Clear for saving
                self.__initialized = False
                self.save()
                self.__lmp = tmp_lmp  # Restore LAMMPS object
                
        except Exception as e:
            LOGGER.error(f"Error in MD-HRMC simulation: {str(e)}")
            raise
            
    def plot_pdf_comparison(self, save_path=None, show=False):
        """
        绘制PDF对比图
        
        参数:
            save_path (str): 保存路径
            show (bool): 是否显示图形
        """
        if not hasattr(self, 'pdf_constraint') or self.pdf_constraint is None:
            LOGGER.error("未找到PDF约束，无法绘图")
            return
            
        try:
            # 获取约束
            pdf_constraint = self.pdf_constraint
            
            # 获取实验数据
            exp_data = pdf_constraint.experimentalData
            exp_r = exp_data[:, 0]
            exp_gr = exp_data[:, 1]
            
            # 计算模拟数据
            hist = pdf_constraint.data
            rho0 = len(self.allElements) / self.volume
            sim_gr = pdf_constraint._PairDistributionConstraint__get_total_Gr(
                data=hist,
                rho0=rho0
            )
            
            # 绘制实验数据和模型数据对比
            plt.figure(figsize=(10, 6))
            
            # 绘制实验数据
            plt.plot(exp_r, exp_gr, 'o', color='blue', markersize=3, label='Experimental')
            
            # 绘制模型数据
            plt.plot(exp_r, sim_gr, '-', color='red', linewidth=2, label='Simulated')
            
            # 添加图例和标签
            plt.legend(loc='best')
            plt.xlabel('r (Å)')
            plt.ylabel('g(r)')
            
            # 添加标题
            plt.title(f'g(r) Comparison - Cycle {self.current_cycle}')
            
            # 保存或显示图片
            if save_path:
                plt.savefig(save_path, dpi=300)
                LOGGER.info(f"PDF对比图已保存至: {save_path}")
            
            if show:
                plt.show()
            else:
                plt.close()
                
        except Exception as e:
            LOGGER.error(f"绘制PDF对比图时出错: {str(e)}")
