"""
HRMCEngine implements energy minimization using ReaxFF and HRMC optimization.
"""
from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from fullrmc.Engine import Engine 
from fullrmc.Globals import LOGGER
from fullrmc.Core.boundary_conditions_collection import transform_coordinates
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.DistanceConstraints import InterMolecularDistanceConstraint
from fullrmc.Constraints.EnergyConstraints import EnergyConstraint

# Element masses (in g/mol)
ELEMENT_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.380
}    

class HRMCEngine(Engine):
    def __init__(self, path, freshStart=False):
        super(HRMCEngine, self).__init__(path=path, freshStart=freshStart)
            
        # Constants - Use proper units for LAMMPS real units
        self.KB = 0.001987204258  # Boltzmann constant in kcal/mol/K
        
        # HRMC parameters
        self.hrmc_initial_temp = None
        self.hrmc_final_temp = None 
        self.hrmc_steps = None
        self.cycles_hrmc = None
        self.current_cycle = 1

        # Counters for statistics
        self._tried = 1  # Start at 1
        self._accepted = 0

        # Store constraints
        self.pdf_constraint = None
        self.dis_constraint = None
        self.energy_constraint = None

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
                'HRMC_energy.txt': "# Step Energy DeltaU"
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
        state = super(HRMCEngine, self).__getstate__()
        state['_HRMCEngine__last_energy'] = self.__last_energy
        return state
        
    def __setstate__(self, state):
        """Handle deserialization."""
        self.__dict__.update(state)

    def cleanup(self):
        """Clean up resources."""
        pass

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
        super(HRMCEngine, self).set_pdb(pdb, *args, **kwargs)
        
        # Skip omega calculation if no structure loaded yet
        if pdb is None:
            return
            
        # Calculate atom counts and omega
        elements = self.allElements
        self.N_total = len(elements)
        self.N_nonH = sum(1 for elem in elements if elem != 'H')
        
        if self.N_nonH == 0:
            raise ValueError("Structure must contain non-hydrogen atoms")
            
        # Calculate omega (parameter used for energy weighting)
        self.omega = (self.N_total**2 / self.N_nonH) * 7 * 23.0609
        
        LOGGER.info(f"Structure loaded with {self.N_total} total atoms ({self.N_nonH} non-hydrogen)")
        LOGGER.info(f"Calculated omega: {self.omega}")

    def add_constraints(self, constraints):
        """Add constraints with proper initialization."""
        try:
            for constraint in constraints:
                if isinstance(constraint, PairDistributionConstraint):
                    self.pdf_constraint = constraint
                elif isinstance(constraint, InterMolecularDistanceConstraint):
                    self.dis_constraint = constraint
                elif isinstance(constraint, EnergyConstraint):
                    self.energy_constraint = constraint
            
            # Call parent's add_constraints AFTER storing references
            super(HRMCEngine, self).add_constraints(constraints)
            
        except Exception as e:
            LOGGER.error(f"Error adding constraints: {str(e)}")
            raise
                
    def initialize_used_constraints(self):
        """Initialize constraints safely."""
        try:
            # 先初始化所有约束数据
            super(HRMCEngine, self).initialize_used_constraints()
        
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

    def setup_hrmc(self, initial_temp, final_temp, steps, cycles):
        """Setup HRMC parameters."""
        self.hrmc_initial_temp = initial_temp
        self.hrmc_final_temp = final_temp
        self.hrmc_steps = steps
        self.cycles_hrmc = cycles
        
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

        # 🔧 KOKKOS性能监控（可选）
        kokkos_energy_times = []
        kokkos_error_count = 0

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
                
                # 安全计算约束 - 捕获能量约束错误
                try:
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
                
                    # 正常检查是否应该接受移动
                    reject_move = self.should_step_get_rejected(real_indexes, relative_indexes)
                
                except Exception as e:
                    error_msg = str(e)
                    # 🔧 增强的KOKKOS错误处理
                    if "bondchk failed" in error_msg or "malformed" in error_msg:
                        LOGGER.warn(f"ReaxFF结构错误 - 原子距离过近或重叠: {error_msg}")
                        reject_move = True
                        self.realCoordinates[relative_indexes] = old_positions
                    elif "kokkos" in error_msg.lower():
                        LOGGER.warn(f"KOKKOS计算错误: {error_msg}")
                        reject_move = True
                        self.realCoordinates[relative_indexes] = old_positions
                    elif "reaxff/kk" in error_msg.lower():
                        LOGGER.warn(f"KOKKOS ReaxFF错误: {error_msg}")
                        reject_move = True
                        self.realCoordinates[relative_indexes] = old_positions
                    else:
                        # 其他未知错误
                        LOGGER.error(f"计算约束时发生未知错误: {error_msg}")
                        reject_move = True
                        self.realCoordinates[relative_indexes] = old_positions
                        reject_move = True
            
                if reject_move:
                    # 拒绝移动
                    self.realCoordinates[relative_indexes] = old_positions
                    self.groupSelector.move_rejected(group_index)
                    try:
                        # 尝试调用约束的reject_move方法，但忽略可能的错误
                        for c in self.constraints:
                            try:
                                c.reject_move(
                                    realIndexes=real_indexes,
                                    relativeIndexes=relative_indexes
                                )
                            except:
                                pass  # 忽略错误，确保继续处理
                    except:
                        pass  # 忽略任何异常
                else:
                    # 接受移动
                    self.increment_accepted()
                    self.groupSelector.move_accepted(group_index)
                    try:
                        for c in self.constraints:
                            c.accept_move(
                                realIndexes=real_indexes,
                                relativeIndexes=relative_indexes
                            )
                    except Exception as e:
                        LOGGER.error(f"接受移动时出错: {str(e)}")
                        # 如果accept_move失败，仍然保持新坐标，但记录错误
                
                # 递增尝试次数
                self.increment_tried()
                
        except Exception as e:
            LOGGER.error(f"Error in HRMC step: {str(e)}")
            raise
        
    def should_step_get_rejected(self, realIndexes, relativeIndexes):
        """
        实现HRMC接受准则，直接使用能量约束的误差值。
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
            
        :Returns:
            #. reject (bool): 是否应该拒绝移动
        """
        try:
            print(f"[HRMC DEBUG] energy_current={self.energy_constraint.data}, energy_last={self.energy_constraint._EnergyConstraint__last_accepted_energy}")

            # 第一阶段：检查分子间距离约束（刚性约束）
            if self.dis_constraint:
                # 检查距离约束是否被违反 - 使用activeAtomsData来判断
                if (self.dis_constraint.activeAtomsDataAfterMove is not None and 
                    self.dis_constraint.activeAtomsDataBeforeMove is not None):
            
                    # 检查数量是否增加 - 表示有更多原子违反了约束
                    after_number = self.dis_constraint.activeAtomsDataAfterMove.get("number", 0)
                    before_number = self.dis_constraint.activeAtomsDataBeforeMove.get("number", 0)
            
                    # 如果违反约束的数量增加，拒绝移动
                    if np.any(after_number > before_number):
                        # 记录违反刚性约束的情况
                        moved_elem = self.allElements[realIndexes[0]]
                        step = getattr(self, 'tried', 0)
                        self._write_data('HRMC_acceptance.txt', 
                            f"{step}    {moved_elem}    REJECTED_RIGID_CONSTRAINT")
                
                        # 违反了刚性约束，立即拒绝
                        return True

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
            T_chi = self.T_chi_0 * (self.cooling_rate**self.current_cycle) * 10

            delta = (0.001 * delta_chi_sq + delta_U / self.omega) / T_chi
            
            # 处理数值溢出
            if delta > 700:
                self._write_data('HRMC_acceptance.txt',
                    f"{self._tried} 0.0 {delta_chi_sq:.6f} {delta_U:.6f} {T_chi:.6f}")
                self._write_data('HRMC_energy.txt',
                    f"{self._tried} {current_energy:.6f} {delta_U:.6f}")
                return True
                
            if delta < -700:
                self._write_data('HRMC_acceptance.txt', 
                    f"{self._tried} 1.0 {delta_U:.6f} {T_chi:.6f}")
                self._write_data('HRMC_energy.txt',
                    f"{self._tried} {current_energy:.6f} {delta_U:.6f}")
                return False
                            
            # 计算接受概率
            prob = min(1.0, np.exp(-delta))

            # 记录数据
            self._write_data('HRMC_acceptance.txt', 
                f"{self._tried} {prob:.6f} {delta:.6f}  {delta_chi_sq:.6f}  {delta_U:.6f}  {T_chi:.6f}")
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

    def run_hrmc(self):
        """Run complete HRMC simulation."""
        try:
            for cycle in range(1, self.cycles_hrmc + 1):
                self.current_cycle = cycle
                LOGGER.info(f"\nStarting HRMC cycle {cycle}/{self.cycles_hrmc}")

                # 运行HRMC
                LOGGER.info("运行HRMC模拟...")
                self.run(numberOfSteps=self.hrmc_steps)    

                # 每个循环后，导出结构
                self.export_pdb(
                    path=os.path.join(self._output_dir, f"cycle_{cycle}.pdb")
                )
                # 导出周期结束PDF对比图
                # 导出周期结束PDF对比图
                plot_pdf_path = os.path.join(self._output_dir, f"pdf_cycle_{cycle}.png")
                self.plot_pdf_comparison(save_path=plot_pdf_path, show=False)
                # ===== 添加绘图 =====
                plot_energy_path = os.path.join(self._output_dir, f"energy_cycle_{cycle}.png")
                self.plot_energy_profile(save_path=plot_energy_path, show=False)

                self.save()
                
        except Exception as e:
            LOGGER.error(f"Error in HRMC simulation: {str(e)}")
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
            
    def plot_energy_profile(self, save_path=None, show=False):
        """
        绘制能量变化曲线
        
        参数:
            save_path (str): 保存路径
            show (bool): 是否显示图形
        """
        try:
            # 读取能量数据
            energy_file = os.path.join(self._output_dir, 'HRMC_energy.txt')
            if not os.path.exists(energy_file):
                LOGGER.error(f"能量文件不存在: {energy_file}")
                return
                
            data = np.loadtxt(energy_file, skiprows=1)
            
            if data.size == 0:
                LOGGER.error("能量数据为空")
                return
                
            # 提取步骤和能量
            steps = data[:, 0]
            energies = data[:, 1]
            
            # 绘制能量变化
            plt.figure(figsize=(10, 6))
            plt.plot(steps, energies, '-', color='blue', linewidth=1.5)
            
            # 添加标签和标题
            plt.xlabel('step')
            plt.ylabel('energy (kcal/mol)')
            plt.title('HRMC energy')
            
            # 添加网格
            plt.grid(True, linestyle='--', alpha=0.7)
            
            # 保存或显示图片
            if save_path:
                plt.savefig(save_path, dpi=300)
                LOGGER.info(f"能量变化图已保存至: {save_path}")
            
            if show:
                plt.show()
            else:
                plt.close()
                
        except Exception as e:
            LOGGER.error(f"绘制能量变化图时出错: {str(e)}")
