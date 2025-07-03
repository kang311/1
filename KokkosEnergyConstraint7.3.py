"""
KokkosEnergyConstraint - KOKKOS包装的ReaxFF能量约束实现
提供动态内存管理、错误恢复和并行计算优化
"""

# Standard libraries imports
from __future__ import print_function
from collections import OrderedDict
import numpy as np
import os
import gc
import time
from contextlib import contextmanager

# Fullrmc imports
from fullrmc.Core.Constraint import Constraint
from fullrmc.Core.Collection import is_number, reset_if_collected_out_of_date, get_path
from fullrmc.Globals import LOGGER, FLOAT_TYPE
from fullrmc.Core.boundary_conditions_collection import transform_coordinates

# Custom exceptions
class LAMMPSError(Exception):
    """LAMMPS-specific error."""
    pass

class KokkosError(Exception):
    """KOKKOS-specific error."""
    pass

class BoxError(Exception):
    """Box definition error."""
    pass

class EnergyComputationError(Exception):
    """Energy computation specific error."""
    pass

# Try to import LAMMPS with KOKKOS support
try:
    from lammps import lammps
    LAMMPS_IMPORT = True
    LOGGER.info("LAMMPS Python interface loaded successfully")
except ImportError:
    LAMMPS_IMPORT = False
    LOGGER.error("LAMMPS Python interface not found")

# Element masses (in g/mol)
ELEMENT_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078
}

# 元素间最小距离定义 (Angstrom)
ELEMENT_MIN_DISTANCES = {
    ('H', 'H'): 0.8,   ('H', 'C'): 1.0,   ('H', 'N'): 0.9,   ('H', 'O'): 0.9,   ('H', 'S'): 1.2,
    ('C', 'C'): 1.2,   ('C', 'N'): 1.1,   ('C', 'O'): 1.1,   ('C', 'S'): 1.5,
    ('N', 'N'): 1.1,   ('N', 'O'): 1.1,   ('N', 'S'): 1.5,
    ('O', 'O'): 1.2,   ('O', 'S'): 1.4,
    ('S', 'S'): 1.8,
}

class KokkosMemoryManager:
    """KOKKOS内存管理器 - 提供动态内存分配和释放"""
    
    def __init__(self):
        self.allocated_buffers = {}
        self.buffer_sizes = {}
        self.allocation_count = 0
        self.deallocation_count = 0
        self.peak_memory_usage = 0
        self.current_memory_usage = 0
        
    def allocate_buffer(self, name, size, dtype=np.float64):
        """动态分配内存缓冲区"""
        try:
            # 如果缓冲区已存在且大小合适，直接复用
            if name in self.allocated_buffers and self.buffer_sizes[name] >= size:
                LOGGER.debug(f"Reusing existing buffer {name} (size: {self.buffer_sizes[name]})")
                return self.allocated_buffers[name]
            
            # 释放旧缓冲区
            if name in self.allocated_buffers:
                self.deallocate_buffer(name)
            
            # 分配新缓冲区
            if isinstance(size, tuple):
                buffer = np.zeros(size, dtype=dtype)
                buffer_size = np.prod(size) * np.dtype(dtype).itemsize
            else:
                buffer = np.zeros(size, dtype=dtype)
                buffer_size = size * np.dtype(dtype).itemsize
                
            self.allocated_buffers[name] = buffer
            self.buffer_sizes[name] = buffer_size
            self.current_memory_usage += buffer_size
            self.peak_memory_usage = max(self.peak_memory_usage, self.current_memory_usage)
            self.allocation_count += 1
            
            LOGGER.debug(f"Allocated buffer {name}: {buffer_size} bytes")
            return buffer
            
        except Exception as e:
            LOGGER.error(f"Failed to allocate buffer {name}: {str(e)}")
            raise KokkosError(f"Memory allocation failed: {str(e)}")
    
    def deallocate_buffer(self, name):
        """释放指定的内存缓冲区"""
        if name in self.allocated_buffers:
            buffer_size = self.buffer_sizes[name]
            del self.allocated_buffers[name]
            del self.buffer_sizes[name]
            self.current_memory_usage -= buffer_size
            self.deallocation_count += 1
            LOGGER.debug(f"Deallocated buffer {name}: {buffer_size} bytes")
            
    def cleanup_all(self):
        """清理所有分配的内存"""
        for name in list(self.allocated_buffers.keys()):
            self.deallocate_buffer(name)
        gc.collect()
        LOGGER.info(f"Memory cleanup complete. Allocations: {self.allocation_count}, "
                   f"Deallocations: {self.deallocation_count}")
    
    def get_memory_stats(self):
        """获取内存使用统计"""
        return {
            'current_usage': self.current_memory_usage,
            'peak_usage': self.peak_memory_usage,
            'allocation_count': self.allocation_count,
            'deallocation_count': self.deallocation_count,
            'active_buffers': len(self.allocated_buffers)
        }

class KokkosDistanceChecker:
    """KOKKOS并行化的原子距离检查器"""
    
    def __init__(self, memory_manager):
        self.memory_manager = memory_manager
        
    def check_atomic_distances(self, positions, elements, min_distance_threshold=0.8):
        """
        并行检查原子间距离
        
        Parameters:
            positions: 原子坐标 (N, 3)
            elements: 原子元素类型列表
            min_distance_threshold: 默认最小距离阈值
            
        Returns:
            violations: 违反距离要求的原子对信息列表
            is_valid: 是否所有距离都合法
        """
        n_atoms = len(positions)
        violations = []
        
        # 分配距离矩阵缓冲区
        distance_matrix = self.memory_manager.allocate_buffer(
            'distance_matrix', (n_atoms, n_atoms), np.float64
        )
        
        try:
            # 计算所有原子对的距离 - 可以用KOKKOS并行化
            self._compute_distance_matrix_parallel(positions, distance_matrix)
            
            # 检查违反情况
            for i in range(n_atoms):
                for j in range(i+1, n_atoms):
                    distance = distance_matrix[i, j]
                    elem_i, elem_j = elements[i], elements[j]
                    pair_key = tuple(sorted([elem_i, elem_j]))
                    
                    min_required = ELEMENT_MIN_DISTANCES.get(pair_key, min_distance_threshold)
                    
                    if distance < min_required:
                        violation_severity = (min_required - distance) / min_required
                        violations.append({
                            'atom_i': i,
                            'atom_j': j,
                            'elements': (elem_i, elem_j),
                            'distance': distance,
                            'min_required': min_required,
                            'violation_magnitude': min_required - distance,
                            'severity': violation_severity
                        })
            
            # 按严重程度排序
            violations.sort(key=lambda x: x['severity'], reverse=True)
            
            return violations, len(violations) == 0
            
        except Exception as e:
            LOGGER.error(f"Distance checking failed: {str(e)}")
            return [], False
    
    def _compute_distance_matrix_parallel(self, positions, distance_matrix):
        """并行计算距离矩阵 - 模拟KOKKOS并行计算"""
        n_atoms = positions.shape[0]
        
        # 在真实的KOKKOS实现中，这里会使用parallel_for
        for i in range(n_atoms):
            for j in range(i, n_atoms):
                if i == j:
                    distance_matrix[i, j] = 0.0
                else:
                    distance = np.linalg.norm(positions[i] - positions[j])
                    distance_matrix[i, j] = distance
                    distance_matrix[j, i] = distance
    
    def adjust_violating_positions(self, positions, violations, adjustment_factor=1.2):
        """调整违反距离要求的原子位置"""
        if not violations:
            return positions
            
        adjusted_positions = positions.copy()
        adjustment_buffer = self.memory_manager.allocate_buffer(
            'position_adjustments', positions.shape, np.float64
        )
        adjustment_buffer.fill(0.0)
        
        # 按严重程度处理违反情况
        for violation in violations[:10]:  # 只处理最严重的10个违反
            i, j = violation['atom_i'], violation['atom_j']
            required_distance = violation['min_required'] * adjustment_factor
            
            # 计算调整向量
            vector = adjusted_positions[j] - adjusted_positions[i]
            current_distance = np.linalg.norm(vector)
            
            if current_distance > 0:
                unit_vector = vector / current_distance
                adjustment_distance = required_distance - current_distance
                
                if adjustment_distance > 0:
                    adjustment = unit_vector * adjustment_distance * 0.5
                    adjustment_buffer[i] -= adjustment
                    adjustment_buffer[j] += adjustment
        
        # 应用调整
        adjusted_positions += adjustment_buffer
        
        return adjusted_positions

class KokkosReaxFFWrapper:
    """KOKKOS包装的ReaxFF计算器"""
    
    def __init__(self, memory_manager, reaxff_params, control_params):
        self.memory_manager = memory_manager
        self.reaxff_params = reaxff_params
        self.control_params = control_params
        self.distance_checker = KokkosDistanceChecker(memory_manager)
        self.computation_cache = {}
        self.error_recovery_attempts = 0
        self.max_recovery_attempts = 3
        
    @contextmanager
    def kokkos_computation_context(self):
        """KOKKOS计算上下文管理器"""
        start_time = time.time()
        initial_memory = self.memory_manager.current_memory_usage
        
        try:
            LOGGER.debug("Entering KOKKOS computation context")
            yield
        except Exception as e:
            LOGGER.error(f"Error in KOKKOS computation context: {str(e)}")
            raise
        finally:
            computation_time = time.time() - start_time
            final_memory = self.memory_manager.current_memory_usage
            memory_delta = final_memory - initial_memory
            LOGGER.debug(f"KOKKOS computation completed in {computation_time:.4f}s, "
                        f"memory delta: {memory_delta} bytes")
    
    def safe_energy_computation(self, lmp, positions, elements, box_lengths):
        """安全的ReaxFF能量计算，包含完整的错误恢复"""
        
        with self.kokkos_computation_context():
            # 第一阶段：距离预检查和调整
            violations, is_valid = self.distance_checker.check_atomic_distances(positions, elements)
            
            if not is_valid:
                LOGGER.warn(f"发现 {len(violations)} 个原子距离违反，尝试调整")
                positions = self.distance_checker.adjust_violating_positions(positions, violations)
                
                # 重新检查
                violations, is_valid = self.distance_checker.check_atomic_distances(positions, elements)
                if not is_valid:
                    LOGGER.error(f"调整后仍有 {len(violations)} 个违反，使用渐进式计算")
                    return self._gradual_energy_computation(lmp, positions, elements, box_lengths)
            
            # 第二阶段：尝试标准能量计算
            try:
                return self._standard_energy_computation(lmp, positions, elements, box_lengths)
            except Exception as e:
                if "bondchk failed" in str(e) or "malformed" in str(e):
                    LOGGER.warn(f"标准计算失败: {str(e)}，尝试恢复策略")
                    return self._error_recovery_computation(lmp, positions, elements, box_lengths, str(e))
                else:
                    raise
    
    def _standard_energy_computation(self, lmp, positions, elements, box_lengths):
        """标准的ReaxFF能量计算"""
        n_atoms = len(positions)
        
        # 分配计算缓冲区
        position_buffer = self.memory_manager.allocate_buffer('positions', (n_atoms, 3))
        charge_buffer = self.memory_manager.allocate_buffer('charges', n_atoms)
        
        # 复制数据到缓冲区
        position_buffer[:] = positions
        
        # 保存旧电荷
        try:
            old_charges = lmp.extract_atom("q")
            if old_charges is not None:
                old_charges_array = np.array(old_charges)
                old_charges_array = np.nan_to_num(old_charges_array, nan=0.0, posinf=0.0, neginf=0.0)
                charge_buffer[:len(old_charges_array)] = old_charges_array
        except:
            charge_buffer.fill(0.0)
        
        # 清理现有fix
        self._cleanup_lammps_fixes(lmp)
        
        # 重建原子系统
        lmp.command("delete_atoms group all")
        lmp.command(f"change_box all x final 0 {box_lengths[0]} y final 0 {box_lengths[1]} z final 0 {box_lengths[2]} units box")
        
        # 坐标包装
        wrapped_positions = self._wrap_positions(position_buffer, box_lengths)
        
        # 创建原子
        self._create_atoms_safe(lmp, wrapped_positions, elements, charge_buffer)
        
        # 验证原子数量
        natoms = lmp.get_natoms()
        if natoms != n_atoms:
            raise LAMMPSError(f"Atom count mismatch: {natoms} != {n_atoms}")
        
        # 配置KOKKOS优化的计算设置
        self._setup_kokkos_computation(lmp)
        
        # QEQ计算
        self._perform_qeq_calculation(lmp)
        
        # 能量最小化
        self._perform_minimization(lmp)
        
        # 提取结果
        return self._extract_computation_results(lmp, n_atoms)
    
    def _gradual_energy_computation(self, lmp, positions, elements, box_lengths, steps=5):
        """渐进式能量计算，避免突然的构型变化"""
        LOGGER.info(f"执行渐进式能量计算，{steps}步")
        
        # 生成安全的目标位置
        violations, _ = self.distance_checker.check_atomic_distances(positions, elements)
        safe_positions = self.distance_checker.adjust_violating_positions(positions, violations, 1.5)
        
        # 分配插值缓冲区
        n_atoms = len(positions)
        interpolated_buffer = self.memory_manager.allocate_buffer('interpolated_positions', (n_atoms, 3))
        
        current_positions = positions.copy()
        
        for step in range(steps):
            alpha = (step + 1) / steps
            
            # 线性插值到目标位置
            interpolated_buffer[:] = (1 - alpha) * current_positions + alpha * safe_positions
            
            try:
                energy, final_positions = self._standard_energy_computation(
                    lmp, interpolated_buffer, elements, box_lengths
                )
                
                LOGGER.debug(f"渐进步骤 {step+1}/{steps}: 能量 = {energy:.6f}")
                
                # 更新当前位置为计算结果
                current_positions = final_positions
                
                if step == steps - 1:
                    return energy, final_positions
                    
            except Exception as e:
                if step == 0:
                    # 第一步就失败，尝试更保守的方法
                    LOGGER.warn("第一步渐进计算失败，使用保守模式")
                    return self._conservative_energy_computation(lmp, positions, elements, box_lengths)
                else:
                    # 返回上一步的结果
                    LOGGER.warn(f"渐进计算在步骤 {step+1} 失败，返回步骤 {step} 的结果")
                    return self._standard_energy_computation(lmp, current_positions, elements, box_lengths)
        
        raise EnergyComputationError("渐进式计算未能完成")
    
    def _conservative_energy_computation(self, lmp, positions, elements, box_lengths):
        """保守的能量计算，使用最小化的参数"""
        LOGGER.info("执行保守模式能量计算")
        
        try:
            # 使用标准流程但更保守的参数
            energy, final_positions = self._standard_energy_computation(lmp, positions, elements, box_lengths)
            
            # 替换最小化设置为更保守的参数
            lmp.command("min_style cg")
            lmp.command("min_modify dmax 0.01 line backtrack")  # 更小的步长
            lmp.command("minimize 1.0e-6 1.0e-8 100 1000")     # 更少的步数，更严格的收敛
            
            # 重新提取结果
            return self._extract_computation_results(lmp, len(positions))
            
        except Exception as e:
            LOGGER.error(f"保守计算也失败: {str(e)}")
            # 最后的尝试：只计算单点能量，不做优化
            return self._single_point_energy_computation(lmp, positions, elements, box_lengths)
    
    def _single_point_energy_computation(self, lmp, positions, elements, box_lengths):
        """单点能量计算，不进行几何优化"""
        LOGGER.info("执行单点能量计算")
        
        n_atoms = len(positions)
        
        # 重建基本系统
        lmp.command("delete_atoms group all")
        lmp.command(f"change_box all x final 0 {box_lengths[0]} y final 0 {box_lengths[1]} z final 0 {box_lengths[2]} units box")
        
        # 创建原子（不进行距离调整）
        self._create_atoms_basic(lmp, positions, elements)
        
        # 只执行QEQ，不做几何优化
        try:
            lmp.command("fix QEQ all qeq/reaxff 1 0.0 10.0 1e-6 reaxff maxiter 200")
            lmp.command("run 0")
            lmp.command("unfix QEQ")
        except:
            # QEQ也失败的话，直接计算能量
            pass
        
        # 提取能量（不做最小化）
        try:
            energy = lmp.extract_compute("thermo_pe", 0, 0)
            return float(energy), positions  # 返回原始坐标
        except:
            # 如果连能量都无法提取，返回高惩罚值
            LOGGER.error("无法提取能量，返回惩罚值")
            return 1e6, positions
    
    def _error_recovery_computation(self, lmp, positions, elements, box_lengths, error_msg):
        """错误恢复计算策略"""
        self.error_recovery_attempts += 1
        
        if self.error_recovery_attempts > self.max_recovery_attempts:
            LOGGER.error("错误恢复尝试次数超限，返回惩罚能量")
            return 1e6, positions
        
        LOGGER.warn(f"执行错误恢复 (尝试 {self.error_recovery_attempts}/{self.max_recovery_attempts})")
        
        # 根据错误类型选择恢复策略
        if "bondchk failed" in error_msg:
            # bondchk错误 - 距离问题
            return self._handle_bondchk_error(lmp, positions, elements, box_lengths)
        elif "malformed" in error_msg:
            # 构型问题
            return self._handle_malformed_error(lmp, positions, elements, box_lengths)
        else:
            # 其他错误
            return self._handle_generic_error(lmp, positions, elements, box_lengths)
    
    def _handle_bondchk_error(self, lmp, positions, elements, box_lengths):
        """处理bondchk失败错误"""
        LOGGER.debug("处理bondchk错误")
        
        # 更激进地调整原子位置
        violations, _ = self.distance_checker.check_atomic_distances(positions, elements, 1.2)
        adjusted_positions = self.distance_checker.adjust_violating_positions(positions, violations, 1.8)
        
        # 尝试渐进式计算
        return self._gradual_energy_computation(lmp, adjusted_positions, elements, box_lengths, steps=3)
    
    def _handle_malformed_error(self, lmp, positions, elements, box_lengths):
        """处理构型错误"""
        LOGGER.debug("处理构型错误")
        
        # 对整个系统进行轻微的随机扰动
        perturbation = np.random.normal(0, 0.1, positions.shape)
        perturbed_positions = positions + perturbation
        
        return self._conservative_energy_computation(lmp, perturbed_positions, elements, box_lengths)
    
    def _handle_generic_error(self, lmp, positions, elements, box_lengths):
        """处理一般错误"""
        LOGGER.debug("处理一般计算错误")
        return self._single_point_energy_computation(lmp, positions, elements, box_lengths)
    
    def _cleanup_lammps_fixes(self, lmp):
        """清理LAMMPS中的fix"""
        fixes_to_remove = ["QEQ", "halt_errors", "halt_bondchk"]
        for fix_name in fixes_to_remove:
            try:
                lmp.command(f"unfix {fix_name}")
            except:
                pass
    
    def _wrap_positions(self, positions, box_lengths):
        """坐标包装到模拟盒子内"""
        wrapped = positions.copy()
        wrapped[:, 0] = positions[:, 0] % box_lengths[0]
        wrapped[:, 1] = positions[:, 1] % box_lengths[1]
        wrapped[:, 2] = positions[:, 2] % box_lengths[2]
        return wrapped
    
    def _create_atoms_safe(self, lmp, positions, elements, charges):
        """安全地创建原子"""
        element_types = {elem: idx+1 for idx, elem in enumerate(sorted(set(elements)))}
        
        for idx, (pos, elem) in enumerate(zip(positions, elements), start=1):
            atom_type = element_types[elem]
            if not np.all(np.isfinite(pos)):
                raise LAMMPSError(f"Invalid position for atom {idx}: {pos}")
            
            lmp.command(f"create_atoms {atom_type} single {pos[0]} {pos[1]} {pos[2]}")
            charge = charges[idx-1] if len(charges) >= idx else 0.0
            lmp.command(f"set atom {idx} charge {charge}")
    
    def _create_atoms_basic(self, lmp, positions, elements):
        """基本的原子创建（用于单点能量计算）"""
        element_types = {elem: idx+1 for idx, elem in enumerate(sorted(set(elements)))}
        
        for idx, (pos, elem) in enumerate(zip(positions, elements), start=1):
            atom_type = element_types[elem]
            lmp.command(f"create_atoms {atom_type} single {pos[0]} {pos[1]} {pos[2]}")
            lmp.command(f"set atom {idx} charge 0.0")
    
    def _setup_kokkos_computation(self, lmp):
        """设置KOKKOS优化的计算参数"""
        # 设置邻居列表
        lmp.command("neighbor 3.0 bin")
        lmp.command("neigh_modify every 1 delay 0 check yes")
        
        # 设置错误处理
        lmp.command("fix halt_errors all halt 1 error continue")
    
    def _perform_qeq_calculation(self, lmp):
        """执行QEQ电荷平衡计算"""
        try:
            lmp.command("suffix off")  # 对QEQ禁用加速
            lmp.command("fix QEQ all qeq/reaxff 1 0.0 10.0 1e-6 reaxff maxiter 1000")
            lmp.command("run 0")
            lmp.command("unfix QEQ")
        except Exception as e:
            LOGGER.warn(f"QEQ计算失败: {str(e)}，继续进行能量计算")
    
    def _perform_minimization(self, lmp):
        """执行能量最小化"""
        try:
            lmp.command("suffix omp")  # 启用OpenMP加速
            lmp.command("min_style cg")
            lmp.command("min_modify dmax 0.1 line quadratic")
            lmp.command("minimize 0.0e-4 1.0e-4 0 0")  # 0 0 表示只计算能量不优化几何
        except Exception as e:
            LOGGER.warn(f"最小化失败: {str(e)}，使用当前构型")
    
    def _extract_computation_results(self, lmp, n_atoms):
        """提取计算结果"""
        try:
            # 提取最终坐标
            x = lmp.extract_atom("x", 3)
            final_positions = np.array([[x[i][0], x[i][1], x[i][2]] for i in range(n_atoms)])
            
            # 提取能量
            energy = lmp.extract_compute("thermo_pe", 0, 0)
            
            return float(energy), final_positions
            
        except Exception as e:
            LOGGER.error(f"结果提取失败: {str(e)}")
            raise EnergyComputationError(f"无法提取计算结果: {str(e)}")

class KokkosEnergyConstraint(Constraint):
    """
    KOKKOS增强的ReaxFF能量约束实现
    
    提供动态内存管理、错误恢复和并行计算优化
    """
    
    def __init__(self, reaxff_params="ffield.reax", control_params="control.reax", 
                 compute_forces=False, weight=1.0, use_kokkos=True):
        # Initialize base class
        super(KokkosEnergyConstraint, self).__init__()

        # Check LAMMPS availability
        if not LAMMPS_IMPORT:
            raise ImportError(LOGGER.error("LAMMPS Python interface not found"))

        # Initialize core components
        self.__lmp = None
        self.__initialized = False
        self.__element_types = {}
        self.__energy_cache = {}
        self.__cached_box = None
        self.__last_accepted_energy = None
        
        # KOKKOS specific attributes
        self.use_kokkos = use_kokkos
        self.memory_manager = KokkosMemoryManager()
        self.reaxff_wrapper = None
        
        # Parameters
        self.__reaxff_params = None
        self.__control_params = None
        self.__compute_forces = None
        
        # Set parameters
        self.set_reaxff_params(reaxff_params)
        self.set_control_params(control_params)
        self.set_compute_forces(compute_forces)

        # Initialize flags for move computation
        self._constraint_before_move_reset = False
        self._constraint_after_move_reset = False

        # Setup frame data
        frame_data = list(self.FRAME_DATA)
        frame_data.extend([
            '_KokkosEnergyConstraint__element_types',
            '_KokkosEnergyConstraint__reaxff_params',
            '_KokkosEnergyConstraint__control_params',
        ])

        runtime_data = list(self.RUNTIME_DATA)

        # Set the class attributes
        object.__setattr__(self, 'FRAME_DATA', tuple(frame_data))
        object.__setattr__(self, 'RUNTIME_DATA', tuple(runtime_data))

    def set_reaxff_params(self, reaxff_params):
        """Set ReaxFF parameter file path."""
        if not os.path.exists(reaxff_params):
            raise ValueError(f"ReaxFF parameter file not found: {reaxff_params}")
        self.__reaxff_params = reaxff_params
        LOGGER.info(f"ReaxFF parameters set to: {reaxff_params}")

    def set_control_params(self, control_params):
        """Set ReaxFF control parameter file path."""
        if not os.path.exists(control_params):
            raise ValueError(f"Control parameter file not found: {control_params}")
        self.__control_params = control_params
        LOGGER.info(f"Control parameters set to: {control_params}")

    def set_compute_forces(self, compute_forces):
        """Set whether to compute atomic forces."""
        self.__compute_forces = bool(compute_forces)

    def _get_box_vectors(self):
        """Get simulation box vectors from engine."""
        try:
            if self.engine.isPBC:
                box_vectors = self.engine.boundaryConditions.get_vectors()
                if box_vectors is None:
                    raise BoxError("No boundary conditions defined in engine")
                
                box_lengths = np.array([np.linalg.norm(v) for v in box_vectors])
                
                if np.any(box_lengths <= 0):
                    raise BoxError("Invalid box dimensions")
                    
                return box_lengths
            else:
                coords = self.engine.realCoordinates
                min_coords = np.min(coords, axis=0) 
                max_coords = np.max(coords, axis=0)
                box_lengths = (max_coords - min_coords) * 1.1
                return box_lengths
                
        except Exception as e:
            LOGGER.error(f"Box vector calculation failed: {str(e)}")
            raise BoxError(f"Failed to get box vectors: {str(e)}")

    def _initialize_lammps_with_kokkos(self):
        """使用KOKKOS支持初始化LAMMPS"""
        try:
            if self.__initialized and self.__lmp is not None:
                return

            # Create LAMMPS instance with KOKKOS support
            logfile = os.path.join(self.engine.path, "lammps_kokkos.log")
            
            if self.use_kokkos:
                # KOKKOS enabled LAMMPS
                cmdargs = [
                    "-log", "lammps.log", 
                    "-screen", logfile,
                    "-kokkos", "on",
                    "-pk", "kokkos", "newton", "on", "neigh", "half"
                ]
            else:
                # Standard LAMMPS
                cmdargs = ["-log", "lammps.log", "-screen", logfile]
                
            self.__lmp = lammps(cmdargs=cmdargs)
            lmp = self.__lmp

            # KOKKOS package configuration
            if self.use_kokkos:
                lmp.command("package kokkos newton on neigh half")
                lmp.command("package kokkos gpu/aware off")
                LOGGER.info("KOKKOS package configured")
            else:
                lmp.command("package omp 0")
                lmp.command("suffix off")
            
            # Basic settings
            lmp.command("units real")
            lmp.command("atom_style charge")
            lmp.command("boundary p p p")
            
            # Get elements and create type mapping
            elements = self.engine.allElements
            unique_elements = sorted(set(elements))
            LOGGER.info(f"Found elements: {unique_elements}")
            self.__element_types = {elem: idx+1 for idx, elem in enumerate(unique_elements)}

            # Get box dimensions
            box_lengths = self._get_box_vectors()
            LOGGER.info(f"Creating simulation box: {box_lengths}")
            
            # Create simulation box
            lmp.command(f"region box block 0 {box_lengths[0]} 0 {box_lengths[1]} 0 {box_lengths[2]} units box")
            lmp.command(f"create_box {len(unique_elements)} box")

            # Output controls
            lmp.command("echo none")
            lmp.command("log none")
            lmp.command("thermo_style custom step pe")
            lmp.command("thermo_modify flush no norm no")
            lmp.command("thermo 0")

            # Set element masses
            LOGGER.info("Setting atomic masses...")
            for elem, type_id in self.__element_types.items():
                mass = ELEMENT_MASSES.get(elem, 0.0)
                lmp.command(f"mass {type_id} {mass}")
            
            # Setup ReaxFF with KOKKOS optimization
            LOGGER.info("Setting up ReaxFF with KOKKOS...")
            if self.use_kokkos:
                lmp.command(f"pair_style reaxff/kk {self.__control_params} safezone 30 mincap 1000")
                lmp.command("pair_modify pair reaxff/kk compute yes")
            else:
                lmp.command(f"pair_style reaxff {self.__control_params} safezone 30 mincap 1000")
                
            lmp.command("neighbor 3.0 bin")
            lmp.command("neigh_modify delay 0 every 1 check yes")

            # Setup pair coefficients
            elem_string = " ".join(unique_elements)
            lmp.command(f"pair_coeff * * {self.__reaxff_params} {elem_string}")

            # Initialize ReaxFF wrapper
            self.reaxff_wrapper = KokkosReaxFFWrapper(
                self.memory_manager, self.__reaxff_params, self.__control_params
            )

            self.__initialized = True
            self.__cached_box = box_lengths.copy()
            
            LOGGER.info("LAMMPS with KOKKOS initialized successfully")

        except Exception as e:
            self.__initialized = False
            self.__lmp = None
            LOGGER.error(f"KOKKOS-LAMMPS initialization failed: {str(e)}")
            raise LAMMPSError(str(e))

    def compute_energy(self):
        """Compute system energy using KOKKOS-enhanced ReaxFF."""
        try:
            # Initialize LAMMPS if needed
            if not self.__initialized:
                self._initialize_lammps_with_kokkos()
                
            # Get current coordinates and box
            positions = self.engine.realCoordinates
            elements = self.engine.allElements
            box_lengths = self._get_box_vectors()

            # Use KOKKOS wrapper for safe computation
            energy, final_positions = self.reaxff_wrapper.safe_energy_computation(
                self.__lmp, positions, elements, box_lengths
            )
        
            # Update engine coordinates only if computation succeeded
            if energy < 1e5:  # If not a penalty energy
                self.engine.realCoordinates[:] = final_positions

                # Transform to box coordinates
                if self.engine.isPBC:
                    final_positions_f32 = np.array(final_positions, dtype=np.float32)
                    recip_basis_f32 = np.array(self.engine.reciprocalBasisVectors, dtype=np.float32)
                    
                    box_coords = transform_coordinates(
                        transMatrix=recip_basis_f32,
                        coords=final_positions_f32
                    )
                else:
                    box_coords = final_positions            
                    
                self.engine.boxCoordinates[:] = np.array(box_coords, dtype=np.float32)

            # Update last accepted energy for HRMC
            if hasattr(self, '_KokkosEnergyConstraint__last_accepted_energy'):
                if self.__last_accepted_energy is None or energy < 1e5:
                    self.__last_accepted_energy = energy

            # Log memory statistics periodically
            if hasattr(self, '_step_count'):
                self._step_count += 1
            else:
                self._step_count = 1
                
            if self._step_count % 100 == 0:
                memory_stats = self.memory_manager.get_memory_stats()
                LOGGER.info(f"Memory stats after {self._step_count} steps: {memory_stats}")

            return energy
                
        except Exception as e:
            LOGGER.error(f"KOKKOS energy calculation failed: {str(e)}")
            # Clean up memory on error
            self.memory_manager.cleanup_all()
            raise LAMMPSError(str(e))

    def accept_move(self, realIndexes, relativeIndexes):
        """Accept the move and update last accepted energy."""
        super(KokkosEnergyConstraint, self).accept_move(realIndexes, relativeIndexes)
        if hasattr(self, 'data') and self.data is not None:
            self.__last_accepted_energy = self.data

    def compute_before_move(self, realIndexes, relativeIndexes):
        """Compute constraint before move."""
        self.compute_energy()

    def compute_after_move(self, realIndexes, relativeIndexes, movedBoxCoordinates):
        """Compute constraint after move.""" 
        self.compute_energy()

    def cleanup(self):
        """Clean up LAMMPS instance and KOKKOS memory."""
        try:
            if self.__lmp is not None:
                self.__lmp.close()
                self.__lmp = None
                
            self.memory_manager.cleanup_all()
            self.__initialized = False
            
            LOGGER.info("KOKKOS energy constraint cleanup completed")
            
        except Exception as e:
            LOGGER.error(f"Cleanup failed: {str(e)}")

    def __del__(self):
        """Destructor to ensure cleanup."""
        try:
            self.cleanup()
        except:
            pass