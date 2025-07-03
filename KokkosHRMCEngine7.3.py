"""
KokkosHRMCEngine - KOKKOS增强的HRMC引擎
集成KOKKOS内存管理和错误恢复机制
"""

import os
import numpy as np
import time
import traceback
from datetime import datetime

# Fullrmc imports
from fullrmc.Engine import Engine
from fullrmc.Core.Collection import is_number, reset_if_collected_out_of_date, get_path
from fullrmc.Core.boundary_conditions_collection import transform_coordinates
from fullrmc.Globals import LOGGER, FLOAT_TYPE

# KOKKOS imports
from EnergyConstraintsKokkos import KokkosEnergyConstraint, KokkosMemoryManager, KokkosError

class KokkosHRMCEngine(Engine):
    """
    KOKKOS增强的HRMC引擎
    
    提供以下增强功能:
    - KOKKOS内存管理
    - 智能错误恢复
    - 性能监控
    - 动态计算策略调整
    """
    
    def __init__(self, path=None, freshStart=True, useKokkos=True):
        # Initialize parent class
        super(KokkosHRMCEngine, self).__init__(path=path, freshStart=freshStart)
        
        # KOKKOS specific attributes
        self.use_kokkos = useKokkos
        self.kokkos_memory_manager = KokkosMemoryManager()
        self.computation_stats = {
            'total_steps': 0,
            'successful_steps': 0,
            'failed_steps': 0,
            'energy_computations': 0,
            'error_recoveries': 0,
            'memory_cleanups': 0
        }
        
        # Performance tracking
        self.performance_data = {
            'step_times': [],
            'energy_times': [],
            'memory_usage': [],
            'error_types': {}
        }
        
        # HRMC specific parameters
        self.hrmc_steps = 50000
        self.initial_temp = 2000
        self.final_temp = 300
        self.cooling_rate = 0.98
        self.current_cycle = 1
        self.cycles_hrmc = 5
        self.T_chi_0 = 300
        self.omega = 1000
        
        # Constraint references
        self.energy_constraint = None
        self.pdf_constraint = None
        self.dis_constraint = None
        
        # Output management
        self._output_dir = None
        self._tried = 0
        self._accepted = 0
        
        # Error handling
        self.max_consecutive_failures = 50
        self.consecutive_failures = 0
        self.adaptive_error_handling = True
        
        # Memory management settings
        self.memory_cleanup_interval = 1000
        self.force_cleanup_threshold = 500 * 1024 * 1024  # 500MB
        
        LOGGER.info(f"KokkosHRMCEngine initialized (KOKKOS: {self.use_kokkos})")

    def setup_output_files(self, output_dir):
        """设置输出文件和目录"""
        self._output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # 创建性能监控文件
        self.performance_log = os.path.join(output_dir, "kokkos_performance.log")
        self.error_log = os.path.join(output_dir, "kokkos_errors.log")
        self.memory_log = os.path.join(output_dir, "kokkos_memory.log")
        
        # 初始化日志文件
        with open(self.performance_log, 'w') as f:
            f.write("# KOKKOS Performance Log\n")
            f.write("# Step, Duration(s), Energy, Memory(MB), Status\n")
            
        with open(self.error_log, 'w') as f:
            f.write("# KOKKOS Error Log\n")
            f.write("# Timestamp, Step, Error_Type, Error_Message, Recovery_Strategy\n")
            
        with open(self.memory_log, 'w') as f:
            f.write("# KOKKOS Memory Log\n")
            f.write("# Step, Current_Usage(MB), Peak_Usage(MB), Active_Buffers, Allocations, Deallocations\n")

    def setup_kokkos_energy_constraint(self, reaxff_params, control_params, weight=1.0):
        """设置KOKKOS能量约束"""
        try:
            self.energy_constraint = KokkosEnergyConstraint(
                reaxff_params=reaxff_params,
                control_params=control_params,
                compute_forces=False,
                weight=weight,
                use_kokkos=self.use_kokkos
            )
            
            self.add_constraints([self.energy_constraint])
            LOGGER.info("KOKKOS energy constraint added successfully")
            
        except Exception as e:
            LOGGER.error(f"Failed to setup KOKKOS energy constraint: {str(e)}")
            raise

    def setup_hrmc_parameters(self, initial_temp=2000, final_temp=300, steps=50000, 
                            cycles=5, cooling_rate=0.98, omega=1000):
        """设置HRMC参数"""
        self.initial_temp = initial_temp
        self.final_temp = final_temp
        self.hrmc_steps = steps
        self.cycles_hrmc = cycles
        self.cooling_rate = cooling_rate
        self.omega = omega
        self.T_chi_0 = initial_temp
        
        LOGGER.info(f"HRMC parameters: T_initial={initial_temp}K, T_final={final_temp}K, "
                   f"steps={steps}, cycles={cycles}")

    def run_kokkos_hrmc_step(self):
        """执行单个KOKKOS增强的HRMC步骤"""
        step_start_time = time.time()
        
        try:
            # 选择原子组进行移动
            group_index = self.groupSelector.select()
            if group_index is None:
                self.computation_stats['failed_steps'] += 1
                return False
                
            group = self.groups[group_index]
            relative_indexes = group.indexes
            real_indexes = self.get_relative_indexes(relative_indexes, group_index)
            
            # 存储原始坐标
            old_positions = self.realCoordinates[relative_indexes].copy()
            
            # 应用移动
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
            
            # KOKKOS增强的约束计算
            success = self._compute_constraints_with_kokkos_recovery(
                real_indexes, relative_indexes, moved_positions, old_positions
            )
            
            if success:
                # 检查移动接受
                reject_move = self.should_step_get_rejected(real_indexes, relative_indexes)
                
                if reject_move:
                    self._reject_move(real_indexes, relative_indexes, old_positions, group_index)
                else:
                    self._accept_move(real_indexes, relative_indexes, group_index)
                    self.consecutive_failures = 0
                    
            else:
                # 约束计算失败，拒绝移动
                self._reject_move(real_indexes, relative_indexes, old_positions, group_index)
                self.consecutive_failures += 1
                
            self.computation_stats['total_steps'] += 1
            self.increment_tried()
            
            # 记录性能数据
            step_duration = time.time() - step_start_time
            self._log_step_performance(step_duration, success)
            
            # 定期内存清理
            if self._tried % self.memory_cleanup_interval == 0:
                self._perform_memory_cleanup()
                
            # 检查连续失败次数
            if self.consecutive_failures >= self.max_consecutive_failures:
                LOGGER.warn(f"连续{self.consecutive_failures}次失败，触发系统重置")
                self._emergency_system_reset()
                
            return success
            
        except Exception as e:
            self._handle_step_error(e, step_start_time)
            return False

    def _compute_constraints_with_kokkos_recovery(self, real_indexes, relative_indexes, 
                                                moved_positions, old_positions):
        """使用KOKKOS错误恢复机制计算约束"""
        energy_start_time = time.time()
        
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
            
            # 记录能量计算时间
            energy_duration = time.time() - energy_start_time
            self.performance_data['energy_times'].append(energy_duration)
            self.computation_stats['energy_computations'] += 1
            
            return True
            
        except Exception as e:
            error_msg = str(e)
            
            # 恢复原始坐标
            self.realCoordinates[relative_indexes] = old_positions
            
            # 分类错误处理
            if "bondchk failed" in error_msg or "malformed" in error_msg:
                LOGGER.warn(f"ReaxFF结构错误: {error_msg}")
                self._log_error("ReaxFF_Structure_Error", error_msg, "Position_Restoration")
                
            elif isinstance(e, KokkosError):
                LOGGER.warn(f"KOKKOS计算错误: {error_msg}")
                self._log_error("KOKKOS_Computation_Error", error_msg, "Memory_Cleanup")
                self._perform_emergency_memory_cleanup()
                
            else:
                LOGGER.error(f"未知约束计算错误: {error_msg}")
                self._log_error("Unknown_Constraint_Error", error_msg, "System_Reset")
                
            self.computation_stats['error_recoveries'] += 1
            return False

    def _accept_move(self, real_indexes, relative_indexes, group_index):
        """接受移动"""
        self.increment_accepted()
        self.groupSelector.move_accepted(group_index)
        
        for c in self.constraints:
            try:
                c.accept_move(
                    realIndexes=real_indexes,
                    relativeIndexes=relative_indexes
                )
            except Exception as e:
                LOGGER.warn(f"约束接受移动失败: {str(e)}")
                
        self.computation_stats['successful_steps'] += 1

    def _reject_move(self, real_indexes, relative_indexes, old_positions, group_index):
        """拒绝移动"""
        self.realCoordinates[relative_indexes] = old_positions
        self.groupSelector.move_rejected(group_index)
        
        for c in self.constraints:
            try:
                c.reject_move(
                    realIndexes=real_indexes,
                    relativeIndexes=relative_indexes
                )
            except Exception as e:
                # 忽略拒绝移动时的错误，确保系统稳定性
                pass
                
        self.computation_stats['failed_steps'] += 1

    def should_step_get_rejected(self, realIndexes, relativeIndexes):
        """KOKKOS增强的HRMC接受准则"""
        try:
            # 获取当前能量
            current_energy = self.energy_constraint.data
            last_energy = getattr(self.energy_constraint, 
                                '_KokkosEnergyConstraint__last_accepted_energy', None)
            
            if last_energy is None:
                last_energy = current_energy
                
            # 计算能量差
            delta_U = current_energy - last_energy
            
            # 处理PDF约束（如果存在）
            delta_chi_sq = 0.0
            if self.pdf_constraint:
                try:
                    old_hist = self.pdf_constraint.data
                    new_hist, _ = self.pdf_constraint.compute_data(update=False)
                    
                    if old_hist is not None and new_hist is not None:
                        # 计算χ²差值
                        exp_data = self.pdf_constraint.experimentalData
                        exp_gr = exp_data[:, 1]
                        
                        def compute_chi_square(hist):
                            return np.sum((hist - exp_gr)**2)
                        
                        chi_sq_new = compute_chi_square(new_hist)
                        chi_sq_old = compute_chi_square(old_hist)
                        delta_chi_sq = chi_sq_new - chi_sq_old
                        
                except Exception as e:
                    LOGGER.warn(f"PDF约束计算失败: {str(e)}")
                    delta_chi_sq = 0.0
            
            # 计算当前温度
            T_chi = self.T_chi_0 * (self.cooling_rate ** self.current_cycle)
            
            # 计算接受概率
            moved_elem = self.allElements[realIndexes[0]]
            is_hydrogen = (moved_elem == 'H')
            
            if is_hydrogen:
                # 氢原子：只考虑能量项
                delta = (delta_U / self.omega) / T_chi
            else:
                # 非氢原子：考虑结构和能量项
                delta = (delta_chi_sq + delta_U / self.omega) / T_chi
            
            # 处理数值溢出
            if delta > 700:
                self._write_acceptance_data(0.0, delta_chi_sq, delta_U, T_chi)
                return True
                
            if delta < -700:
                self._write_acceptance_data(1.0, delta_chi_sq, delta_U, T_chi)
                return False
            
            # 计算接受概率
            prob = min(1.0, np.exp(-delta))
            
            # 记录数据
            self._write_acceptance_data(prob, delta_chi_sq, delta_U, T_chi)
            
            # 做出决定
            accepted = np.random.random() <= prob
            
            if accepted:
                # 更新最后接受的能量
                if hasattr(self.energy_constraint, '_KokkosEnergyConstraint__last_accepted_energy'):
                    self.energy_constraint._KokkosEnergyConstraint__last_accepted_energy = current_energy
                    
            return not accepted
            
        except Exception as e:
            LOGGER.error(f"HRMC接受准则计算失败: {str(e)}")
            return True  # 出错时拒绝移动

    def _write_acceptance_data(self, prob, delta_chi_sq, delta_U, T_chi):
        """写入接受数据"""
        try:
            if hasattr(self, '_output_dir') and self._output_dir:
                acceptance_file = os.path.join(self._output_dir, 'HRMC_acceptance.txt')
                energy_file = os.path.join(self._output_dir, 'HRMC_energy.txt')
                
                with open(acceptance_file, 'a') as f:
                    f.write(f"{self._tried} {prob:.6f} {delta_chi_sq:.6f} {delta_U:.6f} {T_chi:.6f}\n")
                
                if hasattr(self.energy_constraint, 'data'):
                    with open(energy_file, 'a') as f:
                        f.write(f"{self._tried} {self.energy_constraint.data:.6f} {delta_U:.6f}\n")
                        
        except Exception as e:
            LOGGER.warn(f"写入接受数据失败: {str(e)}")

    def _perform_memory_cleanup(self):
        """执行内存清理"""
        try:
            # 检查内存使用情况
            memory_stats = self.kokkos_memory_manager.get_memory_stats()
            current_memory = memory_stats['current_usage']
            
            if current_memory > self.force_cleanup_threshold:
                LOGGER.info(f"内存使用过高 ({current_memory/1024/1024:.1f}MB)，强制清理")
                self._perform_emergency_memory_cleanup()
            else:
                # 常规清理：只清理不常用的缓冲区
                self.kokkos_memory_manager.cleanup_all()
                
            # 记录内存状态
            self._log_memory_stats(memory_stats)
            self.computation_stats['memory_cleanups'] += 1
            
        except Exception as e:
            LOGGER.error(f"内存清理失败: {str(e)}")

    def _perform_emergency_memory_cleanup(self):
        """紧急内存清理"""
        try:
            LOGGER.warn("执行紧急内存清理")
            
            # 清理KOKKOS内存管理器
            self.kokkos_memory_manager.cleanup_all()
            
            # 重新初始化能量约束（如果需要）
            if self.energy_constraint and hasattr(self.energy_constraint, 'cleanup'):
                self.energy_constraint.cleanup()
                
            # 强制垃圾回收
            import gc
            gc.collect()
            
            LOGGER.info("紧急内存清理完成")
            
        except Exception as e:
            LOGGER.error(f"紧急内存清理失败: {str(e)}")

    def _emergency_system_reset(self):
        """紧急系统重置"""
        try:
            LOGGER.warn("执行紧急系统重置")
            
            # 重置连续失败计数
            self.consecutive_failures = 0
            
            # 执行内存清理
            self._perform_emergency_memory_cleanup()
            
            # 重新初始化能量约束
            if self.energy_constraint:
                try:
                    self.energy_constraint.cleanup()
                    # 重新初始化
                    self.energy_constraint._KokkosEnergyConstraint__initialized = False
                    LOGGER.info("能量约束重新初始化完成")
                except Exception as e:
                    LOGGER.error(f"能量约束重新初始化失败: {str(e)}")
            
            LOGGER.info("紧急系统重置完成")
            
        except Exception as e:
            LOGGER.error(f"紧急系统重置失败: {str(e)}")

    def _log_step_performance(self, duration, success):
        """记录步骤性能"""
        self.performance_data['step_times'].append(duration)
        
        # 获取内存使用情况
        if hasattr(self.kokkos_memory_manager, 'current_memory_usage'):
            memory_mb = self.kokkos_memory_manager.current_memory_usage / (1024 * 1024)
            self.performance_data['memory_usage'].append(memory_mb)
        
        # 写入性能日志
        try:
            energy = getattr(self.energy_constraint, 'data', 0.0) if self.energy_constraint else 0.0
            status = "SUCCESS" if success else "FAILED"
            memory_mb = memory_mb if 'memory_mb' in locals() else 0.0
            
            with open(self.performance_log, 'a') as f:
                f.write(f"{self._tried} {duration:.6f} {energy:.6f} {memory_mb:.2f} {status}\n")
                
        except Exception as e:
            LOGGER.warn(f"性能日志写入失败: {str(e)}")

    def _log_error(self, error_type, error_message, recovery_strategy):
        """记录错误信息"""
        try:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # 更新错误统计
            if error_type not in self.performance_data['error_types']:
                self.performance_data['error_types'][error_type] = 0
            self.performance_data['error_types'][error_type] += 1
            
            # 写入错误日志
            with open(self.error_log, 'a') as f:
                f.write(f"{timestamp} {self._tried} {error_type} \"{error_message}\" {recovery_strategy}\n")
                
        except Exception as e:
            LOGGER.warn(f"错误日志写入失败: {str(e)}")

    def _log_memory_stats(self, memory_stats):
        """记录内存统计"""
        try:
            current_mb = memory_stats['current_usage'] / (1024 * 1024)
            peak_mb = memory_stats['peak_usage'] / (1024 * 1024)
            
            with open(self.memory_log, 'a') as f:
                f.write(f"{self._tried} {current_mb:.2f} {peak_mb:.2f} "
                       f"{memory_stats['active_buffers']} "
                       f"{memory_stats['allocation_count']} "
                       f"{memory_stats['deallocation_count']}\n")
                
        except Exception as e:
            LOGGER.warn(f"内存日志写入失败: {str(e)}")

    def _handle_step_error(self, error, step_start_time):
        """处理步骤错误"""
        error_msg = str(error)
        step_duration = time.time() - step_start_time
        
        LOGGER.error(f"HRMC步骤执行失败: {error_msg}")
        LOGGER.error(f"错误堆栈: {traceback.format_exc()}")
        
        # 记录错误
        self._log_error("Step_Execution_Error", error_msg, "Continue")
        
        # 更新统计
        self.computation_stats['failed_steps'] += 1
        self.consecutive_failures += 1
        
        # 记录性能（失败情况）
        self._log_step_performance(step_duration, False)

    def run_kokkos_hrmc_cycle(self):
        """执行完整的KOKKOS HRMC循环"""
        try:
            LOGGER.info(f"开始KOKKOS HRMC循环 {self.current_cycle}/{self.cycles_hrmc}")
            cycle_start_time = time.time()
            
            # 重置统计
            cycle_stats = {
                'accepted': 0,
                'tried': 0,
                'errors': 0
            }
            
            # 执行HRMC步骤
            for step in range(self.hrmc_steps):
                success = self.run_kokkos_hrmc_step()
                
                cycle_stats['tried'] += 1
                if success and hasattr(self, '_accepted'):
                    cycle_stats['accepted'] = self._accepted
                if not success:
                    cycle_stats['errors'] += 1
                
                # 定期报告进度
                if step % 1000 == 0:
                    acceptance_rate = (cycle_stats['accepted'] / cycle_stats['tried']) * 100 if cycle_stats['tried'] > 0 else 0
                    LOGGER.info(f"循环 {self.current_cycle}, 步骤 {step}/{self.hrmc_steps}, "
                               f"接受率: {acceptance_rate:.1f}%, 错误: {cycle_stats['errors']}")
            
            cycle_duration = time.time() - cycle_start_time
            
            # 循环结束统计
            final_acceptance_rate = (cycle_stats['accepted'] / cycle_stats['tried']) * 100 if cycle_stats['tried'] > 0 else 0
            LOGGER.info(f"循环 {self.current_cycle} 完成: 时间={cycle_duration:.2f}s, "
                       f"接受率={final_acceptance_rate:.1f}%, 错误={cycle_stats['errors']}")
            
            # 导出结构
            if self._output_dir:
                cycle_pdb = os.path.join(self._output_dir, f"kokkos_cycle_{self.current_cycle}.pdb")
                self.export_pdb(cycle_pdb)
                
            # 生成性能报告
            self._generate_cycle_performance_report()
            
        except Exception as e:
            LOGGER.error(f"KOKKOS HRMC循环失败: {str(e)}")
            raise

    def run_kokkos_hrmc_simulation(self):
        """运行完整的KOKKOS HRMC模拟"""
        try:
            simulation_start_time = time.time()
            LOGGER.info("开始KOKKOS HRMC模拟")
            
            for cycle in range(1, self.cycles_hrmc + 1):
                self.current_cycle = cycle
                self.run_kokkos_hrmc_cycle()
                
                # 循环间内存清理
                self._perform_memory_cleanup()
            
            simulation_duration = time.time() - simulation_start_time
            
            # 生成最终报告
            self._generate_final_simulation_report(simulation_duration)
            
            LOGGER.info(f"KOKKOS HRMC模拟完成，总时间: {simulation_duration:.2f}秒")
            
        except Exception as e:
            LOGGER.error(f"KOKKOS HRMC模拟失败: {str(e)}")
            raise
        finally:
            # 确保清理
            self.cleanup()

    def _generate_cycle_performance_report(self):
        """生成循环性能报告"""
        try:
            if not self._output_dir:
                return
                
            report_file = os.path.join(self._output_dir, f"cycle_{self.current_cycle}_performance.txt")
            
            with open(report_file, 'w') as f:
                f.write(f"=== KOKKOS HRMC 循环 {self.current_cycle} 性能报告 ===\n\n")
                
                # 基本统计
                f.write("基本统计:\n")
                f.write(f"  总步数: {self.computation_stats['total_steps']}\n")
                f.write(f"  成功步数: {self.computation_stats['successful_steps']}\n")
                f.write(f"  失败步数: {self.computation_stats['failed_steps']}\n")
                f.write(f"  能量计算次数: {self.computation_stats['energy_computations']}\n")
                f.write(f"  错误恢复次数: {self.computation_stats['error_recoveries']}\n")
                f.write(f"  内存清理次数: {self.computation_stats['memory_cleanups']}\n\n")
                
                # 性能统计
                if self.performance_data['step_times']:
                    step_times = self.performance_data['step_times']
                    f.write("步骤时间统计:\n")
                    f.write(f"  平均时间: {np.mean(step_times):.6f} 秒\n")
                    f.write(f"  最小时间: {np.min(step_times):.6f} 秒\n")
                    f.write(f"  最大时间: {np.max(step_times):.6f} 秒\n")
                    f.write(f"  标准差: {np.std(step_times):.6f} 秒\n\n")
                
                if self.performance_data['energy_times']:
                    energy_times = self.performance_data['energy_times']
                    f.write("能量计算时间统计:\n")
                    f.write(f"  平均时间: {np.mean(energy_times):.6f} 秒\n")
                    f.write(f"  最小时间: {np.min(energy_times):.6f} 秒\n")
                    f.write(f"  最大时间: {np.max(energy_times):.6f} 秒\n\n")
                
                # 错误统计
                if self.performance_data['error_types']:
                    f.write("错误类型统计:\n")
                    for error_type, count in self.performance_data['error_types'].items():
                        f.write(f"  {error_type}: {count}\n")
                    f.write("\n")
                
                # 内存统计
                memory_stats = self.kokkos_memory_manager.get_memory_stats()
                f.write("内存使用统计:\n")
                f.write(f"  当前使用: {memory_stats['current_usage']/(1024*1024):.2f} MB\n")
                f.write(f"  峰值使用: {memory_stats['peak_usage']/(1024*1024):.2f} MB\n")
                f.write(f"  活动缓冲区: {memory_stats['active_buffers']}\n")
                f.write(f"  分配次数: {memory_stats['allocation_count']}\n")
                f.write(f"  释放次数: {memory_stats['deallocation_count']}\n")
                
        except Exception as e:
            LOGGER.warn(f"性能报告生成失败: {str(e)}")

    def _generate_final_simulation_report(self, total_duration):
        """生成最终模拟报告"""
        try:
            if not self._output_dir:
                return
                
            report_file = os.path.join(self._output_dir, "kokkos_final_report.txt")
            
            with open(report_file, 'w') as f:
                f.write("=== KOKKOS HRMC 最终模拟报告 ===\n\n")
                
                f.write(f"模拟参数:\n")
                f.write(f"  使用KOKKOS: {self.use_kokkos}\n")
                f.write(f"  总循环数: {self.cycles_hrmc}\n")
                f.write(f"  每循环步数: {self.hrmc_steps}\n")
                f.write(f"  初始温度: {self.initial_temp} K\n")
                f.write(f"  最终温度: {self.final_temp} K\n")
                f.write(f"  冷却率: {self.cooling_rate}\n\n")
                
                f.write(f"性能统计:\n")
                f.write(f"  总模拟时间: {total_duration:.2f} 秒\n")
                f.write(f"  平均每步时间: {total_duration/max(1, self.computation_stats['total_steps']):.6f} 秒\n")
                
                if self.computation_stats['total_steps'] > 0:
                    success_rate = (self.computation_stats['successful_steps'] / self.computation_stats['total_steps']) * 100
                    f.write(f"  成功率: {success_rate:.2f}%\n")
                
                f.write(f"  总错误恢复次数: {self.computation_stats['error_recoveries']}\n")
                f.write(f"  总内存清理次数: {self.computation_stats['memory_cleanups']}\n\n")
                
                # 最终内存状态
                final_memory_stats = self.kokkos_memory_manager.get_memory_stats()
                f.write(f"最终内存状态:\n")
                f.write(f"  峰值内存使用: {final_memory_stats['peak_usage']/(1024*1024):.2f} MB\n")
                f.write(f"  最终内存使用: {final_memory_stats['current_usage']/(1024*1024):.2f} MB\n")
                f.write(f"  内存分配效率: {final_memory_stats['deallocation_count']/max(1, final_memory_stats['allocation_count'])*100:.1f}%\n")
                
        except Exception as e:
            LOGGER.warn(f"最终报告生成失败: {str(e)}")

    def cleanup(self):
        """清理资源"""
        try:
            LOGGER.info("开始KOKKOS HRMC引擎清理")
            
            # 清理能量约束
            if self.energy_constraint and hasattr(self.energy_constraint, 'cleanup'):
                self.energy_constraint.cleanup()
                
            # 清理KOKKOS内存管理器
            self.kokkos_memory_manager.cleanup_all()
            
            # 调用父类清理
            if hasattr(super(), 'cleanup'):
                super().cleanup()
                
            LOGGER.info("KOKKOS HRMC引擎清理完成")
            
        except Exception as e:
            LOGGER.error(f"清理失败: {str(e)}")

    def __del__(self):
        """析构函数"""
        try:
            self.cleanup()
        except:
            pass