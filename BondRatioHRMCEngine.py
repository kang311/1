"""
BondRatioHRMCEngine implements a purely structure-driven HRMC model 
combining PDF constraints and bond ratio constraints.

This engine focuses on structural optimization without energy terms:
- Pair Distribution Function (G(r)) fitting
- Bond ratio distributions fitting (akin to bond constraints)
"""
from __future__ import print_function
import os
import numpy as np
import time
import traceback
from datetime import datetime
from collections import defaultdict

from fullrmc.Engine import Engine 
from fullrmc.Globals import LOGGER, FLOAT_TYPE
from fullrmc.Core.boundary_conditions_collection import transform_coordinates
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.BondRatioConstraints import BondRatioConstraint

class BondRatioHRMCEngine(Engine):
    """
    结合PDF约束和键比例约束的HRMC引擎。
    该引擎不使用能量约束，而是通过结构特征驱动模拟。
    
    :Parameters:
        #. path (str): 引擎保存路径
        #. freshStart (bool): 是否重新开始模拟
    """
    def __init__(self, path, freshStart=False):
        # 初始化父类
        super(BondRatioHRMCEngine, self).__init__(path=path, freshStart=freshStart)
        
        # Constants - 使用与MDHRMCEngine相同的常量
        self.KB = 0.001987204258  # Boltzmann constant in kcal/mol/K
        
        # 退火参数 - 参考run.py中的设置
        self.initial_temperature = 800.0    # 初始退火温度 (K)
        self.final_temperature = 300.0      # 最终退火温度 (K)
        self.current_temperature = self.initial_temperature  # 当前温度
        self.cooling_rate = 0.98            # 冷却率
        self.current_cycle = 1              # 当前周期
        
        # 约束权重参数 - 参考run.py中的设置
        self.w1 = 0.1  # G(r)权重
        self.w2 = 1000  # 键比例约束权重
        
        # 温度转换系数
        self.T_chi_0 = 1.0  # 初始归一化温度
        
        # 约束引用
        self.pdf_constraint = None
        self.bond_constraint = None
        
        # 统计计数器
        self._tried = 1  # 从1开始
        self._accepted = 0
        
        # 初始化输出目录和文件
        self._output_dir = os.path.dirname(path) if path else None 
        self._files = {}
        
        # 设置输出文件（如果路径已提供）
        if self._output_dir:
            self.setup_output_files(self._output_dir)
            
        # 缓存数据
        self._cached_data = {}
        self._lastBondCounts = None
        self._lastTotalBonds = None
        
        # 记录日志
        LOGGER.info(f"BondRatioHRMCEngine 初始化完成，输出目录: {self._output_dir}")
    
    def __getstate__(self):
        """自定义序列化。"""
        state = super(BondRatioHRMCEngine, self).__getstate__()
        # 不需要序列化的临时数据
        state['_cached_data'] = {}
        state['_lastBondCounts'] = None
        state['_lastTotalBonds'] = None
        return state
        
    def __setstate__(self, state):
        """处理反序列化。"""
        self.__dict__.update(state)
        # 初始化缓存数据
        if not hasattr(self, '_cached_data'):
            self._cached_data = {}
        # 恢复后其他可能需要重新初始化的属性
        if not hasattr(self, '_lastBondCounts'):
            self._lastBondCounts = None
        if not hasattr(self, '_lastTotalBonds'):
            self._lastTotalBonds = None

    def cleanup(self):
        """清理资源。"""
        # 关闭可能打开的文件
        for filename in self._files.keys():
            filepath = os.path.join(self._output_dir, filename) if self._output_dir else None
            if filepath and os.path.exists(filepath):
                try:
                    LOGGER.debug(f"清理资源: {filepath}")
                except:
                    pass
        
        # 清理缓存数据
        self._cached_data = {}
        if hasattr(self, '_lastBondCounts'):
            delattr(self, '_lastBondCounts')
        if hasattr(self, '_lastTotalBonds'):
            delattr(self, '_lastTotalBonds')
        
        LOGGER.info("BondRatioHRMCEngine资源已清理")
            
    def setup_output_files(self, output_dir):
        """设置输出目录和文件。"""
        try:
            self._output_dir = output_dir
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # 定义文件和表头 - 与MDHRMCEngine的输出文件格式类似
            self._files = {
                'HRMC_acceptance.txt': "# Step Probability DeltaGr DeltaBond T_chi",
                'HRMC_gr.txt': "# Step StandardError GrDiff",
                'HRMC_bond_ratios.txt': "# Step StandardError DeltaBond"
            }

            # 创建带表头的文件
            for filename, header in self._files.items():
                filepath = os.path.join(output_dir, filename)
                with open(filepath, 'w') as f:
                    f.write(header + '\n')
                    
            LOGGER.info(f"输出文件已设置在 {output_dir}")
                    
        except Exception as e:
            LOGGER.error(f"设置输出文件失败: {str(e)}")
            raise

    def _write_data(self, filename, data):
        """使用上下文管理器将数据写入文件。"""
        try:
            if self._output_dir is None:
                LOGGER.error("输出目录未设置")
                return
                
            filepath = os.path.join(self._output_dir, filename)
            LOGGER.debug(f"尝试写入到 {filepath}")

            # 验证文件是否存在
            if not os.path.exists(filepath):
                with open(filepath, 'w') as f:
                    f.write(self._files.get(filename, "# Data") + '\n')
                    
            # 附加数据
            with open(filepath, 'a') as f:
                f.write(data + '\n')

            LOGGER.debug(f"成功将数据写入 {filepath}")    
                
        except Exception as e:
            LOGGER.error(f"写入 {filename} 时出错: {str(e)}")
            
    def increment_tried(self):
        """增加尝试移动次数，从1开始。"""
        if self._tried == 0:  # 如果是0，从1开始
            self._tried = 1    # 这是第一次移动
        else:
            self._tried += 1   # 增加尝试次数

    def increment_accepted(self):
        """增加接受移动次数。"""  
        self._accepted += 1
            
    def add_constraints(self, constraints):
        """
        添加并正确初始化约束。
        
        :Parameters:
            #. constraints (list): 约束对象列表
        """
        try:
            for constraint in constraints:
                if isinstance(constraint, PairDistributionConstraint):
                    self.pdf_constraint = constraint
                    LOGGER.info(f"已添加PDF约束: {constraint.__class__.__name__}")
                elif isinstance(constraint, BondRatioConstraint):
                    self.bond_constraint = constraint
                    LOGGER.info(f"已添加键比例约束: {constraint.__class__.__name__}")
            
            # 在存储引用后调用父类的add_constraints方法
            super(BondRatioHRMCEngine, self).add_constraints(constraints)
            
        except Exception as e:
            LOGGER.error(f"添加约束出错: {str(e)}")
            raise
                
    def initialize_used_constraints(self):
        """安全地初始化约束。"""
        try:
            # 初始化基本约束
            super(BondRatioHRMCEngine, self).initialize_used_constraints()
            
            # 验证我们需要的特定约束
            if not hasattr(self, 'pdf_constraint') or self.pdf_constraint is None:
                raise ValueError("PDF约束未找到")
            if not hasattr(self, 'bond_constraint') or self.bond_constraint is None:
                raise ValueError("键比例约束(bond)未找到")
                
            # 计算温度转换系数
            elements = self.allElements
            N_total = len(elements)
            N_nonH = sum(1 for elem in elements if elem != 'H')
            
            if N_nonH == 0:
                raise ValueError("结构中必须包含非氢原子")                
                                  
            LOGGER.info(f"结构加载: 总原子数={N_total}, 非氢原子数={N_nonH}")
               
        except Exception as e:
            LOGGER.error(f"初始化约束出错: {str(e)}")
            raise
            
    def setup_hrmc(self, initial_temp, final_temp, steps,warmup_steps=200,target_accept_rate=0.4):
        """
        设置HRMC参数，并根据初期扰动采样计算归一初始温度 T_chi_0。

        Parameters:
            initial_temp (float): 初始物理温度 (K)
            final_temp (float): 最终物理温度 (K)
            steps (int): 总模拟步数
            warmup_steps (int): 前期采样步骤数，用于估算 Δχ²
            target_accept_rate (float): 初始接受率目标，默认 0.4
        """
        self.initial_temperature = initial_temp
        self.final_temperature = final_temp
        self.current_temperature = initial_temp
        self.warmup_steps = warmup_steps
        self.target_accept_rate = target_accept_rate
        self.steps = steps

        # 采样 Δχ²_typ
        delta_list = []
        for _ in range(warmup_steps):
            realIndexes, relativeIndexes = self.random_select_atoms()
        
            # 记录移动前误差
            gr_before = self.pdf_constraint.standardError
            bond_before = self.bond_constraint.standardError

            # 执行移动
            self._apply_random_move(realIndexes, relativeIndexes)

            # 获取移动后误差
            gr_after = self.pdf_constraint.standardError
            bond_after = self.bond_constraint.standardError

            # 初始误差（避免除0）
            gr0 = max(self.pdf_constraint.initialError, 1e-8)
            bond0 = max(self.bond_constraint.initialError, 1e-8)

            delta_gr = (gr_after - gr_before) / gr0
            delta_bond = (bond_after - bond_before) / bond0

            delta_total = self.w1 * delta_gr + self.w2 * delta_bond
            delta_list.append(abs(delta_total))

            # 撤销移动
            self.undo_move()

        delta_typ = np.mean(delta_list)
        self.T_chi_0 = delta_typ / np.log(1.0 / target_accept_rate)

        self._update_generators_temperature(self.current_temperature)

        LOGGER.info(f"HRMC参数已设置：")
        LOGGER.info(f"  初始温度 T_chi_0 = {self.T_chi_0:.6f}")
        LOGGER.info(f"  总步数 = {steps}, 冷却率 = {self.cooling_rate:.6f}")
   
    def set_constraint_weights(self, pdf_weight=1, bond_weight=1):
        """
        设置约束权重
        
        :Parameters:
            #. pdf_weight (float): PDF约束权重
            #. bond_weight (float): 键比例约束权重
        """
        if pdf_weight < 0 or bond_weight < 0:
            raise ValueError("约束权重必须为非负值")
            
        self.w1 = pdf_weight
        self.w2 = bond_weight
        
        LOGGER.info(f"已设置约束权重: PDF={pdf_weight}, 键比例={bond_weight}")
        
    def _update_generators_temperature(self, temperature):
        """更新所有移动生成器的温度。"""
        try:
            for group in self.groups:
                if hasattr(group, 'moveGenerator') and hasattr(group.moveGenerator, 'set_temperature'):
                    group.moveGenerator.set_temperature(temperature)
                    
        except Exception as e:
            LOGGER.warning(f"更新移动生成器温度时出错: {str(e)}")
        
    def run(self, numberOfSteps, saveFrequency=10000):
        """
        运行HRMC步骤。
        
        :Parameters:
            #. numberOfSteps (int): 模拟步数
            #. saveFrequency (int): 保存频率
        """
        try:
            # 记录开始时间
            start_time = time.time()
            LOGGER.info(f"开始HRMC模拟，步数: {numberOfSteps}")
            
            # 确保约束已初始化
            if not self._are_constraints_initialized():
                LOGGER.info("重新初始化约束...")
                self.initialize_used_constraints()
            
            # 计算初始误差值
            initial_pdf_error = self.pdf_constraint.standardError
            initial_bond_error = self.bond_constraint.standardError
            
            LOGGER.info(f"初始PDF误差: {initial_pdf_error:.6f}")
            LOGGER.info(f"初始键比例误差: {initial_bond_error:.6f}")
            
            # 统计接受率
            local_tried = 0
            local_accepted = 0
            
            # 开始模拟循环
            for step in range(1, numberOfSteps + 1):
                # 选择要移动的组
                group_index = self.groupSelector.select_index()
                if group_index is None:
                    continue
                group = self.groups[group_index]
                
                # 获取组信息以便日志
                group_size = len(group.indexes)
                group_elements = [self.allElements[idx] for idx in group.indexes]
                element_counts = {}
                for elem in group_elements:
                    element_counts[elem] = element_counts.get(elem, 0) + 1
                
                # 转换索引为numpy数组
                real_indexes = np.array(group.indexes, dtype=np.int32)
                relative_indexes = np.array([
                    self._atomsCollector.get_relative_index(idx) 
                    for idx in group.indexes
                ], dtype=np.int32)
                
                # 存储旧位置
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
                        
                # 计算移动前的约束
                for c in self.constraints:
                    c.compute_before_move(
                        realIndexes=real_indexes,
                        relativeIndexes=relative_indexes
                    )
                
                # 更新坐标
                self.realCoordinates[relative_indexes] = moved_positions
                
                # 计算移动后的约束
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
                
                # 检查是否接受移动
                rejected = self.should_step_get_rejected(real_indexes, relative_indexes)
                
                # 更新局部统计
                local_tried += 1
                
                if rejected:
                    # 拒绝移动
                    self.realCoordinates[relative_indexes] = old_positions
                    self.groupSelector.move_rejected(group_index)
                    for c in self.constraints:
                        c.reject_move(
                            realIndexes=real_indexes,
                            relativeIndexes=relative_indexes
                        )
                        
                    if step % 500 == 0:  # 减少日志频率
                        LOGGER.debug(f"步骤 {step}: 移动被拒绝 (组 {group_index}，元素 {element_counts})")
                else:
                    # 接受移动
                    self.increment_accepted()
                    local_accepted += 1
                    self.groupSelector.move_accepted(group_index)
                    for c in self.constraints:
                        c.accept_move(
                            realIndexes=real_indexes,
                            relativeIndexes=relative_indexes
                        )
                        
                    if step % 500 == 0:  # 减少日志频率
                        LOGGER.debug(f"步骤 {step}: 移动被接受 (组 {group_index}，元素 {element_counts})")
                
                self.increment_tried()
                
                if step % 100 == 0:
                    self.T_chi = self.T_chi_0 * (self.cooling_rate ** self.current_cycle)
                    self._update_generators_temperature(self.T_chi)
                    self.current_cycle += 1  # 更新冷却轮数
                
                # 定期保存结果
                if step % saveFrequency == 0 or step == numberOfSteps:
                    # 计算接受率
                    if local_tried > 0:
                        local_rate = local_accepted / local_tried
                    else:
                        local_rate = 0
                        
                    # 全局接受率
                    global_rate = self.accepted / max(1, self.tried)
                        
                    # 重置局部统计
                    local_tried = 0
                    local_accepted = 0
                    
                    # 计时信息
                    elapsed = time.time() - start_time
                    remaining = (elapsed / step) * (numberOfSteps - step)
                    
                    # 保存状态
                    save_path = self.path if self.path else "engine_checkpoint.rmc"
                    self.save(path=save_path)
                    
                    # 计算当前误差
                    pdf_error = self.pdf_constraint.standardError
                    bond_error = self.bond_constraint.standardError
                    
                    # 计算改进率
                    pdf_improvement = (initial_pdf_error - pdf_error) / max(initial_pdf_error, 1e-10) * 100
                    bond_improvement = (initial_bond_error - bond_error) / max(initial_bond_error, 1e-10) * 100
                    
                    LOGGER.info(f"步骤 {step}/{numberOfSteps} ({step/numberOfSteps*100:.1f}%)")
                    LOGGER.info(f"- 接受率: 局部={local_rate:.2%}, 全局={global_rate:.2%}")
                    LOGGER.info(f"- PDF误差: {pdf_error:.6f} (改进: {pdf_improvement:.1f}%)")
                    LOGGER.info(f"- 键比例误差: {bond_error:.6f} (改进: {bond_improvement:.1f}%)")
                    LOGGER.info(f"- 当前归一温度 T_chi: {self.T_chi:.6f}")
                    LOGGER.info(f"- 已用时间: {elapsed:.1f}秒, 剩余: {remaining:.1f}秒")
                    
                    # 导出当前PDB
                    pdb_path = os.path.join(self._output_dir, f"structure_step_{step}.pdb")
                    self.export_pdb(path=pdb_path)
                    
                    # 导出当前键比例和G(r)图
                    self.export_analysis_plots(step)
                
            # 完成模拟
            total_time = time.time() - start_time
            overall_rate = self.accepted / max(1, self.tried)
            
            LOGGER.info(f"模拟完成，总步数: {numberOfSteps}，总用时: {total_time:.1f}秒")
            LOGGER.info(f"总体接受率: {overall_rate:.2%} ({self.accepted}/{self.tried})")
            LOGGER.info(f"最终PDF误差: {self.pdf_constraint.standardError:.6f}")
            LOGGER.info(f"最终键比例误差: {self.bond_constraint.standardError:.6f}")
            
            # 保存最终结果
            self.save()
            
            # 导出最终PDB
            final_pdb = os.path.join(self._output_dir, "final_structure.pdb")
            self.export_pdb(path=final_pdb)
            LOGGER.info(f"最终结构已导出到: {final_pdb}")
            
            # 生成最终报告
            report_path = os.path.join(self._output_dir, "final_report.txt")
            self.generate_summary_report(filename=report_path)
            LOGGER.info(f"最终报告已生成: {report_path}")
            
            # 导出最终分析图
            self.export_analysis_plots("final")
            
            return True
            
        except Exception as e:
            LOGGER.error(f"HRMC模拟运行出错: {str(e)}")
            LOGGER.error(traceback.format_exc())
            
            # 尝试保存检查点
            try:
                checkpoint_path = os.path.join(self._output_dir, "error_checkpoint.rmc")
                self.save(path=checkpoint_path)
                LOGGER.info(f"错误检查点已保存: {checkpoint_path}")
            except:
                pass
                
            return False
    
    def export_analysis_plots(self, step_label):
        """
        导出当前分析图表
        
        :Parameters:
            #. step_label (str): 步骤标签用于文件命名
        """
        try:
            # 导出G(r)对比图
            gr_path = os.path.join(self._output_dir, f"gr_plot_{step_label}.png")
            self.pdf_constraint.plot(filename=gr_path)
            
            # 导出键比例对比图
            bond_path = os.path.join(self._output_dir, f"bond_ratio_plot_{step_label}.png")
            self.bond_constraint.plot(filename=bond_path)
            
            LOGGER.info(f"分析图表已导出: G(r)={gr_path}, 键比例={bond_path}")
            
        except Exception as e:
            LOGGER.error(f"导出分析图表失败: {str(e)}")
    
    def export_gr(self, path=None):
        """
        导出当前G(r)数据
        
        :Parameters:
            #. path (str): 输出文件路径
        """
        try:
            if path is None:
                path = os.path.join(self._output_dir, "current_gr.txt")
                
            # 获取实验G(r)数据
            exp_data = self.pdf_constraint.experimentalData
            exp_r = exp_data[:, 0]  
            exp_gr = exp_data[:, 1]
            
            # 计算当前G(r)
            constraint_value = self.pdf_constraint.get_constraint_value()
            sim_gr = constraint_value["total"]
            
            # 计算差异
            diff_gr = sim_gr - exp_gr
            
            # 写入文件
            with open(path, 'w') as f:
                f.write("# r experimental_gr simulated_gr difference\n")
                for i, r in enumerate(exp_r):
                    f.write(f"{r:.4f} {exp_gr[i]:.8e} {sim_gr[i]:.8e} {diff_gr[i]:.8e}\n")
                    
            return path
            
        except Exception as e:
            LOGGER.error(f"导出G(r)数据失败: {str(e)}")
            return None

    def _are_constraints_initialized(self):
        """检查约束是否已初始化。"""
        if not self.constraints:
            return False
            
        for constraint in self.constraints:
            if constraint.data is None:
                return False
                
        return True
        
    def should_step_get_rejected(self, realIndexes, relativeIndexes):
        """
        实现基于结构约束的接受准则。
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
            
        :Returns:
            #. rejected (bool): 如果应该拒绝步骤则返回True
        """
        try:
            # 获取移动原子信息
            moved_elem = self.allElements[realIndexes[0]]
            is_hydrogen = (moved_elem == 'H')
            
            # 获取移动前后的PDF误差
            gr_before = self.pdf_constraint.standardError  # 移动前误差
            gr_after = self.pdf_constraint.afterMoveStandardError  # 移动后误差
            delta_gr = gr_after - gr_before  # 误差变化
            gr_0 = self.pdf_constraint.initialError  # 初始误差 χ²_PDF,0
        
            # 获取移动前后的键比例误差
            bond_before = self.bond_constraint.standardError  # 移动前误差
            bond_after = self.bond_constraint.afterMoveStandardError  # 移动后误差
            delta_bond = bond_after - bond_before  # 误差变化
            bond_0 = self.bond_constraint.initialError  # 初始误差 χ²_bond,0
        
            # 防止除零
            eps = 1e-4
            gr_0 = max(gr_0, eps)
            bond_0 = max(bond_0, eps)

            delta_gr_norm = delta_gr / gr_0
            delta_bond_norm = delta_bond / bond_0

            # 组合两种约束的误差变化，使用权重
            if is_hydrogen:
                # 对于氢原子，主要考虑键比例约束
                delta_total = self.w2 * delta_bond
            else:
                # 对于非氢原子，按指定权重组合
                delta_total = self.w1 * delta_gr + self.w2 * delta_bond
            
            # 计算当前温度下的接受概率
            T_chi = self.T_chi_0 * (self.cooling_rate ** self.current_cycle)
                                        
            # 计算接受概率 P = min[1, exp(-delta/T_chi)]
            ratio = -delta_total / T_chi
            ratio_clipped = np.clip(ratio,-700,700)
            prob = min(1.0,np.exp(ratio_clipped)) 
            
            # 记录数据
            self._write_data('HRMC_acceptance.txt', 
                f"{self._tried} {prob:.6f} {delta_gr:.6f} {delta_bond:.6f} {T_chi:.6f}")
            self._write_data('HRMC_bond_ratios.txt',
                f"{self._tried} {self.bond_constraint.standardError:.6f} {delta_bond:.6f}")
            self._write_data('HRMC_gr.txt',
                f"{self._tried} {self.pdf_constraint.standardError:.6f} {delta_gr:.6f}")
            
            # 做出接受决定
            accepted = np.random.random() <= prob
            return not accepted
                
        except Exception as e:
            LOGGER.error(f"HRMC接受判断出错: {str(e)}")
            LOGGER.error(traceback.format_exc())
            return True  # 出错时默认拒绝
    
    def generate_summary_report(self, filename=None):
        """生成模拟结果摘要报告"""
        try:
            if filename is None:
                filename = os.path.join(self._output_dir, "simulation_report.txt")
                
            with open(filename, 'w') as f:
                f.write("# BondRatioHRMC 模拟结果报告\n")
                f.write(f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"# 生成用户: kang311\n")
                f.write(f"# 生成日期: {datetime.now().strftime('%Y-%m-%d')}\n\n")
                
                f.write("## 模拟参数\n")
                f.write(f"初始温度: {self.initial_temperature} K\n")
                f.write(f"最终温度: {self.current_temperature} K\n")
                f.write(f"权重(G(r)/键比例): {self.w1:.2f}/{self.w2:.2f}\n")
                f.write(f"原子数: {len(self.allElements)}\n\n")
                
                f.write("## 结果统计\n")
                acceptance_rate = self._accepted / max(1, self._tried) * 100
                f.write(f"总步数: {self._tried}\n")
                f.write(f"接受步数: {self._accepted}\n")
                f.write(f"接受率: {acceptance_rate:.2f}%\n\n")
                
                f.write("## 约束误差\n")
                f.write(f"G(r)误差: {self.pdf_constraint.standardError:.6f}\n")
                f.write(f"键比例误差: {self.bond_constraint.standardError:.6f}\n\n")
                
                # 获取键比例详细信息
                if hasattr(self.bond_constraint, 'get_constraint_value'):
                    bond_data = self.bond_constraint.get_constraint_value()
                    target_ratios = self.bond_constraint.targetBondRatios
                    
                    f.write("## 键比例详情\n")
                    f.write("键类型       | 目标比例 | 当前比例 | 差异\n")
                    f.write("------------|---------|---------|----------\n")
                    
                    # 获取键计数
                    bond_counts, total_bonds = None, 0
                    if hasattr(self, 'get_bond_counts'):
                        bond_counts, total_bonds = self.get_bond_counts()
                    
                    for bond_type, target in target_ratios.items():
                        current = bond_data.get(bond_type, 0.0)
                        diff = current - target
                        diff_percent = (diff / max(target, 1e-10)) * 100
                        
                        # 如果有键计数，显示计数信息
                        count_info = ""
                        if bond_counts and bond_type in bond_counts:
                            count = bond_counts[bond_type]
                            percent = (count / total_bonds * 100) if total_bonds > 0 else 0
                            count_info = f" ({count}/{total_bonds}, {percent:.1f}%)"
                            
                        f.write(f"{bond_type:12} | {target:.5f} | {current:.5f} | {diff:+.5f} ({diff_percent:+.1f}%){count_info}\n")
                    
                    # 计算符合率
                    if hasattr(self.bond_constraint, 'tolerances'):
                        tolerances = self.bond_constraint.tolerances
                        compliant_count = 0
                        for bond_type, target in target_ratios.items():
                            current = bond_data.get(bond_type, 0.0)
                            tolerance = tolerances.get(bond_type, 0.05)
                            if abs(current - target) <= tolerance:
                                compliant_count += 1
                                
                        compliance_rate = compliant_count / len(target_ratios) * 100
                        f.write(f"\n符合率: {compliant_count}/{len(target_ratios)} ({compliance_rate:.1f}%)\n")
                
            LOGGER.info(f"报告已生成: {filename}")
            return filename
            
        except Exception as e:
            LOGGER.error(f"生成报告失败: {str(e)}")
            return None
            
    def get_bond_counts(self):
        """
        获取各键类型计数和总键数
        
        :Returns:
            #. bondCounts (dict): 键类型及其计数
            #. totalBonds (int): 键总数
        """
        try:
            # 若缓存有效则直接返回
            if hasattr(self, '_lastBondCounts') and hasattr(self, '_lastTotalBonds'):
                if self._lastBondCounts is not None and self._lastTotalBonds is not None:
                    return self._lastBondCounts, self._lastTotalBonds
                    
            # 若没有有效缓存，则从约束中获取
            if self.bond_constraint and hasattr(self.bond_constraint, "get_bond_counts"):
                bond_counts, total_bonds = self.bond_constraint.get_bond_counts()
                # 缓存结果
                self._lastBondCounts = bond_counts
                self._lastTotalBonds = total_bonds
                return bond_counts, total_bonds
            else:
                LOGGER.warning("无法从约束中获取键计数")
                return {}, 0
                
        except Exception as e:
            LOGGER.error(f"获取键计数出错: {str(e)}")
            return {}, 0
