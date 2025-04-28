"""
GRHRMCEngine implements a structure-driven HRMC model focusing solely on PDF constraints.

This engine focuses on structural optimization using Pair Distribution Function (G(r)) fitting
without temperature annealing and atom swaps.
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
from fullrmc.Constraints.DistanceConstraints import InterMolecularDistanceConstraint
from fullrmc.Selectors.RandomSelectors import RandomSelector
from fullrmc.Core.GroupSelector import RecursiveGroupSelector

class GRHRMCEngine(Engine):
    """
    HRMC引擎专注于PDF约束的结构优化。
    该引擎使用Metropolis蒙特卡洛算法通过G(r)拟合驱动模拟，不使用温度退火。
    
    :Parameters:
        #. path (str): 引擎保存路径
        #. freshStart (bool): 是否重新开始模拟
    """
    def __init__(self, path, freshStart=False):
        # 初始化父类
        super(GRHRMCEngine, self).__init__(path=path, freshStart=freshStart)
        
        # 设置基本参数
        self._chi_parameter = 0.1  # 固定参数，代替之前的温度参数
        self._current_cycle = 1    # 当前周期
        
        # 约束引用
        self.pdf_constraint = None
        self.distance_constraint = None
        
        # 初始化输出目录和文件
        self._output_dir = os.path.dirname(path) if path else None 
        self._files = {}
        
        # 设置输出文件（如果路径已提供）
        if self._output_dir:
            self.setup_output_files(self._output_dir)
            
        # 缓存数据
        self._cache = {
            # 计算结果缓存
            'gr': {},       # G(r)计算结果缓存
            'sij': {},      # S(Q)计算结果缓存
            'atom_ff': {},  # 原子形式因子缓存
            'lorch': {},    # Lorch窗函数缓存
            'shape_function': {}  # 形状函数缓存
        }
        
        # 形状函数参数
        self._shape_function_params = {
            'rmin': 0.0, 
            'rmax': None,
            'dr': 0.1,
            'qmin': 0.0001, 
            'qmax': 30.0, 
            'dq': 0.01,
            'updateFreq': 1000
        }
        
        # 标准误差历史
        self.standardErrorHistory = []
        
        # 记录日志
        LOGGER.info(f"GRHRMCEngine 初始化完成，输出目录: {self._output_dir}")

    def increment_accepted(self):
        """增加接受计数器"""
        self._Engine__accepted += 1

    def increment_tried(self):
        """增加尝试计数器"""
        self._Engine__tried += 1
    
    def __getstate__(self):
        """自定义序列化，明确区分哪些是需要持久化的状态。"""
        state = super(GRHRMCEngine, self).__getstate__()
        
        # 确保历史误差被保存
        state['standardErrorHistory'] = getattr(self, 'standardErrorHistory', [])
        
        # 清除不需要持久化的临时缓存数据
        if '_cache' in state:
            # 保留形状函数参数
            shape_function_params = state['_cache'].get('shape_function', {})
            state['_cache'] = {
                'shape_function': shape_function_params
            }
        
        return state
        
    def __setstate__(self, state):
        """处理反序列化，确保所有必要的状态都被恢复。"""
        self.__dict__.update(state)
        
        # 初始化缓存结构，确保所有必要的缓存键存在
        if not hasattr(self, '_cache'):
            self._cache = {
                'gr': {},
                'sij': {},
                'atom_ff': {},
                'lorch': {},
                'shape_function': {}
            }
        elif 'gr' not in self._cache:
            self._cache.update({
                'gr': {},
                'sij': {},
                'atom_ff': {},
                'lorch': {}
            })
            
        # 确保标准误差历史存在
        if not hasattr(self, 'standardErrorHistory'):
            self.standardErrorHistory = []

    def cleanup(self):
        """清理资源。"""
        # 调用父类的清理方法（如果存在）
        if hasattr(super(GRHRMCEngine, self), 'cleanup'):
            super(GRHRMCEngine, self).cleanup()
        
        # 关闭可能打开的文件
        for filename in self._files.keys():
            filepath = os.path.join(self._output_dir, filename) if self._output_dir else None
            if filepath and os.path.exists(filepath):
                LOGGER.debug(f"清理资源: {filepath}")
        
        # 清理缓存数据
        self._cache['gr'] = {}
        self._cache['sij'] = {}
        self._cache['atom_ff'] = {}
        
        LOGGER.info("GRHRMCEngine资源已清理")
            
    def setup_output_files(self, output_dir):
        """设置输出目录和文件。"""
        self._output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # 定义文件和表头
        self._files = {
            'HRMC_acceptance.txt': "# Step  Element  Probability  DeltaGr  Chi",
            'HRMC_gr.txt': "# Step  gr_new  GrDiff",
        }

        # 创建带表头的文件
        for filename, header in self._files.items():
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                f.write(header + '\n')
                
        LOGGER.info(f"输出文件已设置在 {output_dir}")

    def _write_data(self, filename, data):
        """使用上下文管理器将数据写入文件。"""
        if self._output_dir is None:
            LOGGER.error("输出目录未设置")
            return
            
        filepath = os.path.join(self._output_dir, filename)
        
        # 验证文件是否存在
        if not os.path.exists(filepath):
            with open(filepath, 'w') as f:
                f.write(self._files.get(filename, "# Data") + '\n')
                
        # 附加数据
        with open(filepath, 'a') as f:
            f.write(data + '\n')
            
    def add_constraints(self, constraints):
        """
        添加并正确初始化约束。
        
        :Parameters:
            #. constraints (list): 约束对象列表
        """
        # 在存储引用后调用父类的add_constraints方法
        super(GRHRMCEngine, self).add_constraints(constraints)
        
        # 保存约束引用
        for constraint in constraints:
            if isinstance(constraint, PairDistributionConstraint):
                self.pdf_constraint = constraint
                LOGGER.info(f"已添加PDF约束: {constraint.__class__.__name__}")
            elif isinstance(constraint, InterMolecularDistanceConstraint):
                self.distance_constraint = constraint
                LOGGER.info(f"已添加分子间距离约束: {constraint.__class__.__name__}")
                
    def initialize_used_constraints(self, force=False, sortConstraints=False):
        """安全地初始化约束。"""
        # 初始化基本约束
        usedConstraints, constraints, rigidConstraints = super(GRHRMCEngine, self).initialize_used_constraints(force=force, sortConstraints=sortConstraints)
        
        # 验证我们需要的特定约束
        if not hasattr(self, 'pdf_constraint') or self.pdf_constraint is None:
            for constraint in constraints:
                if isinstance(constraint, PairDistributionConstraint):
                    self.pdf_constraint = constraint
                    break
            
            if self.pdf_constraint is None:
                LOGGER.error("PDF约束未找到")
                raise ValueError("PDF约束未找到，无法继续")
                
        # 检查分子间距离约束
        if not hasattr(self, 'distance_constraint') or self.distance_constraint is None:
            for constraint in constraints:
                if isinstance(constraint, InterMolecularDistanceConstraint):
                    self.distance_constraint = constraint
                    break
                    
        # 输出找到的约束
        LOGGER.info(f"已初始化PDF约束: {self.pdf_constraint.__class__.__name__ if self.pdf_constraint else 'None'}")
        LOGGER.info(f"已初始化分子间距离约束: {self.distance_constraint.__class__.__name__ if self.distance_constraint else 'None'}")
        
        # 计算约束系数
        elements = self.allElements
        N_total = len(elements)
        N_nonH = sum(1 for elem in elements if elem != 'H')
        
        if N_nonH == 0:
            LOGGER.error("结构中必须包含非氢原子")
            raise ValueError("结构中必须包含非氢原子")
                              
        LOGGER.info(f"结构加载: 总原子数={N_total}, 非氢原子数={N_nonH}")
        
        return usedConstraints, constraints, rigidConstraints
    
    def set_shape_function_parameters(self, params):
        """
        设置形状函数参数。
        
        :Parameters:
            #. params (dict): 形状函数参数字典
        """
        if not isinstance(params, dict):
            raise ValueError("参数必须是字典")
            
        # 更新形状函数参数
        self._shape_function_params.update(params)
        
        # 清除旧的缓存
        self._cache['shape_function'] = {}
        
        LOGGER.info(f"形状函数参数已更新: {self._shape_function_params}")
        return True
    
    def verify_constraints_data_updated(self):
        """验证约束数据是否已正确更新。"""
        for c in self.constraints:
            if hasattr(c, 'data') and c.data is None:
                LOGGER.warning(f"约束 {c.__class__.__name__} 的数据为None，尝试重新计算...")
                if hasattr(c, 'compute_data'):
                    c.compute_data(update=True)
            
    def run(self, numberOfSteps, saveFrequency=10000):
        """
        运行HRMC步骤。
        
        :Parameters:
            #. numberOfSteps (int): 模拟步数
            #. saveFrequency (int): 保存频率
        """
        # 记录开始时间
        start_time = time.time()
        LOGGER.info(f"开始HRMC模拟，步数: {numberOfSteps}")
        
        # 确保约束已初始化
        usedConstraints, constraints, rigidConstraints = self.initialize_used_constraints()
        
        # 计算初始误差值
        if not hasattr(self.pdf_constraint, 'initialError') or self.pdf_constraint.initialError is None:
            self.pdf_constraint.initialError = self.pdf_constraint.standardError
        
        LOGGER.info(f"初始PDF误差: {self.pdf_constraint.initialError:.6f}")
                
        # 记录当前步数
        current_step = 0
        
        # 统计局部接受率
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
            element_counts = {}
            for idx in group.indexes:
                elem = self.allElements[idx]
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
                    
            self.increment_tried()
                
            # 检查约束数据一致性（每500步）
            if step % 500 == 0:
                self.verify_constraints_data_updated()

            # 定期保存结果
            if step % saveFrequency == 0 or step == numberOfSteps:
                # 计算接受率
                local_rate = local_accepted / max(1, local_tried)
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
                
                # 添加到历史记录
                self.standardErrorHistory.append(pdf_error)
                
                # 计算改进率
                pdf_improvement = (self.pdf_constraint.initialError - pdf_error) / self.pdf_constraint.initialError * 100
                
                LOGGER.info(f"步骤 {step}/{numberOfSteps} ({step/numberOfSteps*100:.1f}%)")
                LOGGER.info(f"- 接受率: 局部={local_rate:.2%}, 全局={global_rate:.2%}")
                LOGGER.info(f"- PDF误差: {pdf_error:.6f} (改进: {pdf_improvement:.1f}%)")
                LOGGER.info(f"- 已用时间: {elapsed:.1f}秒, 剩余: {remaining:.1f}秒")
                
                # 导出当前PDB
                pdb_path = os.path.join(self._output_dir, f"structure_step_{step}.pdb")
                self.export_pdb(path=pdb_path)
                
                # 导出当前G(r)图
                self.export_analysis_plots(step)
                
                # 清理缓存(保留最近结果)
                self._clear_cache(keep_latest=True)
        
        # 完成模拟
        total_time = time.time() - start_time
        overall_rate = self.accepted / max(1, self.tried)
        
        LOGGER.info(f"模拟完成，总步数: {numberOfSteps}，总用时: {total_time:.1f}秒")
        LOGGER.info(f"总体接受率: {overall_rate:.2%} ({self.accepted}/{self.tried})")
        LOGGER.info(f"最终PDF误差: {self.pdf_constraint.standardError:.6f}")
        
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

    def _clear_cache(self, keep_latest=False):
        """
        清理计算缓存，释放内存
        
        :Parameters:
            #. keep_latest (bool): 是否保留最新计算结果
        """
        if keep_latest and self._cache.get('gr') and 'latest' in self._cache['gr']:
            # 只保留最新结果
            latest = self._cache['gr'].get('latest')
            self._cache['gr'] = {'latest': latest}
        else:
            # 清空所有G(r)缓存
            self._cache['gr'] = {}
            
        # 清理S(Q)缓存
        self._cache['sij'] = {}
    
    def export_analysis_plots(self, step_label):
        """
        导出当前分析图表
        
        :Parameters:
            #. step_label (str): 步骤标签用于文件命名
        """
        # 导出G(r)对比图
        gr_path = os.path.join(self._output_dir, f"gr_plot_{step_label}.png")
        
        self.plot_gr_comparison(save_path=gr_path)
        
        LOGGER.info(f"分析图表已导出: G(r)={gr_path}")
    
    def export_gr(self, path=None):
        """
        导出当前G(r)数据
        
        :Parameters:
            #. path (str): 输出文件路径
        """
        if path is None:
            path = os.path.join(self._output_dir, "current_gr.txt")

        # 获取实验G(r)数据
        exp_data = self.pdf_constraint.experimentalData
        exp_r = exp_data[:, 0]  
        exp_gr = exp_data[:, 1] 
      
        # 使用增量计算获取当前G(r)
        sim_gr, diff_key = self._calculate_cached_gr()
    
        if sim_gr is None:
            LOGGER.error("无法导出G(r)数据：计算失败")
            return None
        
        # 计算差异
        diff_gr = sim_gr - exp_gr
        
        # 写入文件
        with open(path, 'w') as f:
            f.write("# r experimental_gr simulated_gr difference\n")
            for i, r in enumerate(exp_r):
                f.write(f"{r:.4f} {exp_gr[i]:.8e} {sim_gr[i]:.8e} {diff_gr[i]:.8e}\n")
                
        return path

    def _are_constraints_initialized(self):
        """检查约束是否已初始化。"""
        if not self.constraints:
            return False
            
        for constraint in self.constraints:
            if constraint.data is None:
                return False
                
        return True

    def _compute_atom_form_factor(self, elem, q_array):
        """
        计算原子散射因子，包括缓存机制
        
        :Parameters:
            #. elem (str): 元素符号
            #. q_array (numpy.ndarray): Q值数组
        
        :Returns:
            #. f_q (numpy.ndarray): 散射因子值数组
        """
        # Waasmaier-Kirfel参数 - 用于G(r)计算
        _wk_params = {
            "C": ([2.657506, 1.078079, 1.490909, -4.241070, 0.713791],
                 [14.780758, 0.776775, 42.086842, -0.000294, 0.239535],
                 4.297983),
            "H": ([0.413048, 0.294953, 0.187491, 0.080701, 0.023736],
                 [15.569946, 32.398468, 5.711404, 61.889874, 1.334118],
                 0.000049),
            "O": ([2.960427, 2.508818, 0.637853, 0.722838, 1.142756],
                 [14.182259, 5.936858, 0.112726, 34.958481, 0.390240],
                 0.027014),
            "N": ([11.893780, 3.277479, 1.858092, 0.858927, 0.912985],
                 [0.000158, 10.232723, 30.344690, 0.656065, 0.217287],
                 -11.804902),
            "S": ([6.372157, 5.154568, 1.473732, 1.635073, 1.209372],
                 [1.514347, 22.092527, 0.061373, 55.445175, 0.646925],
                 0.154722)
        }
        
        # 生成缓存键
        import hashlib
        q_key = hashlib.md5(q_array.tobytes()).hexdigest()
        cache_key = f"{elem}_{q_key}"
        
        # 检查缓存
        if cache_key in self._cache['atom_ff']:
            return self._cache['atom_ff'][cache_key]
            
        # 不在缓存中，计算散射因子
        if elem in _wk_params:
            a, b, c = _wk_params[elem]
            
            # 使用NumPy向量化操作计算
            f_q = np.zeros_like(q_array)
            q_squared = (q_array / (4.0 * np.pi))**2
            
            for i in range(len(a)):
                f_q += a[i] * np.exp(-b[i] * q_squared)
            f_q += c
        else:
            # 如果没有参数，使用默认值
            f_q = np.ones_like(q_array)
            LOGGER.warning(f"未找到元素 {elem} 的形式因子参数，使用默认值")
            
        # 存入缓存
        self._cache['atom_ff'][cache_key] = f_q
        
        return f_q

    def _calculate_cached_gr(self):
        """
        使用缓存机制计算G(r)，避免重复计算
        
        :Returns:
            #. new_gr (numpy.ndarray): 计算得到的G(r)
            #. cache_key (str): 缓存键
        """
        # 获取数据状态作为缓存键
        import hashlib
        data_state = {
            'intra_hash': hashlib.md5(self.pdf_constraint.data["intra"].tobytes()).hexdigest(),
            'inter_hash': hashlib.md5(self.pdf_constraint.data["inter"].tobytes()).hexdigest(),
            'elements': tuple(self.elements),
            'num_atoms': self.numberOfAtoms
        }
        cache_key = str(hash(frozenset(data_state.items())))
        
        # 检查缓存中是否有结果
        if cache_key in self._cache['gr']:
            LOGGER.debug("使用缓存的G(r)计算结果")
            return self._cache['gr'][cache_key], cache_key
            
        # 计算G(r)
        new_gr = self._calculate_gr_with_shape_function()
        
        # 存入缓存
        self._cache['gr'][cache_key] = new_gr
        self._cache['gr']['latest'] = new_gr  # 同时保存最新结果
        
        return new_gr, cache_key
        
    def _calculate_gr_with_shape_function(self):
        """
        使用形状函数计算G(r)
        
        :Returns:
            #. new_gr (numpy.ndarray): 计算得到的G(r)
        """
        # 获取实验G(r)数据，仅用于确定输出数组形状
        exp_data = self.pdf_constraint.experimentalData
        exp_r = exp_data[:, 0]  
        exp_gr = exp_data[:, 1]
        
        # 获取必要的参数
        rho0 = len(self.allElements) / self.volume  # 系统数密度
        r_array = self.pdf_constraint.shellCenters  # r值数组
        shell_volumes = self.pdf_constraint.shellVolumes  # 壳体积
        total_volume = self.volume  # 系统总体积

        # 从形状函数参数中获取Q值设置
        qmin = self._shape_function_params.get('qmin', 0.0001)
        qmax = self._shape_function_params.get('qmax', 30.0)
        dq = self._shape_function_params.get('dq', 0.01)
        
        # 创建Q值数组用于傅里叶变换
        q_array = np.arange(qmin, qmax, dq)

        # 创建元素浓度字典和原子形式因子字典
        c_dict = {}  # 元素浓度字典
        f_dict = {}  # 原子散射因子字典

        # 计算元素浓度
        for elem in set(self.allElements):
            c_dict[elem] = self.numberOfAtomsPerElement[elem] / len(self.allElements)

        # 计算所有元素的散射因子(使用缓存)
        for elem in c_dict.keys():
            f_dict[elem] = self._compute_atom_form_factor(elem, q_array)

        # 获取 lorch 窗函数（如果没缓存就计算）
        if qmax not in self._cache['lorch']:
            self._cache['lorch'][qmax] = np.sinc(q_array / qmax)

        lorch_window = self._cache['lorch'][qmax]

        # 初始化总S(Q)和散射因子平方和
        S_q_total = np.zeros_like(q_array)
        f_avg_sq = np.zeros_like(q_array)
        
        # 计算分母 <f>^2 (使用NumPy向量化)
        for elem_i, c_i in c_dict.items():
            f_i = f_dict[elem_i]
            for elem_j, c_j in c_dict.items():
                f_j = f_dict[elem_j]
                f_avg_sq += c_i * c_j * f_i * f_j
        
        # 计算所有元素对的g_ij(r)和S_ij(Q)
        for pair in self.pdf_constraint.elementsPairs:
            elem_i, elem_j = pair

            # 获取元素索引
            idi = self.elements.index(elem_i)
            idj = self.elements.index(elem_j)
    
            # 获取原子数
            ni = self.numberOfAtomsPerElement[elem_i]
            nj = self.numberOfAtomsPerElement[elem_j]

            # 获取直方图数据
            if idi == idj:
                nij_hist = self.pdf_constraint.data["intra"][idi, idj, :] + self.pdf_constraint.data["inter"][idi, idj, :]
                normalization = shell_volumes * ni * (ni-1) / total_volume
            else:
                nij_hist = (self.pdf_constraint.data["intra"][idi, idj, :] + 
                            self.pdf_constraint.data["intra"][idj, idi, :] + 
                            self.pdf_constraint.data["inter"][idi, idj, :] + 
                            self.pdf_constraint.data["inter"][idj, idi, :])
                normalization = shell_volumes * ni * nj / total_volume

            # 计算g_ij(r) - 使用NumPy向量化
            mask = normalization > 0
            gij = np.zeros_like(nij_hist, dtype=FLOAT_TYPE)
            gij[mask] = nij_hist[mask] / normalization[mask]
            
            # 计算S_ij(Q) - 使用矢量化积分
            sij_q = np.zeros((len(q_array),), dtype=FLOAT_TYPE)
            
            # 为每个q值计算积分
            for i, qi in enumerate(q_array):
                # 使用矢量化计算
                integrand = r_array * (gij - 1.0) * np.sin(qi * r_array)
                # 使用梯形积分
                dr = r_array[1] - r_array[0]
                integral = np.trapz(integrand, dx=dr)
                
                sij_q[i] = 1.0 + 4.0 * np.pi * rho0 / qi * integral
                
            # 权重计算
            c_i = c_dict[elem_i]
            c_j = c_dict[elem_j]
            f_i = f_dict[elem_i]
            f_j = f_dict[elem_j]
            
            # 将加权S_ij加到总S(Q)
            weight = c_i * c_j * f_i * f_j
            S_q_total += weight * sij_q
        
        # 计算总S(Q)
        s_q_normalized = S_q_total / f_avg_sq

        # === 应用 Lorch 窗函数和形状函数 ===
        lorch_window = np.sinc(q_array / qmax)  # sinc(x) = sin(pi x)/(pi x)
        
        # 应用形状函数（根据实验装置有限尺寸的修正）
        shape_function = self._compute_shape_function(q_array)
        
        # 应用窗函数和形状函数
        fq = q_array * (s_q_normalized - 1.0) * lorch_window * shape_function
        
        # 通过傅里叶变换从S(Q)计算G(r)
        gr = np.zeros((len(r_array),), dtype=FLOAT_TYPE)
        
        # 使用矢量化计算傅里叶变换
        for i, ri in enumerate(r_array):
            # 计算积分
            integrand = fq * np.sin(q_array * ri)
            integral = np.trapz(integrand, dx=dq)
            
            gr[i] = 2.0 / np.pi * integral
        
        # 确保结果长度与实验数据匹配
        if len(gr) != len(exp_gr):
            # 使用简单插值
            from scipy.interpolate import interp1d
            interp_func = interp1d(r_array, gr, bounds_error=False, fill_value="extrapolate")
            gr = interp_func(exp_r)
        
        return gr
        
    def _compute_shape_function(self, q_array):
        """
        计算形状函数
        
        :Parameters:
            #. q_array (numpy.ndarray): Q值数组
            
        :Returns:
            #. shape_function (numpy.ndarray): 形状函数数组
        """
        # 为相同的q_array检查缓存
        import hashlib
        q_key = hashlib.md5(q_array.tobytes()).hexdigest()
        
        if q_key in self._cache['shape_function']:
            return self._cache['shape_function'][q_key]
        
        # 计算高斯形状函数
        # 这是一个简化的形状函数，实际应用中可能需要更复杂的函数
        sigma = 0.5  # 形状函数的宽度参数
        shape_function = np.exp(-(q_array**2)/(2*sigma**2))
        
        # 存入缓存
        self._cache['shape_function'][q_key] = shape_function
        
        return shape_function

    def should_step_get_rejected(self, realIndexes, relativeIndexes):
        """
        实现两阶段约束接受准则：先检查刚性约束，再使用G(r)约束计算Metropolis接受概率。
    
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
        
        :Returns:
            #. rejected (bool): 如果应该拒绝步骤则返回True
        """
        # 第一阶段：检查分子间距离约束（刚性约束）
        if self.distance_constraint:
            # 检查距离约束是否被违反 - 使用activeAtomsData来判断
            if (self.distance_constraint.activeAtomsDataAfterMove is not None and 
                self.distance_constraint.activeAtomsDataBeforeMove is not None):
            
                # 检查数量是否增加 - 表示有更多原子违反了约束
                after_number = self.distance_constraint.activeAtomsDataAfterMove.get("number", 0)
                before_number = self.distance_constraint.activeAtomsDataBeforeMove.get("number", 0)
            
                # 如果违反约束的数量增加，拒绝移动
                if np.any(after_number > before_number):
                    # 记录违反刚性约束的情况
                    moved_elem = self.allElements[realIndexes[0]]
                    step = getattr(self, 'tried', 0)
                    self._write_data('HRMC_rejection.txt', 
                        f"{step}    {moved_elem}    REJECTED_RIGID_CONSTRAINT")
                
                    # 违反了刚性约束，立即拒绝
                    return True
    
        # 第二阶段：基于G(r)约束的Metropolis接受准则
        # 获取移动原子信息
        moved_elem = self.allElements[realIndexes[0]]

        # 计算新旧构型的G(r)差异
        # 获取当前G(r)，使用缓存机制
        old_gr, old_key = self._calculate_cached_gr()
    
        # 获取实验G(r)数据
        exp_data = self.pdf_constraint.experimentalData
        exp_gr = exp_data[:, 1]
    
        # 使用增量方式计算移动后的数据
        dataIntraBefore = self.pdf_constraint.data["intra"]
        dataInterBefore = self.pdf_constraint.data["inter"]

        # 使用增量方式计算移动后的数据
        dataIntraAfter = dataIntraBefore - self.pdf_constraint.activeAtomsDataBeforeMove["intra"] + \
                      self.pdf_constraint.activeAtomsDataAfterMove["intra"]
        dataInterAfter = dataInterBefore - self.pdf_constraint.activeAtomsDataBeforeMove["inter"] + \
                      self.pdf_constraint.activeAtomsDataAfterMove["inter"]
    
        # 临时保存原始数据
        temp_intra = self.pdf_constraint.data["intra"].copy()
        temp_inter = self.pdf_constraint.data["inter"].copy()
    
        # 临时替换数据以计算新G(r)
        self.pdf_constraint.data["intra"] = dataIntraAfter
        self.pdf_constraint.data["inter"] = dataInterAfter
    
        # 计算新G(r)
        new_gr, new_key = self._calculate_cached_gr()
    
        # 恢复原始数据
        self.pdf_constraint.data["intra"] = temp_intra
        self.pdf_constraint.data["inter"] = temp_inter
    
        # 计算χ² = Σ(sim-exp)²/Σ(exp²)
        def compute_chi_gr(sim_gr, exp_gr):
            diff_gr = (sim_gr - exp_gr)**2  
            exp_gr_sum = np.sum(exp_gr**2)
            return np.sum(diff_gr) / exp_gr_sum

        # 计算grΔχ² = χ_j² - χ_i²
        chi_gr_new = compute_chi_gr(new_gr, exp_gr)
        chi_gr_old = compute_chi_gr(old_gr, exp_gr)
        delta_gr = chi_gr_new - chi_gr_old
    
        # 仅使用G(r)约束的误差变化来计算接受概率
        delta_total = delta_gr
    
        # 固定的"温度"参数，控制接受概率
        chi = self._chi_parameter
                            
        # 计算接受概率 P = min[1, exp(-delta_total / chi)]
        prob = min(1, np.exp(-delta_total / chi))
    
        # 记录数据
        step = getattr(self, 'tried', 0)  # 使用父类计数器
        self._write_data('HRMC_acceptance.txt', 
            f"{step}    {moved_elem}    {prob:.6f}   {delta_gr:.6f}    {chi:.6f}")
        self._write_data('HRMC_gr.txt',
            f"{step} {chi_gr_new:.6f}  {delta_gr:.6f}")
    
        # 做出接受决定
        accepted = np.random.random() <= prob
    
        return not accepted

    def generate_summary_report(self, filename=None):
        """生成模拟结果摘要报告"""
        if filename is None:
            filename = os.path.join(self._output_dir, "simulation_report.txt")
            
        with open(filename, 'w') as f:
            f.write("# GRHRMCEngine 模拟结果报告\n")
            f.write(f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# 生成用户: kang311\n")
            f.write(f"# 生成日期: {datetime.now().strftime('%Y-%m-%d')}\n\n")
            
            f.write("## 模拟参数\n")
            f.write(f"固定Chi值: {self._chi_parameter}\n")
            f.write(f"原子数: {len(self.allElements)}\n")
            f.write(f"是否使用PBC: {self.isPBC}\n\n")
            
            f.write("## 结果统计\n")
            acceptance_rate = self.accepted / max(1, self.tried) * 100
            f.write(f"总步数: {self.tried}\n")
            f.write(f"接受步数: {self.accepted}\n")
            f.write(f"接受率: {acceptance_rate:.2f}%\n\n")
            
            f.write("## 约束误差\n")
            f.write(f"G(r)误差: {self.pdf_constraint.standardError:.6f}\n")
            
            if self.distance_constraint:
                f.write(f"分子间距离约束误差: {self.distance_constraint.standardError:.6f}\n")
                
        LOGGER.info(f"报告已生成: {filename}")
        return filename

    def plot_gr_comparison(self, save_path=None):
        """
        绘制G(r)对比图
        
        :Parameters:
            #. save_path (str): 保存路径
        """
        try:
            import matplotlib.pyplot as plt
            
            # 使用缓存方法计算G(r)
            sim_gr, _ = self._calculate_cached_gr()
            
            # 获取实验数据
            exp_data = self.pdf_constraint.experimentalData
            exp_r = exp_data[:, 0]  
            exp_gr = exp_data[:, 1]

            # 绘图
            plt.figure(figsize=(10, 6))
            plt.plot(exp_r, exp_gr, label='Experimental G(r)', linestyle='--', color='blue')
            plt.plot(exp_r, sim_gr, label='Simulated G(r)', linestyle='-', color='red')
            plt.xlabel("r (Å)", fontsize=14)
            plt.ylabel("G(r)", fontsize=14)
            plt.legend(fontsize=12)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300)
                plt.close()
                LOGGER.info(f"G(r) 对比图保存至: {save_path}")
            else:
                plt.show()

        except Exception as e:
            LOGGER.error(f"生成G(r)对比图失败: {str(e)}")
            if 'plt' in locals():
                plt.close()
