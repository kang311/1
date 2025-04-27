"""
GRHRMCEngine implements a structure-driven HRMC model focusing solely on PDF constraints.

This engine focuses on structural optimization using only Pair Distribution Function (G(r)) fitting
without bond ratio constraints.
"""
from __future__ import print_function
import os
import numpy as np
import time
import traceback
import random
import hashlib
from datetime import datetime
from collections import defaultdict

from fullrmc.Engine import Engine 
from fullrmc.Globals import LOGGER, FLOAT_TYPE
from fullrmc.Core.boundary_conditions_collection import transform_coordinates
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Generators.Swaps import SwapPositionsGenerator

class GRHRMCEngine(Engine):
    """
    HRMC引擎专注于PDF约束的结构优化。
    该引擎不使用能量约束，而是仅通过G(r)特征驱动模拟。
    
    :Parameters:
        #. path (str): 引擎保存路径
        #. freshStart (bool): 是否重新开始模拟
    """
    def __init__(self, path, freshStart=False):
        # 初始化父类
        super(GRHRMCEngine, self).__init__(path=path, freshStart=freshStart)
        
        # 退火参数
        self._temperature = {
            'initial': 1000.0,    # 初始退火温度 (K)
            'current': 1000.0,    # 当前温度
            'final': 100.0,       # 最终退火温度 (K)
            'cooling_rate': 0.95, # 冷却率
            'chi_0': 0.1,         # 初始归一化温度
            'chi': 1.0            # 当前归一化温度
        }
        self._current_cycle = 1   # 当前周期
        
        # 约束引用
        self.pdf_constraint = None
        
        # 统计计数器 - 统一使用父类计数器
        
        # 初始化输出目录和文件
        self._output_dir = os.path.dirname(path) if path else None 
        self._files = {}
        
        # 设置输出文件（如果路径已提供）
        if self._output_dir:
            self.setup_output_files(self._output_dir)
            
        # 缓存数据 - 明确区分哪些是临时状态
        self._cache = {
            # 计算结果缓存
            'gr': {},       # G(r)计算结果缓存
            'sij': {},      # S(Q)计算结果缓存
            'atom_ff': {},  # 原子形式因子缓存
            'lorch': {},    
            # 持久化状态
            'simulation': {
                'current_stage': 0,  # 当前阶段
                'swap_count': 0      # 交换次数
            }
        }
        
        # 阶段控制设置
        self._stages = {
            'current': 0,                   # 当前阶段
            'steps': [],                    # 每个阶段的步数
            'swap_frequencies': [],         # 每个阶段的交换频率
            'temperature_update_freq': 1000 # 温度更新频率
        }
        
        # 原子交换生成器 - 统一管理
        self._swap_generators = {}           # 元素 -> 交换生成器的映射
        self._available_swap_elements = []   # 可用于交换的元素列表
        
        # Waasmaier-Kirfel参数 - 用于G(r)计算
        self._wk_params = {
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
        
        # 记录日志
        LOGGER.info(f"GRHRMCEngine 初始化完成，输出目录: {self._output_dir}")

    def increment_accepted(self):
        self._Engine__accepted += 1

    def increment_tried(self):
        self._Engine__tried += 1
    
    def __getstate__(self):
        """自定义序列化，明确区分哪些是需要持久化的状态。"""
        state = super(GRHRMCEngine, self).__getstate__()
        
        # 清除不需要持久化的临时缓存数据
        if '_cache' in state:
            # 只保留需要持久化的状态
            state['_cache'] = {
                'simulation': state['_cache'].get('simulation', {})
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
                'simulation': {
                    'current_stage': 0,
                    'swap_count': 0
                }
            }
        elif 'gr' not in self._cache:
            self._cache.update({
                'gr': {},
                'sij': {},
                'atom_ff': {}
            })
        
        # 确保交换生成器字典存在
        if not hasattr(self, '_swap_generators'):
            self._swap_generators = {}
        
        # 确保可用元素列表存在
        if not hasattr(self, '_available_swap_elements'):
            self._available_swap_elements = []
            
            
        # 确保阶段控制设置存在
        if not hasattr(self, '_stages'):
            self._stages = {
                'current': 0,
                'steps': [],
                'swap_frequencies': [],
                'temperature_update_freq': 1000
            }

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

        # 定义文件和表头 - 简化为只关注G(r)
        self._files = {
            'HRMC_acceptance.txt': "# Step Element Probability DeltaGr T_chi",
            'HRMC_gr.txt': "# Step StandardError GrDiff",
            'HRMC_swaps.txt': "# Step SwapType FromElement ToElement Success"  # 记录交换操作
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
        
        # 保存PDF约束引用
        for constraint in constraints:
            if isinstance(constraint, PairDistributionConstraint):
                self.pdf_constraint = constraint
                LOGGER.info(f"已添加PDF约束: {constraint.__class__.__name__}")
                
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
                
        # 计算温度转换系数
        elements = self.allElements
        N_total = len(elements)
        N_nonH = sum(1 for elem in elements if elem != 'H')
        
        if N_nonH == 0:
            LOGGER.error("结构中必须包含非氢原子")
            raise ValueError("结构中必须包含非氢原子")
                              
        LOGGER.info(f"结构加载: 总原子数={N_total}, 非氢原子数={N_nonH}")
        
        return usedConstraints, constraints, rigidConstraints
            
    def setup_atom_swaps(self):
        """
        设置原子交换生成器，统一管理交换机制。
        
        :Returns:
            #. success (bool): 是否成功设置交换生成器
        """
        # 获取所有原子元素
        allElements = self.allElements
        unique_elements = set(allElements)
    
        # 清空旧的交换生成器
        self._swap_generators.clear()
        self._available_swap_elements = []
    
        # 为每种元素创建交换生成器
        for elem in unique_elements:
            # 获取该元素的所有原子索引
            elem_indexes = [[idx] for idx in range(len(allElements)) if allElements[idx] == elem]
            
            # 只为有多个原子的元素创建交换生成器
            if len(elem_indexes) > 1:
                generator = SwapPositionsGenerator(swapList=elem_indexes)
                self._swap_generators[elem] = generator
                self._available_swap_elements.append(elem)
    
    
        # 重置交换计数
        self._cache['simulation']['swap_count'] = 0
    
        LOGGER.info(f"原子交换生成器已为以下元素设置: {self._available_swap_elements}")
        return len(self._available_swap_elements) > 0

    def setup_hrmc(self, initial_temp, final_temp, steps, 
                   swap_frequencies=None, temp_update_freq=1000,
                   stages=None, target_accept_rate=0.4):
        """
        设置HRMC参数，支持多阶段模拟配置。

        Parameters:
            initial_temp (float): 初始物理温度 (K)
            final_temp (float): 最终物理温度 (K)
            steps (int): 总模拟步数
            swap_frequencies (list, int): 每个阶段的交换频率或单一频率
            temp_update_freq (int): 温度更新频率
            stages (list): 多阶段模拟的配置列表，每个元素是包含步数的字典
            target_accept_rate (float): 目标接受率
        """
        # 设置温度参数
        self._temperature['initial'] = initial_temp
        self._temperature['final'] = final_temp
        self._temperature['current'] = initial_temp
        self._temperature['cooling_rate'] = 0.98
        
        # 设置温度更新频率
        self._stages['temperature_update_freq'] = temp_update_freq
        
        # 设置多阶段模拟
        if stages is not None:
            # 设置多阶段运行参数
            stage_steps = [stage.get('steps', 0) for stage in stages]
            stage_swap_freqs = [stage.get('swap_frequency', 50) for stage in stages]
            
            self._stages['steps'] = stage_steps
            self._stages['swap_frequencies'] = stage_swap_freqs
            self._stages['current'] = 0
        else:
            # 单一阶段运行
            self._stages['steps'] = [steps]
            
            # 设置交换频率
            if swap_frequencies is None:
                swap_frequencies = 50
                
            if isinstance(swap_frequencies, (list, tuple)):
                self._stages['swap_frequencies'] = swap_frequencies
            else:
                self._stages['swap_frequencies'] = [swap_frequencies]
                
        # 确保交换频率与阶段数量匹配
        if len(self._stages['swap_frequencies']) < len(self._stages['steps']):
            # 用最后一个交换频率填充
            last_freq = self._stages['swap_frequencies'][-1] if self._stages['swap_frequencies'] else 50
            self._stages['swap_frequencies'].extend([last_freq] * (len(self._stages['steps']) - len(self._stages['swap_frequencies'])))
            
        # 更新移动生成器的温度
        self._update_generators_temperature(self._temperature['current'])
        
        # 重置周期计数
        self._current_cycle = 1
        
        # 记录配置信息
        LOGGER.info(f"HRMC参数已设置：")
        LOGGER.info(f"  初始温度 = {self._temperature['initial']:.2f} K")
        LOGGER.info(f"  最终温度 = {self._temperature['final']:.2f} K")
        LOGGER.info(f"  温度更新频率 = {self._stages['temperature_update_freq']}")
        LOGGER.info(f"  总阶段数 = {len(self._stages['steps'])}")
        
        for i, (steps, swap_freq) in enumerate(zip(self._stages['steps'], self._stages['swap_frequencies'])):
            LOGGER.info(f"    阶段{i+1}: 步数={steps}, 交换频率={swap_freq}")
   
    def _update_generators_temperature(self, temperature):
        """更新所有移动生成器的温度。"""
        for group in self.groups:
            if hasattr(group, 'moveGenerator') and hasattr(group.moveGenerator, 'set_temperature'):
                group.moveGenerator.set_temperature(temperature)
        
    def verify_constraints_data_updated(self):
        """验证约束数据是否已正确更新。"""
        for c in self.constraints:
            if hasattr(c, 'data') and c.data is None:
                LOGGER.warning(f"约束 {c.__class__.__name__} 的数据为None，尝试重新计算...")
                if hasattr(c, 'compute_data'):
                    c.compute_data(update=True)
            
    def perform_atom_swap(self, element=None, target_element=None):
        """
        执行原子位置交换
    
        Parameters:
            element (str): 源交换元素。如果为None，随机选择一种元素
            target_element (str): 目标交换元素。如果为None，与源元素相同(同种元素交换)
                                  如果指定，则执行不同元素间交换
    
        Returns:
            bool: 交换是否成功
        """
        # 检查是否已设置交换生成器
        if not self._available_swap_elements:
            LOGGER.warning("未设置原子交换生成器或没有可用的元素进行交换")
            return False
    
        # 如果没有指定源元素，随机选择一个可用元素
        if element is None:
            element = np.random.choice(self._available_swap_elements)
    
        # 判断交换类型: 同元素交换或异元素交换
        if target_element is None or target_element == element:
            # 执行同元素内部交换
            return self._perform_same_element_swap(element)
        else:
            # 执行不同元素之间的交换
            return self._perform_different_element_swap(element, target_element)

    def _perform_same_element_swap(self, element):
        """执行同一元素内的原子交换"""
        # 获取该元素的交换生成器
        if element not in self._swap_generators:
            LOGGER.warning(f"找不到元素'{element}'的交换生成器")
            return False
    
        # 获取该元素的所有原子组
        element_groups = [g for g in self.groups if len(g.indexes) == 1 and 
                         self.allElements[g.indexes[0]] == element]
    
        if len(element_groups) < 2:
            LOGGER.warning(f"元素'{element}'的原子数量不足，无法进行交换")
            return False
    
        # 随机选择两个不同的原子组
        group1_idx = np.random.randint(0, len(element_groups))
        group2_idx = group1_idx
        while group2_idx == group1_idx:
            group2_idx = np.random.randint(0, len(element_groups))
    
        group1 = element_groups[group1_idx]
        group2 = element_groups[group2_idx]
    
        # 获取原子索引
        atom1_idx = group1.indexes[0]
        atom2_idx = group2.indexes[0]
        rel1 = self._atomsCollector.get_relative_index(atom1_idx)
        rel2 = self._atomsCollector.get_relative_index(atom2_idx)
    
        # 备份坐标
        old_pos1 = self.realCoordinates[rel1].copy()
        old_pos2 = self.realCoordinates[rel2].copy()
    
        real_indexes = np.array([atom1_idx, atom2_idx], dtype=np.int32)
        relative_indexes = np.array([rel1, rel2], dtype=np.int32)
    
        # 计算移动前约束
        for c in self.constraints:
            c.compute_before_move(realIndexes=real_indexes, relativeIndexes=relative_indexes)
    
        # 交换坐标
        self.realCoordinates[rel1] = old_pos2
        self.realCoordinates[rel2] = old_pos1
    
        # PBC 处理
        if self.isPBC:
            moved_positions = self.realCoordinates[relative_indexes]
            box_coords = transform_coordinates(self.reciprocalBasisVectors, moved_positions) % 1.0
            moved_positions = transform_coordinates(self.basisVectors, box_coords)
            self.realCoordinates[relative_indexes] = moved_positions
    
        # 计算移动后约束
        box_coords = transform_coordinates(self.reciprocalBasisVectors, self.realCoordinates[relative_indexes])
        for c in self.constraints:
            c.compute_after_move(realIndexes=real_indexes, relativeIndexes=relative_indexes, movedBoxCoordinates=box_coords)
    
        # 判断是否接受
        rejected = self.should_step_get_rejected(real_indexes, relative_indexes)
    
        if rejected:
            self.realCoordinates[rel1] = old_pos1
            self.realCoordinates[rel2] = old_pos2
            for c in self.constraints:
                c.reject_move(realIndexes=real_indexes, relativeIndexes=relative_indexes)
            LOGGER.info(f"交换 {element}({atom1_idx}) <-> {element}({atom2_idx}) 被拒绝")
            success = False
        else:
            for c in self.constraints:
                c.accept_move(realIndexes=real_indexes, relativeIndexes=relative_indexes)
            # 使用父类方法增加计数器
            self.increment_accepted()
            LOGGER.info(f"交换 {element}({atom1_idx}) <-> {element}({atom2_idx}) 被接受")
            success = True
    
        # 使用父类方法增加计数器
        self.increment_tried()
        self._cache['simulation']['swap_count'] += 1
    
    
        # 记录交换信息
        step = getattr(self, '_tried', 0)  # 使用父类的计数器
        self._write_data('HRMC_swaps.txt', f"{step} {element}_{element}_SWAP {atom1_idx} {atom2_idx} {success}")
        return success

    def _perform_different_element_swap(self, element1, element2):
        """执行不同元素之间的原子交换"""
        # 获取两种元素的原子组
        elem1_groups = [g for g in self.groups if len(g.indexes) == 1 and 
                       self.allElements[g.indexes[0]] == element1]
        elem2_groups = [g for g in self.groups if len(g.indexes) == 1 and 
                       self.allElements[g.indexes[0]] == element2]
    
        if not elem1_groups or not elem2_groups:
            LOGGER.warning(f"无法进行交换：{element1}或{element2}原子数量不足")
            return False
    
        # 随机选择一对元素组
        group1 = np.random.choice(elem1_groups)
        group2 = np.random.choice(elem2_groups)
    
        # 获取原子索引
        atom1_idx = group1.indexes[0]
        atom2_idx = group2.indexes[0]
        rel1 = self._atomsCollector.get_relative_index(atom1_idx)
        rel2 = self._atomsCollector.get_relative_index(atom2_idx)
    
        # 备份坐标
        old_pos1 = self.realCoordinates[rel1].copy()
        old_pos2 = self.realCoordinates[rel2].copy()
    
        real_indexes = np.array([atom1_idx, atom2_idx], dtype=np.int32)
        relative_indexes = np.array([rel1, rel2], dtype=np.int32)
    
        # 计算移动前约束
        for c in self.constraints:
            c.compute_before_move(realIndexes=real_indexes, relativeIndexes=relative_indexes)
    
        # 交换坐标
        self.realCoordinates[rel1] = old_pos2
        self.realCoordinates[rel2] = old_pos1
    
        # PBC 处理
        if self.isPBC:
            moved_positions = self.realCoordinates[relative_indexes]
            box_coords = transform_coordinates(self.reciprocalBasisVectors, moved_positions) % 1.0
            moved_positions = transform_coordinates(self.basisVectors, box_coords)
            self.realCoordinates[relative_indexes] = moved_positions
    
        # 计算移动后约束
        box_coords = transform_coordinates(self.reciprocalBasisVectors, self.realCoordinates[relative_indexes])
        for c in self.constraints:
            c.compute_after_move(realIndexes=real_indexes, relativeIndexes=relative_indexes, movedBoxCoordinates=box_coords)
    
        # 判断是否接受
        rejected = self.should_step_get_rejected(real_indexes, relative_indexes)
    
        if rejected:
            self.realCoordinates[rel1] = old_pos1
            self.realCoordinates[rel2] = old_pos2
            for c in self.constraints:
                c.reject_move(realIndexes=real_indexes, relativeIndexes=relative_indexes)
            LOGGER.info(f"交换 {element1}({atom1_idx}) <-> {element2}({atom2_idx}) 被拒绝")
            success = False
        else:
            for c in self.constraints:
                c.accept_move(realIndexes=real_indexes, relativeIndexes=relative_indexes)
            # 使用父类方法增加计数器
            self.increment_accepted()
            LOGGER.info(f"交换 {element1}({atom1_idx}) <-> {element2}({atom2_idx}) 被接受")
            success = True
    
        # 使用父类方法增加计数器
        self.increment_tried()
        self._cache['simulation']['swap_count'] += 1
    
    
        # 记录交换信息
        step = getattr(self, '_tried', 0)  # 使用父类的计数器
        self._write_data('HRMC_swaps.txt', f"{step} {element1}_{element2}_SWAP {atom1_idx} {atom2_idx} {success}")
        return success

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
        
        # 设置当前阶段参数
        current_stage = self._stages['current']
        swap_frequency = self._stages['swap_frequencies'][current_stage] if self._stages['swap_frequencies'] else 0
        temp_update_freq = self._stages['temperature_update_freq']
        
        LOGGER.info(f"设置原子交换频率: 每{swap_frequency}步尝试一次交换")

        # 如果未设置阶段步数，使用输入的步数
        if not self._stages['steps']:
            self._stages['steps'] = [numberOfSteps]
                
        # 确保交换生成器已初始化
        if not self._swap_generators:
            self.setup_atom_swaps()
                        
        # 记录当前步数
        current_step = 0
        
        # 统计局部接受率
        local_tried = 0
        local_accepted = 0
        
        # 开始模拟循环
        for step in range(1, numberOfSteps + 1):
            # 执行原子交换（如果到达交换频率）
            if swap_frequency > 0 and step % swap_frequency == 0:
                # 随机选择交换类型: 50%概率同元素交换，50%概率异元素交换
                if len(self._available_swap_elements) >= 2 and random.random() < 0.5:
                    # 执行异元素交换
                    elements = random.sample(self._available_swap_elements, 2)
                    self.perform_atom_swap(element=elements[0], target_element=elements[1])
                else:
                    # 执行同元素交换
                    element = random.choice(self._available_swap_elements)
                    self.perform_atom_swap(element=element)

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
            
            # 更新温度
            if step % temp_update_freq == 0:
                self._temperature['chi'] = self._temperature['chi_0'] * (self._temperature['cooling_rate'] ** self._current_cycle)
                self._temperature['current'] = self._temperature['initial'] * (self._temperature['cooling_rate'] ** self._current_cycle)
                self._update_generators_temperature(self._temperature['current'])
                self._current_cycle += 1  # 更新冷却轮数
                
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
                
                # 计算改进率
                pdf_improvement = (self.pdf_constraint.initialError - pdf_error) / self.pdf_constraint.initialError * 100
                
                LOGGER.info(f"步骤 {step}/{numberOfSteps} ({step/numberOfSteps*100:.1f}%)")
                LOGGER.info(f"- 接受率: 局部={local_rate:.2%}, 全局={global_rate:.2%}")
                LOGGER.info(f"- PDF误差: {pdf_error:.6f} (改进: {pdf_improvement:.1f}%)")
                LOGGER.info(f"- 当前归一温度 T_chi: {self._temperature['chi']:.6f}")
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
        
        # 保留原子形式因子缓存(计算代价高但大小小)
    
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
        # 生成缓存键
        q_key = hashlib.md5(q_array.tobytes()).hexdigest()
        cache_key = f"{elem}_{q_key}"
        
        # 检查缓存
        if cache_key in self._cache['atom_ff']:
            return self._cache['atom_ff'][cache_key]
            
        # 不在缓存中，计算散射因子
        if elem in self._wk_params:
            a, b, c = self._wk_params[elem]
            
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
        new_gr = self._calculate_gr_incremental()
        
        # 存入缓存
        self._cache['gr'][cache_key] = new_gr
        self._cache['gr']['latest'] = new_gr  # 同时保存最新结果
        
        return new_gr, cache_key
        
    def _calculate_gr_incremental(self):
        """
        使用增量计算方法计算G(r)，比原始方法更高效
        
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

        # 创建Q值数组用于傅里叶变换 (减少Q点数以加速计算)
        q_min, q_max, dq = 0.01, 30.0, 0.01
        q_array = np.arange(q_min, q_max, dq)

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
        if q_max not in self._cache['lorch']:
            self._cache['lorch'][q_max] = np.sinc(q_array / q_max)

        lorch_window = self._cache['lorch'][q_max]

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

        # === 应用 Lorch 窗函数 ===
        lorch_window = np.sinc(q_array / q_max)  # sinc(x) = sin(pi x)/(pi x)
        fq = q_array * (s_q_normalized - 1.0) * lorch_window  # 直接嵌入窗函数
        
        # 通过傅里叶变换从S(Q)计算G(r)
        gr = np.zeros((len(r_array),), dtype=FLOAT_TYPE)
        
        # 使用矢量化计算傅里叶变换
        for i, ri in enumerate(r_array):
            # 计算积分
            integrand = fq  * np.sin(q_array * ri)
            integral = np.trapz(integrand, dx=dq)
            
            gr[i] = 2.0 / np.pi * integral
        
        # 确保结果长度与实验数据匹配
        if len(gr) != len(exp_gr):
            # 使用简单插值
            from scipy.interpolate import interp1d
            interp_func = interp1d(r_array, gr, bounds_error=False, fill_value="extrapolate")
            gr = interp_func(exp_r)
        
        return gr

    def should_step_get_rejected(self, realIndexes, relativeIndexes):
        """
        实现基于G(r)约束的接受准则。
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
            
        :Returns:
            #. rejected (bool): 如果应该拒绝步骤则返回True
        """
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
        
        # 使用G(r)约束的误差变化
        delta_total = delta_gr

        # 计算当前温度下的接受概率
        T_chi = self._temperature['chi']
                                
        # 计算接受概率 P = min[1, exp(-delta_total / T_chi)]
        prob = min(1, np.exp(-delta_total / T_chi))
        
        # 记录数据
        step = getattr(self, 'tried', 0)  # 使用父类计数器
        self._write_data('HRMC_acceptance.txt', 
            f"{step}    {moved_elem}    {prob:.6f}   {delta_gr:.6f}    {T_chi:.6f}")
        self._write_data('HRMC_gr.txt',
            f"{step} {self.pdf_constraint.standardError:.6f}  {delta_gr:.6f}")
        
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
            f.write(f"初始温度: {self._temperature['initial']} K\n")
            f.write(f"最终温度: {self._temperature['current']} K\n")
            f.write(f"原子数: {len(self.allElements)}\n")
            f.write(f"交换次数: {self._cache['simulation']['swap_count']}\n\n")
            
            f.write("## 结果统计\n")
            acceptance_rate = self.accepted / max(1, self.tried) * 100
            f.write(f"总步数: {self.tried}\n")
            f.write(f"接受步数: {self.accepted}\n")
            f.write(f"接受率: {acceptance_rate:.2f}%\n\n")
            
            f.write("## 约束误差\n")
            f.write(f"G(r)误差: {self.pdf_constraint.standardError:.6f}\n")
            
                
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
