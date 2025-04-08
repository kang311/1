import os 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from fullrmc.Engine import Engine
from fullrmc.Constraints.NMRDistanceConstraints import NMRDistanceConstraint
from fullrmc.Constraints.DistanceConstraints import IntraMolecularDistanceConstraint
from fullrmc.Constraints.AtomicCoordinationConstraints import AtomicCoordinationNumberConstraint
from fullrmc.Core.MoveGenerator import MoveGenerator
from fullrmc.Generators.Rotations import RotationGenerator
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Core.Group import Group
from fullrmc.Core.GroupSelector import GroupSelector
from fullrmc.Globals import FLOAT_TYPE, LOGGER
import logging
from datetime import datetime
import random
import traceback
import math


# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("simulation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def reset_engine_constraints_statistics(engine):
    """重置引擎中所有约束的统计数据"""
    try:
        for constraint in engine.constraints:
            # 使用正确的方法重置统计数据
            if hasattr(constraint, '_reset_statistics'):
                constraint._reset_statistics()
            elif hasattr(constraint, 'increment_tried'):
                # 使用替代方法重置
                tried_attr = '_Constraint__tried'
                accepted_attr = '_Constraint__accepted'
                
                if hasattr(constraint, tried_attr) and hasattr(constraint, accepted_attr):
                    object.__setattr__(constraint, tried_attr, 0)
                    object.__setattr__(constraint, accepted_attr, 0)
            else:
                logger.warning(f"无法重置约束类型 {type(constraint).__name__} 的统计数据")
                
        logger.info(f"已重置{len(engine.constraints)}个约束的统计数据")
        return True
    except Exception as e:
        logger.error(f"重置约束统计数据失败: {str(e)}")
        logger.error(traceback.format_exc())
        return False


class MDAwareMoveGenerator(MoveGenerator):
    """
    模仿分子动力学行为的移动生成器
    """
    
    def __init__(self, amplitude=0.2, 
                 elementParameters=None, 
                 temperature=800.0):
        # 初始化基类
        super(MDAwareMoveGenerator, self).__init__()
        
        # 默认振幅
        self._amplitude = FLOAT_TYPE(amplitude)
        
        # 每个元素的特定参数
        self._elementParameters = elementParameters or {
            'C': {'amplitude': 0.15, 'weight': 1.0},
            'H': {'amplitude': 0.25, 'weight': 0.8},
            'O': {'amplitude': 0.18, 'weight': 1.2},
            'N': {'amplitude': 0.18, 'weight': 1.2},
            'S': {'amplitude': 0.15, 'weight': 1.1},
        }
        
        # 模拟温度(K)影响移动的随机性
        self._temperature = temperature
        
        # 存储上次移动方向，实现运动惯性
        self._last_moves = {}
        
        # 惯性因子(0-1)，越大则保持方向的趋势越强
        self._momentum = 0.3

    def check_group(self, group):
        """验证组是否适合该移动生成器"""
        if isinstance(group, Group):
            return True, ""
        return False, "Group must be an instance of Group class"
    
    def set_temperature(self, temperature):
        """设置模拟温度"""
        self._temperature = max(10.0, temperature)  # 至少10K
        # 根据温度调整惯性
        self._momentum = max(0.1, min(0.5, 0.5 * (300.0 / self._temperature)))
        
    def set_amplitude(self, amplitude):
        """设置移动振幅"""
        self._amplitude = FLOAT_TYPE(amplitude)
    
    def get_amplitude(self, elements=None):
        """根据元素类型获取振幅"""
        if elements is None or len(elements) == 0:
            return self._amplitude
            
        # 根据元素类型加权平均振幅
        total_weight = 0
        weighted_amplitude = 0
        
        for elem in elements:
            params = self._elementParameters.get(elem, {'amplitude': self._amplitude, 'weight': 1.0})
            weighted_amplitude += params['amplitude'] * params['weight']
            total_weight += params['weight']
            
        if total_weight > 0:
            return weighted_amplitude / total_weight
        else:
            return self._amplitude
    
    def move(self, coordinates):
        """生成移动，包含惯性效应"""
        if coordinates is None:
            return coordinates
            
        moved = np.array(coordinates, dtype=FLOAT_TYPE)
        num_atoms = len(coordinates)
        
        # 温度因子影响移动的随机性 - 增大温度影响系数
        # 从平方根关系改为0.4次方关系，增大高温扰动范围
        temp_factor = (self._temperature / 300.0)**0.4  
        
        # 为每个原子生成移动
        for i in range(num_atoms):
            # 生成随机移动
            random_move = (np.random.random(3) - 0.5) * 2 * self._amplitude * temp_factor
            
            # 获取上次移动方向(如果有)
            if i in self._last_moves:
                last_move = self._last_moves[i]
                # 融合上次方向(惯性)
                current_move = (1 - self._momentum) * random_move + self._momentum * last_move
                # 保持移动大小一致
                move_magnitude = np.linalg.norm(random_move)
                if move_magnitude > 0:
                    current_move = current_move / np.linalg.norm(current_move) * move_magnitude
            else:
                current_move = random_move
            
            # 保存本次移动用于下次的惯性计算
            self._last_moves[i] = current_move
            
            # 应用移动
            moved[i] += current_move
            
        return moved


class CompositeGenerator(MoveGenerator):
    """
    组合多种移动生成器的复合生成器
    可以同时使用平移和旋转扰动
    """
    
    def __init__(self, generators=None, weights=None):
        """
        初始化复合移动生成器
        
        :Parameters:
            #. generators (list): 移动生成器列表
            #. weights (list): 每个生成器的权重列表，用于控制选择概率
        """
        # 初始化基类
        super(CompositeGenerator, self).__init__()
        
        # 保存生成器列表
        self._generators = generators or []
        
        # 设置权重，默认均等权重
        if weights is None:
            self._weights = [1.0] * len(self._generators)
        else:
            assert len(weights) == len(self._generators), "权重数量必须与生成器数量匹配"
            self._weights = [float(w) for w in weights]
        
        # 规范化权重
        total_weight = sum(self._weights)
        if total_weight > 0:
            self._weights = [w/total_weight for w in self._weights]
    
    def check_group(self, group):
        """验证组是否适合该移动生成器"""
        # 只要有一个生成器支持该组，就认为是有效的
        for gen in self._generators:
            is_valid, _ = gen.check_group(group)
            if is_valid:
                return True, ""
        return False, "Group is not compatible with any of the composite generators"
    
    def add_generator(self, generator, weight=1.0):
        """添加新的移动生成器到复合生成器"""
        self._generators.append(generator)
        self._weights.append(float(weight))
        
        # 重新规范化权重
        total_weight = sum(self._weights)
        self._weights = [w/total_weight for w in self._weights]
    
    def set_temperature(self, temperature):
        """设置所有生成器的模拟温度"""
        for gen in self._generators:
            if hasattr(gen, 'set_temperature'):
                gen.set_temperature(temperature)
    
    def set_amplitude(self, amplitude):
        """设置所有生成器的移动振幅"""
        for gen in self._generators:
            if hasattr(gen, 'set_amplitude'):
                gen.set_amplitude(amplitude)
    
    def move(self, coordinates):
        """随机选择一个生成器进行移动"""
        if coordinates is None:
            return coordinates
            
        if not self._generators:
            return coordinates
        
        # 随机选择一个生成器进行扰动
        selected_generator = np.random.choice(self._generators, p=self._weights)
        return selected_generator.move(coordinates)


class ElementAwareGroupSelector(GroupSelector):
    """基于元素类型的原子组选择器"""
    
    def __init__(self, engine=None):
        super(ElementAwareGroupSelector, self).__init__(engine=engine)
        self._elementWeights = {'C': 1.0, 'H': 0.7, 'O': 1.2, 'N': 1.2, 'S': 1.1}
        # 跟踪每个组的选择和接受计数
        self._groupSelectionCount = {}
        self._groupAcceptCount = {}
        self._lastSelectedGroup = None
    
    def set_element_weights(self, weights):
        """设置元素权重"""
        if isinstance(weights, dict):
            for elem, weight in weights.items():
                if not isinstance(weight, (int, float)):
                    logger.warning(f"元素 {elem} 的权重必须是数字，使用默认权重 1.0")
                    weights[elem] = 1.0
            self._elementWeights = weights
    
    def select_index(self):
        """选择一个原子组索引，考虑元素类型"""
        if self.engine is None or not hasattr(self.engine, "groups"):
            return 0  # 默认返回第一个组
        
        # 获取组列表
        groups = self.engine.groups
        
        if not groups:
            logger.warning("没有可用的原子组，无法选择")
            return 0  # 返回默认值，后续会检查
        
        # 计算权重
        weights = []
        for idx, group in enumerate(groups):
            # 基础权重
            base_weight = 1.0
            
            # 应用元素权重
            if hasattr(group, 'indexes') and hasattr(self.engine, 'allElements'):
                elems = [self.engine.allElements[i] for i in group.indexes]
                elem_weight = 0
                for e in elems:
                    elem_weight += self._elementWeights.get(e, 1.0)
                elem_weight /= max(1, len(elems))
                base_weight *= elem_weight
                
            # 应用选择频率平衡
            selection_count = self._groupSelectionCount.get(idx, 0) + 1  # 避免除零
            acceptance_rate = self._groupAcceptCount.get(idx, 0) / selection_count
            
            # 接受率低的组增加权重，提高被选中机会
            if acceptance_rate < 0.2:
                base_weight *= 1.2
            elif acceptance_rate > 0.8:
                base_weight *= 0.9
            
            weights.append(base_weight)
            
        # 归一化权重
        if sum(weights) > 0:
            weights = [w / sum(weights) for w in weights]
        else:
            weights = [1.0 / len(groups)] * len(groups)
            
        # 随机选择
        chosen_idx = np.random.choice(range(len(groups)), p=weights)
        
        # 更新统计
        self._groupSelectionCount[chosen_idx] = self._groupSelectionCount.get(chosen_idx, 0) + 1
        self._lastSelectedGroup = chosen_idx
        
        return chosen_idx
    
    def move_accepted(self, index):
        """当移动被接受时调用"""
        self._groupAcceptCount[index] = self._groupAcceptCount.get(index, 0) + 1
        
    def move_rejected(self, index):
        """当移动被拒绝时调用"""
        # 当前不需要特殊处理
        pass


class SimulationManager:
    """模拟管理类，处理模拟运行和结果分析"""
    
    def __init__(self, dir_path, fresh_start=True, repeat_runs=20):
        self.dir_path = dir_path
        self.fresh_start = fresh_start
        self.repeat_runs = repeat_runs
        self.engine = None
        
        # 文件路径
        self.pdb_file = os.path.join(dir_path, "smoke_coal.pdb")
        self.engine_file = os.path.join(dir_path, "coal_engine.rmc")
        
        # 模拟参数 - 修改为建议值
        self.initial_temperature = 800.0    # 提高初始退火温度到 800K
        self.final_temperature = 100.0      # 最终退火温度 (K)
        self.current_temperature = self.initial_temperature  # 当前温度
        self.cooling_rate = 0.98            # 降低冷却率至0.98
        
        # 元素特定参数
        self.element_parameters = {
            'C': {'amplitude': 0.55, 'weight': 1.0},
            'H': {'amplitude': 0.55, 'weight': 0.7},
            'O': {'amplitude': 0.58, 'weight': 1.2},
            'N': {'amplitude': 0.58, 'weight': 1.2},
            'S': {'amplitude': 0.55, 'weight': 1.1},
        }
        
        # 元素分组大小 - 根据建议设置不同元素的分组大小
        self.element_group_sizes = {
            "C": 12,
            "H": 8,
            "O": 10, 
            "N": 8,
            "S": 7,
        }
        
        # 元素共价半径，用于构建DistanceConstraint
        self.covalent_radii = {
            "H": 0.31,
            "C": 0.76,
            "O": 0.66,
            "N": 0.71,
            "S": 1.05,
        }
        
        # 验证文件
        self._validate_files()
        
    def _validate_files(self):
        """验证所需文件是否存在"""
        for f in [self.pdb_file]:  # 只检查PDB文件
            if not os.path.isfile(f):
                raise FileNotFoundError(f"文件 {f} 不存在")
            logger.info(f"文件 {f} 已找到")
    
    def initialize_engine(self):
        """初始化和配置模拟引擎 - 增强错误处理和状态验证"""
        try:
            logger.info(f"尝试初始化引擎，fresh_start={self.fresh_start}")
        
            # 检查引擎文件是否存在
            engine_exists = os.path.isfile(self.engine_file)
            logger.info(f"引擎文件存在: {engine_exists}")
        
            if not engine_exists or self.fresh_start:
                logger.info("创建新的引擎...")
                self.engine = Engine(path=self.engine_file, freshStart=True)
            
                # 设置PDB
                logger.info(f"设置PDB文件: {self.pdb_file}")
                self.engine.set_pdb(self.pdb_file)
            
                # 验证引擎状态
                if not self._verify_engine_ready():
                    logger.error("引擎初始化失败，请检查PDB文件格式")
                    return False
                
                # 配置定制分组
                logger.info("设置原子分组...")
                self._setup_atom_groups()
                
                # 配置约束
                logger.info("开始设置约束...")
                self._setup_constraints()
            
                # 验证约束是否正确添加
                if len(self.engine.constraints) == 0:
                    logger.error("约束添加失败")
                    return False
                
                logger.info(f"成功添加了 {len(self.engine.constraints)} 个约束")
            
                # 保存引擎状态
                logger.info("保存引擎状态...")
                self.engine.save()
                logger.info("引擎初始化完成并保存")
                return True
                
            else:
                logger.info(f"加载现有引擎文件: {self.engine_file}")
                self.engine = Engine(path=self.engine_file)
                self.engine = self.engine.load(self.engine_file)
                
                # 验证加载后的引擎状态
                if not self._verify_engine_ready():
                    logger.error("引擎加载失败，请检查引擎文件")
                    return False
                    
                # 如果需要，重新初始化约束
                if len(self.engine.constraints) == 0:
                    logger.warning("已加载的引擎没有约束，添加约束...")
                    self._setup_constraints()
                else:
                    # 重新计算约束数据确保有效
                    logger.info("重新计算现有约束数据...")
                    for constraint in self.engine.constraints:
                        constraint.compute_data(update=True)
                
                logger.info(f"引擎已加载，包含 {len(self.engine.constraints)} 个约束")
                return True
                
        except Exception as e:
            logger.error(f"引擎初始化失败: {e}")
            logger.error(traceback.format_exc())
            return False

    def _verify_engine_ready(self):
        """验证引擎是否准备就绪"""
        try:
            # 检查关键属性
            if self.engine is None:
                logger.error("引擎为空")
                return False
                
            # 检查坐标是否可用
            if self.engine.realCoordinates is None or self.engine.realCoordinates.size == 0:
                logger.error("引擎实际坐标不可用")
                return False
                
            # 检查原子信息
            if len(self.engine.allElements) == 0:
                logger.error("引擎元素信息不可用")
                return False
                
            # 打印基本信息
            elements = self.engine.allElements
            element_counts = {}
            for elem in elements:
                element_counts[elem] = element_counts.get(elem, 0) + 1
            
            logger.info(f"引擎信息: 总原子数: {len(elements)}")
            for elem, count in element_counts.items():
                logger.info(f"  - {elem}: {count}原子 ({count/len(elements)*100:.1f}%)")
            
            # 检查周期性边界条件
            if self.engine.isPBC:
                logger.info("周期性边界条件: 已启用")
                box_vectors = self.engine.boundaryConditions.get_vectors()
                box_lengths = [np.linalg.norm(v) for v in box_vectors]
                logger.info(f"模拟盒子尺寸: [{box_lengths[0]:.3f}, {box_lengths[1]:.3f}, {box_lengths[2]:.3f}]")
            else:
                logger.info("周期性边界条件: 未启用")
                
            return True
            
        except Exception as e:
            logger.error(f"引擎验证错误: {e}")
            return False            
    
    def bond_range(self, elem1, elem2, epsilon=0.15):
        """
        计算两元素之间的共价键可能范围
        
        :Parameters:
            #. elem1 (str): 第一个元素符号
            #. elem2 (str): 第二个元素符号
            #. epsilon (float): 距离容差(Å)
            
        :Returns:
            #. lower_bound (float): 键距下限
            #. upper_bound (float): 键距上限
        """
        if elem1 not in self.covalent_radii or elem2 not in self.covalent_radii:
            # 对于未定义的元素，使用默认值
            logger.warning(f"元素{elem1}或{elem2}的共价半径未定义，使用默认值")
            r = 1.5  # 默认键距
            return (round(r - epsilon, 2), round(r + epsilon, 2))
            
        # 计算共价键长度
        r = self.covalent_radii[elem1] + self.covalent_radii[elem2]
        
        # 特殊处理某些键类型
        if (elem1 == 'C' and elem2 == 'C') or (elem1 == 'C' and elem2 == 'N'):
            # C-C和C-N键可能有单键、双键或三键，扩大范围
            return (round(r - epsilon - 0.1, 2), round(r + epsilon + 0.1, 2))
        
        # 返回计算的键距范围
        return (round(r - epsilon, 2), round(r + epsilon, 2))
    
    def _setup_constraints(self):
        """
        设置和配置约束
        
        优化后的三阶段约束策略：
        1. Distance约束：提供基础几何诱导，使原子对靠近，为形成键创造物理空间
        2. Coordination约束：提升靠近原子的配位需求，确保结构紧密组装
        3. NMR约束：负责对成键后的结构进行细节拟合
        """
        try:
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未定义，无法添加约束")
                return False
                
            logger.info("设置约束系统...")
            constraints = []
            
            # 1. 建立Distance约束
            logger.info("创建Distance约束...")
            
            # 定义需要约束的元素对
            pairs = [
                ('C','C'), ('C','H'), ('C','O'), ('C','N'), 
                ('H','H'), ('H','O'), ('H','N'), ('O','O'), ('N','N')
            ]
            
            # 对于含有S的系统，添加S相关的约束
            if 'S' in self.engine.allElements:
                pairs.extend([('C','S'), ('H','S'), ('O','S'), ('N','S')])
            
            # 创建键距定义
            pair_defs = [(a, b, self.bond_range(a, b)[1]) for a, b in pairs]
            
            # 创建Distance约束
            dist_constraint = IntraMolecularDistanceConstraint()
            dist_constraint.set_type_definition("element", pair_defs)
            constraints.append(dist_constraint)
            logger.info(f"Distance约束已创建，包含{len(pair_defs)}个元素对定义")
            
            # 2. 建立Coordination约束
            logger.info("创建Coordination约束...")
            
            # 创建配位数定义
            # (元素1, 元素2, 最小距离, 最大距离, 最小配位数, 最大配位数)
            coord_defs = [
                # C原子的配位约束
                ("C", "C", 1.2, 1.7, 0, 4),   # C周围最多4个C
                ("C", "H", 0.9, 1.2, 0, 3),   # C周围最多3个H
                ("C", "O", 1.2, 1.6, 0, 2),   # C周围最多2个O
                ("C", "N", 1.2, 1.6, 0, 2),   # C周围最多2个N
                
                # H原子的配位约束
                ("H", "C", 0.9, 1.2, 0, 1),   # H周围最多1个C (通常H只连接到一个其他原子)
                ("H", "H", 1.2, 1.8, 0, 0),   # H周围不连接H (避免H-H键)
                
                # O原子的配位约束
                ("O", "C", 1.2, 1.6, 0, 2),   # O周围最多2个C
                ("O", "H", 0.9, 1.1, 0, 2),   # O周围最多2个H (如水)
                
                # N原子的配位约束
                ("N", "C", 1.2, 1.6, 0, 3),   # N周围最多3个C
                ("N", "H", 0.9, 1.1, 0, 2)    # N周围最多2个H
            ]
            
            # 如果存在S原子，添加S的配位约束
            if 'S' in self.engine.allElements:
                coord_defs.extend([
                    ("S", "C", 1.6, 1.9, 0, 2),   # S周围最多2个C
                    ("S", "H", 1.2, 1.4, 0, 1)    # S周围最多1个H
                ])
            
            # 创建Coordination约束
            coord_constraint = AtomicCoordinationNumberConstraint()
            coord_constraint.set_coordination_number_definition(coord_defs)
            constraints.append(coord_constraint)
            logger.info(f"Coordination约束已创建，包含{len(coord_defs)}个配位定义")
            
            # 3. 建立NMR约束
            logger.info("创建NMR约束...")
            
            # 设置初始键比例目标 - 针对煤结构的通用设置
            target_ratios = {
                'C-C_SINGLE': 0.40, 
                'C-C_DOUBLE': 0.20,
                'C-C_TRIPLE': 0.05,
                'C-H_SINGLE': 0.20,
                'C-O_SINGLE': 0.05,
                'C-N_SINGLE': 0.05,
                'C-O_DOUBLE': 0.05
            }
            
            # 设置权重
            weights = {
                'C-C_SINGLE': 1.0,
                'C-C_DOUBLE': 1.2,
                'C-C_TRIPLE': 1.0,
                'C-H_SINGLE': 1.1,
                'C-O_SINGLE': 1.1,
                'C-N_SINGLE': 1.1,
                'C-O_DOUBLE': 1.1
            }
            
            # 设置相对宽松的容差
            tolerances = {
                'C-C_SINGLE': 0.05,
                'C-C_DOUBLE': 0.07,
                'C-C_TRIPLE': 0.05,
                'C-H_SINGLE': 0.08,
                'C-O_SINGLE': 0.07,
                'C-N_SINGLE': 0.07,
                'C-O_DOUBLE': 0.07
            }
            
            # 创建NMR约束实例
            nmr_constraint = NMRDistanceConstraint(
                targetBondRatios=target_ratios,
                weights=weights,
                tolerances=tolerances,
                varianceSquared=0.001,
                rejectProbability=0.25  # 初始拒绝概率设置为25%
            )
            
            constraints.append(nmr_constraint)
            logger.info("NMR约束已创建")
            
            # 添加所有约束到引擎
            logger.info("添加约束到引擎...")
            self.engine.add_constraints(constraints)
            
            # 验证约束添加成功并计算初始数据
            if len(self.engine.constraints) > 0:
                logger.info(f"成功添加了{len(self.engine.constraints)}个约束")
                
                # 计算并验证每个约束的初始数据
                for constraint in self.engine.constraints:
                    logger.info(f"计算{constraint.__class__.__name__}约束初始数据...")
                    data, error = constraint.compute_data(update=True)
                    
                    # 对于NMR约束，额外打印键比例信息
                    if isinstance(constraint, NMRDistanceConstraint):
                        logger.info("初始键比例:")
                        for bond_type, ratio in data.items():
                            logger.info(f"  - {bond_type}: {ratio:.6f}")
                        logger.info(f"初始NMR约束标准误差: {error}")
                        
                return True
            else:
                logger.error("没有约束被成功添加")
                return False
    
        except Exception as e:
            logger.error(f"约束设置失败: {e}")
            logger.error(traceback.format_exc())
            return False

    def _sanitize_tolerances(self, raw_tolerances):
        """确保容差都是 float 类型"""
        cleaned = {}
        for k, v in raw_tolerances.items():
            try:
                cleaned[k] = float(v)
            except Exception:
                logger.error(f"容差 '{k}' 的值无法转为 float: {v}")
                raise
        return cleaned

    def get_bond_counts_from_constraint(self, constraint):
        """获取NMR约束中的键数量统计（绝对数值）"""
        try:
            # 确保是NMR约束
            if not isinstance(constraint, NMRDistanceConstraint):
                logger.warning("无法获取键计数：提供的约束不是NMR约束")
                return {}, 0

            # 从引擎获取坐标和元素信息
            coords = self.engine.realCoordinates
            elements = self.engine.allElements

            # 使用OpenBabel计算键统计
            from openbabel import openbabel as ob
            mol = ob.OBMol()
        
            # 添加原子
            for i, coord in enumerate(coords):
                atom = ob.OBAtom()
                atom.SetAtomicNum(ob.GetAtomicNum(elements[i]))
                atom.SetVector(float(coord[0]), float(coord[1]), float(coord[2]))
                mol.AddAtom(atom)
        
            # 连接和感知键
            mol.ConnectTheDots()
            
            # 增强芳香性检测
            aromatic_typer = ob.OBAromaticTyper()
            aromatic_typer.AssignAromaticFlags(mol)
            
            # 感知键序
            mol.PerceiveBondOrders()
        
            # 初始化计数
            bond_counts = {bond_type: 0 for bond_type in constraint.targetBondRatios.keys()}
            total_bonds = 0
        
            # 统计键类型
            for bond in ob.OBMolBondIter(mol):
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
            
                # 获取元素符号
                begin_idx = begin_atom.GetIdx() - 1  # OpenBabel索引从1开始
                end_idx = end_atom.GetIdx() - 1
            
                if begin_idx >= len(elements) or end_idx >= len(elements):
                    continue
            
                begin_element = elements[begin_idx]
                end_element = elements[end_idx]
            
                # 处理键类型
                bond_type = None
            
                if bond.IsAromatic():
                    bond_type = 'AROMATIC'
                else:
                    order = bond.GetBondOrder()
                    if order in [1, 2, 3]:
                        # 按字母顺序排序元素
                        if begin_element > end_element:
                            begin_element, end_element = end_element, begin_element
                    
                        type_str = {1: '_SINGLE', 2: '_DOUBLE', 3: '_TRIPLE'}[order]
                        bond_type = f"{begin_element}-{end_element}{type_str}"
            
                # 更新计数
                if bond_type in bond_counts:
                    bond_counts[bond_type] += 1
                    total_bonds += 1
        
            return bond_counts, total_bonds
    
        except Exception as e:
            logger.error(f"获取键计数出错: {e}")
            logger.error(traceback.format_exc())
            return {}, 0

    def print_bond_statistics(self, stage_name=""):
        """打印键统计信息"""
        try:
            logger.info(f"\n{'='*30}\n{stage_name}阶段键统计信息\n{'='*30}")
        
            # 获取NMR约束
            nmr_constraint = None
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    nmr_constraint = constraint
                    break
                
            if nmr_constraint is None:
                logger.warning("未找到NMR约束，无法显示键统计信息")
                return
        
            # 获取键比例数据（已存在于约束中）
            ratio_data = nmr_constraint.get_constraint_value()
        
            # 获取键数量（绝对值）
            bond_counts, total_bonds = self.get_bond_counts_from_constraint(nmr_constraint)
        
            # 获取目标比例
            target_ratios = nmr_constraint.targetBondRatios
        
            # 打印标题
            logger.info("键类型          |  目标比例  |  实际计数  |  实际比例  |  差异")
            logger.info("-" * 70)
        
            # 计算并打印每种键类型的统计
            for bond_type in sorted(target_ratios.keys()):
                target = target_ratios[bond_type]
                count = bond_counts.get(bond_type, 0)
                actual_ratio = ratio_data.get(bond_type, 0)
                diff = actual_ratio - target
                diff_percent = (diff / max(target, 1e-10)) * 100
            
                logger.info(f"{bond_type:15} | {target:10.4f} | {count:10d} | {actual_ratio:10.4f} | {diff:+.4f} ({diff_percent:+.1f}%)")
        
            # 打印总键数
            logger.info("-" * 70)
            logger.info(f"总键数: {total_bonds}")
        
            # 计算符合容差的键类型数量
            tolerances = nmr_constraint.tolerances
            compliant_count = 0
            total_keys = len(target_ratios)
        
            for bond_type, target in target_ratios.items():
                actual = ratio_data.get(bond_type, 0)
                tolerance = tolerances.get(bond_type, 0.05)
                diff = abs(actual - target)
                if diff <= tolerance:
                    compliant_count += 1
                
            compliance_rate = compliant_count / total_keys * 100
            logger.info(f"符合目标比例的键类型: {compliant_count}/{total_keys} ({compliance_rate:.1f}%)")
            logger.info(f"当前拒绝概率: {nmr_constraint.rejectProbability:.4f}")
            logger.info(f"当前模拟温度: {self.current_temperature:.1f}K")
            logger.info("=" * 70)
        
        except Exception as e:
            logger.error(f"打印键统计信息出错: {e}")
            logger.error(traceback.format_exc())
    
    def _adjust_simulation_parameters(self, nmr_constraint, batch, total_batches, 
                                      accept_rate, batch_error, initial_error):
        """动态调整模拟参数"""
        
        # 计算进度百分比
        progress = batch / max(1, total_batches)
        
        # 调整退火温度
        # 指数退火策略 - 使用新的冷却率0.98
        self.current_temperature = self.initial_temperature * math.pow(
            self.final_temperature / self.initial_temperature,
            progress
        )
        
        # 动态设置惯性因子（根据温度）
        momentum = min(0.5, max(0.0, 1 - self.current_temperature / self.initial_temperature))

        # 基于接受率调整拒绝概率
        current_reject = nmr_constraint.rejectProbability
        target_accept = 0.5 - 0.3 * progress  # 随着模拟进行，目标接受率降低
        
        if accept_rate > target_accept + 0.1:
            # 接受率过高，减小拒绝概率
            new_reject = max(0.05, current_reject * 0.9)
            logger.info(f"接受率({accept_rate:.2f})过高，减小拒绝概率: {current_reject:.2f} -> {new_reject:.2f}")
            nmr_constraint.set_reject_probability(new_reject)
        elif accept_rate < target_accept - 0.1:
            # 接受率过低，增加拒绝概率
            new_reject = min(0.6, current_reject * 1.1)
            logger.info(f"接受率({accept_rate:.2f})过低，增加拒绝概率: {current_reject:.2f} -> {new_reject:.2f}")
            nmr_constraint.set_reject_probability(new_reject)
        
        # 调整容差 - 随着模拟进行，容差逐渐减小
        # 基于当前的误差与初始误差的比例进行调整
        if initial_error > 0:
            error_ratio = batch_error / initial_error
            
            # 如果当前误差显著小于初始误差，则逐步缩小容差
            if error_ratio < 0.6:
                current_tolerances = nmr_constraint.tolerances
                new_tolerances = {}
                
                reduction_factor = 0.9 + 0.1 * error_ratio  # 0.9-1.0的缩减系数
                
                for bond_type, tol in current_tolerances.items():
                    # 根据进度和误差比例调整容差，确保不会小于0.02
                    min_tol = 0.02 + 0.01 * (1 - progress)  # 随着进度增加，允许更小的容差
                    new_tol = max(min_tol, tol * reduction_factor)
                    new_tolerances[bond_type] = new_tol
                
                # 更新容差
                nmr_constraint.set_tolerances(new_tolerances)
                logger.info(f"误差已减小到初始值的{error_ratio:.2f}，调整容差")
        
        # 更新移动生成器的温度
        for group in self.engine.groups:
            if hasattr(group, 'moveGenerator') and hasattr(group.moveGenerator, 'set_temperature'):
                group.moveGenerator.set_temperature(self.current_temperature)
        
        return {
            'temperature': self.current_temperature,
            'reject_probability': nmr_constraint.rejectProbability,
            'tolerances_updated': True
        }
            
    def _setup_atom_groups(self):
        """设置原子分组 - 优化后的方法"""
        try:
            if self.engine is None:
                logger.error("引擎未初始化，无法设置原子分组")
                return False

            logger.info("开始设置原子分组...")

            all_indexes = list(range(len(self.engine.allElements)))
            elements = self.engine.allElements

            # 按元素类型分组
            element_groups = {}
            for idx, elem in enumerate(elements):
                if elem not in element_groups:
                    element_groups[elem] = []
                element_groups[elem].append(idx)

            logger.info(f"发现元素类型: {', '.join(element_groups.keys())}")

            # 创建具有不同移动生成器的分组
            groups = []
            group_names = []  # 用于存储组名
            group_id = 0  # 用于全局命名唯一性

            def create_sub_groups(element, indexes):
                nonlocal group_id  # 确保 group_id 是全局唯一的
                # 获取该元素的分组大小
                sub_group_size = self.element_group_sizes.get(element, 10)  # 默认分组大小为10
                
                # 对小元素特殊处理
                if element in ['N', 'S'] and len(indexes) < 8:
                    # 如果是N或S元素且原子数量少于8，只生成1~2个分组
                    max_groups = 1 if len(indexes) < 4 else 2
                    avg_size = max(1, len(indexes) // max_groups)
                    sub_groups = [indexes[i:i + avg_size] for i in range(0, len(indexes), avg_size)]
                else:
                    # 正常分组
                    sub_groups = [indexes[i:i + sub_group_size] for i in range(0, len(indexes), sub_group_size)]
                
                # 处理每个子组
                for i, sub_group in enumerate(sub_groups):
                    # 跳过空组
                    if not sub_group:
                        continue
                        
                    # 创建组名
                    group_name = f"{element}_{group_id}"
                    group = Group(indexes=sub_group, name=group_name)
                    
                    # 根据组大小和元素类型选择移动生成器
                    if len(sub_group) == 1:
                        # 单原子组使用平移生成器
                        generator = TranslationGenerator(
                            amplitude=0.30 if element == 'H' else 0.20
                        )
                        group.set_move_generator(generator)
                        logger.info(f"创建{element}单原子组{group_name}，使用平移生成器")
                    elif element == 'H':
                        # 氢原子组只使用平移生成器
                        generator = TranslationGenerator(amplitude=0.30)
                        group.set_move_generator(generator)
                        logger.info(f"创建{element}原子子组{group_name}，包含{len(sub_group)}个原子，使用平移生成器")
                    else:
                        # C、O、N等使用复合生成器(平移+旋转)
                        md_generator = MDAwareMoveGenerator(
                            amplitude=0.20,
                            elementParameters=self.element_parameters,
                            temperature=self.initial_temperature
                        )
                        translation_generator = TranslationGenerator(amplitude=0.15)
                        rotation_generator = RotationGenerator(amplitude=10.0)
                        
                        # 创建复合生成器
                        composite_generator = CompositeGenerator(
                            generators=[md_generator, translation_generator, rotation_generator],
                            weights=[0.6, 0.3, 0.1]  # MD生成器权重最高
                        )
                        group.set_move_generator(composite_generator)
                        logger.info(f"创建{element}原子子组{group_name}，包含{len(sub_group)}个原子，使用复合生成器")
                    
                    groups.append(group)
                    group_names.append(group_name)  # 添加组名到列表
                    group_id += 1  # 增加组ID

            # 创建各元素的子组
            for element in element_groups.keys():
                create_sub_groups(element, element_groups[element])           

            # 设置原子组到引擎
            self.engine.set_groups(groups=groups, name=None)  # 添加 name=None

            # 设置元素感知的组选择器
            selector = ElementAwareGroupSelector(engine=self.engine)
            selector.set_element_weights({
                'C': 1.0,
                'H': 0.7,
                'O': 1.2,
                'N': 1.2,
                'S': 1.1,
            })

            self.engine.set_group_selector(selector=selector)

            logger.info(f"成功创建 {len(groups)} 个原子分组")
            return True

        except Exception as e:
            logger.error(f"设置原子分组失败: {str(e)}")
            logger.error(traceback.format_exc())
            return False

    def _get_indexes_by_element(self):
        """获取按元素分类的原子索引"""
        # 实现获取原子索引的方法
        pass

    def _get_move_generator_for_element(self, element):
        """根据元素类型获取相应的移动生成器"""
        if element == 'H':
            return TranslationGenerator(amplitude=0.4)
        else:
            return RotationGenerator(amplitude=10.0)

    def run_simulation(self, nsteps=500000, saveFrequency=10000):
        """
        改进的模拟运行方法 - 采用优化后的三阶段策略
        
        1. 初始阶段: 促成键、压缩结构 (使用Distance + Coordination + 平移/旋转)
        2. 中间阶段: 优化局部结构 & 键型比例拟合 (使用NMR + 复合扰动（MD主导）)
        3. 后期阶段: 收敛到NMR比例并精细优化 (使用纯MDAwareMoveGenerator)
        """
        try:
            logger.info("开始优化后的三阶段模拟...")

            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未初始化，无法运行模拟")
                return False

            # 获取各种约束
            dist_constraint = None
            coord_constraint = None
            nmr_constraint = None
            
            for constraint in self.engine.constraints:
                if isinstance(constraint, IntraMolecularDistanceConstraint):
                    dist_constraint = constraint
                elif isinstance(constraint, AtomicCoordinationNumberConstraint):
                    coord_constraint = constraint
                elif isinstance(constraint, NMRDistanceConstraint):
                    nmr_constraint = constraint
            
            if nmr_constraint is None:
                logger.error("未找到NMR约束，模拟终止")
                return False
                
            if dist_constraint is None:
                logger.warning("未找到Distance约束，优化效果可能受限")
                
            if coord_constraint is None:
                logger.warning("未找到Coordination约束，优化效果可能受限")

            # 强制初始化约束
            logger.info("重新计算约束数据...")
            for constraint in self.engine.constraints:
                data, error = constraint.compute_data(update=True)
                logger.info(f"{constraint.__class__.__name__} 初始误差: {error}")

            # 保存NMR初始误差用于后续参考
            _, initial_nmr_error = nmr_constraint.compute_data(update=True)
            
            # 验证NMR误差
            if initial_nmr_error is None or initial_nmr_error < 0.0001:
                initial_nmr_error = 0.2
                logger.warning(f"NMR约束计算返回无效误差，使用默认值: {initial_nmr_error}")
                nmr_constraint.set_standard_error(initial_nmr_error)

            # 准备移动生成器
            translation_generator = TranslationGenerator(amplitude=0.5)
            rotation_generator = RotationGenerator(amplitude=15.0)  # 增加旋转角度
            md_aware_generator = MDAwareMoveGenerator(
                amplitude=0.6, 
                elementParameters=self.element_parameters, 
                temperature=self.initial_temperature
            )

            #--------------------------------------------------------------------------
            # 第1阶段: 初期键形成阶段 - 主要使用Distance约束和Coordination约束
            #--------------------------------------------------------------------------
            logger.info("\n第1阶段: 初始键形成与结构压缩")
            total_steps = 0
            phase1_steps = 30000  # 增加初始阶段步数

            # 设置初期参数 - 使用TranslationGenerator和RotationGenerator
            for group in self.engine.groups:
                # 应用移动生成器
                if hasattr(group, 'name'):
                    element = group.name.split('_')[0]  # 获取组名中的元素部分
                    
                    if element == 'H':
                        # 氢原子只用平移，振幅较大
                        trans_gen = TranslationGenerator(amplitude=0.5)
                        group.set_move_generator(trans_gen)
                    elif len(group.indexes) == 1:
                        # 单原子组使用平移
                        trans_gen = TranslationGenerator(amplitude=0.40)
                        group.set_move_generator(trans_gen)
                    else:
                        # 多原子组在初期阶段主要使用旋转+平移
                        # 创建复合生成器，旋转权重较大
                        composite = CompositeGenerator(
                            generators=[rotation_generator, translation_generator],
                            weights=[0.7, 0.3]  # 初期以旋转为主
                        )
                        group.set_move_generator(composite)

            # 设置温度和容差
            self.current_temperature = self.initial_temperature
            
            # 设置NMR约束参数 - 初期阶段NMR约束权重较低
            if nmr_constraint:
                init_tolerances = {k: 0.12 for k, v in nmr_constraint.tolerances.items()}
                nmr_constraint.set_tolerances(init_tolerances)
                nmr_constraint.set_reject_probability(0.3)
                
                # 暂时降低NMR约束权重，让Distance和Coordination约束主导
                nmr_constraint._varianceSquared = 0.01  # 降低NMR约束影响
            
            # 运行第一阶段模拟
            logger.info("运行第1阶段模拟...")
            self.engine.run(numberOfSteps=phase1_steps, saveFrequency=1000)
            total_steps += phase1_steps

            # 打印第一阶段结果
            accept_rate = self.engine.accepted / max(1, self.engine.tried)
            logger.info(f"第1阶段接受率: {accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")
            
            # 重置统计
            reset_engine_constraints_statistics(self.engine)
            
            # 重新计算所有约束数据
            for constraint in self.engine.constraints:
                data, error = constraint.compute_data(update=True)
                logger.info(f"{constraint.__class__.__name__} 第1阶段后误差: {error}")
                
            # 打印键统计
            self.print_bond_statistics(stage_name="第1阶段")
            
            # 保存阶段结果
            phase1_path = os.path.join(self.dir_path, "phase1_result.rmc")
            self.engine.save(path=phase1_path)
            logger.info(f"第1阶段结果已保存: {phase1_path}")

            #--------------------------------------------------------------------------
            # 第2阶段: 中期局部优化 - NMR约束开始主导，使用复合扰动（MD为主）
            #--------------------------------------------------------------------------
            logger.info("\n第2阶段: 中期局部结构与键型比例优化")
            phase2_steps = 40000  # 增加中期阶段步数
            
            # 恢复NMR约束影响
            if nmr_constraint:
                nmr_constraint._varianceSquared = 0.001  # 提高NMR约束影响
                nmr_constraint.set_reject_probability(0.2)
                
                # 适当收紧容差
                phase2_tolerances = {k: 0.07 for k, v in nmr_constraint.tolerances.items()}
                nmr_constraint.set_tolerances(phase2_tolerances)
            
            # 对各原子组设置优化的移动生成器
            for group in self.engine.groups:
                if hasattr(group, 'name'):
                    element = group.name.split('_')[0]
                    
                    if element == 'H':
                        # 氢原子继续使用平移，但振幅减小
                        trans_gen = TranslationGenerator(amplitude=0.45)
                        group.set_move_generator(trans_gen)
                    elif len(group.indexes) == 1:
                        # 单原子组使用平移，振幅适中
                        trans_gen = TranslationGenerator(amplitude=0.4)
                        group.set_move_generator(trans_gen)
                    else:
                        # 多原子组使用MD主导的复合生成器
                        md_gen = MDAwareMoveGenerator(
                            amplitude=0.43,
                            elementParameters=self.element_parameters,
                            temperature=self.current_temperature * 0.9
                        )
                        
                        # 创建复合生成器，MD权重增大，旋转权重减小
                        composite = CompositeGenerator(
                            generators=[md_gen, translation_generator, rotation_generator],
                            weights=[0.7, 0.2, 0.1]  # 中期以MD为主
                        )
                        group.set_move_generator(composite)

            # 设置第二阶段温度参数
            self.current_temperature = self.initial_temperature * 0.7  # 560K左右

            # 运行第二阶段模拟
            logger.info("运行第2阶段模拟...")
            self.engine.run(numberOfSteps=phase2_steps, saveFrequency=1000)
            total_steps += phase2_steps

            # 打印第二阶段结果
            phase2_accept_rate = self.engine.accepted / max(1, self.engine.tried)
            logger.info(f"第2阶段接受率: {phase2_accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")
            
            # 重置统计
            reset_engine_constraints_statistics(self.engine)
            
            # 重新计算所有约束数据
            for constraint in self.engine.constraints:
                data, error = constraint.compute_data(update=True)
                logger.info(f"{constraint.__class__.__name__} 第2阶段后误差: {error}")
            
            # 打印键统计
            self.print_bond_statistics(stage_name="第2阶段")
            
            # 保存阶段结果
            phase2_path = os.path.join(self.dir_path, "phase2_result.rmc")
            self.engine.save(path=phase2_path)
            logger.info(f"第2阶段结果已保存: {phase2_path}")

            #--------------------------------------------------------------------------
            # 第3阶段: 后期精细优化 - 纯MD生成器收尾，专注于NMR约束满足
            #--------------------------------------------------------------------------
            logger.info("\n第3阶段: 后期精细优化与NMR约束收敛")
            phase3_steps = 30000
            
            # 提升NMR约束的影响力
            if nmr_constraint:
                nmr_constraint._varianceSquared = 0.0005  # 进一步提高NMR约束影响
                
                # 进一步收紧容差
                phase3_tolerances = {k: 0.05 for k, v in nmr_constraint.tolerances.items()}
                nmr_constraint.set_tolerances(phase3_tolerances)
                
                # 降低拒绝概率以提高接受率
                nmr_constraint.set_reject_probability(0.1)
            
            # 切换到纯MD生成器
            for group in self.engine.groups:
                if hasattr(group, 'name'):
                    element = group.name.split('_')[0]
                    
                    if element == 'H':
                        # 氢原子使用小振幅平移
                        trans_gen = TranslationGenerator(amplitude=0.45)
                        group.set_move_generator(trans_gen)
                    elif len(group.indexes) == 1:
                        # 单原子组使用平移，振幅较小
                        trans_gen = TranslationGenerator(amplitude=0.42)
                        group.set_move_generator(trans_gen)
                    else:
                        # 所有多原子组使用纯MD生成器
                        pure_md = MDAwareMoveGenerator(
                            amplitude=0.15 if element == 'C' else 0.12,
                            elementParameters=self.element_parameters,
                            temperature=self.current_temperature * 0.8
                        )
                        group.set_move_generator(pure_md)

            # 设置第三阶段温度
            self.current_temperature = self.initial_temperature * 0.35  # 280K左右
            
            # 运行第三阶段模拟
            logger.info("运行第3阶段模拟...")
            self.engine.run(numberOfSteps=phase3_steps, saveFrequency=1000)
            total_steps += phase3_steps

            # 打印第三阶段结果
            phase3_accept_rate = self.engine.accepted / max(1, self.engine.tried)
            logger.info(f"第3阶段接受率: {phase3_accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")
            
            # 重新计算所有约束数据
            for constraint in self.engine.constraints:
                data, error = constraint.compute_data(update=True)
                logger.info(f"{constraint.__class__.__name__} 第3阶段后误差: {error}")
            
            # 打印键统计
            self.print_bond_statistics(stage_name="第3阶段")

            # 保存最终结果
            final_path = os.path.join(self.dir_path, "final_result.rmc")
            self.engine.save(path=final_path)

            final_pdb = os.path.join(self.dir_path, "final_structure.pdb")
            self.engine.export_pdb(path=final_pdb)

            self.summarize_simulation_results()

            logger.info(f"模拟完成，总步数: {total_steps}，最终结果保存在: {final_path}")
            return True

        except Exception as e:
            logger.error(f"模拟出错: {e}")
            logger.error(traceback.format_exc())
            try:
                if self.engine is not None:
                    self.engine.save(path=os.path.join(self.dir_path, "error_checkpoint.rmc"))
                    logger.info("已保存错误检查点")
            except:
                pass
            return False

    def export_final_pdb(self, output_file=None):
        """导出最终构型的PDB文件"""
        if filename is None:
            output_file = os.path.join(self.dir_path, "final_structure.pdb")
        
        try:
            # 导出当前构型为PDB文件
            self.engine.export_pdb(path=output_file)
            logger.info(f"最终构型已导出到: {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"导出PDB文件失败: {str(e)}")
            logger.error(traceback.format_exc())
            return None

    def analyze_nmr_constraints(self, output_file=None):
        """分析NMR约束"""
    # 如果没有提供文件名，则使用默认名称
        if filename is None:
            output_file = os.path.join(self.dir_path, "nmr_analysis.txt")
        
        try:
            # 查找NMR约束
            constraint_found = False
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    constraint_found = True
                    
                    # 强制重新计算数据以确保最新结果
                    logger.info("重新计算NMR约束数据...")
                    data, error = constraint.compute_data(update=True)
                    
                    # 验证计算结果
                    if data is None:
                        logger.error("NMR约束数据计算失败，返回空结果")
                        continue
                    
                    # 获取目标键比例和当前模型的键比例
                    target_ratios = constraint.targetBondRatios
                    current_ratios = data
                    
                    # 准备输出内容
                    output_content = ["# NMR约束分析报告\n", 
                                    f"# 生成时间: {datetime.now()}\n",
                                    f"# 标准误差: {constraint.standardError}\n",
                                    f"# 接受率: {constraint.accepted}/{constraint.tried} " +
                                    f"({constraint.accepted/max(1, constraint.tried):.2%})\n\n",
                                    "键类型\t目标比例\t当前比例\t差异\t容差\t是否符合\n"
                                    "-------\t--------\t--------\t----\t----\t--------\n"]
                    
                    # 统计符合目标的键类型数量
                    compliant_count = 0
                    total_keys = len(target_ratios)
                    
                    # 添加每种键类型的信息
                    for bond_type, target in target_ratios.items():
                        current = current_ratios.get(bond_type, 0.0)
                        tolerance = constraint.tolerances.get(bond_type, 0.05)
                        diff = abs(current - target)
                        diff_percent = (diff / max(target, 1e-10)) * 100
                        
                        # 判断是否在容差范围内
                        is_compliant = diff <= tolerance
                        if is_compliant:
                            compliant_count += 1
                            
                        compliance_text = "✓" if is_compliant else "✗"
                        
                        output_line = f"{bond_type}\t{target:.4f}\t{current:.4f}\t{diff:.4f} ({diff_percent:.1f}%)\t{tolerance:.3f}\t{compliance_text}\n"
                        output_content.append(output_line)
                    
                    # 添加整体符合率
                    compliance_rate = compliant_count / total_keys * 100
                    output_content.append(f"\n整体符合率: {compliant_count}/{total_keys} ({compliance_rate:.1f}%)\n")
                    
                    # 写入文件
                    with open(output_file, 'w') as f:
                        f.writelines(output_content)
                    
                    logger.info(f"NMR约束分析已保存到: {output_file}")
                    logger.info(f"键类型符合率: {compliance_rate:.1f}%")
                    return output_file
            
            if not constraint_found:
                logger.error("未找到NMR约束")
                return None
        
        except Exception as e:
            logger.error(f"分析NMR约束出错: {str(e)}")
            logger.error(traceback.format_exc())
            return None

    def plot_nmr_comparison(self, output_file_path='nmr_comparison.png'):
        """绘制NMR键类型比例比较图"""

        # 如果没有提供文件名，则使用默认名称
        if filename is None:
            output_file_path = os.path.join(self.dir_path, "nmr_comparison.png")

        try:
            # 设置matplotlib后端
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            # 查找NMR约束
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    # 确保数据是最新的
                    data, _ = constraint.compute_data(update=True)
                    
                    # 获取目标键比例和当前模型的键比例
                    target_ratios = constraint.targetBondRatios
                    current_ratios = data
                    
                    # 准备绘图数据
                    bond_types = list(target_ratios.keys())
                    target_values = [target_ratios[bond] for bond in bond_types]
                    current_values = [current_ratios.get(bond, 0.0) for bond in bond_types]
                    
                    # 计算符合容差的键类型
                    compliant = []
                    for i, bond in enumerate(bond_types):
                        tolerance = constraint.tolerances.get(bond, 0.05)
                        is_compliant = abs(current_values[i] - target_values[i]) <= tolerance
                        compliant.append(is_compliant)
                    
                    # 创建图表
                    fig, ax = plt.subplots(figsize=(12, 8))
                    
                    # 设置X轴位置
                    x = np.arange(len(bond_types))
                    width = 0.35
                    
                    # 绘制条形图
                    bar1 = ax.bar(x - width/2, target_values, width, label='Target', color='lightcoral')
                    
                    # 根据是否符合容差使用不同颜色
                    colors = ['limegreen' if c else 'skyblue' for c in compliant]
                    bar2 = ax.bar(x + width/2, current_values, width, label='Current', color=colors)
                    
                    # 添加误差条表示容差范围
                    tolerances = [constraint.tolerances.get(bond, 0.05) for bond in bond_types]
                    ax.errorbar(x - width/2, target_values, yerr=tolerances, fmt='none', ecolor='darkred', capsize=5)
                    
                    # 添加图表元素
                    ax.set_xlabel('Key Type', fontsize=12)
                    ax.set_ylabel('Key Ratio', fontsize=12)
                    
                    # 计算符合率
                    compliant_count = sum(compliant)
                    compliance_rate = compliant_count / len(bond_types) * 100
                    
                    # 设置标题，包含标准误差和符合率
                    title = f"NMR Key Ratio(erro: {constraint.standardError:.4f}, 符合率: {compliance_rate:.1f}%)"
                    ax.set_title(title, fontsize=14)
                    
                    ax.set_xticks(x)
                    ax.set_xticklabels(bond_types, rotation=45, ha='right', fontsize=10)
                    
                    # 自定义图例
                    from matplotlib.patches import Patch
                    legend_elements = [
                        Patch(facecolor='lightcoral', label='Target ratio'),
                        Patch(facecolor='limegreen', label='Current ratio(within tolerance)'),
                        Patch(facecolor='skyblue', label='Current ratio(above tolerance)')
                    ]
                    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
                    
                    # 添加接受率信息
                    accept_info = f"Acceptance rate: {constraint.accepted/max(1, constraint.tried):.2%} ({constraint.accepted}/{constraint.tried})"
                    ax.annotate(accept_info, xy=(0.02, 0.02), xycoords='axes fraction', 
                              bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.3),
                              fontsize=10)
                    
                    # 调整布局
                    plt.tight_layout()
                    
                    # 保存图像
                    plt.savefig(output_file_path, dpi=300)
                    plt.close()
                    
                    logger.info(f"NMR比较图已生成: {output_file_path}")
                    return output_file_path
            
            logger.error("未找到NMR约束")
            return None
                
        except Exception as e:
            logger.error(f"绘制NMR比较图出错: {str(e)}")
            logger.error(traceback.format_exc())
            return None
    
    def summarize_simulation_results(self, filename=None):
        """总结模拟结果并生成详细报告"""
        try:
            # 确定报告文件路径
            if filename is None:
                report_path = os.path.join(self.dir_path, "simulation_summary.txt")

            # 准备报告内容
            report_lines = []
            report_lines.append("="*50)
            report_lines.append("煤体系模拟结果总结报告")
            report_lines.append("="*50)
            report_lines.append(f"生成时间: {datetime.now()}")
            report_lines.append("")
            
            # 添加引擎基本信息
            report_lines.append("一、基本系统信息")
            report_lines.append("-"*30)
            
            # 计算元素统计
            elements = self.engine.allElements
            element_counts = {}
            for elem in elements:
                element_counts[elem] = element_counts.get(elem, 0) + 1
            
            total_atoms = len(elements)
            report_lines.append(f"总原子数: {total_atoms}")
            
            report_lines.append("元素组成:")
            for elem, count in sorted(element_counts.items(), key=lambda x: x[1], reverse=True):
                report_lines.append(f"  - {elem}: {count}原子 ({count/total_atoms*100:.1f}%)")
            
            # 周期性边界信息
            if self.engine.isPBC:
                box_vectors = self.engine.boundaryConditions.get_vectors()
                box_lengths = [np.linalg.norm(v) for v in box_vectors]
                volume = np.abs(np.dot(np.cross(box_vectors[0], box_vectors[1]), box_vectors[2]))
                density = total_atoms / volume
                
                report_lines.append(f"模拟盒子尺寸: [{box_lengths[0]:.3f}, {box_lengths[1]:.3f}, {box_lengths[2]:.3f}]")
                report_lines.append(f"盒子体积: {volume:.3f}")
                report_lines.append(f"原子密度: {density:.5f} 原子/埃^3")
            else:
                report_lines.append("非周期性系统")
            
            report_lines.append("")
            
            # 添加约束信息部分
            report_lines.append("二、约束系统信息")
            report_lines.append("-"*30)
            report_lines.append(f"约束总数: {len(self.engine.constraints)}")
            
            # 分类并统计约束
            constraint_types = {}
            for constraint in self.engine.constraints:
                c_type = constraint.__class__.__name__
                if c_type in constraint_types:
                    constraint_types[c_type] += 1
                else:
                    constraint_types[c_type] = 1
            
            report_lines.append("约束类型分布:")
            for c_type, count in constraint_types.items():
                report_lines.append(f"  - {c_type}: {count}个")
            
            report_lines.append("")
            
            # 获取NMR约束信息
            nmr_constraint = None
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    nmr_constraint = constraint
                    break
            
            if nmr_constraint:
                # 添加NMR约束信息
                report_lines.append("三、NMR约束信息")
                report_lines.append("-"*30)
                report_lines.append(f"标准误差: {nmr_constraint.standardError:.6f}")
                report_lines.append(f"总尝试次数: {nmr_constraint.tried}")
                report_lines.append(f"接受次数: {nmr_constraint.accepted}")
                accept_rate = nmr_constraint.accepted / max(1, nmr_constraint.tried)
                report_lines.append(f"接受率: {accept_rate:.2%}")
                report_lines.append("")
                
                # 获取键信息
                data = nmr_constraint.get_constraint_value()
                target_ratios = nmr_constraint.targetBondRatios
                tolerances = nmr_constraint.tolerances
                
                report_lines.append("四、键比例分析")
                report_lines.append("-"*30)
                
                # 获取键绝对数量
                bond_counts, total_bonds = self.get_bond_counts_from_constraint(nmr_constraint)
                
                report_lines.append(f"总键数: {total_bonds}")
                report_lines.append("")
                report_lines.append("键类型详细分析:")
                report_lines.append(f"{'键类型':<15} | {'目标比例':>10} | {'实际比例':>10} | {'实际计数':>10} | {'差异':>10} | {'容差':>8} | {'状态'}")
                report_lines.append("-"*90)
            
                # 统计符合容差的键类型数量
                compliant_count = 0
                total_keys = len(target_ratios)
            
                # 添加每种键类型的信息
                for bond_type in sorted(target_ratios.keys()):
                    target = target_ratios[bond_type]
                    current = data.get(bond_type, 0.0)
                    count = bond_counts.get(bond_type, 0)
                    tolerance = tolerances.get(bond_type, 0.05)
                    diff = current - target
                    diff_percent = (diff / max(target, 1e-10)) * 100
                
                    # 判断是否在容差范围内
                    is_compliant = abs(diff) <= tolerance
                    if is_compliant:
                        compliant_count += 1
                        status = "符合"
                    else:
                        status = "不符合"
                
                    # 格式化输出
                    line = f"{bond_type:<15} | {target:10.4f} | {current:10.4f} | {count:10d} | {diff:+10.4f} | {tolerance:8.4f} | {status}"
                    report_lines.append(line)
            
                # 添加整体符合率
                compliance_rate = compliant_count / total_keys * 100
                report_lines.append("-"*90)
                report_lines.append(f"符合目标比例的键类型: {compliant_count}/{total_keys} ({compliance_rate:.1f}%)")
            
                report_lines.append("")
            
                # 添加键类型排序（按差异大小）
                report_lines.append("键类型排序（按差异大小）:")
                report_lines.append("-"*30)
            
                # 计算差异并排序
                bond_diffs = []
                for bond_type in target_ratios.keys():
                    target = target_ratios[bond_type]
                    current = data.get(bond_type, 0.0)
                    abs_diff = abs(current - target)
                    diff_percent = (abs_diff / max(target, 1e-10)) * 100
                    bond_diffs.append((bond_type, abs_diff, diff_percent))
            
                # 按差异大小排序
                bond_diffs.sort(key=lambda x: x[1], reverse=True)
            
                report_lines.append(f"{'键类型':<15} | {'目标比例':>10} | {'实际比例':>10} | {'绝对差异':>10} | {'相对差异':>10}")
                report_lines.append("-"*70)
            
                for bond_type, abs_diff, diff_percent in bond_diffs:
                    target = target_ratios[bond_type]
                    current = data.get(bond_type, 0.0)
                    line = f"{bond_type:<15} | {target:10.4f} | {current:10.4f} | {abs_diff:10.4f} | {diff_percent:+10.1f}%"
                    report_lines.append(line)
            
                report_lines.append("")
        
            # 添加模拟参数信息
            report_lines.append("五、模拟参数")
            report_lines.append("-"*30)
            report_lines.append(f"初始退火温度: {self.initial_temperature}K")
            report_lines.append(f"最终退火温度: {self.final_temperature}K")
            report_lines.append(f"当前温度: {self.current_temperature:.1f}K")
            report_lines.append(f"冷却率: {self.cooling_rate}")
            if nmr_constraint:
                report_lines.append(f"最终拒绝概率: {nmr_constraint.rejectProbability:.4f}")
            report_lines.append("")
        
            # 添加文件信息
            report_lines.append("六、相关文件")
            report_lines.append("-"*30)
            report_lines.append(f"输入PDB文件: {os.path.basename(self.pdb_file)}")
            report_lines.append(f"引擎保存文件: {os.path.basename(self.engine_file)}")
            report_lines.append(f"阶段结果文件:")
            report_lines.append(f"  - 第1阶段: phase1_result.rmc")
            report_lines.append(f"  - 第2阶段: phase2_result.rmc")
            report_lines.append(f"最终结构文件:")
            report_lines.append(f"  - 引擎文件: final_result.rmc")
            report_lines.append(f"  - PDB文件: final_structure.pdb")
            report_lines.append("")
        
            # 添加结论
            report_lines.append("七、结论")
            report_lines.append("-"*30)
            
            if nmr_constraint:
                # 根据键比例符合率给出结论
                if compliance_rate >= 90:
                    report_lines.append("模拟结果优异，键类型比例完全符合目标要求。")
                elif compliance_rate >= 70:
                    report_lines.append("模拟结果良好，键类型比例大部分符合目标要求，少量键类型需要进一步优化。")
                elif compliance_rate >= 50:
                    report_lines.append("模拟结果一般，约半数键类型比例符合目标要求，仍需优化。")
                else:
                    report_lines.append("模拟结果不理想，大部分键类型比例未达到目标要求，需要重新审视模拟参数和策略。")
            
            report_lines.append("")
        
            # 写入报告文件
            with open(report_path, 'w') as f:
                f.write('\n'.join(report_lines))
        
            logger.info(f"模拟结果总结报告已生成: {report_path}")
        
            # 同时生成最终的比较图
            graph_filename = os.path.join(self.dir_path,os.path.splitext(os.path.basename(report_path))[0] + "_bond_comparison.png"
            self.plot_nmr_comparison(output_file_path=graph_file_path)
        
            return report_path
        
        except Exception as e:
            logger.error(f"生成模拟结果总结报告失败: {str(e)}")
            logger.error(traceback.format_exc())
            return None

def get_file_paths(dir_path):
    """返回文件路径字典"""
    return {
        'pdb': os.path.join(dir_path, "smoke_coal.pdb"),
        'engine': os.path.join(dir_path, "coal_engine.rmc")
    }

# 主程序
if __name__ == "__main__":
    # 记录当前时间
    current_time = datetime.now()
    logger.info(f"模拟开始时间: {current_time}")
    
    DIR_PATH = os.path.dirname(os.path.realpath(__file__))
    file_paths = get_file_paths(DIR_PATH)
    
    # 根据需要选择是否重新开始
    fresh_start = False  # 设置为True创建新引擎, False加载现有引擎
    
    # 定义重复运行次数
    repeat_runs = 20
    
    try:
        # 创建模拟管理器
        simulation_manager = SimulationManager(DIR_PATH, fresh_start=fresh_start)
        
        # 初始化引擎
        if not simulation_manager.initialize_engine():
            logger.error("引擎初始化失败，程序终止")
            exit(1)
            
        logger.info("引擎初始化成功")
        
        # 执行模拟
        nsteps = 100000  # 总步数
        save_frequency = 10000  # 保存频率
        
        logger.info(f"开始模拟，总步数: {nsteps}，重复运行次数: {repeat_runs}")
        
        # 重复运行指定次数
        for run in range(1, repeat_runs + 1):
            logger.info(f"开始第 {run}/{repeat_runs} 次运行...")
            
            # 运行模拟
            success = simulation_manager.run_simulation(
                nsteps=nsteps, 
                saveFrequency=save_frequency
            )
            
            if success:
                logger.info(f"第 {run}/{repeat_runs} 次运行成功完成")
                
                # 每次运行后导出结构，添加运行编号
                run_pdb = simulation_manager.export_final_pdb(output_file=f"final_structure_run_{run}.pdb")
                logger.info(f"第 {run} 次运行结构已导出到: {run_pdb}")
                
                # 生成中间报告
                if run < repeat_runs:  # 不是最后一次运行
                    summary_file = simulation_manager.summarize_simulation_results(output_file=f"summary_run_{run}.txt")
                    if summary_file:
                        logger.info(f"第 {run} 次运行总结报告已生成：{summary_file}")
            else:
                logger.error(f"第 {run}/{repeat_runs} 次运行失败")
                break  # 如果一次运行失败，退出循环
        
        # 仅在所有运行都完成后生成最终报告
        if success:
            final_pdb = simulation_manager.export_final_pdb()
            logger.info(f"最终结构已导出到: {final_pdb}")
            
            summary_file = simulation_manager.summarize_simulation_results()
            if summary_file:
                logger.info(f"最终总结报告已生成：{summary_file}")
        
        # 记录结束时间和总耗时
        end_time = datetime.now()
        duration = end_time - current_time
        logger.info(f"模拟结束时间: {end_time}")
        logger.info(f"总耗时: {duration}")
        logger.info("=" * 50)
        logger.info(f"完成了 {run}/{repeat_runs} 次模拟运行")
        logger.info("=" * 50)
        
    except Exception as e:
        logger.error(f"程序执行出错: {e}")
        logger.error(traceback.format_exc())
        
        # 尝试保存当前状态
        try:
            if 'simulation_manager' in locals() and simulation_manager.engine is not None:
                simulation_manager.engine.save(path=os.path.join(DIR_PATH, "error_checkpoint.rmc"))
                logger.info("已保存错误检查点")
        except:
            pass
