import os 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from fullrmc.Engine import Engine
from fullrmc.Constraints.NMRDistanceConstraints import NMRDistanceConstraint
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
    参考了MDHRMCEngine的移动策略
    """
    
    def __init__(self, amplitude=0.2, 
                 elementParameters=None, 
                 temperature=300.0):
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
        
        # 温度因子影响移动的随机性
        temp_factor = np.sqrt(self._temperature / 300.0)
        
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


class ElementAwareGroupSelector(GroupSelector):
    """基于元素类型的原子组选择器"""
    
    def __init__(self, engine=None):
        super(ElementAwareGroupSelector, self).__init__(engine=engine)
        self._elementWeights = {'C': 1.0, 'H': 1.0, 'O': 1.0, 'N': 1.0}
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
    """模拟管理类，处理模拟运行和结果分析 - 仅使用NMR约束"""
    
    def __init__(self, dir_path, fresh_start=True, repeat_runs=20):
        self.dir_path = dir_path
        self.fresh_start = fresh_start
        self.repeat_runs = repeat_runs
        self.engine = None
        
        # 文件路径
        self.pdb_file = os.path.join(dir_path, "smoke_coal.pdb")
        self.engine_file = os.path.join(dir_path, "coal_engine.rmc")
        
        # 模拟参数
        self.initial_temperature = 500.0    # 初始退火温度 (K)
        self.final_temperature = 100.0      # 最终退火温度 (K)
        self.current_temperature = self.initial_temperature  # 当前温度
        self.cooling_rate = 0.95            # 冷却率
        
        # 元素特定参数
        self.element_parameters = {
            'C': {'amplitude': 0.15, 'weight': 1.0},
            'H': {'amplitude': 0.25, 'weight': 0.8},
            'O': {'amplitude': 0.18, 'weight': 1.2},
            'N': {'amplitude': 0.18, 'weight': 1.2},
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

    def _setup_atom_groups(self):
        """设置原子分组"""
        try:
            if self.engine is None:
                logger.error("引擎未初始化，无法设置原子分组")
                return False
            
            logger.info("开始设置原子分组...")
        
            # 获取所有原子索引和元素
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
        
            # 创建碳原子分组
            if 'C' in element_groups:
                carbon_indexes = element_groups['C']
                # 碳原子使用MD感知移动生成器
                carbon_generator = MDAwareMoveGenerator(
                    amplitude=0.15,
                    elementParameters=self.element_parameters,
                    temperature=self.initial_temperature
                )
                carbon_group = Group(indexes=carbon_indexes, name="Carbon_Group")
                carbon_group.set_move_generator(carbon_generator)
                groups.append(carbon_group)
                logger.info(f"创建碳原子组，包含 {len(carbon_indexes)} 个原子")
        
            # 创建氢原子分组
            if 'H' in element_groups:
                hydrogen_indexes = element_groups['H']
                # 氢原子使用MD感知移动生成器，但振幅更小
                hydrogen_generator = MDAwareMoveGenerator(
                    amplitude=0.1,
                    elementParameters=self.element_parameters,
                    temperature=self.initial_temperature
                )
                hydrogen_group = Group(indexes=hydrogen_indexes, name="Hydrogen_Group")
                hydrogen_group.set_move_generator(hydrogen_generator)
                groups.append(hydrogen_group)
                logger.info(f"创建氢原子组，包含 {len(hydrogen_indexes)} 个原子")
        
            # 创建氧原子分组
            if 'O' in element_groups:
                oxygen_indexes = element_groups['O']
                # 氧原子使用MD感知移动生成器
                oxygen_generator = MDAwareMoveGenerator(
                    amplitude=0.12,
                    elementParameters=self.element_parameters,
                    temperature=self.initial_temperature
                )
                oxygen_group = Group(indexes=oxygen_indexes, name="Oxygen_Group")
                oxygen_group.set_move_generator(oxygen_generator)
                groups.append(oxygen_group)
                logger.info(f"创建氧原子组，包含 {len(oxygen_indexes)} 个原子")
        
            # 创建氮原子分组
            if 'N' in element_groups:
                nitrogen_indexes = element_groups['N']
                # 氮原子使用MD感知移动生成器
                nitrogen_generator = MDAwareMoveGenerator(
                    amplitude=0.12,
                    elementParameters=self.element_parameters,
                    temperature=self.initial_temperature
                )
                nitrogen_group = Group(indexes=nitrogen_indexes, name="Nitrogen_Group")
                nitrogen_group.set_move_generator(nitrogen_generator)
                groups.append(nitrogen_group)
                logger.info(f"创建氮原子组，包含 {len(nitrogen_indexes)} 个原子")
        
            # 创建其他原子的分组
            other_indexes = []
            for elem, indexes in element_groups.items():
                if elem not in ['C', 'H', 'O', 'N']:
                    other_indexes.extend(indexes)
        
            if other_indexes:
                # 其他原子使用通用移动生成器
                other_generator = MDAwareMoveGenerator(
                    amplitude=0.1,
                    elementParameters=self.element_parameters,
                    temperature=self.initial_temperature
                )
                other_group = Group(indexes=other_indexes, name="Other_Group")
                other_group.set_move_generator(other_generator)
                groups.append(other_group)
                logger.info(f"创建其他元素原子组，包含 {len(other_indexes)} 个原子")
        
            # 设置原子组到引擎
            self.engine.set_groups(groups=groups)
        
            # 设置元素感知的组选择器
            selector = ElementAwareGroupSelector(engine=self.engine)
            selector.set_element_weights({
                'C': 1.0,
                'H': 0.7,
                'O': 1.2,
                'N': 1.2,
            })
        
            self.engine.set_group_selector(selector=selector)
        
            logger.info(f"成功创建 {len(groups)} 个原子分组")
            return True
            
        except Exception as e:
            logger.error(f"设置原子分组失败: {str(e)}")
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
    
    def _setup_constraints(self):
        """设置和配置所有约束 - 只添加NMR约束"""
        try:
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未定义，无法添加约束")
                return False
                
            # 创建约束实例 - 使用NMRDistanceConstraint
            logger.info("创建NMR约束实例...")
            
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
            
            # 添加约束到引擎
            logger.info("添加约束到引擎...")
            self.engine.add_constraints([nmr_constraint])
            
            # 验证约束添加成功并计算初始数据
            if len(self.engine.constraints) > 0:
                for constraint in self.engine.constraints:
                    if isinstance(constraint, NMRDistanceConstraint):
                        logger.info("计算约束初始数据...")
                        data, error = constraint.compute_data(update=True)
                    
                        # 输出所有键类型的比例
                        logger.info("初始键比例:")
                        for bond_type, ratio in data.items():
                            logger.info(f"  - {bond_type}: {ratio:.6f}")
                    
                        logger.info(f"初始约束标准误差: {error}")
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
        # 指数退火策略
        self.current_temperature = self.initial_temperature * math.pow(
            self.final_temperature / self.initial_temperature,
            progress
        )
        
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
            
    def run_simulation(self, nsteps=500000, saveFrequency=10000, xyzFrequency=10000):
        """改进的模拟运行方法 - 采用分阶段策略以提高有效性"""
        try:
            logger.info("开始分阶段模拟...")
            
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未初始化，无法运行模拟")
                return False
                
            # 获取NMR约束
            nmr_constraint = None
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    nmr_constraint = constraint
                    break
            
            if nmr_constraint is None:
                logger.error("未找到NMR约束，模拟终止")
                return False
                
            # 强制初始化约束
            logger.info("重新计算约束数据...")
            data, initial_error = nmr_constraint.compute_data(update=True)
            
            # 验证计算结果
            logger.info(f"初始标准误差: {initial_error}")
            if initial_error is None or initial_error < 0.0001:
                initial_error = 0.2
                logger.warning(f"约束计算返回无效误差，使用默认值: {initial_error}")
                nmr_constraint.set_standard_error(initial_error)
            
            # 第1阶段: 小幅度热身，目的是实现高接受率
            logger.info("模拟第1阶段: 小幅度热身 (50000步)")
            total_steps = 0
            
            # 设置小幅度
            for group in self.engine.groups:
                if hasattr(group.moveGenerator, 'set_amplitude'):
                    if 'Hydrogen' in group.name:
                        group.moveGenerator.set_amplitude(0.15)
                    else:
                        group.moveGenerator.set_amplitude(0.10)
                    
                # 设置初始温度
                if hasattr(group.moveGenerator, 'set_temperature'):
                    group.moveGenerator.set_temperature(self.initial_temperature)
            
            # 设置初始温度
            self.current_temperature = self.initial_temperature
            
            # 设置初始阶段较大的容差和温和的拒绝概率
            init_tolerances = {k: min(0.1, 2*v) for k, v in nmr_constraint.tolerances.items()}
            nmr_constraint.set_tolerances(init_tolerances)
            nmr_constraint.set_reject_probability(0.2)
            
            # 运行初始阶段
            self.engine.run(numberOfSteps=100000, saveFrequency=10000)
            total_steps += 100000
            
            # 检查第一阶段的接受率
            accept_rate = self.engine.accepted / max(1, self.engine.tried)
            logger.info(f"第1阶段接受率: {accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")
            
            # 重新计算约束数据
            data, phase1_error = nmr_constraint.compute_data(update=True)

            # 打印第一阶段的键统计信息
            self.print_bond_statistics(stage_name="第1阶段")
            
            # 调整拒绝概率
            if accept_rate > 0.7:
                logger.warning("第1阶段接受率过高，减少拒绝概率")
                nmr_constraint.set_reject_probability(0.15)
            elif accept_rate < 0.3:
                logger.warning("第1阶段接受率过低，增加拒绝概率")
                nmr_constraint.set_reject_probability(0.3)
            
            # 第2阶段: 中等幅度
            logger.info("模拟第2阶段: 中等幅度 (200000步)")
            
            # 微调容差 - 基于第一阶段的结果
            phase2_tolerances = {}
            for bond_type, tol in nmr_constraint.tolerances.items():
                # 获取当前比例和目标比例
                current = data.get(bond_type, 0)
                target = nmr_constraint.targetBondRatios.get(bond_type, 0)
                diff = abs(current - target)
                
                # 如果差异大，保持较大容差；如果差异小，适当减小容差
                if diff > 0.05:
                    phase2_tolerances[bond_type] = tol  # 保持容差
                else:
                    phase2_tolerances[bond_type] = max(0.04, tol * 0.8)  # 减小容差
                    
            nmr_constraint.set_tolerances(phase2_tolerances)
            
            # 调整温度
            self.current_temperature *= 0.8  # 降低20%
            
            # 设置中等幅度
            for group in self.engine.groups:
                if hasattr(group.moveGenerator, 'set_amplitude'):
                    if 'Hydrogen' in group.name:
                        group.moveGenerator.set_amplitude(0.25)
                    else:
                        group.moveGenerator.set_amplitude(0.18)
                
                # 更新温度
                if hasattr(group.moveGenerator, 'set_temperature'):
                    group.moveGenerator.set_temperature(self.current_temperature)
                
            # 保存第1阶段状态
            self.engine.save(path=os.path.join(self.dir_path, "stage1_checkpoint.rmc"))
            
            # 运行第2阶段
            self.engine.run(numberOfSteps=200000, saveFrequency=1000)
            total_steps += 200000
            
            # 检查第二阶段的接受率和误差
            phase2_accept_rate = self.engine.accepted / max(1, self.engine.tried)
            data, phase2_error = nmr_constraint.compute_data(update=True)
            logger.info(f"第2阶段接受率: {phase2_accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")
            logger.info(f"第2阶段误差: {phase2_error:.6f} (初始误差的 {phase2_error/initial_error:.2f})")

            # 打印第二阶段的键统计信息
            self.print_bond_statistics(stage_name="第2阶段")
            
            # 动态调整参数
            self._adjust_simulation_parameters(
                nmr_constraint, 
                batch=2, 
                total_batches=10, 
                accept_rate=phase2_accept_rate,
                batch_error=phase2_error,
                initial_error=initial_error
            )
            
            # 保存第2阶段状态
            self.engine.save(path=os.path.join(self.dir_path, "stage2_checkpoint.rmc"))
            
            # 第3阶段: 正常幅度，分批运行以便监控进度
            remaining_steps = nsteps - total_steps
            logger.info(f"模拟第3阶段: 正常幅度 (共{remaining_steps}步)")
            
            # 设置正常幅度
            for group in self.engine.groups:
                if hasattr(group.moveGenerator, 'set_amplitude'):
                    amplitude = 0.25
                    # 有针对性地调整不同类型原子的振幅
                    if hasattr(group, 'name'):
                        if 'Carbon' in group.name:
                            amplitude = 0.20
                        elif 'Hydrogen' in group.name:
                            amplitude = 0.30
                        elif 'Oxygen' in group.name:
                            amplitude = 0.22
                        elif 'Nitrogen' in group.name:
                            amplitude = 0.22
                    group.moveGenerator.set_amplitude(amplitude)
                
            # 批量运行，每批次100000步
            batch_size = 100000
            num_batches = max(1, remaining_steps // batch_size)
            
            for batch in range(num_batches):
                batch_steps = min(batch_size, remaining_steps - batch * batch_size)
                if batch_steps <= 0:
                    break
                    
                logger.info(f"运行第3阶段批次 {batch+1}/{num_batches} ({batch_steps}步)")
                batch_start_time = datetime.now()
                
                # 运行一批次
                self.engine.run(numberOfSteps=batch_steps, 
                                saveFrequency=saveFrequency,
                                xyzFrequency=xyzFrequency)
                
                batch_end_time = datetime.now()
                batch_duration = batch_end_time - batch_start_time
                logger.info(f"批次 {batch+1} 运行时间: {batch_duration}")
                
                # 更新总步数
                total_steps += batch_steps
                
                # 检查接受率和误差
                batch_accept_rate = self.engine.accepted / max(1, self.engine.tried)
                data, batch_error = nmr_constraint.compute_data(update=True)
                logger.info(f"当前累计接受率: {batch_accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")
                logger.info(f"当前标准误差: {batch_error:.6f} (初始误差的 {batch_error/initial_error:.2f})")

                # 打印当前批次的键统计信息
                self.print_bond_statistics(stage_name=f"第3阶段批次{batch+1}")
                
                # 动态调整参数
                self._adjust_simulation_parameters(
                    nmr_constraint, 
                    batch=batch+3,  # +3因为已经完成了2个阶段
                    total_batches=num_batches+2, 
                    accept_rate=batch_accept_rate,
                    batch_error=batch_error,
                    initial_error=initial_error
                )
                
                # 每个批次后保存一次
                checkpoint_path = os.path.join(self.dir_path, f"stage3_batch{batch+1}_checkpoint.rmc")
                self.engine.save(path=checkpoint_path)
            
            # 保存最终结果
            final_path = os.path.join(self.dir_path, "final_result.rmc")
            self.engine.save(path=final_path)
            
            # 导出最终结构为PDB和XYZ格式
            final_pdb = os.path.join(self.dir_path, "final_structure.pdb")
            final_xyz = os.path.join(self.dir_path, "final_structure.xyz")
            self.engine.export_pdb(path=final_pdb)
            self.engine.export_xyz(path=final_xyz)
            
            # 总结模拟结果
            self.summarize_simulation_results()
            
            logger.info(f"模拟完成，总步数: {total_steps}，最终结果保存在: {final_path}")
            return True
    
        except Exception as e:
            logger.error(f"模拟出错: {e}")
            logger.error(traceback.format_exc())
            # 尝试保存当前状态
            try:
                if self.engine is not None:
                    self.engine.save(path=os.path.join(self.dir_path, "error_checkpoint.rmc"))
                    logger.info("已保存错误检查点")
            except:
                pass
            return False

    def export_final_pdb(self, output_file=None):
        """导出最终构型的PDB文件"""
        if output_file is None:
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
        if output_file is None:
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
                    bar1 = ax.bar(x - width/2, target_values, width, label='目标比例', color='lightcoral')
                    
                    # 根据是否符合容差使用不同颜色
                    colors = ['limegreen' if c else 'skyblue' for c in compliant]
                    bar2 = ax.bar(x + width/2, current_values, width, label='当前比例', color=colors)
                    
                    # 添加误差条表示容差范围
                    tolerances = [constraint.tolerances.get(bond, 0.05) for bond in bond_types]
                    ax.errorbar(x - width/2, target_values, yerr=tolerances, fmt='none', ecolor='darkred', capsize=5)
                    
                    # 添加图表元素
                    ax.set_xlabel('键类型', fontsize=12)
                    ax.set_ylabel('比例', fontsize=12)
                    
                    # 计算符合率
                    compliant_count = sum(compliant)
                    compliance_rate = compliant_count / len(bond_types) * 100
                    
                    # 设置标题，包含标准误差和符合率
                    title = f"NMR键比例(误差: {constraint.standardError:.4f}, 符合率: {compliance_rate:.1f}%)"
                    ax.set_title(title, fontsize=14)
                    
                    ax.set_xticks(x)
                    ax.set_xticklabels(bond_types, rotation=45, ha='right', fontsize=10)
                    
                    # 自定义图例
                    from matplotlib.patches import Patch
                    legend_elements = [
                        Patch(facecolor='lightcoral', label='目标比例'),
                        Patch(facecolor='limegreen', label='当前比例(在容差内)'),
                        Patch(facecolor='skyblue', label='当前比例(超出容差)')
                    ]
                    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
                    
                    # 添加接受率信息
                    accept_info = f"接受率: {constraint.accepted/max(1, constraint.tried):.2%} ({constraint.accepted}/{constraint.tried})"
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
    
    def summarize_simulation_results(self):
        """总结模拟结果并生成详细报告"""
        try:
            # 创建报告文件
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
            
            # 获取NMR约束信息
            nmr_constraint = None
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    nmr_constraint = constraint
                    break
            
            if nmr_constraint:
                # 添加约束信息
                report_lines.append("二、NMR约束信息")
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
                
                report_lines.append("三、键比例分析")
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
                report_lines.append("四、键类型排序（按差异大小）")
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
            if nmr_constraint:
                report_lines.append(f"最终拒绝概率: {nmr_constraint.rejectProbability:.4f}")
            report_lines.append("")
        
            # 添加文件信息
            report_lines.append("六、相关文件")
            report_lines.append("-"*30)
            report_lines.append(f"输入PDB文件: {os.path.basename(self.pdb_file)}")
            report_lines.append(f"引擎保存文件: {os.path.basename(self.engine_file)}")
            report_lines.append(f"最终结构PDB: final_structure.pdb")
            report_lines.append(f"最终结构XYZ: final_structure.xyz")
            report_lines.append("")
        
            # 写入报告文件
            with open(report_path, 'w') as f:
                f.write('\n'.join(report_lines))
        
            logger.info(f"模拟结果总结报告已生成: {report_path}")
        
            # 同时生成最终的比较图
            self.plot_nmr_comparison(output_file_path=os.path.join(self.dir_path, "final_bond_comparison.png"))
        
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
    fresh_start = True  # 设置为True创建新引擎, False加载现有引擎
    
    try:
        # 创建模拟管理器
        simulation_manager = SimulationManager(DIR_PATH, fresh_start=fresh_start)
        
        # 初始化引擎
        if not simulation_manager.initialize_engine():
            logger.error("引擎初始化失败，程序终止")
            exit(1)
            
        logger.info("引擎初始化成功")
        
        # 执行模拟
        nsteps = 500000  # 总步数
        save_frequency = 10000  # 保存频率
        xyz_frequency = 10000  # 导出XYZ频率
        
        logger.info(f"开始模拟，总步数: {nsteps}")
        
        # 运行模拟
        success = simulation_manager.run_simulation(
            nsteps=nsteps, 
            saveFrequency=save_frequency, 
            xyzFrequency=xyz_frequency
        )
        
        if success:
            logger.info("模拟成功完成")
            
            # 导出最终结构
            final_pdb = simulation_manager.export_final_pdb()
            logger.info(f"最终结构已导出到: {final_pdb}")
            
            # 生成最终报告
            summary_file = simulation_manager.summarize_simulation_results()
            if summary_file:
                logger.info(f"总结报告已生成：{summary_file}")
        else:
            logger.error("模拟过程失败")
        
        # 记录结束时间和总耗时
        end_time = datetime.now()
        duration = end_time - current_time
        logger.info(f"模拟结束时间: {end_time}")
        logger.info(f"总耗时: {duration}")
        logger.info("=" * 50)
        logger.info("模拟任务完成")
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
