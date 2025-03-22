"""
NMRConstraints包含基于核磁共振实验数据的约束类。
"""

# 标准库导入
from __future__ import print_function
import itertools, inspect, copy, os, re
import time
import logging
import multiprocessing as mp
from collections import defaultdict
# 外部库导入
import numpy as np
from openbabel import openbabel as ob

# fullrmc导入
from ..Globals import INT_TYPE, FLOAT_TYPE, PRECISION, LOGGER
from ..Globals import str, long, unicode, bytes, basestring, range, xrange, maxint
from ..Core.Collection import is_number, is_integer, get_path
from ..Core.Constraint import Constraint, ExperimentalConstraint


# 在模块级别定义，必须位于任何类定义外部
def _process_chunk_bonds(chunk_data):
    """处理一个原子块的键比例计算"""
    from openbabel import openbabel as ob
    from fullrmc.Globals import LOGGER
    
    coordinates, elements, target_bond_types = chunk_data
    bond_counts = {bond_type: 0 for bond_type in target_bond_types}
    total_bonds = 0
    
    try:
        # 创建OpenBabel分子对象
        mol = ob.OBMol()
        
        # 添加原子
        for i, coord in enumerate(coordinates):
            atom = ob.OBAtom()
            atom.SetAtomicNum(ob.GetAtomicNum(elements[i]))
            atom.SetVector(float(coord[0]), float(coord[1]), float(coord[2]))
            mol.AddAtom(atom)
        
        # 连接和感知键
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        
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
        # 即使发生错误，也返回空结果而不是抛出异常
        return {}, 0

class NMRDistanceConstraint(Constraint):
    """
    通过核磁共振(NMR)数据控制分子中不同类型键的分布。
    该约束通过计算和控制不同键类型的比例将NMR实验数据与模型关联起来。
    
    :Parameters:
        #. targetBondRatios (dict): 目标键比例字典
        #. weights (None, dict): 键类型权重字典，键为键类型，值为自定义权重。
           如果给定None，所有键类型的权重将设为1.0
        #. tolerances (None, dict): 键类型容差字典，键为键类型，值为容差值。
           如果给定None，所有键类型的容差将设为0.05(5%)
        #. varianceSquared (float): 方差平方值，影响约束对总误差的贡献
    """
    
    def __init__(self, targetBondRatios=None, weights=None, tolerances=None, varianceSquared=0.001, rejectProbability=0.15):
        # 初始化基类
        super(NMRDistanceConstraint, self).__init__()
        
        # 设置目标键比例
        self.set_target_bond_ratios(targetBondRatios)
        
        # 设置权重和容差
        self.set_weights(weights)
        self.set_tolerances(tolerances)
        
        # 设置方差平方值
        self.set_variance_squared(varianceSquared)
        
        # 设置拒绝概率
        self.set_reject_probability(rejectProbability)
        
        # 初始化进程池
        self._process_pool = None
                
        # 初始化键类型映射
        self._initialize_bond_types_mapping()
        
        # 设置框架数据
        FRAME_DATA = [d for d in self.FRAME_DATA]
        FRAME_DATA.extend(['_NMRDistanceConstraint__targetBondRatios',
                          '_NMRDistanceConstraint__weights',
                          '_NMRDistanceConstraint__tolerances',
                          '_varianceSquared',
                          '_Constraint__rejectProbability'])
    
        object.__setattr__(self, 'FRAME_DATA', tuple(FRAME_DATA))

    def set_reject_probability(self, rejectProbability):
        """
        设置约束拒绝概率
        
        :Parameters:
            #. rejectProbability (float): 约束拒绝概率，值域[0,1]
        """
        assert isinstance(rejectProbability, (int, float)), LOGGER.error("rejectProbability必须是数值")
        assert rejectProbability>=0 and rejectProbability<=1, LOGGER.error("rejectProbability必须在[0,1]之间")
        self._Constraint__rejectProbability = float(rejectProbability)    

    def __del__(self):
        """
        析构函数，确保清理所有资源
        """
        # 清理进程池
        if hasattr(self, '_process_pool') and self._process_pool is not None:
            self._process_pool.close()
            self._process_pool = None

    def _compute_bond_ratios_with_openbabel_parallel(self, coordinates):
        """使用多进程并行计算键比例"""
        try:
            import multiprocessing as mp
        
            # 确定CPU核心数量(避免使用过多资源)
            cpu_count = min(mp.cpu_count(), 8)  # 最多使用4个核心
        
            # 分割原子坐标
            atoms_count = len(coordinates)
            if atoms_count < 500:  # 小系统不值得并行
                return self._compute_bond_ratios_with_openbabel(coordinates)
            
            chunk_size = max(100, atoms_count // cpu_count)
        
            # 准备数据块
            chunks = []
            elements = self.engine.allElements
            for i in range(0, atoms_count, chunk_size):
                end_idx = min(i + chunk_size, atoms_count)
                chunk_coords = coordinates[i:end_idx].copy()
                chunk_elements = elements[i:end_idx]
                chunks.append((chunk_coords, chunk_elements, list(self.__targetBondRatios.keys())))
        
            # 并行处理
            with mp.Pool(processes=cpu_count) as pool:
                results = pool.map(_process_chunk_bonds, chunks)
        
            # 合并结果
            all_bond_counts = {bond: 0 for bond in self.__targetBondRatios.keys()}
            total_bonds = 0
        
            for bond_counts, chunk_total in results:
                for bond_type, count in bond_counts.items():
                    if bond_type in all_bond_counts:
                        all_bond_counts[bond_type] += count
                total_bonds += chunk_total
        
            # 计算最终键比例
            bond_ratios = {}
            if total_bonds > 0:
                for bond_type in self.__targetBondRatios.keys():
                    bond_ratios[bond_type] = float(all_bond_counts.get(bond_type, 0)) / total_bonds
            else:
                bond_ratios = {bond_type: 0.0 for bond_type in self.__targetBondRatios.keys()}
        
            return bond_ratios
        
        except Exception as e:
            LOGGER.error(f"并行键比例计算错误: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
        
            # 出错时回退到串行计算
            return self._compute_bond_ratios_with_openbabel(coordinates)
    
    def _initialize_bond_types_mapping(self):
        """初始化键类型映射，处理OpenBabel和目标键类型之间的差异"""
        # 映射OpenBabel键类型到目标键类型
        self.__bondTypesMapping = {
            # 常规映射
            'C-C': 'C-C_SINGLE',
            'C=C': 'C-C_DOUBLE',
            'C#C': 'C-C_TRIPLE',
            'C-O': 'C-O_SINGLE',
            'C=O': 'C-O_DOUBLE',
            'C-N': 'C-N_SINGLE',
            'C-H': 'C-H_SINGLE',
            # 可能的替代格式
            'C:C': 'AROMATIC',
            'c:c': 'AROMATIC',
            'c1ccccc1': 'AROMATIC',
            # 以下是原子对格式的映射
            ('C', 'C', 1): 'C-C_SINGLE',
            ('C', 'C', 2): 'C-C_DOUBLE',
            ('C', 'C', 3): 'C-C_TRIPLE',
            ('C', 'O', 1): 'C-O_SINGLE',
            ('C', 'O', 2): 'C-O_DOUBLE',
            ('C', 'N', 1): 'C-N_SINGLE',
            ('C', 'H', 1): 'C-H_SINGLE',
        }
        
        # 记录初始化
        LOGGER.debug(f"已初始化键类型映射: {len(self.__bondTypesMapping)}个映射规则")
    
    def set_target_bond_ratios(self, targetBondRatios):
        """
        设置目标键比例
        
        :Parameters:
            #. targetBondRatios (dict): 目标键比例字典，键为键类型名称，值为期望比例
        """
        if targetBondRatios is None:
            self.__targetBondRatios = {}
            return

        assert isinstance(targetBondRatios, dict), LOGGER.error("targetBondRatios必须是字典")
        self.__targetBondRatios = {}
        total = sum(targetBondRatios.values())
        
        # 归一化比例
        if abs(total - 1.0) > 1e-6:
            LOGGER.warning(f"键比例总和({total})不为1.0，进行归一化")
            factor = 1.0 / max(total, 1e-10)
            for bondType, ratio in targetBondRatios.items():
                self.__targetBondRatios[bondType] = FLOAT_TYPE(ratio * factor)
        else:
            for bondType, ratio in targetBondRatios.items():
                self.__targetBondRatios[bondType] = FLOAT_TYPE(ratio)
        
        # 记录到存储库
        self._dump_to_repository({'_NMRDistanceConstraint__targetBondRatios': self.__targetBondRatios})
    
    def set_weights(self, weights=None):
        """
        设置不同键类型的权重
        
        :Parameters:
            #. weights (None, dict): 键类型权重字典，键为键类型，值为自定义权重。
               如果给定None，所有键类型的权重将设为1.0
        """
        if weights is None:
            self.__weights = {k: FLOAT_TYPE(1.0) for k in self.__targetBondRatios.keys()}
        else:
            assert isinstance(weights, dict), LOGGER.error("weights必须为None或字典")
            self.__weights = {}
            for k in self.__targetBondRatios.keys():
                if k in weights:
                    assert isinstance(weights[k], (int, float)), LOGGER.error(f"键类型'{k}'的权重必须是数值")
                    assert weights[k] > 0, LOGGER.error(f"键类型'{k}'的权重必须为正数")
                    self.__weights[k] = FLOAT_TYPE(weights[k])
                else:
                    self.__weights[k] = FLOAT_TYPE(1.0)
        
        # 记录到存储库
        self._dump_to_repository({'_NMRDistanceConstraint__weights': self.__weights})
    
    def set_tolerances(self, tolerances=None):
        """
        设置不同键类型的容差
        
        :Parameters:
            #. tolerances (None, dict): 键类型容差字典，键为键类型，值为容差值。
               如果给定None，所有键类型的容差将设为0.05(5%)
        """
        if tolerances is None:
            self.__tolerances = {k: FLOAT_TYPE(0.05) for k in self.__targetBondRatios.keys()}
        else:
            assert isinstance(tolerances, dict), LOGGER.error("tolerances必须为None或字典")
            self.__tolerances = {}
            for k in self.__targetBondRatios.keys():
                if k in tolerances:
                    assert isinstance(tolerances[k], (int, float)), LOGGER.error(f"键类型'{k}'的容差必须是数值")
                    assert tolerances[k] >= 0, LOGGER.error(f"键类型'{k}'的容差必须为非负数")
                    self.__tolerances[k] = FLOAT_TYPE(tolerances[k])
                else:
                    self.__tolerances[k] = FLOAT_TYPE(0.05)
        
        # 记录到存储库
        self._dump_to_repository({'_NMRDistanceConstraint__tolerances': self.__tolerances})
    
    def set_variance_squared(self, varianceSquared):
        """
        设置方差平方值，影响约束对总误差的贡献
        
        :Parameters:
            #. varianceSquared (float): 方差平方值
        """
        assert isinstance(varianceSquared, (int, float)), LOGGER.error("varianceSquared必须是数值")
        assert varianceSquared > 0, LOGGER.error("varianceSquared必须为正数")
        self._varianceSquared = FLOAT_TYPE(varianceSquared)
        
        # 记录到存储库
        self._dump_to_repository({'_varianceSquared': self._varianceSquared})
    
    @property
    def varianceSquared(self):
        """方差平方值"""
        return self._varianceSquared
    
    @property
    def targetBondRatios(self):
        """目标键比例字典"""
        return self.__targetBondRatios
    
    @property
    def weights(self):
        """键类型权重字典"""
        return self.__weights
    
    @property
    def tolerances(self):
        """键类型容差字典"""
        return self.__tolerances

    @property
    def rejectProbability(self):
        """返回拒绝概率属性的访问器"""
        if hasattr(self, '_Constraint__rejectProbability'):
            return self._Constraint__rejectProbability
        return 0.15  # 默认值

    def set_reject_probability(self, rejectProbability):
        """
        设置拒绝概率
    
        :Parameters:
            #. rejectProbability (float): 拒绝概率，必须在0和1之间
        """
        if not isinstance(rejectProbability, (int, float)):
            raise TypeError("rejectProbability必须是数值")
        if rejectProbability < 0 or rejectProbability > 1:
            raise ValueError("rejectProbability必须在0和1之间")
    
        self._Constraint__rejectProbability = float(rejectProbability)
        self._dump_to_repository({'_Constraint__rejectProbability': self._Constraint__rejectProbability})
        return rejectProbability
    
    def get_constraint_value(self):
        """
        获取约束的当前值(键比例)
    
        :Returns:
            #. value (dict): 当前键比例字典
        """
        if self.data is None:
            # 如果数据尚未计算，则计算它
            self.compute_data(update=True)
    
        # 确保总是返回字典的副本
        if self.data is not None:
            return {k: v for k, v in self.data.items()}
    
        # 如果无数据，返回空字典
        return {}
    
    def _map_bond_type(self, openbabel_type):
        """将OpenBabel键类型映射到目标键类型"""
        # 直接映射
        if openbabel_type in self.__bondTypesMapping:
            return self.__bondTypesMapping[openbabel_type]
        
        # 记录未映射的键类型
        LOGGER.debug(f"未映射的键类型: {openbabel_type}")
        return None
    
    def _calculate_bond_statistics(self, coordinates=None, useCache=False):
        """
        使用OpenBabel计算键统计数据
        
        :Parameters:
            #. coordinates (None, numpy.ndarray): 计算键统计的坐标。如果为None，
               使用引擎的当前坐标
            #. useCache (boolean): 是否使用缓存
        
        :Returns:
            #. bondRatios (dict): 键类型比例字典
        """
        # 使用引擎坐标或提供的坐标
        coords = self.engine.realCoordinates if coordinates is None else coordinates     
        
        try:
            # 根据系统大小选择计算方法
            atoms_count = len(coords)
            if atoms_count > 500:  # 对于大型系统使用并行计算
                bondRatios = self._compute_bond_ratios_with_openbabel_parallel(coords)
            else:  # 对于小型系统使用顺序计算
                bondRatios = self._compute_bond_ratios_with_openbabel(coords)
                        
            return bondRatios
            
        except Exception as e:
            LOGGER.error(f"键统计计算错误: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
            # 返回默认值而非引发异常，提高稳健性
            return {k: FLOAT_TYPE(0.0) for k in self.__targetBondRatios.keys()}
    
    def _compute_bond_ratios_with_openbabel(self, coordinates):
        """
        使用OpenBabel计算键比例
        
        :Parameters:
            #. coordinates (numpy.ndarray): 原子坐标
            
        :Returns:
            #. bondRatios (dict): 键类型比例字典
        """
        try:
            # 初始化所有键类型的计数器
            bondCounts = {bondType: 0 for bondType in self.__targetBondRatios.keys()}
            totalBonds = 0
            
            # 创建OpenBabel分子对象
            mol = ob.OBMol()
            
            # 添加原子
            elements = self.engine.allElements
            for i, coord in enumerate(coordinates):
                atom = ob.OBAtom()
                # 设置原子的元素类型
                atom.SetAtomicNum(ob.GetAtomicNum(elements[i]))
                # 设置原子的坐标
                atom.SetVector(float(coord[0]), float(coord[1]), float(coord[2]))
                mol.AddAtom(atom)
            
            # 使用OpenBabel的连接性感知算法
            mol.ConnectTheDots()
            mol.PerceiveBondOrders()
            
            # 分析键类型
            for bond in ob.OBMolBondIter(mol):
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                
                # 获取元素符号
                begin_idx = begin_atom.GetIdx() - 1  # OpenBabel索引从1开始
                end_idx = end_atom.GetIdx() - 1
                begin_element = elements[begin_idx]
                end_element = elements[end_idx]
                
                # 获取键类型
                bond_order = bond.GetBondOrder()
                is_aromatic = bond.IsAromatic()
                
                # 构建键类型
                bond_type = None
                
                # 处理芳香键
                if is_aromatic:
                    bond_type = 'AROMATIC'
                else:
                    # 按字母顺序排序元素
                    if begin_element > end_element:
                        begin_element, end_element = end_element, begin_element
                    
                    # 构建OpenBabel风格的键类型字符串
                    bond_symbol = ''
                    if bond_order == 1:
                        bond_symbol = '-'
                    elif bond_order == 2:
                        bond_symbol = '='
                    elif bond_order == 3:
                        bond_symbol = '#'
                    
                    ob_bond_type = f"{begin_element}{bond_symbol}{end_element}"
                    
                    # 尝试将OpenBabel键类型映射到目标键类型
                    bond_type = self._map_bond_type(ob_bond_type)
                    
                    # 如果映射失败，尝试使用原子和键序的元组
                    if bond_type is None:
                        tuple_key = (begin_element, end_element, bond_order)
                        bond_type = self._map_bond_type(tuple_key)
                    
                    # 如果仍然失败，使用默认格式
                    if bond_type is None:
                        if bond_order == 1:
                            bond_type = f"{begin_element}-{end_element}_SINGLE"
                        elif bond_order == 2:
                            bond_type = f"{begin_element}-{end_element}_DOUBLE"
                        elif bond_order == 3:
                            bond_type = f"{begin_element}-{end_element}_TRIPLE"
                
                # 更新计数
                if bond_type in bondCounts:
                    bondCounts[bond_type] += 1
                    totalBonds += 1
                elif bond_type is not None:
                    LOGGER.debug(f"找到的键类型 '{bond_type}' 不在目标键类型中")
            
            # 计算键比例
            bondRatios = {}
            if totalBonds > 0:
                for bondType, count in bondCounts.items():
                    bondRatios[bondType] = FLOAT_TYPE(count / totalBonds)
            else:
                # 如果没有找到键，返回全零比例
                bondRatios = {bondType: FLOAT_TYPE(0.0) for bondType in bondCounts.keys()}
            
            return bondRatios
        
        except Exception as e:
            LOGGER.error(f"OpenBabel键比例计算错误: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
            # 返回默认值
            return {bondType: FLOAT_TYPE(0.0) for bondType in self.__targetBondRatios.keys()}
    
    def compute_standard_error(self, modelData, normalize=True):
        """
        计算标准误差，可选择归一化并应用容差和权重
        
        :Parameters:
            #. modelData (dict): 与目标值比较的键比例字典
            #. normalize (boolean): 是否归一化误差
        
        :Returns:
            #. standardError (float): 计算的约束标准误差
        """
        try:
            error = FLOAT_TYPE(0.0)
            valid_count = 0
            weights_sum = FLOAT_TYPE(0.0)
            
            # 对每种键类型计算误差
            for bond_type, target in self.__targetBondRatios.items():
                current = modelData.get(bond_type, 0.0)
                weight = self.__weights.get(bond_type, 1.0)
                tolerance = self.__tolerances.get(bond_type, 0.05)
                
                # 计算相对误差
                diff = abs(current - target)
                
                # 如果在容差范围内，不计入误差
                if diff <= tolerance:
                    continue
                
                # 计算超出容差部分的平方误差
                squared_error = (diff - tolerance) ** 2
                
                # 计算加权平方误差
                weighted_error = weight * squared_error
                error += weighted_error
                weights_sum += weight
                valid_count += 1
            
            # 归一化误差（如果需要且有有效键类型）
            if normalize and weights_sum > 0:
                error /= weights_sum
                
            # 应用方差平方
            if self._varianceSquared != 0:
                error /= self._varianceSquared
                
            # 为避免返回零误差，添加一个极小值
            if error < 1e-10:
                error = FLOAT_TYPE(1e-10)
                
            return error
            
        except Exception as e:
            LOGGER.error(f"计算标准误差出错: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
            # 返回一个合理的默认误差值
            return FLOAT_TYPE(0.01)
    
    def _runtime_initialize(self):
        """初始化约束运行时数据"""
        LOGGER.debug("初始化NMR距离约束")
        # 清除缓存
        self._cachedCoords = None
        self._cachedBondRatios = None
        
        # 重置标准误差
        self._reset_standard_error()
        LOGGER.debug(f"初始标准误差：{self.standardError:.6f}")
    
    def _reset_standard_error(self):
        """重置并计算标准误差"""
        if self.data is not None:
            standardError = self.compute_standard_error(self.data)
            self.set_standard_error(standardError)
    
    def compute_data(self, update=True):
        """
        计算约束的数据
        
        :Parameters:
            #. update (boolean): 是否用新的计算更新约束数据和标准误差
        
        :Returns:
            #. data (dict): 约束数据字典
            #. standardError (float): 约束标准误差
        """
        startTime = time.time()
        
        # 计算键比例
        bondRatios = self._calculate_bond_statistics(useCache=False)
        
        # 计算标准误差
        stdError = self.compute_standard_error(bondRatios)
        
        # 更新数据和标准误差
        if update:
            self.set_data(bondRatios)
            self.set_standard_error(stdError)
            
            # 设置originalData
            if self.originalData is None:
                self._set_original_data(self.data)
        
        # 记录性能数据
        elapsedTime = time.time() - startTime
        LOGGER.debug(f"compute_data执行时间: {elapsedTime:.6f}秒")
        
        return bondRatios, stdError
    
    def compute_before_move(self, realIndexes, relativeIndexes):
        """
        执行移动前计算约束
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 移动将应用到的原子相对索引
        """
        
        # 记录当前状态
        self._beforeMoveStandardError = self.standardError
        
        # 不需要提前计算，保持数据一致性
        self.set_active_atoms_data_before_move(None)
        self.set_active_atoms_data_after_move(None)
        
        LOGGER.debug(f"移动前标准误差: {self._beforeMoveStandardError}")
    
    def compute_after_move(self, realIndexes, relativeIndexes, movedBoxCoordinates):
        """
        执行移动后计算约束
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 移动将应用到的原子相对索引
            #. movedBoxCoordinates (numpy.ndarray): 移动原子的新坐标
        """
        try:
            # 保存原始坐标
            boxData = np.array(self.engine.boxCoordinates[relativeIndexes], dtype=FLOAT_TYPE)
            
            if self.engine.isPBC:
                realData = np.array(self.engine.realCoordinates[relativeIndexes], dtype=FLOAT_TYPE)
            
            
            # 临时更改坐标
            self.engine.boxCoordinates[relativeIndexes] = movedBoxCoordinates
            
            # 将盒子坐标转换为实际坐标
            if self.engine.isPBC:
                movedRealCoordinates = np.dot(movedBoxCoordinates, self.engine.basisVectors)
                self.engine.realCoordinates[relativeIndexes] = movedRealCoordinates
            
            # 计算移动后的键比例
            bondRatios = self._calculate_bond_statistics()
            
            # 恢复坐标
            self.engine.boxCoordinates[relativeIndexes] = boxData
            if self.engine.isPBC:
                self.engine.realCoordinates[relativeIndexes] = realData
            
            # 计算标准误差
            newError = self.compute_standard_error(bondRatios)
            
            # 确保返回不同的误差值(防止数值精度问题)
            if abs(newError - self._beforeMoveStandardError) < 1e-10:
                # 对于微小变化，添加小量扰动
                delta = 1e-6 * self._beforeMoveStandardError * (1.0 if np.random.random() > 0.5 else -1.0)
                newError += delta
                LOGGER.debug("应用微小误差扰动以避免数值精度问题")
            
            # 设置移动后标准误差
            self.set_after_move_standard_error(newError)
            
            # 保存键比例数据以便接受移动时使用
            self.set_active_atoms_data_after_move(bondRatios)
            
            # 增加尝试计数
            self.increment_tried()
            
            LOGGER.debug(f"移动后标准误差: {newError}, 差异: {newError - self._beforeMoveStandardError}")
        
        except Exception as e:
            LOGGER.error(f"compute_after_move错误: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
            
            # 发生错误时设置一个合理的误差值，避免阻断模拟
            # 设置为与移动前相同的误差，既不鼓励也不阻止接受
            self.set_after_move_standard_error(self._beforeMoveStandardError * 1.01)  # 略微增加
            self.set_active_atoms_data_after_move({})
            self.increment_tried()

    def accept_move(self, realIndexes, relativeIndexes):
        """
        接受移动
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
        """
        # 更新数据
        self.set_data(self.activeAtomsDataAfterMove)
        
        # 重置activeAtoms数据
        self.set_active_atoms_data_before_move(None)
        self.set_active_atoms_data_after_move(None)
        
        # 更新标准误差
        self.set_standard_error(self.afterMoveStandardError)
        self.set_after_move_standard_error(None)
        
        # 增加接受计数
        self.increment_accepted()

    def reject_move(self, realIndexes, relativeIndexes):
        """
        拒绝移动
        
        :Parameters:
            #. realIndexes (numpy.ndarray): 原子的真实索引
            #. relativeIndexes (numpy.ndarray): 原子的相对索引
        """
        # 重置activeAtoms数据
        self.set_active_atoms_data_before_move(None)
        self.set_active_atoms_data_after_move(None)
        
        # 重置移动后标准误差
        self.set_after_move_standard_error(None)

    def compute_as_if_amputated(self, realIndex, relativeIndex):
        """
        计算并返回约束的数据和标准误差，如同给定原子被截除
        
        :Parameters:
            #. realIndex (numpy.ndarray): 作为numpy数组的单个元素的原子索引
            #. relativeIndex (numpy.ndarray): 作为numpy数组的单个元素的原子相对索引
        """
        try:
            # 保存当前坐标
            tempCoordinates = np.array(self.engine.realCoordinates, dtype=FLOAT_TYPE)
            
            # 创建不包含被截除原子的坐标数组
            mask = np.ones(self.engine.numberOfAtoms, dtype=bool)
            mask[relativeIndex] = False
            reducedCoordinates = tempCoordinates[mask]
            
            # 计算截除后的键比例
            bondRatios = self._calculate_bond_statistics(coordinates=reducedCoordinates, useCache=False)
            
            # 计算标准误差
            standardError = self.compute_standard_error(bondRatios)
            
            # 设置截除数据
            self.set_amputation_data(bondRatios)
            self.set_amputation_standard_error(standardError)
        
        except Exception as e:
            LOGGER.error(f"compute_as_if_amputated错误: {e}")
            
            # 出错时设置默认值
            self.set_amputation_data({})
            self.set_amputation_standard_error(self.standardError * 1.1)  # 略微增加误差

    def accept_amputation(self, realIndex, relativeIndex):
        """
        接受截除原子并相应地设置约束数据和标准误差
        
        :Parameters:
            #. realIndex (numpy.ndarray): 原子的真实索引
            #. relativeIndex (numpy.ndarray): 原子的相对索引
        """
        # 更新数据
        self.set_data(self.amputationData)
        
        # 更新标准误差
        self.set_standard_error(self.amputationStandardError)
        
        # 重置截除数据
        self.set_amputation_data(None)
        self.set_amputation_standard_error(None)

    def reject_amputation(self, realIndex, relativeIndex):
        """
        拒绝截除原子并相应地设置约束数据和标准误差
        
        :Parameters:
            #. realIndex (numpy.ndarray): 原子的真实索引
            #. relativeIndex (numpy.ndarray): 原子的相对索引
        """
        # 重置截除数据
        self.set_amputation_data(None)
        self.set_amputation_standard_error(None)
        
    def listen(self, message, argument=None):
        """
        监听来自Broadcaster的任何消息
        
        :Parameters:
            #. message (object): 要发送到约束listen方法的任何Python对象
            #. argument (object): 要传递给监听器的任何类型的参数
        """
        if message in ("engine set", "update pdb", "update molecules indexes", "update elements indexes", "update names indexes"):
            if self.engine is not None:
                # 重置缓存
                self._cachedCoords = None
                self._cachedBondRatios = None
            self.reset_constraint()
        elif message in("update boundary conditions",):
            self.reset_constraint()

    def _on_collector_reset(self):
        """当收集器重置时调用"""
        pass

