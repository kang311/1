"""
NMRConstraints包含基于核磁共振实验数据的约束类。
"""

# 标准库导入
from __future__ import print_function
import itertools, inspect, copy, os, re
import time
import logging
import traceback
from collections import defaultdict
# 外部库导入
import numpy as np
from openbabel import openbabel as ob
import matplotlib
matplotlib.use('Agg')  # 非交互式后端
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# fullrmc导入
from ..Globals import INT_TYPE, FLOAT_TYPE, PRECISION, LOGGER
from ..Globals import str, long, unicode, bytes, basestring, range, xrange, maxint
from ..Core.Collection import is_number, is_integer, get_path
from ..Core.Constraint import Constraint, ExperimentalConstraint

class BondRatioConstraint(Constraint):
    """
    通过核磁共振(NMR)数据控制分子中不同类型键的分布。
    该约束通过计算和控制不同键类型的比例将NMR实验数据与模型关联起来。
    
    :Parameters:
        #. targetBondRatios (dict): 目标键比例字典
    """
    
    def __init__(self, targetBondRatios=None):
        # 初始化基类
        super(BondRatioConstraint, self).__init__()
        
        # 添加类名属性以便于引擎自动识别
        self.name = "BondRatioConstraint"
        
        # 设置目标键比例
        self.set_target_bond_ratios(targetBondRatios)
                               
        # 初始化键类型映射
        self._initialize_bond_types_mapping()
        
        # 设置框架数据
        FRAME_DATA = [d for d in self.FRAME_DATA]
        FRAME_DATA.extend(['_BondRatioConstraint__targetBondRatios'])
    
        object.__setattr__(self, 'FRAME_DATA', tuple(FRAME_DATA))
    
    def __del__(self):
        """
        析构函数，确保清理所有资源
        """
        pass
    
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

        # 添加更多芳香键的映射
        aromatic_mappings = {
            'c:c': 'AROMATIC',
            'c1ccccc1': 'AROMATIC',
            'c1cccc1': 'AROMATIC',
            'c1ccc1': 'AROMATIC',
            'C:C': 'AROMATIC',
            'C:N': 'AROMATIC',
            'N:C': 'AROMATIC',
            'c:n': 'AROMATIC',
            'n:c': 'AROMATIC',
            # 添加所有可能的芳香环组合
        }
    
        # 更新映射字典
        self.__bondTypesMapping.update(aromatic_mappings)
        
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
        
        # 目标键比例总和应精确等于1.0
        total = sum(targetBondRatios.values())
        if abs(total - 1.0) > 1e-6:
            LOGGER.warning(f"键比例总和({total})不为1.0，这可能导致不精确结果")
        
        self.__targetBondRatios = {k: FLOAT_TYPE(v) for k, v in targetBondRatios.items()}
        
        # 记录到存储库
        self._dump_to_repository({'_BondRatioConstraint__targetBondRatios': self.__targetBondRatios})
        
    @property
    def targetBondRatios(self):
        """目标键比例字典"""
        return self.__targetBondRatios   
    
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
    
    def _calculate_bond_statistics(self, coordinates=None):
        """
        使用OpenBabel计算键统计数据
        
        :Parameters:
            #. coordinates (None, numpy.ndarray): 计算键统计的坐标。如果为None，
               使用引擎的当前坐标
        
        :Returns:
            #. bondRatios (dict): 键类型比例字典
        """
        # 使用引擎坐标或提供的坐标
        coords = self.engine.realCoordinates if coordinates is None else coordinates     
        
        try:
            # 计算键比例和计数
            bondRatios, bondCounts, totalBonds = self._compute_bond_ratios_with_openbabel(coords)
            
            # 存储计数信息供后续使用
            self._lastBondCounts = bondCounts
            self._lastTotalBonds = totalBonds
                        
            return bondRatios
            
        except Exception as e:
            LOGGER.error(f"键统计计算错误: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
            # 返回默认值而非引发异常，提高稳健性
            return {k: FLOAT_TYPE(0.0) for k in self.__targetBondRatios.keys()}
    
    def _compute_bond_ratios_with_openbabel(self, coordinates):
        """
        使用OpenBabel计算键比例（忽略kekulize警告）
        
        :Parameters:
            #. coordinates (numpy.ndarray): 原子坐标
            
        :Returns:
            #. bondRatios (dict): 键类型比例字典
            #. bondCounts (dict): 键类型计数字典
            #. totalBonds (int): 键总数
        """
        try:
            # 初始化所有键类型的计数器
            bondCounts = {bondType: 0 for bondType in self.__targetBondRatios.keys()}
            totalBonds = 0
            
            # 临时禁用OpenBabel警告输出
            old_warn_level = ob.obErrorLog.GetOutputLevel()
            ob.obErrorLog.SetOutputLevel(ob.obError)
            
            try:
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
                
                # 执行键序感知
                mol.PerceiveBondOrders()
                
                # 执行芳香性检测
                aromatic_typer = ob.OBAromaticTyper()
                aromatic_typer.AssignAromaticFlags(mol)
                
                # 分析键类型
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
                    
                    # 获取键类型
                    bond_order = bond.GetBondOrder()
                    is_aromatic = bond.IsAromatic()
                    
                    # 构建键类型
                    bond_type = None
                    
                    # 处理芳香键
                    if is_aromatic:
                        # 检查是否支持芳香键类型
                        if 'AROMATIC' in self.__targetBondRatios:
                            bond_type = 'AROMATIC'
                        else:
                            # 如果不支持AROMATIC类型，根据元素类型映射到其他键类型
                            # 按字母顺序排序元素
                            if begin_element > end_element:
                                begin_element, end_element = end_element, begin_element
                            
                            # 对于芳香C-C键，映射到C-C_DOUBLE
                            if begin_element == 'C' and end_element == 'C' and 'C-C_DOUBLE' in self.__targetBondRatios:
                                bond_type = 'C-C_DOUBLE'
                            # 对于其他芳香键，使用单键作为默认值
                            elif f"{begin_element}-{end_element}_SINGLE" in self.__targetBondRatios:
                                bond_type = f"{begin_element}-{end_element}_SINGLE"
                    else:
                        # 按字母顺序排序元素
                        if begin_element > end_element:
                            begin_element, end_element = end_element, begin_element
                        
                        # 直接构建标准键类型
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
                
                # 计算键比例
                bondRatios = {}
                if totalBonds > 0:
                    for bondType, count in bondCounts.items():
                        bondRatios[bondType] = FLOAT_TYPE(count / totalBonds)
                else:
                    # 如果没有找到键，返回全零比例
                    bondRatios = {bondType: FLOAT_TYPE(0.0) for bondType in bondCounts.keys()}
                
            finally:
                # 恢复原始警告级别
                ob.obErrorLog.SetOutputLevel(old_warn_level)
            
            # 返回比例、计数和总数
            return bondRatios, bondCounts, totalBonds
        
        except Exception as e:
            LOGGER.error(f"OpenBabel键比例计算错误: {e}")
            import traceback
            LOGGER.debug(traceback.format_exc())
            # 返回默认值
            empty_counts = {bondType: 0 for bondType in self.__targetBondRatios.keys()}
            return {bondType: FLOAT_TYPE(0.0) for bondType in self.__targetBondRatios.keys()}, empty_counts, 0
    
    def compute_standard_error(self, modelData):
        """
        计算标准误差
        
        :Parameters:
            #. modelData (dict): 与目标值比较的键比例字典
        
        :Returns:
            #. standardError (float): 计算的约束标准误差
        """
        # 初始化计算所需变量
        numerator = FLOAT_TYPE(0.0)     # 加权平方差
    
        # 对每种键类型计算误差
        for bond_type, target in self.__targetBondRatios.items():
            # 获取当前模拟的键比例
            current = modelData.get(bond_type, 0.0)
                
            # 累加加权平方差
            numerator += (current - target) ** 2       
    
        return numerator 
    
    def _runtime_initialize(self):
        """初始化约束运行时数据"""
        LOGGER.debug("初始化键比例约束")
        
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
        bondRatios = self._calculate_bond_statistics()
        
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
        
        try:
            # 记录当前状态
            self._beforeMoveStandardError = self.standardError
        
            # 计算当前键比例并存储
            current_ratios = self._calculate_bond_statistics()
        
            # 设置移动前数据
            self.set_active_atoms_data_before_move(current_ratios)
            self.set_active_atoms_data_after_move(None)
        
            LOGGER.debug(f"移动前标准误差: {self._beforeMoveStandardError}")
    
        except Exception as e:
            LOGGER.error(f"compute_before_move错误: {str(e)}")
            LOGGER.debug(traceback.format_exc())
    
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
        
        # 清除键计数缓存，确保下次获取时重新计算
        if hasattr(self, '_lastBondCounts'):
            delattr(self, '_lastBondCounts')
        if hasattr(self, '_lastTotalBonds'):
            delattr(self, '_lastTotalBonds')

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
        
    def listen(self, message, argument=None):
        """
        监听来自Broadcaster的任何消息
        
        :Parameters:
            #. message (object): 要发送到约束listen方法的任何Python对象
            #. argument (object): 要传递给监听器的任何类型的参数
        """
        if message in ("engine set", "update pdb", "update molecules indexes", "update elements indexes", "update names indexes"):
            if self.engine is not None:
                pass
            self.reset_constraint()
        elif message in("update boundary conditions",):
            self.reset_constraint()

    def _on_collector_reset(self):
        """当收集器重置时调用"""
        pass

    def get_bond_counts(self):
        """
        获取每种键类型的绝对计数
        
        :Returns:
            #. bondCounts (dict): 键类型及其计数的字典
            #. totalBonds (int): 键的总数
        """
        try:
            # 如果已经计算过，直接返回缓存的结果
            if hasattr(self, '_lastBondCounts') and hasattr(self, '_lastTotalBonds'):
                return self._lastBondCounts, self._lastTotalBonds
                
            # 如果没有缓存的结果，重新计算
            _, bondCounts, totalBonds = self._compute_bond_ratios_with_openbabel(self.engine.realCoordinates)
            
            # 缓存计算结果
            self._lastBondCounts = bondCounts
            self._lastTotalBonds = totalBonds
            
            return bondCounts, totalBonds
            
        except Exception as e:
            LOGGER.error(f"获取键计数出错: {e}")
            return {}, 0

    def plot(self, filename=None, figsize=(14, 10), dpi=300):
        """
        绘制键比例比较图，同时展示键的绝对计数
        
        :Parameters:
            #. filename (None, string): 输出文件路径。如果为None，则使用默认名称
            #. figsize (tuple): 图像大小 (宽, 高)
            #. dpi (int): 图像分辨率
            
        :Returns:
            #. filePath (string): 保存的文件路径
        """
        # 如果没有提供文件名，则使用默认名称
        if filename is None:
            filename = "bond_ratio_comparison.png"
        
        try:
            # 确保数据是最新的
            data = self.get_constraint_value()
            
            # 获取目标键比例和当前模型的键比例
            target_ratios = self.targetBondRatios
            current_ratios = data
            
            # 获取绝对键计数
            bond_counts, total_bonds = self.get_bond_counts()
            
            # 准备绘图数据
            bond_types = list(target_ratios.keys())
            target_values = [target_ratios[bond] for bond in bond_types]
            current_values = [current_ratios.get(bond, 0.0) for bond in bond_types]
            count_values = [bond_counts.get(bond, 0) for bond in bond_types]
                        
            # 创建子图布局: 两行一列
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios':[3, 1]})
            
            # 第一个子图: 比例比较
            x = np.arange(len(bond_types))
            width = 0.35

            colors = plt.cm.Set2(np.linspace(0,1,len(target_values)))
            
            # 绘制条形图 (比例)
            bar1 = ax1.bar(x - width/2, target_values, width, label='Target Ratio', color='lightcoral')            
                    
            bar2 = ax1.bar(x + width/2, current_values, width, label='Current Ratio', color=colors)
            
            # 添加图表元素
            ax1.set_ylabel('Bond Ratio', fontsize=12)     
   
            # 设置标题，包含标准误差和符合率
            title = f"NMR Bond Analysis (StdErr: {self.standardError:.4f})"
            ax1.set_title(title, fontsize=14)
            
            # 设置x轴标签 (比例子图)
            ax1.set_xticks(x)
            ax1.set_xticklabels(bond_types, rotation=45, ha='right', fontsize=10)
            
            # 添加数据标签 (比例)
            for i, value in enumerate(current_values):
                ax1.text(i + width/2, value + 0.01, f"{value:.3f}", 
                       ha='center', va='bottom', fontsize=9)
        
            # 添加目标值标签
            for i, value in enumerate(target_values):
                ax1.text(i - width/2, value + 0.01, f"{value:.3f}", 
                       ha='center', va='bottom', fontsize=9)
            
            # 第二个子图: 绝对计数
            bar3 = ax2.bar(x, count_values, width, color='goldenrod')
            
            # 设置第二个子图的标签
            ax2.set_xlabel('Bond Type', fontsize=12)
            ax2.set_ylabel('Absolute Count', fontsize=12)
            ax2.set_title(f"Bond Counts (Total: {total_bonds} bonds)", fontsize=12)
            
            # 设置x轴标签 (计数子图)
            ax2.set_xticks(x)
            ax2.set_xticklabels(bond_types, rotation=45, ha='right', fontsize=10)
            
            # 添加数据标签 (计数)
            for i, v in enumerate(count_values):
                percentage = (v / total_bonds * 100) if total_bonds > 0 else 0
                label = f"{v}\n({percentage:.1f}%)"
                ax2.text(i, v, label, ha='center', va='bottom', fontsize=8)
            
            # 自定义图例
            legend_elements = [
                Patch(facecolor='lightcoral', label='Target ratio'),
                Patch(facecolor='limegreen', label='Current (diff < 5%)'),
                Patch(facecolor='skyblue', label='Current (diff 5-10%)'),
                Patch(facecolor='salmon', label='Current (diff > 10%)')
            ]
            ax1.legend(handles=legend_elements, loc='upper right', fontsize=10)
            
            # 添加接受率信息
            accept_info = f"Acceptance: {self.accepted/max(1, self.tried):.2%} ({self.accepted}/{self.tried})"
            ax1.annotate(accept_info, xy=(0.02, 0.02), xycoords='axes fraction', 
                      bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.3),
                      fontsize=10)
            
            # 添加详细的键计数表格
            table_data = []
            headers = ["Bond Type", "Count", "Percentage", "Current Ratio", "Target Ratio", "Difference"]
            
            for i, bond in enumerate(bond_types):
                count = count_values[i]
                percentage = (count / total_bonds * 100) if total_bonds > 0 else 0
                current = current_values[i]
                target = target_values[i]
                diff = current - target
                diff_percent = (diff / max(target, 1e-10)) * 100
                
                table_data.append([
                    bond,
                    f"{count}",
                    f"{percentage:.1f}%",
                    f"{current:.3f}",
                    f"{target:.3f}",
                    f"{diff:.3f} ({diff_percent:+.1f}%)"
                ])
            
            # 创建表格文本内容
            table_text = "Bond Analysis Summary:\n\n"
            table_text += f"Total Bonds: {total_bonds}\n"
            table_text += f"Standard Error: {self.standardError:.6f}\n"
            
            # 创建表格格式
            col_widths = [max(len(row[i]) for row in [headers] + table_data) for i in range(len(headers))]
            
            # 添加表头
            header_line = "| "
            for i, header in enumerate(headers):
                header_line += header.ljust(col_widths[i]) + " | "
            table_text += header_line + "\n"
            
            # 添加分隔线
            separator = "| "
            for width in col_widths:
                separator += "-" * width + " | "
            table_text += separator + "\n"
            
            # 添加数据行
            for row in table_data:
                line = "| "
                for i, cell in enumerate(row):
                    line += cell.ljust(col_widths[i]) + " | "
                table_text += line + "\n"
            
            # 保存表格到文本文件
            table_filename = os.path.splitext(filename)[0] + "_table.txt"
            with open(table_filename, "w") as f:
                f.write(table_text)
            
            # 调整布局
            plt.tight_layout()
            
            # 保存图像
            plt.savefig(filename, dpi=dpi)
            plt.close()
            
            LOGGER.info(f"NMR比较图已保存: {filename}")
            LOGGER.info(f"键分析表格已保存: {table_filename}")
            return filename
            
        except Exception as e:
            LOGGER.error(f"绘制NMR比较图出错: {e}")
            LOGGER.error(traceback.format_exc())
            return None
