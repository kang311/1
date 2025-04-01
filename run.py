import os 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from fullrmc.Engine import Engine
from fullrmc.Constraints.NMRDistanceConstraints import NMRDistanceConstraint  # 导入优化后的约束
from fullrmc.Generators.Swaps import SwapPositionsGenerator
from fullrmc.Globals import FLOAT_TYPE
import logging
from datetime import datetime
import random
import traceback


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
                # 注意：这是一种变通方法，应尽可能使用官方API
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

class SimulationManager:
    """模拟管理类，处理模拟运行和结果分析 - 仅使用NMR约束"""
    
    def __init__(self, dir_path, fresh_start=True, repeat_runs=20):
        self.dir_path = dir_path
        self.fresh_start = fresh_start
        self.repeat_runs = repeat_runs
        self.engine = None
        
        # 文件路径
        self.pdb_file = os.path.join(dir_path, "smoke_coal.pdb")
        self.gr_file = os.path.join(dir_path, "smoke_coal.gr")
        self.engine_file = os.path.join(dir_path, "coal_engine.rmc")
        
        # 验证文件
        self._validate_files()
        
    def _validate_files(self):
        """验证所需文件是否存在"""
        for f in [self.pdb_file]:  # 只检查PDB文件，不需要gr文件
            if not os.path.isfile(f):
                raise FileNotFoundError(f"文件 {f} 不存在")
            logger.info(f"文件 {f} 已找到")

    def preprocess_pdb_for_aromaticity(self):
        """预处理PDB文件以增强芳香性检测"""
        try:
            logger.info("开始预处理PDB文件以增强芳香性检测...")
        
            # 检查PDB文件是否存在
            if not os.path.isfile(self.pdb_file):
                logger.error(f"PDB文件不存在: {self.pdb_file}")
                return False
        
            # 创建预处理后的临时文件路径
            processed_pdb = os.path.join(self.dir_path, "processed_smoke_coal.pdb")
        
            from openbabel import openbabel as ob
        
            # 创建OpenBabel转换器
            obConversion = ob.OBConversion()
            obConversion.SetInAndOutFormats("pdb", "pdb")
        
            # 创建分子对象
            mol = ob.OBMol()
        
            # 读取PDB文件
            if not obConversion.ReadFile(mol, self.pdb_file):
                logger.error(f"无法读取PDB文件: {self.pdb_file}")
                return False
        
            logger.info(f"PDB原始分子: {mol.NumAtoms()}原子, {mol.NumBonds()}键")
        
            # 使用增强的芳香性感知
            aromatic_typer = ob.OBAromaticTyper()
        
            # 先重置原子和键的芳香标志
            for atom in ob.OBMolAtomIter(mol):
                atom.SetAromatic(False)
        
            for bond in ob.OBMolBondIter(mol):
                bond.SetAromatic(False)
        
            # 增强芳香性感知
            aromatic_typer.BeginMIDO(20)  # 使用更高的MIDO阈值提高灵敏度
            aromatic_typer.AssignAromaticFlags(mol)
        
            # 感知键序
            mol.PerceiveBondOrders()
        
            # 计算芳香键数量
            aromatic_bonds = 0
            for bond in ob.OBMolBondIter(mol):
                if bond.IsAromatic():
                    aromatic_bonds += 1
        
            logger.info(f"增强芳香性感知后: 检测到{aromatic_bonds}个芳香键")
        
            # 使用力场优化结构
            ff = ob.OBForceField.FindForceField("MMFF94")
            if ff:
                logger.info("应用MMFF94力场优化分子构型...")
                ff.Setup(mol)
                ff.SteepestDescent(500)  # 执行500步最陡下降优化
                ff.GetCoordinates(mol)
                logger.info("力场优化完成")
            else:
                logger.warning("找不到MMFF94力场，跳过力场优化")
        
            # 再次进行芳香性感知
            aromatic_typer.AssignAromaticFlags(mol)
            mol.PerceiveBondOrders()
        
            # 再次计算芳香键数量
            aromatic_bonds = 0
            for bond in ob.OBMolBondIter(mol):
                if bond.IsAromatic():
                    aromatic_bonds += 1
        
            logger.info(f"优化后再次检测: {aromatic_bonds}个芳香键")
        
            # 保存预处理后的PDB文件
            obConversion.WriteFile(mol, processed_pdb)
            logger.info(f"预处理后的PDB文件已保存到: {processed_pdb}")
        
            # 更新PDB文件路径为预处理后的文件
            self.pdb_file = processed_pdb
        
            return True
    
        except Exception as e:
            logger.error(f"预处理PDB文件失败: {e}")
            logger.error(traceback.format_exc())
            return False

    
    def initialize_engine(self):
        """初始化和配置模拟引擎 - 增强错误处理和状态验证"""
        try:
            logger.info(f"尝试初始化引擎，fresh_start={self.fresh_start}")
        
            # 检查引擎文件是否存在
            engine_exists = os.path.isfile(self.engine_file)
            logger.info(f"引擎文件存在: {engine_exists}")
        
            if not engine_exists or self.fresh_start:
                # 预处理PDB文件增强芳香性检测
                self.preprocess_pdb_for_aromaticity()
            
                logger.info("创建新的引擎...")
                self.engine = Engine(path=self.engine_file, freshStart=True)
            
                # 设置PDB
                logger.info(f"设置PDB文件: {self.pdb_file}")
                self.engine.set_pdb(self.pdb_file)
            
                # 验证引擎状态
                if not self._verify_engine_ready():
                    logger.error("引擎初始化失败，请检查PDB文件格式")
                    return False
                
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
            logger.info(f"引擎信息: 原子数: {len(self.engine.allElements)}")
            
            return True
        except Exception as e:
            logger.error(f"引擎验证错误: {e}")
            return False
    
    def _setup_constraints(self):
        """设置和配置所有约束 - 增强错误处理"""
        try:
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未定义，无法添加约束")
                return False
                
            # 创建约束实例 - 使用优化后的NMRDistanceConstraint
            logger.info("创建NMR约束实例...")
            
            # 设置初始键比例目标
            target_ratios = {
                'C-C_SINGLE': 0.45, 
                'C-C_DOUBLE': 0.25,
                'C-C_TRIPLE': 0.05,
                'C-H_SINGLE': 0.20,
                'C-O_SINGLE': 0.05,
                'C-O_DOUBLE': 0.05
            }
            
            # 设置权重和容差
            weights = {
                'C-C_SINGLE': 1.0,
                'C-C_DOUBLE': 1.2,
                'C-C_TRIPLE': 1.3,
                'C-H_SINGLE': 1.1,
                'C-O_SINGLE': 1.1,
                'C-O_DOUBLE': 1.1
            }
            
            tolerances = {
                'C-C_SINGLE': 0.05,
                'C-C_DOUBLE': 0.03,
                'C-C_TRIPLE': 0.02,
                'C-H_SINGLE': 0.03,
                'C-O_SINGLE': 0.03,
                'C-O_DOUBLE': 0.03
            }
            
            # 创建NMR约束实例
            nmr_constraint = NMRDistanceConstraint(
                targetBondRatios=target_ratios,
                weights=weights,
                tolerances=tolerances,
                varianceSquared=0.001
            )
            
            # 添加约束到引擎 - 修正这里的约束添加方式
            logger.info("添加约束到引擎...")
            self.engine.add_constraints([nmr_constraint])  # 直接添加约束实例列表
            
            # 验证约束添加成功并计算初始数据
            if len(self.engine.constraints) > 0:
                for constraint in self.engine.constraints:
                    if isinstance(constraint, NMRDistanceConstraint):
                        logger.info("计算约束初始数据...")
                        data, error = constraint.compute_data(update=True)
                    
                        # 检查芳香键比例
                        if 'AROMATIC' in data:
                            logger.info(f"初始芳香键比例: {data['AROMATIC']:.6f}")
                        else:
                            logger.warning("计算结果中不包含芳香键比例")
                    
                        # 输出所有键类型的比例
                        for bond_type, ratio in data.items():
                            logger.info(f"键类型 {bond_type}: {ratio:.6f}")
                    
                        logger.info(f"约束标准误差: {error}")
                return True
            else:
                logger.error("没有约束被成功添加")
                return False
    
        except Exception as e:
            logger.error(f"约束设置失败: {e}")
            logger.error(traceback.format_exc())
            return False

    def get_bond_counts_from_constraint(self, constraint):
        """获取NMR约束中的键数量统计（绝对数值）
    
        :Parameters:
            #. constraint (NMRDistanceConstraint): NMR约束实例
        
        :Returns:
            #. bond_counts (dict): 键类型及其数量
            #. total_bonds (int): 总键数
        """
        try:
            # 确保是NMR约束
            if not isinstance(constraint, NMRDistanceConstraint):
                logger.warning("无法获取键计数：提供的约束不是NMR约束")
                return {}, 0

            # 从引擎获取坐标和元素信息
            coords = self.engine.realCoordinates
            elements = self.engine.allElements

            # 使用OpenBabel计算键统计
            # 创建OpenBabel分子对象
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
            logger.info("=" * 70)
        
        except Exception as e:
            logger.error(f"打印键统计信息出错: {e}")
            logger.error(traceback.format_exc())
    
    def run_simulation(self, nsteps=250000, saveFrequency=10000, xyzFrequency=5000):
        """改进的模拟运行方法 - 采用分阶段策略以提高有效性"""
        try:
            logger.info("开始分阶段模拟...")
            
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未初始化，无法运行模拟")
                return False
                
            # 强制初始化所有约束
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    
                    # 计算数据并更新
                    logger.info("重新计算约束数据...")
                    data, error = constraint.compute_data(update=True)
                    
                    # 验证计算结果
                    logger.info(f"计算的标准误差: {error}")
                    if error is None or error < 0.0001:
                        error = 0.2
                        logger.warning(f"约束计算返回无效误差，使用默认值: {error}")
                        constraint.set_standard_error(error)
                    
                    # 测试一次移动计算
                    try:
                        logger.info("测试约束移动计算...")
                        indexes = np.array([0], dtype=np.int)
                        constraint.compute_before_move(indexes, indexes)
                        coords = self.engine.realCoordinates[indexes].copy()
                        coords[0,0] += 0.1  # 微调第一个原子的x坐标
                        constraint.compute_after_move(indexes, indexes, coords)
                        
                        before_err = constraint.standardError
                        after_err = constraint.afterMoveStandardError
                        logger.info(f"测试移动结果: 移动前={before_err}, 移动后={after_err}")
                        
                        if abs(after_err - before_err) < 1e-10:
                            logger.warning("约束移动计算可能无效!")
                    except Exception as test_e:
                        logger.error(f"测试约束移动失败: {test_e}")

            logger.info("尝试使用增强的芳香性检测方法...")
            try:
                coords = self.engine.realCoordinates
                elements = self.engine.allElements
    
                # 使用增强的芳香性检测
                aromatic_bonds, bond_types = enhance_aromaticity.enhanced_aromatic_detection(coords, elements)
    
                if aromatic_bonds:
                    logger.info(f"增强检测发现 {len(aromatic_bonds)} 个芳香键")
        
                    # 获取NMR约束
                    for constraint in self.engine.constraints:
                        if isinstance(constraint, NMRDistanceConstraint):
                            # 计算当前数据
                            current_data = constraint.get_constraint_value()
                
                            # 如果没有芳香键，手动添加
                            if 'AROMATIC' not in current_data or current_data['AROMATIC'] < 0.01:
                                # 计算总键数
                                total_bonds = sum(current_data.values()) if current_data else 0
                                # 添加芳香键比例
                                aromatic_ratio = len(aromatic_bonds) / (total_bonds + len(aromatic_bonds))
                    
                                # 创建新数据
                                new_data = current_data.copy() if current_data else {}
                                new_data['AROMATIC'] = aromatic_ratio
                    
                                # 归一化其他键比例
                                total = sum(new_data.values())
                                if total > 0:
                                    for key in new_data:
                                        new_data[key] /= total
                    
                                # 更新约束数据
                                constraint.set_data(new_data)
                                logger.info(f"手动设置芳香键比例为: {new_data['AROMATIC']:.6f}")
                    
                                # 重新计算标准误差
                                new_error = constraint.compute_standard_error(new_data)
                                constraint.set_standard_error(new_error)
                                logger.info(f"更新后的标准误差: {new_error}")
            except Exception as e:
                logger.error(f"增强芳香性检测失败: {e}")
                logger.error(traceback.format_exc())
            
            # 第1阶段: 非常小的幅度，重点强调接受率
            logger.info("模拟第1阶段: 极小幅度热身 (5000步)")
            total_steps = 0
            
            # 设置小幅度
            for group in self.engine.groups:
                group.moveGenerator.set_amplitude(0.05)
                
            # 运行初始阶段
            self.engine.run(numberOfSteps=50000, saveFrequency=1000)
            total_steps += 50000
            
            # 检查第一阶段的接受率
            accept_rate = self.engine.accepted / max(1, self.engine.tried)
            logger.info(f"第1阶段接受率: {accept_rate:.4f} ({self.engine.accepted}/{self.engine.tried})")

            # 打印第一阶段的键统计信息
            self.print_bond_statistics(stage_name="第1阶段")
            
            # 如果接受率过高，增加拒绝概率
            if accept_rate > 0.9:
                logger.warning("第1阶段接受率过高，增加拒绝概率")
                for constraint in self.engine.constraints:
                    if isinstance(constraint, NMRDistanceConstraint):
                        constraint.set_reject_probability(0.4)
            
            # 第2阶段: 中等幅度
            logger.info("模拟第2阶段: 中等幅度 (15000步)")
            for group in self.engine.groups:
                group.moveGenerator.set_amplitude(0.15)
                
            # 保存第1阶段状态
            self.engine.save(path=os.path.join(self.dir_path, "stage1_checkpoint.rmc"))
            
            # 运行第2阶段
            self.engine.run(numberOfSteps=50000, saveFrequency=5000)
            total_steps += 50000
            
            # 检查第二阶段的接受率
            stage2_accepted = self.engine.accepted - self.engine.constraints[0].accepted
            stage2_tried = self.engine.tried - self.engine.constraints[0].tried
            accept_rate = stage2_accepted / max(1, stage2_tried)
            logger.info(f"第2阶段接受率: {accept_rate:.4f} ({stage2_accepted}/{stage2_tried})")

            # 打印第二阶段的键统计信息
            self.print_bond_statistics(stage_name="第2阶段")
            
            # 调整拒绝概率
            target_accept = 0.4  # 目标接受率
            if accept_rate > 0.7:
                # 接受率过高，增加拒绝概率
                for constraint in self.engine.constraints:
                    if isinstance(constraint, NMRDistanceConstraint):
                        current_reject = constraint.rejectProbability  # 使用属性访问器
                        new_reject = min(0.7, current_reject * 1.5)
                        logger.info(f"增加拒绝概率: {current_reject:.2f} -> {new_reject:.2f}")
                        constraint.set_reject_probability(new_reject)
            elif accept_rate < 0.2:
                # 接受率过低，减少拒绝概率
                for constraint in self.engine.constraints:
                    if isinstance(constraint, NMRDistanceConstraint):
                        current_reject = constraint.rejectProbability  # 使用属性访问器
                        new_reject = max(0.1, current_reject * 0.8)
                        logger.info(f"减少拒绝概率: {current_reject:.2f} -> {new_reject:.2f}")
                        constraint.set_reject_probability(new_reject)
            
            # 保存第2阶段状态
            self.engine.save(path=os.path.join(self.dir_path, "stage2_checkpoint.rmc"))
            
            # 第3阶段: 正常幅度，分批运行以便监控进度
            remaining_steps = nsteps - total_steps
            logger.info(f"模拟第3阶段: 正常幅度 (共{remaining_steps}步)")
            
            # 设置正常幅度
            for group in self.engine.groups:
                group.moveGenerator.set_amplitude(0.25)
                
            # 批量运行，每批次5000步
            batch_size = 50000
            num_batches = max(1, remaining_steps // batch_size)
            
            for batch in range(num_batches):
                batch_steps = min(batch_size, remaining_steps - batch * batch_size)
                if batch_steps <= 0:
                    break
                    
                logger.info(f"运行第3阶段批次 {batch+1}/{num_batches} ({batch_steps}步)")
                
                # 运行一批次
                self.engine.run(numberOfSteps=batch_steps, 
                                saveFrequency=saveFrequency,
                                xyzFrequency=xyzFrequency)
                
                # 更新总步数
                total_steps += batch_steps
                
                # 检查接受率
                batch_accept = self.engine.accepted / max(1, self.engine.tried)
                logger.info(f"当前累计接受率: {batch_accept:.4f} ({self.engine.accepted}/{self.engine.tried})")

                # 打印当前批次的键统计信息
                self.print_bond_statistics(stage_name=f"第3阶段批次{batch+1}")
                
                # 动态调整拒绝概率
                if batch < num_batches - 1:  # 除了最后一个批次
                    for constraint in self.engine.constraints:
                        if isinstance(constraint, NMRDistanceConstraint):
                            current_reject = constraint.rejectProbability  # 使用属性访问器
                            if batch_accept > 0.7:
                                # 接受率过高，稍微增加拒绝概率
                                new_reject = min(0.7, current_reject * 1.1)
                                logger.info(f"调整拒绝概率: {current_reject:.2f} -> {new_reject:.2f}")
                                constraint.set_reject_probability(new_reject)
                            elif batch_accept < 0.3:
                                # 接受率过低，稍微减少拒绝概率
                                new_reject = max(0.1, current_reject * 0.9)
                                logger.info(f"调整拒绝概率: {current_reject:.2f} -> {new_reject:.2f}")
                                constraint.set_reject_probability(new_reject)
                
                # 每个批次后保存一次
                checkpoint_path = os.path.join(self.dir_path, f"stage3_batch{batch+1}_checkpoint.rmc")
                self.engine.save(path=checkpoint_path)
                
                # 分析当前约束状态
                self.analyze_nmr_constraints(output_file=os.path.join(self.dir_path, f"nmr_analysis_batch{batch+1}.txt"))
                
                # 生成比较图
                self.plot_nmr_comparison(output_file_path=os.path.join(self.dir_path, f'nmr_comparison_batch{batch+1}.png'))
            
            # 保存最终结果
            final_path = os.path.join(self.dir_path, "final_result.rmc")
            self.engine.save(path=final_path)
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
        """改进的NMR约束分析方法 - 适用于优化后的约束"""
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
                    target_ratios = constraint.targetBondRatios  # 使用属性访问
                    current_ratios = data  # 数据直接存储在data中
                    
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
                        tolerance = constraint.tolerances.get(bond_type, 0.05)  # 使用属性访问
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
        """绘制改进的NMR键类型比例比较图 - 适用于优化后的约束"""
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
                    target_ratios = constraint.targetBondRatios  # 使用属性访问
                    current_ratios = data  # 数据直接存储在data中
                    
                    # 准备绘图数据
                    bond_types = list(target_ratios.keys())
                    target_values = [target_ratios[bond] for bond in bond_types]
                    current_values = [current_ratios.get(bond, 0.0) for bond in bond_types]
                    
                    # 计算符合容差的键类型
                    compliant = []
                    for i, bond in enumerate(bond_types):
                        tolerance = constraint.tolerances.get(bond, 0.05)  # 使用属性访问
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
                    tolerances = [constraint.tolerances.get(bond, 0.05) for bond in bond_types]  # 使用属性访问
                    ax.errorbar(x - width/2, target_values, yerr=tolerances, fmt='none', ecolor='darkred', capsize=5)
                    
                    # 添加图表元素
                    ax.set_xlabel('Key Type', fontsize=12)
                    ax.set_ylabel('Key Ratio', fontsize=12)
                    
                    # 计算符合率
                    compliant_count = sum(compliant)
                    compliance_rate = compliant_count / len(bond_types) * 100
                    
                    # 设置标题，包含标准误差和符合率
                    title = f"NMR Key Ratio(err: {constraint.standardError:.4f}, Pass rate: {compliance_rate:.1f}%)"
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
            
    def analyze_simulation_performance(self):
        """分析模拟性能"""
        try:
            # 提取约束标准误差
            errors = {}
            for constraint in self.engine.constraints:
                class_name = constraint.__class__.__name__
                errors[class_name] = constraint.standardError
                
            # 计算接受率
            total_tried = self.engine.tried
            total_accepted = self.engine.accepted
            acceptance_rate = float(total_accepted) / max(1, total_tried)
            
            logger.info("=" * 50)
            logger.info("模拟性能分析")
            logger.info("=" * 50)
            logger.info(f"模拟总步数：{total_tried}")
            logger.info(f"接受的步数：{total_accepted}")
            logger.info(f"总接受率：{acceptance_rate:.4f} ({acceptance_rate*100:.2f}%)")
            
            # 每个约束的接受/拒绝统计
            for constraint in self.engine.constraints:
                class_name = constraint.__class__.__name__
                constraint_tried = constraint.tried
                constraint_accepted = constraint.accepted
                constraint_rate = float(constraint_accepted) / max(1, constraint_tried) if constraint_tried > 0 else 0
                logger.info(f"{class_name} - 尝试/接受/率：{constraint_tried}/{constraint_accepted}/{constraint_rate:.4f}")
                
                # 更详细的约束状态
                if isinstance(constraint, NMRDistanceConstraint):
                    target_ratios = constraint.targetBondRatios  # 使用属性访问
                    current_data = constraint.get_constraint_value()  # 使用公共方法获取数据
                    
                    # 计算符合目标的键类型数量
                    compliant_count = 0
                    bond_infos = []
                    
                    for bond_type, target in target_ratios.items():
                        current = current_data.get(bond_type, 0.0)
                        tolerance = constraint.tolerances.get(bond_type, 0.05)  # 使用属性访问
                        diff = abs(current - target)
                        is_compliant = diff <= tolerance
                        
                        if is_compliant:
                            compliant_count += 1
                        
                        bond_infos.append({
                            'type': bond_type,
                            'target': target,
                            'current': current,
                            'diff': diff,
                            'compliant': is_compliant
                        })
                    
                    # 计算总体符合率
                    total_types = len(target_ratios)
                    compliance_rate = compliant_count / total_types * 100
                    logger.info(f"键类型符合率: {compliant_count}/{total_types} ({compliance_rate:.1f}%)")
                    
                    # 显示最大偏差的前三个键类型
                    sorted_bonds = sorted(bond_infos, key=lambda x: x['diff'], reverse=True)
                    logger.info("偏差最大的键类型:")
                    for i, bond in enumerate(sorted_bonds[:3]):
                        logger.info(f"  {i+1}. {bond['type']}: 目标={bond['target']:.4f}, 当前={bond['current']:.4f}, 差异={bond['diff']:.4f}")
            
            # 生成性能报告文件
            report_file = os.path.join(self.dir_path, "simulation_performance.txt")
            with open(report_file, 'w') as f:
                f.write("=" * 50 + "\n")
                f.write("模拟性能分析\n")
                f.write("=" * 50 + "\n")
                f.write(f"模拟总步数：{total_tried}\n")
                f.write(f"接受的步数：{total_accepted}\n")
                f.write(f"总接受率：{acceptance_rate:.4f} ({acceptance_rate*100:.2f}%)\n\n")
                
                for constraint in self.engine.constraints:
                    class_name = constraint.__class__.__name__
                    constraint_tried = constraint.tried
                    constraint_accepted = constraint.accepted
                    constraint_rate = float(constraint_accepted) / max(1, constraint_tried) if constraint_tried > 0 else 0
                    f.write(f"{class_name}:\n")
                    f.write(f"  尝试次数: {constraint_tried}\n")
                    f.write(f"  接受次数: {constraint_accepted}\n")
                    f.write(f"  接受率: {constraint_rate:.4f} ({constraint_rate*100:.2f}%)\n")
                    f.write(f"  标准误差: {constraint.standardError}\n\n")
            
            logger.info(f"性能报告已保存到: {report_file}")
            return {
                'total_tried': total_tried,
                'total_accepted': total_accepted,
                'acceptance_rate': acceptance_rate,
                'errors': errors,
                'report_file': report_file
            }
            
        except Exception as e:
            logger.error(f"分析模拟性能出错: {str(e)}")
            logger.error(traceback.format_exc())
            return None

def get_file_paths(dir_path):
    """返回文件路径字典"""
    return {
        'pdb': os.path.join(dir_path, "smoke_coal.pdb"),
        'gr': os.path.join(dir_path, "smoke_coal.gr"),
        'engine': os.path.join(dir_path, "coal_engine.rmc")
    }

# 主程序
if __name__ == "__main__":
    # 记录当前时间
    current_time = datetime.now()
    logger.info(f"模拟开始时间: {current_time}")
    logger.info(f"用户: kang311")
    logger.info(f"日期: {current_time.strftime('%Y-%m-%d')}")
    
    DIR_PATH = os.path.dirname(os.path.realpath(__file__))
    file_paths = get_file_paths(DIR_PATH)
    
    # 根据需要选择是否重新开始
    fresh_start = True  # 设置为True创建新引擎
    
    try:
        simulation_manager = SimulationManager(DIR_PATH, fresh_start=fresh_start)
        
        # 初始化引擎
        if not simulation_manager.initialize_engine():
            logger.error("引擎初始化失败，程序终止")
            exit(1)
            
        logger.info("引擎初始化成功")
        
        # 执行多次模拟
        runs = simulation_manager.repeat_runs
        logger.info(f"计划执行 {runs} 次模拟")
        
        for run_idx in range(runs):
            if run_idx > 0:
                # 重置统计数据以便追踪每个运行的接受率
                reset_engine_constraints_statistics(simulation_manager.engine)
                logger.info(f"已重置约束统计数据，准备进行第 {run_idx + 1} 次模拟")
    
            logger.info(f"开始第 {run_idx + 1}/{runs} 次模拟")
            
            # 运行模拟
            success = simulation_manager.run_simulation(nsteps=250000, 
                                                     saveFrequency=10000, 
                                                     xyzFrequency=5000)
            
            if not success:
                logger.error(f"第 {run_idx + 1} 次模拟失败")
                continue
                
            # 每次运行后分析结果
            logger.info(f"第 {run_idx + 1} 次模拟完成，正在分析结果...")
            
            # 分析约束
            simulation_manager.analyze_nmr_constraints(
                output_file=os.path.join(DIR_PATH, f"nmr_analysis_run{run_idx+1}.txt"))
                
            # 绘制比较图
            simulation_manager.plot_nmr_comparison(
                output_file_path=os.path.join(DIR_PATH, f"nmr_comparison_run{run_idx+1}.png"))
                
            # 导出当前结构
            export_path = os.path.join(DIR_PATH, f"structure_run{run_idx+1}.pdb")
            simulation_manager.export_final_pdb(output_file=export_path)
            
            # 分析性能
            simulation_manager.analyze_simulation_performance()
            
        # 所有运行完成后创建总结
        logger.info("所有模拟完成，创建总结报告")
        summary_file = simulation_manager.create_simulation_summary()
        
        if summary_file:
            logger.info(f"总结报告已生成：{summary_file}")
            
        # 导出最终结构
        final_pdb = simulation_manager.export_final_pdb()
        
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
