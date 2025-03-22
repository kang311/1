import os 
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from fullrmc.Engine import Engine
from fullrmc.Constraints.NMRDistanceConstraints import NMRDistanceConstraint
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.DistanceConstraints import InterMolecularDistanceConstraint
from fullrmc.Generators.Swaps import SwapPositionsGenerator
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Generators.Rotations import RotationGenerator
from fullrmc.Core.MoveGenerator import MoveGeneratorCollector
from fullrmc.Core.GroupSelector import RecursiveGroupSelector
from fullrmc.Selectors.RandomSelectors import RandomSelector
from fullrmc.Globals import FLOAT_TYPE, LOGGER
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
    """模拟管理类，处理模拟运行和结果分析"""
    
    def __init__(self, dir_path, fresh_start=True, repeat_runs=20, debug_level=1):
        self.dir_path = dir_path
        self.fresh_start = fresh_start
        self.repeat_runs = repeat_runs
        self.engine = None
        self._debug_level = debug_level  # 0=无调试输出，1=基本输出，2=详细输出
        
        # 文件路径
        self.pdb_file = os.path.join(dir_path, "smoke_coal.pdb")
        self.gr_file = os.path.join(dir_path, "smoke_coal.gr")
        self.engine_file = os.path.join(dir_path, "coal_engine.rmc")
        
        # 约束实例存储
        self.constraints = {}
        
        # 验证文件
        self._validate_files()
        
    def _debug(self, message, level=1):
        """调试日志输出"""
        if self._debug_level >= level:
            logger.info(f"[DEBUG] {message}")
        
    def _validate_files(self):
        """验证所需文件是否存在"""
        required_files = [self.pdb_file, self.gr_file]
        for f in required_files:
            if not os.path.isfile(f):
                raise FileNotFoundError(f"文件 {f} 不存在")
            self._debug(f"文件 {f} 已找到", level=1)
    
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
                    # 解包约束以便后续使用
                    self._unpack_constraints()
                    
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

    def _unpack_constraints(self):
        """从引擎中解包约束实例"""
        try:
            self.constraints = {}
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    self.constraints['nmr'] = constraint
                elif isinstance(constraint, PairDistributionConstraint):
                    self.constraints['pdf'] = constraint
                elif isinstance(constraint, InterMolecularDistanceConstraint):
                    self.constraints['dist'] = constraint
            
            logger.info(f"从引擎解包出 {len(self.constraints)} 个约束")
        except Exception as e:
            logger.error(f"约束解包失败: {e}")
            raise

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
        """设置和配置所有约束 - 添加PDF约束"""
        try:
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未定义，无法添加约束")
                return False
            
            self.constraints = {}  # 清空现有约束字典
                
            # 1. 创建NMR约束实例
            logger.info("创建NMR约束实例...")
            
            # 设置初始键比例目标
            target_ratios = {
                'C-C_SINGLE': 0.32, 
                'C-C_DOUBLE': 0.15,
                'C-C_TRIPLE': 0.08,
                'C-O_SINGLE': 0.05,
                'C-O_DOUBLE': 0.05,
                'AROMATIC': 0.35  # 芳香键比例
            }
            
            # 设置权重和容差
            weights = {
                'C-C_SINGLE': 1.0,
                'C-C_DOUBLE': 1.2,
                'C-C_TRIPLE': 1.3,
                'C-O_SINGLE': 1.1,
                'C-O_DOUBLE': 1.1,
                'AROMATIC': 3.0  # 给芳香键更高的权重
            }
            
            tolerances = {
                'C-C_SINGLE': 0.05,
                'C-C_DOUBLE': 0.03,
                'C-C_TRIPLE': 0.02,
                'C-O_SINGLE': 0.03,
                'C-O_DOUBLE': 0.03,
                'AROMATIC': 0.10  # 给芳香键较小的容差
            }
            
            # 创建NMR约束实例
            nmr_constraint = NMRDistanceConstraint(
                targetBondRatios=target_ratios,
                weights=weights,
                tolerances=tolerances,
                varianceSquared=0.01
            )
            
            # 2. 创建PDF约束实例
            logger.info("创建PDF约束实例...")
            
            # 检查gr文件路径
            if not os.path.isfile(self.gr_file):
                logger.error(f"PDF实验数据文件不存在: {self.gr_file}")
                return False
                
            # 创建PDF约束
            pdf_constraint = PairDistributionConstraint(
                experimentalData=self.gr_file,
                weighting="atomicNumber"  # 使用原子序数加权
            )
            
            # 初始阶段设置PDF约束范围为短程(0-8Å)
            pdf_constraint.set_limits((None, 8.0))
            
            # 3. 创建元素间距离约束 (可选)
            logger.info("创建元素间距离约束...")
            dist_constraint = InterMolecularDistanceConstraint(
                defaultDistance=0.9,  # 默认最小原子间距
                flexible=True        # 允许灵活调整
            )
            
            # 存储约束实例
            self.constraints['nmr'] = nmr_constraint
            self.constraints['pdf'] = pdf_constraint
            self.constraints['dist'] = dist_constraint
            
            # 调整不同约束之间的相对权重
            nmr_constraint.set_variance_squared(0.1)    # NMR约束初始权重略低
            pdf_constraint.set_variance_squared(1.0)    # PDF约束作为基准权重
            dist_constraint.set_variance_squared(1.2)   # 距离约束权重略高，确保物理合理性
                
            # 添加约束到引擎
            logger.info("添加约束到引擎...")
            self.engine.add_constraints(list(self.constraints.values()))
            
            # 初始阶段禁用NMR约束，先用PDF约束优化短程结构
            nmr_constraint.set_used(False)
            pdf_constraint.set_used(True)
            dist_constraint.set_used(True)
                
            # 验证约束添加成功并计算初始数据
            if len(self.engine.constraints) > 0:
                for constraint in self.engine.constraints:
                    constraint.compute_data(update=True)
                    
                logger.info(f"已添加并初始化 {len(self.engine.constraints)} 个约束")
                return True
            else:
                logger.error("没有约束被成功添加")
                return False
    
        except Exception as e:
            logger.error(f"约束设置失败: {e}")
            logger.error(traceback.format_exc())
            return False
    
    def run_simulation(self, nsteps=500000, saveFrequency=10000, xyzFrequency=5000):
        """运行多阶段优化模拟"""
        try:
            logger.info("开始多阶段优化模拟...")
            
            # 验证引擎状态
            if self.engine is None:
                logger.error("引擎未初始化，无法运行模拟")
                return False
                
            # 验证约束是否已解包
            if not self.constraints:
                self._unpack_constraints()
                
            # 执行多阶段优化
            success = self.run_multi_stage_optimization()
            
            if not success:
                logger.error("多阶段优化失败")
                return False
                
            logger.info("多阶段优化完成")
            
            # 保存最终结果
            final_path = os.path.join(self.dir_path, "final_result.rmc")
            self.engine.save(path=final_path)
            logger.info(f"最终结果保存在: {final_path}")
            
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

    def run_multi_stage_optimization(self):
        """实现从短程到长程的递进式多阶段优化策略"""
        try:
            logger.info("开始执行多阶段优化策略...")
            
            # 确保约束已解包
            if not self.constraints:
                self._unpack_constraints()
                
            if 'nmr' not in self.constraints or 'pdf' not in self.constraints:
                logger.error("缺少必要的约束(NMR/PDF)，无法执行多阶段优化")
                return False
                
            # 开始多阶段优化流程
            # 阶段1: 短程PDF优化 (初始结构调整)
            self._run_stage_1_short_range_pdf()
            
            # 阶段2: 中程PDF优化 (引入NMR约束)
            self._run_stage_2_mid_range_mixed()
            
            # 阶段3: 全程结构优化 (平衡所有约束)
            self._run_stage_3_full_range_balanced()
            
            # 阶段4: 元素特异性优化
            self._run_stage_4_element_specific()
            
            # 阶段5: 芳香结构优化
            self._run_stage_5_aromatic_optimization()
            
            # 阶段6: 最终精修
            self._run_stage_6_final_refinement()
            
            # 保存最终结果
            self.engine.save(path=os.path.join(self.dir_path, "optimized_model.rmc"))
            logger.info("多阶段优化完成")
            
            # 导出最终PDB结构
            final_pdb = os.path.join(self.dir_path, "optimized_model.pdb")
            self.engine.export_pdb(final_pdb)
            logger.info(f"最终模型已导出到: {final_pdb}")
            
            return True
            
        except Exception as e:
            logger.error(f"多阶段优化失败: {e}")
            logger.error(traceback.format_exc())
            return False
            
    def _run_stage_1_short_range_pdf(self):
        """阶段1: 短程PDF优化 - 调整基本结构而不考虑键比例"""
        logger.info("阶段1: 短程PDF优化 (调整基本结构)")
        
        # 设置约束使用状态
        self.constraints['pdf'].set_used(True)
        self.constraints['nmr'].set_used(False)
        self.constraints['dist'].set_used(True)
        
        # 设置PDF约束范围为短程(0-5Å)
        self.constraints['pdf'].set_limits((None, 5.0))
        
        # 设置简单的原子级运动
        self.engine.set_groups_as_atoms()
        
        # 为所有原子设置简单平移运动生成器
        for group in self.engine.groups:
            group.set_move_generator(TranslationGenerator(amplitude=0.2))
        
        # 设置随机选择器
        self.engine.set_group_selector(RandomSelector(self.engine))
        
        # 运行优化 - 短程结构调整
        logger.info("运行短程PDF优化...")
        self.engine.run(numberOfSteps=200000, saveFrequency=5000)
        
        # 保存阶段1结果
        self.engine.save(path=os.path.join(self.dir_path, "stage1_short_range.rmc"))
        logger.info("阶段1: 短程PDF优化完成")
        
        # 分析当前PDF拟合质量
        self._analyze_pdf_fit(stage="stage1")
    
    def _run_stage_2_mid_range_mixed(self):
        """阶段2: 中程PDF优化 (引入NMR约束)"""
        logger.info("阶段2: 中程PDF优化 (引入NMR约束)")
    
        # 获取NMR约束
        if 'nmr' in self.constraints:
            # 获取原始目标键比例
            original_targets = dict(self.constraints['nmr'].targetBondRatios)  # 创建副本而非直接引用
        
            # 计算当前键比例
            self.constraints['nmr'].compute_data(update=True)
            current_data = self.constraints['nmr'].get_constraint_value()
        
            # 创建放宽的目标 (当前值和目标值的中间点)
            relaxed_targets = {}
            for bond_type, target in original_targets.items():
                current = current_data.get(bond_type, 0.0)
                relaxed = current + 0.7 * (target - current)  # 移动70%的距离
                relaxed_targets[bond_type] = relaxed
        
            # 正确设置目标键比例 - 使用合适的方法
            try:
                # 方法1：尝试使用set_target_bond_ratios方法（如果存在）
                if hasattr(self.constraints['nmr'], 'set_target_bond_ratios'):
                    self.constraints['nmr'].set_target_bond_ratios(relaxed_targets)
                # 方法2：尝试直接更新实例字典
                elif hasattr(self.constraints['nmr'], '__dict__'):
                    self.constraints['nmr'].__dict__['targetBondRatios'] = relaxed_targets
                # 方法3：使用constraint_value属性
                elif hasattr(self.constraints['nmr'], 'set_constraint_value'):
                    self.constraints['nmr'].set_constraint_value(relaxed_targets)
                else:
                    logger.warning("无法设置NMR键比例目标，使用原始目标继续")
            except Exception as e:
                logger.warning(f"设置NMR目标时出错: {e}")
                logger.warning("使用原始目标继续")
    
        # 设置中程PDF限制 (增加到10Å)
        if 'pdf' in self.constraints:
            self.constraints['pdf'].set_limits((None, 10.0))
    
        # 调整约束权重 - PDF主导，NMR辅助
        if 'pdf' in self.constraints:
            self.constraints['pdf'].set_variance_squared(0.8)    # PDF约束强调
        if 'nmr' in self.constraints:
            self.constraints['nmr'].set_variance_squared(0.01)    # NMR约束权重降低
            # 启用NMR约束
            self.constraints['nmr'].set_used(True)
    
        # 设置递归选择器
        gs = RecursiveGroupSelector(RandomSelector(self.engine), recur=15)
        self.engine.set_group_selector(gs)
    
        # 运行优化
        self.engine.run(numberOfSteps=400000, saveFrequency=5000)
    
        # 保存阶段2结果
        self.engine.save(path=os.path.join(self.dir_path, "stage2_mid_range.rmc"))
        logger.info("阶段2: 中程PDF优化完成")
    
        # 分析PDF拟合
        self._analyze_pdf_fit(stage="stage2")
        
    def _run_stage_3_full_range_balanced(self):
        """阶段3: 全程结构优化 - 平衡所有约束"""
        logger.info("阶段3: 全程结构优化 (平衡所有约束)")
        
        # 设置PDF约束范围为全程
        self.constraints['pdf'].set_limits((None, None))
        
        # 平衡约束权重
        self.constraints['pdf'].set_variance_squared(1.0)
        self.constraints['nmr'].set_variance_squared(0.01)
        self.constraints['dist'].set_variance_squared(1.0)
        
        # 设置混合的运动生成器 - 包括平移和旋转
        for group in self.engine.groups:
            group.set_move_generator(MoveGeneratorCollector(
                collection=[
                    TranslationGenerator(amplitude=0.25),
                    RotationGenerator(amplitude=5.0)
                ],
                randomize=True
            ))
        
        # 运行优化 - 全程结构平衡
        logger.info("运行全程优化...")
        self.engine.run(numberOfSteps=1000000, saveFrequency=10000)
        
        # 保存阶段3结果
        self.engine.save(path=os.path.join(self.dir_path, "stage3_full_range.rmc"))
        logger.info("阶段3: 全程结构优化完成")
        
        # 分析当前拟合状态
        self._analyze_pdf_fit(stage="stage3")
        self._analyze_nmr_constraints(output_file=os.path.join(self.dir_path, "nmr_analysis_stage3.txt"))
    
    def _run_stage_4_element_specific(self):
        """阶段4: 元素特异性优化 - 针对不同元素进行特定优化"""
        logger.info("阶段4: 元素特异性优化")
        
        # 元素列表及其相应的振幅和递归参数
        elements = [
            {'symbol': 'C', 'amp': 0.25, 'recur': 10, 'steps': 25000},
            {'symbol': 'H', 'amp': 0.20, 'recur': 8, 'steps': 15000},
            {'symbol': 'O', 'amp': 0.20, 'recur': 8, 'steps': 15000},
        ]
        
        # 针对每种元素进行优化
        for elem in elements:
            symbol = elem['symbol']
            logger.info(f"优化元素: {symbol}")
            
            # 查找指定元素的原子
            all_elements = self.engine.allElements
            groups = []
            for idx, el in enumerate(all_elements):
                if el.lower() == symbol.lower():
                    groups.append([idx])
            
            if not groups:
                logger.warning(f"没有找到{symbol}元素的原子，跳过")
                continue
                
            logger.info(f"找到 {len(groups)} 个 {symbol} 原子")
            
            # 设置组
            self.engine.set_groups(groups)
            
            # 设置运动生成器
            for group in self.engine.groups:
                group.set_move_generator(TranslationGenerator(amplitude=elem['amp']))
            
            # 设置递归选择器
            gs = RecursiveGroupSelector(RandomSelector(self.engine), recur=elem['recur'], explore=True)
            self.engine.set_group_selector(gs)
            
            # 运行优化
            logger.info(f"运行{symbol}元素优化... ({elem['steps']}步)")
            self.engine.run(numberOfSteps=elem['steps'], saveFrequency=5000)
            
            # 保存每个元素优化后的结果
            self.engine.save(path=os.path.join(self.dir_path, f"stage4_{symbol}_optimized.rmc"))
            
            # 分析拟合状态
            self._analyze_pdf_fit(stage=f"stage4_{symbol}")
            self._analyze_nmr_constraints(output_file=os.path.join(self.dir_path, f"nmr_analysis_stage4_{symbol}.txt"))
        
        logger.info("阶段4: 元素特异性优化完成")
    
    def _run_stage_5_aromatic_optimization(self):
        """阶段5: 芳香结构优化 - 专注于优化芳香环结构"""
        logger.info("阶段5: 芳香结构优化")
        
        # 1. 识别可能的芳香环结构
        aromatic_clusters = self._identify_aromatic_clusters()
        
        if not aromatic_clusters or len(aromatic_clusters) == 0:
            logger.warning("未能识别出芳香环结构，跳过此阶段")
            return
            
        logger.info(f"识别出 {len(aromatic_clusters)} 个可能的芳香环结构")
        
        # 2. 临时增强芳香键约束
        original_targets = self.constraints['nmr'].targetBondRatios.copy()
        aromatic_enhanced = original_targets.copy()
        
        if 'AROMATIC' in aromatic_enhanced:
            aromatic_enhanced['AROMATIC'] *= 1.3  # 增强芳香键权重
            self.constraints['nmr'].targetBondRatios = aromatic_enhanced
            logger.info("临时增强芳香键权重")
        
        # 3. 设置识别出的芳香环为组，并应用特殊运动
        groups = []
        
        # 添加芳香环作为组
        for cluster in aromatic_clusters:
            if len(cluster) >= 3:  # 至少3个原子的才算一个有效组
                groups.append(cluster)
        
        # 添加剩余原子为单原子组
        all_in_clusters = set()
        for cluster in groups:
            all_in_clusters.update(cluster)
            
        for idx in range(len(self.engine.allElements)):
            if idx not in all_in_clusters:
                groups.append([idx])
        
        # 设置组
        self.engine.set_groups(groups)
        
        # 为不同大小的组设置不同的运动生成器
        for i, group in enumerate(self.engine.groups):
            if len(group.indexes) >= 3:  # 芳香环组
                # 使用组合运动生成器 - 小幅度平移和旋转
                group.set_move_generator(MoveGeneratorCollector(
                    collection=[
                        TranslationGenerator(amplitude=0.15),  # 小幅度平移
                        RotationGenerator(amplitude=3.0)       # 小角度旋转
                    ],
                    randomize=True
                ))
            else:  # 单原子组
                group.set_move_generator(TranslationGenerator(amplitude=0.25))
        
        # 运行优化
        logger.info("运行芳香结构优化...")
        self.engine.run(numberOfSteps=300000, saveFrequency=5000)
        
        # 恢复原始NMR目标
        self.constraints['nmr'].targetBondRatios = original_targets
        
        # 保存阶段5结果
        self.engine.save(path=os.path.join(self.dir_path, "stage5_aromatic_optimized.rmc"))
        logger.info("阶段5: 芳香结构优化完成")
        
        # 分析拟合状态
        self._analyze_pdf_fit(stage="stage5")
        self._analyze_nmr_constraints(output_file=os.path.join(self.dir_path, "nmr_analysis_stage5.txt"))
    
    def _run_stage_6_final_refinement(self):
        """阶段6: 最终精修 - 微调模型以达到最佳拟合"""
        logger.info("阶段6: 最终精修")
        
        # 设置更严格的容差
        if 'nmr' in self.constraints:
            strict_tolerances = {}
            for bond_type, tolerance in self.constraints['nmr'].tolerances.items():
                strict_tolerances[bond_type] = tolerance * 0.8  # 缩小容差
            
            # 应用更严格的容差
            original_tolerances = self.constraints['nmr'].tolerances.copy()
            self.constraints['nmr'].tolerances = strict_tolerances
            logger.info("应用更严格的NMR容差进行精修")
        
        # 设置所有约束为活动状态并平衡权重
        for constraint_type, constraint in self.constraints.items():
            constraint.set_used(True)
            constraint.set_variance_squared(1.0)
        
        # 设置原子级运动 - 小振幅
        self.engine.set_groups_as_atoms()
        for group in self.engine.groups:
            group.set_move_generator(TranslationGenerator(amplitude=0.1))  # 小振幅精修
        
        # 设置精修选择器
        gs = RecursiveGroupSelector(RandomSelector(self.engine), recur=10, refine=True)
        self.engine.set_group_selector(gs)
        
        # 运行精修
        logger.info("运行最终精修...")
        self.engine.run(numberOfSteps=400000, saveFrequency=5000)
        
        # 恢复原始容差
        if 'nmr' in self.constraints:
            self.constraints['nmr'].tolerances = original_tolerances
        
        # 保存阶段6结果
        self.engine.save(path=os.path.join(self.dir_path, "stage6_final_refined.rmc"))
        logger.info("阶段6: 最终精修完成")
        
        # 分析最终拟合状态
        self._analyze_pdf_fit(stage="final")
        self._analyze_nmr_constraints(output_file=os.path.join(self.dir_path, "nmr_analysis_final.txt"))
    
    def _identify_aromatic_clusters(self, max_distance=1.6):
        """识别可能的芳香环结构"""
        try:
            logger.info("识别可能的芳香族集群...")
            
            coords = self.engine.realCoordinates
            all_elements = self.engine.allElements
            carbon_indices = [idx for idx, el in enumerate(all_elements) if el.lower() == 'c']
            
            logger.info(f"找到 {len(carbon_indices)} 个碳原子")
            
            # 寻找相互距离小于max_distance的碳原子集群
            clusters = []
            used = set()
            
            for idx in carbon_indices:
                if idx in used:
                    continue
                
                cluster = [idx]
                used.add(idx)
                stack = [idx]
                
                while stack:
                    current = stack.pop()
                    for other_idx in carbon_indices:
                        if other_idx not in used:
                            dist = np.sqrt(np.sum((coords[current] - coords[other_idx])**2))
                            if dist < max_distance:
                                cluster.append(other_idx)
                                used.add(other_idx)
                                stack.append(other_idx)
                
                if len(cluster) >= 3:  # 至少3个碳原子才视为一个集群
                    clusters.append(cluster)
            
            logger.info(f"识别出 {len(clusters)} 个可能的芳香族集群")
            
            # 为每个集群添加连接的氢原子
            enriched_clusters = []
            for cluster in clusters:
                # 在集群中寻找连接的氢原子
                connected_h = []
                for c_idx in cluster:
                    for h_idx, el in enumerate(all_elements):
                        if el.lower() == 'h':
                            dist = np.sqrt(np.sum((coords[c_idx] - coords[h_idx])**2))
                            if dist < 1.2:  # C-H键典型距离
                                connected_h.append(h_idx)
                
                # 将集群作为一个组
                full_cluster = cluster + connected_h
                if len(full_cluster) > 1:
                    enriched_clusters.append(np.array(full_cluster, dtype=np.int32))
            
            return enriched_clusters
        
        except Exception as e:
            logger.error(f"识别芳香族集群失败: {e}")
            logger.error(traceback.format_exc())
            return []
    
    def _analyze_pdf_fit(self, stage="unknown"):
        """分析PDF拟合质量"""
        try:
            # 查找PDF约束
            pdf_constraint = None
            for constraint in self.engine.constraints:
                if isinstance(constraint, PairDistributionConstraint):
                    pdf_constraint = constraint
                    break
            
            if pdf_constraint is None:
                logger.warning("未找到PDF约束，无法分析拟合质量")
                return
            
            # 确保数据是最新的
            pdf_constraint.compute_data(update=True)
            
            # 记录标准误差
            error = pdf_constraint.standardError
            logger.info(f"阶段 {stage} PDF标准误差: {error:.6f}")
            
            # 绘制拟合曲线
            output_file = os.path.join(self.dir_path, f"pdf_fit_{stage}.png")
            self._plot_pdf_fit(pdf_constraint, output_file)
            
            return error
        
        except Exception as e:
            logger.error(f"分析PDF拟合质量失败: {e}")
            logger.error(traceback.format_exc())
            return None
    
    def _plot_pdf_fit(self, pdf_constraint, output_file):
        """绘制PDF拟合曲线"""
        try:
            # 使用PDF约束的内置绘图方法
            pdf_constraint.plot(
                # 设置x轴标签参数
                xlabelParams={'xlabel':'r($\AA$)', 'size':10},
                # 设置y轴标签参数
                ylabelParams={'ylabel':'g(r)', 'size':10},
                # 设置标题参数
                titleParams={'label':f'Pair Distribution Function (err: {pdf_constraint.standardError:.6f})', 'size':12},
                # 设置图例参数
                legendParams={'loc':'upper right', 'fontsize':8},
                # 调整布局参数
                subplots_adjust={'left':0.15, 'bottom':0.15},
                # 设置是否显示
                show=False
            )
        
            # 保存图像
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"PDF拟合图已保存到: {output_file}")
            
            return output_file
        
        except Exception as e:
            logger.error(f"绘制PDF拟合曲线失败: {e}")
            logger.error(traceback.format_exc())
            return None

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
        """分析NMR约束状态"""
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
                    ax.set_xlabel('Bond Type', fontsize=12)
                    ax.set_ylabel('Bond Ratio', fontsize=12)
                    
                    # 计算符合率
                    compliant_count = sum(compliant)
                    compliance_rate = compliant_count / len(bond_types) * 100
                    
                    # 设置标题，包含标准误差和符合率
                    title = f"NMR Bond Ratio Comparison (err: {constraint.standardError:.4f}, Pass rate: {compliance_rate:.1f}%)"
                    ax.set_title(title, fontsize=14)
                    
                    ax.set_xticks(x)
                    ax.set_xticklabels(bond_types, rotation=45, ha='right', fontsize=10)
                    
                    # 自定义图例
                    from matplotlib.patches import Patch
                    legend_elements = [
                        Patch(facecolor='lightcoral', label='Target ratio'),
                        Patch(facecolor='limegreen', label='Current ratio (within tolerance)'),
                        Patch(facecolor='skyblue', label='Current ratio (above tolerance)')
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

    def validate_coal_model(self):
        """验证烟煤模型的合理性"""
        logger.info("开始验证烟煤模型合理性...")
    
        # 1. 原子接触检查
        close_contacts = self._check_close_contacts()
        if close_contacts:
            logger.warning(f"发现{len(close_contacts)}对过近的原子接触")
        
        # 2. 计算芳香度
        aromaticity = self.calculate_aromaticity()
        logger.info(f"模型芳香度估计: {aromaticity:.4f}")
    
        # 3. 分析键分布
        bond_distribution = self._analyze_bond_distribution()
    
        # 4. 生成验证报告
        report_file = os.path.join(self.dir_path, "model_validation.txt")
        try:
            with open(report_file, 'w') as f:
                f.write("烟煤模型验证报告\n")
                f.write("=" * 50 + "\n\n")
            
                # 原子接触
                f.write("1. 原子接触检查\n")
                f.write("-" * 25 + "\n")
                if close_contacts:
                    f.write(f"发现{len(close_contacts)}对过近的原子接触:\n")
                    for i, (atom1, atom2, dist) in enumerate(close_contacts[:10]):  # 仅显示前10对
                        f.write(f"  {i+1}. {self.engine.allElements[atom1]}{atom1}-{self.engine.allElements[atom2]}{atom2}: {dist:.3f}Å\n")
                    if len(close_contacts) > 10:
                        f.write(f"  ...以及{len(close_contacts)-10}对其他过近接触\n")
                else:
                    f.write("未发现原子过近接触，模型通过距离检查\n")
            
                # 芳香度
                f.write("\n2. 芳香度分析\n")
                f.write("-" * 25 + "\n")
                f.write(f"模型芳香度估计: {aromaticity:.4f}\n")
                # 根据芳香度评价结果
                if aromaticity > 0.35:
                    f.write("芳香度水平较高，符合烟煤特性\n")
                elif aromaticity > 0.25:
                    f.write("芳香度水平一般，基本符合烟煤特性\n")
                else:
                    f.write("芳香度水平偏低，可能需要进一步优化芳香结构\n")
            
                # 键分布
                f.write("\n3. 键分布分析\n")
                f.write("-" * 25 + "\n")
                for bond_type, count in bond_distribution.items():
                    f.write(f"  {bond_type}: {count} 键\n")
                
                # 总体评价
                f.write("\n4. 总体评价\n")
                f.write("-" * 25 + "\n")
            
                # 计算各个指标的评分
                contact_score = 1.0 if not close_contacts else max(0, 1.0 - len(close_contacts)/100)
                aromatic_score = min(1.0, aromaticity / 0.35)
            
                # NMR约束评分
                nmr_score = 0.5  # 默认中等分数
                for constraint in self.engine.constraints:
                    if isinstance(constraint, NMRDistanceConstraint):
                        # 计算当前键比例与目标的差异
                        current_data = constraint.get_constraint_value()
                        target_ratios = constraint.targetBondRatios
                    
                        differences = []
                        for bond_type, target in target_ratios.items():
                            current = current_data.get(bond_type, 0.0)
                            diff = abs(current - target)
                            differences.append(diff)
                    
                        # 平均差异作为得分基础
                        avg_diff = sum(differences) / len(differences)
                        nmr_score = max(0, 1.0 - avg_diff * 10)  # 差异越小分数越高
            
                # PDF约束评分
                pdf_score = 0.5  # 默认中等分数
                for constraint in self.engine.constraints:
                    if isinstance(constraint, PairDistributionConstraint):
                        # 根据标准误差评分
                        error = constraint.standardError
                        pdf_score = max(0, 1.0 - error * 5)  # 误差越小分数越高
            
                # 计算总分 (按重要性加权)
                total_score = (contact_score * 0.3 + 
                              aromatic_score * 0.2 + 
                              nmr_score * 0.25 + 
                              pdf_score * 0.25) * 100
            
                # 写入评分
                f.write(f"总体得分: {total_score:.1f}/100\n")
                f.write(f"  - 原子接触: {contact_score:.2f}/1.0\n")
                f.write(f"  - 芳香度表现: {aromatic_score:.2f}/1.0\n")
                f.write(f"  - NMR约束拟合: {nmr_score:.2f}/1.0\n")
                f.write(f"  - PDF约束拟合: {pdf_score:.2f}/1.0\n\n")
            
                # 结论
                if total_score >= 80:
                    f.write("结论: 模型表现优秀，高度符合烟煤实验特性\n")
                elif total_score >= 65:
                    f.write("结论: 模型表现良好，基本符合烟煤实验特性\n")
                elif total_score >= 50:
                    f.write("结论: 模型表现一般，可以接受但有改进空间\n")
                else:
                    f.write("结论: 模型需要进一步优化\n")
        
            logger.info(f"模型验证报告已保存到: {report_file}")
            return report_file
        
        except Exception as e:
            logger.error(f"生成模型验证报告失败: {str(e)}")
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
                    current_data = constraint.get_constraint_value()  # 使用方法获取数据
                    
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
                
    def run_atoms(self, rang=5, recur=10):
        """运行原子级别优化"""
        logger.info(f"运行原子级别优化 (rang={rang}, recur={recur})...")
        
        # 设置原子组
        self.engine.set_groups_as_atoms()
        
        # 设置递归选择器
        gs = RecursiveGroupSelector(RandomSelector(self.engine), recur=recur)
        self.engine.set_group_selector(gs)
        
        # 设置运动生成器
        for group in self.engine.groups:
            group.set_move_generator(TranslationGenerator(amplitude=0.2))
        
        # 运行模拟
        self.engine.run(numberOfSteps=10000, saveFrequency=2000)
        logger.info("原子级别优化完成")
    
    def run_element_type(self, element='C', rang=8, recur=15):
        """根据元素类型运行特异性优化"""
        logger.info(f"运行{element}元素特异性优化 (rang={rang}, recur={recur})...")
        
        # 查找指定元素的原子
        element = element.lower()
        all_elements = self.engine.allElements
        groups = []
        
        for idx, el in enumerate(all_elements):
            if el.lower() == element:
                groups.append([idx])
        
        if not groups:
            logger.warning(f"未找到{element}元素的原子，跳过优化")
            return
        
        logger.info(f"找到{len(groups)}个{element}元素原子")
        
        # 设置组
        self.engine.set_groups(groups)
        
        # 设置递归选择器
        gs = RecursiveGroupSelector(RandomSelector(self.engine), recur=recur)
        self.engine.set_group_selector(gs)
        
        # 设置运动生成器 - 为不同元素使用特定振幅
        amplitude = 0.2
        if element == 'h':
            amplitude = 0.15  # 氢原子使用较小振幅
        elif element == 'o':
            amplitude = 0.18  # 氧原子使用中等振幅
            
        for group in self.engine.groups:
            group.set_move_generator(TranslationGenerator(amplitude=amplitude))
        
        # 运行模拟
        self.engine.run(numberOfSteps=15000, saveFrequency=3000)
        logger.info(f"{element}元素特异性优化完成")
    
    def run_recurring_atoms(self, rang=10, recur=15, explore=True, refine=False):
        """运行递归原子优化"""
        logger.info(f"运行递归原子优化 (rang={rang}, recur={recur}, explore={explore}, refine={refine})...")
        
        # 设置原子组
        self.engine.set_groups_as_atoms()
        
        # 设置递归选择器
        gs = RecursiveGroupSelector(RandomSelector(self.engine), 
                                   recur=recur, 
                                   explore=explore, 
                                   refine=refine)
        self.engine.set_group_selector(gs)
        
        # 设置运动生成器
        amplitude = 0.25 if explore else 0.15
        for group in self.engine.groups:
            group.set_move_generator(TranslationGenerator(amplitude=amplitude))
        
        # 运行模拟
        steps = 25000 if explore else 15000
        self.engine.run(numberOfSteps=steps, saveFrequency=5000)
        
        mode = "探索性" if explore else "精修性"
        logger.info(f"{mode}递归原子优化完成")
    
    def optimize_coal_aromatics(self):
        """优化烟煤中的芳香结构"""
        logger.info("开始优化芳香结构...")
        
        # 识别可能的芳香环结构
        aromatic_clusters = self._identify_aromatic_clusters()
        
        if not aromatic_clusters:
            logger.warning("未能识别出芳香环结构，使用替代优化方法")
            # 使用元素特异性优化替代
            self.run_element_type(element='C', rang=10, recur=20)
            return
        
        logger.info(f"识别出{len(aromatic_clusters)}个可能的芳香结构")
        
        # 临时增强芳香键约束权重
        if 'nmr' in self.constraints:
            original_weights = self.constraints['nmr'].weights.copy()
            enhanced_weights = original_weights.copy()
            
            if 'AROMATIC' in enhanced_weights:
                enhanced_weights['AROMATIC'] *= 1.5  # 增强芳香键权重
                self.constraints['nmr'].weights = enhanced_weights
                logger.info("临时增强芳香键权重")
        
        # 将芳香族集群设置为组
        groups = []
        for cluster in aromatic_clusters:
            if len(cluster) >= 3:  # 至少3个原子的集群
                groups.append(cluster)
        
        # 添加非集群中的原子作为单独组
        all_in_clusters = set()
        for cluster in groups:
            all_in_clusters.update(cluster)
            
        for idx in range(len(self.engine.allElements)):
            if idx not in all_in_clusters:
                groups.append([idx])
        
        # 设置组
        self.engine.set_groups(groups)
        
        # 设置特定运动生成器 - 芳香结构使用特殊的集体运动
        for i, group in enumerate(self.engine.groups):
            if len(group.indexes) >= 3:  # 芳香环组
                # 使用组合运动生成器
                group.set_move_generator(MoveGeneratorCollector(
                    collection=[
                        TranslationGenerator(amplitude=0.15),  # 整体平移
                        RotationGenerator(amplitude=3.0)       # 轻微旋转
                    ],
                    randomize=True
                ))
            else:  # 单原子组
                group.set_move_generator(TranslationGenerator(amplitude=0.2))
        
        # 运行优化
        self.engine.run(numberOfSteps=30000, saveFrequency=5000)
        
        # 恢复原始NMR权重
        if 'nmr' in self.constraints:
            self.constraints['nmr'].weights = original_weights
            logger.info("已恢复原始NMR约束权重")
        
        logger.info("芳香结构优化完成")
    
    def dynamic_pdf_constraint_adjustment(self, initial_range=5.0, final_range=20.0, steps=3):
        """动态调整PDF约束范围"""
        logger.info(f"开始动态调整PDF约束范围 (从{initial_range}Å到{final_range}Å，分{steps}步)")
        
        if 'pdf' not in self.constraints:
            logger.error("未找到PDF约束，无法调整")
            return False
        
        # 计算每一步的范围增量
        step_increment = (final_range - initial_range) / steps
        current_range = initial_range
        
        for step in range(steps):
            next_range = current_range + step_increment
            logger.info(f"PDF约束范围调整步骤 {step+1}/{steps}: {current_range:.1f}Å -> {next_range:.1f}Å")
            
            # 设置新范围
            self.constraints['pdf'].set_limits((None, next_range))
            
            # 运行优化
            steps_to_run = 20000 if step == 0 else 10000  # 第一步运行更多步数
            self.engine.run(numberOfSteps=steps_to_run, saveFrequency=5000)
            
            # 分析当前PDF拟合
            self._analyze_pdf_fit(stage=f"pdf_range_{next_range:.1f}")
            
            # 更新当前范围
            current_range = next_range
        
        logger.info(f"PDF约束范围已动态调整至{final_range:.1f}Å")
        return True
    
    def apply_progressive_nmr_constraint(self, steps=3):
        """逐步应用NMR约束"""
        logger.info(f"开始逐步应用NMR约束 ({steps}个阶段)")
        
        if 'nmr' not in self.constraints:
            logger.error("未找到NMR约束，无法应用")
            return False
        
        # 获取原始目标和权重
        original_targets = self.constraints['nmr'].targetBondRatios.copy()
        original_weights = self.constraints['nmr'].weights.copy()
        
        # 计算当前键比例作为起点
        self.constraints['nmr'].compute_data(update=True)
        current_data = self.constraints['nmr'].get_constraint_value()
        
        for step in range(steps):
            # 计算这一步的目标值 - 从当前值逐渐靠近目标值
            progress = (step + 1) / steps
            intermediate_targets = {}
            
            for bond_type, target in original_targets.items():
                current = current_data.get(bond_type, 0.0)
                intermediate = current + progress * (target - current)
                intermediate_targets[bond_type] = intermediate
            
            # 设置中间目标
            self.constraints['nmr'].targetBondRatios = intermediate_targets
            
            # 计算这一步的权重 - 逐步增加权重
            intermediate_weights = {}
            for bond_type, weight in original_weights.items():
                # 前期降低权重，后期接近原始权重
                factor = 0.3 + 0.7 * progress
                intermediate_weights[bond_type] = weight * factor
            
            # 设置中间权重
            self.constraints['nmr'].weights = intermediate_weights
            
            logger.info(f"NMR约束应用阶段 {step+1}/{steps} (进度: {progress*100:.0f}%)")
            
            # 运行优化
            self.engine.run(numberOfSteps=15000, saveFrequency=5000)
            
            # 分析当前NMR拟合
            self._analyze_nmr_constraints(output_file=os.path.join(
                self.dir_path, f"nmr_analysis_progressive_{step+1}.txt"))
            
            # 更新当前数据
            self.constraints['nmr'].compute_data(update=True)
            current_data = self.constraints['nmr'].get_constraint_value()
        
        # 恢复原始目标和权重
        self.constraints['nmr'].targetBondRatios = original_targets
        self.constraints['nmr'].weights = original_weights
        logger.info("已恢复原始NMR约束目标和权重")
        
        logger.info("NMR约束已逐步应用完成")
        return True
        
    def monitor_nmr_constraint_evolution(self, steps=3):
        """监控NMR约束演化"""
        logger.info(f"开始监控NMR约束演化 ({steps}个阶段)")
        
        if 'nmr' not in self.constraints:
            logger.error("未找到NMR约束，无法监控")
            return False
        
        # 初始化数据存储
        bond_types = list(self.constraints['nmr'].targetBondRatios.keys())
        evolution_data = {bond_type: [] for bond_type in bond_types}
        evolution_data['step'] = []
        evolution_data['error'] = []
        
        # 计算初始状态
        self.constraints['nmr'].compute_data(update=True)
        initial_data = self.constraints['nmr'].get_constraint_value()
        initial_error = self.constraints['nmr'].standardError
        
        # 记录初始状态
        evolution_data['step'].append(0)
        evolution_data['error'].append(initial_error)
        for bond_type in bond_types:
            evolution_data[bond_type].append(initial_data.get(bond_type, 0.0))
        
        # 运行多个阶段并记录数据
        for step in range(steps):
            logger.info(f"NMR约束监控阶段 {step+1}/{steps}")
            
            # 运行优化
            self.engine.run(numberOfSteps=10000, saveFrequency=5000)
            
            # 计算当前状态
            self.constraints['nmr'].compute_data(update=True)
            current_data = self.constraints['nmr'].get_constraint_value()
            current_error = self.constraints['nmr'].standardError
            
            # 记录当前状态
            evolution_data['step'].append(step + 1)
            evolution_data['error'].append(current_error)
            for bond_type in bond_types:
                evolution_data[bond_type].append(current_data.get(bond_type, 0.0))
        
        # 绘制演化图
        self._plot_nmr_evolution(evolution_data)
        
        logger.info("NMR约束演化监控完成")
        return evolution_data
    
    def _plot_nmr_evolution(self, evolution_data):
        """绘制NMR约束演化图"""
        try:
            import matplotlib.pyplot as plt
            
            # 创建图形
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
            
            # 绘制键比例演化
            steps = evolution_data['step']
            bond_types = [k for k in evolution_data.keys() if k not in ['step', 'error']]
            
            for bond_type in bond_types:
                ax1.plot(steps, evolution_data[bond_type], 'o-', label=bond_type)
            
            # 添加目标线
            target_ratios = self.constraints['nmr'].targetBondRatios
            for bond_type, target in target_ratios.items():
                ax1.axhline(y=target, color='gray', linestyle='--', alpha=0.5)
                # 在右侧标注目标值
                ax1.text(max(steps), target, f'{bond_type}={target:.2f}', 
                        va='center', ha='left', fontsize=8, alpha=0.7)
            
            # 设置图1属性
            ax1.set_ylabel('Bond Ratio')
            ax1.set_title('NMR Bond Ratio Evolution')
            ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
            ax1.grid(True, alpha=0.3)
            
            # 绘制标准误差演化
            ax2.plot(steps, evolution_data['error'], 'ro-', linewidth=2)
            
            # 设置图2属性
            ax2.set_xlabel('Optimization Stage')
            ax2.set_ylabel('Standard Error')
            ax2.set_title('NMR Constraint Standard Error Evolution')
            ax2.grid(True, alpha=0.3)
            
            # 添加注释 - 最终误差
            final_error = evolution_data['error'][-1]
            ax2.annotate(f'Final Error: {final_error:.6f}', 
                        xy=(steps[-1], final_error),
                        xytext=(steps[-1]-0.5, final_error*1.2),
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
            
            # 整体布局
            plt.tight_layout()
            
            # 保存图形
            output_file = os.path.join(self.dir_path, "nmr_evolution.png")
            plt.savefig(output_file, dpi=300)
            plt.close()
            
            logger.info(f"NMR约束演化图已保存到: {output_file}")
            return output_file
            
        except Exception as e:
            logger.error(f"绘制NMR约束演化图失败: {str(e)}")
            logger.error(traceback.format_exc())
            return None
    
    def _check_close_contacts(self, threshold=0.7):
        """检查过近的原子接触"""
        try:
            coords = self.engine.realCoordinates
            all_elements = self.engine.allElements
            close_pairs = []
            
            # 获取所有原子对之间的距离
            for i in range(len(all_elements)):
                for j in range(i+1, len(all_elements)):
                    dist = np.sqrt(np.sum((coords[i] - coords[j])**2))
                    min_allowed = threshold
                    
                    # 根据元素类型调整允许的最小距离
                    if all_elements[i].lower() == 'h' and all_elements[j].lower() == 'h':
                        min_allowed = 0.6  # 氢-氢可以更近
                    
                    if dist < min_allowed:
                        close_pairs.append((i, j, dist))
            
            return close_pairs
            
        except Exception as e:
            logger.error(f"检查原子接触失败: {str(e)}")
            logger.error(traceback.format_exc())
            return []
    
    def _analyze_bond_distribution(self):
        """分析模型中的键分布"""
        try:
            # 使用NMR约束的键检测功能
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    # 在这里直接调用NMR约束的键分析功能
                    return constraint.get_constraint_value()
            
            logger.warning("未找到NMR约束，无法分析键分布")
            return {}
            
        except Exception as e:
            logger.error(f"分析键分布失败: {str(e)}")
            logger.error(traceback.format_exc())
            return {}
    
    def calculate_aromaticity(self):
        """计算模型的芳香度"""
        try:
            # 查找NMR约束
            for constraint in self.engine.constraints:
                if isinstance(constraint, NMRDistanceConstraint):
                    # 计算最新数据
                    data, _ = constraint.compute_data(update=True)
                    
                    # 使用芳香键比例作为芳香度指标
                    if 'AROMATIC' in data:
                        return data['AROMATIC']
                    
                    # 如果没有直接的芳香键数据，尝试使用sp2杂化碳的比例估算
                    if 'C-C_DOUBLE' in data and 'C-C_SINGLE' in data:
                        # 假设sp2杂化的碳倾向于形成双键和芳香键
                        # 这是一个简化的估算，可能需要根据具体情况调整
                        return data.get('C-C_DOUBLE', 0.0) * 1.5
            
            logger.warning("未找到NMR约束，无法计算芳香度")
            return 0.0
            
        except Exception as e:
            logger.error(f"计算芳香度失败: {str(e)}")
            logger.error(traceback.format_exc())
            return 0.0

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
    logger.info(f"日期: {current_time.strftime('%Y-%m-%d')}")
    
    DIR_PATH = os.path.dirname(os.path.realpath(__file__))
    
    # 直接使用多阶段优化
    fresh_start = True  # 默认重新开始
    
    try:
        # 创建并初始化模拟管理器
        simulation_manager = SimulationManager(DIR_PATH, fresh_start=fresh_start)
        simulation_manager.initialize_engine()
        
        # 直接运行多阶段优化流程
        simulation_manager.run_multi_stage_optimization()
        
        # 分析模型
        simulation_manager.calculate_aromaticity()
        #simulation_manager.validate_coal_model()
        
        # 监控NMR约束演化
        simulation_manager.monitor_nmr_constraint_evolution(steps=3)
        
        # 绘制结果
        simulation_manager.plot_nmr_comparison(output_file_path=os.path.join(DIR_PATH, 'output_nmr_optimized.png'))
        
        # 导出最终模型
        final_pdb = os.path.join(DIR_PATH, "final_coal_model.pdb")
        simulation_manager.engine.export_pdb(final_pdb)
        logger.info(f"最终模型已导出到 {final_pdb}")
                
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
            logger.error("无法保存错误检查点")
            pass
