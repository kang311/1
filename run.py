"""
BondRatioHRMCrun.py
纯结构约束驱动的HRMC模拟主运行脚本。
该模拟同时满足实验结构函数G(r)和键类型比例约束，
"""

import os 
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import logging
import traceback
import random
import time

# fullrmc导入
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.BondRatioConstraints import BondRatioConstraint
from fullrmc.Core.Group import Group
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Generators.Rotations import RotationGenerator
from fullrmc.Core.MoveGenerator import MoveGenerator
from fullrmc.Selectors.RandomSelectors import RandomSelector
from fullrmc.Core.GroupSelector import RecursiveGroupSelector
# 导入自定义引擎
from BondRatioHRMCEngine import BondRatioHRMCEngine

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

# 获取当前目录并设置路径
USER_NAME = "user"
CURRENT_DIR = os.getcwd()
DATE_TAG = datetime.now().strftime('%Y%m%d_%H%M%S')
OUTPUT_DIR = os.path.join(CURRENT_DIR, f"bondratio_{USER_NAME}_{DATE_TAG}")
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 输入文件
STRUCTURE_FILE = os.path.join(CURRENT_DIR, "coal_model.pdb")  # PDB文件的初始结构
PDF_DATA = os.path.join(CURRENT_DIR, "coal_model.gr")  # 实验G(r)数据文件

# 引擎文件
ENGINE_PATH = os.path.join(OUTPUT_DIR, "bondratio_engine.rmc")  # 引擎文件保存/加载引擎状态

# 模拟参数
HRMC_STEPS = 100000  # HRMC模拟步数
HRMC_INITIAL_TEMPERATURE = 800  # K，HRMC模拟初始温度
HRMC_FINAL_TEMPERATURE = 100  # K，HRMC模拟最终温度
COOLING_RATE = 0.98  # 温度冷却率

# 约束权重
PDF_WEIGHT = 0.6  # PDF约束权重
BOND_RATIO_WEIGHT = 0.4  # 键比例约束权重

# 可移动的元素列表
MOBILE_ELEMENTS = ['C', 'H', 'O', 'N', 'S']

# 目标键比例设置
TARGET_BOND_RATIOS = {
    'C-C_SINGLE': 0.11, 
    'C-C_DOUBLE': 0.06,
    'C-H_SINGLE': 0.66,
    'C-O_SINGLE': 0.07,
    'C-O_DOUBLE': 0.02,
    'C-N_SINGLE': 0.03,
    'O-H_SINGLE': 0.03,
    'N-H_SINGLE': 0.02
}

# 键比例权重
BOND_WEIGHTS = {
    'C-C_SINGLE': 1.0, 
    'C-C_DOUBLE': 1.2,
    'C-H_SINGLE': 1.0,
    'C-O_SINGLE': 1.1,
    'C-O_DOUBLE':1.1,
    'C-N_SINGLE': 1.1,
    'O-H_SINGLE': 1.1,
    'N-H_SINGLE': 1.1
}

# 初始化BondRatioHRMCEngine
logger.info(f"初始化BondRatioHRMCEngine，输出目录: {OUTPUT_DIR}")
ENGINE = BondRatioHRMCEngine(path=ENGINE_PATH, freshStart=True)

# 设置输出文件
ENGINE.setup_output_files(OUTPUT_DIR)

# 加载PDB结构
logger.info(f"加载结构文件: {STRUCTURE_FILE}")
ENGINE.set_pdb(STRUCTURE_FILE)

# 验证结构已加载
if len(ENGINE.allElements) == 0:
    logger.error("结构加载失败，无法找到原子元素")
    exit(1)

# 获取盒子向量
if ENGINE.isPBC:
    box_vectors = ENGINE.boundaryConditions.get_vectors()
    box_lengths = np.array([np.linalg.norm(v) for v in box_vectors])
    logger.info("系统盒子向量:")
    for i, v in enumerate(box_vectors):
        logger.info(f"向量{i+1}: {v}")
    logger.info(f"盒子尺寸: {box_lengths}")
else:
    logger.info("系统未使用周期性边界条件")

# 设置约束权重
ENGINE.set_constraint_weights(pdf_weight=PDF_WEIGHT, bond_weight=BOND_RATIO_WEIGHT)
logger.info(f"约束权重设置: PDF={PDF_WEIGHT}, 键比例={BOND_RATIO_WEIGHT}")

# 设置PDF约束
logger.info("设置PDF约束...")
PDF_CONSTRAINT = PairDistributionConstraint(
    experimentalData=PDF_DATA,
    weighting="atomicNumber",
    scaleFactor=1.0,
    adjustScaleFactor=[100, 0.1, 10.0]
)

# 设置键比例约束
logger.info("设置键比例约束...")
BOND_RATIO_CONSTRAINT = BondRatioConstraint(
    targetBondRatios=TARGET_BOND_RATIOS,
    weights=BOND_WEIGHTS
)

# 添加约束
try:
    ENGINE.add_constraints([PDF_CONSTRAINT])  # 先添加PDF约束
    ENGINE.add_constraints([BOND_RATIO_CONSTRAINT])  # 然后添加键比例约束
    logger.info("约束已成功添加到引擎")
except Exception as e:
    logger.error(f"添加约束错误: {str(e)}")
    raise

logger.info("设置原子组和移动生成器...")


# 设置组和移动生成器
box_lengths = np.array([np.linalg.norm(v) for v in box_vectors])
max_translation = np.min(box_lengths) * 0.02

print(f"盒子长度: {box_lengths}")
print(f"计算的最大平移距离: {max_translation}")

# 创建平移生成器
translation = TranslationGenerator(amplitude=max_translation)

# 对原子进行分组对于使原子集群（残基、分子等）一起演化和移动至关重要。
# 单个原子索引的组可用于使单个原子与其他原子分开移动。
# 创建组
groups = [[idx] for idx, el in enumerate( ENGINE.allElements ) if (el in MOBILE_ELEMENTS)]
ENGINE.set_groups(groups)    

# 为ENGINE中的每个组设置移动生成器
for group in ENGINE.groups:
    group.set_move_generator(translation)

# 设置组选择器
selector = RandomSelector(ENGINE)
group_selector = RecursiveGroupSelector(
    selector=selector,
    recur=10,     # 减少递归
    refine=False,  # 启用细化
    explore=True
)
ENGINE.set_group_selector(group_selector)

# 保存初始状态
logger.info("保存引擎初始状态...")
ENGINE.save()

# 通过引擎初始化约束
logger.info("初始化约束...")
try:
    ENGINE.initialize_used_constraints()
except Exception as e:  
    logger.error(f"初始化约束错误: {str(e)}")
    raise

# 设置HRMC参数
ENGINE.setup_hrmc(
    initial_temp=HRMC_INITIAL_TEMPERATURE,
    final_temp=HRMC_FINAL_TEMPERATURE,
    steps=HRMC_STEPS,
    cooling_rate=COOLING_RATE
)

logger.info(f"\n开始HRMC模拟，步数: {HRMC_STEPS}, 初始温度: {HRMC_INITIAL_TEMPERATURE}K, 最终温度: {HRMC_FINAL_TEMPERATURE}K")

# 运行模拟
start_time = time.time()
success = ENGINE.run(numberOfSteps=HRMC_STEPS, saveFrequency=5000)
elapsed_time = time.time() - start_time

if success:
    logger.info(f"模拟成功完成，总用时: {elapsed_time:.1f}秒")
    
    # 生成最终报告
    report = ENGINE.generate_summary_report()
    logger.info(f"生成模拟摘要报告: {report}")
    
    # 导出最终PDB
    final_pdb = os.path.join(OUTPUT_DIR, "final_structure.pdb")
    ENGINE.export_pdb(path=final_pdb)