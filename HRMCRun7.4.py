#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HRMCrun.py - 杂化反向蒙特卡罗模拟运行脚本

export OMP_NUM_THREADS=16          # 设置为你机器的物理核心数
export KOKKOS_NUM_THREADS=16       # 通常设置为一样即可
export OMP_PROC_BIND=spread        # 绑定策略：线程均匀分布在物理核心上
export OMP_PLACES=threads 
使用杂化反向蒙特卡罗(HRMC)来优化原子结构，
主要基于能量约束计算系统能量。
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import logging
import argparse
import time
from pdbparser.pdbparser import pdbparser
import multiprocessing

# fullrmc导入
from fullrmc.Constraints.EnergyConstraints import EnergyConstraint
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.DistanceConstraints import InterMolecularDistanceConstraint
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Selectors.RandomSelectors import RandomSelector
from fullrmc.Core.GroupSelector import RecursiveGroupSelector
from fullrmc.Core.MoveGenerator import MoveGeneratorCollector
from fullrmc.Core.Group import Group
from fullrmc.HRMCEngine import HRMCEngine

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

# 命令行参数解析
parser = argparse.ArgumentParser(description='运行HRMC模拟')
parser.add_argument('--structure', type=str, default='coal_model.pdb', help='输入结构文件')
parser.add_argument('--gr_data', type=str, default='coal_model.gr', help='实验g(r)数据文件')
parser.add_argument('--reaxff_params', type=str, default='coal-HCONSB.ff', help='ReaxFF参数文件')
parser.add_argument('--control_params', type=str, default='lmp_control', help='LAMMPS控制文件')
parser.add_argument('--output_dir', type=str, default='', help='输出目录，默认使用时间戳')
parser.add_argument('--restart', action='store_true', help='从现有引擎重启模拟')
parser.add_argument('--hrmc_steps', type=int, default=10000, help='每个周期的HRMC步数')
parser.add_argument('--hrmc_initial_temp', type=float, default=1000, help='HRMC初始温度(K)')
parser.add_argument('--hrmc_final_temp', type=float, default=300, help='HRMC最终温度(K)')
parser.add_argument('--cycles', type=int, default=100, help='HRMC循环周期数')
args = parser.parse_args()

# 设置KOKKOS环境变量
os.environ['OMP_NUM_THREADS'] = str(multiprocessing.cpu_count())
os.environ['KOKKOS_NUM_THREADS'] = str(multiprocessing.cpu_count())
# 设置绑定策略
os.environ['OMP_PROC_BIND'] = 'spread'
os.environ['OMP_PLACES'] = 'threads'
# 添加KOKKOS开关
USE_KOKKOS = os.environ.get('USE_KOKKOS', 'true').lower() == 'true'

# 设置路径
CURRENT_DIR = os.getcwd()
DATE_TAG = datetime.now().strftime('%Y%m%d_%H%M%S')
OUTPUT_DIR = args.output_dir if args.output_dir else os.path.join(CURRENT_DIR, f"hrmc_{DATE_TAG}")
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 输入输出文件
STRUCTURE_FILE = os.path.join(CURRENT_DIR, args.structure) 
PDF_DATA = os.path.join(CURRENT_DIR, args.gr_data)
REAXFF_PARAMS = os.path.join(CURRENT_DIR, args.reaxff_params)
CONTROL_PARAMS = os.path.join(CURRENT_DIR, args.control_params)
ENGINE_PATH = os.path.join(OUTPUT_DIR, "hrmc_engine.rmc")
if not os.path.exists(ENGINE_PATH):
    os.makedirs(ENGINE_PATH)

def main():
    """主模拟流程"""
    logger.info(f"启动HRMC模拟")
    logger.info(f"当前时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"输出目录: {OUTPUT_DIR}")

    # 初始化引擎
    if os.path.exists(ENGINE_PATH) and args.restart:
        # 加载现有引擎
        ENGINE = HRMCEngine(path=ENGINE_PATH)
        ENGINE = ENGINE.load(ENGINE_PATH)
        logger.info("已从保存的引擎状态加载")

        n_atoms = len(ENGINE.allElements)
        if len(np.unique(ENGINE.moleculesIndex)) != n_atoms:
            ENGINE.moleculesIndex[:] = np.arange(n_atoms, dtype=np.int32)
            logger.info(f"[restart] 已重新分配 {n_atoms} 个唯一 molecule ID")
        else:
            logger.info("[restart] molecule ID 唯一，无需修改")
    else:
        # 创建新引擎
        ENGINE = HRMCEngine(path=ENGINE_PATH, freshStart=True)
        logger.info("创建新的引擎实例")
    
    # 如果是新引擎，加载结构和设置约束
    if not args.restart or not os.path.exists(ENGINE_PATH):
        # 设置输出文件
        ENGINE.setup_output_files(OUTPUT_DIR)
        
        # 加载PDB结构
        logger.info(f"加载结构文件: {STRUCTURE_FILE}")
        ENGINE.set_pdb(STRUCTURE_FILE)

        box_vectors = ENGINE.boundaryConditions.get_vectors()
        print("System box vectors:")
        print(box_vectors)

        # ==== 让每个原子拥有唯一的 molecule ID
        n_atoms = len(ENGINE.allElements)
        ENGINE.moleculesIndex[:] = np.arange(n_atoms, dtype=np.int32)
        logger.info(f"已为 {n_atoms} 个原子分配唯一 molecule ID（0 ~ {n_atoms-1}）")

        # 创建约束，指定Q值范围
        q_range = {
            'qmin': 0.0001,  # Å^-1
            'qmax': 30.0,    # Å^-1
            'qstep': 0.01    # Å^-1
        }

        # 创建g(r)约束 - 注意此处使用g(r)而不是G(r)
        logger.info("设置PDF约束...")
        PDF_CONSTRAINT = PairDistributionConstraint(
            experimentalData=PDF_DATA, 
            weighting="xray_scattering",  # 使用X射线散射因子
            qRange=q_range,               # 指定Q值范围
            shapeFuncParams=None,         # 不使用形状函数
            windowFunction=None,
            adjustScaleFactor=(100, 0.8, 1.2)  # 每100步调整一次比例因子
        )               
        
        # 设置分子间距离约束
        logger.info("设置分子间距离约束...")
        dis_CONSTRAINT = InterMolecularDistanceConstraint(
            defaultDistance=1.0,  # 默认最小距离
            typeDefinition='element',
            pairsDistanceDefinition=[
                ('C', 'C', 1.2),   # C-C最小距离
                ('C', 'H', 1.0),   # C-H最小距离
                ('C', 'O', 1.1),   # C-O最小距离
                ('C', 'N', 1.1),   # C-N最小距离
                ('C', 'S', 1.5),   # C-S最小距离
                ('H', 'H', 0.8),   # H-H最小距离
                ('H', 'O', 0.9),   # H-O最小距离
                ('H', 'N', 0.9),   # H-N最小距离
                ('H', 'S', 1.2),   # H-S最小距离
                ('O', 'O', 1.2),   # O-O最小距离
                ('O', 'N', 1.1),   # O-N最小距离
                ('O', 'S', 1.4),   # O-S最小距离
                ('N', 'N', 1.1),   # N-N最小距离
                ('N', 'S', 1.5),   # N-S最小距离
                ('S', 'S', 1.8),   # S-S最小距离
            ],
            flexible=False  # 使约束更严格
        )

        dis_CONSTRAINT.set_type_definition("element")

        # 创建能量约束
        logger.info("设置能量约束...")
        ENERGY_CONSTRAINT = EnergyConstraint(
            reaxff_params=REAXFF_PARAMS,
            control_params=CONTROL_PARAMS,
            compute_forces=False
        )

        # 添加约束 - 首先添加PDF约束和距离约束，然后添加能量约束
        ENGINE.add_constraints([PDF_CONSTRAINT,dis_CONSTRAINT])
        ENGINE.add_constraints([ENERGY_CONSTRAINT])
        logger.info("约束已成功添加到引擎")
       
        
        # 所有要调整的元素
        elements_constraint = ['C', 'H', 'O', 'S', 'N']
        
        # 创建组，只包括指定元素的原子
        groups = [[idx] for idx, el in enumerate(ENGINE.allElements) if el in elements_constraint]
        ENGINE.set_groups(groups)
        logger.info(f"已创建{len(groups)}个原子组")
        
        # 创建平移生成器
        translation = TranslationGenerator(amplitude=0.1)
        
        max_translation = 0.1

        # 为所有组设置移动生成器为TranslationGenerator
        for group in ENGINE.groups:
            group.set_move_generator(translation)
        logger.info(f"已为所有组设置平移移动生成器，振幅: {max_translation}")
        
        # 设置组选择器 - 使用RecursiveGroupSelector提高接受率
        selector = RandomSelector(ENGINE)
        group_selector = RecursiveGroupSelector(
            selector=selector,
            recur=10,     # 调整递归深度
            refine=False, # 关闭精细化
            explore=True  # 启用探索模式
        )
        ENGINE.set_group_selector(group_selector)
        logger.info("已设置递归组选择器")
        
        # 初始化约束
        logger.info("初始化约束...")
        ENGINE.initialize_used_constraints()

        # 记录初始误差
        initial_pdf_error = ENGINE.pdf_constraint.standardError
        initial_energy_error = ENGINE.energy_constraint.data
        logger.info(f"初始PDF误差: {initial_pdf_error:.6f}")
        logger.info(f"初始能量误差: {initial_energy_error:.6f}")

        # 保存初始状态
        logger.info("保存初始引擎状态...")
        ENGINE.save()
        
    # 开始计时
    start_time = time.time()
    

    # 设置HRMC参数
    ENGINE.setup_hrmc(
        initial_temp=args.hrmc_initial_temp,
        final_temp=args.hrmc_final_temp, 
        steps=args.hrmc_steps,
        cycles=args.cycles
    )
    
    # 运行HRMC模拟
    try:
        logger.info("开始HRMC模拟...")
        ENGINE.run_hrmc()
    except Exception as e:
        logger.error(f"HRMC模拟失败: {str(e)}")
        raise
        
    # 计算总用时
    total_time = time.time() - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    # 输出结果摘要
    if hasattr(ENGINE, '_tried') and hasattr(ENGINE, '_accepted'):
        acceptance_rate = 100.0 * ENGINE._accepted / max(1, ENGINE._tried)
        logger.info(f"\n模拟完成! 总用时: {int(hours)}小时 {int(minutes)}分钟 {seconds:.1f}秒")
        logger.info(f"尝试步数: {ENGINE._tried}, 接受步数: {ENGINE._accepted}, 接受率: {acceptance_rate:.2f}%")
    
    # 获取最终能量
    energy_value = None
    if hasattr(ENGINE, 'energy_constraint') and ENGINE.energy_constraint:
        energy_value = ENGINE.energy_constraint.data
        logger.info(f"最终能量: {energy_value:.6f}")
    
    
    # 绘制最终PDF对比图
    final_pdf_plot = os.path.join(OUTPUT_DIR, "final_pdf_comparison.png")
    ENGINE.plot_pdf_comparison(save_path=final_pdf_plot)
    logger.info(f"最终PDF对比图已保存至: {final_pdf_plot}")

    # 绘制能量变化图
    energy_plot = os.path.join(OUTPUT_DIR, "energy_profile.png")
    ENGINE.plot_energy_profile(save_path=energy_plot)
    logger.info(f"能量变化图已保存至: {energy_plot}")
    
    # 保存最终结构
    final_pdb = os.path.join(OUTPUT_DIR, "final_structure.pdb")
    ENGINE.export_pdb(path=final_pdb)
    logger.info(f"最终结构已保存至: {final_pdb}")
    
    logger.info(f"所有结果已保存到: {OUTPUT_DIR}")
    
    # 清理资源
    ENGINE.cleanup()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.exception(f"模拟过程中发生错误: {str(e)}")     
