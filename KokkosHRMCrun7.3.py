"""
KokkosHRMCrun.py
KOKKOS增强的HRMC模拟主运行脚本
集成先进的内存管理和错误恢复机制
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import logging
import traceback
import random
import time
import psutil

# Fullrmc imports
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.DistanceConstraints import InterMolecularDistanceConstraint
from fullrmc.Core.Group import Group
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Generators.Rotations import RotationGenerator
from fullrmc.Core.MoveGenerator import MoveGenerator
from fullrmc.Selectors.RandomSelectors import RandomSelector
from fullrmc.Core.GroupSelector import RecursiveGroupSelector

# KOKKOS HRMC imports
from HRMCEngineKokkos import KokkosHRMCEngine
from EnergyConstraintsKokkos import KokkosEnergyConstraint

# 配置增强日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('kokkos_hrmc_simulation.log'),
        logging.StreamHandler()
    ]
)

##########################################################################################
################################### 模拟参数配置 ########################################

# 基本路径配置
ENGINE_PATH = "./kokkos_hrmc_simulation"
OUTPUT_DIR = "./kokkos_hrmc_results"
STRUCTURE_FILE = "initial_structure.pdb"
PDF_DATA = "experimental_gr.dat"

# ReaxFF力场参数
REAXFF_PARAMS = "ffield.reax"
CONTROL_PARAMS = "control.reax"

# KOKKOS配置
USE_KOKKOS = True
KOKKOS_DEVICE = "OpenMP"  # 或 "Serial", "Cuda"
MEMORY_CLEANUP_INTERVAL = 500  # 内存清理间隔
MAX_MEMORY_USAGE_MB = 1000     # 最大内存使用限制

# HRMC模拟参数
HRMC_STEPS = 50000
HRMC_CYCLES = 5
HRMC_INITIAL_TEMPERATURE = 2000  # K
HRMC_FINAL_TEMPERATURE = 300     # K
COOLING_RATE = 0.98
OMEGA = 1000  # 能量权重因子

# 约束参数
PDF_WEIGHT = 1.0
ENERGY_WEIGHT = 1.0
DISTANCE_CONSTRAINT_WEIGHT = 10000

# 可移动元素
MOBILE_ELEMENTS = ['C', 'H', 'O', 'N', 'S']

# 错误处理配置
MAX_CONSECUTIVE_FAILURES = 100
ADAPTIVE_ERROR_HANDLING = True
EMERGENCY_RESTART_THRESHOLD = 10

##########################################################################################
################################### 性能监控类 ##########################################

class KokkosSimulationMonitor:
    """KOKKOS模拟性能监控器"""
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.start_time = time.time()
        self.system_stats = []
        self.memory_stats = []
        self.energy_stats = []
        self.error_counts = {}
        
        # 创建监控日志文件
        self.monitor_log = os.path.join(output_dir, "system_monitor.log")
        with open(self.monitor_log, 'w') as f:
            f.write("# System Performance Monitor\n")
            f.write("# Timestamp, CPU_Percent, Memory_MB, GPU_Memory_MB, Disk_IO_MB\n")
    
    def log_system_stats(self):
        """记录系统统计"""
        try:
            timestamp = time.time() - self.start_time
            
            # CPU使用率
            cpu_percent = psutil.cpu_percent(interval=1)
            
            # 内存使用
            memory = psutil.virtual_memory()
            memory_mb = memory.used / (1024 * 1024)
            
            # 磁盘IO（简化）
            disk_io = psutil.disk_io_counters()
            disk_io_mb = (disk_io.read_bytes + disk_io.write_bytes) / (1024 * 1024) if disk_io else 0
            
            # GPU内存（如果可用）
            gpu_memory_mb = 0  # 需要GPU监控库如nvidia-ml-py
            
            # 记录数据
            stats = {
                'timestamp': timestamp,
                'cpu_percent': cpu_percent,
                'memory_mb': memory_mb,
                'gpu_memory_mb': gpu_memory_mb,
                'disk_io_mb': disk_io_mb
            }
            
            self.system_stats.append(stats)
            
            # 写入日志
            with open(self.monitor_log, 'a') as f:
                f.write(f"{timestamp:.2f} {cpu_percent:.1f} {memory_mb:.1f} {gpu_memory_mb:.1f} {disk_io_mb:.1f}\n")
                
        except Exception as e:
            logging.warning(f"系统统计记录失败: {str(e)}")
    
    def generate_performance_plots(self):
        """生成性能监控图表"""
        try:
            if not self.system_stats:
                return
                
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            
            timestamps = [s['timestamp'] for s in self.system_stats]
            
            # CPU使用率
            cpu_data = [s['cpu_percent'] for s in self.system_stats]
            axes[0, 0].plot(timestamps, cpu_data, 'b-', linewidth=1)
            axes[0, 0].set_title('CPU Usage (%)')
            axes[0, 0].set_xlabel('Time (s)')
            axes[0, 0].set_ylabel('CPU Percent')
            axes[0, 0].grid(True)
            
            # 内存使用
            memory_data = [s['memory_mb'] for s in self.system_stats]
            axes[0, 1].plot(timestamps, memory_data, 'r-', linewidth=1)
            axes[0, 1].set_title('Memory Usage (MB)')
            axes[0, 1].set_xlabel('Time (s)')
            axes[0, 1].set_ylabel('Memory (MB)')
            axes[0, 1].grid(True)
            
            # 能量演化（如果有数据）
            if self.energy_stats:
                energy_times, energies = zip(*self.energy_stats)
                axes[1, 0].plot(energy_times, energies, 'g-', linewidth=1)
                axes[1, 0].set_title('Energy Evolution')
                axes[1, 0].set_xlabel('Time (s)')
                axes[1, 0].set_ylabel('Energy (kcal/mol)')
                axes[1, 0].grid(True)
            
            # 错误统计
            if self.error_counts:
                error_types, counts = zip(*self.error_counts.items())
                axes[1, 1].bar(range(len(error_types)), counts)
                axes[1, 1].set_xticks(range(len(error_types)))
                axes[1, 1].set_xticklabels(error_types, rotation=45)
                axes[1, 1].set_title('Error Statistics')
                axes[1, 1].set_ylabel('Count')
            
            plt.tight_layout()
            plot_path = os.path.join(self.output_dir, "kokkos_performance_monitor.png")
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"性能监控图表已保存: {plot_path}")
            
        except Exception as e:
            logging.error(f"性能图表生成失败: {str(e)}")

##########################################################################################
################################### 初始化设置 ##########################################

def setup_kokkos_hrmc_engine():
    """设置KOKKOS HRMC引擎"""
    try:
        logging.info("初始化KOKKOS HRMC引擎...")
        
        # 创建输出目录
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        # 初始化引擎
        ENGINE = KokkosHRMCEngine(
            path=ENGINE_PATH, 
            freshStart=True,
            useKokkos=USE_KOKKOS
        )
        
        # 设置输出文件
        ENGINE.setup_output_files(OUTPUT_DIR)
        
        # 加载PDB结构
        if not os.path.exists(STRUCTURE_FILE):
            raise FileNotFoundError(f"结构文件不存在: {STRUCTURE_FILE}")
            
        ENGINE.set_pdb(STRUCTURE_FILE)
        
        # 显示系统信息
        box_vectors = ENGINE.boundaryConditions.get_vectors()
        logging.info("系统盒子向量:")
        for i, vector in enumerate(box_vectors):
            logging.info(f"  {i}: {vector}")
            
        logging.info(f"总原子数: {ENGINE.numberOfAtoms}")
        logging.info(f"元素类型: {set(ENGINE.allElements)}")
        
        return ENGINE
        
    except Exception as e:
        logging.error(f"KOKKOS HRMC引擎初始化失败: {str(e)}")
        raise

def setup_constraints(ENGINE):
    """设置约束条件"""
    try:
        logging.info("设置KOKKOS增强约束...")
        
        constraints = []
        
        # 1. KOKKOS能量约束
        if not os.path.exists(REAXFF_PARAMS):
            raise FileNotFoundError(f"ReaxFF参数文件不存在: {REAXFF_PARAMS}")
        if not os.path.exists(CONTROL_PARAMS):
            raise FileNotFoundError(f"控制参数文件不存在: {CONTROL_PARAMS}")
            
        ENGINE.setup_kokkos_energy_constraint(
            reaxff_params=REAXFF_PARAMS,
            control_params=CONTROL_PARAMS,
            weight=ENERGY_WEIGHT
        )
        
        # 2. PDF约束（如果存在实验数据）
        if os.path.exists(PDF_DATA):
            pdf_constraint = PairDistributionConstraint(
                experimentalData=PDF_DATA,
                weighting="atomicNumber",
                atomsWeight={"H": 0.5},  # 降低氢原子权重
                rmin=1.0,
                rmax=15.0,
                dr=0.1,
                scaleFactor=1.0,
                adjustScaleFactor=True
            )
            pdf_constraint.set_used(True)
            pdf_constraint.set_weight(PDF_WEIGHT)
            constraints.append(pdf_constraint)
            ENGINE.pdf_constraint = pdf_constraint
            logging.info("PDF约束已添加")
        else:
            logging.warning(f"PDF数据文件不存在: {PDF_DATA}")
        
        # 3. 分子间距离约束
        distance_constraint = InterMolecularDistanceConstraint(
            defaultDistance=1.5,
            pairsDistanceDefinition="element",
            flexible=True,
            rejectProbability=0.95
        )
        distance_constraint.set_used(True)
        distance_constraint.set_weight(DISTANCE_CONSTRAINT_WEIGHT)
        constraints.append(distance_constraint)
        ENGINE.dis_constraint = distance_constraint
        
        # 添加约束到引擎
        if constraints:
            ENGINE.add_constraints(constraints)
            logging.info(f"已添加 {len(constraints)} 个额外约束")
        
        return ENGINE
        
    except Exception as e:
        logging.error(f"约束设置失败: {str(e)}")
        raise

def setup_groups_and_generators(ENGINE):
    """设置原子组和移动生成器"""
    try:
        logging.info("设置原子组和KOKKOS优化的移动生成器...")
        
        # 获取可移动原子的索引
        mobile_indexes = []
        for idx, element in enumerate(ENGINE.allElements):
            if element in MOBILE_ELEMENTS:
                mobile_indexes.append(idx)
        
        if not mobile_indexes:
            raise ValueError("没有找到可移动的原子")
            
        logging.info(f"可移动原子数量: {len(mobile_indexes)}")
        
        # 创建原子组
        groups = []
        group_size = min(10, len(mobile_indexes))  # 每组最多10个原子
        
        for i in range(0, len(mobile_indexes), group_size):
            group_indexes = mobile_indexes[i:i+group_size]
            
            # 创建移动生成器
            # 平移生成器（针对KOKKOS优化）
            translation_generator = TranslationGenerator(
                amplitude=0.1,  # 较小的步长，适合KOKKOS精确计算
                direction=None,
                angle=None,
                axis=None
            )
            
            # 旋转生成器（对分子）
            rotation_generator = RotationGenerator(
                amplitude=5.0,  # 度
                direction=None,
                angle=None,
                axis=None
            )
            
            # 组合移动生成器
            move_generator = MoveGenerator(
                collection=[translation_generator, rotation_generator],
                weights=[0.8, 0.2],  # 80%平移，20%旋转
                explore=True
            )
            
            # 创建组
            group = Group(
                indexes=group_indexes,
                moveGenerator=move_generator
            )
            
            groups.append(group)
        
        # 设置组选择器
        group_selector = RecursiveGroupSelector(
            RandomSelector(ENGINE),
            recur=10,
            refine=True,
            explore=True
        )
        
        # 添加到引擎
        ENGINE.set_groups(groups)
        ENGINE.set_group_selector(group_selector)
        
        logging.info(f"已创建 {len(groups)} 个原子组")
        return ENGINE
        
    except Exception as e:
        logging.error(f"组和生成器设置失败: {str(e)}")
        raise

##########################################################################################
################################### 主运行函数 ##########################################

def run_kokkos_hrmc_simulation():
    """运行完整的KOKKOS HRMC模拟"""
    
    simulation_start_time = time.time()
    monitor = None
    
    try:
        logging.info(f"开始KOKKOS HRMC模拟 - {datetime.now()}")
        logging.info(f"KOKKOS模式: {USE_KOKKOS}")
        logging.info(f"目标设备: {KOKKOS_DEVICE}")
        
        # 初始化性能监控
        monitor = KokkosSimulationMonitor(OUTPUT_DIR)
        
        # 设置引擎
        ENGINE = setup_kokkos_hrmc_engine()
        
        # 设置约束
        ENGINE = setup_constraints(ENGINE)
        
        # 设置原子组和移动生成器
        ENGINE = setup_groups_and_generators(ENGINE)
        
        # 设置HRMC参数
        ENGINE.setup_hrmc_parameters(
            initial_temp=HRMC_INITIAL_TEMPERATURE,
            final_temp=HRMC_FINAL_TEMPERATURE,
            steps=HRMC_STEPS,
            cycles=HRMC_CYCLES,
            cooling_rate=COOLING_RATE,
            omega=OMEGA
        )
        
        # 设置内存管理参数
        ENGINE.memory_cleanup_interval = MEMORY_CLEANUP_INTERVAL
        ENGINE.force_cleanup_threshold = MAX_MEMORY_USAGE_MB * 1024 * 1024
        ENGINE.max_consecutive_failures = MAX_CONSECUTIVE_FAILURES
        ENGINE.adaptive_error_handling = ADAPTIVE_ERROR_HANDLING
        
        # 计算初始能量
        logging.info("计算初始配置能量...")
        try:
            initial_energy = ENGINE.energy_constraint.compute_energy()
            logging.info(f"初始能量: {initial_energy:.6f} kcal/mol")
            monitor.energy_stats.append((0, initial_energy))
        except Exception as e:
            logging.error(f"初始能量计算失败: {str(e)}")
            initial_energy = None
        
        # 记录初始内存状态
        initial_memory_stats = ENGINE.kokkos_memory_manager.get_memory_stats()
        logging.info(f"初始内存状态: {initial_memory_stats}")
        
        # 启动系统监控
        monitor.log_system_stats()
        
        # 执行HRMC模拟
        logging.info("开始KOKKOS HRMC主循环...")
        
        restart_count = 0
        max_restarts = EMERGENCY_RESTART_THRESHOLD
        
        while restart_count <= max_restarts:
            try:
                ENGINE.run_kokkos_hrmc_simulation()
                break  # 成功完成，退出重启循环
                
            except Exception as e:
                restart_count += 1
                error_msg = str(e)
                
                logging.error(f"模拟失败 (重启 {restart_count}/{max_restarts}): {error_msg}")
                
                if "memory" in error_msg.lower() or "allocation" in error_msg.lower():
                    logging.warn("检测到内存问题，执行紧急内存清理")
                    ENGINE._perform_emergency_memory_cleanup()
                    monitor.error_counts["Memory_Error"] = monitor.error_counts.get("Memory_Error", 0) + 1
                    
                elif "bondchk failed" in error_msg or "malformed" in error_msg:
                    logging.warn("检测到ReaxFF结构问题，重新初始化")
                    ENGINE._emergency_system_reset()
                    monitor.error_counts["ReaxFF_Error"] = monitor.error_counts.get("ReaxFF_Error", 0) + 1
                    
                else:
                    logging.error(f"未知错误: {error_msg}")
                    monitor.error_counts["Unknown_Error"] = monitor.error_counts.get("Unknown_Error", 0) + 1
                
                if restart_count > max_restarts:
                    logging.error("超过最大重启次数，模拟终止")
                    raise
                else:
                    logging.info(f"将在10秒后尝试重启模拟 ({restart_count}/{max_restarts})")
                    time.sleep(10)
        
        # 记录最终统计
        simulation_duration = time.time() - simulation_start_time
        final_memory_stats = ENGINE.kokkos_memory_manager.get_memory_stats()
        
        logging.info("=== KOKKOS HRMC模拟完成 ===")
        logging.info(f"总模拟时间: {simulation_duration:.2f} 秒")
        logging.info(f"重启次数: {restart_count}")
        logging.info(f"最终内存状态: {final_memory_stats}")
        
        # 保存最终结构
        final_pdb = os.path.join(OUTPUT_DIR, "kokkos_final_structure.pdb")
        ENGINE.export_pdb(final_pdb)
        logging.info(f"最终结构已保存: {final_pdb}")
        
        # 计算最终能量
        try:
            final_energy = ENGINE.energy_constraint.compute_energy()
            logging.info(f"最终能量: {final_energy:.6f} kcal/mol")
            
            if initial_energy is not None:
                energy_change = final_energy - initial_energy
                logging.info(f"能量变化: {energy_change:.6f} kcal/mol")
                
            monitor.energy_stats.append((simulation_duration, final_energy))
            
        except Exception as e:
            logging.error(f"最终能量计算失败: {str(e)}")
        
        # 生成性能报告
        try:
            monitor.log_system_stats()
            monitor.generate_performance_plots()
        except Exception as e:
            logging.warning(f"性能报告生成失败: {str(e)}")
        
        # 清理资源
        ENGINE.cleanup()
        
        return True
        
    except Exception as e:
        logging.error(f"KOKKOS HRMC模拟失败: {str(e)}")
        logging.error(f"错误堆栈: {traceback.format_exc()}")
        return False
        
    finally:
        # 确保生成监控报告
        if monitor:
            try:
                monitor.generate_performance_plots()
            except:
                pass

##########################################################################################
################################### 实用函数 ############################################

def validate_input_files():
    """验证输入文件"""
    required_files = [STRUCTURE_FILE, REAXFF_PARAMS, CONTROL_PARAMS]
    optional_files = [PDF_DATA]
    
    missing_required = []
    missing_optional = []
    
    for file in required_files:
        if not os.path.exists(file):
            missing_required.append(file)
    
    for file in optional_files:
        if not os.path.exists(file):
            missing_optional.append(file)
    
    if missing_required:
        raise FileNotFoundError(f"缺少必需文件: {missing_required}")
    
    if missing_optional:
        logging.warning(f"缺少可选文件: {missing_optional}")
    
    logging.info("输入文件验证通过")

def setup_environment():
    """设置环境变量和KOKKOS配置"""
    try:
        # 设置OpenMP线程数
        os.environ['OMP_NUM_THREADS'] = str(psutil.cpu_count())
        
        # 设置KOKKOS设备
        if USE_KOKKOS:
            os.environ['KOKKOS_DEVICE'] = KOKKOS_DEVICE
            
        # 设置内存限制
        import resource
        memory_limit = MAX_MEMORY_USAGE_MB * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (memory_limit, memory_limit))
        
        logging.info(f"环境配置完成: OMP_THREADS={os.environ.get('OMP_NUM_THREADS')}, "
                    f"KOKKOS_DEVICE={os.environ.get('KOKKOS_DEVICE', 'N/A')}")
        
    except Exception as e:
        logging.warning(f"环境设置部分失败: {str(e)}")

def print_simulation_summary():
    """打印模拟参数摘要"""
    print("\n" + "="*80)
    print("KOKKOS HRMC 模拟参数摘要")
    print("="*80)
    print(f"KOKKOS模式:           {USE_KOKKOS}")
    print(f"KOKKOS设备:           {KOKKOS_DEVICE}")
    print(f"结构文件:             {STRUCTURE_FILE}")
    print(f"ReaxFF参数:           {REAXFF_PARAMS}")
    print(f"HRMC步数:             {HRMC_STEPS}")
    print(f"HRMC循环:             {HRMC_CYCLES}")
    print(f"初始温度:             {HRMC_INITIAL_TEMPERATURE} K")
    print(f"最终温度:             {HRMC_FINAL_TEMPERATURE} K")
    print(f"内存清理间隔:         {MEMORY_CLEANUP_INTERVAL}")
    print(f"最大内存使用:         {MAX_MEMORY_USAGE_MB} MB")
    print(f"输出目录:             {OUTPUT_DIR}")
    print("="*80 + "\n")

##########################################################################################
################################### 主程序入口 ##########################################

if __name__ == "__main__":
    try:
        # 打印模拟摘要
        print_simulation_summary()
        
        # 验证输入文件
        validate_input_files()
        
        # 设置环境
        setup_environment()
        
        # 运行模拟
        success = run_kokkos_hrmc_simulation()
        
        if success:
            print("\n🎉 KOKKOS HRMC模拟成功完成！")
            print(f"📁 结果保存在: {OUTPUT_DIR}")
        else:
            print("\n❌ KOKKOS HRMC模拟失败")
            exit(1)
            
    except KeyboardInterrupt:
        print("\n⏹️  模拟被用户中断")
        exit(0)
    except Exception as e:
        print(f"\n💥 程序执行失败: {str(e)}")
        logging.error(f"程序执行失败: {traceback.format_exc()}")
        exit(1)