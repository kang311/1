#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
StructureCleanup.py - Main script for fullrmc structure preprocessing

This script implements a specialized structure preprocessing pipeline using
fullrmc framework geometric constraints to clean non-physical structures
in molecular models, particularly for coal structures before PDF fitting.

Usage:
    python StructureCleanup.py --input structure.pdb [options]
    
Example:
    python StructureCleanup.py --input coal_model.pdb --cycles 50 --output cleaned_coal.pdb
"""

import os
import sys
import argparse
import logging
import time
from datetime import datetime

# Add current directory to path for local imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from CleanupEngine import CleanupEngine
    import cleanup_config as config
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure CleanupEngine.py and cleanup_config.py are in the same directory")
    sys.exit(1)

def setup_logging(log_level=None, log_file=None):
    """Setup logging configuration."""
    
    if log_level is None:
        log_level = getattr(logging, config.LOG_LEVEL.upper(), logging.INFO)
    
    if log_file is None:
        log_file = config.OUTPUT_LOG_FILE
        
    # Create formatters
    formatter = logging.Formatter(config.LOG_FORMAT)
    
    # Setup root logger
    logger = logging.getLogger()
    logger.setLevel(log_level)
    
    # Clear existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler
    try:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except Exception as e:
        print(f"Warning: Could not setup file logging: {e}")
    
    return logger

def validate_input_file(input_file):
    """Validate input structure file."""
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
        
    if not input_file.lower().endswith('.pdb'):
        raise ValueError(f"Input file must be a PDB file: {input_file}")
        
    # Check if file is readable and has content
    try:
        with open(input_file, 'r') as f:
            content = f.read().strip()
            if not content:
                raise ValueError(f"Input file is empty: {input_file}")
            if 'ATOM' not in content and 'HETATM' not in content:
                raise ValueError(f"Input file does not appear to contain atomic coordinates: {input_file}")
    except Exception as e:
        raise ValueError(f"Error reading input file: {e}")

def create_engine_directory(output_dir, engine_name="cleanup_engine"):
    """Create and return engine directory path."""
    
    engine_dir = os.path.join(output_dir, engine_name)
    
    # Create directory if it doesn't exist
    os.makedirs(engine_dir, exist_ok=True)
    
    return engine_dir

def print_header():
    """Print script header information."""
    
    print("="*80)
    print("fullrmc Structure Cleanup Preprocessing Script")
    print("="*80)
    print(f"Version: 1.0")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    print("This script uses fullrmc geometric constraints to clean non-physical")
    print("structures in molecular models, preparing them for PDF fitting.")
    print("")
    print("Key features:")
    print("- Distance constraints to eliminate abnormal bonds")
    print("- Coordination number constraints for reasonable bonding")
    print("- Bond angle constraints for proper geometry")
    print("- Statistical analysis and progress tracking")
    print("="*80)
    print("")

def print_configuration_summary(args):
    """Print summary of configuration parameters."""
    
    print("CONFIGURATION SUMMARY:")
    print("-"*30)
    print(f"Input file: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print(f"Maximum cycles: {args.cycles}")
    print(f"Initial temperature: {args.temperature}")
    print(f"Steps per cycle: {config.STEPS_PER_CYCLE}")
    print("")
    
    print("DISTANCE CONSTRAINTS:")
    print("-"*20)
    for (elem1, elem2), (min_dist, max_dist, reject_prob) in list(config.DISTANCE_CONSTRAINTS.items())[:5]:
        print(f"  {elem1}-{elem2}: min={min_dist:.1f}Å, reject_prob={reject_prob:.2f}")
    if len(config.DISTANCE_CONSTRAINTS) > 5:
        print(f"  ... and {len(config.DISTANCE_CONSTRAINTS) - 5} more pairs")
    print("")
    
    print("TARGET VALUES:")
    print("-"*15)
    for key, value in config.TARGET_VALUES.items():
        print(f"  {key}: {value}")
    print("")

def run_structure_analysis(input_file, logger):
    """Run initial structure analysis."""
    
    logger.info("Performing initial structure analysis...")
    
    try:
        # Basic file analysis
        with open(input_file, 'r') as f:
            lines = f.readlines()
            
        # Count atoms by element
        atom_counts = {}
        total_atoms = 0
        
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    # Extract element from PDB line
                    element = line[76:78].strip()
                    if not element:
                        # Fallback: extract from atom name
                        atom_name = line[12:16].strip()
                        element = atom_name[0] if atom_name else 'X'
                    
                    atom_counts[element] = atom_counts.get(element, 0) + 1
                    total_atoms += 1
                except:
                    continue
                    
        logger.info(f"Structure contains {total_atoms} atoms:")
        for element, count in sorted(atom_counts.items()):
            logger.info(f"  {element}: {count} atoms")
            
        # Basic validation
        if total_atoms == 0:
            raise ValueError("No atoms found in structure file")
            
        return atom_counts, total_atoms
        
    except Exception as e:
        logger.error(f"Error analyzing structure: {e}")
        raise

def main():
    """Main function for structure cleanup."""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Clean molecular structures using fullrmc geometric constraints',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python StructureCleanup.py --input coal_model.pdb
  python StructureCleanup.py --input structure.pdb --cycles 50 --temperature 1.5
  python StructureCleanup.py --input model.pdb --output-dir results/ --verbose
        """
    )
    
    # Required arguments
    parser.add_argument('--input', '-i', 
                       type=str, required=True,
                       help='Input PDB structure file to clean')
    
    # Optional arguments
    parser.add_argument('--output-dir', '-o',
                       type=str, default='cleanup_results',
                       help='Output directory for results (default: cleanup_results)')
    
    parser.add_argument('--cycles', '-c',
                       type=int, default=config.MAX_CYCLES,
                       help=f'Maximum number of cleanup cycles (default: {config.MAX_CYCLES})')
    
    parser.add_argument('--temperature', '-t',
                       type=float, default=config.INITIAL_TEMPERATURE,
                       help=f'Initial temperature (default: {config.INITIAL_TEMPERATURE})')
    
    parser.add_argument('--engine-name',
                       type=str, default='cleanup_engine',
                       help='Name for engine directory (default: cleanup_engine)')
    
    parser.add_argument('--fresh-start',
                       action='store_true',
                       help='Start fresh simulation (ignore existing engine data)')
    
    parser.add_argument('--verbose', '-v',
                       action='store_true',
                       help='Enable verbose logging')
    
    parser.add_argument('--debug',
                       action='store_true', 
                       help='Enable debug logging')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Determine log level
    log_level = logging.DEBUG if args.debug else (logging.INFO if args.verbose else logging.WARNING)
    
    # Print header
    print_header()
    
    try:
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Setup logging
        log_file = os.path.join(args.output_dir, config.OUTPUT_LOG_FILE)
        logger = setup_logging(log_level, log_file)
        
        logger.info("="*60)
        logger.info("STRUCTURE CLEANUP STARTED")
        logger.info("="*60)
        
        # Print configuration summary
        print_configuration_summary(args)
        
        # Validate input file
        logger.info("Validating input file...")
        validate_input_file(args.input)
        logger.info(f"Input file validation passed: {args.input}")
        
        # Analyze input structure
        atom_counts, total_atoms = run_structure_analysis(args.input, logger)
        
        # Create engine directory
        engine_dir = create_engine_directory(args.output_dir, args.engine_name)
        logger.info(f"Engine directory: {engine_dir}")
        
        # Initialize cleanup engine
        logger.info("Initializing cleanup engine...")
        engine = CleanupEngine(
            path=engine_dir,
            freshStart=args.fresh_start,
            temperature=args.temperature
        )
        
        # Run structure cleanup
        logger.info("Starting structure cleanup process...")
        start_time = time.time()
        
        engine.run_cleanup(
            structure_file=args.input,
            max_cycles=args.cycles
        )
        
        end_time = time.time()
        cleanup_time = end_time - start_time
        
        # Summary
        logger.info("="*60)
        logger.info("STRUCTURE CLEANUP COMPLETED")
        logger.info("="*60)
        logger.info(f"Total time: {cleanup_time:.1f} seconds")
        logger.info(f"Cycles completed: {args.cycles}")
        logger.info(f"Output structure: {os.path.join(args.output_dir, config.OUTPUT_STRUCTURE_FILE)}")
        logger.info(f"Statistics report: {os.path.join(args.output_dir, config.OUTPUT_STATS_FILE)}")
        logger.info(f"Full log: {log_file}")
        
        # Print file summary
        print("\nOUTPUT FILES GENERATED:")
        print("-"*30)
        output_files = [
            (config.OUTPUT_STRUCTURE_FILE, "Cleaned structure"),
            (config.OUTPUT_STATS_FILE, "Statistics report"),
            (config.OUTPUT_LOG_FILE, "Detailed log"),
            (config.OUTPUT_CURVES_FILE, "Progress curves (if matplotlib available)")
        ]
        
        for filename, description in output_files:
            full_path = os.path.join(args.output_dir, filename)
            if os.path.exists(full_path):
                size = os.path.getsize(full_path)
                print(f"✓ {filename} ({size} bytes) - {description}")
            else:
                print(f"✗ {filename} - {description} (not generated)")
        
        print(f"\nStructure cleanup completed successfully!")
        print(f"Check {args.output_dir} for all output files.")
        
        return 0
        
    except KeyboardInterrupt:
        print("\nStructure cleanup interrupted by user")
        logger.info("Structure cleanup interrupted by user")
        return 1
        
    except Exception as e:
        error_msg = f"Structure cleanup failed: {e}"
        print(f"\nERROR: {error_msg}")
        if 'logger' in locals():
            logger.error(error_msg)
            import traceback
            logger.error(traceback.format_exc())
        return 1

if __name__ == '__main__':
    sys.exit(main())