#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
example_usage.py - Example usage of the structure cleanup script

This script demonstrates different ways to use the StructureCleanup.py script
for various types of molecular structures.
"""

import os
import subprocess
import sys

def run_cleanup(input_file, output_dir, **kwargs):
    """Run structure cleanup with specified parameters."""
    
    cmd = [
        sys.executable, 'StructureCleanup.py',
        '--input', input_file,
        '--output-dir', output_dir
    ]
    
    # Add optional parameters
    if 'cycles' in kwargs:
        cmd.extend(['--cycles', str(kwargs['cycles'])])
    if 'temperature' in kwargs:
        cmd.extend(['--temperature', str(kwargs['temperature'])])
    if kwargs.get('verbose', False):
        cmd.append('--verbose')
    if kwargs.get('debug', False):
        cmd.append('--debug')
        
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            print("✓ Cleanup completed successfully")
            return True
        else:
            print(f"✗ Cleanup failed: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ Cleanup timed out")
        return False
    except Exception as e:
        print(f"✗ Error running cleanup: {e}")
        return False

def example_basic_cleanup():
    """Example 1: Basic cleanup with default parameters."""
    
    print("\n" + "="*60)
    print("EXAMPLE 1: Basic Structure Cleanup")
    print("="*60)
    
    # Use our test structure
    input_file = "test_coal_structure.pdb"
    output_dir = "example_basic_cleanup"
    
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print("Parameters: Default (100 cycles, temperature 1.0)")
    
    success = run_cleanup(input_file, output_dir, verbose=True)
    
    if success:
        print(f"\nResults saved to: {output_dir}/")
        print("Key files:")
        print(f"  - {output_dir}/cleaned_structure.pdb")
        print(f"  - {output_dir}/cleanup_statistics.txt")
        print(f"  - {output_dir}/cleanup_progress.png")

def example_quick_cleanup():
    """Example 2: Quick cleanup with fewer cycles."""
    
    print("\n" + "="*60)
    print("EXAMPLE 2: Quick Cleanup (Fewer Cycles)")
    print("="*60)
    
    input_file = "test_coal_structure.pdb"
    output_dir = "example_quick_cleanup"
    
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print("Parameters: 20 cycles for quick testing")
    
    success = run_cleanup(input_file, output_dir, cycles=20, verbose=True)
    
    if success:
        print(f"\nQuick cleanup completed in: {output_dir}/")

def example_high_temperature():
    """Example 3: High temperature for aggressive cleanup."""
    
    print("\n" + "="*60)
    print("EXAMPLE 3: High Temperature Cleanup")
    print("="*60)
    
    input_file = "test_coal_structure.pdb"
    output_dir = "example_high_temp_cleanup"
    
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print("Parameters: High temperature (2.0) for aggressive restructuring")
    
    success = run_cleanup(input_file, output_dir, 
                         cycles=50, temperature=2.0, verbose=True)
    
    if success:
        print(f"\nHigh temperature cleanup completed in: {output_dir}/")

def example_debug_mode():
    """Example 4: Debug mode for detailed analysis."""
    
    print("\n" + "="*60)
    print("EXAMPLE 4: Debug Mode Analysis")
    print("="*60)
    
    input_file = "test_coal_structure.pdb"
    output_dir = "example_debug_cleanup"
    
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print("Parameters: Debug mode with detailed logging")
    
    success = run_cleanup(input_file, output_dir, 
                         cycles=10, debug=True)
    
    if success:
        print(f"\nDebug cleanup completed in: {output_dir}/")
        print("Check the detailed log for analysis insights")

def analyze_results(output_dir):
    """Analyze cleanup results."""
    
    stats_file = os.path.join(output_dir, "cleanup_statistics.txt")
    
    if not os.path.exists(stats_file):
        print(f"Statistics file not found: {stats_file}")
        return
        
    print(f"\nAnalyzing results from: {output_dir}")
    print("-" * 40)
    
    with open(stats_file, 'r') as f:
        content = f.read()
        
    # Extract key statistics
    lines = content.split('\n')
    
    print("KEY IMPROVEMENTS:")
    for line in lines:
        if 'h_h_bonds:' in line and '->' in line:
            print(f"  H-H bonds: {line.split(':')[1].strip()}")
        elif 'avg_c_c_distance:' in line and '->' in line:
            print(f"  C-C distance: {line.split(':')[1].strip()}")
        elif 'avg_c_h_distance:' in line and '->' in line:
            print(f"  C-H distance: {line.split(':')[1].strip()}")
            
    print("\nTARGET VALIDATION:")
    for line in lines:
        if '(target:' in line and ('PASS' in line or 'FAIL' in line):
            parts = line.split()
            metric = parts[0]
            status = parts[-1]
            print(f"  {metric} {status}")

def compare_methods():
    """Compare different cleanup methods."""
    
    print("\n" + "="*60)
    print("COMPARING CLEANUP METHODS")
    print("="*60)
    
    methods = [
        ("Basic", {"cycles": 50, "temperature": 1.0}),
        ("Quick", {"cycles": 20, "temperature": 1.0}),
        ("Aggressive", {"cycles": 50, "temperature": 2.0}),
    ]
    
    input_file = "test_coal_structure.pdb"
    
    results = {}
    
    for method_name, params in methods:
        output_dir = f"comparison_{method_name.lower()}"
        print(f"\nRunning {method_name} method...")
        
        success = run_cleanup(input_file, output_dir, **params)
        
        if success:
            results[method_name] = output_dir
            analyze_results(output_dir)
        else:
            results[method_name] = None
            
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)
    
    for method_name, output_dir in results.items():
        if output_dir:
            print(f"✓ {method_name}: Results in {output_dir}/")
        else:
            print(f"✗ {method_name}: Failed")

def main():
    """Main function to run examples."""
    
    print("fullrmc Structure Cleanup - Usage Examples")
    print("="*60)
    
    # Check if test structure exists
    if not os.path.exists("test_coal_structure.pdb"):
        print("Error: test_coal_structure.pdb not found")
        print("Run the main script first to generate the test structure")
        return
        
    print("\nThis script demonstrates different usage patterns for structure cleanup.")
    print("Each example will create a separate output directory.")
    
    # Run examples
    try:
        example_basic_cleanup()
        example_quick_cleanup()
        example_high_temperature()
        example_debug_mode()
        
        # Compare methods
        compare_methods()
        
        print("\n" + "="*60)
        print("ALL EXAMPLES COMPLETED")
        print("="*60)
        print("\nGenerated directories:")
        
        directories = [
            "example_basic_cleanup",
            "example_quick_cleanup", 
            "example_high_temp_cleanup",
            "example_debug_cleanup",
            "comparison_basic",
            "comparison_quick",
            "comparison_aggressive"
        ]
        
        for directory in directories:
            if os.path.exists(directory):
                print(f"  ✓ {directory}/")
            else:
                print(f"  ✗ {directory}/ (not created)")
                
        print("\nTip: Compare the cleanup_statistics.txt files to see the differences")
        print("     between different parameter settings.")
        
    except KeyboardInterrupt:
        print("\n\nExamples interrupted by user")
    except Exception as e:
        print(f"\nError running examples: {e}")

if __name__ == '__main__':
    main()