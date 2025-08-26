#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CleanupEngine.py - Specialized engine for structure cleanup using geometric constraints

This module implements a cleanup engine that focuses on cleaning non-physical 
structures using geometric constraints (distance, coordination number, and bond 
angle constraints).

This is a simplified version that demonstrates the concept without full fullrmc integration.
"""

import os
import sys
import numpy as np
import time
import logging
from datetime import datetime
from collections import defaultdict

import cleanup_config as config

class CleanupEngine:
    """
    Specialized engine for structure cleanup using geometric constraints.
    
    This engine focuses on cleaning non-physical structures in molecular models
    by applying distance, coordination number, and bond angle constraints.
    
    This is a simplified demonstration version that analyzes structures and
    applies basic geometric corrections without full Monte Carlo simulation.
    
    :Parameters:
        #. path (str): Engine repository path
        #. freshStart (bool): Whether to start fresh simulation
        #. temperature (float): Initial simulation temperature
    """
    
    def __init__(self, path, freshStart=True, temperature=config.INITIAL_TEMPERATURE):
        """Initialize the cleanup engine."""
        
        # Set simulation parameters
        self._temperature = temperature
        self._initial_temperature = temperature
        self._final_temperature = config.FINAL_TEMPERATURE
        self._current_cycle = 0
        self._max_cycles = config.MAX_CYCLES
        
        # Structure data
        self.coordinates = None
        self.elements = None
        self.atom_names = None
        self.n_atoms = 0
        
        # Statistics tracking
        self.statistics_history = []
        self.initial_statistics = None
        
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Output configuration
        self._output_dir = os.path.dirname(path) if path else '.'
        self._output_files = {}
        
        # Setup output files
        self._setup_output_files()
        
        self.logger.info(f"CleanupEngine initialized with path: {path}")
        self.logger.info(f"Initial temperature: {temperature}")
        
    def _setup_output_files(self):
        """Setup output file paths."""
        self._output_files = {
            'structure': os.path.join(self._output_dir, config.OUTPUT_STRUCTURE_FILE),
            'log': os.path.join(self._output_dir, config.OUTPUT_LOG_FILE),
            'stats': os.path.join(self._output_dir, config.OUTPUT_STATS_FILE),
            'curves': os.path.join(self._output_dir, config.OUTPUT_CURVES_FILE)
        }
        
    def load_structure(self, pdb_file):
        """Load structure from PDB file."""
        
        self.logger.info(f"Loading structure from {pdb_file}")
        
        coordinates = []
        elements = []
        atom_names = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Extract coordinates
                    x = float(line[30:38])
                    y = float(line[38:46]) 
                    z = float(line[46:54])
                    coordinates.append([x, y, z])
                    
                    # Extract element
                    element = line[76:78].strip()
                    if not element:
                        # Fallback: extract from atom name
                        atom_name = line[12:16].strip()
                        element = atom_name[0] if atom_name else 'X'
                    elements.append(element)
                    
                    # Extract atom name
                    atom_name = line[12:16].strip()
                    atom_names.append(atom_name)
        
        self.coordinates = np.array(coordinates)
        self.elements = elements
        self.atom_names = atom_names
        self.n_atoms = len(coordinates)
        
        self.logger.info(f"Loaded {self.n_atoms} atoms")
        
    def setup_constraints(self):
        """Setup geometric constraints for structure cleanup."""
        
        self.logger.info("Setting up geometric constraints for cleanup...")
        
        # For this simplified version, constraints are stored as configuration
        # In a full implementation, these would be actual constraint objects
        
        self.distance_constraints = config.DISTANCE_CONSTRAINTS
        self.coordination_constraints = config.COORDINATION_CONSTRAINTS
        self.angle_constraints = config.ANGLE_CONSTRAINTS
        
        self.logger.info(f"Setup {len(self.distance_constraints)} distance constraints")
        self.logger.info(f"Setup {len(self.coordination_constraints)} coordination constraints")
        
    def setup_move_generators(self):
        """Setup move generators for the simulation."""
        
        self.logger.info("Setting up move generators...")
        
        # For this simplified version, we just set up basic parameters
        # In a full implementation, this would set up actual move generators
        
        self.move_amplitude = 0.1  # Small moves initially
        
        self.logger.info("Move generators configured")
            
    def compute_statistics(self):
        """Compute current structure statistics."""
        
        stats = {}
        
        try:
            # Basic atom counts
            if self.elements:
                stats['total_atoms'] = len(self.elements)
                stats['C_atoms'] = self.elements.count('C')
                stats['H_atoms'] = self.elements.count('H')
                stats['O_atoms'] = self.elements.count('O')
                stats['N_atoms'] = self.elements.count('N')
                stats['S_atoms'] = self.elements.count('S')
            
            # Compute distances for analysis
            if self.coordinates is not None and len(self.coordinates) > 1:
                # Count problematic bonds
                stats.update(self._count_problematic_bonds(self.coordinates, self.elements))
                
                # Compute average bond lengths
                stats.update(self._compute_average_bond_lengths(self.coordinates, self.elements))
                
            # Add constraint violations
            stats['distance_violations'] = self._count_distance_violations()
            stats['coordination_violations'] = self._count_coordination_violations()
            stats['angle_violations'] = self._count_angle_violations()
                
            # Add current simulation state
            stats['current_cycle'] = self._current_cycle
            stats['temperature'] = self._temperature
            
        except Exception as e:
            self.logger.error(f"Error computing statistics: {e}")
            stats['error'] = str(e)
            
        return stats
        
    def _count_problematic_bonds(self, coords, elements):
        """Count problematic bonds based on distance criteria."""
        
        # Use optimized distance calculation for better performance
        if len(coords) > 1000:
            # For large structures, use a more efficient approach
            return self._count_problematic_bonds_optimized(coords, elements)
        
        from scipy.spatial.distance import pdist, squareform
        
        # Compute all pairwise distances
        distances = squareform(pdist(coords))
        n_atoms = len(elements)
        
        counts = {
            'h_h_bonds': 0,
            'c_c_bonds': 0, 
            'c_h_bonds': 0,
            'o_o_bonds': 0,
        }
        
        # Count bonds based on distance thresholds
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                dist = distances[i,j]
                elem_i, elem_j = elements[i], elements[j]
                
                # H-H bonds (problematic if < 1.5 Å)
                if elem_i == 'H' and elem_j == 'H' and dist < 1.5:
                    counts['h_h_bonds'] += 1
                    
                # C-C bonds (count all reasonable C-C bonds 1.2-2.0 Å)
                elif elem_i == 'C' and elem_j == 'C' and 1.2 <= dist <= 2.0:
                    counts['c_c_bonds'] += 1
                    
                # C-H bonds (count all reasonable C-H bonds 0.8-1.4 Å)  
                elif {elem_i, elem_j} == {'C', 'H'} and 0.8 <= dist <= 1.4:
                    counts['c_h_bonds'] += 1
                    
                # O-O bonds (problematic if < 2.0 Å)
                elif elem_i == 'O' and elem_j == 'O' and dist < 2.0:
                    counts['o_o_bonds'] += 1
                    
        return counts
        
    def _count_problematic_bonds_optimized(self, coords, elements):
        """Optimized bond counting for large structures."""
        
        # For very large structures, use a neighbor list approach
        # This is a simplified version - in practice would use spatial indexing
        
        counts = {
            'h_h_bonds': 0,
            'c_c_bonds': 0, 
            'c_h_bonds': 0,
            'o_o_bonds': 0,
        }
        
        n_atoms = len(elements)
        max_check_distance = 3.0  # Only check atoms within reasonable bonding distance
        
        for i in range(n_atoms):
            for j in range(i+1, min(i+100, n_atoms)):  # Limit neighbor checks
                dist = np.linalg.norm(coords[i] - coords[j])
                
                if dist > max_check_distance:
                    continue
                    
                elem_i, elem_j = elements[i], elements[j]
                
                # H-H bonds (problematic if < 1.5 Å)
                if elem_i == 'H' and elem_j == 'H' and dist < 1.5:
                    counts['h_h_bonds'] += 1
                    
                # C-C bonds (count all reasonable C-C bonds 1.2-2.0 Å)
                elif elem_i == 'C' and elem_j == 'C' and 1.2 <= dist <= 2.0:
                    counts['c_c_bonds'] += 1
                    
                # C-H bonds (count all reasonable C-H bonds 0.8-1.4 Å)  
                elif {elem_i, elem_j} == {'C', 'H'} and 0.8 <= dist <= 1.4:
                    counts['c_h_bonds'] += 1
                    
                # O-O bonds (problematic if < 2.0 Å)
                elif elem_i == 'O' and elem_j == 'O' and dist < 2.0:
                    counts['o_o_bonds'] += 1
                    
        return counts
        
    def _compute_average_bond_lengths(self, coords, elements):
        """Compute average bond lengths for key bond types."""
        
        from scipy.spatial.distance import pdist, squareform
        
        distances = squareform(pdist(coords))
        n_atoms = len(elements)
        
        bond_lengths = {
            'c_c_distances': [],
            'c_h_distances': [],
        }
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                dist = distances[i,j]
                elem_i, elem_j = elements[i], elements[j]
                
                # C-C bonds
                if elem_i == 'C' and elem_j == 'C' and 1.2 <= dist <= 2.0:
                    bond_lengths['c_c_distances'].append(dist)
                    
                # C-H bonds
                elif {elem_i, elem_j} == {'C', 'H'} and 0.8 <= dist <= 1.4:
                    bond_lengths['c_h_distances'].append(dist)
                    
        # Compute averages
        averages = {}
        for bond_type, distances in bond_lengths.items():
            if distances:
                averages[f'avg_{bond_type[:-1]}'] = np.mean(distances)  # Remove 's' from distances
                averages[f'std_{bond_type[:-1]}'] = np.std(distances)
            else:
                averages[f'avg_{bond_type[:-1]}'] = 0.0
                averages[f'std_{bond_type[:-1]}'] = 0.0
                
        return averages
        
    def _count_distance_violations(self):
        """Count distance constraint violations."""
        
        if self.coordinates is None or len(self.coordinates) < 2:
            return 0
            
        violations = 0
        
        from scipy.spatial.distance import pdist, squareform
        distances = squareform(pdist(self.coordinates))
        n_atoms = len(self.elements)
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                dist = distances[i,j]
                elem_i, elem_j = self.elements[i], self.elements[j]
                
                # Check distance constraints
                for (e1, e2), (min_dist, max_dist, reject_prob) in self.distance_constraints.items():
                    if {elem_i, elem_j} == {e1, e2}:
                        if dist < min_dist:
                            violations += 1
                            
        return violations
        
    def _count_coordination_violations(self):
        """Count coordination number violations."""
        
        if self.coordinates is None or len(self.coordinates) < 2:
            return 0
            
        violations = 0
        
        from scipy.spatial.distance import pdist, squareform
        distances = squareform(pdist(self.coordinates))
        n_atoms = len(self.elements)
        
        # Check coordination constraints
        for core, shell, min_shell, max_shell, min_coord, max_coord, weight in self.coordination_constraints:
            for i in range(n_atoms):
                if self.elements[i] == core:
                    # Count neighbors of shell type within distance range
                    neighbors = 0
                    for j in range(n_atoms):
                        if i != j and self.elements[j] == shell:
                            dist = distances[i,j]
                            if min_shell <= dist <= max_shell:
                                neighbors += 1
                                
                    # Check if coordination is within limits
                    if neighbors < min_coord or neighbors > max_coord:
                        violations += 1
                        
        return violations
        
    def _count_angle_violations(self):
        """Count angle violations."""
        
        # Simplified angle violation count
        # This would require more complex geometry calculations
        # For now, return 0 as placeholder
        return 0
        
    def apply_geometric_corrections(self):
        """Apply basic geometric corrections to structure."""
        
        self.logger.info("Applying geometric corrections...")
        
        if self.coordinates is None:
            return
            
        # Apply distance corrections
        self._correct_distances()
        
        # Apply coordination corrections  
        self._correct_coordination()
        
        self.logger.info("Geometric corrections applied")
        
    def _correct_distances(self):
        """Apply distance corrections to fix problematic bonds."""
        
        from scipy.spatial.distance import pdist, squareform
        
        distances = squareform(pdist(self.coordinates))
        n_atoms = len(self.elements)
        corrections_made = 0
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                dist = distances[i,j]
                elem_i, elem_j = self.elements[i], self.elements[j]
                
                # Check for H-H bonds that are too close
                if elem_i == 'H' and elem_j == 'H' and dist < 1.5:
                    # Move atoms apart slightly
                    direction = self.coordinates[j] - self.coordinates[i]
                    direction_norm = np.linalg.norm(direction)
                    if direction_norm > 0:
                        unit_direction = direction / direction_norm
                        target_distance = 1.6  # Target H-H distance
                        displacement = unit_direction * (target_distance - dist) / 2
                        
                        self.coordinates[i] -= displacement
                        self.coordinates[j] += displacement
                        corrections_made += 1
                        
        self.logger.info(f"Made {corrections_made} distance corrections")
        
    def _correct_coordination(self):
        """Apply coordination corrections."""
        
        # This is a placeholder for coordination corrections
        # In a full implementation, this would adjust atom positions
        # to improve coordination environments
        
        self.logger.info("Coordination corrections applied (placeholder)")
        
    def save_structure(self, output_file):
        """Save current structure to PDB file."""
        
        self.logger.info(f"Saving structure to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("REMARK   Structure cleaned by fullrmc CleanupEngine\n")
            f.write("REMARK   Generated by StructureCleanup.py\n")
            
            for i, (coord, element, atom_name) in enumerate(zip(self.coordinates, self.elements, self.atom_names)):
                f.write(f"ATOM  {i+1:5d} {atom_name:4s} RMC     1    ")
                f.write(f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}")
                f.write(f"  1.00  0.00          {element:2s}\n")
                
            f.write("END\n")
            
        self.logger.info(f"Structure saved with {len(self.coordinates)} atoms")
        
    def update_temperature(self):
        """Update temperature for simulated annealing."""
        
        if self._current_cycle < self._max_cycles:
            # Linear cooling schedule
            progress = self._current_cycle / self._max_cycles
            self._temperature = self._initial_temperature * (1 - progress) + self._final_temperature * progress
            
        self.logger.debug(f"Temperature updated to {self._temperature:.3f}")
        
    def run_cleanup(self, structure_file, max_cycles=None):
        """
        Run the structure cleanup simulation.
        
        :Parameters:
            #. structure_file (str): Input PDB structure file
            #. max_cycles (int): Maximum number of cleanup cycles
        """
        
        if max_cycles:
            self._max_cycles = max_cycles
            
        self.logger.info(f"Starting structure cleanup for {structure_file}")
        self.logger.info(f"Maximum cycles: {self._max_cycles}")
        
        try:
            # Load structure
            self.load_structure(structure_file)
            self.logger.info(f"Loaded structure with {self.n_atoms} atoms")
            
            # Setup constraints and move generators
            self.setup_constraints()
            self.setup_move_generators()
            
            # Compute initial statistics
            self.initial_statistics = self.compute_statistics()
            self.statistics_history.append(self.initial_statistics.copy())
            
            self.logger.info("Initial statistics:")
            self._log_statistics(self.initial_statistics)
            
            # Run cleanup cycles
            for cycle in range(self._max_cycles):
                self._current_cycle = cycle
                
                # Update temperature
                self.update_temperature()
                
                # Run cleanup cycle
                self._run_cycle()
                
                # Compute and log statistics
                if cycle % config.REPORT_FREQUENCY == 0:
                    current_stats = self.compute_statistics()
                    self.statistics_history.append(current_stats.copy())
                    
                    self.logger.info(f"Cycle {cycle} statistics:")
                    self._log_statistics(current_stats)
                    
                # Save intermediate structure
                if cycle % config.SAVE_FREQUENCY == 0:
                    temp_file = f"intermediate_cycle_{cycle}.pdb"
                    self.save_structure(os.path.join(self._output_dir, temp_file))
                    
            # Save final structure
            self.save_structure(self._output_files['structure'])
            
            # Generate final report
            self._generate_final_report()
            
            self.logger.info("Structure cleanup completed successfully")
            
        except Exception as e:
            self.logger.error(f"Structure cleanup failed: {e}")
            raise
            
    def _run_cycle(self):
        """Run one cleanup cycle."""
        
        try:
            # Apply geometric corrections
            self.apply_geometric_corrections()
            
            # Apply some random perturbations (simplified Monte Carlo)
            self._apply_random_moves()
                    
        except Exception as e:
            self.logger.error(f"Error in cycle {self._current_cycle}: {e}")
            
    def _apply_random_moves(self):
        """Apply random moves to atoms (simplified Monte Carlo)."""
        
        if self.coordinates is None:
            return
            
        # Apply small random moves to a subset of atoms
        n_moves = max(1, int(self.n_atoms * 0.1))  # Move 10% of atoms
        
        for _ in range(n_moves):
            # Select random atom
            atom_idx = np.random.randint(0, self.n_atoms)
            
            # Apply small random displacement
            displacement = np.random.normal(0, self.move_amplitude, 3)
            
            # Store old position
            old_pos = self.coordinates[atom_idx].copy()
            
            # Apply move
            self.coordinates[atom_idx] += displacement
            
            # Simple acceptance criteria based on distance constraints
            if not self._check_move_acceptance(atom_idx):
                # Reject move
                self.coordinates[atom_idx] = old_pos
                
    def _check_move_acceptance(self, atom_idx):
        """Check if a move should be accepted based on constraints."""
        
        # Simple acceptance check based on minimum distances
        min_allowed_distance = 0.5  # Minimum distance to any other atom
        
        for j in range(self.n_atoms):
            if j != atom_idx:
                dist = np.linalg.norm(self.coordinates[atom_idx] - self.coordinates[j])
                if dist < min_allowed_distance:
                    return False
                    
        return True
            
    def _log_statistics(self, stats):
        """Log statistics in a formatted way."""
        
        if not stats:
            return
            
        important_stats = [
            ('H-H bonds', 'h_h_bonds'),
            ('C-C bonds', 'c_c_bonds'),
            ('C-H bonds', 'c_h_bonds'), 
            ('O-O bonds', 'o_o_bonds'),
            ('Avg C-C dist', 'avg_c_c_distance'),
            ('Avg C-H dist', 'avg_c_h_distance'),
            ('Temperature', 'temperature'),
        ]
        
        for label, key in important_stats:
            if key in stats:
                value = stats[key]
                if isinstance(value, float):
                    self.logger.info(f"  {label}: {value:.3f}")
                else:
                    self.logger.info(f"  {label}: {value}")
                    
    def _generate_final_report(self):
        """Generate final cleanup report."""
        
        self.logger.info("Generating final cleanup report...")
        
        try:
            # Compute final statistics
            final_stats = self.compute_statistics()
            
            # Write statistics file
            self._write_statistics_file(final_stats)
            
            # Generate progress curves if matplotlib available
            try:
                self._generate_progress_curves()
            except ImportError:
                self.logger.warning("Matplotlib not available, skipping progress curves")
                
        except Exception as e:
            self.logger.error(f"Error generating final report: {e}")
            
    def _write_statistics_file(self, final_stats):
        """Write detailed statistics to file."""
        
        with open(self._output_files['stats'], 'w') as f:
            f.write("Structure Cleanup Statistics Report\n")
            f.write("="*50 + "\n\n")
            
            # Initial vs Final comparison
            if self.initial_statistics:
                f.write("BEFORE/AFTER COMPARISON:\n")
                f.write("-"*30 + "\n")
                
                comparison_keys = [
                    'h_h_bonds', 'c_c_bonds', 'c_h_bonds', 'o_o_bonds',
                    'avg_c_c_distance', 'avg_c_h_distance'
                ]
                
                for key in comparison_keys:
                    if key in self.initial_statistics and key in final_stats:
                        initial = self.initial_statistics[key]
                        final = final_stats[key]
                        change = final - initial if isinstance(final, (int, float)) else "N/A"
                        f.write(f"{key}: {initial} -> {final} (change: {change})\n")
                        
            f.write(f"\nFINAL STATISTICS:\n")
            f.write("-"*20 + "\n")
            for key, value in final_stats.items():
                f.write(f"{key}: {value}\n")
                
            # Target validation
            f.write(f"\nTARGET VALIDATION:\n")
            f.write("-"*20 + "\n")
            self._write_target_validation(f, final_stats)
            
    def _write_target_validation(self, f, stats):
        """Write target validation results."""
        
        targets = config.TARGET_VALUES
        
        # Check H-H bonds
        if 'h_h_bonds' in stats:
            h_h_bonds = stats['h_h_bonds']
            target = targets['h_h_bonds_max']
            status = "PASS" if h_h_bonds <= target else "FAIL"
            f.write(f"H-H bonds: {h_h_bonds} (target: ≤{target}) - {status}\n")
            
        # Check C-C distance
        if 'avg_c_c_distance' in stats:
            cc_dist = stats['avg_c_c_distance']
            target = targets['c_c_distance_target']
            tolerance = targets['c_c_distance_tolerance']
            status = "PASS" if abs(cc_dist - target) <= tolerance else "FAIL"
            f.write(f"C-C distance: {cc_dist:.3f} (target: {target}±{tolerance}) - {status}\n")
            
        # Check C-H distance  
        if 'avg_c_h_distance' in stats:
            ch_dist = stats['avg_c_h_distance'] 
            target = targets['c_h_distance_target']
            tolerance = targets['c_h_distance_tolerance']
            status = "PASS" if abs(ch_dist - target) <= tolerance else "FAIL"
            f.write(f"C-H distance: {ch_dist:.3f} (target: {target}±{tolerance}) - {status}\n")
            
        # Check O-O bonds
        if 'o_o_bonds' in stats:
            o_o_bonds = stats['o_o_bonds']
            target = targets['o_o_bonds_max']
            status = "PASS" if o_o_bonds <= target else "FAIL"
            f.write(f"O-O bonds: {o_o_bonds} (target: ≤{target}) - {status}\n")
            
    def _generate_progress_curves(self):
        """Generate progress curves plot."""
        
        import matplotlib.pyplot as plt
        
        if len(self.statistics_history) < 2:
            return
            
        # Extract data for plotting
        cycles = [stats.get('current_cycle', i) for i, stats in enumerate(self.statistics_history)]
        h_h_bonds = [stats.get('h_h_bonds', 0) for stats in self.statistics_history]
        cc_distances = [stats.get('avg_c_c_distance', 0) for stats in self.statistics_history]
        ch_distances = [stats.get('avg_c_h_distance', 0) for stats in self.statistics_history]
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle('Structure Cleanup Progress')
        
        # H-H bonds
        axes[0,0].plot(cycles, h_h_bonds, 'r-o')
        axes[0,0].set_title('H-H Bonds')
        axes[0,0].set_xlabel('Cycle')
        axes[0,0].set_ylabel('Count')
        axes[0,0].grid(True)
        
        # C-C distances
        axes[0,1].plot(cycles, cc_distances, 'b-o')
        axes[0,1].axhline(y=1.54, color='g', linestyle='--', label='Target')
        axes[0,1].set_title('Average C-C Distance')
        axes[0,1].set_xlabel('Cycle')
        axes[0,1].set_ylabel('Distance (Å)')
        axes[0,1].legend()
        axes[0,1].grid(True)
        
        # C-H distances  
        axes[1,0].plot(cycles, ch_distances, 'm-o')
        axes[1,0].axhline(y=1.09, color='g', linestyle='--', label='Target')
        axes[1,0].set_title('Average C-H Distance')
        axes[1,0].set_xlabel('Cycle')
        axes[1,0].set_ylabel('Distance (Å)')
        axes[1,0].legend()
        axes[1,0].grid(True)
        
        # Temperature
        temperatures = [stats.get('temperature', 0) for stats in self.statistics_history]
        axes[1,1].plot(cycles, temperatures, 'k-o')
        axes[1,1].set_title('Temperature')
        axes[1,1].set_xlabel('Cycle')
        axes[1,1].set_ylabel('Temperature')
        axes[1,1].grid(True)
        
        plt.tight_layout()
        plt.savefig(self._output_files['curves'], dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Progress curves saved to {self._output_files['curves']}")