#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cleanup_config.py - Configuration parameters for structure cleanup

This module contains all the configuration parameters needed for the fullrmc
structure cleanup preprocessing script.
"""

import numpy as np

# ==================== GENERAL SIMULATION PARAMETERS ====================

# Temperature parameters
INITIAL_TEMPERATURE = 1.0  # High initial temperature for large rearrangements
FINAL_TEMPERATURE = 0.1    # Final temperature after cooling
TEMPERATURE_STEPS = 10     # Number of temperature steps for cooling

# Simulation cycles
MAX_CYCLES = 100          # Maximum number of simulation cycles
STEPS_PER_CYCLE = 1000    # Steps per cycle

# ==================== DISTANCE CONSTRAINT PARAMETERS ====================

# Distance constraint parameters with high rejection probabilities
DISTANCE_CONSTRAINTS = {
    # Element pair: (min_distance, max_distance, reject_probability)
    # High rejection probability means moves that violate constraints are strongly rejected
    
    # H-H bonds - strongly eliminate abnormal H-H bonds
    ('H', 'H'): (1.5, 10.0, 0.98),   # < 1.5 Å strongly rejected
    
    # C-C bonds - control to reasonable range around 1.54 Å
    ('C', 'C'): (1.2, 2.0, 0.9),     # Target ~1.54 Å
    
    # C-H bonds - control to reasonable range around 1.09 Å  
    ('C', 'H'): (0.8, 1.4, 0.9),     # Target ~1.09 Å
    
    # O-O bonds - eliminate abnormal O-O bonds
    ('O', 'O'): (2.0, 10.0, 0.95),   # < 2.0 Å strongly rejected
    
    # Other important bonds
    ('C', 'O'): (1.0, 1.8, 0.8),     # C-O single/double bonds
    ('C', 'N'): (1.0, 1.8, 0.8),     # C-N bonds
    ('C', 'S'): (1.4, 2.2, 0.8),     # C-S bonds
    ('O', 'H'): (0.8, 1.2, 0.8),     # O-H bonds
    ('N', 'H'): (0.8, 1.2, 0.8),     # N-H bonds
    ('S', 'H'): (1.2, 1.5, 0.8),     # S-H bonds
}

# Default distance constraint parameters
DEFAULT_DISTANCE = 0.8           # Default minimum distance for unlisted pairs
DEFAULT_REJECT_PROBABILITY = 0.8  # Default rejection probability

# ==================== COORDINATION NUMBER CONSTRAINTS ====================

# Coordination number definitions: (core_element, shell_element, min_shell, max_shell, min_coord, max_coord, weight)
COORDINATION_CONSTRAINTS = [
    # Carbon coordination (2-4 neighbors typical for sp2/sp3)
    ('C', 'C', 1.0, 2.0, 2, 4, 1.0),   # C-C coordination
    ('C', 'H', 0.8, 1.4, 1, 4, 1.0),   # C-H coordination
    ('C', 'O', 1.0, 1.8, 0, 2, 1.0),   # C-O coordination
    ('C', 'N', 1.0, 1.8, 0, 2, 1.0),   # C-N coordination
    ('C', 'S', 1.4, 2.2, 0, 2, 1.0),   # C-S coordination
    
    # Hydrogen coordination (max 1 neighbor)
    ('H', 'C', 0.8, 1.4, 1, 1, 1.0),   # H-C coordination
    ('H', 'O', 0.8, 1.2, 0, 1, 1.0),   # H-O coordination
    ('H', 'N', 0.8, 1.2, 0, 1, 1.0),   # H-N coordination
    ('H', 'S', 1.2, 1.5, 0, 1, 1.0),   # H-S coordination
    ('H', 'H', 1.5, 10.0, 0, 0, 2.0),  # Prevent H-H bonds (high weight)
    
    # Oxygen coordination
    ('O', 'C', 1.0, 1.8, 1, 2, 1.0),   # O-C coordination
    ('O', 'H', 0.8, 1.2, 0, 2, 1.0),   # O-H coordination
    ('O', 'O', 2.0, 10.0, 0, 0, 2.0),  # Prevent O-O bonds (high weight)
    
    # Nitrogen coordination
    ('N', 'C', 1.0, 1.8, 1, 3, 1.0),   # N-C coordination
    ('N', 'H', 0.8, 1.2, 0, 3, 1.0),   # N-H coordination
    
    # Sulfur coordination
    ('S', 'C', 1.4, 2.2, 1, 2, 1.0),   # S-C coordination
    ('S', 'H', 1.2, 1.5, 0, 1, 1.0),   # S-H coordination
]

# Coordination constraint rejection probability
COORDINATION_REJECT_PROBABILITY = 0.9

# ==================== BOND ANGLE CONSTRAINTS ====================

# Bond angle definitions for different molecular environments
# Format: (central_atom, left_atom, right_atom, min_angle, max_angle)
ANGLE_CONSTRAINTS = {
    # sp3 carbon angles (tetrahedral ~109.5°)
    'sp3_carbon': [
        ('C', 'C', 'C', 90, 130),    # C-C-C angles
        ('C', 'C', 'H', 90, 130),    # C-C-H angles  
        ('C', 'H', 'H', 90, 130),    # H-C-H angles
        ('C', 'C', 'O', 90, 130),    # C-C-O angles
        ('C', 'C', 'N', 90, 130),    # C-C-N angles
    ],
    
    # sp2 carbon angles (planar ~120°)
    'sp2_carbon': [
        ('C', 'C', 'C', 110, 140),   # C-C-C angles in aromatic
        ('C', 'C', 'H', 110, 140),   # C-C-H angles in aromatic
    ],
    
    # Prevent very small angles (< 90°)
    'min_angles': [
        ('C', 'C', 'C', 90, 180),    # Minimum C-C-C angle
        ('C', 'C', 'H', 90, 180),    # Minimum C-C-H angle
        ('C', 'H', 'H', 90, 180),    # Minimum H-C-H angle
        ('C', 'C', 'O', 90, 180),    # Minimum C-C-O angle
        ('C', 'C', 'N', 90, 180),    # Minimum C-C-N angle
        ('C', 'C', 'S', 90, 180),    # Minimum C-C-S angle
    ]
}

# Angle constraint rejection probability
ANGLE_REJECT_PROBABILITY = 0.85

# ==================== OUTPUT CONFIGURATION ====================

# Output file names
OUTPUT_STRUCTURE_FILE = "cleaned_structure.pdb"
OUTPUT_LOG_FILE = "cleanup.log"
OUTPUT_STATS_FILE = "cleanup_statistics.txt"
OUTPUT_CURVES_FILE = "cleanup_progress.png"

# Statistics to track
STATISTICS_TO_TRACK = [
    'total_atoms',
    'h_h_bonds',
    'c_c_bonds', 
    'c_h_bonds',
    'o_o_bonds',
    'avg_c_c_distance',
    'avg_c_h_distance',
    'coordination_violations',
    'angle_violations',
    'total_energy'
]

# Reporting frequency
REPORT_FREQUENCY = 10  # Report statistics every N cycles
SAVE_FREQUENCY = 50   # Save intermediate structures every N cycles

# ==================== VALIDATION THRESHOLDS ====================

# Target values for successful cleanup
TARGET_VALUES = {
    'h_h_bonds_max': 5,           # Maximum allowed H-H bonds after cleanup
    'c_c_distance_target': 1.54,  # Target C-C bond length (Å)
    'c_c_distance_tolerance': 0.1, # Tolerance around target (±Å)
    'c_h_distance_target': 1.09,  # Target C-H bond length (Å) 
    'c_h_distance_tolerance': 0.1, # Tolerance around target (±Å)
    'o_o_bonds_max': 2,           # Maximum allowed O-O bonds after cleanup
}

# ==================== LOGGING CONFIGURATION ====================

# Logging level and format
LOG_LEVEL = "INFO"
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# Enable/disable different types of logging
ENABLE_PROGRESS_LOGGING = True
ENABLE_CONSTRAINT_LOGGING = True  
ENABLE_STATISTICS_LOGGING = True
ENABLE_DEBUG_LOGGING = False