# fullrmc Structure Cleanup Preprocessing Script

## Overview

This project implements a specialized structure preprocessing pipeline using the fullrmc framework to clean non-physical structures in molecular models, particularly for coal structures before PDF fitting.

## Key Features

- **Distance Constraints**: Eliminate abnormal bonds (H-H, O-O) and optimize bond lengths (C-C, C-H)
- **Coordination Number Constraints**: Ensure reasonable atomic coordination environments  
- **Bond Angle Constraints**: Establish proper molecular geometry
- **Statistical Analysis**: Comprehensive before/after comparison and progress tracking
- **Modular Design**: Easy parameter configuration and customization

## Problem Statement

Initial coal models often contain serious non-physical structures:
- 215 abnormal H-H bonds (average 0.817 Å)
- C-C bond lengths too long (1.604 Å, should be 1.54 Å)
- C-H bond lengths too long (1.205 Å, should be 1.09 Å)  
- 18 abnormal O-O bonds
- Unreasonable coordination numbers and bond angle distributions

## Installation

1. Install required Python packages:
```bash
pip install numpy scipy matplotlib pdbparser
```

2. Clone this repository:
```bash
git clone <repository-url>
cd <repository-directory>
```

## Usage

### Basic Usage

```bash
# Clean a structure with default parameters
python StructureCleanup.py --input structure.pdb

# Specify output directory and number of cycles
python StructureCleanup.py --input coal_model.pdb --output-dir results/ --cycles 50

# Enable verbose logging
python StructureCleanup.py --input structure.pdb --verbose
```

### Command Line Options

```
--input, -i          Input PDB structure file (required)
--output-dir, -o     Output directory (default: cleanup_results)
--cycles, -c         Maximum cleanup cycles (default: 100) 
--temperature, -t    Initial temperature (default: 1.0)
--engine-name        Engine directory name (default: cleanup_engine)
--fresh-start        Start fresh simulation
--verbose, -v        Enable verbose logging
--debug              Enable debug logging
```

### Example Output

```
Structure cleanup completed successfully!

OUTPUT FILES GENERATED:
✓ cleaned_structure.pdb - Cleaned structure
✓ cleanup_statistics.txt - Statistics report  
✓ cleanup.log - Detailed log
✓ cleanup_progress.png - Progress curves
```

## File Structure

```
StructureCleanup.py    # Main preprocessing script
CleanupEngine.py       # Cleanup engine class with geometric constraints
cleanup_config.py      # Configuration parameters
```

## Configuration

### Distance Constraints

The `cleanup_config.py` file contains distance constraint definitions:

```python
DISTANCE_CONSTRAINTS = {
    ('H', 'H'): (1.5, 10.0, 0.98),   # Eliminate H-H bonds < 1.5 Å
    ('C', 'C'): (1.2, 2.0, 0.9),     # Control C-C bonds ~1.54 Å
    ('C', 'H'): (0.8, 1.4, 0.9),     # Control C-H bonds ~1.09 Å
    ('O', 'O'): (2.0, 10.0, 0.95),   # Eliminate O-O bonds < 2.0 Å
    # ... more pairs
}
```

### Coordination Constraints

```python
COORDINATION_CONSTRAINTS = [
    ('C', 'C', 1.0, 2.0, 2, 4, 1.0),   # C-C coordination (2-4 neighbors)
    ('H', 'C', 0.8, 1.4, 1, 1, 1.0),   # H-C coordination (max 1 neighbor)
    ('H', 'H', 1.5, 10.0, 0, 0, 2.0),  # Prevent H-H bonds
    # ... more definitions
]
```

### Target Values

```python
TARGET_VALUES = {
    'h_h_bonds_max': 5,              # Max allowed H-H bonds
    'c_c_distance_target': 1.54,     # Target C-C distance (Å)
    'c_h_distance_target': 1.09,     # Target C-H distance (Å)  
    'o_o_bonds_max': 2,              # Max allowed O-O bonds
}
```

## Results Interpretation

### Statistics Report

The cleanup generates a detailed statistics report showing:

```
BEFORE/AFTER COMPARISON:
h_h_bonds: 6 -> 0 (change: -6)
c_c_bonds: 3 -> 3 (change: 0)
c_h_bonds: 8 -> 6 (change: -2)
avg_c_c_distance: 1.8 -> 1.84 (change: 0.04)
avg_c_h_distance: 1.18 -> 1.17 (change: -0.01)

TARGET VALIDATION:
H-H bonds: 0 (target: ≤5) - PASS
C-C distance: 1.840 (target: 1.54±0.1) - FAIL  
C-H distance: 1.166 (target: 1.09±0.1) - PASS
O-O bonds: 1 (target: ≤2) - PASS
```

### Success Criteria

- ✅ H-H bonds reduced to near zero
- ✅ C-H bond lengths approaching 1.09 Å
- ⚠️ C-C bond lengths need further optimization
- ✅ O-O bonds controlled

## Performance Notes

- Small structures (< 100 atoms): < 5 seconds
- Medium structures (100-1000 atoms): 10-60 seconds  
- Large structures (> 1000 atoms): May take several minutes

For large structures, consider:
- Reducing the number of cycles
- Using lower temperature values
- Processing in smaller fragments

## Limitations

This is a demonstration implementation with simplified constraints. For production use:

1. **Full fullrmc Integration**: Current version uses simplified geometric corrections instead of full Monte Carlo simulation
2. **Advanced Constraints**: Bond angle constraints are basic placeholders
3. **Performance**: Not optimized for very large structures (> 10,000 atoms)

## Expected Results

After cleanup, structures should show:
- H-H bond count near zero
- C-C bond lengths closer to 1.54 Å  
- C-H bond lengths closer to 1.09 Å
- More reasonable coordination number distributions
- Improved bond angle distributions

## Troubleshooting

### Common Issues

1. **Import Error**: Ensure all dependencies are installed
```bash
pip install numpy scipy matplotlib pdbparser
```

2. **Memory Error**: For large structures, reduce cycles or split the structure

3. **Slow Performance**: Use fewer cycles or smaller move amplitudes

### Debug Mode

Enable debug logging for detailed information:
```bash
python StructureCleanup.py --input structure.pdb --debug
```

## Contributing

This preprocessing script provides a foundation for structure cleanup. Enhancements could include:

- Full fullrmc Monte Carlo integration
- Advanced angle constraint implementations  
- Parallel processing for large structures
- Additional constraint types (torsion angles, etc.)
- Web interface for easy usage

## License

[Add appropriate license information]

## Citation

If you use this script in your research, please cite:
[Add citation information]