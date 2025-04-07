import os, sys
import numpy as np
import pickle
from datetime import datetime
from pdbparser.pdbparser import pdbparser

# fullrmc imports
from fullrmc.Engine import Engine
from fullrmc.Core.Group import Group
from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint
from fullrmc.Constraints.EnergyConstraints import EnergyConstraint
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Generators.Rotations import RotationGenerator
from fullrmc.Core.MoveGenerator import MoveGeneratorCollector
from fullrmc.Core.GroupSelector import RecursiveGroupSelector
from fullrmc.Selectors.RandomSelectors import RandomSelector
from MDHRMCEngine import MDHRMCEngine

##########################################################################################
#################################### USER PARAMETERS #########################################

# Get current directory and setup paths
CURRENT_DIR = os.getcwd()
USER_NAME = "zhang"
DATE_TAG = datetime.utcnow().strftime('%Y%m%d_%H%M%S')
OUTPUT_DIR = os.path.join(CURRENT_DIR, f"simulation_{USER_NAME}_{DATE_TAG}")
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Input files
STRUCTURE_FILE = os.path.join(CURRENT_DIR, "small_coal_model.pdb") # PDB file with initial structure
PDF_DATA = os.path.join(CURRENT_DIR, "Gr-qmax30.gr") # Experimental G(r) data file
REAXFF_PARAMS = os.path.join(CURRENT_DIR, "coal-HCONSB.ff") # ReaxFF force field parameters
CONTROL_PARAMS = os.path.join(CURRENT_DIR, "lmp_control")   # LAMMPS ReaxFF control file

# Engine file
ENGINE_PATH = os.path.join(OUTPUT_DIR, "engine.rmc") # Engine file to save/load engine state

# Simulation parameters
CYCLES_MD_HRMC = 3 # Number of MD-HRMC cycles
MD_STEPS = 10 # MD steps in each MD-HRMC cycle
HRMC_STEPS = 10 # HRMC steps in each MD-HRMC cycle

MD_TEMPERATURE = 2000 # K, temperature for MD simulation

HRMC_INITIAL_TEMPERATURE = 2000 # K, temperature for HRMC simulation
HRMC_FINAL_TEMPERATURE = 300 # K, temperature for HRMC simulation

# All elements whose positions are to be adjusted within HRMC simulation section in STRUCTURE_FILE.
# Atoms not included in any group are held fixed during HRMC simulation, 
# but retain mobility during MD simulation.
elements_constraint = ['C', 'H', 'O', 'N', 'S'] # ['C', 'O', 'N', 'S']  

##########################################################################################
#################################### BASIC SETUP ##########################################

# Initialize MD-HRMC engine
ENGINE = MDHRMCEngine(path=ENGINE_PATH, freshStart=True)

# Setup output files
ENGINE.setup_output_files(OUTPUT_DIR)

# Load PDB structure
ENGINE.set_pdb(STRUCTURE_FILE)

box_vectors = ENGINE.boundaryConditions.get_vectors()
print("System box vectors:")
print(box_vectors)

##########################################################################################
################################# SETUP CONSTRAINTS #####################################

# Setup PDF constraint
PDF_CONSTRAINT = PairDistributionConstraint(
    experimentalData=PDF_DATA,
    weighting="atomicNumber",
    scaleFactor=1.0,
    adjustScaleFactor=[100, 0.1, 10.0]
)

# Setup Energy constraint
ENERGY_CONSTRAINT = EnergyConstraint(
    reaxff_params=REAXFF_PARAMS,
    control_params=CONTROL_PARAMS,
    compute_forces=False
)

# Add constraints
try:
    ENGINE.add_constraints([PDF_CONSTRAINT])  # Add PDF constraint first
    ENGINE.add_constraints([ENERGY_CONSTRAINT])  # Then add Energy constraint
except Exception as e:
    print(f"Error adding constraints: {str(e)}")
    raise

# Initialize explicitly
ENERGY_CONSTRAINT.initialize() # Initialize energy constraint

##########################################################################################
################################# SETUP GROUPS AND GENERATORS ###########################

# Setup groups and move generators
box_lengths = np.array([np.linalg.norm(v) for v in box_vectors])
max_translation = np.min(box_lengths) * 0.02

print(f"Box lengths: {box_lengths}")
print(f"Calculated max translation: {max_translation}")

# Create translation generator
translation = TranslationGenerator(amplitude=max_translation)

# Grouping atoms is essential to make clusters of atoms (residues, molecules, etc) evolve and move together.
# A group of a single atom index can be used to make a single atom move separately from the others.
# Create groups
groups = [[idx] for idx, el in enumerate( ENGINE.allElements ) if (el in elements_constraint)]
ENGINE.set_groups(groups)    

# Set move generator for each group in the ENGINE
for group in ENGINE.groups:
    group.set_move_generator(translation)

# Setup group selector
selector = RandomSelector(ENGINE)
group_selector = RecursiveGroupSelector(
    selector=selector,
    recur=10,     # Reduce recursion
    refine=False,  # Enable refinement
    explore=True
)
ENGINE.set_group_selector(group_selector)

# Save initial state before LAMMPS initialization
print("Saving initial engine state...")
ENGINE.save()

# Initialize constraints through engine only
print("Initializing constraints...")
try:
    ENGINE.initialize_used_constraints()
except Exception as e:  
    print(f"Error initializing constraints: {str(e)}")
    raise

##########################################################################################
#################################### RUN SIMULATION ######################################
    
# Run simulation loop
print(f"\nStarting MD-HRMC simulation at {datetime.utcnow()}")
    
# Setup MD parameters
ENGINE.setup_md(
    temperature=MD_TEMPERATURE,
    steps=MD_STEPS,
    reaxff_params=REAXFF_PARAMS,
    control_params=CONTROL_PARAMS
)

# Setup HRMC parameters  
ENGINE.setup_hrmc(
    initial_temp=HRMC_INITIAL_TEMPERATURE,
    final_temp=HRMC_FINAL_TEMPERATURE, 
    steps=HRMC_STEPS,
    cycles=CYCLES_MD_HRMC
)

# Run MD-HRMC simulation
ENGINE.run_md_hrmc()

# Try to save final state
ENGINE.save()
ENGINE.export_pdb(os.path.join(OUTPUT_DIR, "final_structure.pdb"))

ENGINE.cleanup()
