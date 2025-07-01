"""
EnergyConstraint contains the ReaxFF energy constraint implementation.
This constraint calculates system energy using ReaxFF force field via LAMMPS.
"""

# Standard libraries imports
from __future__ import print_function
from collections import OrderedDict
import numpy as np
import os

# Fullrmc imports
from fullrmc.Core.Constraint import Constraint
from fullrmc.Core.Collection import is_number, reset_if_collected_out_of_date, get_path
from fullrmc.Globals import LOGGER, FLOAT_TYPE
from fullrmc.Core.boundary_conditions_collection import transform_coordinates


# Custom exceptions
class LAMMPSError(Exception):
    """LAMMPS-specific error."""
    pass

class BoxError(Exception):
    """Box definition error."""
    pass

# Try to import LAMMPS 
try:
    from lammps import lammps
    LAMMPS_IMPORT = True
except ImportError:
    LAMMPS_IMPORT = False

# Try to detect KOKKOS availability
KOKKOS_AVAILABLE = False
try:
    # Test if KOKKOS package is available by checking LAMMPS packages
    import subprocess
    result = subprocess.run(['lmp', '-help'], capture_output=True, text=True, timeout=5)
    if 'KOKKOS' in result.stdout:
        KOKKOS_AVAILABLE = True
except:
    pass

# Element masses (in g/mol)
ELEMENT_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.380
}

class KokkosReaxFFWrapper:
    """
    KOKKOS-enhanced wrapper for ReaxFF calculations with improved safety and error handling.
    
    This wrapper provides:
    - Atomic distance pre-checking to prevent bondchk failures
    - Gradual scaling mechanism to avoid sudden energy jumps
    - Enhanced error recovery with KOKKOS memory management
    - Performance optimizations using KOKKOS parallel computing
    """
    
    def __init__(self, lammps_instance):
        """Initialize KOKKOS wrapper with LAMMPS instance."""
        self.lmp = lammps_instance
        self.kokkos_enabled = False
        self.min_bond_distance = 0.5  # Minimum allowed bond distance in Angstroms
        self.scaling_factor = 0.8     # Factor for gradual position scaling
        self.max_scaling_steps = 5    # Maximum scaling attempts
        
    def enable_kokkos(self, device='cpu', threads=None):
        """Enable KOKKOS with specified device and thread configuration."""
        try:
            if not KOKKOS_AVAILABLE:
                LOGGER.warning("KOKKOS package not available, using standard LAMMPS")
                return False
                
            # Configure KOKKOS package
            kokkos_cmd = "package kokkos"
            if device:
                kokkos_cmd += f" {device}"
            if threads:
                kokkos_cmd += f" threads {threads}"
                
            self.lmp.command(kokkos_cmd)
            self.kokkos_enabled = True
            LOGGER.info(f"KOKKOS enabled with device: {device}")
            return True
            
        except Exception as e:
            LOGGER.warning(f"Failed to enable KOKKOS: {str(e)}")
            self.kokkos_enabled = False
            return False
    
    def check_atomic_distances(self, positions, elements):
        """
        Pre-check atomic distances to prevent bondchk failures.
        
        Parameters:
            positions (numpy.ndarray): Atomic positions
            elements (list): Element types
            
        Returns:
            bool: True if distances are safe, False otherwise
            list: Pairs of atoms that are too close (if any)
        """
        n_atoms = len(positions)
        close_pairs = []
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                distance = np.linalg.norm(positions[i] - positions[j])
                
                # Check minimum distance based on element types
                elem_i, elem_j = elements[i], elements[j]
                min_dist = self._get_minimum_distance(elem_i, elem_j)
                
                if distance < min_dist:
                    close_pairs.append((i, j, distance, min_dist))
        
        is_safe = len(close_pairs) == 0
        return is_safe, close_pairs
    
    def _get_minimum_distance(self, elem1, elem2):
        """Get minimum allowed distance between two elements."""
        # Covalent radii in Angstroms (approximate values)
        covalent_radii = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
            'P': 1.07, 'F': 0.57, 'Cl': 0.99, 'Fe': 1.32, 'Cu': 1.32
        }
        
        r1 = covalent_radii.get(elem1, 1.0)
        r2 = covalent_radii.get(elem2, 1.0)
        
        # Minimum distance is sum of covalent radii * safety factor
        return (r1 + r2) * 0.7  # 70% of sum for safety
    
    def apply_gradual_scaling(self, positions, elements, box_lengths):
        """
        Apply gradual position scaling to resolve atomic overlaps.
        
        Parameters:
            positions (numpy.ndarray): Original positions
            elements (list): Element types  
            box_lengths (numpy.ndarray): Box dimensions
            
        Returns:
            numpy.ndarray: Scaled positions
            bool: True if scaling successful, False otherwise
        """
        scaled_positions = positions.copy()
        
        for step in range(self.max_scaling_steps):
            is_safe, close_pairs = self.check_atomic_distances(scaled_positions, elements)
            
            if is_safe:
                LOGGER.info(f"Gradual scaling successful after {step} steps")
                return scaled_positions, True
            
            # Scale positions away from center of mass
            com = np.mean(scaled_positions, axis=0)
            
            for i, j, distance, min_dist in close_pairs:
                # Move atoms apart along their connecting vector
                vector = scaled_positions[j] - scaled_positions[i]
                vector_norm = np.linalg.norm(vector)
                
                if vector_norm > 0:
                    # Calculate required separation
                    required_separation = min_dist - distance + 0.1  # Add small buffer
                    separation_vector = vector / vector_norm * required_separation * 0.5
                    
                    # Move both atoms apart
                    scaled_positions[i] -= separation_vector
                    scaled_positions[j] += separation_vector
                    
                    # Ensure positions stay within box bounds
                    if box_lengths is not None:
                        scaled_positions[i] = np.mod(scaled_positions[i], box_lengths)
                        scaled_positions[j] = np.mod(scaled_positions[j], box_lengths)
        
        LOGGER.warning(f"Gradual scaling failed to resolve overlaps after {self.max_scaling_steps} steps")
        return scaled_positions, False
    
    def safe_reaxff_computation(self, positions, elements, box_lengths):
        """
        Safely compute ReaxFF energy with pre-checking and error recovery.
        
        Parameters:
            positions (numpy.ndarray): Atomic positions
            elements (list): Element types
            box_lengths (numpy.ndarray): Box dimensions
            
        Returns:
            tuple: (success, energy, final_positions, error_message)
        """
        try:
            # Step 1: Pre-check atomic distances
            is_safe, close_pairs = self.check_atomic_distances(positions, elements)
            
            if not is_safe:
                LOGGER.warning(f"Found {len(close_pairs)} close atomic pairs, applying gradual scaling")
                scaled_positions, scaling_success = self.apply_gradual_scaling(positions, elements, box_lengths)
                
                if not scaling_success:
                    return False, None, positions, "Failed to resolve atomic overlaps through scaling"
                
                # Use scaled positions for computation
                computation_positions = scaled_positions
            else:
                computation_positions = positions
            
            # Step 2: Return positions for actual LAMMPS computation
            # The actual energy computation will be done by the calling method
            return True, None, computation_positions, None
            
        except Exception as e:
            error_msg = f"Safe ReaxFF pre-computation failed: {str(e)}"
            LOGGER.error(error_msg)
            return False, None, positions, error_msg

class EnergyConstraint(Constraint):
    """
    ReaxFF energy constraint implementation.
    
    Parameters
    ----------
    reaxff_params : str
        Path to ReaxFF force field parameter file
    control_params : str  
        Path to ReaxFF control parameter file
    compute_forces : bool 
        Whether to compute atomic forces
    weight : float
        Constraint weight for standard error calculation
    """
    def __init__(self, reaxff_params="ffield.reax", control_params="control.reax", compute_forces=False, weight=1.0):
        # Initialize base class
        super(EnergyConstraint, self).__init__()

        # Check LAMMPS availability
        if not LAMMPS_IMPORT:
            raise ImportError(LOGGER.error("LAMMPS Python interface not found"))

        # Initialize attributes
        self.__lmp = None
        self.__initialized = False
        self.__element_types = {}
        self.__energy_cache = {}
        self.__cached_box = None
        self.__reaxff_params = None
        self.__control_params = None
        self.__compute_forces = None
        self.__last_accepted_energy = None
        self.__kokkos_wrapper = None  # KOKKOS wrapper instance

        # Set parameters
        self.set_reaxff_params(reaxff_params)
        self.set_control_params(control_params)
        self.set_compute_forces(compute_forces)

        # Initialize flags for move computation
        self._constraint_before_move_reset = False
        self._constraint_after_move_reset = False

        # Setup frame data by extending parent class tuples
        frame_data = list(self.FRAME_DATA)
        frame_data.extend([
            '_EnergyConstraint__element_types',
            '_EnergyConstraint__reaxff_params',
            '_EnergyConstraint__control_params',
        ])

        runtime_data = list(self.RUNTIME_DATA)

        # Set the class attributes using object.__setattr__
        object.__setattr__(self, 'FRAME_DATA', tuple(frame_data))
        object.__setattr__(self, 'RUNTIME_DATA', tuple(runtime_data))


    def _get_box_vectors(self):
        """Get simulation box vectors from engine."""
        try:
            # Get box vectors directly from engine boundary conditions 
            if self.engine.isPBC:
                box_vectors = self.engine.boundaryConditions.get_vectors()
                if box_vectors is None:
                    raise BoxError("No boundary conditions defined in engine")
                
                # Convert to box lengths
                box_lengths = np.array([np.linalg.norm(v) for v in box_vectors])
                
                # Validate box
                if np.any(box_lengths <= 0):
                    raise BoxError("Invalid box dimensions")
                    
                return box_lengths
            else:
                # For non-periodic systems, estimate box size from coordinates
                coords = self.engine.realCoordinates
                min_coords = np.min(coords, axis=0) 
                max_coords = np.max(coords, axis=0)
                box_lengths = (max_coords - min_coords) * 1.1 # Add 10% padding
                return box_lengths
                
        except Exception as e:
            LOGGER.error(f"Box vector calculation failed: {str(e)}")
            raise BoxError(f"Failed to get box vectors: {str(e)}")
    
    def _initialize_lammps(self):
        """Initialize LAMMPS with proper box dimensions and initial atoms."""
        try:
            if self.__initialized and self.__lmp is not None:
                return

            # Create LAMMPS instance
            logfile = os.path.join(self.engine.path, "lammps_screen.log")
            self.__lmp = lammps(cmdargs=["-log", "lammps.log", "-screen", logfile])
            lmp = self.__lmp

            # Initialize KOKKOS wrapper and try to enable KOKKOS
            self.__kokkos_wrapper = KokkosReaxFFWrapper(lmp)
            kokkos_enabled = self.__kokkos_wrapper.enable_kokkos(device='cpu')
            
            if kokkos_enabled:
                LOGGER.info("KOKKOS package enabled for enhanced ReaxFF safety")
            else:
                LOGGER.info("Using standard LAMMPS without KOKKOS")

            lmp.command("package omp 0")         # 使用 OMP，但默认不启用
            lmp.command("suffix off")            # 默认关闭所有并行加速
            
            # Basic settings
            lmp.command("units real")
            lmp.command("atom_style charge")
            lmp.command("boundary p p p")
            
            # Get elements and create type mapping
            elements = self.engine.allElements
            unique_elements = sorted(set(elements))
            LOGGER.info(f"Found elements: {unique_elements}")
            self.__element_types = {elem: idx+1 for idx, elem in enumerate(unique_elements)}

            # Get box dimensions from engine
            box_lengths = self._get_box_vectors()
            LOGGER.info(f"Creating simulation box with dimensions: {box_lengths}")
            
            # Create simulation box
            lmp.command(f"region box block 0 {box_lengths[0]} 0 {box_lengths[1]} 0 {box_lengths[2]} units box")
            lmp.command(f"create_box {len(unique_elements)} box")


            # NOW we can set output controls after box is defined
            lmp.command("echo none")     # Disable echo
            lmp.command("log none")      # Disable logging  
            lmp.command("thermo_style custom step pe") # Set thermo output
            lmp.command("thermo_modify flush no norm no") # Minimize output
            lmp.command("thermo 0")      # Disable thermo output

            # Set element masses
            LOGGER.info("Setting atomic masses...")
            for elem, type_id in self.__element_types.items():
                mass = ELEMENT_MASSES.get(elem, 0.0)
                lmp.command(f"mass {type_id} {mass}")
            
            # Setup ReaxFF pair style and neighbor settings            
            LOGGER.info("Setting up ReaxFF...")
            lmp.command(f"pair_style reaxff {self.__control_params} safezone 20 mincap 800")
            lmp.command("neighbor 3.0 bin")
            lmp.command("neigh_modify delay 0 every 1 check yes")

            # Setup pair coeffs
            elem_string = " ".join(unique_elements)
            lmp.command(f"pair_coeff * * {self.__reaxff_params} {elem_string}")

            self.__initialized = True
            self.__cached_box = box_lengths.copy()
            
            LOGGER.info("LAMMPS initialization completed successfully")

        except Exception as e:
            self.__initialized = False
            self.__lmp = None
            LOGGER.error(f"LAMMPS initialization failed: {str(e)}")
            raise LAMMPSError(str(e))
    
    def _verify_atom_count(self, positions, elements):
        """Verify atom counts match."""
        n_pos = len(positions)
        n_elem = len(elements)
        if n_pos != n_elem:
            LOGGER.error(f"Mismatch: {n_pos} positions vs {n_elem} elements")
            raise LAMMPSError(f"Atom count mismatch: {n_pos} != {n_elem}")
        return n_pos

    def _compute_configuration_energy(self, positions, box_lengths):
        """Compute ReaxFF energy for configuration with enhanced KOKKOS safety."""
        try:
            lmp = self.__lmp
            elements = self.engine.allElements

            # Verify atom counts before creation
            n_atoms = self._verify_atom_count(positions, elements)

            # Verify positions changed from previous
            if hasattr(self, '_last_positions'):
                delta = np.abs(positions - self._last_positions).max()
                LOGGER.info(f"Max position change: {delta}")
            self._last_positions = positions.copy()

            # KOKKOS-enhanced safety check: Pre-validate atomic distances
            safe_positions = positions.copy()
            if self.__kokkos_wrapper is not None:
                success, _, pre_checked_positions, error_msg = self.__kokkos_wrapper.safe_reaxff_computation(
                    positions, elements, box_lengths
                )
                
                if not success:
                    LOGGER.warning(f"KOKKOS pre-check failed: {error_msg}")
                    # Try with fallback error handling
                    raise LAMMPSError(f"Atomic distance pre-check failed: {error_msg}")
                else:
                    safe_positions = pre_checked_positions
                    if not np.array_equal(positions, safe_positions):
                        LOGGER.info("Applied KOKKOS gradual scaling to resolve atomic overlaps")

            # Try to preserve charges from previous state
            try:
                old_charges = lmp.extract_atom("q")
                # Verify charges are valid numbers
                if old_charges is not None:
                    old_charges = np.array(old_charges)
                    # Replace any NaN or infinite values with 0.0
                    old_charges = np.nan_to_num(old_charges, nan=0.0, posinf=0.0, neginf=0.0)
            except:
                old_charges = None          
            
            try:
                # Try to remove existing QEQ fix first if it exists
                lmp.command("unfix QEQ")
            except:
                # Ignore if fix doesn't exist
                pass            
            
            # Delete existing atoms and reset system
            lmp.command("delete_atoms group all")      

            # Reset box
            lmp.command(f"change_box all x final 0 {box_lengths[0]} y final 0 {box_lengths[1]} z final 0 {box_lengths[2]} units box")

            # Ensure positions are within box bounds before creating atoms
            # Use the safe_positions from KOKKOS pre-check
            wrapped_positions = safe_positions.copy()
            if self.engine.isPBC:
                # Wrap coordinates into primary box
                wrapped_positions[:, 0] = safe_positions[:, 0] % box_lengths[0]
                wrapped_positions[:, 1] = safe_positions[:, 1] % box_lengths[1] 
                wrapped_positions[:, 2] = safe_positions[:, 2] % box_lengths[2]

            # Create atoms with wrapped positions
            for idx, (pos, elem) in enumerate(zip(wrapped_positions, elements), start=1):
                atom_type = self.__element_types[elem]
                if not np.all(np.isfinite(pos)):
                    raise LAMMPSError(f"Invalid position for atom {idx}: {pos}")
                lmp.command(f"create_atoms {atom_type} single {pos[0]} {pos[1]} {pos[2]}")
                charge = old_charges[idx-1] if old_charges is not None else 0.0
                lmp.command(f"set atom {idx} charge {charge}")

            # Verify all atoms were created
            natoms = lmp.get_natoms()
            if natoms != n_atoms:
                raise LAMMPSError(f"Atom count mismatch after creation: {natoms} != {n_atoms}")

            # Update neighbor settings with KOKKOS-aware configuration
            if self.__kokkos_wrapper and self.__kokkos_wrapper.kokkos_enabled:
                lmp.command("neighbor 2.5 bin")  # Slightly smaller neighbor distance for KOKKOS
                lmp.command("neigh_modify every 1 delay 0 check yes page 100000")
            else:
                lmp.command("neighbor 3.0 bin")
                lmp.command("neigh_modify every 1 delay 0 check yes")

            # Enhanced error handling for ReaxFF operations
            try:
                # Apply QEQ with enhanced error recovery
                lmp.command("suffix off") # 禁用 OMP
                lmp.command("fix QEQ all qeq/reaxff 1 0.0 10.0 1e-6 reaxff maxiter 1000")
                lmp.command("run 0")  # Run QEQ without moving atoms
                lmp.command("unfix QEQ")

                # Apply minimize with gradual approach
                if self.__kokkos_wrapper and self.__kokkos_wrapper.kokkos_enabled:
                    lmp.command("suffix kokkos")  # Use KOKKOS acceleration if available
                else:
                    lmp.command("suffix omp")     # 启用 OMP
                    
                lmp.command("min_style cg")  # Use conjugate gradient
                lmp.command("min_modify dmax 0.1 line quadratic")            
                lmp.command("minimize 0.0e-4 1.0e-4 0 0") #0 0 is equivalent to not performing minimization, but only obtaining the energy
           
            except Exception as reaxff_error:
                # Enhanced error handling for known ReaxFF issues
                error_str = str(reaxff_error)
                if "bondchk failed" in error_str or "malformed" in error_str:
                    LOGGER.warning("ReaxFF bondchk error detected, attempting recovery")
                    # Apply more conservative neighbor settings
                    lmp.command("neighbor 4.0 bin")
                    lmp.command("neigh_modify every 1 delay 0 check yes page 200000")
                    # Retry with softer settings
                    try:
                        lmp.command("fix QEQ all qeq/reaxff 1 0.0 10.0 1e-5 reaxff maxiter 500")
                        lmp.command("run 0")
                        lmp.command("unfix QEQ")
                        lmp.command("minimize 0.0e-3 1.0e-3 0 0")
                    except:
                        raise LAMMPSError(f"ReaxFF computation failed even with recovery: {error_str}")
                else:
                    raise LAMMPSError(f"ReaxFF computation failed: {error_str}")
       
            # Get final positions and energy 
            x = lmp.extract_atom("x", 3)
            final_positions = np.array([[x[i][0], x[i][1], x[i][2]] 
                                    for i in range(len(positions))])

            # Get energy
            energy = lmp.extract_compute("thermo_pe", 0, 0)
        
            # Return both energy and minimized positions
            return float(energy), final_positions
        
        except Exception as e:
            LOGGER.error(f"Energy computation failed: {str(e)}")
            raise LAMMPSError(str(e))

    def _export_pdb(self, positions, elements, box_lengths, filename="debug.pdb"):
        """Export configuration as PDB with verification."""
        try:
            n_atoms = len(positions)
            with open(filename, 'w') as f:
                # Write header with counts
                f.write(f"REMARK   Generated by EnergyConstraint\n")
                f.write(f"REMARK   Total atoms: {n_atoms}\n")
                f.write(f"CRYST1 {box_lengths[0]:8.3f} {box_lengths[1]:8.3f} "
                    f"{box_lengths[2]:8.3f}  90.00  90.00  90.00 P 1\n")
                
                # Write atoms with validation
                atoms_written = 0
                for i, (pos, elem) in enumerate(zip(positions, elements), 1):
                    f.write(f"ATOM  {i:5d}  {elem:<3} MOL     1    "
                        f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                        f"  1.00  0.00          {elem:>2}\n")
                    atoms_written += 1
                    
                f.write("END\n")
                
            if atoms_written != n_atoms:
                raise LAMMPSError(f"PDB export mismatch: wrote {atoms_written} of {n_atoms}")
                
            LOGGER.info(f"Exported {atoms_written} atoms to {filename}")
            
        except Exception as e:
            LOGGER.error(f"Failed to export PDB: {str(e)}")
            raise

    def compute_energy(self):
        """Compute system energy using ReaxFF."""
        try:
            # Initialize LAMMPS if needed
            if not self.__initialized:
                self._initialize_lammps()
                
            # Always recompute for new coordinates
            positions = self.engine.realCoordinates
            box_lengths = self._get_box_vectors()

            # Get energy and final positions
            energy, final_positions = self._compute_configuration_energy(positions, box_lengths)
        
            # Update engine coordinates with minimized positions
            self.engine.realCoordinates[:] = final_positions

            # Transform to box coordinates using reciprocal basis vectors
            if self.engine.isPBC:
                # Convert final_positions to float32 before transformation
                final_positions_f32 = np.array(final_positions, dtype=np.float32)
                recip_basis_f32 = np.array(self.engine.reciprocalBasisVectors, dtype=np.float32)
                
                box_coords = transform_coordinates(
                    transMatrix=recip_basis_f32,
                    coords=final_positions_f32
                )
            else:
                box_coords = final_positions            
                
            # Ensure box coordinates are float32        
            self.engine.boxCoordinates[:] = np.array(box_coords, dtype=np.float32)

            # Return just the energy value
            return energy
                
        except Exception as e:
            LOGGER.error(f"Energy calculation failed: {str(e)}")
            raise LAMMPSError(str(e))

    def __getstate__(self):
        """Custom pickle handling"""
        state = self.__dict__.copy()
        # 移除不可序列化的LAMMPS实例和KOKKOS wrapper
        state['_EnergyConstraint__lmp'] = None
        state['_EnergyConstraint__kokkos_wrapper'] = None
        state['_EnergyConstraint__initialized'] = False
        state['_EnergyConstraint__energy_cache'] = {}
        # 确保保存重要属性
        state['_EnergyConstraint__last_accepted_energy'] = self.__last_accepted_energy
        return state
        
    def __setstate__(self, state):
        """Custom unpickle handling"""
        self.__dict__.update(state)
        self.__lmp = None
        self.__kokkos_wrapper = None
        self.__initialized = False
        self.__energy_cache = {}

    def initialize(self):
        """Initialize constraint before simulation."""
        try:
            # Compute initial energy
            energy = self.compute_energy()
            
            # Set initial data
            self.set_data(energy)
            self.referenceData = energy
            self.__last_accepted_energy = energy
            
            # Initialize move flags
            self._constraint_before_move_reset = False
            self._constraint_after_move_reset = False
            
            LOGGER.info(f"Energy constraint initialized: energy={energy}")
            
            return self

        except Exception as e:
            LOGGER.error(f"Energy constraint initialization failed: {str(e)}")
            raise LAMMPSError(f"Initialization failed: {str(e)}")
    
    def compute_data(self):
        """Compute energy data."""
        energy = self.compute_energy()
        self.set_data(energy)
        return energy

    def _on_collector_reset(self):
        """Reset constraint data when collector is reset."""
        self.__energy_cache = {}
        self.__cached_box = None
        if hasattr(self, '_current_step_computed'):
            del self._current_step_computed

    def _runtime_initialize(self):
        """Initialize runtime data."""
        pass

    def _runtime_on_step(self):
        """Called at each step during run."""
        pass

    def set_reaxff_params(self, params):
        """Set ReaxFF parameter file path."""
        if not isinstance(params, str):
            raise TypeError(LOGGER.error("reaxff_params must be a string"))
        if not os.path.exists(params):
            raise ValueError(LOGGER.error(f"ReaxFF parameter file not found: {params}"))
        self.__reaxff_params = params
        self.__initialized = False

    def set_control_params(self, params):
        """Set ReaxFF control parameter file path."""
        if not isinstance(params, str):
            raise TypeError(LOGGER.error("control_params must be a string"))
        if not os.path.exists(params):
            raise ValueError(LOGGER.error(f"ReaxFF control file not found: {params}"))
        self.__control_params = params
        self.__initialized = False

    def set_compute_forces(self, compute_forces):
        """Set force computation flag."""
        assert isinstance(compute_forces, bool)
        self.__compute_forces = compute_forces

    def compute_before_move(self, realIndexes, relativeIndexes):
        """
        Compute energy before move and store as reference.
        
        Parameters:
            realIndexes (numpy.ndarray): Real atom indexes
            relativeIndexes (numpy.ndarray): Relative atom indexes
        """
        try:
            # Only compute if not already done for this move
            if not self._constraint_before_move_reset:
                # Compute current energy
                energy = self.compute_energy()
                
                # Store as reference data
                self.set_data(energy)
                self.referenceData = energy
                
                # Reset flag
                self._constraint_before_move_reset = True
                
                LOGGER.info(f"Before move: energy={energy}")
                
        except Exception as e:
            LOGGER.error(f"Before move computation failed: {str(e)}")
            raise

    def compute_after_move(self, realIndexes, relativeIndexes, movedBoxCoordinates):
        """
        Compute energy after move.
        
        Parameters:
            realIndexes (numpy.ndarray): Real atom indexes
            relativeIndexes (numpy.ndarray): Relative atom indexes
            movedBoxCoordinates (numpy.ndarray): New box coordinates after move
        """
        try:
            # Compute new energy
            energy = self.compute_energy()
            
            # Store new energy
            self.set_data(energy)
            
            # Reset flag
            self._constraint_after_move_reset = True
            
            LOGGER.info(f"After move: energy={energy}")
            
            return True
            
        except Exception as e:
            LOGGER.error(f"After move computation failed: {str(e)}")
            return False

    def accept_move(self, realIndexes, relativeIndexes):
        """
        Accept the current move.
        
        Parameters:
            realIndexes (numpy.ndarray): Real atom indexes
            relativeIndexes (numpy.ndarray): Relative atom indexes
        """
        if self.data is not None:
            self.__last_accepted_energy = self.data
            self.increment_accepted()
        else:
            LOGGER.warning("接受移动时，能量数据为None")
            # 可能需要重新计算能量
            energy = self.compute_energy()
            self.__last_accepted_energy = energy
            self.set_data(energy)
        
        # 重置标志
        self._constraint_before_move_reset = False
        self._constraint_after_move_reset = False

    def reject_move(self, realIndexes, relativeIndexes):
        """
        Reject the current move.
        
        Parameters:
            realIndexes (numpy.ndarray): Real atom indexes
            relativeIndexes (numpy.ndarray): Relative atom indexes
        """
        # Reset flags
        self._constraint_before_move_reset = False 
        self._constraint_after_move_reset = False
