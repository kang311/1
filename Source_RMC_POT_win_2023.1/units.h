//header: units.h
//Last changed 01.03.2023

//define the global constants used throughout the program


#define PI 3.14159265359
#define SQRPI 1.7724539 //square root of Pi
#define INVPI 0.3183099 //inverse of Pi
#define SQRT3 1.7320508075688772 //sqrt(3) (maximum value in sftable)
#define k_BOLTZMANN 1.38064852e-23 //J/K
#define K_SHIFT_CONST 0.262468423789 //2*m_e/h_bar/h_bar [1/A/A/eV] for EXAFS shift calculation
//the tolerable difference coming from the different number representation in the
//binary and decimal number system (if the numbers are represented by 15 digits
//after the decimal point)
#define TOLERANCE 1.0e-15
//tolerable difference at some comparisions
#define TOLERANCE2 1.0e-13
//this is used to ensure the accuracy of the bin->dr conversion
#define GRID_TOL 1.0e-14
//this is a safety increase for array dimensions in NeighbouList object
#define SAFE_ADD 20;
//used during the load of CoordNumbConst and AvCoordConst
#define LOAD_TOL 1.0e-8
//defining the confidence interval used for the calculation of "negative" cosine distribution of bond angles constraints
#define CONF_INT 3
//size of a line for the line buffer for some file processing
#define LINE_SIZE 200
//The CACHE related things are architecture dependent, the given values and caching concept is for the Intel64 architecture.
#define NUMBER_OF_CACHE_LINES_TO_FETCH 1//Number of cache lines to cache in the same time (prefetch)
#define CACHE_LINE_SIZE 64 //Byte .Size of the L1 cache, needed in some cases to optimise cache usage. False sharing between 
						//threads has to be prevented. This can happen, when although different threads are writing different	
						//memory addresses, but the addresses are so close to each other, that they would be cached together
						//into the same cache line (are inside the same CACHE_ALIGNMENT block). Because of this, if one part
						//of a cache line is modified, the whole cache line is written back to memory, so different threads may
						//want to write the same part of the memory holding back each other causing if this happens too often
						//to slow the performance down.
#define CACHE_PADDING	NUMBER_OF_CACHE_LINES_TO_FETCH*CACHE_LINE_SIZE //the thread segments of some arrays have to be separated
						//at least with CACHE_PADDING-size(data_type) amount of bytes have to be kept between the threads segment data

// The size of the file names
#define FILE_NAME_SIZE 50
// The size of names
#define NAME_SIZE 100
//Max number of items to read from topology
#define INPUT_SIZE 50
//number of cosine distr method
#define N_COS_METHOD 4

//These are needed for the topology in case of MD-like molecules
//number of recognized GROMACS directives
#define N_GR_DIR 12
//number of recognized preprocessor directives
#define N_COMP_DIR 5
//number of force field types
#define N_FORCEF 7
//default value for the number of molecule types
#define N_MOLTYPE 5
//Maximum number of segments in the preprocessor arrays (number of topology and include top. files)
#define N_TOPFILE 5
//Maximum number of active (embedded) ifdef statements in a file in the same time 
#define N_ACT_IFDEF 5
//number of different potential types (harmonic bond and angle and Ryckaert-Bellemans dihedral)
#define N_POT_TYPE 3
//number of different dihedral functions
#define N_DIH_FUNCT
//for the Coulomb interaction (1/4/pi/epsilon0) (kJ*A/mol/e2)
#define f_Coulomb 1389.35485

//if atoms (or virtual sites) are too close, the potential can be extremely high, wiping out the meaningful digits
//during updating the potential.This is checked during the updating of the exclusions, where this is most likely to happen.
#define POT_WARNING 1.0e15
//to collect the potentials 
#define POT_SPLIT_HIGH 1.0e10
#define POT_SPLIT_LOW 1.0e-10

//AENET@echo MAC=0		_mac	for _CODE_WARRIOR_MAC
#define WRITE_ATOMIC_ENERGY_DEF 0 // do not write energy for all the atoms at each save
#define RELAX_DEF 0 //whether to relax the atoms inside cutoff
#define AENET_STEP_DEF 5//calculate ANN in each AENET_SET_DEF steps by default

//Default values
#define LAST_MOVE_DEF 10000//the run will end at ngenerated or naccepted = LAST_MOVE_DEF  by default
#define RUN_MODE 0//the run will be controlled by time limit
#define WRITE_LOG_DEF 0 //creating a logfile
#define DEBUG_DEF 0 //write debug info
#define LOG_NONLIN_STEPS_DEF 0 //write the iteration steps to logfile
#define WRITE_EXAFS_COEFFS_DEF 0 //write the EXAFS coeffs into *.expt, only used if E0 shift applied, as it can be a lot of matrix elements
#define R_SWITCH_POWER_DEF 2//if _R_SWITCH_GR_MIN_1 is on, this value is used 
#define SFACTOR_SIZE_DEF 501 //number of elements in the surface factor arrays

#ifdef _LOCAL_INV
	//step for resetting the thread boundaries for the neighbour atoms, if distance based calculation is used the 
	#define LOC_LOAD_BALANCE_STEP 1000
#endif

//Default values for free format data files
#define AUTO_CUTOFF_DEF 1
#define RSPACING_DEF  0.1
#define MAX_MOVES_DEF 0.1
#define MOVEOUT_DEF 1
#define CFG_COLL_DEF 0
#define CFG_COLL_FREQ_DEF 1
#define PRINTSTEP_DEF 10000
#define TIMELIM_DEF 300.0
#define TIMESAVE_DEF 2.0
#define TOOCLOSE_FRACTION_DEF 0.5  //the fraction the moved atoms are chosen from among the 'tooclose' atoms rather than from all the atoms
#define FNC_DEF 0
#define BINSHIFT_DEF 0.0
#define XMAX_DEF 1.0
#define RELOAD_DEF 0
#define HIST_BUFF_DEF 0
#define HIST_STEP_DEF 1
#define CUSTMOVE_DEF 0
#define MAX_GRIDATOM_DEF 5
#define NTHREADS_DEF 1
#define OLD_OUT_DEF 0
#define SUM_PPCF_DEF 0
#define NMOVED_DEF 1
#define LEAD_SERIES_IND_DEF 1
#define SWAP_FRACTION_DEF 0.0
#define E0_SHIFTSTEP_DEF 1000
#define AXS_SHIFTSTEP_DEF 1000
#define IQ_BACKG_CORR_STEP_DEF 1000

#define NO_IVALUE_DEF -999 //no value was supplied at optional data reading
#define NO_DVALUE_DEF -999.0 //no value was supplied at optional data reading
#define RENORM_DEF 0 //used for all the renormalization parameters, no renormalization by default
#define SUBTRACT_DEF 0.0
#define USE_RFACTOR_DEF 0
#define USE_CUBIC_DEF 0
#define EXAFS_CHI2_POWER_DEF 3
#define MAX_E0_SHIFT_DEF 0.0 //this means no shift, even if the DELTAE0_NGRID line is given 
#define EXAFS_NGRID_DEF 0  //this means no shift, even if the DELTAE0_NGRID line is given 
#define IQFIT_DEF 0//whether to fit I(Q) in case of x-ray data sets
#define IQ_ALPHA_DEF 0.0 //for the I(Q) fitting, initial guess
#define READ_COEFFS_DEF 0//whether to read the X-ray coefficients, from free format version #003, for #002 it is 1 (read coeffs) to be able to use the old #002 *.dat without modification
#define NONLIN_MAX_NIT_DEF 1000 //maximum number of iterations for the non-linear regression in case of I(Q) fit
#define NONLIN_EPSILON_DEF 0.001//stop iterating, if delta_par smaller than this 
#define NONLIN_LAMBDA_DEF 0.001 //lambda for nonlin regression
#define NONLIN_FACTOR_DEF 10.0 // multiply or divide lambda by this during iteration
#define NONLIN_TERMINATE_DEF 10//terminate the program after this many failed regression, counted after the last successful regression 
#define IQBACKGCORR_DEF 0//no backgroud correction for I(Q) by default
#define IQ_MU_DEF 0.7//initial correction factor for I(Q) background correction
#define IQ_DMU_MAX_DEF 0.2//initial correction factor for I(Q) background correction
#define COMPTON_DEF 1//whether to use Compton scattering contribution
#define BVS_B_DEF 0.37//default value for the BVS b parameter



#define POT_TYPE_DEF 1
#define POT_CHI_LOW_LIM_FRACTION_DEF 2.0 //this means that it is not set by default
#define POTENTIAL_DEF 0
#define LEAD_SERIES_IND2_DEF 1
#define NB_WEIGHT_MODE_DEF 2
#define VDW_COMB_RULE_DEF 0
#define CUTOFF_DEF 1.0//in reduced units
#define FUDGE14_DEF 1.0
#define LJ_REP_N_DEF 12

#define TEMPERATURE_DEF 298.0

#define WRITE_CNC_DETAIL_DEF 0

#define DCOSTHETA_DEF 0.05
#define DISTRIB_TYPE_DEF "GAUSSIAN"

#ifdef _LOCAL_INV
	#define LOC_CHI2_MODE_DEF 1
#endif
#ifdef _NO_PERIODIC
#define RECENTRE_FLAG_DEF 0
#endif



