//header global.h
//Last changed 24.01.2023

//Contains the file names and extensions to be available for all the source files
#include "units.h" 

#ifdef _OLD_HEADER 
	#include<iostream.h>
	#include<fstream.h>
	#include<iomanip.h>
#else
	#include<iostream>
	#include<fstream>
	#include<iomanip>
	#include <set>
	#include <unordered_set>
	#include <map>
	#include <string>
	#include <thread>
	#include <mutex>
	#include <condition_variable>
	#include<vector>
	#include<cctype>

	#include<math.h>
	#include<stdlib.h>
	#include<stdio.h>
	#include<string.h>
	#include<time.h>
	#include<sys/timeb.h>
	#include <chrono>
	#include <sstream>

	//These are needed for the new style header file, list cannot be used, as it can ve mixed with using declaration
	using std::ifstream;
	using std::ofstream;
	using std::istream;
	using std::ostream;
	using std::stringstream;
	using std::cin;
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::ios;
	using std::string;

#endif


struct tag_table_type
{
	std::string multiplicity;//how many times it can occure in a file
	std::string description;
};
struct key_table_type
{
	std::set<string> alternative_key;
	std::string status;//mandatory, optional
	std::set<std::string> tag_names;
	std::unordered_set<std::string> values;
	std::string default_value;
	std::string description;
};
//these are in global.cpp, as initializing the tables were not accepted here
extern const std::map<std::string, tag_table_type> tag_table;
extern const std::map<std::string, key_table_type> key_table;
extern std::string sigma_text[3];
extern std::string alternative_sigma_text[3];
extern std::string potential_types[11];




//files names
#ifdef _DEF_FILES


	//extern int aenet_nnb_max;
	extern char buf[LINE_SIZE];//for the topology reading
	extern std::streamoff next_line_pos;//for ReadThisLine, position of the next line had to be preserved in case an optional data will be read from the line, but its missing and no comment or anything else follows the mandatory data

	extern char filename[FILE_NAME_SIZE];
	extern char datfilename[FILE_NAME_SIZE+5];	
	extern char cfgfilename[FILE_NAME_SIZE+5], bincfgfilename[FILE_NAME_SIZE+5];
	extern char fncfilename[FILE_NAME_SIZE+5], gridfilename[FILE_NAME_SIZE+5];
	extern char outfilename[FILE_NAME_SIZE+5], hgmfilename[FILE_NAME_SIZE+5], hstfilename[FILE_NAME_SIZE+5], movefilename[FILE_NAME_SIZE+5];
	extern char tcafilename[FILE_NAME_SIZE + 5], psqfilename[FILE_NAME_SIZE + 5];
	extern char cncfilename[FILE_NAME_SIZE+5], cncdfilename[FILE_NAME_SIZE + 5], acnfilename[FILE_NAME_SIZE+5], pgrfilename[FILE_NAME_SIZE+5];
	extern char pfqfilename[FILE_NAME_SIZE+5], pfgfilename[FILE_NAME_SIZE + 5], pekfilename[FILE_NAME_SIZE+5], ppcffilename[FILE_NAME_SIZE+5];
	extern char cosfilename[FILE_NAME_SIZE+5], fitfilename[FILE_NAME_SIZE+5], chifilename[FILE_NAME_SIZE+5];
	extern char tabpotfilename[FILE_NAME_SIZE + 5],  logfilename[FILE_NAME_SIZE + 5];
	extern char sffilename[FILE_NAME_SIZE + 5];
	extern char logfilename1[FILE_NAME_SIZE + 5];
	extern char tempfilename[FILE_NAME_SIZE + 10], cfgcollectfilename[FILE_NAME_SIZE + 10], statefilename[FILE_NAME_SIZE + 10];
	extern char freerfilename[FILE_NAME_SIZE + 10], freersfilename[FILE_NAME_SIZE + 10], freefilename[FILE_NAME_SIZE + 10];
	
#ifdef _USE_LOCAL_INV
	extern char lhgmfilename[FILE_NAME_SIZE+5];
#endif
	extern char topfilename[FILE_NAME_SIZE+5];//topology file for molecules
	extern char exclfilename[FILE_NAME_SIZE+5];//file for exclusions from the vdW cacluation
	extern char potfilename[FILE_NAME_SIZE+5];//binary file for potentials
#ifdef _NEI
	extern char neighfilename[FILE_NAME_SIZE+5];//name of the neighbourlist file
#endif
#ifdef _NO_PERIODIC
	extern char posbinfilename[FILE_NAME_SIZE+5];//name of the central bin index for each atom file
#endif
#ifdef _AV_MOVE
	extern char avmfilename[FILE_NAME_SIZE+5];//name of the average move file
#endif
#ifdef _ADVANCED_GEOM_CONST
	extern char concfilename[FILE_NAME_SIZE + 5];//common neighbour constraint
	extern char sncfilename[FILE_NAME_SIZE + 5];//second neighbour constraint
	extern char bvsfilename[FILE_NAME_SIZE + 5];//bond valence sum constraint
	extern char bbvfilename[FILE_NAME_SIZE + 5];//bond valence sum constraint binary
#endif
#ifdef _AENET
	extern char energyfilename[FILE_NAME_SIZE + 5];//name of the file to save E_i
	extern const char *enerext;//aenet energy
#endif
	extern const char *datext;//text type data file
	extern const char *freeext;//text type free format default data file for output
	extern const char *freerext;//text type free reduced format data file for output
	extern const char *freersext;//text type free reduced  format data file (after scalable sigma calculated) for output
	extern const char *cfgext;//text type coordinates
	extern const char *bcfext;//binary coordinates
	extern const char *fncext;//FNC
	extern const char *outext;//RMCA format output file
	extern const char *hstext;//history
	extern const char *hgmext;//histogram
	extern const char *movext;//cus file for molecular moves
	extern const char *pgrext;//partial g(r)-s
	extern const char *psqext;//partial S(Q)-s
	extern const char *pfqext;//partial F(Q)-s
	extern const char *pfgext;//partial F(g)-s
	extern const char *pekext;//partial E(k)-s
	extern const char *tcaext;//too close atoms
	extern const char *cncext;//coordination constraint
	extern const char *cncdext;//detailed coordination constraint
	extern const char *acnext;//average coordination constraint
	extern const char *gridext;//grid
	extern const char *ppcfext;//ppcf
	extern const char *cosext;//cosine distribution of bond angles
	extern const char *fitext;//new format output file (calculated and experimental totals)
	extern const char *chiext;//chi2 file for start and in case of _WRITE_CHI2_DETAIL for rejected moves	
	extern const char *stateext;//sigma values outputed for the data sets and the constraints
	extern const char *logext;//for the log file
#ifdef _USE_LOCAL_INV
	extern const char *lhgmext;//local histogram
#endif
#ifdef _NEI
	extern const char *neiext;//neighbourlist file
#endif
#ifdef _NO_PERIODIC
	extern const char *binext;//the central bin index for each atom
#endif
	extern const char *cos_method[N_COS_METHOD];//cos dist method 
	extern const char *tabpotext;//writing tabulated potential
	extern const char *topext;//topology file for molecules
	extern const char *exclext;//exclusion file in case of vdW calculation
	extern const char *potext;//binary format potentials
	extern const char *directives[N_GR_DIR+1];
	extern const char *comp_dir[N_COMP_DIR];
	extern const char *force_field[N_FORCEF];
	extern const char *interaction_types[N_POT_TYPE];
#ifdef _AV_MOVE
	extern const char *avmext;//the average moves per each atom
#endif
#ifdef _ADVANCED_GEOM_CONST
	extern const char *concext;//common neighbour constraint
	extern const char *sncext;//second neighbour constraint
	extern const char *bvsext;//bond valence sum constraint
	extern const char *bbvext;//binary bond valence sum constraint
#endif
	extern std::ofstream logfile;
	extern std::ofstream logfile0;
	extern std::ofstream logfile1;
	extern std::ofstream logfile2;
	
	extern int debug;//write debud info
	extern int write_log;//creating a logfile
	extern int log_nlr_steps;//write the iteration steps during non-linear regression
	extern int write_exafs_coeffs;//write the coeffs matrices into *.expt in case of E0 shift
	extern int r_switch_power;//for the  _R_SWITCH_GR_MIN_1 option
	extern int note;//for the notes
	extern int warn;//for the warnings
	extern int end_flag;//0:linux "\n" dec: 10, 1: win "\r\n" dec: 13 10, 2:mac "\r" dec: 13
#else
	//int aenet_nnb_max;
	char buf[LINE_SIZE];//for reding the topology
	std::streamoff next_line_pos=-1;//for ReadThisLine, position of the next line had to be preserved in case an optional data will be read from the line, but its missing and no comment or anything else follows the mandatory data
    std::ofstream logfile;
	std::ofstream logfile0;
	std::ofstream logfile1;
	std::ofstream logfile2;
	char filename[FILE_NAME_SIZE];
	char datfilename[FILE_NAME_SIZE+5];
	char cfgfilename[FILE_NAME_SIZE+5],bincfgfilename[FILE_NAME_SIZE+5];
	char fncfilename[FILE_NAME_SIZE+5], gridfilename[FILE_NAME_SIZE+5];
	char outfilename[FILE_NAME_SIZE+5], hgmfilename[FILE_NAME_SIZE+5], hstfilename[FILE_NAME_SIZE+5],movefilename[FILE_NAME_SIZE+5];
	char cfgcollectfilename[FILE_NAME_SIZE+10],tcafilename[FILE_NAME_SIZE+5],statefilename[FILE_NAME_SIZE+10], psqfilename[FILE_NAME_SIZE + 5];
	char cncfilename[FILE_NAME_SIZE+5], cncdfilename[FILE_NAME_SIZE + 5], acnfilename[FILE_NAME_SIZE+5],pgrfilename[FILE_NAME_SIZE+5];
	char pfqfilename[FILE_NAME_SIZE+5], pfgfilename[FILE_NAME_SIZE + 5], pekfilename[FILE_NAME_SIZE+5],ppcffilename[FILE_NAME_SIZE+5];
	char cosfilename[FILE_NAME_SIZE+5], fitfilename[FILE_NAME_SIZE+5],chifilename[FILE_NAME_SIZE+5];
	char tabpotfilename[FILE_NAME_SIZE + 5], tempfilename[FILE_NAME_SIZE+10], logfilename[FILE_NAME_SIZE + 5];
	char freerfilename[FILE_NAME_SIZE + 10], freersfilename[FILE_NAME_SIZE + 10], freefilename[FILE_NAME_SIZE + 10];
	char sffilename[FILE_NAME_SIZE + 5]="";
	char logfilename1[FILE_NAME_SIZE + 5];
#ifdef _USE_LOCAL_INV
	char lhgmfilename[FILE_NAME_SIZE+5];
#endif
	char topfilename[FILE_NAME_SIZE+5];//topology file for molecules
	char exclfilename[FILE_NAME_SIZE+5];//file for exclusions from the vdW cacluation
	char potfilename[FILE_NAME_SIZE+5];//binary file for potentials
#ifdef _NEI
	char neighfilename[FILE_NAME_SIZE+5];//name of the neighbourlist file
#endif
#ifdef _NO_PERIODIC
	char posbinfilename[FILE_NAME_SIZE+5];//name of the central bin index for each atom file
#endif
#ifdef _AV_MOVE
	char avmfilename[FILE_NAME_SIZE+5];
#endif 
#ifdef _ADVANCED_GEOM_CONST
	char concfilename[FILE_NAME_SIZE + 5];//common neighbour constraint
	char sncfilename[FILE_NAME_SIZE + 5];//second neighbour constraint
	char bvsfilename[FILE_NAME_SIZE + 5];//bond valence sum constraint
	char bbvfilename[FILE_NAME_SIZE + 5];//binary bond valence sum constraint
#endif
#ifdef _AENET
	char energyfilename[FILE_NAME_SIZE + 5];//name of the file to save E_i
	const char *enerext=".en";//aenet energy
#endif
	const char *datext=".dat";//text type data file
	const char *freeext=".free";//text type free format complete data file for output
	const char *freerext=".freer";//text type free reduced format data file for output
	const char *freersext=".freers";//text type free reduced  format data file (after scalable sigma calculated) for output
	const char *cfgext=".cfg";//text type coordinates
	const char *bcfext=".bcf";//binary coordinates
	const char *fncext=".fnc";//FNC
	const char *outext=".out";//RMCA format output file
	const char *hstext=".hst";//history
	const char *hgmext=".hgm";//histogram
	const char *movext=".cus";//cus file for molecular moves
	const char *pgrext=".pgr";//partial g(r)-s
	const char *psqext=".psq";//partial S(Q)-s
	const char *pfqext=".pfq";//partial F(Q)-s
	const char *pfgext=".pfg";//partial F(g)-s
	const char *pekext=".pek";//partial E(k)-s
	const char *tcaext=".tca";//too close atoms
	const char *cncext=".cnc";//coordination constraint
	const char *cncdext = ".cncd";//coordination constraint detail
	const char *acnext=".acn";//average coordination constraint
	const char *gridext=".grid";//grid
	const char *ppcfext=".ppcf";//ppcf
	const char *cosext=".cos";//cosine distribution of bond angles
	const char *fitext=".fit";//new format output file (calculated and experimental totals)
	const char *chiext=".chi";//chi2 file for start and in case of _WRITE_CHI2_DETAIL for rejected moves	
	const char *stateext=".state";//sigma values outputed for the data sets and he constraints
	const char *logext=".log";//for the log file
#ifdef _USE_LOCAL_INV
	const char *lhgmext=".lhgm";//local histogram
#endif
#ifdef _NEI
	const char *neiext=".nei";//neighbourlist file
#endif
#ifdef _NO_PERIODIC
	const char *binext=".bin";//the central bin index for each atom
#endif
	const char *cos_method[N_COS_METHOD] = { "UNIFORM","GAUSSIAN","ABSENT","file"};//cos dist method 
	const char *tabpotext=".tabp";//writing tabulated potential
	const char *topext=".top";//topology file for molecules
	const char *exclext=".excl";//exclusion file in case of vdW calculation
	const char *potext=".pot";//binary format potentials
	const char *directives[N_GR_DIR+1]={"system", "molecules", "moleculetype", "atoms", "bonds", "angles", "dihedrals", "pairs", "exclusions", "virtual_sites2", "virtual_sites3", "virtual_sites4","invalid"};
	//Has to follow in the same order, as in the enumerated type T_comp_dir
	const char *comp_dir[N_COMP_DIR]={"include", "ifdef", "else", "endif","ifndef"}; 
	//Forcefields used in GROMACS
	const char *force_field[N_FORCEF]={"oplsaa", "encads","encadv", "charmm","amber","UNKNOWN_GENERATE_PAIRS","UNKNOWN_READ_PAIRS" };
	//Interaction types, has to match the T_Interaction_type
	const char *interaction_types[N_POT_TYPE]={"BOND","ANGLE","DIHEDRAL"};
#ifdef _AV_MOVE
	const char *avmext=".avm";//the average moves per each atom
#endif
#ifdef _ADVANCED_GEOM_CONST
	const char *concext=".conc";//saving the common neighbour constraint
	const char *sncext = ".snc";//saving the second neighbour constraint
	const char *bvsext = ".bvs";//saving the bond valence sum constraint
	const char *bbvext = ".bbv";//saving the binary bond valence sum constraint
#endif
	int debug = DEBUG_DEF;//write debud info
	int write_log = WRITE_LOG_DEF;//creating a logfile
	int log_nlr_steps = LOG_NONLIN_STEPS_DEF;//write the iteration steps during non-linear regression
	int write_exafs_coeffs = WRITE_EXAFS_COEFFS_DEF;//write the coeffs matrices into *.expt in case of E0 shift
	int r_switch_power=R_SWITCH_POWER_DEF;//for the  _R_SWITCH_GR_MIN_1 option
	int note = 0;//for the notes
	int warn = 0;//for the warnings
	int end_flag=-2;//0:linux "\n" dec: 10, 1: win "\r\n" dec: 13 10, 2:mac "\r" dec: 13
#define _DEF_FILES
	
	

#endif


	


	