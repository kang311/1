//header Topology.h
//Last changed 12.11.2021

//containing the molecular topoloy following the GROMACS convention

#include "utilities.h"
#include <list>
//use it only in case of molecular RMC with flexible, MD-like molecules
	//The molecules are defined following the GROMACS convention, the GROMACS topology file can 
	//be used with some additional parameters specified at the line ends after the text qualifier ;
	//to bind the GROMACS topology to the RMC configuration

//These numbers has to match the values in the enum lists
//forward declarations
class FNC_POT;
class ChiSquared;

//The virtual sites will be added to the end of the RMC configuration.
//Each virtual site declared in [virtual_sitesX ] will create a new virtual type block with as many virtual atoms,
//as the number of molecules in the given type. The virtual site block will be in the order the virtual sites are declared in the 
//topology [ virtualX ] section, which is not necessarily is the same, as their order in the atom section, so the virtual sites has to be
//identified by their .atom_index field!
//for example a configuration with 10 atoms, that has 3 atom types (A, B, C) and 2 virtual sites (V, W) with stochiometry A2BC3VW:
//type1 atoms:			20 lines
//type2 atoms:			10 lines
//type3 atoms:			30 lines
//vtype1 virtual sites:	10 lines
//vtype2 virtual sites:	10 lines
struct virtual_site_struct
{
	int tnumb;//number of atoms making the site
	int type;//gromacs type of the virtual site
	int atom_index;//gromacs index of the virtual site
	int host_index;//index of the host atom, whose type will decide which partial the potential contributes to
	int indices[4];//gromacs indices of the atoms making up the virtual site
	double params[3];//parameters for the virtual site
	short int assigned;//atom index of the virtual site in [atoms], for which the virtual site from [virtual_sitesX] was assigned to
	
};
//This is the class to hold the molecular topology
class Topology
{
	public:
	Topology();

	//regarding the topology
	static int nexclusion_max;//the largest exclusion number for all the molecule types
	static int nmoltype;//number of molecule types according to the number of [ moleculetype ] directive
	static int nsysmoltype;//number of molecule types given after the [ molecules ] directive
	static int n_gr_atoms;//total number of different atoms for gromacs in each different molecule
	static int virtual_count;//total number of virtual sites in the [atoms ] section in each different molecule
	static int n_gr_virtuals;//total number of virtual sites in [ virtual_siteX ] directive 
	static int nGRtypes;//number of different GROMACS atom types for vdW
	static int ncharge_types;//number of different charge types
	static int ncharge_groups;//totla number of different charge groups in the Toplogy files
	static int nRMC_charge_centre;//number of RMC charge group    
	static int npairs;//total number of 1-4 pairs
	static int nbonds;//total number of bonds
	static int nangles;//total number of angles
	static int ndihedrals;//total number of dihedrals
	static int nperdihedrals;//total number of GROMACS type (1) periodic dihedrals
	static int nharmdihedrals;//total number of GROMACS type (2) harmonic dihedrals 
	static int nRBdihedrals;//total number of RBdihedrals GROMACS type (3)
	static int npair_types;//number of different 1-4 pair types
	static int nexclusionlines;//number of exclusion lines in topology
	static int nbond_types;//number of different bond types
	static int nangle_types;//number of different angle types
	static int nperdihedral_types;//number of different periodic dihedral types
	static int nharmdihedral_types;//number of different harmonic dihedral types
	static int nRBdihedral_types;//number of different RBdihedral types
	static int bond_weight_mode;//whether the same weight parameter is used for all the bond types, output is collapsed to one entry
	static int angle_weight_mode;//whether the same weight parameter is used for all the angle types, output is collapsed to one entry
	static int perdih_weight_mode;//whether the same weight parameter is used for all the perdihedral types, output is collapsed to one entry
	static int harmdih_weight_mode;//whether the same weight parameter is used for all the harmdihedral types, output is collapsed to one entry
	static int RBdih_weight_mode;//whether the same weight parameter is used for all the RBdihedral types, output is collapsed to one entry
	static int last_charge_group;//the last charge group number
	static int last_directive;//the last GROMACS directive
	static int gen_pairs;//whether generate pairs based on the combination rule
	static int *nexclusions_per_type;//number of exclusions per moleculetype
	static int *nmol_per_type;//number of molecules per molecula type
	static int *ncharge_group_per_type;//number of charge groups per type
	static int *natoms_per_type;//number of atoms+virtual sites in one molecule /molecule type
	static int *cumul_atoms;//cumulative number of atoms+virtual sites in one molecule /molecule type
	static int *nvirtuals_per_type;//number of virtual sites in one molecule /molecule type
	static int *cumul_virtuals;//cumulative number of virtual sites in one molecule /molecule type
	static int *first_index;//atom index of the given atom in the first molecule of the given residue (1-natoms) (starting with 1)
	static int *second_index;//atom index of the given atom in the first molecule of the given residue (1-natoms) (starting with 1)
	static int *delta_index;//difference between second and first_index
	static int *RMCindex_offset;//offset for the RMC_index of the molecules, useful if it is not the first in the system
	static int *npairs_per_type;//number of 1-4 pairs per molecule type
	static int *cumul_pairs;//cumulative number of 1-4 pairs in one molecule /type
	static int *pair_index;//indices of the pairs in pairs
	static int *pair_type;//type of the pair
	static int *nexclusionlines_per_type;//number of exlusion lines per molecule type
	static int *cumul_exclusionlines;//cumulative number of exclusion lines in a molecule type/molecule type
	static int *exclusion_centre;//index of the central atom for the exclusion
	static int *nexcluded_neighs;//number of excluded atoms for the central atom 
	static int *cumul_nexcluded_neighs;//cumulative number of excluded atoms for the central atom, for each exclusion line
	static int *exclusions;//indices of the excluded atoms for an atom in the exclusion_centre array;
	//there can be deletions from exclusions array, and the empty places filled with non-emty ones form the same exclusion lines,
	//and nexcluded_neighs will be updated, but cumul_nexcluded_neighs array will not be changed, use always explicit indices in exclusions array using cumul_nexcluded_neighs,
	//and not go through with a pointer
	static int *nbonds_per_type;//number of bonds per molecule type
	static int *nbondtypes_per_type;//number of bond types per molecule type
	static int *cumul_bonds;//cumulative number of bonds in a molecule type/molecule type
	static int *cumul_bondtypes;//cumulative number of bond types in a molecule /molecule type
	static int *bond_index;//indices of the bonds in pairs
	static int *bond_type;//type of the bond
	static int *zero_bond_sigma;//size:nmoltype+1, [0] element is 0: no zero sigma at all, 1: there is zero sigma, then index of the first non-zero sigma interaction for molecule type
	static int *nangles_per_type;//number of angles per molecule type
	static int *nangletypes_per_type;//number of angle types per molecule type
	static int *cumul_angles;//cumulative number of angles in a molecule type/molecule type
	static int *cumul_angletypes;//cumulative number of angle types in a molecule /molecule type
	static int *angle_index;//indices of the angles in triplets (neighbour-central-neighbour)
	static int *angle_type;//type of the angle
	static int  zero_angle_sigma;//0: no zero sigma at all, 1: there is zero sigma
	static int *dihedral_index;//indices of the dihedrals in quartetts (1-2-3-4), for all the three dihedral type
	static int *dihedral_type;//type of all the dihedral similarly to GROMACS: 
				//(1) for peiodic GROMOS, (2) harmonic, (3) RB
	static int *ndihedrals_per_type;//total number of dihedrals per molecule type
	static int *nperdihedrals_per_type;//number of perdihedrals per molecule type
	static int *nperdihtypes_per_type;//number of periodic dihedral types per molecule type
	static int *cumul_perdihedrals;//cumulative number of perdihedrals in a molecule type/molecule type
	static int *cumul_perdihtypes;//cumulative number of periodic dihedral types in a molecule /molecule type
	static int *perdihedral_type;//type of the perdihedral
	static int  zero_perdih_sigma;//0: no zero sigma at all, 1: there is zero sigma
	static int *nharmdihedrals_per_type;//number of harmonic dihedrals per molecule type
	static int *nharmdihtypes_per_type;//number of harm dihedral types per molecule type
	static int *cumul_harmdihedrals;//cumulative number of harmonic dihedrals in a molecule type/molecule type
	static int *cumul_harmdihtypes;//cumulative number of harmomic dihedral types in a molecule /molecule type
	static int *harmdihedral_type;//type of the harmonic dihedral among the harmonic dihedrals
	static int  zero_harmdih_sigma;//0: no zero sigma at all, 1: there is zero sigma
	static int *nRBdihedrals_per_type;//number of RBdihedrals per molecule type
	static int *nRBdihtypes_per_type;//number of RB dihedral types per molecule type
	static int *cumul_RBdihedrals;//cumulative number of RBdihedrals in a molecule type/molecule type
	static int *cumul_RBdihtypes;//cumulative number of RB dihedral types in a molecule /molecule type
	static int *RBdihedral_type;//type of the RBdihedral
				//for GROMOS-like force fields (1) periodic proper, (2)harmonic improper (3) RB		
				//for OPLS (1) periodic improper, (3) RB (normal dihedral)
	static int	   zero_RBdih_sigma;//0: no zero sigma at all, 1: there is zero sigma
	static int    *perdihedral_multiplicity;//multiplicity for the periodic dihedrals
	static int	  *GROMACS_type;//the GROMACS atom type of the atom for each GROMACS atom
	static int	  *charge_type;//the charge type of the atom for each GROMACS atom
	static int	  *charge_group;//the charge group of the atom for each GROMACS atom
	static int	  *natoms_per_charge_group;//number of atoms for each charge group
	
	static double bond_r_max;//maximum bond length
	static double *vdW14_pot1;//LJ_sigma in A or C(6)=4*epsilon*(sigma^6) in kJ/mol*A^6, or Buckingham A
	static double *vdW14_pot2;//LJ_epsilon in kJ/mol or C(N)=4*epsilon*(sigma^N) in kJ/mol*A^N or Buckingham B
	static double *vdW14_pot1_partial;//the vdW14_pot1 rendered to GROMACS type
	static double *vdW14_pot2_partial;//the vdW14_pot2 rendered to GROMACS type
	static double *bond_r;//equilibrium value for bonds in Angstrom
	static double *bond_k;//force constant for bonds (read as kJ/mol/nm2 and converted to kJ/mol/A2)
	static double *bond_sigma;//sigma values for the bond types
	static double *angle_ang;//equilibrium value for the angles in degrees 
	static double *angle_k;//force constant for angles (kJ/mol/rad2)
	static double *angle_sigma;//sigma values for the angle types
	static double *perdihedral_ang;//equilibrium value for the (GR type 1 periodic) dihedral angle in degrees 
	static double *perdihedral_k;//force constant for (GR type 1 periodic) dihedrals (kJ/mol/rad2)
	static double *harmdihedral_ang;//equilibrium value for the GR type 2 harmonic dihedral angle in degrees 
	static double *harmdihedral_k;//force constant for GR type 2 harmonic dihedral (kJ/mol/rad2)
	static double *RBdihedral_C0;//C0 force constant for Ryckaert-Bellemans RBdihedrals (kJ/mol)
	static double *RBdihedral_C1;//C1 force constant for Ryckaert-Bellemans RBdihedrals (kJ/mol)
	static double *RBdihedral_C2;//C2 force constant for Ryckaert-Bellemans RBdihedrals (kJ/mol)
	static double *RBdihedral_C3;//C3 force constant for Ryckaert-Bellemans RBdihedrals (kJ/mol)
	static double *RBdihedral_C4;//C4 force constant for Ryckaert-Bellemans RBdihedrals (kJ/mol)
	static double *RBdihedral_C5;//C5 force constant for Ryckaert-Bellemans RBdihedrals (kJ/mol)
	static double *perdihedral_sigma;//sigma values for the periodic dihedral types
	static double *harmdihedral_sigma;//sigma values for the harmonic dihedral types
	static double *RBdihedral_sigma;//sigma values for the RB dihedral types
	static double *charge;//charge of the atoms/per GROMACS atom types
	static double *virtual_params;//the parameters for the virtual sites
	static char (*molname)[NAME_SIZE];//name of the molecules
	static char (*atom_name)[10];//name of the atoms
	static char (*atom_type)[20];//type of the atoms
	static char *currentfilename;//pointer to the name of the file currently handled
	static char *systemname;//name of the GROMACS system
	static char *force_field_name;//name of the force field
	
	static const char *virtual_types[6];//the virtual site types
	static T_force_field force_field_type;//enumerated type of the force_field
	static bool RB_used;//indicator, that RB dhihedral is used. It turned out, that in case of RB, the dihedral potential can
				//be negative for certain dihedral types, which means, that using flexible molecules without NB are not entirely safe. In this case 
				//the bonded potentials contribute to the normal chi2. The proram will check the sign of the total bonded potential
				//at every chi2 calculation, and if it is negative, it will terminate...
				//As far as I know the original OPLS dihedral parameters were set the way, it is above zero, it may be a 
				//GROMACS specific feature, when they transformed the parameters...

	static virtual_site_struct *virtual_site;//array containing the information about the virtul sites
	//for handling the preprocessor directives
	static int nsegment;//number of segments in the preprocessor dir arrays
	static int *nifdef_ifndef;//current number of ifdef and ifndef preproc directives for each segment
	//there will be a segment for each file in this arrays
	static int *pifdef_ifndef;//indicator array of ifdef and ifndef preprocessor directives present in the topology ind itp file
	static int *pelse;//indicator array of else preprocessor directives present in the topology ind itp file
		
	
	static T_interaction_type interaction_type;//the type of the interaction between particles

	static void GetParams();//getting the static members
	static void SetParams();//setting the static members
	static bool IncludeLine(int nfiles);
	static int FindDirective(ifstream &file );//determine the compiler option or directives
	static int GetDirective(ifstream &file, longint file_pos, int n_dir, const char **dir);//get the option or directive
	static int SearchCOD(char *buf, longint &offset);//determining, whether the line contains compiler option, or GROMACS directive
	static int CheckForcefield(ifstream &file,char *name);
	static int CheckPreprocDir(ifstream &file,int nfiles, int flag);
	static int FindEmtyLine(char *buf);
	static int ExtractSemiColon(ifstream &file, int mode=1);//extracting the textqualifier, exit, if it is not found in case of mode=1, no exit, if it is something else
	static void Reset(int nfiles);//resetting some arrays after a file is closed
	static void PrintParams(ostream &file);//writing the parameters
	static void CheckRMCType();//check, whether an atom will have the same molecule type in all the instances of a molecule
	static void CreatevdW14Params();//rendering the vdW14 parameters to GROMACS types
	static void CheckSigma(int *weight_mode, int *cumul_bondtypes, int *ninttypes_per_type, int *zero_sigma, double *sigma, const char *interaction_name, int *leading_sigma);
	static void CheckRepetition(int *array, int count, int line, int moltype, char *directive);//check, whether all the entries are different
	static void CheckVirtualSites();//check, whether all the virtual sites of the [ atoms ] section are defined																							   //checking, if there is nonzero sigma, if there is 0 sigma for an interaction type of a molecule type
	static void CheckExclusions();//Check, whether the exclusions given in the [ exclusions ] directive are selfconsistent
	static double GetSigma(ifstream &file, const char *varname, const char *routine_name, int mode=1);//do not terminate, if sigma is not found, and mode=0
	
	//reading function for the different line types, they will be invoked throuh the read_function pointer array
	static int GetSystem(ifstream &file);
	static int GetMolecules(ifstream &file);
	static int GetMoleculetype(ifstream &file);
	static int GetAtoms(ifstream &file);
	static int GetBonds(ifstream &file);
	static int GetAngles(ifstream &file);
	static int GetDihedrals(ifstream &file);
	static int GetPairs(ifstream &file);
	static int GetExclusions(ifstream &file);
	static int GetVirtual2(ifstream &file);
	static int GetVirtual3(ifstream &file);
	static int GetVirtual4(ifstream &file);
	
	//Routines used to restore sigma parameters read from files (which are present according to RMC_POT 1.1 specification)
	//static void GetIncludes(const std::list<string> &def_opts, const string &filename, ifstream &file, std::list<string> &includes);//Returns with valid include filenames (appearing once) if you provide options and an already open file and its name
	static void GetBADSigmas(const std::list<string> &def_opts, const string &filename, ifstream &file, std::vector<double> &bond,
				 std::vector<string> &bond_tags, std::vector<double> &angle, std::vector<string> &angle_tags, std::vector<double> &dih,
			  std::vector<string> &dih_tags);//Extracts bond, angle and dih(edral angle) sigmas and tags (atomic no.'s) from the already opened
							      //"filename" "file" stream after preprocessing directives processed
	static int GetNoOfAtoms(const std::list<string> &def_opts, const string &filename);//Returns with the number of different atoms defined in topologies
		
	static int (*read_function[N_GR_DIR])(ifstream &file);//will hold the pointer to the above declared functions
	static void ResizeVirtualSiteStruct(int *max, int new_max, virtual_site_struct **array, const char *array_name, const char *routine_name);//resize the virtual site structure array
	static void O_Sort(int length, virtual_site_struct *array);//sorting the virtual sites acoording to their atom index
	static void CopyVirtual(virtual_site_struct *dest, virtual_site_struct *source);//for the sort
	//determining the number of interactions the atoms are involved in for the given interaction type 
	void CreateList(int pottype, FNC_POT &pfnc);
	void CheckAtomIndex(int i, int imol, int imoltype, int iatom);//checking whether the atomindex is not greater than the total number of atoms
	//reading the sigma value for the interaction
	
 // private:
	//static void createFilteredList(const std::list<string> &def_opts, const bool omit_force_field, const string &filename, ifstream &file,
			//	       std::list<string> &filtered_list);
	//Reads the already opened and checked "file" of "filename" and write its content into a list. It filters out preprocessing directives according to "def_opts" and
	//force field include (looks for '"ff' or '.ff' string) definition, if omit_force_field parameter is on
};
