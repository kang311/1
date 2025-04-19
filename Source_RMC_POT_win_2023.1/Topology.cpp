//source Topology.cpp
//Last changed 12.12.2022

//containing the molecular topoloy following the GROMACS convention
#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "classes2.h"//"Topology.h"

//use it only in case of molecular RMC with flexible, MD-like molecules
	//The molecules are defined following the GROMACS convention, the GROMACS topology file can 
	//be used with some additional parameters specified at the line ends after the text qualifier ;
	//to bind the GROMACS topology to the RMC configuration
int  Topology::nexclusion_max=0;//the largest exclusion number for all the molecule types
int  Topology::nmoltype=0;//number of molecule types according to the number of [ moleculetype ] directive
int  Topology::nsysmoltype=0;//number of molecule types given after the [ molecules ] directive
int	 Topology::n_gr_atoms=0;//total number of different atoms for gromacs in each different molecule
int	 Topology::n_gr_virtuals = 0;//total number of different virtual sites for gromacs in each different molecule
int	 Topology::nGRtypes=0;//number of different GROMACS atom types (for vdW)
int  Topology::ncharge_types=0;//number of different charge types
int  Topology::ncharge_groups=0;//number of different charge groups
int  Topology::nRMC_charge_centre=0;//number of RMC charge group    
int  Topology::virtual_count=0;//total number of virtual sites
int  Topology::npairs=0;//total number of 1-4 pairs
int  Topology::nbonds=0;//total number of bonds
int  Topology::nangles=0;//total number of angles
int  Topology::ndihedrals=0;//total number of dihedrals
int  Topology::nperdihedrals=0;//total number of GROMACS type (1) periodic dihedrals
int  Topology::nharmdihedrals=0;//total number of GROMACS type (2) harmonic dihedrals 
int  Topology::nRBdihedrals=0;//total number of RBdihedrals GROMACS type (3)
int  Topology::npair_types=0;//number of different 1-4 pair types
int  Topology::nexclusionlines;//number of exclusion lines in topology
int  Topology::nbond_types=0;//number of different bond types
int  Topology::nangle_types=0;//number of different angle types
int  Topology::nperdihedral_types=0;//number of different periodic dihedral types
int  Topology::nharmdihedral_types=0;//number of different harmonic dihedral types
int  Topology::nRBdihedral_types=0;//number of different RBdihedral types
int  Topology::nsegment=0;//number of segments in the preprocessor dir arrays
int  Topology::bond_weight_mode=1;//whether the same weight parameter is used for all the bond types, output is collapsed to one entry
int  Topology::angle_weight_mode=1;//whether the same weight parameter is used for all the angle types, output is collapsed to one entry
int  Topology::perdih_weight_mode=1;//whether the same weight parameter is used for all the perdihedral types, output is collapsed to one entry
int  Topology::harmdih_weight_mode=1;//whether the same weight parameter is used for all the harmdihedral types, output is collapsed to one entry
int  Topology::RBdih_weight_mode=1;//whether the same weight parameter is used for all the RBdihedral types, output is collapsed to one en
int  Topology::last_charge_group=-9999;//the last charge group number
int  Topology::last_directive = -1;//last directive
int  Topology::gen_pairs=0;//whether generate pairs based on the combination rule
int *Topology::nexclusions_per_type=0;//number of exclusions per moleculetype
int *Topology::nifdef_ifndef=0;//current number of ifdef and ifndef preproc directives for each segment
int *Topology::pifdef_ifndef=0;//indicator array of ifdef  and ifndef preprocessor directives present in the topology ind itp file
int *Topology::pelse=0;//indicator array of else preprocessor directives present in the topology ind itp file
int *Topology::nmol_per_type=NULL;
int *Topology::natoms_per_type=NULL;
int *Topology::nvirtuals_per_type=NULL;//number of virtual sites in one molecule /molecule type
int *Topology::ncharge_group_per_type=0;//number of charge groups per type
int *Topology::cumul_atoms=NULL;
int *Topology::cumul_virtuals=NULL;//cumulative number of virtual sites in one molecule /molecule type
int *Topology::first_index=0;//atom index of the given atom in the first molecule of the given residue (1-natoms)
int *Topology::second_index=0;//atom index of the given atom in the first molecule of the given residue (1-natoms)
int *Topology::delta_index=0;//difference between second and first_index
int *Topology::RMCindex_offset=0;//offset for the RMC_index of the molecules, useful if it is not the first in the system
int *Topology::npairs_per_type=0;//number of 1-4 pairs per molecule type
int *Topology::cumul_pairs=0;//cumulative number of 1-4 pairs in one molecule /type
int *Topology::pair_index=0;//indices of the pairs in pairs
int *Topology::pair_type=0;//type of the pair
int *Topology::nexclusionlines_per_type;//number of exclusion lines per molecule type
int *Topology::cumul_exclusionlines;//cumulative number of exclusion lines in a molecule type/molecule type
int *Topology::exclusion_centre;//index of the central atom for the exclusion
int *Topology::nexcluded_neighs;//number of excluded atoms for the central atom 
int *Topology::cumul_nexcluded_neighs;//cumulative number of excluded atoms for the central atom 
int *Topology::exclusions;//indices of the excluded atoms for an atom in the exclusion_centre array;
int *Topology::nbonds_per_type=0;//number of bonds per molecule type
int *Topology::nbondtypes_per_type=0;//number of bond types per molecule type
int *Topology::cumul_bonds=0;//cumulative number of bonds in a molecule type/ molecule type
int *Topology::cumul_bondtypes=0;//cumulative number of bond types in a molecule /molecule type
int *Topology::bond_index=0;//indices of the bonds in pairs
int *Topology::bond_type=0;//type of the bond
int *Topology::zero_bond_sigma=0;//size:nmoltype+1, [0] element is 0: no zero sigma at all, 1: there is zero sigma, then index of the first non-zero sigma interaction for molecule type
int *Topology::nangles_per_type=0;//number of angles per molecule type
int *Topology::nangletypes_per_type=0;//number of angle types per molecule type
int *Topology::cumul_angles=0;//cumulative number of angles in a molecule type/ molecule type
int *Topology::cumul_angletypes=0;//cumulative number of angle types in a molecule /molecule type
int *Topology::angle_index=0;//indices of the angles in triplets (neighbour-central-neighbour)
int *Topology::angle_type=0;//type of the angle
int  Topology::zero_angle_sigma=0;//0: no zero sigma at all, 1: there is zero sigma
int *Topology::ndihedrals_per_type=0;//total number of dihedrals per molecule type
int *Topology::nperdihedrals_per_type=0;//number of periodic dihedrals per molecule type
int *Topology::nperdihtypes_per_type=0;//number of periodic dihedral types per molecule type
int *Topology::perdihedral_type=0;//type of the periodic dihedral
int *Topology::cumul_perdihedrals=0;//cumulative number of periodic dihedrals in a molecule type/ molecule type
int *Topology::cumul_perdihtypes=0;//cumulative number of periodic dihedral types in a molecule /molecule type
int  Topology::zero_perdih_sigma=0;//0: no zero sigma at all, 1: there is zero sigma
int *Topology::nharmdihedrals_per_type=0;//number of harmonic dihedrals per molecule type
int *Topology::nharmdihtypes_per_type=0;//number of harmonic dihedral types per molecule type
int *Topology::harmdihedral_type=0;//type of the harmonic dihedral
int *Topology::cumul_harmdihedrals=0;//cumulative number of harmonic dihedrals in a molecule type/ molecule type
int *Topology::cumul_harmdihtypes=0;//cumulative number of harmonic dihedral types in a molecule /molecule type
int  Topology::zero_harmdih_sigma=0;//0: no zero sigma at all, 1: there is zero sigma
int *Topology::nRBdihedrals_per_type=0;//number of RBdihedrals per molecule type
int *Topology::nRBdihtypes_per_type=0;//number of RB dihedral types per molecule type
int *Topology::cumul_RBdihedrals=0;//cumulative number of RBdihedrals in a molecule type/ molecule type
int *Topology::cumul_RBdihtypes=0;//cumulative number of RB dihedral types in a molecule /molecule type
int *Topology::RBdihedral_type=0;//type of the RBdihedral
int  Topology::zero_RBdih_sigma=0;//0: no zero sigma at all, 1: there is zero sigma
int *Topology::dihedral_index=0;//indices of the dihedrals in quartetts (1-2-3-4), for all the three dihedral type
int *Topology::dihedral_type=0;//type of the dihedral similatly to GROMACS: 
int *Topology::perdihedral_multiplicity=0;//multiplicity for the periodic dihedrals
int	*Topology::GROMACS_type=0;//the GROMACS atom type of the atom for each GROMACS atom
int	*Topology::charge_type=0;//the charge type of the atom for each GROMACS atom
int	*Topology::charge_group=0;//the charge group of the atom for each GROMACS atom
int	*Topology::natoms_per_charge_group=0;//number of atoms for each charge group
double  Topology::bond_r_max=0.0;//maximum bond length
double *Topology::vdW14_pot1=0;//LJ_sigma in A or C(6)=4*epsilon*(sigma^6) in kJ/mol*A^6, or Buckingham A
double *Topology::vdW14_pot2=0;//LJ_epsilon in kJ/mol or C(N)=4*epsilon*(sigma^N) in kJ/mol*A^N or Buckingham B
double *Topology::vdW14_pot1_partial=0;//the vdW14_pot1 rendered to GROMACS type
double *Topology::vdW14_pot2_partial=0;//the vdW14_pot2 rendered to GROMACS type
double *Topology::bond_r=0;//equilibrium value for bonds in Angstrom
double *Topology::bond_k=0;//force constant for bonds (read as kJ/mol/nm2 and converted to kJ/mol/A2)
double *Topology::bond_sigma=0;//sigma values for the bond types
double *Topology::angle_ang=0;//equilibrium value for the angles in degrees 
double *Topology::angle_k=0;//force constant for angles (kJ/mol/rad2)
double *Topology::angle_sigma=0;//sigma values for the angle types
double *Topology::perdihedral_ang=0;//equilibrium value for the (periodic GR type 1 ) dihedral angle in degrees 
double *Topology::perdihedral_k=0;//force constant for (GR type 1 periodic) dihedrals (kJ/mol/rad2)
double *Topology::harmdihedral_ang=0;//equilibrium value for the harmonic, GR 2 type dihedral angle in degrees 
double *Topology::harmdihedral_k=0;//force constant for harmonic, GR 2 type dihedral (kJ/mol/rad2)
double *Topology::RBdihedral_C0=0;//C0 force constant for Ryckaert-Bellemans dihedrals (kJ/mol)
double *Topology::RBdihedral_C1=0;//C1 force constant for Ryckaert-Bellemans dihedrals (kJ/mol)
double *Topology::RBdihedral_C2=0;//C2 force constant for Ryckaert-Bellemans dihedrals (kJ/mol)
double *Topology::RBdihedral_C3=0;//C3 force constant for Ryckaert-Bellemans dihedrals (kJ/mol)
double *Topology::RBdihedral_C4=0;//C4 force constant for Ryckaert-Bellemans dihedrals (kJ/mol)
double *Topology::RBdihedral_C5=0;//C5 force constant for Ryckaert-Bellemans dihedrals (kJ/mol)
double *Topology::perdihedral_sigma=0;//sigma values for the periodic dihedral types
double *Topology::harmdihedral_sigma=0;//sigma values for the harmonic dihedral types
double *Topology::RBdihedral_sigma=0;//sigma values for the RB dihedral types
double *Topology::charge=0;//charge of the atoms/per GROMACS atom types
double *Topology::virtual_params;//the parameters for the virtual sites
bool	Topology::RB_used=false;//indicator, that RB dhihedral is used

char (*Topology::molname)[NAME_SIZE];//name of the molecules
char (*Topology::atom_name)[10];//name of the atoms
char (*Topology::atom_type)[20];//type of the atoms
char *Topology::currentfilename;//pointer to the name of the file currently handled
char *Topology::systemname;//name of the GROMACS system
char *Topology::force_field_name=NULL;//name of the force field

const char *Topology::virtual_types[6] = { "2","3","3fd","3fad","3out", "4fdn" };//the virtual site types
virtual_site_struct *Topology::virtual_site=NULL;//array containing the information about the virtul sites

//Initializing the pointer array to the reading functions for the different line types
int (*Topology::read_function[N_GR_DIR])(ifstream &file) = {Topology::GetSystem,Topology::GetMolecules,Topology::GetMoleculetype,Topology::GetAtoms,Topology::GetBonds, Topology::GetAngles,Topology::GetDihedrals,Topology::GetPairs,Topology::GetExclusions,Topology::GetVirtual2,Topology::GetVirtual3, Topology::GetVirtual4};

T_interaction_type Topology::interaction_type;//the type of the interaction between particles
T_force_field Topology::force_field_type;//enumerated type of the force_field

Topology::Topology()
{

	//cout<<"Topology object is created"<<endl;
};

void Topology::SetParams()
{
	int i;
	
	SetArraysize(&systemname,2*NAME_SIZE,"systemname","Topology::SetParams");

	//creating the static arrays for preprocessor directive processing
	SetArraysize(&nifdef_ifndef,N_TOPFILE,"nifdef_ifndef","Topology::SetParams");
	SetArraysize(&pifdef_ifndef,N_TOPFILE*N_ACT_IFDEF,"pifdef_ifndef","Topology::SetParams");
	SetArraysize(&pelse,N_TOPFILE*N_ACT_IFDEF,"pelse","Topology::SetParams");
	
	for (i=0;i<N_TOPFILE;i++)
		nifdef_ifndef[i] =0;

	for (i=0;i<N_TOPFILE*N_ACT_IFDEF;i++)
	{
		pifdef_ifndef[i]=0;
		pelse[i]=0;
	}
	

	
}
void Topology::GetParams()
{
	int i;
	int nfiles;//the number of currently open topology+itp files
	//int last_option=-1;
	
	int present;
	bool newfile,includeline;
	std::streamoff file_pos;
	char *buffer;
	char *itpfilenames;
	ifstream *file,*files;
	last_directive = -1;


	files =new ifstream[N_TOPFILE];//pointers to file variables
	SafeOpenTextFile(files[0],topfilename);
	if (CheckFileState(files[0],"Topology::GetParams",topfilename)==0)
	{
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}



	SetArraysize(&buffer,LINE_SIZE,"buffer","Topology::GetParams");
	SetArraysize(&itpfilenames,4*FILE_NAME_SIZE,"itpfilenames","Topology::GetParams");
	
	nfiles=1;//topology file is open
	file=files;//set the file pointer to the current file to read
	currentfilename=topfilename;
	//going through the file, searching for itp files and directives
	newfile=false;//no new itp file yet
	
	includeline = true;//no preprocessor directive to prevent line procession
	
	do
	{
		file_pos=(*file).tellg();
		present=FindDirective(*file);
		
		//preprocessor directives has to be processed, even if we are in a section which is not included, as the end of the section has to be found
		if (includeline && present>-1 && present<100)//new GROMACS directive was found
		{
			last_directive=present;

		}
		else
		{
			
			//preprocessor directive is found or line should be processed if includeline indicates it
			//if includeline is false, should check for else or endif
			switch (present - 100)
			{
				case 0://include found
				{
					if (!includeline)
						break;
					if (CheckForcefield(*file, buffer) == 2)
					{
						mystrncpy(itpfilenames, 4 * FILE_NAME_SIZE, buffer + 1, strlen(buffer) - 1);
						itpfilenames[strlen(buffer) - 2] = '\0';
						//It is assumed,that it is an itp file
						SafeOpenTextFile(files[nfiles], &itpfilenames[nfiles - 1]);
						nfiles++;

						if (CheckFileState(file[nfiles - 1], "Topology::GetParams", &itpfilenames[nfiles - 2]) == 0)
						{
							cout << "Cannot run this way, exiting..." << endl;
							CleanExit();
						}
						//go to the next line in the old file
						(*file).seekg(file_pos, ios::beg);
						SkipLine(*file, currentfilename);
						newfile = true;//do not skip a line in the new file
						file = &files[nfiles - 1];//this will be the current file to read
						currentfilename = &itpfilenames[nfiles - 2];
						cout << "\nReading parameters from the include topology file " << currentfilename << endl;


					}
					break;

				}
				case 1://ifdef found
				{
					
					if (CheckPreprocDir(*file, nfiles, 1) > 0)//the given token was passed to the program from the *.dat, see, if it is in an active section, if there are nestled directives 
					{
						includeline = IncludeLine(nfiles);//this is the actual state considering all the nestled preprocessor directives
						if (::debug)
						{
							if (includeline)
								cout << "\tcompiler option was found for ifdef preprocessor directive, proceeding lines will be included in file " << currentfilename << endl;
							else
								cout << "\tcompiler option was found for ifdef preprocessor directive, but proceeding lines will not be included due to nestling in file processed " << currentfilename << endl;
						}
					}
					else
					{
						//the token was not passed it is inactive regadless the nestled directives
						includeline = false;
						if (::debug)
							cout << "\tinactive ifdef preprocessor directive found, proceeding lines will not be included from file " << currentfilename << endl;
					}
					break;
				}
				case 2://else found
				{
					pelse[(nfiles - 1) * N_ACT_IFDEF + nifdef_ifndef[nfiles - 1] - 1] = 1;//found the else for the latest ifdef
					includeline = IncludeLine(nfiles);//this is the actual state considering all the nestled preprocessor directives
					if (::debug)
					{
						if (includeline)
							cout << "\tlines will be included from the else section of the latest proprocessor directive from file " << currentfilename << endl;
						else
							cout << "\tlines will not be included from the else section of the latest proprocessor directive from file " << currentfilename << endl;
					}
					break;
				}
				case 3://endif found
				{
					//end of the latest ifdef, clear the array elements, and decrease the number of active ifdefs
					pifdef_ifndef[(nfiles - 1) * N_ACT_IFDEF + nifdef_ifndef[nfiles - 1] - 1] = 0;
					pelse[(nfiles - 1) * N_ACT_IFDEF + nifdef_ifndef[nfiles - 1] - 1] = 0;
					nifdef_ifndef[nfiles - 1]--;
					includeline = IncludeLine(nfiles);
					if (::debug)
					{
						cout << "End of the latest preprocessor directive from file " << currentfilename << endl;
						if (includeline)
							cout << "\tlines will be included from file " << currentfilename << endl;
						else
							cout << "\tlines will not be included from file " << currentfilename << endl;
					}
					break;
				}
				case 4://ifndef found
				{
					if (!includeline)
						break;
					if (CheckPreprocDir(*file, nfiles, 2) > 0)
					{
						includeline = false;//this is the actual state based on this, no need to consider all the nestled preprocessor directives
						if (::debug)
							cout << "\tcompiler option was found for ifndef preprocessor directive, proceeding lines will not be included from file " << currentfilename << endl;
							
					}
					else//the token was not given,so this would be an active section, see according to the nestled directives
					{
						includeline = IncludeLine(nfiles);
						if (::debug)
						{
							if (includeline)
								cout << "\tcompiler option was not found for preprocessor directive ifndef, proceeding lines will be included from file " << currentfilename << endl;
							else
								cout << "\tcompiler option was not found for preprocessor directive ifndef, but proceeding lines will not be included due to nestling from file " << currentfilename << endl;
						}
					}
					break;
				}
				case 100://only comment or empty line, no further processing
				{
					break;
				}
				case -102://unimplemented GROMACS directive
				{
					if (!includeline)
						break;
					last_directive = -2;
					break;
				}
				default:
				{
					if (last_directive != -2 && includeline)//see what is in the line, (last_directive=-2 is an unimplemented GROMACS directive, where nothing should be done
					{
						//this line should be processed according to the last directive,if preprocesssor directive allows it
						//includeline is set, when a new preprocessor directive is found
						(*read_function[last_directive])(*file);
					}//end if last_directive==-2
				}
			}//end switch
			
		}

		if (newfile)
			newfile=false;
		else if (!(*file).rdstate()) //If earlier something error happened skip this part (this condition eof bit is kept)
		{
			(*file).seekg(file_pos,ios::beg);
			SkipLine(*file,currentfilename);
		}

		if ((*file).rdstate() & (*file).eofbit)//end of this file or error during at the reading, see, if there is a previos file open
		{
			(*file).clear(ios::goodbit);
			(*file).close();
			Reset(nfiles);//resetting the array element in the pifdef and pelse arrays 
			nfiles--;
			if (nfiles>0)//this was an itp file, close it, and continue with the previous file
			{
				file=&files[nfiles-1];//set the pointer to the previous file
				currentfilename=(nfiles==1 ? topfilename : &itpfilenames[nfiles-2]);
				includeline=IncludeLine(nfiles);
				
			}
		}	
		//cout<<nfiles<<endl;
	} while (nfiles>0);
	
	for (i = 0; i < nmoltype; i++)
	{
		if (nmol_per_type[i] == -1)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The number of molecule types based on the different [ moleculetype ] directive";
			if (nmoltype == 1)
				cout << " is ";
			else
				cout << "s are ";
			cout << nmoltype << endl;
			for (int imol = 0; imol < nmoltype; imol++)
			{
				cout << "\t" << imol + 1 << ". " << molname[imol] << "\t";
				if (nmol_per_type[imol] == -1)
					cout << "not given" << endl;
				else
					cout << nmol_per_type[imol] << endl;
			}

			cout << "Not all molecule types defined with  [ moleculetype ] directive was given after the [ molecules ] directive" << endl;
			cout<<"specifying the number of molecules!" << endl;
			cout << "Correct the topology and try again, exiting..." << endl;
			CleanExit();
		}
	}
	//check the number of exclusions, and set it to zero, if no molecules defined
	for (i=0;i<nmoltype;i++)
	{
		if (nbonds==0 && nangles==0 && ndihedrals==0 && nexclusions_per_type[i]!=0)//no molecules defined
			nexclusions_per_type[i]=0;
		if (nexclusions_per_type[i]>nexclusion_max)
			nexclusion_max=nexclusions_per_type[i];
	}
	int leading_sigma=-1;
	//Checking, whether there is non-zero sigma value for each molecule types interaction types, where zero sigma interaction is also given
	if (nbonds>0)
		CheckSigma(&bond_weight_mode, cumul_bondtypes, nbondtypes_per_type, zero_bond_sigma, bond_sigma, "bond",&leading_sigma);
	if (nangles>0)
		CheckSigma(&angle_weight_mode, cumul_angletypes, nangletypes_per_type, &zero_angle_sigma, angle_sigma, "angle", &leading_sigma);
	if (nperdihedrals>0)
		CheckSigma(&perdih_weight_mode, cumul_perdihtypes, nperdihtypes_per_type, &zero_perdih_sigma, perdihedral_sigma, "periodic dihedral", &leading_sigma);
	if (nharmdihedrals>0)
		CheckSigma(&harmdih_weight_mode, cumul_harmdihtypes, nharmdihtypes_per_type, &zero_harmdih_sigma, harmdihedral_sigma, "harmonic dihedral", &leading_sigma);
	if (nRBdihedrals>0)
		CheckSigma(&RBdih_weight_mode, cumul_RBdihtypes, nRBdihtypes_per_type,  &zero_RBdih_sigma, RBdihedral_sigma, "RB dihedral", &leading_sigma);
	
	
	//check  forcefiled name and type and set it, if necessary
	if (force_field_name == NULL)
	{
		cout << "\n*****ERROR*****"<<endl;
		cout<<"No valid force field name was given in the topology file!" << endl;
		cout << "Presently RMC can only run with the following force fields:" << endl;
		for (i = 0; i <= N_FORCEF-1; i++)
			cout << "\t" << force_field[i] << endl;
		cout << "Choose the appropriate one and try again!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();


	}

	if (::debug)
		cout<<"Finished reading the topology parameters"<<endl;

};

int Topology::FindDirective(ifstream &file )
{
	int pos=-1;
	int option;
	
	longint first_pos;//,last_pos;
	longint offset;
	char *directive;
	
	SetArraysize(&directive,20,"directive","Topology::FindDirective");


	
	first_pos=(longint)file.tellg();
	file.getline(buf,LINE_SIZE);

	if (::debug)
		cout << buf << endl;

	if(!FindEmtyLine(buf) && strlen(buf)>0)
	{

		option=SearchCOD(buf,offset);
		switch (option)	//there is probably a preprocessor or GROMACS directive, find it
		{
			case 1://preprocessor directive
			{
				pos=GetDirective(file, first_pos+ offset,N_COMP_DIR,comp_dir);
				pos+=100;//offset the preprocessor directive
				break;
			}
			case 2://GROMACS directive
			{
				//gives the index of the directive in the array or -2 if unimplemented
				pos=GetDirective(file, first_pos+offset,N_GR_DIR,directives);
				if (::debug)
					if (pos >= 0)
						cout << "GROMACS directive "<<directives[pos] << " was found" << endl;
				break;
			}
			case 3://only comment
			{
				pos=200;//no further processing needed
				break;
			}
			default:
			{
				//line possibly has to be processed, put the pointer back to the beginning of the line
				file.seekg(first_pos,ios::beg);//to set the file position after the line reading, as might be reset					
			}
		}
	}
	else
		pos=200;//no further processing needed


	return pos;
}

//determining, whether the line has a preprocessor or GROMACS directive in it
int Topology::SearchCOD(char *buf, longint &offset)
{
	int i,beg_dir;
	int offset1,offset2,offset3;
	int found=0;
	char *pchar;

	pchar=buf;

	//There can be more, than one control character on a line
	beg_dir='#';
	pchar=strchr(buf,beg_dir);
	offset1=(int)(pchar-buf)+1;
	found=0;
	if (pchar!=NULL && offset1<=(int)strlen(buf))
	{
		found+=1;//found compiler option
		offset=offset1;
	}
	
	beg_dir='[';
	pchar=strchr(buf,beg_dir);
	offset2=int(pchar-buf)+1;
	if (pchar!=NULL && offset2<=(int)strlen(buf))
	{
		found+=2;//found GROMACS directive
		offset=offset2;
	}
	
	beg_dir=';';
	pchar=strchr(buf,beg_dir);
	offset3=int(pchar-buf)+1;
	if (pchar!=NULL && offset3<=(int)strlen(buf))
	{
		found+=4;//found text qualifier
		for (i=0;i<offset3-1;i++)
		{
			if (!(buf[i]==' ' || buf[i]=='\t'))
			{
				found+=4;//there is something else before it, line should be processed
				break;
			}
		}
	}
	
	//see, which control character has precedence
	//# or [ cannot legally be in the same line, if it is, the first one will be processed
	//if the text qualifier is before # or [, than they will be ignored, as they are commented out
	switch (found)
	{
		case 3:
		{
			//both # [ found, which one is the first
			if (offset1<offset2)
			{
				found=1;//#
				offset=offset1;
			}
			else
			{
				offset=offset2;
				found=2;//[
			}
			break;
		}
		case 4:
		{
			//; was found at the beginnig of the line, comment
			found=3;
			offset=offset3;
			break;
		}
		case 5:
		{
			//both # ; found, ; is the first
		}
		case 6:
		{
			//both [ ; found, ; is the first
		}
		case 7:
		{
			//; # and [ ; is the first
			found=3;//; is the first
			offset=offset3;
			break;
		}
		case 8:
		{
			found=4;//some data is before the text qualifier
			offset=offset3;
			break;
		}
		case 9:
		{
			//both ; and # found, there is something before ;, but not necessarily #
			if (offset1<offset3)
			{
				found=1;//#
				offset=offset1;
			}
			else
			{
				found=4;//there is data to process before ;
				offset=offset3;
			}
			break;
		}
		case 10:
		{
			//both ; and [ found, there is something before ;, but not necessarily [
			if (offset2<offset3)
			{
				found=2;//[
				offset=offset2;
			}
			else
			{
				found=4;//there is data to process before ;
				offset=offset3;
			}
			break;
		}
		case 11:
		{
			//; # and [ found, there is something before ;, but not necessarily # and or [
			if (offset1<offset3)
			{
				offset=offset2;
				if (offset1<offset2)
				{
					found=1;//#
					offset=offset1;
				}
				else
					found=2;//[
			}
			else
			{
				if (offset2<offset3)
					found=2;//[
				else
				{
					found=4;//there is data to process before ;
					offset=offset3;
				}
			}
			break;
		}

	}
	return found;
};

//determine, whether the line is empty (only space or tab)
int Topology::FindEmtyLine(char *buf)
{
	int i;
	for (i=0;i<(int)strlen(buf);i++)
	{
		if (!(buf[i]==' ' || buf[i]=='\t' || buf[i]=='\r'))
			return 0;
	}

	return 1;
};

//determine the compiler option or directive
int Topology::GetDirective(ifstream &file, longint file_pos, int n_dir, const char **dir)
{	
	int i,terminator,offset;
	int pos=-1;
	char *directive,*pchar;

	SetArraysize(&directive,20,"directive","Topology::GetDirective");
	file.seekg(file_pos,ios::beg);
	file>>directive;
	//see, if by mistake there is no space between the directive text and the ending bracket
	terminator = ']';
	pchar = strchr(directive, terminator);
	offset = (int)(pchar - directive);
	if (pchar!=NULL && offset < (int)strlen(directive))
		directive[ offset] = '\0';
	
	for (i=0;i<n_dir;i++)
	{
		if (strcmp(directive,dir[i])==0)
		{

			pos=i;//the directive was found
			break;
		}
	}
	if (pos == -1)//the directive or the compiler option is not on the list
	{
		
		if (dir == comp_dir)//compiler option
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The #" << directive << " preprocessor directive is not implemented! in line:" << endl;
			cout << buf << endl;
			cout<<"Correct the topology file " << currentfilename << " and try again. Exiting..." << endl;
			CleanExit();
		}
		else//this is a GROMACS directive not implemented
		{
			cout << "\nWARNING(" << ++warn << "): The [ " << directive << " ] GROMACS directive is not implemented in line:" << endl;
			cout <<"\t"<< buf << endl;
			cout<<"\tThe lines belonging to it will not be processed in file "<<currentfilename<<"!" << endl;
			cout << "\tMake sure, that the topology will be correct next time anyhow! Trying to continue...." << endl;
			pos = -2;//not to process
		}
	}
	return pos;
}

int Topology::GetSystem(ifstream &file)
{
	file>>systemname;
	
	return 1;
};

int Topology::GetMolecules(ifstream &file)
{
	int i,found=0;
	char *temp;
	SetArraysize(&temp,NAME_SIZE,"temp","Topology::GetMolecules");

	if (nsysmoltype==0)
		nsysmoltype++;
	else
		nsysmoltype++;

	if (nsysmoltype>nmoltype)
	{
		cout << "\n****ERROR****"<<endl;
		cout << "Number of molecule types given after the [ molecules ] directive in line:" << endl;
		cout<<buf<<endl;
		cout << "is larger than the number of molecule types defined with the [ moleculetype ] directive!" << endl;
		cout<<"Correct the topology, exiting..."<<endl;
		CleanExit();
	}
	file>>temp;
	for (i=0;i<nmoltype;i++)
	{
		if (strcmp(molname[i],temp)==0)
		{
			found=1;
			file>>nmol_per_type[i];
			break;
		}
		
		
	}
	if (found==0)
	{
		cout << "\n****ERROR****"<<endl;
		cout<<"No definition can be found for the molecule "<<temp<<" given after the [ molecules ] directive in line:"<<endl;
		cout << buf << endl;
		cout<<"Check the topology file "<<currentfilename<<" maybe the include topology file is missing from the topology!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (nsysmoltype==nmoltype)
	{
		//set the number of charge group centres
		for (i=0;i<nmoltype;i++)
			nRMC_charge_centre+=nmol_per_type[i]*ncharge_group_per_type[i];
		
	}

	return found;
};
int Topology::GetMoleculetype(ifstream &file)
{
	//Create the array or resize it
	if (nmoltype==0)
	{
		nmoltype++;
		SetArraysize(&nmol_per_type,nmoltype,"nmol_per_type","Topology::GetMoleculetype");
		SetArraysize(&natoms_per_type,nmoltype,"natoms_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_atoms,nmoltype+1,"cumul_atoms","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nvirtuals_per_type, nmoltype, "nvirtuals_per_type", "Topology::GetMoleculetype");
		SetArraysize(&cumul_virtuals, nmoltype + 1, "cumul_virtuals", "Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&npairs_per_type,nmoltype,"npairs_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_pairs,nmoltype+1,"cumul_pairs","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nexclusionlines_per_type, nmoltype, "nexclusionlines_per_type", "Topology::GetMoleculetype");
		SetArraysize(&cumul_exclusionlines, nmoltype + 1, "cumul_exclusionlines", "Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nbonds_per_type,nmoltype,"nbonds_per_type","Topology::GetMoleculetype");
		SetArraysize(&nbondtypes_per_type,nmoltype,"nbondtypes_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_bonds,nmoltype+1,"cumul_bonds","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&cumul_bondtypes,nmoltype+1,"cumul_bondtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nangles_per_type,nmoltype,"nangles_per_type","Topology::GetMoleculetype");
		SetArraysize(&nangletypes_per_type,nmoltype,"nangletypes_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_angles,nmoltype+1,"cumul_angles","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&cumul_angletypes,nmoltype+1,"cumul_angletypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&ndihedrals_per_type,nmoltype,"ndihedrals_per_type","Topology::GetMoleculetype");
		SetArraysize(&nperdihedrals_per_type,nmoltype,"nperdihedrals_per_type","Topology::GetMoleculetype");
		SetArraysize(&nperdihtypes_per_type,nmoltype,"nperdihtypes_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_perdihedrals,nmoltype+1,"cumul_perdihedrals","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&cumul_perdihtypes,nmoltype+1,"cumul_perdihtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nharmdihedrals_per_type,nmoltype,"nharmdihedrals_per_type","Topology::GetMoleculetype");
		SetArraysize(&nharmdihtypes_per_type,nmoltype,"nharmdihtypes_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_harmdihedrals,nmoltype+1,"cumul_harmdihedrals","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&cumul_harmdihtypes,nmoltype+1,"cumul_harmdihtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nRBdihedrals_per_type,nmoltype,"nRBdihedrals_per_type","Topology::GetMoleculetype");
		SetArraysize(&nRBdihtypes_per_type,nmoltype,"nRBdihtypes_per_type","Topology::GetMoleculetype");
		SetArraysize(&cumul_RBdihedrals,nmoltype+1,"cumul_RBdihedrals","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&cumul_RBdihtypes,nmoltype+1,"cumul_RBdihtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		SetArraysize(&nexclusions_per_type,nmoltype,"nexclusion","Topology::GetMoleculetype");
		SetArraysize(&RMCindex_offset,nmoltype,"RMCindex_offset","Topology::GetMoleculetype");
		SetArraysize(&ncharge_group_per_type,nmoltype,"ncharge_group_per_type","Topology::GetMoleculetype");
		SetArraysize(&zero_bond_sigma,nmoltype+1,"zero_bond_sigma","Topology::GetMoleculetype");//it will be nmoltype +1 long
		molname=new char[nmoltype][NAME_SIZE];
		if (molname==NULL)
			NoArray("Topology::GetMoleculetype","molname");
		cumul_atoms[0]=0;
		cumul_virtuals[0] = 0;
		cumul_pairs[0]=0;
		cumul_exclusionlines[0] = 0;
		cumul_bonds[0]=0;
		cumul_angles[0]=0;
		cumul_perdihedrals[0]=0;
		cumul_harmdihedrals[0]=0;
		cumul_RBdihedrals[0]=0;
		cumul_bondtypes[0]=0;
		cumul_angletypes[0]=0;
		cumul_perdihtypes[0]=0;
		cumul_harmdihtypes[0]=0;
		cumul_RBdihtypes[0]=0;
		nmol_per_type[0] = -1;//not initialized
		zero_bond_sigma[0]=0;
	}
	else
	{
		ResizeArray(&nmoltype,nmoltype+1,&nmol_per_type,"nmol_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		nmol_per_type[nmoltype] = -1;//not initialized
		ResizeArray(&nmoltype,nmoltype+1,&natoms_per_type,"natoms_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_atoms,"cumul_atoms","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype, nmoltype + 1, &nvirtuals_per_type, "nvirtual_per_type", "Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype, nmoltype + 1, &cumul_virtuals, "cumul_virtuals", "Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype -= 2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,(char**)&molname,"molname","Topology::GetMoleculetype",NAME_SIZE);
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&npairs_per_type,"npairs_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_pairs,"cumul_pairs","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype, nmoltype + 1, &nexclusionlines_per_type, "nexclusionlines_per_type", "Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype, nmoltype + 1, &cumul_exclusionlines, "cumul_exclusionlines", "Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype -= 2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nbonds_per_type,"nbonds_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nbondtypes_per_type,"nbondtypes_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_bonds,"cumul_bonds","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&cumul_bondtypes,"cumul_bondtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nangles_per_type,"nangles_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nangletypes_per_type,"nangletypes_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_angles,"cumul_angles","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&cumul_angletypes,"cumul_angletypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&ndihedrals_per_type,"ndihedrals_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nperdihedrals_per_type,"nperdihedrals_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nperdihtypes_per_type,"nperdihtypes_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_perdihedrals,"cumul_perdihedrals","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&cumul_perdihtypes,"cumul_perdihtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nharmdihedrals_per_type,"nharmdihedrals_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nharmdihtypes_per_type,"nharmdihtypes_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_harmdihedrals,"cumul_harmdihedrals","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&cumul_harmdihtypes,"cumul_harmdihtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nRBdihedrals_per_type,"nRBdihedrals_per_type","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&nRBdihtypes_per_type,"nRBdihtypes_per_type","Topology::GetMoleculetype");
		//no setting back to the original, as this array is always 1 larger,than the others
		ResizeArray(&nmoltype,nmoltype+1,&cumul_RBdihedrals,"cumul_RBdihedrals","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&cumul_RBdihtypes,"cumul_RBrdihtypes","Topology::GetMoleculetype");//it will be nmoltype +1 long
		nmoltype-=2;//setting back to the new size, as cumul is one longer, than needed
		ResizeArray(&nmoltype,nmoltype+1,&nexclusions_per_type,"nexclusion","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&RMCindex_offset,"RMCindex_offset","Topology::GetMoleculetype");
		nmoltype--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nmoltype,nmoltype+1,&ncharge_group_per_type,"ncharge_group_per_type","Topology::GetMoleculetype");
		ResizeArray(&nmoltype,nmoltype+1,&zero_bond_sigma,"zero_bond_sigma","Topology::GetMoleculetype");//this has to be nmoltypes+1
		nmoltype--;
				
	}
	//setting the starting values for the arrays, where the number will be increased, if new entry is found
	natoms_per_type[nmoltype-1]=0;
	cumul_atoms[nmoltype] = cumul_atoms[nmoltype - 1];
	nvirtuals_per_type[nmoltype-1] = 0;
	cumul_virtuals[nmoltype]=cumul_virtuals[nmoltype-1];
	npairs_per_type[nmoltype-1]=0;
	cumul_pairs[nmoltype]=cumul_pairs[nmoltype-1];
	nexclusionlines_per_type[nmoltype - 1] = 0;
	cumul_exclusionlines[nmoltype] = cumul_exclusionlines[nmoltype - 1];
	nbonds_per_type[nmoltype-1]=0;
	nbondtypes_per_type[nmoltype-1]=0;
	cumul_bonds[nmoltype]=cumul_bonds[nmoltype-1];
	cumul_bondtypes[nmoltype]=cumul_bondtypes[nmoltype-1];
	nangles_per_type[nmoltype-1]=0;
	nangletypes_per_type[nmoltype-1]=0;
	cumul_angles[nmoltype]=cumul_angles[nmoltype-1];
	cumul_angletypes[nmoltype]=cumul_angletypes[nmoltype-1];
	ndihedrals_per_type[nmoltype-1]=0;
	nperdihedrals_per_type[nmoltype-1]=0;
	nperdihtypes_per_type[nmoltype-1]=0;
	cumul_perdihedrals[nmoltype]=cumul_perdihedrals[nmoltype-1];
	cumul_perdihtypes[nmoltype]=cumul_perdihtypes[nmoltype-1];
	nharmdihedrals_per_type[nmoltype-1]=0;
	nharmdihtypes_per_type[nmoltype-1]=0;
	cumul_harmdihedrals[nmoltype]=cumul_harmdihedrals[nmoltype-1];
	cumul_harmdihtypes[nmoltype]=cumul_harmdihtypes[nmoltype-1];
	nRBdihedrals_per_type[nmoltype-1]=0;
	nRBdihtypes_per_type[nmoltype-1]=0;
	cumul_RBdihedrals[nmoltype]=cumul_RBdihedrals[nmoltype-1];
	cumul_RBdihtypes[nmoltype]=cumul_RBdihtypes[nmoltype-1];
	ncharge_group_per_type[nmoltype-1]=0;

	file>>molname[nmoltype-1];//reading the molecule name
	nexclusions_per_type[nmoltype-1]=ReadThisLine(file,3,3,"nexclusion","Topology::GetMoleculetype");
	ExtractSemiColon(file,0);//extract, if there is a semicolon
	RMCindex_offset[nmoltype-1]=ReadThisLine(file,0,0,"RMCindex_offset","Topology::GetMoleculetype");
	last_charge_group = -9999;//to reset the last charge group to start afresh
	return 1;
};

//it will read the atoms
int Topology::GetAtoms(ifstream &file)
{
	int i;
	int charge_group_number;
	int new_type;//new GROMACS type
	longint offset;
	std::streamoff first_pos;//,last_pos;
	double temp_charge;
	char *buf;
	SetArraysize(&buf,LINE_SIZE,"buf","Topology::GetAtoms");

	//Create the array or resize it
	if (n_gr_atoms==0)
	{
		n_gr_atoms++;		
		nGRtypes++;
		ncharge_groups++;
		ncharge_types++;//number of different charge types, can be larger, than nGRtypes
		SetArraysize(&first_index,n_gr_atoms,"first_index","Topology::GetAtoms");
		SetArraysize(&second_index,n_gr_atoms,"second_index","Topology::GetAtoms");
		SetArraysize(&delta_index,n_gr_atoms,"delta_index","Topology::GetAtoms");
		SetArraysize(&GROMACS_type,n_gr_atoms,"GROMACS_type","Topology::GetAtoms");
		SetArraysize(&charge_type,n_gr_atoms,"charge_type","Topology::GetAtoms");
		SetArraysize(&charge_group,n_gr_atoms,"charge_group","Topology::GetAtoms");
		SetArraysize(&atom_name,n_gr_atoms*10,"atom_name","Topology::GetAtoms");
		SetArraysize(&atom_type,nGRtypes*10,"atom_type","Topology::GetAtoms");
		SetArraysize(&charge,ncharge_types,"charge","Topology::GetAtoms");
		SetArraysize(&natoms_per_charge_group,ncharge_groups,"natoms_per_charge_group","Topology::GetAtoms");
		
		
	}
	else
	{
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,&first_index,"first_index","Topology::GetAtoms");
		n_gr_atoms--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,&second_index,"second_index","Topology::GetAtoms");
		n_gr_atoms--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,&delta_index,"delta_index","Topology::GetAtoms");
		n_gr_atoms--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,&GROMACS_type,"GROMACS_type","Topology::GetAtoms");
		n_gr_atoms--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,&charge_type,"charge_type","Topology::GetAtoms");
		n_gr_atoms--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,&charge_group,"charge_group","Topology::GetAtoms");
		n_gr_atoms--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&n_gr_atoms,n_gr_atoms+1,(char**)&atom_name,"atom_name","Topology::GetAtoms",10);
	}

	file >> i;//this is the atom serial number, should be consecutive and starting with 1, it will be checked
	if (cumul_atoms[nmoltype-1]+i != n_gr_atoms)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The serial index for the " << n_gr_atoms<< ". atom in the [ atoms ] section of the " << currentfilename << " file is " <<i<<","<< endl;
		cout << "and it should be " << n_gr_atoms << "! The atoms in a molecule should be numbered consecutuvely and sterting with 1!" << endl;
		cout << "Correct your topology file and do not forget to adjust the indices in the other sections  (bonds, angles...) of the topology as well!" << endl;
		cout << "Cannot run this way, exiting!" << endl;
		CleanExit();
	}
	file>>buf;//the 1. will be the type of the atom
	if (n_gr_atoms==1)//this is the first atom
		mystrcpy(atom_type[nGRtypes-1],20,buf);
	//Check if it is a new type
	new_type=-1;//new type by default
	for (i=0;i<nGRtypes;i++)
	{
		if (strcmp(buf,atom_type[i])==0)
		{
			//not a new type
			new_type=i;//GROMACS type index
			break;
		}
	}
	if (new_type==-1)//this is a new type
	{
		ResizeArray(&nGRtypes,nGRtypes+1,(char**)&atom_type,"atom_type","Topology::GetAtoms",20);
		mystrcpy(atom_type[nGRtypes-1],20,buf);//store the new type in the array
		//ResizeArray(&ncharge_types,ncharge_types+1,&charge,"charge","Topology::GetAtoms");
		new_type=nGRtypes-1;

	}

	GROMACS_type[n_gr_atoms-1]=new_type;//storing the GROMACS atom type for this atom, will be needed for LJ

	for (i=0;i<3;i++)
		file>>atom_name[n_gr_atoms-1];//the 2. will be the name of the atom
	
	//Not the absolute charge group number counts, only consecutive atoms in the topology having the same charge group number 
	//will belong to the same charge group, the charge groups will be indexed consecutively starting with 0
	charge_group_number=ReadThisLine(file,3,1,"charge group","Topology::GetAtoms",currentfilename);
	if (n_gr_atoms==1)//this is the first atom, this will have the first charge group
	{
		charge_group[n_gr_atoms-1]=0;
		natoms_per_charge_group[0]=1;
		ncharge_group_per_type[0]=1;
		last_charge_group=charge_group_number;//to preserve it to compare with the next atoms'
	}
	else
	{
		//see, if it is a new charge group
		if (charge_group_number==last_charge_group)
		{
			//not a new charge group
			natoms_per_charge_group[ncharge_groups-1]++;
			charge_group[n_gr_atoms-1]=charge_group[n_gr_atoms-2];//same as the previous
		}
		else
		{
			//new charge group
			ResizeArray(&ncharge_groups,ncharge_groups+1,&natoms_per_charge_group,"natoms_per_charge_group","Topology::GetAtoms");// this will increase the ncharge_groups
			natoms_per_charge_group[ncharge_groups-1]=1;
			ncharge_group_per_type[nmoltype-1]++;
			charge_group[n_gr_atoms-1]=charge_group[n_gr_atoms-2]+1;//next index
			last_charge_group=charge_group_number;//to preserve it to compare with the next atoms'
		}			
			
	}

	//read the charge
	temp_charge=ReadThisLine(file,3,1.0,"charge","Topology::GetAtoms",currentfilename);
	if (n_gr_atoms==1)//this is the first atom
		charge[ncharge_types-1]=temp_charge;
	//Check if it is a new charge type
	new_type=-1;//new type by default
	for (i=0;i<ncharge_types;i++)
	{
		if (fabs(temp_charge-charge[i])<TOLERANCE)
		{
			//not a new type
			new_type=i;//charge type index
			break;
		}
	}
	if (new_type==-1)//this is a new type
	{
		ResizeArray(&ncharge_types,ncharge_types+1,&charge,"charge","Topology::GetAtoms");
		charge[ncharge_types-1]=temp_charge;
		new_type=ncharge_types-1;

	}
	charge_type[n_gr_atoms-1]=new_type;

	first_pos=file.tellg();
	file.getline(buf,LINE_SIZE);
	//last_pos=file.tellg();
	if (!(SearchCOD(buf,offset)==3 || SearchCOD(buf,offset)==4))//there is mass or nothing before the text qualifier
	{
		cout << "\n*****ERROR*****"<<endl;
		cout<<"This is not a valid atom definition line:"<<endl;
		cout << buf << endl;
		CleanExit();
	}
	file.seekg(first_pos+(longint)offset,ios::beg);
	first_index[n_gr_atoms-1]=ReadThisLine(file,3,1,"first_index","Topology::GetAtoms",currentfilename);//index start with 1
	second_index[n_gr_atoms-1]=ReadThisLine(file,1,1,"second_index","Topology::GetAtoms",currentfilename);
			
	if (first_index[n_gr_atoms - 1] <0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The RMC index of the " << n_gr_atoms - 1 << ". atom in the first molecule of this type (" << first_index[n_gr_atoms - 1] << ")" << endl;
		cout << "is smaller than 0!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	if (second_index[n_gr_atoms - 1] < 0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The RMC index of the " << n_gr_atoms - 1 << ". atom in the second molecule of this type (" << second_index[n_gr_atoms - 1] << ")" << endl;
		cout << "is smaller than 0!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	if (first_index[n_gr_atoms-1]>SimpleCfg::ntotal)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"The RMC index of the "<<n_gr_atoms-1<<". atom in the first molecule of this type ("<<first_index[n_gr_atoms-1]<<")"<<endl;
		cout<<"is larger, than the total number of atoms ("<<SimpleCfg::ntotal<<")!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (second_index[n_gr_atoms-1]>SimpleCfg::ntotal)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"The RMC index of the "<<n_gr_atoms-1<<". atom in the second molecule of this type ("<<second_index[n_gr_atoms-1]<<")"<<endl;
		cout<<"is larger, than the total number of atoms ("<<SimpleCfg::ntotal<<")!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();

	}

	if (first_index[n_gr_atoms - 1] == 0)//virtual site
	{
		delta_index[n_gr_atoms - 1] = 1;
		second_index[n_gr_atoms - 1] = 0;
	}
	else
		delta_index[n_gr_atoms-1]=second_index[n_gr_atoms-1]-first_index[n_gr_atoms-1];
	natoms_per_type[nmoltype-1]++;//increasing the number of atoms for the moleculetype
	cumul_atoms[nmoltype]++;
	return 1;
};

int Topology::GetBonds(ifstream &file)
{
	if (RunParams::fnc!=4)
		return(1);//do not read, if flexible molecules are not needed
	int i,start,new_type;
	double temp_r,temp_k;//for reading the force constant and the distance
	//Create the array or resize it
	if (nbonds==0)
	{
		nbonds++;	
		nbond_types++;
		SetArraysize(&bond_index,nbonds*2,"bond_index","Topology::GetBonds");
		SetArraysize(&bond_type,nbonds,"bond_type","Topology::GetBonds");
		SetArraysize(&bond_r,nbond_types,"bond_r","Topology::GetBonds");
		SetArraysize(&bond_k,nbond_types,"bond_k","Topology::GetBonds");
		SetArraysize(&bond_sigma,nbond_types,"bond_sigma","Topology::GetBonds");
	}
	else
	{
		ResizeArray(&nbonds,nbonds+1,&bond_index,"bond_index","Topology::GetBonds",2);
		nbonds--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nbonds,nbonds+1,&bond_type,"bond_type","Topology::GetBonds");
	}
	//indices will start with 1!
	bond_index[(nbonds-1)*2]=ReadThisLine(file,3,1,"bond_index 1","Topology::GetBonds");//index of the first atom
	bond_index[(nbonds-1)*2+1]=ReadThisLine(file,3,1,"bond_index 2","Topology::GetBonds");//index of the second atom
	if (bond_index[(nbonds - 1) * 2] == bond_index[(nbonds - 1) * 2 + 1])
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The " << nbonds << ". bond index the " << nmoltype << ". molecule type contains the same atoms twice (index " << bond_index[(nbonds - 1) * 2]<<")!" << endl;
		cout << "This does not make sense, check the topology file " << currentfilename << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();

	}
	i=ReadThisLine(file,3,1,"bond force function type","Topology::GetBonds");//the force type
	nbonds_per_type[nmoltype-1]++;
	cumul_bonds[nmoltype]++;
	//some cheking
	if (bond_index[(nbonds-1)*2]>natoms_per_type[nmoltype-1])
	{
		cout << "\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The first atom index in the "<<nbonds_per_type[nmoltype-1]<<" bond of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<bond_index[(nbonds-1)*2]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (bond_index[(nbonds-1)*2+1]>natoms_per_type[nmoltype-1])
	{
		cout << "\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The second atom index in the "<<nbonds_per_type[nmoltype-1]<<" bond of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<bond_index[(nbonds-1)*2+1]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (i!=1)//the topology file does not use harmonic bond potential
	{
		cout << "\nWARNING(" << ++warn << "): The bond type (" << i << ") for the " << nbonds_per_type[nmoltype - 1] << " bond of the " << nmoltype << endl;
		cout<<"\tmolecule type is not harmonic (1) in the " << currentfilename << " file!" << endl;
		cout<<"\tThe force constant after the text qualifier will be assumed to be for harmonic potential in kJ/mol/nm2!"<<endl;
		ExtractSemiColon(file);
		temp_r=ReadThisLine(file,3,1.0,"bond_r","Topology::GetBonds");//equilibruium distance in nm
		temp_k=ReadThisLine(file,3,1.0,"bond_k","Topology::GetBonds");//force constant in kJ/mol/nm2
	}
	else
	{
		temp_r=ReadThisLine(file,4,-1.0,"bond_r","Topology::GetBonds");//equilibruium distance in nm
		if (fabs(temp_r-(-1.0))<LOAD_TOL)//default value is given back, there was no data
		{
			ExtractSemiColon(file);//maybe it was given after the text qualifier
			temp_r=ReadThisLine(file,3,1.0,"bond_r","Topology::GetBonds");//equilibruium distance in nm
			temp_k=ReadThisLine(file,3,1.0,"bond_k","Topology::GetBonds");//force constant in kJ/mol/nm2
		}
		else
		{
			temp_k=ReadThisLine(file,4,-1.0,"bond_k","Topology::GetBonds");//force constant in kJ/mol/nm2
			if (fabs(temp_k-(-1.0))<LOAD_TOL)//default value is given back, there was no data
			{
				ExtractSemiColon(file);//maybe it was given after the text qualifier
				temp_k=ReadThisLine(file,3,1.0,"bond_k","Topology::GetBonds");//force constant in kJ/mol/nm2
			}
		}
	}
	temp_r*=10.0;//converting it to A
	temp_k/=100.0;//converting it to kJ/mol/A2
	if (nbonds_per_type[nmoltype-1]==1)//this is the first bond for this molecule type
	{
		new_type=1;		
	}
	else
	{
		new_type=1;
		
		//determine, whether it is a new bond type for this molecule type. new types has to be handled separately
		//for each molecule type in order to be able to set sigma independently for each moleculetype, if the same bond type
		//occur for different moleculetypes.
		start=(nmoltype==1 ? 0 : nbondtypes_per_type[nmoltype-2]);
		
		for (i=start;i<nbond_types;i++)
		{
			if (fabs(bond_r[i]-temp_r)<LOAD_TOL && fabs(bond_k[i]-temp_k)<LOAD_TOL)//it is the same, as was already given
			{
				bond_type[nbonds-1]=i;
				new_type=0;
				break;
			}
		}
	}
	if (new_type)
	{
		//it is new
		if (nbonds>1)//resize, if it is not the first bond
		{
			ResizeArray(&nbond_types,nbond_types+1,&bond_r,"bond_r","Topology::GetBonds");
			nbond_types--;//setting back to the original, as ResizeArray increases it to the new value
			ResizeArray(&nbond_types,nbond_types+1,&bond_k,"bond_k","Topology::GetBonds");
			nbond_types--;//setting back to the original, as ResizeArray increases it to the new value
			ResizeArray(&nbond_types,nbond_types+1,&bond_sigma,"bond_sigma","Topology::GetBonds");
		}
		nbondtypes_per_type[nmoltype-1]++;//increasing the number of bond types for this molecule type
		cumul_bondtypes[nmoltype]++;//increasing the cumulative number of bond types for this molecule type
		bond_r[nbond_types-1]=temp_r;
		bond_k[nbond_types-1]=temp_k;
		bond_type[nbonds-1]=nbond_types-1;
		if (fabs(RunParams::NB_weight_mode)!=2)
			bond_sigma[nbond_types-1]=GetSigma(file,"bond sigma","Topology::GetBonds");
		else
		{
			//the sigma vdW_weight[0] will be used for all bonded interactions, no need to read
			if (RunParams::vdW_weight[0]>0)
				bond_sigma[nbond_types-1]=RunParams::vdW_weight[0];
			else
				bond_sigma[nbond_types-1]=-1e6;//it will be set later, when it is calculated
			SkipLine(file);
		}
		if (bond_sigma[nbond_types-1]<0 || fabs(bond_sigma[nbond_types-1]-0.0)<TOLERANCE )
			ChiSquared::calc_sigma=1;//if negative or zero
		//bond weight mode will be determined in CheckSigma
		
	}
	
	if (bond_r_max<bond_r[nbond_types-1])
		bond_r_max=bond_r[nbond_types-1];//finding the largest bond (A)
	return 1;
};

int Topology::GetAngles(ifstream &file)
{
	if (RunParams::fnc!=4)
		return(1);//do not read, if flexible molecules are not needed
	int i,new_type,start;
	double temp_ang,temp_k;

	if (nangles==0)
	{
		nangles++;
		nangle_types++;
		SetArraysize(&angle_index,nangles*3,"angle_index","Topology::GetAngles");
		SetArraysize(&angle_type,nangles,"angle_type","Topology::GetAngles");
		SetArraysize(&angle_ang,nangle_types,"angle_ang","Topology::GetAngles");
		SetArraysize(&angle_k,nangle_types,"angle_k","Topology::GetAngles");
		SetArraysize(&angle_sigma,nangle_types,"angle_sigma","Topology::GetAngles");
	}
	else
	{
		ResizeArray(&nangles,nangles+1,&angle_index,"angle_index","Topology::GetAngles",3);
		nangles--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nangles,nangles+1,&angle_type,"angle_type","Topology::GetAngles");			
	}
	//the indices will remain to start with 1!
	angle_index[(nangles-1)*3]=ReadThisLine(file,3,1,"angle_index 1","Topology::GetAngles");//index of the first atom
	angle_index[(nangles-1)*3+1]=ReadThisLine(file,3,1,"angle_index 2","Topology::GetAngles");//index of the second atom
	angle_index[(nangles-1)*3+2]=ReadThisLine(file,3,1,"angle_index 3","Topology::GetAngles");//index of the third atom
	if ((angle_index[(nangles - 1) * 3] == angle_index[(nangles - 1) * 3 + 1]) || (angle_index[(nangles - 1) * 3] == angle_index[(nangles - 1) * 3 + 2]) || (angle_index[(nangles - 1) * 3+1] == angle_index[(nangles - 1) * 3 + 2]))
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The indices (" << angle_index[(nangles - 1) * 3]<<" "<< angle_index[(nangles - 1) * 3+1]<<" "<< angle_index[(nangles - 1) * 3+2]<<") in the "<<nangles << ". angle of the " << nmoltype << ". molecule type contains the same atom index at least twice!" <<  endl;
		cout << "This does not make sense, check the topology file " << currentfilename << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();

	}
	i=ReadThisLine(file,3,1,"angle force function type","Topology::GetAngles");//the force type
	nangles_per_type[nmoltype-1]++;
	cumul_angles[nmoltype]++;
	//some cheking
	if (angle_index[(nangles-1)*3]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The first atom index in the "<<nangles_per_type[nmoltype-1]<<" angle of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<angle_index[(nangles-1)*3]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (angle_index[(nangles-1)*3+1]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The second atom index in the "<<nangles_per_type[nmoltype-1]<<" angle of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<angle_index[(nangles-1)*3+1]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (angle_index[(nangles-1)*3+2]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The third atom index in the "<<nangles_per_type[nmoltype-1]<<" angle of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<angle_index[(nangles-1)*3+2]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (i!=1)//the topology file does not use harmonic angle potential
	{
		cout << "\nWARNING(" << ++warn << "): The angle type (" << i << ") for the " << nangles_per_type[nmoltype - 1] << " angle of the " << nmoltype << endl;
		cout<<"\tmolecule type is not harmonic (1) in the " << currentfilename << " file!" << endl;
		cout<<"\tThe force constant after the text qualifier will be assumed to be for harmoinic potential in kJ/mol/nm2!"<<endl;
		ExtractSemiColon(file);
		temp_ang=ReadThisLine(file,3,1.0,"angle_ang","Topology::GetAngles");//equilibruium angle in degree
		temp_k=ReadThisLine(file,3,1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
	}
	else
	{
		temp_ang=ReadThisLine(file,4,-1.0,"angle_ang","Topology::GetAngles");//equilibruium angle in degree
		if (fabs(temp_ang-(-1.0))<LOAD_TOL)//default value is given back, there was no data
		{
			ExtractSemiColon(file);//maybe it was given after the text qualifier
			temp_ang=ReadThisLine(file,3,1.0,"angle_r","Topology::GetAngles");//equilibruium angle in degree
			temp_k=ReadThisLine(file,3,1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
		}
		else
		{
			temp_k=ReadThisLine(file,4,-1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
			if (fabs(temp_k-(-1.0))<LOAD_TOL)//default value is given back, there was no data
			{
				ExtractSemiColon(file);//maybe it was given after the text qualifier
				temp_k=ReadThisLine(file,3,1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
			}
		}
	}
	if (nangles_per_type[nmoltype-1]==1)//this is the first angle type for this molecule
	{	
		new_type=1;
	}
	else
	{
		new_type=1;
		//determine, whether it is a new angle type
		//for each molecule type in order to be able to set sigma independently for each moleculetype, if the same angle type
		//occur for different moleculetypes.
		start=(nmoltype==1 ? 0 : nangletypes_per_type[nmoltype-2]);
		
		for (i=start;i<nangle_types;i++)
		{
			if (fabs(angle_ang[i]-temp_ang)<LOAD_TOL && fabs(angle_k[i]-temp_k)<LOAD_TOL)//it is the same, as was already given
			{
				angle_type[nangles-1]=i;
				new_type=0;
				break;
			}
		}
	}
	if (new_type)
	{
		  //it is new
		if (nangles>1)
		{
			ResizeArray(&nangle_types,nangle_types+1,&angle_ang,"angle_ang","Topology::GetAngles");
			nangle_types--;//setting back to the original, as ResizeArray increases it to the new value
			ResizeArray(&nangle_types,nangle_types+1,&angle_k,"angle_k","Topology::GetAngles");
			nangle_types--;//setting back to the original, as ResizeArray increases it to the new value
			ResizeArray(&nangle_types,nangle_types+1,&angle_sigma,"angle_sigma","Topology::GetAngles");
		}
		nangletypes_per_type[nmoltype-1]++;//increasing the number of angle types for this molecule type
		cumul_angletypes[nmoltype]++;//increasing the number of angle types for this molecule type
		angle_ang[nangle_types-1]=temp_ang;
		angle_k[nangle_types-1]=temp_k;
		angle_type[nangles-1]=nangle_types-1;
		if (fabs(RunParams::NB_weight_mode)!=2)
			angle_sigma[nangle_types-1]=GetSigma(file,"angle_sigma","Topology::GetAngles");
		else
		{
			//the sigma vdW_weight[0] will be used for all bonded interactions, no need to read
			if (RunParams::vdW_weight[0]>0)
				angle_sigma[nangle_types-1]=RunParams::vdW_weight[0];
			else
				angle_sigma[nangle_types-1]=-1e6;//it will be set later, when it is calculated
			SkipLine(file);
		}
		
		if (angle_sigma[nangle_types-1]<0 || fabs(angle_sigma[nangle_types-1]-0.0)<TOLERANCE)
			ChiSquared::calc_sigma=1;
	}
		
	return 1;
};

int Topology::GetDihedrals(ifstream &file)
{
	if (RunParams::fnc!=4)
		return(1);//do not read, if flexible molecules are not needed
	int i,new_type,*p_ndih_per_type=NULL,temp_mul=0,start;
	double temp_C0,temp_C1,temp_C2,temp_C3,temp_C4,temp_C5,temp_ang,temp_k;
	if (ndihedrals==0)
	{
		ndihedrals++;	
		SetArraysize(&dihedral_index,ndihedrals*4,"dihedral_index","Topology::GetDihedrals");
		SetArraysize(&dihedral_type,ndihedrals,"Dihedral_type","Topology::GetDihedrals");
	}
	else
	{
		ResizeArray(&ndihedrals,ndihedrals+1,&dihedral_index,"dihedral_index","Topology::GetDihedrals",4);
		ndihedrals--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&ndihedrals,ndihedrals+1,&dihedral_type,"dihedral_type","Topology::GetDihedrals");
	}
	//the indices will remain to start with 0!
	dihedral_index[(ndihedrals-1)*4]  =ReadThisLine(file,3,1,"Dihedral_index 1","Topology::GetDihedrals");//index of the first atom
	dihedral_index[(ndihedrals-1)*4+1]=ReadThisLine(file,3,1,"Dihedral_index 2","Topology::GetDihedrals");//index of the second atom
	dihedral_index[(ndihedrals-1)*4+2]=ReadThisLine(file,3,1,"Dihedral_index 3","Topology::GetDihedrals");//index of the third atom
	dihedral_index[(ndihedrals-1)*4+3]=ReadThisLine(file,3,1,"Dihedral_index 4","Topology::GetDihedrals");//index of the 4th atom
	ndihedrals_per_type[nmoltype-1]++;
	//read the type of the dihedral function
	dihedral_type[ndihedrals-1]=ReadThisLine(file,3,1,"dihedral force function type","Topology::GetDihedrals");//the force type
	switch (dihedral_type[ndihedrals-1])
	{
		case 1:
		{
			//this is a periodic dihedral, in case of GROMOS-like force field a proper, in case of OPLS improper
			nperdihedrals_per_type[nmoltype-1]++;
			cumul_perdihedrals[nmoltype]++;
			p_ndih_per_type=&nperdihedrals_per_type[nmoltype-1];
			break;
		}
		case 2:
		{
			//this is a harmonic improper dihedral
			nharmdihedrals_per_type[nmoltype-1]++;
			cumul_harmdihedrals[nmoltype]++;
			p_ndih_per_type=&nharmdihedrals_per_type[nmoltype-1];
			break;
		}
		case 3:
		{
			//this is an RB dihedral
			RB_used = true;//indicating, that RB dihedrals are used
			nRBdihedrals_per_type[nmoltype-1]++;
			cumul_RBdihedrals[nmoltype]++;
			p_ndih_per_type=&nRBdihedrals_per_type[nmoltype-1];
			break;
		}
		default:
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Dihedral type " << dihedral_type[ndihedrals - 1] << " was given in topology in line:" << endl;
			cout << buf << endl;
			cout<<"It can only be 1, 2 or 3."<<endl;
			cout<<"Correct the topology, and try again! Exiting..."<<endl;
			CleanExit();
		}
	}
	//some cheking
	if (dihedral_index[(ndihedrals-1)*4]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The first atom index in the "<<*p_ndih_per_type<<" dihedral of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<dihedral_index[(ndihedrals-1)*4]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (dihedral_index[(ndihedrals-1)*4+1]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The second atom index in the "<<*p_ndih_per_type<<" dihedral of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<dihedral_index[(ndihedrals-1)*4+1]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (dihedral_index[(ndihedrals-1)*4+2]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The third atom index in the "<<*p_ndih_per_type<<" dihedral of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<dihedral_index[(ndihedrals-1)*4+2]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (dihedral_index[(ndihedrals-1)*4+3]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The fourth atom index in the "<<*p_ndih_per_type<<" dihedral of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<dihedral_index[(ndihedrals-1)*4+3]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	
	switch (dihedral_type[ndihedrals-1])
	{
		case 1:
		{
			//this is a periodic dihedral, in case of GROMOS-like force field a proper, i case of OPLS improper
			if (nperdihedrals==0)
			{
				nperdihedrals++;
				nperdihedral_types++;
				SetArraysize(&perdihedral_k,nperdihedral_types,"perdihedral_k","Topology::GetDihedrals");
				SetArraysize(&perdihedral_ang,nperdihedral_types,"perdihedral_ang","Topology::GetDihedrals");
				SetArraysize(&perdihedral_multiplicity,nperdihedral_types,"perdihedral_multiplicity","Topology::GetDihedrals");
				SetArraysize(&perdihedral_sigma,nperdihedral_types,"perdihedral_sigma","Topology::GetDihedrals");
				SetArraysize(&perdihedral_type,nperdihedrals,"perdihedral_type","Topology::GetDihedrals");
			}
			else
				ResizeArray(&nperdihedrals,nperdihedrals+1,&perdihedral_type,"perdihedral_type","Topology::GetDihedrals");
					

			temp_ang=ReadThisLine(file,4,-1.0,"angle_ang","Topology::GetDihedrals");//equilibruium angle in degree
			if (fabs(temp_ang-(-1.0))<LOAD_TOL)//default value is given back, there was no data
			{
				ExtractSemiColon(file);//maybe it was given after the text qualifier
				temp_ang=ReadThisLine(file,3,1.0,"perdihedral_r","Topology::GetDihedrals");//equilibruium angle in degree
				temp_k=ReadThisLine(file,3,1.0,"perdihedral_k","Topology::GetDihedrals");//force constant in kJ/mol/rad2
				temp_mul=ReadThisLine(file,3,2,"perdihedral_multiplicity","Topology::GetDihedrals");//multiplicity
			}
			else
			{
				temp_k=ReadThisLine(file,4,-1.0,"perdihedral_k","Topology::GetDihedrals");//force constant in kJ/mol/rad2
				if (fabs(temp_k-(-1.0))<LOAD_TOL)//default value is given back, there was no data
				{
					ExtractSemiColon(file);//maybe it was given after the text qualifier
					temp_k=ReadThisLine(file,3,1.0,"perdihedral_k","Topology::GetDihedrals");//force constant in kJ/mol/rad2
				}
				else
				{
					temp_mul=ReadThisLine(file,4,-1,"perdihedral_multiplicity","Topology::GetDihedrals");//force constant in kJ/mol/rad2
					if (fabs(temp_k-(-1.0))<LOAD_TOL)//default value is given back, there was no data
					{
						ExtractSemiColon(file);//maybe it was given after the text qualifier
						temp_mul=ReadThisLine(file,3,2,"perdihedral_multiplicity","Topology::GetDihedrals");//force constant in kJ/mol/rad2
					}
				}
			}
			if (force_field_type==OPLSAA && temp_mul!=2)
			{
				cout << "\nWARNING(" << ++warn << "): In case of OPLS-AA force field the improper dihedrals" << endl;
				cout<<"\tare handled as periodic GROMACS type(1) dihedrals, "<<endl;
				cout<<"\tbut the multiplicity has to be 2! It was set to 2, and continuing, but amend the topology later!"<<endl;
				temp_mul=2;
			}


			if (nperdihedrals_per_type[nmoltype-1]==1)//this is the first periodic dihedral for this moleculetype
			{	
				new_type=1;
			}
			else
			{
				new_type=1;
				//determine, whether it is a new per. dih. type. New types has to be handled separately for each molecule type in order to be able to set sigma independently for each molecule type, if the 
				//same perdih type occur in different molecule types.
				start=(nmoltype==1 ? 0 : nperdihtypes_per_type[nmoltype-2]);
				for (i=start;i<nperdihedral_types;i++)
				{
					if (fabs(perdihedral_ang[i]-temp_ang)<LOAD_TOL && fabs(perdihedral_k[i]-temp_k)<LOAD_TOL &&perdihedral_multiplicity[i]==temp_mul)//it is the same, as was already given
					{
						perdihedral_type[nperdihedrals-1]=i;
						new_type=0;
						break;
					}
				}
			}
			if (new_type)
			{
				  //it is new
				if (nperdihedrals>1)//resize, if it is not the first
				{
					ResizeArray(&nperdihedral_types,nperdihedral_types+1,&perdihedral_ang,"perdihedral_ang","Topology::GetDihedrals");
					nperdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nperdihedral_types,nperdihedral_types+1,&perdihedral_k,"perdihedral_k","Topology::GetDihedrals");
					nperdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nperdihedral_types,nperdihedral_types+1,&perdihedral_multiplicity,"perdihedral_multiplicity","Topology::GetDihedrals");
					nperdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nperdihedral_types,nperdihedral_types+1,&perdihedral_sigma,"perdihedral_sigma","Topology::GetDihedrals");
				}
				nperdihtypes_per_type[nmoltype-1]++;//increasing the number of perdihedral types for this molecule type
				cumul_perdihtypes[nmoltype]++;//increasing the cumulative number of perdihedral types for this molecule type
				perdihedral_ang[nperdihedral_types-1]=temp_ang;
				perdihedral_k[nperdihedral_types-1]=temp_k;
				perdihedral_multiplicity[nperdihedral_types-1]=temp_mul;
				perdihedral_type[nperdihedrals-1]=nperdihedral_types-1;
				if (fabs(RunParams::NB_weight_mode) != 2)
					perdihedral_sigma[nperdihedral_types-1]=GetSigma(file,"perdihedral_sigma","Topology::GetDihedrals");
				else
				{
					//the sigma vdW_weight[0] will be used for all bonded interactions, no need to read
					if (RunParams::vdW_weight[0]>0)
						perdihedral_sigma[nperdihedral_types-1]=RunParams::vdW_weight[0];
					else
						perdihedral_sigma[nperdihedral_types-1]=-1e6;//it will be set later, when it is calculated
					SkipLine(file);
				}
				if (perdihedral_sigma[nperdihedral_types-1]<0 || fabs(perdihedral_sigma[nperdihedral_types-1]-0.0)<TOLERANCE )
					ChiSquared::calc_sigma=1;
			}
			break;
		}
		case 2:
		{
			if (force_field_type==OPLSAA)
			{
				cout<<"\n****ERROR****"<<endl;
				cout << "In line: " << endl;
				cout << buf << endl;
				cout<<"Improper dihedrals GROMACS type (2) cannot be used with OPLS-AA force field!"<<endl;
				cout<<"Proper dihedral GROMACS type (1) has to be used instead, check the ffoplsaabon.itp "<<endl;
				cout<<"for the  available types, the central atom of the planar group has to be the second or third index!"<<endl; 
				cout<<"Fix the topology file, and try again, exiting..."<<endl;
				CleanExit();
			}
			//this is a harmonic improper dihedral
			if (nharmdihedrals==0)
			{
				nharmdihedrals++;
				nharmdihedral_types++;
				SetArraysize(&harmdihedral_k,nharmdihedral_types,"harmdihedral_k","Topology::GetDihedrals");
				SetArraysize(&harmdihedral_ang,nharmdihedral_types,"harmdihedral_ang","Topology::GetDihedrals");
				SetArraysize(&harmdihedral_sigma,nharmdihedral_types,"harmdihedral_sigma","Topology::GetDihedrals");
				SetArraysize(&harmdihedral_type,nharmdihedrals,"harmdihedral_type","Topology::GetDihedrals");
			}
			else
				ResizeArray(&nharmdihedrals,nharmdihedrals+1,&harmdihedral_type,"harmdihedral_type","Topology::GetDihedrals");
					

			temp_ang=ReadThisLine(file,4,-1.0,"angle_ang","Topology::GetAngles");//equilibruium angle in degree
			if (fabs(temp_ang-(-1.0))<LOAD_TOL)//default value is given back, there was no data
			{
				ExtractSemiColon(file);//maybe it was given after the text qualifier
				temp_ang=ReadThisLine(file,3,1.0,"angle_r","Topology::GetAngles");//equilibruium angle in degree
				temp_k=ReadThisLine(file,3,1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
			}
			else
			{
				temp_k=ReadThisLine(file,4,-1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
				if (fabs(temp_k-(-1.0))<LOAD_TOL)//default value is given back, there was no data
				{
					ExtractSemiColon(file);//maybe it was given after the text qualifier
					temp_k=ReadThisLine(file,3,1.0,"angle_k","Topology::GetAngles");//force constant in kJ/mol/rad2
				}
			}
			if (nharmdihedrals_per_type[nmoltype-1]==1)//this is the first harmonic dihedral for this moleculetype
			{	
				new_type=1;
			}
			else
			{
				new_type=1;
				//determine, whether it is a new per. dih. type. New types has to be handled separately for each molecule type in order to be able to set sigma independently for each molecule type, if the 
				//same perdih type occur in different molecule types.
				start=(nmoltype==1 ? 0 : nharmdihtypes_per_type[nmoltype-2]);
				for (i=start;i<nharmdihedral_types;i++)
				{
					if (fabs(harmdihedral_ang[i]-temp_ang)<LOAD_TOL && fabs(harmdihedral_k[i]-temp_k)<LOAD_TOL)//it is the same, as was already given
					{
						harmdihedral_type[nharmdihedrals-1]=i;
						new_type=0;
						break;
					}
				}
			}
			if (new_type)
			{
				  //it is new
				if (nharmdihedrals>1)//resize, if it is not the first
				{
					ResizeArray(&nharmdihedral_types,nharmdihedral_types+1,&harmdihedral_ang,"harmdihedral_ang","Topology::GetDihedrals");
					nharmdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nharmdihedral_types,nharmdihedral_types+1,&harmdihedral_k,"harmdihedral_k","Topology::GetDihedrals");
					nharmdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nharmdihedral_types,nharmdihedral_types+1,&harmdihedral_sigma,"harmdihedral_sigma","Topology::GetDihedrals");
				}
				nharmdihtypes_per_type[nmoltype-1]++;//increasing the number of harmonic dihedral types for this molecule type
				cumul_harmdihtypes[nmoltype]++;//increasing the cumulative number of perdihedral types for this molecule type
				harmdihedral_ang[nharmdihedral_types-1]=temp_ang;
				harmdihedral_k[nharmdihedral_types-1]=temp_k;
				harmdihedral_type[nharmdihedrals-1]=nharmdihedral_types-1;
				
				if (fabs(RunParams::NB_weight_mode) != 2)
					harmdihedral_sigma[nharmdihedral_types-1]=GetSigma(file,"harmdihedral_sigma","Topology::GetDihedrals");
				else
				{
					//the sigma vdW_weight[0] will be used for all bonded interactions, no need to read
					if (RunParams::vdW_weight[0]>0)
						harmdihedral_sigma[nharmdihedral_types-1]=RunParams::vdW_weight[0];
					else
						harmdihedral_sigma[nharmdihedral_types-1]=-1e6;//it will be set later, when it is calculated
					SkipLine(file);
				}
				if (harmdihedral_sigma[nharmdihedral_types-1]<0 || fabs(harmdihedral_sigma[nharmdihedral_types-1]-0.0)<TOLERANCE )
					ChiSquared::calc_sigma=1;
			}
			break;
		}
		case 3:
		{
			//this is an RB dihedral
			if (nRBdihedrals==0)
			{
				nRBdihedrals++;
				nRBdihedral_types++;
				SetArraysize(&RBdihedral_C0,nRBdihedral_types,"RBdihedral_C0","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_C1,nRBdihedral_types,"RBdihedral_C1","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_C2,nRBdihedral_types,"RBdihedral_C2","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_C3,nRBdihedral_types,"RBdihedral_C3","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_C4,nRBdihedral_types,"RBdihedral_C4","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_C5,nRBdihedral_types,"RBdihedral_C5","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_type,nRBdihedral_types,"RBdihedral_type","Topology::GetDihedrals");
				SetArraysize(&RBdihedral_sigma,nRBdihedral_types,"RBdihedral_sigma","Topology::GetDihedrals");
			}
			else
				ResizeArray(&nRBdihedrals,nRBdihedrals+1,&RBdihedral_type,"RBdihedral_type","Topology::GetDihedrals");
			
			//-10000.0 is well out of the usual range of constant parameters, so it is used to indicate, that the reading was not successful
			temp_C0=ReadThisLine(file,4,-10000.0,"RBdihedral_C0","Topology::GetDihedrals");//C0 constant in kJ/mol
			if (fabs(temp_C0-(-10000.0))<LOAD_TOL)//default value is given back, there was no data
			{
				ExtractSemiColon(file);
				temp_C0=ReadThisLine(file,3,0.0,"RBdihedral_C0","Topology::GetDihedrals");//C0 constant in kJ/mol
				temp_C1=ReadThisLine(file,2,0.0,"RBdihedral_C1","Topology::GetDihedrals");//C1 constant in kJ/mol
				temp_C2=ReadThisLine(file,2,0.0,"RBdihedral_C2","Topology::GetDihedrals");//C2 constant in kJ/mol
				temp_C3=ReadThisLine(file,2,0.0,"RBdihedral_C3","Topology::GetDihedrals");//C3 constant in kJ/mol
				temp_C4=ReadThisLine(file,2,0.0,"RBdihedral_C4","Topology::GetDihedrals");//C4 constant in kJ/mol
				temp_C5=ReadThisLine(file,2,0.0,"RBdihedral_C5","Topology::GetDihedrals");//C5 constant in kJ/mol
			}
			else
			{
				temp_C1=ReadThisLine(file,2,0.0,"RBdihedral_C1","Topology::GetDihedrals");//C1 constant in kJ/mol
				temp_C2=ReadThisLine(file,2,0.0,"RBdihedral_C2","Topology::GetDihedrals");//C2 constant in kJ/mol
				temp_C3=ReadThisLine(file,2,0.0,"RBdihedral_C3","Topology::GetDihedrals");//C3 constant in kJ/mol
				temp_C4=ReadThisLine(file,2,0.0,"RBdihedral_C4","Topology::GetDihedrals");//C4 constant in kJ/mol
				temp_C5=ReadThisLine(file,2,0.0,"RBdihedral_C5","Topology::GetDihedrals");//C5 constant in kJ/mol
			}
			if (nRBdihedrals_per_type[nmoltype-1]==1)//this is the first RB dihedral for this moleculetype
			{	
				new_type=1;
			}
			else
			{
				new_type=1;
				//determine, whether it is a new RB dih. type. New types has to be handled separately for each molecule type in order to be able to set sigma independently for each molecule type, if the 
				//same RB dih type occur in different molecule types.
				start=(nmoltype==1 ? 0 : nRBdihtypes_per_type[nmoltype-2]);
				for (i=start;i<nRBdihedral_types;i++)
				{
					if (fabs(RBdihedral_C0[i]-temp_C0)<LOAD_TOL && fabs(RBdihedral_C1[i]-temp_C1)<LOAD_TOL && \
						fabs(RBdihedral_C2[i]-temp_C2)<LOAD_TOL && fabs(RBdihedral_C3[i]-temp_C3)<LOAD_TOL && \
						fabs(RBdihedral_C4[i]-temp_C4)<LOAD_TOL && fabs(RBdihedral_C5[i]-temp_C5)<LOAD_TOL)//it is the same, as was already given
					{
						RBdihedral_type[nRBdihedrals-1]=i;
						new_type=0;
						break;
					}
				}
			}
			if (new_type)
			{
				  //it is new
				if (nRBdihedrals>1)//resize, if it is not the first
				{
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_C0,"RBdihedral_C0","Topology::GetDihedrals");
					nRBdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_C1,"RBdihedral_C1","Topology::GetDihedrals");
					nRBdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_C2,"RBdihedral_C2","Topology::GetDihedrals");
					nRBdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_C3,"RBdihedral_C3","Topology::GetDihedrals");
					nRBdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_C4,"RBdihedral_C4","Topology::GetDihedrals");
					nRBdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_C5,"RBdihedral_C5","Topology::GetDihedrals");
					nRBdihedral_types--;//setting back to the original, as ResizeArray increases it to the new value
					ResizeArray(&nRBdihedral_types,nRBdihedral_types+1,&RBdihedral_sigma,"RBdihedral_sigma","Topology::GetDihedrals");
				}				
				nRBdihtypes_per_type[nmoltype-1]++;//increasing the number of RB dihedral types for this molecule type
				cumul_RBdihtypes[nmoltype]++;//increasing the cumulative number of RBdihedral types for this molecule type
				RBdihedral_C0[nRBdihedral_types-1]=temp_C0;
				RBdihedral_C1[nRBdihedral_types-1]=temp_C1;
				RBdihedral_C2[nRBdihedral_types-1]=temp_C2;
				RBdihedral_C3[nRBdihedral_types-1]=temp_C3;
				RBdihedral_C4[nRBdihedral_types-1]=temp_C4;
				RBdihedral_C5[nRBdihedral_types-1]=temp_C5;
				RBdihedral_type[nRBdihedrals-1]=nRBdihedral_types-1;
				
				if (fabs(RunParams::NB_weight_mode) != 2)
					RBdihedral_sigma[nRBdihedral_types-1]=GetSigma(file,"RBdihedral_sigma","Topology::GetDihedrals");
				else
				{
					//the sigma vdW_weight[0] will be used for all bonded interactions, no need to read
					if (RunParams::vdW_weight[0]>0)
						RBdihedral_sigma[nRBdihedral_types-1]=RunParams::vdW_weight[0];
					else
						RBdihedral_sigma[nRBdihedral_types-1]=-1e6;//it will be set later, when it is calculated
					SkipLine(file);
				}
				if (RBdihedral_sigma[nRBdihedral_types-1]<0 || fabs(RBdihedral_sigma[nRBdihedral_types-1]-0.0)<TOLERANCE )
					ChiSquared::calc_sigma=1;
			}
			break;
		}
	}
	
	return 1;
};

//reading the atom pairs defining the 1-4 interaction pairs (only type 1 is handled)
//No sigma for RMC is necessary, as the same values given for the normal vdW and Coulomb interaction will be used  
int Topology::GetPairs(ifstream &file)
{
	int i,new_type, ind1, ind2;
	double temp_pot1,temp_pot2;//for reading the force constant and the distance

	//if no potential is used, skip
	if (RunParams::potential==0)
		return 0;
	//decide, whether it is a (1) type pair interaction, as the others are not handled
	ind1=ReadThisLine(file,3,1,"pair_index 1","Topology::GetPairs");//index of the first atom
	ind2=ReadThisLine(file,3,1,"pair_index 2","Topology::GetPairs");//index of the second atom
	i=ReadThisLine(file,3,1,"pair force function type","Topology::GetPairs");//the force type
	if (i!=1)//the topology file does not use (1) type 1-4 interaction pair potential
		return 0;//this line will not be read
	
	//Check, whether the 1-4 pairs are excluded from the neighbourlist
	if (nexclusions_per_type[nmoltype-1]<3)
	{
		cout << "\nWARNING(" << ++warn << "): The exclusion number is " << nexclusions_per_type[nmoltype - 1] << " for molecule type " << molname[nmoltype - 1] << endl;
		cout<<"\t and pairs are given for 1-4 interactin calculation. This way the 1-4 intercation would be calculated twice,"<<endl;
		cout<<"\tfirst as normal NB interaction, than as 1-4 interaction. Setting nexclusion to 3!!!"<<endl;
		nexclusions_per_type[nmoltype-1]=3;
	}
	//Create the array or resize it
	if (npairs==0)
	{
		npairs++;	
		npair_types++;
		SetArraysize(&pair_index,npairs*2,"pair_index","Topology::GetPairs");
		SetArraysize(&pair_type,npairs,"pair_type","Topology::GetPairs");
		SetArraysize(&vdW14_pot1,npair_types,"vdW14_pot1","Topology::GetPairs");
		SetArraysize(&vdW14_pot2,npair_types,"vdW14_pot2","Topology::GetPairs");
	}
	else
	{
		ResizeArray(&npairs,npairs+1,&pair_index,"pair_index","Topology::GetPairs",2);
		npairs--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&npairs,npairs+1,&pair_type,"pair_type","Topology::GetPairs");
	}
	//indices will start with 1!
	pair_index[(npairs-1)*2]=ind1;
	pair_index[(npairs-1)*2+1]=ind2;
	
	npairs_per_type[nmoltype-1]++;
	cumul_pairs[nmoltype]++;
	//some cheking
	if (pair_index[(npairs-1)*2]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout<<"The first atom index in the "<<npairs_per_type[nmoltype-1]<<" 1-4 pair of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<pair_index[(npairs-1)*2]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (pair_index[(npairs-1)*2+1]>natoms_per_type[nmoltype-1])
	{
		cout<<"\n****ERROR****"<<endl;
		cout << "The line: " << endl;
		cout << buf << endl;
		cout << "was processed with [ pairs ] directive!" << endl;
		cout<<"The second atom index in the "<<npairs_per_type[nmoltype-1]<<" 1-4 pair  of the "<<nmoltype<<" molecule type is "<<endl;
		cout<<pair_index[(npairs-1)*2+1]<<", which is greater, than the number of atoms in the molecule ("<<natoms_per_type[nmoltype-1]<<")!"<<endl;
		cout<<"Check the file "<<currentfilename<<"! Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	//reading the 1-4 interaction parameters, if necessary
	if (gen_pairs==0)//pairs will not be generated, have to be read
	{
		temp_pot1=ReadThisLine(file,4,-1.0,"vdW14_pot1","Topology::GetPairs");//equilibruium distance in nm
		if (fabs(temp_pot1-(-1.0))<LOAD_TOL)//default value is given back, there was no data
		{
			ExtractSemiColon(file);//maybe it was given after the text qualifier
			temp_pot1=ReadThisLine(file,3,1.0,"vdW14_pot1","Topology::GetPairs");// pot parameter1 according to the vdW type and comb rule
			temp_pot2=ReadThisLine(file,3,1.0,"vdW14_pot2","Topology::GetPairs");// pot parameter2 according to the vdW type and comb rule
		}
		else
		{
			temp_pot2=ReadThisLine(file,4,-1.0,"vdW14_pot2","Topology::GetPairs");// pot parameter2 according to the vdW type and comb rule
			if (fabs(temp_pot2-(-1.0))<LOAD_TOL)//default value is given back, there was no data
			{
				ExtractSemiColon(file);//maybe it was given after the text qualifier
				temp_pot2=ReadThisLine(file,3,1.0,"vdW14_pot2","Topology::GetPairs");// pot parameter2 according to the vdW type and comb rule
			}
		}
		if (RunParams::potential==1)
		{
			if (RunParams::vdW_comb_rule == 0 || RunParams::vdW_comb_rule==1)//LJ, C6, C12 is given (nm->A conversion necessary
			{
				temp_pot1*=1.0e6;
				temp_pot2*=1.0e12;
			}
			else
			{
				//sigma and epsilon is given
				temp_pot1*=10.0;//converting it to A (it is in GROMACS unit, nm 
			}
		}
	}
	else
	{
		temp_pot1=-1;//for some forcefields (OPLS, AMBER, CHARMM..) it will be calculated with scaling the normal interaction parameters
		temp_pot2=-1;
	}
	

	if (npairs==1)
	{	
		vdW14_pot1[0]=temp_pot1;
		vdW14_pot2[0]=temp_pot2;
		pair_type[0]=0;
	}
	else
	{
		new_type=1;
		//determine, whether it is a new pair type
		if (gen_pairs)
		{
			//there will be only 1 pair type, as all the 1-4 interactions are calculated by scaling 
			new_type=0;
			pair_type[npairs-1]=pair_type[npairs-2];
		}
		else
		{
			//new type will be decided based on the pot1 and pot2
			for (i=0;i<npair_types;i++)
			{
				if (fabs(vdW14_pot1[i]-temp_pot1)<LOAD_TOL && fabs(vdW14_pot2[i]-temp_pot2)<LOAD_TOL )
				{
					pair_type[npairs-1]=i;
					new_type=0;
					break;
				}
			}
		}
		if (new_type)
		{
			//it is new
			ResizeArray(&npair_types,npair_types+1,&vdW14_pot1,"vdW14_pot1","Topology::GetPairs");
			npair_types--;//setting back to the original, as ResizeArray increases it to the new value
			ResizeArray(&npair_types,npair_types+1,&vdW14_pot2,"vdW14_pot2","Topology::GetPairs");
			vdW14_pot1[npair_types-1]=temp_pot1;
			vdW14_pot2[npair_types-1]=temp_pot2;
			pair_type[npairs-1]=npair_types-1;
		}
	}
	return 1;
};

//reading the exclusions
int Topology::GetExclusions(ifstream &file)
{
	int i,*temp,il,found,*excl;
	
	SetArraysize(&temp, INPUT_SIZE, "temp", "Topology::GetExclusions");
	//there can be unknown number of exclusion in a line
	
	if (nexclusionlines == 0)
	{
		nexclusionlines++;
		SetArraysize(&nexcluded_neighs, nexclusionlines, "nexcluded_neighs","Topology::GetExclusions");
		SetArraysize(&cumul_nexcluded_neighs, nexclusionlines+1, "cumul_nexcluded_neighs", "Topology::GetExclusions");
		SetArraysize(&exclusion_centre, nexclusionlines, "exclusion_centre", "Topology::GetExclusions");
		cumul_nexcluded_neighs[0] = 0;
	}
	else
	{
		ResizeArray(&nexclusionlines, nexclusionlines + 1, &nexcluded_neighs, "nexcluded_neighs", "Topology::GetExclusions");
		nexclusionlines--;//setting back to the original, as ResizeArray increases it to the new value
		ResizeArray(&nexclusionlines, nexclusionlines + 1, &exclusion_centre, "exclusion_centre", "Topology::GetExclusions");
		ResizeArray(&nexclusionlines, nexclusionlines+1, &cumul_nexcluded_neighs, "cumul_nexcluded_neighs", "Topology::GetExclusions");
		nexclusionlines--;//setting back to the original, as ResizeArray increases it to the new value
	}
	
	//indices will start with 1!
	exclusion_centre[nexclusionlines - 1] = ReadThisLine(file, 3, 1, "index of atom with exclusions", "Topology::Exclusions");//index of the atom for which the exclusions are given
	nexcluded_neighs[nexclusionlines - 1] = 0;
	do
	{
		if (nexcluded_neighs[nexclusionlines - 1] == INPUT_SIZE)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "In line: " << endl;
			cout << buf << endl;
			cout << "The number of excluded neighbours is larger than INPUT_SIZE=" << INPUT_SIZE << " in Topology::Exclusions" << endl;
			cout << "Cannot read all the exclusion indices in " << currentfilename << " file " << nexclusionlines << ". exclusions line" << endl;
			cout << "Incease the value of INPUT_SIZE in units.h and try again!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();

		}
		temp[nexcluded_neighs[nexclusionlines - 1]] = ReadThisLine(file, 4, -100, "index of excluded atom", "Topology::Exclusions");//index of excluded atom 
		found = 0;
		if (temp[nexcluded_neighs[nexclusionlines - 1]] == exclusion_centre[nexclusionlines - 1])
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "In line: " << endl;
			cout << buf << endl;
			cout << "The central index was given as its own exclusion! Most probably the indices were mixed up in file " << currentfilename << endl;
			cout << "Correct the topology and try again!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}
		//check, whether the same excluded atom is not given twice for the same exclusion centre
		//there can be more exclusion lines for one exclusion centre
		if (temp[nexcluded_neighs[nexclusionlines - 1]] != -100)//end of line
		{
			for (il = 0; il < nexclusionlines_per_type[nmoltype - 1] + 1; il++)
			{
				if (exclusion_centre[nexclusionlines - 1] == exclusion_centre[cumul_exclusionlines[nmoltype - 1] + il])//this is for the same exclusion centre
				{
					if (il == nexclusionlines_per_type[nmoltype - 1])//this is the present line
						excl = temp;
					else
						excl = exclusions + cumul_nexcluded_neighs[cumul_exclusionlines[nmoltype - 1] + il];
					for (i = 0; i < nexcluded_neighs[cumul_exclusionlines[nmoltype - 1] + il]; i++)
					{
						if (excl[i] == temp[nexcluded_neighs[nexclusionlines - 1]])
						{
							found = 1;
							cout << "\nWARNING(" << ++warn << "): In line: " << endl;
							cout << "\t"<<buf << endl;
							cout << "\tThe excluded atom index " << temp[nexcluded_neighs[nexclusionlines - 1]] << " was given previously in the ";
							if (nexclusionlines_per_type[nmoltype - 1] == il)
								cout << "same";
							else
								cout << il + 1<<".";
							cout << " line" << endl;
							cout << "\tfor the " << nmoltype << " molecule type in file " << currentfilename << "!" << endl;
							cout << "\tIt will be ignored and continuing, but correct the topology!" << endl;
						}
					}
				}
			}
		}
		if (found==0)
			nexcluded_neighs[nexclusionlines - 1]++;
	} while (temp[nexcluded_neighs[nexclusionlines - 1]-1] != -100);
	nexcluded_neighs[nexclusionlines - 1]--;
	if (nexclusionlines == 1)//first exclusion line
	{
		SetArraysize(&exclusions, nexcluded_neighs[nexclusionlines - 1], "exclusions", "Topology::GetExclusions");
		cumul_nexcluded_neighs[nexclusionlines] = nexcluded_neighs[nexclusionlines - 1];
	}
	else
	{
		//this will set the cumul_nexcluded_neighs+nexclusionlines[nexclusionlines]
		cumul_nexcluded_neighs[nexclusionlines] = cumul_nexcluded_neighs[nexclusionlines - 1];
		ResizeArray(cumul_nexcluded_neighs + nexclusionlines, cumul_nexcluded_neighs[nexclusionlines - 1] + nexcluded_neighs[nexclusionlines - 1], &exclusions, "exclusions", "Topology::GetExclusion");
	}
	for (i = 0; i < nexcluded_neighs[nexclusionlines - 1]; i++)
		exclusions[cumul_nexcluded_neighs[nexclusionlines-1]+i] = temp[i];
	nexclusionlines_per_type[nmoltype - 1]++;
	cumul_exclusionlines[nmoltype]++;
	delete[] temp;

	return 1;

};
//reading the virtual_sites2, the virtual site located on the line between the atoms from a fraction from the first atom 
int Topology::GetVirtual2(ifstream &file)
{
	int i;
	char *name1, *name2,*name3,*conv_numb=NULL;

	SetArraysize(&name1, NAME_SIZE , "name1", "Topology::GetVirtual2");
	SetArraysize(&name2, NAME_SIZE, "name2", "Topology::GetVirtual2");
	SetArraysize(&name3, NAME_SIZE, "name3", "Topology::GetVirtual2");
	
	

	//Create the array or resize it
	if (n_gr_virtuals == 0)
	{
		n_gr_virtuals++;
		SetArraysize(&virtual_site, n_gr_virtuals, "virtual_site", "Topology::GetVirtual2");
	}
	else
	{
		ResizeVirtualSiteStruct(&n_gr_virtuals, n_gr_virtuals + 1, &virtual_site, "virtual_site", "Topology::GetVirtual2");
		
	}
	   	
	mystrcpy(name1, NAME_SIZE, "virtual_site2[");
	IntToStr(&conv_numb, n_gr_virtuals);

	mystrcat(name1, NAME_SIZE, conv_numb);
	mystrcat(name1, NAME_SIZE, "] field: ");
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "index of the virtual site");

	virtual_site[n_gr_virtuals - 1].atom_index=ReadThisLine(file,3,0,name2,"Topology::GetVirtual2", currentfilename);//read the atom indices making up the sites
	nvirtuals_per_type[nmoltype - 1]++;
	cumul_virtuals[nmoltype]++;
	for (i = 0; i < 2; i++)
	{
		mystrcpy(name2, NAME_SIZE, name1);
		mystrcat(name2, NAME_SIZE, "index of the ");
		IntToStr(&conv_numb, i + 1);
		mystrcat(name2, NAME_SIZE, conv_numb);
		mystrcat(name2, NAME_SIZE, ". atom creating the virtual site");
		virtual_site[n_gr_virtuals - 1].indices[i]=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual2", currentfilename);//reading the indices

		if (virtual_site[n_gr_virtuals - 1].indices[i] > natoms_per_type[nmoltype - 1])
		{
			cout << "\n****ERROR****" << endl;
			cout << "In line: " << endl;
			cout << buf << endl;
			cout << "The " << i + 1 << ". atom index in the " << nvirtuals_per_type[nmoltype - 1] << ". 2-atoms virtual site of " << nmoltype << ". molecule type is " << endl;
			cout << virtual_site[n_gr_virtuals - 1].indices[i] << ", which is greater, than the number of atoms in the molecule (" << natoms_per_type[nmoltype - 1] << ")!" << endl;
			cout << "Check the file " << currentfilename << "!" << endl;
			cout<<"Cannot run this way, exiting..." << endl;
			CleanExit();
		}
	}
	mystrcpy(name2, NAME_SIZE, "[ virtual_sites2 ]");
	CheckRepetition(virtual_site[n_gr_virtuals - 1].indices, 2, nvirtuals_per_type[nmoltype - 1], nmoltype, name2);
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "function type");
	virtual_site[n_gr_virtuals - 1].type=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual2", currentfilename); //reading the function type
	if (virtual_site[n_gr_virtuals - 1].type !=1)
	{
		cout << "\n****ERROR****" << endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout << "The function type in the " << nvirtuals_per_type[nmoltype - 1] << ". 2-atoms virtual site of " << nmoltype << ". molecule type is " << endl;
		cout << virtual_site[n_gr_virtuals - 1].type << ", and it can only be 1!" << endl;
		cout << "Check the file " << currentfilename << "!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	if (virtual_site[n_gr_virtuals - 1].type!= 1)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout << "For virtual_sites2 only function type 1 can be given!"<<endl;
		cout << "Check your parameter file! Exiting..." << endl;
		CleanExit();
	}
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "parameter a");
	virtual_site[n_gr_virtuals - 1].params[0]=ReadThisLine(file, 3, 0.0, name2, "Topology::GetVirtual2", currentfilename);// read the parameter a
	virtual_site[n_gr_virtuals - 1].assigned= -1;
	virtual_site[n_gr_virtuals - 1].tnumb = 2;
	virtual_site[n_gr_virtuals - 1].params[1] = -1;//not used
	virtual_site[n_gr_virtuals - 1].params[2] = -1;//not used
	virtual_site[n_gr_virtuals - 1].indices[2] = -1;//not used
	virtual_site[n_gr_virtuals - 1].indices[3] = -1;//not used
	if (ExtractSemiColon(file, 0))//semicolon was give, try to read
	{
		virtual_site[n_gr_virtuals - 1].host_index = ReadThisLine(file, 2, -1, "host_index", "Topology::GetVirtual2", currentfilename);// read host atom index
		if (virtual_site[n_gr_virtuals - 1].host_index==-1)
		{
			cout << "\nWARNING(" << ++warn << "): Altough semicolon was given, but no host index was specified for [ virtual_sites2 ] in line:" << endl;
			cout <<"\t"<< buf << endl;
			cout << "\tThe host index, will be the index of the first atom (" << virtual_site[n_gr_virtuals - 1].indices[0] << ") building the virtual site" << endl;
			virtual_site[n_gr_virtuals - 1].host_index = virtual_site[n_gr_virtuals - 1].indices[0];
		}
	}
	else
	{
		cout << "\nWARNING(" << ++warn << "): No host index was specified for [ virtual_sites2 ] in line:" << endl;
		cout <<"\t"<< buf << endl;
		cout << "\tThe host index, will be the index of the first atom (" << virtual_site[n_gr_virtuals - 1].indices[0] << ") building the virtual site" << endl;
		virtual_site[n_gr_virtuals - 1].host_index = virtual_site[n_gr_virtuals - 1].indices[0];
	}
	delete[] name1;
	delete[] name2;
	delete[] name3;
	if (conv_numb != NULL)
		delete[] conv_numb;
	return 1;
};

//reading the virtual_sites3, the virtual site located on the line between the atoms from a fraction from the first atom 
int Topology::GetVirtual3(ifstream &file)
{
	int i,k;
	char *name1, *name2, *name3,*conv_numb=NULL;

	SetArraysize(&name1, NAME_SIZE, "name1", "Topology::GetVirtual3");
	SetArraysize(&name2, NAME_SIZE, "name2", "Topology::GetVirtual3");
	SetArraysize(&name3, NAME_SIZE, "name3", "Topology::GetVirtual3");
	
	//Create the array or resize it
	if (n_gr_virtuals == 0)
	{
		n_gr_virtuals++;
		SetArraysize(&virtual_site, n_gr_virtuals, "virtual_site", "Topology::GetVirtual3");
	}
	else
	{
		ResizeVirtualSiteStruct(&n_gr_virtuals, n_gr_virtuals + 1, &virtual_site, "virtual_site", "Topology::GetVirtual3");

	}

	mystrcpy(name1, NAME_SIZE, "virtual_site3[");
	IntToStr(&conv_numb, n_gr_virtuals);

	mystrcat(name1, NAME_SIZE, conv_numb);
	mystrcat(name1, NAME_SIZE, "] field: ");
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "index of the virtual site");

	virtual_site[n_gr_virtuals - 1].atom_index=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual3", currentfilename);//read the atom indices making up the sites
	nvirtuals_per_type[nmoltype - 1]++;
	cumul_virtuals[nmoltype]++;
	
	for (i = 0; i < 3; i++)
	{
		mystrcpy(name2, NAME_SIZE, name1);
		mystrcat(name2, NAME_SIZE, "index of the ");
		IntToStr(&conv_numb, i + 1);
		mystrcat(name2, NAME_SIZE, conv_numb);
		mystrcat(name2, NAME_SIZE, ". atom creating the virtual site");
		virtual_site[n_gr_virtuals - 1].indices[i]=ReadThisLine(file, 3,0, name2, "Topology::GetVirtual3", currentfilename);//reading the indices


		if (virtual_site[n_gr_virtuals - 1].indices[i] > natoms_per_type[nmoltype - 1])
		{
			cout << "\n****ERROR****" << endl;
			cout << "In line: " << endl;
			cout << buf << endl;
			cout << "The " << i + 1 << ". atom index in the " << nvirtuals_per_type[nmoltype - 1] << ". line of [ virtual_sites3 ] section of " << nmoltype << ". molecule type is " << endl;
			cout << virtual_site[n_gr_virtuals - 1].indices[i] << ", which is greater, than the number of atoms in the molecule (" << natoms_per_type[nmoltype - 1] << ")!" << endl;
			cout << "Check the file " << currentfilename << "!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}
	}
	mystrcpy(name2, NAME_SIZE, "[ virtual_sites3 ]");
	CheckRepetition(virtual_site[n_gr_virtuals - 1].indices, 3, nvirtuals_per_type[nmoltype - 1], nmoltype,name2);
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "function type");
	virtual_site[n_gr_virtuals - 1].type=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual3", currentfilename); //reading the function type
	if (virtual_site[n_gr_virtuals - 1].type < 1 || virtual_site[n_gr_virtuals - 1].type>4)
	{
		cout << "\n****ERROR****" << endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout << "The function type in the " << nvirtuals_per_type[nmoltype - 1] << ". 3-atoms virtual site of " << nmoltype << ". molecule type is " << endl;
		cout << virtual_site[n_gr_virtuals - 1].type << ", and it can only be 1, 2, 3 or 4!" << endl;
		cout << "Check the file " << currentfilename << "!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	for (i = 0; i < 3; i++)
		virtual_site[n_gr_virtuals - 1].params[i] = -1;//not used
	k = 2;
	if (virtual_site[n_gr_virtuals - 1].type==4)
		k = 3;
	
	for (i = 0; i < k; i++)
	{
		mystrcpy(name2, NAME_SIZE, name1);
		switch (i)
		{
		case 0:
			if (virtual_site[n_gr_virtuals - 1].type == 3)
				mystrcat(name2, NAME_SIZE, "parameter theta");
			else
				mystrcat(name2, NAME_SIZE, "parameter a");
			break;

		case 1:
			if (virtual_site[n_gr_virtuals - 1].type == 2 || virtual_site[n_gr_virtuals - 1].type == 3)
				mystrcat(name2, NAME_SIZE, "parameter d");
			else
				mystrcat(name2, NAME_SIZE, "parameter b");
			break;
		default:
			mystrcat(name2, NAME_SIZE, "parameter c");
			break;

		}

		virtual_site[n_gr_virtuals - 1].params[i]=ReadThisLine(file, 3, 0.0, name2, "Topology::GetVirtual4", currentfilename);// read the parameter a
		//The GROMACS part of the topology uses nm instead of Angstrom as in RMC, convert the parameters if necessary to Ang or 1/Ang
		if (i == 1)//second param  
		{
			if (virtual_site[n_gr_virtuals - 1].type == 2 || virtual_site[n_gr_virtuals - 1].type == 3)
				virtual_site[n_gr_virtuals - 1].params[i] *= 10.0/SimpleCfg::boxedge;//the parameter was in nm, now in reduced unit
		}
		else
		{
			if (i==2 )
				virtual_site[n_gr_virtuals - 1].params[i] *=  SimpleCfg::boxedge/10.0;//the lenth of the cross product scales with distance2. 
					//the equation, which defines the position is given in nm, this conversion is needed, if the virtual site coordinate is calculated in reduced unit
					
		}
		
	}
	if (ExtractSemiColon(file, 0))//semicolon was give, try to read
	{
		virtual_site[n_gr_virtuals - 1].host_index = ReadThisLine(file, 2, -1, "host_index", "Topology::GetVirtual3", currentfilename);// read host atom index
		if (virtual_site[n_gr_virtuals - 1].host_index == -1)
		{
			cout << "\nWARNING(" << ++warn << "): Altough semicolon was given, but no host index was specified for [ virtual_sites3 ] in line:" << endl;
			cout << "\t"<<buf << endl;
			cout << "\tThe host index, will be the index of the first atom (" << virtual_site[n_gr_virtuals - 1].indices[0] << ") building the virtual site" << endl;
			virtual_site[n_gr_virtuals - 1].host_index = virtual_site[n_gr_virtuals - 1].indices[0];
		}
	}
	else
	{
		cout << "\nWARNING(" << ++warn << "): No host index was specified for [ virtual_sites3 ] in line:" << endl;
		cout <<"\t"<< buf << endl;
		cout << "\tThe host index, will be the index of the first atom (" << virtual_site[n_gr_virtuals - 1].indices[0] << ") building the virtual site" << endl;
		virtual_site[n_gr_virtuals - 1].host_index = virtual_site[n_gr_virtuals - 1].indices[0];
	}
	
	virtual_site[n_gr_virtuals - 1].assigned = -1;
	virtual_site[n_gr_virtuals - 1].tnumb = 3;
	virtual_site[n_gr_virtuals - 1].indices[3] = -1;//not used

	delete[] name1;
	delete[] name2;
	delete[] name3;
	if (conv_numb != NULL)
		delete[] conv_numb;
	return 1;
}
//reading the virtual_sites4, the virtual site located on the line between the atoms from a fraction from the first atom 
int Topology::GetVirtual4(ifstream & file)
{
	int i;
	char *name1, *name2, *name3,*conv_numb=NULL;

	SetArraysize(&name1, NAME_SIZE, "name1", "Topology::GetVirtual4");
	SetArraysize(&name2, NAME_SIZE, "name2", "Topology::GetVirtual4");
	SetArraysize(&name3, NAME_SIZE, "name3", "Topology::GetVirtual4");

	//Create the array or resize it
	if (n_gr_virtuals == 0)
	{
		n_gr_virtuals++;
		SetArraysize(&virtual_site, n_gr_virtuals, "virtual_site", "Topology::GetVirtual4");
	}
	else
	{
		ResizeVirtualSiteStruct(&n_gr_virtuals, n_gr_virtuals + 1, &virtual_site, "virtual_site", "Topology::GetVirtual4");

	}

	mystrcpy(name1, NAME_SIZE, "virtual_site4[");
	IntToStr(&conv_numb, n_gr_virtuals);

	mystrcat(name1, NAME_SIZE, conv_numb);
	mystrcat(name1, NAME_SIZE, "] field: ");
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "index of the virtual site");

	virtual_site[n_gr_virtuals - 1].atom_index=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual4", currentfilename);//read the atom indices making up the sites
	nvirtuals_per_type[nmoltype - 1]++;
	cumul_virtuals[nmoltype]++;
	for (i = 0; i < 4; i++)
	{
		mystrcpy(name2, NAME_SIZE, name1);
		mystrcat(name2, NAME_SIZE, "index of the ");
		IntToStr(&conv_numb, i + 1);
		mystrcat(name2, NAME_SIZE, conv_numb);
		mystrcat(name2, NAME_SIZE, ". atom creating the virtual site");
		virtual_site[n_gr_virtuals - 1].indices[i]=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual4", currentfilename);//reading the indices


		if (virtual_site[n_gr_virtuals - 1].indices[i] > natoms_per_type[nmoltype - 1])
		{
			cout << "\n****ERROR****" << endl;
			cout << "In line: " << endl;
			cout << buf << endl;
			cout << "The " << i + 1 << ". atom index in the " << nvirtuals_per_type[nmoltype - 1] << ". 4-atoms virtual site of " << nmoltype << ". molecule type is " << endl;
			cout << virtual_site[n_gr_virtuals - 1].indices[i] << ", which is greater, than the number of atoms in the molecule (" << natoms_per_type[nmoltype - 1] << ")!" << endl;
			cout << "Check the file " << currentfilename << "!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}
	}
	mystrcpy(name2, NAME_SIZE, "[ virtual_sites4 ]");
	CheckRepetition(virtual_site[n_gr_virtuals - 1].indices, 4, nvirtuals_per_type[nmoltype - 1], nmoltype, name2);
	mystrcpy(name2, NAME_SIZE, name1);
	mystrcat(name2, NAME_SIZE, "function type");
	virtual_site[n_gr_virtuals - 1].type=ReadThisLine(file, 3, 0, name2, "Topology::GetVirtual4", currentfilename); //reading the function type

	if (virtual_site[n_gr_virtuals - 1].type == 1)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "In line: " << endl;
		cout << buf << endl;
		cout << "The unstable virtual site 4fd is not supported by the program, as it is not longer supported by the newer GROMACS versions." << endl;
		cout << "Use 4fdn type instead!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	else
	{
		if (virtual_site[n_gr_virtuals - 1].type != 2)
		{
			cout << "\n****ERROR****" << endl;
			cout << "In line: " << endl;
			cout << buf << endl;
			cout << "The function type in the " << nvirtuals_per_type[nmoltype - 1] << ". 4-atoms virtual site of " << nmoltype << ". molecule type is " << endl;
			cout << virtual_site[n_gr_virtuals - 1].type << ", and it can only be 2!" << endl;
			cout << "Check the file " << currentfilename << "!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}
	}
	
	for (i = 0; i < 3; i++)
	{
		mystrcpy(name2, NAME_SIZE, name1);
		switch (i)
		{
		case 0:
			mystrcat(name2, NAME_SIZE, "parameter a");
			break;

		case 1:
			mystrcat(name2, NAME_SIZE, "parameter b");
			break;
		default: 
			mystrcat(name2, NAME_SIZE, "parameter c");
			break;

		}
		
		virtual_site[n_gr_virtuals - 1].params[i]=ReadThisLine(file, 3, 0.0, name2, "Topology::GetVirtual4", currentfilename);// read the parameter a
		if (i == 2)
			virtual_site[n_gr_virtuals - 1].params[i] *= 10.0 / SimpleCfg::boxedge;//to have it in reduced unit
	}
	if (ExtractSemiColon(file, 0))//semicolon was give, try to read
	{
		virtual_site[n_gr_virtuals - 1].host_index = ReadThisLine(file, 2, -1, "host_index", "Topology::GetVirtual4", currentfilename);// read host atom index
		if (virtual_site[n_gr_virtuals - 1].host_index == -1)
		{
			cout << "\nWARNING(" << ++warn << "): Altough semicolon was given, but no host index was specified for [ virtual_sites4 ] in line:" << endl;
			cout << "\t"<<buf << endl;
			cout << "\tThe host index, will be the index of the first atom (" << virtual_site[n_gr_virtuals - 1].indices[0] << ") building the virtual site" << endl;
			virtual_site[n_gr_virtuals - 1].host_index = virtual_site[n_gr_virtuals - 1].indices[0];
		}
	}
	else
	{
		cout << "\nWARNING(" << ++warn << "): No host index was specified for [ virtual_sites4 ] in line:" << endl;
		cout << "\t"<<buf << endl;
		cout << "\tThe host index, will be the index of the first atom (" << virtual_site[n_gr_virtuals - 1].indices[0] << ") building the virtual site" << endl;
		virtual_site[n_gr_virtuals - 1].host_index = virtual_site[n_gr_virtuals - 1].indices[0];
	}
	virtual_site[n_gr_virtuals - 1].assigned = -1;
	virtual_site[n_gr_virtuals - 1].tnumb = 4;
	delete[] name1;
	delete[] name2;
	delete[] name3;
	if (conv_numb != NULL)
		delete[] conv_numb;
	return 1;
}
//Some reading and writing depends on the force field. Presently the following force fields can be used
//name				combination_rule		generate_pairs
//OPLSAA					3						yes
//Encads, Encadv			3						yes
//CHARMM*					2						yes
//AMBER						2						yes
//UNKNOWN_GENERATE_PAIRS	X						yes
//UNKNOWN_READ_PAIRS		X						no

//the last two is to be able to handle the others
int Topology::CheckForcefield(ifstream &file,char *name)
{
	int i,size;
	char *include_name, *pchar,*sep=NULL;
	int good_rule = 0;//no problem with combination rule
	SetArraysize(&include_name, LINE_SIZE, "include_name", "Topology::CheckForcefield");
	char ffname[NAME_SIZE], upper_name[NAME_SIZE];
	
	int result=0;
	file>>name;
	pchar = name + 1;//cut " from the beginning
	if (strncmp(name + 1, "ff", 2) == 0)
		pchar = name + 3;//old format, forcefiled name preceeded with ff
	sep = strchr(pchar, '.');
	
	SkipLine(file,currentfilename,1);
	for (i=0;i<N_FORCEF;i++)
	{
		size =(int)( sep - pchar);
		//the new format charmm and amber force fields are not listed separately, only check for the beginning of the name
		if (strncmp(force_field[i], "amber", 5) == 0)
			size= 5;
		if (strncmp(force_field[i], "charmm", 6) == 0)
			size = 6;
		
		for (int j = 0; j < size; j++)
		{
			ffname[j] = toupper(*(force_field[i] + j));
			upper_name[j] = toupper(pchar[j]);
		}
		if (strncmp(upper_name,ffname,size)==0)
		{
			result=1;
			SetArraysize(&force_field_name,NAME_SIZE,"force_field","Topology::CheckForcefield");
			mystrncpy(force_field_name, NAME_SIZE, pchar, sep-pchar);
			
			force_field_type=(T_force_field)i;
			
			//setting the 1-4 handling, and checking the combination rule
			//the supported force fields should use their own combination  rule
			switch (force_field_type)
			{
				case OPLSAA:
				case Encads:
				case Encadv: 
					if (RunParams::vdW_comb_rule != 3)
						good_rule = 3;
					gen_pairs = 1; 
					break;
				
				case CHARMM:
				case AMBER:
					if (RunParams::vdW_comb_rule != 2)
						good_rule = 2;
					gen_pairs = 1;
					break;

				case UNKNOWN_GENERATE_PAIRS:
					gen_pairs = 1;
					break;
				case UNKNOWN_READ_PAIRS:
					gen_pairs = 0;
					break;
			}
			if (good_rule>0 && (RunParams::potential>0 && RunParams::potential<10))
			{
				cout << force_field_name << " forcefield should use combination rule "<<good_rule<<", change it in the " << datfilename << "!" << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}
			break;
		}
		
	}
	if (result==0)
	{
		//determinig, whether it is an itp file
		if (strncmp(name+strlen(name)-5,".itp",4)==0)
			result =2;
	}
	return result;
};
//Checking, whether the option defined by the ifdef or ifndef was passed to the topology file
//flag =1 for ifdef, 2 for ifndef
int Topology::CheckPreprocDir(ifstream &file,int nfiles, int flag)
{
	int i,result=0;
	char *name;

	if (!(flag == 1 || flag == 2))
	{
		cout << "Invalid function call flag argument " << flag << " for Topology::CheckPreprocDir, exiting..." << endl;
		CleanExit();
	}
	nifdef_ifndef[nfiles-1]++;//increasing the number of ifdefs and ifndefs for this file
	if (nifdef_ifndef[nfiles - 1] > N_ACT_IFDEF)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The number of nestled ifdef - ifndef preprocessor directives is larger than the maximum (N_ACT_IFDEF=" << N_ACT_IFDEF << ")";
		cout << "in topology file " << currentfilename << "!" << endl;
		cout << "Either reorganize the topology file, or compile the code after increasing the N_ACT_IFDEF value in units.h!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}

	SetArraysize(&name,NAME_SIZE+2,"name","Topology::CheckPreprocDir");
	file>>name;
	for (i=0;i<RunParams::n_top_def_option;i++)
	{
		if (strncmp(name,RunParams::top_def_option[i]+2,strlen(RunParams::top_def_option[i])-2)==0 && (strlen(RunParams::top_def_option[i]) - 2)==strlen(name))
		{
			pifdef_ifndef[(nfiles-1)*N_ACT_IFDEF+nifdef_ifndef[nfiles-1]-1]=flag;
			if (::debug)
				cout<<name;
			result = 1;
			return result;
		}
	}
	pifdef_ifndef[(nfiles - 1) * N_ACT_IFDEF + nifdef_ifndef[nfiles - 1] - 1] = -flag;
	return result;
};

//Decide, whether this line should be included according to the active peprocessor directives 
bool Topology::IncludeLine(int nfiles)
{
	int i,offset;
	bool result=true;
	offset=(nfiles-1)*N_ACT_IFDEF;//offset in the array
	//going through the active ifdef or ifndef directives 
	for (i=0;i<nifdef_ifndef[nfiles-1];i++)
	{
		if (pifdef_ifndef[offset + i] == 1 && pelse[offset + i] == 1)//already in the else part of a valid ifdef define statement, no processing
		{
			result = false;
			break;//no need to search futher
		}
		if (pifdef_ifndef[offset+i]==-1 && pelse[offset+i]==0)//in the ifdef part of an invalid ifdef statement, no processing
		{
			result = false;
			break;//no need to search futher
		}
		if (pifdef_ifndef[offset + i] == 2 && pelse[offset + i] == 0)//in the ifndef part of a valid ifndef statement, no processing
		{
			result = false;
			break;//no need to search futher
		}
		if (pifdef_ifndef[offset + i] == -2 && pelse[offset + i] == 1)//already in the else part of an invalid ifndef statement, no processing
		{
			result = false;
			break;//no need to search futher
		}
	}
	return result;
};

//Extracting the text qualifier
int Topology::ExtractSemiColon(ifstream &file,int mode)
{
	int code;
	std::streamoff pos;

	pos=file.tellg();
	do
	{
		code=file.get();
		if (code==EOF || code ==(int)'\n')
		{
			if (mode)
				cout << "\n*****ERROR*****" << endl;
			else
				cout << "\nWARNING(" << ++warn << "): ";
			if (code==EOF )
				cout<<"End of file "<<currentfilename<<" was reached during text qualifier extraction!"<<endl; 
			if (code == (int)'\n')
			{
				cout << "End of line " << endl; 
				cout <<"\t"<< buf << endl;
				cout<<"\twas reached in file "<<currentfilename<<" during text qualifier extraction!"<<endl; 
				cout << "\tThe line was processed with the [ " << directives[last_directive] << " ] directive" << endl;
			}
			if (mode)
			{
				cout<<"There is not enough data, terminating program..."<<endl;
				CleanExit();
			}
			else
			{
				//semicilon was not found, put the file pointer back, where it was and return
				file.seekg(pos,ios::beg);
				return 0;
			}
		}
		
	}while (code!=(int)';');
	return(1);
};
//reset the array elements, if the file is closed
void Topology::Reset(int nfiles)
{
	int i;
	nifdef_ifndef[nfiles-1]=0;
	//zeroing the element of this files segment
	for (i=(nfiles-1)*N_ACT_IFDEF;i<nfiles*N_ACT_IFDEF;i++)
	{
		pifdef_ifndef[i]=0;
		pelse[i]=0;
	}
};

//writing the parameters
void Topology::PrintParams(ostream &file)
{
	int i,j,k,ind=0,iexcl;
	int sum_dih, cumul_dih;
	int prev_charge_gr=0;//index of the last charge group in the previous molecule type

	file<<"\nMolecule binding is created based on the following topology data:"<<endl;
	file << "Force field name: " << force_field_name << endl;
	file<<"Number of different molecule types:\t"<<nmoltype<<endl;
	file<<"Number of different GROMACS types: \t"<<nGRtypes<<endl;
	for (i=0;i<nmoltype;i++)
		file<<J10<<nmol_per_type[i]<<J10<<molname[i]<<J20<<" molecules"<<endl;

	for (i=0;i<nmoltype;i++)
	{
		file<<"\nParameters for molecule "<<molname[i]<<endl;
		file << "The charge groups are only meaningful inside a molecule, and the numbering can restart for each molecule type!" << endl;
		file<<J10<<" "<<J10<<" "<<J10<<" "<<J10<<" "<<J16<<"  "<<J20<<"RMC index"<<endl;
		file<<J10<<"index"<<J10<<"Atom name"<<J10<<"Charge"<<J16<<"Charge group"<<J10<<"1.mol"<<J10<<"2.mol"<<J13<<"Delta index"<<J20<<"GROMACS-type ind"<<endl;
		file.precision(3);
		file.setf(ios::floatfield, ios::fixed);
		for (j=0;j<natoms_per_type[i];j++)
		{
			ind=cumul_atoms[i]+j;
			file<<J10<<j+1<<J10<<atom_name[ind]<<J10<<charge[charge_type[ind]]<<J16<<charge_group[ind]+1-prev_charge_gr<<J10<<first_index[ind]<<J10<<second_index[ind]<<J13<<delta_index[ind]<<J20<<GROMACS_type[ind]+1<<endl;
		}
		prev_charge_gr=charge_group[ind]+1;//to start the charge group index with 0 for each molecule type
		file.precision(4);
		file.unsetf(ios::floatfield);
		if (npairs_per_type[i]>0)
		{
			file<<"\t"<<"\n1-4 pairs "<<endl;
			
			if (RunParams::potential==1)//LJ
			{
				if (Topology::gen_pairs==0)//for OPLSAA, encads... parameters will be calculated by scaling, not given here
				{
					file << J10 << "i" << J10 << "j" << J20;
					if (RunParams::vdW_comb_rule == 0 || RunParams::vdW_comb_rule==1)//C6, C12
						file<<"C6 (kJ/mol*A^6)"<<J20<<"C12 (kJ/mol*A^12)";
					else
						file<<"sigma (A)"<<J20<<"epsilon (kJ/mol)";//LJsigma, epsilon
				}
				else 
				{
					file<<"1-4 parameters are calculated by scaling."<<endl;
					file<<"Fudge vDW:      "<<RunParams::vdW14_fudge<<endl;
					file<<"Fudge Coulomb:  "<<RunParams::Coulomb14_fudge<<endl;
					file << J10 << "i" << J10 << "j" << J20;
				}
				
				
			}
			/*else
			{
				if (RunParams::potential==2)//Buckingham, not implemented yet
					file<<"A (kJ/mol)"<<J20<<"B (1/A)"<<J20<<"C (kJ*A^6/mol)"<<endl;
			}*/
			file<<endl;
		
		}
		for (j=0;j<npairs_per_type[i];j++)
		{
			file<<J10<<pair_index[(cumul_pairs[i]+j)*2]<<J10<<pair_index[(cumul_pairs[i]+j)*2+1];
			if (Topology::gen_pairs==0)
			{
				file.precision(4);
				file.setf(ios::fixed, ios::floatfield);
				file.setf(ios::right, ios::adjustfield);
				file<<J20<<vdW14_pot1[pair_type[cumul_pairs[i]+j]]<<J20<<vdW14_pot2[bond_type[cumul_pairs[i]+j]];
				file.unsetf(ios::fixed);
				file.unsetf(ios::right);
			}
			file<<endl;		
			file.precision(6);
		}
		if (nexclusionlines_per_type[i] > 0)
		{
			file << "\t" << "\nExclusions only defined by [ exclusions ] directive " << endl;
			file << J10 << "central" << J10 << "indices" << endl;
		}
		for (j = 0; j < nexclusionlines_per_type[i]; j++)
		{
			if (exclusion_centre[cumul_exclusionlines[i] + j] > -1)//could have been removed
			{
				file << J10 << exclusion_centre[cumul_exclusionlines[i] + j];
				iexcl = cumul_nexcluded_neighs[cumul_exclusionlines[i] + j];

				for (k = 0; k < nexcluded_neighs[cumul_exclusionlines[i] + j]; k++)
				{
					file << J10 << exclusions[iexcl];
					iexcl++;
				}
				file << endl;
			}
			
		}
		if (nvirtuals_per_type[i] > 0)
		{
			file << "\t" << "\nVirtual sites " << endl;
			
		}
		for (j = 0; j < nvirtuals_per_type[i]; j++)
		{
			file << J20 << "virtual site type: ";
			switch (virtual_site[cumul_virtuals[i] + j].tnumb)
			{
			case 2:
				file << J7 << 2;
				break;
			case 3:
				switch (virtual_site[cumul_virtuals[i] + j].type)
				{
				case 1:
					file << J7 << 3;
					break;
				case 2:
					file << J7 << "3fd";
					break;
				case 3:
					file << J7 << "3fad";
					break;
				case 4:
					file << J7 << "3out";
					break;

				}
				break;
			case 4:
				file << J7 << 4;
				break;
			}
			file << " for virtual atom index " << J7 << virtual_site[cumul_virtuals[i] + j].atom_index << endl;
			file << "\tbuilt from atoms: ";
			for (k = 0; k < virtual_site[cumul_virtuals[i] + j].tnumb; k++)
				file << J7 << virtual_site[cumul_virtuals[i] + j].indices[k];
			file << endl;
			file << "\twith parameters:" << endl;
			switch (virtual_site[cumul_virtuals[i] + j].tnumb)
			{
			case 2:
				file << "\t\ta = " << virtual_site[cumul_virtuals[i] + j].params[0] << endl;;
				break;
			case 3:
				switch (virtual_site[cumul_virtuals[i] + j].type)
				{
				case 1:
					file << J7 << "\t\ta = " << virtual_site[cumul_virtuals[i] + j].params[0]<<endl;
					file << J7 << "\t\tb = " << virtual_site[cumul_virtuals[i] + j].params[1] << endl;
					break;
				case 2:
					file << J7 << "\t\ta = " << virtual_site[cumul_virtuals[i] + j].params[0] << endl;
					file << J7 << "\t\td = " << virtual_site[cumul_virtuals[i] + j].params[1]/10.0 * SimpleCfg::boxedge <<" nm" << endl;
					break;
				case 3:
					file << J7 << "\t\ttheta = " << virtual_site[cumul_virtuals[i] + j].params[0] << " deg"<<endl;
					file << J7 << "\t\td = " << virtual_site[cumul_virtuals[i] + j].params[1]/10.0 * SimpleCfg::boxedge << " nm" << endl;
					break;
				case 4:
					file << J7 << "\t\ta = " << virtual_site[cumul_virtuals[i] + j].params[0] << endl;
					file << J7 << "\t\tb = " << virtual_site[cumul_virtuals[i] + j].params[1] << endl;
					file << J7 << "\t\tc = " << virtual_site[cumul_virtuals[i] + j].params[2]*10.0 / SimpleCfg::boxedge << " nm-1"<<endl;
					break;

				}
				break;
			case 4:
				file << J7 << "\t\ta = " << virtual_site[cumul_virtuals[i] + j].params[0] << endl;
				file << J7 << "\t\tb = " << virtual_site[cumul_virtuals[i] + j].params[1] << endl;
				file << J7 << "\t\tc = " << virtual_site[cumul_virtuals[i] + j].params[2]/10.0 * SimpleCfg::boxedge << " nm" << endl;
				break;
			}
			file << "\tHost atom type for the potential partials: "<< SimpleCfg::virtual_host_type[cumul_virtuals[i] + j]+1<< endl;

		}

		if (nbonds_per_type[i]>0)
		{
			file<<"\t"<<"\nBonds between atoms "<<endl;
			file<<J10<<"i"<<J10<<"j"<<J20<<"r (A)"<<J20<<"k (kJ/mol/A2)"<< J10 << "sigma" << endl;
		}
		for (j=0;j<nbonds_per_type[i];j++)
		{
			file<<J10<<bond_index[(cumul_bonds[i]+j)*2]<<J10<<bond_index[(cumul_bonds[i]+j)*2+1];
			file.precision(4);
			file.setf(ios::fixed, ios::floatfield);
			file.setf(ios::right, ios::adjustfield);
			file<<J20<<bond_r[bond_type[cumul_bonds[i]+j]]<<J20<<bond_k[bond_type[cumul_bonds[i]+j]];
			file.setf(ios::scientific, ios::floatfield);
			file<<J20<<bond_sigma[bond_type[cumul_bonds[i]+j]]<<endl;
			file.precision(6);
			file.unsetf(ios::scientific);
			file.unsetf(ios::right);
		}
		if (nangles_per_type[i]>0)	

		{
			file<<"\t"<<"\nAngles between atoms"<<endl;
			file<<J10<<"i"<<J10<<"j"<<J10<<"k"<<J20<<"angle (degree)"<<J20<<"k (kJ/mol/rad2)"<< J10 << "sigma" << endl;
		}
		for (j=0;j<nangles_per_type[i];j++)
		{
			file<<J10<<angle_index[(cumul_angles[i]+j)*3]<<J10<<angle_index[(cumul_angles[i]+j)*3+1]<<J10<<angle_index[(cumul_angles[i]+j)*3+2];
			file.precision(4);
			file.setf(ios::fixed, ios::floatfield);
			file.setf(ios::right, ios::adjustfield);
			file<<J20<<angle_ang[angle_type[cumul_angles[i]+j]]<<J20<<angle_k[angle_type[cumul_angles[i]+j]];
			file.setf(ios::scientific, ios::floatfield);
			file<<J20<<angle_sigma[angle_type[cumul_angles[i]+j]]<<endl;
			file.precision(6);
			file.unsetf(ios::scientific);
			file.unsetf(ios::right);
		}
		//The indices of  all the dihedrals are kept in one array, as it will make the craetion of one index list for FNCD_MOL easier
		if (nperdihedrals_per_type[i]>0)
		{
			if (force_field_type==OPLSAA)
				file<<"\t"<<"\nImproper (periodic) dihedrals between atoms"<<endl;

			else
				file<<"\t"<<"\nProper (periodic) dihedrals between atoms"<<endl;
		file<<J10<<"i"<<J10<<"j"<<J10<<"k"<<J10<<"l"<<J10<<"phi (degree)"<<J10<<"k (kJ/mol)"<<J10<<"multiplicity"<<J10<<"sigma"<<endl;	
		}
		k=0;
		sum_dih=nperdihedrals_per_type[i]+nharmdihedrals_per_type[i]+nRBdihedrals_per_type[i];
		cumul_dih=cumul_perdihedrals[i]+cumul_harmdihedrals[i]+cumul_RBdihedrals[i];
		for (j=0;j<sum_dih;j++)
		{
			if (dihedral_type[cumul_dih+j]==1)
			{
				file<<J10<<dihedral_index[(cumul_dih+j)*4]<<J10<<dihedral_index[(cumul_dih+j)*4+1]<<J10<<dihedral_index[(cumul_dih+j)*4+2]<<J10<<dihedral_index[(cumul_dih+j)*4+3];
				file.precision(4);
				file.setf(ios::fixed, ios::floatfield);
				file.setf(ios::right, ios::adjustfield);
				file<<J10<<perdihedral_ang[perdihedral_type[cumul_perdihedrals[i]+k]]<<J10<<perdihedral_k[perdihedral_type[cumul_perdihedrals[i]+k]];
				file<<J10<<perdihedral_multiplicity[perdihedral_type[cumul_perdihedrals[i]+k]];
				file.setf(ios::scientific, ios::floatfield);
				file<<J20<<perdihedral_sigma[perdihedral_type[cumul_perdihedrals[i]+k]]<<endl;
				file.precision(6);
				file.unsetf(ios::scientific);
				file.unsetf(ios::right);
				k++;
			}
		}
		if (nharmdihedrals_per_type[i]>0)
		{
			file<<"\t"<<"\nImproper (harmonic) dihedrals between atoms"<<endl;
			file<<J10<<"i"<<J10<<"j"<<J10<<"k"<<J10<<"l"<<J10<<"phi (degree)"<<J10<<"k (kJ/mol)"<<J10<<"sigma"<<endl;	
		}
		k=0;
		sum_dih=nperdihedrals_per_type[i]+nharmdihedrals_per_type[i]+nRBdihedrals_per_type[i];
		cumul_dih=cumul_perdihedrals[i]+cumul_harmdihedrals[i]+cumul_RBdihedrals[i];
		for (j=0;j<sum_dih;j++)
		{
			if (dihedral_type[cumul_dih+j]==2)
			{
				file<<J10<<dihedral_index[(cumul_dih+j)*4]<<J10<<dihedral_index[(cumul_dih+j)*4+1]<<J10<<dihedral_index[(cumul_dih+j)*4+2]<<J10<<dihedral_index[(cumul_dih+j)*4+3];
				file.precision(4);
				file.setf(ios::fixed, ios::floatfield);
				file.setf(ios::right, ios::adjustfield);
				file<<J10<<harmdihedral_ang[harmdihedral_type[cumul_harmdihedrals[i]+k]]<<J10<<harmdihedral_k[harmdihedral_type[cumul_harmdihedrals[i]+k]];
				file.setf(ios::scientific, ios::floatfield);
				file<<J20<<harmdihedral_sigma[harmdihedral_type[cumul_harmdihedrals[i]+k]]<<endl;
				file.precision(6);
				file.unsetf(ios::scientific);
				file.unsetf(ios::right);
				k++;
			}
		}
		if (nRBdihedrals_per_type[i]>0)
		{
			if (force_field_type==OPLSAA)
				file<<"\t"<<"\nProper (Ryckaert-Bellemans) dihedrals between atoms"<<endl;
			else
				file<<"\t"<<"\nRyckaert-Bellemans dihedrals between atoms"<<endl;
			file<<J10<<"i"<<J10<<"j"<<J10<<"k"<<J10<<"l"<<J10<<"C0"<<J10<<"C1"<<J10<<"C2"<<J10<<"C3"<<J10<<"C4"<<J10<<"C5"<<" (kJ/mol)"<<J10<<"sigma"<<endl;	
		}
		k=0;
		sum_dih=nperdihedrals_per_type[i]+nharmdihedrals_per_type[i]+nRBdihedrals_per_type[i];
		cumul_dih=cumul_perdihedrals[i]+cumul_harmdihedrals[i]+cumul_RBdihedrals[i];
		for (j=0;j<sum_dih;j++)
		{
			if (dihedral_type[cumul_dih+j]==3)
			{
				file<<J10<<dihedral_index[(cumul_dih+j)*4]<<J10<<dihedral_index[(cumul_dih+j)*4+1]<<J10<<dihedral_index[(cumul_dih+j)*4+2]<<J10<<dihedral_index[(cumul_dih+j)*4+3];
				file.precision(4);
				file.setf(ios::fixed, ios::floatfield);
				file.setf(ios::right, ios::adjustfield);
				file<<J10<<RBdihedral_C0[RBdihedral_type[cumul_RBdihedrals[i]+k]]<<J10<<RBdihedral_C1[RBdihedral_type[cumul_RBdihedrals[i]+k]];
				file<<J10<<RBdihedral_C2[RBdihedral_type[cumul_RBdihedrals[i]+k]]<<J10<<RBdihedral_C3[RBdihedral_type[cumul_RBdihedrals[i]+k]];
				file<<J10<<RBdihedral_C4[RBdihedral_type[cumul_RBdihedrals[i]+k]]<<J10<<RBdihedral_C5[RBdihedral_type[cumul_RBdihedrals[i]+k]];
				file.setf(ios::scientific, ios::floatfield);
				file<<J20<<RBdihedral_sigma[RBdihedral_type[cumul_RBdihedrals[i]+k]]<<endl;
				file.precision(6);
				file.unsetf(ios::scientific);
				file.unsetf(ios::right);
				k++;
			}
		}
		
		
	}

};
//determining the number of interactions the atoms are involved in for the given interaction type 
void Topology::CreateList(int pottype, FNC_POT &pfnc)
{
	int imoltype,imol,iint,iat,ivirt;
	int charge_c_offset=0;
	int ind1,ind2=0,ind3,ind4,iper,iharm,iRB;
	int this_charge_group,my_ind;
	int prev_charge_group=0;
	int *pind=NULL, *pper_type=NULL,*ppertype=NULL,*pharmtype=NULL,*pRBtype=NULL,*count=NULL,*conv=NULL,*ptype=NULL,*pdihtype=NULL;
	int natoms_involved=0,sign=0;

	switch (pottype)
	{
		case BOND:
		{
			pind=bond_index;
			pper_type=nbonds_per_type;
			natoms_involved=2;
			break;
		}
		case ANGLE:
		{
			pind=angle_index;
			pper_type=nangles_per_type;
			natoms_involved=3;
			break;
		}
		case DIHEDRAL:
		{
			pind=dihedral_index;
			pper_type=ndihedrals_per_type;
			natoms_involved=4;
			break;
		}
	}

	
	
	if (RunParams::fnc==4)//only if flexible molecules present
	{
		//going through all the molecule types, and counting the neighbours for each atom
		for (imoltype=0;imoltype<nmoltype;imoltype++)
		{
			//going through all the interactions of the given type for this molecule type
			for (iint=0;iint<pper_type[imoltype];iint++)
			{
				for (iat=0;iat<natoms_involved;iat++)
				{
					for (imol=0;imol<nmol_per_type[imoltype];imol++)
					{
						//*pind is the GROMACS index of the atom in the molecule, it start with 1, so 1 has to be subtracted
						//first_index[cumul_atom[imoltype]+*pind]is the RMC_index of the first occurence of this atomtype, if this molecule
						//	type would be the first in the system
						//RMCindex_offset[imoltype] is the offset of the RMC index to take into account that might not be the first moleculetype
						//imol*delta_index[cumul_atom[imoltype]+*pind] is the offset of this occurence's index relative to the first
						pfnc.numbneigh[first_index[cumul_atoms[imoltype]+*pind-1]-1+RMCindex_offset[imoltype]+imol*delta_index[cumul_atoms[imoltype]+*pind-1]]++;
//	cout<<"iat  "<<iat<<"  "<<imol<<"  "<<first_index[cumul_atoms[imoltype]+*pind]+imol*delta_index[cumul_atoms[imoltype]+*pind]<<endl;
				
					}
					pind++;
				}
				
			}
		}
		//initializing the converter to find the first array element in the neighbour and type array for each atom
		pfnc.converter[0]=0;
		pfnc.ndistances=pfnc.numbneigh[0];//counting all the neighbours, needed for array dimensions
		for (iat=1;iat<pfnc.natoms;iat++)
		{
			pfnc.converter[iat]=pfnc.converter[iat-1]+pfnc.numbneigh[iat-1];
			pfnc.ndistances+=pfnc.numbneigh[iat];
		}
		//resizing the arrays to the right size
		pfnc.ResizeNeighArrays();
		conv=pfnc.converter;//to make it shorter
	}
	//creating the index lists
	//BONDS: 1---2 (both atoms are 'end' aroms) (no special order!)
	//ANGLE: 1---2---3 (1 and 3 are 'end' atoms, 2 is 'middle')	
	//DIHEDRAL: 1---2---3---4 (1 and 4 are 'end' atoms, 2 and 3 are 'middle')
	//the sign of the index in the neighbours array will show, whether the central atom is an end (-) or a middle (+) one
	//so the order of the atoms can be decided
	//the indices in the arrays will be the following
	//ANGLE:    central atom		 1	 2	 3		DIHEDRAL:	central atom		 1	 2	 3	 4
	//			array neighbour		-2	 1	-2					array neighbour		-2	 1	 4	-3
	//			array neighbour2	 3	 3	 1					array neighbour2	 3	 3	 2	 2
	//															array neighbour3	 4	 4	 1	 1

	switch (pottype)//resetting the pointer to the beginnig of the right array
	{
		case BOND:
		{
			pind=bond_index;
			ptype=bond_type;
			sign=1;
			
			//setting up the list with the charge type indices, it take less space, than storing the charges for each atom
			//setting the charge group and the relative position of the given atom in the charge group
			if (RunParams::potential==1)//Coulomb has to be calculated
			{
				this_charge_group=0;//charge groups in the charge group arry are numbered consecutively regardless the molecule type starting with 0
				ivirt = 0;
				for (imoltype=0;imoltype<nmoltype;imoltype++)
				{
					my_ind=-1;
					for (iat=0;iat<natoms_per_type[imoltype];iat++)
					{
						//iat is the GROMACS index of the atom in the molecule, it start with 0
						//first_index[cumul_atom[imoltype]+iat] is the RMC_index of the first occurence of this atomtype, if this molecule
						//	type would be the first in the system
						//RMCindex_offset[imoltype] is the offset of the RMC index to take onto account that might not be the first moleculetype
						//imol*delta_index[cumul_atom[imoltype]+iat] is the offset of the imol occurence's index relative to the first
						
						//RMC_index of the atom in the first instance, has to be calculated differently for virtial sites
						if (first_index[cumul_atoms[imoltype] + iat] == 0)
						{
							ind1 = SimpleCfg::cumul[SimpleCfg::ntypes + ivirt];
							ivirt++;
						}
						else
							ind1=first_index[cumul_atoms[imoltype]+iat]-1+RMCindex_offset[imoltype];
						if (RunParams::fnc==4)
						{
							if (this_charge_group==charge_group[cumul_atoms[imoltype]+iat])
								my_ind++;
							else
							{
								this_charge_group=charge_group[cumul_atoms[imoltype]+iat];
								my_ind=0;
							}
						}
						for (imol=0;imol<nmol_per_type[imoltype];imol++)
						{
							ind2=ind1+imol*delta_index[cumul_atoms[imoltype]+iat];
							pfnc.charge_type[ind2]=charge_type[cumul_atoms[imoltype]+iat];
							if (RunParams::fnc==4)
							{
								pfnc.charge_group[ind2]=charge_group[cumul_atoms[imoltype]+iat];
								pfnc.charge_centre[ind2]=charge_c_offset+imol*ncharge_group_per_type[imoltype]+charge_group[cumul_atoms[imoltype]+iat]-prev_charge_group;
							}
						}//end of imol cycle
					}//end of iat cycle
					if (RunParams::fnc==4)
					{
						charge_c_offset=pfnc.charge_centre[ind2]+1;//to continue with the next available charge centre 
						prev_charge_group+=ncharge_group_per_type[imoltype];//to restart the charge group index in the next molecule type with 0
					}
				}//end of imoltype cycle
				
			}//end of if there is a potential
			break;
		}

		case ANGLE:
		{
			pind=angle_index;
			ptype=angle_type;
			sign=-1;//to multiply the index in the neighbours array for central 1 to show whether the central atom is an 'end' one  
			break;
		}
		case DIHEDRAL:
		{
			iper=0;
			iharm=0;
			iRB=0;
			pind=dihedral_index;
			pdihtype=dihedral_type;
			//ptype will be assigned later, when we know which is the actual GROMACS dihedral type
			ppertype=perdihedral_type;
			pharmtype=harmdihedral_type;
			pRBtype=RBdihedral_type;
			sign=-1;//to multiply the index in the neighbours array for central 1 and 4 to show whether the atom is an 'end' one  
			break;
		}
	}
	SetArraysize(&count,pfnc.natoms,"count","Topology::CreateList");//temporary array to count the already found neighbours
	for (iat=0;iat<pfnc.natoms;iat++)
		count[iat]=0;
	
	if (RunParams::fnc==4)//only if flexible molecules present
	{
		for (imoltype=0;imoltype<nmoltype;imoltype++)
		{
			//going through all the interactions of the given type for this molecule type
			for (iint=0;iint<pper_type[imoltype];iint++)
			{
				if (pottype==DIHEDRAL)
				{
					switch (*pdihtype)
					{
						case 1:
						{

							ptype=ppertype+iper;//to set it to the next interaction's type
							iper++;//advance the interaction type index for this interaction
							break;
						}
						case 2:
						{
							ptype=pharmtype+iharm;//to set it to the next interaction's type
							iharm++;//advance the interaction type index for this interaction
							break;
						}
						case 3:
						{
							ptype=pRBtype+iRB;//to set it to the next interaction's type
							iRB++;//advance the interaction type index for this interaction
							break;
						}
					}
				}
				for (imol=0;imol<nmol_per_type[imoltype];imol++)
				{
					//*pind is the GROMACS index of the atom in the molecule, it start with 1, so 1 has to be subtracted
					//first_index[cumul_atom[imoltype]+*pind] is the RMC_index of the first occurence of this atomtype, if this molecule
					//	type would be the first in the system
					//RMCindex_offset[imoltype] is the offset of the RMC index to take onto account that might not be the first moleculetype
					//imol*delta_index[cumul_atom[imoltype]+*pind] is the offset of this occurence's index relative to the first
					//the first atom is the central,putting the index of the second atom into its neighbours array
					//RMC_index of the first atom
					ind1=first_index[cumul_atoms[imoltype]+*pind-1]+RMCindex_offset[imoltype]+imol*delta_index[cumul_atoms[imoltype]+*pind-1];
					CheckAtomIndex(*pind-1,imol,imoltype,ind1);//checking, whether the index is smaller, than the total number of atoms
					//RMC_index of the second atom
					ind2=first_index[cumul_atoms[imoltype]+*(pind+1)-1]+RMCindex_offset[imoltype]+imol*delta_index[cumul_atoms[imoltype]+*(pind+1)-1];
					CheckAtomIndex(*(pind+1)-1,imol,imoltype,ind2);//checking, whether the index is smaller, than the total number of atoms
					//first atom is the central, putting the index of the second atom into its neighbours and types array
					pfnc.neighbours[conv[ind1-1]+count[ind1-1]]= sign*ind2;
					pfnc.consttypes[conv[ind1-1]+count[ind1-1]]=*ptype;
						
					//second atom is the central, putting the index of the first atom into its neighbours and types array
					pfnc.neighbours[conv[ind2-1]+count[ind2-1]]=ind1;
					pfnc.consttypes[conv[ind2-1]+count[ind2-1]]=*ptype;
					switch (pottype)//resetting the pointer to the beginnig of the right array
					{
						case DIHEDRAL://there is a 3rd and 4th atom
						{
							//RMC_index of the 3rd atom
							ind3=first_index[cumul_atoms[imoltype]+*(pind+2)-1]+RMCindex_offset[imoltype]+imol*delta_index[cumul_atoms[imoltype]+*(pind+2)-1];
							CheckAtomIndex(*(pind+2)-1,imol,imoltype,ind3);//checking, whether the index is smaller, than the total number of atoms
							//RMC_index of the 4th atom
							ind4=first_index[cumul_atoms[imoltype]+*(pind+3)-1]+RMCindex_offset[imoltype]+imol*delta_index[cumul_atoms[imoltype]+*(pind+3)-1];
							CheckAtomIndex(*(pind+3)-1,imol,imoltype,ind4);//checking, whether the index is smaller, than the total number of atoms
							
							//first atom is the central, putting the index of the 4th atom into its neighbours3 array
							pfnc.neighbours3[conv[ind1-1]+count[ind1-1]]= ind4;

							//second atom is the central, putting the index of the 4th atom into its neighbours3 array
							pfnc.neighbours3[conv[ind2-1]+count[ind2-1]]= ind4;
						
							//third atom is the central, putting the index of the 4th atom into its neighbours array
							//this has to be done here, as it is not the same as for the ANGLE
							pfnc.neighbours[conv[ind3-1]+count[ind3-1]]= ind4;

							//third atom is the central, putting the index of the 2nd atom into its neighbours2 array
							//this has to be done here, as it is not the same as for the ANGLE
							pfnc.neighbours2[conv[ind3-1]+count[ind3-1]]= ind2;

							//third atom is the central, putting the index of the 1rst atom into its neighbours3 array
							pfnc.neighbours3[conv[ind3-1]+count[ind3-1]]= ind1;

							//fourth atom is the central, putting the index of the 3rd atom into its neighbours and type array
							pfnc.neighbours[conv[ind4-1]+count[ind4-1]]= sign*ind3;
						
							pfnc.consttypes[conv[ind4-1]+count[ind4-1]]=*ptype;

							//fourth atom is the central, putting the index of the 2nd atom into its neighbours2 array
							pfnc.neighbours2[conv[ind4-1]+count[ind4-1]]= ind2;

							//fourth atom is the central, putting the index of the 1st atom into its neighbours3 array
							pfnc.neighbours3[conv[ind4-1]+count[ind4-1]]= ind1;

							//dih_type array
							pfnc.dih_type[conv[ind1-1]+count[ind1-1]]= *pdihtype;
							pfnc.dih_type[conv[ind2-1]+count[ind2-1]]= *pdihtype;
							pfnc.dih_type[conv[ind3-1]+count[ind3-1]]= *pdihtype;
							pfnc.dih_type[conv[ind4-1]+count[ind4-1]]= *pdihtype;

							count[ind4-1]++;
						
						}
						case ANGLE://there is a 3rd atom
						{
							//RMC_index of the 3rd atom
							ind3=first_index[cumul_atoms[imoltype]+*(pind+2)-1]+RMCindex_offset[imoltype]+imol*delta_index[cumul_atoms[imoltype]+*(pind+2)-1];
							CheckAtomIndex(*(pind+2)-1,imol,imoltype,ind3);//checking, whether the index is smaller, than the total number of atoms
							//first atom is the central, putting the index of the 3rd atom into its neighbours2 array
							pfnc.neighbours2[conv[ind1-1]+count[ind1-1]]= ind3;

							//second atom is the central, putting the index of the 3rd atom into its neighbours2 array
							pfnc.neighbours2[conv[ind2-1]+count[ind2-1]]=ind3;
							
							if (pottype==ANGLE)//this is different for the DIHEDRAL, it was already set for it
							{
								//third atom is the central, putting the index of the second atom into its neighbours array
								pfnc.neighbours[conv[ind3-1]+count[ind3-1]]= -ind2;
							
								//third atom is the central, putting the index of the first atom into its neighbours2 array
								pfnc.neighbours2[conv[ind3-1]+count[ind3-1]]= ind1;
							}
							//third atom is the central, putting the type into its consttypes array
							pfnc.consttypes[conv[ind3-1]+count[ind3-1]]=*ptype;
							count[ind3-1]++;
							break;
						}
					}
					count[ind1-1]++;
					count[ind2-1]++;
										
				}//end of imol cycle
				if (pottype==DIHEDRAL)
					pdihtype++;
				pind+=natoms_involved;
				ptype++;
				
				
			}//end of iint cycle
		}//end of imoltype cycle
	}//end if
	delete [] count;
};

//checking, whether the index is smaller, than the total number of atoms
void Topology::CheckAtomIndex(int i, int imol, int imoltype, int iatom)
{
	if (iatom>SimpleCfg::ntotal)
	{
		cout << "\n****ERROR****"<<endl;
		cout<<"The RMC index of the "<<i+1<<". atom of the "<<imoltype+1<<". moleculetype's "<<imol+1<<". molecule "<<endl;
		cout<<"is "<<iatom<<", which is greater, than the total number of atoms ("<<SimpleCfg::ntotal<<")! Check the topology file!"<<endl;
		cout<<"Cannot run this way,exiting..."<<endl;
		CleanExit();
	}
};
//reading the sigma value for the interaction
double Topology::GetSigma(ifstream &file, const char *varname, const char *routine_name, int mode)
{
	double temp;
	temp=ReadThisLine(file,4,-100000.0,varname,routine_name);
	if (fabs(temp-(-100000.0))<LOAD_TOL)//default value is given back, there was no data
	{
		ExtractSemiColon(file,0);//maybe it was given after the text qualifier
		temp=ReadThisLine(file,mode,-100000.0,varname,routine_name);
	}
	return temp;
};

//check, whether an atom will have the same molecule type in all the instances of a molecule
void Topology::CheckRMCType()
{
	int i,imoltype,iat,ind1,ind2,RMC_type1,RMC_type2;

	//although checking, whether the number of GROMACS types is the same as from the *.dat file
	if (RunParams::potential==1 && RunParams::nGRtypes!=nGRtypes)
	{
		cout<<"Number of GROMACS types is "<<RunParams::nGRtypes<<" from the "<<datfilename<<" and "<<nGRtypes<<endl;
		cout<<"from the topology!"<<endl;
		cout << "The GROMACS types from the topology are the following:" << endl;
		for (i = 0; i < nGRtypes; i++)
			cout << atom_type[i] << endl;
		cout<<"Cannot run this way,exiting..."<<endl;
		CleanExit();
	}
	//ivirt = 0;
	for (imoltype=0;imoltype<nmoltype;imoltype++)
	{
		for (iat=0;iat<natoms_per_type[imoltype];iat++)
		{
			//iat is the GROMACS index of the atom in the molecule, it start with 0
			//first_index[cumul_atom[imoltype]+iat] is the RMC_index of the first occurence of this atomtype, if this molecule
			//	type would be the first in the system
			//RMCindex_offset[imoltype] is the offset of the RMC index to take onto account that might not be the first moleculetype
			//(nmol_per_type-1)*delta_index[cumul_atom[imoltype]+iat] is the offset of the last occurence's index relative to the first
			//first make some further check concerning the first index, which could not have been done in GetAtoms
			/*if (first_index[cumul_atoms[imoltype] + iat] > nmol_per_type[imoltype] * natoms_per_type[imoltype])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The first index (" << first_index[cumul_atoms[imoltype] + iat] << ") of the " << iat + 1 << ". atom in the " << imoltype + 1 << ". molecule type is " << endl;
				cout << "larger than than the total number of atoms in this moleculetype, which is " << nmol_per_type[imoltype] * natoms_per_type[imoltype] <<"."<< endl;
				cout << "The first index has to be between 1 and the total number of atoms in all the molecules of this molecule type, even if there are other" << endl;
				cout << "molecule types in the configuration!" << endl;
				cout << "Correct the topology, and try again! Cannot run this way, exiting..." << endl;
				CleanExit();
			}
			if (second_index[cumul_atoms[imoltype] + iat] > nmol_per_type[imoltype] * natoms_per_type[imoltype] && nmol_per_type[imoltype] * natoms_per_type[imoltype]>1)
			{//if there is just one instance of an 1 atom molecule, tthen second index does not matter
				cout << "\n*****ERROR*****" << endl;
				cout << "The second index (" << second_index[cumul_atoms[imoltype] + iat] << ") of the " << iat + 1 << ". atom in the " << imoltype + 1 << ". molecule type is " << endl;
				cout << "larger than than the total number of atoms in this moleculetype, which is " << nmol_per_type[imoltype] * natoms_per_type[imoltype] << "." << endl;
				cout << "The second index has to be between 2 and the total number of atoms in all the molecules of this molecule type, even if there are other" << endl;
				cout << "molecule types in the configuration!" << endl;
				cout << "Correct the topology, and try again! Cannot run this way, exiting..." << endl;
				CleanExit();
			}
			//check last index
			if (first_index[cumul_atoms[imoltype] + iat]+(nmol_per_type[imoltype]-1)* delta_index[cumul_atoms[imoltype] + iat] > nmol_per_type[imoltype] * natoms_per_type[imoltype] && nmol_per_type[imoltype] * natoms_per_type[imoltype] > 1)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The last index (" << first_index[cumul_atoms[imoltype] + iat] + (nmol_per_type[imoltype] - 1) * delta_index[cumul_atoms[imoltype] + iat] << ") of the " << iat + 1 << ". atom in the " << imoltype + 1 << ". molecule type is " << endl;
				cout << "larger than than the total number of atoms in this moleculetype, which is " << nmol_per_type[imoltype] * natoms_per_type[imoltype] << "." << endl;
				cout << "The last index has to be between 1 and the total number of atoms in all the molecules of this molecule type," << endl;
				cout << "even if there are other molecule types in the configuration!" << endl;
				cout << "Correct the topology, and try again! Cannot run this way, exiting..." << endl;
				CleanExit();
			}*/
			
			//RMC_index of the atom in the first instance
			if (first_index[cumul_atoms[imoltype] + iat] == 0)//virtual sites
			{
				continue;//cumul is not available for virtual sites yet
				//ind1 = SimpleCfg::cumul[SimpleCfg::ntypes+ ivirt];
				//ivirt++;
			}
			else
				ind1=first_index[cumul_atoms[imoltype]+iat]-1+RMCindex_offset[imoltype];
			ind2=ind1+(nmol_per_type[imoltype]-1)*delta_index[cumul_atoms[imoltype]+iat];
			//find the type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
			RMC_type1=0;
			while (ind1>=SimpleCfg::cumul[(RMC_type1) +1])
				(RMC_type1)++;
			//RMC_type1 now points to the type of the atom
			//find the type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
			RMC_type2=0;
			while (ind2>=SimpleCfg::cumul[(RMC_type2) +1])
				(RMC_type2)++;
			//RMC_type2 now points to the type of the atom
			if (RMC_type1!=RMC_type2)
			{
				cout << "\n****ERROR****"<<endl;
				cout<<"Molecule type: "<<imoltype+1<<", GROMACS atom index: "<<iat+1<<endl;
				cout<<"The RMC type= "<<RMC_type1+1<<", RMC atom index= "<<ind1+1<<" for the first instance of this molecule type,"<<endl;  
				cout<<"The RMC type= "<<RMC_type2+1<<", RMC atom index= "<<ind2+1<<" for the last instance of this molecule type"<<endl;  
				cout<<"They are not the same, something is wrong with the topology-configuration consistency!"<<endl;
				cout<<"An atom has to have the same RMC type in all the molecules of a given molecule type!"<<endl;
				cout<<"Correct it, and try again, exiting..."<<endl;
				CleanExit();
			}
		}
	}
				
};

//check, whether all the virtual sites of the [ atoms ] section are defined
//chack, whether in the bonds, angle, dihedral and pair section only real atoms are given
void Topology::CheckVirtualSites()
{
	int imoltype, iat, i,j,k,ind1,ind2,ind3,ind4,my_ch_group;
	int found;
	int *pind;
	bool rearrange=false;
	
	iat = 0;
	//counting the number of virtual sites in the [ atoms ] section
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (i = 0; i < natoms_per_type[imoltype]; i++)
		{
			if (first_index[iat] == 0)
				virtual_count ++;//increasing the number of virtual sites
			iat++;
		}
	}
	if (virtual_count != n_gr_virtuals)//not all the virtual sites are defined, or multiple definition
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The total number of virtual sites given in the [ atoms ] section of the topology file(s) is " << virtual_count <<","<< endl;
		cout << "while total " << n_gr_virtuals << " virtual sites were declared in the [ virtual_siteX ] sections!" << endl;
		cout << "Correct the topology file(s), and try again!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	if (n_gr_virtuals == 0)
		return;
	i = 0;

	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (iat = 0; iat < natoms_per_type[imoltype]; iat++)
		{
			found = -1;
			if (first_index[i] == 0)
			{
				for (j = 0; j < nvirtuals_per_type[imoltype]; j++)
				{
					if (virtual_site[cumul_virtuals[imoltype] + j].atom_index == iat + 1)//this virtual site is defined in the [virtal_sitesX] section
					{
						if (found>-1)//this virtual site was already described in [ virtual_sitesX ]section
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The " << iat + 1 << ". virtual site from the [ atoms ] section of the topology of the " << imoltype + 1 << ". molecule type" << endl;
							cout << "is defined more than once in the [ virtual_siteX ] section(s). The first definition is the " << found + 1 << endl;
							cout << "virtual site in the [ virtual_siteX ] section(s), the second is the " << j + 1 << "." << endl;
							cout << "Correct the topology file, and try again! Cannot run this way,exiting..." << endl;
							CleanExit();
						}
						else
						{
							virtual_site[cumul_virtuals[imoltype] + j].assigned = iat;
							found = j;
						}
					}
				}//j cycle through the virtual site structures defined in [ virtual_sitesX ]
			}//if virtual site
			i++;
		}//iat cycle
	}//imoltype cyle
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (j = 0; j < nvirtuals_per_type[imoltype]; j++)
		{
			if (virtual_site[cumul_virtuals[imoltype] + j].assigned == -1)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The " <<j + 1 << ". virtual site from the [ virtual_sitesX ] section(s) of the topology of the " << imoltype + 1 << ". molecule type" << endl;
				cout << "does not define a virtual site in the [ atoms ] section." << endl;
				cout << "Correct the topology file, and try again! Cannot run this way,exiting..." << endl;
				CleanExit();
			}
		}
	}
	//check, whether only real atoms are involved in the building of virtual sites, as RMC can only handle this
	//To make calculation easier the all the atoms building a virtual site should belong to the same charge group, 
	//and the virtual site as well. Check it!
	//Check, whether the host index is valid!
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (j = 0; j < nvirtuals_per_type[imoltype]; j++)
		{
			my_ch_group = charge_group[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + j].atom_index - 1];
			for (k = 0; k < virtual_site[cumul_virtuals[imoltype] + j].tnumb; k++)//go through the atoms building the virtual site
			{
				ind1 = virtual_site[cumul_virtuals[imoltype] + j].indices[k];//GROMACS index in the molecule
				if (first_index[cumul_atoms[imoltype] +ind1- 1] == 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The atom index " << ind1 << " in the " << j + 1 << ". virtual site" << endl;
					cout << "of the " << imoltype + 1 << ". molecule type refers not to a real atom, but a virtual site!" << endl;
					cout << "This is not allowed in RMC, correct the topology and try again! Exiting..." << endl;
					CleanExit();
				}
				if (my_ch_group != charge_group[cumul_atoms[imoltype] + ind1 - 1])
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The charge group of the atom index " << virtual_site[cumul_virtuals[imoltype] + j].atom_index << " virtual site" << endl;
					cout << "of the " << imoltype + 1 << ". molecule type is " << my_ch_group+1 << ", which is not the same as the charge group "<< charge_group[cumul_atoms[imoltype] + ind1 - 1]+1<<" of the " << endl;
					cout << k + 1 << ". atom with atom index " << ind1 << " building the site!" << endl;
					cout << "(The charge group indices given here are not necessarily the same as in the topology, they are referring to" << endl;
					cout << "all the charge groups in the configuration!)"<<endl;
					cout << "In RMC_POT all the atoms building a virtual site and the virtual site itself have to belong to the same charge group!" << endl;
					cout << "Correct the topology and try again! Exiting..." << endl;
					CleanExit();
				}

			}//end of k cycle for the building atoms
			if (virtual_site[cumul_virtuals[imoltype] + j].host_index > natoms_per_type[imoltype])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The atom index for the host atom " << virtual_site[cumul_virtuals[imoltype] + j].host_index << " in the " << j + 1 << ". virtual site" << endl;
				cout << "of the " << imoltype + 1 << ". moleculetype is larger than the number of atoms in the molecule ("<< natoms_per_type[imoltype]<<")!" << endl;
				cout << "Correct the topology and try again! Exiting..." << endl;
				CleanExit();
			}
			else
			{
				if (first_index[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + j].host_index - 1] == 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The atom index for the host atom " << virtual_site[cumul_virtuals[imoltype] + j].host_index << " in the " << j + 1 << ". virtual site" << endl;
					cout << "of the " << imoltype + 1 << ". moleculetype refers not to a real atom, but a virtual site!" << endl;
					cout << "Correct the topology and try again! Exiting..." << endl;
					CleanExit();
				}
			}
		}//j cycle through the virtual site structures defined in [ virtual_sitesX ]
	}//imoltype cyle
	//now check the bonded interactions, whether no virtial atoms are given by mistake
	pind = pair_index;
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (j = 0; j < npairs_per_type[imoltype]; j++)
		{
			ind1 = *pind++;
			ind2 = *pind++;

			if (first_index[cumul_atoms[imoltype] + ind1 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The first index (" << ind1 << ") in the " << j + 1 << ". pair of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind2 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The second index (" << ind2 << ") in the " << j + 1 << ". pair of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
		}
	}
	pind = bond_index;
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (j = 0; j < nbonds_per_type[imoltype]; j++)
		{
			ind1 = *pind++;
			ind2 = *pind++;

			if (first_index[cumul_atoms[imoltype] + ind1 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The first index (" << ind1 << ") in the " << j + 1 << ". bond of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind2 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The second index (" << ind2 << ") in the " << j + 1 << ". bond of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
		}
	}
	pind = angle_index;
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (j = 0; j < nangles_per_type[imoltype]; j++)
		{
			ind1 = *pind++;
			ind2 = *pind++;
			ind3 = *pind++;

			if (first_index[cumul_atoms[imoltype] + ind1 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The first index (" << ind1 << ") in the " << j + 1 << ". angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind2 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The second index (" << ind2 << ") in the " << j + 1 << ". angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind3 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The third index (" << ind3 << ") in the " << j + 1 << ". angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
		}
	}
	
	pind = dihedral_index;
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (j = 0; j < ndihedrals_per_type[imoltype]; j++)
		{
			ind1 = *pind++;
			ind2 = *pind++;
			ind3 = *pind++;
			ind4 = *pind++;

			if (first_index[cumul_atoms[imoltype] + ind1 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The first index (" << ind1 << ") in the " << j + 1 << ". dihedral angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind2 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The second index (" << ind2 << ") in the " << j + 1 << ". dihedral angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind3 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The third index (" << ind3 << ") in the " << j + 1 << ". dihedral angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
			if (first_index[cumul_atoms[imoltype] + ind4 - 1] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The fourth index (" << ind4 << ") in the " << j + 1 << ". dihedral angle of the " << imoltype + 1 << ". moleculetype " << endl;
				cout << "does not belong to a real atom, but to a virtual site, as it has 0 for its RMC index!" << endl;
				cout << "Cannot run this way, correct the topology and try again!" << endl;
				CleanExit();

			}
		}
	}
	//now check, whether the virtual sites are given in the same order, as in the atom section, if not eorder them
	//accoording to increasing atom index for each moltype 

	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (i = 0; i < nvirtuals_per_type[imoltype] - 1; i++)
		{
			for (j = i + 1; j < nvirtuals_per_type[imoltype]; j++)
			{
				if (virtual_site[cumul_virtuals[imoltype] + i].atom_index > virtual_site[cumul_virtuals[imoltype] + j].atom_index)
				{
					rearrange = true;
					break;

				}
			}
			if (rearrange)
				break;
		}
	
		O_Sort(nvirtuals_per_type[imoltype], virtual_site+cumul_virtuals[imoltype]);
	}


	
};

//giving the number of virtual sites the RMC_ind atom is involved in,
//optionally giving back the RMC indices of the virtual sites the RMC_ind atoms is involved in,
//and the RMC indices of atoms building the virtual sites, both pointers point to arrays supposed to have sufficent size to store the values
/*int Topology::GetVirtual(int RMC_ind, int *virtual_indices, int *RMCind_source)
{
	

	
	imoltype = 0;

	for (imoltype = 0; imoltype < nmoltype - 1; imoltype++)
	{
		if (RMC_ind >= RMCindex_offset[imoltype] && RMC_ind < RMCindex_offset[imoltype])
			break;
	}
	//imoltype will be the GROMACS molecule type of the RMC_ind atom
	for (vi = 0; vi < nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
	{
		found = -1;
		for (j = 0; j < virtual_site[cumul_virtuals[imoltype] + vi].tnumb; j++)//going through the indices of the atoms making the site 
		{
			//RMC_index of the atom building the site in the first instance, and the incement for the following ones
			ind1 = first_index[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + vi].indices[j] - 1] - 1 + RMCindex_offset[imoltype];
			delta = delta_index[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + vi].indices[j] - 1];
			if (RMC_ind >= ind1 && RMC_ind <= (nmol_per_type[imoltype] - 1) * delta)
			{
				//the index is in the right range

				temp = (double)(RMC_ind - ind1) / (double)delta;
				if (fabs(temp - (int)temp) < TOLERANCE)//RMC_ind is one of the ind1+imol*delta series, and belongs to the this virtual site
				{
					imol = (int)temp;
					found = 1;
				}
			}
			if (found)
			{
				if (RMCind_source != NULL)
				{
					for (i = 0; i < virtual_site[cumul_virtuals[imoltype] + vi].tnumb; i++)
					{
						//The index of the first instance and the increment  
						ind1 = first_index[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + vi].indices[j] - 1] - 1 + RMCindex_offset[imoltype];
						delta = delta_index[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + vi].indices[i] - 1];
						RMCind_source[count * 4 + i] = ind1 + imol * delta;
					}
				}
				if (virtual_indices != NULL)//normally max_diff_virtual was already set
				{
					//for a virtual site first_index is 0, has to get RMC index differently, for the ith virtual site simply as the values in the SimpleCFg::cumul[ntypes+i] array
					ind1 = SimpleCfg::cumul[SimpleCfg::ntypes + cumul_virtuals[imoltype] + vi];
					delta = delta_index[cumul_atoms[imoltype] + virtual_site[cumul_virtuals[imoltype] + vi].atom_index - 1];
					virtual_indices[count] = ind1 + imol * delta;
				}
				count++;
				break;//j cycle
			}
		}//j cycle
	}//vi cyle
	return count;
};*/

//rendering the vdW14 parameters to GROMACS types
void Topology::CreatevdW14Params()
{
	int i,ipair,ind1,ind2,type1,type2,partialind;

	if (Topology::npairs>0)
	{
		SetArraysize(&RunParams::used_14partials,SimpleCfg::npartials,"RunParams::used_14partials","Topology::CreatevdW14Params");
		SetArraysize(&vdW14_pot1_partial,nGRtypes*(nGRtypes+1)/2,"vdW14_pot1_partial","Topology::CreatevdW14Params");
		SetArraysize(&vdW14_pot2_partial,nGRtypes*(nGRtypes+1)/2,"vdW14_pot2_partial","Topology::CreatevdW14Params");
	
		for (i=0;i<SimpleCfg::npartials;i++)
			RunParams::used_14partials[i]=0;

		for (i=0;i<nGRtypes*(nGRtypes+1)/2;i++)
		{
			vdW14_pot1_partial[i]=-1.0;
			vdW14_pot2_partial[i]=-1.0;
		}

		for (ipair=0;ipair<npairs;ipair++)
		{
			ind1=pair_index[ipair*2];
			ind2=pair_index[ipair*2+1];
			type1=GROMACS_type[ind1-1];
			type2=GROMACS_type[ind2-1];
			partialind=(type1<=type2 ? (type1*nGRtypes-(type1*(type1+1)/2)+type2) : (type2*nGRtypes-(type2*(type2+1)/2)+type1));
			
			if (fabs(vdW14_pot1_partial[partialind]- -1.0)<TOLERANCE)//this was not set yet
				vdW14_pot1_partial[partialind]=vdW14_pot1[pair_type[ipair]];
			else
			{
				if (fabs(vdW14_pot1_partial[partialind] - vdW14_pot1[pair_type[ipair]]) > TOLERANCE)
				{

					cout << "\n*****ERROR*****" << endl;
					cout << "There are multiple 1-4 first potential parameter entries for the GROMACS partial " << partialind + 1 << endl;
					cout << "The multiple values are " << vdW14_pot1_partial[partialind] << " and " << vdW14_pot1[pair_type[ipair]] << " and the repetition with wrong parameter occures " << endl;
					cout<<"at pair "<<ind1<<" - "<<ind2<<"! (The length dimension is Angstrom, if applies)" << endl;
					cout << "There are something wrong with the entries, check the topology..." << endl;
					CleanExit();
				}
			}
			if (fabs(vdW14_pot2_partial[partialind]- -1.0)<TOLERANCE)//this was not set yet
				vdW14_pot2_partial[partialind]=vdW14_pot2[pair_type[ipair]];
			else
			{
				if (fabs(vdW14_pot2_partial[partialind] - vdW14_pot2[pair_type[ipair]]) > TOLERANCE)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "There are multiple 1-4 second potential parameter entries for the GROMACS partial " << partialind + 1 << endl;
					cout << "The multiple values are " << vdW14_pot2_partial[partialind] << " and " << vdW14_pot2[pair_type[ipair]] << " and the repetition with wrong parameter occures " << endl;
					cout<<"at pair " << ind1 << " - " << ind2 << "! (The length dimension is Angstrom, if applies)" << endl;
					cout << "There are something wrong with the entries, check the topology..." << endl;
					CleanExit();
				}
			}
		}
	}
};



//checking, whether there is non-zero sigma, if 0 sigma is given for an interaction type of each molecule type 
//If 0 is given for an interaction type's sigma in a molecule type, then this interaction type's sigma will be 
//identical to the first interaction type of this molecule with non-zero sigma, or if there is not non-zero sigma
//for any bond of this molecule type, than the first overall bondtype with non-zero sigma regardless the molecule type
//(keeping the potential ratio the same, as in the MD)
//It will be called for all the existing bonded interaction type of the molecules (bond, angles dihedrals)
void Topology::CheckSigma(int *weight_mode, int *cumul_inttypes, int *ninttypes_per_type, int *zero_sigma, double *sigma, const char *interaction_name, int *leading_sigma)
{
	int imoltype,inttype;
	int is_zero, all_zero=1, neg_sigma_ind=0;
	int pos_sigma_same=1;//all the positive sigma is the same by default
	int flag=0;//indicating, whether to search for the non-zero sigma, only in case of the bond
	int nonzero_sigma = 0;//indocator, that for at least one bond type of at least one mol type there is a non-zero sigma, which can be used to set the others
	double pos_sigma=-1;
	if (strcmp(interaction_name, "bond") == 0)
		flag = 1;
	
	*zero_sigma=0;
	for (imoltype=0;imoltype<nmoltype;imoltype++)
	{
		
		is_zero=0;
		//first find the first non-zero sigma for this molecule type
		if (flag)
			zero_sigma[imoltype+1]=-1;
		for (inttype=0;inttype<ninttypes_per_type[imoltype];inttype++)
		{
			if (sigma[cumul_inttypes[imoltype]+inttype]>0 || sigma[cumul_inttypes[imoltype]+inttype]<0  )
			{
				//find the first non-zero sigma for the bonds, this will be used to set the zero sigmas
				if (flag && zero_sigma[imoltype + 1] == -1)
				{
					zero_sigma[imoltype + 1] = cumul_inttypes[imoltype] + inttype;//index of the first non-zero sigma for this molecule
					nonzero_sigma++;
					if (*leading_sigma == -1)
						*leading_sigma = cumul_inttypes[imoltype] + inttype;
				}
				
				//check, whether the positive sigmas are the same
				if (sigma[cumul_inttypes[imoltype]+inttype]>0)
				{
					if (pos_sigma<0 )
						pos_sigma=sigma[cumul_inttypes[imoltype]+inttype];//this is the first positive sigma
					else
					{
						//check, if this is the same as the first positive 
						if (fabs(sigma[cumul_inttypes[imoltype]+inttype]-pos_sigma)>TOLERANCE)
							pos_sigma_same=0;//not all the positive sigma is the same
					}

				}
				
				if (sigma[cumul_inttypes[imoltype]+inttype]<0)
					neg_sigma_ind++;
				
			}
			//see, if it is zero
			if (fabs(sigma[cumul_inttypes[imoltype]+inttype]-0.0)<TOLERANCE)
				is_zero=1;
			else
				all_zero=0;//this is not zero, so not all of it is
		}
		
		if (is_zero )
			*zero_sigma=1;
	}
	
	//if there is zero sigma, and there is no nonzero sigma for any bond type of any mol type, 
	//or if it is not a bond, check, if there is zero sigma, wheteher there is nonzero bond sigma somewhere, if not give error
	//if it not bond, it is necessary to check, as it is possible, that the topology is incorrect, and no bond was given at all
	if ((flag && *zero_sigma && nonzero_sigma==0) || (flag==0 && *zero_sigma && *leading_sigma==-1))
	{
		
		cout << "\n***** ERROR *****" << endl;
		cout << "There is no bond interaction with non-zero sigma value for any molecule type";
		if (flag)
			cout << "!" << endl;
		else
			cout << ",\nand there is zero sigma for at least one " << interaction_name << " interaction!" << endl;
		cout << "The sigma values for the interaction types with zero sigma cannot be set! Only use zero sigma for an" << endl;
		cout << "interaction type in any molecule type, if there is at least one bond interaction type with positive (absolute)" << endl;
		cout << "or negative (scalable) sigma parameter for at least one of the molecule types!" << endl;
		cout << "Correct the toplogy or include topology file, and try again!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	//for bonds:there is zero sigma, and there is not nonzero sigma for all the moltypes, 
	//set the leading_sigma for those, which does not have non-zero sigma
	if (flag && *zero_sigma && nonzero_sigma > 0 && nonzero_sigma != nmoltype)
	{
		for (imoltype = 0; imoltype < nmoltype; imoltype++)
		{
			if (zero_sigma[imoltype + 1] == -1)
				zero_sigma[imoltype + 1] = *leading_sigma;
		}
	}
	//set the weight mode, 
	//0: The chi2 calculation for the interaction type can be done together, and only one chi2 value is displayed for the interaction
	//1: the chi2 calculation of the different types of the given bonded interaction has to be done and displayed separately,
	//		as not all of them have the same sigma
	if (fabs(RunParams::NB_weight_mode)==2)
		*weight_mode=0;//the same weight (vdW_weight[0]) will be used for all the potential contributions, regardless it fix or scalable
	else
	{
		if (flag)//for bonds
		{
			//all the positive sigmas are the same, no negative scalable (can be zero)
			if (pos_sigma>0 && pos_sigma_same && neg_sigma_ind==0)
				*weight_mode=0;//can be displayed together, as the zero sigma will be set to the non-zero value
			if (nmoltype==1)
			{
				if (ninttypes_per_type[0]==1 || (neg_sigma_ind==1 && pos_sigma<0))//only one bond type or one negative sigma and the others 0
					*weight_mode=0;//can be displayed together, as the zero sigma will be set to the non-zero value

			}
		}
		else//for angles and dihedral angles
		{
			//all the positive sigmas are the same, no negative scalable or zero; or all zero, and the bond weight mode is collapsed
			if ((pos_sigma>0 && pos_sigma_same && neg_sigma_ind==0 && *zero_sigma==0) || (all_zero==1 && bond_weight_mode==0))
				*weight_mode=0;
			if (nmoltype==1)
			{
				if (ninttypes_per_type[0]==1 || (all_zero==1 ))//only one bond type or all zero
					*weight_mode=0;//can be displayed together, as the zero sigma will be set to the non-zero value

			}
		}
	}
};

//resize the virtual site structure array
void Topology::ResizeVirtualSiteStruct(int *max, int new_max, virtual_site_struct **array, const char *array_name, const char *routine_name)
{
	int i,j;
	bool isarray = true;
	virtual_site_struct *temp_array, *newp, *oldp;

	if (*array == NULL)//was not created
		isarray = false;
	if (isarray)
	{
		temp_array = new virtual_site_struct[*max ];//the original array will be copied here
		if (temp_array == NULL)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Resizing of " << array_name << " array in ResizeVirtualSiteStruct called from " << routine_name << " routine was not successful," << endl;
			cout << "as temporary array cannot be created, most probably because of the shortage of memory." << endl;
			CleanExit();
		}

		newp = temp_array;
		oldp = *array;
		for (i = 0; i < *max; i++)
		{
			(*newp).atom_index = (*oldp).atom_index;//copy the original to temp
			(*newp).host_index = (*oldp).host_index;
			(*newp).tnumb = (*oldp).tnumb;
			for (j = 0; j < 4; j++)
				(*newp).indices[j] = (*oldp).indices[j];
			for (j = 0; j < 3; j++)
				(*newp).params[j] = (*oldp).params[j];
			(*newp).type= (*oldp).type;
			(*newp).assigned = (*oldp).assigned;
			newp++;
			oldp++;


		}
		delete[] * array;
	}

	*array = new virtual_site_struct[new_max ];//creating the larger new array
	if (*array == NULL)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Resizing of " << array_name << " array in ResizeVirtalSiteStruct called from " << routine_name << " routine was not successful," << endl;
		cout << "as the new, larger array cannot be created, most probably because of the shortage of memory." << endl;
		CleanExit();
	}
	if (isarray)
	{
		newp = *array;
		oldp = temp_array;
		for (i = 0; i < *max; i++)
		{
			(*newp).atom_index = (*oldp).atom_index;//copy the original to temp
			(*newp).host_index = (*oldp).host_index;
			(*newp).tnumb = (*oldp).tnumb;
			for (j = 0; j < 4; j++)
				(*newp).indices[j] = (*oldp).indices[j];

			for (j = 0; j < 3; j++)
				(*newp).params[j] = (*oldp).params[j];
			(*newp).type = (*oldp).type;
			(*newp).assigned = (*oldp).assigned;
			newp++;
			oldp++;

		}

		delete[] temp_array;
	}
	*max = new_max;
};

//check, whether all the entries are different, line, moltype starts  with 1
void Topology::CheckRepetition(int *array, int count, int line, int moltype, char *directive)
{
	int i,j;
	for (i = 0; i < count-1; i++)
	{
		for (j=i+1;j<count;j++)
			if (array[i] == array[j])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The " << i+1 << ". " << " and " << j+1 << ". " << " atom indices are the same, value "<<array[i]<<"."<<endl;
				cout<<"This occurs for the "<< moltype<<". molecule type's "<<directive<<" directive "<< line << ". line." << endl;
				cout << "The atom indices has to be different, correct the topology file and try again!" << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}
	}
};

//Check, whether the exclusions given in the [ exclusions ] directive are selfconsistent
//it is important to make sure, that if A atom is given as an exclusion for B, then B should be given for A as well!
//Otherwise the histogram calculation and update will have errors
void Topology::CheckExclusions()
{
	int imoltype,il1,il2,iex1,iex2,iloffset1,iloffset2;
	int found;
	int excl_ind1,excl_ind2;
	bool contract = false;
	
	for (imoltype = 0; imoltype < nmoltype; imoltype++)
	{
		for (il1 = 0; il1 < nexclusionlines_per_type[imoltype]; il1++)
		{
			iloffset1= Topology::cumul_exclusionlines[imoltype] + il1;//the index of the exclusion line
			if (exclusion_centre[iloffset1] == -1)//it was removed already go to next
				continue;
			//atom A is exclusion_centre[iloffset1]
			for (iex1 = 0; iex1 < nexcluded_neighs[iloffset1]; iex1++)//going through the exclusions in iloffset1 exclusion line
			{
				excl_ind1 = exclusions[cumul_nexcluded_neighs[iloffset1] + iex1];//atom B index of the exclusion in iloffset1 exclusion line
				found = 0;//not found by default
				for (il2 = 0; il2 < nexclusionlines_per_type[imoltype]; il2++)
				{
					iloffset2 = Topology::cumul_exclusionlines[imoltype] + il2;//the index of the exclusion line
					if (il1 == il2 || exclusion_centre[iloffset2] == -1)//it is the same line or was removed already go to next
						continue;//this is the same line, go to next
					if (excl_ind1 == exclusion_centre[iloffset2])//the exclusion has its own exclusion line
					{
						for (iex2 = 0; iex2 < nexcluded_neighs[iloffset2]; iex2++)//going through the exclusions in iloffset2 exclusion line
						{
							excl_ind2 = exclusions[cumul_nexcluded_neighs[iloffset2] + iex2];//index of the exclusion in iloffset2 exclusion line
							if (excl_ind2 == exclusion_centre[iloffset1])
							{
								if (found > 0)//it was given already, delete it
								{
									cout << "\nWARNING(" << ++warn << "): The excluded atom index " << excl_ind1 << " for exclusion centre " << endl;
									cout<<"\t"<<exclusion_centre[iloffset1] << " for molecule type " << imoltype + 1 << endl;
									cout << "\t(" << il1 + 1 << ". line in [ exclusions ] section) is given for exclusion centre " << exclusion_centre[iloffset2] << " both in exclusion line " << endl;
									cout <<"\t"<< found << " and " << il2 + 1 << "! The second occurrence will be eliminated and continuing, but correct the topology!" << endl;
									exclusions[cumul_nexcluded_neighs[iloffset2] + iex2] = -1;//the array will be contracted later
									contract = true;
								}
								else
								{
									found = il2 + 1;//given in both ways, mark it with the index of the exclusion line index in its own file i2 STARTING WITH 1
									break;//iex2 cycle, no need to search this line any more
								}
							}
						}//end of iex2 cycle for the exclusions
						if (found <1)
							found = -1;//indicating that atom B has its own line, but A was not given on it, there can be an other line theoretically for atom B, keep searching
					}
				}//end of il2 cycle for the exclusion lines
				if (found < 1)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The excluded atom index " << excl_ind1 << " for exclusion centre " << exclusion_centre[iloffset1] << " for molecule type " << imoltype + 1 << endl;
					cout << "can be found in the " << il1 + 1 << ", line of the [ exclusions ] section." << endl;
					
					if (found == -1)
						cout << "Athough exclusion line is given for atom " << excl_ind1 << " no exclusion " << exclusion_centre[iloffset1]<<" can be found!"<<endl;
					else
						cout << "No exclusion line is given for atom " << excl_ind1 << "!"<<endl;
					cout << "Check the topology and correct it! Exclusions in the exclusion section has to be declared in both way!" << endl;
					cout << "Cannot run this way, exiting..." << endl;
					CleanExit();
				}
			}//end of iex1 cycle for the exclusions
		}//end of il1 cycle for the exclusion lines
	}//end of imoltype cycle
	//there was deletion, contract the exclusion lines by filling the emptyplaces and update the number of exclusions in nexcluded_neighs array
	//cumul_nexcluded_neighs array will not be changed, use always explicit indices in exclusions array using cumul_nexcluded_neighs, and not go through with a pointer
	if (contract)
	{
		for (il1 = 0; il1 < nexclusionlines; il1++)
		{
			iex2 = 0;
			for (iex1 = 0; iex1 < Topology::nexcluded_neighs[il1]; iex1++)
			{
				if (Topology::exclusions[Topology::cumul_nexcluded_neighs[il1] + iex1] != -1)
				{
					Topology::exclusions[Topology::cumul_nexcluded_neighs[il1] + iex2] = Topology::exclusions[Topology::cumul_nexcluded_neighs[il1] + iex1];
					iex2++;
				}
			}
			Topology::nexcluded_neighs[il1] = iex2;
			if (Topology::nexcluded_neighs[il1] == 0)
				Topology::exclusion_centre[il1] = -1;//no exclusions for this atom, remove it
			
		}
	}
};

// Sorting is done with the Oosterwal algorithm for the virtual site structure array 
void Topology::O_Sort(int length, virtual_site_struct *array)
{
	int i, step;
	double SngFib, SngPhi;
	virtual_site_struct temp;

	SngPhi = 0.78;//Define phi value
	SngFib = length * SngPhi;//Set initial real step size
	step = int(SngFib);//Set initial integer step size

	while (step > 0)
	{

		for (i = 0; i < length - step; i++)//Range of lower cell
		{
			if (array[i].atom_index > array[i + step].atom_index) //Compare cells
			{
				CopyVirtual(&temp, &array[i]);
				CopyVirtual(&array[i],&array[i + step]);//Swap Cells
				CopyVirtual(&array[i + step],&temp);
			}
		}

		SngFib = SngFib * SngPhi;//Decrease the Real step size
		step = int(SngFib);//Set the integer step value

	}
}
void Topology::CopyVirtual(virtual_site_struct *dest, virtual_site_struct *source)//for the sort
{
	int j;

	dest->atom_index = source->atom_index;
	dest->host_index = source->host_index;
	dest->tnumb = source->tnumb;
	for (j = 0; j < 4; j++)
		dest->indices[j] = source->indices[j];

	for (j = 0; j < 3; j++)
		dest->params[j] = source->params[j];
	dest->type = source->type;
	dest->assigned = source->assigned;
}
/*
void Topology::ReadTopology(const std::list<string> defs, const string &topology_name)
//This function reads the topology file in
{
  std::list<string> topfiles;
  topfiles.push_back(topfilename);
  //Discover the full list of topology files
  for(std::list<string>::iterator it=includes.begin(); it!=includes.end(); it++)
  {
    //Open the file for read
    ifstream inp(it->c_str());
    if(!inp)
    {
      cout<<"\nERROR: the \""<<it->c_str()<<"\"topology file does not exist!"<<endl;
      CleanExit();
    } 
    Topology::GetIncludes(defs,*it,inp,includes);
    inp.close();
  }
}
*/


