//source FNC_POT.cpp
//Last changed 24.01.2023

//It contains the Fixed Neighbour Constraint-related things, 
//or in case of POTENTIAL: RunParams::fnc=4 the flexible molecule handling,
//						   RunParams::potential>0 non-bonded interactions, sofar 1:LJ, 10:tabulated


//The working of the FNC constraint is changed a bit from the combined rmc version 1.0. Now there are several option possibilities for FNC:
//If FNC option is given in the *.dat file as before, but now it can be 0, 1, 2, 3
//Option 0: no FNC (as before).
//Option 1: FNC is applied as before: the minimum and maximum distances provided in the *.fnc file strictly applied, if
//			the starting configuration does not satisfy the FNC constraint, than the simulation cannot be started.
//Option 2: FNC is applied, the FNC range is widened. The minimum and maximum FNC distances provided in the *.fnc file are reset in case the FNC is not satisfied for
//			the starting configuration by the smallest and largest distances occuring in the starting configuration for the given constraint.
//			Later on the FNC is applied as before, only distances inside the constrained range are accepted. Each FNC constraint
//			is handled separately, only those are reset, which necessary.
//Option 3:	FNC is applied, the FNC range is narrowed. If the starting configuration does not satisfy the FNC, only warning is given and the simulation continues. If a new
//			FNC constrained distances is closer to the FNC range after a move than before, than it is accepted, hopefully guiding the simulation to satisfy a smaller range
//			FNC constraint, than could be found for the starting configuration. 
//Option 4: flexible, MD-like molecules, kept together by forces
//			No *.fnc file is needed, instead GROMACS type *.top file with few additional parameter will define the molecular structure
//fnc_conflict in case of automatic cutoff will be set first to 1 to go into fnc neighbour determination during histogram calculation,
//and when the cutoff not including fnc neighbours are determined, fnc_conflict is determined again
//it has to be noted, that the same fnc constraint can be used for atom pairs with different type combination (belonging to different partials!) 

#define _DEF_FILES //not redefine the file names included through global.h
#include "Threads.h"//classes1.h is included through this

//Defining static members
int  FNC_POT::nconstraints;//the number of FNC_POT constraints
int  FNC_POT::natoms;//number of FNC_POT central atoms
int  FNC_POT::nvirtuals;//number of virtual sites
int  FNC_POT::ndist;//the total number of constrained distances define the size of the arrays, normal FNC
int  FNC_POT::sum_out_of_range;//total number of pairs out of the FNC range (for option 3)
int *FNC_POT::nout_of_range;//number of pairs/fnc type out of the FNC range (for option 3)
int *FNC_POT::dec_out_of_range;//decrement the number of out of range atoms /fnc type, if the move is accepted
int  FNC_POT::index_offset;//1 for BOND and 0 for normal fnc
int  FNC_POT::nexclusions=0;//total number of excluded atoms
int  FNC_POT::check_exclusions=0;//boolean whether to check the exclusions during potential correction due to the move of the charge group centres
int  FNC_POT::pot_dim=0;//dimension of the potential array
int  FNC_POT::nGRpartials=0;//number of GR partials for the NB parameter arrays
int  FNC_POT::nthreads=0;//number of threads
int  FNC_POT::vdW_warning = 0;//0 no warning, index of the partial (start with 1) where larger than POT_WARNING absolute value for the exclusion potential last occured
int  FNC_POT::Coul_warning = 0;//0 no warning, index of the partial (start with 1) where larger than POT_WARNING absolute value for the exclusion potential last occured
int	 FNC_POT::nused_potpartials;//number of partials to use for tabulated potential (number of potential files)
int *FNC_POT::npot_points;//array of number of data points in the tabulated potential files
int *FNC_POT::tabpot_offset;//offset for each used partial for r_pot and u_pot
int *FNC_POT::tabpot_index;//index of the partial to use the tabulated potential for, for each used partial
double  FNC_POT::bond_tot_pot=0.0;//total bond potential
double  FNC_POT::angle_tot_pot=0.0;//total angle potential
double  FNC_POT::perdihedral_tot_pot=0.0;//total periodic dihedral potential
double  FNC_POT::harmdihedral_tot_pot=0.0;//total harmonic dihedral potential
double  FNC_POT::RBdihedral_tot_pot=0.0;//total RB dihedral potential
double  FNC_POT::vdW_tot_pot=0.0;//total vdW potential
double  FNC_POT::Coulomb_tot_pot=0.0;//total Coulomb potential
double  FNC_POT::vdW14_tot_pot=0.0;//total 1-4 vdW potential
double  FNC_POT::Coulomb14_tot_pot=0.0;//total 1-4 Coulomb potential
double  FNC_POT::f_Coulomb_A=0.0;//f_Coulomb/boxedge to transform the distance from reduced unit to A
double *FNC_POT::bond_pot=0;//potential energy of the bonds 
double *FNC_POT::angle_pot=0;//potential energy of the angles 
double *FNC_POT::perdihedral_pot=0;//potential energy of the periodic dihedrals
double *FNC_POT::harmdihedral_pot=0;//potential energy of the harmonic dihedrals
double *FNC_POT::RBdihedral_pot=0;//potential energy of the RB dihedrals 
double *FNC_POT::vdW_pot=0;//non-bonded potential energy  
double *FNC_POT::vdW_pot_s = 0;//non-bonded potential energy for the samll values
double *FNC_POT::vdW_pot_l = 0;//non-bonded potential energy for the large values
double *FNC_POT::vdW14_pot=0;//non-bonded potential energy for 1-4 interaction 
double *FNC_POT::Coulomb_pot=0;//Coulomb potential energy
double *FNC_POT::Coulomb_pot_s = 0;//Coulomb potential energy for the small values
double *FNC_POT::Coulomb_pot_l = 0;//Coulomb potential energy for the large values
double *FNC_POT::Coulomb14_pot=0;//Coulomb potential energy for 1-4 interaction 
double *FNC_POT::bond_pot_old=0;//potential energy of the bonds for the previous move
double *FNC_POT::angle_pot_old=0;//potential energy of the angles for the previous move
double *FNC_POT::perdihedral_pot_old=0;//potential energy of the periodic dihedrals for the previous move
double *FNC_POT::harmdihedral_pot_old=0;//potential energy of the harmonic dihedrals for the previous move
double *FNC_POT::RBdihedral_pot_old=0;//potential energy of the RB dihedrals for the previous move
double *FNC_POT::vdW_pot_old=0;//non-bonded potential energy for the previous move
double *FNC_POT::vdW14_pot_old=0;//non-bonded potential energy for 1-4 interaction for the previous move
double *FNC_POT::Coulomb_pot_old=0;//Coulomb potential energy for the previous move
double *FNC_POT::Coulomb14_pot_old=0;//Coulomb potential energy for 1-4 interaction for the previous move
double *FNC_POT::r_pot;//r values for the tabulated potential files
double *FNC_POT::u_pot;//potential values for the tabulated potential files
double *FNC_POT::tab_cutoff;//cutoff for each used partial tabulated potential
double *FNC_POT::tabpot_dr;//reduced spacing of the tabulated potential for each used partial;
char  (*FNC_POT::tabpot_filename)[FILE_NAME_SIZE];//name of the file(s) containing the tabulated potential(s)
FNC_POT::TCalcvdW FNC_POT::CalcvdW=0;//pointer to the actual vdW function
FNC_POT::TCalcvdW FNC_POT::CalcvdW14=0;//pointer to the actual 1-4 vdW function
Threads *FNC_POT::thread_object=0;
std::vector<std::vector<int>> FNC_POT::fnc_partial_ind;//in case of automatic cutoff and fnc it will indicate for each partial in which fnc constraint it is involved in 


SimpleCfg *FNC_POT::config;
RunParams *FNC_POT::rundat;

//-----default constructor--------------------
FNC_POT::FNC_POT()
{
	int i;
	//Checking the main parameters 
	
	switch (RunParams::fnc)
	{
		
		case 1:
		case 2:
		case 3:
		{
			ndistances=ndist;//normal fnc
			break;
		}
		case 4:
		{
			ndistances=0;//will be set later for flexible molecules
			break;
		}
		default:
		{
			ndistances=0;
			break;
		}
	}
		
	//These two will be created with 0 length for all the instances, as later they will be resized, if needed
	SetArraysize(&neighbours2,0,"neighbours2","FNC_POT::FNC_POT");//array of neighbours2 indices
								//\nWARNING! the first atom has index 0 in the
								//program (see remark above)
	SetArraysize(&neighbours3,0,"neighbours3","FNC_POT::FNC_POT");//array of neighbours3 indices
								//\nWARNING! the first atom has index 0 in the
								//program (see remark above)
	SetArraysize(&dih_type,0,"dih_type","FNC_POT::FNC_POT");//array to the GROMACS types of the dihedrals
	//to show,that it was not yet created, and will not be created for all the instances
	//where needed, they will be created in GenerateExclusion or Init
	vdW_excl_numb=NULL;
	vdW_exclusions=NULL;
	vdW_conv=NULL;
	LJ_C6=NULL;
	LJ_CN=NULL;
	LJ14_C6=NULL;
	LJ14_CN=NULL;
	bond_r=NULL;
	bond_k=NULL;
	RBdih_C=NULL;
	charge_type=NULL;
	charge_group=NULL;
	charge_centre=NULL;
	atom_for_charge_gr=NULL;
	atom_for_charge_gr_finder=NULL;	
	//array of the number of neighbours	of each atom

	SetArraysize(&numbneigh,natoms,"numbneigh","FNC_POT::FNC_POT");
	
	//array of positions of atoms in the neighbours and consttypes array
	SetArraysize(&converter,natoms,"converter","FNC_POT::FNC_POT");
	
	//array of neighbours indices, \nWARNING! the first atom has index 0 in the program (see remark above)
	SetArraysize(&neighbours,ndistances,"neighbours","FNC_POT::FNC_POT");
	
	//array of constraint types, \nWARNING the constraint types start at 0 in the program instead of 1 in the fnc file
	SetArraysize(&consttypes,ndistances,"consttypes","FNC_POT::FNC_POT");
	
	//array of minimal distances 
	SetArraysize(&dmin,nconstraints,"dmin","FNC_POT::FNC_POT");
	//array of maximum distances
	SetArraysize(&dmax,nconstraints,"dmax","FNC_POT::FNC_POT");
	SetArraysize(&dminsq,nconstraints,"dminsq","FNC_POT::FNC_POT");
	SetArraysize(&dmaxsq,nconstraints,"dmaxsq","FNC_POT::FNC_POT");
	//setting all the array elements to 'zero' values
	for(i=0;i<nconstraints;i++)
	{
		*(dmin+i)=0.0;
		*(dmax+i)=0.0;
		*(dminsq+i)=0.0;
		*(dmaxsq+i)=0.0;
	}

};

void FNC_POT::Init(T_interaction_type pottype)
{

	pot_type=pottype;//type of the interaction (bond, angle dihedral)
		
	int i;
	nGRpartials=Topology::nGRtypes*(Topology::nGRtypes+1)/2;

	for (i=0;i<natoms;i++)
		numbneigh[i]=0;
		
	

	if (pot_type==BOND)
	{
		SetArraysize(&bond_k,Topology::nbond_types,"bond_k","FNC_POT::Init");
		SetArraysize(&bond_r,Topology::nbond_types,"bond_r","FNC_POT::Init");
		for (i=0;i<Topology::nbond_types;i++)
		{
			bond_r[i]=Topology::bond_r[i]/SimpleCfg::boxedge;//to get in in reduced units
			bond_k[i]=Topology::bond_k[i]*SimpleCfg::boxedge*SimpleCfg::boxedge;//to get the potential in kJ/mol, even if the distance is in reduced unit
		}
		if (RunParams::potential == 1)
		{
			LJ_pow=&RunParams::LJ_rep_N;
			SetArraysize(&LJ_C6,nGRpartials,"LJ_C6","FNC_POT::Init");
			SetArraysize(&LJ_CN,nGRpartials,"LJ_CN","FNC_POT::Init");
			SetArraysize(&LJ14_C6,nGRpartials,"LJ14_C6","FNC_POT::Init");
			SetArraysize(&LJ14_CN,nGRpartials,"LJ14_CN","FNC_POT::Init");
			for (i = 0 ; i < nGRpartials; i++)
			{
				LJ14_C6[i] = 0;
				LJ14_CN[i] = 0; 
			}
			//setting the arrays connected to the charges
			SetArraysize(&charge_type,SimpleCfg::ntotal+SimpleCfg::ntotvirtual,"charge_type","FNC_POT::Init");
			SetArraysize(&charge_group,SimpleCfg::ntotal + SimpleCfg::ntotvirtual,"charge_group","FNC_POT::Init");
			SetArraysize(&charge_centre,SimpleCfg::ntotal + SimpleCfg::ntotvirtual,"charge_centre","FNC_POT::Init");
			SetArraysize(&atom_for_charge_gr,SimpleCfg::ntotal + SimpleCfg::ntotvirtual,"atom_for charge_gr","FNC_POT::Init");
			SetArraysize(&atom_for_charge_gr_finder,Topology::nRMC_charge_centre,"atom_for charge_gr_finder","FNC_POT::Init");
			

			//This array, atom_for_charge_gr and its finder will be set by SimpleCfg::SetChargeGroupCentre
			for (i=0;i<SimpleCfg::ntotal + SimpleCfg::ntotvirtual;i++)
				atom_for_charge_gr[i]=-1;//not yet assigned
		
			//LJ_sigma has to be in reduced units, as distance
			if (RunParams::vdW_comb_rule == 0 || RunParams::vdW_comb_rule==1)//C(6) and C(N) is given, only reduction is necessary
			{
				for (i=0;i<nGRpartials;i++)
				{
					LJ_C6[i]=RunParams::vdW_pot1[i]/pow(SimpleCfg::boxedge,6.0);
					LJ_CN[i]=RunParams::vdW_pot2[i]/pow(SimpleCfg::boxedge,(double)(*LJ_pow));
					
				}
			}
			else
			{
				//sigma and epsion is given
				for (i=0;i<nGRpartials;i++)
				{
					LJ_C6[i]=4*RunParams::vdW_pot2[i]*pow(RunParams::vdW_pot1[i]/SimpleCfg::boxedge,6.0);
					LJ_CN[i]=4*RunParams::vdW_pot2[i]*pow(RunParams::vdW_pot1[i]/SimpleCfg::boxedge,(double)(*LJ_pow));
				}
			}
			//now the 1-4 interaction parameters only handles the LJ
			if (Topology::npairs > 0)
			{
				if (RunParams::vdW_comb_rule == 0 || RunParams::vdW_comb_rule == 1)//C(6) and C(N) is given, 
				{
					for (i = 0; i < nGRpartials; i++)
					{
						if (Topology::gen_pairs)//no known forcefield, but can appear with unknown
						{
							LJ14_C6[i] = LJ_C6[i] * RunParams::vdW14_fudge;//LJ 1-4 will be scaled from LJ, but it is already reduced
							LJ14_CN[i] = LJ_CN[i] * RunParams::vdW14_fudge;
						}
						else
						{
							//there is no implemented force field, only unknown
							//Only conversion to reduced units is needed
							LJ14_C6[i] = Topology::vdW14_pot1_partial[i] / pow(SimpleCfg::boxedge, 6.0);
							LJ14_CN[i] = Topology::vdW14_pot2_partial[i] / pow(SimpleCfg::boxedge, (double)(*LJ_pow));
						}

					}
				}
				else
				{
					//combination rule 2,3
					//LJ sigma and epsilon was given
					for (i = 0; i < nGRpartials; i++)
					{

						if (Topology::gen_pairs)//typically OPLSAA, Encads, Encadv, AMBER, CHARMM
						{
							LJ14_C6[i] = LJ_C6[i] * RunParams::vdW14_fudge;//LJ 1-4 will be scaled from LJ, but it is already reduced
							LJ14_CN[i] = LJ_CN[i] * RunParams::vdW14_fudge;
						}
						else
						{
							//sigma and epsilon is given, and reduction is needed, no known forcefield, can appear as unknown
							LJ14_C6[i] = 4 * Topology::vdW14_pot2_partial[i] * pow(Topology::vdW14_pot1_partial[i] / SimpleCfg::boxedge, 6.0);//this has to be redue
							LJ14_CN[i] = 4 * Topology::vdW14_pot2_partial[i] * pow(Topology::vdW14_pot1_partial[i] / SimpleCfg::boxedge, (double)(*LJ_pow));
						}
					}
				}
			}
		}
	}
	if (pot_type==ANGLE)
	{
		if (Topology::nangles>0)
		{
			SetArraysize(&angle_ang_rad,Topology::nangle_types,"angle_ang_rad","FNC_POT::Init");
			for (i=0;i<Topology::nangle_types;i++)
				angle_ang_rad[i]=Topology::angle_ang[i]/180*PI;
		}
	}
	if (pot_type==DIHEDRAL)
	{
		if (Topology::nperdihedrals>0)
		{
			SetArraysize(&perdih_ang_rad,Topology::nperdihedral_types,"perdih_ang_rad","FNC_POT::Init");
			for (i=0;i<Topology::nperdihedral_types;i++)
				perdih_ang_rad[i]=Topology::perdihedral_ang[i]/180*PI;
		}
		if (Topology::nharmdihedrals>0)
		{
			SetArraysize(&harmdih_ang_rad,Topology::nharmdihedral_types,"harmdih_ang_rad","FNC_POT::Init");
			for (i=0;i<Topology::nharmdihedral_types;i++)
				harmdih_ang_rad[i]=Topology::harmdihedral_ang[i]/180*PI;
		}
		if (Topology::nRBdihedrals>0)
		{
			SetArraysize(&RBdih_C,6,"RBdih_C","FNC_POT::Init");
			//setting the pointer array to the diherdral coefficients of the Topology object, to make the usage easier
			RBdih_C[0]=Topology::RBdihedral_C0;
			RBdih_C[1]=Topology::RBdihedral_C1;
			RBdih_C[2]=Topology::RBdihedral_C2;
			RBdih_C[3]=Topology::RBdihedral_C3;
			RBdih_C[4]=Topology::RBdihedral_C4;
			RBdih_C[5]=Topology::RBdihedral_C5;
		}

	}
	
}

	
//-----------reads the FNC_POT parameters from .fnc file----------------------------------
void FNC_POT::GetFNCParams()
{
	int i,k, dumb;
	ifstream file;

	SafeOpenTextFile(file,fncfilename);
	if (CheckFileState(file,"FNC_POT::Load",fncfilename)==0)
	{
		cout<<"Cannot run this way,exiting..."<<endl;
		CleanExit();
	};
	cout<<"\nLoading the FNC parameters from the "<<fncfilename<<" file"<<endl;

	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	nconstraints=ReadThisLine(file,1,1,"nconstraints", "FNC_POT::GetFNCParams",fncfilename);//number of FNC constraints
	SetArraysize(&nout_of_range,nconstraints,"nout_of_range","FNC_POT::GetFNCParams");
	SetArraysize(&dec_out_of_range,nconstraints,"dec_out_of_range","FNC_POT::GetFNCParams");
	for (i=0;i<nconstraints;i++)
	{
		nout_of_range[i]=0;
		dec_out_of_range[i]=0;
	}
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	natoms=ReadThisLine(file,1,1,"natoms", "FNC_POT::GetFNCParams",fncfilename);//number of central FNC atoms in the file
	SkipLine(file,fncfilename,1);
	k=0;//used to compute the total number of constrained distances
	for(i=1;i<=natoms;i++)//for each constrained atom except the last
	{
		dumb=ReadThisLine(file,3,1,"central atom index", "FNC_POT::GetFNCParams");//central atom index
		if (i<natoms)
		{
			dumb=ReadThisLine(file,1,1,"number of neighbours", "FNC_POT::GetFNCParams",fncfilename);//that's the number of neighbours

		
			if(dumb>0)
			{
				SkipLine(file,fncfilename,1);//neighbour(s) index(es)
				SkipLine(file,fncfilename,1);//constraint(s) type(s)
			}
		}
		else
			dumb=ReadThisLine(file,1,1,"number of neighbours", "FNC_POT::GetFNCParams");//that's the number of neighbours
		k+=dumb;
	}
	
	ndist=k;

	//Checking, whether loading was successful
	if (!CheckReadFileState(file,"FNC_POT::GetFNCParams",fncfilename))
		CleanExit();
	next_line_pos = -1;
	file.close();
};

//-------save-------
void FNC_POT::Save(const char *file_name) const
{
	//saves the FNC data to the file given in argument
	//the format file mimicks the actual FNC format
	//but it adds the number of distances and some comments to
	//improve readability: this is the 'save' format
	int i,j,k,offset=0;
	int *p_neigh, *p_type;//pointers to number of neighbours
							//and constraint type	
	int *p_neigh2,*p_neigh3,*p_GRtype;//pointers to number of neighbours2 and 3

	ofstream file;
	OpenFile(file,file_name,"FNC_POT::Save",0);//open file, check, whether it was successfully opened

	file<<"This is an object of type FNC_POT. File created by FNC_POT::save "<<endl;
	if (RunParams::fnc==4)
	{
		file<<"\n Number of "<<interaction_types[pot_type]<<"S"<<endl;
		offset=1;//all the indices will have to start with 1, otherwise the - sign will not have effect!
		switch (pot_type)
		{
			case BOND:
			{
				file<<J7<<Topology::nbond_types<<endl;
				for(i=0;i<Topology::nbond_types;i++)
					file<<J13<<Topology::bond_r[i];
				file<<"\t <-equilibrium distances (A)"<<endl;
				for(i=0;i<Topology::nbond_types;i++)
					file<<J13<<Topology::bond_k[i];
				file<<"\t <-force constant (kJ/mol/A2)"<<endl;
				break;
			}
			case ANGLE:
			{
				file<<J7<<Topology::nangle_types<<endl;
				for(i=0;i<Topology::nangle_types;i++)
					file<<J13<<*(Topology::angle_ang+i);
				file<<"\t <-equilibrium angle (degree)"<<endl;
				for(i=0;i<Topology::nangle_types;i++)
					file<<J13<<*(Topology::angle_k+i);
				file<<"\t <-force constant (kJ/mol/rad2)"<<endl;
				break;
			}
			case DIHEDRAL:
			{
				file<<J7<<Topology::nperdihedral_types<<endl;
				for(i=0;i<Topology::nperdihedral_types;i++)
					file<<J13<<*(Topology::perdihedral_ang+i);
				file<<"\t <-equilibrium angle (degree)"<<endl;
				for(i=0;i<Topology::nperdihedral_types;i++)
					file<<J13<<*(Topology::perdihedral_k+i);
				file<<"\t <-force constant (kJ/mol/rad2)"<<endl;
					
				file<<J7<<Topology::nharmdihedral_types<<endl;
				for(i=0;i<Topology::nharmdihedral_types;i++)
					file<<J13<<*(Topology::harmdihedral_ang+i);
				file<<"\t <-equilibrium angle (degree)"<<endl;
				for(i=0;i<Topology::nharmdihedral_types;i++)
					file<<J13<<*(Topology::harmdihedral_k+i);
				file<<"\t <-force constant (kJ/mol/rad2)"<<endl;

				file<<J7<<Topology::nRBdihedral_types<<endl;
				for(i=0;i<Topology::nRBdihedral_types;i++)
					file<<J13<<Topology::RBdihedral_C0[i];
				file<<"\t <-C0 coefficient (kJ/mol)"<<endl;
				for(i=0;i<Topology::nRBdihedral_types;i++)
					file<<J13<<Topology::RBdihedral_C1[i];
				file<<"\t <-C1 coefficient (kJ/mol)"<<endl;
				for(i=0;i<Topology::nRBdihedral_types;i++)
					file<<J13<<Topology::RBdihedral_C2[i];
				file<<"\t <-C2 coefficient (kJ/mol)"<<endl;
				for(i=0;i<Topology::nRBdihedral_types;i++)
					file<<J13<<Topology::RBdihedral_C3[i];
				file<<"\t <-C3 coefficient (kJ/mol)"<<endl;
				for(i=0;i<Topology::nRBdihedral_types;i++)
					file<<J13<<Topology::RBdihedral_C4[i];
				file<<"\t <-C4 coefficient (kJ/mol)"<<endl;
				for(i=0;i<Topology::nRBdihedral_types;i++)
					file<<J13<<Topology::RBdihedral_C5[i];
				file<<"\t <-C5 coefficient (kJ/mol)"<<endl;
				break;
			}
			file<<endl;

		}
		
	}
	else
	{
		file<<"\n Number of FNC constraints; FNC intervals "<<endl;
		file<<J7<<nconstraints<<endl;
		for(i=0;i<nconstraints;i++)
			file<<J13<<*(dmin+i);
		file<<"\t <-minimum distances"<<endl;
		for(i=0;i<nconstraints;i++)
			file<<J13<<*(dmax+i);
		file<<"\t <-maximum distances"<<endl;
		file<<endl;
	}
	file<<J7<<natoms<<J7<<ndistances<<\
	"\t <-number of constrained atoms & distances "<<endl;
	file<<endl;
	p_neigh=neighbours;
	p_neigh2=neighbours2;
	p_neigh3=neighbours3;
	p_GRtype=dih_type;
	p_type=consttypes;
	for(i=0;i<natoms;i++)
	{
		file<<J10<<i+offset<<J10<<*(numbneigh+i)<<\
		"\t <-atom index, number of neighbours"<<endl;
		k=*(numbneigh+i);
		for(j=0;j<k;j++)//saving the neighbours indices
		{
			file<<J10<<*p_neigh;
			p_neigh++;
		}
		file<<"\t <-neighbours indices"<<endl;
		if (RunParams::fnc==4 && (pot_type==(int)ANGLE || pot_type==(int)DIHEDRAL))
		{
			for(j=0;j<k;j++)//saving the neighbours indices
			{
				file<<J10<<*p_neigh2;
				p_neigh2++;
			}
			file<<"\t <-neighbours2 indices"<<endl;
		}
		if (RunParams::fnc==4 && pot_type==(int)DIHEDRAL)
		{
			for(j=0;j<k;j++)//saving the neighbours indices
			{
				file<<J10<<*p_neigh3;
				p_neigh3++;
			}
			file<<"\t <-neighbours3 indices"<<endl;
		}
		if (RunParams::fnc==4 && pot_type==(int)DIHEDRAL)
		{
			for(j=0;j<k;j++)//saving the GROMACS dihedral types
			{
				file<<J10<<*p_GRtype;
				p_GRtype++;
			}
			file<<"\t <-GROMACS dihedral types"<<endl;
		}
		for(j=0;j<k;j++)//saving the constraint types
		{
			file<<J10<<*p_type;
			p_type++;
		}
		file<<"\t <-constraint types"<<endl;
	}
	file.close();
}

//---------------rawsave---------
void FNC_POT::RawSave(const char *file_name) const
{//saves the FNC data in the orginal RMCA .fnc format

	int i,j,k;
	int *p_neigh, *p_type;//pointers to number of neighbours
							//and constraint type
	ofstream file;
	OpenFile(file,file_name,"FNC_POT::RawSave",0);//open file, check, whether it was successfully opened
	
	file<<"This is a .fnc file created by FNC_POT::rawsave "<<endl;
	file<<"\n No. of possible rmin-rmax pairs: "<<endl;
	file<<J7<<nconstraints<<endl;
	file.precision(9);
	file.setf(ios::fixed, ios::floatfield); 
	file.setf(ios::right, ios::adjustfield);
	for(i=0;i<nconstraints;i++)
		file<<J13<<*(dmin+i);
	file<<endl;
	for(i=0;i<nconstraints;i++)
		file<<J13<<*(dmax+i);
	file<<"\n "<<endl;
	file<<J7<<natoms<<"\n "<<endl;
	file.precision(6);
	file.unsetf(ios::fixed);
	file.unsetf(ios::right);

	p_neigh=neighbours;
	p_type=consttypes;
	for(i=0;i<natoms;i++)
	{
		file<<J10<<i+1<<J10<<*(numbneigh+i)<<endl;
		k=*(numbneigh+i);
		for(j=0;j<k;j++)//saving the neighbours indices
						//with re-conversion to file type (index starting at 1)
		{
			file<<J10<<*p_neigh+1;
			p_neigh++;
		}	file<<endl;
		for(j=0;j<k;j++)//saving the constraint types
					//with re-conversion to file type (index starting at 1)
		{
			file<<J4<<*p_type+1;
			p_type++;
		}
	file<<endl;
	}
	file.close();
}

//------------load----------
void FNC_POT::Load()
{//loads from a 'save' file (checks for compatibility of 
	//dimensions nconstraints, natoms, ndistances

	ifstream file;

	SafeOpenTextFile(file,fncfilename);
	if (CheckFileState(file,"FNC_POT::Load",fncfilename)==0)
	{
		cerr<<"There is no open file to load the FNC_POT object from, exiting..."<<endl;
		CleanExit();
	}
		
	int i,j,k,temp;
	int *p_neighb, *p_type;//pointers to neighbour index
							// and constraint types
	
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	temp=ReadThisLine(file,1,1,"nconstraints", "FNC_POT::Load");//number of constraints
	if(nconstraints!=temp)//check compatibility
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "FNC_POT::load-> the number of constraints in the " << fncfilename << " file is "<<temp<< endl;
		cout<<"which is not of the same as in the program: "<<nconstraints<<endl;
		cout<<"\tExiting..."<<endl;
		CleanExit();
	} 
	
	for(i=0;i<nconstraints;i++)
	{
		if (i<nconstraints-1)
			*(dmin+i)=ReadThisLine(file,3,1.0,"dmin", "FNC_POT::Load");//minimum distances
		else
			*(dmin+i)=ReadThisLine(file,1,1.0,"dmin", "FNC_POT::Load","fnc3.txt");
		*(dminsq+i)=*(dmin+i)* (*(dmin+i));
	}
	
	for(i=0;i<nconstraints;i++)
	{
		if (i<nconstraints-1)
			*(dmax+i)=ReadThisLine(file,3,1.0,"dmax", "FNC_POT::Load");//maximum distances
		else
			*(dmax+i)=ReadThisLine(file,1,1.0,"dmax", "FNC_POT::Load","fnc3.txt");
		*(dmaxsq+i)=*(dmax+i)*(*(dmax+i));
		if (*(dmax+i)< *(dmin+i))
		{
			cout << "\n*****ERROR*****" << endl;
			cout<<"The FNC dmax "<<*(dmax+i)<<" A is smaller than the FNC dmin "<<*(dmin+i)<<" A for the "<<i+1<<". constraint!"<<endl;
			cout<<"Correct it and start again! Cannot run this way, CleanExit!"<<endl;
			CleanExit();
		}
	}
	SkipLine(file,fncfilename,1);
	temp=ReadThisLine(file,3,1,"natoms", "FNC_POT::Load",fncfilename);;//number of constrained atoms indicated
	
	if(natoms!=temp)//check compatibility
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "FNC_POT::Load-> the number of atoms in the " << fncfilename << " file is " << temp << endl;
		cout << "which is not of the same as in the program: "<<natoms<<endl;
		cout<<"Exiting..."<<endl;
		CleanExit();
	} 
	temp=ReadThisLine(file,1,1,"ndistances", "FNC_POT::Load","fnc3.txt");;//number of constrained distances indicated 
	//in the FNC file
	
	if(ndistances!=temp)//check compatibility
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "FNC_POT::load-> the number of constrained distances in the " << fncfilename << " file is " << temp << endl;
		cout << "which is not of the same as in the program: " <<ndistances<<endl;
		cout<<"Exiting..."<<endl;
		CleanExit();
	} 
			
	SkipLine(file,fncfilename,1);

	k=0;//that is the converter array current value
	p_neighb=neighbours;
	p_type=consttypes;//positioning the pointer
	for(i=0;i<natoms;i++)//for all atoms in the .FNC file
	{
		*(converter+i)=k;
		temp=ReadThisLine(file,3,1,"central index", "FNC_POT::Load");//atom cfg index in the fnc file
		temp=ReadThisLine(file,1,1,"number of neighbours", "FNC_POT::Load","fnc3.txt");//number of neighbours
		k+=temp;//update the converter content
		*(numbneigh+i)=temp;//fill the array for number of neighbours
		
		for(j=0;j<temp;j++)
		{
			//filling the neighb. indices array
			if (j<temp-1)
				*p_neighb=ReadThisLine(file,3,1,"neighbour index", "FNC_POT::Load");
			else
				*p_neighb=ReadThisLine(file,1,1,"neighbour index", "FNC_POT::Load","fnc3.txt");

			p_neighb++;
		}
		
		for(j=0;j<temp;j++)
		{	
			//filling the constr. type array
			if (j<temp-1)
				*p_type=ReadThisLine(file,3,1,"neighbour type", "FNC_POT::Load");
			else
			{
				if (i<natoms-1)
					*p_type=ReadThisLine(file,1,1,"neighbour type", "FNC_POT::Load","fnc3.txt");
				else
					*p_type=ReadThisLine(file,1,1,"neighbour type", "FNC_POT::Load");
			}
			p_type++;
		}
	}
	if (!CheckReadFileState(file,"FNC_POT::Load",fncfilename))
		CleanExit();
	next_line_pos = -1;
	file.close();
}

//--------------Rawload--------------
void FNC_POT::RawLoad()
{
	//loads from a .fnc format file
	//\nWARNING ! compatibility of the number of distances
	// is NOT checked
	
	ifstream file;

	SafeOpenTextFile(file,fncfilename);
	if (CheckFileState(file,"FNC_POT::RawLoad",fncfilename)==0)
	{
		cerr<<"There is no open file to load the FNC_POT object from, exiting... "<<endl;
		CleanExit();
	}
	
	int i,j,k,temp;
	int *p_neighb, *p_type;//pointers to neighbour index
							// and constraint types
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	SkipLine(file,fncfilename,1);
	temp=ReadThisLine(file,1,1,"nconstraints", "FNC_POT::Rawload");//number of constraints
	if(nconstraints!=temp)//check compatibility
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "FNC_POT::RawLoad-> the number of constraints in the " << fncfilename << " file is " << temp << endl;
		cout << "which is not of the same as in the program: " << nconstraints << endl;
		cout << "\tExiting..." << endl;
		CleanExit();
	} 
	

	for(i=0;i<nconstraints;i++)
	{
		if (i<nconstraints-1)
			*(dmin+i)=ReadThisLine(file,3,1.0,"dmin", "FNC_POT::Rawload");//minimum distances
		else
			*(dmin+i)=ReadThisLine(file,1,1.0,"dmin", "FNC_POT::Rawload",fncfilename);
		*(dminsq+i)=*(dmin+i)* (*(dmin+i));
	}
	
	for(i=0;i<nconstraints;i++)
	{
		if (i<nconstraints-1)
			*(dmax+i)=ReadThisLine(file,3,1.0,"dmax", "FNC_POT::Rawload");//maximum distances
		else
			*(dmax+i)=ReadThisLine(file,1,1.0,"dmax", "FNC_POT::Rawload",fncfilename);
		*(dmaxsq+i)=*(dmax+i)*(*(dmax+i));
		if (*(dmax+i)< *(dmin+i))
		{
			cout << "\n*****ERROR*****" << endl;
			cout<<"The FNC dmax "<<*(dmax+i)<<" A is smaller than the FNC dmin "<<*(dmin+i)<<" A for the "<<i+1<<". constraint!"<<endl;
			cout<<"Correct it and start again! Cannot run this way, CleanExit!"<<endl;
			CleanExit();
		}
	}
	SkipLine(file,fncfilename,1);
	temp=ReadThisLine(file,1,1,"natoms", "FNC_POT::Rawload",fncfilename);;//number of constrained atoms indicated

	//in the FNC file
	if(natoms!=temp)//check compatibility
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "FNC_POT::Load-> the number of atoms in the " << fncfilename << " file is " << temp << endl;
		cout << "which is not of the same as in the program: " << natoms << endl;
		cout << "Exiting..." << endl;
		CleanExit();
	} 
					
	SkipLine(file,fncfilename,1);
	
	k=0;//that is the converter current content
	p_neighb=neighbours;
	p_type=consttypes;//positioning the pointer

	for(i=0;i<natoms;i++)//for all atoms in the .FNC file

	{
		*(converter+i)=k;
		temp=ReadThisLine(file,3,1,"central index", "FNC_POT::Rawload");//atom cfg index in the fnc file
		temp=ReadThisLine(file,1,1,"number of neighbours", "FNC_POT::Rawload",fncfilename);;//number of neighbours
		
		k+=temp;//update the converter content
		*(numbneigh+i)=temp;//fill the array for number of neighbours


		for(j=0;j<temp;j++)
		{
			//filling the neighb. indices array

			if (j<temp-1)
				*p_neighb=ReadThisLine(file,3,1,"neighbour index", "FNC_POT::Rawload")-1;
			else
				*p_neighb=ReadThisLine(file,1,1,"neighbour index", "FNC_POT::Rawload",fncfilename)-1;
			p_neighb++;
		}

		
		for(j=0;j<temp;j++)
		{

			//filling the constr. type array
			if (j<temp-1)
				*p_type++=ReadThisLine(file,3,1,"neighbour type", "FNC_POT::Rawload")-1;
			else
			{
				if (i<natoms-1)
					*p_type++=ReadThisLine(file,1,1,"neighbour type", "FNC_POT::Rawload",fncfilename)-1;
				else
					*p_type=ReadThisLine(file,1,1,"neighbour type", "FNC_POT::Rawload")-1;
			}
			
		}

	}

	//finally unscaling the FNC distances squared
	for(i=0;i<nconstraints;i++)
	{
		*(dmaxsq+i)/=(RunParams::boxedge*RunParams::boxedge);
		*(dminsq+i)/=(RunParams::boxedge*RunParams::boxedge);
	}
	//Checking, whether loading was successful
	if (!CheckReadFileState(file,"FNC_POT::RawLoad",fncfilename))
		CleanExit();
	next_line_pos = -1;
	file.close();
}

//--------------Check the satisfaction of the FNC constraint-------------
//Determine, whether there is a conflict between the cutoffs and the FNC constraint

int FNC_POT::CheckFNC(double *cutsq)//return value is the conflict indicator
{
	int type1,type2,i,j,partialind;//cycle variable 
	int ntypes;
	int array_ind,const_type,neighbour_ind;
	int conflict;//1 if the one of the FNC constraint min. distances is below cutoff, 0 otherwise
	double x_central,y_central,z_central;//coordinates of the central atom
	double dx,dy,dz;//differences in the coordinates
	double *dmin_conf, *dmax_conf;//minimum, maximum distances occured for the constraints
	double dcutsq;//squared, reduced cutoff
	double dsquare;//squared, reduced distance of two atoms
	double tol=0;//3.464102e-7;//this is the tolerance for reduced, squared distances
	SimpleCfg &conf=*config;
	//Initialisation
	conflict=0;//Sets the conflict indicator to default no conflict
	if (RunParams::auto_cutoff&1)//there is autocutoff datermination
	{
		for (i = 0; i < SimpleCfg::npartials; i++)
			fnc_partial_ind.push_back({});
	}
	ntypes=conf.ntypes;
	if (RunParams::fnc==2)
		i=nconstraints;//to collect the smallest and largest values out of range
	else
		i=0;
	SetArraysize(&dmin_conf,i,"dmin_conf","FNC_POT::CheckFNC");
	SetArraysize(&dmax_conf,i,"dmax_conf","FNC_POT::CheckFNC");
	if (RunParams::fnc==2)
	{
		for (i=0;i<nconstraints;i++)
		{
			dmin_conf[i]=1000000;//default, no smaller value than the one given in the fnc file for rmin occured
			dmax_conf[i]=-1;//default, no larger value than the one given in the fnc file for rmax occured
		}
	}

	for (type1=0;type1<ntypes;type1++)//for all the atom types
	{
		for (i=conf.cumul[type1];i<conf.cumul[type1+1];i++)//for all the central atoms of type1
		{
			//the index of the first element of the neighbours and consttypes arrays
			//for the i-th atom
			array_ind=converter[i];
			x_central=conf.positions[3*i];//x ccordinate of the central atom
			y_central=conf.positions[3*i+1];//y ccordinate of the central atom
			z_central=conf.positions[3*i+2];//z ccordinate of the central atom
			for (j=0;j<numbneigh[i];j++)//for all the neighbours of the central atom 
			{
				const_type=consttypes[array_ind+j];//constraint type
				neighbour_ind=neighbours[array_ind+j];//index of the neighbour

				if (conflict==0)//if it is still necessary to check for conflict
				{
					//find the type of the neighbour atom
					type2=0;
					while (neighbour_ind>=conf.cumul[type2+1])
						type2++;
					if (type1 <= type2)//type1 is the central atom
						partialind=type1 * ntypes - (type1 * (type1 + 1) / 2) + type2;
					else//type2 is the central atom
						partialind = type2 * ntypes - (type2 * (type2 + 1) / 2) + type1;
					int found = 0;
					if (RunParams::auto_cutoff & 1)//confilct stays 0, set fnc
					{
						//first see, if this fnc constrint type is already present among the values of this partial
						for (unsigned int iv = 0; iv < fnc_partial_ind[partialind].size(); iv++)
							if (fnc_partial_ind[partialind][iv] == const_type)
								found = 1;
						if (found == 0)
							fnc_partial_ind[partialind].push_back(const_type);
					}
					else
					{
						dcutsq = cutsq[partialind];
						if (dminsq[const_type] < dcutsq)//FNC minimum distance is below cutoff
						{
							conflict = 1;//conflict between FNC and cutoff
							cout << "\nNOTE(" << ++note << "):FNC_POT::CheckFNC" << endl;
							cout << "\tFNC minimum distance is below cutoff at least for pair type1=" << type1 << " type2=" << type2 << endl;

						}
					}
				}//end of checking for conflict

				//Calculating the distance of the central and neighbour atoms
				dx=x_central-conf.positions[3*neighbour_ind];
				dy=y_central-conf.positions[3*neighbour_ind+1];
				dz=z_central-conf.positions[3*neighbour_ind+2];

				//Take into account the periodical boundary conditions
				if (dx>1) dx-=2;
				else if (dx<-1) dx+=2;
				if (dy>1) dy-=2;
				else if (dy<-1) dy+=2;
				if (dz>1) dz-=2;
				else if (dz<-1) dz+=2;

				dsquare=dx*dx+dy*dy+dz*dz;//squared, reduced distance of the atoms

				//Check, whether the FNC constraint is satisfied
				if (dsquare+tol<dminsq[const_type] || dsquare>dmaxsq[const_type]+tol)
				{
					//The given atoms do not satisfy the FNC
					sum_out_of_range++;
					nout_of_range[const_type]++;
					if (RunParams::fnc!=3)//do not write, as most probably a lot will be out of range
					{
						cout << "\nWARNING(" << ++warn << "): FNC_POT::CheckFNC"<<endl;
						cout<<"\tAtoms "<<i+1<<" and "<<neighbour_ind+1<<" do not satisfy the FNC constraint number "<<const_type+1<<" !"<<endl;
						cout<<"\t(Indexes start at 1!)"<<endl;
					}

					if (dsquare>dmaxsq[const_type])
					{
						if (RunParams::fnc!=3)//do not write, as most probably a lot will be out of range
							cout << "\nWARNING(" << ++warn << "): Distance: "<<sqrt(dsquare)*conf.boxedge<<\
							" greater, than maximum distance: "<<dmax[const_type]<<" Angstrom"<<endl;
						if (RunParams::fnc==2 && dsquare>dmax_conf[const_type])//in case of fnc option 2 the constraint range has to be widened
							dmax_conf[const_type]=dsquare;//this is larger than the sofar largest distance
					}
					else
					{
						if (RunParams::fnc!=3)//do not write, as most probably a lot will be out of range
							cout << "\nWARNING(" << ++warn << "): Distance: "<<sqrt(dsquare)*conf.boxedge<<\
							" smaller, than minimum distance: "<<dmin[const_type]<<" Angstrom"<<endl;
						if (RunParams::fnc==2 && dsquare<dmin_conf[const_type])//in case of fnc option 2 the constraint range has to be widened
							dmin_conf[const_type]=dsquare;//this is smaller, than the sofar smallest distance
					}
					cout.precision(6);
					

				}
			}//end of cycling through the neighbours
		}//end of cycling through the central atoms of type1
	}//end of cycle type1

	if (sum_out_of_range)
		cout<<"-------------------------------------------------------------"<<endl;

	for (i=0;i<nconstraints;i++)
	{
		if (nout_of_range[i]&1)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Something is wrong with the FNC, as odd number of FNC pairs are out of range for the " << i + 1 << ". constraint," << endl;
			cout<<"it should be an even number as a pairs should be considered twice, both being central atom. Check the fnc file!"<<endl;
			CleanExit();
		}
	
		nout_of_range[i]/=2;
		if (nout_of_range[i]>0)
		{
		
			if (RunParams::fnc == 1)
				cout << "\n*****ERROR*****" << endl;
			else
				cout << "\nWARNING(" << ++warn << "): ";
			cout<<"FNC constraint "<<i+1<<" is not satisfied for "<<nout_of_range[i]<<" atom pair(s)!"<<endl;
			if (RunParams::fnc==1)
			{
				cout<<"Cannot run this way, exiting..."<<endl;
				CleanExit();
			}
		}
	}
	if (sum_out_of_range)
	{
		if (RunParams::fnc==2)//reset the fnc ranges, if necessary
		{
			cout << "\nWARNING(" << ++warn << "): FNC option=2 is applied, so if there are distances out of the FNC range," << endl;
			cout<<"\tthan the minimum and maximum fnc distance will be reset to accomodate the out-of range values"<<endl;
			for (i=0;i<nconstraints;i++)
			{ 
				if (fabs(dmin_conf[i]-1000000)>TOLERANCE)//there was smaller value, than the minimum value given in the *.fnc file
				{
					cout<<"\tThe minimum distance for FNC constraint "<<i+1<<" (index starting with 1) was reset from "<<dmin[i]<<" to "<<sqrt(dmin_conf[i]-tol)*RunParams::boxedge<<" A!"<<endl; 
					dmin[i]=sqrt(dmin_conf[i]-tol)*RunParams::boxedge;
					dminsq[i]=dmin_conf[i]-tol;

				}
				if (dmax_conf[i]>0)//there was larger value, than the maximum value given in the *.fnc file
				{
					cout<<"\the maximum distance for FNC constraint "<<i+1<<" (index starting with 1) was reset from "<<dmax[i]<<" to "<<sqrt(dmax_conf[i]+tol)*RunParams::boxedge<<" A!"<<endl; 
					dmax[i]=sqrt(dmax_conf[i]+tol)*RunParams::boxedge;
					dmaxsq[i]=dmax_conf[i]+tol;
				}
			}

		}
		if (RunParams::fnc==3)
		{
			cout << "\nNOTE(" << ++note << "): FNC option=3 is applied, so RMC will try to eliminate the distances out of the FNC range"<<endl;
		}

	}
	if (RunParams::auto_cutoff & 1)
		conflict = 1;//cutoff is not known, set conflict to 1, so during cutoff determination in HistoSet it will be checked, if a pair is an fnc pair
			//as those does not count during cutoff determination
	return(conflict);//returns the conflict indicator
};

//----------Check, whether the FNC is satisfied after each move----------------
int FNC_POT::CheckFNCChange(Move &move)
{
	int imoved,array_ind,nnb,acc,ind,i;
	int *p_const_type,*p_neighbour_ind,*p_ind;
	double dx,dy,dz,d1,d2,dsq,dsq_old;
	double *p_newpos,*p_oldpos,*p_pos_neigh;
	SimpleCfg &conf=*config;

	acc=1;//set the indicator of acceptance to true
	p_newpos=move.newpos;//points to the (new) coordinates of the central atom
	p_oldpos=move.oldpos;//points to the (old) coordinates of the central atom

	for (i=0;i<nconstraints;i++)
		dec_out_of_range[i]=0;
	for (imoved=0;imoved<move.nmoved;imoved++)//for all atoms to be moved in the Move (no swaps are alowed in case of FNC)
	{
		//index of the moved atom
		array_ind=converter[move.indices[imoved]];
		p_const_type=consttypes+array_ind;//constraint type
		p_neighbour_ind=neighbours+array_ind;//index of the neighbour
		
		for (nnb=0;nnb<numbneigh[move.indices[imoved]];nnb++)//for each neighbour
		{
			d1=dminsq[*p_const_type];//reduced minimum distance squared
			d2=dmaxsq[*p_const_type];//reduced maximum distance squared

			p_pos_neigh=conf.positions+*p_neighbour_ind*3;//totconf.positions contains the new coordinates of the moved atoms
			
			dx=*p_pos_neigh++ -*p_newpos;
			dy=*p_pos_neigh++ -*(p_newpos+1);
			dz=*p_pos_neigh   -*(p_newpos+2);
			
			//take into account the toroidal structure of the cell
			if (dx>=1)
				dx-=2;
			else
				if(dx<-1)
					dx+=2;
			if (dy>=1)
				dy-=2;
			else
				if(dy<-1)
					dy+=2;
			if (dz>=1)
				dz-=2;
			else
				if (dz<-1)
					dz+=2;
			
			dsq=dx*dx+dy*dy+dz*dz;//distance squared
			
			if (nout_of_range[*p_const_type]>0)//there are still FNC pairs out of range, see, if this is one
			{
				//calulating the old distance
				p_pos_neigh-=2;//put the neighbour coordinate to x
				p_ind=move.indices;
				//see if the neighbour is moved as well
				for (ind=0;ind<move.nmoved;ind++)
				{
					if (*p_ind==*p_neighbour_ind)
					{
						p_pos_neigh=p_oldpos+3* ind;
						break;
					}
				}
					
				dx=*p_pos_neigh++ -*p_oldpos;
				dy=*p_pos_neigh++ -*(p_oldpos+1);
				dz=*p_pos_neigh   -*(p_oldpos+2);
		
				//take into account the toroidal structure of the cell
				if (dx>=1)
					dx-=2;
				else
					if(dx<-1)
						dx+=2;
				if (dy>=1)
					dy-=2;
				else
					if(dy<-1)
						dy+=2;
				if (dz>=1)
					dz-=2;
				else
					if (dz<-1)
						dz+=2;
				
				dsq_old=dx*dx+dy*dy+dz*dz;//old distance squared
				
				if (dsq_old>d2 || dsq_old<d1)//the old distance is out of range
				{
					if (dsq>d2 || dsq<d1)//the new is out of range
					{
						if ((dsq_old>d2 && dsq>d2 && dsq_old<dsq) || (dsq_old<d1 && dsq<d1 && dsq_old>dsq))
						{
							//the new distance is out of range, further than before
							acc=0;
							break;
						}
					}
					else
					{
						dec_out_of_range[*p_const_type]++;//this pair is not out of range any more, decrease the counter
						if (nout_of_range[*p_const_type]==0)
							cout << "\nNOTE(" << ++note << "): There are no fnc pairs out of range any more!!!"<<endl;
					}
				}
				else
				{
					if (dsq>d2 || dsq<d1)//this is a normal fnc pair, out of range cannot be accepted
					{	
						acc=0;
						break;
					}				
				}

			}
			else //there are no fnc pairs out of range, if the new distance is out of range, cannot be accepted
			{
				if (dsq>d2 || dsq<d1)
				{
					acc=0;
					break;
				}				
			}
			p_neighbour_ind++;//go to next neighbour
			p_const_type++;//go to next constraint type
		}//next neighbour
		if (acc==0)
			break;
		p_newpos+=3;//points to the new x-coordinate of the next central atom
		p_oldpos+=3;//points to the old x-coordinate of the next central atom
	}//next central atom
	return (acc);//return, whether the FNC is satisfied, and the move is acceptable
};


//As when the constructor is called, the number of neighbours is not known, these arrays are created with 0 length
//now they are recreated 
void FNC_POT::ResizeNeighArrays()
{
	int i=ndistances;
	
	ndistances=0;//reset it to 0 to have the same value, when the arays were originally created
	ResizeArray(&ndistances,i,&neighbours,"neighbours","FNC_POT::ResizeNeighArrays");
	ndistances=0;//reset it to 0 to have the same value, when the arays were originally created
	ResizeArray(&ndistances,i,&consttypes,"consttypes","FNC_POT::ResizeNeighArrays");
	if (RunParams::fnc==4)
	{
		if (Topology::nangles>0)
		{
			ndistances=0;//reset it to 0 to have the same value, when the arays were originally created
			ResizeArray(&ndistances,i,&neighbours2,"neighbours2","FNC_POT::ResizeNeighArrays");
		}
		if (Topology::ndihedrals>0)
		{
			ndistances=0;//reset it to 0 to have the same value, when the arays were originally created
			ResizeArray(&ndistances,i,&neighbours3,"neighbours3","FNC_POT::ResizeNeighArrays");
			ndistances=0;//reset it to 0 to have the same value, when the arays were originally created
			ResizeArray(&ndistances,i,&dih_type,"dih_type","FNC_POT::ResizeNeighArrays");
		}
	}
};

//The non-bonded interactions are handled through the exclusion list, which tells, that the vdW should not be calculated for the excluded pairs,
//if fudge for the given pair is 1.0. 1-4 interactions are handled by it the way, that 1-4 pairs will be included in the exclusion list, but with
//0<fudge<1.0, and only fudge*pot is subtracted, leaving (1-fudge)*pot for the vdW interaction (so fudge here will be 1-fudge_GROMACS!!!
//The nexclusions has to be 3 to exclude the 1-4 pairs from the list, and the 1-4 pairs to calculate for should be given at the [pairs] section of the topology!
//The list has to be atom-based, and it is contained in the pfnc[0] object (which also handles the BONDS), so for each atom the number of the 
//excluded atoms are contained in the vdW_excl_num list, and the indices of the exluded atoms is the vdW_exclusions array, and the vdW_Conv array 
//shows the first element belonging to a given atom in the vdW_exclusions array.
//This way each excluded pair is represented twice, once in the list of both participating atoms. This is necessary
//in RMC, as during the loop only the interactions of the moved atom(s) has to be updated, so atom based lists has to be used.
//Exclusion list will be cheked during the initial histogram calculation, and vdW will only be calculated for the pairs not on the exclusion std::list. 
//However, during the histogram change calculation vdW interaction energy will be calculated for all the atom pairs (this way we avoid the continuous checking), 
//and then the interactions for the excluded pairs are recalculated and subtracted, as during the loop there will be only a few excluded pairs for the
//moved atoms.
//To create a positive list (containing which atom pairs the vdW should be calculated for) seemed unfeasible, as this list would be huge.
//This routine generates the exclusion list for the atom pairs exluded from the non-bonded interaction calculation.
//All the atoms being nexclusions bonds away from each other will be excluded from each others's list, and additionally atoms from the 
// [ exclusion ] directive will be added, if there are any, which is not coming from the grapph search. This is necessary because of the virtual sites! 
//The exclusion will contain the indices of the atoms beginning with 1 to be able to have a sign (similarly to neighbour2 array)
void FNC_POT::GenerateExclusion(Topology &topology)
{
	int i,il,vi,iex1,iex2,imoltype, iat, imol,ipair,pair14=0,iexcl,pair_offset,iloffset;
	int ind1=0;//RMC_index of the atom in the first instance of the molecule type
	int ind2;//offset of the RMC index compared to the first molecule of this molecule type
	int ind_excl_centre;//RMC index of the exclusion centre
	int ind_excluded=0;//RMC index of the excluded atom
	int count;//number of exclusion from graph search
	int max_exclusions=0;//maximum number of exclusions for an atom
	int type1, type2, partialind;//RMC types and partial index for the 1-4 pairs
	int *temp_excl;//temporary array to hold the exclusions for an atom, as the found exclusions has to be checked, as in case of ring the same atom can be reached both way
	int *temp_delta;//to store temporarily the delta indices for the exclusions of an atom

	//The atoms are according to the RMC order in the following two arrays, as appear in the *.cfg file!
	SetArraysize(&temp_excl, INPUT_SIZE, "temp_excl", "FNC_POT::GenerateExclusion");//the number of exlusions for an the atom, cannot be more than INPUT_SIZE
	SetArraysize(&vdW_excl_numb,natoms+nvirtuals,"vdW_excl_numb","FNC_POT::GenerateExclusion");//the number of exlusions for all the atoms
	SetArraysize(&vdW_conv,natoms+nvirtuals,"vdW_conv","FNC_POT::GenerateExclusion");//the converter array for the exclusions
	//Determining the number of exclusions for each atom of each molecule type
	for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
	{
		for (iat=0;iat<Topology::natoms_per_type[imoltype];iat++)
		{
			//Topology::first_index[Topology::cumul_atom[imoltype]+iat]is the RMC_index of the first occurence of this atomtype, in this molecule
			//	if this type would be the first in the system
			//Topology::RMCindex_offset[imoltype] is the offset of the RMC index to take into account that might not be the first moleculetype
			if (Topology::first_index[Topology::cumul_atoms[imoltype] + iat]!=0)
				ind1 = Topology::first_index[Topology::cumul_atoms[imoltype] + iat] - 1 + Topology::RMCindex_offset[imoltype];
			else
			{
				//for virtual site has to be done differently
				for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
				{
					if (iat+1 == Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].atom_index)
					{
						ind1 = SimpleCfg::cumul[SimpleCfg::ntypes + Topology::cumul_virtuals[imoltype] + vi];
						break;//vi cycle
					}
				}
			}
			
			//graph search only needed for real atoms
			if (Topology::nexclusions_per_type[imoltype]>0 && Topology::first_index[Topology::cumul_atoms[imoltype] + iat] != 0)
				vdW_excl_numb[ind1]= GetExclusions(ind1,-1,Topology::nexclusions_per_type[imoltype]-1,temp_excl,temp_excl);
			else
				vdW_excl_numb[ind1]=0;
			
			//now determine the additional number of exclusions for this atom which are given in the exclusion section
			//only the atoms not found by the graph search will be determined, the others will be discarded
			for (il = 0; il < Topology::nexclusionlines_per_type[imoltype]; il++)//go through the exclusion lines of this moltype
			{
				iloffset = Topology::cumul_exclusionlines[imoltype] + il;//the index of the exclusion line
				ind_excl_centre = Topology::exclusion_centre[iloffset];//GROMACS index
				if (ind_excl_centre == -1)
					continue;//already removed during exclusions check
				if (iat + 1 == ind_excl_centre)
				{
					for (iex1 = 0; iex1 < Topology::nexcluded_neighs[iloffset]; iex1++)
					{
						//determine the RMC_index +1 (index start with 1)
						//for virtual sites first_index cannot be used
						if (Topology::first_index[Topology::cumul_atoms[imoltype] + Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] - 1] == 0)
						{
							for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
							{
								if (Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1]  == Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].atom_index)
								{
									ind_excluded = SimpleCfg::cumul[SimpleCfg::ntypes + Topology::cumul_virtuals[imoltype] + vi]+1;
									break;//vi cycle
								}
							}
						}
						else
							ind_excluded = Topology::first_index[Topology::cumul_atoms[imoltype] + Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] - 1] + Topology::RMCindex_offset[imoltype];
						for (iex2 = 0; iex2 < vdW_excl_numb[ind1]; iex2++)//these are the exclusion found by the graph search
						{
							if (temp_excl[iex2] == ind_excluded)
							{
								//remove it from the excluded list in the exlusion line by setting it to -1
								Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] = -1;
								break;//no need to search for more in temp_excl array for ind_excluded
							}
						}
					}
					//now contract the array by filling the empty places and update nexcluded_neighs
					i = 0;
					for (iex1 = 0; iex1 < Topology::nexcluded_neighs[iloffset]; iex1++)
					{
						if (Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] != -1)
						{
							Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + i ]= Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1];
							i++;
						}
					}
					Topology::nexcluded_neighs[iloffset] = i;
					if (Topology::nexcluded_neighs[iloffset] == 0)
						Topology::exclusion_centre[iloffset] = -1;//no exclusions for this atom, remove it
					//the empty exclusion line will remain, as it would be difficult to remove it 
					vdW_excl_numb[ind1] += Topology::nexcluded_neighs[Topology::cumul_exclusionlines[imoltype] + il];
				}//this exclusion line is for the ind1 atom
			}
				
			for (imol=1;imol<Topology::nmol_per_type[imoltype];imol++)//for the other molecules of the same type
			{
				//imol*delta_index[cumul_atom[imoltype]+iat] is the offset of this occurence's index relative to the first
				ind2=imol*Topology::delta_index[Topology::cumul_atoms[imoltype]+iat];
				vdW_excl_numb[ind1+ind2]= vdW_excl_numb[ind1];//number of exclusion for this atom
		
			}
			nexclusions+=vdW_excl_numb[ind1]*Topology::nmol_per_type[imoltype];//to determine the total  number of exclusions
			if (vdW_excl_numb[ind1] > max_exclusions)
				max_exclusions = vdW_excl_numb[ind1];
		}

	}
	delete [] temp_excl;
	//Now we can create the vdW_Conv array to find the first element for each atom
	SetArraysize(&vdW_exclusions,nexclusions,"vdW_exclusions","FNC_POT::GenerateExclusion");//the number of exlusions for all the atoms
	
	i=0;
	for (iat=0;iat<natoms+nvirtuals;iat++)//the Gromacs index of the atom in the molecule -1
	{
		vdW_conv[iat]=i;//number of exclusion for this atom
		i+=vdW_excl_numb[iat];
			
	}
	SetArraysize(&temp_delta,max_exclusions,"temp_delta", "FNC_POT::GenerateExclusion");//this will hold the delta indices for an atom's excluded atoms
	//Now filling the vdW_exclusion array with the indices of the excluded atoms
	for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
	{
		if (Topology::nexclusions_per_type[imoltype]>0 || Topology::nexclusionlines_per_type[imoltype]>0)
		{
			for (iat=0;iat<Topology::natoms_per_type[imoltype];iat++)
			{
				//Topology::first_index[Topology::cumul_atom[imoltype]+iat]is the RMC_index of the first occurence of this atomtype, in this molecule
				//	if this type would be the first in the system
				//Topology::RMCindex_offset[imoltype] is the offset of the RMC index to take into account that might not be the first moleculetype
				//first see the exculsion for the first instance of the molecule, to determine the pairs involved in 1-4 interaction
				if (Topology::first_index[Topology::cumul_atoms[imoltype] + iat]!=0)
					ind1=Topology::first_index[Topology::cumul_atoms[imoltype]+iat]-1+ Topology::RMCindex_offset[imoltype];
				else
				{
					//for virtualsite has to be done differently
					for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
					{
						if (iat + 1 == Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].atom_index)
						{
							ind1 = SimpleCfg::cumul[SimpleCfg::ntypes + Topology::cumul_virtuals[imoltype] + vi];
							break;//vi cycle
						}
					}
				}
				//graph search only needed for real atoms
				if (Topology::first_index[Topology::cumul_atoms[imoltype] + iat] != 0)
					count = GetExclusions(ind1, -1, Topology::nexclusions_per_type[imoltype] - 1, &vdW_exclusions[vdW_conv[ind1]], &vdW_exclusions[vdW_conv[ind1]]);
				else
					count = 0;
				for (iex1 = 0; iex1 < count; iex1++)
					temp_delta[iex1] = -1;//this will not be used
				//now add the excluded atoms from the exclusion lines
				for (il = 0; il < Topology::nexclusionlines_per_type[imoltype]; il++)//go through the exclusion lines of this moltype
				{
					iloffset= Topology::cumul_exclusionlines[imoltype] + il;//the index of the exclusion line
					ind_excl_centre = Topology::exclusion_centre[iloffset] ;//GROMACS index
					if (ind_excl_centre == -1)
						continue;//already removed
					if (ind_excl_centre == iat+1)
					{
						iex2 = 0;
						for (iex1 = 0; iex1 < Topology::nexcluded_neighs[iloffset]; iex1++)
						{
							//for virtual sites first_index cannot be used
							if (Topology::first_index[Topology::cumul_atoms[imoltype] + Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] - 1] == 0)
							{
								for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
								{
									if (Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] == Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].atom_index)
									{
										//RMC index +1
										ind_excluded = SimpleCfg::cumul[SimpleCfg::ntypes + Topology::cumul_virtuals[imoltype] + vi]+1;
										break;//vi cycle
									}
								}
							}
							else
								ind_excluded = Topology::first_index[Topology::cumul_atoms[imoltype] + Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] - 1] + Topology::RMCindex_offset[imoltype];
							vdW_exclusions[vdW_conv[ind1]+count+iex2] = ind_excluded;
							temp_delta[count + iex2] = Topology::delta_index[Topology::cumul_atoms[imoltype] + Topology::exclusions[Topology::cumul_nexcluded_neighs[iloffset] + iex1] - 1];
							iex2++;
						}
					}
				}
				//Determining, whether this atom and any of its exclusion is involved in a pair, and set the excluded pair's index to -x, if it is
				if (Topology::first_index[Topology::cumul_atoms[imoltype] + iat] != 0)//for real atoms
				{
					pair14 = 0;//no 1-4 interactions pairs as default
					for (ipair = 0; ipair < Topology::npairs_per_type[imoltype]; ipair++)
					{
						pair_offset = -1;
						if (Topology::pair_index[(Topology::cumul_pairs[imoltype] + ipair) * 2] - 1 == iat)//this  atom is involved in an exclusion, find its pair
							pair_offset = 1;
						else
						{
							if (Topology::pair_index[(Topology::cumul_pairs[imoltype] + ipair) * 2 + 1] - 1 == iat)
								pair_offset = 0;
						}
						if (pair_offset > -1)
						{
							for (iexcl = 0; iexcl < count; iexcl++)
							{
								//ind1 and ind2 start with 0, vdE_exclusions start with 1
								ind2 = Topology::first_index[Topology::cumul_atoms[imoltype] + Topology::pair_index[(Topology::cumul_pairs[imoltype] + ipair) * 2 + pair_offset] - 1] - 1 + Topology::RMCindex_offset[imoltype];;
								if (ind2 == vdW_exclusions[vdW_conv[ind1] + iexcl] - 1)
								{
									pair14 = 1;
									vdW_exclusions[vdW_conv[ind1] + iexcl] *= -1;//this excluded pair is on the 1-4 list
									//determine the type of the RMC partial to set the used_14partial indicator array
									//find the type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
									type1 = 0;
									while (ind1 >= SimpleCfg::cumul[(type1)+1])
										(type1)++;
									//type1 now points to the RMC type of ind1 
									//find the type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
									type2 = 0;
									while (ind2 >= SimpleCfg::cumul[(type2)+1])
										(type2)++;
									//type2 now points to the RMC type of ind1 
									partialind = (type1 <= type2 ? (type1 * SimpleCfg::ntypes - (type1 * (type1 + 1) / 2) + type2) : (type2 * SimpleCfg::ntypes - (type2 * (type2 + 1) / 2) + type1));
									if (RunParams::used_14partials[partialind] == 0)
									{
										RunParams::used_14partials[partialind] = 1;
										RunParams::tot_used14part++;
									}
									break;
								}
							}
						}
					}
				}
				for (imol=1;imol<Topology::nmol_per_type[imoltype];imol++)//for the other molecules of the same type
				{
					//imol*delta_index[cumul_atom[imoltype]+iat] is the offset of this occurence's index relative to the first
					ind2=imol*Topology::delta_index[Topology::cumul_atoms[imoltype]+iat];
					if (Topology::first_index[Topology::cumul_atoms[imoltype] + iat] != 0)//for the real atoms
					{
						GetExclusions(ind1 + ind2, -1, Topology::nexclusions_per_type[imoltype] - 1, &vdW_exclusions[vdW_conv[ind1 + ind2]], &vdW_exclusions[vdW_conv[ind1 + ind2]]);
						
						if (pair14)
						{
							for (iexcl = 0; iexcl < count; iexcl++)
							{
								if (vdW_exclusions[vdW_conv[ind1] + iexcl] < 0)
									vdW_exclusions[vdW_conv[ind1 + ind2] + iexcl] *= -1;
							}
						}
					}
					//fill in the virtual sites, if there is any, or fill all the exclusions for a virtual site
					for (iexcl = 0; iexcl < vdW_excl_numb[ind1] - count; iexcl++)
						vdW_exclusions[vdW_conv[ind1 + ind2] + count + iexcl] = imol*temp_delta[count+iexcl]+vdW_exclusions[vdW_conv[ind1] + count + iexcl];
												   							
				}//end of imol
			}//end of iat
		}//end of if there are exclusions
	}//end of cycle imoltype
}
  
//determing the number of exclusions for an atom
/*int FNC_POT::GetNExclusion(int iatom, int iprev,int iexcl, int *array_start)
{
	int inumb=0,ineigh;

	for (ineigh=0;ineigh<numbneigh[iatom];ineigh++)
	{
		if ((neighbours[converter[iatom]+ineigh]-1)==iprev)
			continue;
		inumb++;
		if (iexcl>0)
			inumb+=GetNExclusion(neighbours[converter[iatom]+ineigh]-1,iatom,iexcl-1);
	}
	return inumb;
};*/

//getting the number and indices of the excluded atom pairs for an atom
int FNC_POT::GetExclusions(int iatom, int iprev,int iexcl, int *array_start, int *excl_array, int checksize )
{
	int inumb=0,ineigh,i,found=0;
	if (iexcl < 0)
		return 0;
	for (ineigh=0;ineigh<numbneigh[iatom];ineigh++)
	{
		//to prevent to turn back
		if ((neighbours[converter[iatom]+ineigh]-1)==iprev)
			continue;
		
		//to check, whether the neighbour was not added already, as it can happen for 6-member rings on case of nexcl=3
		for (i = 0; i < excl_array - array_start; i++)
		{
			if (array_start[i] == neighbours[converter[iatom] + ineigh])
			{
				found = 1;
				break;
			}
		}
		if (found == 0)
		{
			if (checksize > 0 && excl_array - array_start+inumb==INPUT_SIZE)//only needs checking, when the number of exclusions are determined
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The number of exclusions found in the graph search of the molecular topology exceeds the array size defined by INPUT_SIZE=" << INPUT_SIZE << endl;
				cout << "Increase the INPUT_SIZE in units.h, compile it and try again!" << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}
			excl_array[inumb] = neighbours[converter[iatom] + ineigh];//index of the excluded atom starting with 1)
			inumb++;

			if (iexcl > 0 && iatom <natoms)
				inumb += GetExclusions(neighbours[converter[iatom] + ineigh] - 1, iatom, iexcl - 1, array_start, excl_array + inumb);
		}
		found = 0;
	}
	return inumb;
};

//save the exclusions for each atom
void FNC_POT::SaveExclusions(const char *file_name) const
{
	int i,j,k=0,*p_excl;//pointers to vdW_exclusion array
	ofstream file;

	OpenFile(file,file_name,"FNC_POT::SaveExclusions",0);//open file, check, whether it was successfully opened

	file<<"These are the exluded atoms for each atom and virtual sites from the non-bonded interaction calculation created by FNC_POT::SaveExclusions"<<endl;

	file<<J7<<natoms+nvirtuals<<J7<<nexclusions<<\
	"\t <-number of atoms & exclusions "<<endl;
	file<<endl;
	p_excl=vdW_exclusions;

	for(i=0;i<natoms+nvirtuals;i++)
	{
		file<<J10<<i+1<<J10<<*(vdW_excl_numb+i)<<\
		"\t <-atom index, number of exclusions"<<endl;
		
		for(j=0;j<vdW_excl_numb[i];j++)//saving the neighbours indices
		{
			file<<J10<<(*p_excl);
			p_excl++;
			k++;
		}
		file<<"\t <-excusion indices"<<endl;
	}
	file.close();
};


//to avoid an if inside the function there will be separate functions for the normal and 1-4 interaction

double FNC_POT::CalcLJ(int ipartial, double dsquare)
{
	double r6,rN;
	r6=pow(dsquare,3.0);
	rN=pow(dsquare,(double)(*LJ_pow/2));
	return(LJ_CN[ipartial]/rN-LJ_C6[ipartial]/r6);
};

double FNC_POT::CalcLJ14(int ipartial, double dsquare)
{
	double r6,rN;
	r6=pow(dsquare,3.0);
	rN=pow(dsquare,(double)(*LJ_pow/2));
	return(LJ14_CN[ipartial]/rN-LJ14_C6[ipartial]/r6);
};

double FNC_POT::CalcCoulomb(int ind1, int ind2, double dsquare, double scale)
{
	return(scale*f_Coulomb_A*Topology::charge[charge_type[ind1]]*Topology::charge[charge_type[ind2]]/sqrt(dsquare));
};

//calculate the NB interaction for old and new coordinates of the exclusions to subtract them, as during the histogram update NB
//is calculated for all the pairs
void FNC_POT::CalcNBExclusion(Move &move)
{
	int i,j,imoved,jmoved,icoord,nnb;
	int array_ind,index,partialind,pot_partialind,GRpartialind=0;
	int type1, type2,GRtype1,GRtype2;
	int skip_cycle;
	int *p_neighbour_ind;
	double *pmoved_chargegr,*pneigh_chargegr;//pointer to the charge group centres
	double *p_oldchargegr;//pointer to the old charge group centre coordinates of the moved atom
	double dcharge_c_sq,dcharge_c_sq_old;
	double dr[3],dsq,dsq_old,temp;
	double *p_newpos,*p_oldpos,*p_pos_neigh,*p_oldpos_neigh;
	SimpleCfg &conf=*config;

	p_newpos=move.newpos;//points to the (new) coordinates of the central atom
	p_oldpos=move.oldpos;//points to the (old) coordinates of the central atom
	p_oldchargegr=move.old_charge_c_pos;//sets the pointer to the old coordinates of the first moved atom's charge group centre coordinates

	for (imoved=0;imoved<move.nmoved+move.tot_moved_virtuals;imoved++)//for all atoms to be moved in the Move (no swaps are alowed in case of molecules)
		//but there can be virtual sites affected by the move
	{
		index=move.indices[imoved];//index of the moved atom 
		array_ind=vdW_conv[index];//index in the exclusion array
		p_neighbour_ind=vdW_exclusions+array_ind;//index of the excluded pair (starting with 1)
		pmoved_chargegr=conf.charge_gr_centre+3*charge_centre[index];//position of the charge group centre for the moved atom
		//find the type (0 to ntypes-1) of the moved atom
		type1=0;
		while (index>=conf.cumul[type1+1])
			(type1)++;
		//type1 now points to the type of the moved atom
		//find the GROMACS type (0 to nsep_GRtypes-1) of the atom 
		GRtype1=0;
		while (index>=SimpleCfg::cumul_GR[(GRtype1) +1])
			(GRtype1)++;
		//GRtype1 now points to the type index of the seperate GROMACS type segment, where the atom is located (the atoms index is continuously increasing)
		//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
		GRtype1=SimpleCfg::GRsep_GR_type[GRtype1];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters

		for (nnb=0;nnb<vdW_excl_numb[index];nnb++)//for each excluded pair
		{
			//first to check, whether this is not a moved-moved pair already calculated, 
			//among them the central and neighbour can be the same
			skip_cycle=0;//flag to indicate, whether the rest of the nnb cycle has to be skipped
			for (jmoved=0;jmoved<imoved;jmoved++)
			{
				if (abs(*p_neighbour_ind)-1==move.indices[jmoved])//this pair was aleady calculated 
				{
					//this is a moved-moved pair, the rest of the j cycle is skipped
					skip_cycle=1;
					p_neighbour_ind++;//skip the index of the neighbour in its own type
					break;
				}
			}
			if (skip_cycle)
				continue;//continue with the next moved atom
			//find the type (0 to ntypes-1) of the excluded pair atom
			type2=0;
			while (abs(*p_neighbour_ind)-1>=conf.cumul[(type2) +1])
				type2++;
			//type2 now points to the type of the excluded pair atom
			//Determining the index of the partial 
			partialind = (type1 <= type2 ? (type1 * SimpleCfg::ntypes - (type1 * (type1 + 1) / 2) + type2) : \
				(type2 * SimpleCfg::ntypes - (type2 * (type2 + 1) / 2) + type1));
			
			if (type1 < SimpleCfg::ntypes && type2 < SimpleCfg::ntypes)//normal atoms
			{
				pot_partialind = partialind;
			}
			else
			{
				//for the virtual sites the host atom's type will decide the partialind

				if (type1 >= SimpleCfg::ntypes)
					i = SimpleCfg::virtual_host_type[type1 - SimpleCfg::ntypes];
				else
					i = type1;
				if (type2 >= SimpleCfg::ntypes)
					j = SimpleCfg::virtual_host_type[type2 - SimpleCfg::ntypes];
				else
					j = type2;

				pot_partialind = (i <= j ? (i * SimpleCfg::ntypes - (i * (i + 1) / 2) + j) : (j * SimpleCfg::ntypes - (j * (j + 1) / 2) + i));
			}
			
			
			//calculating the new distance
			p_pos_neigh=conf.positions+(abs(*p_neighbour_ind)-1)*3;//totconf.positions contains the new coordinates of the moved atoms
			
			for (icoord=0;icoord<3;icoord++)
				dr[icoord]=*(p_pos_neigh+icoord) -*(p_newpos+icoord);
			
			GetMinImage(dr);
			
			dsq=0;
			for (icoord=0;icoord<3;icoord++)
				dsq+=dr[icoord]*dr[icoord];//distance squared
			
			//calculate the old distance, check, if the neighbour is a moved atom too, following this moved atom in the indices array
			p_oldpos_neigh=p_pos_neigh;//default, the neighbour was not moved
			pneigh_chargegr=conf.charge_gr_centre+3*charge_centre[abs(*p_neighbour_ind)-1];//for the charge group centre position, default, not moved

			for (jmoved=imoved+1;jmoved<move.tot_moved_atoms+move.tot_moved_virtuals;jmoved++)
			{
				if (abs(*p_neighbour_ind)-1==move.indices[jmoved])//this is a moved atom as well, old coordinates have to come from the oldpos array!
				{
					p_oldpos_neigh=move.oldpos+3*jmoved;
					pneigh_chargegr=move.old_charge_c_pos+3*jmoved;//for the charge group centre position

					break;
				}

			}
			
			for (icoord=0;icoord<3;icoord++)
				dr[icoord]=*(p_oldpos_neigh+icoord) -*(p_oldpos+icoord);
			
			GetMinImage(dr);
			
			dsq_old=0;
			for (icoord=0;icoord<3;icoord++)
				dsq_old+=dr[icoord]*dr[icoord];//distance squared

			//calculating the old distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
			for (icoord=0;icoord<3;icoord++)
				dr[icoord]=*(p_oldchargegr+icoord)- *pneigh_chargegr++;
			GetMinImage(dr);
			dcharge_c_sq_old=0;
			for (icoord=0;icoord<3;icoord++)
				dcharge_c_sq_old+=pow(dr[icoord],2.0);

			//calculating the distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
			pneigh_chargegr=conf.charge_gr_centre+3*charge_centre[abs(*p_neighbour_ind)-1];
			for (icoord=0;icoord<3;icoord++)
				dr[icoord]=*(pmoved_chargegr+icoord)- *pneigh_chargegr++;
			GetMinImage(dr);
			dcharge_c_sq=0;
			for (icoord=0;icoord<3;icoord++)
				dcharge_c_sq+=pow(dr[icoord],2.0);	

			if (dcharge_c_sq<=rundat->vdW_cutoff_sq || dcharge_c_sq_old<=rundat->vdW_cutoff_sq)//check the cutoff
			{
				//find the GROMACS type (0 to nsep_GRtypes-1) of the atom
				GRtype2=0;
				while (abs(*p_neighbour_ind)-1>=SimpleCfg::cumul_GR[(GRtype2) +1])
					(GRtype2)++;
				//GRtype2 now points to the type index of the seperate GROMACS type segment where the atom is located (the atoms index is continuously increasing)

				//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
				GRtype2=SimpleCfg::GRsep_GR_type[GRtype2];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters

				//GRpartialind is GROMACS type based partial index, needed to have the vdW parameters, but the histogram is split into RMC-type partials
				GRpartialind=(GRtype1<=GRtype2 ? (GRtype1*Topology::nGRtypes-(GRtype1*(GRtype1+1)/2)+GRtype2) : \
					(GRtype2*Topology::nGRtypes-(GRtype2*(GRtype2+1)/2)+GRtype1));

			
			}
			//logfile.precision(15);//GO
			//logfile.setf(ios::scientific, ios::floatfield);
			

			//adding the excluded pair's vdW potential, as it was calculated during the histogram change update
			if (dcharge_c_sq_old<=rundat->vdW_cutoff_sq)
			{
				temp = (this->*CalcvdW)(GRpartialind, dsq_old);
				if (fabs(temp) > POT_SPLIT_HIGH)
					vdW_pot_l[pot_partialind] += temp;
				else if (fabs(temp) < POT_SPLIT_LOW)
					vdW_pot_s[pot_partialind] += temp;
				else
					vdW_pot[pot_partialind]+=temp ;
				if (fabs(temp) > POT_WARNING)
				{
					vdW_warning = pot_partialind + 1;//gives warning in the hst file and on screen that unreasonbly large potential was present
					//logfile << "excl+ " << " " << pot_partialind << " " << index << " " << abs(*p_neighbour_ind) - 1 << " " << dsq_old << " " << temp << " " << vdW_pot[pot_partialind] << " l " << vdW_pot_l[pot_partialind] << " s " << vdW_pot_s[pot_partialind] << endl;
				}
				//now remove the contribution of the new distance for 1-4
				if (*p_neighbour_ind<0)//this is an 1-4 pair
					vdW14_pot[pot_partialind]-=(this->*CalcvdW14)(GRpartialind,dsq_old);
			}

			if (dcharge_c_sq <= rundat->vdW_cutoff_sq)
			{

				//subtracting the excluded pair's vdW potential, as it was calculated during the histogram change update
				temp = (this->*CalcvdW)(GRpartialind, dsq);
				if (fabs(temp) > POT_SPLIT_HIGH)
					vdW_pot_l[pot_partialind] -= temp;
				else if (fabs(temp) < POT_SPLIT_LOW)
					vdW_pot_s[pot_partialind] -= temp;
				else
					vdW_pot[pot_partialind] -= temp;
				if (fabs(temp) > POT_WARNING)
				{
					vdW_warning = pot_partialind + 1;//gives warning in the hst file and on screen that unreasonbly large potential was present
					//logfile << "excl- " << " " << pot_partialind << " " << index << " " << abs(*p_neighbour_ind) - 1 << " " << dsq << " " << temp << " " << vdW_pot[pot_partialind] << " l " << vdW_pot_l[pot_partialind] << " s " << vdW_pot_s[pot_partialind] << endl;
				}
				//now add the contribution to 1-4 interaction, if this is an 1-4 pair
				if (*p_neighbour_ind < 0)//this is an 1-4 pair
					vdW14_pot[pot_partialind] += (this->*CalcvdW14)(GRpartialind, dsq);
			}
			

			//adding the Coulomb contribution for the old distance, as it was removed during Histogram change calculation
			if (dcharge_c_sq_old<=rundat->Coulomb_cutoff_sq)//check the cutoff
			{
				temp= CalcCoulomb(index, abs(*p_neighbour_ind) - 1, dsq_old);
				if (fabs(temp) > POT_WARNING)
					Coul_warning = pot_partialind+1;//gives warning in the hst file and on screen that unreasonbly large potential was present
				Coulomb_pot[pot_partialind] += temp;

				//logfile << "excl+ " << " " << pot_partialind << " " << index << " " << *p_neighbour_ind << " " << dsq_old << " " << temp << " " << Coulomb_pot[partialind] << endl;

				if (*p_neighbour_ind<0)//this is an 1-4 pair
					Coulomb14_pot[pot_partialind]-=RunParams::Coulomb14_fudge*temp;
			}
			if (dcharge_c_sq <= rundat->Coulomb_cutoff_sq)//check the cutoff
			{
				//subtracting the Coulomb contribution for the new distance
				temp = CalcCoulomb(index, abs(*p_neighbour_ind) - 1, dsq);
				if (fabs(temp) > POT_WARNING)
					Coul_warning = pot_partialind + 1;//gives warning in the hst file and on screen that unreasonbly large potential was present
				Coulomb_pot[pot_partialind] -= temp;


				//logfile << "excl- " << " " << pot_partialind << " " << index <<" "<<*p_neighbour_ind << " "<< dsq << " " << temp<< " "<<Coulomb_pot[partialind]<<endl;

				//now add the contribution to 1-4 interaction, if this is an 1-4 pair
				if (*p_neighbour_ind < 0)//this is an 1-4 pair
					Coulomb14_pot[pot_partialind] += RunParams::Coulomb14_fudge * temp;
			}
			//logfile.precision(6);//GO
			//logfile.unsetf(ios::scientific);
			p_neighbour_ind++;//go to next pair
		}
		p_newpos+=3;//go to the position of the next moved atom
		p_oldpos+=3;//go to the position of the next moved atom
		p_oldchargegr+=3;//sets the pointer to the old charge group centre coordinates of the next moved atom
	}
};

//Now calculate the interaction potentials for the BONDS
void FNC_POT::CalcBonds()
{
	int i,icoord,ctype,nnb,array_ind;
	double *posi,*posj,vec[3],d,dsq;
	SimpleCfg &conf=*config;
	
	for (i=0;i<natoms;i++)
	{
		array_ind=converter[i];
		for (nnb=0;nnb<numbneigh[i];nnb++)
		{

			posi=conf.positions+3*i;
			posj=conf.positions+3*(neighbours[array_ind]-1);
			
			for (icoord=0;icoord<3;icoord++)
				vec[icoord]=posi[icoord]-posj[icoord];
			
			GetMinImage(vec);
			dsq=0;
			for (icoord=0;icoord<3;icoord++)
				dsq+=vec[icoord]*vec[icoord];
			
			d=sqrt(dsq);

			ctype=consttypes[array_ind];
			//calculating the harmonic bond potential 1/2*k*(r-b)^2
			bond_pot[ctype]+=0.5*bond_k[ctype]*pow(d-bond_r[ctype],2);
			array_ind++;
		}
		
		
	}
	
	//A bond was calculated twice, for i central j neighbour and the other way around, and divide by 2 as part of the formula
	bond_tot_pot=0;
	for (i=0;i<Topology::nbond_types;i++)
	{
		bond_pot[i]=bond_pot[i]/2.0;
		bond_tot_pot+=bond_pot[i];
	}


}

//calculate the change in the BOND interaction 
void FNC_POT::CalcBondsChange(Move &move)
{
	int imoved,icoord,nnb,ind,ctype,i;
	int array_ind,index;
	int do_calc;//flag to indicate, whether to calculate for this pair
	int *p_neighbour_ind,*p_ind;
	double dr[3],dsq,d;
	double *p_newpos,*p_oldpos,*p_pos_neigh;
	SimpleCfg &conf=*config;

	p_newpos=move.newpos;//points to the (new) coordinates of the central atom
	p_oldpos=move.oldpos;//points to the (old) coordinates of the central atom

	//keep in mind, that the neighbours array are holding the atom indices starting with 1, so 
	//1 has to be subtracted
	for (imoved=0;imoved<move.nmoved;imoved++)//for all atoms to be moved in the Move (no swaps are alowed in case of molecules, virtual sites cannot be involved in the move)
	{
		index=move.indices[imoved];//index of the moved atom 
		array_ind=converter[index];//index in the arrays
		p_neighbour_ind=neighbours+array_ind;//index of the bonded pair
		
		for (nnb=0;nnb<numbneigh[index];nnb++)//for each bond
		{
			//see if the neighbour is moved as well, and preceding this moved atom in the moved
			//indices array, as in this case potential change was already calculated
			p_ind=move.indices;
			do_calc=1;
			for (ind=0;ind<imoved;ind++)
			{
				if (*p_ind==*p_neighbour_ind-1)
				{
					do_calc=0;//this pair was already calculated
					break;
				}
				p_ind++;
			}
			if (do_calc)//this angle was not yet calculated
			{
				//calculating the new distance
				p_pos_neigh=conf.positions+(*p_neighbour_ind-1)*3;//totconf.positions contains the new coordinates of the moved atoms
				ctype=consttypes[array_ind];
				
				for (icoord=0;icoord<3;icoord++)
					dr[icoord]=*(p_pos_neigh+icoord) -*(p_newpos+icoord);
				
				GetMinImage(dr);
				dsq=0;
				for (icoord=0;icoord<3;icoord++)
					dsq+=dr[icoord]*dr[icoord];//distance squared
				d=sqrt(dsq);
				//calculating the harmonic bond potential 1/2*k*(r-b)^2
				bond_pot[ctype]+=0.5*bond_k[ctype]*pow(d-bond_r[ctype],2);

				//calculating the old distance
				//now we must check, whether the neighbour is a moved atom following the present one in the list
				//as in this case the old coordinates have to come from the old position array
				p_ind=move.indices+imoved+1;
				for (ind=imoved+1;ind<move.nmoved;ind++)
				{
					if (*p_ind==*p_neighbour_ind-1)//the neighbour is a moved atom too
					{
						p_pos_neigh=move.oldpos+3* ind;
						break;
					}
					p_ind++;
				}	
				for (icoord=0;icoord<3;icoord++)
					dr[icoord]=*(p_pos_neigh+icoord) -*(p_oldpos+icoord);
									
				GetMinImage(dr);
				dsq=0;
				for (icoord=0;icoord<3;icoord++)
					dsq+=dr[icoord]*dr[icoord];//distance squared
				d=sqrt(dsq);
				//calculating the harmonic bond potential 1/2*k*(r-b)^2
				bond_pot[ctype]-=0.5*bond_k[ctype]*pow(d-bond_r[ctype],2);
			}
			array_ind++;//go to the next neighbour's value in the consttypes array
			p_neighbour_ind++;//go to the next neighbour
		}
		p_oldpos+=3;
		p_newpos+=3;
	}
	bond_tot_pot=0;
	for (i=0;i<Topology::nbond_types;i++)
		bond_tot_pot+=bond_pot[i];
	
};

//Now calculate the interaction potentials for the ANGLES
void FNC_POT::CalcAngles()
{
	int i,icoord,ctype,array_ind,nnb;
	double *pos_end,*posj,vecji[3],vecjk[3],dji,djk,costheta;
	SimpleCfg &conf=*config;

	for (i=0;i<natoms;i++)
	{

		array_ind=converter[i];
		for (nnb=0;nnb<numbneigh[i];nnb++)
		{
			//all angle will be calculated only ones, when the central atom is the middle one
			if (neighbours[array_ind]>0)//the central atom is the middle one
			{

				posj=conf.positions+3*i;//central (j)
				
				pos_end=conf.positions+3*(neighbours[array_ind]-1);//neighbour (i)
				for (icoord=0;icoord<3;icoord++)
					vecji[icoord]=*pos_end++ - *posj++;
				
				posj-=3;//put the coordinate to x
				pos_end=conf.positions+3*(neighbours2[array_ind]-1); //neighbour (k)
				//calculating the ij vector components
				for (icoord=0;icoord<3;icoord++)
					vecjk[icoord]=*pos_end++ - *posj++;
				
				//taking into account the periodical boundary conditions 
				GetMinImage(vecji);
				GetMinImage(vecjk);
						
				costheta=0;
				for (icoord=0;icoord<3;icoord++)
					costheta+=*(vecji+icoord) * *(vecjk+icoord);
				dji=0;
				djk=0;
				for (icoord=0;icoord<3;icoord++)
				{
					dji+=vecji[icoord]*vecji[icoord];
					djk+=vecjk[icoord]*vecjk[icoord];
				}

				costheta/=sqrt(dji * djk);//dividing by the length of the vectors
				//to avoid rounding errors, put back into the -1 - +1 range
				costheta=(costheta<-1 ? -1 : (costheta>1 ? 1 :costheta));
			
				ctype=consttypes[array_ind];
				//calculating the harmonic angle potential 1/2*k*(theta-theta0)^2

				angle_pot[ctype]+=Topology::angle_k[ctype]*pow(acos(costheta)-angle_ang_rad[ctype],2);
			}
			array_ind++;
		}
	}
	angle_tot_pot=0;
	for (i=0;i<Topology::nangle_types;i++)
	{
		angle_pot[i]=angle_pot[i]/2;//now dividing by 2 according to the formula
		angle_tot_pot+=angle_pot[i];
	}
}

//calculate the change in the interaction 
void FNC_POT::CalcAnglesChange(Move &move)
{
	int imoved,icoord,nnb,nnb2,ind,ctype,i;
	int array_ind,array_ind2,index,index2;
	int mid_ind, end1_ind, end2_ind, old1,old2,neigh_ind;
	int do_calc;//flag to indicate, whether to calculate for this pair
	int *p_neighbour_ind1,*p_neighbour_ind2,*p_nneighbour_ind1,*p_nneighbour_ind2,*p_ind;
	double vecji[3],vecjk[3],dji,djk,costheta;
	double *p_oldpos,*posj,*pos_end;
	SimpleCfg &conf=*config;

	p_oldpos=move.oldpos;//points to the (old) coordinates of the central atom

	//keep in mind, that the neighbours and neighbours2 arrays are holding the atom indices starting with 1, so 
	//1 has to be subtracted
	for (imoved=0;imoved<move.nmoved;imoved++)//for all atoms to be moved in the Move (no swaps are alowed in case of molecules, virtual sites are not involved in bonding)
	{
		index=move.indices[imoved];//index of the moved atom 
		array_ind=converter[index];//index in the exclusion array
		p_neighbour_ind1=neighbours+array_ind;//index of the first neighbour in the angle
		p_neighbour_ind2=neighbours2+array_ind;//index of the second neighbour in the angle			
		
		for (nnb=0;nnb<numbneigh[index];nnb++)//for each bond
		{	
			//see, if one of the neighbours is moved as well, and potential change was already calculated
			p_ind=move.indices;
			do_calc=1;
			for (ind=0;ind<imoved;ind++)
			{
				neigh_ind=0;//the neighbour does not match any previous moved atom
				if (*p_ind==abs(*p_neighbour_ind1)-1)//this previous moved atom is among the angle components of this one
					neigh_ind=*p_neighbour_ind2;//store the value of the third neighbour
				else
				{
					if (*p_ind==*p_neighbour_ind2-1)//this previous moved atom is among the angle components of this one
						neigh_ind=abs(*p_neighbour_ind1);//store the value of the second neighbour
				}
				if (neigh_ind>0)//one neighbour is the same, as one of the previous moved atoms, see, if the third atom is the same, and this triple was already calculated
				{
					index2=move.indices[ind];//index of the moved atom 
					array_ind2=converter[index2];//index in the exclusion array
					p_nneighbour_ind1=neighbours+array_ind2;//index of the first neighbour in the angle
					p_nneighbour_ind2=neighbours2+array_ind2;//index of the second neighbour in the angle			
					for (nnb2=0;nnb2<numbneigh[index2];nnb2++)//going through the angles of this previously moved atom
					{
						if ((index==abs(*p_nneighbour_ind1)-1 && neigh_ind==*p_nneighbour_ind2) || (index==*p_nneighbour_ind2-1 && neigh_ind==abs(*p_nneighbour_ind1)))
						{
							do_calc=0;//this angle was already calculated
							break;
						}
						p_nneighbour_ind1++;//go to the next neighbour
						p_nneighbour_ind2++;//go to the next neighbour
					}
					if (do_calc==0)
						break;//do not search more, this angle was already calculated
				}
				p_ind++;
			}
			if (do_calc)
			{
				//the angle is defined by i--j--k atoms

				//first decide, whether any of the two neighbours are not moved atoms in the moved indices array,
				//as in that case the coordinates have to come from the oldpos array for those as well
				p_ind=move.indices;
				old1=-1;//not a moved atom by default
				old2=-1;
				for (ind=0;ind<move.nmoved;ind++)
				{
					if (ind==imoved)
					{
						p_ind++;
						continue;//this is the same as the actual moved atom, go to next
					}
						
					if (*p_ind==abs(*p_neighbour_ind1)-1)//this neighbour is an other moves atom
						old1=ind;
												
					if (*p_ind==*p_neighbour_ind2-1)//this neighbour is an other moves atom
						old2=ind;
					if (old1>-1 && old2>-1)
						break;//do not search more, both neighbours are moved as well
					p_ind++;
				}

				if (*p_neighbour_ind1>0)//deciding, whether the moved atom is a middle or an end one
				{
					//middle one
					mid_ind=index;
					end1_ind=abs(*p_neighbour_ind1)-1;
					posj=p_oldpos;
					if (old1>-1)
						pos_end=move.oldpos+3*old1;
					else
						pos_end=conf.positions+3*(end1_ind);//neighbour (i)
				}
				else
				{
					//end one
					end1_ind=index;
					mid_ind=abs(*p_neighbour_ind1)-1;
					pos_end=p_oldpos;
					if (old1>-1)
						posj=move.oldpos+3*old1;
					else
						posj=conf.positions+3*mid_ind;//central (j)

				}

				end2_ind=*p_neighbour_ind2-1;
				ctype=consttypes[array_ind];
				
				//now calculated for the old angle first
				for (icoord=0;icoord<3;icoord++)
					vecji[icoord]=*pos_end++ - *posj++;
				
				if (old2>-1)
					pos_end=move.oldpos+3*old2;
				else
					pos_end=conf.positions+3*(end2_ind); //neighbour (k)
				
				posj-=3;//put the coordinate to x
				//calculating the ij vector components
				for (icoord=0;icoord<3;icoord++)
					vecjk[icoord]=*pos_end++ - *posj++;
				
				//taking into account the periodical boundary conditions 
				GetMinImage(vecji);
				GetMinImage(vecjk);
						
				costheta=0;
				for (icoord=0;icoord<3;icoord++)
					costheta+=*(vecji+icoord) * *(vecjk+icoord);
				dji=0;
				djk=0;
				for (icoord=0;icoord<3;icoord++)
				{
					dji+=vecji[icoord]*vecji[icoord];
					djk+=vecjk[icoord]*vecjk[icoord];
				}

				costheta/=sqrt(dji * djk);//dividing by the length of the vectors
				//to avoid rounding errors, put back into the -1 - +1 range
				costheta=(costheta<-1 ? -1 : (costheta>1 ? 1 :costheta));
			
				//calculating the harmonic angle potential 1/2*k*(theta-theta0)^2
				angle_pot[ctype]-=0.5*Topology::angle_k[ctype]*pow(acos(costheta)-angle_ang_rad[ctype],2);

				//now calculate for the new angle
			
				posj=conf.positions+3*mid_ind;//central (j)
			
				pos_end=conf.positions+3*(end1_ind);//neighbour (i)
				for (icoord=0;icoord<3;icoord++)
					vecji[icoord]=*pos_end++ - *posj++;
				
				posj-=3;//put the coordinate to x
				pos_end=conf.positions+3*(end2_ind); //neighbour (k)
				//calculating the ij vector components
				for (icoord=0;icoord<3;icoord++)
					vecjk[icoord]=*pos_end++ - *posj++;
				
				//taking into account the periodical boundary conditions 
				GetMinImage(vecji);
				GetMinImage(vecjk);
						
				costheta=0;
				for (icoord=0;icoord<3;icoord++)
					costheta+=*(vecji+icoord) * *(vecjk+icoord);
				dji=0;
				djk=0;
				for (icoord=0;icoord<3;icoord++)
				{
					dji+=vecji[icoord]*vecji[icoord];
					djk+=vecjk[icoord]*vecjk[icoord];
				}

				costheta/=sqrt(dji * djk);//dividing by the length of the vectors
				//to avoid rounding errors, put back into the -1 - +1 range
				costheta=(costheta<-1 ? -1 : (costheta>1 ? 1 :costheta));
			
				//calculating the harmonic angle potential 1/2*k*(theta-theta0)^2
				angle_pot[ctype]+=0.5*Topology::angle_k[ctype]*pow(acos(costheta)-angle_ang_rad[ctype],2);
			}
			array_ind++;//go to the next neighbour's value in the consttypes array
			p_neighbour_ind1++;//go to the next neighbour
			p_neighbour_ind2++;//go to the next neighbour
		}
		p_oldpos+=3;
	}
	angle_tot_pot=0;
	for (i=0;i<Topology::nangle_types;i++)
		angle_tot_pot+=angle_pot[i];
};

//Now calculate the interaction potentials for the DIHEDRALS
void FNC_POT::CalcDihedrals()
{
	int i,icoord,ctype,nnb,array_ind;
	double *pos_end,*pos_mid,*pos_mid2,phi;
	double vecij[3],vecjk[3],veckl[3];

	SimpleCfg &conf=*config;
	for (i=0;i<natoms;i++)
	{
		array_ind=converter[i];
		for (nnb=0;nnb<numbneigh[i];nnb++)
		{
			//all dihedral will be calculated twice, when the neighbour array sign is negative,
			if (neighbours[array_ind]<0)//the central atom is an end one
			{
				pos_end=conf.positions+3*i;//atom (i)
				pos_mid=conf.positions+3*(abs(neighbours[array_ind])-1);//atom (j)
				//calculating the ij vector components
				for (icoord=0;icoord<3;icoord++)
					vecij[icoord]=*pos_mid++ - *pos_end++;
				
				pos_mid-=3;//put the coordinate to x
				pos_mid2=conf.positions+3*(neighbours2[array_ind]-1);//atom (k)
								
				//calculating the jk vector components
				for (icoord=0;icoord<3;icoord++)
					vecjk[icoord]=*pos_mid2++ - *pos_mid++;
				
				pos_mid2-=3;//put the coordinate to x
				pos_end=conf.positions+3*(neighbours3[array_ind]-1); //atom (l)

				//calculating the kl vector components
				for (icoord=0;icoord<3;icoord++)
					veckl[icoord]=*pos_end++ - *pos_mid2++;
				//taking into account the periodical boundary conditions 
				GetMinImage(vecij);
				GetMinImage(vecjk);
				GetMinImage(veckl);

				phi=CalcDihAngle(vecij,vecjk,veckl);
				ctype=consttypes[array_ind];
			
				//calling the appropriate dihedral function depending on the type of the dihedral
				(this->*DihedralFunction[dih_type[array_ind]-1])(ctype,1,phi);
			}
			array_ind++;
		}
	}
	//divide it by two as both diherdral was calculated twice
	perdihedral_tot_pot=0;
	for (i=0;i<Topology::nperdihedral_types;i++)
	{
		perdihedral_pot[i]=perdihedral_pot[i]/2;
		perdihedral_tot_pot+=perdihedral_pot[i];
	}
	harmdihedral_tot_pot=0;
	for (i=0;i<Topology::nharmdihedral_types;i++)
	{
		harmdihedral_pot[i]=harmdihedral_pot[i]/2;
		harmdihedral_tot_pot+=harmdihedral_pot[i];
	}
	RBdihedral_tot_pot=0;
	for (i=0;i<Topology::nRBdihedral_types;i++)
	{
		RBdihedral_pot[i]=RBdihedral_pot[i]/2;
		RBdihedral_tot_pot+=RBdihedral_pot[i];
	}
}

	//calculate the change in the interaction 
void FNC_POT::CalcDihedralsChange(Move &move)
{
	int imoved,icoord,nnb,nnb2,ind,ctype,i;
	int array_ind,array_ind2,index,index2;
	int neigh_ind1,neigh_ind2;
	int mid_ind,mid2_ind,end1_ind,end2_ind,old1,old2,old3;
	int do_calc;//flag to indicate, whether to calculate for this pair
	int *p_neighbour_ind1,*p_neighbour_ind2,*p_neighbour_ind3;
	int *p_nneighbour_ind1,*p_nneighbour_ind2,*p_nneighbour_ind3,*p_ind;
	double vecij[3],vecjk[3],veckl[3],phi;
	double *p_oldpos,*pos_end,*pos_mid,*pos_mid2;
	SimpleCfg &conf=*config;

	p_oldpos=move.oldpos;//points to the (old) coordinates of the central atom

	//keep in mind, that the neighbours, neighbours2 and neighbours3 arrays are holding the atom indices starting with 1, so 
	//1 has to be subtracted
	for (imoved=0;imoved<move.nmoved;imoved++)//for all atoms to be moved in the Move (no swaps are alowed in case of molecules, virtual sites are not involved in bonding)
	{
		index=move.indices[imoved];//index of the moved atom 
		array_ind=converter[index];//index in the exclusion array
		p_neighbour_ind1=neighbours+array_ind;//index of the first neighbour in the dihedral
		p_neighbour_ind2=neighbours2+array_ind;//index of the second neighbour in the dihedral			
		p_neighbour_ind3=neighbours3+array_ind;//index of the second neighbour in the dihedral	
		
		for (nnb=0;nnb<numbneigh[index];nnb++)//for each bond
		{				
			//see, if one of the neighbours is moved as well, and potential change was already calculated
			p_ind=move.indices;
			do_calc=1;
			for (ind=0;ind<imoved;ind++)
			{
				neigh_ind1=0;//the neighbour does not match any previous moved atom
				neigh_ind2=0;//the neighbour does not match any previous moved atom
				if (*p_ind==abs(*p_neighbour_ind1)-1)
				{ //this previous moved atom is among the dihedral components of this one
						neigh_ind1=*p_neighbour_ind2;//store the value of the third neighbour
						neigh_ind2=*p_neighbour_ind3;//store the value of the 4th neighbour
				}
				else
				{
					if (*p_ind==*p_neighbour_ind2-1)
					{//this previous moved atom is among the dihedral components of this one
						neigh_ind1=abs(*p_neighbour_ind1);//store the value of the second neighbour
						neigh_ind2=*p_neighbour_ind3;//store the value of the 4th neighbour
					}
					else
					{
						if(*p_ind==*p_neighbour_ind3-1)
						{//this previous moved atom is among the dihedral components of this one
							neigh_ind1=abs(*p_neighbour_ind1);//store the value of the second neighbour
							neigh_ind2=*p_neighbour_ind2;//store the value of the 4th neighbour
						}
					}
				}
				if (neigh_ind1>0)//one neighbour is the same, as one of the previous moved atoms, see, if the third atom is the same, and this triple was already calculated
				{
					index2=move.indices[ind];//index of the moved atom 
					array_ind2=converter[index2];//index in the exclusion array
					p_nneighbour_ind1=neighbours+array_ind2;//index of the first neighbour in the dihedral
					p_nneighbour_ind2=neighbours2+array_ind2;//index of the second neighbour in the dihedral			
					p_nneighbour_ind3=neighbours3+array_ind2;//index of the third neighbour in the dihedral			
					for (nnb2=0;nnb2<numbneigh[index2];nnb2++)//going through the dihedrals of this previously moved atom
					{
						if ((index==abs(*p_nneighbour_ind1)-1 && neigh_ind1==*p_nneighbour_ind2 && neigh_ind2==*p_nneighbour_ind3) || \
							(index==abs(*p_nneighbour_ind1)-1 && neigh_ind2==*p_nneighbour_ind2 && neigh_ind1==*p_nneighbour_ind3) || \
							(index==*p_nneighbour_ind2-1 && neigh_ind1==abs(*p_nneighbour_ind1) && neigh_ind2==*p_nneighbour_ind3) || \
							(index==*p_nneighbour_ind2-1 && neigh_ind2==abs(*p_nneighbour_ind1) && neigh_ind1==*p_nneighbour_ind3) || \
							(index==*p_nneighbour_ind3-1 && neigh_ind1==abs(*p_nneighbour_ind1) && neigh_ind2==*p_nneighbour_ind2) || \
							(index==*p_nneighbour_ind3-1 && neigh_ind2==abs(*p_nneighbour_ind1) && neigh_ind1==*p_nneighbour_ind2))
						{
							do_calc=0;//this dihedral was already calculated
							break;
						}
						p_nneighbour_ind1++;//go to the next neighbour
						p_nneighbour_ind2++;//go to the next neighbour
						p_nneighbour_ind3++;//go to the next neighbour
					}
					if (do_calc==0)
						break;//do not search more, this dihedral was already calculated
				}//end if neigh_ind>0
				p_ind++;
			}//end of ind cycle
			if (do_calc)
			{
				//the dihedral is defined by i--j--k--l atoms

				//first decide, whether any of the two neighbours are not moved atoms in the moved indices array,
				//as in that case the coordinates have to come from the oldpos array for those as well
				p_ind=move.indices;
				old1=-1;//not a moved atom by default
				old2=-1;
				old3=-1;
				for (ind=0;ind<move.nmoved;ind++)
				{
					if (ind==imoved)
					{
						p_ind++;
						continue;//this is the same as the actual moved atom, go to next
					}
					if (*p_ind==abs(*p_neighbour_ind1)-1)//this neighbour is an other moved atom
						old1=ind;
					if (*p_ind==*p_neighbour_ind2-1)//this neighbour is an other moved atom
						old2=ind;
					if (*p_ind==*p_neighbour_ind3-1)//this neighbour is an other moved atom
						old3=ind;
					if (old1>-1 && old2>-1  && old3>-1)
						break;//do not search more, all the three neighbours are moved as well
					p_ind++;
				}

				if (*p_neighbour_ind1>0)//deciding, whether the moved atom is a middle or an end one
				{
					//middle one
					mid_ind=index;
					end1_ind=abs(*p_neighbour_ind1)-1;
					pos_mid=p_oldpos;
					if (old1>-1)
						pos_end=move.oldpos+3*old1;
					else
						pos_end=conf.positions+3*end1_ind;//neighbour (i)
				}
				else
				{
					//end one
					end1_ind=index;
					mid_ind=abs(*p_neighbour_ind1)-1;
					pos_end=p_oldpos;
					if (old1>-1)
						pos_mid=move.oldpos+3*old1;
					else
						pos_mid=conf.positions+3*mid_ind;//central (j)

				}

				mid2_ind=*p_neighbour_ind2-1;
				end2_ind=*p_neighbour_ind3-1;
				ctype=consttypes[array_ind];
				
				//now calculating for the old dihedral first
				//calculating the ij vector components
				for (icoord=0;icoord<3;icoord++)
					vecij[icoord]=*pos_mid++ - *pos_end++;
				
				if (old2>-1)
					pos_mid2=move.oldpos+3*old2;
				else
					pos_mid2=conf.positions+3*mid2_ind;//atom (k)

				if (old3>-1)
					pos_end=move.oldpos+3*old3;
				else
					pos_end=conf.positions+3*end2_ind; ////atom (l)
			
				pos_mid-=3;//put the coordinate to x
										
				//calculating the jk vector components
				for (icoord=0;icoord<3;icoord++)
					vecjk[icoord]=*pos_mid2++ - *pos_mid++;
				
				pos_mid2-=3;//put the coordinate to x

				//calculating the kl vector components
				for (icoord=0;icoord<3;icoord++)
					veckl[icoord]=*pos_end++ - *pos_mid2++;
				//taking into account the periodical boundary conditions 
				GetMinImage(vecij);
				GetMinImage(vecjk);
				GetMinImage(veckl);
				phi=CalcDihAngle(vecij,vecjk,veckl);
											
				//calling the appropriate dihedral function depending on the type of the dihedral
				(this->*DihedralFunction[dih_type[array_ind]-1])(ctype,-1,phi);

				//calculating for the new distance
				pos_end=conf.positions+3*end1_ind;//atom (i)
				pos_mid=conf.positions+3*mid_ind;//atom (j)
				//calculating the ij vector components
				for (icoord=0;icoord<3;icoord++)
					vecij[icoord]=*pos_mid++ - *pos_end++;
				
				pos_mid-=3;//put the coordinate to x
				pos_mid2=conf.positions+3*mid2_ind;//atom (k)
								
				//calculating the jk vector components
				for (icoord=0;icoord<3;icoord++)
					vecjk[icoord]=*pos_mid2++ - *pos_mid++;
				
				pos_mid2-=3;//put the coordinate to x
				pos_end=conf.positions+3*(end2_ind); //atom (l)

				//calculating the kl vector components
				for (icoord=0;icoord<3;icoord++)
					veckl[icoord]=*pos_end++ - *pos_mid2++;
				//taking into account the periodical boundary conditions 
				GetMinImage(vecij);
				GetMinImage(vecjk);
				GetMinImage(veckl);
				phi=CalcDihAngle(vecij,vecjk,veckl);
				
				//calling the appropriate dihedral function depending on the type of the dihedral
				(this->*DihedralFunction[dih_type[array_ind]-1])(ctype,1,phi);
		
			}
			array_ind++;//go to the next neighbour's value in the consttypes array
			p_neighbour_ind1++;//go to the next neighbour
			p_neighbour_ind2++;//go to the next neighbour
			p_neighbour_ind3++;//go to the next neighbour
		}
		p_oldpos+=3;
	}
	perdihedral_tot_pot=0;
	for (i=0;i<Topology::nperdihedral_types;i++)
		perdihedral_tot_pot+=perdihedral_pot[i];
	
	harmdihedral_tot_pot=0;
	for (i=0;i<Topology::nharmdihedral_types;i++)
		harmdihedral_tot_pot+=harmdihedral_pot[i];
	
	RBdihedral_tot_pot=0;
	for (i=0;i<Topology::nRBdihedral_types;i++)
		RBdihedral_tot_pot+=RBdihedral_pot[i];
	
};

//calculating the dihedral angle
double FNC_POT::CalcDihAngle(double *vecij, double *vecjk, double *veckl)
{
	int icoord;
	double djk,a,b,prodijk[3],prodjkl[3];
	
	//The dihedral angle will be calculated by:
	//phi=atan2(|vecjk|vecij*[vecjk � veckl],[vecij � vecjk] * [vecjk � veckl]),
	//where  i----->j----->k----->l
	//        vecij  vecjk  veckl
	//* means the scalar or dot or inner product, and � is the cross or vector product
	
	djk=0;
	for (icoord=0;icoord<3;icoord++)
		djk+=vecjk[icoord]*vecjk[icoord];
	sqrt(djk);
	
	VectorProduct(vecij,vecjk,prodijk);
	VectorProduct(vecjk,veckl,prodjkl);
	
	djk=0;
	for (icoord=0;icoord<3;icoord++)
		djk+=vecjk[icoord]*vecjk[icoord];
	djk=sqrt(djk);

	a=0;
	for (icoord=0;icoord<3;icoord++)
		a+=*(vecij+icoord) * *(prodjkl+icoord);

	b=0;
	for (icoord=0;icoord<3;icoord++)
		b+=*(prodijk+icoord) * *(prodjkl+icoord);
					
	return(atan2(djk*a,b));//returning the dihedral angle in radian
}


void FNC_POT::CalcPeriodicDihedrals(int consttype, int sign,double theta)
{
	perdihedral_pot[consttype]+=sign*Topology::perdihedral_k[consttype]*(1+cos(Topology::perdihedral_multiplicity[consttype]*theta-perdih_ang_rad[consttype]));
}

void FNC_POT::CalcHarmonicDihedrals(int consttype, int sign,double theta)
{
	harmdihedral_pot[consttype]+=sign*0.5* Topology::harmdihedral_k[consttype]*pow(theta-harmdih_ang_rad[consttype],2);
}

void FNC_POT::CalcRBDihedrals(int consttype, int sign,double theta)
{
	int j;
	//for RB the angle is defined, that the 0� is the trans
	for (j=0;j<6;j++)
		RBdihedral_pot[consttype]+=sign* *(RBdih_C[j]+consttype)*pow(cos(theta-PI),j);
}

//reading the tabulated potential
// it is assumed, that it is given in r (A) - U(J)
void FNC_POT::ReadTabPot(ifstream &file)
{
	int i,j,ntotpoints=0;
	double r_first, r_prev, dr,*r,*u;
	char *name, *conv_numb = NULL;
	ifstream potfile;
	tabpot_filename = new char[nused_potpartials][FILE_NAME_SIZE];
	SetArraysize(&name,NAME_SIZE,"name", "FNC_POT::ReadTabPot");
	
	SetArraysize(&npot_points, nused_potpartials, "npot_points", "FNC_POT::ReadTabPot");
	SetArraysize(&tabpot_offset, nused_potpartials, "tabpot_offset", "FNC_POT::ReadTabPot");
	SetArraysize(&tabpot_index, nused_potpartials, "tabpot_index", "FNC_POT::ReadTabPot");
	SetArraysize(&tabpot_dr, nused_potpartials, "tabpot_dr", "FNC_POT::ReadTabPot");
	SetArraysize(&tab_cutoff, nused_potpartials, "tab_cutoff", "FNC_POT::ReadTabPot");

	for (i = 0; i < nused_potpartials; i++)
	{
		file >> tabpot_filename[i];
		SkipLine(file);
		SafeOpenTextFile(potfile, tabpot_filename[i]);
		if (CheckFileState(potfile, "FNC_POT::ReadTabPot", tabpot_filename[i]) == 0)
		{
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		};
		cout << "Reading data from " << tabpot_filename[i] << endl;
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "npot_points[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		npot_points[i] = ReadThisLine(potfile, 1, 1, name, "FNC_POT::ReadTabPot");
		
		if (i == 0)
		{
			ntotpoints += npot_points[i];
			SetArraysize(&r_pot, ntotpoints, "r_pot", "FNC_POT::ReadTabPot");
			SetArraysize(&u_pot, ntotpoints, "u_pot", "FNC_POT::ReadTabPot");
		}
		else
		{
			ResizeArray(&ntotpoints, npot_points[i] + ntotpoints, &r_pot, "r_pot", "FNC_POT::ReadTabPot");
			ntotpoints -= npot_points[i];
			ResizeArray(&ntotpoints, npot_points[i] + ntotpoints, &u_pot, "u_pot", "FNC_POT::ReadTabPot");
		}
		tabpot_offset[i] = ntotpoints- npot_points[i];//offset for each partial, can be used for both r_pot and u_pot
		SkipLine(potfile);
		r = r_pot + tabpot_offset[i];
		u = u_pot + tabpot_offset[i];
		for (j = 0; j < npot_points[i]; j++)
		{
			*r = ReadThisLine(potfile, 3, 1.0, "r_pot", "FNC_POT::ReadTabPot");
			*u++ = ReadThisLine(potfile, 1, 1.0, "u_pot", "FNC_POT::ReadTabPot");
			*r /= RunParams::boxedge;
			r++;
		}
		next_line_pos = -1;
		potfile.close();
		//determine, whether there is potential given for  the cutoff range
		tab_cutoff[i]=((RunParams::vdW_cutoff / RunParams::boxedge > r_pot[ntotpoints - 1])? r_pot[ntotpoints - 1]: RunParams::vdW_cutoff / RunParams::boxedge);

		tabpot_index[i] = ReadThisLine(file, 1, 1, "tabpot_index", "FNC_POT::ReadTabPot");
		tabpot_index[i]--;//to start the indexing with 0
		for (j = 0; j < i; j++)//check, whether this partial was not already given
		{
			if (tabpot_index[i] == tabpot_index[j])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The partial index "<<tabpot_index[i]+1<<" for tabulated potential file " << tabpot_filename[i] << " is the same as for " << endl;
				cout << "tabulated potential file " << tabpot_filename[j]<<"!" << endl;
				cout << "This makes no sense! Exiting..." << endl;
				CleanExit();

			}
		}
		RunParams::vdW_weight[i] = ReadThisLine(file, 1, 1.0, "vdW_weight", "FNC_POT::ReadTabPot");
		if (RunParams::vdW_weight[i] < 0)
			ChiSquared::calc_sigma = 1;
		//check, whether the potential is equidistant
		r_prev = r_pot[tabpot_offset[i]+1];
		r_first = r_pot[tabpot_offset[i]];
		dr = r_prev - r_first;
		r = r_pot + tabpot_offset[i];
		for (j = 2; j < npot_points[i]; j++)
		{
			if ((fabs(r[j] - r_prev) - dr) > LOAD_TOL)
			{
				cout.setf(ios::floatfield, ios::scientific);
				cout.precision(8);
				cout << "\n*****ERROR*****" << endl;
				cout << "The spacing of tabulated potential file " << tabpot_filename[i] << " is not equidistant!" << endl;
				cout << "The disfference between the first and second point is " << dr << ", while between point " << j << " and " << j + 1 << " is " << fabs(r[j] - r_prev) << endl;
				cout << "Only equidistant tabulated potential can be used! Exiting..." << endl;
				CleanExit();
			}
			r_prev = r_pot[tabpot_offset[i]+j];
		}
		
		tabpot_dr[i] = r_pot[tabpot_offset[i]+1] - r_pot[tabpot_offset[i]];
	}
	
	delete [] name;
	if (conv_numb != NULL)
		delete [] conv_numb;
	
};

//reading the tabulated potential files, and finishing the initialization
// it is assumed, that it is given in r (A) - U(J)
void FNC_POT::ReadTabPotFree()
{
	int i, j, ntotpoints = 0;
	double r_first, r_prev, dr, *r, *u;
	char *name, *conv_numb = NULL;
	ifstream potfile;
	SetArraysize(&name, NAME_SIZE, "name", "FNC_POT::ReadTabPot");

	SetArraysize(&npot_points, nused_potpartials, "npot_points", "FNC_POT::ReadTabPot");
	SetArraysize(&tabpot_offset, nused_potpartials, "tabpot_offset", "FNC_POT::ReadTabPot");
	SetArraysize(&tabpot_dr, nused_potpartials, "tabpot_dr", "FNC_POT::ReadTabPot");
	SetArraysize(&tab_cutoff, nused_potpartials, "tab_cutoff", "FNC_POT::ReadTabPot");

	for (i = 0; i < nused_potpartials; i++)
	{
		SafeOpenTextFile(potfile, tabpot_filename[i]);
		if (CheckFileState(potfile, "FNC_POT::ReadTabPot", tabpot_filename[i]) == 0)
		{
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		};
		cout << "\nLoading data from " << tabpot_filename[i] << endl;
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "npot_points[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		npot_points[i] = ReadThisLine(potfile, 1, 1, name, "FNC_POT::ReadTabPot");

		if (i == 0)
		{
			ntotpoints += npot_points[i];
			SetArraysize(&r_pot, ntotpoints, "r_pot", "FNC_POT::ReadTabPot");
			SetArraysize(&u_pot, ntotpoints, "u_pot", "FNC_POT::ReadTabPot");
		}
		else
		{
			ResizeArray(&ntotpoints, npot_points[i] + ntotpoints, &r_pot, "r_pot", "FNC_POT::ReadTabPot");
			ntotpoints -= npot_points[i];
			ResizeArray(&ntotpoints, npot_points[i] + ntotpoints, &u_pot, "u_pot", "FNC_POT::ReadTabPot");
		}
		tabpot_offset[i] = ntotpoints - npot_points[i];//offset for each partial, can be used for both r_pot and u_pot
		SkipLine(potfile);
		r = r_pot + tabpot_offset[i];
		u = u_pot + tabpot_offset[i];
		for (j = 0; j < npot_points[i]; j++)
		{
			*r = ReadThisLine(potfile, 3, 1.0, "r_pot", "FNC_POT::ReadTabPot");
			*u++ = ReadThisLine(potfile, 1, 1.0, "u_pot", "FNC_POT::ReadTabPot");
			*r /= RunParams::boxedge;
			r++;
		}
		next_line_pos = -1;
		potfile.close();
		//determine, whether there is potential given for  the cutoff range
		tab_cutoff[i] = ((RunParams::vdW_cutoff / RunParams::boxedge > r_pot[ntotpoints - 1]) ? r_pot[ntotpoints - 1] : RunParams::vdW_cutoff / RunParams::boxedge);

		tabpot_index[i]--;//to start the indexing with 0
		for (j = 0; j < i; j++)//check, whether this partial was not already given
		{
			if (tabpot_index[i] == tabpot_index[j])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The partial index " << tabpot_index[i] + 1 << " for tabulated potential file " << tabpot_filename[i] << " is the same as for " << endl;
				cout << "tabulated potential file " << tabpot_filename[j] << "!" << endl;
				cout << "This makes no sense! Exiting..." << endl;
				CleanExit();

			}
		}
		
		//check, whether the potential is equidistant
		r_prev = r_pot[tabpot_offset[i] + 1];
		r_first = r_pot[tabpot_offset[i]];
		dr = r_prev - r_first;
		r = r_pot + tabpot_offset[i];
		for (j = 2; j < npot_points[i]; j++)
		{
			if ((fabs(r[j] - r_prev) - dr) > LOAD_TOL)
			{
				cout.setf(ios::floatfield, ios::scientific);
				cout.precision(8);
				cout << "\n*****ERROR*****" << endl;
				cout << "The spacing of tabulated potential file " << tabpot_filename[i] << " is not equidistant!" << endl;
				cout << "The disfference between the first and second point is " << dr << ", while between point " << j << " and " << j + 1 << " is " << fabs(r[j] - r_prev) << endl;
				cout << "Only equidistant tabulated potential can be used! Exiting..." << endl;
				CleanExit();
			}
			r_prev = r_pot[tabpot_offset[i] + j];
		}

		tabpot_dr[i] = r_pot[tabpot_offset[i] + 1] - r_pot[tabpot_offset[i]];
	}

	delete [] name;
	if (conv_numb != NULL)
		delete [] conv_numb;

};

//calculate the tabulated potential, it is assumed (it was checked), that the spacing is equdisant
//ipartial will be index of this partial among the used tabulated potential files!
//first find the potential interval, than determine the U(r) with linear interpolaion
double  FNC_POT::CalcTabPot(int ipartial, double dsquare)
{
	int index;
	double *r_first,*u_first,dist=sqrt(dsquare);
	r_first = r_pot +FNC_POT::tabpot_offset[ipartial];
	u_first = u_pot + FNC_POT::tabpot_offset[ipartial];
	index = (int)((dist -*r_first)/ tabpot_dr[ipartial]);

	return((dist- r_first[index])/tabpot_dr[ipartial]*(u_first[index+1]-u_first[index])+u_first[index]);

};

void FNC_POT::VectorProduct(double *a, double *b, double *result)
{
	result[0]=a[1]*b[2]-a[2]*b[1];
	result[1]=a[2]*b[0]-a[0]*b[2];
	result[2]=a[0]*b[1]-a[1]*b[0];
}


//mode 0: Copying the values from the normal array to the old, 1: from the old to the normal
void FNC_POT::UpdatePot(int mode)
{
	int i;
	double *pb_dest,*pb_source,*pa_dest,*pa_source,*ppd_dest=NULL,*ppd_source=NULL,*phd_dest=NULL,*phd_source=NULL,*pRBd_dest=NULL,*pRBd_source=NULL;
	double *pvdW_dest,*pvdW_source,*pvdW14_dest,*pvdW14_source,*pCoulomb_dest,*pCoulomb_source,*pCoulomb14_dest,*pCoulomb14_source;

	if (mode==0)//from new to old
	{
		pvdW_dest = vdW_pot_old;
		pvdW_source = vdW_pot;
		pvdW14_dest = vdW14_pot_old;
		pvdW14_source = vdW14_pot;
		pCoulomb_dest = Coulomb_pot_old;
		pCoulomb_source = Coulomb_pot;
		pCoulomb14_dest = Coulomb14_pot_old;
		pCoulomb14_source = Coulomb14_pot;
		pb_dest = bond_pot_old;
		pb_source = bond_pot;
		pa_dest = angle_pot_old;
		pa_source = angle_pot;
		ppd_dest = perdihedral_pot_old;
		ppd_source = perdihedral_pot;
		phd_dest = harmdihedral_pot_old;
		phd_source = harmdihedral_pot;
		pRBd_dest = RBdihedral_pot_old;
		pRBd_source = RBdihedral_pot;
		
	}
	else
	{//from old to new
		pvdW_dest = vdW_pot;
		pvdW_source = vdW_pot_old;
		pvdW14_dest = vdW14_pot;
		pvdW14_source = vdW14_pot_old;
		pCoulomb_dest = Coulomb_pot;
		pCoulomb_source = Coulomb_pot_old;
		pCoulomb14_dest = Coulomb14_pot;
		pCoulomb14_source = Coulomb14_pot_old;
		pb_dest = bond_pot;
		pb_source = bond_pot_old;
		pa_dest = angle_pot;
		pa_source = angle_pot_old;
		ppd_dest = perdihedral_pot;
		ppd_source = perdihedral_pot_old;
		phd_dest = harmdihedral_pot;
		phd_source = harmdihedral_pot_old;
		pRBd_dest = RBdihedral_pot;
		pRBd_source = RBdihedral_pot_old;
	}
	if (pot_type==BOND)
	{
		for (i=0;i<pot_dim;i++)
		{
			*pvdW_dest++=*pvdW_source++;
			if (RunParams::potential == 1)
			{
				*pvdW14_dest++ = *pvdW14_source++;
				*pCoulomb_dest++ = *pCoulomb_source++;
				*pCoulomb14_dest++ = *pCoulomb14_source++;
			}
		}

		for (i=0;i<Topology::nbond_types;i++)
			*pb_dest++=*pb_source++;
	}
	if (pot_type==ANGLE)
		for (i=0;i<Topology::nangle_types;i++)
			*pa_dest++=*pa_source++;

	if (pot_type==DIHEDRAL)
	{
		for (i=0;i<Topology::nperdihedral_types;i++)
			*ppd_dest++=*ppd_source++;
		for (i=0;i<Topology::nharmdihedral_types;i++)
			*phd_dest++=*phd_source++;
		for (i=0;i<Topology::nRBdihedral_types;i++)
			*pRBd_dest++=*pRBd_source++;
	}
	
};

//saving the potentials (bond, angle, dihedral, vdW, Coulomb, 1-4) in binary format
void FNC_POT::SavePotBinary() const
{
	int i,sum;
	ofstream file;
	
	OpenFile(file,potfilename,"FNC_POT::SavePotBinary",1);//open file, check, whether it was successfully opened

	file.write(reinterpret_cast<const char *>(&Topology::nbond_types), sizeof(Topology::nbond_types));
	if (Topology::nbond_types>0)
		file.write(reinterpret_cast<const char *>(bond_pot),  Topology::nbond_types*sizeof(*bond_pot));

	file.write(reinterpret_cast<const char *>(&Topology::nangle_types), sizeof(Topology::nangle_types));
	if (Topology::nangle_types>0)
		file.write(reinterpret_cast<const char *>(angle_pot),  Topology::nangle_types*sizeof(*angle_pot));

	file.write(reinterpret_cast<const char *>(&Topology::nperdihedral_types), sizeof(Topology::nperdihedral_types));
	if (Topology::nperdihedral_types>0)
		file.write(reinterpret_cast<const char *>(perdihedral_pot),  Topology::nperdihedral_types*sizeof(*perdihedral_pot));

	file.write(reinterpret_cast<const char *>(&Topology::nharmdihedral_types), sizeof(Topology::nharmdihedral_types));
	if (Topology::nharmdihedral_types>0)
		file.write(reinterpret_cast<const char *>(harmdihedral_pot),  Topology::nharmdihedral_types*sizeof(*harmdihedral_pot));

	file.write(reinterpret_cast<const char *>(&Topology::nRBdihedral_types), sizeof(Topology::nRBdihedral_types));
	if (Topology::nRBdihedral_types>0)
		file.write(reinterpret_cast<const char *>(RBdihedral_pot),  Topology::nRBdihedral_types*sizeof(*RBdihedral_pot));

	file.write(reinterpret_cast<const char *>(&RunParams::potential), sizeof(RunParams::potential));
	if (RunParams::potential > 0)
	{
		file.write(reinterpret_cast<const char *>(&pot_dim), sizeof(pot_dim));
		file.write(reinterpret_cast<const char *>(vdW_pot), pot_dim * sizeof(*vdW_pot));
		
	}
	switch (RunParams::potential)
	{
	case 1:
	{
		file.write(reinterpret_cast<const char *>(vdW14_pot), pot_dim * sizeof(*vdW14_pot));
		file.write(reinterpret_cast<const char *>(Coulomb_pot), pot_dim * sizeof(*Coulomb_pot));
		file.write(reinterpret_cast<const char *>(Coulomb14_pot), pot_dim * sizeof(*Coulomb14_pot));

		file.write(reinterpret_cast<const char *>(LJ_pow), sizeof(*LJ_pow));
		file.write(reinterpret_cast<const char *>(&nGRpartials), sizeof(nGRpartials));
		//writing the reduced LJ C(6) and C(12) 
		file.write(reinterpret_cast<const char *>(LJ_C6), nGRpartials * sizeof(*LJ_C6));
		file.write(reinterpret_cast<const char *>(LJ14_C6), nGRpartials * sizeof(*LJ_C6));
		file.write(reinterpret_cast<const char *>(LJ_CN), nGRpartials * sizeof(*LJ14_CN));
		file.write(reinterpret_cast<const char *>(LJ14_CN), nGRpartials * sizeof(*LJ14_CN));
		file.write(reinterpret_cast<const char *>(&RunParams::vdW_cutoff), sizeof(RunParams::vdW_cutoff));
		file.write(reinterpret_cast<const char *>(&RunParams::Coulomb_cutoff), sizeof(RunParams::Coulomb_cutoff));
		break;
	}
	case 10:
	{
	
		file.write(reinterpret_cast<const char *>(&nused_potpartials), sizeof(nused_potpartials));
		file.write(reinterpret_cast<const char *>(&RunParams::temperature), sizeof(RunParams::temperature));
		file.write(reinterpret_cast<const char *>(tabpot_index), nused_potpartials*sizeof(*tabpot_index));
		file.write(reinterpret_cast<const char *>(npot_points), nused_potpartials *sizeof(*npot_points));
		sum = 0;
		for (i = 0; i < nused_potpartials; i++)
			sum += npot_points[i];
		file.write(reinterpret_cast<const char *>(r_pot), sum*sizeof(*r_pot));
		file.write(reinterpret_cast<const char *>(u_pot), sum*sizeof(*u_pot));
		file.write(reinterpret_cast<const char *>(&RunParams::vdW_cutoff), sizeof(RunParams::vdW_cutoff));
		break;
	}
	//ANN potential is saved by Aenet::SavePotBinary, but the same *.pot file is used
	}

	file.close();
};

//loading the potentials (bond, angle, dihedral, vdW) from the binary format *.pot file
int FNC_POT::LoadPotBinary() const
{
	int i,j,temp, *temp_i;
	double *temp_d;
	string warning = "\tHistogram and potentials has to be recalculated!";
	string warning2 = "cutoff could not be read and checked, the *.pot file was craeted with an older version of RMC than 2.4!\nOnly continue, if you are sure, that the same cutoff was used in the original run!";

	ifstream file;

	CleanOpen(file,potfilename,1);//open file
	if (CheckFileState(file,"FNC_POT::LoadPotBinary",potfilename)==0)//check, whether it was successfully opened
	{
		cout << "\nWARNING(" << ++warn << "): There is no open "<<potfilename<<" file to load the potentials from"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);//loading was not successful
	}

	//reading the dimension of the bond_pot array
	if (ReadBin(file, &temp, 1, "Topology::nbond_types", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp!=Topology::nbond_types)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The dimension of the bond potential array is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<Topology::nbond_types<<" based on the topology file!"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);
	}
	if (Topology::nbond_types>0)
		if (ReadBin(file, bond_pot, Topology::nbond_types, "bond_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		
	bond_tot_pot=0;
	for (i=0;i<Topology::nbond_types;i++)
		bond_tot_pot+=bond_pot[i];
	
	//reading the dimension of the angle_pot array
	if (ReadBin(file, &temp, 1, "Topology::nangle_types", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp!=Topology::nangle_types)
	{
		cout << "\nWARNING(" << ++warn << "):I NCONSISTENCY! The dimension of the angle potential array is " << temp << endl;
		cout<<"\tfrom the " << potfilename << " file, and "<<Topology::nangle_types<<" based on the topology file!"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);
	}
	if (Topology::nangle_types>0)
		if (ReadBin(file, angle_pot, Topology::nangle_types, "angle_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }

	angle_tot_pot=0;
	for (i=0;i<Topology::nangle_types;i++)
		angle_tot_pot+=angle_pot[i];

	//reading the dimension of the perdihedral_pot array
	if (ReadBin(file, &temp, 1, "Topology::nperdihedral_types", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp!=Topology::nperdihedral_types)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The dimension of the periodic dihedral potential array is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<Topology::nperdihedral_types<<" based on the topology file!"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);
	}
	if (Topology::nperdihedral_types>0)
		if (ReadBin(file, perdihedral_pot, Topology::nperdihedral_types, "perdihedral_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }

	perdihedral_tot_pot=0;
	for (i=0;i<Topology::nperdihedral_types;i++)
		perdihedral_tot_pot+=perdihedral_pot[i];

	//reading the dimension of the harmdihedral_pot array
	if (ReadBin(file, &temp, 1, "Topology::nharmdihedral_types", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); };
	if (temp!=Topology::nharmdihedral_types)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The dimension of the harmonic dihedral potential array is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<Topology::nharmdihedral_types<<" based on the topology file!"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);
	}
	if (Topology::nharmdihedral_types>0)
		if (ReadBin(file, harmdihedral_pot, Topology::nharmdihedral_types, "harmdihedral_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }

	harmdihedral_tot_pot=0;
	for (i=0;i<Topology::nharmdihedral_types;i++)
		harmdihedral_tot_pot+=harmdihedral_pot[i];
	
	//reading the dimension of the RBdihedral_pot array
	if (ReadBin(file, &temp, 1, "Topology::nRBdihedral_types", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp!=Topology::nRBdihedral_types)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The dimension of the Ryckaert-Bellemans dihedral potential array is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<Topology::nRBdihedral_types<<" based on the topology file!"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);
	}
	if (Topology::nRBdihedral_types>0)
		if (ReadBin(file, RBdihedral_pot, Topology::nRBdihedral_types, "RBdihedral_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }

	RBdihedral_tot_pot=0;
	for (i=0;i<Topology::nRBdihedral_types;i++)
		RBdihedral_tot_pot+=RBdihedral_pot[i];
	

	//reading the RunParams::potential
	if (ReadBin(file, &temp, 1, "RunParams::potential", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp!=RunParams::potential)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The potential switch is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<RunParams::potential<<" from the "<<datfilename<<" file!"<<endl;
		cout<<"\tHistogram and potentials has to be recalculated!"<<endl;
		return(0);
	}
	//reading the dimension of the vdW array
	if (ReadBin(file, &temp, 1, "pot_dim", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (RunParams::potential>0)
	{
		if (temp != pot_dim)
		{
			cout << "\nWARNING(" << ++note << "): ";
			if (RunParams::potential == 1)
				cout << "INCONSISTENCY! The dimension of the non-bonded"; 
			else
				cout << "INCONSISTENCY! The dimension of the tabulated" << endl;;
			cout<<" interaction potential array is " << temp<< endl;
			cout<<"\tfrom the " << potfilename << " file, and " << pot_dim;
			if (RunParams::potential == 1)
				cout<< " based on the configuration file!" << endl;
			else
				cout << " based on the "<<datfilename<<" file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		if (pot_dim > 0)
		{
			if (ReadBin(file,vdW_pot, pot_dim, "vdW_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
			if (RunParams::potential == 1)
			{
				if (ReadBin(file, vdW14_pot, pot_dim, "vdW_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
				if (ReadBin(file, Coulomb_pot, pot_dim, "vdW_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
				if (ReadBin(file, Coulomb14_pot, pot_dim, "vdW_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
			}

		}
	}
	//Calculating the totals
	for (i=0;i<pot_dim;i++)
	{
		vdW_tot_pot+=vdW_pot[i];
		if (RunParams::potential == 1)
		{
			Coulomb_tot_pot += Coulomb_pot[i];
			vdW14_tot_pot += vdW14_pot[i];
			Coulomb14_tot_pot += Coulomb14_pot[i];
		}
	};
	switch (RunParams::potential)
	{
	case 1:
	{
		if (ReadBin(file, &temp, 1, "LJ_pow", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		if (temp != *LJ_pow)
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The power of the LJ repulsive term is " << temp << endl;
			cout <<  "\tfrom the " << potfilename << " file, and " << *LJ_pow << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		if (ReadBin(file, &temp, 1, "nGRpartials", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		if (temp != nGRpartials)
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of the GROMACS partials is " << temp <<  endl;
			cout << "\tfrom the " << potfilename << " file, and " << nGRpartials << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		SetArraysize(&temp_d, nGRpartials, "temp_d", "FNC_POT::LoadPotBinary");
		if (ReadBin(file, temp_d, nGRpartials, "LJ_C6", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		for (i = 0; i < nGRpartials; i++)
		{
			if (fabs(temp_d[i] - LJ_C6[i]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The " << i + 1 << ". reduced LJ(6) parameter is " << temp_d[i] << endl;
				cout << "\tfrom the " << potfilename << " file, and " << LJ_C6[i] << " based on the " << datfilename << " file!" << endl;
				cout << "\tHistogram and potentials has to be recalculated!" << endl;
				delete [] temp_d;
				return(0);
			}

		}
		if (ReadBin(file, temp_d, nGRpartials, "LJ14_C6", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		for (i = 0; i < nGRpartials; i++)
		{
			if (fabs(temp_d[i] - LJ14_C6[i]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The " << i + 1 << ". reduced LJ(6) parameter is " << temp_d[i] << endl;
				cout << "\tfrom the " << potfilename << " file, and " << LJ14_C6[i] << " based on the " << datfilename << " file!" << endl;
				cout << "\tHistogram and potentials has to be recalculated!" << endl;
				delete [] temp_d;
				return(0);
			}
		}
		if (ReadBin(file, temp_d, nGRpartials, "LJ_CN", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		for (i = 0; i < nGRpartials; i++)
		{
			if (fabs(temp_d[i] - LJ_CN[i]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The " << i + 1 << ". reduced LJ_C" << *LJ_pow << "[i] parameter is " << temp_d[i] <<  endl;
				cout << "\tfrom the " << potfilename << " file, and " << LJ_CN[i] << " based on the " << datfilename << " file!" << endl;
				cout << "\tHistogram and potentials has to be recalculated!" << endl;
				delete [] temp_d;
				return(0);
			}
		}
		if (ReadBin(file, temp_d, nGRpartials, "L14_CN", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		for (i = 0; i < nGRpartials; i++)
		{
			if (fabs(temp_d[i] - LJ14_CN[i]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The " << i + 1 << ". reduced LJ14_C" << *LJ_pow << "[i] parameter is " << temp_d[i] <<endl;
				cout <<  "\tfrom the " << potfilename << " file, and " << LJ14_CN[i] << " based on the " << datfilename << " file!" << endl;
				cout << "\tHistogram and potentials has to be recalculated!" << endl;
				delete [] temp_d;
				return(0);
			}
		}
		if (ReadBin(file, temp_d, 1, "vdW_cutoff", "FNC_POT::LoadPotBinary", potfilename,-1) == 0) 
		{
			cout << "\nWARNING(" << ++warn << "): vdW "<<warning2 << endl;
			file.close();
			return(0); 
		}
		if (temp_d[0] != RunParams::vdW_cutoff )
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The vdW cutoff is " << temp_d[0] <<  endl;
			cout << "\tfrom the " << potfilename << " file, and " << RunParams::vdW_cutoff  << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		if (ReadBin(file, temp_d, 1, "Coulomb_cutoff", "FNC_POT::LoadPotBinary", potfilename,-1) == 0) 
		{
			cout << "\nWARNING(" << ++warn << "): Coulomb "<<warning2 << endl;
			file.close();
			return(0); 
		}
		if (temp_d[0] != RunParams::Coulomb_cutoff )
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The Coulomb cutoff is " << temp_d[0] <<  endl;
			cout << "\tfrom the " << potfilename << " file, and " << RunParams::Coulomb_cutoff  << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		delete [] temp_d;
		break;
	}
	case 10:
	{
		if (ReadBin(file, &temp, 1, "nused_potpartials", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		if (temp != nused_potpartials)
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of tabulated potential file(s) is " << temp <<  endl;
			cout << "\tfrom the " << potfilename << " file, and " << nused_potpartials << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		SetArraysize(&temp_d, nused_potpartials, "temp_d", "FNC_POT::LoadPotBinary");
		if (ReadBin(file, temp_d, 1, "temperature", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		if (fabs(temp_d[0] - RunParams::temperature) > LOAD_TOL)
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The temperature is " << temp_d[0] << endl;
			cout << "\tfrom the " << potfilename << " file, and " << RunParams::temperature << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			delete [] temp_d;
			return(0);
		}
		SetArraysize(&temp_i, nused_potpartials, "temp_i", "FNC_POT::LoadPotBinary");
		if (ReadBin(file, temp_i, nused_potpartials, "tabpot_index", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		for (i = 0; i < nused_potpartials; i++)
		{
			if (temp_i[i] != tabpot_index[i])
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The index of the partial the tabulated potential is used for in case of the " << i + 1 << ".  tabulated potential  file is " << temp_i[i] <<endl;
				cout <<"\tfrom the " << potfilename << " file, and " << tabpot_index[i] << " based on the " << datfilename << " file!" << endl;
				cout << "\tHistogram and potentials has to be recalculated!" << endl;
				delete [] temp_d;
				delete [] temp_i;
				return(0);
			}
		}
		if (ReadBin(file, temp_i, nused_potpartials, "npot_points", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		temp = 0;
		for (i = 0; i < nused_potpartials; i++)
		{
			if (temp_i[i] != npot_points[i])
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of points in the " << i + 1 << ". tabulated potential file is " << temp_i[i] << endl;
				cout << "\tfrom the " << potfilename << " file, and " << npot_points[i] << " based on the " << datfilename << " file!" << endl;
				cout << "\tHistogram and potentials has to be recalculated!" << endl;
				delete [] temp_d;
				delete [] temp_i;
				return(0);
			}
			temp += temp_i[i];//total number of points
		}
		delete [] temp_d;
		SetArraysize(&temp_d, temp, "temp_d", "FNC_POT::LoadPotBinary");
		if (ReadBin(file, temp_d, temp, "r_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		temp = 0;
		for (i = 0; i < nused_potpartials; i++)
		{
			for (j = 0; j < npot_points[i]; j++)
			{
				if (fabs(temp_d[temp] - r_pot[temp])> LOAD_TOL)
				{
					cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The "<<j+1<<". r value in the " << i + 1 << ". tabulated potential file is " << temp_d[temp] << endl;
					cout << "\tfrom the " << potfilename << " file, and " << r_pot[temp] << " based on the " << tabpot_filename[i] << " file!" << endl;
					cout << "\tHistogram and potentials has to be recalculated!" << endl;
					delete [] temp_d;
					delete [] temp_i;
					return(0);
				}
				temp++;
			}
		}
		if (ReadBin(file, temp_d, temp, "u_pot", "FNC_POT::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
		temp = 0;
		for (i = 0; i < nused_potpartials; i++)
		{
			for (j = 0; j < npot_points[i]; j++)
			{
				if (fabs(temp_d[temp] - u_pot[temp]) > LOAD_TOL)
				{
					cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The " << j + 1 << ". U value in the " << i + 1 << ". tabulated potential file is " << temp_d[temp] << endl;
					cout << "\tfrom the " << potfilename << " file, and " << u_pot[temp] << " based on the " << tabpot_filename[i] << " file!" << endl;
					cout << "\tHistogram and potentials has to be recalculated!" << endl;
					delete [] temp_d;
					delete [] temp_i;
					return(0);
				}
				temp++;
			}
		}
		if (ReadBin(file, temp_d, 1, "tabulared potential cutoff", "FNC_POT::LoadPotBinary", potfilename,-1) == 0) 
		{
			cout << "\nWARNING(" << ++warn << "): Tabulated "<<warning2 << endl;
			file.close();
			return(0); 
		}
		if (temp_d[0] != RunParams::vdW_cutoff )
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The vdW cutoff is " << temp_d[0] <<  endl;
			cout << "\tfrom the " << potfilename << " file, and " << RunParams::vdW_cutoff  << " based on the " << datfilename << " file!" << endl;
			cout << "\tHistogram and potentials has to be recalculated!" << endl;
			return(0);
		}
		delete [] temp_i;
		delete [] temp_d;
	}
	}
	if (::debug)
		cout<<"Loading of the potential related parameters was successful"<<endl;
	return (1);//loading was successful
};

//setting the static potential arrays
void FNC_POT::SetPotParams()
{
	int padding_offset;
	double min_cutoff;

	if (RunParams::fnc==4)
		index_offset=-1;//to compensate, that for the BONDS, for which the indices are stored in the
	//neighbours array, and are checked by the same code parts as the normal fnc,if not fnc=4 is used, 
	//for normal fnc index_offset=0

	//setting the size of some arrays
	switch (RunParams::potential)
	{
	case 1:
	{
		pot_dim = SimpleCfg::ntypes*(SimpleCfg::ntypes + 1) / 2;
		//determinig, whether to check for exclusions during the potential correction due to the movement of the charge centres between the 
		//not moved atoms of the involved charge centres
		//The largest bond length and exclusion regardless the molecule type will be used and the smaller of the two cutoffs
		min_cutoff = (RunParams::vdW_cutoff > RunParams::Coulomb_cutoff ? RunParams::Coulomb_cutoff : RunParams::vdW_cutoff);

		if (Topology::bond_r_max*(Topology::nexclusion_max + 1) > min_cutoff)
			check_exclusions = 1;
		CalcvdW = &FNC_POT::CalcLJ;
		CalcvdW14 = &FNC_POT::CalcLJ14;
		break;
	}
	case 10:
	{
		pot_dim = nused_potpartials;
		nused_potpartials = RunParams::nused_potpartials;
		break;
	}
	}
	nthreads=RunParams::nthreads;
	f_Coulomb_A=f_Coulomb/SimpleCfg::boxedge;
	//vdW_pot is updated at histogram and its change calculation, so each thread has to have its own segment
	//all the arrays will be created if potential >0, regardless that in case of tabulated potential only vdW is used
	padding_offset=(RunParams::potential>0 ? (int)((CACHE_PADDING-sizeof(*FNC_POT::vdW_pot))/sizeof(*FNC_POT::vdW_pot)) : 0);//number of dummy double elements
	SetArraysize(&vdW_pot, pot_dim*nthreads+(nthreads-1)*padding_offset,"vdW_pot","FNC_POT::SetPotParams");
	SetArraysize(&vdW_pot_s, pot_dim * nthreads + (nthreads - 1) * padding_offset, "vdW_pot_s", "FNC_POT::SetPotParams");
	SetArraysize(&vdW_pot_l, pot_dim * nthreads + (nthreads - 1) * padding_offset, "vdW_pot_l", "FNC_POT::SetPotParams");
	SetArraysize(&vdW_pot_old,SimpleCfg::npartials,"vdW_pot_old","FNC_POT::SetPotParams");
	SetArraysize(&vdW14_pot,pot_dim*nthreads+(nthreads-1)*padding_offset,"vdW14_pot","FNC_POT::SetPotParams");
	SetArraysize(&vdW14_pot_old,SimpleCfg::npartials,"vdW14_pot_old","FNC_POT::SetPotParams");
	SetArraysize(&Coulomb_pot,pot_dim*nthreads+(nthreads-1)*padding_offset,"Coulomb_pot","FNC_POT::SetPotParams");
	SetArraysize(&Coulomb_pot_s, pot_dim * nthreads + (nthreads - 1) * padding_offset, "Coulomb_pot_s", "FNC_POT::SetPotParams");
	SetArraysize(&Coulomb_pot_l, pot_dim * nthreads + (nthreads - 1) * padding_offset, "Coulomb_pot_l", "FNC_POT::SetPotParams");
	SetArraysize(&Coulomb_pot_old,SimpleCfg::npartials,"Coulomb_pot_old","FNC_POT::SetPotParams");
	SetArraysize(&Coulomb14_pot,pot_dim*nthreads+(nthreads-1)*padding_offset,"Coulomb14_pot","FNC_POT::SetPotParams");
	SetArraysize(&Coulomb14_pot_old,SimpleCfg::npartials,"Coulomb14_pot_old","FNC_POT::SetPotParams");

	
	SetArraysize(&bond_pot,Topology::nbond_types,"bond_pot","FNC_POT::SetPotParams");
	SetArraysize(&bond_pot_old,Topology::nbond_types,"bond_pot_old","FNC_POT::SetPotParams");

	SetArraysize(&angle_pot,Topology::nangle_types,"angle_pot","FNC_POT::SetPotParams");
	SetArraysize(&angle_pot_old,Topology::nangle_types,"angle_pot_old","FNC_POT::SetPotParams");
		
	SetArraysize(&perdihedral_pot,Topology::nperdihedral_types,"perdihedral_pot","FNC_POT::SetPotParams");
	SetArraysize(&perdihedral_pot_old,Topology::nperdihedral_types,"perdihedral_pot_old","FNC_POT::SetPotParams");
	
	SetArraysize(&harmdihedral_pot,Topology::nharmdihedral_types,"harmdihedral_pot","FNC_POT::SetPotParams");
	SetArraysize(&harmdihedral_pot_old,Topology::nharmdihedral_types,"harmdihedral_pot_old","FNC_POT::SetPotParams");
	
	SetArraysize(&RBdihedral_pot,Topology::nRBdihedral_types,"RBdihedral_pot","FNC_POT::SetPotParams");
	SetArraysize(&RBdihedral_pot_old,Topology::nRBdihedral_types,"RBdihedral_pot_old","FNC_POT::SetPotParams");
	InitPot();
};

void FNC_POT::InitPot()
{
	int i;
	int padding_offset;

	padding_offset=(RunParams::potential>0 ? (int)((CACHE_PADDING-sizeof(*FNC_POT::vdW_pot))/sizeof(*FNC_POT::vdW_pot)) : 0);//number of dummy double elements
	for (i=0;i<pot_dim*nthreads+(nthreads-1)*padding_offset;i++)
	{
		vdW_pot[i]=0.0;//initializing
		Coulomb_pot[i]=0.0;
		vdW14_pot[i]=0.0;//initializing
		Coulomb14_pot[i]=0.0;
	}
	FNC_POT::vdW_tot_pot=0.0;
	FNC_POT::Coulomb_tot_pot=0.0;
	FNC_POT::vdW14_tot_pot=0.0;
	FNC_POT::Coulomb14_tot_pot=0.0;
	for (i=0;i<Topology::nbond_types;i++)
		bond_pot[i]=0.0;

	for (i=0;i<Topology::nangle_types;i++)
		angle_pot[i]=0.0;
	
	for (i=0;i<Topology::nperdihedral_types;i++)
		perdihedral_pot[i]=0.0;

	for (i=0;i<Topology::nharmdihedral_types;i++)
		harmdihedral_pot[i]=0.0;

	for (i=0;i<Topology::nRBdihedral_types;i++)
		RBdihedral_pot[i]=0.0;
};

//In case of continuation with missing *.hgm file but loaded pot file check, whether the loaded and recalculated
//potential values are in side a tolerable margin of errors
int FNC_POT::CheckPot()
{
	int i;
	double norm;
	char *name;
	SetArraysize(&name,50,"name","FNC_POT::CheckPot");
	
	
	for (i=0;i<Topology::nbond_types;i++)
	{
		if (fabs(bond_pot[i]-0)>TOLERANCE)
			norm=bond_pot[i];
		else
			norm=1.0;
		if (fabs((bond_pot[i]-bond_pot_old[i])/norm)>LOAD_TOL)
		{
			mystrcpy(name,50, "bond");	
			ShowPotError(i,name,bond_pot[i],bond_pot_old[i]);
			return(0);
		}
	}
	
	for (i=0;i<Topology::nangle_types;i++)
	{
		if (fabs(angle_pot[i]-0)>TOLERANCE)
			norm=angle_pot[i];
		else
			norm=1.0;
		if (fabs((angle_pot[i]-angle_pot_old[i])/norm)>LOAD_TOL)
		{
			mystrcpy(name, 50, "angle");
			ShowPotError(i,name,angle_pot[i],angle_pot_old[i]);
			return(0);
		}
	}
	
	for (i=0;i<Topology::nperdihedral_types;i++)
	{
		if (fabs(perdihedral_pot[i]-0)>TOLERANCE)
			norm=perdihedral_pot[i];
		else
			norm=1.0;
		if (fabs((perdihedral_pot[i]-perdihedral_pot_old[i])/norm)>LOAD_TOL)
		{
			if (Topology::force_field_type==OPLSAA)
				mystrcpy(name, 50, "improper (periodic) dihedral");
			else	
				mystrcpy(name, 50, "proper (periodic) dihedral");
			ShowPotError(i,name,perdihedral_pot[i],perdihedral_pot_old[i]);
			return(0);
		}
	}
	
	for (i=0;i<Topology::nharmdihedral_types;i++)
	{
		if (fabs(harmdihedral_pot[i]-0)>TOLERANCE)
			norm=harmdihedral_pot[i];
		else
			norm=1.0;
		if (fabs((harmdihedral_pot[i]-harmdihedral_pot_old[i])/norm)>LOAD_TOL)
		{
			mystrcpy(name, 50, "improper harmonic dihedral");
			ShowPotError(i,name,harmdihedral_pot[i],harmdihedral_pot_old[i]);
			return(0);
		}
	}
	
	for (i=0;i<Topology::nRBdihedral_types;i++)
	{
		if (fabs(RBdihedral_pot[i]-0)>TOLERANCE)
			norm=RBdihedral_pot[i];
		else
			norm=1.0;
		if (fabs((RBdihedral_pot[i]-RBdihedral_pot_old[i])/norm)>LOAD_TOL)
		{
			if (Topology::force_field_type==OPLSAA)
				mystrcpy(name, 50, "proper (RB) dihedral");
			else	
				mystrcpy(name, 50, "Ryckaert-Bellemans dihedral");
		
			ShowPotError(i,name,RBdihedral_pot[i],RBdihedral_pot_old[i]);
			return(0);
		}
	}
		
	if (RunParams::potential>0)
	{
		for (i=0;i<pot_dim;i++)
		{
			
			if (fabs(vdW_pot[i]-0)>TOLERANCE)
				norm=vdW_pot[i];
			else
				norm=1.0;
			if (fabs((vdW_pot[i]-vdW_pot_old[i])/norm)>LOAD_TOL)
			{
				if (RunParams::potential==1)
					mystrcpy(name, 50, "van der Waals");
				else
					mystrcpy(name, 50, "tabulated");//potential=10
				ShowPotError(i,name,vdW_pot[i],vdW_pot_old[i]);
				return(0);
			}
			if (RunParams::potential == 1)
			{
				if (fabs(vdW14_pot[i] - 0) > TOLERANCE)
					norm = vdW14_pot[i];
				else
					norm = 1.0;
				if (fabs((vdW14_pot[i] - vdW14_pot_old[i]) / norm) > LOAD_TOL)
				{
					mystrcpy(name, 50, "van der Waals 1-4");
					ShowPotError(i, name, vdW14_pot[i], vdW14_pot_old[i]);
					return(0);
				}
				if (fabs(Coulomb_pot[i] - 0) > TOLERANCE)
					norm = Coulomb_pot[i];
				else
					norm = 1.0;
				if (fabs((Coulomb_pot[i] - Coulomb_pot_old[i]) / norm) > LOAD_TOL)
				{
					mystrcpy(name, 50, "Coulomb");
					ShowPotError(i, name, Coulomb_pot[i], Coulomb_pot_old[i]);
					return(0);
				}
				if (fabs(Coulomb14_pot[i] - 0) > TOLERANCE)
					norm = Coulomb14_pot[i];
				else
					norm = 1.0;
				if (fabs((Coulomb14_pot[i] - Coulomb14_pot_old[i]) / norm) > LOAD_TOL)
				{
					mystrcpy(name, 50, "Coulomb 1-4");
					ShowPotError(i, name, Coulomb14_pot[i], Coulomb14_pot_old[i]);
					return(0);
				}
			}
		}
		
	}
	return(1);
};

void FNC_POT::ShowPotError(int i,char *potname, double recalc, double loaded)
{
	cout.precision(15);
	cout.setf(ios::right, ios::adjustfield);
	cout.setf(ios::scientific, ios::floatfield);
	cout << "\nWARNING(" << ++warn << "): The loaded ("<<loaded<<") and recalculated ("<<recalc<<") "<<i+1<<". "<<potname<<" potential is outside the tolerable margin of errors!"<<endl;
	cout << "\tMost probably due to a mismatch between the coordination and *.pot file." << endl;
	cout << "\tExact continuation not possible, everything will be recalculated, the run will continue " << endl;
	cout << "\twith the scalable sigmas (if there is any) and move counts from the *.state file!" << endl;
	cout.precision(6);
	cout.unsetf(ios::scientific);
	cout.unsetf(ios::right);
	
};
//updating first the vdW and Coulomb poitential by adding the high and low potential segments (hopefully the high potential values,
//which mostly belonged to the exluded atoms were eliminated), so only small values are added only to the real potetial 
void FNC_POT::UpdateTotPot()
{
	int ipartial;
	vdW_tot_pot=0.0;
	Coulomb_tot_pot=0.0;
	vdW14_tot_pot=0.0;
	Coulomb14_tot_pot=0.0;


	for (ipartial=0;ipartial<pot_dim;ipartial++)
	{
		vdW_pot[ipartial] += vdW_pot_l[ipartial];
		vdW_pot[ipartial] += vdW_pot_s[ipartial];
		
		vdW_tot_pot+=vdW_pot[ipartial];
		if (RunParams::potential == 1)
		{
			Coulomb_pot[ipartial] += Coulomb_pot_l[ipartial];
			Coulomb_pot[ipartial] += Coulomb_pot_s[ipartial];
			Coulomb_tot_pot += Coulomb_pot[ipartial];
			vdW14_tot_pot += vdW14_pot[ipartial];
			Coulomb14_tot_pot += Coulomb14_pot[ipartial];
		}
	}
}

//correcting the non-bonded potential for the atoms belonging to the charge group of the moved atoms
//if the new and the old charge group centre distance to an other charge group minus the Coulomb or vdW cutoff has different
//sign, so the contribution of this charge group for these atom pairs has to be added or removed from the potential
//This is esssentially very similar to the histogram change update (the atoms in the charge group of the moved atoms replacing 
//the moved atoms), first the correction is done regardless the exclusions, then the correction for the exclusions are removed.
//multi-threading is used by splitting the charge groups among the threads 
void FNC_POT::NBPotChargeGrCorr(ThreadArg &thread_arg, Move &move)
{
		int i,j,imoved,iat,jat,jmoved,nnb;
		int igr,jgr;//index of the charge group centres
		int nat_i,nat_j;//number of atoms in the charge groups
		int correction,skip,sign;
		int GRtype1,GRtype2,GRpartialind,partialind,pot_partialind,type1,type2;
		int array_ind,indj;
		int *pindi,*pindj,*p_charge_c_indices,*p_neighbour_ind;
		double temp,*posigr,*posjgr,*posi, *posj,*p_oldchargegr_c1,dcharge_c_sq,dcharge_c_sq_old,d;
		double *my_vdWpot,*my_Coulombpot, *my_s_vdWpot, *my_s_Coulombpot, *my_l_vdWpot, *my_l_Coulombpot;//pointer for this thread potential
		double dx,dy,dz;//coordinate difference
		SimpleCfg &conf=*config;
		
		my_vdWpot=thread_arg.vdW_pot_start;//the starting point of the vdW potential for this thread
		my_Coulombpot=thread_arg.Coulomb_pot_start;//the starting point of the Coulomb potential for this thread
		my_s_vdWpot = thread_arg.vdW_pot_s_start;//the starting point of the vdW potential for this thread for the small contributions
		my_s_Coulombpot = thread_arg.Coulomb_pot_s_start;//the starting point of the Coulomb potential for this thread for the small contributions
		my_l_vdWpot = thread_arg.vdW_pot_l_start;//the starting point of the vdW potential for this thread for the large contributions
		my_l_Coulombpot = thread_arg.Coulomb_pot_l_start;//the starting point of the Coulomb potential for this thread for the large contributions

		//First calculating the correction for all the atom pairs of the atoms in the charge group of the moved atoms
		//It is done by going through the charge groups
		p_charge_c_indices=move.charge_c_indices;//charge group centre index of the moved atom
		p_oldchargegr_c1=move.old_charge_c_pos;//sets the pointer to the first moved atom's charge group centre's old coordinates
		for (imoved=0;imoved<move.tot_moved_atoms;imoved++)//for virtuals not needed, as they are in the same charge group as their moved atom pair
		{
			igr=*p_charge_c_indices;//the charge group centre of this moved atom
			posigr=conf.charge_gr_centre+3*igr;//sets the pointer to the first moved atom's charge group centre's new coordinates
			
			for (jgr=thread_arg.charge_gr_min_ind;jgr<=thread_arg.charge_gr_max_ind;jgr++)//charge group centre index of the neighbour
			{
				if (igr==jgr)
					continue;
				//new distance for the charge group centres
				posjgr=conf.charge_gr_centre+3*jgr;
				dx=*(posigr) - *posjgr++;
				dy=*(posigr+1) - *posjgr++;
				dz=*(posigr+2) - *posjgr++;

				//taking into account the periodical boundary conditions 			
				if (dx>1) dx-=2;
				else if (dx<-1) dx+=2;
				if (dy>1) dy-=2;
				else if (dy<-1) dy+=2;
				if (dz>1) dz-=2;
				else if (dz<-1) dz+=2;
				dcharge_c_sq=dx*dx+dy*dy+dz*dz;

				posjgr-=3;
				//check, whether the jgr charge group does not have a moved atom as well, and consequently its position is changed too
				for (jmoved=0;jmoved<move.nmoved;jmoved++)//here enough to see only the atoms, as the virtual sites are in the same charge group as the moved atoms
				{
					if (move.charge_c_indices[jmoved]==jgr)
					{
						posjgr=move.old_charge_c_pos+3*jmoved;//old position of the charge group
						break;
					}
				}
				//calculating the old distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
				dx=*(p_oldchargegr_c1)- *posjgr++;
				dy=*(p_oldchargegr_c1+1)- *posjgr++;
				dz=*(p_oldchargegr_c1+2)- *posjgr++;

				//taking into account the periodical boundary conditions 			
				if (dx>1) dx-=2;
				else if (dx<-1) dx+=2;
				if (dy>1) dy-=2;
				else if (dy<-1) dy+=2;
				if (dz>1) dz-=2;
				else if (dz<-1) dz+=2;
				dcharge_c_sq_old=dx*dx+dy*dy+dz*dz;
				
				correction=3;
				if ((dcharge_c_sq_old-rundat->Coulomb_cutoff_sq)*(dcharge_c_sq-rundat->Coulomb_cutoff_sq)>0)
					correction--;//no correction necessary, the old and the new charge group distance is either both smaller or larger than cutoff
				if ((dcharge_c_sq_old-rundat->vdW_cutoff_sq)*(dcharge_c_sq-rundat->vdW_cutoff_sq)>0)
					correction-=2;//no correction necessary, the old and the new charge group distance is either both smaller or larger than cutoff
				if (correction==0)
					continue;//go to next charge group

				//correction is necessary, so calculate the distances of the atoms involved
				pindi=atom_for_charge_gr_finder[igr];//pointer to the first atom of this charge group
				
				//natoms_per_charge_group is according to the charge group indices of the topology, not the charge group centres!
				nat_i=Topology::natoms_per_charge_group[charge_group[*pindi]];//number of atoms in this charge group
				for (iat=0;iat<nat_i;iat++)//going through the atoms of the charge group of the moved atom
				{
					//firts check, whether this atom is not a moved atom (the imoved, or an other, if there is more)
					skip=0;
					for (jmoved=0;jmoved<move.tot_moved_atoms+move.tot_moved_virtuals;jmoved++)
					{
						if (*pindi==move.indices[jmoved])
						{
							skip=1;
							pindi++;
							break;//this was already calculated in histoset change update
						}
					}
					if (skip)
						continue;//go to th next atom of the charge group
					
					posi=conf.positions+3* *pindi;

					//find the type (0 to ntypes-1) of the atom
					type1=0;
					while (*pindi>=conf.cumul[type1+1])
						(type1)++;
					//type1 now points to the type of the atom
					//find the GROMACS type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
					GRtype1=0;
					while (*pindi>=SimpleCfg::cumul_GR[(GRtype1) +1])
						(GRtype1)++;
					//GRtype1 now points to the type index of the seperate GROMACS type segment, where the atom is located (the atoms index is continuously increasing)
					//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
					GRtype1=SimpleCfg::GRsep_GR_type[GRtype1];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters

					pindj=atom_for_charge_gr_finder[jgr];//pointer to the first atom of this charge group

					//natoms_per_charge_group is according to the charge group indices of the topology, not the charge group centres!
					nat_j=Topology::natoms_per_charge_group[charge_group[*pindj]];//number of atoms in this charge group
					for (jat=0;jat<nat_j;jat++)//going through the atoms of the other charge group
					{
						//firts check, whether this atom is not a moved atom (the imoved, or an other, if there is more)
						skip=0;
						for (jmoved=0;jmoved<move.tot_moved_atoms + move.tot_moved_virtuals;jmoved++)
						{
							if (*pindj==move.indices[jmoved])
							{
								skip=1;
								pindj++;
								break;//this was already calculated in histoset change update
							}
						}
						if (skip)
							continue;//go to the next atom of the charge group

						//calculate the distance between the atoms, as none of them can be moved it does not changed,  only one distance! 
						posj=conf.positions+3* *pindj;
						dx=*(posi)- *posj++;
						dy=*(posi+1)- *posj++;
						dz=*(posi+2)- *posj++;

						//taking into account the periodical boundary conditions 			
						if (dx>1) dx-=2;
						else if (dx<-1) dx+=2;
						if (dy>1) dy-=2;
						else if (dy<-1) dy+=2;
						if (dz>1) dz-=2;
						else if (dz<-1) dz+=2;
						d=dx*dx+dy*dy+dz*dz;
				
						type2=0;
						while (*pindj>=conf.cumul[type2+1])
						(type2)++;
						//type2 now points to the type of the atom

						//Determining the index of the partial 
						partialind=(type1<=type2 ? (type1*SimpleCfg::ntypes-(type1*(type1+1)/2)+type2) : \
							(type2*SimpleCfg::ntypes-(type2*(type2+1)/2)+type1));
						//Determining the index of the pot_partial, as if virual sites are present, it has to be one of the norma partials
						if (type1 < SimpleCfg::ntypes && type2 < SimpleCfg::ntypes)//normal atoms
							pot_partialind = partialind;
						else
						{
							//for the virtual sites the host atom's type will decide the partialind
							if (type1 >= SimpleCfg::ntypes)
								i = SimpleCfg::virtual_host_type[type1 - SimpleCfg::ntypes];
							else
								i = type1;
							if (type2 >= SimpleCfg::ntypes)
								j = SimpleCfg::virtual_host_type[type2 - SimpleCfg::ntypes];
							else
								j = type2;
							pot_partialind = (i <= j ? (i * SimpleCfg::ntypes - (i * (i + 1) / 2) + j) : (j * SimpleCfg::ntypes - (j * (j + 1) / 2) + i));
						}
						sign=-1;//new distance is above cutoff, old is below, pot correction has to be subtracted
						if (dcharge_c_sq<dcharge_c_sq_old)
							sign=1;//new distance is below cutoff, old is above, pot correction has to be added
						
						if (correction & 1)//correction is 1 or 3, Coulomb has to be updated
						{
							temp = sign * (CalcCoulomb)(*pindi, *pindj, d);//calculating the potential 
							if (fabs(temp) > POT_SPLIT_HIGH)
								*(my_l_Coulombpot + pot_partialind) += temp;
							else if (fabs(temp) < POT_SPLIT_LOW)
								*(my_s_Coulombpot + pot_partialind) += temp;
							else
								*(my_Coulombpot + pot_partialind) += temp;
						}
							
						if (correction>1) //2 (only vdW) or 3 (coulomb and vdW update)
						{
							//find the type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
							GRtype2=0;
							while (*pindj>=SimpleCfg::cumul_GR[(GRtype2) +1])
								(GRtype2)++;
							//GRtype2 now points to the type index of the seperate GROMACS type segmentm where the atom is located (the atoms index is continuously increasing)
							//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
							GRtype2=SimpleCfg::GRsep_GR_type[GRtype2];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters
							//GRpartialind is GROMACS type based partial index, needed to have the vdW parameters, but the histogram is split into RMC-type partials
							GRpartialind=(GRtype1<=GRtype2 ? (GRtype1*Topology::nGRtypes-(GRtype1*(GRtype1+1)/2)+GRtype2) : \
								(GRtype2*Topology::nGRtypes-(GRtype2*(GRtype2+1)/2)+GRtype1));
							temp=sign * (this->*CalcvdW)(GRpartialind, d);//calculating the potential 
							if (fabs(temp) > POT_SPLIT_HIGH)
								*(my_l_vdWpot + pot_partialind) += temp;
							else if (fabs(temp) < POT_SPLIT_LOW)
								*(my_s_vdWpot + pot_partialind) += temp;
							else
								*(my_vdWpot + pot_partialind) += temp;
							
						}
						pindj++;
					}//end of jat cycle
					pindi++;
				}//end of iat cycle
			}//end of charge group cycle jgr
			posigr+=3;
			p_oldchargegr_c1+=3;
			p_charge_c_indices++;
		};//end of imoved cycle for the moved atoms

		//Now go through the exclusions of the atoms in the charge group of the moved atoms
		//This only has to be cheked, if the min(vdW, Coulomb cutoff) < N_exclusions* max(bond length)
		if (check_exclusions)
		{
			p_oldchargegr_c1=move.old_charge_c_pos;//sets the pointer to the first moved atom's charge group centre's old coordinates
			for (imoved=0;imoved<move.nmoved;imoved++)//for all atoms to be moved in the Move (no swaps are alowed in case of molecules)
			{//if there are virtual sites shifted by the move, they belong to the same charge group, so do not needed here as centres
				
				igr=charge_centre[move.indices[imoved]];//the charge group centre of this moved atom
				pindi=atom_for_charge_gr_finder[igr];//pointer to the first atom of this charge group (atom i)
				posigr=conf.charge_gr_centre+3*igr;//sets the pointer to the first moved atom's charge group centre's new coordinates
				nat_i=Topology::natoms_per_charge_group[charge_group[*pindi]];//number of atoms in this charge group
				for (iat=0;iat<nat_i;iat++)//going through the atoms of the charge group of the moved atom
				{
					//firts check, whether this atom is not a moved atom (the imoved, or an other, if there is more
					skip=0;
					for (jmoved=0;jmoved<move.tot_moved_atoms + move.tot_moved_virtuals;jmoved++)//here the virtuals are nedded as well
					{
						if (*pindi==move.indices[jmoved])
						{
							skip=1;
							pindi++;
							break;//this was already calculated in histoset change update
						}
					}
					if (skip)
						continue;//go to th next atom of the charge group

					array_ind=vdW_conv[*pindi];//index in the exclusion array
					p_neighbour_ind=vdW_exclusions+array_ind;//index of the excluded pair (starting with 1)
					posi=conf.positions+3* *pindi;//pointer to the position of the atom
				
					//find the type (0 to ntypes-1) of the atom i
					type1=0;
					while (*pindi>=conf.cumul[type1+1])
						(type1)++;
					//type1 now points to the type of the atom i

					//find the GROMACS type (0 to nsep_GRtypes-1) of the atom
					GRtype1=0;
					while (*pindi>=SimpleCfg::cumul_GR[(GRtype1) +1])
						(GRtype1)++;
					//GRtype1 now points to the type index of the seperate GROMACS type segment, where the atom is located (the atoms index is continuously increasing)
					//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
					GRtype1=SimpleCfg::GRsep_GR_type[GRtype1];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters

					for (nnb=0;nnb<vdW_excl_numb[*pindi];nnb++)//for each excluded pair
					{
						//first to check, whether this is not a not_moved-moved pair already calculated, 
						//among them the central and neighbour can be the same
						skip=0;//flag to indicate, whether the rest of the nnb cycle has to be skipped
						indj=abs(*p_neighbour_ind)-1;//index of the excluded pair;
						pindj=&indj;//just to handle it similarly to the previous section

						for (jmoved=0;jmoved<move.tot_moved_atoms+move.tot_moved_virtuals;jmoved++)
						{
							if (*pindj==move.indices[jmoved])//this pair was aleady calculated 
							{
								//this is a not_moved-moved pair, the rest of the j cycle is skipped
								skip=1;
								p_neighbour_ind++;//skip the index of the neighbour in its own type
								break;
							}
						}
						if (skip)
							continue;//continue with the next moved atom

						jgr=charge_centre[*pindj];//the charge group centre of the exluded pair
					
						if (jgr<thread_arg.charge_gr_min_ind || jgr>thread_arg.charge_gr_max_ind || igr==jgr)
						{
							p_neighbour_ind++;
							continue;//this is the same charge group, no correction needed, or this jgr charge
									//group is handled by another thread
						}

						//calculate the charge group centre distance
						//new distance for the charge group centres
						posjgr=conf.charge_gr_centre+3*jgr;
						dx=*(posigr)- *posjgr++;
						dy=*(posigr+1)- *posjgr++;
						dz=*(posigr+2)- *posjgr++;

						//taking into account the periodical boundary conditions 			
						if (dx>1) dx-=2;
						else if (dx<-1) dx+=2;
						if (dy>1) dy-=2;
						else if (dy<-1) dy+=2;
						if (dz>1) dz-=2;
						else if (dz<-1) dz+=2;
						dcharge_c_sq=dx*dx+dy*dy+dz*dz;
						
						//calculating the old distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
						posjgr-=3;
						//check, whether the jgr charge group does not have a moved atom as well, and consequently its position is changed too
						for (jmoved=0;jmoved<move.nmoved;jmoved++)
						{
							if (move.charge_c_indices[jmoved]==jgr)
								posjgr=move.old_charge_c_pos+3*jmoved;//old position of the charge group
						}
						dx=*(p_oldchargegr_c1)- *posjgr++;
						dy=*(p_oldchargegr_c1+1)- *posjgr++;
						dz=*(p_oldchargegr_c1+2)- *posjgr++;

						//taking into account the periodical boundary conditions 			
						if (dx>1) dx-=2;
						else if (dx<-1) dx+=2;
						if (dy>1) dy-=2;
						else if (dy<-1) dy+=2;
						if (dz>1) dz-=2;
						else if (dz<-1) dz+=2;
						dcharge_c_sq_old=dx*dx+dy*dy+dz*dz;

						correction=3;
						
						if ((dcharge_c_sq_old-rundat->Coulomb_cutoff_sq)*(dcharge_c_sq-rundat->Coulomb_cutoff_sq)>0)
							correction--;//no correction necessary, the old and the new charge group distance is either both smaller or larger than cutoff
						if ((dcharge_c_sq_old-rundat->vdW_cutoff_sq)*(dcharge_c_sq-rundat->vdW_cutoff_sq)>0)
							correction-=2;//no correction necessary, the old and the new charge group distance is either both smaller or larger than cutoff
						if (correction==0)
						{
							p_neighbour_ind++;
							continue;//go to next excluded pair
						}

						//find the type (0 to ntypes-1) of the excluded pair atom
						type2=0;
						while (*pindj>=conf.cumul[(type2) +1])
							type2++;
						//type2 now points to the type of the excluded pair atom

						//Determining the index of the partial 
						partialind=(type1<=type2 ? (type1*SimpleCfg::ntypes-(type1*(type1+1)/2)+type2) : \
						(type2*SimpleCfg::ntypes-(type2*(type2+1)/2)+type1));
						//Determining the index of the pot_partial, as if virual sites are present, it has to be one of the norma partials
						if (type1 < SimpleCfg::ntypes && type2 < SimpleCfg::ntypes)//normal atoms
							pot_partialind = partialind;
						else
						{
							//for the virtual sites the host atom's type will decide the partialind
							if (type1 >= SimpleCfg::ntypes)
								i = SimpleCfg::virtual_host_type[type1 - SimpleCfg::ntypes];
							else
								i = type1;
							if (type2 >= SimpleCfg::ntypes)
								j = SimpleCfg::virtual_host_type[type2 - SimpleCfg::ntypes];
							else
								j = type2;
							pot_partialind = (i <= j ? (i * SimpleCfg::ntypes - (i * (i + 1) / 2) + j) : (j * SimpleCfg::ntypes - (j * (j + 1) / 2) + i));
						}
						//calculating the distance, both atoms not moved
						posj=conf.positions+*pindj*3;//totconf.positions contains the new coordinates of the moved atoms
						
						dx=*(posi)- *posj++;
						dy=*(posi+1)- *posj++;
						dz=*(posi+2)- *posj++;

						//taking into account the periodical boundary conditions 			
						if (dx>1) dx-=2;
						else if (dx<-1) dx+=2;
						if (dy>1) dy-=2;
						else if (dy<-1) dy+=2;
						if (dz>1) dz-=2;
						else if (dz<-1) dz+=2;
						d=dx*dx+dy*dy+dz*dz;
						
						sign=1;//new distance is above cutoff, old is below, pot correction has to be added, as it was subtracted
						if (dcharge_c_sq<dcharge_c_sq_old)
							sign=-1;//new distance is below cutoff, old is above, pot correction has to be subtracted, as it was added
						if (correction & 1)//correction is 1 or 3, Coulomb has to be updated
						{
							temp= sign * (CalcCoulomb)(*pindi, *pindj, d);//calculating the potential 
							if (fabs(temp) > POT_SPLIT_HIGH)
								*(my_l_Coulombpot + pot_partialind) += temp;
							else if (fabs(temp) < POT_SPLIT_LOW)
								*(my_s_Coulombpot + pot_partialind) += temp;
							else
								*(my_Coulombpot + pot_partialind) += temp;
							    
							if (*p_neighbour_ind<0)//1-4 pair
							{
								//as this cannot happen too often, all the threads are updating the same array, and mutex lock is used to
								//to ensure, that only one thread updates the arrays at a time
								std::lock_guard<std::mutex> guard(thread_object->tcl_mutex);//lock mutex
								Coulomb14_pot[pot_partialind]-=RunParams::Coulomb14_fudge*temp;//with opposite sign as the normal Coulomb
							}
						}
						
						if (correction>1) //2 (only vdW) or 3 (coulomb and vdW update)
						{
							//find the type (0 to ntypes-1) of the atom in the first and last instance of the molecule type
							GRtype2=0;
							while (*pindj>=SimpleCfg::cumul_GR[(GRtype2) +1])
								(GRtype2)++;
							//GRtype2 now points to the type index of the seperate GROMACS type segmentm where the atom is located (the atoms index is continuously increasing)
							//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
							GRtype2=SimpleCfg::GRsep_GR_type[GRtype2];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters
							//GRpartialind is GROMACS type based partial index, needed to have the vdW parameters, but the histogram is split into RMC-type partials
							GRpartialind=(GRtype1<=GRtype2 ? (GRtype1*Topology::nGRtypes-(GRtype1*(GRtype1+1)/2)+GRtype2) : \
								(GRtype2*Topology::nGRtypes-(GRtype2*(GRtype2+1)/2)+GRtype1));
							temp= sign * (this->*CalcvdW)(GRpartialind, d);//calculating the potential
							if (fabs(temp) > POT_SPLIT_HIGH)
								*(my_l_vdWpot + pot_partialind) += temp;
							else if (fabs(temp) < POT_SPLIT_LOW)
								*(my_s_vdWpot + pot_partialind) += temp;
							else
								*(my_vdWpot + pot_partialind) += temp;
							if (*p_neighbour_ind<0)//1-4 pair
							{
								//as this cannot happen too often, all the threads are updating the same array, and mutex lock is used to
								//to ensure, that only one thread updates the arrays at a time
								std::lock_guard<std::mutex> guard(thread_object->tcl_mutex);//lock mutex
								vdW14_pot[pot_partialind]-=sign*(this->*CalcvdW14)(GRpartialind,d);
							}
						}
						p_neighbour_ind++;//go to next pair
					}//end of excluded pair cycle nnb
					pindi++;
				}//end of iat cycle for the atoms of the charge group i
				p_oldchargegr_c1+=3;//sets the pointer to the old charge group centre coordinates of the next moved atom
			}//end of imoved cycle for zhe moved atoms
		}//end of if check_exclusion 
	};

	//save the loaded tabulated potential
void FNC_POT::SaveTabulatedPotential(const char *file_name)
{
	int i,j;
	double *r, *u;

	ofstream file;
	
	OpenFile(file, tabpotfilename, "FNC_POT::SaveTabulatedPotential",0);//open file, check, whether it was successfully opened
	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	file.precision(8);
	file << "Tabulated potential is used for " << nused_potpartials << " partial(s)" << endl;
	file << "Teperature: " << RunParams::temperature << " K" << endl;
	
	r = r_pot;
	u = u_pot;
	for (i = 0; i < nused_potpartials; i++)
	{
		file << "\nFor RMC partial " << tabpot_index[i] + 1 << "\n"<<endl;
		file << "Cutoff: " << tab_cutoff[i] * SimpleCfg::boxedge << " A" << endl;
		file << J16 << "r [A]" << J16 << "U(r) [J]" << endl;
		for (j = 0; j < npot_points[i]; j++)
		{
			file << J16 << *r *SimpleCfg::boxedge << J16 << *u++  << endl;
			r++;
		}

	}
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
};

