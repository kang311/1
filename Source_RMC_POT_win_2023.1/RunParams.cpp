//source Runparams.cpp
//Last changed 01.03.2023

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
//extern int aenet_Rc_max;
#include"classes2.h"//classes1.h is included through classes2.h

//Defining the static members
string	RunParams::title;//title
int		RunParams::continuation;//the run should continue based on the configuration, *.dat and *.state file
int		RunParams::runmode;//indicator of program termination
int		RunParams::cfgnumb;//number of configurations to collect after convergence
int		RunParams::coll_frequency;//collect the configuration at each coll_frequency save, if 0, not collected
int		RunParams::printstep;//generated step interval for printing	
int		RunParams::fnc;//this is the fnc switch
int		RunParams::natoms;//the number of atoms in the cell
int		RunParams::ntypes;//the number of atom type (NECESSARY for memory allocation)
int		RunParams::nthreads;//total number of threads to use
int		RunParams::stateversion=3;//to identify the format of the state file, only the same version can be read
					//3: from version 2.3
int		RunParams::ngr;//number of g(r) constraints
int		RunParams::nsq;//number of S(Q) - neutron constraints
int		RunParams::nfq;//number of F(Q) - X-ray constraints
int		RunParams::nfg;//number of F(g) - electron diffraction
int		RunParams::nek;//number of EXAFS constraints
int		RunParams::ntot_datasets;// total number of data sets
int		RunParams::ncosdistr;//number of constraint for the cosine distribution of bond angles
int		RunParams::ncommonneigh;//number of common neighbour constraints
int		RunParams::nsecondneigh;//number of second neighbour constraints
int		RunParams::nbvs;//number of bond valnce sum constraints
int		RunParams::nicoord;//number of individual coordination constraints
int		RunParams::navcoord;//number of average coordination constraints
int		RunParams::ntotal_points;//total number of data points used for output
int		RunParams::potential;//indicator if there is a potential used
int		RunParams::ndiffbin;//number of different bin sizes
int		RunParams::nmoved;//number of atoms to move in a single move
int		RunParams::moveout;//this is a boolean switch indicating whether to use the moveout option
int     RunParams::Rfit;//this is a boolean switch indicating whether to fit the Rfactor instead of chi2 default:0
int		RunParams::exafs_shiftstep=E0_SHIFTSTEP_DEF;//shift the k values in each exafs_shiftstep step with EXAFS data sets with valid MAX-KSHIFT_NGRID
int		RunParams::AXS_shiftstep=AXS_SHIFTSTEP_DEF;//change f' for AXS I(Q) in each AXS_shiftstep. step
int		RunParams::IQ_backg_corr_step = IQ_BACKG_CORR_STEP_DEF;//change the mu for I(Q) background correction for in each IQ_backgcorr_step. step
int		RunParams::reload;//whether to load histogram, and coordination numbers from file
int		RunParams::custmove;//indicator of custom move
int		RunParams::histbuffsize;//size of the history buffer
int		RunParams::histstepratio;//number of saves between each history buffering
int     RunParams::max_gridatom;//maximum number of atoms in a grid cell
int		RunParams::nswap_pairs;//number of mixed partials, where swap can occure
int		RunParams::lead_series_ind;//index of the leading series
int		RunParams::max_option = 10;//maximum number of define options for topology
int		RunParams::n_top_def_option;//number of define options for topology
int		RunParams::offset_lead2 = 0;//indicator how to offset the lead_series_ind2 in case of freeformat reading
int		RunParams::niter_nonlin = NONLIN_MAX_NIT_DEF;//max number of iteration in case of non-linear regression 
int		RunParams::terminate_nonlin = NONLIN_TERMINATE_DEF;//Terminate after this many successive moves ended with failed non-linear regression
int		RunParams::dat_version=1;//the version of the *.dat file
int	   *RunParams::assign_hist;//assigning the different histogram bin sizes to the data sets dim [ntot_datasets]
int    *RunParams::nbins;//the number of histogram bins for each bin size
int    *RunParams::firstbin;//the first histogram bin to use, starts with zero
int    *RunParams::swap_type1;//array containing the type of the first atom of the allowed swap pairs
int    *RunParams::swap_type2;//array containing the type of the second atom of the allowed swap pairs
int    *RunParams::pnatoms = 0;//number of atoms per type	

longint	RunParams::n_generated = 0;//number of generated moves
longint	RunParams::n_swapgen = 0;//number of generated swaps
longint RunParams::n_E0shiftgen = 0;//number of E0 shifts
longint RunParams::n_fprimeshiftgen=0;//number of f' shifts
longint RunParams::n_mucorrgen = 0;//number of generated I(Q) mu correction
longint	RunParams::n_tried = 0;//number of tried moves
longint	RunParams::n_potaccepted = 0;//number of moves accepted based on the chi2_pot 
longint	RunParams::n_accepted = 0;//number of accepted moves
longint	RunParams::n_swapacc = 0;//number of accepted swaps
longint RunParams::n_E0shiftacc = 0;//number of accepted E0 shifts
longint RunParams::n_fprimeshiftacc=0;//number of accepted f' shifts
longint RunParams::n_mucorracc = 0;//number of accepted I(Q) mu correction
longint RunParams::n_goodnonlinreg = 0;//number of tried steps ended with successful non-linear regression
longint	RunParams::last_gen;//number of generated moves at last display,
longint	RunParams::last_tried;//number of tried moves at last display
longint	RunParams::last_accepted;//number of accepted moves at last display,
longint	RunParams::last_swap;//number of generated swaps at last display
longint	RunParams::last_swapacc;//number of accepted swaps at last display
longint RunParams::last_E0shift;//number of generated E0 shifts at last display
longint RunParams::last_E0shiftacc;//number of accepted E0 shifts at last display
longint RunParams::last_mucorracc;//number of accepted mu correction steps at last display
longint RunParams::last_mucorr;//number of generated mu correction steps at last display
longint RunParams::last_fprimeshift;//number of generated f' shifts at last display
longint RunParams::last_fprimeshiftacc;//number of accepted f' shifts at last display

double  RunParams::binshift_def;//default binshift
double  RunParams::epsilon_nonlin= NONLIN_EPSILON_DEF;//stop non-linear iteration, when change in parameter is less than this
double  RunParams::lambda_nonlin= NONLIN_LAMBDA_DEF;//fudge for nonlinera regression
double  RunParams::factor_nonlin= NONLIN_FACTOR_DEF;//factor to change lambda during nonlinear regression
double *RunParams::binshift;//for the calculation of xmin, number of bins to leave out from the calculation at the beginning of the
						//histogram, can be different from the default for the g(r) date sets, where this is automatically detected from the 
						//r point spacing, but only one value is read from the *.dat file
double *RunParams::binshift_diff;//for the different binsize-binshiftcombination, dim [ndiffbin]
double	RunParams::boxedge;//the half length of the simulation cell in real units
double  RunParams::rspacing_def;//the default bin size
double *RunParams::rspacing;//array dim:[ntot_datasets] of the width of the histogram bins as given for each data sets,as  there can be different bin sizes for different data sets
double *RunParams::rspacing_diff;//array dim:[ndiffbin] of the width of the different histogram bins
double *RunParams::gr_rspacing_ori;//for the g(r) data sets, keep the original rpacing values, as it can be reset, if the gr data set's spacing is different, only for pronting run params
double  RunParams::xmax_ori;//the original maximum value for interatomic distances
double *RunParams::xmax;//the maximum value for interatomic distances in reduced units for each bin size, dim [ndiffbin]
double *RunParams::xmin;//the minimum value for interatomic distances in reduced units, dim [ndiffbin]
double  RunParams::too_close_fraction=TOOCLOSE_FRACTION_DEF;//the fraction of moves for too close atoms in case of moveout option
double  RunParams::swap_fraction;//fraction of swaps
double	RunParams::runlimit;//upper running time limit for the run in minutes, or number of accepeted moves to run to, if negative value, it means the number of steps to generate
double	RunParams::timesave;//time interval for saving
double  RunParams::rho;//the number density
string  RunParams::lead_series_name;//name of the lead series
string  RunParams::lead_series_name2;//name of the lead pot series
string *RunParams::chem_symbols=nullptr;//chemical symbol for each type given in *.dat, using the symbols in the Xray_coeffs_W table 
string *RunParams::chem_symbols_standard=nullptr;//standard chemical symbol for each type extracted from chem_symbols, to find them in Compton and N_scat_length

double *RunParams::pcutoff;//the distances of closest approach
double *RunParams::pmaxmove=0;//the maximum moves (1 per atom type)
char   (*RunParams::top_def_option)[FILE_NAME_SIZE];//define option for the topology file, governing which part of it is included
int 	RunParams::auto_cutoff=0;//whether to determine the cutoffs automatically, this is for fixed format input
bool    RunParams::old_out=OLD_OUT_DEF;//Setting writing old out format
bool    RunParams::sum_ppcf=SUM_PPCF_DEF;//Setting summing ppcf's at each save
bool    RunParams::binshift_flag=false;//Indicating, that there are custom bin shift(s)
bool	RunParams::use_custom_sftable=false;//Indicating, that custom sfactortable should be used
bool	RunParams::write_density=false;//if the simulation box was rescaled, write the density in the *.freer(s) files

int		RunParams::nused_potpartials;//number of partials to use for tabulated potential (number of potential files)
int		RunParams::LJ_rep_N=0;//power for LJ repulsion term in 6-N LJ
int		RunParams::nGRtypes=0;//number of different GROMACS types
int		RunParams::NB_weight_mode=0;	//0: the same weighting parameter will be used for all the partials in NB interaction, if there is any, 
									//only one has to be given, 1: different weighting will be used for each partial, npartials weighting parameter has to be given
									//2: the vdW_weight[0] will be used for all potential related contribution, only this is read
int		RunParams::vdW_comb_rule=0;//combination rule for the creation of the mixed partial's potential parameters 
int 	RunParams::tot_used14part=0;//total number of used 1-4 partials
int		RunParams::ntotal_points2=0;//total number of data points used for output in case of non-bonded interactions for the sum of the bonded and non-bonded interactions
int	   *RunParams::used_14partials=0;//array containing 1 for the partials used in 1-4 interactions, 0 otherwise (not to display the unused partials)
double *RunParams::vdW_weight=0;//weighting parameter for the non-bonded interaction 
double *RunParams::Coulomb_weight=0;//weighting parameter for the Coulomb interaction 
double *RunParams::vdW_pot1=0;//LJ_sigma in A or C(6)=4*epsilon*(sigma^6) in kJ/mol*A^6
double *RunParams::vdW_pot2=0;//LJ_epsilon in kJ/mol or C(N)=4*epsilon*(sigma^N) in kJ/mol*A^N
double RunParams::vdW14_fudge=0.0;//this factor will scale the vdW potential, to be able to handle the 1-4 interactions
double RunParams::Coulomb14_fudge=0.0;//this factor will scale the Coulomb potential, to be able to handle the 1-4 interactions
double RunParams::Coulomb_cutoff=0.0;//cutoff value for the Coulomb intercation (A)
double RunParams::vdW_cutoff=0.0;//cutoff value for the vdW intercation (A)
double RunParams::pot_chi2_low_lim_fraction= POT_CHI_LOW_LIM_FRACTION_DEF;//the fraction of the initial total chi2 wich would serve 
									//as the lower limit, which the pot chi2 could not go under
double   RunParams::temperature;//for the tabulated potential, only for output

int RunParams::nbond_types=0;//number of bond types
int RunParams::nangle_types=0;//number of angle types
int RunParams::nperdihedral_types=0;//number of periodic dihedrals types
int RunParams::nharmdihedral_types=0;//number of harmonic dihedrals types
int RunParams::nRBdihedral_types=0;//number of RB dihedrals types	
int RunParams::lead_series_ind2=0;//index of the leading series for the potential sigma scaling, starting with the vdW interaction
char  RunParams::vdW_name[15];//name of the vdW potential


std::list<string> RunParams::cus_entry;//Custom move entries

#ifdef _AENET
	int	RunParams::aenet_step;//perform ANN pot calculation in each aenet_step step
	double RunParams::aenet_cutoff=-1;//cutoff for the ANN potential
	double RunParams::aenet_weight;//weight parameter
	bool RunParams::write_E;//whether write the atomic energies to enerfilename
	bool RunParams::do_relax;//whether do structure relaxation for the atoms inside cutoff around the moved atom
#endif

#ifdef _LOCAL_INV
	int		RunParams::nlocint;//number of local invariance intervals
	int		*RunParams::loc_npoints;//number of data points for local invariance (depend on the calculation mode)
	int		RunParams::loc_chi2_mode;//how to calculate the local invariance chi2, 0 based on bins, 1 based on distance
	double	*RunParams::loc_inv_sigma;//weight for the local invariance
	double *RunParams::loc_at_ratio;//ratio of atoms at the interval bordersfor  the local invariance calculation in case of loc_chi2_mode=1 (nlocint+1)
#endif
#ifdef _USE_LOCAL_INV
	int		RunParams::min_loc_nbins;//first normal histogram bin index to use in the local invariance 
	int		RunParams::max_loc_nbins;//number of bins to use for local invariance calculation 
	double  RunParams::loc_rspacing;//the width of the local invariance histogram bins
	double	RunParams::min_loc_r;//minimum distance for the local invariance calculation in reduced unit
	double	RunParams::max_loc_r;//maximum distance for the local invariance calculation in reduced unit
#endif
	
#ifdef _NO_PERIODIC
	int RunParams::recentre_flag;//indicator, whether recentre the simulation box
	double	RunParams::R0;//the radius of the speric sample inside the box
#endif
#ifdef _VIBR_AMP
	double *RunParams::T_corr_sigma_t;//sigma parameter for the Gaussian convolution function for each type
#endif

RunParams::RunParams()//default constructor
{
	//this constructor allocates space for the non-static run parameters data
	int i;
	int npartials = ntypes * (ntypes + 1) / 2;//for number of cutoffs

	rmax = 0;
	fnc_conflict = 0;//sets it to default no conflict
	SetArraysize(&pcutsq, npartials, "pcutsq", "RunParams::RunParams");//the distances of closest approach squared

	for (i = 0; i < npartials; i++)//initialising the array of cutoffs squared and reduced
		pcutsq[i] = pcutoff[i] * pcutoff[i] / boxedge / boxedge;
	//this is set here for LJ-Coulomb, otherwise it can be reset later
	Coulomb_cutoff_sq = pow(Coulomb_cutoff / boxedge, 2.0);//to get it in reduced unit
	vdW_cutoff_sq = pow(vdW_cutoff / boxedge, 2.0);

}

//------------setting static parameters concerning the configurations
void RunParams::SetParams(int n_atoms,int n_types, double box_edge, int *p_natoms)
{
	int i;
	ntypes=n_types;
	natoms=n_atoms;
	boxedge=box_edge;
	rho=n_atoms/(box_edge*box_edge*box_edge*8);//number density
	top_def_option=new char[max_option][FILE_NAME_SIZE];
	if (top_def_option==NULL)
		NoArray("RunParams::SetParams","top_def_option");
	#ifdef _AENET
		Aenet::SetParams();
	#endif
	

	SetArraysize(&pnatoms,ntypes,"pnatoms","RunParams::SetParams");//number of atoms/type
	
	for (i=0;i<ntypes;i++)
		pnatoms[i]=p_natoms[i];
	
#ifdef _LOCAL_INV	
	SetArraysize(&loc_npoints,nlocint,"loc_npoints","RunParams::SetParams");
#endif
#ifdef _VIBR_AMP
	SetArraysize(&T_corr_sigma_t,ntypes,"T_corr_sigma_t","RunParams::SetParams");
#endif
}

//determines, how many different binsize has to be used, and assign it to the data sets
//The concept is that only as many different bin sizes will be handled, as it is absolutely necessary. 
//As in case of g(r) data sets, the r point spacing and binshift is automatically detected, it can happen, that even in case of the 
//same binsize the binshift is different, so different histograms have to be used! In case of the same binsize but different binshifts
//it is examined, whether it is actually the same histogram with different starting bin. (For example binshift 0, 1, 2... are the same
//and will be considered as binshift 0, or 0.5, 1.5 and 2.5... will be considered as binshift 0.5. So only binshift<1 will occure for the different histograms! 
//So where binshifts with integer difference will share the same histogram, and only the first valid bin for the given data set is determined.
//To make thing easier, the histogram will be calculated from the first possible binshift (so if for example fom the 0, 1 and 2 series actually 1
//is the smallest used, it will be calculated for 0 as well.

//The different histograms will be determined here, and their index will be assigned to the dat sets.
//if possible, the default binsize with default binshift will be the first, if it is not used, which can happen, then the first data set's binsize will be the first
void RunParams::AssignBinsize()
{
	
	int i, j,index = -1, index2=-1,first;
	int newbinsize;
	for (i = 0; i < ntot_datasets; i++)
	{
		if (fabs(rspacing[i] - rspacing_def) < TOLERANCE && fabs(binshift_def-binshift[i])<TOLERANCE)
		{
			index = i;//find the first default
			break;
		}
		else
		{
			if (fabs(rspacing[i] - rspacing_def) < TOLERANCE)
			{
				if (index2==-1)
					index2 = i;//find the first default
			}
		}
	}
	if (index == -1 && index2 > -1)//there is none with both default, but there is one with def binsize, use that
		index = index2;

	SetArraysize(&rspacing_diff, 1, "rspacing_diff", "RunParams::AssignBinsize");
	SetArraysize(&binshift_diff, 1, "binshift_diff", "RunParams::AssignBinsize");
	if (index == -1)//no default binsize with default binshift is used, all the data sets have custom bin sizes or binshits
	{
		first = 0;
		rspacing_diff[0] = rspacing[0];//set the size of the first unique bin size
		binshift_diff[0] = binshift[0]-firstbin[0];
	}
	else
	{
		first = index;//default bin size is used, this is the first occurence
		rspacing_diff[0] = rspacing[index];//set the size of the first unique bin size, which is the default
		binshift_diff[0] = binshift[index] - firstbin[index];
	}
	//count the number of different binsizes
	ndiffbin = 1;
	assign_hist[0] = 0;//in case there is no data sets
	for (i = 0; i < ntot_datasets; i++)
	{
		if (first == i)
		{
			assign_hist[i] = 0;
			continue;
		}
		//see, whether this bin size and binshift  occured before
		newbinsize = 1;
		for (j = 0; j < ndiffbin; j++)
		{
			if (fabs(rspacing[i] - rspacing_diff[j]) < TOLERANCE2)
			{
				if (fabs(binshift[i]-firstbin[i] - binshift_diff[j]) < TOLERANCE2)//this is not a new bin size)
				{
					newbinsize = 0;
					assign_hist[i] = j;//index of the binsize
					break;
				}
			}
		}
		if (newbinsize)//this is new unique binsize
		{
			assign_hist[i] = ndiffbin;//index of the binsize
			ResizeArray( &ndiffbin, ndiffbin + 1, &rspacing_diff, "rspacing_diff", "RunParams::AssignBinsize");
			ndiffbin--;
			ResizeArray(&ndiffbin, ndiffbin + 1, &binshift_diff, "binshift_diff", "RunParams::AssignBinsize");
			rspacing_diff[ndiffbin-1] = rspacing[i];
			binshift_diff[ndiffbin-1] = fabs(binshift[i]-firstbin[i]);
		}
	}
}

//-------------Setting the limits for the histogram calculation--------------
void RunParams::SetXmax()
{
	int i,dim=ndiffbin;//default

	xmax_ori = xmax[0];

	rmax = 0;
#ifdef _NO_PERIODIC
	dim = 1;//there cannot be different bin sizes
#endif
#ifdef _VIBR_AMP
	dim = 1;//there cannot be different bin sizes
#endif

	SetArraysize(&xmin, dim, "xmin", "RunParams::SetXmax");
	
	
	if (dim > 1)
	{
		i = 1;
		ResizeArray(&i, dim, &nbins, "nbins", "RunParams::SetXmax");
		i = 1;
		ResizeArray(&i,dim,&xmax,"xmax", "RunParams::SetXmax");
	}
#ifndef _NO_PERIODIC
	for (i = 0; i < dim; i++)
	{
		xmin[i] = binshift_diff[i] * rspacing_diff[i] / boxedge;//keep in mind, that binshift_diff is always below 1, the lowest binshift for the series
		//of histograms it represetnt. For example for the 0, 1, 2... binshift series binshift=0, for the 0.5, 1.5, 2.5 ... series 0.5
		if (fabs(xmin[i]) < TOLERANCE2)
			xmin[i] = 0.0;															

		//to be able to use the same histogram, if the custom histograms only differ in which first bin to use, the histogram will be calculated from the start
		//binshift[i]-firstbin[i], and for each data set to use the correct first bin firstbin[i] will be addded, when the histogram is positioned 
		nbins[i] = int((xmax_ori - xmin[i])*boxedge / rspacing_diff[i]) + 1;
		//adjusting xmax so that it coincides with the upper limit of the last bin
		xmax[i] = xmin[i] + nbins[i] * rspacing_diff[i] / boxedge;
		while (xmin[i] + (nbins[i] - 1)*rspacing_diff[i] / boxedge > SQRT3)//the last bin should should not be over the sqrt(3), the maximum value in sftable
		{
			nbins--;
			xmax[i] = xmin[i] + nbins[i] * rspacing_diff[i] / boxedge;
		}
		
		rmax = (xmax[i] > rmax ? xmax[i]:rmax);//finding the largest
		cout << "\nNOTE(" << ++note << "): xmax was adjusted to " << xmax[i];
		cout << ", number of histogram bins " << nbins[i] << endl;
	}
	rmax *= boxedge;
	cout << "\nNOTE(" << ++note << "): The largest distance to calculate for, rmax is " << rmax << " Angstroms" << endl;
#else
	//rspacing and number of bins were already determined at reading based on the sample size. It will not be reset here, the whole sample size 
	//will be used. xmax can only be 0.5,it is checked at reading
	xmin[0] = binshift_diff[0] * rspacing_diff[0] / boxedge;
	min_loc_r=xmin[0];
	max_loc_r=xmax[0];
	max_loc_nbins=nbins[0];
#endif

#ifdef _LOCAL_INV
	//first check, that the local inv first and last bin to use are in the normal bin range, and set it, if necessary
	if (min_loc_r>=max_loc_r)
	{
		cout << "\nWARNING(" << ++warn << "): The minimum reduced distance for the local invariance calculation is larger, than the maximum distance. "<<endl;
		cout<<"\tIt wlll be set to the same as the default normal histogram!"<<endl;
		min_loc_r=xmin[0];
	}
	if (min_loc_r<xmin[0])
	{
		cout << "\nWARNING(" << ++warn << "): The minimum reduced distance for the local invariance calculation is smaller,"<<endl;
		cout<<"\tthan minimum distance for the normal histogram calculation! "<<endl;
	}
	if (max_loc_r>xmax[0])
	{
		cout << "\nWARNING(" << ++warn << "): The maximum reduced distance for the local invariance calculation is larger, "<<endl;
		cout<<"\tthan maximum distance for the normal histogram calculation! "<<endl;
	}
		
	min_loc_nbins=RoundNearInt( min_loc_r*boxedge/loc_rspacing);
	min_loc_r=min_loc_nbins*loc_rspacing/boxedge;
	max_loc_nbins=int((max_loc_r-min_loc_r)*boxedge/loc_rspacing)+1;
	max_loc_r=min_loc_r+max_loc_nbins*loc_rspacing/boxedge;
	
	while (min_loc_r+(max_loc_nbins-1)*loc_rspacing/boxedge>SQRT3 )//the bin before the last bin should should not be over the sqrt(3)
	{
		max_loc_nbins--;
		max_loc_r=min_loc_r+max_loc_nbins*loc_rspacing/boxedge;
	}

#endif
}
//set the total number of data points
void RunParams::SetPoints()
{
	//ntotal points will be updated with the bonded interactions by ChiSquared::SetChiSquaredParams
	ntotal_points = ExptsData::rused_tot + ExptsData::sqused_tot + ExptsData::fqused_tot + ExptsData::fgused_tot + ExptsData::ekused_tot + ncosdistr * CosDistrConst::ncos_bin + CoordNumbConst::tot_subconst + navcoord;
	//the ntotal_point will be updated by Chisquare::SetChisquareParams
#ifdef _ADVANCED_GEOM_CONST
	ntotal_points += CommonNeighConst::nconstraints + SecondNeighConst::nconstraints + BondValenceSumConst::nconstraints;
#endif
}
//-----reading logical parameters from file, can be given in new and old format-----
//New format: 0 for false, 1 for true, old format .false. or .true..
//if mode=1, then other integer value than 0 and 1 is read, only for potential
//if mode is 0 or even number, then go to new line, otherwise leave file pointer as it is after reading
//if mode is 3 logical reading (0,1,true, false not to go to new line
void RunParams::ReadLogical(ifstream &file, int &param, const char *param_name, const char *routine_name, int mode)
{
	char dumb[10];
	char *pdumb;
	int repeat; 
	longint length;
	std::streamoff file_pos;

	file_pos = file.tellg();
	file >> dumb;
	length = strlen(dumb);
	if (length > 1 && mode == 1 && *dumb!='.')
	{
		file.seekg(file_pos, ios::beg);//setting back the file pointer, before reading
		param = ReadThisLine(file, 3, 0, "potential", "RunParams::ReadLogical", datfilename);
		return;
	};
	pdumb = dumb;
	do
	{
		repeat=0;//do not repeat by default
		switch (*pdumb)
		{
			case '0'://new format
			case 'f'://old format
			{
				param=0;
				break;
			}
			case '1'://new format
			case 't'://old format
			{
				param=1;
				break;
			}
			default:
			{
				if (*pdumb=='.')
				{
					//most probably old format beginning with ., read another character
					pdumb++;
					repeat=1;
				}
				else
				{
					if (mode==1)
					{
						//this is basically an integer entry, but boolean is accepted for compatibility reasons
						file.seekg(file_pos,ios::beg);//setting back the file pointer, before reading
						param=ReadThisLine(file,3,0,"potential","RunParams::ReadLogical",datfilename);
					}
					else
					{
						if (::debug)
						{
							//Neither new or old format is used
							cout << "\nWARNING(" << ++warn << "): RunParams::ReadLogical does not read a " << param_name << " switch" << " in the expected format! " << endl;
							cout << "\tThe format can be 0 for false or 1 for true! Setting the switch value to 0!" << endl;
						}
						param=0;
					}
				}
			}
		}//end of switch dumb
	}
	while (repeat);
	
	if (mode%2 ==0)
		SkipLine(file);
};
//needed if an optional logical is read, and not given, and there is no comment
//-----reading logical parameters from string stream, can be given in new and old format-----
//New format: 0 for false, 1 for true, old format .false. or .true..
void RunParams::ReadLogicalS(stringstream &ss, int &param, const char *param_name, const char *routine_name)
{
	string dumb;
			
	ss >> dumb;
	if ( dumb == "0")	param = 0;
	else if (dumb == "false")	param = 0;
	else if (dumb == ".false.")	param = 0;
	else if (dumb == "1") param = 1; 
	else if (dumb == "true") param = 1;
	else if (dumb == ".true.") param = 1;
	else
	{
		//Neither new or old format is used
		cout << "\nWARNING(" << ++warn << "): RunParams::ReadLogicalS does not read a " << param_name << " switch" << " in the expected format! " << endl;
		cout << "\tThe format can be 0 for false or 1 for true! Setting the switch value to 0!" << endl;
		param = 0;
		return;
	}

	

};
//------------reading .dat file to get the satic parameters
void RunParams::GetParams()
{
	ifstream file;
	string line;
	SafeOpenTextFile(file,datfilename);
	if (CheckFileState(file,"RunParams::GetParams",datfilename)==0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};
	

	
	SetArraysize(&xmax, 1, "xmax", "RunParams::GetParams");
	SetArraysize(&nbins, 1, "nbins", "RunParams::GetParams");
			
	getline(file,line);
	string datversion = line.substr(0, 4);
	if (datversion == "#002")//this is the older free fromat *.dat, where the X-ray coeff reading was the default
		dat_version = 2;
	else if (datversion == "#003")//this is the newer free fromat *.dat, where the X-ray coeff calculation is the default
		dat_version = 3;
	if (dat_version==1)
	{
		Ltrim(line);
		Rtrim(line);
		title = line;
		GetFixedFormatParams(file);
	}
	else GetFreeFormatParams(file);
	next_line_pos = -1;
	file.close();
}

//This function rescales the boxedge
void RunParams::AdjustBoxedge()
{
	
#ifndef _NO_PERIODIC
	int intbox, temp;
	double diff;
	
	//here the box will be rescaled according to the density value given in the *.dat file
	boxedge=pow(natoms/rho,1.0/3.0)/2;
	diff=fabs(boxedge-SimpleCfg::boxedge);
	//round the value to 6 digits, as this will be saved and further used, if the run continued, and cause discrapencies
	intbox=(int)boxedge;
	temp=(int)((boxedge-intbox)*1e6);
	if ((boxedge-intbox)*1e6 -(double)temp>0.5)
		temp++;
	boxedge=intbox+temp/1e6;
	if (diff>5e-7)
	{
		write_density = true;
		cout.precision(16);
		cout << "\nWARNING(" << ++warn << "): The boxedge is " << SimpleCfg::boxedge << " A in the configuration file and " << boxedge << " A" << endl;\
		cout<<"\tcalculated from the number density given in the *.dat file are not the same!"<<endl;
		cout<<"\tThe value from the *.dat file will be used, which means rescaling of the simulation box!"<<endl;
		SimpleCfg::boxedge=boxedge;
	}
	//Even if the density given in *.dat and calculated from the data from the configuration file agree, the rho will bw recalculated
	//from the data in the configuration file to be the same using the fixed and free format data reading
	rho = SimpleCfg::ntotal / pow(2 * SimpleCfg::boxedge, 3);
#endif
  return;
}


void RunParams::GetFixedFormatParams(ifstream &file)
{
	int i, j, npartials;
	int temp = 0;
	int size = 1;
	int cont;
	char *name;
	char *conv_numb = NULL;
#ifdef _LOCAL_INV
	char *text;

#endif
#ifdef _NO_PERIODIC
	int ineightype;
	double dtemp;
	double v_R0;
#endif

	SetArraysize(&name, NAME_SIZE, "name", "RunParams::GetFixedFormatParams");
	rho = ReadThisLine(file, 1, 0.0, "rho", "RunParams::GetFixedFormatParams", datfilename);//reading the number density from the *.dat file
#ifndef _NO_PERIODIC
	AdjustBoxedge();//temp is not needed in this case
#endif
	//EXAFS E0 shift can only be used in case of free format input, only default value is set
	exafs_shiftstep = E0_SHIFTSTEP_DEF;
	//cutoffs
	if (RunParams::ntypes == 0)
	{
		cout << "\n" << "*****ERROR*****" << endl;
		cout << "RunParams::GetFixedFormatParams " << endl;
		cout << "Number of types is 0 in RunParams object" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	npartials = ntypes * (ntypes + 1) / 2;
	SetArraysize(&pcutoff, npartials, "pcutoff", "RunParams::GetFixedFormatParams");
	string myline;
	getline(file, myline);
	stringstream ss(myline);

	j = 0;
	int k = 0;//index in the row
	for (i = 0; i < npartials; i++)
	{
		ss >> pcutoff[i];

		if (ss.fail())
		{
			if (k != ntypes -j)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Variable cutoff was not given either in one line or in upper triangular form, " << endl;
				cout << npartials << " values are expected for the " << ntypes << " RMC types" << endl;
				CleanExit();
			}
			getline(file, myline);
			j++;
			k = 0;
			ss.str("");
			ss.clear();
			ss << myline;

			ss >> pcutoff[i];
			if (ss.fail())
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Variable cutoff[" << i + 1 << "] could not be read in RunParams::GetExptParams, probably end of line was reached!" << endl;
				CleanExit();
			}
		}
		k++;
	}



	//maximal atomic movements /type
	SetArraysize(&pmaxmove, ntypes, "pmaxmove", "RunParams::GetFixedFormatParams");
	for (i = 0; i < ntypes; i++)//the last value will be read differently, to make sure there is still data in the line
	{
		if (i < ntypes - 1)
			pmaxmove[i] = ReadThisLine(file, 3, MAX_MOVES_DEF, "pmaxmove", "RunParams::GetFixedFormatParams");
		else
			pmaxmove[i] = ReadThisLine(file, 1, MAX_MOVES_DEF, "pmaxmove", "RunParams::GetFixedFormatParams", datfilename);
	}

	i = 1;//go to new line after rspacing reading
#ifdef _LOCAL_INV
	i = 3;//stay in the same line
#endif
	rspacing_def = ReadThisLine(file, i, RSPACING_DEF, "rspacing_def", "RunParams::GetFixedFormatParams", datfilename);//the width of the histogram bins
#ifdef _LOCAL_INV
	loc_rspacing = ReadThisLine(file, 6, rspacing_def, "loc_rspacing", "RunParams::GetFixedFormatParams", datfilename);//the width of the local invariance histogram bins
#else
#ifdef _USE_LOCAL_INV
	loc_rspacing = rspacing_def;
#endif
#endif
	//this is the moveout switch (cannot be read directly, as it can be in the old format .false. or .true.,
	//but it can also be 0 for false and 1 for true)
	ReadLogical(file, moveout, "moveout", "RunParams::GetFixedFormatParams");

	cfgnumb = ReadThisLine(file, 3, CFG_COLL_DEF, "cfgnumb", "RunParams::GetFixedFormatParams");//number of configurations to collect
	coll_frequency = ReadThisLine(file, 6, CFG_COLL_FREQ_DEF, "cfgnumb", "RunParams::GetFixedFormatParams", datfilename);//collect the configuration at each coll_frequency save, if 0, not collected

	if (coll_frequency <= 0)
	{
		cfgnumb = 0;//configurations will not be collected
		coll_frequency = 0;
	}

	printstep = ReadThisLine(file, 1, PRINTSTEP_DEF, "printstep", "RunParams::GetFixedFormatParams", datfilename);//generated step interval for printing
	runlimit = ReadThisLine(file, 3, TIMELIM_DEF, "runlimit", "RunParams::GetFixedFormatParams");//time limit for the run
	if (runlimit < 0)
	{
		runmode = 1;//run to -runlimit ngenerated steps
		runlimit=round(runlimit);
	}
#ifdef _TEST_MODE
	else
	{
		runmode = 2;//to be the same as with freeformat
		runlimit=round(runlimit);
	}
#endif
	timesave = ReadThisLine(file, 1, TIMESAVE_DEF, "timesave", "RunParams::GetFixedFormatParams", datfilename);//time interval for saving

	ngr = ReadThisLine(file, 3, 0, "ngr", "RunParams::GetFixedFormatParams");//number of g(r) data sets
	nsq = ReadThisLine(file, 3, 0, "nsq", "RunParams::GetFixedFormatParams");//number of S(Q) - neutron data sets
	nfq = ReadThisLine(file, 3, 0, "nfq", "RunParams::GetFixedFormatParams");//number of F(Q) - X-ray data sets
#ifdef  _READ_EDIFF //only read it, if compiled this way, to be able to use the old *.dat files otherwise
	nfg = ReadThisLine(file, 3, 0, "nfg", "RunParams::GetFixedFormatParams");//number of F(g) - electron difraction data sets
#endif
	nek = ReadThisLine(file, 3, 0, "nek", "RunParams::GetFixedFormatParams");//number of EXAFS data sets
	lead_series_ind = ReadThisLine(file, 6, LEAD_SERIES_IND_DEF, "lead_series_ind", "RunParams::GetFixedFormatParams", datfilename);//index of the lead series in case of sigma calculation

 // if _VIBR_AMP and/or _NO_PERIODIC is used, then only one bin size is allowed and will be used, but the array is created and filled anyway
	ntot_datasets = ngr + nsq + nfq + nfg + nek;
	if (ntot_datasets > size)
		size = ntot_datasets;
	SetArraysize(&rspacing, size, "rspacing", "RunParams::GetFixedFormatParams");


	for (i = 0; i < size; i++)
		rspacing[i] = rspacing_def;//fill it with the default binsize, it will be reset, where necessary

	SetArraysize(&assign_hist, size, "assign_hist", "RunParams::AssignBinsize");

	//read the parameters for the data constraints
	ExptsData::GetExptParam(file);

	ncosdistr = ReadThisLine(file, 3, 0, "ncosdistr", "RunParams::GetFixedFormatParams");//number of coordination constrains
	if (ncosdistr)
		CosDistrConst::GetCosDistrConst(file);
	else
		SkipLine(file, datfilename, 1);

	nicoord = ReadThisLine(file, 3, 0, "nicoord", "RunParams::GetFixedFormatParams");//number of coordination constrains
	if (nicoord > 0)
		CoordNumbConst::GetCoordConst(file);
	else
		SkipLine(file, datfilename, 1);

	navcoord = ReadThisLine(file, 1, 0, "navcoord", "RunParams::GetFixedFormatParams", datfilename);//number of average coordination constraint
	if (navcoord > 0)
		AvCoordConst::GetAvCoordConst(file);

#ifdef _ADVANCED_GEOM_CONST
	ncommonneigh = ReadThisLine(file, 3, 0, "ncommonneigh", "RunParams::GetFixedFormatParams", datfilename);//number of common neighbor constraint
	if (ncommonneigh > 0)
		CommonNeighConst::GetCommonNeighConst(file);
	else
		SkipLine(file);
	//form version 1.7 the fixed format input is not possible for the new features, only the free format
#endif

	//this is the potential switch (cannot be read directly, as it can be in the old format .false. or .true.,
	//but it can also be 0 for false and 1 for true)

	ReadLogical(file, potential, "potential", "RunParams::GetFixedFormatParams", 1);

	lead_series_ind2 = -1;//default

	if (potential > 0)
		GetPotParams(file);//getting the potential related parameters
	else
		SkipLine(file);

	getline(file, myline);
	ss.str("");
	ss.clear();
	ss << myline;
	ss >> fnc;
	if (ss.fail())
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Variable fnc could not be read from " << datfilename << " file in routine RunParams::GetFixedFormatParams" << endl;
		cout << "Exiting..." << endl;
		CleanExit();
	}


	if (potential > 0 && (fnc > 0 && fnc < 4))
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Non-bonded potential cannot be used with standard FNC constraint!" << endl;
		cout << "Use fnc with flexible molecules and proper topology file, if you want to use non-bonded potential!" << endl;
		cout << "Cannot run this way, exiting!!!" << endl;
		CleanExit();
	}
#ifdef _USE_LOCAL_INV
	if (potential > 0 || fnc == 4)
	{
		cout << "\n*****ERROR*****" << endl;
#ifdef _LOCAL_INV
		cout << "LOCAL INVARIANCE ";
#ifdef _NO_PERIODIC
		cout << "and ";
#endif
#endif

#ifdef _NO_PERIODIC
		cout << "PERIODIC BOUNDARY CONDITIONS ";
#endif
		cout << "cannot be used with non-bonded or bonded potential!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
#endif
    if (fnc==4)
	{
		//trying to read the define optopns which should be passed to the topology file processing
		cont=1;
	 
		do
		{
			  ss>>top_def_option[n_top_def_option];
			  if (strlen(top_def_option[n_top_def_option]) > FILE_NAME_SIZE)
			  {
				  cout << "\n*****ERROR*****" << endl;
				  cout << "The size of the " << n_top_def_option + 1 << ". topology define option (" << top_def_option[n_top_def_option] << ") in the " << endl;
				  cout << datfilename << " file exceeds " << FILE_NAME_SIZE << "!" << endl;
				  cout << "Either use shorter define options or change the value of FILE_NAME_SIZE parameter in units.h and recompile the program!" << endl;
				  cout << "Cannot run this way, exiting..." << endl;
				  CleanExit();
			  }
		  
			  //normally, the comments begins with an ! in the *.dat file, but not necessarily, and a valid option has to begin with -D
			  if (strncmp(top_def_option[n_top_def_option],"!",1)==0 || strncmp(top_def_option[n_top_def_option],"-D",2)!=0)
			  {
				  cont=0;//this is not a valid option, delete it, and finish option reading
				  n_top_def_option--;
			  }
			  else
			  {
				  n_top_def_option++;
				  if (n_top_def_option>max_option-1)
				  {
					  //resize array, it will be resized, as a simple char array
					  //it is done in two steps to prevent a g++ 4.4.5 compiler warning 
					  ResizeArray(&max_option,max_option+3,(char**)(&top_def_option),"top_def_option","RunParams::GetFixedFormatParams",FILE_NAME_SIZE);
					  ResizeArray(&max_option,max_option+2,(char**)(&top_def_option),"top_def_option","RunParams::GetFixedFormatParams",FILE_NAME_SIZE);
				  
				  }
			  }
		  
		}
		while (cont);
		n_top_def_option++;//this is the number of valid options
	  //check, whether the last option does not contain the ! and some othe characters after it, if no space was given between it and the comment, and truncateit if necessary
		char *pos_char;
		pos_char =  strchr(top_def_option[n_top_def_option - 1], '!');
		if (pos_char!=NULL && (pos_char- top_def_option[n_top_def_option-1]<(int)strlen(top_def_option[n_top_def_option - 1])))
			*pos_char='\0';
	}
	
	binshift_def=ReadThisLine(file,1,BINSHIFT_DEF,"binshift", "RunParams::GetFixedFormatParams",datfilename);//binshift (number of bins to leave out of the calculation)
	SetArraysize(&binshift, size, "binshift", "RunParams::GetFixedFormatParams");
	SetArraysize(&firstbin, size, "firstbin", "RunParams::GetFixedFormatParams");
	//not to complicate thing even more, in case of _NO_PERIODIC binshift has to be below 1
	j = 0;
	if (fabs(binshift_def) < LOAD_TOL)
		binshift_def = 0;
	else
	{
		//as due to rounding error in case of binshift 1, 2... the (int)binshift is not necessarily correct, an integer firstbin will be calculated here
		if (fabs(binshift_def - round(binshift_def)) < LOAD_TOL)//this supposed to be integer
			j = (int)round(binshift_def);
		else
			j = (int)binshift_def;
	}
#ifdef _NO_PERIODIC
	if (j> 0)
	{
		cout << "\nWARNING(" << ++warn << "): In case of no periodic boundary conditions binshift cannot be larger than 1.0!" << endl;
		cout << "\tIt will be reset from " << binshift_def << " to " << binshift_def - j <<"!"<< endl;
		binshift_def -= j;
		j = 0;
	}
#endif
	for (i = 0; i < size; i++)
	{
		binshift[i] = binshift_def;
		firstbin[i] = j;
	}
	xmax[0]=ReadThisLine(file,1,XMAX_DEF,"xmax[0]", "RunParams::GetFixedFormatParams",datfilename);//maximum value for histogram calculation (in reduced units)

	ExptsData::CheckgrSets();

	nmoved=ReadThisLine(file,1,NMOVED_DEF,"nmoved", "RunParams::GetFixedFormatParams",datfilename);//number of atoms to move in a single move
	histbuffsize=ReadThisLine(file,1,HIST_BUFF_DEF,"histbuffsize", "RunParams::GetFixedFormatParams",datfilename);//size of the history buffer
	histstepratio=ReadThisLine(file,1,HIST_STEP_DEF,"histstepratio", "RunParams::GetFixedFormatParams",datfilename);//number of saves between each history buffering
	custmove=ReadThisLine(file,1,CUSTMOVE_DEF,"custmove", "RunParams::GetFixedFormatParams",datfilename);//indicator of custom move	
	if (ncommonneigh > 0 && custmove > 0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Custom move cannot be used with common neighbour constraint!" << endl;
		cout << "Exiting..." << endl;
		CleanExit();
	}
	if (ncommonneigh > 0 && nmoved > 1)
	{
		cout << "\nWARNING(" << ++warn << "): Common neighbour constraint can only be used with 1 moved atom!" << endl;
		cout<<"\tThe number of moved atoms will be set to 1!" << endl;
		nmoved = 1;
	}
	reload=ReadThisLine(file,1,RELOAD_DEF,"reload", "RunParams::GetFixedFormatParams",datfilename);//whether to load the histogram (and the coordination numbers,
					//if there is any) from files. (histogram:.hst, CoordNumbConst:.cnc, AvCoordConst:.acn)
	max_gridatom=ReadThisLine(file,1,MAX_GRIDATOM_DEF,"max_gridatom", "RunParams::GetFixedFormatParams",datfilename);//maximum number of atoms in a  grid cell
	
	//ratio of swaps among the moves
	swap_fraction=ReadThisLine(file,3,SWAP_FRACTION_DEF,"swap_fraction", "RunParams::GetFixedFormatParams");

	if (swap_fraction>1.0 || swap_fraction<0)
	{
		cout << "\nWARNING(" << ++warn << "): RunParams::GetFixedFormatParams "<<endl;
		cout<<"\tFraction of swaps related to the moves has to be between 0 (no swaps) and 1 (only swaps, no moves)!"<<endl;
		cout<<"\tSetting swap fraction to 0!"<<endl;
		swap_fraction=0;
	}

	if (swap_fraction>0 && (fnc>0 || ntypes==1))
	{
		cout << "\nWARNING(" << ++warn << "): RunParams::GetFixedFormatParams "<<endl;
		cout<<"\tSwaps cannot be used in case of flexible molecular/fixed neighbour constraint!"<<endl;
		cout<<"\tSetting swap fraction to 0!"<<endl;
		swap_fraction=0;
	}
	
	if (swap_fraction>0)//there will be swaps
		temp=ntypes*(ntypes-1)/2;//number of mixed partials
	else
		temp=0;
	//allocating the arrays to hold the type of atoms pairs allowed to be swapped
	SetArraysize(&swap_type1,temp,"swap_type1","RunParams::GetFixedFormatParams");
	SetArraysize(&swap_type2,temp,"swap_type2","RunParams::GetFixedFormatParams");
	
	if (swap_fraction>0)//there will be swaps
	{
		for (i=0; i<ntypes-1;i++)
		{
			for (j=i+1; j<ntypes;j++)
			{
				//reading the boolean indicator, whether this mixed partial is allowed to be involved in swap
				if (i==ntypes-2 && j==ntypes-1)//go to next line after reading
					temp=ReadThisLine(file,1,0,"swap_type", "RunParams::GetFixedFormatParams",datfilename);
				else//not to go to next line after reading
					temp=ReadThisLine(file,3,0,"swap_type", "RunParams::GetFixedFormatParams");//reading the boolean indicator, whether this mixed partial is allowed to be involved in swap

				if (temp)
				{
					//this partial is allowed to be involved in swap
					swap_type1[nswap_pairs]=i;//type of the first atom of the pair
					swap_type2[nswap_pairs]=j;//type of the second atom of the pair
					nswap_pairs++;
				}
			}
		}
	}
	else
		SkipLine(file,datfilename,1);

	if (swap_fraction>0 && nswap_pairs==0)
	{
		cout << "\nWARNING(" << ++warn << "): RunParams::GetFixedFormatParams "<<endl;
		cout<<"\tThere are no partials, where swap is alllowed! Cannot swap this way!"<<endl;
		cout<<"\tSetting swap fraction to 0!"<<endl;
		swap_fraction=0; 
	}
	//Decide whether this will be the last to read
	i=3;
#ifdef _LOCAL_INV
	i=1;
#endif
#ifdef _NO_PERIODIC
	i=1;
#endif
#ifdef _VIBR_AMP
	i=1;
#endif
	nthreads=ReadThisLine(file,i,1,"nthreads", "RunParams::GetFixedFormatParams");//total number of threads
	if (nthreads<1)
	{
		cout << "\nWARNING(" << ++warn << "): The number of threads has to be positive! Resetting the number of threads to 1"<<endl;
		nthreads=1;
	}
	
#ifdef _LOCAL_INV
	//ntotal_points will be updated by ChiSquared
	SetArraysize(&text, NAME_SIZE,"text","RunParams::GetFixedFormatParams");
	//number of local invariance intervals (for bin_based only one possible)
	nlocint=ReadThisLine(file,3,1,"nlocint", "RunParams::GetFixedFormatParams");
	if (nlocint<1)
	{
		cout << "\nWARNING(" << ++warn << "): There has to be at least 1 local invariance interval, setting it to 1!"<<endl;
		nlocint=1;
	}

	SetArraysize(&loc_inv_sigma,nlocint,"loc_inv_sigma","RunParams::GetFixedFormatParams");

	//local invariance will be calculated
	for (i=0;i<nlocint;i++)
	{
		mystrcpy(name, NAME_SIZE, "loc_inv_sigma[");
		IntToStr(&conv_numb,i + 1);
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		loc_inv_sigma[i]=ReadThisLine(file,3,1.0,name, "RunParams::GetFixedFormatParams");
		if (loc_inv_sigma[i]<0)
			ChiSquared::calc_sigma=1;
	}

	loc_chi2_mode=ReadThisLine(file,3,1,"loc_chi2_mode", "RunParams::GetFixedFormatParams");//calculation mode for local chi2, 0 or 1
	if (loc_chi2_mode<0 || loc_chi2_mode>1)
	{
		cout << "\nWARNING(" << ++warn << "): The calculation mode for the local chi2 can be 0 (bin based) or 1 (distance base). "<<endl;
		cout<<"\tIt will be set to the default "<<LOC_CHI2_MODE_DEF<<"!"<<endl;
		loc_chi2_mode= LOC_CHI2_MODE_DEF;
	}
	if (loc_chi2_mode==0)
	{//for the bin-based calculation a minimum and a maximum distance can be given 
		if (nlocint>1)
		{
			cout << "\nWARNING(" << ++warn << "): The number of local invariance intervals cannot be more than 1 for bin-based calculation!"<<endl;
			cout<<"\tChange the "<<datfilename<<" file, and try again! Exiting..."<<endl;
			CleanExit();
		}
		min_loc_r=ReadThisLine(file,3,0.0,"min_loc_r", "RunParams::GetFixedFormatParams");//minimum dist in reduced units for local inv. calc.
		max_loc_r=ReadThisLine(file,3,xmax[0],"max_loc_r", "RunParams::GetFixedFormatParams");//maximum dist in reduced units for local inv. calc.
		if (min_loc_r<0 || max_loc_r<0)
	
		{
			if (min_loc_r<0)
			{
				mystrcpy(text, NAME_SIZE, "minimum");
				if (max_loc_r<0)
					mystrcat(text, NAME_SIZE, " and maximum distance for the local invariance calculation");
			}
			else
			{
				if (max_loc_r<0)
					mystrcpy(text, NAME_SIZE, "maximum distance for the local invariance calculation");

			}
			
			cout << "\nWARNING(" << ++warn << "): The "<<text<<" in the "<<datfilename<<" file is negative!"<<endl;
			cout<<"\tIt has to be in the range of the normal histogram, between "<<binshift_def*rspacing[0]/boxedge;
			cout<<"\tand "<<xmax[0]<<" in reduced units! Set it right, and try again! Exiting..."<<endl;
			CleanExit();
		}
	}
	else //(loc_chi2_mode==1)
	{
		//distance-based calculation, the local histograms has to be calculated for all the possibale bins, as we do not know how much the local 
		//environments differ
		//percentage of the atoms will decide, how many atoms (between int(loc_at_ration[i]*natom_iytpe) and int(loc_at_ration[i+1]*natom_iytpe))
		//of a type will be involved in the distance based calculation
		max_loc_r=SQRT3;
		SetArraysize(&loc_at_ratio,nlocint+1,"loc_at_ratio","RunParams::GetFixedFormatParams");
	
		j=3;
		for (i=0;i<nlocint+1;i++)
		{
			mystrcpy(name, NAME_SIZE, "loc_at_ratio[");
			IntToStr(&conv_numb,i + 1);
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, "]");
			if (i==nlocint)
				j=1;
			loc_at_ratio[i]=ReadThisLine(file,j,0.0,name, "RunParams::GetFixedFormatParams");//ratio of the atoms for the intervalum borders in loc inv. calc. 
		}
		if (loc_at_ratio[0]<0 || loc_at_ratio[0]>=1.0)
		{
			cout << "\nWARNING(" << ++warn << "): The ratio of atoms to start the distance based calculation of the local invariance has to be greater or equl than 0,"<<endl;
			cout<<"\tsmaller than 1.0! It will be set to 0.0"<<endl;
			loc_at_ratio[0]=0.0;
		}
		for (i=1;i<nlocint+1;i++)
		{
			if (loc_at_ratio[i]<=0 || loc_at_ratio[i]>1.0 || loc_at_ratio[i]<loc_at_ratio[i-1])
			{
				cout << "\n*****ERROR*****"<<endl;
				cout<<"The ratio of atoms for the end of the "<<i<<". interval in the distance based calculation of "<<endl;
				cout<<"the local invariance ("<<loc_at_ratio[i]<<") has to be greater than 0 and <=1.0 "<<endl;
				cout<<"and greater than the previous ratio ("<<loc_at_ratio[i-1]<<")"<<endl;
				cout<<"Change the "<<datfilename<<" file and try again! Exiting..."<<endl;
				CleanExit();
			}
		}

	}
	
#endif



#ifdef _NO_PERIODIC

	if (fabs(xmax[0]-0.5)>TOLERANCE2)
	{
		cout << "\nWARNING(" << ++warn << "):xmax has to be 0.5 in case of non-periodic boundary conditions!" << endl;
		cout<<"\txmax will be reset to 0.5!"<<endl;
		xmax[0]=0.5;

	}
	R0=ReadThisLine(file,3,1.0,"R0", "RunParams::GetFixedFormatParams");//the radius of the spheric sample inside the box
	i=8;
#ifdef _VIBR_AMP
	i=6;
#endif
	recentre_flag=ReadThisLine(file,i,RECENTRE_FLAG_DEF,"recentre_flag", "RunParams::GetFixedFormatParams",datfilename);//the width of the local invariance histogram bins
	//if no periodic boundary conditions are used, then the sample have to be in the middle of the box not exceeding 
	//the quarter of the box edge from the centre in each coordinate axis direction

	//check, if fnc>0 was given
	if (RunParams::fnc>0)
	{
		cout << "\n***** ERROR *****"<<endl;
		cout<<"FNC constraint or flexible molecules cannot be used in case of non-periodic boundary conditions!"<<endl;
		cout<<"Set fnc option to 0, and try again! exiting"<<endl;
		CleanExit();
	}

	v_R0=4.0/3.0*PI*pow(R0,3);
	dtemp=SimpleCfg::ntotal/v_R0;
	if (fabs(dtemp-rho)>5e-7)
	{
		cout << "\nWARNING(" << ++warn << "): The radius of the sample (R0) does not correspond to the density,"<<endl;
		cout<<"\tthe R0 will be reset from "<<R0;
		R0=pow((SimpleCfg::ntotal/rho)*3.0/4.0/PI,1.0/3.0);
		cout<<" to "<<R0<<" A!"<<endl;
	}

	//We have to make sure, that the R0 is divisible by dr to make the volume correction easier
	dtemp=R0/rspacing_def;
	temp=(int)dtemp;
	if (fabs(dtemp-(double)temp-binshift_def)>TOLERANCE)
	{
		cout << "\nWARNING(" << ++warn << "): The radius of the sample is not divisible by the bin width taking " << endl;
		cout<<"\tinto account the binshift as well!"<<endl;
		cout<<"\tThe binwidth will be reset from "<<rspacing_def<<" to ";
		if (dtemp-temp-binshift_def>0.5)
		{
			rspacing_def=R0/(temp+binshift_def+1);//make it a bit smaller
			temp++;
		}
		else
			rspacing_def=R0/(temp+binshift_def);//make it a bit larger
		cout<<rspacing_def<<endl;
		for (i=0;i<ntot_datasets;i++)
			rspacing[i] = rspacing_def;
	}
	nbins[0]=temp;
	if (fabs(SimpleCfg::boxedge/2.0-R0)>1e-7)
	{
		cout << "\nWARNING(" << ++warn << "): The total boxlength has to be four times the sample radius, the half  "<<endl;
		cout<<"\tboxlength will be reset from "<<SimpleCfg::boxedge << " to ";
		SimpleCfg::boxedge=R0*2;
		boxedge=SimpleCfg::boxedge;
		cout<<SimpleCfg::boxedge<<endl;
		for (i=0;i<navcoord;i++)
		{
			//Reset minimum and maximum distances for the constraints squared in reduced units
			AvCoordConst::udminsq[i]=AvCoordConst::dmin[i]*AvCoordConst::dmin[i]/boxedge/boxedge;
			AvCoordConst::udmaxsq[i]=AvCoordConst::dmax[i]*AvCoordConst::dmax[i]/boxedge/boxedge;

			
	}
		for (i=0;i<nicoord;i++)
		{
			//Initialising the remainder of the static members
			//Minimum and maximum distances for the constraints squared in reduced units
			for (ineightype=0;ineightype<CoordNumbConst::n_neightype[i];ineightype++)
			{
				CoordNumbConst::udminsq[CoordNumbConst::cum_n_neightype[i]+ineightype]=	\
					CoordNumbConst::dmin[CoordNumbConst::cum_n_neightype[i]+ineightype]*CoordNumbConst::dmin[CoordNumbConst::cum_n_neightype[i]+ineightype]/boxedge/boxedge;
				CoordNumbConst::udmaxsq[CoordNumbConst::cum_n_neightype[i]+ineightype]=	\
					CoordNumbConst::dmax[CoordNumbConst::cum_n_neightype[i]+ineightype]*CoordNumbConst::dmax[CoordNumbConst::cum_n_neightype[i]+ineightype]/boxedge/boxedge;

			}
		
		}
	}
	//this has to be set, regardless if the boxedge was changed
	for (i=0;i<navcoord;i++)
	{
		AvCoordConst::dav_red[i]=(AvCoordConst::dmax[i]+AvCoordConst::dmin[i])/2/boxedge;//reduced maximum distance
		AvCoordConst::ncentral[i]= static_cast<int>(round(AvCoordConst::ncentral[i]*(pow(R0-(AvCoordConst::dav_red[i]*boxedge),3)/pow(R0,3))));//number of central atoms, which coordination atom is inside the sample
		if (AvCoordConst::ncentral[i] == 0)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The number of available central atoms for the " << i + 1 << ". average coordinaion constraint data set is zero" << endl;
			cout << "Check the constraint parameters in the " << datfilename << " file!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}
	}
	
#endif

#ifdef _VIBR_AMP
	char *name2;
	SetArraysize(&name2,30,"name2","RunParams::GetFixedFormatParams");
	j=3;
	for (i=0;i<ntypes;i++)
	{
		mystrcpy(name2,30,"T_corr_sigma_t[");
		IntToStr(&conv_numb,i + 1);
		mystrcat(name2,30,conv_numb);
		mystrcat(name2,30,"]");
		if (i==ntypes-1)
			j=1;
		T_corr_sigma_t[i]=ReadThisLine(file,j,0.01,name2, "RunParams::GetFixedFormatParams");//the sigma parameter of the Gauss distribution for 
			//thermal correction convolution, connected to the root mean square atomic vibrational ampiltude
		if (T_corr_sigma_t[i]*4<rspacing[0]/2)
		{
			cout << "\n*****ERROR*****" << endl;
			cout<<"The sigma parameter for thermal correction ("<<T_corr_sigma_t[i]<<") is too small compared to the bin size."<<endl;
			cout<<"The 4*T_corr_sigma_t accounting for 99.99% of the values are inside one histgram bin, so"<<endl;
			cout<<"no thermal correction will be performed!!! Either change the bin size or the T_sigma_corr"<<endl;
			cout<<"value! Exiting..."<<endl;
			CleanExit();
		}
	}
	delete [] name2;
	
#endif
	if (conv_numb != NULL)
		delete [] conv_numb;
	delete [] name;
	//Checking, whether loading was successful
	if (CheckReadFileState(file,"RunParams::GetFixedFormatParams",datfilename))
		cout<<"\nLoading the parameters from the fixed format "<<datfilename<<" file was successful!"<<endl;//load was successful
	else
		CleanExit();//loading was not successful

	//Reading keywords
	if(SkipLine(file))
	{
		bool endcyc=false;
		string line;
		getline(file,line);
		do
		{
			std::size_t eqpos=line.find("=");
			if(eqpos!=string::npos)
			{
				stringstream ss(line.substr(eqpos+1,line.length()-eqpos-1));
				int value;
				ss>>value;
				if(ss.fail())
				{
					cout << "\nWARNING(" << ++warn << "): Invalid integer value appeared in the following line:\n"<<line<<endl;
				}
				else
				{
					if((line.substr(0,eqpos)).find("OLD_OUT")!=string::npos) RunParams::old_out=(value!=0);
					else if((line.substr(0,eqpos)).find("SUM_PPCF")!=string::npos) RunParams::sum_ppcf=(value!=0);
					else cout << "\nWARNING(" << ++warn << "): The line does not contain any valid keywords such as OLD_OUT, SUM_PPCF in\n"<<line<<endl;
				}
				getline(file,line);
				endcyc=file.eof();
			}
			else endcyc=true;
		}while(!endcyc);
	}
	
	cout << "\nNOTE(" << ++note << "): The following setup will be used for some infrequently altered options [keyword:actual:default]:\n";
	cout<<"\t[OLD_OUT:"<<RunParams::old_out<<":"<<OLD_OUT_DEF<<"] Old output format "<<(RunParams::old_out ? "will" : "won't")<<" be written."<<endl;
	cout<<"\t[SUM_PPCF:"<<RunParams::sum_ppcf<<":"<<SUM_PPCF_DEF<<"] Summed and averaged ppcf-s "<<(RunParams::sum_ppcf ? "will" : "won't")<<" be written at each save!"<<endl;
	
	CoordNumbConst::nthreads=nthreads;
	AvCoordConst::nthreads=nthreads;

	
	//Checking, whether loading was successful
	if (file.rdstate() & file.failbit)
		cout<<"\nLoading the optional keywords at the end of the fixed format parameter file "<<datfilename<<" file was successful!"<<endl;//load was successful
	else
		file.clear(ios::goodbit);
	
}

//reading the potential related parameters from the *.dat file
void RunParams::GetPotParams(ifstream &file)
{
	int i, j;
	int npartials,nGRpartials;
	char *name,*conv_numb=NULL;
	ifstream potfile;

	SetArraysize(&name,NAME_SIZE,"name","RunParams::GetPotParams");

	npartials=ntypes*(ntypes+1)/2;
	//1:LJ 
	//10:tabulated potential
	
	switch (potential)
	{
		case 1:
		{
			mystrcpy(vdW_name, 15, "LJ");
			//shows, how many weigth parameter (sigma) is needed
			//negative value shows, that a percentage value after the last sigma for LJ will be read as well
			//0:  one weight paramter for vdW and one for Coulomb, kept for compatibility
			//1,-1:  npartials weight parameter for both vdW and for Coulomb
			//2,-2:  only one weight parameter for vDW, and the sigma (calculated based on this if it is negative, scalable) 
			//	  will be used for all the bonded and non-bonded potential contributions, to leave their original ration to each other
			//3,-3:  one weight paramter for vdW and one for Coulomb, same as 0, but had to be introduced, as negative value cannot used otherwise

			NB_weight_mode = ReadThisLine(file, 3, 0, "NB_weight_mode", "RunParams::GetPotParams");

			//combination rule for creating the parameters for the mixed partial, following GROMACS (except option 0)
			//0: sigma and epsilon for all the partials will be given
			//1: C(6) and C(N) is given for each type, for the  mixed pairs C_ij(M)=sqrt(C_i(M)*C_j(M)) M=6 or N, where N is LJ_rep_N
			//2: sigma and epsilon is given, sigma_ij=0.5*(sigma_i+sigma_j), epsilon_ij=sqrt(epsilon_i*epsilon_j)
			//3: sigma and epsilon is given, sigma_ij=sqrt(sigma_i*sigma_j), epsilon_ij=sqrt(epsilon_i*epsilon_j) (like in OPLS)
			vdW_comb_rule = ReadThisLine(file, 3, 0, "vdW_comb_rule", "RunParams::GetPotParams");

			if (vdW_comb_rule < 0 || vdW_comb_rule>3)
			{
				cout << "\nWARNING(" << ++warn << "): The combination rule to create the potential parameters for the mixed partials can be 0-3!" << endl;
				cout << "\t"<<vdW_comb_rule << " was given, setting it to default 0, meaining that trying to read the sigma and epsilon values for all the partials!" << endl;
				vdW_comb_rule = VDW_COMB_RULE_DEF;
			}

			vdW14_fudge = ReadThisLine(file, 3, FUDGE14_DEF, "vdW14_fudge", "RunParams::GetPotParams");
			Coulomb14_fudge = ReadThisLine(file, 3, FUDGE14_DEF, "Coulomb14_fudge", "RunParams::GetPotParams");
			lead_series_ind2 = ReadThisLine(file, 3, LEAD_SERIES_IND2_DEF, "lead_series_ind2", "RunParams::GetPotParams");
			//power of the repulsion term, in 6-N LJ 
			LJ_rep_N = ReadThisLine(file, 8, LJ_REP_N_DEF, "LJ_rep_N", "RunParams::GetPotParams");
			SkipLine(file, "RunParams::GetPotParams", 1);
			vdW_cutoff = ReadThisLine(file, 3, CUTOFF_DEF, "vdW_cutoff", "RunParams::GetPotParams");
			Coulomb_cutoff = ReadThisLine(file, 1, CUTOFF_DEF, "Coulomb_cutoff", "RunParams::GetPotParams");
			nGRtypes = ReadThisLine(file, 1, 1, "nGRtypes", "RunParams::GetPotParams");
			nGRpartials = nGRtypes * (nGRtypes + 1) / 2;
			SetArraysize(&vdW_pot1, nGRpartials, "vdW_pot1", "RunParams::GetPotParams");//if LJ, this is the sigma parameter of the potential, describing where the pot crosses the x-axis
			SetArraysize(&vdW_pot2, nGRpartials, "vdW_pot2", "RunParams::GetPotParams");//if LJ, epsilon is the depth of the potential valley

			if (potential == 1)//there can be the same reading before this for other LJ-similar potential
			{
				if (vdW_comb_rule == 0)//parameters for all the partials will be supplied
				{
					for (i = 0; i < nGRpartials; i++)
					{
						if (i < nGRpartials - 1)
							vdW_pot1[i] = ReadThisLine(file, 3, 1.0, "LJ_sigma", "RunParams::GetPotParams");
						else
							vdW_pot1[i] = ReadThisLine(file, 1, 1.0, "LJ_sigma", "RunParams::GetPotParams", datfilename);//go to next line
					}
					for (i = 0; i < nGRpartials; i++)
					{
						if (i < nGRpartials - 1)
							vdW_pot2[i] = ReadThisLine(file, 3, 1.0, "LJ_epsilon", "RunParams::GetPotParams");
						else
							vdW_pot2[i] = ReadThisLine(file, 1, 1.0, "LJ_epsilon", "RunParams::GetPotParams", datfilename);//go to next line
					}
				}
				else
				{
					//parameters for the different types will be given, and the parameters for the mixed partials will be calculated according to the 
					//combination rule
					for (i = 0; i < nGRtypes; i++)
					{
						for (j = i; j < nGRtypes; j++)
						{
							if (i == j)//only read into the ii partials
							{
								if (i*nGRtypes + j - (i*(i + 1) / 2) < nGRpartials - 1)
									vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)] = ReadThisLine(file, 3, 1.0, "vdW_pot1", "RunParams::GetPotParams");
								else
									vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)] = ReadThisLine(file, 1, 1.0, "vdW_pot1", "RunParams::GetPotParams", datfilename);//go to next line
							}
						}
					}
					for (i = 0; i < nGRtypes; i++)
					{
						for (j = i; j < nGRtypes; j++)
						{
							if (i == j)//only read into the ii partials
							{
								if (i*nGRtypes + j - (i*(i + 1) / 2) < nGRpartials - 1)
									vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)] = ReadThisLine(file, 3, 1.0, "vdW_pot2", "RunParams::GetPotParams");
								else
									vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)] = ReadThisLine(file, 1, 1.0, "vdW_pot2", "RunParams::GetPotParams", datfilename);//go to next line
							}
						}
					}

					//Now create the parameters for the mixed partials according to the combination rule
					for (i = 0; i < nGRtypes; i++)
					{
						for (j = i; j < nGRtypes; j++)
						{
							if (i != j)//mixed partials
							{
								switch (vdW_comb_rule)
								{
								case 1://C6 and CN is given for LJ
								case 3://sigma and epsilon is given, geometric average
								{
									
									vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)] = sqrt(vdW_pot1[i*nGRtypes + i - (i*(i + 1) / 2)] * vdW_pot1[j*nGRtypes + j - (j*(j + 1) / 2)]);
									vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)] = sqrt(vdW_pot2[i*nGRtypes + i - (i*(i + 1) / 2)] * vdW_pot2[j*nGRtypes + j - (j*(j + 1) / 2)]);
									break;
								}
								case 2://sigma and epsilon is given, arithmetic average
								{
									
									vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)] = 0.5*(vdW_pot1[i*nGRtypes + i - (i*(i + 1) / 2)] + vdW_pot1[j*nGRtypes + j - (j*(j + 1) / 2)]);
									vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)] = sqrt(vdW_pot2[i*nGRtypes + i - (i*(i + 1) / 2)] * vdW_pot2[j*nGRtypes + j - (j*(j + 1) / 2)]);
								}
								}
							}
						}//end of j cycle
					}//end of i cycle
				}//end of vdW_comb_rule!=0
			}//end of LJ	

			//The weighting parameters for the NB interactions
			SetArraysize(&vdW_weight, npartials, "vdW_weight", "RunParams::GetPotParams");
			SetArraysize(&Coulomb_weight, npartials, "Coulomb_weight", "RunParams::GetPotParams");

			if (NB_weight_mode == 0 || fabs((double)NB_weight_mode) == 2 || fabs((double)NB_weight_mode) == 3)//only one parameter is read
			{
				if (NB_weight_mode < 0)
					i = 3;
				else
					i = 1;
				vdW_weight[0] = ReadThisLine(file, i, 1.0, "vdW_weight", "RunParams::GetPotParams", datfilename);
				CheckSigma(vdW_weight, "vdw_weight[1] parameter");//check, if it not zero
				if (vdW_weight[0] < 0)
					ChiSquared::calc_sigma = 1;
				if (NB_weight_mode < 0)//read pot_chi2_low_lim_fraction
					pot_chi2_low_lim_fraction = ReadThisLine(file, 1, POT_CHI_LOW_LIM_FRACTION_DEF, "pot_chi2_low_lim_fraction", "RunParams::GetPotParams", datfilename);
			}
			else
			{
				if (NB_weight_mode < 0)
					j = 3;
				else
					j = 1;
				for (i = 0; i < npartials; i++)
				{
					if (i < npartials - 1)
						vdW_weight[i] = ReadThisLine(file, 3, 1.0, "vdW_weight", "RunParams::GetPotParams");
					else
					{

						vdW_weight[i] = ReadThisLine(file, j, 1.0, "vdW_weight", "RunParams::GetPotParams", datfilename);
					}
					mystrcpy(name, NAME_SIZE, "vdW_weight[");
					IntToStr(&conv_numb, i + 1);
					mystrcat(name, NAME_SIZE, conv_numb);
					mystrcat(name, NAME_SIZE, "] parameter never");
					CheckSigma(vdW_weight + i, name);//check, if it not zero
					if (vdW_weight[i] < 0)
						ChiSquared::calc_sigma = 1;
				}
				if (NB_weight_mode < 0)//read pot_chi2_low_lim_fraction
					pot_chi2_low_lim_fraction = ReadThisLine(file, 1, POT_CHI_LOW_LIM_FRACTION_DEF, "pot_chi2_low_lim_fraction", "RunParams::GetPotParams", datfilename);


			}

			if (NB_weight_mode == 0)//only one parameter is read
			{
				Coulomb_weight[0] = ReadThisLine(file, 1, 1.0, "Coulomb_weight", "RunParams::GetPotParams", datfilename);
				if (Coulomb_weight[0] < 0 || fabs(Coulomb_weight[0] - 0.0) < TOLERANCE)//this can be zero as well
					ChiSquared::calc_sigma = 1;
				if (fabs(Coulomb_weight[0] - 0.0) < TOLERANCE)
				{
					if (vdW_weight[0] > 0)
						Coulomb_weight[0] = vdW_weight[0];//fixed weight, can be set now
				}

			}
			else
			{
				if (fabs(NB_weight_mode) == 1)
				{
					for (i = 0; i < npartials; i++)
					{
						if (i < npartials - 1)
							Coulomb_weight[i] = ReadThisLine(file, 3, 1.0, "Coulomb_weight", "RunParams::GetPotParams");
						else
							Coulomb_weight[i] = ReadThisLine(file, 1, 1.0, "Coulomb_weight", "RunParams::GetPotParams", datfilename);

						mystrcpy(name, NAME_SIZE, "Coulomb_weight[");
						IntToStr(&conv_numb, i + 1);
						mystrcat(name, NAME_SIZE, conv_numb);
						mystrcat(name, NAME_SIZE, "] parameter together with the non-bonded interaction weight mode=1");
						CheckSigma(Coulomb_weight + i, name);//check, if it not zero
						if (Coulomb_weight[i] < 0)
							ChiSquared::calc_sigma = 1;
					}
				}
				else
				{
					SkipLine(file);//weight mode is presumably 2, the vdW sigma will be used, so skip the line
					if (vdW_weight[0] > 0)
						Coulomb_weight[0] = vdW_weight[0];//fixed weight, can be set now
					else
						Coulomb_weight[0] = -1e6;//scalable, will be set later
				}

			}
			if (fabs(NB_weight_mode) != 1)
			{
				//fill the remaining array elements with the same value, as for the sake of simplicity the potentilas will be always split into partials
				for (i = 1; i < npartials; i++)
				{
					vdW_weight[i] = vdW_weight[0];
					Coulomb_weight[i] = Coulomb_weight[0];
				}
			}
			if (fabs(NB_weight_mode) == 2 && lead_series_ind2 != 1 && vdW_weight[0] < 0)
			{
				cout << "\nWARNING(" << ++warn << "): In case of NB weight mode=2 (all the potential related contributions " << endl;
				cout<<"\twill be weighted with the sigma of the van der Waals interaction)" << endl;
				cout << "\tand scalable sigma for this the leading potential series index has to be 1!" << endl;
				cout << "\tIt is reset accordingly." << endl;
				lead_series_ind2 = 1;
			}
			break;

		}
		case 10:  //tabulated potential read from file(s), there can be as many as the number of RMC partials, one file for a given partial
		{
			mystrcpy(vdW_name, 15, "tabulated");
			lead_series_ind2 = ReadThisLine(file, 6, 1, "lead_series_ind2", "RunParams::GetPotParams");
			vdW_cutoff = ReadThisLine(file, 1, 1.0, "vdW_cutoff", "RunParams::GetPotParams");//Angstrom
			temperature= ReadThisLine(file, 1, 1.0, "temperature", "RunParams::GetPotParams");
			nused_potpartials = ReadThisLine(file, 1, 1, "nused_potpartials", "RunParams::GetPotParams");
			FNC_POT::nused_potpartials = nused_potpartials;
			if (nused_potpartials <= 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Number of tabulated partial potential data is " << nused_potpartials << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}
		
			SetArraysize(&vdW_weight, nused_potpartials, "vdW_weight", "RunParams::GetPotParams");
			FNC_POT::ReadTabPot(file);
			
			NB_weight_mode = 1;//different sigmas for the different partials
			break;
		}
		
		default:
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Presently potential parameter can only be 1 (LJ) and 10 (tabulated potential)! " << endl;
			cout << "Cannor run this way, exiting..." << endl;
			CleanExit();
			

		}
	}//end of switch potential
	
	
	if (conv_numb != NULL)
		delete [] conv_numb;
	//ntotal_points2 will be set by ChiSquared::SetChiSquaredParams
}
	

//checking the leading series index, as it could not be done at reading, as the bonded interactions were not available then
//here lead_series_ind and lead_series_ind2 are starting with 1 at the beginning, but set it to satrt with 0 at the end!
void RunParams::CheckLeadSerInd()
{
	int tot;
	
	tot=ngr+nsq+nfq+nfg+nek+ncosdistr+navcoord+CoordNumbConst::tot_subconst;
#ifdef _ADVANCED_GEOM_CONST
	tot += CommonNeighConst::nconstraints + SecondNeighConst::nconstraints + BondValenceSumConst::nconstraints;
#endif
#ifdef _LOCAL_INV
	tot+=nlocint;
#endif


	int reset=0;
	int index,i;
	if ((potential==0) && (fnc==4))
		tot+=ChiSquared::nbonds+ChiSquared::nangles+ChiSquared::nperdihs+ChiSquared::nharmdihs+ChiSquared::nRBdihs;
	if (offset_lead2 > 0)
		OffsetPotLead();
	if (lead_series_ind<1)
	{
		cout << "\nWARNING(" << ++warn << "): The lead series index cannot be zero. Setting it to 1"<<endl;
		lead_series_ind=1;
	}
	if (lead_series_ind>tot)
	{
		cout << "\nWARNING(" << ++warn << "): The leading series index specified in the *.dat file is ("<<lead_series_ind<<" is greater,"<<endl;
		cout<<"\tthan the total number of series and constraints, which are "<<endl;
		cout<<"\tngr=                 "<<ngr<<endl;
		cout<<"\tnsq=                 "<<nsq<<endl;
		cout<<"\tnfq=                 "<<nfq<<endl;
		cout<<"\tnfg=                 "<<nfg<< endl;
		cout<<"\tnexafs=              "<<nek<<endl;
		cout<<"\tncosine=             "<<ncosdistr<<endl;
		cout<<"\tncoord=              "<<CoordNumbConst::tot_subconst<<endl;
		cout<<"\tnavcoord=            "<<navcoord<<endl;
#ifdef _ADVANCED_GEOM_CONST
		cout << "\tncommonneighbour=  " << ncommonneigh << endl;
		cout << "\tnsecondneighbour=  " << nsecondneigh << endl;
		cout << "\tnbvs=			  " << BondValenceSumConst::nconstraints << endl;
#endif
		
#ifdef _LOCAL_INV
		cout<<"\tLocal invariance   "<<nlocint<<endl;
#endif
		
		if ((potential==0)&&(fnc==4))
		{
			cout<<"\tbonds:             "<<ChiSquared::nbonds<<endl;
			cout<<"\tangles:            "<<ChiSquared::nangles<<endl;
			cout<<"\tper. dihedrals:    "<<ChiSquared::nperdihs<<endl;
			cout<<"\tharm. dihs:        "<<ChiSquared::nharmdihs<<endl;
			cout<<"\tRB. dihs:          "<<ChiSquared::nRBdihs<<endl;
		}
		cout<<"\tTotal=    "<<tot<<endl;
		cout<<"\tThe index will be set to 1!"<<endl;
		lead_series_ind=1;
	}
	index = 0;
#ifdef _LOCAL_INV
	index = nlocint;
#endif
	int *cum_series;
	//number of data sets before this
	SetArraysize(&cum_series, 20, "cum_series", "RunParams::CheckLeadSerInd");
	cum_series[0] = 0;
	cum_series[1] = cum_series[0]+ngr;
	cum_series[2] = cum_series[1] +nsq;
	cum_series[3] = cum_series[2] + nfq;
	cum_series[4] = cum_series[3] + nfg;
	cum_series[5] = cum_series[4] + nek;
	cum_series[6] = cum_series[5] + ncosdistr;
	cum_series[7] = cum_series[6] + CoordNumbConst::tot_subconst;
	cum_series[8] = cum_series[7] + navcoord;
	cum_series[9] = cum_series[8] + ncommonneigh;
	cum_series[10] = cum_series[9] + nsecondneigh;
	cum_series[11] = cum_series[10] + nbvs;
	cum_series[12] = cum_series[11] + index;
	cum_series[13] = cum_series[12] + ChiSquared::nbonds;
	cum_series[14] = cum_series[13] + ChiSquared::nangles;
	cum_series[15] = cum_series[14] + ChiSquared::nperdihs;
	cum_series[16] = cum_series[15] + ChiSquared::nharmdihs;
	cum_series[17] = cum_series[16] + ChiSquared::nRBdihs;
		

	if (lead_series_ind <= cum_series[1]) lead_series_name = std::to_string(lead_series_ind) + ". series, which is the "+ std::to_string(lead_series_ind-cum_series[0])+". g(r)  data set";
	else if (lead_series_ind <= cum_series[2]) lead_series_name = std::to_string(lead_series_ind) + ". series, which is the " + std::to_string(lead_series_ind - cum_series[1]) + ". neutron data set";
	else if (lead_series_ind <= cum_series[3]) lead_series_name = std::to_string(lead_series_ind) + ". series, which is the " + std::to_string(lead_series_ind - cum_series[2]) + ". X-ray data set";
	else if (lead_series_ind <= cum_series[4]) lead_series_name = std::to_string(lead_series_ind) + ". series, which is the " + std::to_string(lead_series_ind - cum_series[3]) + ". electron diffraction data set";
	else if (lead_series_ind <= cum_series[5]) lead_series_name = std::to_string(lead_series_ind) + ". series, which is the " + std::to_string(lead_series_ind - cum_series[4]) + ". EXAFS data set";
	else if (lead_series_ind <= cum_series[6]) lead_series_name = std::to_string(lead_series_ind) + ". constraint, which is the " + std::to_string(lead_series_ind - cum_series[5]) + ". cosine distribution of bond angles";
	else if (lead_series_ind <= cum_series[7]) lead_series_name = std::to_string(lead_series_ind) + ". constraint, which is the " + std::to_string(lead_series_ind - cum_series[6]) + ". coordination constraint";
	else if (lead_series_ind <= cum_series[8]) lead_series_name = std::to_string(lead_series_ind) + ". constraint, which is the " + std::to_string(lead_series_ind - cum_series[7]) +". average coordination constraint";
#ifdef  _ADVANCED_GEOM_CONST	
	else if (lead_series_ind <= cum_series[9]) lead_series_name = std::to_string(lead_series_ind) + ". constraint, which is the " + std::to_string(lead_series_ind - cum_series[8]) +". common neighbour constraint";
	else if (lead_series_ind <= cum_series[10]) lead_series_name = std::to_string(lead_series_ind) + ". constraint, which is the " + std::to_string(lead_series_ind - cum_series[9]) +". second neighbour constraint";
	else if (lead_series_ind <= cum_series[11]) lead_series_name = std::to_string(lead_series_ind) + ". constraint, which is the " + std::to_string(lead_series_ind - cum_series[10]) + ". bond valence sum constraint";
#endif
#ifdef _LOCAL_INV
	else if (lead_series_ind <= cum_series[12]) lead_series_name = std::to_string(lead_series_ind) +", which is the " + std::to_string(lead_series_ind - cum_series[11]) + ". local invariance interval";
#endif
	else
	{
		if ((potential == 0) && (fnc == 4))
		{
			if (lead_series_ind <= cum_series[13]) lead_series_name = std::to_string(lead_series_ind) + ", which is the " + std::to_string(lead_series_ind - cum_series[12]) + ". bond interaction";
			else if (lead_series_ind <= cum_series[14]) lead_series_name = std::to_string(lead_series_ind) + ", which is the " + std::to_string(lead_series_ind - cum_series[13]) + ". angle interaction";
			else if (lead_series_ind <= cum_series[15]) lead_series_name = std::to_string(lead_series_ind) + ", which is the " + std::to_string(lead_series_ind - cum_series[14]) + ". periodic dihedral interaction";
			else if (lead_series_ind <= cum_series[16]) lead_series_name = std::to_string(lead_series_ind) + ", which is the " + std::to_string(lead_series_ind - cum_series[15]) + ". harmonic dihedral interaction";
			else if (lead_series_ind <= cum_series[17]) lead_series_name = std::to_string(lead_series_ind) + ", which is the " + std::to_string(lead_series_ind - cum_series[16]) + ". RB dihedral interaction";

		}
	}
	delete[] cum_series;
	lead_series_ind--;
	
	if ((potential>0 && lead_series_ind2<1) || potential==0)
		return;//there is no scaling, no checking is necessary
	if (potential == 10)
	{
		if (lead_series_ind2 > ChiSquared::nNB)
		{
			cout << "\nWARNING(" << ++warn << "): The leading potential series index specified in the *.dat file is (" << lead_series_ind2 << endl;
			cout << "\tis greater than the total number of used tabulated partial potential interactions which is " << ChiSquared::nNB <<endl;
			cout << "\tThe index will be set to 1, to the (first partial of the) tabulated potential interaction!" << endl;
			lead_series_ind2 = 1;
			
		}

		lead_series_name2 = std::to_string(lead_series_ind2) + ". tabulated potential";
		lead_series_ind2--;
		return;
	}
	if (potential==20)
	{
		if (lead_series_ind2 > 1)//this cannot happen,but anyway
		{
			cout << "\nWARNING(" << ++warn << "): The leading potential series index specified in the *.dat file is (" << lead_series_ind2 << endl;
			cout << "\tis greater than 1, and in case of ANN potential ther can be only this potential term!" <<endl;
			cout << "\tThe index will be set to 1!" << endl;
			lead_series_ind2 = 1;
			
		}

		lead_series_name2 = "ANN potential";
		lead_series_ind2--;
		return;
	}
	index=2*ChiSquared::nNB;
	if (potential==1 && Topology::npair_types>0)
	{
		index+=2*ChiSquared::nNB;
		if (lead_series_ind2>2*ChiSquared::nNB && lead_series_ind2<4*ChiSquared::nNB+1)//the lead series index is on a 1-4 interaction, it cannot be, there is no separate sigma for it
		{
			cout << "\nWARNING(" << ++warn << "): It makes no sense to put the leading potential-based interaction index on an 1-4 interaction," << endl;
			cout<<"\tit is reset to the corresponding";
			if (lead_series_ind2<3*ChiSquared::nNB+1)
				cout<<"wdW partial"<<endl;
				
			else
				cout<<"Coulomb partial"<<endl;
			lead_series_ind2-=2*ChiSquared::nNB;
			cout<<"\tSo the lead series index for the potential-based interactions is "<<lead_series_ind2<<endl;
			reset=1;
		
		}
		if (abs(NB_weight_mode)==2 && fabs(Coulomb_weight[0]-0.0)<TOLERANCE && lead_series_ind2==2)//Coulomb would be the leading pot series
		{
			cout << "\nWARNING(" << ++warn << "): If the Coulomb potential's weight is zero " << endl;
			cout<<"\t(it will be set to the same as the vdW potential's weight)"<<endl;
			cout<<"\tit cannot be the leading potential series, the vdW potantial will be the leading potential series instead!"<<endl;
			lead_series_ind2=1;	//vdW will be
		}
	}
	//have to check, whether the leading potential series index is not on a bonded interaction with zero sigma (non-bonded was already checked) 
	index++;//the first bond interaction
	if (potential == 1)
	{
		for (i = 0; i < ChiSquared::nbonds; i++)
		{
			if (fabs(Topology::bond_sigma[i] - 0.0) < TOLERANCE && index == lead_series_ind2)
				NoLeadZeroSigma(i, "bond");
			index++;
		}
		for (i = 0; i < ChiSquared::nangles; i++)
		{
			if (fabs(Topology::angle_sigma[i] - 0.0) < TOLERANCE && index == lead_series_ind2)
				NoLeadZeroSigma(i, "angle");
			index++;
		}
		for (i = 0; i < ChiSquared::nperdihs; i++)
		{
			if (fabs(Topology::perdihedral_sigma[i] - 0.0) < TOLERANCE && index == lead_series_ind2)
				NoLeadZeroSigma(i, "periodic dihedral");
			index++;
		}
		for (i = 0; i < ChiSquared::nharmdihs; i++)
		{
			if (fabs(Topology::harmdihedral_sigma[i] - 0.0) < TOLERANCE && index == lead_series_ind2)
				NoLeadZeroSigma(i, "harmonic dihedral");
			index++;
		}
		for (i = 0; i < ChiSquared::nRBdihs; i++)
		{
			if (fabs(Topology::RBdihedral_sigma[i] - 0.0) < TOLERANCE && index == lead_series_ind2)
				NoLeadZeroSigma(i, "RB dihedral");
			index++;
		}


		tot = 2 * ChiSquared::nNB + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs + ChiSquared::nRBdihs;
		if (Topology::npair_types > 0)
			tot += 2 * ChiSquared::nNB;
		if (!reset && potential > 0 && lead_series_ind2 > tot)
		{
			cout << "\nWARNING(" << ++warn << "): The leading potential series index specified in the *.dat file is (" << lead_series_ind2 << " is greater," << endl;
			cout << "\tthan the total number of separately displayed potential-based interactions which are " << endl;
			cout << "\tvDW:         " << ChiSquared::nNB << endl;
			cout << "\tCoulomb:     " << ChiSquared::nNB << endl;
			cout << "\t1-4 vDW:     " << ChiSquared::nNB << endl;
			cout << "\t1-4 Coulomb: " << ChiSquared::nNB << endl;
			cout << "\tbonds:           " << ChiSquared::nbonds << endl;
			cout << "\tangles:          " << ChiSquared::nangles << endl;
			cout << "\tper. dihedrals:  " << ChiSquared::nperdihs << endl;
			cout << "\tharm. dihs:      " << ChiSquared::nharmdihs << endl;
			cout << "\tRB. dihs:        " << ChiSquared::nRBdihs << endl;
			cout << "\tTotal:           " << tot << endl;
			cout << "\tThe index will be set to 1, to the (first partial of the)  vdW interaction!" << endl;
			lead_series_ind2 = 1;
		}
		index = 2;
		if (Topology::npair_types > 0)
			index = 4;
		if (lead_series_ind2 <= ChiSquared::nNB) lead_series_name2 = std::to_string(lead_series_ind2) + ". vdW interaction";
		else if (lead_series_ind2 <= 2 * ChiSquared::nNB) lead_series_name2 = std::to_string(lead_series_ind2) + ", which is the " + std::to_string(lead_series_ind2 - ChiSquared::nNB) + ". Coulomb interaction";
		else if (lead_series_ind2 <= index * ChiSquared::nNB + ChiSquared::nbonds) lead_series_name2 = std::to_string(lead_series_ind2) + ", which is the " + std::to_string(lead_series_ind2 - index * ChiSquared::nNB) + ". bond interaction";
		else if (lead_series_ind2 <= index * ChiSquared::nNB + ChiSquared::nbonds + ChiSquared::nangles) lead_series_name2 = std::to_string(lead_series_ind2) + ", which is the " + std::to_string(lead_series_ind2 - index * ChiSquared::nNB - ChiSquared::nbonds) + ". angle interaction";
		else if (lead_series_ind2 <= index * ChiSquared::nNB + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs) lead_series_name2 = std::to_string(lead_series_ind2) + ", which is the " + std::to_string(lead_series_ind2 - index * ChiSquared::nNB - ChiSquared::nbonds - ChiSquared::nangles) + ". periodic dihedral interaction";
		else if (lead_series_ind2 <= index * ChiSquared::nNB + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs) lead_series_name2 = std::to_string(lead_series_ind2) + ", which is the " + std::to_string(lead_series_ind2 - index * ChiSquared::nNB - ChiSquared::nbonds - ChiSquared::nangles - ChiSquared::nperdihs) + ". harmonic dihedral interaction";
		else if (lead_series_ind2 <= index * ChiSquared::nNB + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs + ChiSquared::nRBdihs) lead_series_name2 = std::to_string(lead_series_ind2) + ", which is the " + std::to_string(lead_series_ind2 - index * ChiSquared::nNB - ChiSquared::nbonds - ChiSquared::nangles - ChiSquared::nperdihs - ChiSquared::nharmdihs) + ". RB dihedral interaction";
	}
	lead_series_ind2--;
}

//check, whether it is not zero, and stop, if it is
void RunParams::CheckSigma(double *sigma, const char *name)
{
	if (fabs(*sigma-0.0)<TOLERANCE)
	{
		cout<<"\n*****ERROR*****"<<endl;
		cout<<"The value of the "<<name<<" cannot be zero!"<<endl;
		cout<<"Check the parameter file, and try again!)"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();

	}
};
//----writing the status of the program at regular intervals-----
#ifdef _ADVANCED_GEOM_CONST
	void RunParams::PrintStatus(ExptsData &edata, ChiSquared &chi, CoordNumbConst &icc, AvCoordConst &avcc, CommonNeighConst &conc, SecondNeighConst &snc, BondValenceSumConst &bvs, SimpleCfg &config, int status)
#else
	void RunParams::PrintStatus(ExptsData &edata, ChiSquared &chi, CoordNumbConst &icc,AvCoordConst &avcc, SimpleCfg &config, int status)
#endif
{
	int i,j,index,isubconst;
	int precision_s,precision_l;
	double genstep=(double)n_generated-last_gen;//number of moves generated since the last display
	double triedstep= (double)n_tried-last_tried;//number of moves tried since the last display
	double accstep = (double)n_accepted - last_accepted;//number of moves accepted since the last display
	double *px,*py,*a,*b,*c,*d,*e,sum,background;//for Rw calculation
	cout.precision(2);
	
	cout << "\nTOTAL NUMBER OF MOVES:" << endl;
	cout << J20 << "GENERATED" << J13 << "TRIED" << J13 << "ACCEPTED";
	if (RunParams::swap_fraction > 0)
		cout << J20 << "SWAP GENERATED" << J20 << "SWAP ACCEPTED";
	if (RunParams::potential > 0)
		cout << J20 << "POTENTIAL ACCEPETED";
	cout << endl;
	cout << J20 << n_generated << J13 << n_tried << J13 << n_accepted;
	if (RunParams::swap_fraction > 0)
		cout << J20 <<n_swapgen << J20 << n_swapacc;
	if (RunParams::potential > 0)
		cout << J20 << n_potaccepted;
	cout << endl;

	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << "I(Q) RELATED:  ";
	if (ExptsData::is_IQ == 2)
		cout << J23 << "COVERGED NON-LIN REG.";
	if (ExptsData::is_IQmucorr)
		cout << J20 << "MU CORR. GENERATED" << J20 << "MU CORR. ACCEPTED";
	if (ExptsData::is_AXS > 0)
		cout << J20 << "f'-SHIFT GENERATED" << J20 << "f'-SHIFT ACCEPTED";
	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << endl;

	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << "               ";
	if (ExptsData::is_IQ == 2)
		cout << J23 << n_goodnonlinreg;
	if (ExptsData::is_IQmucorr)
		cout << J20 << n_mucorrgen << J20 << n_mucorracc;
	if (ExptsData::is_AXS > 0)
		cout << J20 << n_fprimeshiftgen << J20 << n_fprimeshiftacc;
	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << endl;

	if (ExptsData::is_E0shift)
	{
		cout << "EXAfS RELATED: " << J23 << "E0-SHIFT GENERATED" << J20 << "E0-SHIFT ACCEPTED" << endl;
		cout << "               " << J23 << n_E0shiftgen << J20 << n_E0shiftacc << endl;
	}

	if (triedstep>0 && genstep>0)
	{
		cout<<"AVERAGE RATIOS:"<<J23<<"TRIED/GENERATED"<<J20<<"ACCEPTED/TRIED"<<"  OVER LAST "<<(int)genstep<<" MOVES:"<<endl;
		if (genstep>0 && triedstep>0)
			cout<<"               "<<J23<<triedstep/genstep<<J20<<accstep/triedstep<<endl;
	}
	
	precision_s=2;
	precision_l=6;
	cout.precision(precision_l);
	cout.setf(ios::scientific, ios::floatfield);

	if(ntotal_points>0)
	{
		cout<<"Global chi-squared value (per data point):\t\t\t\t\t\t"<<J10<<chi.total/ntotal_points<<endl;
	}
	if (potential > 0)
	{
		switch (potential)
		{
			case 1:
				cout << "Global chi-squared value for the bonded and non-bonded interactions (per data point):\t" << J10 << chi.total2 / ntotal_points2 << endl;
				break;
			case 10:
				cout << "Global chi-squared value for the tabulated potential interactions (per data point):\t" << J10 << chi.total2 / ntotal_points2 << endl;
				break;
			case 20:
				cout<<"Global chi-squared value for the ANN potential interaction:\t\t\t\t"<<J10 <<chi.total2 / ntotal_points2 << endl;
				break;
		}
	}

	cout<<"Details per data sets "<<endl;
	a=chi.a;//multiplicative factor
	b=chi.b;//constant
	c=chi.c;//linear
	d=chi.d;//quadratic
	e=chi.e;//cubic

	index=0;
	for(i=0;i<ngr;i++)
	{	
		py = edata.gr_gfinder[i];
		px = edata.gr_rfinder[i];
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout<<"g(r) set #"<<i+1<<"/"<<ngr;
		cout.setf(ios::scientific, ios::floatfield);
		//the experimental data has to be multiplied with this renormalization factor to fit the calculated data
		cout<<"\tRenorm. a="<<*a;
		if (ExptsData::use_cubic[i])
		{
			cout<<"\tConstant b="<<*b<<endl;
			cout<<"\tLinear c="<<*c;
			cout<<"\tQuadratic d="<<*d<<endl;
			cout<<"\tCubic e="<<*e;
		}
		cout<<endl;
		cout.precision(precision_l);
		if (ExptsData::gruseR[i])
		{
			//cout << "\tRw/(sigma*npoints)=\t\t\t\t\t" << J10 << *(chi.chicomp + i) / (*(ExptsData::grused + i)) << endl;
			if (fabs(chi.Rw[index]+1000)>TOLERANCE)
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			else
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << "cannot be calculted, denominator zero" << endl;
		}
		else
		{
			cout << "\tChi-square/npoints=\t\t\t\t\t\t\t\t" << J10 << *(chi.chicomp + i) / (*(ExptsData::grused + i)) << endl;

			//caclculating the Rw
			sum = 0.0;
			if (ExptsData::use_cubic[i])
			{
				for (j = 0; j < *(ExptsData::grused + index); j++)
				{
					//background consists of the constant, linear, quadratic and cubic part 
					background = (*b + *c * (*px) + *d * (*px) * (*px) + *e * (*px) * (*px) * (*px));
					sum += pow(*a * *py + background, 2);
					px++;
					py++;
				}
			}
			else
			{
				for (j = 0; j < *(ExptsData::grused + index); j++)
				{
					sum += *py * *py;
					py++;
				}
			}
			sum *= pow(*a, 2);
			if (sum > TOLERANCE)
			{
				chi.Rw[index] = sqrt(*(chi.chicomp + index)) * ExptsData::grsigma[index] / sqrt(sum) * 100;
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			}
		}
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\tBin size for this dataset: " << rspacing[index] << " A" << endl;
			cout << "\tBin shift for this dataset: " << binshift[index] << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
		a++;//go to next
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next 
		index++;
	}
	
	for(i=0;i<nsq;i++)
	{
		px = edata.sq_qfinder[i];
		py = edata.sq_sfinder[i];
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout<<"S(Q) set #"<<i+1<<"/"<<nsq;
		cout.setf(ios::scientific, ios::floatfield);
		//the experimental data-background has to be multiplied with this renormalization factor to fit the calculated data
		cout<<"\tRenorm. a="<<*a;
		//background consists of the constant, linear, quadratic and if use_cubic was given then a cubic part:
		//this constant has to be added to the expt data before renormalization for comparison with RMC data		
		cout<<"\tConstant b="<<*b<<endl;
		//this linear coefficient*Q has to be added to the expt data before renormalization for comparison with RMC data
		cout<<"\tLinear c="<<*c;
		//this quadratic coefficient*Q*Q has to be added to the expt data before renormalization for comparison with RMC data
		cout<<"\tQuadratic d="<<*d<<endl;
		//this cubic coefficient*Q*Q*Q has to be added to the expt data before renormalization for comparison with RMC data
		if (ExptsData::use_cubic[index])
		{
			//this cubic coefficient*Q*Q*Q has to be added to the expt data before renormalization for comparison with RMC data
			cout<<"\tCubic e="<<*e<<endl;
		}
		cout.precision(precision_l);
		if (ExptsData::squseR[i])
		{
			//cout << "\tRw/(sigma*npoints)=\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::sqused + i)) << endl;
			if (fabs(chi.Rw[index] + 1000) > TOLERANCE)
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			else
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << "cannot be calculted, denominator zero" << endl;
		}
		else
		{
			cout<<"\tChi-square/npoints=\t\t\t\t\t\t\t\t"<<J10<<*(chi.chicomp+index)/(*(ExptsData::sqused+i))<<endl;
			sum=0.0;
			background=0;
			for (j=0;j<*(ExptsData::sqused+i);j++)
			{

				background=(*b + *c * (*px)+ *d * (*px) * (*px)+ *e * (*px) * (*px) * (*px));
				sum+=pow(*a * *py+background,2);
				px++;
				py++;
			}
				
			if (sum > TOLERANCE)
			{
				chi.Rw[index] = sqrt(*(chi.chicomp + index)) * ExptsData::sqsigma[i] / sqrt(sum) * 100;
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			}
			
		}
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\tBin size for this dataset: " << rspacing[index] << " A" << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
		a++;//go to next
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next 
		index++;

	}
	
	for(i=0;i<nfq;i++)
	{
		px = edata.fq_qfinder[i];
		py = edata.fq_ffinder[i];
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout<<"F(Q) set #"<<i+1<<"/"<<nfq;
		cout.setf(ios::scientific, ios::floatfield);
		//the experimental data-background hast to be multiplied with this renormalization factor to fit the calculated data
		cout<<"\tRenorm. a="<<*a;
		//background consists of the constant, linear, quadratic and if use_cubic was given then a cubic part
		//this constant has to be added to the expt data before renormalization for comparison with RMC data		
		cout<<"\tConstant b="<<*b<<endl;
		if (ExptsData::fqfitIQ[i])
		{
			cout << "\tCompton alpha=" << chi.alpha[i] << endl;
			if (ExptsData::fqIQbackgcorr[i])
				cout << "\tmu=" << edata.fqmuact[i] << endl;
			if (ExptsData::fqAXS[i] > 0)
				cout << "\tactual f'=" << ExptsData::fqfprimeact[i * ntypes + ExptsData::fqAXS[i] - 1]<<" for the " << ExptsData::fqAXS[i] << ". type  "<< endl;;
		}
		else
		{
			//this linear coefficient*Q has to be added to the expt data before renormalization for comparison with RMC data
			cout << "\tLinear c=" << *c;
			//this quadratic coefficient*Q*Q has to be added to the expt data before renormalization for comparison with RMC data
			cout << "\tQuadratic d=" << *d << endl;
			//this cubic coefficient*Q*Q*Q has to be added to the expt data before renormalization for comparison with RMC data
			if (ExptsData::use_cubic[index])
			{
				//this cubic coefficient*Q*Q*Q has to be added to the expt data before renormalization for comparison with RMC data
				cout << "\tCubic e=" << *e << endl;
			}
		}
		cout.precision(precision_l);
		if (ExptsData::fquseR[i])
		{
			//cout << "\tRw/(sigma*npoints)=\t\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::fqused + i)) << endl;
			if (fabs(chi.Rw[index] + 1000) > TOLERANCE)
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			else
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << "cannot be calculted, denominator zero" << endl;
		}
		else
		{
			cout << "\tChi-square/npoints=\t\t\t\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::fqused + i)) << endl;
			//caclculating the Rw
			sum = 0.0;
			background = 0.0;
			for (j = 0; j < *(ExptsData::fqused + i); j++)
			{
				background = (*b + *c * (*px) + *d * (*px) * (*px) + *e * (*px) * (*px) * (*px));
				sum += pow(*a * *py + background, 2);
				px++;
				py++;
			}
			if (sum > TOLERANCE)
			{
				chi.Rw[index] = sqrt(*(chi.chicomp + index)) * ExptsData::fqsigma[i] / sqrt(sum) * 100;
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			}
		}
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\tBin size for this dataset: " << rspacing[index] << " A" << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
		if (chi.fit_index[index] > 31)//nonlin regression
		{
			cout << "\tNon-linear regression last accepted chi2_end/chi2_start:\t\t\t" << chi.last_nlr_chi2fr[i] << endl;
			cout << "\tNon-linear regression largest decrease in chi2_end/chi2_start:\t\t\t" << chi.maxdec_nlr_chi2fr[i] << endl;

		}
		a++;//go to next
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next 
		index++;

	}
	
	for (i = 0; i < nfg; i++)
	{
		px = edata.fg_gfinder[i];
		py = edata.fg_ffinder[i];
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout << "F(g) set #" << i + 1 << "/" << nfg;
		cout.setf(ios::scientific, ios::floatfield);
		//the experimental data-background hast to be multiplied with this renormalization factor to fit the calculated data
		cout << "\tRenorm. a=" << *a;
		//background consists of the constant, linear, quadratic and if use_cubic was given then a cubic part
		//this constant has to be added to the expt data before renormalization for comparison with RMC data		
		cout << "\tConstant b=" << *b << endl;
		//this linear coefficient*g has to be added to the expt data before renormalization for comparison with RMC data
		cout << "\tLinear c=" << *c;
		//this quadratic coefficient*g*g has to be added to the expt data before renormalization for comparison with RMC data
		cout << "\tQuadratic d=" << *d << endl;
		//this cubic coefficient*g*g*g has to be added to the expt data before renormalization for comparison with RMC data
		if (ExptsData::use_cubic[index])
		{
			//this cubic coefficient*g*g*g has to be added to the expt data before renormalization for comparison with RMC data
			cout << "\tCubic e=" << *e << endl;
		}
		cout.precision(precision_l);
		if (ExptsData::fguseR[i])
		{
			//cout << "\tRw/(sigma*npoints)=\t\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::fgused + i)) << endl;
			if (fabs(chi.Rw[index] + 1000) > TOLERANCE)
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			else
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << "cannot be calculted, denominator zero" << endl;
		}
		else
		{
			cout << "\tChi-square/npoints=\t\t\t\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::fgused + i)) << endl;
			//caclculating the Rw
			sum = 0.0;
			background = 0.0;
			for (j = 0; j < *(ExptsData::fgused + i); j++)
			{
				background = (*b + *c * (*px) + *d * (*px) * (*px) + *e * (*px) * (*px) * (*px));
				sum += pow(*a * *py + background, 2);
				px++;
				py++;
			}
			if (sum > TOLERANCE)
			{
				chi.Rw[index] = sqrt(*(chi.chicomp + index)) * ExptsData::fgsigma[i] / sqrt(sum) * 100;
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			}
		}
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\tBin size for this dataset: " << rspacing[index] << " A" << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
		a++;//go to next
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next 
		index++;

	}
		
	for(i=0;i<nek;i++)
	{
		py = edata.ek_efinder[i];
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout<<"E(k) set #"<<i+1<<"/"<<nek;
		cout.setf(ios::scientific, ios::floatfield);
		//the experimental data-constant hast to be multiplied with this renormalization factor to fit the calculated data
		cout<<"\tRenorm. a="<<*a;
		//this constant has to be added to the expt data before renormalization for comparison with RMC data		
		cout<<"\tConstant b="<<*b<<endl;
		cout.precision(precision_l);
		if (ExptsData::ekuseR[i])
		{
			//cout << "\tRw/(sigma*npoints)=\t\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::ekused + i)) << endl;
			if (fabs(chi.Rw[index] + 1000) > TOLERANCE)
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			else
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << "cannot be calculted, denominator zero" << endl;
		}
		else
		{
			cout << "\tChi-square/npoints=\t\t\t\t\t\t\t\t" << J10 << *(chi.chicomp + index) / (*(ExptsData::ekused + i)) << endl;
			sum = 0.0;

			background = *b;

			for (j = 0; j < *(ExptsData::ekused + i); j++)
			{
				sum += pow(*a * *py + background, 2);
				py++;
			}
			if (sum > TOLERANCE)
			{
				chi.Rw[index] = sqrt(*(chi.chicomp + index)) * ExptsData::eksigma[i] / sqrt(sum) * 100;
				cout << "\tRw for the set =\t\t\t\t\t\t\t\t" << J10 << chi.Rw[index] << " %" << endl;
			}
		}
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\tBin size for this dataset: " << rspacing[index] << " A" << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
		if (ExptsData::ek_ngrid_in[i]>0)
		{//the k axis and E0 exis is going in opposite directions
			cout << "\tE0 shift grid index:\t\t\t\t\t\t\t\t" << ExptsData::ek_gridind[i]+ExptsData::ek_gridfrom[i] << endl;
			cout<<"\tdE0 actual:\t\t\t\t\t\t\t\t\t" << -(ExptsData::ek_gridind[i]+ExptsData::ek_gridfrom[i])* ExptsData::ek_dE0[i] / ExptsData::ek_ngrid_in[i] << " eV"<<endl;
		}
		a++;//go to next
		b++;//go to next
		index++;

	}
	for(i=0;i<ncosdistr;i++)
	{
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout<<"Cosine distr. #"<<i+1<<"/"<<ncosdistr;
		cout.precision(precision_l);
		cout.setf(ios::scientific, ios::floatfield);
		cout<<"\tChi-square/npoints=\t\t\t\t\t\t"<<J10<<*(chi.chicomp+index)/CosDistrConst::ncos_bin<<endl;
		index++;
	}
	//Individual coordination constraints information
	for(i=0;i<nicoord;i++)
	{
		for (isubconst=0;isubconst<icc.n_subconst[i];isubconst++)
		{
			cout.precision(precision_s);
			cout.setf(ios::fixed, ios::floatfield);
			cout.setf(ios::right, ios::adjustfield);
			cout<<"Coordination constraint "<<i+1<<"/"<<icc.nconstraints<<" subconstraint "<<isubconst+1<<endl;
			cout<<"\tFraction="<<*(icc.nsatisfy+icc.cum_n_subconst[i]+isubconst)<<"/"<<*(icc.ncentral+i)<<"="\
			<<double(*(icc.nsatisfy+icc.cum_n_subconst[i]+isubconst))/double(*(icc.ncentral+i))*100.0<<\
			"%, target="<<100* *(icc.fraction+icc.cum_n_subconst[i]+isubconst)<<"%"<<endl;
			cout.precision(precision_l);
			cout.setf(ios::scientific, ios::floatfield);
			cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
			index++;
		}
	}

	//Average coordination constraints information
	for(i=0;i<navcoord;i++)
	{
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout<<"Average coordination constraint "<<i+1<<"/"<<avcc.nconstraints<<endl;
		cout<<"\tAv. Coord. number is="<<double(*(avcc.neighbourcount+i))/double(*(avcc.ncentral+i))<<\
		", target="<< *(avcc.acnreq+i)<<endl;
		cout.precision(precision_l);
		cout.setf(ios::scientific, ios::floatfield);
		cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
		index++;
	}
#ifdef _ADVANCED_GEOM_CONST
	for (i=0;i<ncommonneigh;i++)
	{
		double temp;
		if (double(*(conc.nprimary + i)) == 0)
			temp = 0;
		else
			temp = double(*(conc.nsatisfy + i)) / double(*(conc.nprimary + i)) * 100.0;
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout << "Common neighbour constraint " << i + 1 << "/" << ncommonneigh << endl;
		cout << "\tFraction=" << *(conc.nsatisfy + i) << "/" << *(conc.nprimary + i) << "="\
			<< temp<< "%, target=" << 100 * *(conc.fraction + i) << "%" << endl;
		cout.precision(precision_l);
		cout.setf(ios::scientific, ios::floatfield);
		cout << "\tChi-square add-on:\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
		index++;
	}
	for (i = 0; i < nsecondneigh; i++)
	{
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout << "Second neighbour constraint " << i + 1 << "/" << nsecondneigh << endl;
		cout << "\tFraction=" << *(snc.nsatisfy + i) << "/" << *(snc.ncentral + i) << "="\
			<< double(*(snc.nsatisfy + i)) / double(*(snc.ncentral + i))*100.0 << \
			"%, target=" << 100 * *(snc.fraction + i) << "%" << endl;
		cout.precision(precision_l);
		cout.setf(ios::scientific, ios::floatfield);
		cout << "\tChi-square add-on:\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
		index++;
	}
	for (i = 0; i < nbvs; i++)
	{
		cout.precision(precision_s);
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);
		cout << "BVS constraint " << i + 1 << "/" << nbvs << endl;
		cout << "\tTarget valence="<<bvs.target_valence[i]<< endl;
		cout.precision(precision_l);
		cout.setf(ios::scientific, ios::floatfield);
		cout << "\tChi-square add-on:\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
		index++;
	}
#endif

#ifdef _LOCAL_INV
	int itype,jtype;
	for (i=0;i<nlocint;i++)
	{
		cout<<"Local invariance interval "<<i+1<<"/npoints "<<endl;
		cout.precision(precision_l);
		cout.setf(ios::scientific, ios::floatfield);
		if (loc_npoints[i]>0)
			cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)/loc_npoints[i]<<endl;
		else
			cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
		if (RunParams::loc_chi2_mode)
		{
			cout.precision(precision_s);
			cout.setf(ios::fixed, ios::floatfield);
			cout.setf(ios::right, ios::adjustfield);
			cout<<"\taverage distances of for the first and last atom of the intervals"<<endl;
			for (itype=0;itype<ntypes;itype++)
			{
				for (jtype=0;jtype<ntypes;jtype++)
				{
					cout<<"\t "<<itype+1<<"-"<<jtype+1<<". partial from "<<chi.min_loc_r_now[itype*ntypes*nlocint+jtype*nlocint+i]<<" to "<<chi.max_loc_r_now[itype*ntypes*nlocint+jtype*nlocint+i]<<" A"<<endl;              
				}
			}
		}
		index++;
	}
#endif

	//Potential
	char *potname,*potname2,*name_ext=NULL;
	SetArraysize(&potname,50,"potname","RunParams::PrintStatus");
	SetArraysize(&potname2,50,"potname2","RunParams::PrintStatus");
	
	if (potential>0)
	{
		cout.precision(precision_l);
		cout.setf(ios::scientific);
		if (fabs(NB_weight_mode)!=1)
		{

			cout<<vdW_name<<" potential\t\t\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
			index++;
		}
		else
		{
			for (i=0;i<ChiSquared::nNB;i++)
			{
				switch (potential)
				{
					case 1:
					{
						cout << vdW_name << " potential partial  " << i + 1 << ".\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
						break;
					}
					case 10:
					{
						cout << vdW_name << " potential partial  " << FNC_POT::tabpot_index[i] + 1 << ".\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
						cout << "E/kT: " << FNC_POT::vdW_pot[i] * 1000 / k_BOLTZMANN / temperature << endl;//using kJ for U in the program
					}
					
				}
				index++;
			}
		}
		if (potential == 1)
		{
			if (fabs(NB_weight_mode) != 1)
			{
				cout << "Coulomb potential\t\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
				index++;
			}
			else
			{
				for (i = 0; i < ChiSquared::nNB; i++)
				{
					cout << "Coulomb potential partial " << i + 1 << ".\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
					index++;
				}
			}

			if (Topology::npair_types > 0)
			{
				if (fabs(NB_weight_mode) != 1)
				{

					cout << "1-4 " << vdW_name << " potential\t\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
					index++;
				}
				else
				{
					for (i = 0; i < ChiSquared::nNB; i++)
					{
						if (RunParams::used_14partials[i])//only display, if it was used
							cout << "1-4 " << vdW_name << " potential partial " << i + 1 << ".\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
						index++;
					}
				}

				if (fabs(NB_weight_mode) != 1)
				{

					cout << "1-4 Coulomb potential\t\t\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
					index++;
				}
				else
				{
					for (i = 0; i < ChiSquared::nNB; i++)
					{
						if (RunParams::used_14partials[i])
							cout << "1-4 Coulomb potential partial " << i + 1 << ".\t\t\t\t\t\t\t" << *(chi.chicomp + index) << endl;
						index++;
					}
				}
			}
		}
		cout.unsetf(ios::scientific);

	}

	if (nbond_types>0)
	{
		if (Topology::bond_weight_mode==0)
		{
			cout.precision(precision_l);
			cout.setf(ios::scientific);
			cout<<"Bond chi-square add-on:\t\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
			cout.unsetf(ios::scientific);
			index++;
		}
		else
		{
			for (i=0; i<ChiSquared::nbonds;i++)
			{
				cout.precision(precision_l);
				cout.setf(ios::fixed, ios::floatfield);
				cout.setf(ios::right, ios::adjustfield);
				cout<<"Bond type "<<i+1<<"/"<<nbond_types<<endl;
				cout<<"\tEquilibrium dist.: "<<Topology::bond_r[i]<<" A, force const.: "<< Topology::bond_k[i]<<" kj/mol/A^2"<<endl;
				cout.setf(ios::scientific, ios::floatfield);
				cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
				index++;
			}
		}
	}

	if (nangle_types>0)
	{
		if (Topology::angle_weight_mode==0)
		{
			cout.precision(precision_l);
			cout.setf(ios::scientific);
			cout<<"Angle chi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
			cout.unsetf(ios::scientific);
			index++;
		}
		else
		{
			for (i=0; i<ChiSquared::nangles;i++)
			{
				cout.precision(precision_l);
				cout.setf(ios::fixed, ios::floatfield);
				cout.setf(ios::right, ios::adjustfield);
				cout<<"Angle type "<<i+1<<"/"<<nangle_types<<endl;
				cout<<"\tEquilibrium angle : "<<Topology::angle_ang[i]<<" degree, force const.: "<< Topology::angle_k[i]<<" kj/mol/rad^2"<<endl;
				cout.setf(ios::scientific, ios::floatfield);
				cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
				index++;
			}
		}
	}

	if (nperdihedral_types>0)
	{
		if (Topology::force_field_type==OPLSAA)
				mystrcpy(potname,50,"Improper dihedral");
			else
				mystrcpy(potname, 50, "Proper dihedral");
		if (Topology::perdih_weight_mode==0)
		{
			cout.precision(precision_l);
			cout.setf(ios::scientific);
			cout<<potname<<" chi-square add-on:\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
			cout.unsetf(ios::scientific);
			index++;
		}
		else
		{
			for (i=0; i<ChiSquared::nperdihs;i++)
			{
				mystrcpy(potname2, 50, potname);
				
				IntToStr(&name_ext,i+1);
				mystrcat(potname2, 50, name_ext);
				cout<<potname2<<"/"<<nperdihedral_types<<endl;
				cout.precision(precision_l);
				cout.setf(ios::fixed, ios::floatfield);
				cout.setf(ios::right, ios::adjustfield);
				cout<<"\tEquilibrium angle : "<<Topology::perdihedral_ang[i]<<" degree, force const.: "<< Topology::perdihedral_k[i]<<" kj/mol"<<endl;
				cout.setf(ios::scientific, ios::floatfield);
				cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
				index++;
			}
		}
	}

	if (nharmdihedral_types>0)
	{
		if (Topology::harmdih_weight_mode==0)
		{
			cout.precision(precision_l);
			cout.setf(ios::scientific);
			cout<<"Improper dihedral chi-square add-on:\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
			cout.unsetf(ios::scientific);
			index++;
		}
		else
		{
			for (i=0; i<ChiSquared::nharmdihs;i++)
			{
				cout.precision(precision_l);
				cout.setf(ios::fixed, ios::floatfield);
				cout.setf(ios::right, ios::adjustfield);
				cout<<"Improper dihedral type "<<i+1<<"/"<<nharmdihedral_types<<endl;
				cout<<"\tEquilibrium angle : "<<Topology::harmdihedral_ang[i]<<" degree, force const.: "<< Topology::harmdihedral_k[i]<<" kj/mol/rad^2"<<endl;
				cout.setf(ios::scientific, ios::floatfield);
				cout<<"\tChi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
				index++;
			}
		}
	}

	if (nRBdihedral_types>0)
	{
		int count;
		if (Topology::force_field_type == OPLSAA)
		{
			mystrcpy(potname, 50, "Proper dihedral");
			count = 7;
		}
		else
		{
			mystrcpy(potname, 50, "RB dihedral");
			count = 8;
		}
		if (Topology::RBdih_weight_mode==0)
		{
			cout.precision(precision_l);
			cout.setf(ios::scientific);
			cout << potname << " chi-square add-on:";
			for (int ii = 0; ii < count; ii++)
				cout << "\t";
			cout<<*(chi.chicomp+index)<<endl;
			cout.unsetf(ios::scientific);
			index++;
		}
		else
		{
			for (i=0; i<ChiSquared::nRBdihs;i++)
			{
				mystrcpy(potname2, 50, potname);
				IntToStr(&name_ext,i+1);
				mystrcat(potname2, 50, name_ext);
				cout.precision(precision_l);
				cout.setf(ios::fixed, ios::floatfield);
				cout.setf(ios::right, ios::adjustfield);
				cout<<potname2<<"/"<<nRBdihedral_types<<endl;
				cout<<"\tC0 : "<<Topology::RBdihedral_C0[i]<<" C1 : "<<Topology::RBdihedral_C1[i]<<" C2 : "<<Topology::RBdihedral_C2[i]<<endl;
				cout<<"\tC3 : "<<Topology::RBdihedral_C3[i]<<" C4 : "<<Topology::RBdihedral_C4[i]<<" C5 : "<<Topology::RBdihedral_C5[i]<<" kj/mol"<<endl;
				cout.setf(ios::scientific, ios::floatfield);
				cout<<potname<<" chi-square add-on:\t\t\t\t\t\t\t\t"<<*(chi.chicomp+index)<<endl;
				index++;
			}
		}
	}


	cout<<endl;
	cout.unsetf(ios::floatfield);
	cout.unsetf(ios::right);
	cout.precision(6);	

#ifdef _AV_MOVE
	//Displaying the average moved distance of the atoms
	if (status)
		cout<<"The average distance an atom moved:	"<<config.CalcAvMove()<<" Angstrom"<<endl;
#endif
	if (FNC_POT::sum_out_of_range>0)
	{
		for (i=0;i<FNC_POT::nconstraints;i++)
			cout<<"The number of FNC pairs still out of FNC range for constraint "<<i+1<<" : "<<FNC_POT::nout_of_range[i]<<endl;
	}
	if (name_ext != NULL)
		delete [] name_ext;
	delete [] potname;
	delete [] potname2;
}

//-----------print----------------
//displays the Run parameters on screen
void RunParams::PrintRunParams(ExptsData &edata) const
{
	int i,j,k,m,ineightype,isubconst;
	//int l;
	int npartials=ntypes*(ntypes+1)/2;//number of cutoffs
	cout<<endl;
	cout<<"*************************************************************************"<<endl;
	cout<<"****    THE RMC++ PROGRAM WILL RUN WITH THE FOLLOWING PARAMETERS:    ****"<<endl;
	cout<<"*************************************************************************"<<"\n"<<endl;

	cout<<"Title: "<<title<<endl;//title
	if (custmove)
		cout<<nthreads<<"-threaded molecular RMC++ with "<<nmoved<<" moved atom(s)."<<endl;
	else
		cout<<nthreads<<"-threaded atomic RMC++ with "<<nmoved<<" moved atom(s)."<<endl;
	if (fnc>0)
		cout<<"FNC will be used."<<endl;
	else
	{
		if (ntypes>1)
		{
			cout<<"Fraction of swaps: "<<swap_fraction<<endl;
			cout<<"Number of swap pair types: "<<nswap_pairs<<endl;
			cout<<"Type of atoms in allowed swap pairs: "<<endl;
			for (i=0;i<nswap_pairs;i++)
				cout<<"\t"<<i+1<<". swap pair type "<<swap_type1[i]+1<<" - "<<swap_type2[i]+1<<endl;
		}
	}

	if (reload)
	{
		if (CoordNumbConst::nconstraints>0 || AvCoordConst::nconstraints>0)
			cout<<"Histogram and coordination numbers will be loaded from files!"<<endl;
		else
			cout<<"Histogram will be loaded from file!"<<endl;
	}
	if (nmoved==1)
	{
		if (moveout)
			cout<<"Moveout option will be used!"<<endl;
		else
			cout<<"Moveout option will not be used!"<<"\n"<<endl;
	}
	if (ExptsData::is_IQmucorr)
		cout<<"Scalable I(Q) background correction will be performed in each "<<RunParams::IQ_backg_corr_step<<". step"<<endl;
	if (ExptsData::is_AXS)
		cout << "Anomalous X-ray sacttering f' will be shifted in each "<<RunParams::AXS_shiftstep << ". step" << endl;
	if (ExptsData::is_E0shift)
		cout << "EXAFS coefficient shift will be performed in each " << RunParams::exafs_shiftstep << ". step" << endl;
	cout.precision(6);
	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::right, ios::adjustfield);;//formatting output
	cout<<"\n************************CONFIGURATION RELATED****************************"<<"\n"<<endl;
	cout<<"The half-length of the simulation box is "<<boxedge<<" Angstrom."<<endl;
#ifdef _NO_PERIODIC
	cout<<"Non-periodic boundary conditions are used, the spherical sample is located at the "<<endl;
	cout<<"middle of the simulation box, with radius smaller than half box lenth /2"<<endl;
	cout<<"Radius of the spherical sample: "<<RunParams::R0<<" Angstrom."<<endl;
	cout<<"The number density is defined as the total number of atoms/volume of the spherical sample"<<endl;
#endif
	
	cout<<"The number density is "<<rho<<" atoms per cubic Angstrom."<<endl;
#ifdef _VIBR_AMP
	cout<<"\n_VIBR_AMP option is ON. The thermal atomic vibrations will be calculated by the "<<endl;
	cout<<"convolution of the histogram with a Gaussian distribution with sigma parameter"<<endl;
	cout<<"for the different partials:"<<endl;
	int ipartial=0;
	for (i=0;i<ntypes;i++)
	{
		for (j=i;j<ntypes;j++)
		{
			cout<<"\t["<<i<<","<<j<<"] partial "<<HistoSet::T_corr_sigma[ipartial]<<" Angstrom "<<endl;
			ipartial++;
		}
	}
#endif

	cout.unsetf(ios::fixed);
	cout.unsetf(ios::right);
	cout<<"\n"<<J10<<natoms<<"\t atoms of "<<ntypes<<" type(s)"<<endl;
	for (i = 0; i < ntypes; i++)
	{
		cout << J10 << pnatoms[i] << "\t atoms of type " << i + 1 << " , maximum movement " << pmaxmove[i] << " Angstrom.";
		if (RunParams::chem_symbols != nullptr)
			cout << " Symbol: " << RunParams::chem_symbols[i];
		cout << endl;
	}
	cout<<"\n"<<J7<<npartials<<" partial(s). The distances of closest approach for [i,j]-th partial:"<<endl;
	k=0;
	//l=
	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::right, ios::adjustfield);//formatting output
	for (i=0;i<ntypes;i++)
		for (j=i;j<ntypes;j++)
			cout<<"\t["<<i+1<<","<<j+1<<"] "<<J10<<*(pcutoff+k++)<<"\t Angstrom."<<endl;

	cout<<"\n**************************HISTOGRAM RELATED******************************"<<"\n"<<endl;
	cout << "The same histograms is used for the same binsize and integer difference in the bin shifts, " << endl;
	cout << "\tcalculated for the smallest binshift in the series!" << endl;
	for (i = 0; i < ndiffbin; i++)
	{
		cout.precision(6);
		cout << i + 1 << ". unique bin size: " << rspacing_diff[i] << " Angstrom";
		cout.precision(2);
		if (binshift_flag)
			cout << " with bin shift: " << binshift_diff[i];
		cout<< endl;
		cout << "\tThe number of bins: " << nbins[i] << endl;
		
		cout << "\tThe histogram will be calculated from " << xmin[i] * boxedge << " up to " << xmax[i] * boxedge << " Angstrom."  << endl;
	}
	cout.precision(8);
	cout.unsetf(ios::fixed);
	cout.unsetf(ios::right);
	cout<<"\n***************************OUTPUT RELATED********************************"<<"\n"<<endl;
#ifdef _TEST_MODE 
	cout << "Running in TEST MODE" << endl;
#endif
	switch (runmode)
	{
	case 1:
		cout << J10 << -runlimit << "\t steps will be generated." << endl;
		break;
	case 2:
		cout << J10 << runlimit << "\t steps will be accepted." << endl;
		break;
	default:
		cout << J10 << runlimit << "\t minutes time limit for the run." << endl;
	}
	cout << J10 << printstep << "\t generated steps for printing." << endl;
	cout<<J10<<timesave<<"\t minutes time interval for saving."<<"\n"<<endl;
	cout<<J10<<cfgnumb<<"\t configuration(s) to collect."<<endl;

	cout<<"\n**************************CONSTRAINS RELATED*****************************"<<"\n"<<endl;
	//g(r) constraints
	cout<<"\n"<<J10<<ngr<<" g(r) DATA SET";
	if (ngr > 1)
		cout << "S";
	cout << endl;
	m=0;//to increment the coefficient index
	for (i=0;i<ngr;i++)
	{
		//l=
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);//formatting output
		cout << "\n" << J4<<i+1<<". g(r) data set: "<<endl;
		cout<<"\t Number of r values in the file: "<<ExptsData::grsize[i];
		cout<<"\t Number of USED r values: "<<ExptsData::grused[i]<<endl;
		cout<<"\t Coherent coefficients for [i,j]-th partial:"<<endl;
	
		for(j=0;j<ntypes;j++)
			for (k=j;k<ntypes;k++)
				cout<<"\t ["<<j+1<<","<<k+1<<"] "<<J10<<*(ExptsData::gr_coeffs+m++)<<endl;
		
		cout<<"\t Constant to substarct from the g(r) data: "<<ExptsData::grsub[i]<<endl;
		if (ExptsData::grrenorm[i])
			cout<<"\t Renormalisation will be used!"<<endl;
		else
			cout<<"\t Renormalisation will not be used!"<<endl;
		if (ExptsData::use_cubic[i])
		{
			if (ExptsData::groffset[i])
				cout<<"\t Offset will be used!"<<endl;
			else
				cout<<"\t Offset will not be used!"<<endl;
			if (ExptsData::grlinear[i])
				cout<<"\t Linear background will be used!"<<endl;
			else
				cout<<"\t Linear background will not be used!"<<endl;
			if (ExptsData::grquadratic[i])
				cout<<"\t Quadratic background will be used!"<<endl;
			else
				cout<<"\t Quadratic background will not be used!"<<endl;
			if (ExptsData::grcubic[i])
				cout<<"\t Cubic background will be used!"<<endl;
			else
				cout<<"\t Cubic background will not be used!"<<endl;
		}
		cout.setf(ios::scientific, ios::floatfield);
		cout<<"\t Sigma factor : "<<ExptsData::grsigma[i]<<endl;
		if (ExptsData::gruseR[i])
			cout << "\t Rw/sigma will be used instead of the normal chisquare"<< endl;
		cout.unsetf(ios::scientific);
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\t Bin size for this dataset: " << rspacing[i] << " A" << endl;
			cout << "\t Bin shift for this dataset: " << binshift[i] << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
	};

	//S(Q) constraints
	cout << "\n" << J10 <<nsq<<" S(Q) DATA SET";
	if (nsq > 1)
		cout << "S";
	cout << endl;
	m=0;//to increment the coefficient index
	for (i=0;i<nsq;i++)
	{
		//l=
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);//formatting output
		cout << "\n" << J4<<i+1<<". S(Q) data set: "<<endl;
		cout<<"\t Number of data points in the file: "<<ExptsData::sqsize[i];
		cout<<"\t Number of USED data points: "<<ExptsData::sqused[i]<<endl;
		if (ExptsData::sqreadcoeffs[i] == 0)
		{
			cout << "\t Coefficients were calculated for " << endl;
			for (j = 0; j < ntypes; j++)
			{
				cout << "\t "<<j + 1 << ". type " << chem_symbols_standard[j] << endl;
				if (ExptsData::isotope_count[i*ntypes+j].count == 0 && chem_symbols_standard[j].compare("dummy")!=0)
					cout << "\t\tnatural abundance" << endl;
				else
					for (k = 0; k < ExptsData::isotope_count[i * ntypes + j].count; k++)
						cout << "\t\t" << ExptsData::isotope_count[i * ntypes + j].isotopes[k].symbol << " ratio: " << ExptsData::isotope_count[i * ntypes + j].isotopes[k].ratio << endl;
			}
		}
		cout<<"\t Coherent coefficients for [i,j]-th partial:"<<endl;
	
		for(j=0;j<ntypes;j++)
			for (k=j;k<ntypes;k++)
				cout<<"\t ["<<j+1<<","<<k+1<<"] "<<J10<<*(ExptsData::sq_coeffs+m++)<<endl;
		
		cout<<"\t Constant to substarct from the S(Q) data: "<<ExptsData::sqsub[i]<<endl;
		if (ExptsData::sqrenorm[i])
			cout<<"\t Renormalisation will be used!"<<endl;
		else
			cout<<"\t Renormalisation will not be used!"<<endl;
		if (ExptsData::sqoffset[i])
			cout<<"\t Offset will be used!"<<endl;
		else
			cout<<"\t Offset will not be used!"<<endl;
		if (ExptsData::sqlinear[i])
			cout<<"\t Linear background will be used!"<<endl;
		else
			cout<<"\t Linear background will not be used!"<<endl;
		if (ExptsData::sqquadratic[i])
			cout<<"\t Quadratic background will be used!"<<endl;
		else
			cout<<"\t Quadratic background will not be used!"<<endl;
		if (ExptsData::use_cubic[i+ngr])
		{
			if (ExptsData::sqcubic[i])
				cout<<"\t Cubic background will be used!"<<endl;
			else
				cout<<"\t Cubic background will not be used!"<<endl;
		}
		cout.setf(ios::scientific, ios::floatfield);
		cout<<"\t Sigma factor : "<<ExptsData::sqsigma[i]<<endl;
		if (ExptsData::squseR[i])
			cout << "\t Rw/sigma will be used instead of the normal chisquare" <<  endl;
		cout.unsetf(ios::scientific);
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\t Bin size for this dataset: " << rspacing[i+ngr] << " A" << endl;
			cout << "\t Bin shift for this dataset: " << binshift[i+ngr] << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
	};

	//F(Q) constraints
	cout << "\n" << J10 <<nfq<<" F(Q) DATA SET";
	if (nfq > 1)
		cout << "S";
	cout << endl;
	for (i=0;i<nfq;i++)
	{
		//l=
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);//formatting output
		cout << "\n" << J4<<i+1<<". F(Q) data set: "<<endl;
		cout<<"\t Number of data points in the file: "<<ExptsData::fqsize[i];
		cout<<"\t Number of USED data points: "<<ExptsData::fqused[i]<<endl;
		if (!ExptsData::fqIQbackgcorr[i])
			cout<<"\t Constant to substarct from the F(Q) data: "<<ExptsData::fqsub[i]<<endl;
		if (ExptsData::fqrenorm[i])
			cout<<"\t Renormalisation will be used!"<<endl;
		else
			cout<<"\t Renormalisation will not be used!"<<endl;
		if (ExptsData::fqoffset[i])
			cout<<"\t Offset will be used!"<<endl;
		else
			cout<<"\t Offset will not be used!"<<endl;
		if (ExptsData::fqlinear[i])
			cout<<"\t Linear background will be used!"<<endl;
		else
			cout<<"\t Linear background will not be used!"<<endl;
		if (ExptsData::fqquadratic[i])
			cout<<"\t Quadratic background will be used!"<<endl;
		else
			cout<<"\t Quadratic background will not be used!"<<endl;
		if (ExptsData::use_cubic[i+ngr+nsq])
		{
			if (ExptsData::fqcubic[i])
				cout<<"\t Cubic background will be used!"<<endl;
			else
				cout<<"\t Cubic background will not be used!"<<endl;
		}
		
		if (!ExptsData::fqreadcoeffs[i])
		{
			if (ExptsData::fqfitIQ[i])
				cout << "\t f(Q)-s, coefficients and A(Q) were calculated" << endl;
			else
				cout << "\t f(Q)-s and normalized coefficients were calculated" << endl;
		}
		else if (ExptsData::fqfitIQ[i])
			cout << "\t Coefficients were calculated" << endl;
		if (ExptsData::fqfitIQ[i])
		{
			cout << "\t I(Q) is fitted instead of structure factor" << endl;
			if (ExptsData::fqusecompton[i])
			{
				cout << "\t\t Compton term will be used" << endl;
				if (ExptsData::fqrenalpha[i])
					cout << "\t\t Compton alpha parameter renormalization will be used!" << endl;
				else
					cout << "\t\t Compton alpha parameter renormalization will not be used!" << endl;
			}
			else
				cout << "\t\t Compton term will not be used" << endl;
						
			if (ExptsData::fqIQbackgcorr[i])
				cout<< "\t\t Scalable background correction with mu: "<<ExptsData::fqmu[i]<<", dmu_max: "<<ExptsData::fqdmumax[i]<< endl;
			if (ExptsData::fqAXS[i] > 0)
			{
				cout << "\t\t AXS is used for atom type: " << ExptsData::fqAXS[i] << " with f': " << ExptsData::fqfprime[i * ntypes + ExptsData::fqAXS[i]-1];
				cout << " with fraction to vary: " << ExptsData::fqfprimefactor[i] << endl;
			}	
		}
		if ((ExptsData::fqAXS[i] > 0 && ExptsData::fqfprimecount[i] > 1) || (ExptsData::fqAXS[i] == 0 && ExptsData::fqfprimecount[i] > 0))//there are fixed f' using sets
		{
			cout << "\t Fixed f' using atom types:" << endl;
			for (int itype = 0; itype < ntypes; itype++)
			{
				if ((itype!= ExptsData::fqAXS[i]-1) && ExptsData::fqfprimeindex[i] & (int)pow(2, itype))//this type has fixed f' on it is not AXS
					cout << "\t\t" << itype + 1 << ". type with f': " << ExptsData::fqfprime[i * ntypes + itype] << endl;
			}			
		}
		cout.setf(ios::scientific, ios::floatfield);
		cout<<"\t Sigma factor : "<<ExptsData::fqsigma[i]<<endl;
		if (ExptsData::fquseR[i])
			cout << "\t Rw/sigma will be used instead of the normal chisquare" <<  endl;
		cout.unsetf(ios::scientific);
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\t Bin size for this dataset: " << rspacing[i + ngr+nsq] << " A" << endl;
			cout << "\t Bin shift for this dataset: " << binshift[i + ngr+nsq] << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
	};
	//F(g) constraints
	cout << "\n" << J10 << nfg << " F(g) DATA SET" ;
	if (nfg > 1)
		cout << "S";
	cout << endl;
	for (i = 0; i < nfg; i++)
	{
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);//formatting output
		cout << "\n" << J4 << i + 1 << ". F(g) data set: " << endl;
		cout << "\t Number of data points in the file: " << ExptsData::fgsize[i];
		cout << "\t Number of USED data points: " << ExptsData::fgused[i] << endl;

		cout << "\t Constant to substarct from the F(g) data: " << ExptsData::fgsub[i] << endl;
		if (ExptsData::fgrenorm[i])
			cout << "\t Renormalisation will be used!" << endl;
		else
			cout << "\t Renormalisation will not be used!" << endl;
		if (ExptsData::fgoffset[i])
			cout << "\t Offset will be used!" << endl;
		else
			cout << "\t Offset will not be used!" << endl;
		if (ExptsData::fglinear[i])
			cout << "\t Linear background will be used!" << endl;
		else
			cout << "\t Linear background will not be used!" << endl;
		if (ExptsData::fgquadratic[i])
			cout << "\t Quadratic background will be used!" << endl;
		else
			cout << "\t Quadratic background will not be used!" << endl;
		if (ExptsData::use_cubic[i + ngr + nsq + nfq])
		{
			if (ExptsData::fgcubic[i])
				cout << "\t Cubic background will be used!" << endl;
			else
				cout << "\t Cubic background will not be used!" << endl;
		}
		cout.setf(ios::scientific, ios::floatfield);
		cout << "\t Sigma factor : " << ExptsData::fgsigma[i] << endl;
		if (ExptsData::fguseR[i])
			cout << "\t Rw/sigma will be used instead of the normal chisquare" << endl;
		cout.unsetf(ios::scientific);
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\t Bin size for this dataset: " << rspacing[i + ngr+nsq+nfq] << " A" << endl;
			cout << "\t Bin shift for this dataset: " << binshift[i + ngr + nsq + nfq] << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
	};

	//E(k) constraints
	cout << "\n" << J10 <<nek<<" E(k) DATA SET";
	if (nek > 1)
		cout << "S";
	cout << endl;
	for (i=0;i<nek;i++)
	{
		//l=
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::right, ios::adjustfield);//formatting output
		cout << "\n" << J4<<i+1<<". E(k) data set: "<<endl;
		cout<<"\t Type of the absorbing particle: "<<ExptsData::ek_abstype[i]+1<<endl;
		cout<<"\t Number of data points in the file: "<<ExptsData::eksize[i];
		cout<<"\t Number of USED data points: "<<ExptsData::ekused[i]<<endl;
		
		if (ExptsData::ekrenorm[i])
			cout<<"\t Renormalisation will be used!"<<endl;
		else
			cout<<"\t Renormalisation will not be used!"<<endl;
		if (ExptsData::ekoffset[i])
			cout<<"\t Offset will be used!"<<endl;
		else
			cout<<"\t Offset will not be used!"<<endl;
		if (ExptsData::ek_ngrid_in[i] > 0)
		{
			cout << "\t E0 shift is applied with dE0:	           " << ExptsData::ek_dE0[i] << " eV" << endl;
			cout << "\t Number of grid points in one direction : " << ExptsData::ek_ngrid_in[i] << endl;
			cout << "\t Grid points used from " << ExptsData::ek_gridfrom[i] << " to " << ExptsData::ek_gridto[i] << endl;
			cout << "\t Starting with grid index " << ExptsData::ek_gridind[i]+ ExptsData::ek_gridfrom[i] << endl;
		}
		cout.setf(ios::scientific, ios::floatfield);
		cout << "\t E(k)*(Q^weight) factor : " << ExptsData::ekchipower[i] << endl;
		cout<<"\t Sigma factor : "<<ExptsData::eksigma[i]<<endl;
		if (ExptsData::ekuseR[i])
			cout << "\t Rw/sigma will be used instead of the normal chisquare" << endl;
		if (ndiffbin > 1)
		{
			cout.unsetf(ios::scientific);
			cout.setf(ios::fixed, ios::floatfield);
			cout << "\t Bin size for this dataset: " << rspacing[i + ngr + nsq + nfq+nfg] << " A" << endl;
			cout << "\t Bin shift for this dataset: " << binshift[i + ngr + nsq + nfq+nfg] << endl;
			cout.unsetf(ios::fixed);
			cout.setf(ios::scientific, ios::floatfield);
		}
		
		
		cout.unsetf(ios::scientific);
	};

	//Cosine distribution of bond angle constraints
	cout << "\n" << J10 <<ncosdistr<<" COSINE DISTRIBUTION OF BOND ANGLES CONSTRAINT";
	if (ncosdistr > 1)
		cout << "S";
	cout << endl;
	for (i=0;i<ncosdistr;i++)
	{
		cout << "\n" << J4<<i+1<<". constraint: "<<endl;
		//The type has to start with one at the output
		cout<<"\t Central atom type: "<<CosDistrConst::central[i]+1;
		cout<<"\t Neighbour atom types: ";
		cout<<"\t "<<CosDistrConst::neighbour1[i]+1<<", "<<CosDistrConst::neighbour2[i]+1<<endl;
		cout<<"\t For 1-st neighbour from "<<CosDistrConst::dmin1[i]<<" A to "<<CosDistrConst::dmax1[i]<<" A"<<endl;
		cout<<"\t For 2-nd neighbour from "<<CosDistrConst::dmin2[i]<<" A to "<<CosDistrConst::dmax2[i]<<" A"<<endl;
		switch (CosDistrConst::method[i])
		{
			case 0:
			{
				cout<<"\t Calculation method is a step function, desired range is from "<<\
				CosDistrConst::angle[i]-CosDistrConst::wcontrol[i]<<" to "<<CosDistrConst::angle[i]+CosDistrConst::wcontrol[i]<<" degree"<<endl;
				break;
			}
			case 1:
			{
				cout<<"\t Calculation method is normal distribution with maximum at "<<CosDistrConst::angle[i]<<" degree and sigma "<<CosDistrConst::wcontrol[i]<<endl;
				break;
			}
			case 2:
			{
				cout<<"\t No angles are desired around "<<CosDistrConst::angle[i]<<" degree with normal distribution sigma "<<CosDistrConst::wcontrol[i]<<endl;
				break;
			}
			case 3:
			{
				cout<<"\t Experimental cosine distribution is loaded from file: "<<CosDistrConst::datafilename[i]<<endl;
			}
		}
		cout<<"\t Weight factor: "<<CosDistrConst::weights[i]<<endl;
	};

	//Coordination number constraints
	cout << "\n" << J10 <<nicoord<<" COORDINATION NUMBER CONSTRAINT";
	if (nicoord > 1)
		cout << "S";
	cout << endl;
	for (i=0;i<nicoord;i++)
	{
		cout << "\n" << J4<<i+1<<". constraint: "<<endl;
		//The type has to start with one at the output
		cout<<"\t Central atom type: "<<CoordNumbConst::central[i]+1;
		cout<<"\t Number of neighbour types: "<<CoordNumbConst::n_neightype[i]<<endl;
		
		for (ineightype=0;ineightype<CoordNumbConst::n_neightype[i];ineightype++)
		{
			cout<<"\t "<<ineightype+1<<". Neighbour atom type: "<<CoordNumbConst::neighbours[CoordNumbConst::cum_n_neightype[i]+ineightype]+1;
			cout<<" from "<<CoordNumbConst::dmin[CoordNumbConst::cum_n_neightype[i]+ineightype];
			cout<<" to "<<CoordNumbConst::dmax[CoordNumbConst::cum_n_neightype[i]+ineightype]<<" A "<<endl;
		}		
		for (isubconst=0;isubconst<CoordNumbConst::n_subconst[i];isubconst++)
		{
			cout<<J10<<isubconst+1<<". Subconstraint: "<<endl;
			cout<<"\t\t Target coordination number: "<<CoordNumbConst::target_coord[CoordNumbConst::cum_n_subconst[i]+isubconst];
		
			cout<<"\t\t for "<<CoordNumbConst::fraction[CoordNumbConst::cum_n_subconst[i]+isubconst]*100<<" % of the central atoms"<<endl; 
			cout<<"\t\t Weight factor: "<<CoordNumbConst::weights[CoordNumbConst::cum_n_subconst[i]+isubconst]<<endl;
		}
		if (CoordNumbConst::write_detail[i])
			cout << "\t Details of the neighbours are written into the *.cncd file" << endl;
	};

	//Average coordination number constraints
	cout << "\n" << J10 <<navcoord<<" 'AVERAGE' COORDINATION CONSTRAINT";
	if (navcoord > 1)
		cout << "S";
	cout << endl;
	for (i=0;i<navcoord;i++)
	{
		cout << "\n" << J4<<i+1<<". constraint: "<<endl;
		//The type has to start with one at the output
		cout<<"\t Central atom type: "<<AvCoordConst::central[i]+1;
		cout<<"\t Neighbour atom type: "<<AvCoordConst::neighbours[i]+1<<endl;
		cout<<"\t Target coordination number: "<<AvCoordConst::acnreq[i]<<endl;
		cout<<"\t From "<<AvCoordConst::dmin[i]<<" A to "<<AvCoordConst::dmax[i]<<" A"<<endl;
		cout<<"\t Weight factor: "<<AvCoordConst::weights[i]<<endl;
	};

#ifdef _ADVANCED_GEOM_CONST

	//Common neighbour constraints
	cout << "\n" << J10 << ncommonneigh << " COMMOM NEIGHBOUR CONSTRAINT";
	if (ncommonneigh > 1)
		cout << "S";
	cout<< endl;
	for (i = 0; i < ncommonneigh; i++)
	{
		cout << "\n" << J4 << i + 1 << ". constraint: " << endl;
		//The type has to start with one at the output
		cout << "\t First primary atom type:    " << CommonNeighConst::primary1[i] + 1;
		cout << "\t Second primary atom type:   " << CommonNeighConst::primary2[i] + 1;
		cout << "\t From " << CommonNeighConst::dmin[i] << " A to " << CommonNeighConst::dmax[i] << " A" << endl;
		for (j = 0; j < CommonNeighConst::nsec_neigh[i]; j++)
		{
			cout << "\t "<<j+1<<". Secondary atom type: " << CommonNeighConst::secondary[CommonNeighConst::cumul_sec[i]+j] + 1 << endl;
			cout << "\t\t Distance range with primary1 from " << CommonNeighConst::dmin1[CommonNeighConst::cumul_sec[i] + j] << " A to " << CommonNeighConst::dmax1[CommonNeighConst::cumul_sec[i] + j] << " A" << endl;
			cout << "\t\t Distance range with primary2 from " << CommonNeighConst::dmin2[CommonNeighConst::cumul_sec[i] + j] << " A to " << CommonNeighConst::dmax2[CommonNeighConst::cumul_sec[i] + j] << " A" << endl;
		}
		cout << "\t Target coordination number: " << CommonNeighConst::target_coord[i] << endl;
		cout << "\t Fraction required:          " << CommonNeighConst::fraction[i] << endl;
		cout << "\t Weight factor:              " << CommonNeighConst::weights[i] << endl;
	}
	//Second neighbour constraints
	cout << "\n" << J10 << nsecondneigh << " SECOND NEIGHBOUR CONSTRAINT";
	if (nsecondneigh > 1)
		cout << "S";
	cout << endl;
	for (i = 0; i < nsecondneigh; i++)
	{
		cout << "\n" << J4 << i + 1 << ". constraint: " << endl;
		//The type has to start with one at the output
		cout << "\t Central atom type        :   " << SecondNeighConst::central[i] + 1;
		cout << "\t First neighbour atom type:   " << SecondNeighConst::fneighbour[i] + 1;
		cout << "\t From " << SecondNeighConst::dmin1[i] << " A to " << SecondNeighConst::dmax1[i] << " A" << endl;
		for (j = 0; j < SecondNeighConst::nsec_type[i]; j++)
		{
			cout << "\t " << j + 1 << ". Second neighbour atom type: " << SecondNeighConst::sneighbours[SecondNeighConst::nsec_type_cum[i] + j] + 1 << endl;
			cout << "\t\t Distance range with first neighbour from " << SecondNeighConst::dmin2[SecondNeighConst::nsec_type_cum[i] + j] << " A to " << SecondNeighConst::dmax2[SecondNeighConst::nsec_type_cum[i] + j] << " A" << endl;
		}
		cout << "\t Target coordination number: " << SecondNeighConst::target_coord[i] << endl;
		cout << "\t Fraction required:          " << SecondNeighConst::fraction[i] << endl;
		cout << "\t Weight factor:              " << SecondNeighConst::weights[i] <<  endl;
	}
	//Bond valence sum constraints
	cout << "\n" << J10 << nbvs << " BOND VALENCE SUM CONSTRAINT";
	if (nbvs > 1)
		cout << "S";
	cout << endl;
	for (i = 0; i < nbvs; i++)
	{
		cout << "\n" << J4 << i + 1 << ". constraint: " << endl;
		//The type has to start with one at the output
		cout << "\t Central atom type:    " << BondValenceSumConst::central[i] + 1;
		if (BondValenceSumConst::central_charge[i] !=NO_IVALUE_DEF)
			cout<< " oxidation state: " << BondValenceSumConst::central_charge[i];
		cout<<endl;
		
		for (ineightype = 0; ineightype < BondValenceSumConst::n_neightype[i]; ineightype++)
		{
			cout << "\t " << ineightype + 1 << ". Neighbour atom type: "<< BondValenceSumConst::neighbours[BondValenceSumConst::cum_n_neightype[i] + ineightype] + 1;
			cout << " to " << BondValenceSumConst::dmax[BondValenceSumConst::cum_n_neightype[i] + ineightype]<<" A "<<endl;
			cout << "\t\tR0: " << BondValenceSumConst::R0[BondValenceSumConst::cum_n_neightype[i] + ineightype] << " A" << endl;
			cout << "\t\tb: " << BondValenceSumConst::b[BondValenceSumConst::cum_n_neightype[i] + ineightype] << " A" << endl;
			if (BondValenceSumConst::neigh_charge[BondValenceSumConst::cum_n_neightype[i] + ineightype] != NO_IVALUE_DEF)
				cout << "\t\toxidation state: " << BondValenceSumConst::neigh_charge[BondValenceSumConst::cum_n_neightype[i] + ineightype] << endl;
		}
		
		cout << "\t Target valence:             " << BondValenceSumConst::target_valence[i] << endl;
		cout << "\t Weight factor:              " << BondValenceSumConst::weights[i] << endl;
	}
#endif
#ifdef _LOCAL_INV
	//Local invariance
	cout << "\n" << "         LOCAL INVARIANCE \n"<<endl;
	
	if (RunParams::loc_chi2_mode==0)
	{
		cout<<"         Local histogram is calculated from "<<min_loc_r<<" to "<<max_loc_r<<" in reduced units"<<endl;
		cout<<"\t corresponding to "<<min_loc_r*boxedge<<" to "<<max_loc_r*boxedge<<" in A."<<endl;
		cout<<"\t Sigma factor: "<<loc_inv_sigma[0]<<endl;
	}
	else
	{
		cout<<"         There are "<<nlocint<<" local invariance intervals, interval borders are :"<<endl;
		for (i=0;i<nlocint;i++)
		{
			cout<<"\t Interval "<<i+1<<": "<<loc_at_ratio[i]<<" to "<<loc_at_ratio[i+1]<<" part of the neighbour atoms"<<endl;
			cout<<"\t\tstarting with index for "<<endl;
			for (j=0;j<ntypes;j++)
				cout<<"\t\t"<<j+1<<".atom type "<<ChiSquared::loc_natoms_min[j*(RunParams::nlocint+1)+i]+1<<endl;
			cout<<"\t Sigma factor: "<<loc_inv_sigma[i]<<endl;
		}
		cout<<"\t This corresponds to the following total number of atoms per type"<<endl;
		for (j=0;j<ntypes;j++)
			cout<<"\t\t"<<j+1<<". type "<<ChiSquared::loc_natoms[j]<<endl;
	}
	cout<<"\t\tThe bin size is "<<loc_rspacing<<" A"<<endl;
#endif

	switch (potential)
	{
		case 0:
		{
			cout<<"\n         POTENTIAL will not be used"<<endl;
			break;
		}
		case 1:
		{
			cout<<"\n\t 6-"<<LJ_rep_N<<" LENNARD-JONES POTENTIAL will be used"<<endl;
			cout<<"\t\t Force field name: " << Topology::force_field_name << endl;
			cout<<"\t\t LJ cutoff:  "<<vdW_cutoff<<" A"<<endl;
			cout<<"\t\t Coulomb cutoff:  "<<Coulomb_cutoff<<" A"<<endl;;
			cout<<"\t\t Scaling factor for 1-4 LJ interactions: "<<vdW14_fudge<<endl;
			cout<<"\t\t Scaling factor for 1-4 Coulomb interactions: "<<Coulomb14_fudge<<endl;
			cout<<"\t\t For each different GROMACS partial:"<<endl;
			cout<<"\t\t "<<J13<<"Atom type1"<<J13<<"Atom type2"<<J16<<"C6 (kj/mol*A^6)"<<J4<<"C"<<LJ_rep_N<<" (kj/mol*A^"<<LJ_rep_N<<endl;
			cout.setf(ios::fixed, ios::floatfield);
			cout.setf(ios::right, ios::adjustfield);

			k=0;
			for (i=0;i<nGRtypes;i++)
			{
				for (j=i;j<nGRtypes;j++)
				{					
					if (vdW_comb_rule==1)
						cout<<"\t\t "<<J13<<Topology::atom_type[i]<<J13<<Topology::atom_type[j]<<J16<<vdW_pot1[k]<<J20<<vdW_pot2[k]<<endl;
					else
						cout<<"\t\t "<<J13<<Topology::atom_type[i]<<J13<<Topology::atom_type[j]<<J16<<4*vdW_pot2[k]*pow(vdW_pot1[k],6.0)<<J20<<4*vdW_pot2[k]*pow(vdW_pot1[k],LJ_rep_N)<<endl;
					k++;
				}
			}
			cout<< "\t\t " << J16 << "   "<<"Weights for the LJ "<<J16<<" Coulomb interaction:"<<endl;
			if (fabs(NB_weight_mode) != 1)
				cout<<"\t\t "<<J16<<"   "<<J16<<vdW_weight[0]<<J16<<Coulomb_weight[0]<<endl;
			else
			{
				for (i=0;i<npartials;i++)
					cout<<"\t\t "<<J2<<i+1<<J16<<". RMC partial   "<<J16<<vdW_weight[i]<<J16<<Coulomb_weight[i]<<endl;
			}

			cout.unsetf(ios::fixed);
			cout.unsetf(ios::right);
			break;
		}
		case 10:
		{
			cout << "\n         TABULATED POTENTIAL will be used for " << nused_potpartials << " partial";
			if (nused_potpartials > 1)
				cout << "s";
			cout<< endl;
			
			cout << "\t\t Temperature: " << temperature << " K "<<endl;
			cout << "\t\t For partials:" << endl;
			for (i = 0; i < nused_potpartials; i++)
			{
				cout << FNC_POT::tabpot_index[i] + 1 << ". partial" << endl;
				cout << "\t\t\t Cutoff: " << FNC_POT::tab_cutoff[i] * SimpleCfg::boxedge << " weight parameter : " << vdW_weight[i] << endl;
			}
			break;
		}
#ifdef _AENET
		case 20:
			cout<<"\n         ANN POTENTIAL calculated by AENET will be used"<<endl;
			cout<<"\t Cutoff: "<<Aenet::cutoff<<" A"<<endl;
			cout<<"\t ANN potential will be caluclated at each "<<aenet_step<<". step"<<endl;
			if (write_E)
				cout<<"\t Energy for each atom will be saved in "<<energyfilename<<endl;
			cout<<"\t Weight parameter : "<< aenet_weight<<endl;
			break;
#endif
	}
		
}

//-----------print----------------
/*void RunParams::Print() const
{//displays on screen
	int i;
	int npartials=ntypes*(ntypes+1)/2;//number of cutoffs

	
	cout<<"This is an object of class RunParams "<<endl;
	cout<<title<<endl;//title
	cout<<cfgnumb<<"\t number of configurations to collect "<<endl;
	cout<<printstep<<"\t step for printing (in number of moves)"<<endl;	
	if (runlimit<0)
		cout<<"Program will run till "<<J20<<-runlimit<<" generated step."<<endl;
	else
		cout<<J10<<runlimit<<"\t minutes time limit for the run."<<endl;
	cout<<timesave<<"\t time interval for saving "<<endl;
	cout<<rho<<"\t number density (atoms per cubic Angstrom)"<<endl;
	for(i=1;i<=npartials;i++)
		cout<<J7<<*(pcutoff+i-1);
	cout<<"\t distances of closest approach "<<endl;
	for(i=1;i<=npartials;i++)
		cout<<J7<<*(pcutsq+i-1);
	cout<<"\t distances of closest approach squared"<<endl;
	for(i=1;i<=ntypes;i++)
		cout<<J7<<*(pmaxmove+i-1);
	cout<<"\t maximum moves "<<endl;
	cout<< rspacing<<"\t width of the histogram bins (r spacing) "<<endl;
	cout<<moveout<<"\t moveout option switch "<<endl;
	cout<<fnc<<"\t fnc option switch "<<endl;
	cout<<rmax<<"\t maximum value for r (eqv. cell size)"<<endl;
	cout<<natoms<<"\t total number of atoms in the cell "<<endl;
	cout<<ntypes<<"\t number of atom types "<<endl;
	for(i=1;i<=ntypes;i++)
		cout<<J7<<*(pnatoms+i-1);
	cout<<"\t number of atoms per type "<<endl;
	cout<<ngr<<"\t number of g(r) constraints "<<endl;
	cout<<nsq<<"\t number of S(Q) - neutron data - constraints "<<endl;
	cout<<nfq<<"\t number of F(Q) - x-ray data - constraints "<<endl;
	cout << nfg << "\t number of F(g) - electron diffraction data - constraints " << endl;
	cout<<nek<<"\t number of EXAFS constraints "<<endl;
	cout<<nicoord<<"\t number of 'individual' coordination constraints "<<endl;
	cout<<navcoord<<"\t number of 'average' coordination constraints "<<endl;
	cout<<potential<<"\t use of potential"<<endl;
}
*/
//saving the sigma values in binary format, important for continuation of the run, if sigma was calculated
//Each data set and constraint types will have an ID, to see, if the data set and constraint types match 
//the actual data sets of the run, when the sigma is read from file. Some additional parameters are saved as well
//IDs:
//	 0: g(r)
//	 1: S(Q)
//	 2: F(Q)
//  14: F(g) it was added later, the original id order was not changed, but it will be written out after F(Q)
//	 3: E(k)
//	 4: cosine distribution of bond angles constraint
//	 5: Coordination number constraint
//	 6: Average coordination constraint
//  16: Common neighbour constraint
//  17: Second neighbour constraint
//  26: Bond valence sum constraint
//	 7: non-bonded potential
//  15: Tabulated potential
//	25: ANN potential
//	 8: bond
//	 9: angle
//	10: periodic dihedral
//	11: harmonic dihedral
//	12: RB dihedral
//	13: local invariance
//---------------------
//	18: ngrid_in for EXAFS E0 shift
//	19: gridind for EXAFS E0 shift
//  20: alpha, a,b for I(Q) sets in case of nonlinreg
//  21: mu for I(Q) sets
//  22: f' for I(Q) set with AXS
//  23: cutoff
//  24: gridfrom, gridto for E0shift state version 3

void RunParams::SaveState(ChiSquared &chi)
{
	int i,id;
	ofstream file;
	char c = '#';

	OpenFile(file,statefilename,"RunParams::SaveState",1);//open file, check, whether it was successfully opened

	//to identify the format of the stat file, DO NOT modify the statefile size without changing the the version in the header
	file.write(reinterpret_cast<const char *>(&c), sizeof(c));
	file.write(reinterpret_cast<const char *>(&stateversion), sizeof(stateversion));
	//writing the generated, tried and accepted moves
	file.write(reinterpret_cast<const char *>(&n_generated), sizeof(n_generated));
	file.write(reinterpret_cast<const char *>(&n_tried), sizeof(n_tried));
	file.write(reinterpret_cast<const char *>(&n_accepted), sizeof(n_accepted));
	file.write(reinterpret_cast<const char *>(&n_potaccepted), sizeof(n_potaccepted));
	file.write(reinterpret_cast<const char *>(&n_swapgen), sizeof(n_swapgen));
	file.write(reinterpret_cast<const char *>(&n_swapacc), sizeof(n_swapacc));
	file.write(reinterpret_cast<const char *>(&n_E0shiftgen), sizeof(n_E0shiftgen));
	file.write(reinterpret_cast<const char *>(&n_E0shiftacc), sizeof(n_E0shiftacc));
	file.write(reinterpret_cast<const char *>(&n_mucorrgen), sizeof(n_mucorrgen));
	file.write(reinterpret_cast<const char *>(&n_mucorracc), sizeof(n_mucorracc));
	file.write(reinterpret_cast<const char *>(&n_fprimeshiftgen), sizeof(n_fprimeshiftgen));
	file.write(reinterpret_cast<const char *>(&n_fprimeshiftacc), sizeof(n_fprimeshiftacc));
	file.write(reinterpret_cast<const char *>(&n_goodnonlinreg), sizeof(n_goodnonlinreg));
	

	//the g(r) components
	id=0;
	for(i=0;i<ngr;i++)//for each g(r) data set
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&ExptsData::grsigma[i]), sizeof(*ExptsData::grsigma));
	}
	//the S(Q) components
	id=1;
	for(i=0;i<nsq;i++)//for each S(Q) data set
	{

		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&ExptsData::sqsigma[i]), sizeof(*ExptsData::sqsigma));
	}
	//the F(Q) components
	id=2;
	for(i=0;i<nfq;i++)//for each F(Q) data set
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&ExptsData::fqsigma[i]), sizeof(*ExptsData::fqsigma));
	}
	//the F(g) components it was added later to the code, the original ids will not be changed
	//but it will be written out after the X-ray data
	id = 14;
	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&ExptsData::fgsigma[i]), sizeof(*ExptsData::fgsigma));
	}
	//the E(k) components
	id=3;
	for(i=0;i<nek;i++)//for each E(k) data set
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&ExptsData::eksigma[i]), sizeof(*ExptsData::eksigma));
	}
	
	//the Cosine disrtribution of bond angles constraint
	id=4;
	for(i=0;i<ncosdistr;i++)//for each CosDistr constraint
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&CosDistrConst::weights[i]), sizeof(*CosDistrConst::weights));
	}
	//the Coordination number constraint
	id=5;
	for(i=0;i<nicoord;i++)//for each Coordination number constraint
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&CoordNumbConst::weights[i]), sizeof(*CoordNumbConst::weights));
	}
	//the Average coordination constraint
	id=6;
	for(i=0;i<navcoord;i++)//for each Average coordination constraint
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&AvCoordConst::weights[i]), sizeof(*AvCoordConst::weights));
	}
#ifdef _ADVANCED_GEOM_CONST
	//the Common neighbour constraint
	id = 16;
	for (i = 0; i < ncommonneigh; i++)//for each Common neighbour constraint
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&CommonNeighConst::weights[i]), sizeof(*CommonNeighConst::weights));
	}
	//the Second neighbour constraint
	id = 17;
	for (i = 0; i < nsecondneigh; i++)//for each Second neighbour constraint
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&SecondNeighConst::weights[i]), sizeof(*SecondNeighConst::weights));
	}
	//the Bond valence sumr constraint
	id = 26;
	for (i = 0; i < nbvs; i++)//for each Bond Valence sum constraint
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&BondValenceSumConst::weights[i]), sizeof(*BondValenceSumConst::weights));
	}
#endif
	//LJ potential
	id=7;
	if (potential==1)
	{
		for (i=0;i<SimpleCfg::npartials;i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&vdW_weight[i]), sizeof(*vdW_weight));
			file.write(reinterpret_cast<const char *>(&Coulomb_weight[i]), sizeof(*Coulomb_weight));
		}
	}
	//Tabulated potential
	id = 15;
	if (potential ==10)
	{
		for (i = 0; i < nused_potpartials; i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&vdW_weight[i]), sizeof(*vdW_weight));
		}
	}
#ifdef _AENET
	//ANN potential
	id = 25;
	if (potential ==20)
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&aenet_weight), sizeof(aenet_weight));
	}
#endif
	//bonds
	id=8;
	if (nbond_types>0)
	{
		for (i=0;i<nbond_types;i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&Topology::bond_sigma[i]), sizeof(*Topology::bond_sigma));
		}
	}
	//angles
	id=9;
	if (nangle_types>0)
	{
		for (i=0;i<nangle_types;i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&Topology::angle_sigma[i]), sizeof(*Topology::angle_sigma));
		}
	}
	//periodic dihedrals
	id=10;
	if (nperdihedral_types>0)
	{
		for (i=0;i<nperdihedral_types;i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&Topology::perdihedral_sigma[i]), sizeof(*Topology::perdihedral_sigma));
		}
	}
	//periodic dihedrals
	id=11;
	if (nharmdihedral_types>0)
	{
		for (i=0;i<nharmdihedral_types;i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&Topology::harmdihedral_sigma[i]), sizeof(*Topology::harmdihedral_sigma));
		}
	}
	//RB dihedrals
	id=12;
	if (nRBdihedral_types>0)
	{
		for (i=0;i<nRBdihedral_types;i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&Topology::RBdihedral_sigma[i]), sizeof(*Topology::RBdihedral_sigma));
		}
	}
#ifdef _LOCAL_INV
	id=13;
	for (i=0;i<nlocint;i++)
	{
		file.write(reinterpret_cast<const char *>(&id), sizeof(id));
		file.write(reinterpret_cast<const char *>(&RunParams::loc_inv_sigma[i]), sizeof(*RunParams::loc_inv_sigma));
	}

#endif
	if (ExptsData::is_E0shift)//from version 1.8
	{
		//the number of grid points in one direction for EXAFS E0 shift grid and the last grid index is saved
		for (i = 0; i < nek; i++)
		{
			id = 18;
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&ExptsData::ek_ngrid_in[i]), sizeof(*ExptsData::ek_ngrid_in));
			id = 19;
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&ExptsData::ek_gridind[i]), sizeof(*ExptsData::ek_gridind));
			id = 24;
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&ExptsData::ek_gridfrom[i]), sizeof(*ExptsData::ek_gridfrom));
			file.write(reinterpret_cast<const char *>(&ExptsData::ek_gridto[i]), sizeof(*ExptsData::ek_gridto));
		}
	}
	if (ExptsData::is_IQ)
	{
		id = 20;
		for (i = 0; i < nfq; i++)
		{
	
			if (ExptsData::fqrenalpha[i])
			{
				file.write(reinterpret_cast<const char *>(&id), sizeof(id));
				file.write(reinterpret_cast<const char *>(&chi.alpha[i]), sizeof(*chi.alpha));
				file.write(reinterpret_cast<const char *>(&chi.a[i + ngr + nsq]), sizeof(*chi.a));
				file.write(reinterpret_cast<const char *>(&chi.b[i + ngr + nsq]), sizeof(*chi.b));
			}
		}
	}
	if (ExptsData::is_IQmucorr)
	{
		id = 21;
		for (i = 0; i < nfq; i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&ExptsData::fqmuact[i]), sizeof(*ExptsData::fqmuact));
		}
	}
	if (ExptsData::is_AXS)
	{
		id = 22;
		for (i = 0; i < nfq; i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&ExptsData::fqfprimeact[i*ntypes+ExptsData::fqAXS[i]-1]), sizeof(*ExptsData::fqfprimeact));
		}
	}
	if (RunParams::auto_cutoff !=0)
	{
		id = 23;
		for (i = 0; i < SimpleCfg::npartials; i++)
		{
			file.write(reinterpret_cast<const char *>(&id), sizeof(id));
			file.write(reinterpret_cast<const char *>(&RunParams::pcutoff[i]), sizeof(*RunParams::pcutoff));
		}
	}
	file.close();
}

void RunParams::LoadState()
{
	int i, id, temp_i;
	longint temp_l;
	char *datatype;
	char temp_c;

	SetArraysize(&datatype, 70, "datatype", "RunParams::LoadState");


	ifstream file;

	CleanOpen(file, statefilename, 1);//open file
	if (CheckFileState(file, "RunParams::LoadState", statefilename) == 0)//check, whether it was successfully opened
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "There is no open " << statefilename << " file to load the sigma and chi2 values" << endl;
		cout << "Exact continuation of the run is not possible, exiting!" << endl;
		CleanExit();//loading was not successful
	}
	//checking the version of the state file, only the same version can be loaded
	ReadBin(file, &temp_c, 1, "state_file_version", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_c != '#')
	{
		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! This state file is not compatible with this RMC++ version," << endl;
		cout << "excact continuation is not possible! Start the program without the cont option" << endl;
	}
	ReadBin(file, &temp_i, 1, "stateversion", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_i != stateversion)
	{
		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! This state file is " << temp_i << ", which is not compatible with the state file ";
		cout << "version " << stateversion << " of this RMC++ code!" << endl;
		cout << "Exact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	//reading the generated, tried and accepted ....moves
	ReadBin(file, &temp_l, 1, "n_generated", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_generated)
	{
		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of generated moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[0] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_tried", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_tried)
	{
		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of tried moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[1] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_accepted", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_accepted)
	{
		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of accepted moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[2] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_potaccepted", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_potaccepted)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of accepted moves based on the potential is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[3] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_swap_gen", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_swapgen)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of swap generated moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[4] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_swapacc", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_swapacc)
	{
		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of swap accepted moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[5] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_E0shiftgen", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_E0shiftgen)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of E0 shift generated moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[6] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_E0shiftacc", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_E0shiftacc)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of E0 shift accepted moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[7] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_mucorrgen", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_mucorrgen)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of I(Q) mu correction generated moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[8] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_mucorracc", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_mucorracc)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of I(Q) mu correction accepted moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[9] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_fprimeshiftgen", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_fprimeshiftgen)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of f' shift generated moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[10] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_fprimeshiftacc", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_fprimeshiftacc)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of f' shift accepted moves is " << temp_l << " from the " << statefilename << " file," << endl;
		cout << "and " << *SimpleCfg::nmoves[11] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}
	ReadBin(file, &temp_l, 1, "n_goodnonlinreg", "SimpleCfg::LoadBinary", statefilename, 2);
	if (temp_l != n_goodnonlinreg)
	{

		cout << "\n****ERROR****" << endl;
		cout << "INCONSISTENCY! The number of tried steps ended with successful nonlinear regression " << temp_l << " from the " << endl;
		cout << statefilename << " file," << "and " << *SimpleCfg::nmoves[12] << " from the configuration file!" << endl;
		cout << "Excact continuation is not possible! Start the program without the cont option" << endl;
		CleanExit();
	}

	//the g(r) components
	id = 0;
	for (i = 0; i < ngr; i++)//for each g(r) data set
	{
		ReadBin(file, &temp_i, 1, "gr_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". g(r) series a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &ExptsData::grsigma[i], 1, "ExptsData::grsigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the S(Q) components
	id = 1;
	for (i = 0; i < nsq; i++)//for each S(Q) data set
	{
		ReadBin(file, &temp_i, 1, "sq_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". S(Q) series a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &ExptsData::sqsigma[i], 1, "ExptsData::sqsigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the F(Q) components
	id = 2;
	for (i = 0; i < nfq; i++)//for each F(Q) data set
	{
		ReadBin(file, &temp_i, 1, "fq_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". F(Q) series a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &ExptsData::fqsigma[i], 1, "ExptsData::fqsigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the F(g) components
	id = 14;
	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		ReadBin(file, &temp_i, 1, "fg_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". F(g) series a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &ExptsData::fgsigma[i], 1, "ExptsData::fgsigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the E(k) components
	id = 3;
	for (i = 0; i < nek; i++)//for each E(k) data set
	{
		ReadBin(file, &temp_i, 1, "exafs_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". E(k) series a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &ExptsData::eksigma[i], 1, "ExptsData::eksigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}

	//the Cosine distribution of bond angles constraint
	id = 4;
	for (i = 0; i < ncosdistr; i++)//for each CosDistr data set
	{
		ReadBin(file, &temp_i, 1, "cos_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". cosine distribution of bond angles constraint a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &CosDistrConst::weights[i], 1, "CosDistrConst::weights", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the Coordination number constraint
	id = 5;
	for (i = 0; i < nicoord; i++)//for each Cooord numb constraint
	{
		ReadBin(file, &temp_i, 1, "cnc_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". coordination number constraint a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &CoordNumbConst::weights[i], 1, "CoordNumbConst::weights", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the Average coordination constraint
	id = 6;
	for (i = 0; i < navcoord; i++)//for each Avreage coord constraint
	{
		ReadBin(file, &temp_i, 1, "avcnc_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". average coordination constraint a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &AvCoordConst::weights[i], 1, "AvCoordConst::weights", "SimpleCfg::LoadBinary", statefilename, 2);
	}
#ifdef _ADVANCED_GEOM_CONST
	//the Common neighbour constraint
	id = 16;
	for (i = 0; i < ncommonneigh; i++)//for each Common neighbour constraint
	{
		ReadBin(file, &temp_i, 1, "conc_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". common neighbour constraint a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &CommonNeighConst::weights[i], 1, "CommonNeighConst::weights", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the Second neighbour constraint
	id = 17;
	for (i = 0; i < nsecondneigh; i++)//for each Second neighbour constraint
	{
		ReadBin(file, &temp_i, 1, "snc_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". second neighbour constraint a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &SecondNeighConst::weights[i], 1, "SecondNeighConst::weights", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the Bond valence sum constraint
	id = 26;
	for (i = 0; i < nbvs; i++)//for each Bond valence sum constraint
	{
		ReadBin(file, &temp_i, 1, "bvs_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". bond valence sum constraint a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &BondValenceSumConst::weights[i], 1, "BondValenceSumConst::weights", "SimpleCfg::LoadBinary", statefilename, 2);
	}
#endif
	//the LJ potential
	id = 7;
	if (RunParams::potential == 1)
	{
		for (i = 0; i < SimpleCfg::npartials; i++)//for each partial
		{
			ReadBin(file, &temp_i, 1, "LJ_id", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". potential partial a " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				CleanExit();
			}
			ReadBin(file, &vdW_weight[i], 1, "vdW_weight", "SimpleCfg::LoadBinary", statefilename, 2); 
			ReadBin(file, &Coulomb_weight[i], 1, "Coulomb_weight", "SimpleCfg::LoadBinary", statefilename, 2);
		}
	}
	//the tabulated potential
	id = 15;
	if (RunParams::potential == 10)
	{
		for (i = 0; i < nused_potpartials; i++)//for each used partial
		{
			ReadBin(file, &temp_i, 1, "tabpot_id", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". potential partial a " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				CleanExit();
			}
			ReadBin(file, &vdW_weight[i], 1, "vdW_weight", "SimpleCfg::LoadBinary", statefilename, 2);
		}
	}
#ifdef _AENET
	//the ANN potential
	id = 25;
	if (RunParams::potential == 20)
	{
		ReadBin(file, &temp_i, 1, "aenet_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the ANN potential a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &aenet_weight, 1, "aenet_weight", "SimpleCfg::LoadBinary", statefilename, 2);
	}
#endif
	//the bonds
	id = 8;
	for (i = 0; i < nbond_types; i++)//for each bond type
	{
		ReadBin(file, &temp_i, 1, "bond_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". bond type a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &Topology::bond_sigma[i], 1, "Topology::bond_sigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the angles
	id = 9;
	for (i = 0; i < nangle_types; i++)//for each angle type
	{
		ReadBin(file, &temp_i, 1, "angle_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". angle type a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &Topology::angle_sigma[i], 1, "Topology::angle_sigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the perdihedral
	id = 10;
	for (i = 0; i < nperdihedral_types; i++)//for each perdihedral type
	{
		ReadBin(file, &temp_i, 1, "perdih_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". periodic dihedral type a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &Topology::perdihedral_sigma[i], 1, "Topology::perdihedral_sigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the harmdihedral
	id = 11;
	for (i = 0; i < nharmdihedral_types; i++)//for each harmdihedral type
	{
		ReadBin(file, &temp_i, 1, "harmdih_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". harmonic dihedral type a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &Topology::harmdihedral_sigma[i], 1, "Topology::harmdihedral_sigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}
	//the RBdihedral
	id = 12;
	for (i = 0; i < nRBdihedral_types; i++)//for each RBdihedral type
	{
		ReadBin(file, &temp_i, 1, "RBdih_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". RB dihedral type a " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &Topology::RBdihedral_sigma[i], 1, "Topology::RBdihedral_sigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}

#ifdef _LOCAL_INV
	id = 13;
	for (i = 0; i < nlocint; i++)
	{
		ReadBin(file, &temp_i, 1, "locinv_id", "SimpleCfg::LoadBinary", statefilename, 2);
		if (temp_i != id)
		{
			GetDataType(temp_i, datatype);//getting the data type from the file
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! Instead of the " << i + 1 << ". local invariance interval type " << datatype << " can be " << endl;
			cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
			cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
			CleanExit();
		}
		ReadBin(file, &RunParams::loc_inv_sigma[i], 1, "RunParams::loc_inv_sigma", "SimpleCfg::LoadBinary", statefilename, 2);
	}

#endif
	//----------------------- this is not sigma, and only handled from version 1.8
	if (ExptsData::is_E0shift)
	{
		for (i = 0; i < nek; i++)
		{
			ReadBin(file, &temp_i, 1, "E0shift_gridpoints_id", "SimpleCfg::LoadBinary", statefilename, 2);
			id = 18;
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". EXAFS E0 shift ngrid type " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				cout << "Exiting..." << endl;
				CleanExit();
			}
			ReadBin(file, &temp_i, 1, "ExptsData::ek_ngrid_in", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != ExptsData::ek_ngrid_in[i])
			{
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! The number of EXAFS shift grid points in one direction is " << ExptsData::ek_ngrid_in[i] << " from the " << datfilename << " file" << endl;
				cout << "and " << temp_i << " from the " << statefilename << " file for the  " << i + 1 << ". EXAFS data!" << endl;
				cout << "Exact continuation is not possible, exiting..." << endl;
				CleanExit();
			}
			id = 19;
			ReadBin(file, &temp_i, 1, "E0shift_gridind_id", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". EXAFS E0 shift grid index type " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				CleanExit();
			}
			ReadBin(file, &ExptsData::ek_gridind[i], 1, "ExptsData::ek_gridind", "SimpleCfg::LoadBinary", statefilename, 2);
			id = 24;
			ReadBin(file, &temp_i, 1, "E0shift_gridrange_id", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". EXAFS E0 shift grid from - to indices type " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				CleanExit();
			}
			ReadBin(file, &temp_i, 1, "ExptsData::ek_gridfrom", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != ExptsData::ek_gridfrom[i])
			{
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! The grid index of the beginning (FROM) of the used EXAFS shift grid point range is " << ExptsData::ek_gridfrom[i] << " from the " << datfilename << " file" << endl;
				cout << "and " << temp_i << " from the " << statefilename << " file for the  " << i + 1 << ". EXAFS data!" << endl;
				cout << "Exact continuation is not possible, exiting..." << endl;
				CleanExit();
			}
			ReadBin(file, &temp_i, 1, "ExptsData::ek_gridto", "SimpleCfg::LoadBinary", statefilename, 2);
			if (temp_i != ExptsData::ek_gridto[i])
			{
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! The The grid index of the end (TO) of the used EXAFS shift grid point range is " << ExptsData::ek_gridto[i] << " from the " << datfilename << " file" << endl;
				cout << "and " << temp_i << " from the " << statefilename << " file for the  " << i + 1 << ". EXAFS data!" << endl;
				cout << "Exact continuation is not possible, exiting..." << endl;
				CleanExit();
			}
			ExptsData::ek_gridstart[i] = ExptsData::ek_gridind[i]+ExptsData::ek_gridfrom[i];
		}
	}
	if (ExptsData::is_IQ)
	{
		for (i = 0; i < nfq; i++)
		{
			if (ExptsData::fqrenalpha[i])
			{
				ReadBin(file, &temp_i, 1, "IQfit_id", "SimpleCfg::LoadBinary", statefilename, 2);
				id = 20;
				if (temp_i != id)
				{
					GetDataType(temp_i, datatype);//getting the data type from the file
					cout << "\n****ERROR****" << endl;
					cout << "INCONSISTENCY! Instead of the " << i + 1 << ". I(Q) a, b and alpha value data type " << datatype << " can be " << endl;
					cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
					cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
					cout << "Exiting..." << endl;
					CleanExit();
				}
				ReadBin(file, &ExptsData::fqalpha[i], 1, "ExptsData::fqalpha", "SimpleCfg::LoadBinary", statefilename, 2);
				ReadBin(file, &ExptsData::fqa[i], 1, "ExptsData::fqa", "SimpleCfg::LoadBinary", statefilename, 2);
				ReadBin(file, &ExptsData::fqb[i], 1, "ExptsData::fqb", "SimpleCfg::LoadBinary", statefilename, 2);
			}
		}
	}

	if (ExptsData::is_IQmucorr)
	{
		for (i = 0; i < nfq; i++)
		{
			ReadBin(file, &temp_i, 1, "IQmucorr_id", "SimpleCfg::LoadBinary", statefilename, 2);
			id = 21;
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". I(Q) actual mu value data type " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				cout << "Exiting..." << endl;
				CleanExit();
			}
			ReadBin(file, &ExptsData::fqmuact[i], 1, "ExptsData::fqmuact", "SimpleCfg::LoadBinary", statefilename, 2);
		}
	}
	if (ExptsData::is_AXS)
	{
		for (i = 0; i < nfq; i++)
		{
			ReadBin(file, &temp_i, 1, "IQAXS_id", "SimpleCfg::LoadBinary", statefilename, 2);
			id = 22;
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". anomalous I(Q)-fitting X-ray data set's actual f' value data type " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				cout << "Exiting..." << endl;
				CleanExit();
			}
			ReadBin(file, &ExptsData::fqfprimeact[i * ntypes + ExptsData::fqAXS[i] - 1], 1, "ExptsData::fqfprimeact", "SimpleCfg::LoadBinary", statefilename, 2);
		}
	}
	if (abs(RunParams::auto_cutoff) != 0)
	{
		for (i = 0; i < SimpleCfg::npartials; i++)
		{
			ReadBin(file, &temp_i, 1, "cutoff_id", "SimpleCfg::LoadBinary", statefilename, 2);
			id = 23;
			if (temp_i != id)
			{
				GetDataType(temp_i, datatype);//getting the data type from the file
				cout << "\n****ERROR****" << endl;
				cout << "INCONSISTENCY! Instead of the " << i + 1 << ". cutoff value data type " << datatype << " can be " << endl;
				cout << "found in its place in the " << statefilename << " file! The " << datfilename << " file is not " << endl;
				cout << "consistent with the " << statefilename << " file, cannot perform exact continuation!" << endl;
				cout << "Exiting..." << endl;
				CleanExit();
			}
			ReadBin(file, &RunParams::pcutoff[i], 1, "RunParams::pcutoff[i]", "SimpleCfg::LoadBinary", statefilename, 2);
		}
	}
	if (!CheckReadFileState(file,"RunParams::LoadState",statefilename))
	{
		cout<<"Loading of the "<<statefilename<<" file has failed!"<<endl;
		cout<<"Exact continuation cannot be performed, start again without the continuation option!"<<endl;
		CleanExit();
	}
	file.close();
}

void RunParams::GetDataType(int i, const char *data_type)
{
	switch (i)
	{
		case 0:
		{
			data_type="g(r) sigma";
			break;
		}
		case 1:
		{
			data_type="S(Q) sigma";
			break;
		}
		case 2:
		{
			data_type="F(Q) sigma";
			break;
		}
		case 3:
		{
			data_type="E(k) sigma";
			break;
		}
		case 4:
		{
			data_type="Cosine distribution of bond angles constraint sigma";
			break;
		}

		case 5:
		{
			data_type="coordination number constraint sigma";
			break;
		}
		case 6:
		{
			data_type="averaga coordination constraint sigma";
			break;
		}
		case 7:
		{
			data_type="non-bonded potential sigma";
			break;
		}
		case 8:
		{
			data_type="bond sigma";
			break;
		}
		case 9:
		{
			data_type="angle sigma";
			break;
		}
		case 10:
		{
			data_type="periodic dihedral sigma";
			break;
		}
		case 11:
		{
			data_type="harmionic dihedral sigma";
			break;
		}
		case 12:
		{
			data_type="RB dihedral sigma";
			break;
		}
		case 13:
		{
			data_type="local invariance sigma";
			break;
		}
		case 14:
		{
			data_type = "F(g) sigma";
			break;
		}
		case 15:
		{
			data_type = "tabulated potential sigma";
			break;
		}
		case 16:
		{
			data_type = "common neighbour constraint sigma";
			break;
		}
		case 17:
		{
			data_type = "second neighbour constraint sigma";
			break;
		}
		case 18:
		{
			data_type = "EXAFS E0 shift number of grid points in one direction";
			break;
		}
		case 19:
		{
			data_type = "EXAFS E0 shift grid index";
			break;
		}
		case 20:
		{
			data_type = "I(Q) a, b, alpha for non-linear regression";
			break;
		}
		case 21:
		{
			data_type = "I(Q) actual mu";
			break;
		}
		case 22:
		{
			data_type = "AXS I(Q) actual f'";
			break;
		}
		case 23:
		{
			data_type = "cutoff'";
			break;
		}
		case 24:
		{
			data_type = "gridfrom, gridto for E0shift'";
			break;
		}
		case 25:
		{
			data_type = "ANN potential sigma'";
			break;
		}

	}
}

//setting the parameters for the bonded interactions
void RunParams::SetBondedParams()
{
	nbond_types=Topology::nbond_types;
	nangle_types=Topology::nangle_types;
	nperdihedral_types=Topology::nperdihedral_types;
	nharmdihedral_types=Topology::nharmdihedral_types;
	nRBdihedral_types=Topology::nRBdihedral_types;
};

//giving error message in case of the leading potential series index is on a bonded interaction with zero sigma
void RunParams::NoLeadZeroSigma(int i, const char *name)
{
	cout<<"\n*****ERROR*****"<<endl;
	cout<<"The leading potential series index is "<<lead_series_ind2<<", and refers to the "<<i+1<<". "<<name<<endl;
	cout<<"interaction, which has zero sigma. Choose an other leading series with non-zero sigma!"<<endl;
	cout<<"Cannot run this way, exiting..."<<endl;
	CleanExit();
};


//Check whether the given type read from input is in the range
void RunParams::CheckType(int type, char *name, const char *routine_name)
{
	if (type > SimpleCfg::ntypes)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The type for " << name << " in routine " << routine_name << " is " << type << ", which is larger than " << SimpleCfg::ntypes <<","<< endl;
		cout << "the number of types given in " << cfgfilename << " file! " << endl;
		cout << "Cannot run this way, exiting...." << endl;
		CleanExit();
	}
};


void RunParams::GetFreeFormatParams(ifstream &file)
{
	string line;
	bool usedebug=false;//as the keys are not processed yet, if DEBUG is given it will be assumed, that it is on
	bool continuesw = true;
	using std::list;
	std::list<string> gener;//Stores general string
	std::list<std::list<string> > exper;//Experiment strings
	std::list<std::list<string> > cocon;//Coordination constraint strings
	std::list<std::list<string> > avcocon;//Average coordination constraint strings
	std::list<std::list<string> > coscon;//Cosine distribution constraint strings
	//The data for options only available if the program complied with the proper compiler option will be read if present, but it will only be processed, if the program was compiled for it 
	std::list<std::list<string> > conc;//Common neighbour constraint strings
	std::list<std::list<string> > snc;//Second neighbour constraint strings
	std::list<std::list<string> > bvs;//BondValenceSum constraint strings

	std::list<string>  bpot;//Stores bonded potential related strings
	std::list<string> swp;//Stores swap-moves related strings

	std::list<string> loc;//Local invariance strings

	std::list<string> vibramp;//Vibrational amplitude  strings

	std::list<string> noper;//no periodic boundary conditions strings

	std::list<string> nbpot;//Stores non-bonded potential related strings
	std::list<string> aenet;//stores anete potential related strings
	std::list<string> *apoi = 0;//Pointer to string structures

	//Read keywords in, categorizing and wrapping each line
	bool fst = true;//To avoid reading any lines before the first tag item
	bool conn = false;//To avoid qualifying the first non-keyword containing wrong line after tag item as "continue line"
	fnc = FNC_DEF;//fnc can be set in this routine, if BPOT tag is read

	while (continuesw&&getline(file, line))
	{
		//Stripping the comments and capitalize
		std::size_t le = line.find_first_of('!');
		if (le != string::npos) line = line.substr(0, le);
		line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
		line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
		string origline(line);
		std::transform(line.begin(), line.end(), line.begin(), [](char ch) {if ((ch == '\r') || (ch == '\t')) ch = ' '; return ch; });
		std::transform(line.begin(), line.end(), line.begin(), ::toupper);
		if (line.rfind("DEBUG") != std::string::npos)
			usedebug = true;
		if ((line.find('[') != string::npos) && (line.find(']') != string::npos))//tags
		{
			//A tag item has been found remove blank spaces and the out []
			fst = false;
			conn = false;
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			line = line.substr(line.find_first_of('[') + 1, line.find_last_of(']') - line.find_first_of('[') - 1);
			if (line == CheckTag("END")) continuesw = false;
			else if (line == CheckTag("GENERAL")) apoi = &gener;
			else if (line == CheckTag("EXP"))
			{
				std::list<string> newexp;
				exper.push_back(newexp);
				apoi = &exper.back();
			}
			else if (line == CheckTag("COS"))
			{
				std::list<string> newexp;
				coscon.push_back(newexp);
				apoi = &coscon.back();
			}
			else if (line == CheckTag("COORD"))
			{
				std::list<string> newexp;
				cocon.push_back(newexp);
				apoi = &cocon.back();
			}
			else if (line == CheckTag("AVCOORD"))
			{
				std::list<string> newexp;
				avcocon.push_back(newexp);
				apoi = &avcocon.back();
			}
			else if (line == CheckTag("CONC"))
			{
#ifndef _ADVANCED_GEOM_CONST
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("CONC") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _ADVANCED_GEOM_CONST compiler option!" << endl;
#endif
				std::list<string> newexp;
				conc.push_back(newexp);
				apoi = &conc.back();
			}
			else if (line == CheckTag("SNC"))
			{
#ifndef _ADVANCED_GEOM_CONST
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("SNC") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _ADVANCED_GEOM_CONST compiler option!" << endl;
#endif
				std::list<string> newexp;
				snc.push_back(newexp);
				apoi = &snc.back();
			}
			else if (line == CheckTag("BVS"))
			{
#ifndef _ADVANCED_GEOM_CONST
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("BVS") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _ADVANCED_GEOM_CONST compiler option!" << endl;
#endif
				std::list<string> newexp;
				bvs.push_back(newexp);
				apoi = &bvs.back();
			}
			else if (line == CheckTag("VIBRAMP"))
			{
#ifndef _VIBR_AMP
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("VIBRAMP") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _VIBR_AMP compiler option!" << endl;
#endif
				apoi = &vibramp;
			}
			else if (line == CheckTag("NOPER"))
			{
#ifndef _NO_PERIODIC	
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("NOPER") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _NO_PERIODIC compiler option!" << endl;
#endif
				apoi = &noper;
			}
			else if (line == CheckTag("SWAP")) apoi = &swp;
			else if (line == CheckTag("LOCINV"))
			{
#ifndef _LOCAL_INV
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("LOCINV") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _LOCAL_INV compiler option!" << endl;
#endif
				apoi = &loc;
			}
			else if (line == CheckTag("CUSTMOVE")) apoi = &cus_entry;
			else if (line == CheckTag("NBPOT"))	apoi = &nbpot;
			//else if ((line == "[NB]") && ((apoi == &pothead) || (apoi == &nbpot) || HaveReference(apoi, bpot))) apoi = &nbpot;
			else if (line == CheckTag("AENET"))
			{
#ifndef _AENET
				cout << "\nWARNING(" << ++warn << "): [ " << CheckTag("AENET") << " ] tag was found, but it will be ignored, as the program was not compiled with" << endl;
				cout << "\tthe _AENET compiler option!" << endl;
#endif
				apoi = &aenet;
			}
			else if (line == CheckTag("BPOT"))
			{
				apoi = &bpot;
				fnc = -4;//it is possible that there will be no keywords in this section, the tag in itself indicate, that flexible molecules should be used.
				//the fnc switch can be set in the general section as well, so this indicate, that BPOT tag was found. It will be set to 4 in the beginning of the ReadFreeBpot
			}
			else
			{
				cout << "\nERROR: The \"" << line << "\" tag is invalid!" << endl;
				CleanExit();
			}
		}
		else if ((!fst) && (line.length() > 0) && (line.find_first_not_of(' ') != string::npos))
		{
			if  ((std::find_if(line.begin(), line.end(), [](char c) { return std::isalpha(c); }) != line.end()) && (!((line.find('=') == string::npos) && ((apoi == &nbpot) ))))
			{
				if (line.find('=') == string::npos)
				{
					//Missing '=' sign after the keyword -- ERROR
					std::cerr << "\nERROR: Missing \"=\" sign in the following line: " << endl << "\"" << line << "\"" << endl;
					file.close();
					CleanExit();
				}
				else
				{
					
					//Add the result to a new line
					if ((line.find(CheckKey("TITLE")) != string::npos) || (line.find("FILE") != string::npos) || (line.find(CheckKey("CUSTOM-SFACTORTABLE")) != string::npos) || (line.find(CheckKey("CHEMICAL-SYMBOLS")) != string::npos) || (line.find(CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS")) != string::npos) || ((line.find("UNIFORM") == string::npos) && (line.find("GAUSSIAN") == string::npos) && (line.find("ABSENT") == string::npos) && (line.find("DISTRIB-TYPE") != string::npos)))
					{
						Ltrim(origline);
						std::transform(origline.begin(), origline.begin()+origline.find_first_of('='), origline.begin(), ::toupper);
						apoi->push_back(origline);
					}
					else
					{
						Ltrim(line);
						if ((line.find(CheckKey("LOC-MODE")) != string::npos)|| (line.find(CheckKey("NGR-TYPES")) != string::npos))
							apoi->push_front(line);
						else
							apoi->push_back(line);
					}
						
					conn = true;
				}
			}
			else if (conn)
			{
				//This line might contain additional data continuing the last. Append it
				string base = apoi->back() + line;
				apoi->pop_back();
				apoi->push_back(base);
			}
			else cout << "\nWARNING(" << ++warn << "): The following line does not contain any useful information: " << endl << "\"" << line << "\"" << endl;
		}
	}
#ifdef _ADVANCED_GEOM_CONST
	if (( !bpot.empty()) && !bvs.empty())//BPOT and BVS is specified, cannot use both
	{
		cout << "\nERROR:  Both BPOT (using bonding potential) and BVS (Bond valence sum constraint) " << endl;
		cout << "\tis specified, which is not allowed in the same time!" << endl;
		cout << "\tCorrect the " << datfilename << " file, and try again!" << endl;
		CleanExit();
	}
#endif
#ifdef _AENET
	if ((!nbpot.empty() || !bpot.empty()) && !aenet.empty())//NBPOT and AENET is specified, cannot use both
	{
		cout << "\nWARNING(" << ++warn << "): NBPOT and/or BPOT (using classical potential) and AENET "<<endl;
		cout<<"\t(using artificial neural network potential) is specified, which is not allowed in the same time!"<<endl;
		cout<<"\tAs the code is compiled explicitly for AENET usage, AENET will be used, and NBPOT/BPOT is ignored!"<<endl;
		cout<<"Continuing, but correct the "<<datfilename <<" file!"<<endl;
		bpot.clear();
		nbpot.clear();
		fnc=0;
	} 	
#else
	if (!aenet.empty())
		aenet.clear();
#endif
	if (::debug ||usedebug)
	{
		if (!gener.empty())
			cout << "\nGENERALS:" << endl;
		for (std::list<string>::iterator it = gener.begin(); it != gener.end(); it++) cout << " \"" << *it << "\"" << endl;
		if (!exper.empty())
			cout << "\nEXPERIMENTS:" << endl;
		int num = 0;
		for (std::list<std::list<string> >::iterator it2 = exper.begin(); it2 != exper.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
		if (!cocon.empty())
			cout << "\nCOORDINATION CONSTRAINTS:" << endl;
		num = 0;
		for (std::list<std::list<string> >::iterator it2 = cocon.begin(); it2 != cocon.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
		if (!avcocon.empty())
			cout << "\nAVERAGE COORDINATION CONSTRAINTS:" << endl;
		num = 0;
		for (std::list<std::list<string> >::iterator it2 = avcocon.begin(); it2 != avcocon.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
		if (!coscon.empty())
			cout << "\nCOSINE DISTRIBUTION CONSTRAINTS:" << endl;
		num = 0;
		for (std::list<std::list<string> >::iterator it2 = coscon.begin(); it2 != coscon.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
#ifdef _ADVANCED_GEOM_CONST
		if (!conc.empty())
			cout << "\nCOMMON NEIGHBOUR CONSTRAINTS:" << endl;
		num = 0;
		for (std::list<std::list<string> >::iterator it2 = conc.begin(); it2 != conc.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
		if (!snc.empty())
			cout << "\nSECOND NEIGHBOUR CONSTRAINTS:" << endl;
		num = 0;
		for (std::list<std::list<string> >::iterator it2 = snc.begin(); it2 != snc.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
		if (!bvs.empty())
			cout << "\nBOND VALENCE SUM CONSTRAINTS:" << endl;
		num = 0;
		for (std::list<std::list<string> >::iterator it2 = bvs.begin(); it2 != bvs.end(); it2++)
		{
			cout << "NO " << ++num << endl;
			for (std::list<string>::iterator it = it2->begin(); it != it2->end(); it++)
				cout << " \"" << *it << "\"" << endl;
		}
#endif
		if (!cus_entry.empty())
			cout << "\nCUSTOM MOVES:" << endl;
		for (std::list<string>::iterator it = cus_entry.begin(); it != cus_entry.end(); it++) cout << " \"" << *it << "\"" << endl;
		if (!swp.empty())
			cout << "\nSWAP MOVES:" << endl;
		for (std::list<string>::iterator it = swp.begin(); it != swp.end(); it++) cout << " \"" << *it << "\"" << endl;
#ifdef _LOCAL_INV
		if (!loc.empty())
			cout << "\nLOCAL INVARIANCE:" << endl;
		for (std::list<string>::iterator it = loc.begin(); it != loc.end(); it++) cout << " \"" << *it << "\"" << endl;
#endif
#ifdef _VIBR_AMP
		if (!vibramp.empty())
			cout << "\nVIBRATIONAL AMPLITUDE:" << endl;
		for (std::list<string>::iterator it = vibramp.begin(); it != vibramp.end(); it++) cout << " \"" << *it << "\"" << endl;
#endif
#ifdef _NO_PERIODIC
		if (!noper.empty())
			cout << "\nNO PERIODIC BOUNDARY CONDITIONS:" << endl;
		for (std::list<string>::iterator it = noper.begin(); it != noper.end(); it++) cout << " \"" << *it << "\"" << endl;
#endif
		//AENET and NBPOT - BPOT is mutually exclusive, but the unnecessary tags were already removed
		if (!aenet.empty() )
			cout << "\nAENET RELATED::" << endl;
		for (std::list<string>::iterator it = aenet.begin(); it != aenet.end(); it++) cout << " \"" << *it << "\"" << endl;
		if (!nbpot.empty())
			cout << "\nNON-BONDED POTENTIAL RELATED::" << endl;
		for (std::list<string>::iterator it = nbpot.begin(); it != nbpot.end(); it++) cout << " \"" << *it << "\"" << endl;

		if (!bpot.empty())
			cout << "\nBONDED POTENTIAL RELATED:" << endl;
		for (std::list<string>::iterator it = bpot.begin(); it != bpot.end(); it++) cout << " \"" << *it << "\"" << endl;
	}

	lead_series_ind = -1;//Initialization of leading series index
	lead_series_ind2 = -1;//Initialization of leading series index
	//set some defaults, in case there are reading errors, which does not prevent running the simulatio	swap_fraction = SWAP_FRACTION_DEF;
	bool result = ReadFreeGeneral(gener ) && ReadFreeExp(exper ) && ReadFreeCos(coscon ) && ReadFreeCoord(cocon ) && ReadFreeAcoord(avcocon ) && ReadFreeCustom(cus_entry ) && ReadFreeSwap(swp ) && ReadFreeNbpot(nbpot) && ReadFreeBpot(bpot);
#ifdef _AENET
	result=result && (ReadFreeAenet(aenet));
#endif
#ifdef _ADVANCED_GEOM_CONST
	result = result && (ReadFreeConc(conc ));
	result = result && (ReadFreeSnc(snc ));
	result = result && (ReadFreeBvs(bvs));
#endif
#ifdef _LOCAL_INV
	result = result && (ReadFreeLoc(loc ));
#endif
#ifdef _VIBR_AMP
	result= result && ReadFreeVibramp(vibramp );
#endif
#ifdef _NO_PERIODIC
	result = result && ReadFreeNoper(noper );
#endif
#ifdef _USE_LOCAL_INV
	if (potential > 0 || fnc == 4)
	{
		cout << "\n*****ERROR*****" << endl;
#ifdef _LOCAL_INV
		cout << "LOCAL INVARIANCE ";
#ifdef _NO_PERIODIC
		cout << "and ";
#endif
#endif

#ifdef _NO_PERIODIC
		cout << "PERIODIC BOUNDARY CONDITIONS ";
#endif
		cout << "cannot be used with non-bonded or bonded potential!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
#endif
	if (!result)
	{
		cout << "Cannot continue, exiting..." << endl;
		CleanExit();//Process the datasets and exit in the case of error
	}

	//setting the lead series indices
	
	if (potential == 0)//no non-bonded potential, all the sets, constraints and interaction contribute to the same chi2
	{
		if (fnc == 4)//there are flexible molecules, which contribute to normal chi2 and can be leading series index
		{
			if (lead_series_ind > 0 && lead_series_ind2 > 0)
				lead_series_ind = LEAD_SERIES_IND_DEF;//lead series was set multiple times, set default
			else
			{
				if (lead_series_ind2 > 0)
				{
					lead_series_ind = lead_series_ind2;//bonded interaction is the lead
					lead_series_ind2 = -1;
				}
			}
			// if only normal constraint and data set is the lead, that was set already
		}
		else//normal RMC without any potential
		{
			if (lead_series_ind == -2)//was set incorrectly, set default
				lead_series_ind = LEAD_SERIES_IND_DEF;
		}
		//cout << endl << "\nWARNING(" << ++warn << "):lead potential series index was set in BPOT section and no non-bonding potential is used.\nIt will be ignored, as the bonding potential chi2 contributes to the normal chi2!" << endl;
	}
	else//there is non-bonding potential, the potentila interaction contribute to chi_pot and have their own leading series
	{
		if (lead_series_ind2 <0 && potential > 0)//was not set correctly neither by NBPOT nor by BPOT
		{
			cout << "Default: " << LEAD_SERIES_IND2_DEF << " will be set" << endl;
			lead_series_ind2 = 1;
		}
	}

	cout << "\n******************************************" << endl;
	cout << "**** Number of warning messages: " << J4 << warn << " ****" << endl;
	cout << "**** Number of notes:            " << J4 << note << " ****" << endl;
	cout << "******************************************\n" << endl;

	CoordNumbConst::nthreads = nthreads;
	AvCoordConst::nthreads = nthreads;
#ifdef _ADVANCED_GEOM_CONST
	BondValenceSumConst::nthreads = nthreads;
#endif
	
	return;
}

bool RunParams::ReadFreeGeneral(std::list<string> &pool)
//Process free format strings related to general setup. Returns with false, if incoming datasets are insufficient.
{
	
	if (ntypes == 0)
	{
		cout << "\n" << "****ERROR!****" << endl;
		cout << "RunParams::GetParams " << endl;
		cout << "Number of types is 0 in RunParams object" << endl;
		cout << "Please, check the configuration file!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		return false;
	}
	int npartials = ntypes * (ntypes + 1) / 2;
	bool set_histstep_flag = false;
	bool set_cfg_coll_flag = false;
#ifndef _TEST_MODE
	int runtime = 0;
#endif
	int last=0;
	char *name;

	SetArraysize(&name, NAME_SIZE, "name", "RunParams::ReadFreeGeneral");
	SetArraysize(&pcutoff, npartials, "cutoff array", "RunParams::ReadFreeGeneral");
	SetArraysize(&pmaxmove, ntypes, "maximal moves array", "RunParams::ReadFreeGeneral");

	bool retval = true;
	//Initialize by default values
	//if(std::find_if(pool.begin(),pool.end(), [](string c) { return c.find("TITLE")!=string::npos; }) == pool.end()) title="";
	if (std::find_if(pool.begin(), pool.end(),
		[](string c) { return (c.find(CheckKey("NDENS")) != string::npos) || (c.find(CheckKey("HALFBOX")) != string::npos); }) == pool.end())
	{
#ifdef _NO_PERIODIC
		rho = natoms/(4.0 / 3.0 * pow(SimpleCfg::boxedge/2.0, 3) * PI);
#else
		rho = natoms / pow(2 * SimpleCfg::boxedge, 3);
#endif
		cout << "\nNOTE(" << ++note << "): The number density (" << rho;
		cout << " A^-3) has been set according to the information from the " << cfgfilename << " file" << endl;
#ifdef _NO_PERIODIC
		cout << "\tusing the total number of atoms ( " << SimpleCfg::ntotal << ") and the radius of the sample (" << SimpleCfg::boxedge/2.0 << "A)" << endl;
#else
		cout << "\tusing the total number of atoms ( " << SimpleCfg::ntotal << ") and the the halfbox of the cell (" << SimpleCfg::boxedge << "A)" << endl;
#endif
		
	}
	if (std::find_if(pool.begin(), pool.end(),
		[](string c) { return c.find(CheckKey("MAX-MOVES")) != string::npos; }) == pool.end())
	{
		for (int i = 0; i < ntypes; i++) pmaxmove[i] = 0.1;
	}
	
#ifdef _TEST_MODE
	pool.remove_if([](string c) { return c.find(CheckKey("RUN-TIME")) != string::npos; });
#else
	for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
	{
		if ((*it).find(CheckKey("RUN-TIME")) != string::npos)
			runtime++;//present
	}
#endif
	
	for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
	{
		if ((*it).find(CheckKey("LAST-ACC")) != string::npos)
			last++;//present
	}
	for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
	{
		if ((*it).find(CheckKey("LAST-GEN")) != string::npos)
			last++;//present
	}
		
	
#ifdef _TEST_MODE
	if (last > 1)
	{
		for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
		{
			if ((*it).find(CheckKey("LAST-GEN")) != string::npos || (*it).find(CheckKey("LAST-ACC")) != string::npos)
			{
				pool.erase(it);
				last--;
				if (last == 1)
					break;
			}
		}
	}
#else
	if (runtime > 1)
	{

		for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
		{
			if ((*it).find(CheckKey("RUN-TIME")) != string::npos)//the last is kept, if there is multiple entries, RUN-TIME should be kept
			{
				pool.erase(it);
				runtime--;
				if (runtime == 1)
					break;
			}
		}
	}
	if (runtime+last>1)
	{
		for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
		{
			if ((*it).find(CheckKey("LAST-GEN")) != string::npos || (*it).find(CheckKey("LAST-ACC")) != string::npos)
			{
				pool.erase(it);
				last--;
				if (runtime+last == 1)
					break;
			}
		}
	}
#endif

	//most of it already got the default value at the definition anyhow
	too_close_fraction = TOOCLOSE_FRACTION_DEF;
	rspacing_def = RSPACING_DEF;
	moveout = MOVEOUT_DEF;
	printstep = PRINTSTEP_DEF;
#ifdef _TEST_MODE
	runlimit = -LAST_MOVE_DEF;
	runmode = 1;
#else
	runlimit = TIMELIM_DEF;
	runmode = RUN_MODE;
#endif
	
	timesave = TIMESAVE_DEF;//SAVE-TIME
	histbuffsize = HIST_BUFF_DEF;
	histstepratio = HIST_STEP_DEF;
	cfgnumb = CFG_COLL_DEF; 
	coll_frequency = CFG_COLL_FREQ_DEF; 
	binshift_def = BINSHIFT_DEF;
	xmax[0] = XMAX_DEF;
	reload = RELOAD_DEF;
	max_gridatom = MAX_GRIDATOM_DEF;
	exafs_shiftstep = E0_SHIFTSTEP_DEF;
	nthreads = NTHREADS_DEF;
	old_out = OLD_OUT_DEF;
	sum_ppcf = SUM_PPCF_DEF;
	niter_nonlin = NONLIN_MAX_NIT_DEF;
	epsilon_nonlin = NONLIN_EPSILON_DEF;
	lambda_nonlin = NONLIN_LAMBDA_DEF;
	factor_nonlin = NONLIN_FACTOR_DEF;
	terminate_nonlin = NONLIN_TERMINATE_DEF;
	//some parameters are initialized at definition
	if (AUTO_CUTOFF_DEF == 1)
	{
		if (abs(continuation) == 1)//set the default, if the CUT-OFF key world is not given
			auto_cutoff = -1;
		else
			auto_cutoff = 1;
	}
	else
		auto_cutoff = 0;
	for (int i = 0; i < npartials; i++)
		pcutoff[i] = 0.0;//set the cutoff to zero for the time being to avoid FNC-conflict
	
	for (std::list<string>::iterator it = pool.begin(); (it != pool.end()) && retval; it++)
	{
		string key, params;
		WrapFree(*it, key, params);

		stringstream ss(params);

		if (key == CheckKey("TITLE"))
		{
			title=params;
		//	std::copy(params.begin(), params.end(), title);
		}
		else if ((key == CheckKey("NDENS")) || (key == CheckKey("HALFBOX")))
		{
			if (!(ss >> rho)) return ErrLackValue(*it);
			if (key == CheckKey("HALFBOX")) rho = natoms / pow(2 * rho, 3);
			AdjustBoxedge();
		}
		else if (key == CheckKey("CUT-OFF"))
		{
			bool ret = true;
			string first;
			if (!(ss >> first)) return ErrLackValue(*it);
			if (first.compare("MIN-ALL") == 0)
			{
				if (abs(continuation) == 1)
					auto_cutoff = -1;

				else
					auto_cutoff = 1;
			}
			else if (first.compare("MIN") == 0)
			{

				if (abs(continuation) == 1)
					auto_cutoff = -3;
				else
					auto_cutoff = 3;
			}
			else
			{
				if ((std::find_if(first.begin(), first.end(), [](char c) { return std::isalpha(c); }) != first.end()))
				{
					cout << "\nERROR: You should give MIN or " << npartials << " values in the following line: " << endl << "\t" << *it << endl;
					return false;
				}
				auto_cutoff = 0;
				pcutoff[0] = std::stod(first);

				for (int i = 1; (i < npartials) && ret; i++) ret = static_cast<bool>(ss >> pcutoff[i]);
				if (!ret)
				{
					cout << "\nERROR: too few parameters. You should provide " << npartials << " values in the following line: " << endl << "\t" << *it << endl;
					return false;
				}
			}
		}
		else if (key == CheckKey("MAX-MOVES"))
		{
			int last = -1;
			bool rv = true;
			for (int i = 0; i < ntypes; i++)
			{
				if (rv) rv = static_cast<bool>(ss >> pmaxmove[i]);
				if (!rv)
				{
					if (last == -1)
					{
						cout << "\nERROR: too few parameters. You should provide at least one value in the following line: " << endl<<"\t"<<*it << endl;
						return false;
					}
					else pmaxmove[i] = pmaxmove[last];
				}
				else last = i;
			}
			if ((last != ntypes - 1) && ::debug)
				cout << "\nNOTE(" << ++note << "): maximal atomic moves from type " << last + 2 << " to type" << ntypes << " are/is adjusted to " << pmaxmove[last] << endl;
		}
		else if (key == CheckKey("R-SPACING"))
		{
			if (!(ss >> rspacing_def)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("MOVEOUT"))
		{
			if (!(ss >> moveout)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "MOVEOUT");
			Check0_1(moveout, name, "RunParams::ReadFreeGeneral");
		}
		else if (key == CheckKey("TOOCLOSE-FRACTION"))
		{
			if (!(ss >> too_close_fraction)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("PRINT-STEP"))
		{
			if (!(ss >> printstep)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("RUN-TIME"))
		{
			if (!(ss >> runlimit)) return ErrLackValue(*it);
			if (runlimit < 0)
			{
				cout << "\nWARNING(" << ++warn << "): Negative value was given for " << CheckKey("RUN-TIME") << ", it will be considered as number of moves to generate" << endl;
				cout << "\tas in case of the fixed format input! Next time use " << CheckKey("LAST-GEN") << " with positive value instead!" << endl;
				runmode = 1;
				runlimit = round(runlimit);

			}
			
		}
		else if (key == CheckKey("LAST-GEN"))//now it can be used by normal mode as well as an alternative to RUN-TIME and LAST-ACC
		{
			if (!(ss >> runlimit)) return ErrLackValue(*it);
			if (runlimit > 0)
				runlimit *= -1.0;//to be compatible with fixed format
			runmode = 1;
			runlimit = round(runlimit);
		}
		else if (key == CheckKey("LAST-ACC"))//it is an alternative to RUN-TIME (and LAST-GEN) both in normal mode as well as in _TEST_MODE
		{
			if (!(ss >> runlimit)) return ErrLackValue(*it);
			runmode = 2;
			if (runlimit < 0)
			{
				cout << "\nWARNING(" << ++warn << "): Negative value was given for " << CheckKey("LAST-GEN") << ", absolute value will be used!" << endl;
				runlimit *= -1.0;//to be positive
			}
			runlimit = round(runlimit);
		}
		else if (key == CheckKey("SAVE-TIME"))
		{
			if (!(ss >> timesave)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("HST_BUFFSIZE") || key == CheckKey("HST-BUFFSIZE"))
		{
			if (!(ss >> histbuffsize)) return ErrLackValue(*it);
		}
		else if ((key == CheckKey("HST-STEP-FACTOR")) || (key == CheckKey("HST-SAVE-TIME")))
		{
			if (!(ss >> histstepratio)) return ErrLackValue(*it);
			if (key == CheckKey("HST-SAVE-TIME"))
				set_histstep_flag = true;//not sure, that timesave was already given
		}
		else if ((key == CheckKey("COLL-STEP-FACTOR")) || (key == CheckKey("COLL-SAVE-TIME")))
		{
			if (!(ss >> coll_frequency)) return ErrLackValue(*it);
			if (key == CheckKey("COLL-SAVE-TIME"))
				set_cfg_coll_flag = true;//not sure, that timesave was already given
		}
		else if (key == CheckKey("COLL-NUMBER"))
		{
			if (!(ss >> cfgnumb)) return ErrLackValue(*it);
			if (cfgnumb <= 0)
			{
				coll_frequency = 0;
				cfgnumb = CFG_COLL_DEF;
			}
		}
		else if (key == CheckKey("FNC-TYPE"))
		{
			if (fnc != -4)//this was set as BPOT tag was given in the dat file, indicating, that flexible molecules will be used
			{
				params.erase(std::remove_if(params.begin(), params.end(), [](char c) { return c == ' '; }), params.end());
				if (params.length() == 0) return  ErrLackValue(*it);
				if (params == "NONE") fnc = 0;
				else if (params == "NORMAL") fnc = 1;
				else if (params == "ADJUST") fnc = 2;
				else if (params == "MOVE_IN") fnc = 3;
				else
				{
					cout << "\nERROR: Wrong entry in the following line: " << endl<<"\t"<<*it << endl;
					return false;
				}
			}
		}
		else if (key == CheckKey("BIN-SHIFT"))
		{
			if (!(ss >> binshift_def)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("XMAX-FACTOR"))
		{
			if (!(ss >> xmax[0])) return ErrLackValue(*it);
		}
		else if (key == CheckKey("RELOAD"))
		{
			if (!(ss >> reload)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "RELOAD");
			Check0_1(reload, name, "RunParams::ReadFreeGeneral");
		}
		else if (key == CheckKey("ATOMS-IN-GRIDCELL"))
		{
			if (!(ss >> max_gridatom)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("CUSTOM-SFACTORTABLE"))
		{
			if (!(ss >> sffilename))	return ErrLackValue(*it);
			else use_custom_sftable = true;
		}
		else if (key == CheckKey("EPSILON-NONLIN"))
		{
			if (!(ss >> epsilon_nonlin)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("MAX-NITER-NONLIN"))
		{
			if (!(ss >> niter_nonlin)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("LAMBDA-NONLIN"))
		{
			if (!(ss >> lambda_nonlin)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("FACTOR-NONLIN"))
		{
			if (!(ss >> factor_nonlin)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("TERMINATE-NONLIN"))
		{
			if (!(ss >> terminate_nonlin)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("AXS-SHIFTSTEP"))
		{
			if (!(ss >> AXS_shiftstep)) return ErrLackValue(*it);
			else
			{
				if (AXS_shiftstep <= 0)
				{
					cout << "\nWARNING(" << ++warn << "): For " << CheckKey("AXS-SHIFTSTEP") << " only positive values can be given, " << AXS_shiftstep << endl;
					cout << "\twas read from file " << datfilename << "! It is set to the default value " << AXS_SHIFTSTEP_DEF << endl;
					AXS_shiftstep = E0_SHIFTSTEP_DEF;
				}
			}
		}
		else if (key == CheckKey("I(Q)BACKG-STEP"))
		{
			if (!(ss >> IQ_backg_corr_step)) return ErrLackValue(*it);
			else
			{
				if (IQ_backg_corr_step < 0)
				{
					cout << "\nWARNING(" << ++warn << "): For " << CheckKey("I(Q)BACKG-STEP") << " only positive values or 0 can be given, (in case of 0 fixed background correction is used), " <<endl;
					cout << "\t"<< IQ_backg_corr_step << " was read from file " << datfilename << "! It is set to the default value "<<IQ_BACKG_CORR_STEP_DEF<< endl;
					IQ_backg_corr_step = IQ_BACKG_CORR_STEP_DEF;
				}
			}
		}
		else if (key == CheckKey("EXAFS-SHIFTSTEP"))
		{
			if (!(ss >> exafs_shiftstep)) return ErrLackValue(*it);
			else
			{
				if (exafs_shiftstep <= 0)
				{
					cout << "\nWARNING(" << ++warn << "): For " << CheckKey("EXAFS-SHIFTSTEP") << " only positive values can be given, " << exafs_shiftstep << endl;
					cout << "\twas read from file " << datfilename << "! It is set to the default value " << E0_SHIFTSTEP_DEF << endl;
					exafs_shiftstep = E0_SHIFTSTEP_DEF;
				}
			}
		}
		else if (key == CheckKey("WRITE-LOG"))
		{
			if (!(ss >> write_log)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "WRITE-LOG");
			Check0_1(write_log, name, "RunParams::ReadFreeGeneral");
		}
		else if (key == CheckKey("DEBUG"))
		{
			if (!(ss >> debug)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "DEBUG");
			Check0_1(debug, name, "RunParams::ReadFreeGeneral");
		}
		else if (key == CheckKey("LOG-NONLIN-STEPS"))
		{
			if (!(ss >> log_nlr_steps)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "LOG-NONLIN-STEPS");
			Check0_1(log_nlr_steps, name, "RunParams::ReadFreeGeneral");
		}
		
		else if (key == CheckKey("WRITE-EXAFS-COEFFS"))
		{
			if (!(ss >> write_exafs_coeffs)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "WRITE-EXAFS-COEFFS");
			Check0_1(write_exafs_coeffs, name, "RunParams::ReadFreeGeneral");
		}
		else if (key == CheckKey("R-SWITCH-POWER"))
		{
			if (!(ss >> r_switch_power)) return ErrLackValue(*it);
		}
		else if (key == CheckKey("CHEMICAL-SYMBOLS"))
		{
			int i;
			if (chem_symbols == nullptr)
				SetArraysize(&chem_symbols, ntypes, "chem_symbols", "RunParams::ReadFreeGeneral");
			if (chem_symbols_standard == nullptr)
				SetArraysize(&chem_symbols_standard, ntypes, "chem_symbols_standard", "RunParams::ReadFreeGeneral");
			
			for (i = 0; i < ntypes; i++)
				if (!(ss >> chem_symbols[i])) return ErrLackValue(*it);
			for (i = 0; i < ntypes; i++)
				chem_symbols_standard[i] = ExtractChemSymbol(chem_symbols[i]);
			
		}
		else if (key == CheckKey("THREADS"))
		{
			if (!(ss >> nthreads)) return ErrLackValue(*it);
			if (nthreads < 1)
			{
				cout << "\nWARNING(" << ++warn << "): The number of threads has to be positive! Resetting the number of threads to 1" << endl;
				nthreads = 1;
			}
		}
		else if (key == CheckKey("CREATE-OUT"))
		{
			if (!(ss >> old_out)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "CTREATE-OUT");
			Check0_1(old_out, name, "RunParams::ReadFreeGeneral");
		}
		else if (key == CheckKey("PPCF-AVERAGE"))
		{
			if (!(ss >> sum_ppcf)) return ErrLackValue(*it);
			mystrcpy(name, NAME_SIZE, "PPCF-AVERAGE");
			Check0_1(sum_ppcf, name, "RunParams::ReadFreeGeneral");
			
		}
		else
		{
			CheckKeyTag(key, CheckTag("GENERAL"));//write the appropriate error mesage
			return false;
		}
	}
	if (set_histstep_flag)
	{
		if (histbuffsize == 0)//buffsize was not given, but step or save time was, record history
			histbuffsize = 1;
		histstepratio =(int) round(histstepratio / timesave) - 1;
		if (histstepratio < 0)
			histstepratio = 0;
	}
	if (set_cfg_coll_flag)
	{
		coll_frequency = (int)round(coll_frequency / timesave) - 1;
		if (coll_frequency < 0)
			coll_frequency = 0;
	}
	if (reload && abs(auto_cutoff)&1)
	{
		cout << "\nWARNING(" << ++warn << "): Cannot use RELOAD with automatic cutoff determination, the histogram will be calculated!" << endl;
		reload = 0;
	}
	if (log_nlr_steps)
		write_log = 1;//make sure that the logfile is open
	if (name != NULL)
		delete[] name;
	return retval;
}



bool RunParams::ReadFreeExp(std::list<std::list<string> > &pool)
//Process free format strings related to setup of experimental constraint. Returns with false, if incoming datasets are insufficient.
{
	ngr = 0; nsq = 0; nfq = 0; nfg = 0; nek = 0;
	int j = 0;
	int size=1;
	int i = 1;//EXP entry number
	using std::list;
	for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
	{
		if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("TYPE")) != string::npos; }) == it->end()) ||
			(std::find_if(it->begin(), it->end(), [](string c) { return c.find("DATAFILE") != string::npos; }) == it->end()) ||
			(std::find_if(it->begin(), it->end(), [](string c) { return (c.find("SIGMA") != string::npos); }) == it->end()))
		{
			cout << "\nERROR!!: **TYPE**,**DATAFILE** and **SIGMA(-MASTER)**/**SIGMA-SCALABLE** are mandatory keywords in EXP section." << endl;
			
			cout << " At least, one of them is missing in the " << i << ". "<<CheckTag("EXP")<<" section!" << endl;
			return false;
		}
		if (std::count_if(it->begin(), it->end(), [](string c) { return (c.find(CheckKey("SIGMA")) != string::npos) || (c.find(CheckKey("SIGMA-MASTER")) != string::npos) || (c.find(CheckKey("SIGMA-SCALABLE")) != string::npos); }) != 1)
		{
			cout << "\nWARNING(" << ++warn << "): Multiple declaration of SIGMA/SIGMA-MASTER/SIGMA-SCALABLE has been found in one of the " <<i<<". "<< CheckTag("EXP") << " section." << endl;
			cout<<"\tThe last one will be used!" << endl;
		}
		bool coeff_flag = false;//This is GR 
		bool coeff_flag2 = false;//This is ND
		bool exafsflag = false;//This is EXAFS
		bool ppcf = true;
		bool abst = true;
		bool bsf = true;
		bool isotopeflag = true;
		bool grid = true;
		bool agrid = true;

		for (std::list<string>::iterator it2 = it->begin(); it2 != it->end(); it2++)
		{
			if ((it2->find(CheckKey("TYPE")) != string::npos) && (it2->find(CheckKeyValue("TYPE", "GR")) != string::npos))
			{ 
				coeff_flag = true;
				ngr++; 
			}
			if ((it2->find(CheckKey("TYPE")) != string::npos) && (it2->find(CheckKeyValue("TYPE", "ND")) != string::npos))
			{
				coeff_flag2 = true;
				nsq++;
			}
			if ((it2->find(CheckKey("TYPE")) != string::npos) && (it2->find(CheckKeyValue("TYPE", "XRD")) != string::npos))
				nfq++;
			if ((it2->find(CheckKey("TYPE")) != string::npos) && (it2->find(CheckKeyValue("TYPE", "EDIFF")) != string::npos))
				nfg++;
			if ((it2->find(CheckKey("TYPE")) != string::npos) && (it2->find(CheckKeyValue("TYPE", "EXAFS")) != string::npos))
			{
				exafsflag = true;
				nek++; 
			}
			if (it2->find(CheckKey("PARTIAL-COEFFS")) != string::npos)
				ppcf = false;//present
			if (it2->find(CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS")) != string::npos)
				isotopeflag = false;
			if (it2->find(CheckKey("ABSORBER-TYPE")) != string::npos)
				abst = false;
			if (it2->find(CheckKey("BACKSCATT-FILE")) != string::npos)
				bsf = false;
			if (it2->find(CheckKey("DELTAE0_NGRID")) != string::npos)
				grid = false;
			if (it2->find(CheckKey("E0GRID-FROM_TO_START")) != string::npos)
				agrid = false;
		}
		if (coeff_flag&&ppcf)//GR
		{
			cout << "\nERROR!!: The **"<<CheckKey("PARTIAL-COEFFS")<<"** is mandatory keyword in  " << CheckTag("EXP")<<" section if " << CheckKey("TYPE") << " is GR " << endl;
			cout << "It is missing in the " << i << ". "<<CheckTag("EXP")<<" section!" << endl;
			return false;
		}
		if (coeff_flag2 && (!ppcf && !isotopeflag))//ND, if PARTIAL-COEFFS and ISOTOPE-COUNT_SYMBOLS_RATIOS both given, PARTIAL-COEFFS will be used, remove the others
		{
			cout << "\nWARNING(" << ++warn << "): Both " << CheckKey("PARTIAL-COEFFS") << " and " << CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS") << " was given in " << endl;
			cout << "\t" <<i << ". "<< CheckTag("EXP") << " section, describing a neutron scattering data set " <<"! " << CheckKey("PARTIAL-COEFFS") << " will be used, the line()s) containing " << endl;
			cout << "\t" << CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS") << " will be removed!" << endl;
			it->remove_if([](string c) { return c.find(CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS")) != string::npos; });
		}
		
		if (coeff_flag2 && ppcf && chem_symbols==nullptr)//ND
		{
			cout << "\nERROR!!: The **" << CheckKey("PARTIAL-COEFFS") << "** is missing in the " << i << ". " << CheckTag("EXP") << " section describing an ND set," << endl;
			cout << "and the **" << CheckKey("CHEMICAL-SYMBOLS") << "** key word is not present in the " << CheckTag("GENERAL") << " section!" << endl;
			cout << "One of them is needed for the neutron scattering coefficients!" << endl;
			return false;
		}
		int icount = std::count_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS")) != string::npos; });
		if ( icount> ntypes)
		{

			cout << "\nWARNING(" << ++warn << "): "<< "More than "<<ntypes<<" "<<CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS")<<" can be found for the "<<i <<". "<<CheckTag("EXP") << " section." << endl;
			cout << "\tThe last "<<ntypes<<" will be used!" << endl;
			for (j = 0; j < icount-ntypes; j++)
			{
				std::list<string> ::iterator fit = std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS")) != string::npos; });
				cout << "Removing lines: "<<*fit << endl;
				it->erase(fit);
			}
		}
		if (coeff_flag2 && ppcf)//neutron coeff should be calculated
			ExptsData::is_Ncoeff_calc = 1;
		if (exafsflag && (abst || bsf))
		{
			cout << "\nERROR!!: The **" << CheckKey("ABSORBER-TYPE") << "** and **" << CheckKey("BACKSCATT-FILE") << "** are mandatory keywords in " << CheckTag("EXP")<<" section if " << CheckKey("TYPE")<<" is EXAFS." << endl;
			cout << " At least, one of them is missing in the " << i << ". "<<CheckTag("EXP")<<" section!" << endl;
			return false;
		}
		if (exafsflag && !agrid && grid)
		{
			cout << "\nERROR!!: The **" << CheckKey("E0GRID-FROM_TO_START") << "** is present, but **" << CheckKey("DELTAE0_NGRID") << "** is not in the " <<i << ". " << CheckTag("EXP") << " section" << endl;
			cout << "describing an EXAFS set. "<<CheckKey("DELTAE0_NGRID")<<" has to be present, if " << CheckKey("E0GRID-FROM_TO_START") << " is used!"<<endl;
			return false;
		}
	}
	ntot_datasets = ngr + nsq + nfq + nfg + nek;
	if (ntot_datasets != (int)pool.size())
	{
		cout << "\nERROR!!: Some of TYPE dataset is out of the valid range. It/they should be "<<GiveKeyValues("TYPE")<<". Or maybe you defined TYPE labels more than once for a data set." << endl;
		return false;
	}

	pool.sort(CompareExp);
	if (ntot_datasets > size)
		size = ntot_datasets;
	SetArraysize(&rspacing, size, "rspacing", "RunParams::ReadFreeExp");
	SetArraysize(&binshift, size,"binshift", "RunParams::ReadFreeExp");
	SetArraysize(&firstbin, size, "firstbin", "RunParams::ReadFreeExp");
	if (fabs(binshift_def) < LOAD_TOL)
		binshift_def = 0;
	else
	{
		//as due to rounding error in case of binshift 1, 2... the (int)binshift is not necessarily correct, an integer firstbin will be calculated here
		if (fabs(binshift_def - round(binshift_def)) < LOAD_TOL)//this supposed to be integer
			j = (int)round(binshift_def);
		else
			j = (int)binshift_def;
	}
	for (i = 0; i < size; i++)
	{
		rspacing[i] = rspacing_def;
		binshift[i] = binshift_def;
		firstbin[i] = j;
	}

	SetArraysize(&assign_hist, size, "assign_hist", "RunParams::ReadFreeExp");

	//read the parameters for the data constraints
	return ExptsData::GetExptParamFree(pool);
}

bool RunParams::CompareExp(std::list<string> &fst, std::list<string> &scnd)
//Compare two objects according to gr,Sq,Fq,Fg,Ek order
{
	bool nfound = true;
	int first_index = 10;
	for (std::list<string>::iterator it2 = fst.begin(); (it2 != fst.end()) && nfound; it2++)
	{
		if (it2->find("TYPE") != string::npos)
		{
			nfound = false;
			if (it2->find("GR") != string::npos) first_index = 1;
			if (it2->find("ND") != string::npos) first_index = 2;
			if (it2->find("XRD") != string::npos) first_index = 3;
			if (it2->find("EDIFF") != string::npos) first_index = 4;
			if (it2->find("EXAFS") != string::npos) first_index = 5;
		}
	}
	for (std::list<string>::iterator it2 = scnd.begin(); it2 != scnd.end(); it2++)
	{
		if (it2->find("TYPE") != string::npos)
		{
			if (it2->find("GR") != string::npos) 
				return false;
			if (it2->find("ND") != string::npos) 
				return (first_index < 2);
			if (it2->find("XRD") != string::npos) 
				return (first_index < 3);
			if (it2->find("EDIFF") != string::npos)
				return (first_index < 4);
			if (it2->find("EXAFS") != string::npos) 
				return (first_index < 5);
		}
	}
	return true;
}
bool RunParams::ReadFreeCos(std::list<std::list<string> > &pool )
//Process free format strings related to setup of cosine distribution constraint. Returns with false, if incoming datasets are insufficient.
{
	ncosdistr = (int)pool.size();//number of coordination constrains
	if (ncosdistr > 0)
	{
		//Check the required keywords
		int i = 1;
		for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
		{
			if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("CENT-TYPE")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return (c.find(CheckKey("SIGMA")) != string::npos); }) == it->end()))
			{
				cout << "\nERROR!!: **" << CheckKey("CENT-TYPE") << "**, **" << CheckKey("DISTRIB-TYPE") << "** and one of **" << CheckKey("SIGMA") << "(-MASTER)**/**" << CheckKey("SIGMA-SCALABLE") << "** are mandatorykeywords in " << CheckTag("COS")<<" section." << endl;
				cout << " At least, one of them is missing in the " << i << ". " << CheckTag("COS")<<" section!" << endl;
				return false;
			}
			int neigh_counter = 0;
			bool dd = true;//DISTRIB-DEGREES
			bool dw = true;//DISTRIB-WIDTH
			bool findicator = true;//DISTRIB-TYPE!=FILE indicator
			for (std::list<string>::iterator it2 = it->begin(); it2 != it->end(); it2++)
			{
				if (it2->find(CheckKey("NEIGH-TYPE_FROM_TO")) != string::npos) neigh_counter++;
				if (it2->find(CheckKey("DISTRIB-DEGREES")) != string::npos) dd = false;
				if (it2->find(CheckKey("DISTRIB-WIDTH")) != string::npos) dw = false;
				if ((it2->find(CheckKey("DISTRIB-TYPE")) != string::npos) && (it2->find("UNIFORM") == string::npos) && (it2->find("GAUSSIAN") == string::npos) && (it2->find("ABSENT") == string::npos)) findicator = false;
			}
			if (neigh_counter != 2)
			{
				cout << "\nERROR!!: Two **" << CheckKey("NEIGH-TYPE_FROM_TO") << "** are mandatory keywords in " << CheckTag("COS")<<" section." << endl;
				cout << " You provided less or more in the " << i << ". " << CheckTag("COS")<<" section!" << endl;
				return false;
			}
			if (findicator && (dd || dw))
			{
				cout << "\nERROR!!: The **" << CheckKey("DISTRIB-DEGREES") << "** and **" << CheckKey("DISTRIB-WIDTH") << "** are mandatory keywords in " << CheckTag("COS")<<" section if "<<CheckKey("DISTRIB-TYPE")<<" is other than filename." << endl;
				cout << " At least, one of them is missing in the " << i << ". " << CheckTag("COS")<<" section!" << endl;
				return false;
			}
		}
		return CosDistrConst::GetCosDistrConstFree(pool );
	}
	else CosDistrConst::nconstraints = 0;
	return true;
}

bool RunParams::ReadFreeCoord(std::list<std::list<string> > &pool )
//Process free format strings related to setup of coordination constraint. Returns with false, if incoming datasets are insufficient.
{
	nicoord = (int)pool.size();
	if (nicoord > 0)
	{
		//Check the required keywords
		int i = 1;
		for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
		{
			if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("CENT-TYPE")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("NEIGH-TYPE_FROM_TO")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return (c.find(CheckKey("COORDNUM_FRACT_SIGMA")) != string::npos); }) == it->end()))
			{
				cout << "\nERROR!!: **" << CheckKey("CENT-TYPE") << "**,**" << CheckKey("NEIGH-TYPE_FROM_TO") << "** and one of **" << CheckKey("COORDNUM_FRACT_SIGMA") << "(_MASTER)**/**" << CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE") << "**) are mandatory keywords in " << CheckTag("COORD")<<" section." << endl;
				
				cout << " At least one of them is missing in the " << i << ". " << CheckTag("COORD")<<" section!" << endl;
				return false;
			}
		}
		return CoordNumbConst::GetCoordConstFree(pool );
	}
	else
	{
		//Initialize some important variables asked by other functions
		CoordNumbConst::tot_neightype = 0;
		CoordNumbConst::tot_subconst = 0;
		CoordNumbConst::ncentral_tot = 0;
		CoordNumbConst::nconstraints = 0;
	}
	return true;
}


bool RunParams::ReadFreeAcoord(std::list<std::list<string> > &pool )
//Process free format strings related to setup of average coordination constraint. Returns with false, if incoming datasets are insufficient.
{
	navcoord = (int)pool.size();
	if (navcoord > 0)
	{
		//Check the required keywords
		int i = 1;
		for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
		{
			if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("CENT-TYPE")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("NEIGH-TYPE_FROM_TO")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return (c.find(CheckKey("COORDNUM_SIGMA")) != string::npos); }) == it->end()))
			{
				cout << "\nERROR!!: **" << CheckKey("CENT-TYPE") << "**,**" << CheckKey("NEIGH-TYPE_FROM_TO") << "** and **" << CheckKey("COORDNUM_SIGMA") << "(-MASTER)**/**" << CheckKey("COORDNUM_SIGMA-SCALABLE") << "** are mandatory keywords in " << CheckTag("AVCOORD")<<" section." << endl;
				cout << " At least one of them is missing in the " << i << ". "<<CheckTag("AVCOORD")<<" section!" << endl;

				return false;
			}
		}
		return AvCoordConst::GetAvCoordConstFree(pool );
	}
	else AvCoordConst::nconstraints = 0;
	return true;
}
#ifdef _ADVANCED_GEOM_CONST
bool RunParams::ReadFreeConc(std::list<std::list<string> > &pool )
//Process free format strings related to setup of common neighbour constraint. Returns with false, if incoming datasets are insufficient.
{
	ncommonneigh = pool.size();
	if (ncommonneigh > 0)
	{
		//Check the required keywords
		int i = 1;
		for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
		{
			if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("PTYPE1_PTYPE2_FROM_TO")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("STYPE_FROM1_TO1_FROM2_TO2")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return (c.find(CheckKey("COORDNUM_FRACT_SIGMA")) != string::npos); }) == it->end()))
			{
				cout << "\nERROR!!: **"<<CheckKey("PTYPE1_PTYPE2_FROM_TO")<<"**,**" << CheckKey("STYPE_FROM1_TO1_FROM2_TO2")<<"** and one of **"<<CheckKey("COORDNUM_FRACT_SIGMA")<<"(_MASTER)**/**"<<CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE")<<"**) are mandatory keywords in "<<CheckTag("CONC")<<" section." << endl;
				cout << " At least one of them is missing in the " << i << ". "<<CheckTag("CONC")<<" section!" << endl;
				return false;
			}
		}
		return CommonNeighConst::GetCommonNeighConstFree(pool );
	}
	else
	{
		//Initialize some important variables asked by other functions
		CommonNeighConst::nconstraints = 0;
		CommonNeighConst::tot_sec = 0;
	}
	return true;
}

bool RunParams::ReadFreeSnc(std::list<std::list<string> > &pool )
//Process free format strings related to setup of second neighbour constraint. Returns with false, if incoming datasets are insufficient.
{
	nsecondneigh = pool.size();
	if (nsecondneigh > 0)
	{
		//Check the required keywords
		int i = 1;
		for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
		{
			if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("CTYPE_FTYPE_FROM_TO")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("STYPE_FROM_TO")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return (c.find(CheckKey("COORDNUM_FRACT_SIGMA")) != string::npos); }) == it->end()))
			{
				cout << "\nERROR!!: **" << CheckKey("CTYPE_FTYPE_FROM_TO") << "**,**" << CheckKey("STYPE_FROM_TO") << "** and one of **" << CheckKey("COORDNUM_FRACT_SIGMA") << "(_MASTER)**/**" << CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE") << "**) are mandatory keywords in " << CheckTag("SNC") << " section." << endl;
				cout << " At least one of them is missing in the " << i << ". " << CheckTag("SNC") << " section!" << endl;
				return false;
			}
		}
		return SecondNeighConst::GetSecondNeighConstFree(pool );
	}
	else
	{
		//Initialize some important variables asked by other functions
		SecondNeighConst::nconstraints = 0;
		SecondNeighConst::tot_sec = 0;
	}
	return true;
}
bool RunParams::ReadFreeBvs(std::list<std::list<string> > &pool)
//Process free format strings related to setup of bond valence sum constraint. Returns with false, if incoming datasets are insufficient.
{
	nbvs = (int)pool.size();
	if (nbvs > 0)
	{
		//Check the required keywords
		int i = 1;
		for (std::list<std::list<string> >::iterator it = pool.begin(); it != pool.end(); it++, i++)
		{
			if ((std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("CENT-TYPE_CHARGE")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("NEIGH-TYPE_RMAX_R0_B_CHARGE")) != string::npos; }) == it->end()) ||
				(std::find_if(it->begin(), it->end(), [](string c) { return c.find(CheckKey("VALENCE_SIGMA")) != string::npos; }) == it->end()))
			{
				cout << "\nERROR!!: **" << CheckKey("CENT-TYPE_CHARGE") << "**,**" << CheckKey("NEIGH-TYPE_RMAX_R0_B_CHARGE") << "** and one of **" << CheckKey("VALENCE_SIGMA") << ", " << CheckKey("VALENCE_SIGMA-MASTER")<<", " << endl;
				cout<< CheckKey("VALENCE_SIGMA-SCALABLE") << "** are mandatory keywords in " << CheckTag("BVS") << " section." << endl;

				cout << " At least one of them is missing in the " << i << ". " << CheckTag("BVS") << " section!" << endl;
				return false;
			}
		}
		return BondValenceSumConst::GetBondValenceSumConstFree(pool);
	}
	else
	{
		//Initialize some important variables asked by other functions
		BondValenceSumConst::tot_neightype = 0;
		BondValenceSumConst::ncentral_tot = 0;
		BondValenceSumConst::nconstraints = 0;
	}
	return true;
}

#endif

bool RunParams::ReadFreeCustom(std::list<string> &pool )
//Process free format strings related to custom moves setup. Returns with false, if incoming datasets are insufficient.
{
	custmove = (pool.size() > 0 ? 1 : 0);
	if (custmove)
	{
		//Check the required keywords
		if (!CheckKeyInPool(pool, CheckKey("NMOVED-ATOMS")))
		{
			cout << "\nERROR!!: ****" << CheckKey("NMOVED-ATOMS") << "**** is mandatory item in CUSTMOVE section. Please, provide it to allow the program to run!" << endl;
			return false;
		}
		int ma = 0;

		for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);
			if (key == CheckKey("NMOVED-ATOMS"))
			{
				if (ma++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("NMOVED-ATOMS") << " definition in the " << CheckTag("CUSTMOVE") << " section!" << endl;
					cout << "\tThe last value will be used." << endl;
				}
				if (!GetInt(ss, *it, &nmoved,"RunParams::ReadFreeCustom")) return ErrLackValue(*it);
				if (nmoved < 1)
				{
					cout << "\nERROR: The actual value (" << nmoved << ") of number of moved atoms in a single move in " << CheckTag("CUSTMOVE")<<" section is out of valid range: >0." << endl;
					return false;
				}
			}
		}
	}
	else nmoved = NMOVED_DEF;
	if (ncommonneigh > 0 && custmove > 0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Custom move cannot be used with common neighbour constraint!" << endl;
		cout << "Exiting..." << endl;
		CleanExit();
	}
	if (ncommonneigh > 0 && nmoved > 1)
	{
		cout << "\nWARNING(" << ++warn << "): Common neighbour constraint can only be used with 1 moved atom! The number of moved atoms will be set to 1!" << endl;
		nmoved = 1;
	}
	return true;
}

bool RunParams::ReadFreeSwap(std::list<string> &pool )
//Process free format strings related to swap moves setup. Returns with false, if incoming datasets are insufficient.
{
	if ((pool.size() > 0) && (ntypes > 1))
	{
		//Check the required keywords
		if (!(CheckKeyInPool(pool, CheckKey("FRACT")) && CheckKeyInPool(pool, CheckKey("PAIRS"))))
		{
			cout << "\nWARNING(" << ++warn << ") ****" << CheckKey("FRACT")<<"**** and ****" << CheckKey("PAIRS") << "**** are mandatory items in " << CheckTag("SWAP") << " section, and they are missing. No swaps will be performed!" << endl;
			return true;
		}
		else std::iter_swap(pool.begin(), std::find_if(pool.begin(), pool.end(), [](string c) { return c.find(CheckKey("FRACT")) != string::npos; }));//Insert FRACT to head
		int fc = 0, pc = 0;
		for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("FRACT"))
			{
				if (fc++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple "<<CheckKey("FRACT")<<" definition in the "<<CheckTag("SWAP")<<" section!" << endl << "\tThe last value will be set used." << endl;
				if (!(ss >> swap_fraction)) return ErrLackValue(*it);
				if (swap_fraction > 1.0 || swap_fraction <= 0)
				{
					cout << "\nERROR: The actual value (" << swap_fraction << ") of fraction of swaps is out from valid range: 0<..<=1." << endl;
					return false;
				}
				if (fnc > 0)
				{
					cout << "\nERROR: There is not possible to use " << CheckTag("SWAP")<<" together with " << CheckKey("FNC-TYPE")<<" other than NONE or with flexible molecular constraint" << endl;
					return false;
				}
				
			}
			else if (key == CheckKey("PAIRS"))
			{
				if (pc++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("PAIRS")<<" definition in the " << CheckTag("SWAP")<<" section!" << endl << "\tThe last value will be used." << endl;
				//Allocating an array to store the swap pairs
				std::vector<bool> swset(ntypes*ntypes, false);
				bool nfound = true;
				if (params.find(CheckKeyValue("PAIRS", "ALL")) != string::npos)
				{
					swset.flip();
					nfound = false;
				}
				else
				{
					//Reading the input line
					string p;
					while (ss >> p)
					{
						if (p.find("-ANY") != string::npos)
						{
							try
							{
								int number = std::stoi(p.substr(0, p.find("-ANY")));
								for (int i = 0; i < ntypes; i++)
								{
									swset[(number - 1)*ntypes + i] = true;
									swset[(number - 1) + i * ntypes] = true;
								}
							}
							catch (const std::invalid_argument& )
							{							
								cout << "\nERROR: Invalid argument error in the following line: " << endl<<"\t"<<*it << endl;
								return false;
							}
							nfound = false;
						}
						else
						{
							std::transform(p.begin(), p.end(), p.begin(), [](char ch) {if (ch == '-') ch = ' '; return ch; });
							stringstream ss2(p);
							int num1 = 0, num2 = 0;
							if (!(ss2 >> num1 >> num2)) return ErrLackValue(*it);
							swset[(num1 - 1)*ntypes + num2 - 1] = true;
							swset[(num2 - 1)*ntypes + num1 - 1] = true;
							if (num1 != num2) nfound = false;
						}
					}
				}
				if (nfound)
				{
					cout << "\nERROR: Although  " << CheckKey("FRACT")<<">0 at " << CheckTag("SWAP")<<", but you did not provided any valid statement under  " << CheckKey("PAIRS")<<" keyword in the following line:" << endl<<"\t"<<*it << endl;
					return false;
				}
				//Writing the swap datasets into array
				//allocating the arrays to hold the type of atoms pairs allowed to be swapped
				if (swap_type1 == 0) SetArraysize(&swap_type1, ntypes*(ntypes - 1) / 2, "swap_type1", "RunParams::GetParams");
				if (swap_type2 == 0) SetArraysize(&swap_type2, ntypes*(ntypes - 1) / 2, "swap_type2", "RunParams::GetParams");
				nswap_pairs = 0;
				for (int i = 0; i < ntypes - 1; i++)
				{
					for (int j = i + 1; j < ntypes; j++)
					{
						if (swset[i*ntypes + j])
						{
							swap_type1[nswap_pairs] = i;
							swap_type2[nswap_pairs++] = j;
						}
					}
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("SWAP"));//write the appropriate error mesage
				return false;
			}
		}
	}
	else
	{
		if ((ntypes == 1) && (pool.size() != 0)) cout << "\nWARNING(" << ++warn << "): There is only one type of atom. Therefore " << CheckTag("SWAP")<<" is not meaningful!" << endl;
		swap_fraction = 0.0;
	}
	return true;
}

#ifdef _LOCAL_INV
bool RunParams::ReadFreeLoc(std::list<string> &pool )
//Process free format strings related to local invariance setup. Returns with false, if incoming datasets are insufficient.
{
	char  *name=NULL;
	int lm = 0, nr = 0, nt=0;
	SetArraysize(&name, NAME_SIZE, "name", "RunParams::ReadFreeLoc");
	nlocint = 0;
	if ((pool.size() > 0) )
	{
		//Check the required keywords
		if (!(CheckKeyInPool(pool, "LOC-INT-FROM_SIGMA")) && (CheckKeyInPool(pool, "LOC-INT-FROM_SIGMA-SCALABLE")) && (CheckKeyInPool(pool, "LOC-INT-FROM_SIGMA-MASTER")))
		{
			cout<<"\nERROR!!: **** One from " <<  CheckKey("LOC-INT-FROM_SIGMA")<<", "<< CheckKey("LOC-INT-FROM_SIGMA-MASTER")<<" or " << CheckKey("LOC-INT-FROM_SIGMA-SCALABLE")<<" is a mandatory item in " << CheckTag("LOCINV") << " section, and it is missing!" << endl;
			return false;
		}
				
		if (!(CheckKeyInPool(pool, CheckKey("LOC-MODE"))))
		{
			//No LOC-MODE is given, set the default
			loc_chi2_mode = LOC_CHI2_MODE_DEF;
			
		}
		
		loc_rspacing = rspacing_def;//default
		for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("LOC-MODE"))
			{
				string temp;
				if (lm++ > 0)
				{
					cout << "\nERROR: Multiple " << CheckKey("LOC-MODE") << " definition in the " << CheckTag("LOCINV") << " section, correct the *.dat file!" << endl;
					return (0);
				}
				if (!(ss >> temp)) return ErrLackValue(*it);
				if (temp == CheckKeyValue("LOC-MODE","BIN-BASED"))
					loc_chi2_mode = 0;
				else
					loc_chi2_mode = 1;
					
				if (loc_chi2_mode < 0 || loc_chi2_mode>1)
				{

					cout << "\nWARNING(" << ++warn << "): The calculation mode for the local chi2 can be 0 (bin based) or 1 (distance base). " << endl;
					cout << "\tIt will be set to the default "<< LOC_CHI2_MODE_DEF<<"!" << endl;
					loc_chi2_mode = LOC_CHI2_MODE_DEF;
				}
				
			}
			else if (key == CheckKey("LOC-R-SPACING"))
			{
				if (nr++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("LOC-R-SPACING") << " definition in the " << CheckTag("LOCINV") << " section!" << endl << "\tThe last value will be used." << endl;
				if (!(ss >>loc_rspacing)) return ErrLackValue(*it);
				if (loc_rspacing < 0)
				{
					cout << "\nWARNING(" << ++warn << "): The local invariance histogram spacing is negative!" << endl;
					cout << "\tIt will be set to default value, which is the default spacing of the normal histogram!" << endl;
					loc_rspacing = rspacing_def;
				}
			}
			else if ((key == CheckKey("LOC-INT-FROM_SIGMA")) || (key == CheckKey("LOC-INT-FROM_SIGMA-MASTER")) || (key == CheckKey("LOC-INT-FROM_SIGMA-SCALABLE")))
			{
				if (nlocint==0)
					SetArraysize(&loc_inv_sigma, 1, "loc_inv_sigma", "RunParams::ReadFreeLoc");
				
				if (loc_chi2_mode == 0)//for the bin-based calculation a minimum distance is given 
				{
					if (nlocint > 0)
					{
						cout << "\nWARNING(" << ++warn << "): Multiple LOC-INT-FROM_SIGMA... definition in the " << CheckTag("LOCINV") << " section!" << endl << "\tThe last value will be used." << endl;
						if (lead_series_ind == ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord + ncommonneigh + nsecondneigh + nlocint + 1)
							lead_series_ind = -1;
						nlocint = 0;//set it back
					}
					if (!(ss >> min_loc_r >> loc_inv_sigma[nlocint])) return ErrLackValue(*it);
				}
				else //(loc_chi2_mode==1)
				{
					//distance-based calculation, the local histograms has to be calculated for all the possibale bins, as we do not know how much the local 
					//environments differ
					//percentage of the atoms will decide, how many atoms (between int(loc_at_ration[i]*natom_iytpe) and int(loc_at_ration[i+1]*natom_iytpe))
					//of a type will be involved in the distance based calculation
					max_loc_r = SQRT3;
					if (nlocint==0)//make sure, that it is created only in case of this is the first interval
						SetArraysize(&loc_at_ratio, 2, "loc_at_ratio", "RunParams::readfree:loc");
					else
					{
						//resize arrays
						ResizeArray(&nlocint,nlocint+1,&loc_inv_sigma, "loc_inv_sigma", "RunParams::ReadFreeLoc");
						ResizeArray(&nlocint, nlocint + 1, &loc_at_ratio, "loc_at_ratio", "RunParams::ReadFreeLoc");//dim: nlocint+1
						nlocint-=2;//to be the array dimension

					}
					if (!(ss >> loc_at_ratio[nlocint] >> loc_inv_sigma[nlocint])) return ErrLackValue(*it);
				}
				
				if (key == CheckKey("LOC-INT-FROM_SIGMA") || key == CheckKey("LOC-INT-FROM_SIGMA-MASTER"))
					loc_inv_sigma[nlocint] = fabs(loc_inv_sigma[nlocint]);
				if (key == CheckKey("LOC-INT-FROM_SIGMA-MASTER"))
				{
					if (lead_series_ind == -1)
						lead_series_ind = ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord + ncommonneigh +nsecondneigh + nlocint+1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry was declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
						cout << ", " << CheckTag("LOCINV");
						cout << " section!" << endl;
						lead_series_ind = -2;//incorrectly set
					}
				}
				else if (key == CheckKey("LOC-INT-FROM_SIGMA-SCALABLE"))
				{
					loc_inv_sigma[nlocint] = -fabs(loc_inv_sigma[nlocint]);
					ChiSquared::calc_sigma = 1;
				}
				nlocint++;
			}
			else if (key == CheckKey("LOC-INT-TO"))
			{
				if (nt++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("LOC-INT-TO") << " definition in the " << CheckTag("LOCINV") << " section!" << endl << "\tThe last value will be used." << endl;
				if (loc_chi2_mode == 0)
				{
					if (!(ss >> max_loc_r)) return ErrLackValue(*it);
				}
				else
					if (!(ss >> loc_at_ratio[nlocint])) return ErrLackValue(*it);
			}
			else
			{
				CheckKeyTag(key, CheckTag("LOCINV"));//write the appropriate error mesage
				return false;
			}
		}
	}
	else
	{
	cout << "\n*****ERROR*****" << endl;
	cout << "The code was compiled with the _LOCAL_INV compiler option to calculate locale invariance," << endl;
	cout << "but no [ LOCINV ] section was found in the " << datfilename << " file!" << endl;
	cout << "Use the normal code, if do not want to calculate local invariance, cannot run this way, exiting..." << endl;
	CleanExit();
	}
	//make some cheks
	if (loc_chi2_mode == 0)
	{//for the bin-based calculation a minimum and a maximum distance can be given 

		if (min_loc_r < 0 || max_loc_r < 0)
		{
			if (min_loc_r < 0)
			{
				mystrcpy(name, NAME_SIZE, "minimum");
				if (max_loc_r < 0)
					mystrcat(name, NAME_SIZE, " and maximum");
				mystrcat(name, NAME_SIZE, " distance for the local invariance calculation");
			}
			else
			{
				if (max_loc_r < 0)
					mystrcpy(name, NAME_SIZE, "maximum distance for the local invariance calculation");

			}

			cout << "\n*****ERROR*****" << endl;
			cout << "The " << name << " in the " << datfilename << " file is negative!" << endl;
			cout << "It has to be in the range of the normal histogram, between " << binshift_def * rspacing[0] / boxedge;
			cout << " and " << xmax[0] << " in reduced units! Set it right, and try again! Exiting..." << endl;
			CleanExit();
		}
	}
	else //(loc_chi2_mode==1)
	{
		//distance-based calculation, the local histograms has to be calculated for all the possibale bins, as we do not know how much the local 
		//environments differ
		//percentage of the atoms will decide, how many atoms (between int(loc_at_ration[i]*natom_iytpe) and int(loc_at_ration[i+1]*natom_iytpe))
		//of a type will be involved in the distance based calculation
		max_loc_r = SQRT3;

		if (loc_at_ratio[0] < 0 || loc_at_ratio[0] >= 1.0)
		{
			cout << "\nWARNING(" << ++warn << "): The ratio of atoms to start the distance based calculation of the local invariance has to be greater or equl than 0," << endl;
			cout << "\tsmaller than 1.0! It will be set to 0.0" << endl;
			loc_at_ratio[0] = 0.0;
		}
		for (int i = 1; i < nlocint + 1; i++)
		{
			if (loc_at_ratio[i] <= 0 || loc_at_ratio[i] > 1.0 || loc_at_ratio[i] < loc_at_ratio[i - 1])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The ratio of atoms for the end of the " << i << ". interval in the distance based calculation of " << endl;
				cout << "the local invariance (" << loc_at_ratio[i] << ") has to be greater than 0 and <=1.0 " << endl;
				cout << "and greater than the previous ratio (" << loc_at_ratio[i - 1] << ")" << endl;
				cout << "Change the " << datfilename << " file and try again! Exiting..." << endl;
				CleanExit();
			}
		}
	}
	if (name != NULL)
		delete[] name;

	return true;
}
#endif

#ifdef _VIBR_AMP
bool RunParams::ReadFreeVibramp(std::list<string> &pool )
//Process free format string related to vibrational amplitude setup. Returns with false, if incoming datasets are insufficient.
{
	char  *name = NULL;
	SetArraysize(&name, NAME_SIZE, "name", "RunParams::ReadFreevibramp");

	if ((pool.size() > 0))
	{
		//Check the required keywords
		if (!(CheckKeyInPool(pool, CheckKey("GAUSSIAN-SIGMA"))))
		{
			cout << "\nERROR!!: **** " << CheckKey("GAUSSIAN-SIGMA") << " is a mandatory item in " << CheckTag("VIBRAMP") << " section, and it is missing!" << endl;
			return false;
		}

		if (std::count_if(pool.begin(), pool.end(), [](string c) { return (c.find(CheckKey("GAUSSIAN-SIGMA")) != string::npos); }) != 1)
		{
			cout << "\nWARNING(" << ++warn << "): Multiple declaration of "<< CheckKey("GAUSSIAN-SIGMA")<<" has been found in " << CheckTag("VIBRAMP") << " section." << endl << "\tThe last one will be used!" << endl;
		}

		SetArraysize(&T_corr_sigma_t, ntypes, "T_corr_sigma_t", "RunParams::ReadFreevibramp");
		for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("GAUSSIAN-SIGMA"))
			{
				for (int i = 0; i < ntypes; i++)
				{
					if (!(ss >> T_corr_sigma_t[i])) return ErrLackValue(*it);

					if (T_corr_sigma_t[i] * 4 < rspacing[0] / 2)
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The sigma parameter for thermal correction (" << T_corr_sigma_t[i] << ") is too small compared to the bin size." << endl;
						cout << "The 4*T_corr_sigma_t accounting for 99.99% of the values are inside one histgram bin, so" << endl;
						cout << "no thermal correction will be performed!!! Either change the bin size or the T_sigma_corr" << endl;
						cout << "value! Exiting..." << endl;
						CleanExit();
					}
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("VIBRAMP"));//write the appropriate error mesage
				return false;
			}
		}
	}
	else
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The code was compiled with the _VIBR_AMP compiler option to include the vibrational motion of the atoms," << endl;
		cout << "but no [ VIBRAMP ] section was found in the " << datfilename << " file!" << endl;
		cout << "Use the normal code, if do not want vibrational motion, cannot run this way, exiting..." << endl;
		CleanExit();
	}
	return true;
}
#endif

#ifdef _NO_PERIODIC
bool RunParams::ReadFreeNoper(std::list<string> &pool )
//Process free format string related to no periodic boundary condition setup. Returns with false, if incoming datasets are insufficient.
{
	char  *name = NULL;
	int i;
	SetArraysize(&name, NAME_SIZE, "name", "RunParams::ReadFreenoper");

	if ((pool.size() > 0))
	{
		//Check the required keywords
		if (!(CheckKeyInPool(pool, CheckKey("SAMPLE-RADIUS"))))
		{
			cout << "\nERROR!!: **** " << CheckKey("SAMPLE-RADIUS") << " is a mandatory item in " << CheckTag("NOPER") << " section, and it is missing!" << endl;
			return false;
		}

		if (std::count_if(pool.begin(), pool.end(), [](string c) { return (c.find(CheckKey("SAMPLE-RADIUS")) != string::npos); }) != 1)
		{
			cout << "\nWARNING(" << ++warn << "): Multiple declaration of "<<CheckKey("SAMPLE-RADIUS")<<" has been found in " << CheckTag("NOPER") << " section." << endl << "\tThe last one will be used!" << endl;
		}

		for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("SAMPLE-RADIUS"))
			{
				if (!(ss >> R0)) return ErrLackValue(*it);

				if (fabs(xmax[0]- 0.5)<TOLERANCE2)
				{
					cout << "\nWARNING(" << ++warn << "): xmax has to be 0.5 in case of non-periodic boundary conditions!" << endl;
					cout << "\txmax will be reset to 0.5!" << endl;
					xmax[0] = 0.5;

				}
					
				//check, if fnc>0 was given
				if (RunParams::fnc > 0)
				{
					cout << "\n***** ERROR *****" << endl;
					cout << "FNC constraint or flexible molecules cannot be used in case of non-periodic boundary conditions!" << endl;
					cout << "Set fnc option to 0, and try again! exiting" << endl;
					CleanExit();
				}

				double v_R0 = 4.0 / 3.0*PI*pow(R0, 3);
				double dtemp = SimpleCfg::ntotal / v_R0;
				if (fabs(dtemp - rho) > 5e-7)
				{
					cout << "\nWARNING(" << ++warn << "): The radius of the sample (R0) does not correspond to the density, the R0 will be reset" << endl;
					cout << "\tfrom " << R0;
					R0 = pow((SimpleCfg::ntotal / rho)*3.0 / 4.0 / PI, 1.0 / 3.0);
					cout << " to " << R0 << " A!" << endl;
				}

				//We have to make sure, that the R0 is divisible by dr to make the volume correction easier
				dtemp = R0 / rspacing_def;
				int temp = (int)dtemp;
				if (fabs(dtemp - (double)temp - binshift_def) > TOLERANCE)
				{
					cout << "\nWARNING(" << ++warn << "): The radius of the sample is not divisible by the bin width" << endl;
					cout<<"\ttaking into account the binshift as well!" << endl;
					cout << "\tThe binwidth will be reset from " << rspacing_def << " to ";
					if (dtemp - temp - binshift_def > 0.5)
					{
						rspacing_def = R0 / (temp + binshift_def + 1);//make it a bit smaller
						temp++;
					}
					else
						rspacing_def = R0 / (temp + binshift_def);//make it a bit larger
					cout << rspacing_def << endl;
					for (i = 0; i < ntot_datasets; i++)
						rspacing[i] = rspacing_def;
				}
				nbins[0] = temp;
				if (fabs(SimpleCfg::boxedge / 2.0 - R0) > 1e-7)
				{
					cout << "\nWARNING(" << ++warn << "): The total boxlength has to be four times the sample radius, the half  " << endl;
					cout << "\tboxlength will be reset from " << SimpleCfg::boxedge << " to ";
					SimpleCfg::boxedge = R0 * 2;
					boxedge = SimpleCfg::boxedge;
					cout << SimpleCfg::boxedge << endl;
					for (i = 0; i < navcoord; i++)
					{
						//Reset minimum and maximum distances for the constraints squared in reduced units
						AvCoordConst::udminsq[i] = AvCoordConst::dmin[i] * AvCoordConst::dmin[i] / boxedge / boxedge;
						AvCoordConst::udmaxsq[i] = AvCoordConst::dmax[i] * AvCoordConst::dmax[i] / boxedge / boxedge;


					}
					for (i = 0; i < nicoord; i++)
					{
						//Initialising the remainder of the static members
						//Minimum and maximum distances for the constraints squared in reduced units
						for (int ineightype = 0; ineightype < CoordNumbConst::n_neightype[i]; ineightype++)
						{
							CoordNumbConst::udminsq[CoordNumbConst::cum_n_neightype[i] + ineightype] = \
								CoordNumbConst::dmin[CoordNumbConst::cum_n_neightype[i] + ineightype] * CoordNumbConst::dmin[CoordNumbConst::cum_n_neightype[i] + ineightype] / boxedge / boxedge;
							CoordNumbConst::udmaxsq[CoordNumbConst::cum_n_neightype[i] + ineightype] = \
								CoordNumbConst::dmax[CoordNumbConst::cum_n_neightype[i] + ineightype] * CoordNumbConst::dmax[CoordNumbConst::cum_n_neightype[i] + ineightype] / boxedge / boxedge;

						}

					}
				}
				//this has to be set, regardless if the boxedge was changed
				for (i = 0; i < navcoord; i++)
				{
					AvCoordConst::dav_red[i] = (AvCoordConst::dmax[i] + AvCoordConst::dmin[i]) / 2 / boxedge;//reduced maximum distance
					AvCoordConst::ncentral[i] = static_cast<int>(round(AvCoordConst::ncentral[i] * (pow(R0 - (AvCoordConst::dav_red[i] * boxedge), 3) / pow(R0, 3))));//number of central atoms, which coordination atom is inside the sample
					if (AvCoordConst::ncentral[i] == 0)
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The number of available central atoms for the " << i + 1 << ". average coordinaion constraint data set is zero" << endl;
						cout << "Check the constraint parameters in the " << datfilename << " file!" << endl;
						cout << "Cannot run this way, exiting..." << endl;
						CleanExit();
					}
				}
			}
			else if (key == CheckKey("RECENTRE-SPHERE"))
			{
				if (!(ss >> recentre_flag)) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "RECENTRE-SPHERE");
				Check0_1(recentre_flag, name, "RunParams::ReadFreenoper");
				//if no periodic boundary conditions are used, then the sample have to be in the middle of the box not exceeding 
				//the quarter of the box edge from the centre in each coordinate axis direction
			}
			else
			{
				CheckKeyTag(key, CheckTag("NOPER"));//write the appropriate error mesage
				return false;
			}
		}
	}
	return true;
}
#endif

bool RunParams::ReadFreeNbpot(std::list<string> &pool )
//Process free format strings related to the general potential setup. Returns with false, if incoming datasets are insufficient.
{
	int npartials = ntypes * (ntypes + 1) / 2;
	int nGRpartials=0;
	int my_vdW_weight_mode=-10;//needed if both VDW-SIGMA and COUL-SIGMA is read
	int my_Coul_weight_mode = -10;
	string key, params;
	bool noCoulomb = true;
	if (pool.size() > 0)
	{
		//int sop=0;
		NB_weight_mode = -10;//was not set yet
		lead_series_ind2 = -1;//was not set yet, default will not be set in this routine, if none is given, only in GetFreeFormatParams, there can be lead bpot
		//Set the default values for the optional parameters, if not given
		if (!CheckKeyInPool(pool, CheckKey("NB-TYPE")))
		{
			potential = POT_TYPE_DEF;
		}
		else
		{
			int pott = 0;
			for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
			{
				
				WrapFree(*it, key, params);
				stringstream ss(params);


				if (key == CheckKey("NB-TYPE"))//this has to be known to check for the others
				{
					string temp;
					if (pott++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("NB-TYPE") << " definitions in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> temp)) return ErrLackValue(*it);
					if (temp == CheckKeyValue("NB-TYPE","LJ")) potential = 1;
					else if (temp == CheckKeyValue("NB-TYPE", "TABULATED")) potential = 10;
					else
					{
						cout << "\nERROR: only the following potential types are implemented: " << GiveKeyValues("NB-TYPE") << endl;
						return false;
					}
					
				}
			}
		}
		if (potential == 1)
		{
			mystrcpy(vdW_name, 15, "LJ");
			if (!CheckKeyInPool(pool, CheckKey("NGR-TYPES")))
			{
				nGRtypes = ntypes;
				nGRpartials = nGRtypes * (nGRtypes + 1) / 2;
				SetArraysize(&vdW_pot1, nGRpartials, "vdW_pot1", "RunParams::GetPotParams");//if LJ, this is the sigma parameter of the potential, describing where the pot crosses the x-axis
				SetArraysize(&vdW_pot2, nGRpartials, "vdW_pot2", "RunParams::GetPotParams");//if LJ, epsilon is the depth of the potential valley
			}
			if (!CheckKeyInPool(pool, CheckKey("CHI-LOWLIMIT-RATIO"))) pot_chi2_low_lim_fraction = POT_CHI_LOW_LIM_FRACTION_DEF;
			if (!CheckKeyInPool(pool, CheckKey("LJ-POW"))) LJ_rep_N = LJ_REP_N_DEF;
			if (!CheckKeyInPool(pool, CheckKey("1-4-SCALE-VDW_COUL")))
			{
				vdW14_fudge = FUDGE14_DEF;
				Coulomb14_fudge = FUDGE14_DEF;
			}
			if (!CheckKeyInPool(pool, CheckKey("CUTOFF-VDW_COUL")))
			{
				vdW_cutoff = CUTOFF_DEF * SimpleCfg::boxedge;
				Coulomb_cutoff = vdW_cutoff;
			}
			int lc = 0, lp = 0, v14 = 0, co = 0, ng = 0;

			//Check the presence of required keywords connected to the weight 
			int vdwc = 0, csc = 0, soc = 0;
		
			for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
			{
				if (it->find("VDW-SIGMA") != string::npos) vdwc++;
				if (it->find("COUL-SIGMA") != string::npos)
				{
					csc++;
					noCoulomb = false;
				}
				if (it->find("SIGMA-OVERALLPOT") != string::npos) soc++;
			}
			
			if ((vdwc > 1) || (csc > 1) || (soc > 1))
			{
				cout << "\nERROR!!: Multiple definition of ****" << CheckKey("VDW-SIGMA") << "**** or ****" << CheckKey("COUL-SIGMA") << "**** or ****" << CheckKey("SIGMA-OVERALLPOT") << " in " << CheckTag("NBPOT") << " section!." << endl;
				cout << "Please provide one " << CheckKey("SIGMA-OVERALLPOT") << " or a pair of " << CheckKey("VDW-SIGMA") << " and " << CheckKey("COUL-SIGMA") << "!" << endl;
				return false;
			}
			if (!(((vdwc == 1) && (csc == 1) && (soc == 0)) || ((vdwc == 1) && (csc ==0) && (soc == 0)) || ((vdwc == 0) && (csc == 0) && (soc == 1))))
			{
				cout << "\nERROR!!: ****" << CheckKey("VDW-SIGMA") << "**** and ****" << CheckKey("COUL-SIGMA") << "****, or ****" << CheckKey("SIGMA-OVERALLPOT") << "*** are mandatory items in " << CheckTag("NBPOT") << " section!." << endl;
				cout << "Please provide one " << CheckKey("SIGMA-OVERALLPOT") << ",or one "<< CheckKey("VDW-SIGMA")<<" or\n a pair of " << CheckKey("VDW-SIGMA") << " and " << CheckKey("COUL-SIGMA") << "!" << endl;
				return false;
			}
			
			//hopefully all the incorrect possibilities were covered

			//check the LJ-parameters
			int sige = 0, c6 = 0, siga = 0, sigg = 0, epse = 0, cpow = 0, epsg = 0;
			for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
			{
				if (it->find("LJ-C6-EACH") != string::npos) sige++;
				if (it->find("LJ-C6-GMEAN") != string::npos) c6++;
				if (it->find("LJ-SIG-AMEAN") != string::npos) siga++;
				if (it->find("LJ-SIG-GMEAN") != string::npos) sigg++;
				if (it->find("LJ-CPOW-EACH") != string::npos) epse++;
				if (it->find("LJ-CPOW-GMEAN") != string::npos) cpow++;
				if (it->find("LJ-EPS-GMEAN") != string::npos) epsg++;
			}

			if (!(((sige == 1) && (epse == 1)) || ((c6 == 1) && (cpow == 1)) || ((siga == 1) && (epsg == 1)) || ((sigg == 1) && (epsg == 1))))
			{
				cout << "\nERROR!!: Only the following combinations are accepted for the LJ parameters in the " << CheckTag("NBPOT") << " section:" << endl;
				cout << " One keyword of each " << CheckKey("LJ-C6-EACH") << " and " << CheckKey("LJ-CPOW-EACH") << " for combination rule: 0" << endl;
				cout << " One keyword of each " << CheckKey("LJ-C6-GMEAN") << " and " << CheckKey("LJ-CPOW-GMEAN") << " for combination rule: 1" << endl;
				cout << " One keyword of each " << CheckKey("LJ-SIG-AMEAN") << " and " << CheckKey("LJ-EPS-GMEAN") << " for combination rule: 2" << endl;
				cout << " One keyword of each " << CheckKey("LJ-SIG-GMEAN") << " and " << CheckKey("LJ-EPS-GMEAN") << " for combination rule: 3" << endl;
				cout << " Check the manual and correct the " << datfilename << " file!" << endl;
				return false;

			}
			//now only a valid option is possible, set the combination rule
			if (CheckKeyInPool(pool, CheckKey("LJ-C6-EACH")) && CheckKeyInPool(pool, CheckKey("LJ-CPOW-EACH"))) vdW_comb_rule = 0;
			if (CheckKeyInPool(pool, CheckKey("LJ-C6-GMEAN")) && CheckKeyInPool(pool, CheckKey("LJ-CPOW-GMEAN")))	vdW_comb_rule = 1;
			if (CheckKeyInPool(pool, CheckKey("LJ-SIG-AMEAN")) && CheckKeyInPool(pool, CheckKey("LJ-EPS-GMEAN"))) vdW_comb_rule = 2;
			if (CheckKeyInPool(pool, CheckKey("LJ-SIG-GMEAN")) && CheckKeyInPool(pool, CheckKey("LJ-EPS-GMEAN"))) vdW_comb_rule = 3;

			SetArraysize(&vdW_weight, npartials, "vdW_weight", "RunParams::GetPotParams");
			SetArraysize(&Coulomb_weight, npartials, "Coulomb_weight", "RunParams::GetPotParams");
		
			for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
			{
				string key, params;
				WrapFree(*it, key, params);
				stringstream ss(params);

				if (key == CheckKey("NB-TYPE")) {}//this was cheked, has to be included, to show that this is a valid key
				else if (key == CheckKey("SIGMA-OVERALLPOT"))//NB_weight_mode = 2
				{
					lead_series_ind2 = RegularizeSigmaEntry(params);
					if (lead_series_ind2 < 1)
					{
						cout << "\nWARNING(" << ++warn << ") : multiple definition of \"m\" or \"M\" or no leading potential index at all was found in the following line:" << endl<<"\t"<<*it << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					ss.str(""); ss.clear();
					ss << params;
					if (!(ss >> vdW_weight[0])) return ErrLackValue(*it);
					CheckSigma(&vdW_weight[0], "SIGMA parameter for non-bonded potential interactions");//check, if it not zero
					if (vdW_weight[0] < 0)
					{
						ChiSquared::calc_sigma = 1;
						Coulomb_weight[0] = -1e6;//scalable, will be set later
					}
					else Coulomb_weight[0] = vdW_weight[0];//fixed weight, can be set now

					for (int i = 1; i < npartials; i++)
					{
						vdW_weight[i] = vdW_weight[0];
						Coulomb_weight[i] = Coulomb_weight[0];
					}
					NB_weight_mode = 2;
				}
				else if (key == CheckKey("VDW-SIGMA"))//NB_weight_mode = 0 or 1, depending on how many value is given
				{
					int mylead = RegularizeSigmaEntry(params);
					if (lead_series_ind2 > -1 && mylead > 0)//already a valid index was found among the Coulomb sigmas
					{
						cout << "\nWARNING(" << ++warn << ") : multiple leading series index definitions in the " << CheckKey("VDW-SIGMA") << " and " << CheckKey("COUL-SIGMA") << " lines!" << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					else
					{
						if (mylead < 0)
						{
							cout << "\nWARNING(" << ++warn << ") : multiple definition of \"m\" or \"M\" of leading potential index was found in the following line:" << endl<<"\t"<<*it << endl;
							lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
						}
						else
						{
							if (mylead > 0)
								lead_series_ind2 = mylead;//this is a valid lead series index
						}
					}

					ss.str(""); ss.clear();
					ss << params;
					int last = -1;
					my_vdW_weight_mode = 0;//deafult only one parameter is needed, if more than one given, even if not npartials, then for the missing ones the last valid value is used
					bool rv = true;
					for (int i = 0; i < npartials; i++)
					{
						if (rv) rv = static_cast<bool>(ss >> vdW_weight[i]);
						if (!rv)
						{
							if (last == -1)
							{
								cout << "\nERROR: too few parameters. You should provide at least one value in the following line: " << endl<<"\t"<<*it << endl;
								return false;
							}
							else
								vdW_weight[i] = vdW_weight[last];

						}
						else
						{
							CheckSigma(vdW_weight + i, "VDW-SIGMA parameter");//check, if it is not zero
							if (vdW_weight[i] < 0) ChiSquared::calc_sigma = 1;
							if (i > 0) my_vdW_weight_mode = 1;
							last = i;
						}
					}

					if ((last != npartials - 1) && (my_vdW_weight_mode == 1) && ::debug)
						cout << "\nNOTE(" << ++note << "): VDW-SIGMA parameter from partial " << last + 2 << " to partial " << npartials << " are/is adjusted to " << vdW_weight[last] << endl;
				}
				else if (key == CheckKey("COUL-SIGMA"))//NB_weight_mode = 0 or 1, depending on how many value is given
				{
					int mylead = RegularizeSigmaEntry(params);
					if (lead_series_ind2 > -1 && mylead > 0)//already a valid index was found among the vdW sigmas
					{
						cout << "\nWARNING(" << ++warn << ") : multiple leading series index definitions in the " << CheckKey("VDW-SIGMA") << " and " << CheckKey("COUL-SIGMA") << " lines!" << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					else
					{
						if (mylead < 0)
						{
							cout << "\nWARNING(" << ++warn << ") : multiple definition of \"m\" or \"M\" of leading potential was found in the following line:" << endl<<"\t"<<*it << endl;
							lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
						}
						else
						{
							if (mylead > 0)
								lead_series_ind2 = mylead;//this is a valid lead series index
						}
					}

					ss.str(""); ss.clear();
					ss << params;
					int last = -1;
					my_Coul_weight_mode = 0;//deafult only one parameter is needed, if more than one given, even if not npartials, then for the missing ones the last valid value is used
					bool rv = true;
					for (int i = 0; i < npartials; i++)
					{
						if (rv) rv = static_cast<bool>(ss >> Coulomb_weight[i]);
						if (!rv)
						{
							if (last == -1)
							{
								cout << "\nERROR: too few parameters. You should provide at least one value in the following line: " << endl<<"\t"<<*it << endl;
								return false;
							}
							else
								Coulomb_weight[i] = Coulomb_weight[last];

						}
						else
						{
							CheckSigma(vdW_weight + i, "COUL-SIGMA parameter");//check, if it is not zero
							if (Coulomb_weight[i] < 0) ChiSquared::calc_sigma = 1;
							if (i > 0) my_Coul_weight_mode = 1;
							last = i;
						}
					}

					if ((last != npartials - 1) && (my_Coul_weight_mode == 1) && ::debug)
						cout << "\nNOTE(" << ++note << "): COUL-SIGMA parameter from partial " << last + 2 << " to partial " << npartials << " are/is adjusted to " << vdW_weight[last] << endl;
				}
				else if (key == CheckKey("NGR-TYPES"))//if it was given, this is the first in the list, if not given set and the arrays created at the beginning of the routine
				{
					if (ng++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("NGR-TYPES") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> nGRtypes)) return ErrLackValue(*it);
					nGRpartials = nGRtypes * (nGRtypes + 1) / 2;
					SetArraysize(&vdW_pot1, nGRpartials, "vdW_pot1", "RunParams::GetPotParams");//if LJ, this is the sigma parameter of the potential, describing where the pot crosses the x-axis
					SetArraysize(&vdW_pot2, nGRpartials, "vdW_pot2", "RunParams::GetPotParams");//if LJ, epsilon is the depth of the potential valley

				}
				else if (key == CheckKey("LJ-C6-EACH"))
				{
					for (int i = 0; i < nGRpartials; i++)
						if (!(ss >> vdW_pot1[i])) return ErrLackValue(*it);
				}
				else if (key == CheckKey("LJ-CPOW-EACH"))
				{
					for (int i = 0; i < nGRpartials; i++)
						if (!(ss >> vdW_pot2[i])) return ErrLackValue(*it);
				}
				else if ((key == CheckKey("LJ-C6-GMEAN")) || (key == CheckKey("LJ-SIG-AMEAN")) || (key == CheckKey("LJ-SIG-GMEAN")))
				{
					//parameters for the different types will be given, and the parameters for the mixed partials will be calculated according to the 
					//combination rule
					for (int i = 0; i < nGRtypes; i++)
					{
						for (int j = i; j < nGRtypes; j++)
						{
							if (i == j)//only read into the ii partials
							{
								if (!(ss >> vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)])) return ErrLackValue(*it);
							}
						}
					}
				}
				else if ((key == CheckKey("LJ-CPOW-GMEAN")) || (key == CheckKey("LJ-EPS-GMEAN")))
				{
					//parameters for the different types will be given, and the parameters for the mixed partials will be calculated according to the 
					//combination rule
					for (int i = 0; i < nGRtypes; i++)
					{
						for (int j = i; j < nGRtypes; j++)
						{
							if (i == j)//only read into the ii partials
							{
								if (!(ss >> vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)])) return ErrLackValue(*it);
							}
						}
					}
				}
				else if (key == CheckKey("CUTOFF-VDW_COUL"))
				{
					if (co++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("CUTOFF-VDW_COUL") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> vdW_cutoff)) return ErrLackValue(*it);
					bool rv = static_cast<bool>(ss >> Coulomb_cutoff);
					if (!rv)
					{
						if (::debug)
						{
							cout << "\nWARNING(" << ++warn << " cutoff limit for the Coulomb interactions was not given in line: " << endl<<"\t"<<*it << endl;
							cout << "\nThe vdW cutoff be used!" << endl;
						}
						Coulomb_cutoff = vdW_cutoff;
					}
				}
				else if (key == CheckKey("CHI-LOWLIMIT-RATIO"))
				{
					if (lc++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("CHI-LOWLIMIT-RATIO") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (params.find("NONE") == std::string::npos)
						if (!(ss >> pot_chi2_low_lim_fraction)) return ErrLackValue(*it);
				}
				else if (key == CheckKey("LJ-POW"))
				{
					if (lp++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("LJ-POW") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> LJ_rep_N)) return ErrLackValue(*it);
				}
				else if (key == CheckKey("1-4-SCALE-VDW_COUL"))
				{
					if (v14++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("1-4-SCALE-VDW_COUL") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> vdW14_fudge)) return ErrLackValue(*it);
					bool rv = static_cast<bool>(ss >> Coulomb14_fudge);
					if (!rv)
					{
						if (::debug)
						{
							cout << "\nWARNING(" << ++warn << " scaling factor for 1-4 Coulomb interactions was not given in line: " << endl<<"\t"<<*it << endl;
							cout << "\nThe scaling factor for the 1-4 vdW interactions will be used!" << endl;
						}
						Coulomb14_fudge = vdW14_fudge;
					}
				}
				else
				{
					CheckKeyTag(key, CheckTag("NBPOT"), CheckKeyValue("NB-TYPE", "LJ"));//write the appropriate error mesage
					return false;
				}

			}
		}
		else if (potential==10)
		{
		mystrcpy(vdW_name, 15, "Tabulated");
			string myfile_names;
			NB_weight_mode = 1;//different sigma for all the partials

			if (!CheckKeyInPool(pool, CheckKey("TEMPERATURE"))) temperature = TEMPERATURE_DEF;
			if (!CheckKeyInPool(pool, CheckKey("CUTOFF-TAB"))) vdW_cutoff = CUTOFF_DEF * SimpleCfg::boxedge;
			int t = 0; int co = 0;
			for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
			{
				string key, params, mysigma,dummy;
				WrapFree(*it, key, params);
				stringstream ss(params);
				

				if (key == CheckKey("NB-TYPE")) {}//this was cheked, has to be included, to show that this is a valid key
				else if (key == CheckKey("TEMPERATURE"))
				{
					if (t++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("TEMPERATURE") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> temperature)) return ErrLackValue(*it);
				}
				else if (key == CheckKey("CUTOFF-TAB"))
				{
					if (co++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("CUTOFF-TAB") << " definition in the " << CheckTag("NBPOT") << " section!" << endl << "\tThe last value will be used." << endl;
					if (!(ss >> vdW_cutoff)) return ErrLackValue(*it);
				}
				else if ((key == CheckKey("IND_SIGMA_FILE")) || (key == CheckKey("IND_SIGMA-MASTER_FILE")) || (key == CheckKey("IND_SIGMA-SCALABLE_FILE")))
				{
					if (nused_potpartials == 0)
					{
						SetArraysize(&vdW_weight, 1, "vdW_weight", "RunParams::ReadFreeNbpot");
						SetArraysize(&FNC_POT::tabpot_index, 1, "tabpot_index", "RunParams::ReadFreeNbpot");
					}
					else
					{
						ResizeArray(&nused_potpartials, nused_potpartials+1, &vdW_weight, "vdW_weight", "RunParams::ReadFreeNbpot");
						nused_potpartials--;
						ResizeArray(&nused_potpartials, nused_potpartials+1, &FNC_POT::tabpot_index, "FNC_POT::tabpot_index", "RunParams::ReadFreeNbpot");
						nused_potpartials--;
					}
					if (!GetInt(ss, *it, FNC_POT::tabpot_index+nused_potpartials,"RunParams::ReadFreeNbpot")) return ErrLackValue(*it);
					ss >> mysigma;
					int mylead = RegularizeSigmaEntry(mysigma);
					if (mylead > 0)
						mylead = nused_potpartials + 1;//it has to be set to the index of this tabulated potential, as the regularize... routine is not exactly for this purpose
					if (lead_series_ind2 > -1 && mylead > 0)//already a valid index was found among the other tabulated potentials
					{
						cout << "\nWARNING(" << ++warn << ") : multiple leading series index definitions in the " << CheckKey("IND_SIGMA_FILE") << " and " << CheckKey("IND_SIGMA-MASTER_FILE") << " and " << CheckKey("IND_SIGMA-SCALABLE_FILE") << " lines!" << endl;
						//tabulated potential is not recommeneded to be used with lexible molecules, but anyway
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					else
					{
						if (mylead < 0)//this should not happen here
						{
							cout << "\nWARNING(" << ++warn << ") : multiple definition of \"m\" or \"M\" of leading potential index was found in the following line:" << endl<<"\t"<<*it << endl;
							lead_series_ind2 = -2;//Indicating, tat it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
						}
						else
						{
							if (mylead > 0)
								lead_series_ind2 = mylead;//this is a valid lead series index
						}
					}
					//first finish this line, read the filename
					if (!(ss >>dummy)) return ErrLackValue(*it);
					myfile_names += " ";
					myfile_names+=dummy;
					//get the sigma
					ss.str(""); ss.clear();
					ss << mysigma;
					if (!(ss >> vdW_weight[nused_potpartials])) return ErrLackValue(*it);
					CheckSigma(&vdW_weight[nused_potpartials], "SIGMA parameter for non-bonded potential interactions");//check, if it not zero
					if (vdW_weight[nused_potpartials] < 0)
						ChiSquared::calc_sigma = 1;
					nused_potpartials++;
				}
				else
				{
					CheckKeyTag(key, CheckTag("NBPOT"), CheckKeyValue("NB-TYPE", "TABULATED"));//write the appropriate error mesage
					return false;
				}

			}
			FNC_POT::tabpot_filename = new char[nused_potpartials][FILE_NAME_SIZE];
			stringstream ss(myfile_names);
			for (int i = 0;i<nused_potpartials;i++)
			ss>>FNC_POT::tabpot_filename[i];
		}	//end of if tabpot
	}
	else
	{
		potential = POTENTIAL_DEF;
		lead_series_ind2 = -1;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
		NB_weight_mode = NB_WEIGHT_MODE_DEF;
	}
	
	//Set weight mode
	if (potential == 1)
	{
		if (fabs(NB_weight_mode) != 2)
		{
			if ((my_vdW_weight_mode > -10))//it was set
			{
				if (my_Coul_weight_mode > -10)//both was set
				{
					if (my_vdW_weight_mode == my_Coul_weight_mode)//the same weight mode was found for both
						NB_weight_mode = my_vdW_weight_mode;
					else
					{
						//not the same, weight mode 3 (same as 0) will be set
						NB_weight_mode = 3;//this can be negated
						for (int i = 0; i < npartials; i++)
						{
							vdW_weight[i] = vdW_weight[0];
							Coulomb_weight[i] = Coulomb_weight[0];
						}
					}
				}
				else//only vdW was set, set Coulomb accordingly
				{
					NB_weight_mode = my_vdW_weight_mode;
					for (int i = 0; i < npartials; i++)
						Coulomb_weight[i] = vdW_weight[i];
				}
			}
		}
		if (noCoulomb)
		{
			for (int i = 0; i < npartials; i++)
				Coulomb_weight[i] = vdW_weight[i];
		}
		if (pot_chi2_low_lim_fraction<2.0)
			NB_weight_mode *= -1;//low limit should be used
		if (vdW_comb_rule!=0)
		{
			//Now create the parameters for the mixed partials according to the combination rule
			for (int i = 0; i < nGRtypes; i++)
			{
				for (int j = i; j < nGRtypes; j++)
				{
					if (i != j)//mixed partials
					{
						switch (vdW_comb_rule)
						{
						case 1://C6 and CN is given for LJ
						case 3:
						{
							//sigma and epsilon is given, geometric average
							vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)] = sqrt(vdW_pot1[i*nGRtypes + i - (i*(i + 1) / 2)] * vdW_pot1[j*nGRtypes + j - (j*(j + 1) / 2)]);
							vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)] = sqrt(vdW_pot2[i*nGRtypes + i - (i*(i + 1) / 2)] * vdW_pot2[j*nGRtypes + j - (j*(j + 1) / 2)]);
							break;
						}
						case 2:
						{
							//sigma nad epsilon is given, arithmetic average
							vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)] = 0.5*(vdW_pot1[i*nGRtypes + i - (i*(i + 1) / 2)] + vdW_pot1[j*nGRtypes + j - (j*(j + 1) / 2)]);
							vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)] = sqrt(vdW_pot2[i*nGRtypes + i - (i*(i + 1) / 2)] * vdW_pot2[j*nGRtypes + j - (j*(j + 1) / 2)]);
						}
						}
					}
				}//end of j cycle
			}//end of i cycle
		}//end of vdW_comb_rule!=0
	}
	else if (potential==10)
	{
		FNC_POT::nused_potpartials = nused_potpartials;
		FNC_POT::ReadTabPotFree();
	}
	
	return true;
}


bool RunParams::ReadFreeBpot(std::list<string> &pool )
//Process free format strings related to the bonded potential setup. Returns with false, if incoming datasets are insufficient.
{
	n_top_def_option = 0;
	if (fnc == -4)
		fnc = 4;//BPOT tag was given, flexible molecules should be used
	//lead_series_ind2 will not be set to default only in GetFreeFormatParams, as NBPOT can be present
	//[BPOT] tag was used, then fnc was set to for already, and flexible molecules will be attemted to be used, regardless if there is
	//no additional keywords here, as these are all optional!
	if (pool.size() > 0)
	{
		int count = 0;
		if (CheckKeyInPool(pool, CheckKey("LEAD-BOND"))) count++;
		else if (CheckKeyInPool(pool, CheckKey("LEAD-ANGLE"))) count++;
		else if (CheckKeyInPool(pool, CheckKey("LEAD-PERDIH"))) count++;
		else if (CheckKeyInPool(pool, CheckKey("LEAD-HARMDIH"))) count++;
		else if (CheckKeyInPool(pool, CheckKey("LEAD-RBDIH"))) count++;
		if (count > 1)
		{
			cout << "\nWARNING(" << ++warn << ") : Multiple definition of leading potential series index in section " << CheckTag("BPOT") << "!" << endl;
			lead_series_ind2 = -2;//was not set correctly
		}
		for (std::list<std::string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			std::string key, params;
			WrapFree(*it, key, params);
			std::stringstream ss(params);

			if (key == CheckKey("COMP-OPTION"))
			{
				if (params.find("NONE") == std::string::npos)
				{
					std::string input;
					while (ss >> input)
					{
						input.insert(0, "-D");
						if (n_top_def_option > max_option - 1) ResizeArray(&max_option, max_option + 2, (char**)(&top_def_option), "top_def_option", "RunParams::ReadFreeBpot", FILE_NAME_SIZE);
						if (input.length()>FILE_NAME_SIZE-1)
						{
							cout << "\n*****ERROR*****" << endl;
				  			cout << "The size of the " << n_top_def_option + 1 << ". topology define option (" << input << ") in the " << endl;
				  			cout << datfilename << " file exceeds " << FILE_NAME_SIZE << "!" << endl;
				  			cout << "Either use shorter define options or change the value of FILE_NAME_SIZE parameter in units.h and recompile the program!" << endl;
				  			return false;
						}
						std::copy(input.begin(), input.end(), top_def_option[n_top_def_option]);
						top_def_option[n_top_def_option++][input.length()] = 0;
					}
				}
			}
			else if (key == CheckKey("LEAD-BOND"))
			{
				if (lead_series_ind2 != -2)//if -2 multiple definition alraedy, no need to process
				{
					if (lead_series_ind2 > -1)//it was set
					{
						cout << "\nWARNING(" << ++warn << ") : redefinition of leading potential series index in line:" << endl<<"\t"<<*it << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					if (!(ss >> lead_series_ind2)) return ErrLackValue(*it);
					//it has to be offset later, wih the number of potential related interactions with different sigma, not known now
					offset_lead2 += 1;
					}
			}
			else if (key == CheckKey("LEAD-ANGLE"))
			{
				if (lead_series_ind2 != -2)//if -2 multiple definition alraedy, no need to process
				{
					if (lead_series_ind2 > -1)//it was set
					{
						cout << "\nWARNING(" << ++warn << ") : redefinition of leading potential series index in line:" << endl<<"\t"<<*it << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					if (!(ss >> lead_series_ind2)) return ErrLackValue(*it);
					//it has to be offset later, wih the number of potential related interactions with different sigma, not known now
					offset_lead2 += 2;
				}
			}
			else if (key == CheckKey("LEAD-PERDIH"))
			{
				if (lead_series_ind2 != -2)//if -2 multiple definition alraedy, no need to process
				{
					if (lead_series_ind2 > -1)//it was set
					{
						cout << "\nWARNING(" << ++warn << ") : redefinition of leading potential series index in line:" << endl<<"\t"<<*it << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					if (!(ss >> lead_series_ind2)) return ErrLackValue(*it);
					//it has to be offset later, wih the number of potential related interactions with different sigma, not known now
					offset_lead2 += 4;
				}
			}
			else if (key == CheckKey("LEAD-HARMDIH"))
			{
				if (lead_series_ind2 != -2)//if -2 multiple definition alraedy, no need to process
				{
					if (lead_series_ind2 > -1)//it was set
					{
						cout << "\nWARNING(" << ++warn << ") : redefinition of leading potential series index in line:" << endl<<"\t"<<*it << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					if (!(ss >> lead_series_ind2)) return ErrLackValue(*it);
					//it has to be offset later, wih the number of potential related interactions with different sigma, not known now
					offset_lead2 += 8;
				}
			}
			else if (key == CheckKey("LEAD-RBDIH"))
			{
				if (lead_series_ind2 != -2)//if -2 multiple definition alraedy, no need to process
				{
					if (lead_series_ind2 > -1)//it was set
					{
						cout << "\nWARNING(" << ++warn << ") : redefinition of leading potential series index in line:" << endl<<"\t"<<*it << endl;
						lead_series_ind2 = -2;//Indicating, that it was not set yet correctly, default value only set after BPOT is processed as well, as if both present, BPOT can be the lead
					}
					if (!(ss >> lead_series_ind2)) return ErrLackValue(*it);
					//it has to be offset later, wih the number of potential related interactions with different sigma, not known now
					offset_lead2 += 16;
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("BPOT"));//write the appropriate error mesage
				return false;
			}
		}
		if (fnc != 4)
		{
			cout << "\nWARNING(" << ++warn << ") : fnc switch is currently " << fnc << ", for flexible molecules it has to be, fnc=4, it will be set to this value!" << endl;
			fnc = 4;
		}
	}
	return true;
}

#ifdef _AENET
bool RunParams::ReadFreeAenet(std::list<string> &pool)
{
	char  *name = NULL;
	SetArraysize(&name, NAME_SIZE, "name", "RunParams::ReadFreeAenet");
	
	//initialize the defaults
	lead_series_ind2=1;//only one potential contribution exist
	write_E=WRITE_ATOMIC_ENERGY_DEF;
	do_relax=RELAX_DEF;
	aenet_step=AENET_STEP_DEF;
	
	if ((pool.size() > 0))
	{
		
		//Check the required keywords
		if (std::find_if(pool.begin(), pool.end(), [](string c) { return (c.find("SIGMA") != string::npos); }) == pool.end())
		{
			cout << "\nERROR!!: **SIGMA(-MASTER)**/**SIGMA-SCALABLE** are mandatory keywords in AENET section." << endl;
			cout << "Specify one of them in the "<<CheckTag("AENET")<<" section!" << endl;
			return false;
		}
		int i=std::count_if(pool.begin(), pool.end(), [](string c) { return (c.find(CheckKey("IND_ANN-FILE_TYPE-NAME")) != string::npos); });
		if (i!=ntypes) 
		{
			cout << "\nERROR!!: There has to be one " <<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" for each atom type in " << CheckTag("AENET") << " section."  << endl;
			cout<<i<<" "<<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" has been found instead of "<<ntypes<<"."<<endl;
			return false;
		}

		if (std::count_if(pool.begin(), pool.end(), [](string c) { return (c.find(CheckKey("SIGMA")) != string::npos); }) != 1)
			cout << "\nWARNING(" << ++warn << "): Multiple declaration of "<< CheckKey("SIGMA")<<" has been found in " << CheckTag("AENET") << " section." << endl << "\tThe last one will be used!" << endl;
		
		if (chem_symbols==nullptr) 
		{
			cout << "\nERROR!!: ANN potential usage needs "<<CheckKey("CHEMICAL-SYMBOLS")<<" in " << CheckTag("GENERAL") << " section."  << endl;
			cout<<"Correct the "<<datfilename<<" and try again!"<<endl;
			return false;
		}

		for (int i=0;i<ntypes;i++)
		{
			std::copy(chem_symbols_standard[i].begin(), chem_symbols_standard[i].end(), Aenet::aenet_type_names[i]);
			Aenet::aenet_type_names[i][chem_symbols_standard[i].length()]='\0';
		}

		for (std::list<string>::iterator it = pool.begin(); it != pool.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (!(ss >> aenet_weight)) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA")||key == CheckKey("SIGMA-MASTER"))
					aenet_weight = fabs(aenet_weight);
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					aenet_weight = -fabs(aenet_weight);
					ChiSquared::calc_sigma = 1;
				}
			}
			else if (key == CheckKey("IND_ANN-FILE_TYPE-NAME"))
			{
				int index;
				if (!GetInt(ss, *it, &index,"RunParams::ReadFreeAenet")) return ErrLackValue(*it);
				if (index<0 || index>ntypes)
				{
					cout << "\nERROR!!: The index in line "<<endl;
					cout<<*it<<endl;
					cout<<"in the " << CheckTag("AENET") << " section is out of range, has to 1->" <<ntypes<< endl << "\tThe last one will be used!" << endl;
					return false;
				}
				string input;
				ss >> input;
				if (input.length()==0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"No file name was specified for the "<<index<<". atom type after the  "<<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" keyword!"<<endl;
					return false;
				}
				
				if (input.length()>FILE_NAME_SIZE-1)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The length of the ANN file name cannot exceed "<<FILE_NAME_SIZE-1<<" character!"<<endl;
					cout<<input<<" was given for the "<<index<<". atom type after the  "<<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" keyword!"<<endl;
					return false; 
				}
				if (strlen(Aenet::ANN_filenames[index-1])>0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The "<<index<<". ANN potential file name was already given! One file name has to be specified for"<<endl;
					cout<<"each atom type with the "<<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" keyword in the "<<CheckTag("AENET")<<" section!"<<endl;
					return false;
				}
				std::copy(input.begin(), input.end(),Aenet::ANN_filenames[index-1]);
				Aenet::ANN_filenames[index-1][input.length()] = 0;
				input.clear();
				ss >> input;
				if (input.length()>0)
				{
					if (input.length()>2)
					{
						if (chem_symbols_standard!=nullptr)
							cout << "\nWARNING(" << ++warn << "): ";
						else
							cout << "\n*****ERROR*****" << endl;
						cout<<"The length of the ANN atom type name cannot be longer than 2 characters!"<<endl; 
						cout<<input<<" was given for the "<<index<<". atom type after the  "<<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" keyword!"<<endl;
						if (chem_symbols_standard!=nullptr)
						{
							cout<<"Element name "<<chem_symbols_standard[index-1]<<" will be used!"<<endl;
							//it is copied already as a default
						}
						else
						 return false;
					}
					else
					{
						std::copy(input.begin(), input.end(), Aenet::aenet_type_names[index-1]);
						Aenet::aenet_type_names[i][chem_symbols_standard[index-1].length()]='\0';
					}
				}
			}
			else if (key == CheckKey("WRITE-ATOMIC-ENERGY"))
			{
				if (!(ss >> write_E)) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "WRITE-ATOMIC-ENERGY");
				Check0_1(write_E, name, "RunParams::ReadFreeAenet");
			}
			else if (key == CheckKey("CUTOFF-AENET"))
			{
				if (!(ss >> aenet_cutoff)) return ErrLackValue(*it);
				if (aenet_cutoff<0)
				{
					cout << "\nWARNING(" << ++warn << "): AENET cutoff cannot be negativ, default (same as during potential construction) will be used!"<<endl;;
					aenet_cutoff=-1;
				}
			}
			else if (key == CheckKey("AENET-STEP"))
			{
				if (!(ss >> aenet_step)) return ErrLackValue(*it);
				if (aenet_step<0)
				{
					cout << "\nWARNING(" << ++warn << "): Calculation step for ANN potential, AENET-STEP cannot be negative, default ("<<AENET_STEP_DEF <<") is set!" <<endl;;
					aenet_step=AENET_STEP_DEF;
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("AENET"));//write the appropriate error mesage
				return false;
			}
		}
	}
	else
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"The code was compiled with the _AENET compiler option to use ANN potentials," << endl;
		cout << "but no [ AENET ] section was found in the " << datfilename << " file!" << endl;
		cout << "Use the normal code, if do not want to use ANN potential!" << endl;
		cout<<"Exiting..."<<endl;
		CleanExit();
	}
	//check if there is a potential file name for each type
	for (int i=0;i<ntypes;i++)
	{
		if (strlen(Aenet::ANN_filenames[i])==0)
		{
			cout << "\n*****ERROR*****" << endl;
			cout<<"The "<<i+1<<". ANN potential file name is missing! A file name has to be specified for"<<endl;
			cout<<"each atom type with the "<<CheckKey("IND_ANN-FILE_TYPE-NAME")<<" keyword in the "<<CheckTag("AENET")<<" section!"<<endl;
			return false;
		}
	}
	if (name != NULL)
		delete[] name;
	potential=20;//ANN is used
	NB_weight_mode=0;//only one value
	mystrcpy(vdW_name, 15, "ANN");
	return true;
};
#endif

int RunParams::RegularizeSigmaEntry(string &parameters)
//This function regularizes the sigma entry as writes "-" signs instead of "s" in parameters and 
//returns with the master index (start from 1) if it has been found,
//0 if not
//-1 if more than one master tags have been found
{
	int retval = 0;
	//Looking for the master tag(s)
	if (parameters.find_first_of("Mm") != string::npos)
	{
		std::size_t pos = parameters.find_first_of("Mm");
		if (parameters.find_first_of("Mm", pos + 1) != string::npos) return -1;//More than one master tags have been found
		//Determine the no and erase the item from the original line
		stringstream ss(parameters);
		string str;
		retval = 1;
		while ((ss >> str) && (str.find_first_of("Mm") == string::npos)) retval++;
		parameters.erase(parameters.find_first_of("Mm"), 1);
	}
	std::size_t position = 0;
	while ((position = parameters.find_first_of("Ss", position)) != string::npos) parameters.replace(position, 1, "-");
	return retval;
}


bool RunParams::HaveReference(std::list<string> *pointer, std::list<std::list<string> > &target)
//Returns with true, if pointer points to an element of the target
{
	if (target.size() == 0) return false;
	bool retval = false;
	for (std::list<std::list<string> >::iterator it = target.begin(); (it != target.end()) && (!retval); it++) retval = (pointer == &(*it));
	return retval;
}

void RunParams::CreateFreeFormat(const int write_reduced_free)
//This method creates .free free format file with all the parameters and reduced .freer and .freers free format files
//the difference between .freer and .freers is that if scalable sigma applied or automatic cutoff dtermination is used, then the freer contains the original scaling
//fractions for sigma and the MIN for CUTOFF, while the .freers the actual scaled sigma values and cutoffs
{
	char *tempfilename=NULL;
	ofstream outdef;
	
	switch (write_reduced_free)
	{
	case 0: tempfilename = freefilename; break;
	case 1: tempfilename = freerfilename; break;
	case 2: tempfilename = freersfilename;
	}
	OpenFile(outdef, tempfilename, "RunParams::CreateFreeFormat",0);


	CreateFreeGeneral(outdef, write_reduced_free);
	for (int i = 0; i < ngr + nsq + nfq + nfg+nek; i++) CreateFreeExp(i, outdef, write_reduced_free);
	for (int i = 0; i < ncosdistr; i++) CreateFreeCos(i, outdef, write_reduced_free);
	for (int i = 0; i < nicoord; i++) CreateFreeCoord(i, outdef, false, write_reduced_free);
	for (int i = 0; i < navcoord; i++) CreateFreeCoord(i, outdef, true, write_reduced_free);
#ifdef _ADVANCED_GEOM_CONST
	for (int i = 0; i < ncommonneigh; i++) CreateFreeConc(i, outdef, write_reduced_free);
	for (int i = 0; i < nsecondneigh; i++) CreateFreeSnc(i, outdef, write_reduced_free);
	for (int i = 0; i < nbvs; i++) CreateFreeBvs(i, outdef, write_reduced_free);
#endif
	if (custmove > 0) CreateFreeCustom(outdef);
	if (swap_fraction > 0.0) CreateFreeSwap(outdef);

	if (potential ==1 || potential == 10) CreateFreeNbpot(outdef, write_reduced_free);
#ifdef _AENET
	if (potential ==20) CreateFreeAenet(outdef,write_reduced_free);
#endif
	if (fnc == 4) CreateFreeBpot(outdef, write_reduced_free);
#ifdef _LOCAL_INV
	if (nlocint > 0) CreateFreeLoc(outdef, write_reduced_free);
#endif
#ifdef _VIBR_AMP
	CreateFreeVibramp(outdef);
#endif
#ifdef _NO_PERIODIC
	CreateFreeNoper(outdef, write_reduced_free);
#endif
	outdef << endl << "[ "<<CheckTag("END")<<" ]" << endl;
	outdef.close();
}

void RunParams::CreateFreeGeneral(std::ofstream &out, const bool onlyndef)
//This method writes general tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	out << "#003" << endl << "[ "<<CheckTag("GENERAL")<<" ]" << endl;
	string titl(RunParams::title);
	titl.erase(std::remove(titl.begin(), titl.end(), '\n'), titl.end());
	titl.erase(std::remove(titl.begin(), titl.end(), '\r'), titl.end());
	string titl2 = titl;
	titl2.erase(std::remove(titl2.begin(), titl2.end(), ' '), titl2.end());
	if ((!onlyndef) || (onlyndef && (titl2 != ""))) out << " " << CheckKey("TITLE")<<" = " << titl << "\t! --Title of the run (default: blank)" << endl;
	std::streamsize ss = out.precision();
	out.precision(15);
	out.setf(ios::scientific, ios::floatfield);
	if (!onlyndef || write_density) out << " " << CheckKey("NDENS") << " = " << rho << "\t\t! --Number density A^{-3} (A=angstrom) units" << endl << "\t\t\t! alternative: HALFBOX = " << SimpleCfg::boxedge << " ! Half box length in A units;" << endl << "\t\t\t! (default: determine it from the .cfg file)" << endl;
	out.unsetf(ios::floatfield);
	out.setf(ios::fixed, ios::floatfield);
	out.precision(ss);
	if ((!onlyndef) || (onlyndef && ((abs(auto_cutoff)&1) != (AUTO_CUTOFF_DEF&1))))//1, -1 and 3, -3 all considered def
	{
		out << " " << CheckKey("CUT-OFF") << " =";
		if (abs(auto_cutoff)==1)
			out << " MIN-ALL";
		else if (abs(auto_cutoff) == 3)
		{
			out << " MIN";
		}
		else
		{
			double *p = pcutoff;
			if (ntypes > 1) out << endl;
			for (int i = 0; i < ntypes; i++)
			{
				out << " ";
				for (int j = 0; j < ntypes - i; j++, p++) out << " " << *p;
				if (i != ntypes - 1) out << endl;
			}
		}
		out << "\t\t! --Distance of closest approach for each partials in A (in 11,12,...,22,23.. order);" << endl << "\t\t\t! Give as many values as the number of partials, or MIN for automatic cutoff determination" << endl;
	}
	bool found = !onlyndef;
	for (int i = 0; (i < ntypes) && (!found); i++) if ((!found) && (pmaxmove[i] != 0.1)) found = true;
	if (found)
	{
		out << " " << CheckKey("MAX-MOVES") << " =";
		int last = 0;
		if ((onlyndef) && (ntypes > 1))
		{
			for (int i = ntypes - 2; (i >= 0) && (pmaxmove[i] == pmaxmove[i + 1]); i--) last++;
		}
		for (int i = 0; i < ntypes - last; i++) out << " " << pmaxmove[i];
		out << "\t! --Maximal moves for each atomtypes in A units; (default: "<<MAX_MOVES_DEF<<" for each);" << endl << "\t\t\t! when you provide less values than required, the last value passed to the remaining ones" << endl;
	}
	//*.free			*.freer or *.freers
	if ((!onlyndef) || (onlyndef && (rspacing_def != RSPACING_DEF))) out << " " << CheckKey("R-SPACING") << " = " << rspacing_def << "\t! --R-spacing in A units (default: "<<RSPACING_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (moveout !=MOVEOUT_DEF))) out << " " << CheckKey("MOVEOUT") << " = " << moveout << "\t\t! --Whether to use moveout option (default: "<<MOVEOUT_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (too_close_fraction != TOOCLOSE_FRACTION_DEF))) out << " " << CheckKey("TOOCLOSE-FRACTION") << " = " << too_close_fraction << "\t\t! --The fraction of moves to choose from the too close atoms in case of moveout option (default: " << TOOCLOSE_FRACTION_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (printstep != PRINTSTEP_DEF))) out << " " << CheckKey("PRINT-STEP") << " = " << printstep << "\t! --Generated number of steps for printing (default: "<<PRINTSTEP_DEF<<")" << endl;
#ifndef _TEST_MODE	
	if (runmode == 0)
	{
		if ((!onlyndef) || (onlyndef && runlimit != TIMELIM_DEF))
		{
			out << " " << CheckKey("RUN-TIME") << " = " << runlimit << "\t! --Time limit of the simulation in minutes (default: " << TIMELIM_DEF << ")" << endl;
			out << "\t\t\t! alternatives: " << CheckKey("LAST-GEN") << " = " << LAST_MOVE_DEF << ", or " << CheckKey("LAST-ACC") << " = " << LAST_MOVE_DEF << ", (default " << CheckKey("RUN-TIME") << " = " << TIMELIM_DEF << ")" << endl;
		}
	}
	if (runmode == 1)
	{
		out << " " << CheckKey("LAST-GEN") << " = " << (int)-runlimit << "\t! --Number of moves to generate during the simulation " << endl;
		out << "\t\t\t! alternatives: " << CheckKey("RUN-TIME") << " = " << TIMELIM_DEF << ", or " << CheckKey("LAST-ACC") << " = " << LAST_MOVE_DEF << ", (default " << CheckKey("RUN-TIME") << " = " << TIMELIM_DEF << ")" << endl;
	}
	if (runmode == 2)
	{
		out << " " << CheckKey("LAST-ACC") << " = " << (int)runlimit << "\t! --Number of accepted moves during the simulation " << endl;
		out << "\t\t\t! alternatives: " << CheckKey("RUN-TIME") << " = " << TIMELIM_DEF << ", or " << CheckKey("LAST-GEN") << " = " << LAST_MOVE_DEF << ", (default " << CheckKey("RUN-TIME") << " = " << TIMELIM_DEF << ")" << endl;
	}
#else
	if (runmode == 1)
	{
		if ((!onlyndef) || (onlyndef  && -runlimit != LAST_MOVE_DEF))
		{
			out << " " << CheckKey("LAST-GEN") << " = " << (int)-runlimit << "\t! --Number of moves to generate during the simulation " << endl;
			out << "\t\t\t! alternatives in _TEST_MODE: " << CheckKey("LAST-ACC") << " = " << LAST_MOVE_DEF << ", (default " << CheckKey("LAST-GEN") << " = " << LAST_MOVE_DEF << ")" << endl;
		}
	}
	else if (runmode == 2)
	{
		out << " " << CheckKey("LAST-ACC") << " = " << (int)runlimit << "\t! --Number of accepted moves during the simulation " << endl;
		out << "\t\t\t! alternatives in _TEST_MODE: " << CheckKey("LAST-GEN") << " = " << LAST_MOVE_DEF << ", (default " << CheckKey("LAST-GEN") << " = " << LAST_MOVE_DEF << ")" << endl;
	}
#endif
	if ((!onlyndef) || (onlyndef && (timesave != TIMESAVE_DEF))) out << " " << CheckKey("SAVE-TIME") << " = " << timesave << "\t\t! --Time step for saving in minutes (default: "<<TIMESAVE_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (histbuffsize > 0)))
	{
		if ((!onlyndef) || (onlyndef && (histbuffsize != HIST_BUFF_DEF))) out << " " << CheckKey("HST-BUFFSIZE") << " = " << histbuffsize << "\t\t! --Size of the history buffer in lines (default: " << HIST_BUFF_DEF << ", no history)" << endl;
		if ((!onlyndef) || (onlyndef && (histbuffsize != HIST_BUFF_DEF && histstepratio != HIST_STEP_DEF)))
		{
			out << " " << CheckKey("HST-STEP-FACTOR") << " = " << histstepratio << "\t! --Number of savings between each history buffering for the .hst file (default: " << HIST_STEP_DEF << " -- means history saved at every " << HIST_STEP_DEF << ". saving)." << endl;
			out << "\t\t\t! alternative: " << CheckKey("HST-SAVE-TIME")<<" = " << (histstepratio + 1)*timesave << " ! --History save time to *.hst file in minutes" << endl;
			out << "\t\t\t! this value is adjusted internally to be an integer multiple of SAVE-TIME." << endl;
		}
	}
	if ((!onlyndef) || (onlyndef && (cfgnumb != CFG_COLL_DEF))) out << " " << CheckKey("COLL-NUMBER") << " = " << cfgnumb << "\t! --Number of configuration to collect after RUN-TIME has been elapsed (default: "<<CFG_COLL_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (cfgnumb != CFG_COLL_DEF) && (coll_frequency != CFG_COLL_FREQ_DEF)))
	{
		out << " " << CheckKey("COLL-STEP-FACTOR") << " = " << coll_frequency << "\t ! --Collection step, configuration is collected at every COLL-STEP_FACTOR. saving (default: "<<CFG_COLL_FREQ_DEF<<")"<< endl;
		out << "\t\t\t! alternative: " << CheckKey("COLL-SAVE-TIME")<<" = " << coll_frequency * timesave << " ! --Time step for saving collected configuration in minutes (default: equals to SAVE-TIME)," << endl;
		out << "\t\t\t! this value is adjusted internally to be an integer multiple of SAVE-TIME." << endl;
		
	}
	
	if ((!onlyndef) || (onlyndef && (fnc != FNC_DEF) && (fnc != 4)))
	{
		out << " " << CheckKey("FNC-TYPE") << " = ";
		switch (fnc)
		{
		case 1: out << "NORMAL";
			break;
		case 2: out << "ADJUST";
			break;
		case 3: out << "MOVE_IN";
			break;
		default: out << "NONE";
		}
		out << "\t! --FNC switch: NONE(default) or NORMAL or ADJUST or MOVE_IN" << endl;
	}
	if (chem_symbols != nullptr)
	{
		out << " " << CheckKey("CHEMICAL-SYMBOLS") << " = ";
		if (chem_symbols != nullptr)
			for (int i= 0; i < ntypes; i++)
				out << chem_symbols[i] << " ";
		out << "\t\t! --Chemical symbols in the order of the atom types in the configuration (default: not given)" << endl;
	}
	if ((!onlyndef) || (onlyndef && (binshift_def != BINSHIFT_DEF))) out << " " << CheckKey("BIN-SHIFT") << " = " << binshift_def << "\t\t! --Initial bin shift (default: "<<BINSHIFT_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (xmax_ori != XMAX_DEF))) out << " " << CheckKey("XMAX-FACTOR") << " = " << xmax_ori << "\t! --Xmax= FACTOR*halbox used in the run (default: "<<XMAX_DEF<<")" << endl;
	if (strcmp(sffilename,"")!=0) out << " " << CheckKey("CUSTOM-SFACTORTABLE") << " = " << sffilename << "\t! --Custom surface factor table is used (default: built-in is used)" << endl;
	if ((!onlyndef) || (onlyndef && (reload != RELOAD_DEF))) out << " " << CheckKey("RELOAD") << " = " << reload << "\t\t! --Whether to load the histogram from file if possible (default: "<<RELOAD_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (max_gridatom != MAX_GRIDATOM_DEF))) out << " " << CheckKey("ATOMS-IN-GRIDCELL") << " = " << max_gridatom << "\t! --Maximum number of atoms in a gridcell (default: "<< MAX_GRIDATOM_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (exafs_shiftstep != E0_SHIFTSTEP_DEF))) out << " " << CheckKey("EXAFS-SHIFTSTEP") << " = " << exafs_shiftstep << "\t! --Shift E0 for selected EXAFS data sets in each " << CheckKey("EXAFS-SHIFTSTEP") << ". simulation step (default: " << E0_SHIFTSTEP_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (lambda_nonlin != NONLIN_LAMBDA_DEF))) out << " " << CheckKey("LAMBDA-NONLIN") << " = " << lambda_nonlin << "\t! --Lambda fudge factor for Levenberg-Marquardt non-linear regression (default: " << NONLIN_LAMBDA_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (factor_nonlin != NONLIN_FACTOR_DEF))) out << " " << CheckKey("FACTOR-NONLIN") << " = " << factor_nonlin << "\t! --Factor to change fudge for Levenberg-Marquardt non-linear regression (default: " << NONLIN_FACTOR_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (terminate_nonlin != NONLIN_TERMINATE_DEF))) out << " " << CheckKey("TERMINATE-NONLIN") << " = " << terminate_nonlin << "\t! --Terminate after this many successive moves ended with failed non-linear regression (default: " << NONLIN_TERMINATE_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (niter_nonlin != NONLIN_MAX_NIT_DEF))) out << " " << CheckKey("MAX-NITER-NONLIN") << " = " << niter_nonlin << "\t! --Maximum number of iteration in case of I(Q) fitting with non-linear regression (default: " << NONLIN_MAX_NIT_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (epsilon_nonlin != NONLIN_EPSILON_DEF))) out << " " << CheckKey("EPSILON-NONLIN") << " = " << epsilon_nonlin << "\t! --Iteration stops, if the fitted parameters change less than this (default: " << NONLIN_EPSILON_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (IQ_backg_corr_step != IQ_BACKG_CORR_STEP_DEF))) out << " " << CheckKey("I(Q)BACKG-STEP") << " = " << IQ_backg_corr_step << "\t\t! --Change mu background correction factor for I(Q) fitting in each " << CheckKey("I(Q)BACKG-STEP") << ". simulation step, (default: " << IQ_BACKG_CORR_STEP_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (AXS_shiftstep != AXS_SHIFTSTEP_DEF))) out << " " << CheckKey("AXS-SHIFTSTEP") << " = " << AXS_shiftstep << "\t\t! --Change the f' for AXS I(Q) fitting in each " << CheckKey("AXS-SHIFTSTEP") << ". simulation step, (default: " << AXS_SHIFTSTEP_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (debug != DEBUG_DEF))) out << " " << CheckKey("DEBUG") << " = " << debug << "\t\t! --Write debug information (default: " << DEBUG_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (write_log != WRITE_LOG_DEF))) out << " " << CheckKey("WRITE-LOG") << " = " << write_log << "\t\t! --Whether to create a log file (default: " << WRITE_LOG_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (log_nlr_steps != LOG_NONLIN_STEPS_DEF))) out << " " << CheckKey("LOG-NONLIN-STEPS") << " = " << log_nlr_steps << "\t! --Write non-linear regression steps, if there is any (default: " << LOG_NONLIN_STEPS_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (write_exafs_coeffs != WRITE_EXAFS_COEFFS_DEF))) out << " " << CheckKey("WRITE-EXAFS-COEFFS") << " = " << write_exafs_coeffs << "\t\t! --Whether to write the EXAFS coefficients into *.expt (only active if dE0 shift is applied (default: " << WRITE_EXAFS_COEFFS_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (r_switch_power != R_SWITCH_POWER_DEF))) out << " " << CheckKey("R-SWITCH-POWER") << " = " << r_switch_power << "\t\t! --Modify the default value for the _R_SWITCH_GR_MIN_1 compiler option, if it is on, then\n\t\t\t! experimental g(r) data sets will be interpreted as (r^r_switch_power)*[g(r)-1], and this will be calculated and fitted, default: " << R_SWITCH_POWER_DEF << ")" << endl;
	if ((!onlyndef) || (onlyndef && (nthreads != NTHREADS_DEF))) out << " " << CheckKey("THREADS") << " = " << nthreads << "\t\t! --Total number of threads to use (default: "<<NTHREADS_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (old_out != OLD_OUT_DEF))) out << " " << CheckKey("CREATE-OUT") << " = " << old_out << "\t\t! --Create old .out output format containing ppcf's, psq's and agreement to experimental data in single file (default: "<<OLD_OUT_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (sum_ppcf != SUM_PPCF_DEF))) out << " " << CheckKey("PPCF-AVERAGE") << " = " << sum_ppcf << "\t! --Append the average of ppcf's of each saved configuration to the .ppcf file (default: "<<SUM_PPCF_DEF<<")" << endl;

	return;
}

void RunParams::CreateFreeExp(const int expno, std::ofstream &out, const bool onlyndef)
//This method writes the expno'th experiment related tags into free control files defined the stream "out"
//If onlyndef is true, it writes only the non-default values to out
{
	int i;
	int from, to;
	int datano,size;
	int myindex;//index of the data set in its own type
	double datafromq, datatoq, ctrlfromq, ctrltoq;
	string line;
	enum Exptype { GR, ND, XRD, EDIFF, EXAFS };
	Exptype et;

	out << endl << "[ "<< CheckTag("EXP")<<" ]" << endl;
	out << " TYPE = ";
	
	
	if (ngr > expno)
	{
		et = GR;
		myindex = expno;
	}
	else if (ngr + nsq > expno)
	{
		et = ND;
		myindex = expno - ngr;
	}
	else if (ngr + nsq + nfq > expno)
	{
		et = XRD;
		myindex = expno - ngr - nsq;
	}
	else if (ngr + nsq + nfq + nfg > expno)
	{
		et = EDIFF;
		myindex = expno - ngr - nsq - nfq;
	}
	else
	{
		et = EXAFS;
		myindex = expno - ngr - nsq - nfq - nfg;
	}
	switch (et)
	{
	case GR: out << CheckKeyValue("TYPE","GR")<<" "; break;
	case ND: out << CheckKeyValue("TYPE", "ND")<<" "; break;
	case XRD: out << CheckKeyValue("TYPE", "XRD") << " "; break;
	case EDIFF: out << CheckKeyValue("TYPE", "EDIFF") << " "; break;
	case EXAFS:
	default:out << CheckKeyValue("TYPE", "EXAFS") << " ";
	}
	out << "\t\t! --Type of the experiment, **MANDATORY ITEM** (it should be "<<GiveKeyValues("TYPE")<< endl;
	out << " "<< CheckKey("DATAFILE")<<" = " << ExptsData::datafilename[expno + ((et == EXAFS) ? (myindex) : 0)] << "\t! --File name containing the experimental data, **MANDATORY ITEM**" << endl;
	
	switch (et)
	{
	case GR: from = ExptsData::grmin[expno]; to = ExptsData::grmax[expno]; size = ExptsData::grsize[expno]; break;
	case ND: from = ExptsData::sqmin[myindex]; to = ExptsData::sqmax[myindex]; size = ExptsData::sqsize[myindex]; break;
	case XRD: from = ExptsData::fqmin[myindex]; to = ExptsData::fqmax[myindex]; size = ExptsData::fqsize[myindex]; break;
	case EDIFF: from = ExptsData::fgmin[myindex]; to = ExptsData::fgmax[myindex]; size = ExptsData::fgsize[myindex]; break;
	default: from = ExptsData::ekmin_ori[myindex]; to = ExptsData::ekmax_ori[myindex]; size = ExptsData::eksize[myindex];//to write the original range in case of E0 shift
	}
	ifstream infile(ExptsData::datafilename[((et == EXAFS) ? (ngr+nsq+nfq+nfg+2*(myindex)) : expno)]);
	
	infile >> datano;//might not be correct, the files were already processed, the corrected values had to be used

	getline(infile, line);
	getline(infile, line);
	
	for (i = 0; i < size; i++)
	{
		double data;
		infile >> data;
		getline(infile, line);
		if (i == 0)
			datafromq = data;
		else
			if (i == size - 1)
				datatoq = data;
		if (i == from - 1)
			ctrlfromq = data;
		if (i == to - 1)
			ctrltoq = data;
	}
	infile.close();
	if ((!onlyndef) || (onlyndef && ((from != 1) || (to != size))))
	{
		string text = "Range";
		string text2 = "applied";
		switch (et)
		{
		case GR: out << " "<< CheckKey("R-RANGE")<<" = "; break;
		case ND:
		case XRD: out << " " << CheckKey("Q-RANGE") << " = "; break;
		case EDIFF: out << " " << CheckKey("G-RANGE") << " = "; break;
		default: 
			out << " " << CheckKey("K-RANGE") << " = "; 
			if (ExptsData::ek_ngrid_in[myindex] > 0)
			{
				text.assign("Original range");
				text2.assign("original");
			}
		}
		out << ctrlfromq << " " << ctrltoq << "\t! --"<<text<<" used for fit (default: " << datafromq << ".." << datatoq << ", the whole range of the data)" << endl;
		out << "\t\t\t! alternative: "<< CheckKey("POINT-RANGE")<<" = " << from << " " << to << " ! The "<< text2 << " point range to fit (default: 1.." << size << ", the whole range of the data)" << endl;
		if (et==EXAFS && ExptsData::ek_ngrid_in[myindex]>0)
			out << "\t\t\t! actual range used due to E0 shit interpolation data point loss: " << ExptsData::ekmin[myindex] << " " << ExptsData::ekmax[myindex] << endl;

	}
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
	if ((!onlyndef) || (onlyndef && (rspacing[expno] != RSPACING_DEF))) out << " " << CheckKey("R-SPACING") << " = " << rspacing[expno] << "\t! --R-spacing in A units (default: " << RSPACING_DEF << ")" << endl;
#endif
#endif
	double myvalue;
	switch (et)
	{
	case GR: myvalue = ExptsData::grsub[expno]; break;
	case ND: myvalue = ExptsData::sqsub[myindex]; break;
	case XRD: myvalue = ExptsData::fqsub[myindex]; break;
	case EDIFF: myvalue = ExptsData::fgsub[myindex]; break;
	default: myvalue = 0;
	}
	if (((!onlyndef) && (et != EXAFS)) || (onlyndef && (myvalue != RENORM_DEF))) out << " " << CheckKey("CONST-SUBTRACT") << " = " << myvalue << "\t! --Constant to subtract (default: " << RENORM_DEF << ")" << endl;
	double *par = 0;
	switch (et)
	{
	case GR: par = ExptsData::gr_coeffs + expno * (ntypes*(ntypes + 1)) / 2; break;
	case ND: if (ExptsData::sqreadcoeffs[myindex] || ExptsData::is_Ncoeff_calc==2) par = ExptsData::sq_coeffs + (myindex)*(ntypes*(ntypes + 1)) / 2; break;
	//for neutron the coeffs are saved in the *.freers file, after they were calculated
	default: par = 0;
	}
	if (par != 0)
	{
		out.precision(16);
		out << " " << CheckKey("PARTIAL-COEFFS") << " = ";
		if (ntypes > 1) out << endl;
		for (int i = 0; i < ntypes; i++)
		{
			out << " ";
			for (int j = 0; j < ntypes - i; j++, par++) out << " " << *par;
			if (i != ntypes - 1) out << endl;
		}
		out << "\t\t! --Weight of each partials (in 11,12,...,22,23.. order);" << endl << "\t\t\t! **MANDATORY ITEM**: you should provide as many values as each partial of your system have" << endl;
		out.precision(8);
	}
	if (et == ND)
	{
		for (int j = 0; j < ntypes; j++)
		{
			if (ExptsData::isotope_count[myindex * ntypes + j].count > 0 && ExptsData::is_Ncoeff_calc == 1)//only for *.free and freer
			{
				out << " " << CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS") << " = " << ExptsData::isotope_count[myindex * ntypes + j].count;
				for (int k = 0; k < ExptsData::isotope_count[myindex * ntypes + j].count; k++)
					out << " " << ExptsData::isotope_count[myindex * ntypes + j].isotopes[k].symbol;
				for (int k = 0; k < ExptsData::isotope_count[myindex * ntypes + j].count; k++)
					out << " " << ExptsData::isotope_count[myindex * ntypes + j].isotopes[k].ratio;
				out << "\t\t! --Number of isotopes for an element, symbols for each isotope in Cl35 format, ratios, sum of ratios should be 1.0!"<<endl;
			}
		}
	}
	switch (et)
	{
	case GR: myvalue = ExptsData::grsigma[expno]; break;
	case ND: myvalue = ExptsData::sqsigma[myindex]; break;
	case XRD: myvalue = ExptsData::fqsigma[myindex]; break;
	case EDIFF: myvalue = ExptsData::fgsigma[myindex]; break;
	default: myvalue = ExptsData::eksigma[myindex];
	}
	std::streamsize ss = out.precision();
	if (fabs(myvalue) < pow(10, -ss))//would not fit into the frame
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::scientific, ios::floatfield);
	}
	i = 0;
	if (expno == lead_series_ind - 1)
	{
		out << " " << CheckKey("SIGMA-MASTER") << " = " << myvalue;
		i = 1;
	}
	else
	{
		if (myvalue < 0.0)
		{
			out << " " << CheckKey("SIGMA-SCALABLE") << " = " << -myvalue;
			i = 2;
		}
		else
			out << " " << CheckKey("SIGMA") << " = " << myvalue;
	}
	if (fabs(myvalue) < pow(10, -ss))//set it back to fixed
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::fixed, ios::floatfield);
	}
	//The sigma description texts are stored beginning with lower case, as they are used for the combined sigmas as well, transform it here
	string mytext= sigma_text[i];
	std::transform(mytext.begin(), ++mytext.begin(), mytext.begin(), ::toupper);
	out<< "\t! --" <<mytext << endl << alternative_sigma_text[i] << endl;

	switch (et)
	{
		case GR: myvalue = ExptsData::gruseR[expno]; break;
		case ND: myvalue = ExptsData::squseR[myindex]; break;
		case XRD: myvalue = ExptsData::fquseR[myindex]; break;
		case EDIFF: myvalue = ExptsData::fguseR[myindex]; break;
		default: myvalue = ExptsData::ekuseR[myindex];
	}
	if ((!onlyndef) || (onlyndef && (myvalue != USE_RFACTOR_DEF))) out << " " << CheckKey("USE-RFACTOR") << " = " << myvalue << "\t! --Whether to use R-factor (1) instead of normal chi2 for this set (default: " << USE_RFACTOR_DEF << ")" << endl;
	
	if (et == XRD)
	{
		if ((!onlyndef) || (onlyndef && (ExptsData::fqreadcoeffs[myindex] != READ_COEFFS_DEF)))//whether to read X-ray coeffs writing, if free or reduced,but not def
		{
			out << " " << CheckKey("READ-COEFFS") << " = " << ExptsData::fqreadcoeffs[myindex];
			out << "\t! --Whether to read the X-ray coefficients based on the Waasmaier-Kirfel table (default: " << READ_COEFFS_DEF << ")" << endl;
		}
		if ((!onlyndef) || (onlyndef && (ExptsData::fqfitIQ[myindex] != IQFIT_DEF)))//I(Q) fit writing, if free or reduced,but not def
		{
			out << " " << CheckKey("I(Q)_A_B_ALPHA") << " = " << ExptsData::fqfitIQ[myindex] << " ";
			if (ExptsData::fqfitIQ[myindex])//there is I(Q) fit, write the parameters regardless whether they are default or not, default will be given for the fitted
			{
				if (ExptsData::fqrenorm[myindex])
					out << 1 << " ";
				else
					out << ExptsData::fqa[myindex] << " ";

				if (ExptsData::fqoffset[myindex])
					out << 0 << " ";
				else
					out << ExptsData::fqb[myindex] << " ";

				if (ExptsData::fqrenalpha[myindex])
					out << 0;
				else
				{
					if (ExptsData::fqoffset[myindex])//fitindex 3:a, b fitted, alpha=0
						out << 0;
					else// fitindex 1: a fitted, b and alpha custom
						out << ExptsData::fqalpha[myindex];
				}
			}
			out << "\t! --Whether to fit I(Q) instead of structure factor (default: " << IQFIT_DEF << ")" << endl;
			out << "\t\t\t! If custom values for a, b or alpha are needed, should be supplied here. If only b and alpha are needed, give something for a as well!" << endl;
		}
		if ((!onlyndef) || (onlyndef && (ExptsData::fqusecompton[myindex] != COMPTON_DEF)))//using Compton contributionfor I(Q), if free or reduced,but not def
		{
			out << " " << CheckKey("COMPTON") << " = " << ExptsData::fqusecompton[myindex] << " ";
			out << "\t! --Whether to use Compton contribution during I(Q) calculation (default: " << COMPTON_DEF << ")" << endl;
		}
		if ((!onlyndef) || (onlyndef && (ExptsData::fqIQbackgcorr[myindex] != IQBACKGCORR_DEF)))
		{
			out << " " << CheckKey("I(Q)BACKG_MU_DMU-MAX") << " = " << ExptsData::fqIQbackgcorr[myindex];
			if (ExptsData::fqIQbackgcorr[myindex])
				out << " " << ExptsData::fqmu[myindex] << " " << ExptsData::fqdmumax[myindex];

			out << "\t! --Whether to use background correction in case of I(Q) data fit (default: " << IQBACKGCORR_DEF << ")" << endl;
			out << "\t\t\t! If custom values mu and dmu_max are needed, should be supplied here. If only dmu_max is needed, mu has to be given as well!" << endl;
			out << "\t\t\t! Default values are for mu: " << IQ_MU_DEF << " dmu_max: " << IQ_DMU_MAX_DEF << endl;
		}
		

		if ((!onlyndef) || (onlyndef && (ExptsData::fqAXS[myindex] > 0)))//there is AXS
		{
			out << " " << CheckKey("AXS-INDEX_FPRIME_FRACTION") << " = " << ExptsData::fqAXS[myindex] << " ";
			if (ExptsData::fqAXS[myindex] > 0)
				out << ExptsData::fqfprimeact[myindex * ntypes + ExptsData::fqAXS[myindex] - 1] << " " << ExptsData::fqfprimefactor[myindex];
			out << "\t! --Index of the atom type to use to use AXS for, f' and fraction\n\t\t\t! to vary f' between f'-f'*fraction->f'+f'*fraction (default: 0, not to use it)" << endl;
		}
		if ((!onlyndef) || (onlyndef && (ExptsData::fqfprimeindex[myindex] - (int)pow(2, ExptsData::fqAXS[myindex] - 1) > 0)))//there is  f' usage beside AXS
		{
			if ((ExptsData::fqAXS[myindex] > 0 && ExptsData::fqfprimecount[myindex] > 1) || (ExptsData::fqAXS[myindex] == 0 && ExptsData::fqfprimecount[myindex] > 0))//there are fixed f' using sets
			{
				for (int itype = 0; itype < ntypes; itype++)
				{
					if ((itype != ExptsData::fqAXS[myindex] - 1) && ExptsData::fqfprimeindex[myindex] & (int)pow(2, itype))//this type has fixed f' on it is not AXS
					{
						out << " " << CheckKey("FPRIME-INDEX_FPRIME") << " = " << itype + 1 << " ";
						out << ExptsData::fqfprime[myindex * ntypes + itype];
						out << "\t! --Index of the atom type to use f' for and f' (default: 0, not to use it)" << endl;
					}
				}

			}
			else
			{
				out << " " << CheckKey("FPRIME-INDEX_FPRIME") << " = " << 0 << " ";
				out << "\t! --Index of the atom type to use f' for and f' (default: 0, not to use it)" << endl;
			}

		}
		
	}

	switch (et)
	{
		case GR: to = ExptsData::grrenorm[expno]; break;
		case ND: to = ExptsData::sqrenorm[myindex]; break;
		case XRD: to = ExptsData::fqrenorm[myindex]; break;
		case EDIFF: to = ExptsData::fgrenorm[myindex]; break;
		default: to = ExptsData::ekrenorm[myindex];
	}
	if ((!onlyndef) || (onlyndef && (to != 0))) out << " " << CheckKey("RENORM") << " = " << to << "\t\t! --Whether to vary amplitudes/renormalize (default: "<<RENORM_DEF<<")" << endl;
	int var[4];
	switch (et)
	{
	case GR: var[0] = ExptsData::groffset[expno]; var[1] = ExptsData::grlinear[expno]; var[2] = ExptsData::grquadratic[expno]; var[3] = ExptsData::grcubic[expno]; break;
	case ND: var[0] = ExptsData::sqoffset[myindex]; var[1] = ExptsData::sqlinear[myindex]; var[2] = ExptsData::sqquadratic[myindex]; var[3] = ExptsData::sqcubic[myindex]; break;
	case XRD: var[0] = ExptsData::fqoffset[myindex]; var[1] = ExptsData::fqlinear[myindex]; var[2] = ExptsData::fqquadratic[myindex]; var[3] = ExptsData::fqcubic[myindex]; break;	
	case EDIFF: var[0] = ExptsData::fgoffset[myindex]; var[1] = ExptsData::fglinear[myindex]; var[2] = ExptsData::fgquadratic[myindex]; var[3] = ExptsData::fgcubic[myindex]; break;
	default: var[0] = ExptsData::ekoffset[myindex]; var[1] = var[2] = var[3] = 0;
	}
	if ((!onlyndef) || (onlyndef && (var[0] || var[1] || var[2] || var[3])))
	{
		int last = 3;
		for (int i = last; (i >= 0) && (!var[last]); i--, last--);
		out << " " << CheckKey("POLY-BACK-FLAGS") << " =";
		for (int i = 0; i <= (onlyndef ? last : 3); i++) out << " " << var[i];
		out << "\t! --Whether to vary polynomial terms: 0th,1st,.. order (default: " << RENORM_DEF << " for each)" << endl;
	}
	if (et==XRD)
		if ((!onlyndef) || (onlyndef && (ExptsData::fqrenalpha[myindex] != RENORM_DEF))) out << " " << CheckKey("COMPTON-SUPR-COEFF") << " = " << ExptsData::fqrenalpha[myindex] << "\t! --Whether to fit alpha parameter in I(Q) (default: " << RENORM_DEF << ")" << endl;
	if (et == EXAFS)
	{
		out << " " << CheckKey("ABSORBER-TYPE") << " = " << ExptsData::ek_abstype[myindex] + 1 << "\t! --Type of absorbing atom, **MANDATORY ITEM**" << endl;
		if ((!onlyndef) || (onlyndef && (ExptsData::ekchipower[myindex] != EXAFS_CHI2_POWER_DEF))) out << " " << CheckKey("CHIK-POWER")<<" = " << ExptsData::ekchipower[myindex] << "\t\t! --Power of chi(k), (default: "<<EXAFS_CHI2_POWER_DEF<<")" << endl;
		out << " " << CheckKey("BACKSCATT-FILE") << " = " << ExptsData::datafilename[expno +myindex+1] << "\t! --File name containing the E(k,r) coefficients, **MANDATORY ITEM**" << endl;
		//Get the r-range
		std::vector<int> mn(ntypes, -1), mx(ntypes, -1);
		std::vector<double> mnv(ntypes, -1.0), mxv(ntypes, -1.0);
		if (!ExptsData::FindLimitsExafsBackscattering(ngr + nsq + nfq + nfg+myindex,ExptsData::datafilename[expno + myindex + 1], &mn[0], &mx[0], &mnv[0], &mxv[0], ntypes)) CleanExit();
		bool found = !onlyndef;
		for (int i = 0; (i < ntypes) && (!found); i++) if (((!found) && (ExptsData::ek_rmin[(myindex)*ntypes + i] != mn[i])) || ((!found) && (ExptsData::ek_rmax[(myindex)*ntypes + i] != mx[i]))) found = true;
		if (found)
		{
			out << " " << CheckKey("R-RANGE") << " =";
			int last = 0;
			if ((onlyndef) && (ntypes > 1))
			{
				for (int i = ntypes - 2; (i >= 0) && (ExptsData::ek_rmin[(myindex)*ntypes + i] == ExptsData::ek_rmin[(myindex)*ntypes + i + 1]) && (ExptsData::ek_rmax[(myindex)*ntypes + i] == ExptsData::ek_rmax[(myindex)*ntypes + i + 1]); i--) last++;
			}
			for (int i = 0; i < ntypes - last; i++)
				out << endl << " " << rspacing[ngr + nsq + nfq + nfg + myindex] * (ExptsData::ek_rmin[(myindex)*ntypes + i] - 0.5 + binshift[ngr + nsq + nfq + nfg + myindex]) << " " << rspacing[ngr + nsq + nfq + nfg + myindex] * (ExptsData::ek_rmax[(myindex)*ntypes + i] - 0.5 + binshift[ngr + nsq + nfq + nfg + myindex]);
			out << "\t\t! --R-range in A (centre of bins) for the coeffs to use in the fit for each partials (default: the whole r range of data -- for the 1st partial: " << mnv[0] << ".." << mxv[0] << ");" << endl << "\t\t\t! when you provide less values than required, the last value passed to the remaining ones" << endl;
			out << "\t\t\t! alternative (example for the first partial only): " << CheckKey("R-POINT-RANGE")<<" = " << ExptsData::ek_rmin[(myindex)*ntypes] << " " << ExptsData::ek_rmax[(myindex)*ntypes] << " ! The applied point range to fit (default: 1.." << mx[0] << ", the whole range of the data)" << endl;
		}
		if ((!onlyndef) || (onlyndef && ((ExptsData::ek_dE0[myindex] != MAX_E0_SHIFT_DEF || ExptsData::ek_ngrid_in[myindex] != EXAFS_NGRID_DEF)))) out << " " << CheckKey("DELTAE0_NGRID") << " = " << ExptsData::ek_dE0[myindex] << " "<< ExptsData::ek_ngrid_in[myindex]<<"\t\t! --Maximum E0 shift in eV, number of grid points in one direction, (default: " << MAX_E0_SHIFT_DEF<<", "<< EXAFS_NGRID_DEF << ", no shift is performed)." << endl;
	}
	return;
}

void RunParams::CreateFreeCos(const int cosno, std::ofstream &out, const bool onlyndef)
//This method writes the cosno'th cosine angle distribution constraint related tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	int deftype=0;
	for (int i = 0; i < N_COS_METHOD; i++)
	{
		if (strcmp(cos_method[i], DISTRIB_TYPE_DEF) == 0)
		{
			deftype = i;
			break;
		}
	}
	out << endl << "[ " << CheckTag("COS") << " ]" << endl;
	if ((!onlyndef) || (onlyndef && (CosDistrConst::method[cosno] != deftype)))
	{
		out << " " << CheckKey("DISTRIB-TYPE") << " = ";
		if (CosDistrConst::method[cosno] == 3)
			out << CosDistrConst::datafilename[cosno];
		else
			out << cos_method[CosDistrConst::method[cosno]];
			
		out << "\t\t! --Target distribution: " << CheckKeyValue("DISTRIB-TYPE", DISTRIB_TYPE_DEF) << "(default) or " << CheckKeyValue("DISTRIB-TYPE", "UNIFORM") << " or " << CheckKeyValue("DISTRIB-TYPE", "ABSENT") << " or filename(to be read from)" << endl;
	}
	if (CosDistrConst::method[cosno] != 3)
	{
		out << " " << CheckKey("DISTRIB-DEGREES") << " = " << CosDistrConst::angle[cosno] << "\t\t\t! --Mean angle in degrees, **MANDATORY ITEM**" << endl;
		out << " " << CheckKey("DISTRIB-WIDTH") << " = " << CosDistrConst::wcontrol[cosno] << "\t\t\t! --Width parameter, **MANDATORY ITEM**, interpreted as" << endl;
		out << "\t\t\t\t\t! half-width in angles for " << CheckKeyValue("DISTRIB-TYPE", "UNIFORM") << " and sigma in cosine for " << CheckKeyValue("DISTRIB-TYPE", "GAUSSIAN") << " and " << CheckKeyValue("DISTRIB-TYPE", "ABSENT") << endl;
		if ((!onlyndef) || (onlyndef && (CosDistrConst::dcos_theta[cosno] != DCOSTHETA_DEF))) out << " DCOSTH = " << CosDistrConst::dcos_theta[cosno] << "\t\t\t\t! --Spacing in cos(theta) space, (default: " << DCOSTHETA_DEF << "). It should be the same for all COS constraints (except in the case of file)." << endl;
	}
	out << " " << CheckKey("CENT-TYPE") << " = " << CosDistrConst::central[cosno] + 1 << "\t\t\t\t! --Central atom type, **MANDATORY ITEM**" << endl;
	out << " " << CheckKey("NEIGH-TYPE_FROM_TO") << " = " << CosDistrConst::neighbour1[cosno] + 1;
	out << " " << CosDistrConst::dmin1[cosno];
	out << " " << CosDistrConst::dmax1[cosno] << "\t\t! --Neighbour type, rmin and rmax (in A units), **MANDATORY ITEM**" << endl;
	out << " " << CheckKey("NEIGH-TYPE_FROM_TO") << " = " << CosDistrConst::neighbour2[cosno] + 1;
	out << " " << CosDistrConst::dmin2[cosno];
	out << " " << CosDistrConst::dmax2[cosno] << "\t\t! --Neighbour type, rmin and rmax (in A units), **MANDATORY ITEM**" << endl;
	std::streamsize ss = out.precision();
	if (fabs(CosDistrConst::weights[cosno]) < pow(10, -ss))//would not fit into the frame
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::scientific, ios::floatfield);
	}
	int i = 0;
	if (cosno + ngr + nsq + nfq + nfg + nek == lead_series_ind - 1)
	{
		out << " " << CheckKey("SIGMA-MASTER") << " = " << CosDistrConst::weights[cosno];
		i = 1;
	}
	else
	{
		if (CosDistrConst::weights[cosno] < 0.0)
		{
			out << " " << CheckKey("SIGMA-SCALABLE") << " = " << -CosDistrConst::weights[cosno];
			i = 2;
		}
		else
			out << " " << CheckKey("SIGMA") << " = " << CosDistrConst::weights[cosno];
	}
	if (fabs(CosDistrConst::weights[cosno]) < pow(10, -ss))//set it back to fixed
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::fixed, ios::floatfield);
	}
	//The sigma description texts are stored beginning with lower case, as they are used for the combined sigmas as well, transform it here
	string mytext = sigma_text[i];
	std::transform(mytext.begin(), ++mytext.begin(), mytext.begin(), ::toupper);
	out << "\t\t! --" << mytext << endl << alternative_sigma_text[i] << endl;
	return;
}

void RunParams::CreateFreeCoord(const int cono, std::ofstream &out, const bool avcc, const bool onlyndef)
//This method writes the cono'th (if average -- avcc=true) coordination constraint related tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	out << endl << "[ " << (avcc ? "AV" : "") << "COORD ]" << endl;
	out << " CENT-TYPE = " << (avcc ? AvCoordConst::central[cono] : CoordNumbConst::central[cono]) + 1 << "\t\t\t\t! --Central atom type, **MANDATORY ITEM**" << endl;
	for (int i = 0; i < (avcc ? 1 : CoordNumbConst::n_neightype[cono]); i++)
	{
		out << " "<<CheckKey("NEIGH-TYPE_FROM_TO")<<" = " << (avcc ? AvCoordConst::neighbours[cono] : CoordNumbConst::neighbours[CoordNumbConst::cum_n_neightype[cono] + i]) + 1;
		out << " " << (avcc ? AvCoordConst::dmin[cono] : CoordNumbConst::dmin[CoordNumbConst::cum_n_neightype[cono] + i]);
		out << " " << (avcc ? AvCoordConst::dmax[cono] : CoordNumbConst::dmax[CoordNumbConst::cum_n_neightype[cono] + i]) << "\t\t! --Neighbour type, rmin and rmax (in A units), **MANDATORY ITEM**" << endl;
	}
	for (int i = 0; i < (avcc ? 1 : CoordNumbConst::n_subconst[cono]); i++)
	{
		out << " COORDNUM_" << (avcc ? "" : "FRACT_");
		int option = 0;
		if (i + (avcc ? CoordNumbConst::tot_subconst + cono : CoordNumbConst::cum_n_subconst[cono]) + ngr + nsq + nfq + nfg + nek + ncosdistr == lead_series_ind - 1)
		{
			out << "SIGMA-MASTER = ";
			option = 1;
		}
		else if ((avcc ? AvCoordConst::weights[cono] : CoordNumbConst::weights[CoordNumbConst::cum_n_subconst[cono] + i]) < 0.0)
		{
			out << "SIGMA-SCALABLE = ";
			option = 2;
		}
		else out << "SIGMA = ";
		///out << (avcc ? AvCoordConst::acnreq[cono] : CoordNumbConst::target_coord[CoordNumbConst::cum_n_subconst[cono] + i]);
		if (avcc) out << AvCoordConst::acnreq[cono];//DO NOT COMBINE DIFFERENT tYPES IN ONE OUT, it will write double for CoordNumb as well
		else out << CoordNumbConst::target_coord[CoordNumbConst::cum_n_subconst[cono] + i];
		if (!avcc) out << " " << CoordNumbConst::fraction[CoordNumbConst::cum_n_subconst[cono] + i];
		double sigma = (avcc ? AvCoordConst::weights[cono] : CoordNumbConst::weights[CoordNumbConst::cum_n_subconst[cono] + i]);
		std::streamsize ss = out.precision();
		if (fabs(sigma) < pow(10, -ss))//would not fit into the frame
		{
			out.unsetf(ios::floatfield);
			out.setf(ios::scientific, ios::floatfield);
		}
		out << " " << fabs(sigma);
		if (fabs(sigma) < pow(10, -ss))//set it back to fixed
		{
			out.unsetf(ios::floatfield);
			out.setf(ios::fixed, ios::floatfield);
		}
		out << (avcc ? "\t" : "") << "\t! --Desired average coordination number, " << (avcc ? "" : "fraction, ");
		
		if (option == 0) out <<sigma_text[0] << endl << alternative_sigma_text[0] << endl;
		else if (option == 1) out <<sigma_text[1] << endl << alternative_sigma_text[1] << endl;
		else out <<sigma_text[2] << endl << alternative_sigma_text[2] << endl;
	}
	if (!avcc)
		if ((!onlyndef) || (onlyndef && (CoordNumbConst::write_detail[cono] != WRITE_CNC_DETAIL_DEF))) out << " WRITE-CNC-DETAIL = " << CoordNumbConst::write_detail[cono] << "\t! --Whether to write detailed information for each central, (default: " << WRITE_CNC_DETAIL_DEF << ")." << endl;
	return;
}

#ifdef _ADVANCED_GEOM_CONST
void RunParams::CreateFreeConc(const int concno, std::ofstream &out, const bool onlyndef)
//This method writes the concno'th common neighbour constraint related tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	out << endl << "[ "<<CheckTag("CONC")<<" ]" << endl;
	
	out << " " << CheckKey("PTYPE1_PTYPE2_FROM_TO") << " = " << CommonNeighConst::primary1[concno] + 1 << " " << CommonNeighConst::primary2[concno] + 1;
	out<< " " << CommonNeighConst::dmin[concno]<<" "<<CommonNeighConst::dmax[concno]<<"\t\t\t\t! --Primary type 1, primary type 2, minimum and maximum distance between them, **MANDATORY ITEM**" << endl;
	for (int i = 0; i < CommonNeighConst::nsec_neigh[concno]; i++)
	{
		int index = CommonNeighConst::cumul_sec[concno] + i;
		out << " " << CheckKey("STYPE_FROM1_TO1_FROM2_TO2") << " = " << CommonNeighConst::secondary[index] + 1;
		out << " " << CommonNeighConst::dmin1[index] << " " << CommonNeighConst::dmax1[index];
		out << " " << CommonNeighConst::dmin2[index] << " " << CommonNeighConst::dmax2[index] << "\t\t! --Secondary type, minimum and maximum distance to primary1, minimum and maximum distance to primary2 (in A units), **MANDATORY ITEM**" << endl;
	}
	int option = 0;
	if (concno + ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord == lead_series_ind - 1)
	{
		option = 1;
		out << " " << CheckKey("COORDNUM_FRACT_SIGMA-MASTER");
	}
	else
	{
		if (CommonNeighConst::weights[concno] < 0.0)
		{
			option = 2;
			out << " " << CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE");
		}
		else
		{
			option = 0;
			out << " " << CheckKey("COORDNUM_FRACT_SIGMA");
		}
	}
	out << " = " << CommonNeighConst::target_coord[concno] << " " << CommonNeighConst::fraction[concno] << " ";
	std::streamsize ss = out.precision();
	if (fabs(CommonNeighConst::weights[concno]) < pow(10, -ss))//would not fit into the frame
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::scientific, ios::floatfield);
	}
	out << fabs(CommonNeighConst::weights[concno]);
	if (fabs(CommonNeighConst::weights[concno]) < pow(10, -ss))//set it back to fixed
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::fixed, ios::floatfield);
	}
	out << "\t! --Desired coordination number, fraction, " << sigma_text[option] << endl << alternative_sigma_text[option] << endl;
	
	return;
}

void RunParams::CreateFreeSnc(const int sncno, std::ofstream &out, const bool onlyndef)
//This method writes the sncno'th second neighbour constraint related tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	out << endl << "[ " << CheckTag("SNC") << " ]" << endl;

	out << " " << CheckKey("CTYPE_FTYPE_FROM_TO") << " = " << SecondNeighConst::central[sncno] + 1 << " " << SecondNeighConst::fneighbour[sncno] + 1;
	out << " " << SecondNeighConst::dmin1[sncno] << " " << SecondNeighConst::dmax1[sncno] << "\t\t\t\t! --Central type, first neighbour type, minimum and maximum distance between them, **MANDATORY ITEM**" << endl;
	for (int i = 0; i < SecondNeighConst::nsec_type[sncno]; i++)
	{
		int index = SecondNeighConst::nsec_type_cum[sncno] + i;
		out << " " << CheckKey("STYPE_FROM_TO") << " = " << SecondNeighConst::sneighbours[index] + 1;
		out << " " << SecondNeighConst::dmin2[index] << " " << SecondNeighConst::dmax2[index] << "\t\t! --Second neighbour type, minimum and maximum distance to first neighbour (in A units), **MANDATORY ITEM**" << endl;
	}
	int option = 0;
	if (sncno + ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord + ncommonneigh == lead_series_ind - 1)
	{
		option = 1;
		out << " " << CheckKey("COORDNUM_FRACT_SIGMA-MASTER");
	}
	else
	{
		if (SecondNeighConst::weights[sncno] < 0.0)
		{
			option = 2;
			out << " " << CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE");

		}
		else
			out << " " << CheckKey("COORDNUM_FRACT_SIGMA");
	}
	out << " = " << SecondNeighConst::target_coord[sncno] << " " << SecondNeighConst::fraction[sncno]<<" ";
	std::streamsize ss = out.precision();
	if (fabs(SecondNeighConst::weights[sncno]) < pow(10, -ss))//would not fit into the frame
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::scientific, ios::floatfield);
	}

	out<< fabs(SecondNeighConst::weights[sncno]);

	if (fabs(SecondNeighConst::weights[sncno]) < pow(10, -ss))//set it back to fixed
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::fixed, ios::floatfield);
	}
	out<< "\t! --Desired coordination number, fraction, " << sigma_text[option] << endl << alternative_sigma_text[option] << endl;
	return;
}
void RunParams::CreateFreeBvs(const int bvsno, std::ofstream &out, const bool onlyndef)
//This method writes the bvs'th bond valence sum constraint related tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	double *myR0, *myb;
	out << endl << "[ " << CheckTag("BVS") << " ]" << endl;
	if (BondValenceSumConst::is_set_def == 2)//freers, write the calculated values, if there is any
	{
		myR0 = BondValenceSumConst::R0;
		myb = BondValenceSumConst::b;
	}
	else
	{
		myR0 = BondValenceSumConst::R0_in;
		myb = BondValenceSumConst::b_in;
	}

	out << " " << CheckKey("CENT-TYPE_CHARGE") << " = " << BondValenceSumConst::central[bvsno] + 1;
	if (BondValenceSumConst::central_charge[bvsno] !=NO_IVALUE_DEF)
		out<< " " << BondValenceSumConst::central_charge[bvsno];
	out << "\t\t\t\t! --Central type, central charge (optional), **MANDATORY ITEM**" << endl;
	for (int i = 0; i < BondValenceSumConst::n_neightype[bvsno]; i++)
	{
		int index = BondValenceSumConst::cum_n_neightype[bvsno] + i;
		out << " " << CheckKey("NEIGH-TYPE_RMAX_R0_B_CHARGE") << " = " << BondValenceSumConst::neighbours[index] + 1;
		out << " " << BondValenceSumConst::dmax[index] << " " << myR0[index] << " " << myb[index];
		if (BondValenceSumConst::neigh_charge[index] != NO_IVALUE_DEF)
			out << " " << BondValenceSumConst::neigh_charge[index];
		out<< "\t\t! --Neighbour type, maximum distance (in A units), bond length R0 (in A, optional), b parameter (in A, optional),"<<endl;
		out<<"\t\t\t! neighbour charge for default value lookup (optional), **MANDATORY ITEM**" << endl;
		out << "\t\t\t! Give -1 for R0 and/or b, if default should be used!" << endl;
	}
	int option = 0;
	if (bvsno + ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord + ncommonneigh +nsecondneigh== lead_series_ind - 1)
	{
		option = 1;
		out << " " << CheckKey("VALENCE_SIGMA-MASTER");
	}
	else
	{
		if (BondValenceSumConst::weights[bvsno] < 0.0)
		{
			option = 2;
			out << " " << CheckKey("VALENCE_SIGMA-SCALABLE");

		}
		else
			out << " " << CheckKey("VALENCE_SIGMA");
	}
	out << " = " << BondValenceSumConst::target_valence[bvsno] << " ";
	std::streamsize ss = out.precision();
	if (fabs(BondValenceSumConst::weights[bvsno]) < pow(10, -ss))//would not fit into the frame
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::scientific, ios::floatfield);
	}

	out << fabs(BondValenceSumConst::weights[bvsno]);

	if (fabs(BondValenceSumConst::weights[bvsno]) < pow(10, -ss))//set it back to fixed
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::fixed, ios::floatfield);
	}
	out << "\t! --Desired valence, " << sigma_text[option] << endl << alternative_sigma_text[option] << endl;
	return;
}
#endif

void RunParams::CreateFreeCustom(std::ofstream &out)
//This method writes the custom moves option related tags into free control files defined the stream "out"
{
	out << endl << "[ "<<CheckTag("CUSTMOVE")<<" ]" << endl;
	out << " "<< CheckKey("NMOVED-ATOMS")<<" = " << nmoved << "\t\t! --Number of atoms moved in a single move, **MANDATORY ITEM**" << endl << "\t\t\t\t! The move parameters should be listed after each CUSTOM-LINE keywords line by line" << endl;
	if (cus_entry.size() > 0)
	{
		//Input from new format -- writing them out as they are, except NMOVED-ATOMS entry
		for (std::list<string>::iterator it = cus_entry.begin(); it != cus_entry.end(); it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			if (key != CheckKey("NMOVED-ATOMS")) out << *it << endl;
		}
	}
	else
	{
		//Input from old format
		ifstream inp(movefilename);
		if (!inp)
		{
			cout << "\nERROR: The custom move file \"" << movefilename << "\" doesn\'t exist!" << endl;
			CleanExit();
		}
		string line;
		bool rv = true;
		rv = static_cast<bool>(getline(inp, line));//Comment line
		if (rv) rv = static_cast<bool>(getline(inp, line));//Number of custom parameters 
		if (!rv)
		{
			cout << "\nERROR: At the reading of the custom move file \"" << movefilename << "\" in the following line: " << endl << line << endl;
			CleanExit();
		}
		stringstream ss(line);
		int no_params;
		if (!(ss >> no_params))
		{
			cout << "\nERROR: At the reading of the custom move file \"" << movefilename << "\" in the following line: " << endl << line << endl;
			CleanExit();
		}
		for (int i = 0; i < no_params; i++)
		{
			if (!getline(inp, line))
			{
				cout << "\nERROR: At the reading of the custom move file \"" << movefilename << "\" in the following line: " << endl << line << endl;
				CleanExit();
			}
			else
			{
				line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
				line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
				out << " "<<CheckKey("CUSTOM-LINE")<<" = " << line << endl;
			}
		}
		inp.close();
	}
	return;
}

void RunParams::CreateFreeSwap(std::ofstream &out)
//This method writes swap related tags into free control files defined the stream "out"
{
	out << endl << "[ "<<CheckTag("SWAP")<<" ]" << endl;
	out << " " << CheckKey("FRACT")<<" = " << swap_fraction << "\t\t\t\t! --Fraction of swap moves, **MANDATORY ITEM**" << endl;
	//Create array to store the setup
	std::vector<bool> swset((unsigned int)pow(ntypes, 2), false);
	for (int i = 0; i < ntypes; i++) swset[i*(ntypes + 1)] = true;
	for (int i = 0; i < nswap_pairs; i++)
	{
		swset[swap_type1[i] * ntypes + swap_type2[i]] = true;
		swset[swap_type2[i] * ntypes + swap_type1[i]] = true;
	}
	out << " " << CheckKey("PAIRS")<<" =";
	if (std::all_of(swset.begin(), swset.end(), [](bool value) {return value; })) out << " " << CheckKeyValue("PAIRS","ALL");
	else
	{
		//Looking for larger groups
		std::vector<int> to_clear;
		for (int i = 0; i < ntypes; i++)
		{
			if (std::all_of(swset.begin() + i * ntypes, swset.end() - (ntypes - 1 - i)*ntypes, [](bool value) {return value; }))
			{
				out << " " << i + 1 << "-ANY";
				to_clear.push_back(i);
			}
		}
		//Clearing the array
		for (std::vector<int>::iterator it = to_clear.begin(); it != to_clear.end(); it++)
		{
			for (int i = 0; i < ntypes; i++) swset[(*it > i) ? i * ntypes + *it : *it*ntypes + i] = false;
		}
		//Writing the remaining one out
		for (int i = 0; i < ntypes - 1; i++)
		{
			for (int j = i + 1; j < ntypes; j++) if (swset[i*ntypes + j]) out << " " << i + 1 << "-" << j + 1;
		}
	}
	out << "\t\t\t\t! --Atom pairs to swap, **MANDATORY ITEM**, syntax: list of keyword values (consecutively or on multiple lines without keywords):" << endl << "\t\t\t\t\t! "<<CheckKeyValue("PAIRS","ALL")<<": all type of atoms to swap; Type-ANY: all pairs containing atom type {Type} are swapped; Type1-Type2: only pairs with atom types Type1 and Type can be swapped" << endl;
	return;
}
void RunParams::CreateFreeBpot(std::ofstream &out, const bool onlyndef)
//This method writes bonded potential related tags (only those, which are normally read from *.dat) into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	out << endl << "[ BPOT ]" << endl;
	if ((n_top_def_option > 0) || (!onlyndef))
	{
		if (n_top_def_option == 0) out << " " << CheckKey("COMP-OPTION") << " = NONE" << endl;
		else
		{
			for (int i = 0; i < n_top_def_option; i++)
			{
				out << " "<<CheckKey("COMP-OPTION")<<" = ";
				out << string(top_def_option[i]).substr(2) << " ";
				out << "\t\t! --Compiler option for bonded potential for all of the topology file(s)" << endl;
			}
		}
	}
	int n=0,*p;
	//keep in mind, that lead pot series index start with 1 here
	if (potential == 0)
	{
		n = ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord;
		p = &lead_series_ind;
	}
	else
	{
		if (potential == 10)
			n = 2 * ChiSquared::nNB;
		else
		{
			if (potential == 1)
			{
				n = 2 * ChiSquared::nNB;
				if (Topology::npair_types > 0)
					n += 2 * ChiSquared::nNB;
			}
		}
		p = &lead_series_ind2;
	}
#ifdef _ADVANCED_GEOM_CONST
	n += ncommonneigh + nsecondneigh;
#endif
#ifdef _LOCAL_INV
	n += nlocint;
#endif
	if (*p > n && *p <= n + ChiSquared::nbonds)
	{
		out << " " << CheckKey("LEAD-BOND") << " = ";
		out << *p-n;
		out << "\t\t\t\t! --Index of the leading bond interaction in its own type starting with 1." << endl;
		out << "\t\t\t\t! If there is non-bonding potential this will be used to calculate the leading potential series index," << endl;
		out << "\t\t\t\t! if there is not, then the normal leading series index." << endl;
		out << "\t\t\t\t! alternatives: " << CheckKey("LEAD-ANGLE") << ", " << CheckKey("LEAD-PERDIH") << ", " << CheckKey("LEAD-HARMDIH") << ", " << CheckKey("LEAD-RBDIH") << endl;
		out << "\t\t\t\t! Same concept for the angles, periodic, harmonic and RB dihedrals" << endl;
	}
	else if (*p > n + ChiSquared::nbonds && *p <= n + ChiSquared::nbonds+ChiSquared::nangles)
	{
		out << " " << CheckKey("LEAD-ANGLE") << " = ";
		out << *p - n - ChiSquared::nbonds;
		out << "\t\t\t\t! --Index of the leading angle interaction in its own type starting with 1." << endl;
		out << "\t\t\t\t! If there is non-bonding potential this will be used to calculate the leading potential series index," << endl;
		out << "\t\t\t\t! if there is not, then the normal leading series index." << endl;
		out << "\t\t\t\t! alternatives: " << CheckKey("LEAD-BOND") << ", " << CheckKey("LEAD-PERDIH") << ", " << CheckKey("LEAD-HARMDIH") << ", " << CheckKey("LEAD-RBDIH") << endl;
		out << "\t\t\t\t! Same concept for the bonds, periodic, harmonic and RB dihedrals" << endl;
	}
	else if (*p > n+ ChiSquared::nbonds + ChiSquared::nangles && *p <= n + ChiSquared::nbonds + ChiSquared::nangles+ChiSquared::nperdihs)
	{
		out << " " << CheckKey("LEAD-PERDIH") << " = ";
		out << *p - n - ChiSquared::nbonds - ChiSquared::nangles;
		out << "\t\t\t\t! --Index of the leading periodic dihedral interaction in its own type starting with 1." << endl;
		out << "\t\t\t\t! If there is non-bonding potential this will be used to calculate the leading potential series index," << endl;
		out << "\t\t\t\t! if there is not, then the normal leading series index." << endl;
		out << "\t\t\t\t! alternatives: " << CheckKey("LEAD-BOND") << ", " << CheckKey("LEAD-ANGLE") << ", " << CheckKey("LEAD-HARMDIH") << ", " << CheckKey("LEAD-RBDIH") << endl;
		out << "\t\t\t\t! Same concept for the bonds, angles, harmonic and RB dihedrals" << endl;
	}
	else if (*p > n + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs && *p <= n + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs+ ChiSquared::nharmdihs)
	{
		out << " " << CheckKey("LEAD-HARMDIH") << " = ";
		out << *p - n - ChiSquared::nbonds - ChiSquared::nangles - ChiSquared::nperdihs;
		out << "\t\t\t\t! --Index of the leading harmonic dihedral interaction in its own type starting with 1." << endl;
		out << "\t\t\t\t! If there is non-bonding potential this will be used to calculate the leading potential series index," << endl;
		out << "\t\t\t\t! if there is not, then the normal leading series index." << endl;
		out << "\t\t\t\t! alternatives: " << CheckKey("LEAD-ANGLE") << ", " << CheckKey("LEAD-PERDIH") << ", " << CheckKey("LEAD-HARMDIH") << ", " << CheckKey("LEAD-RBDIH") << endl;
		out << "\t\t\t\t! Same concept for the bonds, angles, periodic and RB dihedrals" << endl;
	}
	else if (*p > n + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs && *p <= n + ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs + ChiSquared::nRBdihs)
	{
		out << " " << CheckKey("LEAD-RBDIH") << " = ";
		out << *p - n - ChiSquared::nbonds - ChiSquared::nangles - ChiSquared::nperdihs - ChiSquared::nharmdihs;
		out << "\t\t\t\t! --Index of the leading RB dihedral interaction in its own type starting with 1." << endl;
		out << "\t\t\t\t! If there is non-bonding potential this will be used to calculate the leading potential series index," << endl;
		out << "\t\t\t\t! if there is not, then the normal leading series index." << endl;
		out << "\t\t\t\t! alternatives: " << CheckKey("LEAD-BOND") << ", " << CheckKey("LEAD-ANGLE") << ", " << CheckKey("LEAD-PERDIH") << ", " << CheckKey("LEAD-HARMDIH")  << endl;
		out << "\t\t\t\t! Same concept for the bonds, angles, periodic and harmonic dihedrals" << endl;
	}
	
}

void RunParams::CreateFreeNbpot(std::ofstream &out, const bool onlyndef)
//This method writes non-bonded potential related tags into free control files defined the stream "out"
//If onlydef is true, it writes only the non-default values to out
{
	string text_S0 = "(combination rule = 0): LJ sigma parameters [in A] for each GROMACS partials";
	string text_S1 = "(combination rule = 1): C6 coefficient [in kJ*A^6/mol] for each GROMACS types, values for the partials calculated by geometric mean";
	string text_S2 = "(combination rule = 2): LJ sigma parameters for each GROMACS types [in A], values for the partials calculated by arithmetic mean";
	string text_S3 = "(combination rule = 3): LJ sigma parameters for each GROMACS types [in A], values for the partials calculated by geometric mean";
	string text_E0 = "(combination rule = 0): LJ epsilon parameters [in kJ/mol] for each GROMACS partials";
	string text_E1 = "(combination rule = 1): CN coefficient [in kJ*A^N/mol] for each GROMACS types, values for the partials calculated by geometric mean";
	string text_E2 = "(combination rule = 2): LJ epsilon parameters for each GROMACS types [in kJ/mol], values for the partials calculated by geometric mean";
	string text_E3 = "(combination rule = 3): LJ epsilon parameters for each GROMACS types [in kJ/mol], values for the partials calculated by geometric mean";
	out << endl << "[ "<<CheckTag("NBPOT")<<" ]" << endl;
	if ((!onlyndef) || (potential != POT_TYPE_DEF)) out << " " << CheckKey("NB-TYPE") << " = " << potential_types[potential] << "\t\t\t\t! -- Number of GROMACS types (default: number of RMC types)," << endl;
	int npartials = (ntypes*(ntypes + 1)) / 2;
	int nGRpartials = (nGRtypes*(nGRtypes + 1)) / 2;

	switch (potential)
	{
	case 1:
		if ((!onlyndef) || (nGRtypes != ntypes)) out << " " << CheckKey("NGR-TYPES") << " = " << nGRtypes << "\t\t\t\t! -- Number of GROMACS types (default: number of RMC types)," << endl;
		
		
		if ((!onlyndef) || (fabs(vdW_cutoff - CUTOFF_DEF > TOLERANCE2) || (fabs(Coulomb_cutoff - CUTOFF_DEF > TOLERANCE2)))) out << " "<<CheckKey("CUTOFF-VDW_COUL")<<" = " << vdW_cutoff;
		if ((!onlyndef) || (vdW_cutoff != Coulomb_cutoff)) out << " " << Coulomb_cutoff;
		if ((!onlyndef) || (fabs(vdW_cutoff - CUTOFF_DEF > TOLERANCE2) || (fabs(Coulomb_cutoff - CUTOFF_DEF > TOLERANCE2)))) out<< "\t\t! -- VdW and Coulomb cutoffs in A, (default: " << CUTOFF_DEF << ", " << CUTOFF_DEF << "), in reduced units, corresponding here to " << CUTOFF_DEF*SimpleCfg::boxedge<<" A"<< endl << "\t\t\t\t\t! if only the vdW cutoff is provided, the same value will be used for the Coulomb term, as well." << endl<< "\t\t\t\t\t! if only the factor for vdW provided, the same value will be applied for Coulomb term, as well." << endl;
		out.precision(12);
		out.setf(ios::scientific, ios::floatfield);
		if (vdW_comb_rule == 0)
		{
			out << " " << CheckKey("LJ-C6-EACH")<<" =";
			for (int i = 0; i < nGRpartials; i++)
				out << " " << vdW_pot1[i];
			
			out << "\t\t! --"<<text_S0<<", **MANDATORY ITEM**" << endl;
			out << "\t\t\t\t\t! alternatives: " << CheckKey("LJ-C6-GMEAN")<<" "<<text_S1 << endl;
			out << "\t\t\t\t\t! or " << CheckKey("LJ-SIG-AMEAN") << " "<<text_S2 << endl;
			out << "\t\t\t\t\t! or " << CheckKey("LJ-SIG-GMEAN") << " "<<text_S3<< endl;
			
			out << " " << CheckKey("LJ-CPOW-EACH") << " =";
			for (int i = 0; i < nGRpartials; i++)
			out << " " << vdW_pot2[i];
			out << "\t\t! --" << text_E0 << ", **MANDATORY ITEM**" << endl;
			out << "\t\t\t\t\t! alternatives: "<<CheckKey("LJ-CPOW-GMEAN") << text_E1 << endl;
			out << "\t\t\t\t\t! or "<<CheckKey("LJ-EPS-GMEAN")<<" " << text_E2 << endl;
			out << "\t\t\t\t\t! or "<< CheckKey("LJ-EPS-GMEAN")<<" " << text_E3 << endl;
		}
		else
		{
			//write sigma parameters
			switch (vdW_comb_rule)
			{
			case 1: out << " " << CheckKey("LJ-C6-GMEAN") << " = "; break;
			case 2: out << " " << CheckKey("LJ-SIG-AMEAN") << " = "; break;
			case 3: out << " " << CheckKey("LJ-SIG-GMEAN") << " = ";
			}
			for (int i = 0; i < nGRtypes; i++)
			{
				for (int j = i; j < nGRtypes; j++)
				{
					if (i == j)//only write the ii partials
						out<<" "<< vdW_pot1[i*nGRtypes + j - (i*(i + 1) / 2)];
				}
			}
			out << "\t\t! --";
			switch (vdW_comb_rule)
			{
			case 1: out << " " << text_S1 << " = "; break;
			case 2: out << " " << text_S2 << " = "; break;
			case 3: out << " " << text_S3 << " = ";
			} 
			out << ", **MANDATORY ITEM**" << endl;
			out << "\t\t\t\t\t! alternatives: " << CheckKey("LJ-C6-EACH") << " " << text_S0 << endl;//this is the first alternative for all
			switch (vdW_comb_rule)
			{
			case 1:
				out << "\t\t\t\t\t! or " << CheckKey("LJ-SIG-AMEAN") << " " << text_S2 << endl;
				out << "\t\t\t\t\t! or " << CheckKey("LJ-SIG-GMEAN") << " " << text_S3 << endl;
				break;
			case 2: 
				out << "\t\t\t\t\t! or " << CheckKey("LJ-C6-GMEAN") << " " << text_S1 << endl;
				out << "\t\t\t\t\t! or " << CheckKey("LJ-SIG-GMEAN") << " " << text_S3 << endl;
				break;
			case 3:
				out << "\t\t\t\t\t! or " << CheckKey("LJ-C6-GMEAN") << " " << text_S1 << endl;
				out << "\t\t\t\t\t! or " << CheckKey("LJ-SIG-AMEAN") << " " << text_S2 << endl;
			}
			
			//write epsilon parameters
			switch (vdW_comb_rule)
			{
			case 1: out << " "<<CheckKey("LJ-CPOW-GMEAN")<<" = "; break;
			case 2:
			case 3: out << " "<<CheckKey("LJ-EPS-GMEAN")<<" = ";
			}
			for (int i = 0; i < nGRtypes; i++)
			{
				for (int j = i; j < nGRtypes; j++)
				{
					if (i == j)//only write the ii partials
						out << " " << vdW_pot2[i*nGRtypes + j - (i*(i + 1) / 2)];
				}
			}
			out << "\t\t! --";
			switch (vdW_comb_rule)
			{
			case 1: out << " " << text_E1 << " = "; break;
			case 2: out << " " << text_E2 << " = "; break;
			case 3: out << " " << text_E3 << " = ";
			}
			out << ", **MANDATORY ITEM**" << endl;
			out << "\t\t\t\t\t! alternatives: " << CheckKey("LJ-CPOW-EACH") << " " << text_E0 << endl;//this is the first alternative for all
			switch (vdW_comb_rule)
			{
			case 1:
				out << "\t\t\t\t\t! or " << CheckKey("LJ-EPS-GMEAN") << " " << text_E2 << endl;
				out << "\t\t\t\t\t! or " << CheckKey("LJ-EPS-GMEAN") << " " << text_E3 << endl;
				break;
			case 2:
				out << "\t\t\t\t\t! or " << CheckKey("LJ-CPOW-GMEAN") << " " << text_E1 << endl;
				out << "\t\t\t\t\t! or " << CheckKey("LJ-EPS-GMEAN") << " " << text_E3 << endl;
				break;
			case 3:
				out << "\t\t\t\t\t! or " << CheckKey("LJ-CPOW-GMEAN") << " " << text_E1 << endl;
				out << "\t\t\t\t\t! or " << CheckKey("LJ-EPS-GMEAN") << " " << text_E2 << endl;
			}
		}
		
		//the weighting parameters depending on the combination mode
		//do not encourage not giving enough data and guessing, write out as many parameters, as it is expected!
		//leading series index will be proceeded with 'm' and scalable sigma with 's'
		if (abs(NB_weight_mode) == 2)
		{
			out << " " << CheckKey("SIGMA-OVERALLPOT") << " = " << ((lead_series_ind2 == 1) ? "m" : "") << ((vdW_weight[0] < 0.0) ? "s" : "") << fabs(vdW_weight[0]) << "\t\t\t! --Weighting parameter for all the non-bonded potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma, **MANDATORY ITEM**" << endl;
			out << "\t\t\t\t\t! alternatives: standard deviation(s) for VdW (" << CheckKey("VDW-SIGMA")<<": mandatory) and Coulomb (" << CheckKey("COUL-SIGMA")<<": optional, if not given, then " << CheckKey("VDW-SIGMA")<<" values(s) are used" << endl;
		}
		

		switch (NB_weight_mode)
		{
			case 3:
			case -3:
			case 0:
				out << " " << CheckKey("VDW-SIGMA") << " = " << ((lead_series_ind2 == 1) ? "m" : "") << ((vdW_weight[0] < 0.0) ? "s" : "") << fabs(vdW_weight[0]) << "\t\t\t! --Weighting parameter for all the LJ potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma, **MANDATORY ITEM**" << endl;
				if (!onlyndef || (vdW_weight[0] != Coulomb_weight[0]) || (lead_series_ind2 == 2))
				{
					out << " " << CheckKey("COUL-SIGMA") << " = " << ((lead_series_ind2 == 2) ? "m" : "") << ((Coulomb_weight[0] < 0.0) ? "s" : "") << fabs(Coulomb_weight[0]) << "\t\t\t! ---Weighting parameter for all the Coulomb potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma, **MANDATORY ITEM**" << endl;
					out << "\t\t\t\t\t! (when you provide less than no. of RMC partials, the last passed to the remaining), (default: same as VdW)" << endl;
				}
				out << "\t\t\t\t\t! alternative: (instead of " << CheckKey("VDW-SIGMA") << " and " << CheckKey("COUL-SIGMA") << ") " << CheckKey("SIGMA-OVERALLPOT") << ": weighting parameter for all the non-bonded potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma" << endl;
				break;
			case -1:
			case 1: out << " " << CheckKey("VDW-SIGMA") << " = ";
			{
				double *p = vdW_weight;
				int index = 1;
				for (int i = 0; i < ntypes; i++)
				{
					for (int j = i; j < ntypes; j++)
					{
						out << (lead_series_ind2 == index ? " m" : " ") << ((*p < 0.0) ? "s" : "") << fabs(*p);
						p++;
						index++;
					}
					if (i != ntypes - 1) out << endl;
					
				}

				out << "\t\t\t\t\t! ---Weighting parameters for npartials LJ potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma, **MANDATORY ITEM**" << endl;
				out << "\t\t\t\t\t! (when you provide less than no. of RMC partials, the last passed to the remaining)" << endl;
				bool same = true;
				for (int i = 0; i < npartials; i++)
					if (fabs(vdW_weight[i] - Coulomb_weight[i]) > TOLERANCE2) same = false;
				if (!onlyndef || !same || (lead_series_ind2 >= npartials))
				{
					out << " " << CheckKey("COUL-SIGMA") << " = ";
					p = Coulomb_weight;
					int index = 1;
					for (int i = 0; i < ntypes; i++)
					{
						for (int j = i; j < ntypes; j++)
						{
							out << (lead_series_ind2 == npartials+index ? " m" : " ") << ((*p < 0.0) ? "s" : "") << fabs(*p);
							p++;
							index++;
						}
						if (i != ntypes - 1) out << endl;
					}
				}
				out << "\t\t\t\t\t! --Weighting parameters for npartials Coulomb potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma, **MANDATORY ITEM**" << endl;
				out << "\t\t\t\t\t! (when you provide less than no. of RMC partials, the last passed to the remaining), if not given at all, the VDW-SIGMA values are used" << endl;
				out << "\t\t\t\t\t! alternative: (instead of " << CheckKey("VDW-SIGMA") << " and " << CheckKey("COUL-SIGMA") << ") " << CheckKey("SIGMA-OVERALLPOT") << ": weighting parameter for all the non-bonded potential terms, no prefix: standard deviation; 'm': leading potential series index; 's': scalable sigma" << endl;
			}
		}
		out.unsetf(ios::scientific);
		out.setf(ios::fixed, ios::floatfield);
		out.precision(6);
		if ((!onlyndef) || (pot_chi2_low_lim_fraction < POT_CHI_LOW_LIM_FRACTION_DEF))//it was set
		{
			out << " " << CheckKey("CHI-LOWLIMIT-RATIO") << " = ";
			if (pot_chi2_low_lim_fraction < POT_CHI_LOW_LIM_FRACTION_DEF)
				out << pot_chi2_low_lim_fraction;
			else out << "NONE";
			out << "\t\t! -- The chi^2 belonging to the potential term does not decrease below the limit given by this ratio of the initial\n\t\t\t\t\t! chi^2 of the potential term (default: NONE) -- this effects only, if non-bonded interactions are present" << endl;
		}
		if ((!onlyndef) || (LJ_rep_N != LJ_REP_N_DEF)) out << " "<<CheckKey("LJ-POW")<<" = " << LJ_rep_N << "\t\t\t\t! --Power of the repulsion term in 6-{POW} LJ (default: "<<LJ_REP_N_DEF<<")" << endl;
			
		if ((!onlyndef) || (vdW14_fudge != FUDGE14_DEF) || (Coulomb14_fudge != FUDGE14_DEF)) out << " " << CheckKey("1-4-SCALE-VDW_COUL")<<" = " << vdW14_fudge;
		if ((!onlyndef) || (vdW14_fudge != Coulomb14_fudge)) out << " " << Coulomb14_fudge;
		if ((!onlyndef) || (vdW14_fudge != FUDGE14_DEF) || (Coulomb14_fudge != FUDGE14_DEF)) out << "\t! --Factor(s) for scaling the  1-4 VdW and Coulombic parameters (default: "<<FUDGE14_DEF<<", "<<FUDGE14_DEF<<")," << endl << "\t\t\t\t\t! if only the factor for vdW provided, the same value will be applied for Coulomb term, as well." << endl;
		break;
	case 10:
		if ((!onlyndef) || (temperature != TEMPERATURE_DEF)) out << " " << CheckKey("TEMPERATURE") << " = " << temperature << "\t\t\t\t! --Temperature in K, no algorithmic meaning, for tabulated potential (default: "<<TEMPERATURE_DEF<<")," << endl;
		if ((!onlyndef) || (fabs(vdW_cutoff - CUTOFF_DEF > TOLERANCE2))) out << " " << CheckKey("CUTOFF-TAB") << " = " << vdW_cutoff<< "\t\t\t\t! --Cutoff for the tabulated potential(s) in A, (default: " << CUTOFF_DEF << "), in reduced units, corresponding here to " << CUTOFF_DEF * SimpleCfg::boxedge << " A" << endl ;
		
		
		int  j;
		for (int i=0;i<nused_potpartials;i++)
		{
			if (i == lead_series_ind2 - 1)
			{
				out << " " << CheckKey("IND_SIGMA-MASTER_FILE") << " = " ;
				j = 1;
			}
			else
			{
				if (vdW_weight[i] < 0.0)
				{
					out << " " << CheckKey("IND_SIGMA-SCALABLE_FILE") << " = ";
					j = 2;
				}
				else
				{
					out << " " << CheckKey("IND_SIGMA_FILE") << " = ";
					j = 0;
				}
			}
			out<<FNC_POT::tabpot_index[i]+1<<" "<< ((lead_series_ind2 == i+1) ? "m" : "") << ((vdW_weight[i] < 0.0) ? "s" : "") << fabs(vdW_weight[i])<<" "<< FNC_POT::tabpot_filename[i]<<"\t! --Index of the partial for this tabulated potential, "<<sigma_text[j]<<" file name for the tabulated potential"<<endl;
			
			out << "\t\t\t\t"<< alternative_sigma_text[j] << endl;
			
		}
	
	}
	return;
}

#ifdef _AENET
void RunParams::CreateFreeAenet(std::ofstream &out, const bool onlyndef)
//This method writes aenet related tags into free control files defined the stream "out"
{
	out << endl << "[ "<<CheckTag("AENET")<<" ]" << endl;
	for (int i=0;i<ntypes;i++)
		out << " " << CheckKey("IND_ANN-FILE_TYPE-NAME")<<" = " << i+1 << " "<<Aenet::ANN_filenames[i]<<" "<<Aenet::aenet_type_names[i] <<"\t\t\t\t! --Atom type index, ANN potential file, optionally the atom type name for AENET, **MANDATORY ITEM**" << endl;
	
	if ((!onlyndef) || (onlyndef && (write_E != WRITE_ATOMIC_ENERGY_DEF))) out << " " << CheckKey("WRITE-ATOMIC-ENERGY") << " = " << write_E << "\t! --Whether to write the energy for each atom into the *.en file (default: "<<WRITE_ATOMIC_ENERGY_DEF<<")" << endl;
	if ((!onlyndef) || (onlyndef && (!(Aenet::is_def_cutoff)))) out << " " << CheckKey("CUTOFF-AENET") << " = " << Aenet::cutoff << "\t! --The cutoff in A to used for AENET energy calculations (default: value used for the potential creation, "<<aenet_Rc_max<<" A)" << endl;
	if ((!onlyndef) || (onlyndef && (!(aenet_step!=AENET_STEP_DEF)))) out << " " << CheckKey("AENET-STEP") << " = " << RunParams::aenet_step << "\t! --The ANN potential is calculated by the AENET library in each AENET-STEP step (default: "<<AENET_STEP_DEF <<")"<< endl;
	std::streamsize ss = out.precision();
	if (fabs(aenet_weight) < pow(10, -ss))//would not fit into the frame
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::scientific, ios::floatfield);
	}
	int i=0;
	//cannot be any other potential conponenet, either master or scalable
	if (aenet_weight>0)
	{
		out << " " << CheckKey("SIGMA-MASTER") << " = " << aenet_weight;
		i = 1;
	}
	else
	{
		out << " " << CheckKey("SIGMA-SCALABLE") << " = " << -aenet_weight;
		i = 2;
	}
	if (fabs(aenet_weight) < pow(10, -ss))//set it back to fixed
	{
		out.unsetf(ios::floatfield);
		out.setf(ios::fixed, ios::floatfield);
	}
	//The sigma description texts are stored beginning with lower case, as they are used for the combined sigmas as well, transform it here
	string mytext= sigma_text[i];
	std::transform(mytext.begin(), ++mytext.begin(), mytext.begin(), ::toupper);
	out<< "\t! --" <<mytext << endl << alternative_sigma_text[i] << endl;



	return;
}
#endif
//in case of free format reading the topology related parameters are not known, so the lead bonded pot series index has to be calculated here
void RunParams::OffsetPotLead()
{
	//NB offset is needed for all the cases
	if (potential > 0)
	{
		if (potential == 1)
		{
			lead_series_ind2 += 2 * ChiSquared::nNB;
			if (Topology::npair_types > 0)
				lead_series_ind2 += 2 * ChiSquared::nNB;//there are npartials items for the LJ and Coulomb
		}
		else
		{
			if (potential == 10)
				lead_series_ind2 += ChiSquared::nNB;
		}
		switch (offset_lead2)
		{
		case 2:
			lead_series_ind2 += ChiSquared::nbonds;
			break;
		case 4:
			lead_series_ind2 += ChiSquared::nbonds + ChiSquared::nangles;
			break;
		case 8:
			lead_series_ind2 += ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs;
			break;
		case 16:
			lead_series_ind2 += ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs;
			break;
		}
	}
	else
	{
		lead_series_ind += ngr + nsq + nfq + nfg + nek + ncosdistr + navcoord + CoordNumbConst::tot_subconst;
#ifdef _ADVANCED_GEOM_CONST
		lead_series_ind += CommonNeighConst::nconstraints + SecondNeighConst::nconstraints;
#endif
#ifdef _LOCAL_INV
		lead_series_ind += nlocint;
#endif
		switch (offset_lead2)
		{
		case 2:
			lead_series_ind2 += ChiSquared::nbonds;
			break;
		case 4:
			lead_series_ind2 += ChiSquared::nbonds + ChiSquared::nangles;
			break;
		case 8:
			lead_series_ind2 += ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs;
			break;
		case 16:
			lead_series_ind2 += ChiSquared::nbonds + ChiSquared::nangles + ChiSquared::nperdihs + ChiSquared::nharmdihs;
			break;
		}
	}
}




#ifdef _LOCAL_INV
void RunParams::CreateFreeLoc(std::ofstream &out, const bool onlyndef)
//This method writes the local invariance related tags into free control files defined by the stream "out"
{
	int sigind;//serial number in comb_sigma and comb_alt_sigma arrays for explanation texts
	double *temp;
	string text, text2;
	//the number of data sets and constraint sofar
	int n = ngr + nsq + nfq + nfg + nek + ncosdistr + CoordNumbConst::tot_subconst + navcoord;
#ifdef _ADVANCED_GEOM_CONST
	n += ncommonneigh + nsecondneigh;
#endif

	out << endl << "[ " << CheckTag("LOCINV") << " ]" << endl;

	if ((!onlyndef) || (onlyndef && (loc_chi2_mode != LOC_CHI2_MODE_DEF))) out << " " << CheckKey("LOC-MODE") << " = " << (loc_chi2_mode == 0 ? CheckKeyValue("LOC-MODE", "BIN-BASED") : CheckKeyValue("LOC-MODE", "DISTANCE-BASED")) << "\t\t! --Local invariance calculation mode, can be "<<\
		CheckKeyValue("LOC-MODE", "BIN-BASED")<<" or "<< CheckKeyValue("LOC-MODE", "DISTANCE-BASED")<<" (default)" << endl;
	if ((!onlyndef) || (onlyndef && (loc_rspacing != rspacing_def))) out << " " << CheckKey("LOC-R-SPACING") << " = "<< loc_rspacing<< "\t\t! --Local invariance histogram bin size in A (default: first normal histogram bin size)" << endl;
	for (int i = 0; i < nlocint; i++)
	{
		if (n + i == lead_series_ind - 1)
		{
			out << " " << CheckKey("LOC-INT-FROM_SIGMA-MASTER") << " = ";
			sigind = 1;
		}
		else
		{
			if (loc_inv_sigma[i] < 0)
			{
				out << " " << CheckKey("LOC-INT-FROM_SIGMA-SCALABLE") << " = ";
				sigind = 2;
			}
			else
			{
				out << " " << CheckKey("LOC-INT-FROM_SIGMA") << " = ";
				sigind = 0;
			}
		}
		
		if (RunParams::loc_chi2_mode == 0)
			out << " " << min_loc_r;
		else
			out << " " << loc_at_ratio[i];
		std::streamsize ss = out.precision();
		if (fabs(loc_inv_sigma[i]) < pow(10, -ss))//would not fit into the frame
		{
			out.unsetf(ios::floatfield);
			out.setf(ios::scientific, ios::floatfield);
		}
		out<< " " << fabs(loc_inv_sigma[i])<<"\t\t! --Local invariance interval beginning, BIN-BASED method: minimum reduced distance for local invariance calculation,\n\t\t\t! (DISTANCE-BASED method: the fraction of neighbour atoms for the local invariance calculation at the beginning of the current interval);\n\t\t\t! " << sigma_text[sigind] << endl << alternative_sigma_text[sigind] << endl;
		if (fabs(loc_inv_sigma[i]) < pow(10, -ss))//set it back to fixed
		{
			out.unsetf(ios::floatfield);
			out.setf(ios::fixed, ios::floatfield);
		}
	}
	if (RunParams::loc_chi2_mode == 0)
		temp = &max_loc_r;
	else
		temp = loc_at_ratio+nlocint;
	if ((!onlyndef) || (onlyndef && (fabs(*temp - xmax_ori) > TOLERANCE2)))
		out << " " << CheckKey("LOC-INT-TO") << " = " << *temp << "\t\t! --Local invariance interval end, BIN-BASED method: minimum reduced distance for local invariance calculation,\n\t\t\t! (DISTANCE-BASED method: the fraction of neighbour atoms for the local invariance calculation at the end of the final interval)" << endl;
	return;
}
#endif

#ifdef _VIBR_AMP
void RunParams::CreateFreeVibramp(std::ofstream &out)
//This method writes the vibrationa amplitude related tags into free control files defined by the stream "out"
{
	out << endl << "[ " << CheckTag("VIBRAMP") << " ]" << endl;
	out << " " << CheckKey("GAUSSIAN-SIGMA") << " =";
	for (int i = 0; i < ntypes; i++)
		out << " " << T_corr_sigma_t[i];
	out << "\t\t! --Sigma parameters in A controlling the width of the Gaussian distribution for each atom type, **MANDATORY ITEM**" << endl;
	
	return;
}
#endif

#ifdef _NO_PERIODIC
void RunParams::CreateFreeNoper(std::ofstream &out, const bool onlyndef)
//This method writes the no periodic boundary condition related tags into free control files defined by the stream "out"
{
	out << endl << "[ " << CheckTag("NOPER") << " ]" << endl;
	out << " " << CheckKey("SAMPLE-RADIUS") << " = "<<R0<< "\t\t! --Radius of the spherical sample in A, **MANDATORY ITEM**" << endl;
		if ((!onlyndef) || (onlyndef && (recentre_flag != RECENTRE_FLAG_DEF))) out << " " << CheckKey("RECENTRE-SPHERE") << " = " << recentre_flag << "\t\t! --Whether to recentre the spherical sample in the box, 0: no, 1: yes, (default:0)" << endl;
	return;
}
#endif

