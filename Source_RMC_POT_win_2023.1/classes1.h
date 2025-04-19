//header file for 'primary' classes
//Last changed 06.02.2023

#ifdef _AENET
	#include "Aenet.h"
#else
	#include "Topology.h"
#endif
#include <algorithm>
#include <list>
#include <sstream>
#include <string>

//forward declarations
class RunParams;
class ChiSquared;
class Move;
class NeighbourList;
class FNC_POT;
class Threads;
struct ThreadArg;
#ifdef _NO_PERIODIC
	class DataMat;
#endif


//=======================================================================
class ExptsData
{
	//this class is used to store all the NECESSARY (for a defined set of
	//RMC run parameters given in a .dat file) experimental data
	//in memory for subsequent RMC use. In particular it allows the use
	//of arrays of data-constraints objects. The pointers to data information
	//in these objects (GofR,SofQ, etc...) should refer to data stored in
	//the ExptsData object of the RMC run
	friend class ChiSquared;
	public:
	ExptsData(RunParams &rundat);
	
	~ExptsData(); //destructor
	
	static int ntypes;//number of atom types
	static int ngr;//number of g(r) data files
	static int nsq;//number of S(Q) data files
	static int nfq;//number of F(Q) data files
	static int nfg;//number of F(g) data files
	static int nek;//number EXAFS data files
	static int ntot_datasets;//total number of data series
	static int ndiffbin;//number of different bin sizes
	static int *assign_hist;//assigning the different histogram bin sizes to the data sets dim [ntot_datasets]
	static int *use_cubic;//array to indicate to use cubic renormalization for the g(r), S(Q) and F(Q) data 
	static int *gruseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
	static double *rspacing;//pointer to Runparams::rspacing
	static char (*datafilename)[FILE_NAME_SIZE];//names of the experimental data files
	static int is_IQ;//indicator, whether I(Q) fitting is used at all
	static int is_E0shift;//indicator whether any EXAFS set uses E0_shift
	static int is_Xcoeff_calc;//indicator, whether Xray coeff calculation is used at all 
	static int is_IQmucorr;//indicator, whether there is a mu correction in any of the I(Q) fitting sets
	static int is_AXS;//indicator, whether AXS is used
	static int is_Ncoeff_calc;//indicator, whether neutron coeff calculation is used at all, 0 if not, 1 if yes before calculated, 2 after calculated

	//g(r) fitting
	static int  rused_tot;//total number of used r points
	static int *grsize;//number of data points in g(r) files
	static int *grused;//number of used data points in g(r) files
		//(i.e. we don't store data before nmix or after the nmax indices)
		//there is one size integer for each datafile
	static int *grmin;//indices of first data point to be used (START AT 1!)
	static int *grmax;//indices of last data point to be used
	static int *grrenorm;//renormalization switches for g(r) files
	static int *groffset;//offset switches for g(r) files
	static int *grlinear;//linear background switches for g(r) files
	static int *grquadratic;//quadratic background switches for g(r) files
	static int *grcubic;//cubic background switches for g(r) files
	static int *grused_cum;//cumulative number of USED g(r) data points in previous data sets; size:ngr+1
	static double *grsub;//constants to be substracted from data for g(r) files
	static double *grsigma;//sigma values  for g(r) files
	static double *gr_coeffs;//partials coefficients of g(r) files
	static int *bin_flag;//to indicate, whether the binsize and shift was reset based on the dr spacing in the experimental gr data set
	
	//S(Q) - neutron fitting
	static int sqused_tot;//total number of used S(Q) points
	static int *sqsize;//number of data points in S(Q) files
	static int *sqused;//number of USED data points in S(Q) files
	static int *sqmin;//indices of first data point to be used (START AT 1!)
	static int *sqmax;//indices of last data point to be used (START AT 1!)
	static int *sqrenorm;//renormalization switch for S(Q) files
	static int *sqoffset;//offset switches for S(Q) files
	static int *sqlinear;//linear background switches for S(Q) files
	static int *sqquadratic;//quadratic background switches for S(Q) files
	static int *sqcubic;//cubic background switches for S(Q) files
	static int *sqused_cum;//cumulative number of USED S(Q) data points in previous data sets; size:nsq+1
	static int *squseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
	static int *sqreadcoeffs;//wherther to read the coefficients instead of calculating them
	static double *sqsub;//constants to be substracted from data for S(Q) files
	static double *sqsigma;//sigma values for S(Q) files
	static double *sq_coeffs;//partials coefficients of S(Q) files
	static isotope_count_type *isotope_count;//[nfq*ntypes] number of isotopes for each data set's each component

	//F(Q) - X-ray fitting
	static int fqused_tot;//total number of used F(Q) points
	static int *fqsize;//number of data points in F(Q) files
	static int *fqused;//number of used data points in F(Q) files
	static int *fqmin;//indices of first data point to be used (START AT 1!)
	static int *fqmax;//indices of last data point to be used (START AT 1!)
	static int *fqrenorm;//renormalization switch for F(Q) files
	static int *fqoffset;//offset switches for F(Q) files
	static int *fqlinear;//linear background switches for F(Q) files
	static int *fqquadratic;//quadratic background switches for F(Q) files
	static int *fqcubic;//cubic background switches for F(Q) files
	static int *fqrenalpha;//compton alpha correction switches for F(Q) files
	static int *fqused_cum;//cumulative number of USED F(Q) data points in previous data sets; size:nfq+1
	static int *fqreadcoeffs;//wherther to read the coefficients instead of calculating them
	static int *fquseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
	static int *fqAXS;//0: no AXS, > 0 index(STARTING WITH 1) of the anomalous X - ray scattering type / data set
	static int *fqfprimeindex;//sum_fprime_itype(2**itype) for each set, indicating which types are using f'
	static int *fqfprimecount;//sum of how many types using f' for the given set
	static double *fqsub;//constants to be subsstracted from data for F(Q) files
	static double *fqsigma;//sigma values for F(Q) files
	static double *fqa;//a renormalization parameter in I(Q) fit
	static double *fqb;//a renormalization parameter in I(Q) fit
	static double *fqalpha;//alpha renormalization parameter in I(Q) fit
	static double *fqmu;//parameter for I(Q) background correction
	static double *fqdmumax;//max value to vary mu, between fqmu - fqdmumax -> fqmu + fqdmumax
	static double *fqmuact;//actual value of the mu parameterm between fqmu - fqdmumax -> fqmu + fqdmumax
	static double *fqfprime;//real part of the original f'(E) for each atom type
	static double *fqfprimeact;//real part of the actual f'(E) for each atom type
	static double *fqfprimefactor;//for varying f', values are in the range f'-f'*fqpfactor->f'+f'*fqpfactor
	static bool *fqfitIQ;//Fit I(Q) instead of F(Q), for each data set
	static bool *fqIQbackgcorr;//use background correction only in case of I(Q) fitting
	static bool *fqusecompton;//whether to add the Compton contribution during the IQ) calculation for the given data set
	

	//F(g) - electron diffraction fitting
	static int fgused_tot;//total number of used F(g) points
	static int *fgsize;//number of data points in F(g) files
	static int *fgused;//number of used data points in F(g) files
	static int *fgmin;//indices of first data point to be used (START AT 1!)
	static int *fgmax;//indices of last data point to be used (START AT 1!)
	static int *fgrenorm;//renormalization switch for F(g) files
	static int *fgoffset;//offset switches for F(g) files
	static int *fglinear;//linear background switches for F(g) files
	static int *fgquadratic;//quadratic background switches for F(g) files
	static int *fgcubic;//cubic background switches for F(g) files
	static int *fgused_cum;//cumulative number of USED F(g) data points in previous data sets; size:nfg+1
	static int *fguseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
	static double *fgsub;//constants to be subsstracted from data for F(g) files
	static double *fgsigma;//sigma values for F(g) files

	//EXAFS fitting//
	static int ekused_tot;//total number of used E(k) points
	static int ek_rused_tot;//total number of r points used in E(k) fitting
	static int *eksize;//number of data points in E(k) files
	static int *ekused;//number of used data points in E(k) files
	static int *ekmin;//indices of first data point to be used (START AT 1!)
	static int *ekmax;//indices of last data point to be used (START AT 1!)
	static int *ekmin_ori;//indices of first data point to be used originally (needed in case of E0 shift) (START AT 1!)
	static int *ekmax_ori;//indices of last data point to be used originally (needed in case of E0 shift) (START AT 1!)
	static int *ekrenorm;//renormalization switch for E(k) files
	static int *ekoffset;//offset switches for E(k) files
	static int *ekused_cum;//cumulative number of USED E(k) data points in previous data sets; size:nek+1
	static int *ekchipower;//weighting factor (E(k)*(k power of ekchipower)) for the experimental data points
	static int *ek_rused;//number of r points to use
	static int *ek_rmin;//index of first r point (histogram bin) to use
	static int *ek_rmax;//index of last r point (histogram bin) to use
	static int *ek_abstype;//type of the absorbing particle
	static int *ekuseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
	static int *ek_ngrid_in;//number of grid points to read from *.dat including the maximum for k-shift (0-: no shift, 1:only max_shift...)
	static int *ek_ngrid;//total number of grid points, ek_gridto-ek_gridfrom+1
	static int *ek_gridind;//the index of the actually chosen grid point during the E0 shift for each E(k) data series 
	static int *ek_gridstart;//the index of the starting grid point, default is zero (no shift)
	static int *ek_gridfrom;//index of the first grid point to use in case of assymetric grid (between -NGRID -> +NGRID)
	static int *ek_gridto;//index of the last grid point to use in case of assymetric grid (between -NGRID -> +NGRID)
	static int *ek_cf_cum;//cumulative number of coeff matrices for each data set  
	static int *ek_ngrid_cum;//cumulativ number of gridpoints for each set
	static double *eksigma;//sigma values for E(k) files
	static double *ek_dE0;//maximum E0-shift for both direction
	static double *ek_max_dkl;//maximum k-shift left for each data point
	static double *ek_max_dkr;//maximum k-shift right for each data point
	static double *ek_dE0_grid;//dE0 values of the dE0 grid
	static double *dE0_actual;//the actual dE0 values belonging to each grid point, this is  not equdistantly set 
	
	//data arrays
	double *gr_rvalues;//r values of g(r) files
	double *gr_gvalues;//g(r) values
	double *sq_qvalues;//Q values of S(Q) files
	double *sq_svalues;//S(Q) values
	double *fq_qvalues;//Q values of F(Q) files
	double *fq_fvalues;//F values or I values of F(Q) files, depending on what is fitted
	double *fq_IQvalues_ori=NULL;//original experimental I(Q) values of F(Q) files, only in case of mu correction
	double *fq_ffromIQrenvalues = NULL;//renormalized exp F(Q) calculated form I(Q), only in case of mu correction
	double *fq_coeffs;//Q-dependent coeffs of F(Q) files, can be read by the constructor, 
					//order Set1,part1:Q1-Qmax,part2:Q1-Qmax....part_last:Q1-Qmax
					//	    Set2,part1:Q1-Qmax,part2:Q1-Qmax....part_last:Q1-Qmax......
	double *fq_scatfactors=nullptr;//Q-dependent scattering factors of F(Q) files, can be read by the constructor, 
					//order Set1,type1:Q1-Qmax,type2:Q1-Qmax....type_last:Q1-Qmax
					//	    Set2,type1:Q1-Qmax,type2:Q1-Qmax....type_last:Q1-Qmax......
	double *fq_A=NULL;//<f^2> parameter for I(Q) fitting
	double *fq_B=NULL;//Compton parameter for I(Q) fitting
	double *fq_IB = NULL;//Q-dependent background for I(Q) fitting
	double *fg_gvalues;//g values of F(g) files
	double *fg_fvalues;//F values of F(g) files
	double *fg_coeffs;//g-dependent coeffs of F(g) files, will be read by the constructor, 
	//order Set1,part1:g1-gmax,part2:g1-gmax....part_last:g1-gmax
	//	    Set2,part1:g1-gmax,part2:g1-gmax....part_last:g1-gmax......		

	double *ek_kvalues_ori;//for the k values to read in
	double *ek_kvalues;//k values of E(k) files, in case of no E0 shift it is the same as ek_kvalues_ori, in case of E0 shift it will
	//contain the whole k-grid, the order is:
	//cycle data set, ngrid, k
	//FIRST SET:  k1[1]-kmax[1] for -ek_grid[1], k1[1]-kmax[1] for -ek_grid[1]+1,....., k1[1]-kmax[1] for -1, original k1[1]-kmax[1], 
	//            k1[1]-kmax[1] for +1,.....,k1[1]-kmax[1] for ek_grid[1]-1, k1[1]-kmax[1] for ek_grid[1],
	//SECOND SET: k1[2] - kmax[2] for - ek_grid[1], k1[2] - kmax[2] for - ek_grid[i] + 1, ....., k1[2] - kmax[2] for - 1, original k1[2]-kmax[2],
	//            k1[2]-kmax[2] for +1,.....,k1[2]-kmax[2] for ek_grid[2]-1, k1[2]-kmax[2] for ek_grid[2],
	//....

	double *pek_kvalues;//shows the beginning of the k1-kmax array part for the actual shift in ek_kvalues for the i-th set
	double *ek_evalues_ori;//original E values of E(k) files read from file
	double *ek_evalues;///E values of E(k) files, in case of no E0 shift it is the same as ek_evalues_ori, in case of E0 shift it will
	//contain the contracted number of datapoints left after interpolation
	//the ek_coeff matrix (matrices) are stored similarly as in case of the k values array, if there is E0 shift, then all the matrices
	//are stored in the same array in the loop order nek, ek_grid[i],ntypes, k,r:
	double *ek_coeffs;//(k,r)-dependent coeffs of E(k) files, if no E0 shift same as ek_coeffs ori, if there is E0 shift the grid matrices will stored here
	double *ek_coeffs_ori;//original (k,r)-dependent coeffs of E(k) files read by the constructor
	double *pek_coeffs;//similarly to pek_kvalues, this points to the beginning of the actually used (k,r) matrix in ek_coeffs array
	//contains the data in the order of cycles: data set, grid point, ntypes,k,r
	//the symmetric E0 shift results in assymetric, k-dependent k-shift 
	double *ek_maxklvalues;//maximum left k values for E(k) in case of E0 shift
	double *ek_maxkrvalues;//maximum right k values for E(k) in case of E0 shift
	
	

			
	//finder pointers
	double **gr_rfinder;//finder of r array of g(r) files
	double **gr_gfinder;//finder of g array of g(r) files
	double **sq_qfinder;//finder of Q array of S(Q) files
	double **sq_sfinder;//finder of s array of S(Q) files
	double **fq_qfinder;//finder of Q array of F(Q) files
	double **fq_ffinder;//finder of f array of F(Q) files
	double **fq_ffinder_ori=NULL;//finder of original f (meaning here I(Q)) array of F(Q) files
	double **fq_IBfinder=NULL;//finder of fq_IB array of F(Q) files in case of mu correction
	double **fq_ffromIQrenfinder = NULL;//finder of fq_ffromIBvalues array of F(Q) files in case of I(Q) fitting
	double **fq_cfinder;//finder of coeffs array of F(Q) files for each partials of each data set
	double **fq_sffinder=nullptr;//finder of scattering factor array of F(Q) files for each type of each data set
	double **fg_gfinder;//finder of g array of F(g) files
	double **fg_ffinder;//finder of f array of F(g) files
	double **fg_cfinder;//finder of coeffs array of F(g) files for each partials of each data set	
	double **ek_kfinder_ori;//finder of original k array of E(k) files, first to ek_kvalues_ori, then points to the beginning
							//of the original k segment for each data set in ek_kvalues array, if the original segment is not used due to E0GRID-FROM_TO...
							//setting, the whole, original ek_kvalues_ori array will be kept, and this pointer will point to the first used k point of the appropriate set
	double **ek_kfinder;//finder of k array of E(k) files, to ek_kvalues
	double **ek_efinder_ori;//finder of e array of E(k) files,  to ek_evalues_ori
	double **ek_efinder;//finder of E array of E(k) files, to ek_evalues
	double **ek_cfinder_ori;//finder of coeffs array of E(k) files, ek_coeffs_ori
	double **ek_cfinder;//finder for the beginning of the first coeffs array in case of E(k) for a k-grid point, ek_coeffs
	
	

	//member functions
	static void GetExptParam(ifstream &file);//getting the main experimental parameters
	void CalcNDCoeffs(int i);//Determining the neutron coefficients for a given set, only call for sets, were it is needed
	void CalcXrayCoeffs(int i, int mode=0);//calculate the X-ray coefficients for the given data set
	void CalcXrayCoeffsChange();//calculating the change in the coeffs  due to the change of f' in case of AXS
	void CalcCompton(int i);//calculate the Compton contribution for the given data set
	void NormEKCoeff();//normalize the EXAFS coefficients with the number of the edge particles
	void SetE0ShiftParams();//set the params for the E0 shift for EXAFS
	static void CheckgrSets();//checking the rspacing in the gr data files and resetting the rspacing in RunParams::rspacing to match the rspacing in the files, if necessary
	void Save() const;//saves on a file
	void CalcIQmucorr();//calculate the actual I(Q)_corr in the loop


	static bool GetExptParamFree(std::list<std::list<string> > pool);//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
	static bool FindLimitsExafsBackscattering(int expno, const char *filename, int *minlim, int *maxlim, double *minval, double *maxval, const int dimension);//Looks for the actual limit values given in filename: can be point range(min/maxlim), datarange (min/maxval) or using the whole data. If not mentioned, it should initialize with -1/-1.0
	
	
    private:
	static bool FindLimits(const char *filename, int &minlim, int &maxlim, double &minval, double &maxval);//Looks for the actual limit values given in filename: can be point range(min/maxlim), datarange (min/maxval) or using the whole data. If not mentioned, it should initialize with -1/-1.0
	

};	

//destructor of the ExptsData class
inline ExptsData::~ExptsData()
{
	if (ngr>0)
	{
		delete [] gr_rvalues;//r values of g(r) files
		delete [] gr_gvalues;//g(r) values
		delete [] gr_rfinder;//finder of r array of g(r) files
		delete [] gr_gfinder;//finder of g array of g(r) files
	}
	if (nsq>0)
	{
		delete [] sq_qvalues;//Q values of S(Q) files
		delete [] sq_svalues;//S(Q) values
		delete [] sq_qfinder;//finder of Q array of S(Q) files
		delete [] sq_sfinder;//finder of s array of S(Q) files
	}
	if (nfq>0)
	{
		delete [] fq_qvalues;//Q values of F(Q) files
		delete [] fq_fvalues;//f values of F(Q) files
		delete [] fq_coeffs;//Q-dependent coeffs of F(Q)
		delete [] fq_qfinder;//finder of Q array of F(Q) files
		delete [] fq_ffinder;//finder of F array of F(Q) files
		delete [] fq_cfinder;//finder of coeffs array of F(Q) files
		if (fq_scatfactors != NULL)// for I(Q), atomic scattering factor
			delete[] fq_scatfactors;
		if (fq_sffinder != NULL)// for I(Q), finder of atomic scattering factor
			delete[] fq_sffinder;
		if (fq_A != NULL)// I(Q) fitting parameter
			delete[] fq_A;
		if (fq_B != NULL)//I(Q) fitting parameter
			delete[] fq_B;
		if (fq_IB != NULL)//I() background correction
			delete[] fq_IB;
		if (fq_IBfinder != NULL)
			delete[] fq_IBfinder;
		if (fq_IQvalues_ori != NULL)
			delete[] fq_IQvalues_ori;
		if (fq_ffinder_ori != NULL)
			delete[] fq_ffinder_ori;
		if (fq_ffromIQrenvalues != NULL)
			delete[] fq_ffromIQrenvalues;
		if (fq_ffromIQrenfinder != NULL)
			delete[] fq_ffromIQrenfinder;
	}
	if (nfg > 0)
	{
		delete [] fg_gvalues;//g values of F(g) files
		delete [] fg_fvalues;//f values of F(g) files
		delete [] fg_coeffs;//g-dependent coeffs of F(g)
		delete [] fg_gfinder;//finder of g array of F(g) files
		delete [] fg_ffinder;//finder of F array of F(g) files
		delete [] fg_cfinder;//finder of coeffs array of F(g) files
	}
	if (nek>0)
	{
		delete [] ek_kvalues;//k values of E(k) files
		delete [] ek_evalues;//E values of E(k) files
		delete [] ek_coeffs;//(k,r)-dependent coeffs of E(k) files
		delete [] ek_kfinder;//finder of k array of E(k) files
		delete [] ek_efinder;//finder of E array of E(k) files
		delete [] ek_cfinder;//finder of coeffs array of E(k) files
	}
	if(::debug) cout<<"ExptsData destructor"<<endl;
};

//========================================================================
//-----class SimpleCfg-------------

//this is the simplified class for a configuration
//i.e. without molecular features and with a square box
class SimpleCfg
{
	public:
		SimpleCfg(); //default constructor
		SimpleCfg(SimpleCfg &source); //copy constructor

		~SimpleCfg(); //destructor
		static int ntotal;//total number of atoms (not including the virtual sites)
		static int ntypes;// number of atom types (not including the virtual sites)
		static int npartials;//number of partials (not including the virtual sites)
		static int ntotvirtual;//total number of virtual sites in the configuration
		static int nvirtualtypes;//number of virtual types in the configuration
		static int cfg_exist;//whether the *.cfg file exist
		static int nmoves_dim;//dim of the nmoves array

		static int *pnatoms;//pointer to number of atoms and virtual sites(per type)
		static int *cumul;//array of cumulative number of atoms and virtual sites per type
		//The Gromacs types will not necessary be continuous, if for example atoms from different molecule types will have
		//the GROMACS type. We have to know for each atom, which RMC and GROMACS type it belongs to.  n_sepGRtype>=nGRtypes,
		//the number of not continuos GROMACS type segement.
		static int n_sepGRtypes;//n_sepGRtypes>=nGRtypes
		static int *pGRatoms;//pointer to number of atoms (per n_sepGRtype)
		static int *cumul_GR;//array of cumulative number of atoms split to GROMACS types (dimension n_sepGRtypes)
		static int *GRsep_GR_type;//array for each to show the GROMACS type for each separate GROMACS type segment (dim: n_sepGRtypes)
		static int max_diff_virtual;//maximum number of virtual site an atom is involved in
		static int *ndiff_virtual;//number of virtual sites a given atom is involved in [ntotal]
		static int *virtual_ind;//the RMC atom indices of the virtual sites for this atom [max_diff_virtual*ntotal]
		static int *virtual_source;//the RMC indices ofthe atoms building the virtual site for each virtual site [ntotvirtual*4]
		static int *virtual_host_type;//RMC type of the host atom for the virtual sites, which decide, which partials the potentials contribute to
		double *charge_gr_centre;//the positions of the centre of the charge groups, needed for proper handling of the Coulomb cutoff
				//i.e (0,ntype1, ntype1 + ntype2,...,ntotal-nlasttype,ntotal))
		static longint **nmoves;//pointer to number of generated, tried and accepeted moves
		static double boxedge;//half-length of the box edge (in Angstrom)
		static double *fractions;//molar fractions
		
		double *positions;//pointer to array of (x,y,z) positions, after a new move is made, the NEW positions of the moved atom(s) are kept here!!!
						//[dimension: (ntotal+ntotvirtual)*3
		double **finder;//array of pointers to the start of the coordinates
					 	//for one atom or virtual type in the positions array
#ifdef _AV_MOVE		
		double *start_pos;//pointer to array containing the starting positions of the atoms, to be able to calculate the moved distances from the original positions
		double *moved_dist;//pointer to the array storing the moved distance of the atoms
#endif
		//member functions
		static int GetParamsCfg();//Reading the static parameters (number of atoms, types, timeframes) from the cfg file
												  //and initializing the static arrays pnatoms, cumul
												  //Setting the static parameters for the SimpleCfg object 
													//gives back 1, if it was succesfull, 0 if failed
		static void GetParamsBcf();//Reading the static parameters (number of atoms, types, timeframes) from the cfg file
												  //and initializing the static arrays pnatoms, cumul
												  //Setting the static parameters for the SimpleCfg object 
													//gives back 1, if it was succesfull, 0 if failed	
		static void SetGRtype();//setting up the arrays connected to the GROMACS types
		static void SetVirtual();//setting up the arrays connected to the virtual sites
		void ResizePosition();//in case of virtual sites resize the positions array and the finders
		void SetChargeGroupCentre(FNC_POT &pfnc);//setting the centre of the charge groups
		void CalcVirtualCoords();//calculate the coordinates of the virtual sites
		void CalcVirtualCoord(int ind);//creating the coordinates of the virtual sites, if there is any
		void PutToBox();
		void ConfigError(double pos, const char *file_name, const char *text);
		void Copy(SimpleCfg &source);//copy to the target
		void Save(const string title,const char *file_name="",int flag=0) const;//save to a file
		void Save_binary() const;//save the CombinedCfg object in binary format
		void Load();//loads from a text file
		int LoadBinary(longint **moves);//loading the data from a binary file,
										//gives back 1, if it was succesfull, 0 if failed
		int Comp(const SimpleCfg &conf) const;//compare the coordinates of the caller object with conf
#ifdef _AV_MOVE
		double CalcAvMove();//calculate the average moved distance of the atoms
		void SaveMovedDist();//saving the avrage move of the atoms to the avmfile
#endif
#ifdef _NO_PERIODIC
		static double boxedge_ori;//keep the original boxedge, it maybe needed for output
		void CheckSample();//check and rescale the sample if necessary, as the reduced coordinates have to be between -0.5 and 0.5
#endif
		
};
	
//destructor of the SimpleCfg class
inline SimpleCfg::~SimpleCfg()
{
	if (positions!=0)
	{
		//to make sure, that the program will not try to delete the array,
		//if it was already deleted 
		delete [] positions;
		positions=0;
	}
	if (finder!=0)
	{
		//to make sure, that the program will not try to delete the array,
		//if it was already deleted 
		delete [] finder;
		finder=0;
	}
#ifdef _AV_MOVE
	if (start_pos!=0)
	{
		//to make sure, that the program will not try to delete the array,
		//if it was already deleted 
		delete [] positions;
		positions=0;
	}
#endif
	if (charge_gr_centre!=NULL)
		delete [] charge_gr_centre;
	charge_gr_centre=0;
	if(::debug)  cout<<"SimpleCfg destructor"<<endl;
};
	
//==========CoordNumbConst==========
//this class implements the coordination number constraints
//i.e. constraints on the proportion of atoms having a given INTEGER
//coordination number for some neighbour at some distance
class CoordNumbConst
{
	public:
	CoordNumbConst(NeighbourList &neigh);//constructor
	CoordNumbConst(NeighbourList &neigh, CoordNumbConst &source);//Copy constructor
	
	~CoordNumbConst();//destructor
	static int nconstraints;//the number of constraints
	static int nthreads;//the total number of threads to use = RunParams::nthreads
	static int ncentral_tot;//total number of central atoms for all the constraints
	//if details are written, a separate neighbour list has to be calculated
	static int max_neigh;
	static int ncells_nlist;
	static int overlap_nlist;
	
	//all the following arrays are of dimension nconstraints
	//(i.e. one quantity per constraint)
	static int tot_neightype;//total number of neighbour types
	static int tot_subconst;//total number of subconstraints
	static int *central;//array of types of central particules
	static int *n_neightype;//number of neighbour types /constraint
	static int *neighbours;//array of types of neighbour particles
	static int *n_subconst;//number of subconstraint/constraint
	static double *dmin;//array of minimal distances 
	static double *dmax;//array of maximum distances 
	static double *weights;//array of weights for each CC 
	static int *target_coord;//array of DESIRED coord. numbers
	static double *fraction;//array of fraction of atoms desired 
					// with the coord. num.
	static double *udminsq;//array of reduced min. dist. squared
	static double *udmaxsq;//array of reduced max. dist. squared
	static double maxdistsq;//the maximum distance among the udmaxsq arrays, this will be the limit of the neighbourlist calculation if details is written
	


	static int *ncentral;//number of central atoms in the cfg
	static int *ncentral_cum;//cumulative number of central atoms before the given constraint for each constraint
	static longint *cctype;//array containing the partial indices in binary format
	static int *cum_n_neightype;//cumulative number of neighbour types
	static int *cum_n_subconst;//number of subconstraint/constraint

	static bool *write_detail;//array of whetehr to write detailed cnc information to cncd file
	static bool is_detail;//whether  there is detail for any constraint
	static int *neigh_list;//neighbour indices for each central of each constraint
	static int *neigh_type;//neighbour types for each central of each constraint
	static double *dist;//central-neighbour distance in A for each central of each constraint
	static int **neigh_finder;//finder of neigh_list for each constraint
	static int **type_finder;//finder of type_list for each constraint
	static double **dist_finder;//finder of dist_list for each constraint

	//The following members can differ from instance to instance
	int *nsatisfy;//number of central atoms satisfying the C.C 
	
	
	int *coordnumbs;//for each constraint, the coordination number for each central atom
	//if details are saved, it will be calculated separately only before saving and not updated, not during histogram calculation at each step
	
	
	int **finder;//finder for the beginnning of each constraint for the main thread in the coordnumbs array
	int **thread_finder;//finder for the beginning of each thread's segment in the coordnumbs array
	

	//handle for object
	//static Move *move;//this cannot be initialized with the constructor, has to be assigned later
	NeighbourList &neighlist;

	//member functions
	static void GetCoordConst(ifstream &file);//gets the parameters for the coordination const. from .dat file, and sets the static members
	static bool GetCoordConstFree(std::list<std::list<string> > pool);//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
	static void SetParams();
	static void ResizeCoordArrays(int **ind_list, int **type_list, double **dist_list);

	int &CoordinationNumb(int iconst, int iatom); //gets or sets the value for the coordnumbs array
	void CopyModified(CoordNumbConst &target, ThreadArg &thread_arg);//copy the modified coordination numbers
	int UpdateCoordnumb(ThreadArg &thread_arg); //Updating the coordnumbs array with the contribution of the threads, if there is any, and recalculating the number of atoms satisfying the constraint
	//iconst-th constraints, iatom-th central atom in the coordnumbs array
	int Load(double *sigma_percentage);//loading the coordination numbers, return 1, if it was successful, 0, if it was not
	void Save(const char *file_name="");//saving the CoordNumbConst object to a file
	void InitCoordnumbs();//sets the elements of array coordnumbs to 0
	void SaveDetail(const char *file_name="");//calculate and save the detials of the neighbours in the *.cncd file
	
	
};	

//destructor of the CoordNumbConst class
inline CoordNumbConst::~CoordNumbConst()
{
	if (nconstraints > 0)
	{
		delete[] nsatisfy;
		delete[] coordnumbs;
		delete[] finder;
		delete[] thread_finder;
		
	}
	if(::debug) cout<<"CoordNumbConst destructor"<<endl;
};

//==========AvCoordConst==========
//this class implements the average coordination number constraints

class AvCoordConst
{
	public:
	AvCoordConst();//default constructor
	AvCoordConst(AvCoordConst &source);//copy constructor
		
	~AvCoordConst();//destructor
	static int nconstraints;//the number of constraints
	static int  nthreads;//the total number of threads to use = RunParams::nthreads

	//all the following arrays are of dimension nconstraints
	//(i.e. one quantity per constraint)
	static int *central;//array of types of central particules
	static int *neighbours;//array of types of neighbour particles
	static int *avcctype;//array containing the partial index
	static double *dmin;//array of minimal distances 
	static double *dmax;//array of maximum distances 
	static double *weights;//array of weights for each CC 
	static double *acnreq;//array of average coordination numbers required
	
	static int *ncentral;//number of central atoms (in case of periodic only those, which are inside the R0-dmax)
	static double *udminsq;//reduced min. dist. squared
	static double *udmaxsq;//reduced max. dist. squared
#ifdef _NO_PERIODIC
	static double *dav_red;//array of average reduced distances 
	
#endif
	
	//the following members differ from instance to instance
	int *neighbourcount;//array of number of neighbours for the configurations per constraint
	
	//handle for object
	//static Move *move;//this cannot be initialized with the constructor, has to be assigned later

	//member functions
	int **finder;
	
	static void GetAvCoordConst(ifstream &file);
	static bool GetAvCoordConstFree(std::list<std::list<string> > pool);
	//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
	void InitNeighbourcount();//sets the elements of array neighbourcount to 0
	void Save(const char *file_name="") const;//save to a file
	int Load(double *sigma_percentage);//loads from a file
	//Copy the neighbourcount values for the constraints affected by the move from source to target
	void CopyModified(AvCoordConst &target);
	
};	

//destructor of the AvCoordConst class
inline AvCoordConst::~AvCoordConst()
{
	if (nconstraints > 0)
	{
		delete[] neighbourcount;
	}
	if(::debug)  cout<<"AvCoordConst destructor"<<endl;
};

class CosDistrConst
{
	//constraint for the bond angle defined by the neighbour1-central-neighbour2 triplett 
	public:
	CosDistrConst(NeighbourList &neigh_list, SimpleCfg &config);//default constructor
	CosDistrConst(CosDistrConst &source, NeighbourList &neigh_list, SimpleCfg &config);//copy constructor
		
	~CosDistrConst();//destructor
	static int nconstraints;//the number of constraints
	static int ncos_bin;//number of bins for method=0, 1, 2 constrains
	static int ntot_bins;//total number of bins
	static int *npoints;//number of bins for each contsraints;
	static int *cumul_bins;//cumulative number of bins for each constraint
	static int *max_neigh;//pointer to the maximum number of neighbopurs for the neighbourlist for an atom
	static longint *consttype;//pointer to NeighbourList::consttype, each bit represents a type, containing a 1 for each type that is involved in a constraint
	static longint cosconsttype;//each bit represents a type, containing a 1 for each type that is involved in a COS constraint

	//parameters for each constraint
	static int *method;//array of switches for the calculation method: 0: step function, 1: normal distribution
	static int *central;//array of types of central particules
	static int *neighbour1;//array of types of neighbour1 particles
	static int *neighbour2;//array of types of neighbour2 particles
	static double *dmin1;//array of minimal distances for the neighbour1
	static double *dmin2;//array of minimal distances for the neighbour2
	static double *dmax1;//array of maximum distances for the neighbour1
	static double *dmax2;//array of maximum distances for the neighbour2
	static double *angle;//array of desired bond angle
	static double *wcontrol;//array of parameters controlling the width of the distribution, its meaning depend on the calculation method
	static double *weights;//array of weights for each CC 

	static int *mod_const;//whteher the constraint was effected by the move
	static int *neigh_partial;//index of the partial calculated from the two neighbour types
	//for method=2 ("negative" constraint: no angle desired in a given region) the constraint is for the given region and not for all the bins, as it could interfere otherwise
	//with a "positive" constraint for the same particle types
	static int *first_bin;//index of the first bin, needed for the negative constraints
	static int *last_bin;//index of the last bin, needed for the negative constraints					
	static double *dcos_theta;//width of the cosine theta bin (full range goes -1 -> +1) for each constraint (can differ from the others for method=3 constraints)
	static double *inv_binwidth;// 1./dcos_theta for each constraint (can differ from the others for method=3 constraints)
	static double *min1sq;//reduced min. dist. squared
	static double *min2sq;//reduced min. dist. squared
	static double *max1sq;//reduced max. dist. squared
	static double *max2sq;//reduced max. dist. squared
	static double maxdistsq;//the maximum distance among the max1sq and max2sq arrays, this will be the limit of the neighbourlist calculation

	static double *theordistr;//array of theoretical or experimental distribiutions

	static char (*datafilename)[FILE_NAME_SIZE];//names of the experimental data files

	//the following members differ from instance to instance
	int *cos_hist;//bin counts for the calculated disribution
	double *cosinedistr;//array of cosine distribution of bond angles
	
	
	//handles to object
	NeighbourList &neighlist;
	SimpleCfg &config;

	//member functions
	static void GetCosDistrConst(ifstream &file);
	static bool GetCosDistrConstFree(std::list<std::list<string> > pool);
	//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
	static void SetParams();//initializing some of the static members
	void CalcTheoretical();//calculate the desired theoretical distribution
	void LoadCosDistr(int i);//loading the cosine distribution from a file, if option 3 is given
	void CalcHist();//calculate the initial histogram
	void UpdateHist(Move &move, int sign_switch);//updating the histogram
	void CalcDistribution();//calculating the cosine distribution of the bond angles
	void InitBincount();//sets the elements of array cos_hist to 0
	void Save(const char *file_name="") const;//save to a file
	//Copy the cos_hist values for the constraints affected by the move from source to target
	void CopyModified(CosDistrConst &target);

	
};	

//destructor of the CosDistrConst class
inline CosDistrConst::~CosDistrConst()
{
	if (nconstraints>0)
	{
		delete [] cos_hist;
		delete [] cosinedistr;
	}
	if(::debug)  cout<<"CordDistConst destructor"<<endl;
};

//====================================================================================
#ifdef _ADVANCED_GEOM_CONST
class CommonNeighConst
{
	//constraint for the common neighbours defined by the primary1-secondary1-primary2 (primary1-secondaryN-primary2) tripletts,
	//there can be nsec_neigh different secondary neighbour for a constraint

public:
	CommonNeighConst(NeighbourList &neigh_list, SimpleCfg &config);//default constructor
	CommonNeighConst(CommonNeighConst &source, NeighbourList &neigh_list, SimpleCfg &config);//copy constructor

	~CommonNeighConst();//destructor
	static int  nconstraints;//the number of constraints
	static int  tot_sec;//total number of secondary neighbours for all the constraints
	static int  tot_primary_pairs;//total number possible primary pairs
	static int  max_hist_size;//number of histogram bins(number of common neighbours+1)
	static int *nsec_neigh;//number of secondary neighbours for a costraint
	static int *cumul_sec;//cumulative number of secondary neighbours for each constraint
	static int *cumul_prim;//cumulative number of primary pairs for each constraint
	static int *max_neigh;//pointer to the maximum number of neighbours for the neighbourlist for an atom
	static longint *consttype;//pointer to NeighbourList::consttype, each bit represents a type, containing a 1 for each type that is involved in a constraint
	static longint comconsttype;//each bit represents a type, containing a 1 for each type that is involved in a CONC constraint


	//parameters for each constraint
	
	static int *primary1;//array of types of primary1 particles [nconstraints]
	static int *primary2;//array of types of primary2 particles [nconstraints]
	static int *secondary;//array of types of secondary neighbour particles  [tot_sec]
	static int *target_coord;//array of DESIRED coord. numbers [nconstraints]
	static double *dmin;//array of minimal distances between primary1 and primary2 [nconstraints]
	static double *dmax;//array of maximum distances between primary1 and primary2 [nconstraints]
	static double *dmin1;//array of minimal distances between primary1 and secondary [tot_sec]
	static double *dmin2;//array of minimal distances between primary2 and secondary [tot_sec]
	static double *dmax1;//array of maximum distances between primary1 and secondary [tot_sec]
	static double *dmax2;//array of maximum distances between primary2 and secondary [tot_sec]
	static double *weights;//array of weights for each constraint
	static double *fraction;//desired fraction for each constraint
	
	static int *mod_const;//whteher the constraint was effected by the move
	static int *primary_partial;//index of the partial calculated from the two primary types
	static double *minsq;//reduced min. dist. squared
	static double *min1sq;//reduced min. dist. squared
	static double *min2sq;//reduced min. dist. squared
	static double *maxsq;//reduced min. dist. squared
	static double *max1sq;//reduced max. dist. squared
	static double *max2sq;//reduced max. dist. squared
	static double maxdistsq;//the maximum distance among the maxsq, max1sq and max2sq arrays, this is used determining the limit of the neighbourlist calculation

	//the following members differ from instance to instance
	int *nprimary;//number of primary atom pairs in the given range for each constraint [nconstraints]
	int *nsatisfy;//number of central atoms satisfying the C.C [nconstraints]
	int *hist;//histogramm for the common neighbours
	
	//handles to object
	NeighbourList &neighlist;
	SimpleCfg &config;

	//member functions
	static void GetCommonNeighConst(ifstream &file);
	static bool GetCommonNeighConstFree(std::list<std::list<string> > pool);
	//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
	static void SetParams();//initializing some of the static members
	void CalcHist();//calculate the histogram for the common neighbours
	void UpdateHist(Move &move, int sign_switch);//updating the histogram 
	void SaveHist(const char *file_name = "") const;//save to a file
	//Copy the coordnumbs values for the constraints affected by the move from source to target
	void CopyModified(CommonNeighConst &target);


};


inline CommonNeighConst::~CommonNeighConst()
{
	if (nconstraints > 0)
	{
		delete [] nprimary;
		delete [] nsatisfy;
		
	}
	if (::debug)  cout << "CommonNeighConst destructor" << endl;
};


//========================================================================
class SecondNeighConst
{
	//constraint for the secondn neighbours defined by the primary1-secondary1-primary2 (primary1-secondaryN-primary2) tripletts,
	//there can be nsec_neigh different secondary neighbour for a constraint

public:
	SecondNeighConst(NeighbourList &neigh_list, SimpleCfg &config);//default constructor
	SecondNeighConst(SecondNeighConst &source, NeighbourList &neigh_list, SimpleCfg &config);//copy constructor

	~SecondNeighConst();//destructor
	static int  nconstraints;//the number of constraints
	static int  ncentral_tot;//total number of central atoms for all the constraints
	static int  nthreads;//the total number of threads to use = RunParams::nthreads
	static int  tot_sec;//total number of second neighbours for all the constraints
	static int *ncentral;//number of central atoms in the cfg
	static int *ncentral_cum;//cumulative number of central atoms before the given constraint for each constraint
	static int *nsec_type;//number of second neighbours for a costraint
	static int *nsec_type_cum;//cumulative number of second neighbours for each constraint
	static int *max_neigh;//pointer to the maximum number of neighbours for the neighbourlist for an atom

	static int *central;//array of types of central particles [nconstraints]
	static int *fneighbour;//array of types of first neighbour particles  [nconstraints]
	static int *sneighbours;//array of types of second neighbour particles  [tot_sec]
	static int *target_coord;//array of DESIRED coord. numbers [nconstraints]
	static double *dmin1;//array of minimal distances between central and first neighbour [nconstraints]
	static double *dmax1;//array of maximum distances between central and first neighbour [nconstraints]
	static double *dmin2;//array of minimal distances between first and second neighbours [tot_sec]
	static double *dmax2;//array of maximum distances between first and second neighbours [tot_sec]
	static double *weights;//array of weights for each constraint [nconstraints]
	static double *fraction;//desired fraction for each constraint [nconstraints]
	
	static double *min1sq;//reduced min. dist. squared
	static double *min2sq;//reduced min. dist. squared
	static double *max1sq;//reduced max. dist. squared
	static double *max2sq;//reduced max. dist. squared
	static double maxdistsq;//the maximum distance among the maxsq, max1sq and max2sq arrays, this is used determining the limit of the neighbourlist calculation

	static int *mod_const;//whteher the constraint was effected by the move
	static int *cp_partial;//index of the partial calculated from the central and first neighbour types

	static int involved_count;//number of involved indices
	static int *involved_indices;//indices of the atoms to recalculate for (the moved atoms and their first and second neighbours)

	static longint *consttype;//pointer to NeighbourList::consttype, each bit represents a type, containing a 1 for each type that is involved in a constraint
	static longint secconsttype;//each bit represents a type, containing a 1 for each type that is involved in a SNC constraint

	int *nsatisfy;//number of central atoms satisfying the C.C [nconstraints]
	int *coordnumbs;//coordination numbers  array (i.e. for each constraint, the coordination number for each central atom)
	
	int **finder;//finder for the beginnning of each constraint in the coordnumbs array

	//handles to object
	NeighbourList &neighlist;
	SimpleCfg &config;

	//member functions
	//static void GetSecondNeighConst(ifstream &file);//form vrsion 1.7 only the free format input is implemeted for the new features
	static bool GetSecondNeighConstFree(std::list<std::list<string> > pool);//get the freeformat parameters
	static void SetParams();//initializing some of the static members

	int &CoordinationNumb(int iconst, int iatom); //gets or sets the value for the coordnumbs array
	void CalcNeigh();//calculate the firts and second neighbours
	void UpdateNeigh(Move &move, int sign_switch);//update the firts and second neighbours
	void Save(const char *file_name = "");//save to a file (cannot be made const, because then CoordinationNumb cannot be used due to type mismatch)
	void InitCoordnumbs();//initialize the coordnumbs array
	void CopyModified(SecondNeighConst &target);//Copying the modified coordnumbs


};
//destructor of the SecondNeighConst class
inline SecondNeighConst::~SecondNeighConst()
{
	if (nconstraints > 0)
	{
		delete[] nsatisfy;
		delete[] coordnumbs;
		delete[] finder;
	}
	if (::debug)  cout << "SecondNeighConst destructor" << endl;
	
};

class BondValenceSumConst
{
	//constraint for the bond valence sum,
	//there can be nbvs_neigh different type if neighbours for a constraint
	//there can be n_subconst different valence states for a constraint

public:
	BondValenceSumConst();//default constructor
	BondValenceSumConst(BondValenceSumConst &source);//copy constructor

	~BondValenceSumConst();//destructor
	static int  nconstraints;//the number of constraints
	static int  nthreads;//the total number of threads to use, set by RunParams::GetParams
	static int  tot_neightype;//total number of neighbour types for all the constraints
			
	static longint *bvsctype;//array containing the partial indices in binary format


	//parameters for each constraint
	static int ncentral_tot;//total number of central atoms for all the constarint
	static int *central;//array of types of central particles [nconstraints]
	static int *neighbours;//array of types of neighbour particles  [nconstraints]
	static int *n_neightype;//number of bvs neighbours types for a costraint
	static int *ncentral;//number of central atoms for each constarint
	static int *ncentral_cum;//cumulative number of central atoms before the given constraint for each constraint
	static int *cum_n_neightype;//cumulative number of neighbour types
		
	static double *dmax;//array of maximum distances between central and neighbour [tot_bvs_neigh]
	static double *R0;//ideal bond length for each neighbour of each constraint [tot_bvs_neigh]
	static double *b;//b parameter for each neighbour of each constraint [tot_bvs_neigh]
	static double *R0_in;//ideal bond length for each neighbour of each constraint [tot_bvs_neigh] to keep it for *.free 
	static double *b_in;//b parameter for each neighbour of each constraint [tot_bvs_neigh] to keep it for *.free 
	static int *neigh_charge;//assumed oxidation state of the neighbour atom [tot_bvs_neigh] 
	
	static double *target_valence;//array of DESIRED valence [nconstraints]
	static double *weights;//array of weights for each subconstraint [nconstraints]
	static int *central_charge;//assumed oxidation state of the central atom [nconstraints]

	static int is_set_def;//set default R0 or b for any constraint
	static int *set_def;//whether to set default R0, b for a neighbour, 1 is added, if R0 should be set, 2 is added if b should be set [tot_bvs_neigh]
	
	static int *mod_const;//whteher the constraint was effected by the move
	
	static double *maxsq;//reduced min. dist. squared
	static double maxdistsq;//the maximum distance in the maxsq, this is used determining the limit of the neighbourlist calculation
	static string *pair_name;//the names with valence states for the atom pairs, if it is given [tot_bvs_neigh]  
	//The following members can differ from instance to instance
	double *valence;//current valence for each central atom of each constraint [nconstraints]

	double **finder;//finder for the beginnning of each constraint for the main thread in the valence array
	double **thread_finder;//finder for the beginning of each thread's segment in the valence array
	   
	//member functions
	static bool GetBondValenceSumConstFree(std::list<std::list<string> > pool);
	//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
	static void GetDefParams();//get the default R0 and b from the table
	static void CheckR0();//check whether the R0 is inside rmax
	void InitValence();//set the valence array elements to zero
	int UpdateValence(ThreadArg &thread_arg); //Updating the valence array with the contribution of the threads, if there is any
	
	//Copy the bond valence sum values for the constraints affected by the move from source to target
	void CopyModified(BondValenceSumConst &target, ThreadArg &thread_arg);
	double &BondValence(int iconst, int iatom);//gets or sets the valence
	int LoadBinary(double *sigma_percentage) const;//loading the valences from binary file, return 1, if it was successful, 0, if it was not
	void Save(const char *file_name = "");//save to a file
	void SaveBinary() const;//save to binary file
	


};
//destructor of the SecondNeighConst class
inline BondValenceSumConst::~BondValenceSumConst()
{
	if (nconstraints > 0)
	{
		delete[] valence;
		delete[] finder;
		delete[] thread_finder;
	}
	if (::debug)  cout << "BondValenceSumConst destructor" << endl;

};

#endif

//========================================================================
//-----class RunParams-------------

//this class contains the parameters needed for
//one RMC run ( contained in the .dat file, but also some information about
//the configuration and cell size, redundent with the SimpleCfg class)
class RunParams
{
	public:
	RunParams(); //constructor
	
	~RunParams(); //destructor
	
	//first the data needed for history records and other logs
	static string title;//title
	static int continuation;//the run should continue based on the configuration, *.dat and *.state file
	static int stateversion;//to identify the format of the state file, only the same version can be read
	static int printstep;//tries interval step for printing	
	static int histbuffsize;//size of the history buffer
	static int histstepratio;//number of saves between each history buffering
	static longint n_generated;//number of generated movesof generated swaps
	static longint n_swapgen;//number of swaps
	static longint n_E0shiftgen;//number of E0 shifts
	static longint n_mucorrgen;//number of generated I(Q) mu correction
	static longint n_fprimeshiftgen;//number of f' shifts
	static longint n_tried;//number of tried moves
	static longint n_accepted;//number of accepted moves
	static longint n_potaccepted;//number of moves accepted based on the chi2_pot 
	static longint n_swapacc;//number of accepted swaps 
	static longint n_E0shiftacc;//number of accepted E0 shifts
	static longint n_mucorracc;//number of accepted I(Q) mu correction
	static longint n_goodnonlinreg;//number of tried steps ended with successful non-linear regression
	static longint n_fprimeshiftacc;//number of accepted f' shifts
	static longint last_gen;//number of generated moves at last display
	static longint last_tried;//number of tried moves at last display
	static longint last_accepted;//number of accepted moves at last display
	static longint last_swap;//number of generated swaps at last display
	static longint last_swapacc;//number of accepted swaps at last display
	static longint last_E0shift;//number of generated E0 shifts at last display
	static longint last_E0shiftacc;//number of accepted E0 shifts at last display
	static longint last_mucorracc;//number of accepted mu correction steps at last display
	static longint last_mucorr;//number of generated mu correction steps at last display
	static longint last_fprimeshift;//number of generated f' shifts at last display
	static longint last_fprimeshiftacc;//number of accepted f' shifts at last display
	
	static double runlimit;//upper running time limit for the run in minutes  or to runlimit accepted steps depending on the runmode, or if negative value, it means the number of steps to generate
	static double timesave;//time interval for saving
	static string lead_series_name;//name of the lead series
	static string lead_series_name2;//name of the lead potential series

	//second the data that play an algorithmic role
	static int dat_version;//the version of the *.dat file
	static int runmode;//indicator of program termination
	//		0: default, timelimit in minutes was given for the run
	//		1: runnig to -runlimit generated steps, either in _TEST_MODE or in normal mode (- sign is for compatibility with fixed format input) 
	//		2: running to runlimit naccepted steps either in _TEST_MODE or in normal mode (in fixed format reading with TEST_MODE positive runlimit will mean running to last acc to be compatible with free format)
	static int nthreads;//total number of threads to use
	static int custmove;//custom move indicator
	static int nmoved;//number of atoms to move in a single move
	static int reload;//boolean switch, whether to load histogram from .hgm file and
						//coordination numbers (if there is coordination constraints) 
						//for CoordNumbConst: .cnc file, AvCoordConst: .acn file
	static int moveout;//this is a boolean switch indicating whether to use the moveout option
	static int Rfit;//this is a boolean switch indicating whether to fit the Rfactor instead of chi2 default:0
	static int fnc;//this is the fnc switch
	static int max_gridatom;//maximum number of atoms in a grid cell
	static int nswap_pairs;//number of mixed partials, where swap can occure
	static int niter_nonlin;//max number of iteration in case of non-linear regression 
	static int terminate_nonlin;//Terminate after this many successive moves ended with failed non-linear regression
	static int IQ_backg_corr_step;//change the mu for I(Q) background correction for in each IQ_backgcorr_step. step
	static int AXS_shiftstep;//change f' for AXS I(Q) in each AXS_shiftstep. step
	static int exafs_shiftstep;//shift the k values in each exafs_shiftstep step with EXAFS data sets with valid MAX-KSHIFT_NGRID
	static int *swap_type1;//array containing the type of the first atom of the allowed swap pairs
	static int *swap_type2;//array containing the type of the second atom of the allowed swap pairs
	static int *firstbin;//the first histogram bin to use, starts with zero
	static double rho;//the number density
	static double rspacing_def;//the default bin size
	static double *rspacing;//array dim:[ntot_datasets] of the width of the histogram bins as given for each data sets,as  there can be different bin sizes for different data sets
	static double *rspacing_diff;//array dim:[ndiffbin] of the width of the different histogram bins
					//vdW: it has a minimum value related to Q range
	static double *gr_rspacing_ori;//for the g(r) data sets, keep the original rpacing values, as it can be reset, if the gr data set's spacing is different, only for pronting run params
	static double too_close_fraction;//the fraction of moves for too close atoms in case of moveout option
	static double swap_fraction;//fraction of swaps
	static double *pcutoff;//the distances of closest approach
	static double *pmaxmove;//the maximum moves (1 per atom type)
	static double epsilon_nonlin;//stop non-linear iteration, when change in parameter is less than this
	static double lambda_nonlin;//fudge for nonlinera regression
	static double factor_nonlin;//factor to change lambda during nonlinear regression
	static string *chem_symbols;//chemical symbol for each type given in *.dat, using the symbols in the Xray_coeffs_W table 
	static string *chem_symbols_standard;//standard chemical symbol for each type extracted from chem_symbols, to find them in Compton and N_scat_length
		
	int fnc_conflict;//indicator between the conflict of the cutoffs and FNC
	static double binshift_def;//default binshift, this is read from the *.datfile
	static double *binshift;//number of bins to leave out from the calculation at the beginning of the histogram for a given data set
							//can differ, if gr data sets applied, and the automatically calculated binshift not the same as the default. 
							//for example , if for a data set binshift=1.5, then it will belong to the binshift_diff=0.5 histogram, and will start
							//with the second bin of it, set by firstbin offset
	static double *binshift_diff;//for the calculation of xmin, for the different binsize-binshift combination, dim [ndiffbin], this is always below 1!
								//always refers to the histogram beginning, not to the actual histogram parts, used for the data sets
	static int 	auto_cutoff;//whether to determine the cutoffs automatically:
							// 0:no
							// 1:yes, for flexible molecules all distances considered for minimum distance dteremination, similarly as for atomic systems 
							//-1:continuation of a 1 autocutoff run, read it from state
							// 3:yes, for flexible molecules bond and 1-3 angle distances not considered at min distance determination
							//-3:continuation of a 3 autocutoff run, read it from state
							// 2:already determined or read for autocutoff run,  
							
	static bool old_out;//True if an old .out file should be created
	static bool sum_ppcf;//Creates summed and averaged ppcf's at configuration collection
	static bool binshift_flag;//Indicating, that there are custom bin shift(s)
	static bool use_custom_sftable;//Indicating, that custom sfactortable should be used
	static bool write_density;//if the simulation box was rescaled, write the density in the *.freer(s) files
	double *pcutsq;//the distances of closest approach SQUARED

	//potential related
	static int potential;//indicator if there is a potential used
	static int offset_lead2;//indicator how to offset the lead_series_ind2 in case of freeformat reading
	static int nused_potpartials;//number of partials to use for tabulated potential (number of potential files)
	static int LJ_rep_N;//power for LJ repulsion term in 6-N LJ
	static int nGRtypes;//number of different GROMACS types
	static int NB_weight_mode;	//0: the same weighting parameter will be used for all the partials in NB interaction, if there is any, 
								//only one has to be given, 1: different weighting will be used for each partial, npartials weighting parameter has to be given
								//2: the vdW_weight[0] will be used for all potential related contribution, only this is read
	static int vdW_comb_rule;//combination rule for the creation of the mixed partial's potential parameters 
	static int max_option;//maximum number of define options for topology
	static int n_top_def_option;//number of define options for topology
	static double vdW14_fudge;//this factor will scale the vdW potential, to be able to handle the 1-4 interactions
						//if (1.0), normal vdW interaction, real exclusion, the whole potential previously calculated has to be subtracted
						//if (0<fudge<1.0) the whole pot is subtracted from the normal potential, and fudge*potential is added to the 14pot
	static double Coulomb14_fudge;//this factor will scale the Coulomb potential, similarly to vdW14_fudge
	static double *vdW_pot1;//LJ_sigma in A or C(6)=4*epsilon*(sigma^6) in kJ/mol*A^6 or Buckingham A
	static double *vdW_pot2;//LJ_epsilon in kJ/mol or C(N)=4*epsilon*(sigma^N) in kJ/mol*A^N or Buckingham B
	static double *vdW_weight;//weighting parameter for the vdW interaction (per partials or in case of NB_weight_mode=0 or 2, there is 1 weight parameter for all the partilas)
							//this is used for tabulated potential as well, always weight mode=1 is assumed
	static double *Coulomb_weight;//weighting paramaters for the Coulomb interaction (per partials or in case of NB_weight_mode=0 or 2, there is 1 weight parameter for all the partilas)
	static double Coulomb_cutoff;//cutoff value for the Coulomb intercation (A)
	static double vdW_cutoff;//cutoff value for the vdW intercation (A)
	static double pot_chi2_low_lim_fraction;//the fraction of the initial total chi2 wich would serve 
										//as the lower limit, which the pot chi2 could not go under
	double Coulomb_cutoff_sq;//squared cutoff value for the Coulomb intercation (reduced)
	double vdW_cutoff_sq;//squared cutoff value for the vdW intercation (reduced)
	static double temperature;//for the tabulated potential, only for output
	static char  vdW_name[15];//name of the vdW potential 
	static char(*top_def_option)[FILE_NAME_SIZE];//define option for the topology file, governing which part of it is included
#ifdef _AENET
	static int	aenet_step;//perform ANN pot calculation in each aenet_step step
	static double aenet_cutoff;//cutoff for the ANN potential
	static double aenet_weight;//weight parameter
	static bool write_E;//whether write the atomic energies to enerfilename
	static bool do_relax;//whether do structure relaxation for the atoms inside cutoff around the moved atom
#endif
#ifdef _LOCAL_INV
	static int nlocint;//number of local invariance intervals
	static int loc_chi2_mode;//how to calculate the local invariance chi2, 0 based on bins, 1 based on distance
	static int *loc_npoints;//number of data points for local invariance (depend on the calculation mode)
	static double *loc_inv_sigma;//weight for the local invariance (intervals)
	static double *loc_at_ratio;//ratio of atoms at the interval bordersfor  the local invariance calculation in case of loc_chi2_mode=1 (nlocint+1)
#endif
#ifdef _USE_LOCAL_INV
	static int min_loc_nbins;//first normal histogram bin index to use in the local invariance 
	static int max_loc_nbins;//number of bins to use for local invariance calculation 
	static double loc_rspacing;//the width of the local invariance histogram bins
	static double min_loc_r;//minimum distance for the local invariance calculation in reduced unit
	static double max_loc_r;//maximum distance for the local invariance calculation in reduced unit
#endif
#ifdef _NO_PERIODIC
	static int recentre_flag;//indicator, whether recentre the simulation box
	static double R0;//the radius of the speric sample inside the box
#endif	
#ifdef _VIBR_AMP
	static double *T_corr_sigma_t;//sigma parameter for the Gaussian convolution function for each types
#endif
	//third some information about the configuration
	//!this appears in the configuration files not in the .dat files!
	static int natoms;//the number of atoms in the cell
	static int ntypes;//the number of atom type (NECESSARY for memory allocation)
	static int cfgnumb;//number of configurations to collect after convergence
	static int coll_frequency;//collect the configuration at each coll_frequency save, if 0, not collected
	static int ndiffbin;//number of different bin sizes
	static int *pnatoms;//number of atoms per type	
	static int *assign_hist;//assigning the different histogram bin sizes to the data sets dim [ntot_datasets]
	static int *nbins;//the number of histogram bins
	static double boxedge;//the half-length of the simulation cell in real units
	static double xmax_ori;//the original maximum value for interatomic distances
	static double *xmax;//array of the maximum value for interatomic distances in reduced units, dim [ndiffbin]

	static double *xmin;//the minimum value for the histogram's interatomic distances in reduced units, dim [ndiffbin]
	double rmax;//the maximum value for interatomic distances
				//related to the size (half edge) of the cell-box.
	
	//fourth information related to the datasets
	static int ntot_datasets;// total number of data sets
	static int ngr;//number of g(r) constraints
	static int nsq;//number of S(Q) - neutron data - constraints
	static int nfq;//number of F(Q) - x-ray data - constraints
	static int nfg;//number of F(g) -electron diffraction data - constraints
	static int nek;//number of EXAFS constraints
	static int ncosdistr;//number of constraint for the cosine distribution of bond angles
	static int ncommonneigh;//number of common neighbour constraints
	static int nsecondneigh;//number of second neighbour constraints
	static int navcoord;//number of average coordination constraints
	static int nicoord;//number of individual coordination constraints
	static int nbvs;//number of bond valence sum constraints
	static int nbond_types;//number of bond types
	static int nangle_types;//number of angle types
	static int nperdihedral_types;//number of periodic dihedrals types
	static int nharmdihedral_types;//number of harmonic dihedrals types
	static int nRBdihedral_types;//number of RB dihedrals types
	static int tot_used14part;//total number of used 1-4 partials
	static int *used_14partials;//array containing 1 for the partials used in 1-4 interactions, 0 otherwise (not to display the unused partials)
	static int ntotal_points2;//total number of data points used for output in case of non-bonded interactions for the sum of the bonded and non-bonded interactions
	static int lead_series_ind2;//index of the leading series for the potential sigma scaling, starting with the vdW interaction
	static int ntotal_points;//total number of data points used for output
	static int lead_series_ind;//index of the leading series
	static std::list<string> cus_entry;//Entries for custom moves

	
	
	//member functions
	static void CheckSigma(double *sigma, const char *name);//check, whether it is not zero, and stop, if it is
	static void ReadLogical(ifstream &file, int &param, const char *param_name, const char *routine_name, int mode=0);//reading logical parameters, can be given in new and old format
	static void ReadLogicalS(stringstream &string, int &param, const char *param_name, const char *routine_name);//reading logical parameters, can be given in new and old format
	static void GetParams();//read data from the .dat file
	static void GetFixedFormatParams(ifstream &file);//read fixed format parameters
	static void GetFreeFormatParams(ifstream &file);//read free format parameters
	static void SetParams(int n_atoms,int n_types,double box_edge, int *p_natoms);//setting the static parameters related to the configurations
	static void SetBondedParams();//setting the parameters for the bonded interactions
	static void AdjustBoxedge();//This function rescales the boxedge
	static void SetPoints();//set the total number of data points
	void SetXmax();//the histogram will be calculated up to xmax in reduced units
	void AssignBinsize();//determines, how many different binsize has to be used, and assign it to the data sets
	void PrintRunParams(ExptsData &edata) const;//displays the parameters on screen
	//writing the status of the program at regular intervals
#ifdef _ADVANCED_GEOM_CONST
	void PrintStatus(ExptsData &edata, ChiSquared &chi, CoordNumbConst &icc,AvCoordConst &avcc, CommonNeighConst &conc, SecondNeighConst &snc, BondValenceSumConst &bvs, SimpleCfg &config, int status);
#else
	void PrintStatus(ExptsData &edata, ChiSquared &chi, CoordNumbConst &icc, AvCoordConst &avcc, SimpleCfg &config, int status);
#endif
	//saving the sigma values and chi2 in binary format, important for continuation of the run, if sigma was calculated
	void SaveState(ChiSquared &chi);
	//loading the sigma and chi2 values from binary format file in case of exact continuation of the run, important, if sigma was calculated
	static void LoadState();
	//get the data type of the seriesis
	static void GetDataType(int i, const char *data_type);
	static void CheckLeadSerInd();
	static void CheckType(int type, char *name, const char *routine_name);//Check whether the given type read from input is in the range
	static void NoLeadZeroSigma(int i, const char *name);//Error message, if the leading pot series have zero sigma
	static void GetPotParams(ifstream &file);//reading the potential related parameters from the *.dat file
	//checking the leading series index, as it could not be done at reading, as the bonded interactions were not available than
	
	static void CreateFreeFormat(const int write_reduced_free=0); //This method creates (0) .free complete, (1) .freer reduced and (2) .freers reduced free format files
																  //.freers only created, if scalable sigma parameters are used, and here the scaled sigmas are written. In case of the .freer the original scaling factors for the scalable sigma can be found
	
	static void OffsetPotLead();//offseting the lead_series_ind2 in case of free format reading
private:
	//In the following if onlydef is true, it writes only the non-default values out
	static void CreateFreeGeneral(std::ofstream &out, const bool onlyndef=false);  //This method writes general tags into free control files defined the stream "out"
	static void CreateFreeExp(const int expno, std::ofstream &out, const bool onlyndef=false);  //This method writes the expno'th experiment related tags into free control files defined the stream "out"
	static void CreateFreeCos(const int cosno, std::ofstream &out, const bool onlyndef = false);//This method writes the cosno'th cosine angle distribution constraint related tags into free control files defined the stream "out"
	static void CreateFreeCoord(const int cono, std::ofstream &out, const bool avcc=false, const bool onlyndef=false);//This method writes the cono'th (if average -- avcc=true) coordination constraint related tags into free control files defined the stream "out"
#ifdef _ADVANCED_GEOM_CONST
	static void CreateFreeConc(const int concno, std::ofstream &out, const bool onlyndef=false);//This method writes the concno'th common neighbourn constraint related tags into free control files defined the stream "out"
	static void CreateFreeSnc(const int sncno, std::ofstream &out, const bool onlyndef = false);//This method writes the sncno'th second neighbourn constraint related tags into free control files defined the stream "out"
	static void CreateFreeBvs(const int bvsno, std::ofstream &out, const bool onlyndef = false);//This method writes the bvsno'th bond valence sum constraint related tags into free control files defined the stream "out"
#endif
	static void CreateFreeNbpot(std::ofstream &out, const bool onlyndef=false);//This method writes non-bonded potential related tags into free control files defined the stream "out", if onlydef is true, it writes only the non-default values to out
	static void CreateFreeBpot(std::ofstream &out, const bool onlyndef=false);//This method writes bonded potential related tags (only those read normally from fixed *.dat) into free control files defined the stream "out", if onlydef is true, it writes only the non-default values to out
	static void CreateFreeSwap(std::ofstream &out);//This method writes swap related tags into free control files defined the stream "out"
	static void CreateFreeCustom(std::ofstream &out);//This method writes the custom moves option related tags into free control files defined the stream "out"
#ifdef _AENET
	static void CreateFreeAenet(std::ofstream &outconst, bool onlyndef);//This method writes aenet related tags into free control files defined the stream "out"
#endif
#ifdef _LOCAL_INV
	static void CreateFreeLoc(std::ofstream &out, const bool onlyndef = false);//This method writes the local invariance related tags into free control files defined the stream "out"
#endif
#ifdef _NO_PERIODIC
	static void CreateFreeNoper(std::ofstream &out, const bool onlyndef = false);//This method writes the no periodic boundary condition related tags into free control files defined the stream "out"
#endif
#ifdef _VIBR_AMP
	static void CreateFreeVibramp(std::ofstream &out);//This method writes vibrational amplitude related tags into free control files defined the stream "out"
#endif

	//These functions are processing the relevant initialization according to keywords provided in pool. They returns with false, if incoming datasets are insufficient.
	//Note and warn are counters required for statistics
	static bool ReadFreeGeneral(std::list<string> &pool);//general parameters
	static bool ReadFreeExp(std::list<std::list<string> > &pool);//experimental data
	static bool ReadFreeCos(std::list<std::list<string> > &pool);//cosine distribution of bond angle constraint
	static bool ReadFreeCoord(std::list<std::list<string> > &pool);//coordination constraint
	static bool ReadFreeAcoord(std::list<std::list<string> > &pool);//average coordination constraint
#ifdef _ADVANCED_GEOM_CONST
	static bool ReadFreeConc(std::list<std::list<string> > &pool);//common neighbour constraint
	static bool ReadFreeSnc(std::list<std::list<string> > &pool);//second neighbour constraint
	static bool ReadFreeBvs(std::list<std::list<string> > &pool);//bond valence sum constraint
#endif
	static bool ReadFreeCustom(std::list<string> &pool);//custom molecular move
	static bool ReadFreeSwap(std::list<string> &pool);//swap
	static bool ReadFreeNbpot(std::list<string> &pool);//non-bonded potential related (only parameters what would be in the fixed datfile
	static bool ReadFreeBpot(std::list<string> &pool);//bonded potential related (only parameters what would be in the fixed datfile
#ifdef _AENET
	static bool ReadFreeAenet(std::list<string> &pool);// AENET related
#endif
#ifdef _LOCAL_INV
	static  bool ReadFreeLoc(std::list<string> &pool);//local invariance
#endif
#ifdef _NO_PERIODIC
	static bool ReadFreeNoper(std::list<string> &pool);//no periodic boundary condition 
#endif
#ifdef _VIBR_AMP
	static  bool ReadFreeVibramp(std::list<string> &pool);//vibrational amplitude
#endif
	static bool CheckKeyInPool(const std::list<string> &pool, const string key) {return std::find_if(pool.begin(),pool.end(), [&key](string c) { return c.find(key)!=string::npos; }) != pool.end();};//Looking for a key string in pool std::list. Return true, if it contains one of it.
	static int RegularizeSigmaEntry(string &parameters);//This function regularizes the sigma entry as writes "-" signs instead of "s" in parameters and returns with the master index (start from 1) if it has been found and 0 if not.
	static bool HaveReference(std::list<string> *pointer, std::list<std::list<string> > &target);//Returns with true, if pointer points to an element of the target
	
	static bool CompareExp(std::list<string> &fst, std::list<string> &scnd);//Compare two objects according to gr,Sq,Fq,Ek order
	
};
	
//destructor of the RunParams class
inline RunParams::~RunParams()
{
	delete [] pcutsq;
	if(::debug)  cout<<"RunParams destructor"<<endl;

};

//==========================================================================
class FNC_POT
{
	//this class implements fixed neighbours constraints
	//if fnc option is 4, and code is compiled with _DEf_FLEX_MOL (flexible molecule definition)
	//than fnc is modified to be able to hold the angle and dihedral constraints
	//In this case 3 instance of FNC is created, one for the bonds, one for the angles and one for the dihedrals
	public:
	FNC_POT();//constructor, cannot figure out how to have member initializer list and create an array of the object in the same time
	

	~FNC_POT();//destructor
	static int nconstraints;//the number of FNC_POT constraints
	static int natoms;//number of FNC_POT central atoms
	static int nvirtuals;//number of virtual sites
	static int vdW_warning;//0 no warning, index of the partial (start with 1) where larger than POT_WARNING absolute value for the exclusion potential last occured
	static int Coul_warning;//0 no warning, index of the partial (start with 1) where larger than POT_WARNING absolute value for the exclusion potential last occured
	static int ndist;//the total number of constrained distances, define the size of the arrays, has to be static
	static int sum_out_of_range;//total number of pairs out of the FNC range (for option 3)
	static int *nout_of_range;//number of pairs/fnc type out of the FNC range (for option 3)
	static int *dec_out_of_range;//decrement the number of out of range atoms, if the move is accepted for each fnc type
	static int index_offset;//1 for BOND and 0 for normal fnc
	static int nexclusions;//total number of excluded atoms
	static int check_exclusions;//boolean whether to check the exclusions during potential correction due to the move of the charge group centres
	static int pot_dim;//dimension of the potential array
	static int nthreads;//number of threads
	static int nGRpartials;//number of GR partials for the NB parameter arrays
	static int	 nused_potpartials;//number of partials to use for tabulated potential (number of potential files)
	static int *npot_points;//array of number of data points in the tabulated potential files
	static int *tabpot_offset;//offset for each used partial for r_pot and u_pot
	static int *tabpot_index;//index of the partial to use the tabulated potential for, for each used partial (indexing start with 0)
	static std::vector<std::vector<int>> fnc_partial_ind;//in case of automatic cutoff and fnc it will indicate for each partial in which fnc constraint it is involved in 
	//these arrays will be static to be able to use them easier
	static double bond_tot_pot;//total bond potential
	static double angle_tot_pot;//total angle potential
	static double perdihedral_tot_pot;//total periodic dihedral potential
	static double harmdihedral_tot_pot;//total harmonic dihedral potential
	static double RBdihedral_tot_pot;//total RB dihedral potential
	static double vdW_tot_pot;//total vdW potential
	static double Coulomb_tot_pot;//total Coulomb potential
	static double vdW14_tot_pot;//total 1-4 vdW potential
	static double Coulomb14_tot_pot;//total 1-4 Coulomb potential
	static double f_Coulomb_A;//f_Coulomb/boxedge to transform the distance from reduced unit to A
	static double *bond_pot;//potential energy of the bonds 
	static double *angle_pot;//potential energy of the angles 
	static double *perdihedral_pot;//potential energy of the periodic dihedrals
	static double *harmdihedral_pot;//potential energy of the harm dihedrals 
	static double *RBdihedral_pot;//potential energy of the RB dihedrals 
	static double *vdW_pot;//non-bonded potential energy  
	static double *vdW_pot_s;//non-bonded potential energy for the samll values
	static double *vdW_pot_l;//non-bonded potential energy for the large values
	static double *vdW14_pot;//non-bonded potential energy  for 1-4 interaction
	static double *Coulomb_pot;//Coulomb potential energy 
	static double *Coulomb_pot_s;//Coulomb potential energy for the samll values
	static double *Coulomb_pot_l;//Coulomb potential energy for the large values 
	static double *Coulomb14_pot;//Coulomb potential energy for 1-4 interaction
	static double *bond_pot_old;//potential energy of the bonds for the previous move
	static double *angle_pot_old;//potential energy of the angles for the previous move
	static double *perdihedral_pot_old;//potential energy of the periodic dihedrals for the previous move
	static double *harmdihedral_pot_old;//potential energy of the harmonic dihedrals for the previous move
	static double *RBdihedral_pot_old;//potential energy of the RB dihedrals for the previous move
	static double *vdW_pot_old;//non-bonded potential energy for the previous move
	static double *Coulomb_pot_old;//Coulomb potential energy for the previous move
	static double *vdW14_pot_old;//non-bonded potential energy  for 1-4 interaction for the previous move
	static double *Coulomb14_pot_old;//Coulomb potential energy for 1-4 interaction for the previous move 
	static double *r_pot;//r values for the tabulated potential files
	static double *u_pot;//potential values for the tabulated potential files
	static double *tab_cutoff;//cutoff for each used partial tabulated potential
	static double *tabpot_dr;//reduced spacing of the tabulated potential for each used partial;
	static char(*tabpot_filename)[FILE_NAME_SIZE];//name of the file(s) containing the tabulated potential(s)
	static Threads *thread_object;

	int ndistances;//same as ndist, if normal FNC, differs for each instance, if molecular run
	
	//arrays of dimensions nconstraints:
	double *dmin;//array of minimal distances (1 per constraint) in Angstrom
	double *dmax;//array of maximum distances (1 per constraint) in Angstrom
	double *dminsq;//array of minimal distances squared in reduced units
	double *dmaxsq;//array of maximum distances squared in reduced units

	//arrays of dimensions natoms:
	int *numbneigh;//array of the number of neighbours of each atom
	 
	int *converter;//array of positions of atoms in the 
					//neighbours and consttypes array 	
	//arrays of dimensions ndistances		
	int *neighbours;//array of neighbour indices
	int *consttypes;//array of constraint types	

	T_interaction_type pot_type;//showing the interaction type
	int *neighbours2;//array of second neighbour indices (for angle and dihedrals)
	int *neighbours3;//array of third neighbour indices (for dihedrals)
	int *dih_type;//array of the types for the dihedrals
	int *vdW_excl_numb;//array for all the atoms with the number of exclusions of each atom for the non-bonded interaction calculation, atoms in RMC_order
	int *vdW_exclusions;//array for all the atoms with the indices of the excluded atoms for the non-bonded interaction calculation, atoms in RMC_order
	int *vdW_conv;//converter for the LJ_* arrays
	int *LJ_pow;//pointer to RunParams' power to the repulsion term of LJ interaction
	//the concept os the same as in GROMACS, the atoms are rendered to charge groups, and the coordinate for the centre of the charge group
	//is used deciding, whether all the atoms are inside or outside the cutoff 
	int *charge_type;//charge type of the atoms (the index of the charge in the Topology::charge array
	int *charge_group;//the charge group for each atom
	int *atom_for_charge_gr;//giving the atoms belonging to each charge group
	int *charge_centre;//the index of the charge centre in the SimpleCfg::charge_gr_centre array
	int **atom_for_charge_gr_finder;//the beginning of the atom indices belonging to the given charge group 
	double *LJ_C6;//LJ 4*epsilon*sigma^6 per partial
	double *LJ_CN;//LJ 4*epsilon*sigma^N per partial
	double *LJ14_C6;//LJ 1-4 4*epsilon*sigma^6 per partial
	double *LJ14_CN;//LJ 1-4 4*epsilon*sigma^N per partial
	double *bond_k;//original bonds_k*boxedge*boxedge to make calculation quicker
	double *bond_r;//reduced, original bonds_r/boxedge
	double **RBdih_C;//pointer to the C coefficients of the dihedral interaction 
	double *angle_ang_rad;//angles in radian for angles
	double *perdih_ang_rad;//angles in radian for the periodic dihedral
	double *harmdih_ang_rad;//angles in radian for the harmonic dihedral
	
	//pointer to object
	static SimpleCfg *config;
	static RunParams *rundat;

	//member functions
	static void GetFNCParams();//reads the FNC parameters from .fnc file
	static void ReadTabPot(ifstream &file);//reading the tabulated potential
	static void ReadTabPotFree();//reading the tabulated potential
	static void SaveTabulatedPotential(const char *file_name = "");//save the loaded tabulated potential
	void Save(const char *file_name="") const;//save to a file
	void RawSave(const char *file_name="") const;//produce a raw .fnc file
	void Load();//loads from a 'save' file
	void RawLoad();//loads from a .fnc file, and reduces the distances squared
	//Check, whether the configurations satisfy the FNC constraints, and whether the FNC 
	//distances are not below cutoff. Second parameter is the squared, reduced cutoff distances
	int CheckFNC(double *cutsq);
	//check, whether the FNC is satisfied after each move
	int CheckFNCChange(Move &movem);//gives back, whether it is acceptable
	void	Init(T_interaction_type type);
	void	ResizeNeighArrays();//resizing the neighbour, neighbour2,neighbour3 arrays
	void	GenerateExclusion( Topology &topology);//generate the exclusion list in case of nordonded interactions
	void	UpdatePot(int mode);//mode 0: Copying the values from the normal array to the old, 1: from the old to the normal
	void	SaveExclusions(const char *file_name="") const;//save the exclusions for each atom
	void	CalcNBExclusion(Move &move);//calculate the NB interaction for the exclusions to subtract them 
	void	CalcBonds();//calculate the potential for the bonds
	void	CalcAngles();//calculate the potential for the angles
	void	CalcDihedrals();//calculate the potential for the dihedrals
	void	CalcBondsChange(Move &move);//calculate the change in the  potential for the bonds
	void	CalcAnglesChange(Move &move);//calculate the change in the potential for the angles
	void	CalcDihedralsChange(Move &move);//calculate the change in the potential for the dihedrals
	void	CalcPeriodicDihedrals(int consttype, int sign,double theta);//calculate periodic dihedral potential
	void	CalcHarmonicDihedrals(int consttype,int sign, double theta);//calculate harmonic dihedral potential
	void	CalcRBDihedrals(int consttype, int sign,double theta);//calculate RB dihedral potential
	void	VectorProduct(double *a, double *b, double *result);//calculate vector product
	void	SavePotBinary() const;//saving the potentials (bond, angle, dihedral, vdW) in//	int		GetNExclusion(int iatom, int iprev, int iexcl);//get the number of exclusions for this atom 
	int		GetExclusions(int iatom, int iprev, int iexcl, int *array_start, int *excl_array, int checksize=0);//get the number of exclusions for this atom 
	int		LoadPotBinary() const;//loading the potentials (bond, angle, dihedral, vdW) from the binary format *.pot file
	double	CalcLJ(int ipartial, double dsquare);//calculate the LJ interaction
	double	CalcLJ14(int ipartial, double dsquare);//calculate the LJ interaction
	double  CalcTabPot(int ipartial, double dsquare);//calculate the tabulated potential
	double  CalcCoulomb(int ind1,int ind2, double dsquare, double scale=1.0);//calculate the Coulomb interaction
	double	CalcDihAngle(double *vecij, double *vecjk, double *veckl);//calculating the dihedral angle
	void    UpdateTotPot();//updating the total potential
	void	NBPotChargeGrCorr(ThreadArg &thread_arg, Move &move);//correcting the non-bonded potential for the atoms belonging to the charge group of the moved atoms
	static void SetPotParams();//setting the static potential arrays	
	static void InitPot();//setting the potential array elements to zero
	static int CheckPot();//check, whether the loaded and reacluated pot is the same
	static void ShowPotError(int i, char *potname, double recalc, double loaded);//error message during CheckPot
	typedef double (FNC_POT::*TCalcvdW)(int,double);	
	static TCalcvdW CalcvdW;//pointer to the actual vdW function
	static TCalcvdW CalcvdW14;//pointer to the actual 1-4 vdW function

};

#ifdef _DEF_INTERACTION_FUNC
	extern void (FNC_POT::*InteractionFunct[N_POT_TYPE])();
	extern void (FNC_POT::*InteractionFunctChange[N_POT_TYPE])(Move &move);
	extern void (FNC_POT::*DihedralFunction[N_DIH_FUNCT])(int consttype,int sign,double theta);
#else
	//will hold the pointer to the above declared functions
	void (FNC_POT::*InteractionFunct[N_POT_TYPE])()={&FNC_POT::CalcBonds,&FNC_POT::CalcAngles,&FNC_POT::CalcDihedrals};
	void (FNC_POT::*InteractionFunctChange[N_POT_TYPE])(Move &move)={&FNC_POT::CalcBondsChange,&FNC_POT::CalcAnglesChange,&FNC_POT::CalcDihedralsChange};
	void (FNC_POT::*DihedralFunction[N_DIH_FUNCT])(int consttype,int sign,double theta)={&FNC_POT::CalcPeriodicDihedrals,&FNC_POT::CalcHarmonicDihedrals,&FNC_POT::CalcRBDihedrals};
	#define _DEF_INTERACTION_FUNC
#endif


//destructor of the FNC_POT class
inline FNC_POT::~FNC_POT()
{
	delete [] dmin;
	delete [] dmax;
	delete [] dminsq;
	delete [] dmaxsq;
	delete [] numbneigh;
	delete [] converter;
	delete [] neighbours;
	delete [] consttypes;

	//to show,that it was not yet created
	if (neighbours2!=NULL)
		delete [] neighbours2;
	if (neighbours3!=NULL)
		delete [] neighbours3;
	if (dih_type!=NULL)
		delete [] dih_type;
	if (vdW_excl_numb!=NULL)
		delete [] vdW_excl_numb;
	if (vdW_exclusions!=NULL)
		delete [] vdW_exclusions;
	if (vdW_conv!=NULL)
		delete [] vdW_conv;
	if (LJ_pow!=NULL)
		delete [] LJ_pow;
	if (LJ_C6!=NULL)
		delete [] LJ_C6;
	if (LJ_CN!=NULL)
		delete [] LJ_CN;
	if (LJ14_C6!=NULL)
		delete [] LJ14_C6;
	if (LJ14_CN!=NULL)
		delete [] LJ14_CN;
	if (bond_r!=NULL)
		delete [] bond_r;
	if (bond_k!=NULL)
		delete [] bond_k;
	if (RBdih_C!=NULL)
		delete [] RBdih_C;
	if (angle_ang_rad!=NULL)
		delete [] angle_ang_rad;
	if (perdih_ang_rad!=NULL)
		delete [] perdih_ang_rad;
	if (harmdih_ang_rad!=NULL)
		delete [] harmdih_ang_rad;
	if (charge_type!=NULL)
		delete [] charge_type;
	if (charge_group!=NULL)
		delete [] charge_group;
	if (charge_centre!=NULL)
		delete [] charge_centre;
	if (atom_for_charge_gr!=NULL)
		delete [] atom_for_charge_gr;
	if (atom_for_charge_gr_finder!=NULL)
		delete [] atom_for_charge_gr_finder;
	if(::debug) cout<<"FNC_POT destructor"<<endl;
	
};
	

//========================================================================
//The concept is that only as many different bin sizes will be handled, as it is absolutely necessary. 
//As in case of g(r) data sets, the r point spacing and binshift is automatically detected, it can happen, that even in case of the 
//same binsize the binshift is different, so different histograms have to be used! In case of the same binsize but different binshifts
//it is examined, whether it is actually the same histogram with different starting bin. (For example binshift 0, 1, 2... are the same
//and will be considered as binshift 0, or 0.5, 1.5 and 2.5... will be considered as binshift 0.5. So only binshift<1 will occure for the different histograms! 
//So where binshifts with integer difference will share the same histogram, and only the first valid bin for the given data set is determined.
//To make thing easier, the histogram will be calculated from the first possible binshift (so if for example fom the 0, 1 and 2 series actually 1
//is the smallest used, it will be claculated for 0 as well.

class HistoSet
{
	public:
	//constructor
	//Bond Valence Sum is implemented as a handle
#ifdef _AENET
HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst,	AvCoordConst &acoordconst,
		 Aenet &aen, FNC_POT* &fnc_const, Threads &thread_obj);//default 
HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst,	AvCoordConst &acoordconst,
		 Aenet &aen, FNC_POT* &fnc_const, Threads &thread_obj, HistoSet &source); //Copy constructor
#else
HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst,\
		AvCoordConst &acoordconst, FNC_POT* &fnc_const, Threads &thread_obj);//default 
HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst,\
		AvCoordConst &acoordconst, FNC_POT* &fnc_const, Threads &thread_obj, HistoSet &source); //Copy constructor
#endif
	~HistoSet(); //destructor
	static longint *ngen;
	static int  ntypes;//number of atom types
	static int  nvirtualtypes;//number of virtual types
	static int  npartials;//number of partials
	static int  nthreads;//the total number of threads to use = RunParams::nthreads
	static int  ndiffbin;//number of different bin sizes
	static int  ntot_datasets;//total number of data series
	static int  ntot_bins;//total number of bins for all the different bin sizes
	static int *assign_hist;//pointer to ExptsData::assign_hist, to assigning the different histogram bin sizes to the data sets, original dim [ntot_datasets]
	static int *nbins;//number of bins (same for each histogram) for each different bin size
	static int *offset;//array elements for each partial and different bin size: total number of bins from the beginning of the given threads beginning in pcounts
	static int *calc_tabpot;//in case of potential=10 indicator for each partial, which partial is used dim [naprtials]
	static double boxedge;//that's the value of the boxedge unit in Angstrom
	static double *xmin;//minimum value for the abcissas (in boxedge units) for all the different unique bin sizes. dim [ndiffbin]
	static double *xmax;//maximum value for the abcissas (in boxedge units) for all the different unique bin sizes. dim [ndiffbin]
	static double *deltax;//bins width (in boxedge units), dim [ndiffbin]
	static double *inv_binwidth;//inverse of the bin width in reduced units or all the different unique bin sizes, dim [ndiffbin]

	static bool write_hist;//write the histogram, even if it was not calculated, the data set"s bin size assigment changed within acceptable limits

	//the following members differ from instance to instance
	int *tclarray;//pointer to the temporary array for counting the aroms below cutoff
#ifdef _TEST_MODE
	int maxtcp;//maximum number of too close pairs for a thread, needed for array dimensions
	int ntcp;//number of too close pairs
	int *tcpairs;//indices of too close pairs for the threads
#endif
		
	longint *pcounts;//pointer to the array of counts in each bin
	longint *ptotal;//pointer to the total number of counts
				//in each individual histogram
	longint **finder;//pointer to one particular histogram
	longint **thread_histcount_finder;//pointer to the beginning of the histogram of each thread 
#ifdef _VIBR_AMP
	longint *pcounts_T;//pointer to the array of counts in each bin
	longint *ptotal_T;//pointer to the total number of counts
				//in each individual histogram
	longint **finder_T;//pointer to one particular histogram
#endif
#ifdef _NO_PERIODIC
	double *periodic_pcounts;//histogram converted to periodic
	double **periodic_finder;//finder for periodic_counts for the partials
#endif
	//creating a pointer to the non-bonded interaction function, this pointer will be used in other part of the code
	//to make it possible to use different interactions without if structure in the loop, choice has to be made only once
	typedef double (FNC_POT::*TCalcvdW)(int,double);	
	static TCalcvdW CalcvdW;//pointer to the actual vdW function
	static TCalcvdW CalcvdW14;//pointer to the actual 1-4 vdW function

#ifdef _USE_LOCAL_INV 
	static int *loc_offset;//array elements for each partial: total number of local bins to the beginning of the given partial
	static int sum_bins;//max_loc_nbins*ntypes
	static int sum_bins_offset;//sum_bins+padding
	static int max_loc_nbins;//last bin compared to the normal histogram to use for local invariance calculation 
	static double min_loc_r;//minimum distance for the local invariance calculation in reduced unit
	static double max_loc_r;//maximum distance for the local invariance calculation in reduced unit
	static double inv_locbinwidth;//inverse of the local invariance bin width in reduced units

	int *local_count;//local histogram counts for each atom
	longint *local_sum_count;//sum of local histogram counts, X-Y and Y-X partials together, as in case of normal histogram 
	longint **local_sum_count_finder;//finder to local_sum_count

	longint **thread_lochistcount_finder;//pointer to the beginning of the histogram of each thread 
	//it would require a lot of memory to have a local count segment for each thread. As only for the moved atoms' local
	//histogram is used by all the threads (as the neighbours are divided among the treads) a delta local count array is written for the 
	//moved atoms, and it added to the local histogram of these atoms, after all the threads have finished
	int *dlocal_count;//the local count change in one step for all the moved atoms
	int **dlocal_count_finder;//finder for the beginning of the threads segment
	double delta_chi2_inv;//contribution to chi2 from local invariance terms
#endif

#ifdef _NO_PERIODIC
	static double R0_red;
	static int *atom_binind;//the index of the bin the atom can be found in
#endif

#ifdef _VIBR_AMP
	static  int  nbins_conv;//the number of bins for the Gaussian convolution funtion
	static	int	 nbins_used;//the number of bins where the ppcf will be good, nbins-3*nbins_conv/5
	static  int *temp_conv_hist;//the convoluted histogram counts for one histogram bin
	static  double *thermal_conv;//convolution function for thermal vibration correction 
	static  double **thermal_conv_finder;//finder for each partial in the thermal conv array
	
	static  double *T_corr_sigma;//sigma parameter for the Gaussian convolution function for each partial
#endif
	//Handles to objects
	SimpleCfg &conf;
	RunParams &rundat;
	CoordNumbConst &coordnc;
	AvCoordConst &acoordc;
#ifdef _ADVANCED_GEOM_CONST
	BondValenceSumConst *bvs;//pointer, it would have been troublesome to put in the member initializer list
#endif
#ifdef _AENET
	Aenet &aenet;
#endif
	FNC_POT* &fnc;
	Threads &thread_object;


	static Move *move;//this cannot be initialized with the constructor, has to be assigned later

	//member functions
	static void SetHistParams(RunParams &rundata);//Sets the static members
	void InitHist();//initializing the histogram to zero
	//calculating the whole histograms
	void HistCalc(int &ntooclose);
	//calculating the change in the histograms due to the move
	void HistCalcChange();
	void Save(const char *file_name="") const;//save to a file
#ifdef _USE_LOCAL_INV
	void SaveLoc(const char *file_name="") const;//save to a file
	void LoadLocCount(ifstream &file);//loading the local histogram counts
	void InitLocHist();//initializing the local histogram to zero
#endif
#ifdef _NO_PERIODIC
	int LoadPosBin();//loading the central bin index of the atoms
	void SavePosBin();//saving the central bin index of the atoms
	void CalcPerHist(ThreadArg &thread_arg, DataMat &datmat);//calculating the histogram, which would correspond to periodic boundary conditions
#endif
	int Load(int mode=0);//loads the init part from a file
	int LoadCount(ifstream &file,int mode, int binsize_index=0);//loading the histogram counts
	void CopyModified(HistoSet &target, ThreadArg &thread_arg);//Copy the modified parts of the histogram to target
	void HistCalcThread(ThreadArg &thread_arg);//the calculation of the histogram performed by the threads
	void HistCalcChangeThread(ThreadArg &thread_arg);//the calculation of the change in the histogram performed by the threads
#ifdef _TEST_MODE
	void SaveTooClosePairs(char *file_name);//saving the too close pairs
#endif
#ifdef _VIBR_AMP
	static void CalcConvFunction();//calculating the convolution function for 
	void ThermalCorr();//calculating the convoluted histogram
#endif
};	

//destructor of the Histoset class
inline HistoSet::~HistoSet()
{
	delete [] pcounts;
	delete [] finder;
	if (thread_histcount_finder!=NULL)
		delete [] thread_histcount_finder;
	delete [] ptotal;
#ifdef _USE_LOCAL_INV
	delete [] local_count;
	delete [] local_sum_count;
	delete [] local_sum_count_finder;
	if (nthreads>1 && dlocal_count!=NULL)
		delete [] dlocal_count; 
	if (nthreads>1 && dlocal_count_finder!=NULL)
		delete [] dlocal_count_finder; 
	if (thread_lochistcount_finder!=NULL)
		delete [] thread_lochistcount_finder;

#endif
#ifdef _NO_PERIODIC
	delete [] periodic_pcounts;
	delete [] periodic_finder;
#endif
#ifdef _VIBR_AMP
	delete [] pcounts_T;
	delete [] ptotal_T;
	delete [] finder_T;
#endif
	if(::debug) cout<<"HistoSet destructor"<<endl;
};
