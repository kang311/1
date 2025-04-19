//header for 'secondary' classes
//Last changed 12.12.2022

#include "NeighbourList.h"//classes1.h" //which includes utilities

#ifdef _ATLAS //ATLAS is used
	#include "rmcblas.h"
#endif
#ifdef _LOCAL_INV
	struct ThreadArg;
#endif
//========================================================================
//-------------class DataMat--------------
//this class implements the tables needed to convert histograms
//to partial PCFs  and PCFs to structure factors
class DataMat
	{//this class implements the matrices that will convert the 
	//histogram counts to partial PPCF and SF
	public:
	DataMat(ExptsData &edata);//constructor
		
	~DataMat();//destructor
	static int  ngr;//number of g(r) data sets
	static int  nsq;//number of S(Q) data sets
	static int  nfq;//number of F(Q) data set
	static int  nfg;//number of F(g) data set
	static int  npartials;//number of partials
	static int  ntypes;//number of atom types
	static int  calcmode;//whether to calculate from the local hist, default is not
	static int  ndiffbin;//number of different bin sizes
	static int  sfactor_size;//number of points in the surface factor
	static int *nbins;//number of bins in the g(r) histograms
	static int *mybins;//the number of used bins, can be unique, but can be smaller than nbins
	static int *assign_hist;//assigning the different histogram bin sizes to the data sets dim [ntot_datasets]
	
#ifdef _VIBR_AMP
	static int nbins_used;//number of bins used in the g(r) histograms
#endif
	static int *pnatoms;//pointer to number of atoms per type
	static double *xmin;//pointer to RunParams::xmin, minimum value for the histogram in reduced units, as many as the number of different unique binsizes, dim [ndiffbin]
	static double *xmax;//pointer to RunParams::xmax, the maximum value for the histogram in reduced units, as many as the number of different unique binsizes, dim [ndiffbin]
	static double boxedge;//size of the simulation cell in Angstrom
	static double dr;//r spacing in Angstrom
	static double volume;//sample volume
	static double rho;//total number density
#ifdef _LOCAL_INV
	static double *dV_loc;//the histogram bin volumes has to be stored in case of locale invariance calculation
#endif
#ifdef _NO_PERIODIC
	static double R0;//the radius of the speric sample inside the box
	static double *dV_in;//the actual volume of the ribbon inside sphere A for all the possible origins of sphere B and C
	static double *dV;//the volume of the normal volume element, needed if exafs or vibr_amp is used for non-periodic
#endif
	int *grsize;//pointer to the array of g(r) data sets sizes
	int *sqsize;//pointer to the array of  S(Q) data sets sizes
					//i.e. number of points
	int *fqsize;//pointer to the array of  F(Q) data sets sizes
					//i.e. number of points	
	int *fgsize;//pointer to the array of  F(g) data sets sizes
					//i.e. number of points	
		
	//tables
	double *ntable;//table of normalisation factors to obtain g(r) values
					//from histograms counts
	double *nfactor;//for calculating the factor to normalize the histogram beside the volume element
	//matrices
	double *grmij;//array of matrices elements for g(r) data;
					//i.e. (r bins->  data points converter)
	// (sine Fourier matrices);
	double *sqmij;//array of matrices elements for S(Q) data
	double *fqmij;//array of matrices elements for F(Q) data;
	double *fgmij;//array of matrices elements for F(g) data
	
	double **grfinder;//pointer to the start of individual matrices
					//in the grmij array for g(r) data sets
	double **sqfinder;//pointer to the start of individual matrices
					//in the sqmij array for S(Q) data sets
	double **fqfinder;//pointer to the start of individual matrices
					//in the fqmij array for F(Q) data sets
	double **fgfinder;//pointer to the start of individual matrices
					//in the fgmij array for F(g) data sets
	double **ntable_finder;//finders for  each partial ntable
					
	static void SetDataMatParams(RunParams &rundat);//sets the static members
	void Save() const;//saves the contents on a file
	void CalcdV(int mode, double *rtable, double *sftable);//calculating the renormalization table in case of mode=0, or storing the dV for mode=1 
};	

//destructor of the DataMat class
inline DataMat::~DataMat()
{
	delete [] grsize;
	delete [] sqsize;
	delete [] fqsize;
	delete [] fgsize;
	delete [] ntable;
	delete [] nfactor;
	delete [] grmij;
	delete [] sqmij;
	delete [] fqmij;
	delete [] fgmij;
	delete [] grfinder;
	delete [] sqfinder;
	delete [] fqfinder;
	delete [] fgfinder;
	delete [] ntable_finder;
	if(::debug)  cout<<"DataMat destructor"<<endl;
};


//========================================================================
class PPCFSet
{
	public:
	PPCFSet(HistoSet &histoset, DataMat &datamat,Move &move_obj); //default constructor
	PPCFSet(PPCFSet &source,HistoSet &histoset, DataMat &datamat,Move &move_obj);//copy constructor
		
	~PPCFSet(); //destructor
	static int  ntypes;//number of atom types
	static int  npartials;//number of partials
	static int  ndiffbin;//number of different bin sizes
	static int *nbins;//number of points (equals number of hist. bins) hor each bin size
	static int *mybins;//the number of used bins for each bin size
#ifdef _VIBR_AMP
	static	int	 nbins_used;//the number of bins where the ppcf will be good, nbins-3*nbins_conv/5
#endif
#ifdef _NO_PERIODIC
	static int calcmode;//whether to calculate from the local hist
#endif
	static int *offset;//array elements for each partial: total number of bins
						//to the beginning of the given partial
	
	static double *deltar;//pointer to RunParams__rspacing_diff, the discretization step (in Angstroms) for each different bin size
	static double *rmin;//array of minimum value for the abcissas (in Angstroms) for each bin size dim [ndiffbin]
	static double *rmax;//array of maximum value for the abcissas (in Angstroms) for each bin size dim [ndiffbin]

	int ncollect;//the number of times, the gvalues were summed
	double *gvalues;//pointer to the array of g(r) values
	double *gvalues_sum;//sum of the gvalues, to calculate the average
	double **finder;//pointer to one particular partial
	double **finder_sum;//pointer to one particular partial in gvalues_sum

#ifdef _ATLAS
	static double *gr_0;//ppcf-1, needed for the ATLAS in case of S(Q) and F(Q) and F(g) sets
#endif 

	//handles to objects
	HistoSet &hist;
	DataMat &datmat;
	Move &move;

	//member functions
	static void SetPPCFParams(RunParams &rundata);//Sets the static members
	void CalcPPCF(ThreadArg &thread_arg);//calculating the PPCF from the histogram
	void CalcModPPCF(ThreadArg &thread_arg);//calculating the PPCF for the modified types from the histogram
	void CopyModified(PPCFSet &target, ThreadArg &thread_arg);//Copy the modified parts of the ppcf to target	
	void Save(const char *file_name="");//save to a file
	void SaveOldOut(ofstream &file, const char *file_name="") const;//save in RMCA format to *.out file
	void AddPPCF();//Add the present ppcf-s to the sum

};	

//destructor of the PPCFSet class
inline PPCFSet::~PPCFSet()
{
	delete [] gvalues;
	if (gvalues_sum!=NULL)
	{
		delete [] gvalues_sum;
		delete [] finder_sum;
	}
	delete [] finder;
	if(::debug)  cout<<"PPCFSet destructor"<<endl;
};

//========================================================================
class CalcPart
{
	//this class is used to store the calculated partials necessary for
	//comparison with experimental data
	//this is used as an algorithmic tool and in principle it is not
	//intended to stand as its own therefore all the structure (size)
	//parameters are not stored but are those of the companion ExptsData
	//object
	public:
	CalcPart(ExptsData &edat, DataMat &datamat, PPCFSet &ppcfset, HistoSet &histog, Move &move_obj);//since the size of the CalcPart object
	//derives mainly from the different data sets used (and the number
	//of partials -i.e. of atom types-) the constructor takes an ExptsData
	//object in argument, the others are initializing the handles
	CalcPart(CalcPart &source, ExptsData &edat, DataMat &datamat, PPCFSet &ppcfset,  HistoSet &histog, Move &move_obj);//Copy constructor
	
	~CalcPart(); //destructor

	static int  ntypes;//number of types
	static int  npartials;//number of partials
	static int *nbins;//pointer to RunParams::nbins number of bins
#ifdef _VIBR_AMP
	static	int	 nbins_used;//the number of bins where the ppcf will be good, nbins-3*nbins_conv/5
#endif

	static int ngr;//number of g(r) data sets
	static int nsq;//number of S(Q) data sets
	static int nfq;//number of F(Q) data set
	static int nfg;//number of F(g) data set
	static int nek;//number of F(Q) data set
	//number of used data points /data set
	static int *grused;//array of number of g(r) points/set
	static int *sqused;//array of number of S(Q) points/set
	static int *fqused;//array of number of F(Q) points/set
	static int *fgused;//array of number of F(g) points/set
	static int *ekused;//array of number of E(k) points/set
	
	//data arrays dimensions [n_sets x npartials x n_data_points]
	double *grvalues;//g(r) values 
	double *sqvalues;//S(Q) values
	double *fqvalues;//F(Q) values 
	double *fgvalues;//F(g) values 
	double *ekvalues;//E(k) values 
	
	//finder pointers
	double **grfinder;//finder of g(r) partials
	double **sqfinder;//finder of S(Q) partials
	double **fqfinder;//finder of F(Q) partials
	double **fgfinder;//finder of F(g) partials
	double **ekfinder;//finder of E(k) partials
	
	//handles to objects
	ExptsData &edata;
	DataMat &datmat;
	PPCFSet &ppcf;
	HistoSet &hist;//for EXAFS
	Move &move;

	//member functions
	static void SetCalcPartParams();//sets the static parameters
	//calculate the partial pair correlation functions for the same data points, as the experimental
	void CalcPartialgr(ThreadArg &thread_arg);
	//calculate the partial F(Q) functions
	void CalcPartialFq(ThreadArg &thread_arg);
	//calculate the partial F(g) functions
	void CalcPartialFg(ThreadArg &thread_arg);
	//calculate the partial S(Q) functions
	void CalcPartialSq(ThreadArg &thread_arg);
	//calculate the partial E(k) functions
	void CalcPartialEk(ThreadArg &thread_arg);
	//calculate the modified part of the partial pair correlation functions for the same data points, as the experimental
	//in case of g(r) fitting, the partial S(Q), F(Q) and E(k) in case of S(Q), F(Q) and E(k) fitting
	void CalcModPartial(ThreadArg &thread_arg);
	void Savegr() const;//saves the partial g(r)-s to a file
	void SaveSQ() const;//saves the partial S(Q)-s to a file
	void SaveFQ() const;//saves the partial F(Q)-s to a file
	void SaveFg() const;//saves the partial F(g)-s to a file
	void SaveEK() const;//saves the partial E(k)-s to a file
	void SaveOldOut(ofstream &file, const char *file_name="") const;//save in RMCA format to *.out file
	void CopyModified(CalcPart &target,ThreadArg &thread_arg);//Copy the modified parts of the partials to target
};

//destructor of the CalcPart class
inline CalcPart::~CalcPart()
{
	delete [] grvalues;
	delete [] sqvalues;
	delete [] fqvalues;
	delete [] fgvalues;
	delete [] ekvalues;
	delete [] grfinder;
	delete [] sqfinder;
	delete [] fqfinder;
	delete [] fgfinder;
	delete [] ekfinder;
	if(::debug)  cout<<"CalcPart destructor"<<endl;
};
	
//=======================================================================
class CalcData
{
	//this class describes calculated data, mimicking the real
	//data given in an ExptsData object. Thhis class
	//cannot be used  without the source ExptsData object
	public:
	CalcData(ExptsData &edat, CalcPart &calcpartial, ChiSquared &chisq); //constructor
	CalcData(CalcData &source, ExptsData &edat, CalcPart &calcpartial, ChiSquared &chisq); //constructor
	
	~CalcData(); //destructor

	static int ntypes;//number of types
	static int npartials;//number of partials
	static int ngr;//number of g(r) data sets
	static int nsq;//number of S(Q) data sets
	static int nfq;//number of F(Q) data set
	static int nfg;//number of F(g) data set
	static int nek;//number of F(Q) data set

	//data arrays
	double *grvalues;//g(r) values
	double *sqvalues;//S(Q) values
	double *fqvalues;//F(Q) values
	double *fgvalues;//F(g) values
	double *ekvalues;//E(k) values
	
	//in case of xray I(Q) fitting
	double *iqvalues;//I(Q) values 
	double *compton;//B(Q)*exp(-alpha*Q*Q) 
		
	//finder pointers
	double **grfinder;//finder of g(r) data
	double **sqfinder;//finder of S(Q) data
	double **fqfinder;//finder of F(Q) data
	double **fgfinder;//finder of F(g) data
	double **ekfinder;//finder of E(k) data
	
	//in case of xray I(Q) fitting
	double **iqfinder;//finder of I(Q) data 
	double **comptfinder;//finder of compton correction

	//handles to objects
	ExptsData &edata;
	CalcPart &calcp;
	ChiSquared &chi;
			
	//member functions
	static void SetCalcDataParams();//sets the static parameters
	void CalcTotal(ThreadArg &thread_arg);
	void RecalcIQ(int i);//recalc I(Q) data sets during non-linear regression iteration steps
	void Copy(CalcData &target);//copying the totals, needed only for TESt MODE
	void Save() const;//save to a file
	void SaveResult(const char *file_name="") const;//saves the calculated and experimental data in a file
	void SaveOldOut(ofstream &file, const char *file_name="") const;//save in RMCA format to *.out file
};	

//destructor of the CalcData class
inline CalcData::~CalcData()
{
	delete [] grvalues;
	delete [] sqvalues;
	delete [] fqvalues;
	delete [] fgvalues;
	delete [] ekvalues;
	delete [] grfinder;
	delete [] sqfinder;
	delete [] fqfinder;
	delete [] fgfinder;
	delete [] ekfinder;
	if (iqvalues != NULL)
		delete[] iqvalues;
	if (iqfinder != NULL)
		delete[] iqfinder;
	if (compton != NULL)
		delete[] compton;
	if (comptfinder != NULL)
		delete[] comptfinder;
	if(::debug)  cout<<"CalcData destructor"<<endl;
};

//========================================================================
class ChiSquared
{
	public:
	#ifdef _AENET
		ChiSquared(ExptsData &edat,CosDistrConst &cosc, HistoSet &histog, Threads &thread_object, Aenet &aen); //constructor
	#else
		ChiSquared(ExptsData &edat,CosDistrConst &cosc, HistoSet &histog, Threads &thread_object); //constructor
	#endif

	~ChiSquared(); //destructor
	static int chisize;//number of elements in the chicomp array
	static int ngr;//number of g(r) data sets
	static int nsq;//number of S(Q) data sets
	static int nfq;//number of F(Q) data sets
	static int nfg;//number of F(g) data sets
	static int nek;//number of EXAFS data sets
	static int nacc;//number of average coordination constraints
	static int nicc;//number of individual coordination constraints
	static int ncos;//number of individual cosine distribution of bond angle constraints
	static int ncommonneigh;//number of individual common neighbour constraints
	static int nsecondneigh;//number of individual second neighbour constraints
	static int nbvs;//total number of bond valence sum constraints
	static int tot_cc_subconst;//total number of coordination subconstraints
	static int calc_sigma;//indicator, whether the sigma should be calculated
	static int nNB;//number of non-bonded interactions with different sigma
	static int npair_types;//number of 1-4 non-bonded interaction types
	static int nbonds;//number of bond types to use in chisquare calculation (if Topology::bond_weight_mode==0 ->1, or nbond_types)
	static int nangles;//number of angle types to use in chisquare calculation(if Topology::angle_weight_mode==0 ->1, or nangle_types)
	static int nperdihs;//number of periodic dihedral types to use in chisquare calculation(if Topology::perdih_weight_mode==0 ->1, or nperdihedral_types)
	static int nharmdihs;//number of harmonic dihedral types to use in chisquare calculation(if Topology::harmdih_weight_mode==0 ->1, or nharmdihedral_types)
	static int nRBdihs;//number of RB dihedral types to use in chisquare calculation(if Topology::RBdih_weight_mode==0 ->1, or nRBdihedral_types)
	static int chi_potstart;//beginning of the first potential related member in the chicomp array
	static int failed_nlr;//number of failed nonlinear regression since the last successful 
	static int *tot_failed_nlr;//total number for failed non-lin regressions for an I(Q) set
	static double *p_pot_chitotal;//pointer to the total chi2 value to be used for the potential related components
	static double total2;//sum of the chi squared potential components if non-bonded/aenet potential option is on, all the bonded and non-bonded
					//potential conribution is added to this. If there is only bonded potential, than its contribution is added to total
	int lead_ser_pot;//the index of the leading series for the bonded interactions  
	int *fit_index;//for each data set, defined by  int *(renorm+i) + *(offset+i) * 2 + *(linear+i) * 4 + *(quadratic+i) * 8 + *(cubic+i) * 16 for non-xrayray sets
		//for x-ray the above + ExptsData::fqfitIQ[i] * 32. This gives a unique index for the renormalization factor combination
	double lead_chi2_pot;//the value of the leading potential-based interactions ' chi2
	double lead_sigma_pot;//the sigma value for the leading potential-based interaction
	static double pot_chi2_low_limit;//if the potential does not have to be
	static double total;//sum of the chi squared components
	double lead_chi2;//the value of the leading series' chi2
	double lead_sigma;//the sigma value for the leading series
	double *chicomp;//array of chi squares components
	double *sigma_percentage;//the starting sigma values

					//one for each data set
	//now the sums each array has a number of elements equal to
	//the number of data sets, X can be r or Q depending on the data set
	double *f;//sums of experimental values
	double *ff;//sums os of experimental values squared
	double *fx;//sums of experimental values multiplyed by X values
	double *fx2;//sums of experimental values multiplyed by X square values
	double *fx3;//sums of experimental values multiplyed by X ^3 values
	double *s;//sums of calculated values
	double *ss;//sums of calculated values squared
	double *sx;//sums of calculated values multiplyed by X values
	double *sx2;//sums of calculated values multiplyed by X square values
	double *sx3;//sums of calculated values multiplyed by X ^3 values
	double *fs;//sums of cross products
	double *x;//sums  X values of data
	double *x2;//sums X square values of data
	double *x3;//sums X cubic values of data
	double *x4;//sums X ^4 values of data
	double *x5;//sums X ^5 values of data
	double *x6;//sums X ^6 values of data
	double *a;//multiplicative correction factor
	double *b;//additive correction factor
	double *c;//linear correction factor
	double *d;//quadratic correction factor
	double *e;//cubic correction factor
	double *alpha;//for I(Q)-fit correction factor
	double *Rw;//measurement of the goodness of the fit
	//for the nonlinear regression in case of I(Q) fitting
	double *C;//sum of B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j)
	double *CC;//sum of (B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j))^2
	double *fC;//sum of (experimental*( B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j))
	double *sC;// sum of(I(Q)_calc * (B_j * Q_j * Q_j * exp(-alpha * Q_j * Q_j))
	double FD1;//first derivate according to a (beta1 in numerical recepies https://e-maxx.ru/bookz/files/numerical_recipes.pdf 15.5.1)
	double FD2;//first derivate according to alpha, beta2
	double FD3;//first derivate according to b, beta3
	double SD11;//second derivate according to a,a (this is alpha in numerical recepies
	double SD12;//second derivate according to a,alpha
	double SD13;//second derivate according to a,b
	double SD22;//second derivate according to alpha,alpha
	double SD23;//second derivate according to alpha,b
	double SD33;//second derivate according to b,b
	double *last_nlr_chi2fr;//for each xray data set, in case of non-linear regression, for the last accepted move chi2_end/chi2_start for the regression
	double *maxdec_nlr_chi2fr;//for each xray data set, in case of non-linear regression,maximum decrease in chi2, chi2_end/chi2_start is the smallest for the regression

#ifdef _LOCAL_INV
	static int nlocint;//number of local inv intervals
	static int sum_loc_natoms_tot;//sum of loc_natoms, needed for some arrays
	static int *sum_loc_natoms;//sum of loc_natoms, needed for some arrays for each interval
	static int *cumalpart_loc_atoms;//cumulative number of loc atoms for each partials, where 1-2<> 2-1 (ntypes^2 partials)
	static int *loc_natoms_min;//the first atom to be used from each type in case of loc_chi2_mode=1 for each loc int (first all int of type 1, then type2...)
	static int *loc_natoms;//X atom to be used from each type in case of loc_chi2_mode=1
	static int *nused_loc;//number of used threads for local inv. per type
	static int *first_loc_ind;//the first histogram bin to use for a central atom/per central atom for each thread
	static int *starting_count;//the starting count for the first used bin/per central atom for each thread
	static int *av_loc_hist_offset;//offset to the begininnig of each type in av_loc_hist for a thread from the beginning of this thread's segment
	static double *av_loc_hist;//the averaged histogram bins according to the distance of the i-th atom of a type
	double *loc_chi;//values for the local invariance chi2 for each thread, padding will be used
	double **loc_chi_finder;//finder for the local invariance chi2 for each thread, needed because of the cache padding
	double loc_bin_width;//the width of the local invariance histogram
	double *min_loc_r_now;//minimum average distance for the local invariance calculation at present for each type
	double *max_loc_r_now;//maximum average distance for the local invariance calculation at present for each type
#endif

	char *name,*name_ext=NULL;//for display in CalcChiSuared, it is declated here and not as a local variable to avoid
						//	creating and deleteing at each subroutine call


	//handle
	ExptsData &edata;
	CosDistrConst &cosconst;
	HistoSet &hist;//only needed in case of local invariance
	Threads &thread_obj;//only needed in case of local invariance
#ifdef _AENET
	Aenet &aenet;
#endif
	
	//member functions
	static void SetChiSquaredParams(RunParams &rundat);//sets the static members
	//Calculating the chisquared
#ifdef _ADVANCED_GEOM_CONST
	longint CalcChiSquared(CalcData &calc, CoordNumbConst &icc, AvCoordConst &avcc, CommonNeighConst &conc, SecondNeighConst &snc, BondValenceSumConst &bvs);
#else
	longint CalcChiSquared(CalcData &calc, CoordNumbConst &icc, AvCoordConst &avcc);
#endif
	void CalcPotChiSquared();//calculating the chi2 for the potential related components
#ifdef _LOCAL_INV
	void CalcLocChiSquaredThread_bin(ThreadArg &thread_arg);//calculating the chisquared contribution based on bins by the threads 
	void CalcLocAvThread_dist(ThreadArg &thread_arg);//calculating the distance_based averega histogram, multi-thraeding cannot be introduced sufficiently
	void CalcLocChiSquaredThread_dist(ThreadArg &thread_arg);//calculating the chisquared contribution based on distance by the threads 
	typedef void (ChiSquared::*TCalcLocChi2Thread)(ThreadArg &);	
	static TCalcLocChi2Thread CalcLocChi2Thread;//pointer to the actual local invariane chi2 function

#endif
	void Save(ofstream &file) const;//save to a file
	void Coefficient(int nset, int first_set, int *renorm, int *offset, int *linear, int *quadratic, int *cubic, int *nused, bool checkIQ);
	void CalcSigma();//Calculate the sigma values, if they were given as a percentage 
	void ZeroChi2Error(const char *name, int mode=0); //error message
	void MissingRenorCombError(int i);//in case of Xray I(Q) fitting only some of the possible combinations is meaningful and therefore implemented, otherwise this message is given

	private:
	void five(double A5, double B5, double C5, double D5, double E5, double F5, double G5, double H5,\
					  double I5, double J5, double K5, double L5, double M5, double N5, double O5, double P5,\
					  double Q5, double R5, double S5, double T5, double &out1, double &out2, double &out3, \
					  double &out4, double &out5);
	void four(double j, double k, double l, double m, double n, double o, double p, double q, double r, double s, double t,\
			  double u, double v, double w, double &out1, double &out2, double &out3, double &out4);
	void three(double A, double B, double C, double D, double E, double F,\
	 			double G, double H, double I, double &out1, double &out2, double &out3);
	void two(double A, double B, double C, double D, double E, double &out1, double &out2);
	int NonLinReg(CalcData &calc);//Finding parameters with non-linear regression in case of I(Q) is fitting, gives back 1 if any of the I(Q) sets did not converge
	void CalcDerivates(CalcData &calc,int i);//calculate the sums for chi2 and the derivates for non-lin regression
};

//destructor of the ChiSquared class
inline ChiSquared::~ChiSquared()
{
	delete [] chicomp;
	delete [] f;
	delete [] ff;
	delete [] fx;
	delete [] fx2;
	delete [] fx3;
	delete [] s;
	delete [] ss;
	delete [] sx;
	delete [] sx2;
	delete [] sx3;
	delete [] fs;
	delete [] x;
	delete [] x2;
	delete [] x3;
	delete [] x4;
	delete [] x5;
	delete [] x6;
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] d;
	delete [] e;
	delete [] Rw;
	delete[] fit_index;
	if (alpha != NULL)
		delete[] alpha;
	if (C!= NULL)
		delete[] C;
	if (CC != NULL)
		delete[] CC;
	if (sC != NULL)
		delete[] sC;
	if (fC != NULL)
		delete[] fC;
	
	if (maxdec_nlr_chi2fr != NULL)
		delete[] maxdec_nlr_chi2fr;
	if (last_nlr_chi2fr != NULL)
		delete[] last_nlr_chi2fr;
		
	delete [] name;
	if (name_ext != NULL)
		delete [] name_ext;
#ifdef _LOCAl_INV
	delete [] loc_chi;
	delete [] loc_chi_finder;
#endif
	
	if(::debug)  cout<<"ChiSquared destructor"<<endl;
};
//=====================History================
class History
{
	public:
	History();//constructor
	//the arguments are 
	//the chi	squared
	//the buffers size
	//and the steps ratio (with respect to the screen updates)
	
	~History();//destructor
	static int configtype;//1, 3 if the coordinates come from the binary .bcf file, 0 if from the text .cfg file 
	static int swap;//boolean, whether there are swaps
	static int buffsize;//buffer size
	static int chisize;//number of chi squared parameters to buffer
	static int stepratio;//how often should values be kept
				//(with respect to saving)
	static int fill;//current number of buffered steps
	static int datapoint_tot;//total number of data points
	static int chisize2;//number of standard, (not potential related) chi squared parameters to save
	static int datapoint_tot2;//total number of data points in case of non-bonded pot present for the bonded and non-bonded interacions
	static int nNB;//number of non-bonded interactions with different sigma
	static int chi_offset1;//first element of 1-4 poptential params in chicomp array
	static int chi_offset2;//first element of 1-4 poptential params in chicomp array
	static longint *generated;//number of generated moves
	static longint *swapgen;//number of generated swaps
	static longint *E0shiftgen;//number of generated E0 shifts
	static longint *fprimeshiftgen;//number of generated f' shifts
	static longint *mucorrgen;//number of generated I(Q) mu corrections
	static longint *tried;//number of tried moves
	static longint *goodnonlinreg;//number of moves ended in good nonlin reg
	static longint *accepted;//number of accepted moves
	static longint *potaccepted;//number of accepted moves
	static longint *swapacc;//number of accepted swaps
	static longint *E0shiftacc;//number of accepted E0 shifts
	static longint *fprimeshiftacc;//number of accepted f' shifts
	static longint *mucorracc;//number of accepted I(Q) mu corrections
	static longint last_gen;//number of generated moves since the last buffering, cannot use rundat's as displaying and history buffering happens at differnt times
	static longint last_tried;//number of tried moves since the last buffering
	static longint last_acc;//number of accepted moves since the last buffering
	static longint calc_method;//indicating the calculation method of chi2
	static double *time;//time values
	static double *ir1;//instant ratio of number of tried/generated moves
				//over the last nstep moves (defined by the run parameters)
	static double *ir2;//instant ratio of number of accepted/tried moves
				//over the last nstep moves (defined by the run parameters)
	static double *chicomp;//chi squared components (one per data set + 1 per CC+
	
#ifdef _LOCAL_INV
	static int nlocint;//number of local invariance intervals
	static int *loc_npoints;//number of data points for local invariance (depend on the calculation mode)
#endif
		
	//member functions
	static void SetHistoryParams(RunParams &rundat);//sets the static members
	void Inifile(ofstream &file, RunParams &rundat, ExptsData &edata);//writes the beginning of the .hst file
	void Save(ofstream &file);//save (appends) to a file
	void Buffline(double &chron, RunParams &rundat, ChiSquared &chi, ExptsData &edata);//buffer one line
								///for the .hst file	
	void EndNotes(ofstream &file, ChiSquared &chi) const;//writing the calculation method at the end of the *.hst file 
	static void PotWarning(ostream &file, int pind, const char *name);//extremely large potential values occured, give warning
};	
	
//destructor of the History class
inline History::~History()
{
	 if(::debug) cout<<"History destructor"<<endl;
};
