//header Threads.h
//Last changed 12.12.2022

//containing the global thread realated variables

#if !defined _MOVE
	#define _MOVE
	#include "Move.h"
#endif


#if !defined _THREAD_STRUCTS
	#define _THREAD_STRUCTS
 
	struct ThreadArg
	{
		int thread_index;//the index of the thread
	
		//During the histogram and its change calculation all threads should update the same histogram counts array, 
		//and the arrays used for the CoordNumbConst and AvCoordConst, which would mean excessive mutex locking decreasing the speed.
		//To avoid that the storage place of the necessary arrays were multiplied by the number of threads, each thread having its
		//own full segment. The order in the original arrays are the following: first the main thread, than the auxiliary threads' segment
		//in the order of their increasing ID
		//for the histogram calculation
	
		longint **hist_count_finder;//finder array to the partials of the threads histogram counts for each different binsize as well, dim [npartials*ndiffbin]
		longint *hist_tot_start;//pointer to the starting place of the threads histogram totals (does not have to be deleted)
		
		//for histogram calculation
		int *min_xx;//first central atom index to calculate the X-X type pure partials for for each type
		int *max_xx;//last central atom index to calculate  the X-X type pure partials for for each type
		int *min_atom_index;//first cetral atom index to calculate for for each type
		int *max_atom_index;//last central atom index to calculate for for each type

		//for the coordination number constraint
		
		int **coord_numb_finder;//finder to the starting point for each constraint in the coordnumbs array for the thread
		int **nsatisfy_finder;//to the beginning of the nsatisfy segment of each thread
		int *cc_central_offset_min;//the offset of the first central for each constraint (counted from the first central of this constraint) for the given thread
		int *cc_central_offset_max;//the offset of the last central for each constraint (counted from the first central of this constraint) for the given thread
								//both used when the centrals of the constraint are divided between the threads for copying
		
		//for the average coordination number constraint
		int *av_coord_start;//pointer to the starting point in the neighbourcount array
#ifdef _ADVANCED_GEOM_CONST
		//for BVS
		double **valence_finder;//finder to the starting point for each constraint in the valence array for the thread
		int *bvs_central_offset_min;//the offset of the first central for each constraint (counted from the first central of this constraint) for the given thread
		int *bvs_central_offset_max;//the offset of the last central for each constraint (counted from the first central of this constraint) for the given thread
								//both used when the centrals of the constraint are divided between the threads for copying
#endif

		//for histogram bin-based split (histogram, ppcf)
		int *min_bin;//array for each different bin size
		int *max_bin;

#ifdef _VIBR_AMP
		int min_bin_used;
		int max_bin_used;
#endif

#ifdef _USE_LOCAL_INV
		//for local histogram bin-based split 
		int loc_min_bin;
		int loc_max_bin;
		int *loc_hist_count;//pointer to the beginning of each threads local histogram segment
		longint **hist_loc_count_finder;//finder array to the partials of the threads histogram counts
#endif
#ifdef _LOCAL_INV
		int loc_hist_offset;//offset to the beginning of each threads part in the ChiSquared::first_loc_ind and start_count arrays in loc chi2 calc
		int *min_loc_atom_index;//first neighbour atom index to calculate for for each type in case of distance-based calculation
		int *max_loc_atom_index;//last neighbour atom index to calculate for for each type in case of distance-based calculation
		double *av_loc_hist_finder;//pointer to the beginning of the thread's segment in chisquare.av_loc_hist array for each thread in av_loc_hist calc
#endif

		//for r-based split in case of g(r) calculation
		int *min_r_index;//index of the starting value for the cycle/data set
		int *max_r_index;//index of the finishing value for the cycle/data set

		//for Q-based split in case of S(Q),F(Q),F(g) calculation
		int *min_Q_g_index;//index of the starting value for the cycle/data set
		int *max_Q_g_index;//index of the finishing value for the cycle/data set

		//for k-based split in case of E(k) calculation
		int *min_k_index;//index of the starting value for the cycle/data set
		int *max_k_index;//index of the finishing value for the cycle/data set

		int charge_gr_min_ind;//first charge groups to check for NB pot corrections
		int charge_gr_max_ind;//last charge groups to check for NB pot corrections

#ifdef _AENET
		//Although the aenet library is parallelized with mpi, it seems it is possible to call it multiple times form
		//different threads during energy calculation  
        //atom-based split, but as it is used as well during the loop for the neighbours of the moved atom, not
		//segments are stored, but index list
		int *aenet_at_index;//indices of the atoms to calculate for 
		int aenet_at_count;//number of atoms in aenet_at_index to calculate for
		int aenet_max_nneigh;//max number of neighbours, has to kept separately for each thread in case of resizing	
		int aenet_nneigh;//actual number of neighbours
		int *aenet_neighlist;//neighbour list for the actual atom
		int *aenet_neightype;//type of the neighbours
		double *aenet_neighcoord;//coordinates of the neighbours
		double *aenet_E_i;//atomic energies for the atoms pf this thread
#endif
		//as there can be large values for the potential if atoms are close (so mainly for the excluded atoms, and as the
		//calculation of the potential is performed the way that during hostogram calculation all the atoms inside cutoff contribute to the potential
		//and the contribution of the exclusions is removed later the large positiv and large negativ contributions will be collected separately
		//otherwise due to floationg point errors the real contributions will be truncated too much
		double *vdW_pot_start;//pointer array to the starting point of the threads vdW  or tabulated potential
		double *Coulomb_pot_start;//pointer array to the starting point of the threads Coulomb potential
		double *vdW_pot_s_start;//pointer array to the starting point of the threads vdW  or tabulated potential for the small values
		double *Coulomb_pot_s_start;//pointer array to the starting point of the threads Coulomb potential for the small values
		double *vdW_pot_l_start;//pointer array to the starting point of the threads vdW  or tabulated potential for the large values
		double *Coulomb_pot_l_start;//pointer array to the starting point of the threads Coulomb potential for the large values
		
		double *vdW14_pot_start;//pointer array to the starting point of the threads vdW14 potential
		double *Coulomb14_pot_start;//pointer array to the starting point of the threads Coulomb14 potential

		int *ek_range_error;//to indicate, that for the given thread there was no E(k) data point with non-zero value, possible the istogram and coeff range does not overlap
		
	};

#endif	


#if !defined _DEF_EXTERN
	#define _DEF_EXTERN
	extern int acceptable;
#endif

class Threads
{
	public:
	Threads();//default constructor
	~Threads();//destructor

	static int nthreads;//number of threads

	//static Pointers to objects
	static SimpleCfg *conf;
	static HistoSet *histnew_p,*histold_p;
	static ExptsData *edata_p;
	static CoordNumbConst *iccnew_p,*iccold_p;
	static AvCoordConst *avccnew_p,*avccold_p;
#ifdef _ADVANCED_GEOM_CONST
	static BondValenceSumConst *bvsnew_p, *bvsold_p;
#endif
	static Move *move_p;
	static CalcData *calcdat_p;
	static PPCFSet *ppcfnew_p,*ppcfold_p;
	static CalcPart *calcpnew_p,*calcpold_p;
	static Threads *thread_p;
#ifdef _LOCAL_INV
	static ChiSquared *chi_p;
	typedef void (ChiSquared::*TCalcLocChi2Thread)(ThreadArg &);	
	static TCalcLocChi2Thread CalcLocChi2Thread;//pointer to the actual local invariane chi2 function
#endif
#ifdef _NO_PERIODIC
	static DataMat *datmat_p;
#endif
#ifdef _AENET
	static Aenet *aenetnew_p,*aenetold_p;
#endif
	//to replace the corresponding parameter of the Move object, so the values can be preserved for the threads
	static int *mod_partial;//array of modified partials

#ifdef _TEST_MODE
	static std::chrono::time_point<std::chrono::high_resolution_clock> *ctime1, *ctime2;//for CoordNumbConst threaded calc
	static std::chrono::time_point<std::chrono::high_resolution_clock> *thtime1, *thtime2;//for ThreadEntry, ThreadLoop
	static std::chrono::time_point<std::chrono::high_resolution_clock> hctime1, hctime2;//for HistoSet and Chisquared threaded calc
	static double dur_00a,dur_00,dur_01,dur_02,dur_03,dur_04,dur_05,dur_06;
	static double *dur_10,*dur_11,*dur_12,*dur_12a,*dur_12b,*dur_13,*dur_13a,*dur_13b,*dur_14,*dur_15,*dur_15b,*dur_16,*dur_17;
	static std::chrono::duration<double, std::milli> *elapsed;
#ifdef _ADVANCED_GEOM_CONST
	static std::chrono::time_point<std::chrono::high_resolution_clock> *btime1, *btime2;//for BondValenceSumConst threaded calc
	static double *dur_18;
#endif
#ifdef _USE_LOCAL_INV
	static double *dur_loc1,*dur_loc2,*dur_loc3,*dur_loc4;
#endif
#ifdef _LOCAL_INV
	static std::chrono::time_point<std::chrono::high_resolution_clock> *lctime1,*lctime2;//localc inv chi2
	static double *dur_chiloc1,*dur_chiloc2,*dur_chiloc3,*dur_chiloc4;
#endif
#ifdef _AENET
	static std::chrono::time_point<std::chrono::high_resolution_clock> *anntime1, *anntime2;
	static double *dur_ann1,*dur_ann2;
#endif
#endif

	static int loop_flag;//indicator, whether the main loop has to be executed
	static int calc_hist_flag;//1 if the histogram has to be calculated, 0 otherwise
	static int start_flag;//indicator, that the calculation can be started
	
	static int start_count;//counter for the threads to start the loop again
	static int complete_count;//counter for the threads finished the actual calculation
	static int complete_count2;//counter for the threads finished the actual calculation
	static int CC_count;//counter for the threads finished the actual calculation
	static int UC_count;//counter for the threads finished the actual calculation
	static int CMP_count;//counter for the threads finished the actual calculation
	static int T_count;//counter for the threads finished the actual calculation

	static int *EXAFS_gridind;//grid index
	static double *IQ_muact;//actual mu values for I(Q) fitting
	static double *IQ_fprimeact;//actual fprime values only one (for the AXS particle) /data set
	static bool is_E0_shift;//indicator mainly for the threads that in this step only E0 shift is done
	static bool is_IQ_mucorr;//indicator, whether there is a change in mu for I(Q) background correction
	static bool is_fprime_shift;//indicator, whether there is a change in the f'
	//static int return_wait;

	//for signaling the thread to start some work, and wait for completion
	static std::mutex start1_mutex,start2_mutex,start3_mutex,count_mutex,tcl_mutex, cut_mutex;//mutex for locking
	static std::mutex CC_count_mutex, UC_count_mutex, CMP_count_mutex, T_count_mutex;
#ifdef _NO_PERIODIC
	static int  NOP_count;//counter for the threads finished the actual calculation
	static std::mutex NOP_count_mutex;
	static std::condition_variable NOP_check_count;
#endif
#ifdef _AENET
	static std::mutex aenet_mutex;//mutex for locking during hist. calc. and update for aenet neighbours 
	static int calc_ANN;//indicator whether to calculate the ANN potential, same as Aenet::calc_ANN
	static int aenetstart_flag;//indicator, that the final ANN calculation can start
	static std::condition_variable check_aenet_start;//conditional variable
		
#endif
#ifdef _USE_LOCAL_INV
	static int LOC_count, LOC_count2;
	static std::mutex inv_mutex, LOC_count_mutex, LOC_count_mutex2;
	static std::condition_variable check_inv_start, LOC_check_count, LOC_check_count2;//conditional variables for synchronization
#endif
	static std::condition_variable check_hist_start,check_cd_start,check_cp_start;//conditional variables for synchronization
	static std::condition_variable check_count, check_loop_start;//conditional variable for synchronization
	static std::condition_variable CC_check_count, UC_check_count, CMP_check_count, T_check_count;
	
	std::thread *threads;//pointer to th earray of thread objects

	ThreadArg *thread_arg;//the arguments /thread (different for each thread)

	static void SetThreadsParams();//set the parameters
	static void ThreadEntry(ThreadArg &thread_arg);//this is called at thread creation
	static void ThreadLoop(ThreadArg &my_arg);//executing the main loop for the threads
	void SplitToThreads();//spliiting the work among the threads
#ifdef _AENET
	void ResizeAenetArrays(int id, int **ind_list, int **type_list, double **coord_list);
	void SplitAenetAtoms(int n, int *ind);//split the neighbour atoms to calculate for among the threads 
	static void ThreadAenetFinal(ThreadArg &my_arg);//final ANN calculation, if necessary 
#endif	

#ifdef _LOCAL_INV
	void ResetSegment();//resetting the segment boundaries in case of loc_chi2_mode to help load balancing
#endif
};

inline Threads::~Threads()
{
	int ithread;
	for (ithread=0;ithread<nthreads-1;ithread++) threads[ithread].join();
	delete [] threads;
	for (ithread=0;ithread<nthreads;ithread++)
	{
		delete [] thread_arg[ithread].min_atom_index;
		delete [] thread_arg[ithread].max_atom_index;
		delete [] thread_arg[ithread].min_xx;
		delete [] thread_arg[ithread].max_xx;
		delete [] thread_arg[ithread].min_bin;
		delete [] thread_arg[ithread].max_bin;
		delete [] thread_arg[ithread].min_r_index;
		delete [] thread_arg[ithread].max_r_index;
		delete [] thread_arg[ithread].min_Q_g_index;
		delete [] thread_arg[ithread].max_Q_g_index;
		delete [] thread_arg[ithread].min_k_index;
		delete [] thread_arg[ithread].max_k_index;
		delete [] thread_arg[ithread].hist_count_finder;
#ifdef _ADVANCED_GEOM_CONST
		delete [] thread_arg[ithread].valence_finder;
		delete[] thread_arg[ithread].bvs_central_offset_min;
		delete[] thread_arg[ithread].bvs_central_offset_max;
#endif
#ifdef _AENET
		delete[] thread_arg[ithread].aenet_at_index;
		delete[] thread_arg[ithread].aenet_neighlist;
		delete[] thread_arg[ithread].aenet_neightype;
		delete[] thread_arg[ithread].aenet_neighcoord;
		delete[] thread_arg[ithread].aenet_E_i;
#endif
#ifdef _USE_LOCAL_INV
		delete [] thread_arg[ithread].hist_loc_count_finder;
#endif
#ifdef _LOCAL_INV
		delete [] thread_arg[ithread].min_loc_atom_index;
		delete [] thread_arg[ithread].max_loc_atom_index;
#endif
		delete [] thread_arg[ithread].coord_numb_finder;
		delete [] thread_arg[ithread].nsatisfy_finder;
		delete [] thread_arg[ithread].cc_central_offset_min;
		delete [] thread_arg[ithread].cc_central_offset_max;
		delete [] thread_arg[ithread].ek_range_error;
	}
	delete [] thread_arg;
};
