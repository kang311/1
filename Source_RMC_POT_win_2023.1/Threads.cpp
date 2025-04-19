//source Threads.cpp
//Last changed 12.12.2022

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h" 

int	 Threads::nthreads;//number of threads
int *Threads::mod_partial;//array to replace the move.modpart, to preserve the values for the threads
int *Threads::EXAFS_gridind;//grid index
double *Threads::IQ_muact;//actual mu values for I(Q) fitting
double *Threads::IQ_fprimeact;//actual fprime values only one (for the AXS particle) /data set
bool Threads::is_E0_shift;//indicator mainly for the threads that in this step only E0 shift is done
bool Threads::is_IQ_mucorr;//indicator, whether there is a change in mu for I(Q) background correction
bool Threads::is_fprime_shift;//indicator, whether there is a change in the f'
#ifdef _AENET
	int Threads::calc_ANN = 1;//indicator whether to calculate the ANN potential
#endif
SimpleCfg *Threads::conf;
HistoSet *Threads::histnew_p,*Threads::histold_p;//pointers to objects
ExptsData *Threads::edata_p;
CoordNumbConst *Threads::iccnew_p,*Threads::iccold_p;
AvCoordConst *Threads::avccnew_p,*Threads::avccold_p;
#ifdef _ADVANCED_GEOM_CONST
BondValenceSumConst *Threads::bvsnew_p, *Threads::bvsold_p;
#endif
Move *Threads::move_p;
CalcData *Threads::calcdat_p;
PPCFSet *Threads::ppcfnew_p,*Threads::ppcfold_p;
CalcPart *Threads::calcpnew_p,*Threads::calcpold_p;
Threads *Threads::thread_p;
#ifdef _LOCAL_INV
ChiSquared *Threads::chi_p;
#endif
#ifdef _NO_PERIODIC
DataMat *Threads::datmat_p;
#endif
#ifdef _AENET
Aenet *Threads::aenetnew_p,*Threads::aenetold_p;
#endif
#ifdef _TEST_MODE
	std::chrono::time_point<std::chrono::high_resolution_clock> *Threads::ctime1, *Threads::ctime2;//for CoordNumbConst threaded calc
	std::chrono::time_point<std::chrono::high_resolution_clock> *Threads::thtime1, *Threads::thtime2;//for ThreadEntry, ThreadLoop
	std::chrono::time_point<std::chrono::high_resolution_clock> Threads::hctime1, Threads::hctime2;//for HistoSet and Chisquared threaded calc
	double Threads::dur_00a, Threads::dur_00,Threads::dur_01,Threads::dur_02,Threads::dur_03,Threads::dur_04,Threads::dur_05,Threads::dur_06;
	double *Threads::dur_10,*Threads::dur_11,*Threads::dur_12,*Threads::dur_12a,*Threads::dur_12b,*Threads::dur_13,*Threads::dur_13a,*Threads::dur_13b;
	double *Threads::dur_14,*Threads::dur_15,*Threads::dur_15b,*Threads::dur_16,*Threads::dur_17;
	std::chrono::duration<double, std::milli> *Threads::elapsed;
#ifdef _ADVANCED_GEOM_CONST	
	std::chrono::time_point<std::chrono::high_resolution_clock> *Threads::btime1, *Threads::btime2;//for BondValenceSumConst  threaded calc
	double *Threads::dur_18;
#endif
#ifdef _USE_LOCAL_INV	
	double *Threads::dur_loc1,*Threads::dur_loc2,*Threads::dur_loc3,*Threads::dur_loc4;
#endif
#ifdef _LOCAL_INV
	double *Threads::dur_chiloc1,*Threads::dur_chiloc2,*Threads::dur_chiloc3,*Threads::dur_chiloc4;
	std::chrono::time_point<std::chrono::high_resolution_clock> *Threads::lctime1,*Threads::lctime2;//localc inv chi2
#endif
#ifdef _AENET
	std::chrono::time_point<std::chrono::high_resolution_clock> *Threads::anntime1, *Threads::anntime2;
	double *Threads::dur_ann1, *Threads::dur_ann2;
#endif
#endif

	extern int acceptable;

	int  Threads::loop_flag=1;//indicator, whether the main loop has to be executed
	int  Threads::calc_hist_flag=0;//1 if the histogram has to be calculated, 0 otherwise
	int  Threads::start_flag=0;//indicator, that the calculation can be started
	int  Threads::start_count=0;//counter for the threads to start the loop again
	int  Threads::complete_count=0;//counter for the threads finished the actual calculation
	int  Threads::complete_count2=0;//counter for the threads finished the actual calculation
	int  Threads::CC_count = 0;//counter for the threads finished the actual calculation
	int  Threads::UC_count = 0;//counter for the threads finished the actual calculation
	int  Threads::CMP_count = 0;//counter for the threads finished the actual calculation
	int  Threads::T_count = 0;//counter for the threads finished the actual calculation
	//int  Threads::return_wait;


	std::mutex Threads::start1_mutex;//mutex for locking during histogram calculation
	std::mutex Threads::start2_mutex;//mutex for locking during calculated data calculation
	std::mutex Threads::start3_mutex;//mutex for locking during data copy
	std::mutex Threads::count_mutex;//mutex for locking
	std::mutex Threads::tcl_mutex;//mutex for locking during hist. calc. and update for too close array update 
	std::mutex Threads::cut_mutex;//mutex for locking during hist. calc for automatic cutoff determination update
	std::mutex Threads::CC_count_mutex;
	std::mutex Threads::UC_count_mutex;
	std::mutex Threads::CMP_count_mutex;
	std::mutex Threads::T_count_mutex;
#ifdef _AENET
	int  Threads::aenetstart_flag = 0;//indicator, that the final ANN calculation can start
	std::mutex Threads::aenet_mutex;//mutex for locking during hist. calc. and update for aenet neighbours 
	std::condition_variable Threads::check_aenet_start;//conditional variable
#endif
#ifdef _NO_PERIODIC
	int  Threads::NOP_count = 0;//counter for the threads finished the actual calculation
	std::mutex Threads::NOP_count_mutex;
	std::condition_variable Threads::NOP_check_count;
#endif
#ifdef _USE_LOCAL_INV
	int Threads::LOC_count=0, Threads::LOC_count2 = 0;
	std::mutex Threads::inv_mutex,Threads::LOC_count_mutex, Threads::LOC_count_mutex2;//mutex for locking during hist calculation the local histogram
	std::condition_variable Threads::check_inv_start, Threads::LOC_check_count, Threads::LOC_check_count2;//conditional variable for the local invariance chisquare synchronization
#endif
	std::condition_variable Threads::check_hist_start;//conditional variable for the histogram synchronization
	std::condition_variable Threads::check_cd_start;//conditional variable for the calcdata synchronization
	std::condition_variable Threads::check_cp_start;//conditional variable for the data copy synchronization
	std::condition_variable Threads::check_count;//conditional variable
	std::condition_variable Threads::check_loop_start;//conditional variable
	std::condition_variable Threads::CC_check_count;
	std::condition_variable Threads::UC_check_count;
	std::condition_variable Threads::CMP_check_count;
	std::condition_variable Threads::T_check_count;
	
#ifdef _LOCAL_INV
	Threads::TCalcLocChi2Thread Threads::CalcLocChi2Thread;//pointer to the actual local invariance chi2 function
#endif

//default constructor
Threads::Threads()
{
	int ithread;
	int size, size2;
	char *name=NULL, *name2=NULL, *conv_numb = NULL;
	SetArraysize(&name, NAME_SIZE, "name", "Threads::Threads");
	SetArraysize(&name2, NAME_SIZE, "name2", "Threads::Threads");

	SetArraysize(&thread_arg,nthreads,"thread_arg","Threads::Threads");
	size = edata_p->nsq + edata_p->nfq + edata_p->nfg;//Q and g-based split for S(Q), F(Q) and F(g), technically F(g) is the same as the other two
	size2=(RunParams::nicoord>0  ? RunParams::nicoord : 1);//to avoid problems at the destructor
#ifdef _ADVANCED_GEOM_CONST
	int size3;
	size3 = (RunParams::nbvs > 0 ? RunParams::nbvs : 1);//to avoid problems at the destructor
#endif
	for (ithread=0;ithread<nthreads;ithread++)
	{
		mystrcpy(name, NAME_SIZE, "thread_arg[");
		IntToStr(&conv_numb, ithread + 1);
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "].");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_atom_index");
		SetArraysize(&thread_arg[ithread].min_atom_index,RunParams::ntypes + SimpleCfg::nvirtualtypes,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_atom_index");
		SetArraysize(&thread_arg[ithread].max_atom_index,RunParams::ntypes + SimpleCfg::nvirtualtypes,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_xx");
		SetArraysize(&thread_arg[ithread].min_xx,RunParams::ntypes+SimpleCfg::nvirtualtypes,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_xx");
		SetArraysize(&thread_arg[ithread].max_xx,RunParams::ntypes + SimpleCfg::nvirtualtypes,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_bin");
		SetArraysize(&thread_arg[ithread].min_bin, RunParams::ndiffbin, name2, "Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_bin");
		SetArraysize(&thread_arg[ithread].max_bin, RunParams::ndiffbin, name2, "Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_r_index");
		SetArraysize(&thread_arg[ithread].min_r_index,edata_p->ngr,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_r_index");
		SetArraysize(&thread_arg[ithread].max_r_index,edata_p->ngr,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_Q_g_index");
		SetArraysize(&thread_arg[ithread].min_Q_g_index,size,name2,"Threads::Threads");//it handles F(g) as well
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_Q_g_index");
		SetArraysize(&thread_arg[ithread].max_Q_g_index,size,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_k_index");
		SetArraysize(&thread_arg[ithread].min_k_index,edata_p->nek,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_k_index");
		SetArraysize(&thread_arg[ithread].max_k_index,edata_p->nek,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "hist_count_finder");
		SetArraysize(&thread_arg[ithread].hist_count_finder,(RunParams::ntypes*(RunParams::ntypes+1)/2)*RunParams::ndiffbin,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "coord_numb_finder");
		SetArraysize(&thread_arg[ithread].coord_numb_finder,size2,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "nsatisfy_finder");
		SetArraysize(&thread_arg[ithread].nsatisfy_finder,size2,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "cc_central_offset_min");
		SetArraysize(&thread_arg[ithread].cc_central_offset_min,size2,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "cc_central_offset_max");
		SetArraysize(&thread_arg[ithread].cc_central_offset_max,size2,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "ek_range_error");
		SetArraysize(&thread_arg[ithread].ek_range_error, RunParams::nek, name2, "Threads::Threads");
#ifdef _ADVANCED_GEOM_CONST
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "valence_finder");
		SetArraysize(&thread_arg[ithread].valence_finder, size3, name2, "Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "bvs_central_offset_min");
		SetArraysize(&thread_arg[ithread].bvs_central_offset_min, size3, name2, "Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "bvs_central_offset_max");
		SetArraysize(&thread_arg[ithread].bvs_central_offset_max, size3, name2, "Threads::Threads");
#endif
#ifdef _USE_LOCAL_INV
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "hist_loc_count_finder");
		SetArraysize(&thread_arg[ithread].hist_loc_count_finder,RunParams::ntypes*(RunParams::ntypes+1)/2,name2,"Threads::Threads");
#endif
#ifdef _LOCAL_INV
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "min_loc_atom_index");
		SetArraysize(&thread_arg[ithread].min_loc_atom_index,RunParams::ntypes,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "max_loc_atom_index");
		SetArraysize(&thread_arg[ithread].max_loc_atom_index,RunParams::ntypes,name2,"Threads::Threads");
#endif
#ifdef _AENET
		
		thread_arg[ithread].aenet_max_nneigh=(int)(SimpleCfg::ntotal/8.*4./3.*pow(Aenet::cutoff/SimpleCfg::boxedge,3)*PI*1.1);//initial value
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "aenet_at_index");
		SetArraysize(&thread_arg[ithread].aenet_at_index,SimpleCfg::ntotal/nthreads+1,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "aenet_neighlist");
		int padding_offset=(int)((CACHE_PADDING-sizeof(*thread_arg[ithread].aenet_neighlist))/sizeof(*thread_arg[ithread].aenet_neighlist));//number of dummy elements
		SetArraysize(&thread_arg[ithread].aenet_neighlist,thread_arg[ithread].aenet_max_nneigh,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "aenet_neightype");
		padding_offset=(int)((CACHE_PADDING-sizeof(*thread_arg[ithread].aenet_neightype))/sizeof(*thread_arg[ithread].aenet_neightype));//number of dummy elements
		SetArraysize(&thread_arg[ithread].aenet_neightype,thread_arg[ithread].aenet_max_nneigh,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "aenet_neighcoord");
		SetArraysize(&thread_arg[ithread].aenet_neighcoord,thread_arg[ithread].aenet_max_nneigh*3,name2,"Threads::Threads");
		mystrcpy(name2, NAME_SIZE, name);
		mystrcat(name2, NAME_SIZE, "aenet_E_i");
		padding_offset=(int)((CACHE_PADDING-sizeof(*thread_arg[ithread].aenet_E_i))/sizeof(*thread_arg[ithread].aenet_E_i));//number of dummy elements
		SetArraysize(&thread_arg[ithread].aenet_E_i,SimpleCfg::ntotal/nthreads+padding_offset,name2,"Threads::Threads");
#endif

	}
	SetArraysize(&threads,(nthreads-1),"threads","Threads::Threads");

	if (::debug)
		cout<<"\nThread object is created "<<this<<endl;
	if (name != NULL)
		delete[] name;
	if (name2 != NULL)
		delete[] name2;
}

//-----------set the static parameters------------------
void Threads::SetThreadsParams()
{
	int i;
	nthreads=RunParams::nthreads;//total number of threads including main
	
	//Used, even in case of consecutive compilation!
	//creating the array for the modified partials
	//this array belonged originally to Move. As copying is done multithreaded, for the new threads the values of these arrays had 
	//to be preserved for the duration of the copying by the HistoSet, PPCFSet and CalcPart objects, these arrays could be set just
	//before the signaling to start the histogram change calculation, because then all of the threads finished copying .
	//Therefore the arrays are relocated to the Thread object, and the Move counterparts are removed.
	SetArraysize(&mod_partial,RunParams::ntypes*(RunParams::ntypes+1)/2,"mod_partial","Threads::SetThreadsParams");
		//setting the mod_partial indicator to 1 for the case of ntypes=1 or custmove=1, as it does not change during the simulation
	//it is recalculated in each loop, if ntypes>1, and custmove=0
	
	//initialising the array of modified partial indicators
	for (i = 0; i < RunParams::ntypes * (RunParams::ntypes + 1) / 2; i++)
		mod_partial[i] = 1;//this is the default, it will be reset, when needed
	
	SetArraysize(&EXAFS_gridind, ExptsData::nek, "EXAFS_gridind", "Threads::SetThreadsParams");
	for (i = 0; i < ExptsData::nek; i++)
		EXAFS_gridind[i] = ExptsData::ek_gridind[i];
	SetArraysize(&IQ_muact, ExptsData::nfq, "IQ_muact", "Threads::SetThreadsParams");
	SetArraysize(&IQ_fprimeact, ExptsData::nfq, "IQ_fprimeact", "Threads::SetThreadsParams");
	for (i = 0; i < ExptsData::nfq; i++)
	{
		IQ_muact[i] = ExptsData::fqmuact[i];
		IQ_fprimeact[i] = ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1];
	}

#ifdef _TEST_MODE
	int ithread;
	SetArraysize(&elapsed,nthreads,"elapsed","Threads::SetThreadsParams");
	SetArraysize(&ctime1,nthreads,"ctime1","Threads::SetThreadsParams");
	SetArraysize(&ctime2,nthreads,"ctime2","Threads::SetThreadsParams");
	SetArraysize(&thtime1,nthreads,"thtime1","Threads::SetThreadsParams");
	SetArraysize(&thtime2,nthreads,"thtime2","Threads::SetThreadsParams");
	SetArraysize(&dur_10,nthreads,"dur_10","Threads::SetThreadsParams");
	SetArraysize(&dur_11,nthreads,"dur_11","Threads::SetThreadsParams");
	SetArraysize(&dur_12,nthreads,"dur_12","Threads::SetThreadsParams");
	SetArraysize(&dur_12a,nthreads,"dur_12a","Threads::SetThreadsParams");
	SetArraysize(&dur_12b,nthreads,"dur_12b","Threads::SetThreadsParams");
	SetArraysize(&dur_13,nthreads,"dur_13","Threads::SetThreadsParams");
	SetArraysize(&dur_13a,nthreads,"dur_13a","Threads::SetThreadsParams");
	SetArraysize(&dur_13b,nthreads,"dur_13b","Threads::SetThreadsParams");
	SetArraysize(&dur_14,nthreads,"dur_14","Threads::SetThreadsParams");
	SetArraysize(&dur_15,nthreads,"dur_15","Threads::SetThreadsParams");
	SetArraysize(&dur_15b,nthreads,"dur_15b","Threads::SetThreadsParams");
	SetArraysize(&dur_16,nthreads,"dur_16","Threads::SetThreadsParams");
	SetArraysize(&dur_17,nthreads,"dur_17","Threads::SetThreadsParams");
#ifdef _ADVANCED_GEOM_CONST	
	SetArraysize(&btime1, nthreads, "btime1", "Threads::SetThreadsParams");
	SetArraysize(&btime2, nthreads, "btime2", "Threads::SetThreadsParams");
	SetArraysize(&dur_18, nthreads, "dur_18", "Threads::SetThreadsParams");
#endif
#ifdef _USE_LOCAL_INV
	SetArraysize(&dur_loc1,nthreads,"dur_loc1","Threads::SetThreadsParams");
	SetArraysize(&dur_loc2,nthreads,"dur_loc2","Threads::SetThreadsParams");
	SetArraysize(&dur_loc3,nthreads,"dur_loc3","Threads::SetThreadsParams");
	SetArraysize(&dur_loc4,nthreads,"dur_loc4","Threads::SetThreadsParams");
#endif
#ifdef _LOCAL_INV
	SetArraysize(&dur_chiloc1,nthreads,"dur_chiloc1","Threads::SetThreadsParams");
	SetArraysize(&dur_chiloc2,nthreads,"dur_chiloc2","Threads::SetThreadsParams");
	SetArraysize(&dur_chiloc3,nthreads,"dur_chiloc3","Threads::SetThreadsParams");
	SetArraysize(&dur_chiloc4,nthreads,"dur_chiloc4","Threads::SetThreadsParams");
	SetArraysize(&lctime1,nthreads,"lctime1","Threads::SetThreadsParams");
	SetArraysize(&lctime2,nthreads,"lctime2","Threads::SetThreadsParams");
#endif
#ifdef _AENET
	SetArraysize(&dur_ann1,nthreads,"dur_ann1","Threads::SetThreadsParams");
	SetArraysize(&dur_ann2,nthreads,"dur_ann2","Threads::SetThreadsParams");
	SetArraysize(&anntime1,nthreads,"anntime1","Threads::SetThreadsParams");
	SetArraysize(&anntime2,nthreads,"anntime2","Threads::SetThreadsParams");
#endif
	for (ithread=0;ithread<nthreads;ithread++)
	{
		dur_10[ithread]=0;
		dur_11[ithread]=0;
		dur_12[ithread]=0;
		dur_12a[ithread]=0;
		dur_12b[ithread]=0;
		dur_13[ithread]=0;
		dur_13a[ithread]=0;
		dur_13b[ithread]=0;
		dur_14[ithread]=0;
		dur_15[ithread]=0;
		dur_15b[ithread]=0;
		dur_16[ithread]=0;
		dur_17[ithread]=0;
#ifdef _ADVANCED_GEOM_CONST
		dur_18[ithread] = 0;
#endif
#ifdef _USE_LOCAL_INV
		dur_loc1[ithread]=0;
		dur_loc2[ithread]=0;
		dur_loc3[ithread]=0;
		dur_loc4[ithread]=0;
#endif
#ifdef _LOCAL_INV
		dur_chiloc1[ithread]=0;
		dur_chiloc2[ithread]=0;
		dur_chiloc3[ithread]=0;
		dur_chiloc4[ithread]=0;
#endif
#ifdef _AENET
		dur_ann1[ithread]=0;
		dur_ann2[ithread]=0;
#endif
	}
#endif
#ifdef _LOCAL_INV
	if (RunParams::loc_chi2_mode==1)
		CalcLocChi2Thread=&ChiSquared::CalcLocChiSquaredThread_dist;
	else
		CalcLocChi2Thread=&ChiSquared::CalcLocChiSquaredThread_bin;
#endif
};

//---------------spliiting the task between the threads-----------------
void Threads::SplitToThreads()
{
	int ib,ithread,iexpt,itype,ipartial,iconst;
	double N_pairs;
	int nexpt,*npoints;
	int step,remainder,offset;//for splitting the work between threads
	const int npartials=RunParams::ntypes*(RunParams::ntypes+1)/2;
	const int ndiffbin = RunParams::ndiffbin;
	double  n_pairs_thread;

	int nused_threads=0;



	//to avoid warning
	npoints=NULL;

	//the histogram's pcounts array will contain the data in the following order
	//for ithread
	//for binsize cycle
	//for partials
	//setting the finder of the thread's histogram count
	for (ithread=0; ithread<nthreads; ithread++)
	{
		thread_arg[ithread].hist_count_finder[0]=histnew_p->thread_histcount_finder[ithread];//first setting the starting point for each thread
#ifdef _USE_LOCAL_INV
		thread_arg[ithread].hist_loc_count_finder[0]=histnew_p->thread_lochistcount_finder[ithread];//first setting the starting point for each thread
#endif
		offset = 0;
		//setting the remaining partials
		for (ib=0;ib<ndiffbin;ib++)
		{
			for (ipartial = 0; ipartial < npartials; ipartial++)
			{
				if (!(ib == 0 && ipartial == 0))//the first of the thread was already set
				{
					thread_arg[ithread].hist_count_finder[ib*npartials+ipartial] = thread_arg[ithread].hist_count_finder[0] + offset;
									
						
#ifdef _USE_LOCAL_INV

					if (ib == 0)
						thread_arg[ithread].hist_loc_count_finder[ipartial] = thread_arg[ithread].hist_loc_count_finder[0] + ipartial * RunParams::max_loc_nbins;
#endif
				}
				offset+= RunParams::nbins[ib];
			}
		}
	}

#ifdef _AENET
//splitting the atoms for ANN potential calculation among the threads for the initial calculation
//as the same arrays are used later in the loop it has to be an index list not segments
	if ((int)SimpleCfg::ntotal>nthreads)
	{
		step=(int)SimpleCfg::ntotal/nthreads;
		remainder=(SimpleCfg::ntotal% nthreads);
		offset=step;//only assign step atoms for the original thread, regardless of remainder
		if (remainder>0)
		{
			//to spread the remainder
			offset++;
			remainder--;
		}
		nused_threads=nthreads;
	}
	else
	{
		step=1;
		offset=step;
		remainder=0;
		nused_threads=(int)SimpleCfg::ntotal;
	}
	//assign for ithread=0 0->offset-1
	for (int j=0;j<offset;j++)
		thread_arg[0].aenet_at_index[j]=j;
	thread_arg[0].aenet_at_count=offset;

	if (nthreads>1)
	{
		int min;
		for (ithread=1;ithread<nused_threads-1;ithread++)
		{
			min=thread_arg[ithread-1].aenet_at_index[thread_arg[ithread-1].aenet_at_count-1]+1;
			offset=step;
			if (remainder>0)
			{
				//to spread the remainder
				offset=step+1;
				remainder--;
			} 
			for (int j=0;j<offset;j++)
				thread_arg[ithread].aenet_at_index[j]=min+j;
			thread_arg[ithread].aenet_at_count=offset;
		}
		min=thread_arg[nused_threads-2].aenet_at_index[thread_arg[nused_threads-2].aenet_at_count-1]+1;
		thread_arg[nused_threads-1].aenet_at_count=SimpleCfg::ntotal-min;
		for (int j=0;j<thread_arg[nused_threads-1].aenet_at_count;j++)
			thread_arg[nused_threads-1].aenet_at_index[j]=min+j;
		
		for (ithread=nused_threads;ithread<nthreads;ithread++)
		{
			thread_arg[ithread].aenet_at_index[0]=-1;//not to do anything
			thread_arg[ithread].aenet_at_count=0;

		}
	}
#endif
	//setting the beginning of the histogram's total count for each thread
	 //Paddinng will be used, to make sure that the segments of different threads will not be cached together.
	offset=(histnew_p->nbins[0]>0 ? (int)((CACHE_PADDING-sizeof(*histnew_p->ptotal))/sizeof(*histnew_p->ptotal)) : 0);//number of dummy longinteger elements
	thread_arg[nthreads-1].hist_tot_start=histnew_p->ptotal;// for the main thread
	
	for (ithread = 0; ithread < nthreads - 1; ithread++)//for the auxiliary threads
		thread_arg[ithread].hist_tot_start = histnew_p->ptotal + (ithread + 1)*(ndiffbin* npartials + offset);

	//setting the pointers for the beginning of each CoordNumbConst in the given threads segment in coordnumbconst array and nsatisfy array
	 //Paddinng will be used, to make sure that the segments of different threads will not be cached together.
	//nsatisfy does not follow the pattern, that the nthreads-1. thread is the main for nsatisfy calculation, but does not matter, as later the
	//sum is collected at the first segment of the array
	offset=(iccnew_p->tot_subconst>0 ? (int)((CACHE_PADDING-sizeof(*iccnew_p->nsatisfy))/sizeof(*iccnew_p->nsatisfy)) : 0);//number of dummy integer elements
	for (ithread=0;ithread<nthreads;ithread++)
	{
		thread_arg[ithread].coord_numb_finder[0]=iccnew_p->thread_finder[ithread];//first setting the starting point for each thread
		thread_arg[ithread].nsatisfy_finder[0]=iccnew_p->nsatisfy+ithread*(CoordNumbConst::tot_subconst+offset);
		for (iconst=1;iconst<CoordNumbConst::nconstraints;iconst++)
		{
			thread_arg[ithread].coord_numb_finder[iconst]=thread_arg[ithread].coord_numb_finder[0]+CoordNumbConst::ncentral_cum[iconst];
			thread_arg[ithread].nsatisfy_finder[iconst]=iccnew_p->nsatisfy+ithread*(CoordNumbConst::tot_subconst+offset)+iccnew_p->cum_n_subconst[iconst];
		}
	}

	//splitting the central atoms of the coordinaton constraints among the threads
	for (iconst=0;iconst<CoordNumbConst::nconstraints;iconst++)
	{
		if ((int)iccnew_p->ncentral[iconst]>nthreads)
		{
			step=(int)iccnew_p->ncentral[iconst]/nthreads;
			remainder=(iccnew_p->ncentral[iconst] % nthreads);
			offset=step;//only assign step centrals for the original thread, regardless of remainder
			if (remainder>0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads=nthreads;
		}
		else
		{
			step=1;
			offset=step;
			remainder=0;
			nused_threads=(int)iccnew_p->ncentral[iconst];
		}
		//assign for ithread=0
		thread_arg[0].cc_central_offset_min[iconst]=0;
		thread_arg[0].cc_central_offset_max[iconst]=offset-1;
		if (nthreads>1)
		{
			for (ithread=1;ithread<nused_threads-1;ithread++)
			{
				thread_arg[ithread].cc_central_offset_min[iconst]=thread_arg[ithread-1].cc_central_offset_max[iconst]+1;
				offset=step;
				if (remainder>0)
				{
					//to spread the remainder
					offset=step+1;
					remainder--;
				} 
				thread_arg[ithread].cc_central_offset_max[iconst]=thread_arg[ithread].cc_central_offset_min[iconst]+offset-1;
			}
			thread_arg[nused_threads-1].cc_central_offset_min[iconst]=thread_arg[nused_threads-2].cc_central_offset_max[iconst]+1;
			thread_arg[nused_threads-1].cc_central_offset_max[iconst]=iccnew_p->ncentral[iconst]-1;//to make sure, that all is assigned, and no more
			for (ithread=nused_threads;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].cc_central_offset_min[iconst]=-1;//not to do anything
				thread_arg[ithread].cc_central_offset_max[iconst]=-2;

			}
		}
	}

	//setting the beginning of the neighbourcount of the AvCoordConst
	offset = (avccnew_p->nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*avccnew_p->neighbourcount)) / sizeof(*avccnew_p->neighbourcount)) : 0);//number of dummy integer elements
	thread_arg[nthreads-1].av_coord_start=avccnew_p->neighbourcount;// for the main thread
	for (ithread=0;ithread<nthreads-1;ithread++)//for the auxiliary threads
		thread_arg[ithread].av_coord_start=avccnew_p->neighbourcount+(ithread+1)*(AvCoordConst::nconstraints+offset);

	//Splitting the work equally among the threads
	//first for the calculation of the histogram and its change:each thread will calculate for each partial of each different binsize.
	//For the initial calculation of the X-Y mixed partials splitting the atoms of type X among the threads,
	//and for the histogram calc. change (neighbour atoms)
	//for each atom-types the atoms will be divided equally among the threads
	//this has to include the virtual sites, if there is any
	//If there are more threads
	for (itype=0;itype<RunParams::ntypes + SimpleCfg::nvirtualtypes;itype++)
	{
		if (conf->pnatoms[itype]>nthreads)
		{
			step=(int)conf->pnatoms[itype]/nthreads;
			remainder=((int)conf->pnatoms[itype] % nthreads);
			offset=step;//only assign step centrals for the original thread, regardless of remainder
			if (remainder>0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads=nthreads;
		}
		else
		{
			step=1;
			offset=step;
			remainder=0;
			nused_threads=(int)conf->pnatoms[itype];
		}
		
		//assign for ithread=0
		thread_arg[0].min_atom_index[itype]=0;
		thread_arg[0].max_atom_index[itype]=offset-1;
		if (nused_threads>1)
		{
			for (ithread=1;ithread<nused_threads-1;ithread++)
			{
				thread_arg[ithread].min_atom_index[itype]=thread_arg[ithread-1].max_atom_index[itype]+1;
				offset=step;
				if (remainder>0)
				{
					//to spread the remainder
					offset=step+1;
					remainder--;
				} 
				thread_arg[ithread].max_atom_index[itype]=thread_arg[ithread].min_atom_index[itype]+offset-1;
			}
			thread_arg[nused_threads-1].min_atom_index[itype]=thread_arg[nused_threads-2].max_atom_index[itype]+1;
			thread_arg[nused_threads-1].max_atom_index[itype]=conf->pnatoms[itype]-1;//to make sure, that all is assigned, and no more
		}
		for (ithread=nused_threads;ithread<nthreads;ithread++)
		{
		  thread_arg[ithread].min_atom_index[itype]=-1;//not to do anything
		  thread_arg[ithread].max_atom_index[itype]=-2;
		}
	}//end of itype cycle

#ifdef _ADVANCED_GEOM_CONST
	//setting the pointers for the beginning of each BVS constraint in the given threads segment in the valence array
	 //Paddinng will be used, to make sure that the segments of different threads will not be cached together.
	offset = (bvsnew_p->nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*bvsnew_p->valence)) / sizeof(*bvsnew_p->valence)) : 0);//number of dummy integer elements
	for (ithread = 0; ithread < nthreads; ithread++)
	{
		thread_arg[ithread].valence_finder[0] = bvsnew_p->thread_finder[ithread];//first setting the starting point for each thread
		for (iconst = 1; iconst < BondValenceSumConst::nconstraints; iconst++)
			thread_arg[ithread].valence_finder[iconst] = thread_arg[ithread].valence_finder[0] + BondValenceSumConst::ncentral_cum[iconst];
	}
	//splitting the central atoms of the BVS constraints among the threads
	for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
	{
		if ((int)bvsnew_p->ncentral[iconst] > nthreads)
		{
			step = (int)bvsnew_p->ncentral[iconst] / nthreads;
			remainder = (bvsnew_p->ncentral[iconst] % nthreads);
			offset = step;//only assign step centrals for the original thread, regardless of remainder
			if (remainder > 0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads = nthreads;
		}
		else
		{
			step = 1;
			offset = step;
			remainder = 0;
			nused_threads = (int)bvsnew_p->ncentral[iconst];
		}
		//assign for ithread=0
		thread_arg[0].bvs_central_offset_min[iconst] = 0;
		thread_arg[0].bvs_central_offset_max[iconst] = offset - 1;
		if (nthreads > 1)
		{
			for (ithread = 1; ithread < nused_threads - 1; ithread++)
			{
				thread_arg[ithread].bvs_central_offset_min[iconst] = thread_arg[ithread - 1].bvs_central_offset_max[iconst] + 1;
				offset = step;
				if (remainder > 0)
				{
					//to spread the remainder
					offset = step + 1;
					remainder--;
				}
				thread_arg[ithread].bvs_central_offset_max[iconst] = thread_arg[ithread].bvs_central_offset_min[iconst] + offset - 1;
			}
			thread_arg[nused_threads - 1].bvs_central_offset_min[iconst] = thread_arg[nused_threads - 2].bvs_central_offset_max[iconst] + 1;
			thread_arg[nused_threads - 1].bvs_central_offset_max[iconst] = bvsnew_p->ncentral[iconst] - 1;//to make sure, that all is assigned, and no more
			for (ithread = nused_threads; ithread < nthreads; ithread++)
			{
				thread_arg[ithread].bvs_central_offset_min[iconst] = -1;//not to do anything
				thread_arg[ithread].bvs_central_offset_max[iconst] = -2;

			}
		}
	}
#endif
	//Initial histogram calculation for X-X type pure partials, each pair is calculated only once, for each central, only the neighbours
	//with smaller index than the central's are considered-> the number of neighbours are increasing with increasing central index
	//this has to include the virtual sites, if there is any
	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::right, ios::adjustfield);
	//it is not likely, that the number of pairs would be smaller, than the number of threads, it is not checked
	for (itype=0;itype<RunParams::ntypes + SimpleCfg::nvirtualtypes;itype++)
	{
		//total number of pairs to calculate for
		N_pairs=(double)conf->pnatoms[itype]*((double)conf->pnatoms[itype]-1)/2;
		n_pairs_thread=(double)N_pairs/nthreads;
		thread_arg[0].min_xx[itype]=0;
		for (ithread=0;ithread<nthreads-1;ithread++)
		{
			step=RoundNearInt((1+sqrt(1+8*(ithread+1)*n_pairs_thread))/2);
			thread_arg[ithread].max_xx[itype]=step-1;
			thread_arg[ithread+1].min_xx[itype]=step;
		}
		thread_arg[nthreads-1].max_xx[itype]=conf->pnatoms[itype]-1;//to make sure, that all is assigned, and no more

	}//end of itype cycle

	//Making the split based on the histogram bins (for the PPCF calculation)
	//Each thread will handle a consecutive segment of bins
	//The split is made based on the total number of bins, if binshift >1 for a data set it has to be handled in CalcPart gr calculation
	for (ib = 0; ib < ndiffbin; ib++)
	{
		if (RunParams::nbins[ib] > nthreads)
		{
			step = RunParams::nbins[ib] / nthreads;
			remainder = (RunParams::nbins[ib] % nthreads);
			offset = step;//only assign step atoms for the original thread, regardless of remainder
			if (remainder > 0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads = nthreads;
		}
		else
		{
			step = 1;
			offset = step;
			remainder = 0;
			nused_threads = RunParams::nbins[ib];

		}
		//assign for ithread=0
		thread_arg[0].min_bin[ib] = 0;
		thread_arg[0].max_bin[ib] = offset - 1;
		if (nthreads > 1)
		{
			for (ithread = 1; ithread < nused_threads - 1; ithread++)
			{
				thread_arg[ithread].min_bin[ib] = thread_arg[ithread - 1].max_bin[ib] + 1;
				offset = step;
				if (remainder > 0)
				{
					//to spread the remainder
					offset = step + 1;
					remainder--;
				}
				thread_arg[ithread].max_bin[ib] = thread_arg[ithread].min_bin[ib] + offset - 1;
			}
			thread_arg[nused_threads - 1].min_bin[ib] = thread_arg[nused_threads - 2].max_bin[ib] + 1;
			thread_arg[nused_threads - 1].max_bin[ib] = RunParams::nbins[ib] - 1;
			for (ithread = nused_threads; ithread < nthreads; ithread++)
			{
				thread_arg[ithread].min_bin[ib] = -1;//not to do anything
				thread_arg[ithread].max_bin[ib] = -2;

			}

		}
	}
#ifdef _VIBR_AMP
	//Making the split based on the used histogram bins (for the PPCF calculation)
	//Each thread will handle a consecutive segment of bins	
	if (HistoSet::nbins_used>nthreads)
	{
		step=HistoSet::nbins_used/nthreads;
		remainder=(HistoSet::nbins_used % nthreads);
		offset=step;//only assign step atoms for the original thread, regardless of remainder
		if (remainder>0)
		{
			//to spread the remainder
			offset++;
			remainder--;
		}
		nused_threads=nthreads;
	}
	else
	{
		step=1;
		offset=step;
		remainder=0;
		nused_threads=HistoSet::nbins_used;

	}
	//assign for ithread=0
	thread_arg[0].min_bin_used=0;
	thread_arg[0].max_bin_used=offset-1;
	if (nthreads>1)
	{
		for (ithread=1;ithread<nused_threads-1;ithread++)
		{
			thread_arg[ithread].min_bin_used=thread_arg[ithread-1].max_bin_used+1;
			offset=step;
			if (remainder>0)
			{
				//to spread the remainder
				offset=step+1;
				remainder--;
			} 
			thread_arg[ithread].max_bin_used=thread_arg[ithread].min_bin_used+offset-1;
		}
		thread_arg[nused_threads-1].min_bin_used=thread_arg[nused_threads-2].max_bin_used+1;
		thread_arg[nused_threads-1].max_bin_used=HistoSet::nbins_used-1;
		for (ithread=nused_threads;ithread<nthreads;ithread++)
		{
			thread_arg[ithread].min_bin_used=-1;//not to do anything
			thread_arg[ithread].max_bin_used=-2;

		}

	}
#endif
	//Making the split based on the r data points (for g(r) fitting)
	//Each thread will handle a consecutive segment of the r point range
	for (iexpt=0;iexpt<RunParams::ngr;iexpt++)
	{
		if ((int)ExptsData::grused[iexpt]>nthreads)
		{
			step=(int)ExptsData::grused[iexpt]/nthreads;
			remainder=(ExptsData::grused[iexpt] % nthreads);
			offset=step;//only assign step atoms for the original thread, regardless of remainder
			if (remainder>0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads=nthreads;
		}
		else
		{
			step=1;
			offset=step;
			remainder=0;
			nused_threads=(int)ExptsData::grused[iexpt];

		}
		//assign for ithread=0
		thread_arg[0].min_r_index[iexpt]=0;
		thread_arg[0].max_r_index[iexpt]=offset-1;
		if (nthreads>1)
		{
			for (ithread=1;ithread<nused_threads-1;ithread++)
			{
				thread_arg[ithread].min_r_index[iexpt]=thread_arg[ithread-1].max_r_index[iexpt]+1;
				offset=step;
				if (remainder>0)
				{
					//to spread the remainder
					offset=step+1;
					remainder--;
				} 
				thread_arg[ithread].max_r_index[iexpt]=thread_arg[ithread].min_r_index[iexpt]+offset-1;
			}
			thread_arg[nused_threads-1].min_r_index[iexpt]=thread_arg[nused_threads-2].max_r_index[iexpt]+1;
			thread_arg[nused_threads-1].max_r_index[iexpt]=ExptsData::grused[iexpt]-1;
			for (ithread=nused_threads;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].min_r_index[iexpt]=-1;//not to do anything
				thread_arg[ithread].max_r_index[iexpt]=-2;

			}
		}
	}//end of iexpt cycle

	//Making the split based on the Q and g data points (for S(Q), F(Q), F(g)), they will be kept in the same two arrays 
	//calculation-wise the electron diffraction F(g) is the same as X-ray diffraction F(Q) and neutron S(Q)
	//Each thread will handle a consecutive segment of the Q point range
	nexpt = ExptsData::nsq + ExptsData::nfq + ExptsData::nfg;
	for (iexpt=0;iexpt<nexpt;iexpt++)
	{
		//setting the pointer for the number of used data points when a new type of data set comes
		if (iexpt==0)
			npoints=ExptsData::sqused;
		if (iexpt==ExptsData::nsq)
			npoints=ExptsData::fqused;
		if (iexpt == (ExptsData::nsq + ExptsData::nfq))
			npoints = ExptsData::fgused;
		
		if (*npoints>nthreads)
		{
			step=int(*npoints/nthreads);//number of data points to handle for each thread
			remainder=(*npoints % nthreads);//reaminder of data points to divide among the threads one by one among the 
			offset=step;//only assign step atoms for the original thread, regardless of remainder
			if (remainder>0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads=nthreads;
		}
		else
		{
			step=1;
			offset=step;
			remainder=0;
			nused_threads=*npoints;

		}
		//assign for ithread=0
		thread_arg[0].min_Q_g_index[iexpt]=0;
		thread_arg[0].max_Q_g_index[iexpt]=offset-1;
		if (nthreads>1)
		{
			for (ithread=1;ithread<nused_threads-1;ithread++)
			{
				thread_arg[ithread].min_Q_g_index[iexpt]=thread_arg[ithread-1].max_Q_g_index[iexpt]+1;
				offset=step;
				if (remainder>0)
				{
					//to spread the remainder
					offset=step+1;
					remainder--;
				} 
				thread_arg[ithread].max_Q_g_index[iexpt]=thread_arg[ithread].min_Q_g_index[iexpt]+offset-1;
			}
			thread_arg[nused_threads-1].min_Q_g_index[iexpt]=thread_arg[nused_threads-2].max_Q_g_index[iexpt]+1;
			thread_arg[nused_threads-1].max_Q_g_index[iexpt]=*npoints-1;
			for (ithread=nused_threads;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].min_Q_g_index[iexpt]=-1;//not to do anything
				thread_arg[ithread].max_Q_g_index[iexpt]=-2;

			}
		}
		npoints++;//go to the number of points of the next data set
	}//end of iexpt cycle

	//Making the split based on the k data points (E(k) fitting)
	//Each thread will handle a consecutive segment of the k point range
	npoints=ExptsData::ekused;
	for (iexpt=0;iexpt<ExptsData::nek; iexpt++)
	{
		if (*npoints>nthreads)
		{
			step=int(*npoints/nthreads);//number of k points to handle for each thread
			remainder=(*npoints % nthreads);//reaminder of k points to divide among the threads one by one among the 
			offset=step;//only assign step atoms for the original thread, regardless of remainder
			if (remainder>0)
			{
				//to spread the remainder
				offset++;
				remainder--;
				nused_threads=nthreads;
			}
		}
		else
		{
			step=1;
			offset=step;
			remainder=0;
			nused_threads=*npoints;

		}
		
		//assign for ithread=0
		thread_arg[0].min_k_index[iexpt]=0;
		thread_arg[0].max_k_index[iexpt]=offset-1;
		if (nthreads>1)
		{
			for (ithread=1;ithread<nthreads-1;ithread++)
			{
				thread_arg[ithread].min_k_index[iexpt]=thread_arg[ithread-1].max_k_index[iexpt]+1;
				offset=step;
				if (remainder>0)
				{
					//to spread the remainder
					offset=step+1;
					remainder--;
				} 
				thread_arg[ithread].max_k_index[iexpt]=thread_arg[ithread].min_k_index[iexpt]+offset-1;
			}
			thread_arg[nthreads-1].min_k_index[iexpt]=thread_arg[nthreads-2].max_k_index[iexpt]+1;
			thread_arg[nthreads-1].max_k_index[iexpt]=*npoints-1;
			for (ithread=nused_threads;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].min_k_index[iexpt]=-1;//not to do anything
				thread_arg[ithread].max_k_index[iexpt]=-2;

			}

		}
		npoints++;//go to the number of points of the next data set
	}//end of iexpt cycle

	//setting the beginning of the non-bonded potential for the threads
	//Paddinng will be used, to make sure that the segments of different threads will not be cached together.
	offset=(RunParams::potential>0 ? (int)((CACHE_PADDING-sizeof(*FNC_POT::vdW_pot))/sizeof(*FNC_POT::vdW_pot)) : 0);//number of dummy double elements
	thread_arg[nthreads - 1].vdW_pot_start=FNC_POT::vdW_pot;// for the main thread
	thread_arg[nthreads - 1].Coulomb_pot_start=FNC_POT::Coulomb_pot;// for the main thread
	thread_arg[nthreads - 1].vdW_pot_s_start = FNC_POT::vdW_pot_s;// for the main thread
	thread_arg[nthreads - 1].Coulomb_pot_s_start = FNC_POT::Coulomb_pot_s;// for the main thread
	thread_arg[nthreads - 1].vdW_pot_l_start = FNC_POT::vdW_pot_l;// for the main thread
	thread_arg[nthreads - 1].Coulomb_pot_l_start = FNC_POT::Coulomb_pot_l;// for the main thread
	thread_arg[nthreads - 1].vdW14_pot_start=FNC_POT::vdW14_pot;// for the main thread
	thread_arg[nthreads - 1].Coulomb14_pot_start=FNC_POT::Coulomb14_pot;// for the main thread
	for (ithread=0; ithread<nthreads-1;ithread++)//for the auxiliary threads
	{
		thread_arg[ithread].vdW_pot_start=FNC_POT::vdW_pot+(ithread+1)*(FNC_POT::pot_dim+offset);
		thread_arg[ithread].Coulomb_pot_start=FNC_POT::Coulomb_pot+(ithread+1)*(FNC_POT::pot_dim +offset);
		thread_arg[ithread].vdW_pot_s_start = FNC_POT::vdW_pot_s + (ithread + 1) * (FNC_POT::pot_dim + offset);
		thread_arg[ithread].Coulomb_pot_s_start = FNC_POT::Coulomb_pot_s + (ithread + 1) * (FNC_POT::pot_dim + offset);
		thread_arg[ithread].vdW_pot_l_start = FNC_POT::vdW_pot_l + (ithread + 1) * (FNC_POT::pot_dim + offset);
		thread_arg[ithread].Coulomb_pot_l_start = FNC_POT::Coulomb_pot_l + (ithread + 1) * (FNC_POT::pot_dim + offset);
		thread_arg[ithread].vdW14_pot_start=FNC_POT::vdW14_pot+(ithread+1)*(FNC_POT::pot_dim +offset);
		thread_arg[ithread].Coulomb14_pot_start=FNC_POT::Coulomb14_pot+(ithread+1)*(FNC_POT::pot_dim +offset);
	}

	//splitting the work during the potential correction due to the charge group centre movement among the threads
	//the charge group centres has to be splitted
	if (Topology::nRMC_charge_centre > 0)
	{
		if (Topology::nRMC_charge_centre > nthreads)
		{
			step = int(Topology::nRMC_charge_centre / nthreads);//number of charge group centres to handle for each thread
			remainder = (Topology::nRMC_charge_centre % nthreads);//remainder of charge group centres to divide among the threads one by one among the 
			offset = step;//only assign step centres for the original thread, regardless of remainder
			if (remainder > 0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads = nthreads;//all the threads will be used
		}
		else
		{
			//there is not enough charge centres for each thread
			step = 1;
			offset = 1;
			remainder = 0;
			nused_threads = Topology::nRMC_charge_centre;
		}
		//assign for ithread=0
		thread_arg[0].charge_gr_min_ind = 0;
		thread_arg[0].charge_gr_max_ind = offset - 1;
		if (nused_threads>1)
		{
			for (ithread = 1; ithread < nused_threads - 1; ithread++)
			{
				thread_arg[ithread].charge_gr_min_ind = thread_arg[ithread - 1].charge_gr_max_ind + 1;
				offset = step;
				if (remainder > 0)
				{
					//to spread the remainder
					offset = step + 1;
					remainder--;
				}
				thread_arg[ithread].charge_gr_max_ind = thread_arg[ithread].charge_gr_min_ind + offset - 1;
			}
			thread_arg[nused_threads - 1].charge_gr_min_ind = thread_arg[nused_threads - 2].charge_gr_max_ind + 1;
			thread_arg[nused_threads - 1].charge_gr_max_ind = Topology::nRMC_charge_centre - 1;
			for (ithread = nused_threads; ithread < nthreads; ithread++)
			{
				thread_arg[ithread].charge_gr_min_ind = -1;//not to do anything
				thread_arg[ithread].charge_gr_max_ind = -2;

			}
		}
	}

#ifdef _USE_LOCAL_INV
	int tot_moved;

	if (RunParams::swap_fraction>0)
		tot_moved=RunParams::nmoved+1;//to have place for the extra atom the (nmoved-1)-th is swapped with
	else
		tot_moved=RunParams::nmoved;//only nmoved atom
	//Paddinng will be used, to make sure that the segments of different threads will not be cached together.
	offset=(tot_moved*histnew_p->sum_bins>0 ? (int)((CACHE_PADDING-sizeof(*histnew_p->local_count))/sizeof(*histnew_p->local_count)) : 0);//number of dummy integer elements
	//setting the pointer to the thread's local_histogram count
	thread_arg[nthreads-1].loc_hist_count=histnew_p->local_count;//main thread will work on the local count array
	for (ithread=0; ithread<nthreads-1; ithread++)
		thread_arg[ithread].loc_hist_count=histnew_p->dlocal_count+ithread*(tot_moved*HistoSet::sum_bins+offset);
	
	//Making the split based on the local histogram bins 
	//Each thread will handle a consecutive segment of bins	
	if (RunParams::max_loc_nbins>nthreads)
	{
		step=RunParams::max_loc_nbins/nthreads;
		remainder=(RunParams::max_loc_nbins % nthreads);
		offset=step;//only assign step atoms for the original thread, regardless of remainder
		if (remainder>0)
		{
			//to spread the remainder
			offset++;
			remainder--;
		}
		nused_threads=nthreads;
	}
	else
	{
		step=1;
		offset=step;
		remainder=0;
		nused_threads=RunParams::max_loc_nbins;

	}
	//assign for ithread=0
	thread_arg[0].loc_min_bin=0;
	thread_arg[0].loc_max_bin=offset-1;
	if (nused_threads>1)
	{
		for (ithread=1;ithread<nused_threads-1;ithread++)
		{
			thread_arg[ithread].loc_min_bin=thread_arg[ithread-1].loc_max_bin+1;
			offset=step;
			if (remainder>0)
			{
				//to spread the remainder
				offset=step+1;
				remainder--;
			} 
			thread_arg[ithread].loc_max_bin=thread_arg[ithread].loc_min_bin+offset-1;
		}
		thread_arg[nused_threads-1].loc_min_bin=thread_arg[nused_threads-2].loc_max_bin+1;
		thread_arg[nused_threads-1].loc_max_bin=RunParams::max_loc_nbins-1;
		for (ithread=nused_threads;ithread<nthreads;ithread++)
		{
			thread_arg[ithread].loc_min_bin=-1;//not to do anything
			thread_arg[ithread].loc_max_bin=-2;

		}

	}
#endif
#ifdef _LOCAL_INV
	if (RunParams::loc_chi2_mode==1)
	{
		//If there are more threads
		
		for (itype=0;itype<RunParams::ntypes;itype++)
		{
			//we have to assign ChiSquared::loc_natoms[itype] atoms equally 
			if (ChiSquared::loc_natoms[itype]>nthreads)
			{
				step=(int)(ChiSquared::loc_natoms[itype]/nthreads);
				remainder=(ChiSquared::loc_natoms[itype] % nthreads);
				offset=step;//only assign step atoms for the original thread, regardless of remainder
				if (remainder>0)
				{
					//to spread the remainder
					offset++;
					remainder--;
				}
				nused_threads=nthreads;
			}
			else
			{
				if (ChiSquared::loc_natoms[itype]>0)
					nused_threads=ChiSquared::loc_natoms[itype];
				else
					nused_threads=0;
				step=1;
				offset=step;
				remainder=0;
			}
			
			//assign for ithread=0
			if (nused_threads>0)
			{
				thread_arg[0].min_loc_atom_index[itype]=ChiSquared::loc_natoms_min[itype*(RunParams::nlocint+1)];//the first index will be ChiSquared::loc_natoms_min[itype*(RunParams::nlocint+1)]
				thread_arg[0].max_loc_atom_index[itype]=offset-1+ChiSquared::loc_natoms_min[itype*(RunParams::nlocint+1)];
			}
			else
			{
				thread_arg[0].min_loc_atom_index[itype]=-1;
				thread_arg[0].max_loc_atom_index[itype]=-2;
			}
			if (nused_threads>1)
			{
				for (ithread=1;ithread<nused_threads-1;ithread++)
				{
					thread_arg[ithread].min_loc_atom_index[itype]=thread_arg[ithread-1].max_loc_atom_index[itype]+1;
					offset=step;
					if (remainder>0)
					{
						//to spread the remainder
						offset=step+1;
						remainder--;
					} 
				
					thread_arg[ithread].max_loc_atom_index[itype]=thread_arg[ithread].min_loc_atom_index[itype]+offset-1;
				}
				if (nused_threads>0)
				{
					thread_arg[nused_threads-1].min_loc_atom_index[itype]=thread_arg[nused_threads-2].max_loc_atom_index[itype]+1;
					thread_arg[nused_threads-1].max_loc_atom_index[itype]=ChiSquared::loc_natoms[itype]+ChiSquared::loc_natoms_min[itype*(RunParams::nlocint+1)]-1;//to make sure, that all is assigned, and no more
				}
	
			}
			for (ithread=nused_threads;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].min_loc_atom_index[itype]=-1;//not to do anything
				thread_arg[ithread].max_loc_atom_index[itype]=-2;
			}
			ChiSquared::nused_loc[itype]=nused_threads;
			if (nused_threads>0 && thread_arg[nused_threads-1].max_loc_atom_index[itype]==thread_arg[nused_threads-1].min_loc_atom_index[itype])
				ChiSquared::nused_loc[itype]*=-1;//this indicates, that for the clean partials the previous thread will be the last
		}//end of itype cycle
		
		//setting the loc_hist_offset and the av_loc_hist_finder
offset=0;//=(RunParams::nthreads>1 ? (int)((CACHE_PADDING-sizeof(*ChiSquared::first_loc_ind))/sizeof(*ChiSquared::first_loc_ind)) : 0);//number of dummy double elements
		//cache padding would not be practical, as different segmentation is used for av_loc_hist and chi2 loc dist calculation
		for (ithread=0;ithread<nthreads;ithread++)
			thread_arg[ithread].loc_hist_offset=ithread*SimpleCfg::ntypes*SimpleCfg::ntotal;

		for (ithread=0;ithread<nthreads;ithread++)
			thread_arg[ithread].av_loc_hist_finder=chi_p->av_loc_hist+ithread*SimpleCfg::ntypes*ChiSquared::sum_loc_natoms_tot;

		
	}
	
#endif

};

#ifdef _LOCAL_INV
//resetting the segment boundaries in case of loc_chi2_mode to help load balancing, as equal splitting according to the 
//neighbour atoms can result load inbalancing due to the unequal distribution of the count in the local histogram bins usually
//resulting, that the thread 0 has to work much more, than the others, as here there are fever counts 
	void Threads::ResetSegment()
	{
		int i,iav,k,itype,jtype,temp_ind;
		double step,*my_av_loc;
		k=0;
		for (itype=0;itype<SimpleCfg::ntypes;itype++)
		{
			for (jtype=0;jtype<SimpleCfg::ntypes;jtype++)
			{
				if (ChiSquared::loc_natoms[jtype]>nthreads)//no need to do it, if there is not enough neighbour
				{
					//calculate the difference between the index of the last and first average local histogram bin for this partial
					step=ChiSquared::av_loc_hist[ChiSquared::av_loc_hist_offset[k+1]-1]-ChiSquared::av_loc_hist[ChiSquared::av_loc_hist_offset[k]];
					if (step>=nthreads)
						step/=nthreads;
					else
						continue;//no need to balance, just a few atoms
					iav=0;
					my_av_loc=ChiSquared::av_loc_hist+ChiSquared::av_loc_hist_offset[k];
					//the 0-th thread min_loc_atom and the last thread max_loc_atom does not change
					for (i=0;i<nthreads-1;i++)
					{
						temp_ind=(int)(*my_av_loc+(i+1)*step);
						//now find the closest larger average index
						while (my_av_loc[iav]<temp_ind)
							iav++;
						thread_arg[i].max_loc_atom_index[jtype]=iav;
						thread_arg[i+1].min_loc_atom_index[jtype]=iav+1;
						
					}
				
				}
				k++;
			}
		}
	};
#endif

//-------this is called at thread creation------------------------
void Threads::ThreadEntry(ThreadArg &thread_arg)
{	
	int i;

	//ThreadArg *p_arg=(struct ThreadArg*) entry_arg;
	//ThreadArg &thread_arg=*p_arg;//to create a reference for the thread arguments
#ifdef _NO_PERIODIC
	DataMat &datmat=*datmat_p;
#endif
#ifdef _TEST_MODE
	int index=thread_arg.thread_index;
#endif
	
	//Calculating the histogram, if necessary
	for (i=0;i<calc_hist_flag;i++)
	{
	  {
		std::unique_lock<std::mutex> guard(start1_mutex);//lock mutex
		if (start_flag) check_hist_start.notify_all();
		else check_hist_start.wait(guard);
	  }
#ifdef _TEST_MODE
		thtime1[index] = std::chrono::high_resolution_clock::now();
#endif
		histnew_p->HistCalcThread(thread_arg);//calculate histogram for the given configuration

#ifdef _TEST_MODE
		thtime2[index] = std::chrono::high_resolution_clock::now();
		elapsed[index] = thtime2[index]-thtime1[index];
		dur_10[index]+=elapsed[index].count();
#endif

		//to wait for all the threads to finish, and signal to continue, if this was the last one
		std::unique_lock<std::mutex> guard(count_mutex);//lock mutex
		complete_count++;
		if (complete_count==nthreads)//time to proceed
		{
			start_flag=0;//reset it to 0
			complete_count=0;//reset it to 0
			check_count.notify_all();
		}
		else check_count.wait(guard);
	}	
	
#ifdef _NO_PERIODIC	
	//wait, till the signal comes to proceed with execution
	{
	  std::unique_lock<std::mutex> guard(start2_mutex);//lock mutex
	  if (start_flag) check_cd_start.notify_all();
	  else check_cd_start.wait(guard);
	}

	if (datmat.calcmode==0)
		histnew_p->CalcPerHist(thread_arg,datmat);//calculate the peiodic histogram, if necessary	

	//Threads have to wait, till all of them finished
	{
	  std::unique_lock<std::mutex> guard(NOP_count_mutex);//lock mutex
	  NOP_count++;
	  if (NOP_count==nthreads)//time to proceed
	  {
		start_flag=0;//reset it to 0
		NOP_count=0;//reset it to 0
		NOP_check_count.notify_all();
	  }
	  else  NOP_check_count.wait(guard);
	}
#endif //_NO_PERIODIC
#ifdef _TEST_MODE
	thtime1[index] = std::chrono::high_resolution_clock::now();
	elapsed[index] = thtime1[index]-thtime2[index];
	dur_11[index]+=elapsed[index].count();		
#endif

	//PPCF calculation
	//wait, till the signal comes to proceed with execution
	{
		std::unique_lock<std::mutex> guard(start2_mutex);//lock mutex
		if (start_flag) check_cd_start.notify_all();
		else check_cd_start.wait(guard);
	}

#ifdef _AENET
	//calculating the initial ANN potential
#ifdef _TEST_MODE
	anntime1[index] = std::chrono::high_resolution_clock::now();
#endif
	if (calc_ANN)
		aenetnew_p->CalcANN(thread_arg);
#ifdef _TEST_MODE
	anntime2[index] = std::chrono::high_resolution_clock::now();
	elapsed[index]=anntime2[index]-anntime1[index]; 
	dur_ann1[index]=elapsed[index].count();
#endif
	
#endif
	//Each thread handles a different histogram bin-segment	
	ppcfnew_p->CalcPPCF(thread_arg);	

	//Threads have to wait, till all of them finished
	{
	  std::unique_lock<std::mutex> guard(count_mutex);//lock mutex
	  complete_count2++;
	  if (complete_count2==nthreads)//time to proceed
	  {
		start_flag=0;//reset it to 0
		complete_count2=0;//reset it to 0
		check_count.notify_all();
	  }
	  else  check_count.wait(guard);
	}

	//calculating the initial calculated data
	//Each thread handles a different segment	
	if (RunParams::ngr>0)//g(r) fitting
		calcpnew_p->CalcPartialgr(thread_arg);//converts the ppcfs from bins to r data point

	if (RunParams::nsq>0)
		calcpnew_p->CalcPartialSq(thread_arg);//calculating the partial S(Q)-s by Fourier transformation
		
	if (RunParams::nfq>0)
		calcpnew_p->CalcPartialFq(thread_arg);//calculating the partial F(Q)-s by Fourier transformation

	if (RunParams::nfg > 0)
		calcpnew_p->CalcPartialFg(thread_arg);//calculating the partial F(g)-s by Fourier transformation

	if (RunParams::nek>0)
		calcpnew_p->CalcPartialEk(thread_arg);//calculating the partial E(k)-s by Fourier transformation

	calcdat_p->CalcTotal(thread_arg);//calculate the total from the partials
	
	//to wait for all the threads to finish, and signal to continue, if this was the last one
	{
	  std::unique_lock<std::mutex> guard(T_count_mutex);//lock mutex
	  T_count++;
	  if (T_count==nthreads)//time to proceed
	  {
		  start_flag=0;//reset it to 0
		  T_count=0;//reset it to 0
		  T_check_count.notify_all();
	  }
	  else T_check_count.wait(guard);
	}
	

#ifdef _LOCAL_INV
	//wait, till the signal comes to proceed with execution
	if (RunParams::loc_chi2_mode==1)
	{
		//first calculate the average local histogram
		{
		  std::unique_lock<std::mutex> guard(inv_mutex);//lock mutex
		  if (start_flag) check_inv_start.notify_all();
		  else check_inv_start.wait(guard);
		}
		//calculate the chisqure contribution coming from this thread
		chi_p->CalcLocAvThread_dist(thread_arg);

		//Threads have to wait, till all of them finished, as each of them will use all the bins
		std::unique_lock<std::mutex> guard(LOC_count_mutex);//lock mutex
		LOC_count++;
		if (LOC_count==nthreads)//time to proceed
		{
			start_flag=0;//reset it to 0
			LOC_count=0;//reset it to 0
			LOC_check_count.notify_all();
		}
		else LOC_check_count.wait(guard);
	}

	{
	  std::unique_lock<std::mutex> guard(inv_mutex);//lock mutex
	  if (start_flag) check_inv_start.notify_all();
	  else check_inv_start.wait(guard);
	}

	//calculate the chisqure contribution coming from this thread
	(chi_p->*CalcLocChi2Thread)(thread_arg);

	//Threads have to wait, till all of them finished, as each of them will use all the bins
	{
	  std::unique_lock<std::mutex> guard(LOC_count_mutex2);//lock mutex
	  LOC_count2++;
	  if (LOC_count2==nthreads)//time to proceed
	  {
		  start_flag=0;//reset it to 0
		  LOC_count2=0;//reset it to 0
		  LOC_check_count2.notify_all();
	  }
	  else LOC_check_count2.wait(guard);
	}
#endif
	
	//Wait, till the main creates the second PPCF and CalcPart and initialize the pointers, chi2 is calculated, History is created...
	{
	  std::unique_lock<std::mutex> guard(start3_mutex);//lock mutex
	  if (start_flag) check_loop_start.notify_all();
	  else check_loop_start.wait(guard);
	}

	//Here starts the main loop
	ThreadLoop(thread_arg);	
	
#ifdef _AENET

	if (RunParams::aenet_step>1)
	{ 
		//Wait, till the loop finishes, and the split is made, if no calculatiom is necessary thread_arg.aenet_at_count=0 
		{
			std::unique_lock<std::mutex> guard(aenet_mutex);//lock mutex
			if (aenetstart_flag) check_loop_start.notify_all();
			else check_aenet_start.wait(guard);
		}
		ThreadAenetFinal(thread_arg);
	}
#endif
	return;
}

//-------executing the main loop for the threads, it is called for each thread separately---------
void Threads::ThreadLoop(ThreadArg &my_arg)
{
	int i,itype;
	int break_flag=0;//this indicates, if the while loop was quit through 
	int ntypes=RunParams::ntypes;
	int npartials=ntypes*(ntypes+1)/2;
	int partialind,*p_type1;
	int my_ngen = 0;

#ifdef _TEST_MODE
int index;
index=my_arg.thread_index;
#endif

	HistoSet &histnew=*histnew_p;
	HistoSet &histold=*histold_p;
	CoordNumbConst &iccnew=*iccnew_p;
	CoordNumbConst &iccold=*iccold_p;

	Move &move=*move_p;
	
	PPCFSet &ppcfnew=*ppcfnew_p;
	PPCFSet &ppcfold=*ppcfold_p;
	CalcPart &calcpnew=*calcpnew_p;
	CalcPart &calcpold=*calcpold_p;
	CalcData &calcdat=*calcdat_p;
#ifdef _LOCAL_INV
	ChiSquared &chi=*chi_p;
#endif
#ifdef _NO_PERIODIC
	DataMat &datmat=*datmat_p;
#endif
#ifdef _AENET
	Aenet &aenetnew=*aenetnew_p;
	Aenet &aenetold=*aenetold_p;
#endif
#ifdef _TH_ORDER
	ofstream *plog;
	if (RunParams::nthreads>1)
	{
		switch (my_arg.thread_index)
		{
		case 0:
			plog=&logfile0;
			break;
		case 1:
			plog=&logfile1;
			break;
		case 2:
			plog=&logfile2;
			break;
		default:
			plog = &logfile;
		}
	}
#endif
	while (loop_flag)
	{
		//wait, till all the threads are ready to go
		{
			std::unique_lock<std::mutex> guard(count_mutex);//lock mutex
			start_count++;
			my_ngen++;
#ifdef _TH_ORDER
			*plog << my_ngen << "HCC "<<start_count<<endl;
#endif
			if (start_count == nthreads)//time to proceed, this was the last to reach this point
			{
#ifdef _TH_ORDER
				switch (index)
				{
				case 0:
					cout << "a" << endl;
					break;
				case 1:
					cout << "A" << endl;
					break;
				case 2:
					cout << "b" << endl;
					break;
				}
#endif
				acceptable = 1;//the move is acceptable, if it does not prove otherwise
				if (Move::nomove)//no atomic move is made
				{
					if (Move::makeE0shift)
						is_E0_shift = 1;//this is needed that the actual value should be preserved for the duration of the thread loop
					else
						is_E0_shift = 0;
					if (Move::makemucorr)
						is_IQ_mucorr = 1;
					else
						is_IQ_mucorr = 0;
									
					if (Move::makefprimeshift)
						is_fprime_shift = 1;
					else
						is_fprime_shift = 0;
#ifdef _TH_ORDER
					*plog << "c " << is_IQ_mucorr << endl;
#endif
#ifdef _AENET
					//if ANN is not calculated in each step, then if there were accepted moves without ANN
					//calculation since the last ANN calculation, then has to calculate ANN
					if (Aenet::last_gen_accepted!=Aenet::last_gen_aenet_accepted && Aenet::calc_ANN==2)
						calc_ANN=Aenet::calc_ANN;
					else
						calc_ANN=0;//no need to calculate if ANN is calculated in each step
#endif
					if (!RunParams::custmove && RunParams::ntypes > 1)//for custmove or ntypes=1 no update is necessary, as all the partials and types are modified
					{
						for (i = 0; i < Move::npartials; i++)
							mod_partial[i] = 0;

						if (Move::makeE0shift)//technically only those partials are modified, which contain the edge particle, but the others are not used
						{
							for (i = 0; i < Move::npartials; i++)
								mod_partial[i] = 1;
						}
						//mucorr and f' shift does not alter the partial F(Q)-s
					}
					
				}
				else //normal step with atomic moves
				{
					is_E0_shift = 0;
					is_IQ_mucorr = 0;
					is_fprime_shift = 0;
#ifdef _AENET
					calc_ANN=Aenet::calc_ANN;
					if (calc_ANN==1)//calculated in every step, use aenet_changed_neigh_ind
						thread_p->SplitAenetAtoms(Aenet::nchanged,Aenet::aenet_changed_neigh_ind);//split the neighbour atoms for which the ANN pot should be calculated among the threads
					else if (calc_ANN==2)
						thread_p->SplitAenetAtoms(Aenet::nstored,Aenet::aenet_stored_ind);//use the stored indices, al the atoms involved in change since last calc
			
#endif
					if (!RunParams::custmove && RunParams::ntypes > 1)//for custmove or ntypes=1 no update is necessary, as all the partials and types are modified
					{

						for (i = 0; i < npartials; i++)
							mod_partial[i] = 0;

						p_type1 = move.types;//points to the type of the first moved atom
						for (i = 0; i < Move::tot_moved_atoms; i++)
						{
							//Determining the index of the partial 
							for (itype = 0; itype < ntypes; itype++)//go through all the partials containing this type
							{
								partialind = (*p_type1 <= itype ? (*p_type1 * ntypes - (*p_type1 * ((*p_type1) + 1) / 2) + itype) : \
									(itype * ntypes - (itype * (itype + 1) / 2) + *p_type1));
								mod_partial[partialind] = 1;//this partial was modified
							}
							p_type1++;
						}
					}
				}
				for (i = 0; i < ExptsData::nek; i++)
					EXAFS_gridind[i] = ExptsData::ek_gridind[i];
				for (i = 0; i < ExptsData::nfq; i++)
					IQ_muact[i] = ExptsData::fqmuact[i];
				for (i = 0; i < ExptsData::nfq; i++)
					IQ_fprimeact[i] = ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1];
				start_flag = 0;//reset it for later use
				start_count = 0;//all the threads finished calculation, reset it to 0
												
				check_count.notify_all();
#ifdef _TH_ORDER
				cout << "<" << endl;
#endif
			}
			else
				check_count.wait(guard);
		}
		
		
		if (!loop_flag)
		{
#ifdef _TH_ORDER
			*plog << "b b H" << endl;
#endif
			break_flag=1;//the loop was abandoned through break
			break;//to exit loop, if time is up, as most probably the thread is in conditional wait already
		}
#ifdef _TEST_MODE
		thtime1[index] = std::chrono::high_resolution_clock::now();
#endif
		if (!(is_E0_shift || is_IQ_mucorr || is_fprime_shift))
		{
			//Calculate the change in the histogram
			histnew.HistCalcChangeThread(my_arg);
		}
		
#ifdef _TEST_MODE
		thtime2[index] = std::chrono::high_resolution_clock::now();
		elapsed[index] = thtime2[index]-thtime1[index];
		dur_10[index]+=elapsed[index].count();
#endif

#ifdef _AENET
#ifdef _TEST_MODE
		anntime1[index] = std::chrono::high_resolution_clock::now();
#endif
		if (calc_ANN)
			aenetnew.CalcANN(my_arg);
#ifdef _TEST_MODE
		anntime2[index] = std::chrono::high_resolution_clock::now();
		elapsed[index]=anntime2[index]-anntime1[index]; 
		dur_ann2[index]+=elapsed[index].count();
#endif
#endif

		//wait, till all the threads finished calculation
		{
			std::unique_lock<std::mutex> guard(CC_count_mutex);//lock mutex
			CC_count++;
#ifdef _TH_ORDER
			*plog << "CC " << CC_count << endl;
			switch (index)
			{
			case 0:
				cout << "c" << endl;
				break;
			case 1:
				cout << "d" << endl;
				break;
			case 2:
				cout << "e" << endl;
				break;
			}
#endif
			if (CC_count==nthreads)//time to proceed
			{
				start_flag=0;//reset it to 0
				CC_count=0;//all the threads finished calculation, reset it to 0
				CC_check_count.notify_all();
#ifdef _TH_ORDER
				cout << ">" << endl;
#endif
			}
			else
			{
				CC_check_count.wait(guard);
#ifdef _TH_ORDER
				cout << "." << endl;
#endif
			}
		}

#ifdef _TEST_MODE
		thtime1[index] = std::chrono::high_resolution_clock::now();
		elapsed[index] = thtime1[index]-thtime2[index];
		dur_11[index]+=elapsed[index].count();
#endif			
		if (!(is_E0_shift || is_IQ_mucorr || is_fprime_shift))//normal move
		{
			//Updating the coordnumbs array and calculating the number of atoms satisfying the coordination constraint
			if (iccnew.nconstraints > 0)//there is/are coordination constraints
				iccnew.UpdateCoordnumb(my_arg);
#ifdef _ADVANCED_GEOM_CONST
			//updating the valence array
			if (bvsnew_p->nconstraints > 0)
				bvsnew_p->UpdateValence(my_arg);
#endif
			if (iccnew.nconstraints > 0 || RunParams::nbvs > 0)
			{
#ifdef _TEST_MODE
				thtime2[index] = std::chrono::high_resolution_clock::now();
				elapsed[index] = thtime2[index]-thtime1[index];
				dur_12[index]+=elapsed[index].count();
#endif

				//wait, till all the threads finished calculation
				std::unique_lock<std::mutex> guard(UC_count_mutex);//lock mutex
				UC_count++;
#ifdef _TH_ORDER
				*plog << "UC " << UC_count << endl;
				switch (index)
				{
				case 0:
					cout << "F" << endl;
					break;
				case 1:
					cout << "f" << endl;
					break;
				case 2:
					cout << "g" << endl;
					break;
				}
#endif
				if (UC_count == nthreads)//time to proceed
				{
					start_flag = 0;//reset it to 0
					UC_count = 0;//all the threads finished calculation, reset it to 0
					UC_check_count.notify_all();
#ifdef _TH_ORDER
					cout << "#" << endl;
#endif
				}
				else 
					UC_check_count.wait(guard);
			};

#ifdef _TEST_MODE
			thtime1[index] = std::chrono::high_resolution_clock::now();
			elapsed[index] = thtime1[index]-thtime2[index];
			dur_12a[index]+=elapsed[index].count();
#endif

			//wait, till the signal comes to proceed with execution
			{
				std::unique_lock<std::mutex> guard(start2_mutex);//lock mutex
				if (start_flag)
				{
#ifdef _TH_ORDER
					*plog << "WACC s1 " << acceptable<<endl;
					switch (index)
					{
					case 0:
						cout << "B" << endl;
						break;
					case 1:
						cout << "C" << endl;
						break;
					case 2:
						cout << "D" << endl;
						break;
					}
					check_cd_start.notify_all();
					cout << "(" << endl;
#endif
				}
				else
				{
#ifdef _TH_ORDER
					*plog << "WACC s0 " << endl;
					switch (index)
					{
					case 0:
						cout << "G" << endl;
						break;
					case 1:
						cout << "H" << endl;
						break;
					case 2:
						cout << "I" << endl;
						break;
					}
#endif
					check_cd_start.wait(guard);
#ifdef _TH_ORDER
					cout << ")" << endl;
#endif
				}
			}

			//we have to check, whether the move is still acceptable on the bases of the potentia-based chi2
			if (!acceptable)
			{
#ifdef _TH_ORDER
				*plog << "BP!" << endl;
				switch (index)
				{
				case 0:
					cout << "J" << endl;
					break;
				case 1:
					cout << "K" << endl;
					break;
				case 2:
					cout << "L" << endl;
					break;
				}
#endif
				iccold.CopyModified(iccnew, my_arg);//Copy the modified coordination numbers 
#ifdef _ADVANCED_GEOM_CONST
				bvsold_p->CopyModified(*bvsnew_p, my_arg);//Copy the modified valences
#endif
				histold.CopyModified(histnew, my_arg);//Copying back the unmodified histogram parts to histnew
#ifdef _AENET
				if (calc_ANN)
					aenetold.Copy(aenetnew,move.indices,my_arg);//Copy the modified energis back
#endif	
#ifdef _TEST_MODE
				thtime2[index] = std::chrono::high_resolution_clock::now();
				elapsed[index] = thtime2[index]-thtime1[index];
				dur_12b[index]+=elapsed[index].count();
#endif
				continue;//start the loop again
			}//if not acc
#ifdef _TEST_MODE
			thtime1[index] = std::chrono::high_resolution_clock::now();
#endif

#ifdef _NO_PERIODIC
			if (datmat.calcmode == 0)
				histnew.CalcPerHist(my_arg, datmat);//calculate the peiodic histogram, if necessary	
			//Threads have to wait, till all of them finished, as each of them will use all the bins
			{
				std::unique_lock<std::mutex> guard(NOP_count_mutex);//lock mutex
				NOP_count++;
				if (NOP_count == nthreads)//time to proceed
				{
					start_flag = 0;//reset it to 0
					NOP_count = 0;//reset it to 0
					NOP_check_count.notify_all();
				}
				else NOP_check_count.wait(guard);
			}


			//wait, till the signal comes to proceed with execution, as the vibrational convp�ution has to be done by main thread
			{
				std::unique_lock<std::mutex> guard(start2_mutex);//lock mutex
				if (start_flag)
				{
					check_cd_start.notify_all();
					//cout << "{" << endl;
				}
				else 
					check_cd_start.wait(guard);
			}
#endif
			//Each thread handles bin segment 
			//calculate the modidified parts of the ppcf
			ppcfnew.CalcModPPCF(my_arg);

#ifdef _TEST_MODE
			thtime2[index] = std::chrono::high_resolution_clock::now();
			elapsed[index] = thtime2[index]-thtime1[index];
			dur_13[index]+=elapsed[index].count();
#endif

			//Threads have to wait, till all of them finished, as each of them will use all the bins
			{
				std::unique_lock<std::mutex> guard(CMP_count_mutex);//lock mutex
				CMP_count++;
#ifdef _TH_ORDER
				*plog << "CMP " << CMP_count << endl;
				switch (index)
				{
				case 0:
					cout << "h" << endl;
					break;
				case 1:
					cout << "i" << endl;
					break;
				case 2:
					cout << "j" << endl;
					break;
				}
#endif
				if (CMP_count == nthreads)//time to proceed
				{
					start_flag = 0;//reset it to 0
					CMP_count = 0;//reset it to 0
					CMP_check_count.notify_all();
#ifdef _TH_ORDER
					cout << "&" << endl;
#endif
				}
				else 
					CMP_check_count.wait(guard);
			}

#ifdef _TEST_MODE
			thtime1[index] = std::chrono::high_resolution_clock::now();
			elapsed[index] = thtime1[index]-thtime2[index];
			dur_13a[index]+=elapsed[index].count();
#endif
		}
		//Each thread handles a different r segment for g(r) fitting and Q segment for S(Q), F(Q), g segment for F(g) and E(k) fitting
		//calculate the modidified partial g(r), S(Q), F(Q), F(g) and  E(k) if necessary
		calcpnew.CalcModPartial(my_arg);

#ifdef _TEST_MODE
		thtime2[index] = std::chrono::high_resolution_clock::now();
		elapsed[index] = thtime2[index]-thtime1[index];
		dur_13b[index]+=elapsed[index].count();
#endif
		
		//calculate the totals
		calcdat.CalcTotal(my_arg);

#ifdef _TEST_MODE
		thtime1[index] = std::chrono::high_resolution_clock::now();
		elapsed[index] = thtime1[index]-thtime2[index];
		dur_14[index]+=elapsed[index].count();
#endif

		//to wait for all the threads to finish, and signal to continue, if this was the last one
		{
			std::unique_lock<std::mutex> guard(T_count_mutex);//lock mutex
			T_count++;
#ifdef _TH_ORDER
			*plog << "T " << T_count << endl;
			switch (index)
			{
			case 0:
				cout << "k" << endl;
				break;
			case 1:
				cout << "l" << endl;
				break;
			case 2:
				cout << "m" << endl;
				break;
			}
#endif
			if (T_count == nthreads)//time to proceed
			{
				start_flag = 0;//reset it to 0
				T_count = 0;//reset it to 0
				T_check_count.notify_all();
#ifdef _TH_ORDER
				cout << "@" << endl;
#endif
			}
			else 
				T_check_count.wait(guard);
		}
#ifdef _TEST_MODE
		thtime2[index] = std::chrono::high_resolution_clock::now();
		elapsed[index] = thtime2[index]-thtime1[index];
		dur_15[index]+=elapsed[index].count();
#endif

		if (!(is_E0_shift || is_IQ_mucorr || is_fprime_shift))
		{
#ifdef _LOCAL_INV
			//wait, till the signal comes to proceed with execution
#ifdef _TEST_MODE
			lctime1[index] = std::chrono::high_resolution_clock::now();

#endif	
			if (RunParams::loc_chi2_mode == 1)
			{
				//first the average local histogram has to be calculated
				{
					std::unique_lock<std::mutex> guard(inv_mutex);//lock mutex
					if (start_flag) check_inv_start.notify_all();
					else check_inv_start.wait(guard);
				}

				//calculate the chisqure contribution coming from this thread
				chi.CalcLocAvThread_dist(my_arg);

				//Threads have to wait, till all of them finished, as each of them will use all the bins
				std::unique_lock<std::mutex> guard(LOC_count_mutex);//lock mutex
				LOC_count++;
				if (LOC_count == nthreads)//time to proceed
				{
					start_flag = 0;//reset it to 0
					LOC_count = 0;//reset it to 0
					LOC_count2 = 0;//reset it to 0
					LOC_check_count.notify_all();

				}
				else LOC_check_count.wait(guard);
			}
#ifdef _TEST_MODE
			lctime2[index] = std::chrono::high_resolution_clock::now();
			elapsed[index] = lctime2[index]-lctime1[index];
			dur_chiloc1[index]+=elapsed[index].count();
#endif
			{
				std::unique_lock<std::mutex> guard(inv_mutex);//lock mutex
				if (start_flag) check_inv_start.notify_all();
				else check_inv_start.wait(guard);
			}

				//calculate the chisqure contribution coming from this thread
			(chi.*CalcLocChi2Thread)(my_arg);

#ifdef _TEST_MODE
			lctime1[index] = std::chrono::high_resolution_clock::now();
			elapsed[index] = lctime1[index]-lctime2[index];
			dur_chiloc3[index]+=elapsed[index].count();
#endif

			//Threads have to wait, till all of them finished, as each of them will use all the bins
			{
				std::unique_lock<std::mutex> guard(LOC_count_mutex2);//lock mutex
				LOC_count2++;
				if (LOC_count2 == nthreads)//time to proceed
				{
					start_flag = 0;//reset it to 0
					LOC_count = 0;//reset it to 0
					LOC_count2 = 0;//reset it to 0
					LOC_check_count2.notify_all();
				}
				else LOC_check_count2.wait(guard);

			}
#ifdef _TEST_MODE
			lctime2[index] = std::chrono::high_resolution_clock::now();
			elapsed[index] = lctime2[index]-lctime1[index];
			dur_chiloc4[index]+=elapsed[index].count();
#endif

#endif
		}//end of if no shift
		//the main thread decides, whether the move is acceptable

		//copying the data depending on whether the move was acceptable or not
		//wait, till the signal comes to proceed with execution
		{
		  std::unique_lock<std::mutex> guard(start3_mutex);//lock mutex
		  if (start_flag)
		  {
			  check_cp_start.notify_all();
#ifdef _TH_ORDER
			  cout << "}" << endl;
#endif
		  }
		  else 
			  check_cp_start.wait(guard);
		}

		if (acceptable)
		{
			if (!(is_E0_shift || is_IQ_mucorr || is_fprime_shift))
			{
				iccnew.CopyModified(iccold, my_arg);//Copy the modified coordination numbers 
#ifdef _ADVANCED_GEOM_CONST
				bvsnew_p->CopyModified(*bvsold_p, my_arg);//Copy the modified valences
#endif
#ifdef _TEST_MODE
				thtime1[index] = std::chrono::high_resolution_clock::now();
#endif
				
				histnew.CopyModified(histold, my_arg);//Copying the modified histogram parts to histold
#ifdef _TEST_MODE
				thtime2[index] = std::chrono::high_resolution_clock::now();
				elapsed[index] = thtime2[index]-thtime1[index];
				dur_15b[index]+=elapsed[index].count();
#endif
				ppcfnew.CopyModified(ppcfold, my_arg);//Copying the modified ppcf parts to ppcfold
			}
#ifdef _AENET
			if (calc_ANN)
				aenetnew.Copy(aenetold,move.indices,my_arg);//Copy the modified energy to aenetold
#endif

			calcpnew.CopyModified(calcpold,my_arg);//Copying the modified g(r) partial parts to calcpold
			
		}
		else
		{
			if (!(is_E0_shift || is_IQ_mucorr || is_fprime_shift))
			{
				iccold.CopyModified(iccnew, my_arg);//Copy the modified coordination numbers 
#ifdef _ADVANCED_GEOM_CONST
				bvsold_p->CopyModified(*bvsnew_p, my_arg);//Copy the modified valences
#endif
#ifdef _TEST_MODE
				thtime1[index] = std::chrono::high_resolution_clock::now();
#endif
				
				histold.CopyModified(histnew, my_arg);//Copying back the unmodified histogram parts to histnew
#ifdef _TEST_MODE
				thtime2[index] = std::chrono::high_resolution_clock::now();
				elapsed[index] = thtime2[index]-thtime1[index];
				dur_15b[index]+=elapsed[index].count();
#endif

				ppcfold.CopyModified(ppcfnew, my_arg);//Copying back the unmodified ppcf parts to ppcfnew
#ifdef _AENET
				if (calc_ANN)
					aenetold.Copy(aenetnew,move.indices,my_arg);//Copy the old energy to aenetmew
#endif
			}
			calcpold.CopyModified(calcpnew,my_arg);//Copying back the unmodified g(r) partial parts to calcpnew
		}
	}//end of if loop_flag

	if (!break_flag)
	{
		//the loop was quit at the while(loop_flag), so the start_count was not increased for this thread
		std::unique_lock<std::mutex> guard(count_mutex);//lock mutex
		start_count++;
#ifdef _TH_ORDER
		*plog << "EL " << start_count << endl;
		switch (index)
		{
		case 0:
			cout << "n" << endl;
			break;
		case 1:
			cout << "o" << endl;
			break;
		case 2:
			cout << "O" << endl;
			break;
		}
#endif
    
		if (start_count==nthreads)//time to proceed, this was the last to reach this point
		{
			start_count=0;//all the threads finished calculation, reset it to 0
			check_count.notify_all();
#ifdef _TH_ORDER
			cout << "=" << endl;
#endif
		}
		else 
		{
			check_count.wait(guard);
#ifdef _TH_ORDER
			cout << "*" << endl;
#endif
		}
	}
#ifdef _TH_ORDER
	cout<<"_"<<endl;
#endif
};
	

#ifdef _AENET
	void Threads::ResizeAenetArrays(int id,  int **ind_list, int **type_list, double **coord_list)
	{
		ResizeArray(&thread_arg[id].aenet_max_nneigh, thread_arg[id].aenet_max_nneigh+5,&thread_arg[id].aenet_neighlist,"thread_arg.aenet_neighlist","Threads::ResizeAenetArrays");
		thread_arg[id].aenet_max_nneigh-=5;
		ResizeArray(&thread_arg[id].aenet_max_nneigh, thread_arg[id].aenet_max_nneigh+5,&thread_arg[id].aenet_neightype,"thread_arg.aenet_neightype","Threads::ResizeAenetArrays");
		thread_arg[id].aenet_max_nneigh-=5;
		ResizeArray(&thread_arg[id].aenet_max_nneigh, thread_arg[id].aenet_max_nneigh+5,&thread_arg[id].aenet_neighcoord,"thread_arg.aenet_neighcoord","Threads::ResizeAenetArrays",3);
		*ind_list=thread_arg[id].aenet_neighlist;
		*type_list=thread_arg[id].aenet_neightype;
		*coord_list=thread_arg[id].aenet_neighcoord;

	};

	//split the neighbour atoms to calculate E for among the threads 
	void Threads::SplitAenetAtoms(int n, int *ind)
	{
		int ithread,step,remainder,offset,nused_threads;
		//not consecutive, index list is used
		if (n==0)
		{
			for (ithread=0;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].aenet_at_index[0]=-1;//not to do anything
				thread_arg[ithread].aenet_at_count=0;
			}	
			return;
		}
		if (n>nthreads)
		{
			step=(int)n/nthreads;
			remainder=(n% nthreads);
			offset=step;//only assign step atoms for the original thread, regardless of remainder
			if (remainder>0)
			{
				//to spread the remainder
				offset++;
				remainder--;
			}
			nused_threads=nthreads;
		}
		else
		{
			step=1;
			offset=step;
			remainder=0;
			nused_threads=n;
		}
		//assign for ithread=0 0->offset-1
		//logfile<<"split offset "<<offset<<endl;
		for (int j=0;j<offset;j++)
		{
			thread_arg[0].aenet_at_index[j]=ind[j];
			//logfile<<j<<" "<<thread_arg[0].aenet_at_index[j]<<endl;
		}
		thread_arg[0].aenet_at_count=offset;

		if (nthreads>1)
		{
			int first=0;
			for (ithread=1;ithread<nused_threads-1;ithread++)
			{
				first+=offset;
				offset=step;
				if (remainder>0)
				{
					//to spread the remainder
					offset=step+1;
					remainder--;
				} 
				for (int j=0;j<offset;j++)
					thread_arg[ithread].aenet_at_index[j]=ind[first+j];
				thread_arg[ithread].aenet_at_count=offset;
			}
			first+=offset;
			thread_arg[nused_threads-1].aenet_at_count=n-first;
			for (int j=0;j<thread_arg[nused_threads-1].aenet_at_count;j++)
				thread_arg[nused_threads-1].aenet_at_index[j]=ind[first+j];
			
			for (ithread=nused_threads;ithread<nthreads;ithread++)
			{
				thread_arg[ithread].aenet_at_index[0]=-1;//not to do anything
				thread_arg[ithread].aenet_at_count=0;

			}
		}
	};
	//final ANN calculation, if necessary 
	void Threads::ThreadAenetFinal(ThreadArg &my_arg)
	{
		/*
		ofstream *plog;
		if (RunParams::nthreads>1)
		{
			plog=&logfile;//for main
			switch (my_arg.thread_index)
			{
			case 0:
				plog=&logfile0;
				break;
			case 1:
				plog=&logfile1;
				break;
			case 2:
				plog=&logfile2;
				break;
			}
		}*/
			
		aenetnew_p->CalcANN(my_arg);
		
		//wait, till all the threads finished calculation
		{
			std::unique_lock<std::mutex> guard(count_mutex);//lock mutex
			complete_count2++;
			if (complete_count2==nthreads)//time to proceed
			{
				aenetstart_flag=0;//reset it to 0
				complete_count2=0;//all the threads finished calculation, reset it to 0
				check_count.notify_all();
			}
			else check_count.wait(guard);		
		}
	};
#endif	
