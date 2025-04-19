//source HistoSet.cpp
//Last changed 02.03.2023

//The concept is that only as many different bin sizes will be handled, as it is absolutely necessary. 
//As in case of g(r) data sets, the r point spacing and binshift is automatically detected, it can happen, that even in case of the 
//same binsize the binshift is different, so different histograms have to be used! In case of the same binsize but different binshifts
//it is examined, whether it is actually the same histogram with different starting bin. (For example binshift 0, 1, 2... are the same
//and will be considered as binshift 0, or 0.5, 1.5 and 2.5... will be considered as binshift 0.5. So only binshift<1 will occure for the different histograms! 
//So where binshifts with integer difference will share the same histogram, and only the first valid bin for the given data set is determined.
//To make thing easier, the histogram will be calculated from the first possible binshift (so if for example fom the 0, 1 and 2 series actually 1
//is the smallest used, it will be calculated for 0 as well.

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//classes1.h is included through this
longint *HistoSet::ngen;
int	 HistoSet::ntypes;//number of atom types
int  HistoSet::nvirtualtypes;//number of virtual types
int  HistoSet::npartials;//number of partials
int  HistoSet::nthreads;//the total number of threads to use
int  HistoSet::ndiffbin;//number of different bin sizes
int  HistoSet::ntot_datasets;//total number of data series
int  HistoSet::ntot_bins;//total number of bons for all the different bin sizes
int *HistoSet::assign_hist;//pointer to ExptsData::assign_hist, to assigning the different histogram bin sizes to the data sets, original dim [ntot_datasets]
int *HistoSet::nbins;//pointer to RunParams::nbin, as there can be different bin sizes, dim [ndiffbin]
int *HistoSet::offset;//array elements for each partial of each different bin size
int *HistoSet::calc_tabpot;//in case of potential=10 indicator for each partial, which partial is used dim [naprtials]

double HistoSet::boxedge;//that's the value of the boxedge unit in Angstrom
double *HistoSet::xmin;//minimum value for the abcissas (in boxedge units) dim [ndiffbin]
double *HistoSet::xmax;//pointer to RunParams::xmax, maximum value for the abcissas (in boxedge units) for each different bin size dim [ndiffbin]
double *HistoSet::deltax;//array of bins width (in reduced units) for each different bin size
double *HistoSet::inv_binwidth;//array of inverse of the bin width in reduced units for each different bin size

bool HistoSet::write_hist;//write the histogram, even if it was not calculated, the data set"s bin size assigment changed within acceptable limits

HistoSet::TCalcvdW HistoSet::CalcvdW;//pointer to the actual vdW function
HistoSet::TCalcvdW HistoSet::CalcvdW14;//pointer to the actual 1-4 vdW function
#ifdef _USE_LOCAL_INV
	int HistoSet::sum_bins;//max_loc_nbins*ntypes
	int HistoSet::sum_bins_offset;//sum_bins+padding
	int	HistoSet::max_loc_nbins;//number of bins to use for local invariance calculation 
	int *HistoSet::loc_offset;//array elements for each partial: total number of local bins to the beginning of the given partial
	double HistoSet::min_loc_r;//minimum distance for the local invariance calculation in reduced unit
	double HistoSet::max_loc_r;//maximum distance for the local invariance calculation in reduced unit
	double HistoSet::inv_locbinwidth;//inverse of the local invariance bin width in reduced units
#endif
#ifdef _NO_PERIODIC
	int *HistoSet::atom_binind;//the index of the bin the atom can be found in
	double HistoSet::R0_red;//reduced samole radius
#endif
#ifdef _VIBR_AMP
	int  HistoSet::nbins_conv;//the number of bins for the Gaussian convolution funtion
	int	 HistoSet::nbins_used;//the number of bins where the ppcf will be good, nbins-3*nbins_conv/5
	int *HistoSet::temp_conv_hist;//the convoluted histogram counts for one histogram bin
	double *HistoSet::thermal_conv;//convolution function for thermal vibration correction 
	double **HistoSet::thermal_conv_finder;//finder for each partial in the thermal conv array
	double *HistoSet::T_corr_sigma;//sigma parameter for the Gaussian convolution function for each partial
#endif
Move *HistoSet::move;//this cannot be initialized with the constructor, has to be assigned later

//-----------------------------default constructor----------------------
#ifdef _AENET
HistoSet::HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst, AvCoordConst &acoordconst,\
					Aenet &aen, FNC_POT* &fnc_const, Threads &thread_obj)
					:conf(config),rundat(rundata),coordnc(coordconst),acoordc(acoordconst), \
					aenet(aen), fnc(fnc_const), thread_object(thread_obj)
#else
HistoSet::HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst, \
						AvCoordConst &acoordconst,FNC_POT* &fnc_const, Threads &thread_obj)
							:conf(config),rundat(rundata),coordnc(coordconst),acoordc(acoordconst), \
					fnc(fnc_const), thread_object(thread_obj)
#endif
{ 
	int i,ib,padding_offset;
	int imax;
	if (nbins[0]<=0)
	{
		cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"HistoSet constructor"<<endl;
		cout<<"Number of bins was not calculated!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};	
	//Each thread will have its own histogram to avoid the use of mutexes during the histogram and its change calculation
	//the same applies for the totals.
	//Padding will be used for ptotal and pcounts to make sure, that the segments of different threads will not be cached together.
	//Presently longint is the type of both variable, if it is changed, than the padding_size calculation has to be changed accordingly
	imax=npartials*ntot_bins*nthreads;//size of storage array
	padding_offset=(imax>0 ? (int)((CACHE_PADDING-sizeof(*ptotal))/sizeof(*ptotal)) : 0);//number of dummy longinteger elements
	
	SetArraysize(&pcounts,imax+(nthreads-1)*padding_offset,"pcounts","HistoSet::HistoSet");
	SetArraysize(&ptotal,npartials*ndiffbin*nthreads+(nthreads-1)*padding_offset,"ptotal","HistoSet::HistoSet");
	SetArraysize(&finder,npartials*ndiffbin,"finder","HistoSet::HistoSet");//for the counts of the main thread only
	SetArraysize(&thread_histcount_finder,nthreads,"thread_histcount_finder","HistoSet::HistoSet");//to the beginning of segment in pcounts array for each thread
	//initializing the finders
	*(thread_histcount_finder+nthreads-1)=pcounts;//the main thread's histogram will be located at the beginning of the array,
					//as this will collect the contributions of the other threads and used in the calculation of the PPCF
	for (i=0;i<nthreads-1;i++)//for each auxiliary thread
		*(thread_histcount_finder+i)=pcounts+(i+1)*(npartials*ntot_bins+padding_offset);

	imax = 0;
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (i = 0; i < npartials; i++)//for each partial of each different bin size
		{
			*(finder + ib * npartials + i) = pcounts + imax;
			imax += nbins[ib];
		}
	}
	
#ifdef _VIBR_AMP
	SetArraysize(&pcounts_T,npartials*(nbins[0]+2),"pcounts_T","HistoSet::HistoSet");//it will be longer by two per partial, to put in those counts wich would go over the array bounds at convolution
	SetArraysize(&ptotal_T,npartials,"ptotal_T","HistoSet::HistoSet");
	SetArraysize(&finder_T,npartials,"finder_T","HistoSet::HistoSet");
	//initializing the finders
		for (i=0;i<npartials;i++)//for each partial
			*(finder_T+i)=pcounts_T+i*(nbins[0]+2);
#endif

#ifdef _NO_PERIODIC
	int size=1;//the finder will be cerated, but will not be used
	imax=1;//the ppcf will be calculated from the local histogram directly, the periodic_hist will not be calculated, so
				//time will be saved by not to copy it back and forth 
#ifdef _VIBR_AMP
	imax=npartials*nbins[0];
#endif
	if (RunParams::nek>0)
		imax=npartials*nbins[0];
	SetArraysize(&periodic_pcounts,imax,"periodic_pcounts","HistoSet::HistoSet");
	for (i=0;i<imax;i++)
		periodic_pcounts[i]=0.0;//initializing to zero

	if (imax>1)
		size=npartials;//the finder is needed
	SetArraysize(&periodic_finder,size,"periodic_finder","HistoSet::HistoSet");
	//initialising the finder
	periodic_finder[0]=periodic_pcounts;//set it anyhow
	if (imax>1)
	{
		for(i=1;i<npartials;i++)//for each partial
			*(periodic_finder+i)=periodic_pcounts+i*nbins[0];
	}

#endif
	
	//Initialization of  the histograms to 0
	InitHist();
	
	//In test mode in case of tooclose atoms the indices of the too close pairs are written in a file 
	// As all the threads can encounter tooclose atoms during the original histogram calculation, each thread
	//will store the indices of its too close pairs in the tcpairs array, and the ntcp will contain for the
	//number of too close pairs 
#ifdef _TEST_MODE
	maxtcp=config.ntotal* (config.ntotal-1)/2/10;//setting the array dim to the tenth of the possible pairs
	SetArraysize(&tcpairs,maxtcp,"tcpairs","HistoSet::HistoSet");
	ntcp=0;
#endif

#ifdef _USE_LOCAL_INV
	int *p1;
	int tot_moved;
	//array dimensions depending on the number of moved atoms, the Move::tot_moved_atoms will be set later
	if (RunParams::swap_fraction>0)
		tot_moved=RunParams::nmoved+1;//to have place for the extra atom the (nmoved-1)-th is swapped with
	else
		tot_moved=RunParams::nmoved;//only nmoved atom
	//local histograms for initial calculation
	imax=(nthreads-1)*tot_moved*sum_bins;//as the main thread works on the local_count array
	padding_offset=(imax>0 ? (int)((CACHE_PADDING-sizeof(*local_count))/sizeof(*local_count)) : 0);//number of dummy integer elements
	sum_bins_offset=sum_bins+(SimpleCfg::ntotal*sum_bins>0 ? (int)((CACHE_PADDING-sizeof(*local_count))/sizeof(*local_count)) : 0);
	//as different threads can write different atom's local hist., apply padding
	SetArraysize(&local_count,SimpleCfg::ntotal*sum_bins_offset,"local_count","HistoSet::HistoSet");
	SetArraysize(&dlocal_count,imax+(nthreads-2)*padding_offset,"dlocal_count","HistoSet::HistoSet");
	SetArraysize(&dlocal_count_finder,(nthreads-1),"dlocal_count_finder","HistoSet::HistoSet");
		
	p1=dlocal_count;
	for (i=0;i<imax+(nthreads-2)*padding_offset;i++)
		*p1++=0;
	//initializing the finders
	for (i=0;i<(nthreads-1);i++)
		dlocal_count_finder[i]=dlocal_count+i*(tot_moved*sum_bins+padding_offset);
	
	padding_offset=(max_loc_nbins*npartials*nthreads>0 ? (int)((CACHE_PADDING-sizeof(*local_sum_count))/sizeof(*local_sum_count)) : 0);//number of dummy longinteger elements
	SetArraysize(&local_sum_count,max_loc_nbins*npartials*nthreads+(nthreads-1)*padding_offset,"local_count","HistoSet::HistoSet");
	SetArraysize(&local_sum_count_finder,npartials,"local_sum_count_finder","HistoSet::HistoSet");
	SetArraysize(&thread_lochistcount_finder,nthreads,"thread_lochistcount_finder","HistoSet::HistoSet");//to the beginning of segment in local_sum_count array for each thread
	
	for (i=0;i<npartials;i++)//for each partial
			*(local_sum_count_finder+i)=local_sum_count+i*max_loc_nbins;
	*(thread_lochistcount_finder+nthreads-1)=local_sum_count;//the main thread's histogram will be located at the beginning of the array,
	for (i=0;i<nthreads-1;i++)//for each auxiliary thread
		*(thread_lochistcount_finder+i)=local_sum_count+(i+1)*(npartials*max_loc_nbins+padding_offset);

	InitLocHist();
	
#endif
	
};	
//--------------------copy constructor, not exact copy, only the parts used by the main thread are created and copied, as the others -------------------
//---------------are not used by the histold object--------------------------------------
#ifdef _AENET
HistoSet::HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst, AvCoordConst &acoordconst,
					Aenet &aen, FNC_POT* &fnc_const, Threads &thread_obj, HistoSet &source)
					:conf(config),rundat(rundata),coordnc(coordconst),acoordc(acoordconst), \
					aenet(aen), fnc(fnc_const), thread_object(thread_obj)
#else
HistoSet::HistoSet(SimpleCfg &config, RunParams &rundata, CoordNumbConst &coordconst, AvCoordConst &acoordconst,
					FNC_POT* &fnc_const, Threads &thread_obj, HistoSet &source)
					:conf(config),rundat(rundata),coordnc(coordconst),acoordc(acoordconst), \
					fnc(fnc_const), thread_object(thread_obj)
#endif
{ 
	longint *pl1,*pl2;
	int i, ib;
	int imax;
			
	if (nbins[0]<=0)
	{
		cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"HistoSet constructor"<<endl;
		cout<<"Number of bins was not calculated!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};
	
	//Presently longint is the type of both variable, if it is changed, than the padding_size calculation has to be changed accordingly
	imax=npartials*ntot_bins;//size of storage array
	
	SetArraysize(&pcounts,imax,"pcounts","HistoSet::HistoSet");
	SetArraysize(&ptotal,npartials*ndiffbin,"ptotal","HistoSet::HistoSet");
	SetArraysize(&finder,npartials*ndiffbin,"finder","HistoSet::HistoSet");//for the counts of the main thread only
	thread_histcount_finder=NULL;//not needed here
	
	//Copying the source histograms to target
	pl1=pcounts;
	pl2=source.pcounts;
	for (i=0;i<imax;i++)
		*pl1++=*pl2++;
	pl1=ptotal;
	pl2=source.ptotal;
	for (i=0;i<npartials*ndiffbin;i++)
		*pl1++=*pl2++;
		
	//initializing the finders
	imax = 0;
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (i = 0; i < npartials; i++)//for each partial of each different bin size
		{
			*(finder + ib * npartials + i) = pcounts + imax;
			imax += nbins[ib];
		}
	}
	
#ifdef _VIBR_AMP
	SetArraysize(&pcounts_T,npartials*(nbins[0]+2),"pcounts_T","HistoSet::HistoSet");//it will be longer by two/partial, to put in those counts which would go over the array bounds at convolution
	SetArraysize(&ptotal_T,npartials,"ptotal_T","HistoSet::HistoSet");
	SetArraysize(&finder_T,npartials,"finder_T","HistoSet::HistoSet");//for the counts of the main thread only

	//Copying the source histograms to target
	pl1=pcounts_T;
	pl2=source.pcounts_T;
	for (i=0;i<npartials*(nbins[0]+2);i++)
		*pl1++=*pl2++;
	pl1=ptotal_T;
	pl2=source.ptotal_T;
	for (i=0;i<npartials;i++)
		*pl1++=*pl2++;
	//initializing the finders
	for (i=0;i<npartials;i++)//for each partial
			*(finder_T+i)=pcounts_T+i*(nbins[0]+2);
#endif

#ifdef _NO_PERIODIC
	int size=1;//the finder will be cerated, but will not be used
	imax=1;//the ppcf will be calculated from the local histogram directly, the periodic_hist will not be calculated, so
				//time will be saved by not to copy it back and forth 
#ifdef _VIBR_AMP
	imax=npartials*nbins[0];
#endif
	if (RunParams::nek>0)
		imax=npartials*nbins[0];
	SetArraysize(&periodic_pcounts,imax,"periodic_pcounts","HistoSet::HistoSet");
	
	for (i=0;i<imax;i++)
		periodic_pcounts[i]=source.periodic_pcounts[i];

	if (imax>1)
		size=npartials;//the finder is needed
	SetArraysize(&periodic_finder,size,"periodic_finder","HistoSet::HistoSet");
	//initialising the finder
	periodic_finder[0]=periodic_pcounts;//set it anyhow
	if (imax>1)
	{
		for(i=1;i<npartials;i++)//for each partial
			*(periodic_finder+i)=periodic_pcounts+i*nbins[0];
	}

#endif

#ifdef _USE_LOCAL_INV
	int *p1,*p2;
	thread_lochistcount_finder=NULL;//not needed here
		
	//local histograms for initial calculation
		SetArraysize(&local_count,SimpleCfg::ntotal*sum_bins_offset,"local_count","HistoSet::HistoSet");
	SetArraysize(&local_sum_count,max_loc_nbins*npartials,"local_count","HistoSet::HistoSet");
	SetArraysize(&local_sum_count_finder,npartials,"local_sum_count_finder","HistoSet::HistoSet");
	
	dlocal_count=NULL;
	dlocal_count_finder=NULL;
	p1=local_count;
	p2=source.local_count;
	for (i=0;i<SimpleCfg::ntotal*sum_bins_offset;i++)
		*p1++=*p2++;

	pl1=local_sum_count;
	pl2=source.local_sum_count;
	for (i=0;i<max_loc_nbins*npartials;i++)
		*pl1++=*pl2++;
	for (i=0;i<npartials;i++)//for each partial
			*(local_sum_count_finder+i)=local_sum_count+i*max_loc_nbins;

#endif
	
};	

//-------------Sets the elements of the histograms and total count array to zero----
void HistoSet::InitHist()
{
	int i,imax,padding_offset;
	imax=npartials*ntot_bins*nthreads;//size of storage array
	padding_offset=(imax>0 ? (int)((CACHE_PADDING-sizeof(*ptotal))/sizeof(*ptotal)) : 0);//number of dummy longinteger elements

	//initialisation of histogram counts to 0
	for (i=0;i<imax+(nthreads-1)*padding_offset;i++)
		pcounts[i]=0;//histogram initialisation
	for (i=0;i< npartials*ndiffbin*nthreads + (nthreads - 1)*padding_offset; i++)
		ptotal[i]=0;//total counts init.
#ifdef _VIBR_AMP
	for (i=0;i<npartials*(nbins[0]+2);i++)
		pcounts_T[i]=0;//histogram initialisation
	for (i=0;i<npartials; i++)	
		ptotal_T[i]=0;//total counts init.
#endif
}	

#ifdef _USE_LOCAL_INV
	//-------------Sets the elements of the local histograms count array to zero----
	void HistoSet::InitLocHist()
	{
		int i,padding_offset;
		for (i=0;i<SimpleCfg::ntotal*sum_bins_offset;i++)
			local_count[i]=0;
		padding_offset=(max_loc_nbins*npartials*nthreads>0 ? (int)((CACHE_PADDING-sizeof(*local_sum_count))/sizeof(*local_sum_count)) : 0);//number of dummy longinteger elements
		for (i=0;i<max_loc_nbins*npartials*nthreads+(nthreads-1)*padding_offset;i++)
			local_sum_count[i]=0;
	
	}
#endif

//-------------Sets the static members------------------------
void HistoSet::SetHistParams(RunParams &rundata)
{
	int i,ib;
	int sum;

	//Checking the parameters
	if (RunParams::ntypes==0)
	{
		cout<<"\n*****ERROR*****"<<endl;
		cout<<"HistoSet::SetHistParams: Number of types is 0!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};
	if (rundata.nbins[0]<=0)
	{
		cout<<"\n*****ERROR*****"<<endl;
		cout<<"HistoSet::SetHistParams: Number of bins was not calculated!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	ntypes=RunParams::ntypes;//number of atom types
	nvirtualtypes = SimpleCfg::nvirtualtypes;//number of virtual types
	ntot_datasets = ExptsData::ntot_datasets;//total number of data sets
	ndiffbin = RunParams::ndiffbin;//number of different bin sizes
	assign_hist = ExptsData::assign_hist;//assigning the data sets to the histogram bin types, needed for output
	npartials=ntypes*(ntypes+1)/2;//number of partials
	nbins=rundata.nbins;//number of bins (same for each histogram)
	nthreads=RunParams::nthreads;//number of threads
	boxedge=RunParams::boxedge;//that's the value of the boxedge unit in Angstrom
	xmin=rundata.xmin;//minimum value for the abcissas (in reduced units)
	xmax=RunParams::xmax;//maximum value for the abcissas (in reduced units)
	SetArraysize(&deltax, ndiffbin, "deltax", "HistoSet::SetHistParams");
	SetArraysize(&inv_binwidth, ndiffbin, "inv_binwidth", "HistoSet::SetHistParams");
	ntot_bins = 0;
	for (i = 0; i < ndiffbin; i++)
	{
		deltax[i] = (xmax[i] - xmin[i]) / nbins[i];//bins width (in reduced units)
		inv_binwidth[i] = 1.0 / deltax[i];//inverse of the bin width in reduced units
		ntot_bins += nbins[i];
	}
	//Memory allocation and initialisation for offset element
	SetArraysize(&offset,npartials*ndiffbin,"offset","HistoSet::SetHistParams");
#ifdef _USE_LOCAL_INV
	SetArraysize(&loc_offset,npartials,"loc_offset","HistoSet::SetHistParams");
	max_loc_nbins=RunParams::max_loc_nbins;//number of bins for local inv calculation
	min_loc_r=RunParams::min_loc_r;//minimum distance for local inv calculation
	max_loc_r=RunParams::max_loc_r;//maximum distance for local inv calculation
	sum_bins=max_loc_nbins*ntypes;
	
	inv_locbinwidth=1.0/(max_loc_r-min_loc_r)*max_loc_nbins;//inverse of the local bin width in reduced units
#endif
	sum = 0;
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (i = 0; i < npartials; i++)
		{
			offset[ib*npartials + i] = sum;
			sum += nbins[ib];
		}
	}
#ifdef _USE_LOCAL_INV
	for (i = 0; i < npartials; i++)
		loc_offset[i]=max_loc_nbins*i;
#endif
	

	switch (RunParams::potential)
	{
		case 1:
		{
			CalcvdW = &FNC_POT::CalcLJ;
			CalcvdW14 = &FNC_POT::CalcLJ14;
			break;
		}
		case 10:
		{
			CalcvdW = &FNC_POT::CalcTabPot;
			SetArraysize(&calc_tabpot,npartials,"calc_tabpot", "HistoSet::SetHistParams");
			for (i = 0; i < SimpleCfg::npartials; i++)
				calc_tabpot[i] = 0;
		
			for (i = 0; i < FNC_POT::nused_potpartials; i++)
				calc_tabpot[FNC_POT::tabpot_index[i]] = 1;
			
			break;
		}

	}

#ifdef _NO_PERIODIC
	SetArraysize(&atom_binind,SimpleCfg::ntotal,"atom_binind","HistoSet::SetHistParams");
	R0_red=RunParams::R0/boxedge;
#endif

#ifdef _VIBR_AMP
	int itype,jtype,ipartial=0;
	SetArraysize(&T_corr_sigma,npartials,"T_corr_sigma","HistoSet::SetHistParams");
	for (itype=0; itype<ntypes; itype++)//for each first atom type
	{
		for(jtype=itype;jtype<ntypes;jtype++)//for each second atom type
		{
			if (itype==jtype)
				T_corr_sigma[ipartial]=RunParams::T_corr_sigma_t[itype];
			else
				T_corr_sigma[ipartial]=(RunParams::T_corr_sigma_t[itype]+RunParams::T_corr_sigma_t[jtype])/2;
			ipartial++;
		}
	}
#endif

}

//----------------save to the *.hgm file--------------------------
void HistoSet::Save( const char *file_name) const
{
	int ib;
	longint *browser1, *browser2;
	int i,count,itype, jtype;
	char *name,*conv_numb=NULL;

	
	SetArraysize(&name,40,"name","HistoSet::Save");
#if defined(_VIBR_AMP) || defined(_NO_PERIODIC)
	int ipartial;
	ipartial=0;
#endif
#ifdef _VIBR_AMP
	mystrcpy(name,40,"Original histogram ");
#else
	mystrcpy(name,40,"");
#endif
	ofstream file;
	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, hgmfilename);
	OpenFile(file,tempfilename,"HistoSet::Save",0);//open file, check, whether it was successfully opened

	
	file.precision(18);
	file<<"This is an object of Class HistoSet containing the unique-pair histogram(s)"<<endl;
	file<<ntypes<<"\t number of atom types"<<endl;
	file << ndiffbin << "\t number of different bin sizes" << endl;
	file<<boxedge<<"\t boxedge (i.e. cell edge) in Angstroms"<<endl;
		
	browser1=pcounts;
	browser2=ptotal;
	for (ib=0;ib<ndiffbin;ib++)//for each different bin type
	{
		
		file << "\nFor data set(s): ";
		count = 0;
		for (i = 0; i < ntot_datasets; i++)
		{
			if (assign_hist[i] == ib)//this binsize is used for this data set
			{
				if (count > 0)
					file << ", ";
				file << i + 1;
				count++;
			}
		}
		file << endl;
		file << nbins[ib] << "\t number of histogram bins " << endl;
		file.setf(ios::scientific, ios::floatfield);
		file.precision(18);
		file << xmin[ib] << "\t xmin value (reduced) " << "\t rmin "; 
		file.precision(5); 
		file<< xmin[ib] * boxedge << " A" << endl;
		file.precision(18);
		file << xmax[ib] << "\t xmax value (reduced) ";
		file.precision(5); 
		file<< "\t rmax " << xmax[ib] * boxedge << " A" << endl;
		file.precision(18);
		file << deltax[ib] << "\t bin width (deltax) ";
		file.precision(5);
		file<< "\t bin width " << RunParams::rspacing_diff[ib] << " A" << endl;
		file << "\n bin number, counts in each bin:" << endl;

		for (itype=1; itype<=ntypes; itype++)//for each first atom type
		{
			for (jtype = itype; jtype <= ntypes; jtype++)//for each second atom type
			{
				
				file << endl << name << "partial " << itype << "-" << jtype << endl;
				file << *browser2++ << "\t total counts " << endl;
				for (i = 1; i <= nbins[ib]; i++)//for each histogram bin
					file << i << "\t" << *browser1++ << endl;

			}//end jytpe
		}//end itype
	}//end ib


#ifdef _VIBR_AMP

	//here always the default, nbins[0] histogram bin is used
	mystrcpy(name,40,"Histogram with thermal correction ");
	ipartial=0;
	browser2=ptotal_T;
	for (itype=1; itype<=ntypes; itype++)//for each first atom type
	{
		for(jtype=itype;jtype<=ntypes;jtype++)//for each second atom type
		{
			browser1=finder_T[ipartial];//as the convoluted histogram has two extra bins for each partial
			file<<endl<<name<<"partial "<<itype<<"-"<<jtype<<endl;
			file<<T_corr_sigma[ipartial]<<"\t the sigma for the Gauss convolution function for thermal vibration correction in Angstrom for this type"<<endl;
			file<<*browser2++<<"\t total counts "<<endl;
			for(i=1;i<=nbins[0];i++)//for each histogram bin
				file<<i<<"\t"<<*browser1++<<endl;
			ipartial++;
		}
	}
#endif
#ifdef _NO_PERIODIC
	//here always the default, nbins[0] histogram bin is used
	double *p;
	if (DataMat::calcmode==0)
	{
		mystrcpy(name,40,"Periodic histogram ");
		ipartial=0;
	
		for (itype=1; itype<=ntypes; itype++)//for each first atom type
		{
			for(jtype=itype;jtype<=ntypes;jtype++)//for each second atom type
			{
				p=periodic_finder[ipartial];//as the convoluted histogram has two extra bins for each partial
				file<<endl<<name<<"partial "<<itype<<"-"<<jtype<<endl;
				
				for(i=1;i<=nbins[0];i++)//for each histogram bin
					file<<i<<"\t"<<*p++<<endl;
				ipartial++;
			}
		}
	}
#endif
	delete [] name;
	if (conv_numb != NULL)
		delete [] conv_numb;
	file.close();

};

#ifdef _USE_LOCAL_INV
	//----------------save to the *.hgm file--------------------------
	void HistoSet::SaveLoc( const char *file_name) const
	{
		int *browser;
		int i,itype, jtype,ic,n_centre,jc;
		int cum_off;
		int iat,n[2];
		longint *psum;
//		longint sum=0; 
		ofstream file;
		if (strlen(file_name)!=0)
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
		else
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, lhgmfilename);
		OpenFile(file,tempfilename,"HistoSet::SaveLoc",0);//open file, check, whether it was successfully opened
		
		file.precision(18);
		file<<"This is an object of Class HistoSet for the local histograms"<<endl;
		file<<ntypes<<"\t number of atom types"<<endl;
		file<<max_loc_nbins<<"\t number of local histogram bins"<<endl;
		file<<boxedge<<"\t boxedge (i.e. cell edge) in Angstroms"<<endl;
		file<<min_loc_r<<"\t minimum value (reduced)"<<endl;
		file<<max_loc_r<<"\t maximum value (reduced)"<<endl;
		file<<1.0/inv_locbinwidth<<"\t bin width (reduced)"<<endl;
		file<<"\n bin number, counts in each bin for each atoms:"<<endl;

		
		psum=local_sum_count;	
		for (itype=0; itype<ntypes; itype++)//for each first atom type
		{
			for(jtype=itype;jtype<ntypes;jtype++)//for each second atom type
			{
				file<<endl<<"partial "<<itype+1<<"-"<<jtype+1;
				n_centre=1;//default, clean partial
				if (itype!=jtype)
				{
					n_centre=2;//mixed partial
					file<<" and partial "<<jtype+1<<"-"<<itype+1;
				}
				file<<" for all atoms:";
				//type of the actual centre atoms
				n[0]=itype;
				n[1]=jtype;
				for (ic=0;ic<n_centre;ic++)		
				{
					for (iat=SimpleCfg::cumul[n[ic]];iat<SimpleCfg::cumul[n[ic]+1];iat++)
						file<<"\t "<<iat+1;
				}
				file<<"\t sum"<<endl;

				file<<endl;//skip the row containing the total count for normal histogram
								
				for(i=0;i<max_loc_nbins;i++)//for each local histogram bin
				{
				
					file<<i+1;//index 

					//calculating the total value for a bin, needed for checking, this can be useful to check against the double of the normal histogram, but files written with this cannot be read back!
/*//
					sum=0;
					for (ic=0;ic<n_centre;ic++)		
					{	
						cum_off=0;
						jc=!(ic&1);
						browser=local_count+SimpleCfg::cumul[n[ic]]*sum_bins_offset+loc_offset[n[jc]]+i;//setting to the beginning of the local hist of the first atom of this type for this bin
						for (iat=SimpleCfg::cumul[n[ic]];iat<SimpleCfg::cumul[n[ic]+1];iat++)
						{
							sum+=*(browser+cum_off);
							cum_off+=sum_bins_offset;
						}
					}

					file<<"\t"<<sum;

*///					

					for (ic=0;ic<n_centre;ic++)		
					{
						cum_off=0;
						jc=!(ic&1);
						browser=local_count+SimpleCfg::cumul[n[ic]]*sum_bins_offset+loc_offset[n[jc]]+i;//setting to the beginning of the local hist of the first atom of this type for this bin
						for (iat=SimpleCfg::cumul[n[ic]];iat<SimpleCfg::cumul[n[ic]+1];iat++)
						{
							file<<"\t"<<*(browser+cum_off);
							cum_off+=sum_bins_offset;//go to the next atom's local hist count for this partial and bin
						}
					}

					file<<"\t"<<*psum++<<endl;
				}
				
			}
		}
#ifdef _LOCAL_INV
		//saving the average histogram in case of distance-based calculation
		if (RunParams::loc_chi2_mode == 1)
		{
			int j, k;
			file << "\nAverage local histogram bins of the neighbours in the whole range where the localc invariance is calculated" << endl;
			file << "in case of distance-based calculation for all the partials" << endl;
			file << "Only for the neighbours to calulate for, so for practical reasons it will be written " << endl;
			file << "NOT in RMC-order, but in order 1-1, 2-1, 3-1 ... " << endl;
			file << "                               1-2, 2-2, 3-2 ... " << endl;
			file << "                               1-3, 2-3, 3-3 ... " << endl;
			for (j = 0; j < ntypes; j++)//neighbour type
			{
				file << "\nneigh_ind";
				for (i = 0; i < ntypes; i++)//central type
				{
					file << "\t " << i + 1 << "-" << j + 1;
					if (i == ntypes - 1)
						file << endl;
				}
				for (k=0;k<ChiSquared::loc_natoms[j];k++)//going through the neighbours
				{
					file << k + 1;
					for (i = 0; i < ntypes; i++)//central type
						file << "\t " << *(ChiSquared::av_loc_hist + ChiSquared::av_loc_hist_offset[i * ntypes + j] + k);
					file << endl;
				}

			}
		}
#endif
		file.close();
		
	}

	void HistoSet::LoadLocCount(ifstream &file)//loads the histogram counts from a file
	{
		int i,iat,ic,jc,itype,jtype;
		int cum_off,n_centre,n[2];
		int *browser;
		longint *psum;
		//the file has to be already open, not cheked

		psum=local_sum_count;
		for (itype=0; itype<ntypes; itype++)//for each first atom type
		{
			for(jtype=itype;jtype<ntypes;jtype++)//for each second atom type
			{
				SkipLine(file,lhgmfilename,1);
				SkipLine(file,lhgmfilename,1);
				SkipLine(file,lhgmfilename,1);

				n_centre=1;//default, clean partial
				if (itype!=jtype)
					n_centre=2;//mixed partial
				
				//type of the actual centre atoms
				n[0]=itype;
				n[1]=jtype;
			
				for(i=0;i<max_loc_nbins;i++)//for each local histogram bin
				{
					ic=ReadThisLine(file,3,1,"bin index", "HistoSet::LoadLocCount",lhgmfilename);//that's the bin index, not needed
					for (ic=0;ic<n_centre;ic++)		
					{
						cum_off=0;
						jc=!(ic&1);
						browser=local_count+SimpleCfg::cumul[n[ic]]*sum_bins_offset+loc_offset[n[jc]]+i;//setting to the beginning of the local hist of the first atom of this type for this bin
						for (iat=SimpleCfg::cumul[n[ic]];iat<SimpleCfg::cumul[n[ic]+1]-1;iat++)
						{
							file>>*(browser+cum_off);//this would take a long time for a big configuration, no check is used =ReadThisLine(file,3,1,"local_count", "HistoSet::LoadLocCount",lhgmfilename);//bin count
							cum_off+=sum_bins_offset;//go to the next atom's local hist count for this partial and bin
						}
						//for the last atom
						*(browser+cum_off)=ReadThisLine(file,3,0,"local_count", "HistoSet::LoadLocCount",lhgmfilename);//bin count

					}
					*psum++=ReadThisLine(file,1, (longint)0,"local_sum_count", "HistoSet::LoadLocCount",lhgmfilename);//average bin count
				}
			}
		}
		return;
	};
#endif

	int HistoSet::Load(int mode)//loads from a file
	//0 normal histogram
	//1 local histogram
	{
		int i,j, dumb;
		int flag = 1;//return status
		std::streamoff current_pos, next_pos;
		double realdumb, diff;
		char *fname, *histname,*text2;
		char *text = NULL;
		
		ifstream file;

		SetArraysize(&histname, FILE_NAME_SIZE + 5, "histname", "HistoSetLoad");
		fname = hgmfilename;//default
		mystrcpy(histname, FILE_NAME_SIZE + 5, "Histogram");
#ifdef _USE_LOCAL_INV
		if (mode)
		{
			fname = lhgmfilename;
			mystrcpy(histname, FILE_NAME_SIZE + 5, "Local histogram");
		}
#endif
		SafeOpenTextFile(file, fname);
		if (CheckFileState(file, "HistoSet::Load", fname) == 0)
		{
			cout << "\nWARNING(" << ++warn << "): Histogram has to be calculated!" << endl;
			return(0);//load was not successful
		};
#ifdef _USE_LOCAL_INV
		if (mode)
			cout << "\nLoading the local histogram!" << endl;
		else
			cout << "\nLoading the histogram!" << endl;
#else
		cout << "\nLoading the histogram!" << endl;
#endif


		SkipLine(file, fname, 1);
		dumb = ReadThisLine(file, 1, 1, "ntypes", "HistoSet::Load", fname);//ntypes:number of atom types
		if (dumb != ntypes)
		{
			cout << "\nWARNING(" << ++warn << "): The number of atom types is " << ntypes << " in the .cfg and " << dumb << \
				" in the " << fname << " file are not the same!" << endl;
			cout <<"\t"<< histname << " has to be calculated!" << endl;
			return(0);//load was not successful
		}
#ifdef _USE_LOCAL_INV
		if (mode)//there can be different binsizes for the normal histogram, they are removed from the header
		{
			dumb = ReadThisLine(file, 1, 1, "max_loc_nbins", "HistoSet::Load", fname);//max_loc_nbins:number of local histogram bins
			if (dumb != max_loc_nbins)
			{
				cout << "\nWARNING(" << ++warn << "): The number of max_loc_bins is " << max_loc_nbins << " calculated according"<<endl;
				cout<<"\tto the parameters in .dat file and " << dumb << " in the " << fname << " file are not the same!" << endl;
				cout << histname << " has to be calculated!" << endl;
				return(0);//load was not successful
			}
		}
		else
		{
#endif
			dumb = ReadThisLine(file, 1, 1, "ndiffbin", "HistoSet::Load", fname);//number of different bin sizes
			if (dumb != ndiffbin)
			{
				cout << "\nWARNING(" << ++warn << "): The number of different histogram bin sizes is " << ndiffbin << " calculated according to the parameters in" ;
				cout<<"\t the "<<datfilename<<" file and " << dumb << " in the " << fname << " file are not the same!" << endl;
				cout << "\t"<<histname << " has to be calculated!" << endl;
				return(0);//load was not successful
			}
#ifdef _USE_LOCAL_INV
		}

#endif

		realdumb = ReadThisLine(file, 1, 1.0, "boxedge", "HistoSet::Load", fname);//boxedge
		diff = fabs(realdumb - boxedge);
		if (diff > 0)
		{
			cout.precision(16);
			cout << "\nWARNING(" << ++warn << "): The boxedge is " << boxedge << " in .cfg file and " << realdumb << \
				" in the " << fname << " file are not the same!" << endl;
			if (diff > 1e-15)
			{
				cout <<"\t"<< histname << " has to be calculated!" << endl;
				return(0);//load was not successful
			}
			else
				cout << "\tThe difference is " << diff << ", the value from the .cfg file will be used!" << endl;
		}
#ifdef _USE_LOCAL_INV
		if (mode)
		{
			realdumb = ReadThisLine(file, 1, 1.0, "min_loc_r", "HistoSet::Load", fname);//min_loc_r value
			diff = fabs(realdumb - min_loc_r);
			if (diff > 0)
			{
				cout.precision(16);
				cout << "\nWARNING(" << ++warn << "): The minimum distance to calculate the local histogram is " << min_loc_r << " calculated " << endl;
				cout<<"\tin the program and	"<< realdumb << " in the " << fname << " file are not the same!" << endl;
				if (diff > 1e-14)
				{
					cout <<"\t"<< histname << " has to be calculated!" << endl;
					return(0);//load was not successful
				}
				else
				{
					cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
				}
			}
			//It is dangerous to use the rounded max_loc_r value, it can cause overflow
			realdumb = ReadThisLine(file, 1, 1.0, "max_loc_r", "HistoSet::Load", fname);//max_loc_r value
			diff = fabs(realdumb - max_loc_r);
			if (diff > 0)
			{
				cout.precision(16);
				cout << "\nWARNING(" << ++warn << "): The maximum distance to calculate the local histogram is " << max_loc_r << " calculated" << endl;
				cout<<"\t in the program and "<< realdumb << " in the " << fname << " file are not the same!" << endl;
				if (diff > 1e-14)
				{
					cout <<"\t"<< histname << " has to be calculated!" << endl;
					return(0);//load was not successful
				}
				else
				{
					cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
				}
			}
			//It is dangerous to use this rounded deltax value, it can produce unforseen effects	
			realdumb = ReadThisLine(file, 1, 1.0, "locbinwidth", "HistoSet::Load", fname);//locbinwidth
			diff = fabs(realdumb - 1.0 / inv_locbinwidth);
			if (diff > 0)
			{
				cout.precision(16);
				cout << "\nWARNING(" << ++warn << "): The bin width of the histogram is " << 1.0 / inv_locbinwidth << " calculated " << endl;
				cout<<"\t"	<< realdumb << " in the " << fname << " file are not the same!" << endl;
				if (diff > 1e-14)
				{
					cout <<"\t"<< histname << " has to be calculated!" << endl;
					return(0);//load was not successful
				}
				else
				{
					cout.setf(ios::right, ios::adjustfield);
					cout.setf(ios::scientific, ios::floatfield);
					cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
					cout.setf(ios::fixed, ios::floatfield);

				}
			}
		}
		else
		{
#endif
			

#ifdef _USE_LOCAL_INV
		}
#endif

		

#ifdef _USE_LOCAL_INV
		if (mode)
		{
			SkipLine(file);//empty
			SkipLine(file);//comment (bin number, counts in each bin:)
			LoadLocCount(file);//loading the local histogram counts
		}
#endif
		if (mode == 0)
		{
			SetArraysize(&text2, NAME_SIZE, "text2", "HistoSet::Load");
			write_hist = false;
			for (i = 0; i < ndiffbin; i++)
			{
				SkipLine(file, fname, 1);//blank
				current_pos = file.tellg();
				SkipLine(file);
				next_pos = file.tellg();
				file.seekg(current_pos);
				for (j = 0; j < 3; j++)
					file >> text2;
				
				do
				{

					dumb = ReadThisLine(file, 4, -1, "data set index", "HistoSet::Load");//
					if (assign_hist[dumb - 1] != i)
						write_hist = true;
					current_pos = file.tellg();
					file >> text2;
				} while (current_pos<next_pos);
				file.seekg(next_pos);//reste it to the beginning of the next line 
				IntToStr(&text, i + 1);
				mystrcpy(text2, NAME_SIZE,"nbins[");
				mystrcat(text2,NAME_SIZE, text);
				mystrcat(text2, NAME_SIZE, "]");
				dumb = ReadThisLine(file, 1, 1, text2, "HistoSet::Load", fname);//nbins:number of histogram bins
				if (dumb != nbins[i])
				{
					cout << "\nWARNING(" << ++warn << "): The number of bins for the " << i + 1 << " bin size is " << nbins[i] <<" calculated according to" <<endl;
					cout<<"\tthe parameters in .dat file and " << dumb << " in the " << fname << " file are not the same!" << endl;
					cout <<"\t"<< histname << " has to be calculated!" << endl;
					return(0);//load was not successful
				}
				mystrcpy(text2, NAME_SIZE, "xmin[");
				mystrcat(text2, NAME_SIZE, text);
				mystrcat(text2, NAME_SIZE, "]");
				realdumb = ReadThisLine(file, 1, 1.0, text2, "HistoSet::Load", fname);//xmin value
				diff = fabs(realdumb - xmin[i]);
				if (diff > 0)
				{
					cout.precision(16);
					cout << "\nWARNING(" << ++warn << "): The minimum distance to calculate the histogram for the " << i + 1 << " bin size is " << xmin[i] << " calculated" << endl;
					cout<<"\tin the program and "<< realdumb << " in the " << fname << " file are not the same!" << endl;
					if (diff > 1e-14)
					{
						cout <<"\t"<< histname << " has to be calculated!" << endl;
						return(0);//load was not successful
					}
					else
					{
						cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
					}
				}
				//It is dangerous to use the rounded xmax value, it can cause overflow
				mystrcpy(text2, NAME_SIZE, "xmax[");
				mystrcat(text2, NAME_SIZE, text);
				mystrcat(text2, NAME_SIZE, "]");
				realdumb = ReadThisLine(file, 1, 1.0, text2, "HistoSet::Load", fname);//xmax value
				diff = fabs(realdumb - xmax[i]);
				if (diff > 0)
				{
					cout.precision(16);
					cout << "\nWARNING(" << ++warn << "): The maximum distance to calculate the histogram for the " << i + 1 << " bin size is " << xmax[i] << " calculated  " << endl;
					cout<<"\tin the program and "<< realdumb << " in the " << fname << " file are not the same!" << endl;
					if (diff > 1e-14)
					{
						cout << "\t"<<histname << " has to be calculated!" << endl;
						return(0);//load was not successful
					}
					else
					{
						cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
					}
				}
				//It is dangerous to use this rounded deltax value, it can produce unforseen effects	
				mystrcpy(text2, NAME_SIZE, "deltax[");
				mystrcat(text2, NAME_SIZE, text);
				mystrcat(text2, NAME_SIZE, "]");
				realdumb = ReadThisLine(file, 1, 1.0, text2, "HistoSet::Load", fname);//deltax
				diff = fabs(realdumb - deltax[i]);
				if (diff > 0)
				{
					cout.precision(16);
					cout << "\nWARNING(" << ++warn << "): The bin width of the histogram is for the " << i + 1 << " bin size is " << deltax[i] << " calculated" << endl;
					cout<<"\tin the program and "<< realdumb << " in the " << fname << " file are not the same!" << endl;
					if (diff > 1e-14)
					{
						cout <<"\t"<< histname << " has to be calculated!" << endl;
						return(0);//load was not successful
					}
					else
					{
						cout.setf(ios::right, ios::adjustfield);
						cout.setf(ios::scientific, ios::floatfield);
						cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
						cout.setf(ios::fixed, ios::floatfield);
						cout.setf(ios::right, ios::adjustfield);

					}
				}
				SkipLine(file);//empty
				SkipLine(file);//comment (bin number, counts in each bin:)
				LoadCount(file,0,i);//loading the histogram counts
			}
		}//is mode 0, normal histogram 
	
#ifndef _USE_LOCAL_INV
#ifdef _VIBR_AMP
		flag=LoadCount(file,1);
#endif
#endif

	//Checking, whether loading was successful
	if (CheckReadFileState(file,"HistoSet::Load",fname))
		return(flag);//load was successful
	else
	{
		cout << "\nWARNING(" << ++warn << "): Loading was not successful, histogram and coordination numbers will be calculated!"<<endl;
		return(flag);//loading was not successful
	}
	next_line_pos = -1;
	file.close();
	if (text != NULL)
		delete [] text;
	delete [] text2;
	return (flag);
};

int HistoSet::LoadCount(ifstream &file,int mode, int binsize_index)//loads the histogram counts from a file
{
	//mode 0: loads the normal histogram, binsize index is needed, as is loads the counts only for the given binsize at a call
	//mode 1: loads the convoluted histogram in case of _VIBR_AMP

	int i,ipartial;//,dumb
	longint *browser1=0, *browser2;
	//the file has to be already open, not cheked
#ifdef _VIBR_AMP
	double realdumb,diff;
	if (mode)
	{
		browser2=ptotal_T;
	}
	else
	{
#endif
		browser1=finder[binsize_index*npartials];
		browser2=ptotal+binsize_index*npartials;
#ifdef _VIBR_AMP
	}
#endif
	for (ipartial = 0; ipartial < npartials; ipartial++)//for each partial
	{
		SkipLine(file, hgmfilename, 1);//empty
		SkipLine(file, hgmfilename, 1);//xx partial
#ifdef _VIBR_AMP
		if (mode == 1)
		{
			browser1 = finder_T[ipartial];//as the convoluted histogram has two extra bins for each partial
			realdumb = ReadThisLine(file, 1, 0.01, "T_corr_sigma[ipartial]", "HistoSet::Load", hgmfilename);//deltax
			diff = fabs(realdumb - T_corr_sigma[ipartial]);
			if (diff > 0)
			{
				cout.precision(16);
				cout << "\nWARNING(" << ++warn << "): The sigma for the Gaussian thermal vibration convolution function " << T_corr_sigma[ipartial] << endl;
				cout<<"\tread in RunParams::GetParams and "	<< realdumb << " in the " << hgmfilename << " file are not the same!" << endl;
				if (diff > 1e-14)
				{
					cout << "\tHistogram has to be calculated!" << endl;
					return(0);//load was not successful
				}
				else
				{
					cout.unsetf(ios::right);
					cout.unsetf(ios::fixed);
					cout.setf(ios::right, ios::scientific);
					cout << "\tThe difference is " << diff << ", the value calculated in the program will be used!" << endl;
					cout.setf(ios::fixed, ios::floatfield);
					cout.setf(ios::right, ios::adjustfield);

				}
			}
		}
#endif
		*browser2 = ReadThisLine(file, 1, (longint)1, "total count", "HistoSet::LoadCount", hgmfilename);//total counts
		browser2++;
		for (i = 1; i <= nbins[binsize_index]; i++)//for each histogram bin
		{
			ReadThisLine(file, 3, 1, "bin index", "HistoSet::LoadCount", hgmfilename);//that's the bin index
			*browser1 = ReadThisLine(file, 1, (longint)1, "count", "HistoSet::LoadCount", hgmfilename);//bin count
			browser1++;
		}
	}
	return(1);
	
};

//-------Copy the modified parts of the histogram to target-------
void HistoSet::CopyModified(HistoSet &target,ThreadArg &thread_arg)
{
	int ipartial,ibin,ib;
	int min_ind,max_ind;//cycle boundaries
	longint *p_source,*p_target;
#ifdef _VIBR_AMP
	longint *p_source_T,*p_target_T;
#endif
#ifdef _NO_PERIODIC
	double *p_source_per=0, *p_target_per=0;
#endif
#ifdef _USE_LOCAL_INV
	int i,*min_loc_ind,*max_loc_ind,itype,jtype,loc_min_ind, loc_max_ind;//cycle boundaries
	int *ploc_source,*ploc_target;
		
	min_loc_ind=thread_arg.min_atom_index;
	max_loc_ind=thread_arg.max_atom_index;
	loc_min_ind=thread_arg.loc_min_bin;
	loc_max_ind=thread_arg.loc_max_bin;
#endif

	for (ib = 0; ib < ndiffbin; ib++)
	{
		min_ind = thread_arg.min_bin[ib];
		max_ind = thread_arg.max_bin[ib];
		for (ipartial = 0; ipartial < npartials; ipartial++)
		{
			if (Threads::mod_partial[ipartial])
			{

				//counts
				p_source = finder[ib*npartials + ipartial] + min_ind;
				p_target = target.finder[ib*npartials + ipartial] + min_ind;
#ifdef _VIBR_AMP
				p_source_T = finder_T[ipartial] + min_ind;
				p_target_T = target.finder_T[ipartial] + min_ind;
#endif

				for (ibin = min_ind; ibin <= max_ind; ibin++)
				{
					*p_target++ = *p_source++;
#ifdef _VIBR_AMP
					*p_target_T++ = *p_source_T++;
#endif
				}
#ifdef _NO_PERIODIC
				if (DataMat::calcmode == 0)//only in this case is periodic used
				{
					p_source_per = periodic_finder[ipartial] + min_ind;
					p_target_per = target.periodic_finder[ipartial] + min_ind;
					for (ibin = min_ind; ibin <= max_ind; ibin++)
						*p_target_per++ = *p_source_per++;

				}
#endif

#ifdef _USE_LOCAL_INV
				//counts calculated using the local histogram bin size 
				p_source = local_sum_count_finder[ipartial] + loc_min_ind;
				p_target = target.local_sum_count_finder[ipartial] + loc_min_ind;
				for (ibin = loc_min_ind; ibin <= loc_max_ind; ibin++)
					*p_target++ = *p_source++;
#endif
				//total counts
				if (min_ind == 0)
				{
					*(target.ptotal + ib * npartials + ipartial) = *(ptotal + ib * npartials + ipartial);//only the thread_id=0 will copy total
#ifdef _VIBR_AMP
					*(target.ptotal_T + ipartial) = *(ptotal_T + ipartial);//only the thread_id=0 will copy total
#endif 
				}

			}//end of if this partial has been modified
		}//end of cycle ipartial
	}//end of cycle diff bin size
#ifdef _USE_LOCAL_INV
	//here only those partials exist for each atom, in which the given atom type is involved, ntypes partial for each atom
	//has to be handled separately from histogram
	for (itype=0;itype<ntypes;itype++)
	{
		for (jtype=0;jtype<ntypes;jtype++)
		{
			ipartial=(itype<=jtype ? (itype*ntypes-(itype*(itype+1)/2)+ jtype) : \
				(jtype*ntypes-(jtype*(jtype+1)/2)+ itype));
			if (Threads::mod_partial[ipartial])
			{
				for (i=min_loc_ind[itype];i<=max_loc_ind[itype];i++)
				{
					//copying the local histograms
					ploc_source=local_count+(SimpleCfg::cumul[itype]+i)*sum_bins_offset+loc_offset[jtype];
					ploc_target=target.local_count+(SimpleCfg::cumul[itype]+i)*sum_bins_offset+loc_offset[jtype];
					for (ibin=0;ibin<max_loc_nbins;ibin++)
						*ploc_target++=*ploc_source++;
				}
			}
		}
	}
#endif

}

//----------------calculating the histogram----------------------------------------------	
void HistoSet::HistCalc( int &ntooclose)
{
	int j,iconst,isubconst;
	int *pcoordn;//pointer to the coordnumbs array
	int ip=0;
	int i,ithread;
	int *pcoordn_aux;
	longint *p_main, *p_aux;
	double *p_maind=0, *p_auxd=0;
	int th_offset;
#ifdef _TEST_MODE
	std::chrono::duration<double, std::milli> elapsed;
#endif
	th_offset=0;//to prevent to cause compiler warning

	thread_object.complete_count=0;//counter indicating how many thread finished with the calculation, reset it to 0
	
	if (nthreads>1)//signal to threads to start histogram calculation 
	{
		std::unique_lock<std::mutex> guard(thread_object.start1_mutex);//lock mutex
		thread_object.start_flag=1;//time to proceed
		thread_object.check_hist_start.notify_all();
	}
	
	cout<<"\nCalculating histogram, determining atoms below cutoffs, and checking coordination constraints, if there is any"<<endl;
	if (RunParams::potential > 0  && RunParams::potential<20)
		cout << "Calculating non-bonding potential" << endl;

#ifdef _TEST_MODE
	Threads::hctime1 = std::chrono::high_resolution_clock::now();
#endif

	HistCalcThread(thread_object.thread_arg[nthreads-1]);//starting histogram calculation for the main thread

#ifdef _TEST_MODE
	Threads::hctime2 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime2-Threads::hctime1; 
	Threads::dur_00=elapsed.count();
#endif
	
	//has to be done regardless it is multi-threaded compilation or not
	if (RunParams::potential>0)
	{
		//summing the main threads contribution
		for (ip=0;ip<FNC_POT::pot_dim;ip++)
		{
			FNC_POT::vdW_tot_pot+=FNC_POT::vdW_pot[ip];
			if (RunParams::potential == 1)
			{
				FNC_POT::Coulomb_tot_pot += FNC_POT::Coulomb_pot[ip];
				FNC_POT::vdW14_tot_pot += FNC_POT::vdW14_pot[ip];
				FNC_POT::Coulomb14_tot_pot += FNC_POT::Coulomb14_pot[ip];
			}
		}
	}

	//waiting for all threads to finish calculation of conf1 before saving the result to file
	if (nthreads>1)
	{
		std::unique_lock<std::mutex> guard(thread_object.count_mutex);//lock mutex
		thread_object.complete_count++;
		if (thread_object.complete_count==nthreads)//time to proceed
		{
			thread_object.start_flag=0;//reset it to 0
			thread_object.complete_count=0;//reset it to 0
			thread_object.check_count.notify_all();
		}
		else  thread_object.check_count.wait(guard);
	};

#ifdef _TEST_MODE
	Threads::hctime1 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime1-Threads::hctime2; 
	Threads::dur_01=elapsed.count();
#endif

	//Now update the histogram by adding the contributions of the auxiliary threads, and reset the auxiliary threads histogram segment 
	for (ithread=0;ithread<nthreads-1;ithread++)
	{
		p_main=pcounts;
		p_aux=thread_histcount_finder[ithread];

		for (i=0;i<ntot_bins*npartials;i++)
		{
			*p_main+=*p_aux;
			*p_aux++=0;//reset it to zero
			p_main++;
		}
		//total counts
		p_main=ptotal;
		p_aux=thread_object.thread_arg[ithread].hist_tot_start;
		for (i=0;i<npartials*ndiffbin;i++)
		{
			*p_main+=*p_aux;
			*p_aux++=0;//reset it to zero
			p_main++;
		}
#ifdef _USE_LOCAL_INV
		p_main=local_sum_count;
		p_aux=thread_lochistcount_finder[ithread];

		for (i=0;i<max_loc_nbins*npartials;i++)
		{
			*p_main+=*p_aux;
			*p_aux++=0;//reset it to zero
			p_main++;
		}
#endif
	}//end of ithread cycle

#ifdef _TEST_MODE
	Threads::hctime2 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime2-Threads::hctime1; 
	Threads::dur_02=elapsed.count();
#endif
	
	//Updating the potential array, if necessary
	if (RunParams::potential > 0)//both for 1 and 10
	{

		for (ithread = 0; ithread < nthreads - 1; ithread++)
		{
			p_maind = FNC_POT::vdW_pot;
			p_auxd = thread_object.thread_arg[ithread].vdW_pot_start;

			for (i = 0; i < FNC_POT::pot_dim; i++)
			{
				*p_maind += *p_auxd;
				FNC_POT::vdW_tot_pot += *p_auxd;//total vdW
				*p_auxd++ = 0;//reset it to zero
				p_maind++;
			}
		}//end of i cycle
	}
	if (RunParams::potential == 1)//only used for LJ
	{
		for (ithread=0;ithread<nthreads-1;ithread++)
		{
			p_maind=FNC_POT::Coulomb_pot;
			p_auxd=thread_object.thread_arg[ithread].Coulomb_pot_start;

			for (i=0;i<FNC_POT::pot_dim;i++)
			{
				*p_maind+=*p_auxd;
				FNC_POT::Coulomb_tot_pot+=*p_auxd;//total Coulomb
				*p_auxd++=0;//reset it to zero
				p_maind++;
			}
		}//end of i cycle
		if (Topology::npairs>0)//there is 1-4 interaction
		{
			for (ithread=0;ithread<nthreads-1;ithread++)
			{
				p_maind=FNC_POT::vdW14_pot;
				p_auxd=thread_object.thread_arg[ithread].vdW14_pot_start;

				for (i=0;i<FNC_POT::pot_dim;i++)
				{
					*p_maind+=*p_auxd;
					FNC_POT::vdW14_tot_pot+=*p_auxd;//total 1-4 vdW
					*p_auxd++=0;//reset it to zero
					p_maind++;
				}
			}//end of i cycle
			for (ithread=0;ithread<nthreads-1;ithread++)
			{
				p_maind=FNC_POT::Coulomb14_pot;
				p_auxd=thread_object.thread_arg[ithread].Coulomb14_pot_start;

				for (i=0;i<FNC_POT::pot_dim;i++)
				{
					*p_maind+=*p_auxd;
					FNC_POT::Coulomb14_tot_pot+=*p_auxd;//total 1-4 Coulomb
					*p_auxd++=0;//reset it to zero
					p_maind++;
				}
			}//end of i cycle
		}
	}

	//Updating the average coordination number constraint 
	int offset = (AvCoordConst::nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*AvCoordConst::neighbourcount)) / sizeof(*AvCoordConst::neighbourcount)) : 0);//number of dummy integer elements
	for (iconst=0;iconst<AvCoordConst::nconstraints;iconst++)
	{
		for (ithread=1;ithread<nthreads;ithread++)
		{
			//number of neighbours
			th_offset=ithread*(AvCoordConst::nconstraints+offset)+iconst;
			*(acoordc.neighbourcount+iconst)+=*(acoordc.neighbourcount+th_offset);
			*(acoordc.neighbourcount+th_offset)=0;
		}//end of if this constraint was affected by the move
	}//end of constraint cycle iconst

#ifdef _TEST_MODE
	Threads::hctime1 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime1-Threads::hctime2; 
	Threads::dur_03=elapsed.count();
#endif

	//Updating the coordnumbs array and calculating the number of atoms satisfying the coordination constraint
	//sets the pointer to the beginning of the coordnumbs array
	for (iconst=0;iconst<CoordNumbConst::nconstraints;iconst++)
	{
		//first add the contribution of the thread's coordnumbs array to the main
		for (ithread=0;ithread<nthreads-1;ithread++)
		{
			pcoordn=coordnc.finder[iconst];//sets the pointer to the first atom of this constraint in the coordnumbs array
			pcoordn_aux=coordnc.thread_finder[ithread]+CoordNumbConst::ncentral_cum[iconst];
			for (j=0;j<CoordNumbConst::ncentral[iconst];j++)
			{
				*pcoordn+=*pcoordn_aux;
				*pcoordn_aux++=0;
				pcoordn++;
			}
		}

		//now calculate the number of atoms satisfying the constraint
		pcoordn=coordnc.finder[iconst];
		for (isubconst=0;isubconst<coordnc.n_subconst[iconst];isubconst++)
			coordnc.nsatisfy[coordnc.cum_n_subconst[iconst]+isubconst]=0;
		for (j=0;j<CoordNumbConst::ncentral[iconst];j++)
		{
			for (isubconst=0;isubconst<coordnc.n_subconst[iconst];isubconst++)
			{
				if (*pcoordn==CoordNumbConst::target_coord[coordnc.cum_n_subconst[iconst]+isubconst])
					coordnc.nsatisfy[coordnc.cum_n_subconst[iconst]+isubconst]++;
			}
			pcoordn++;
		}
	}

#ifdef _ADVANCED_GEOM_CONST
	//Updating the valence array for bvs constraint
	//sets the pointer to the beginning of the valence array
	double *pval,*pval_aux;
	for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
	{
		//first add the contribution of the thread's valence array to the main
		for (ithread = 0; ithread < nthreads - 1; ithread++)
		{
			pval = bvs->finder[iconst];//sets the pointer to the first atom of this constraint in the coordnumbs array
			pval_aux = bvs->thread_finder[ithread] + BondValenceSumConst::ncentral_cum[iconst];
			for (j = 0; j < BondValenceSumConst::ncentral[iconst]; j++)
			{
				*pval += *pval_aux;
				*pval_aux++ = 0;
				pval++;
			}
		}
	} 
#endif
#ifdef _TEST_MODE
	Threads::hctime2 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime2-Threads::hctime1; 
	Threads::dur_04=elapsed.count();
#endif
	if (RunParams::auto_cutoff==0)
	{
		//Counting the atoms below cutoff
		ntooclose = 0;
		for (j = 0; j < conf.ntotal; j++)
		{
			if (tclarray[j] > 0)
				ntooclose++;
		}
		if (ntooclose > 0)
			cout << "\nNOTE(" << ++note << "): Number of atoms below cutoff: " << ntooclose << endl;
	}
	else
	{
		if (RunParams::fnc > 0 && RunParams::fnc != 4)//only for normal fnc
		{
		//the cutoff was determined, see, if there is fnc conflict
			for (i = 0; i < npartials; i++)
			{
				for (j = 0; j < (int)FNC_POT::fnc_partial_ind[i].size(); j++)//going through the fnc constraint this partial is involved in
				{
					if (fnc[0].dminsq[FNC_POT::fnc_partial_ind[i][j]] < rundat.pcutsq[i])
					{
						rundat.fnc_conflict = 1;
						//	cout << i << " " << sqrt(fnc[0].dminsq[FNC_POT::fnc_partial_ind[i][j]]) * SimpleCfg::boxedge << " " << sqrt(rundat.pcutsq[i]) * SimpleCfg::boxedge << endl;
						break;
					}
				}
			}
		}
	}
};


//-----------the calculation of the histogram performed by the threads---------------------
//all the other functions, which depend on the pair distances, as coordination constraints and potential will be done here
void HistoSet::HistCalcThread(ThreadArg &thread_arg)
{
	int type1, type2, i, j, ib, icoord, iconst, iavconst, inb, neightype;//cycle variable
	int jmax;//cycle upper limits
	int cc_neigh_type;//type of the neighbour atom for a CoordNumbConst
	int cc_neigh_offset;//starting point for the given constraint in CoordNumbConst arrays having tot_neightype elements 
	int thread_id;
	int min_ind, max_ind;//cycle limits
	int partialind, c_partind, pot_partialind;//index of the partial;pow(2,index of the partial),if virtual sites are present, then the host atom's type decides the partialind
	int bin_index;//index of the bin to put the distance in
	int tooclose_flag;//(1) means to add the too close atoms to the tooclose list, (0) not to in case of it is an fnc conflict pair
	int *p_fnc_neighind;//for the fnc neighbours or BONDS
	int *p_fnc_neighind2 = 0;//for the atoms in flexible ANGLE 
	int *p_excl_neighind = 0;//for the exclusions
	int  mytabpot_index = 0;//index among the given tabulated potentials for the presen partial
	longint *my_total;//pointer for this thread hist_tot_start
	double *coordi, *coordj;//pointers to the coordinates 
	double *cc_min, *cc_max;//for the actual central neighbour pair for each constarint
	double distance, dsquare;//distance, squared distance
	double dcomp[3];//the differences between the coordinates of two atoms
	int calcNB_flag = 0;//(1) means, that NB interactions should be calculated for the given pair, (0) means not, (-1) means 1-4 NB interaction should be calculated
	int GRtype1 = 0, GRtype2 = 0, GRpartialind = 0;//type of the atom according to the separate GROMACS types arrays
	double *my_vdWpot = 0, *my_Coulombpot = 0, *my_vdW14pot = 0, *my_Coulomb14pot = 0;//pointer for this thread potential
	double tabpot_rminsq = 0;//the square of the first tabulated r value for the given partial
	double tabcutoff_sq = 0;//the squared cutoff for the given data set
	double *pi = 0, *pj = 0, dcharge_c_square = 0;
	double *mypcutsq=nullptr;//for cutoff determination
	bool normal_calc;//if true then for the given partial normal calculation is performed, if not, then virtual site is involved 
					//in the partial, only potential calculation is performed
#ifdef _ADVANCED_GEOM_CONST
	int bvs_neigh_type;//type of the neighbour atom for bvs const
	int bvs_neigh_offset;//starting point for the given constraint in BondValenceSumConst arrays having tot_neightype elements 
	double *bvs_max;//cutoff for the current neighbour of each constraint
	double *bvs_R0, *bvs_b;//R0 and b for the current neighbour of each constraint
	double bv;
#endif	
#ifdef _NO_PERIODIC
	double distance_i=0, distance_j=0;
#endif

	//Each thread calculates for each partial, consecutive central atom segments are assigned to each thread,
	//and the given thread calculates for all the neighbours of its central atoms 
	//Calculate histogram, and determine atoms below cutoff, 
	//and calculate the number of neighbours, if there are coordination constraints

	thread_id = thread_arg.thread_index;//the ID of the thread
	my_total = thread_arg.hist_tot_start;//the starting points of the total counts for this thread
	my_vdWpot = thread_arg.vdW_pot_start;//the starting points of the vdW potential for this thread
	my_Coulombpot = thread_arg.Coulomb_pot_start;//the starting points of the Coulomb potential for this thread
	my_vdW14pot = thread_arg.vdW14_pot_start;//the starting points of the vdW potential for this thread
	my_Coulomb14pot = thread_arg.Coulomb14_pot_start;//the starting points of the Coulomb potential for this thread
	SetArraysize(&cc_min, RunParams::nicoord, "cc_min", "HistoSet::HistCalcThread");
	SetArraysize(&cc_max, RunParams::nicoord, "cc_max", "HistoSet::HistCalcThread");
#ifdef _ADVANCED_GEOM_CONST
	SetArraysize(&bvs_max, BondValenceSumConst::nconstraints, "bvs_max", "HistoSet::HistCalcThread");
	SetArraysize(&bvs_R0, BondValenceSumConst::nconstraints, "bvs_R0", "HistoSet::HistCalcThread");
	SetArraysize(&bvs_b, BondValenceSumConst::nconstraints, "bvs_b", "HistoSet::HistCalcThread");
#endif
#ifdef _USE_LOCAL_INV
	int loc_offset1=0;//starting point of the histogram for the central atom
#endif
	
	partialind = 0;
	c_partind = 1;
	if (RunParams::auto_cutoff & 1)
	{
		SetArraysize(&mypcutsq, npartials, "mypcutsq", "HistoSet::HistCalcThread");
		for (i = 0; i < npartials;i++)
			mypcutsq[i] = rundat.pcutsq[i];
	}
	for (type1 = 0; type1 < ntypes + nvirtualtypes; type1++)//central atom's type
	{
		for (type2 = type1; type2 < ntypes + nvirtualtypes; type2++)//neighbour atom's type 
		{
			if (type1 < ntypes && type2 < ntypes)//normal atoms
			{
				normal_calc = true;
				pot_partialind = partialind;
			}
			else
			{
				normal_calc = false;
				//for the virtual sites the host atom's type will decide the partialind

				cout << thread_id << " \tPartial potential for ";
				if (type1 >= ntypes)
				{
					i = SimpleCfg::virtual_host_type[type1 - ntypes];
					cout << "virtual type: " << type1 - ntypes;
				}
				else
				{
					cout << "type: " << type1;
					i = type1;
				}
				cout << " - ";
				if (type2 >= ntypes)
				{
					j = SimpleCfg::virtual_host_type[type2 - ntypes];
					cout << "virtual type: " << type2 - ntypes;
				}
				else
				{
					j = type2;
					cout << "type: " << type2;
				}

				pot_partialind = (i <= j ? (i * ntypes - (i * (i + 1) / 2) + j) : (j * ntypes - (j * (j + 1) / 2) + i));
				cout << " is calculated, contributing to potential partial " << i << ", " << j << endl;
			}

			if (normal_calc)
			{
				if (RunParams::potential == 10 && calc_tabpot[partialind])//is there is a tabulated potential for this partial
				{

					for (i = 0; i < RunParams::nused_potpartials; i++)
					{
						if (FNC_POT::tabpot_index[i] == partialind)
						{
							tabcutoff_sq = pow(FNC_POT::tab_cutoff[i], 2.0);//reduced, squared
							mytabpot_index = i;
							tabpot_rminsq = pow(FNC_POT::r_pot[i], 2);
							break;
						}
					}
				}

				//setting the actual minimal and maximal distance sqaures for each CoordNumbconst constraint
				for (iconst = 0; iconst < CoordNumbConst::nconstraints; iconst++)
				{
					if (CoordNumbConst::cctype[iconst] & c_partind)//there is a constraint for this type pair
					{
						if (type1 == CoordNumbConst::central[iconst])//first particle is central for the constraint
							cc_neigh_type = type2;//second is the neighbour
						else
							cc_neigh_type = type1;//first is the neighbour
						cc_neigh_offset = CoordNumbConst::cum_n_neightype[iconst];
						for (neightype = 0; neightype < CoordNumbConst::n_neightype[iconst]; neightype++)
						{
							if (cc_neigh_type == CoordNumbConst::neighbours[cc_neigh_offset + neightype])
							{
								//the type of the neighbour atom match the neightype-th neighbour type of this constraint
								cc_min[iconst] = CoordNumbConst::udminsq[cc_neigh_offset + neightype];
								cc_max[iconst] = CoordNumbConst::udmaxsq[cc_neigh_offset + neightype];
								break;
							}
						}
					}
					else
					{
						//the iconst constraint is not for this type pair, sets the minimum and maximum distance squares to unrealistic values
						cc_min[iconst] = 1e6;
						cc_max[iconst] = -1;//this is enough to make it impossible to find a distance square in range
					}
				}
#ifdef _ADVANCED_GEOM_CONST
				//setting the actual maximal distance sqaures and other parameters for each BVS constraint
				for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
				{
					if (BondValenceSumConst::bvsctype[iconst] & c_partind)//there is a constraint for this type pair
					{
						if (type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
							bvs_neigh_type = type2;//second is the neighbour
						else
							bvs_neigh_type = type1;//first is the neighbour
						bvs_neigh_offset = BondValenceSumConst::cum_n_neightype[iconst];
						for (neightype = 0; neightype < BondValenceSumConst::n_neightype[iconst]; neightype++)
						{
							if (bvs_neigh_type == BondValenceSumConst::neighbours[bvs_neigh_offset + neightype])
							{
								//the type of the neighbour atom match the neightype-th neighbour type of this constraint
								bvs_max[iconst] = BondValenceSumConst::maxsq[bvs_neigh_offset + neightype];
								bvs_R0[iconst]= BondValenceSumConst::R0[bvs_neigh_offset + neightype];
								bvs_b[iconst] = BondValenceSumConst::b[bvs_neigh_offset + neightype];
								break;
							}
						}
					}
					else
					{
						//the iconst constraint is not for this type pair, sets the minimum and maximum distance squares to unrealistic values
						bvs_max[iconst] = -1;//this is enough to make it impossible to find a distance square in range

					}
				}
#endif
			}
			if (type1 == type2)
			{
				min_ind = thread_arg.min_xx[type1];//the index of the first atom in its own type			
				max_ind = thread_arg.max_xx[type1];//the index of the last atom in its own type			
			}
			else
			{
				min_ind = thread_arg.min_atom_index[type1];//the index of the first atom in its own type
				max_ind = thread_arg.max_atom_index[type1];//the index of the last atom in its own type
			}
			//Sets the coordi pointer to the coordinates of the first central atom of type1 and conf1	
			coordi = conf.finder[type1] + 3 * min_ind;
			
			for (i = min_ind; i <= max_ind; i++)//Central atom: Cycling through the atoms of the first type
			{
				//sets the j loop upper limit
				if (type1 == type2)
					jmax = i;
				else
					jmax = conf.pnatoms[type2];



				//Sets the coordj pointer to the coordinates of the neighbour atom of type2 and conf2	
				coordj = conf.finder[type2];
				if (RunParams::potential == 1)//LJ
				{
					if (RunParams::fnc == 4)
						pi = conf.charge_gr_centre + 3 * fnc[0].charge_centre[conf.cumul[type1] + i];
					//find the GROMACS type (0 to nsep_GRtypes-1) of the atom 
					GRtype1 = 0;
					while (conf.cumul[type1] + i >= SimpleCfg::cumul_GR[GRtype1 + 1])
						(GRtype1)++;
					//GRtype1 now points to the type index of the seperate GROMACS type segment, where the atom is located (the atoms index is continuously increasing)
					//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
					GRtype1 = SimpleCfg::GRsep_GR_type[GRtype1];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters
				}
				if (normal_calc)
				{
#ifdef _USE_LOCAL_INV
					loc_offset1 = (conf.cumul[type1] + i) * sum_bins_offset + loc_offset[type2];
#endif
				}

#ifdef _NO_PERIODIC
				if (normal_calc)
				{
					dsquare = 0.0;
					for (icoord = 0; icoord < 3; icoord++)
						dsquare += pow(coordi[icoord], 2);
					distance_i = sqrt(dsquare);

					//calculate the bin index, where the central atom is located in the spherical sample 
					bin_index = int((distance_i - xmin[0]) * inv_binwidth[0]);
					atom_binind[conf.cumul[type1] + i] = bin_index;//the bin index centred on the origin of the sample for this atom
				}

#endif

				for (j = 0; j < jmax; j++)//Neighbour atom: Cycling through the atoms of the second type
				{
					//calculating the distance of the two atoms
					for (icoord = 0; icoord < 3; icoord++)
						dcomp[icoord] = coordi[icoord] - coordj[icoord];

					coordj += 3;

					GetMinImage(dcomp);//taking into account the periodical boundary conditions 
					dsquare = 0.0;
					for (icoord = 0; icoord < 3; icoord++)
						dsquare += pow(dcomp[icoord], 2);
					distance = sqrt(dsquare);
					if (RunParams::potential == 1)
					{
						calcNB_flag = 1;
						if (RunParams::fnc == 4)
						{
							//Decide, whether this is an excluded pair or not
							p_excl_neighind = fnc[0].vdW_exclusions + *(fnc[0].vdW_conv + conf.cumul[type1] + i);
							for (inb = 0; inb < *(fnc[0].vdW_excl_numb + conf.cumul[type1] + i); inb++)
							{
								//for each neighbour of this central atom
								if (abs(*p_excl_neighind) - 1 == conf.cumul[type2] + j)
								{
									if (*p_excl_neighind < 0)
									{
										calcNB_flag = -1;//calculate NB 14 interaction
										break;
									}
									else
									{
										calcNB_flag = 0;//not to calculate NB interaction
										break;//(inb cycle)
									}
								}
								p_excl_neighind++;//next neighbour index
							}//next excluded neighbour
						}

						if (abs(calcNB_flag) == 1)
						{
							//calculating the distances of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
							if (RunParams::fnc == 4)
							{
								pj = conf.charge_gr_centre + 3 * fnc[0].charge_centre[conf.cumul[type2] + j];
								for (icoord = 0; icoord < 3; icoord++)
									dcomp[icoord] = pi[icoord] - *pj++;
								GetMinImage(dcomp);	//taking into account the periodical boundary conditions 
								dcharge_c_square = 0;
								for (icoord = 0; icoord < 3; icoord++)
									dcharge_c_square += pow(dcomp[icoord], 2);

							}
							else
								dcharge_c_square = dsquare;//if there are no molecules, ->no charge centres, the normal atomic distance counts

							
							if (dcharge_c_square <= rundat.Coulomb_cutoff_sq)//check the cutoff
							{

								if (calcNB_flag > 0)
									*(my_Coulombpot + pot_partialind) += fnc[0].CalcCoulomb(conf.cumul[type1] + i, conf.cumul[type2] + j, dsquare);//calculating the Coulomb interaction
									
								else
									*(my_Coulomb14pot + pot_partialind) += RunParams::Coulomb14_fudge * fnc[0].CalcCoulomb(conf.cumul[type1] + i, conf.cumul[type2] + j, dsquare);//calculating the Coulomb interaction
							}

							if (dcharge_c_square <= rundat.vdW_cutoff_sq)//check the cutoff
							{
								//find the GROMACS type (0 to nsep_GRtypes-1) of the atom 
								GRtype2 = 0;
								while (conf.cumul[type2] + j >= SimpleCfg::cumul_GR[(GRtype2)+1])
									(GRtype2)++;
								//GRtype2 now points to the type index of the seperate GROMACS type segment where the atom is located (the atoms index is continuously increasing)

								//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
								GRtype2 = SimpleCfg::GRsep_GR_type[GRtype2];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters

								//GRpartialind is GROMACS type based partial index, needed to have the vdW parameters, but the histogram is split into RMC-type partials
								GRpartialind = (GRtype1 <= GRtype2 ? (GRtype1 * Topology::nGRtypes - (GRtype1 * (GRtype1 + 1) / 2) + GRtype2) : \
									(GRtype2 * Topology::nGRtypes - (GRtype2 * (GRtype2 + 1) / 2) + GRtype1));
								
								if (calcNB_flag > 0)
									*(my_vdWpot + pot_partialind) += (fnc[0].*CalcvdW)(GRpartialind, dsquare);//calculating the nonbonded interaction
								else
									*(my_vdW14pot + pot_partialind) += (fnc[0].*CalcvdW14)(GRpartialind, dsquare);//calculating the nonbonded interaction
								
							}
						}
					}//end of if potential==1
					if (!normal_calc)
						continue;//skip j cycle go to next second atom
					//checking, whether the distance is below cutoff
					
					if ((RunParams::auto_cutoff == 0 && dsquare < rundat.pcutsq[partialind]) || (RunParams::auto_cutoff == 3 && dsquare < mypcutsq[partialind]))
					{
						//if there is FNC do not included the FNC pair, for flexible molecules the auto_cutoff=-1 will not consider 1-2 bongd and 1-3 atoms of an angle

						//The index of the k-th atom of the l-th atom type, 
						//if k,l starts with 0: cumulative number of atoms for types 0 -> (l-1) + k

						tooclose_flag = 1;//the too close atoms are added to the too close list by default,
						//or a sofar smallest distance updates the cutoff
						if (rundat.fnc == 4 || rundat.fnc_conflict)//There is a conflict between the cutoffs and FNC

						{
							//In case of FNC conflict check whether this is an FNC pair, as FNC pairs will not be included into the
							//too close list, as we have no intention to move them out!
							//In case of rundat.fnc==4 check, whether these atoms belong to the same bond, as if they are they will not be considered for cutoff datermination
							p_fnc_neighind = fnc[0].neighbours + *(fnc[0].converter + conf.cumul[type1] + i);
							for (inb = 0; inb < *(fnc[0].numbneigh + conf.cumul[type1] + i); inb++)
							{
								//for each neighbour of this central atom
								if (*p_fnc_neighind + FNC_POT::index_offset == conf.cumul[type2] + j)
								{
									tooclose_flag = 0;//not to add to the too close list, do not set as cutoff
									//logfile << "fnc " << partialind<<" "<<conf.cumul[type1] + i << " " << conf.cumul[type2] + j << " "<<sqrt(dsquare)*RunParams::boxedge<<endl;
									break;//(inb cycle)
								}
								p_fnc_neighind++;//next neighbour index
							}//next fnc neighbour
						}//end of if conflict
						//check, if they are 1-3 atoms (two neighbour for the middle atom of the angle) in case of flexible ANGLE 
						//the 1-2 or 2-3 pairs are not checked, as it was checked during BOND checking
						if (rundat.fnc == 4 && tooclose_flag)
						{
							p_fnc_neighind = fnc[1].neighbours + *(fnc[1].converter + conf.cumul[type1] + i);
							p_fnc_neighind2 = fnc[1].neighbours2 + *(fnc[1].converter + conf.cumul[type1] + i);
							for (inb = 0; inb < *(fnc[1].numbneigh + conf.cumul[type1] + i); inb++)//for each neighbour of this central atom
							{
								if (*p_fnc_neighind < 0 && *p_fnc_neighind2+FNC_POT::index_offset == conf.cumul[type2] + j)//- sign in neighbour array shows that that atom is the middle in the angle
								{
									tooclose_flag = 0;//this is 1-3 atom in the angle, the potential should take care of it keeping them at a proper distance
									//logfile << "angle " << conf.cumul[type1] + i << " " << conf.cumul[type2] + j << " " << *p_fnc_neighind<<" "<<sqrt(dsquare)*SimpleCfg::boxedge << endl;
									break;//(inb cycle)
								}
								p_fnc_neighind++;//next neighbour index
								p_fnc_neighind2++;//next neighbour index
							}//next FNC neighbour
						}

						if (tooclose_flag)
						{
							if (RunParams::auto_cutoff == 0)
							{
								//as this cannot happen too often, all the threads are updating the same array, and mutex lock is used to
								//to ensure, that only one thread updates the arrays at a time
								std::unique_lock<std::mutex> guard(thread_object.tcl_mutex);//lock mutex

								//Increase the number of too close pairs of the central atom is involved in the tclarray
								tclarray[conf.cumul[type1] + i]++;
								//Increase the number of too close pairs of the neighbour atom is involved in the tclarray
								tclarray[conf.cumul[type2] + j]++;
#ifdef _TEST_MODE
								//logfile << "tooclose " << " t1 " << type1 << " t2 " << type2 << " " << conf.cumul[type1] + i << " " << conf.cumul[type2] + j << " " << sqrt(dsquare) * SimpleCfg::boxedge << endl;
								if (ntcp + 2 > maxtcp - 1)//would be over the array dimension, there are a lot of tooclose pairs 
									ResizeArray(&maxtcp, maxtcp + 100, &tcpairs, "tcpairs", "HistoSet::HistCalcChangeThread");
								tcpairs[ntcp] = conf.cumul[type1] + i;
								ntcp++;
								tcpairs[ntcp] = conf.cumul[type2] + j;
								ntcp++;
#endif
							}
							else
							{
								if (RunParams::auto_cutoff == 3)
								{
									//auto cutoff determination
									mypcutsq[partialind] = dsquare;
									//logfile << "min " << partialind<<" "<<conf.cumul[type1] + i << " " << conf.cumul[type2] + j <<  " "<<sqrt(dsquare)*SimpleCfg::boxedge<<endl;
								}
							}
						}
					}
					else if (RunParams::auto_cutoff == 1 && dsquare < mypcutsq[partialind])
					{
						//auto cutoff determination
						mypcutsq[partialind] = dsquare;
						//logfile << "minall " << partialind<<" "<<conf.cumul[type1] + i << " " << conf.cumul[type2] + j <<  " "<<sqrt(dsquare)*SimpleCfg::boxedge<<endl;
					}
#ifdef _NO_PERIODIC
					if (CoordNumbConst::nconstraints > 0 || AvCoordConst::nconstraints > 0)
						distance_j = sqrt(*(coordj - 3) * *(coordj - 3) + *(coordj - 2) * *(coordj - 2) + *(coordj - 1) * *(coordj - 1));
#endif		
					//Check coordination constraint, if there is any
					for (iconst = 0; iconst < CoordNumbConst::nconstraints; iconst++)
					{
						if (dsquare >= cc_min[iconst] && dsquare <= cc_max[iconst])
						{
							//the distance is in the range, the number of neighbours
							//has to be increased for central atom 
							if (type1 == type2)//both can be central and neighbour
							{
								*(thread_arg.coord_numb_finder[iconst] + i) += 1;
								*(thread_arg.coord_numb_finder[iconst] + j) += 1;
							}
							else
							{
								if (type1 == CoordNumbConst::central[iconst])//first particle is central for the constraint
									*(thread_arg.coord_numb_finder[iconst] + i) += 1;
								else //second particle is central for the constraint
									*(thread_arg.coord_numb_finder[iconst] + j) += 1;
							}//end of if (type1==type2)
						}//distance is in the range

					}//end of checking and updating coordination number constraints

					//Check average coordination constraint, if there is any
					for (iavconst = 0; iavconst < AvCoordConst::nconstraints; iavconst++)
					{
						if (AvCoordConst::avcctype[iavconst] == partialind)//there is a constraint for this atom pair
							if (dsquare >= AvCoordConst::udminsq[iavconst] && \
								dsquare <= AvCoordConst::udmaxsq[iavconst])
							{
								//the distance is in the range, the counter for this constraint
								//has to be increased for this configuration 
#ifdef _NO_PERIODIC
								if (distance_i + AvCoordConst::dav_red[iavconst] <= R0_red)
#endif
									*(thread_arg.av_coord_start + iavconst) += 1;
								if (type1 == type2)//second atom can be a central two
#ifdef _NO_PERIODIC
									if (distance_j + AvCoordConst::dav_red[iavconst] <= R0_red)
#endif
										*(thread_arg.av_coord_start + iavconst) += 1;

							}
					}//end of checking and updating av coordination number constraints
#ifdef _ADVANCED_GEOM_CONST
					//Check bond valence sum constraint, if there is any
					for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
					{
						if (dsquare <= bvs_max[iconst])
						{
							//the distance is in the range, calculate the bond valence
							bv = exp((bvs_R0[iconst] - distance*boxedge) / bvs_b[iconst]);
							if (type1 == type2)//both can be central and neighbour
							{
								*(thread_arg.valence_finder[iconst] + i) += bv;
								*(thread_arg.valence_finder[iconst] + j) += bv;
							}
							else
							{
								if (type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
									*(thread_arg.valence_finder[iconst] + i) += bv;
								else //second particle is central for the constraint
									*(thread_arg.valence_finder[iconst] + j) += bv;
							}//end of if (type1==type2)*/
						}//distance is in the range
					}//end of checking and updating bvs constraints
					
#endif
					
					if (RunParams::potential == 10 && calc_tabpot[partialind])
					{

						//check whether the distance is in the range of the tabulated potential for this partial
						if (dsquare >= tabpot_rminsq && dsquare < tabcutoff_sq)
						{
							*(my_vdWpot + mytabpot_index) += (fnc[0].*CalcvdW)(mytabpot_index, dsquare);

						}

					}

					
					//calculate the bin index, where this distance should be put
					for (ib = 0; ib < ndiffbin; ib++)
					{
						if (distance >= xmin[ib] && distance < xmax[ib])
						{
							bin_index = int((distance - xmin[ib]) * inv_binwidth[ib]);
							*(thread_arg.hist_count_finder[ib * npartials + partialind] + bin_index) += 1;//count for this bin
							*(my_total + ib * npartials + partialind) += 1;//total count for this partial
						}
					}

#ifdef _USE_LOCAL_INV
					//local histogram has to be stored
					//as this cannot happen too often, all the threads are updating the same array, and mutex lock is used to
					//to ensure, that only one thread updates the arrays at a time
					//only for the update of the neighbour atoms histogram mutex has to be used, as the central atoms are distributed amond the threads
					//to decrease competition among the threads for this mutex lock, a bin-dependent mutex is used
					std::unique_lock<std::mutex> guard(thread_object.inv_mutex);//lock mutex

					if (distance >= min_loc_r && distance <= max_loc_r)
					{
						bin_index = int((distance - min_loc_r) * inv_locbinwidth);
						*(local_count + loc_offset1 + bin_index) += 1;
						*(local_count + (conf.cumul[type2] + j) * sum_bins_offset + loc_offset[type1] + bin_index) += 1;
						*(thread_arg.hist_loc_count_finder[partialind] + bin_index) += 2;//sum. count for this bin
					}
#endif
				}//end of neighbour atom indices of second atom type (cycle j)
				coordi += 3;
			}//end of central atom indices of first atom type (cycle i)
			
			if (normal_calc)
			{
				//write on screen the number of binned distances
				for (ib = 0; ib < ndiffbin; ib++)
					cout << thread_id << " Number of binned distances for " << ib << " bin size, " << type1 << "," << type2 << " partial \t:" << *(my_total + ib * npartials + partialind) << endl;

				//Incrementing the index of the partial 
				partialind++;
				c_partind <<= 1;//for the coordinaion number and bvs constraint
			}

		}//end of neighbour atom type (cycle type2)
	}//end of central atom type (cycle type1)
	if (RunParams::auto_cutoff & 1)
	{
		//update the pcutsq array
		std::unique_lock<std::mutex> guard(thread_object.cut_mutex);//lock mutex
		for (partialind = 0; partialind < npartials; partialind++)
			if (mypcutsq[partialind] < rundat.pcutsq[partialind])
				rundat.pcutsq[partialind] = mypcutsq[partialind];
		delete[] mypcutsq;
	}
	delete[] cc_min;
	delete[] cc_max;
#ifdef _ADVANCED_GEOM_CONST
	delete[] bvs_max;
	delete[] bvs_b;
	delete[] bvs_R0;
#endif
	
}



//----Calculating the change in the histograms due to the move---
//All the checking and the too close list updates were done earlier by NeighboutList::CheckCutoff!
void HistoSet::HistCalcChange()
{
	int i, j, ib, ithread, itype, ipartial;
	int cc_mod=0;
	int partialind = 0, c_partind, pot_partialind;//index of the partial;pow(2,index of the partial),where to put the pot contribution of virtual sites
	int bin_index = 0, bin_index_old = 0;//index of the histogram bin
	int th_offset;
	int cc_neigh_type;//type of the neighbour atom for a CoordNumbConst
	int cc_neigh_offset;//starting point for the given constraint in CoordNumbConst arrays having tot_neightype elements 
	int padding_offset;//number of dummy integer elements in CoordNumbConst::nsatisfy array for cache padding
	int central_index, neighbour_index;//index of the central and neighbour atoms 
	int imoved, jmoved, iconst, iavconst;
	int neightype, *pcoordn, *pcoord_aux;//for calculating the number of atoms satisfying he CoordNumbConst
	int *p_indexm1, *p_indexm2, *p_type1, *p_type2;
	int  mytabpot_index = 0;//index among the given tabulated potentials for the presen partial

	longint *p_main, *p_aux;
	//	longint start,finish;//to check the speed of program parts

	double dsquare, dsquare_old, distance, distance_old;
	double *p_newcoord1, *p_newcoord2, *p_oldcoord1, *p_oldcoord2;
	double dx, dy, dz;//coordinate difference
	double tabpot_rminsq = 0;//the square of the first tabulated r value for the given partial
	double tabcutoff_sq = 0;//the squared cutoff for the given data set
	double *cc_min, *cc_max;//for the actual central neighbour pair for each constarint
	//double duration;//to see time differences
#ifdef _ADVANCED_GEOM_CONST
	int bvs_neigh_type;//type of the neighbour atom for bvs const
	int bvs_neigh_offset;//starting point for the given constraint in BondValenceSumConst arrays having tot_neightype elements 
	double *bvs_max;
	double *bvs_R0, *bvs_b;//R0 and b for the current neighbour of each constraint
	double bv;
#endif
#ifdef _USE_LOCAL_INV
	int loc_offset1;
#endif
	int igr = 0, jgr;//charge group centre indices
	int GRtype1 = 0, GRtype2 = 0, GRpartialind = 0;//type of the atom according to the separate GROMACS types arrays
	double *my_vdWpot = 0, *my_Coulombpot = 0, *my_s_vdWpot = 0, *my_s_Coulombpot = 0, *my_l_vdWpot = 0, *my_l_Coulombpot = 0;//pointer for this thread potential
	double *pi = 0, *pj = 0;//pointer to the charge group centres
	double *p_oldchargegr_c1 = 0, *p_oldchargegr_c2 = 0;//pointers to the old charge group centre coordinates of the moved atoms
	double dcharge_c_sq = 0, dcharge_c_sq_old = 0;
	my_vdWpot = thread_object.thread_arg[nthreads - 1].vdW_pot_start;//the starting point of the vdW potential for this thread
	my_Coulombpot = thread_object.thread_arg[nthreads - 1].Coulomb_pot_start;//the starting point of the Coulomb potential for this thread
	my_s_vdWpot = thread_object.thread_arg[nthreads - 1].vdW_pot_s_start;//the starting point of the vdW potential for this thread
	my_s_Coulombpot = thread_object.thread_arg[nthreads - 1].Coulomb_pot_s_start;//the starting point of the Coulomb potential for this thread
	my_l_vdWpot = thread_object.thread_arg[nthreads - 1].vdW_pot_l_start;//the starting point of the vdW potential for this thread
	my_l_Coulombpot = thread_object.thread_arg[nthreads - 1].Coulomb_pot_l_start;//the starting point of the Coulomb potential for this thread
	bool normal_calc;//true if between atom type, false, if any of it a virtual type
	SetArraysize(&cc_min, RunParams::nicoord, "cc_min", "HistoSet::HistCalcChange");
	SetArraysize(&cc_max, RunParams::nicoord, "cc_max", "HistoSet::HistCalcChange");
#ifdef _ADVANCED_GEOM_CONST
	SetArraysize(&bvs_max, BondValenceSumConst::nconstraints, "bvs_max", "HistoSet::HistCalcChange");
	SetArraysize(&bvs_R0, BondValenceSumConst::nconstraints, "bvs_R0", "HistoSet::HistCalcChange");
	SetArraysize(&bvs_b, BondValenceSumConst::nconstraints, "bvs_b", "HistoSet::HistCalcChange");
#endif
#ifdef _TEST_MODE
	std::chrono::duration<double, std::milli> elapsed;
#endif
#ifdef _NO_PERIODIC
	double distance_i, distance_j, distance_o_i, distance_o_j;
#endif
	if (RunParams::potential > 0)
	{
		for (ithread = 0; ithread < nthreads; ithread++)
		{
			for (ipartial = 0; ipartial < FNC_POT::pot_dim; ipartial++)
			{
				thread_object.thread_arg[ithread].vdW_pot_s_start[ipartial] = 0;
				thread_object.thread_arg[ithread].vdW_pot_l_start[ipartial] = 0;
				thread_object.thread_arg[ithread].Coulomb_pot_s_start[ipartial] = 0;
				thread_object.thread_arg[ithread].Coulomb_pot_l_start[ipartial] = 0;
			}
		}
	}
	//------MULTI THREADING BEGINS-----
	//The main thread will handle the calculation between the moved atoms, and its portion of the calculation between the moved
	//atoms and other atoms
	
	//Determine, whether the cutoffs, and the coordination constraints are satisfied
	//In case of the moveout option, if the moved was among the "too close" atoms, checking
	//whether it is still below cutoff, or it can be removed from the list. The too close atom can be any of the
	//moved atoms!
	//Histogram change is calculated as well

	if (nthreads>1)//signal to threads to start histogram change calculation 
	{
		//wait, if not all the new threads finished with the previous loop
		std::unique_lock<std::mutex> guard(thread_object.count_mutex);//lock mutex
		thread_object.start_count++;
#ifdef _TH_ORDER
		cout << "x" << endl;
		logfile << "HCC " << thread_object.start_count << endl;
#endif
		if (thread_object.start_count==nthreads)//time to proceed, this was the last thread to reach this point
		{
			acceptable=1;//the move is acceptable, if it does not prove otherwise
			Threads::is_E0_shift = 0;
			Threads::is_IQ_mucorr = 0;
			Threads::is_fprime_shift = 0;
#ifdef _AENET
			Threads::calc_ANN=Aenet::calc_ANN;
			//split has to be done here after all the threads finished with updates
			if (Threads::calc_ANN==1)//calculated in every step, use aenet_changed_neigh_ind
				thread_object.SplitAenetAtoms(Aenet::nchanged,Aenet::aenet_changed_neigh_ind);//split the neighbour atoms for which the ANN pot should be calculated among the threads
			else if (Aenet::calc_ANN==2)
				thread_object.SplitAenetAtoms(Aenet::nstored,Aenet::aenet_stored_ind);//use the stored indices, al the atoms involved in change since last calc
			
#endif
			for (i = 0; i < ExptsData::nek; i++)
				Threads::EXAFS_gridind[i] = ExptsData::ek_gridind[i];
			for (i = 0; i < ExptsData::nfq; i++)
			{
				Threads::IQ_muact[i] = ExptsData::fqmuact[i];
				if (ExptsData::fqAXS[i] > 0)
					Threads::IQ_fprimeact[i] = ExptsData::fqfprimeact[i * ntypes + ExptsData::fqAXS[i] - 1];
			}
			thread_object.start_flag=0;//reset it for later use
			thread_object.start_count=0;//reset it to 0
			if (!RunParams::custmove && ntypes>1)//for custmove or ntypes=1 no update is necessary, as all the partials and types are modified
			{
				for(i=0;i<npartials;i++)
					Threads::mod_partial[i]=0;
			
				p_type1=move->types;//points to the type of the first moved atom
				for (i=0;i<Move::tot_moved_atoms;i++)
				{
					//Determining the index of the partial 
					for(itype=0;itype<ntypes;itype++)//go through all the partials containing this type
					{
						partialind=(*p_type1<=itype ? (*p_type1 * ntypes-(*p_type1 *((*p_type1)+1)/2)+ itype) : \
								(itype*ntypes-(itype*(itype+1)/2)+ *p_type1));
						//for (ib=0;ib<ndiffbin;ib++)
							Threads::mod_partial[partialind]=1;//this partial was modified
					}
					p_type1++;
				}	
			}
			thread_object.check_count.notify_all();
#ifdef _TH_ORDER
			cout << ";" << endl;
#endif

		}
		else thread_object.check_count.wait(guard);
	}
	else
	{
		//Compiled as MULTI, but used only with nthreads==1
		acceptable=1;//the move is acceptable, if it does not prove otherwise
		for (i = 0; i < ExptsData::nek; i++)
			Threads::EXAFS_gridind[i] = ExptsData::ek_gridind[i];
		for (i = 0; i < ExptsData::nfq; i++)
		{
			Threads::IQ_muact[i] = ExptsData::fqmuact[i];
			if (ExptsData::fqAXS[i] > 0)
				Threads::IQ_fprimeact[i] = ExptsData::fqfprimeact[i * ntypes + ExptsData::fqAXS[i] - 1];
		}
		Threads::is_E0_shift = 0;
		Threads::is_IQ_mucorr = 0;
		Threads::is_fprime_shift = 0;
#ifdef _AENET
			Threads::calc_ANN=Aenet::calc_ANN;
			//split has to be done here to be compatible with multithreading
			if (Threads::calc_ANN==1)//calculated in every step, use aenet_changed_neigh_ind
				thread_object.SplitAenetAtoms(Aenet::nchanged,Aenet::aenet_changed_neigh_ind);//split the neighbour atoms for which the ANN pot should be calculated among the threads
			else if (Aenet::calc_ANN==2)
				thread_object.SplitAenetAtoms(Aenet::nstored,Aenet::aenet_stored_ind);//use the stored indices, al the atoms involved in change since last calc
			
#endif
		if (!RunParams::custmove && ntypes > 1)//if ntypes=1 or custom move no update is necessary, as all the partials and types are modified
		{

			for(i=0;i<npartials;i++)
				Threads::mod_partial[i]=0;
		
			p_type1=move->types;//points to the type of the first moved atom
			for (i=0;i<Move::tot_moved_atoms;i++)
			{
				//Determining the index of the partial 
				for(itype=0;itype<ntypes;itype++)//go through all the partials containing this type
				{
					partialind=(*p_type1<=itype ? (*p_type1 * ntypes-(*p_type1 *((*p_type1)+1)/2)+ itype) : \
							(itype*ntypes-(itype*(itype+1)/2)+ *p_type1));
					//for (ib=0;ib<ndiffbin;ib++)
						Threads::mod_partial[partialind]=1;//this partial was modified
				}
				p_type1++;
			}	
		}
	}


#ifdef _TEST_MODE
	Threads::hctime1 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime1-Threads::hctime2; 
	Threads::dur_00a+=elapsed.count();
#endif

	//Only by the main thread
	//Between the moved atoms, if there is only one, this is skipped
	p_indexm1=move->indices;//sets the pointer to the index of the first moved atom
	p_type1=move->types;//sets the pointer to the type of the first moved atom
	p_newcoord1=move->newpos;//sets the pointer to the new coordinates of the first moved atom
	p_oldcoord1=move->oldpos;//sets the pointer to the old coordinates of the first moved atom
	p_oldchargegr_c1=move->old_charge_c_pos;//sets the pointer to the old coordinates of the first moved atom's charge group centre coordinates

	for (imoved=0;imoved<Move::tot_moved_atoms+Move::tot_moved_virtuals-1;imoved++)
	{
		//index of the central atom within its own type
		central_index=*p_indexm1-conf.cumul[*p_type1];

		p_indexm2=p_indexm1+1;//sets the pointer to the index of the second moved atom
		p_type2=p_type1+1;//sets the pointer to the type of the second moved atom
		p_newcoord2=p_newcoord1+3;//sets the pointer to the new coordinates of the second moved atom
		p_oldcoord2=p_oldcoord1+3;//sets the pointer to the old coordinates of the second moved atom
		if (RunParams::potential==1)
		{
			if (RunParams::fnc==4)
			{
				p_oldchargegr_c2=p_oldchargegr_c1+3;//sets the pointer to the old coordinates of the first moved atom's charge group centre coordinates
				igr=fnc[0].charge_centre[*p_indexm1];//charge group centre of the moved atom
				pi=conf.charge_gr_centre+3*igr;
			}
			//find the GROMACS type (0 to nsep_GRtypes-1) of the atom 
			GRtype1=0;
			while (*p_indexm1>=SimpleCfg::cumul_GR[(GRtype1) +1])
				(GRtype1)++;
			//GRtype1 now points to the type index of the seperate GROMACS type segment, where the atom is located (the atoms index is continuously increasing)
			//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
			GRtype1=SimpleCfg::GRsep_GR_type[GRtype1];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters
		}
#ifdef _NO_PERIODIC
		distance_i=sqrt( pow(*p_newcoord1,2) + pow(*(p_newcoord1+1),2) + pow(*(p_newcoord1+2),2));
		distance_o_i=sqrt( pow(*p_oldcoord1,2) + pow(*(p_oldcoord1+1),2) + pow(*(p_oldcoord1+2),2));
#endif
		for (jmoved=imoved+1;jmoved<Move::tot_moved_atoms + Move::tot_moved_virtuals;jmoved++)
		{
#ifdef _USE_LOCAL_INV
			loc_offset1=(*p_indexm1)*sum_bins_offset+loc_offset[*p_type2];
#endif
	
			//Determining the index of the partial 
			partialind = (*p_type1 <= *p_type2 ? (*p_type1 * ntypes - (*p_type1 * (*p_type1 + 1) / 2) + *p_type2) : \
				(*p_type2 * ntypes - (*p_type2 * (*p_type2 + 1) / 2) + *p_type1));
			c_partind = (int)pow((double)2, partialind);//for the coordination number and bvsconstraint
			if (*p_type1 < ntypes && *p_type2 < ntypes)//normal atoms
			{
				normal_calc = true;
				pot_partialind = partialind;
			}
			else
			{
				normal_calc = false;
				//for the virtual sites the host atom's type will decide the partialind

				if (*p_type1 >= ntypes)
					i = SimpleCfg::virtual_host_type[*p_type1 - ntypes];
				else
					i = *p_type1;
				if (*p_type2 >= ntypes)
					j = SimpleCfg::virtual_host_type[*p_type2 - ntypes];
				else
					j = *p_type2;
			
				pot_partialind = (i <= j ? (i * ntypes - (i * (i + 1) / 2) + j) : (j * ntypes - (j * (j + 1) / 2) + i));
			}
			
			if (normal_calc)
			{
				//setting the actual minimal and maximal distance sqaures for each CoordNumbconst constraint
				if (RunParams::potential == 10 && calc_tabpot[partialind])//is there is a tabulated potential for this partial
				{
					for (i = 0; i < RunParams::nused_potpartials; i++)
					{
						if (FNC_POT::tabpot_index[i] == partialind)
						{
							tabcutoff_sq = pow(FNC_POT::tab_cutoff[i], 2.0);//reduced, squared
							mytabpot_index = i;
							tabpot_rminsq = pow(FNC_POT::r_pot[i], 2);
							break;
						}
					}
				}
				for (iconst = 0; iconst < CoordNumbConst::nconstraints; iconst++)
				{
					if (CoordNumbConst::cctype[iconst] & c_partind)//there is a constraint for this type pair
					{
						if (*p_type1 == CoordNumbConst::central[iconst])//first moved particle is central for the constraint
							cc_neigh_type = *p_type2;//second is the neighbour
						else
							cc_neigh_type = *p_type1;//first is the neighbour
						cc_neigh_offset = CoordNumbConst::cum_n_neightype[iconst];
						for (neightype = 0; neightype < CoordNumbConst::n_neightype[iconst]; neightype++)
						{
							if (cc_neigh_type == CoordNumbConst::neighbours[cc_neigh_offset + neightype])
							{
								//the type of the neighbour atom match the neightype-th neighbour type of this constraint
								cc_min[iconst] = CoordNumbConst::udminsq[cc_neigh_offset + neightype];
								cc_max[iconst] = CoordNumbConst::udmaxsq[cc_neigh_offset + neightype];
								break;
							}
						}
					}
					else
					{
						//the iconst constraint is not for this type pair, sets the minimum and maximum distance squares to unrealistic values
						cc_min[iconst] = 1e6;
						cc_max[iconst] = -1;//this is enough to make it impossible to find a distance square in range
					}
				}
			}
#ifdef _ADVANCED_GEOM_CONST
			//setting the actual maximal distance sqaures and other parameters for each BVS constraint
			for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
			{
				if (BondValenceSumConst::bvsctype[iconst] & c_partind)//there is a constraint for this type pair
				{
					if (*p_type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
						bvs_neigh_type = *p_type2;//second is the neighbour
					else
						bvs_neigh_type = *p_type1;//first is the neighbour
					bvs_neigh_offset = BondValenceSumConst::cum_n_neightype[iconst];
					for (neightype = 0; neightype < BondValenceSumConst::n_neightype[iconst]; neightype++)
					{
						if (bvs_neigh_type == BondValenceSumConst::neighbours[bvs_neigh_offset + neightype])
						{
							//the type of the neighbour atom match the neightype-th neighbour type of this constraint
							bvs_max[iconst] = BondValenceSumConst::maxsq[bvs_neigh_offset + neightype];
							bvs_R0[iconst] = BondValenceSumConst::R0[bvs_neigh_offset + neightype];
							bvs_b[iconst] = BondValenceSumConst::b[bvs_neigh_offset + neightype];
							break;
						}
					}
				}
				else
				{
					//the iconst constraint is not for this type pair, sets the minimum and maximum distance squares to unrealistic values
					bvs_max[iconst] = -1;//this is enough to make it impossible to find a distance square in range

				}
			}
#endif
			//calculating the distance of the two atoms
			dx=*p_newcoord1-*p_newcoord2++;
			dy=*(p_newcoord1+1)-*p_newcoord2++;
			dz=*(p_newcoord1+2)-*p_newcoord2++;//goes automatically to the next moved atom

			//taking into account the periodical boundary conditions 			
			if (dx>1) dx-=2;
			else if (dx<-1) dx+=2;
			if (dy>1) dy-=2;
			else if (dy<-1) dy+=2;
			if (dz>1) dz-=2;
			else if (dz<-1) dz+=2;
			dsquare=dx*dx+dy*dy+dz*dz;

			//old distance
			dx=*p_oldcoord1-*p_oldcoord2++;
			dy=*(p_oldcoord1+1)-*p_oldcoord2++;
			dz=*(p_oldcoord1+2)-*p_oldcoord2++;//goes automatically to the next moved atom

			//taking into account the periodical boundary conditions 			
			if (dx>1) dx-=2;
			else if (dx<-1) dx+=2;
			if (dy>1) dy-=2;
			else if (dy<-1) dy+=2;
			if (dz>1) dz-=2;
			else if (dz<-1) dz+=2;
			dsquare_old=dx*dx+dy*dy+dz*dz;

#ifdef _NO_PERIODIC
			distance_j=sqrt( pow(*(p_newcoord2-3),2) + pow(*(p_newcoord2-2),2) + pow(*(p_newcoord2-1),2));
			distance_o_j=sqrt( pow(*(p_oldcoord2-3),2) + pow(*(p_oldcoord2-2),2) + pow(*(p_oldcoord2-1),2));
#endif
			
			distance=sqrt(dsquare);
			//updating the LJ potential if it is any
			if (RunParams::potential==1)
			{
				if (RunParams::fnc==4)
				{
					//calculating the distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
					jgr=fnc[0].charge_centre[*p_indexm2];//charge grouip centre of the 2nd moved atom
					dcharge_c_sq=0;
					dcharge_c_sq_old=0;//this will stay, if the two group is the same
					if (igr != jgr)
					{
						pj = conf.charge_gr_centre + 3 * jgr;
						dx = *pi - *pj++;
						dy = *(pi + 1) - *pj++;
						dz = *(pi + 2) - *pj++;

						//taking into account the periodical boundary conditions 			
						if (dx > 1) dx -= 2;
						else if (dx < -1) dx += 2;
						if (dy > 1) dy -= 2;
						else if (dy < -1) dy += 2;
						if (dz > 1) dz -= 2;
						else if (dz < -1) dz += 2;
						dcharge_c_sq = dx * dx + dy * dy + dz * dz;

						//calculating the old distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
						dx = *(p_oldchargegr_c1)-*p_oldchargegr_c2++;
						dy = *(p_oldchargegr_c1 + 1) - *p_oldchargegr_c2++;
						dz = *(p_oldchargegr_c1 + 2) - *p_oldchargegr_c2++;

						//taking into account the periodical boundary conditions 			
						if (dx > 1) dx -= 2;
						else if (dx < -1) dx += 2;
						if (dy > 1) dy -= 2;
						else if (dy < -1) dy += 2;
						if (dz > 1) dz -= 2;
						else if (dz < -1) dz += 2;
						dcharge_c_sq_old = dx * dx + dy * dy + dz * dz;
					}
					else
						p_oldchargegr_c2 += 3;//go to next moved's data
				}
				else
				{
					//no molecules, no charge groups, the atomic distances has to be used
					dcharge_c_sq_old=dsquare_old;
					dcharge_c_sq=dsquare;
				}
				
				if (dcharge_c_sq <= rundat.Coulomb_cutoff_sq)//check the cutoff
				{
					double temp = (fnc[0].CalcCoulomb)(*p_indexm1, *p_indexm2, dsquare);//calculating the potential for the new distance
					if (fabs(temp) > POT_SPLIT_HIGH)
						*(my_l_Coulombpot + pot_partialind) += temp;//collect the large values separately
					else if (fabs(temp) < POT_SPLIT_LOW)
						*(my_s_Coulombpot + pot_partialind) += temp;//collect the small values separately
					else
						*(my_Coulombpot + pot_partialind) += temp;//normal
				}
				
				
				if (dcharge_c_sq_old <= rundat.Coulomb_cutoff_sq)//check the cutoff
				{
					double temp = (fnc[0].CalcCoulomb)(*p_indexm1, *p_indexm2, dsquare_old);//calculating the potential for the old distance
					if (fabs(temp) > POT_SPLIT_HIGH)
						*(my_l_Coulombpot + pot_partialind) -= temp;//collect the large values separately
					else if (fabs(temp) < POT_SPLIT_LOW)
						*(my_s_Coulombpot + pot_partialind) -= temp;//collect the small values separately
					else
						*(my_Coulombpot + pot_partialind) -= temp;//normal
				}
				if (dcharge_c_sq<=rundat.vdW_cutoff_sq || dcharge_c_sq_old<=rundat.vdW_cutoff_sq)//check the cutoff
				{
					//find the GROMACS type (0 to nsep_GRtypes-1) of the atom in the first and last instance of the molecule type
					GRtype2=0;
					while (*p_indexm2>=SimpleCfg::cumul_GR[(GRtype2) +1])
						(GRtype2)++;
					//GRtype2 now points to the type index of the seperate GROMACS type segmentm where the atom is located (the atoms index is continuously increasing)
					//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
					GRtype2=SimpleCfg::GRsep_GR_type[GRtype2];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters

					//GRpartialind is GROMACS type based partial index, needed to have the vdW parameters, but the histogram is split into RMC-type partials
					GRpartialind=(GRtype1<=GRtype2 ? (GRtype1*Topology::nGRtypes-(GRtype1*(GRtype1+1)/2)+GRtype2) : \
						(GRtype2*Topology::nGRtypes-(GRtype2*(GRtype2+1)/2)+GRtype1));
					
					if (dcharge_c_sq <= rundat.vdW_cutoff_sq)
					{
						double temp = (fnc[0].*CalcvdW)(GRpartialind, dsquare);//calculating the potential for the new distance
						if (fabs(temp) >POT_SPLIT_HIGH)
							*(my_l_vdWpot + pot_partialind) += temp;//collect the large values separately
						else if (fabs(temp) <POT_SPLIT_LOW)
							*(my_s_vdWpot + pot_partialind) += temp;//collect the small values separately
						else
							*(my_vdWpot + pot_partialind) += temp;//normal
						//if (fabs(temp) > POT_WARNING)
							//logfile << "HCC+ " << " " << pot_partialind << " " << *p_indexm1<<" "<< *p_indexm2 << " " << dsquare << " " << temp << " " << *(my_vdWpot + pot_partialind) << " l "<< *(my_l_vdWpot + pot_partialind) << " s "<< * (my_s_vdWpot + pot_partialind) <<endl;
					}
					if (dcharge_c_sq_old <= rundat.vdW_cutoff_sq)
					{
						double temp= (fnc[0].*CalcvdW)(GRpartialind, dsquare_old);//calculating the potential for the old distance
						if (fabs(temp) > POT_SPLIT_HIGH)
							*(my_l_vdWpot + pot_partialind) -= temp;//collect the large values separately
						else if (fabs(temp) < POT_SPLIT_LOW)
							*(my_s_vdWpot + pot_partialind) -= temp;//collect the small values separately
						else
							*(my_vdWpot + pot_partialind) -= temp;//normal
						//if (fabs(temp) > POT_WARNING)
							//logfile << "HCC- " << " " << pot_partialind << " " << *p_indexm1 << " " << *p_indexm2 << " " << dsquare_old << " " << temp << " " << *(my_vdWpot + pot_partialind) << " l " << *(my_l_vdWpot + pot_partialind) << " s " << *(my_s_vdWpot + pot_partialind)<< endl;
					}
				}
			}
			if (!normal_calc)
			{
				p_indexm2++;//sets the pointer to the index of the next second moved atom
				p_type2++;//sets the pointer to the type of the next second moved atom
				continue;
			}
			if (RunParams::potential == 10 && calc_tabpot[partialind])
			{
				//check whether the distance is in the range of the tabulated potential for this partial
				if (dsquare >= tabpot_rminsq && dsquare < tabcutoff_sq)
				{
					double temp= (fnc[0].*CalcvdW)(mytabpot_index, dsquare);
					if (fabs(temp) > POT_SPLIT_HIGH)
						*(my_l_vdWpot + mytabpot_index) += temp;
					else if (fabs(temp) < POT_SPLIT_LOW)
						*(my_s_vdWpot + mytabpot_index) += temp;
					else
						*(my_vdWpot + mytabpot_index) += temp;
				}
				if (dsquare_old >= tabpot_rminsq && dsquare_old < tabcutoff_sq)
				{
					double temp = (fnc[0].*CalcvdW)(mytabpot_index, dsquare_old);
					if (fabs(temp) > POT_SPLIT_HIGH)
						*(my_l_vdWpot + mytabpot_index) -= temp;
					else if (fabs(temp) < POT_SPLIT_LOW)
						*(my_s_vdWpot + mytabpot_index) -= temp;
					else
						*(my_vdWpot + mytabpot_index) -= temp;
				}
			}

			//calculate the bin index, where this distance should be put
			for (ib = 0; ib < ndiffbin; ib++)
			{
				if (distance >= xmin[ib] && distance <= xmax[ib])
				{
					bin_index = int((distance - xmin[ib])*inv_binwidth[ib]);
					*(finder[ib*npartials + partialind] + bin_index) += 1;//count for this bin
					*(ptotal + ib * npartials + partialind) += 1;//total count for this partial
				}
			}
			//Removing the contribution of the old distances
			//old coordinates of the first moved atom
			distance_old=sqrt(dsquare_old);

			//calculate the bin index, where this distance should be put
			for (ib = 0; ib < ndiffbin; ib++)
			{
				if (distance_old >= xmin[ib] && distance_old <= xmax[ib])
				{
					bin_index_old = int((distance_old - xmin[ib])*inv_binwidth[ib]);
					*(finder[ib*npartials + partialind] + bin_index_old) -= 1;//count for this bin
					*(ptotal + ib * npartials + partialind) -= 1;//total count for this partial
				}
			}
			
			//update local invariance, if necessary, only with one mutex lock
#ifdef _USE_LOCAL_INV
			//locale histogram has to be stored	
			//as this cannot happen too often, all the threads are updating the same array, and mutex lock is used to
			//to ensure, that only one thread updates the arrays at a time
			{
			  std::unique_lock<std::mutex> guard(thread_object.inv_mutex);//lock mutex
			
			  if (distance>=min_loc_r  && distance<=max_loc_r)
			  {
				  bin_index=int((distance-min_loc_r)*inv_locbinwidth);
				  *(local_count+loc_offset1+bin_index)+=1;
				  *(local_count+(*p_indexm2 *sum_bins_offset+loc_offset[*p_type1]+bin_index))+=1;
				  *(local_sum_count_finder[partialind]+bin_index)+=2;//sum count for this bin
			  }
			  if (distance_old>=min_loc_r  && distance_old<=max_loc_r)
			  {
				  bin_index_old=int((distance_old-min_loc_r)*inv_locbinwidth);
				  *(local_count+loc_offset1+bin_index_old)-=1;
				  *(local_count+(*p_indexm2 *sum_bins_offset+loc_offset[*p_type1]+bin_index_old))-=1;
				  *(local_sum_count_finder[partialind]+bin_index_old)-=2;//av count for this bin
			  }
			}
			
#endif		
			//index of the neighbour atom within its own type
			neighbour_index=*p_indexm2-conf.cumul[*p_type2];

			//Check coordination constraint, if there is any
			for (iconst=0;iconst<CoordNumbConst::nconstraints;iconst++)
			{
				if (dsquare>=cc_min[iconst] && dsquare<=cc_max[iconst])
				{
					//the distance is in the range, the number of neighbours
					//has to be increased for central atom 

					if (*p_type1==*p_type2)//both can be central and neighbour
					{
						*(coordnc.finder[iconst]+central_index)+=1;
						*(coordnc.finder[iconst]+neighbour_index)+=1;
					}
					else
					{
						if (*p_type1==CoordNumbConst::central[iconst])//moved particle is central for the constraint
							*(coordnc.finder[iconst]+central_index)+=1;
						else //second particle is central for the constraint
							*(coordnc.finder[iconst]+neighbour_index)+=1;
					}//end of if *p_type1==*p_type2
				}
				//Remove the contribution for the old distance
				if (dsquare_old>=cc_min[iconst]	&& dsquare_old<=cc_max[iconst])
				{ 
					//the distance is in the range, the number of neighbours
					//has to be decreased for moved atom 

					if (*p_type1==*p_type2)//both can be central and neighbour
					{
						*(coordnc.finder[iconst]+central_index)-=1;
						*(coordnc.finder[iconst]+neighbour_index)-=1;
					}
					else
					{
						if (*p_type1==CoordNumbConst::central[iconst])//moved particle is central for the constraint
							*(coordnc.finder[iconst]+central_index)-=1;
						else //second particle is central for the constraint
							*(coordnc.finder[iconst]+neighbour_index)-=1;
					}//end of if *p_type1==*p_type2
				}//distance is in the range
			}//end of checking and updating coordination number constraints

			//Check average coordination constraint, if there is any
			for (iavconst=0;iavconst<AvCoordConst::nconstraints;iavconst++)
			{
				if (AvCoordConst::avcctype[iavconst]==partialind)//there is a constraint for this atom pair
				{
					//Add the contribution for the new distance
					if (dsquare>=AvCoordConst::udminsq[iavconst] &&\
					   dsquare<=AvCoordConst::udmaxsq[iavconst])
					{
						//the distance is in the range, the counter for this constraint
						//has to be increased for this configuration 
#ifdef _NO_PERIODIC
						if (distance_i+AvCoordConst::dav_red[iavconst]<=R0_red)
#endif
							acoordc.neighbourcount[iavconst]++;
						if (*p_type1==*p_type2)//second atom can be a central two
#ifdef _NO_PERIODIC
							if (distance_j+AvCoordConst::dav_red[iavconst]<=R0_red)
#endif
								acoordc.neighbourcount[iavconst]++;
					}
					//Remove the contribution for the old distance
					if (dsquare_old>=AvCoordConst::udminsq[iavconst] &&\
					   dsquare_old<=AvCoordConst::udmaxsq[iavconst])
					{
						//the distance is in the range, the counter for this constraint
						//has to be decreased for this configuration
#ifdef _NO_PERIODIC
						if (distance_o_i+AvCoordConst::dav_red[iavconst]<=R0_red)
#endif
							acoordc.neighbourcount[iavconst]--;
						if (*p_type1==*p_type2)//second atom can be a central two
#ifdef _NO_PERIODIC
							if (distance_o_j+AvCoordConst::dav_red[iavconst]<=R0_red)
#endif
								acoordc.neighbourcount[iavconst]--;
					}
				}//end of if there is a constraint
			}//end of checking and updating coordination number constraints

#ifdef _ADVANCED_GEOM_CONST
			//Check bond valence sum constraint, if there is any
			for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
			{
				if (dsquare <= bvs_max[iconst])
				{
					//the distance is in the range, calculate the bond valence
					bv = exp((bvs_R0[iconst] - distance*boxedge) / bvs_b[iconst]);
					if (*p_type1 == *p_type2)//both can be central and neighbour
					{
						*(bvs->finder[iconst] + central_index) += bv;
						*(bvs->finder[iconst] + neighbour_index) += bv;
					}
					else
					{
						if (*p_type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
							*(bvs->finder[iconst] + central_index) += bv;
						else //second particle is central for the constraint
							*(bvs->finder[iconst] + neighbour_index) += bv;
					}//end of if (type1==type2)
					
				}
				//remove the old contribution
				if (dsquare_old <= bvs_max[iconst])
				{
					//the old distance is in the range, calculate the bond valence
					bv = exp((bvs_R0[iconst] - distance_old*boxedge) / bvs_b[iconst]);
					if (*p_type1 == *p_type2)//both can be central and neighbour
					{
						*(bvs->finder[iconst] + central_index) -= bv;
						*(bvs->finder[iconst] + neighbour_index) -= bv;
					}
					else
					{
						if (*p_type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
							*(bvs->finder[iconst] + central_index) -= bv;
						else //second particle is central for the constraint
							*(bvs->finder[iconst] + neighbour_index) -= bv;
					}//end of if (type1==type2)
				}//distance is in the range
			}//end of checking and updating bvs constraints
			
#endif
			p_indexm2++;//sets the pointer to the index of the next second moved atom
			p_type2++;//sets the pointer to the type of the next second moved atom
		}//end of second atom cycle jmoved
		p_indexm1++;//sets the pointer to the index of the next first moved atom
		p_type1++;//sets the pointer to the type of the next first moved atom
		p_newcoord1+=3;//sets the pointer to the new coordinates of the next first moved atom
		p_oldcoord1+=3;//sets the pointer to the old coordinates of the next first moved atom
		p_oldchargegr_c1+=3;//sets the pointer to the old charge group centre coordinates of the next first moved atom
	};//end of first moved atom cycle imoved, end of calculation between the moved atoms


#ifdef _NO_PERIODIC
	//Calculates the central bin indices for the moved atoms
	p_indexm1=move->indices;//sets the pointer to the index of the first moved atom
	p_newcoord1=move->newpos;//sets the pointer to the new coordinates of the first moved atom
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		dsquare=pow(*p_newcoord1,2)+pow(*(p_newcoord1+1),2)+pow(*(p_newcoord1+2),2);
		distance=sqrt(dsquare);
		//calculate the bin index, where the central atom is located in the spherical sample 
		bin_index=int((distance-xmin[0])*inv_binwidth[0]);
		atom_binind[*p_indexm1]=bin_index;//the bin index centred on the origin of the sample for this atom
		p_newcoord1+=3;
		p_indexm1++;
	}
		
#endif
	HistCalcChangeThread(thread_object.thread_arg[nthreads-1]);
#ifdef _TEST_MODE
	Threads::hctime2 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime2-Threads::hctime1; 
	Threads::dur_00+=elapsed.count();
#endif

#ifdef _AENET
#ifdef _TEST_MODE
	Threads::anntime1[nthreads-1]= std::chrono::high_resolution_clock::now();
#endif
	//for threading reason the main thread's ANN change will be called from here
	if (Aenet::calc_ANN)
		aenet.CalcANN(thread_object.thread_arg[nthreads-1]);
#ifdef _TEST_MODE
	Threads::anntime2[nthreads-1] = std::chrono::high_resolution_clock::now();
	elapsed=Threads::anntime2[nthreads-1]-Threads::anntime1[nthreads-1]; 
	Threads::dur_ann2[nthreads-1]+=elapsed.count();
#endif
#endif


	//waiting for all threads to finish histogram change calculation 
	if (nthreads>1)
	{
		std::unique_lock<std::mutex> guard(thread_object.CC_count_mutex);//lock mutex
		thread_object.CC_count++;
#ifdef _TH_ORDER
		cout << "p" << endl;
		logfile << "CC " << thread_object.CC_count << endl;
#endif
		if (thread_object.CC_count==nthreads)//time to proceed
		{
			thread_object.start_flag=0;//reset it to 0
			Threads::CC_count=0;//reset it to 0
			thread_object.CC_check_count.notify_all();
#ifdef _TH_ORDER
			cout << ":" << endl;
#endif
		}
		else thread_object.CC_check_count.wait(guard);
	};
	//------END OF MULTI-THREADING------

#ifdef _TEST_MODE
	Threads::hctime1 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime1-Threads::hctime2; 
	Threads::dur_01+=elapsed.count();
#endif

	//Now update the histogram by adding the changes of the auxiliary threads, and reset the auxiliary threads histogram segment 
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (ipartial = 0; ipartial < npartials; ipartial++)
		{
			if (Threads::mod_partial[ipartial])
			{
				for (ithread = 0; ithread < nthreads - 1; ithread++)
				{
					p_main = finder[npartials*ib+ipartial];
					p_aux = thread_histcount_finder[ithread] + offset[npartials*ib + ipartial];

					for (i = 0; i < nbins[ib]; i++)
					{
						*p_main += *p_aux;
						*p_aux++ = 0;//reset it to zero
						p_main++;
					}
					//total counts
					*(ptotal + npartials * ib + ipartial) += *(thread_object.thread_arg[ithread].hist_tot_start + npartials * ib + ipartial);//*(ptotal+th_offset);
					*(thread_object.thread_arg[ithread].hist_tot_start + npartials * ib + ipartial) = 0;
#ifdef _USE_LOCAL_INV
					//update the average histogram calculated on the bin size of the local histogram 
					p_main = local_sum_count_finder[ipartial];
					p_aux = thread_lochistcount_finder[ithread] + loc_offset[ipartial];

					for (i = 0; i < max_loc_nbins; i++)
					{
						*p_main += *p_aux;
						*p_aux++ = 0;//reset it to zero
						p_main++;
					}
#endif
				}//end of if this partial has been modified
			}//end of ithread cycle
		}//end of cycle ipartial
	}//end of cycle different bin size

#ifdef _USE_LOCAL_INV
	int *pi_main,*pi_aux;
	//Now update the local histogram by adding the changes of the auxiliary threads, and reset the auxiliary threads histogram segment
	p_type1=move->types;//type of the moved atoms
	p_indexm1=move->indices;//index of the moved atoms
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		for (itype=0;itype<ntypes;itype++)
		{
			//determinig the partial for this moved atom
			ipartial=(*p_type1<=itype ? (*p_type1*ntypes-(*p_type1*(*p_type1+1)/2)+itype) : \
				(itype*ntypes-(itype*(itype+1)/2)+*p_type1));
			if (Threads::mod_partial[ipartial])
			{
				for (ithread=0;ithread<nthreads-1;ithread++)
				{
					pi_main=local_count+ *p_indexm1*sum_bins_offset+loc_offset[itype];
					pi_aux=dlocal_count_finder[ithread]+imoved*sum_bins+loc_offset[itype];

					for (i=0;i<max_loc_nbins;i++)
					{
						*pi_main+=*pi_aux;
						*pi_aux++=0;//reset it to zero
						pi_main++;
					}
					
				}//end of if this partial has been modified
			}//end of ithread cycle
		}//end of cycle itype
		p_type1++;
		p_indexm1++;
	}//end of imoved cycle
#endif

#ifdef _TEST_MODE
	Threads::hctime2 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime2-Threads::hctime1; 
	Threads::dur_02+=elapsed.count();
#endif

	if (RunParams::potential > 0)
	{
		for (ipartial = 0; ipartial < FNC_POT::pot_dim; ipartial++)
		{
			for (ithread = 0; ithread < nthreads - 1; ithread++)
			{
				//potentials
				FNC_POT::vdW_pot[ipartial] += *(thread_object.thread_arg[ithread].vdW_pot_start + ipartial);
				*(thread_object.thread_arg[ithread].vdW_pot_start + ipartial) = 0;
				FNC_POT::vdW_pot_s[ipartial] += *(thread_object.thread_arg[ithread].vdW_pot_s_start + ipartial);
				*(thread_object.thread_arg[ithread].vdW_pot_s_start + ipartial) = 0;
				FNC_POT::vdW_pot_l[ipartial] += *(thread_object.thread_arg[ithread].vdW_pot_l_start + ipartial);
				*(thread_object.thread_arg[ithread].vdW_pot_l_start + ipartial) = 0;
				
			}
		}
	}
	if (RunParams::potential ==1)//only for LJ
	{
		for (ipartial = 0; ipartial < FNC_POT::pot_dim; ipartial++)
		{
			for (ithread = 0; ithread < nthreads - 1; ithread++)
			{
				FNC_POT::Coulomb_pot[ipartial]+=*(thread_object.thread_arg[ithread].Coulomb_pot_start+ipartial);
				*(thread_object.thread_arg[ithread].Coulomb_pot_start+ipartial)=0;
				FNC_POT::Coulomb_pot_s[ipartial] += *(thread_object.thread_arg[ithread].Coulomb_pot_s_start + ipartial);
				*(thread_object.thread_arg[ithread].Coulomb_pot_s_start + ipartial) = 0;
				FNC_POT::Coulomb_pot_l[ipartial] += *(thread_object.thread_arg[ithread].Coulomb_pot_l_start + ipartial);
				*(thread_object.thread_arg[ithread].Coulomb_pot_l_start + ipartial) = 0;
				FNC_POT::vdW14_pot[ipartial]+=*(thread_object.thread_arg[ithread].vdW14_pot_start+ipartial);
				*(thread_object.thread_arg[ithread].vdW14_pot_start+ipartial)=0;
				FNC_POT::Coulomb14_pot[ipartial]+=*(thread_object.thread_arg[ithread].Coulomb14_pot_start+ipartial);
				*(thread_object.thread_arg[ithread].Coulomb14_pot_start+ipartial)=0;
			}//end of ithread cycle
		}//end of cycle ipartial
	}
	//Updating the average coordination number constraint 
	int offset = (AvCoordConst::nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*AvCoordConst::neighbourcount)) / sizeof(*AvCoordConst::neighbourcount)) : 0);//number of dummy integer elements
	for (iconst=0;iconst<AvCoordConst::nconstraints;iconst++)
	{
		if (Threads::mod_partial[AvCoordConst::avcctype[iconst]])
		{
			for (ithread=1;ithread<nthreads;ithread++)
			{
				//number of neighbours
				th_offset=ithread*(AvCoordConst::nconstraints+offset)+iconst;
				*(acoordc.neighbourcount+iconst)+=*(acoordc.neighbourcount+th_offset);
				*(acoordc.neighbourcount+th_offset)=0;
			}
		}//end of if this constraint was affected by the move
	}//end of constraint cycle iconst

#ifdef _TEST_MODE
	Threads::hctime1 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime1-Threads::hctime2; 
	Threads::dur_03+=elapsed.count();
#endif
	//---------------START MULTI-THREADING------------------
	if (coordnc.nconstraints>0)//there is/are coordination constraints
		//Updating the coordnumbs array and calculating the number of atoms satisfying the coordination constraint
		cc_mod=coordnc.UpdateCoordnumb(thread_object.thread_arg[nthreads-1]);
#ifdef _ADVANCED_GEOM_CONST
	//updating the valence array
	if (bvs->nconstraints > 0)
		bvs->UpdateValence(thread_object.thread_arg[nthreads - 1]);
#endif
	if (coordnc.nconstraints > 0 || RunParams::nbvs > 0)//there was multithreading
	{
#ifdef _TEST_MODE
		Threads::hctime2 = std::chrono::high_resolution_clock::now();
		elapsed = Threads::hctime2 - Threads::hctime1;
		Threads::dur_04 += elapsed.count();
#endif
		//waiting for all the threads to finish
		if (nthreads > 1)
		{
			std::unique_lock<std::mutex> guard(thread_object.UC_count_mutex);//lock mutex
			thread_object.UC_count++;
#ifdef _TH_ORDER
			cout << "q" << endl;
			logfile << "UC " << thread_object.UC_count << endl;
#endif
			if (thread_object.UC_count == nthreads)//time to proceed
			{
				thread_object.start_flag = 0;//reset it to 0
				thread_object.UC_count = 0;//reset it to 0
				thread_object.UC_check_count.notify_all();
			}
			else thread_object.UC_check_count.wait(guard);
		};

		//------END OF MULTI-THREADING------

#ifdef _TEST_MODE
		Threads::hctime1 = std::chrono::high_resolution_clock::now();
		elapsed = Threads::hctime1 - Threads::hctime2;
		Threads::dur_05 += elapsed.count();
#endif
	}
	if (coordnc.nconstraints > 0)
	{
		padding_offset=(coordnc.tot_subconst>0 ? (int)((CACHE_PADDING-sizeof(*coordnc.nsatisfy))/sizeof(*coordnc.nsatisfy)) : 0);//number of dummy integer elements
		if (cc_mod && nthreads>0)//if multi-threading and at least one constraint was modified, update the nsatisfy array for all the constraint by copying the contributions of the tread to the first segment of the array 
		{
			for (ithread=1;ithread<nthreads;ithread++)
			{
				pcoordn=coordnc.nsatisfy;//sets the pointer to the beginning
				pcoord_aux=coordnc.nsatisfy+ithread*(coordnc.tot_subconst+padding_offset);//the beginning of this thread in the nsatisfy array
				for (j=0;j<coordnc.tot_subconst;j++)
				{
					*pcoordn+=*pcoord_aux;
					*pcoord_aux++=0;
					pcoordn++;
				}
			}
		}
	}

#ifdef _TEST_MODE
	Threads::hctime2 = std::chrono::high_resolution_clock::now();
	elapsed=Threads::hctime2-Threads::hctime1; 
	Threads::dur_06+=elapsed.count();
#endif
#ifdef _ADVANCED_GEOM_CONST
	delete[] bvs_max;
	delete[] bvs_b;
	delete[] bvs_R0;
#endif
	delete[] cc_min;
	delete[] cc_max;
};//the calculation of the change in the histogram performed by the threads	


void HistoSet::HistCalcChangeThread(ThreadArg &thread_arg)
{
	int i,ib,icoord;
	int partialind,c_partind=0,pot_partialind;//index of the partial;pow(2,index of the partial),partial ind based on the host atom's type fpr the virtuals
	int min_ind,max_ind;//cycle limits
	int bin_index=0,bin_index_old=0;//index of the histogram bin
	int cc_neigh_type;//type of the neighbour atom for a CoordNumbConst
	int cc_neigh_offset;//starting point for the given constraint in CoordNumbConst arrays having tot_neightype elements 

	//int moveout_flag;//indicator that the moved atom is a "too close" one still below cutoff
	int central_index,neighbour_index;//index of the central and neighbour atoms 
	int skip_cycle;//flag to indicate whether the calculation has to be performed for the given two atoms
	int imoved, jmoved,j,iconst,iavconst,type2,neightype;
	int *p_indexm1,*p_indexm2,*p_type1;
	int  mytabpot_index=0;//index among the given tabulated potentials for the presen partial
    //	longint start,finish;//to check the speed of program parts
	longint *my_total;//pointer for this thread hist_tot_start
	double dsquare,dsquare_old,distance,distance_old;
	double *cc_min, *cc_max;//for the actual central neighbour pair for each constarint
	double *p_newcoord1,*p_oldcoord1,*coordj;
	double dcomp[3];//the differences between the coordinates of two atoms
	int igr=0,jgr=0,GRtype1=0,GRtype2=0, GRpartialind=0;//type of the atom according to the separate GROMACS types arrays
	int *p_charge_c_indices=0;//pointer to the charge group centre index of the first moved atom
	double *my_vdWpot=0,*my_Coulombpot=0, *my_s_vdWpot = 0, *my_s_Coulombpot = 0, *my_l_vdWpot = 0, *my_l_Coulombpot = 0;//pointer for this thread potential
	double *pi=0,*pj=0;//pointer to the charge group centres
	double *p_oldchargegr_c1=0;//pointer to the old charge group centre coordinates of the moved atom
	double dcharge_c_sq=0,dcharge_c_sq_old=0;
	double tabpot_rminsq=0;//the square of the first tabulated r value for the given partial
	double tabcutoff_sq=0;//the squared cutoff for the given data set
	
	my_vdWpot=thread_arg.vdW_pot_start;//the starting point of the vdW potential for this thread
	my_Coulombpot=thread_arg.Coulomb_pot_start;//the starting point of the Coulomb potential for this thread
	my_s_vdWpot = thread_arg.vdW_pot_s_start;//the starting point of the vdW potential for this thread
	my_s_Coulombpot = thread_arg.Coulomb_pot_s_start;//the starting point of the Coulomb potential for this thread
	my_l_vdWpot = thread_arg.vdW_pot_l_start;//the starting point of the vdW potential for this thread
	my_l_Coulombpot = thread_arg.Coulomb_pot_l_start;//the starting point of the Coulomb potential for this thread
#ifdef _ADVANCED_GEOM_CONST
	int bvs_neigh_type;//type of the neighbour atom for bvs const
	int bvs_neigh_offset;//starting point for the given constraint in BondValenceSumConst arrays having tot_neightype elements 
	double *bvs_max;
	double *bvs_R0, *bvs_b;//R0 and b for the current neighbour of each constraint
	double bv;
#endif
#ifdef _USE_LOCAL_INV
	int loc_offset1;
	int *my_local_count,*parray_pos,my_sum_bins;
	longint *partial_loc_finder;//pointer to the local_sum_count 
#endif
#ifdef _NO_PERIODIC
	double distance_i=0,distance_j=0,distance_o_i=0;
#endif
	bool normal_calc;//if virtual is invoved, only potential is calculated
	//double duration;//to see time differences

	//Calculation between the moved atoms and other atoms 
	p_indexm1=move->indices;//sets the pointer to the index of the first moved atom
	p_type1=move->types;//sets the pointer to the type of the first moved atom
	p_newcoord1=move->newpos;//sets the pointer to the new coordinates of the first moved atom
	p_oldcoord1=move->oldpos;//sets the pointer to the old coordinates of the first moved atom
	p_charge_c_indices=move->charge_c_indices;//sets the pointer to the charge group centre index of the first moved atom
	p_oldchargegr_c1=move->old_charge_c_pos;//sets the pointer to the old coordinates of the first moved atom's charge group centre coordinates
	SetArraysize(&cc_min, RunParams::nicoord, "cc_min", "HistoSet::HistCalcChangeThread");
	SetArraysize(&cc_max, RunParams::nicoord, "cc_max", "HistoSet::HistCalcChangeThread");
#ifdef _ADVANCED_GEOM_CONST
	SetArraysize(&bvs_max, BondValenceSumConst::nconstraints, "bvs_max", "HistoSet::HistCalcChangeThread");
	SetArraysize(&bvs_R0, BondValenceSumConst::nconstraints, "bvs_R0", "HistoSet::HistCalcChangeThread");
	SetArraysize(&bvs_b, BondValenceSumConst::nconstraints, "bvs_b", "HistoSet::HistCalcChangeThread");
#endif
	
	my_total=thread_arg.hist_tot_start;//the starting points of the total counts for this thread
#ifdef _USE_LOCAL_INV
	int thread_id;//ID of this thread
	thread_id=thread_arg.thread_index;//the ID of the thread
	my_local_count=thread_arg.loc_hist_count;//this will be used only for the local count update of the moved atoms, not for the neighbours!
	if (thread_id==nthreads-1)
	{
		//for the main thread, it will work on the local_count arary, where sum_bins_offset is used
		parray_pos=p_indexm1;
		my_sum_bins=sum_bins_offset;
	}
	else
	{
		//for the aux threads, they will work on the dlocal_count arary, where sum_bins is used
		parray_pos=&imoved;
		my_sum_bins=sum_bins;
	}

#endif
	
	for (imoved=0;imoved<Move::tot_moved_atoms+Move::tot_moved_virtuals;imoved++)
	{
		//index of the central atom within its own type
		central_index=*p_indexm1-conf.cumul[*p_type1];
		if (RunParams::potential==1)
		{
			if (RunParams::fnc==4)
			{
				igr=*p_charge_c_indices;//charge centre index of the moved atom
				pi=conf.charge_gr_centre+3*igr;
			}

			//find the GROMACS type (0 to n_sepGRtypes-1) of the atom 
			GRtype1=0;
			while (*p_indexm1>=SimpleCfg::cumul_GR[(GRtype1) +1])
				(GRtype1)++;
			//GRtype1 now points to the type index of the seperate GROMACS type segment, where the atom is located (the atoms index is continuously increasing)
			//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
			GRtype1=SimpleCfg::GRsep_GR_type[GRtype1];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters
		}

		for (type2=0;type2<ntypes+nvirtualtypes;type2++)//neighbour atom's type 
		{
			partialind = (*p_type1 <= type2 ? (*p_type1 * ntypes - (*p_type1 * (*p_type1 + 1) / 2) + type2) : \
				(type2 * ntypes - (type2 * (type2 + 1) / 2) + *p_type1));
			//Determining the index of the partial 
			if (*p_type1 < ntypes && type2 < ntypes)//normal atoms
			{
				normal_calc = true;
				pot_partialind = partialind;
			}
			else
			{
				normal_calc = false;
				//for the virtual sites the host atom's type will decide the partialind
				if (*p_type1 >= ntypes)
					i = SimpleCfg::virtual_host_type[*p_type1 - ntypes];
				else
					i = *p_type1;
				if (type2 >= ntypes)
					j = SimpleCfg::virtual_host_type[type2 - ntypes];
				else
					j = type2;
				pot_partialind = (i <= j ? (i * ntypes - (i * (i + 1) / 2) + j) : (j * ntypes - (j * (j + 1) / 2) + i));
			}

			if (normal_calc)
			{
				c_partind = (int)pow((double)2, partialind);//for the coordination number constraint
				if (RunParams::potential == 10 && calc_tabpot[partialind])//is there is a tabulated potential for this partial
				{
					for (i = 0; i < RunParams::nused_potpartials; i++)
					{
						if (FNC_POT::tabpot_index[i] == partialind)
						{
							tabcutoff_sq = pow(FNC_POT::tab_cutoff[i], 2.0);//reduced, squared
							mytabpot_index = i;
							tabpot_rminsq = pow(FNC_POT::r_pot[i], 2);
							break;
						}
					}
				}
				//setting the actual minimal and maximal distance sqaures for each CoordNumbconst constraint
				for (iconst = 0; iconst < CoordNumbConst::nconstraints; iconst++)
				{
					if (CoordNumbConst::cctype[iconst] & c_partind)//there is a constraint for this type pair
					{
						if (*p_type1 == CoordNumbConst::central[iconst])//moved particle is central for the constraint
							cc_neigh_type = type2;//second is the neighbour
						else
							cc_neigh_type = *p_type1;//moved is the neighbour
						cc_neigh_offset = CoordNumbConst::cum_n_neightype[iconst];
						for (neightype = 0; neightype < CoordNumbConst::n_neightype[iconst]; neightype++)
						{
							if (cc_neigh_type == CoordNumbConst::neighbours[cc_neigh_offset + neightype])
							{
								//the type of the neighbour atom match the neightype-th neighbour type of this constraint
								cc_min[iconst] = CoordNumbConst::udminsq[cc_neigh_offset + neightype];
								cc_max[iconst] = CoordNumbConst::udmaxsq[cc_neigh_offset + neightype];
								break;
							}
						}
					}
					else
					{
						//the iconst constraint is not for this type pair, sets the minimum and maximum distance squares to unrealistic values
						cc_min[iconst] = 1e6;
						cc_max[iconst] = -1;//this is enough to make it impossible to find a distance square in range
					}
				}
			}
#ifdef _ADVANCED_GEOM_CONST
			//setting the actual maximal distance sqaures and other parameters for each BVS constraint
			for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
			{
				if (BondValenceSumConst::bvsctype[iconst] & c_partind)//there is a constraint for this type pair
				{
					if (*p_type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
						bvs_neigh_type = type2;//second is the neighbour
					else
						bvs_neigh_type = *p_type1;//first is the neighbour
					bvs_neigh_offset = BondValenceSumConst::cum_n_neightype[iconst];
					for (neightype = 0; neightype < BondValenceSumConst::n_neightype[iconst]; neightype++)
					{
						if (bvs_neigh_type == BondValenceSumConst::neighbours[bvs_neigh_offset + neightype])
						{
							//the type of the neighbour atom match the neightype-th neighbour type of this constraint
							bvs_max[iconst] = BondValenceSumConst::maxsq[bvs_neigh_offset + neightype];
							bvs_R0[iconst] = BondValenceSumConst::R0[bvs_neigh_offset + neightype];
							bvs_b[iconst] = BondValenceSumConst::b[bvs_neigh_offset + neightype];
							break;
						}
					}
				}
				else
				{
					//the iconst constraint is not for this type pair, sets the minimum and maximum distance squares to unrealistic values
					bvs_max[iconst] = -1;//this is enough to make it impossible to find a distance square in range

				}
			}
#endif
			min_ind=conf.cumul[type2]+thread_arg.min_atom_index[type2];//the index of the first atom in the conf
			max_ind=conf.cumul[type2]+thread_arg.max_atom_index[type2];//the index of the last atom in the conf
					
			//Sets the coordj pointer to the coordinates of the first neighbour atom for this thread of type2 	
			coordj=conf.positions+(3*min_ind);

			neighbour_index=thread_arg.min_atom_index[type2];//sets the index of the neighbour in its own type to the first atom for this thread
			if (normal_calc)
			{
#ifdef _USE_LOCAL_INV

				loc_offset1 = (*parray_pos) * my_sum_bins + loc_offset[type2];
				partial_loc_finder = thread_arg.hist_loc_count_finder[partialind];//the first loc_av_counts for this thread and partials

#endif
#ifdef _NO_PERIODIC
				distance_i = sqrt(pow(*p_newcoord1, 2) + pow(*(p_newcoord1 + 1), 2) + pow(*(p_newcoord1 + 2), 2));
				distance_o_i = sqrt(pow(*p_oldcoord1, 2) + pow(*(p_oldcoord1 + 1), 2) + pow(*(p_oldcoord1 + 2), 2));
#endif
			}
			for (j=min_ind;j<=max_ind;j++)//Neighbour atom: Cycling through the atoms of the second type
			{
				//first to check, whether this is not a moved-moved pair, 
				//among them the central and neighbour can be the same
				p_indexm2=move->indices;
				skip_cycle=0;//flag to indicate, whether the rest of the j cycle has to be done
				for (jmoved=0;jmoved<Move::tot_moved_atoms + Move::tot_moved_virtuals;jmoved++)
				{
					if (j==p_indexm2[jmoved])
					{
						//this is a moved-moved pair, the rest of the j cycle is skipped
						skip_cycle=1;
						coordj+=3;//skip the coordinates
						neighbour_index++;//skip the index of the neighbour in its own type
						break;
					}
				}
				if (skip_cycle)
					continue;//continue with the next moved atom
				//the moved atom is not a "too close" one still below cutoff by default
				//moveout_flag=0;

				//calculating the distance of the two atoms
				for (icoord=0;icoord<3;icoord++)
					dcomp[icoord]=p_newcoord1[icoord]-coordj[icoord];
				
				GetMinImage(dcomp);//taking into account the periodical boundary conditions 
				
				dsquare = 0;
				for (icoord = 0; icoord < 3; icoord++)
					dsquare += pow(dcomp[icoord], 2);

				//calculating the old distance of the two atoms
				for (icoord = 0; icoord < 3; icoord++)
					dcomp[icoord]=p_oldcoord1[icoord]-*coordj++;
				
				GetMinImage(dcomp);//taking into account the periodical boundary conditions 
				dsquare_old=0;
				for (icoord = 0; icoord < 3; icoord++)
					dsquare_old += pow(dcomp[icoord], 2);
#ifdef _NO_PERIODIC
			distance_j=sqrt( pow(*(coordj-3),2) + pow(*(coordj-2),2) + pow(*(coordj-1),2));
#endif				
				//updating the potential if it is any

				if (RunParams::potential==1)//LJ
				{
					if (RunParams::fnc==4)
					{
						//calculating the distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
						dcharge_c_sq=0;
						dcharge_c_sq_old=0;
						jgr=fnc[0].charge_centre[j];//charge group centre index of the neighbour
						if (igr!=jgr)//only calculate, if they are not in the same charge group centre
						{
							pj=conf.charge_gr_centre+3*jgr;
							for (icoord = 0; icoord < 3; icoord++)
								dcomp[icoord] = pi[icoord] - pj[icoord];
							
							GetMinImage(dcomp);//taking into account the periodical boundary conditions 
							dcharge_c_sq = 0;
							for (icoord = 0; icoord < 3; icoord++)
								dcharge_c_sq += pow(dcomp[icoord], 2);
							
							//calculating the old distance of the centre of the charge groups, as this has to be used deciding, whether the distance is inside cutoff
							//although the neighbour cannot be a moved atom, but in can be in the charge group of a moved atom, so it has to be checked
							//the neighbour"s charge group does not moved by default
							for (jmoved=0;jmoved<move->tot_moved_atoms;jmoved++)//the virtuals if tgere is any are in the same charge group, not needed here
							{
								if (fnc[0].charge_centre[move->indices[jmoved]]==jgr)//atom j is in the charge group of a moved atom
									pj=move->old_charge_c_pos+3*jmoved;
							}

							for (icoord = 0; icoord < 3; icoord++)
								dcomp[icoord]=p_oldchargegr_c1[icoord] - *pj++;
													
							GetMinImage(dcomp);//taking into account the periodical boundary conditions 
							dcharge_c_sq_old = 0;
							for (icoord = 0; icoord < 3; icoord++)
								dcharge_c_sq_old += pow(dcomp[icoord], 2);
													
						}
					}
					else
					{
						//no molecules, no charge groups, the atomic distances has to be used 
						dcharge_c_sq=dsquare;
						dcharge_c_sq_old=dsquare_old;
					}
					
					if (dcharge_c_sq <= rundat.Coulomb_cutoff_sq)
					{
						double temp = (fnc[0].CalcCoulomb)(*p_indexm1, j, dsquare);//calculating the potential for the new distance
						if (fabs(temp) > POT_SPLIT_HIGH)
							*(my_l_Coulombpot + pot_partialind) += temp;
						else if (fabs(temp) < POT_SPLIT_LOW)
							*(my_s_Coulombpot + pot_partialind) += temp;
						else
							*(my_Coulombpot + pot_partialind) += temp;
					}
					if (dcharge_c_sq_old <= rundat.Coulomb_cutoff_sq)
					{
						double temp = (fnc[0].CalcCoulomb)(*p_indexm1, j, dsquare_old);//calculating the potential for the old distance
						if (fabs(temp) > POT_SPLIT_HIGH)
							*(my_l_Coulombpot + pot_partialind) -= temp;
						else if (fabs(temp) < POT_SPLIT_LOW)
							*(my_s_Coulombpot + pot_partialind) -= temp;
						else
							*(my_Coulombpot + pot_partialind) -= temp;
					}
					if (dcharge_c_sq <= rundat.vdW_cutoff_sq || dcharge_c_sq_old <= rundat.vdW_cutoff_sq)
					{
						//find the GROMACS type (0 to nsepGRtypes-1) of the atom 
						GRtype2 = 0;
						while (j >= SimpleCfg::cumul_GR[(GRtype2)+1])
							(GRtype2)++;
						//GRtype2 now points to the type index of the seperate GROMACS type segmentm where the atom is located (the atoms index is continuously increasing)
						//there can be segments with the same GROMACS type at different places of the array, determine the GROMACS type
						GRtype2 = SimpleCfg::GRsep_GR_type[GRtype2];//this is the real GROMACS type of the atoms, corresponding to the LJ parameters
						//GRpartialind is GROMACS type based partial index, needed to have the vdW parameters, but the histogram is split into RMC-type partials
						GRpartialind = (GRtype1 <= GRtype2 ? (GRtype1 * Topology::nGRtypes - (GRtype1 * (GRtype1 + 1) / 2) + GRtype2) : \
							(GRtype2 * Topology::nGRtypes - (GRtype2 * (GRtype2 + 1) / 2) + GRtype1));

						
						if (dcharge_c_sq <= rundat.vdW_cutoff_sq)
						{
							double temp= (fnc[0].*CalcvdW)(GRpartialind, dsquare);//calculating the potential for the new distance
							if (fabs(temp) > POT_SPLIT_HIGH)
								*(my_l_vdWpot + pot_partialind) += temp;
							else if (fabs(temp)<POT_SPLIT_LOW)
								*(my_s_vdWpot + pot_partialind) += temp;
							else
								*(my_vdWpot + pot_partialind) += temp;
						}
						if (dcharge_c_sq_old <= rundat.vdW_cutoff_sq)
						{
							double temp = (fnc[0].*CalcvdW)(GRpartialind, dsquare_old);//calculating the potential for the old distance
							if (fabs(temp) > POT_SPLIT_HIGH)
								*(my_l_vdWpot + pot_partialind) -= temp;
							else if (fabs(temp) < POT_SPLIT_LOW)
								*(my_s_vdWpot + pot_partialind) -= temp;
							else
								*(my_vdWpot + pot_partialind) -= temp;
						}
					}

				}	
			
				if (!normal_calc)
				{
					neighbour_index++;//index of the neighbour in its own type
					continue;//j cycle
				}
				if (RunParams::potential == 10 && calc_tabpot[partialind])
				{
					//check whether the distance is in the range of the tabulated potential for this partial
					if (dsquare >= tabpot_rminsq && dsquare < tabcutoff_sq)
					{
						double temp = (fnc[0].*CalcvdW)(mytabpot_index, dsquare);
						if (fabs(temp)>POT_SPLIT_HIGH)
							*(my_l_vdWpot + mytabpot_index) += temp;
						else if (fabs(temp) < POT_SPLIT_LOW)
							*(my_s_vdWpot + mytabpot_index) += temp;
						else
							*(my_vdWpot + mytabpot_index) += temp;
					}

					if (dsquare_old >= tabpot_rminsq && dsquare_old < tabcutoff_sq)
					{
						double temp = (fnc[0].*CalcvdW)(mytabpot_index, dsquare_old);
						if (fabs(temp) > POT_SPLIT_HIGH)
							*(my_l_vdWpot + mytabpot_index) -= temp;
						else if (fabs(temp) < POT_SPLIT_LOW)
							*(my_s_vdWpot + mytabpot_index) -= temp;
						else
							*(my_vdWpot + mytabpot_index) -= temp;
					}
				}
		
				//Updating the histogram
				distance=sqrt(dsquare);
				//calculate the bin index, where this distance should be put
				for (ib = 0; ib < ndiffbin; ib++)
				{
					if (distance >= xmin[ib] && distance <= xmax[ib])
					{
						bin_index = int((distance - xmin [ib])*inv_binwidth[ib]);
						*(thread_arg.hist_count_finder[ib*npartials+partialind] + bin_index) += 1;//count for this bin
						*(my_total + ib*npartials+partialind) += 1;//total count for this partial
					
					}
				}
				//Removing the contribution of the old distances
				//old coordinates of the first moved atom
				
				distance_old=sqrt(dsquare_old);

				//calculate the bin index, where this distance should be put
				for (ib = 0; ib < ndiffbin; ib++)
				{
					if (distance_old >= xmin[ib] && distance_old <= xmax[ib])
					{
						bin_index_old = int((distance_old - xmin[ib])*inv_binwidth[ib]);
						*(thread_arg.hist_count_finder[ib*npartials + partialind] + bin_index_old) -= 1;//count for this bin
						*(my_total + ib * npartials + partialind) -= 1;//total count for this partial
						
					}
				}
				//update local invariance, if necessary
#ifdef _USE_LOCAL_INV
				//locale histogram has to be stored, and the average for the same bin size
				if (distance>=min_loc_r  && distance<=max_loc_r)
				{
					bin_index=int((distance-min_loc_r)*inv_locbinwidth);
					*(my_local_count+loc_offset1+bin_index)+=1;
					*(local_count+j *sum_bins_offset+loc_offset[*p_type1]+bin_index)+=1;
					*(partial_loc_finder+bin_index)+=2;//sum count for this bin
				}

				if (distance_old>=min_loc_r  && distance_old<=max_loc_r)
				{
					bin_index_old=int((distance_old-min_loc_r)*inv_locbinwidth);
					*(my_local_count+loc_offset1+bin_index_old)-=1;
					*(local_count+j *sum_bins_offset+loc_offset[*p_type1]+bin_index_old)-=1;
					*(partial_loc_finder+bin_index_old)-=2;//average count for this bin
				}
#endif	

				//Check coordination constraint, if there is any
				for (iconst=0;iconst<CoordNumbConst::nconstraints;iconst++)
				{
					if (dsquare>=cc_min[iconst] && dsquare<=cc_max[iconst])
					{						
						//the distance is in the range, the number of neighbours
						//has to be increased for the moved atom 

						if (*p_type1==type2)//both can be central and neighbour
						{
							*(thread_arg.coord_numb_finder[iconst]+central_index)+=1;
							*(thread_arg.coord_numb_finder[iconst]+neighbour_index)+=1;
						}
						else
						{
							if (*p_type1==CoordNumbConst::central[iconst])//moved particle is central for the constraint
								*(thread_arg.coord_numb_finder[iconst]+central_index)+=1;
							else //second particle is central for the constraint
								*(thread_arg.coord_numb_finder[iconst]+neighbour_index)+=1;
							
						}//end of if *p_type1==type2
					}//end of if dsquare is in the range

					//remove the contribution of the old distance
					if (dsquare_old>=cc_min[iconst] && dsquare_old<=cc_max[iconst])
					{
						//the distance is in the range, the number of neighbours
						//has to be decreased for central atom 

						if (*p_type1==type2)//both can be central and neighbour
						{
							*(thread_arg.coord_numb_finder[iconst]+central_index)-=1;
							*(thread_arg.coord_numb_finder[iconst]+neighbour_index)-=1;
						}
						else
						{
							if (*p_type1==CoordNumbConst::central[iconst])//moved particle is central for the constraint
								*(thread_arg.coord_numb_finder[iconst]+central_index)-=1;
							else //second particle is central for the constraint
								*(thread_arg.coord_numb_finder[iconst]+neighbour_index)-=1;
						}//end of if *p_type1==type2
					}//end of if dsquare_old is in the range
				}//end of checking and updating coordination number constraints

				//Check average coordination constraint, if there is any
				for (iavconst=0;iavconst<AvCoordConst::nconstraints;iavconst++)
				{
					if (AvCoordConst::avcctype[iavconst]==partialind)//there is a constraint for this atom pair
					{
						//Add the contribution of the new distance
						if (dsquare>=AvCoordConst::udminsq[iavconst] &&\
							dsquare<=AvCoordConst::udmaxsq[iavconst])
						{
							//the distance is in the range, the counter for this constraint
							//has to be increased for this configuration 
#ifdef _NO_PERIODIC
							if (distance_i+AvCoordConst::dav_red[iavconst]<=R0_red)//the neighbour's position did not change
#endif
								thread_arg.av_coord_start[iavconst]++;

							if (*p_type1==type2)//second atom can be a central two
#ifdef _NO_PERIODIC
								if (distance_j+AvCoordConst::dav_red[iavconst]<=R0_red)//the neighbour's position did not change
#endif
									thread_arg.av_coord_start[iavconst]++;
						}//end of if dsquare is in the range

						//Remove the contribution of the old distance
						if (dsquare_old>=AvCoordConst::udminsq[iavconst] &&\
						   dsquare_old<=AvCoordConst::udmaxsq[iavconst])
						{
							//the distance is in the range, the counter for this constraint
							//has to be decreased for this configuration 
#ifdef _NO_PERIODIC
							if (distance_o_i+AvCoordConst::dav_red[iavconst]<=R0_red)//the neighbour's position did not change
#endif
								thread_arg.av_coord_start[iavconst]--;
	
							if (*p_type1==type2)//second atom can be a central two
#ifdef _NO_PERIODIC
								if (distance_j+AvCoordConst::dav_red[iavconst]<=R0_red)//the neighbour's position did not change
#endif
									thread_arg.av_coord_start[iavconst]--;
						}//end of if dsquare_old is in the range
					}//end of if there is a constraint
				}//end of checking and updating coordination number constraints
#ifdef _ADVANCED_GEOM_CONST
				//Check bond valence sum constraint, if there is any
				for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
				{
					if (dsquare <= bvs_max[iconst])
					{
						//the distance is in the range, calculate the bond valence
						bv = exp((bvs_R0[iconst] - distance*boxedge) / bvs_b[iconst]);
						if (*p_type1 == type2)//both can be central and neighbour
						{
							*(thread_arg.valence_finder[iconst] + central_index) += bv;
							*(thread_arg.valence_finder[iconst] + neighbour_index) += bv;
						}
						else
						{
							if (*p_type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
								*(thread_arg.valence_finder[iconst] + central_index) += bv;
							else //second particle is central for the constraint
								*(thread_arg.valence_finder[iconst] + neighbour_index) += bv;
						}//end of if (type1==type2)
					}
					//remove the old contribution
					if (dsquare_old <= bvs_max[iconst])
					{
						//the old distance is in the range, calculate the bond valence
						bv = exp((bvs_R0[iconst] - distance_old*boxedge) / bvs_b[iconst]);
						if (*p_type1 == type2)//both can be central and neighbour
						{
							*(thread_arg.valence_finder[iconst] + central_index) -= bv;
							*(thread_arg.valence_finder[iconst] + neighbour_index) -= bv;
						}
						else
						{
							if (*p_type1 == BondValenceSumConst::central[iconst])//first particle is central for the constraint
								*(thread_arg.valence_finder[iconst] + central_index) -= bv;
							else //second particle is central for the constraint
								*(thread_arg.valence_finder[iconst] + neighbour_index) -= bv;
						}//end of if (type1==type2)
					}//distance is in the range
				}//end of checking and updating bvs constraints
#endif
				neighbour_index++;//index of the neighbour in its own type
			}//end of cycle neighbour atom, cycle j
			
		}//end of cycle type2
		p_indexm1++;//sets the pointer to the index of the next moved atom
		p_type1++;//sets the pointer to the type of the next moved atom
		p_newcoord1+=3;//sets the pointer to the new coordinates of the next moved atom
		p_oldcoord1+=3;//sets the pointer to the old coordinates of the next moved atom
#ifdef _USE_LOCAL_INV
		if (thread_id==nthreads-1)
			parray_pos++;//sets the pointer to the next moved atom's index, only needed for the main thread, otherwise the pointer's value is increased
#endif

		p_charge_c_indices++;//go to the charge group centre of the next moved atom
		p_oldchargegr_c1+=3;//sets the pointer to the old charge group centre coordinates of the next moved atom
	}//end of cycle imoved
	
	//Potential correction, if it is any

	if (move->update_charge_c)//if charge groups has to be updated (potenial with molecules)
		fnc[0].NBPotChargeGrCorr(thread_arg,*move);
#ifdef _ADVANCED_GEOM_CONST
	delete[] bvs_max;
	delete[] bvs_b;
	delete[] bvs_R0;
#endif
	delete[] cc_min;
	delete[] cc_max;
};
#ifdef _TEST_MODE

//Saving the too close pairs
void HistoSet::SaveTooClosePairs(char *file_name)
{
	int i;
	ofstream file;
	OpenFile(file,file_name,"HistoSet::SaveTooClosePairs",0);//open file, check, whether it was successfully opened
	file<<"Total number of too close pairs:  "<<ntcp/2<<endl;
	for (i=0;i<ntcp/2;i++)
		file<<tcpairs[i*2]<<"\t"<<tcpairs[i*2+1]<<endl;
	file.close();
	delete [] tcpairs;
}

#endif

#ifdef _VIBR_AMP
	//calculating the convolution function for 
	void HistoSet::CalcConvFunction()
	{
		int i,ipartial;
		double norm,sum,max=0,*pconv;


		norm=1/sqrt(2*PI);
		//find the largest sigma
		for (ipartial=0;ipartial<npartials;ipartial++)
		{
			if (T_corr_sigma[ipartial]>max)
				max=T_corr_sigma[ipartial];

		}
		//the largest T_corr_sigma will be used for the determination of the convolution interval, and then the effctive 
		//length of the convoluted histograms, as all the partials should have the same size in the convoluted histograms array
		//the mean of the Gauss distribution will be in the middle of the histogram bin to be convoluted
		//This will be the 0th convolution bin, and the others will follow up to 3*sigma, so only one side of the 
		//Gauss distribution will be calculated here
		nbins_conv=(int)((5*max-RunParams::rspacing[0])/RunParams::rspacing[0])+1;//to be on the safe side
		SetArraysize(&thermal_conv,nbins_conv*npartials,"thermal_conv","HistoSet::CalcConvFunction");
		SetArraysize(&thermal_conv_finder,npartials,"thermal_conv_finder","HistoSet::CalcConvFunction");
		SetArraysize(&temp_conv_hist,2*nbins_conv-1,"temp_conv_hist","HistoSet::CalcConvFunction");
		
		for (ipartial=0;ipartial<npartials;ipartial++)
		{
			sum=0;
			thermal_conv_finder[ipartial]=thermal_conv+ipartial*nbins_conv;//initializing the finder
			pconv=thermal_conv_finder[ipartial];
			for (i=0;i<nbins_conv;i++)
			{
				*pconv=exp(-pow(i*RunParams::rspacing[0]/T_corr_sigma[ipartial],2)/2)*norm/T_corr_sigma[ipartial];
				sum+=*pconv;
				if (i>0)
					sum+=*pconv;
				pconv++;
			}
			pconv=thermal_conv_finder[ipartial];
			for (i=0;i<nbins_conv;i++)
			{
				*pconv/=sum;
				pconv++;
			}
		}
		nbins_used=nbins[0]-3*nbins_conv/5;//at the end of the convoluted histogram the bin count will be 
										//smaller, than the average, as the non-existing bins after nbins
										//should contribute to the values, so the last few bins corresponding to 
										//3*sigma will be discarded

	}

	//calculating the convoluted histogram taking into account the atomic vibrations 
	void HistoSet::ThermalCorr()
	{
		int ibin,i,j,sign,ipartial,sum,diff;

		longint *phist_T;
		double *pconv;
#ifdef _NO_PERIODIC
		longint int_count;//the count rounded to integer
		double *phist_d;//the periodic histogram is double precision
#endif
		longint *phist;

	
		for (ipartial=0;ipartial<npartials;ipartial++)
		{
			if (Threads::mod_partial[ipartial])
			{
				phist_T=finder_T[ipartial];
				for (i=0;i<(nbins[0]+2);i++)
					*phist_T++=0; //set it to zero
#ifdef _NO_PERIODIC
				phist_d=periodic_finder[ipartial];//the periodic histogram has to be used
#else
				phist=finder[ipartial];
#endif
				phist_T=finder_T[ipartial];
				pconv=thermal_conv_finder[ipartial];

				for (ibin=0;ibin<nbins[0];ibin++)
				{
#ifdef _NO_PERIODIC
					int_count=RoundNearInt(*phist_d);
					phist=&int_count;
					
#endif

					if (*phist==0)
					{
						phist++;
#ifdef _NO_PERIODIC	
						phist_d++;
#endif

						continue;
					}

					sign=-1;//sign for the negative part of the function
					sum=0;
					for (i=0;i<nbins_conv*2-1;i++)//going through the whole convolution function
					{
						j=sign*(i-nbins_conv+1);//index in the thermal_conv array
						temp_conv_hist[i]=RoundNearInt(*(pconv+j)* *phist);
						sum+=temp_conv_hist[i];
						if (j==0)
							sign=1;//set the sign for the positive part of the function
					}
					//now we have to make sure, that sum is equal to the original hist count, and correct it if necessary
					if (sum!=*phist)
					{
						diff=(int)*phist-sum;
						if (diff>0)
							sign=1;
						else
							sign=-1;
						if (diff&1)
						{
							//diff is odd, one will be added or subtracted from the middle bin
							temp_conv_hist[nbins_conv-1]+=sign;
							diff-=sign;//decrease diff
						}
						if (diff>2*(nbins_conv-1))//this cannot happen, just to be on the safe side
						{
							cout<<"\n*****ERROR*****"<<endl;
							cout<<"The remainder is larger than the number of convolution bins,"<<endl;
							cout<<"something is wrong with the vibrational amplitude convolution! Exiting..."<<endl;
							CleanExit();
						}
						//divide the remaining (or remove the excess) symetrically beginning with the bins 
						//closest to the middle, if excess and farthest from the middle if remove
						for (i=1;i<nbins_conv;i++)
						{
							j=i;
							if (diff<0)
							{
								j=nbins_conv-i;//to start from the edge and go toward the middle
								if (temp_conv_hist[nbins_conv-1+j]==0)//find the first with value
									continue;
							}
							if (diff==0)
								break;
							temp_conv_hist[nbins_conv-1+j]+=sign;
							temp_conv_hist[nbins_conv-1-j]+=sign;
							diff-=2*sign;
						}
						//Now make sure, that the middle bin is not smaller, than the two neighbours,
						//can happen, if they were equal, and the odd numner of excess was removed from it
						if (temp_conv_hist[nbins_conv-1]<temp_conv_hist[nbins_conv])//the difference can only be one
						{
							temp_conv_hist[nbins_conv]--;//remove one to make it equal to middle
							temp_conv_hist[nbins_conv-2]--;//remove one to make it equal to middle
							temp_conv_hist[nbins_conv+1]++;//add to the next
							temp_conv_hist[nbins_conv-3]++;//add to the next
						}
						//maybe now this second from the middle became larger than the next
						i=1;
						while (temp_conv_hist[nbins_conv-1+i]<temp_conv_hist[nbins_conv+i]  && i<nbins_conv-2)//make sure, that not going out of the array
						{//roll the excess toward the edge if necessary
							temp_conv_hist[nbins_conv+i]--;
							temp_conv_hist[nbins_conv-2-i]--;
							i++;
							temp_conv_hist[nbins_conv+i]++;
							temp_conv_hist[nbins_conv-2-i]++;
							
						}
					}//end of make corrections if necessary
					//now update the histogram
					for (i=-nbins_conv+1;i<nbins_conv;i++)
					{
						j=ibin+i;
						//if the convolution window would be outside the histogram, then it will only be stored for checking at the end of the histogram for this partial
						if (j<0) // this will not happen possibly because of the cutoff
							j=nbins[0];//put it in the first auxiliary bin after the histogram
						if (j>nbins[0]-1)
							j=nbins[0]+1;//put it in the second auxiliary bin after the histogram
						*(phist_T+j)+=temp_conv_hist[i+nbins_conv-1];
					}
#ifdef _NO_PERIODIC
					phist_d++;
#else
					phist++;
#endif
				}//end of ibin cycle
			
				//calculate the total
				ptotal_T[ipartial]=0;
				phist_T=finder_T[ipartial];
				for (ibin=0;ibin<nbins[0];ibin++)
					ptotal_T[ipartial]+=*phist_T++;
#ifndef _NO_PERIODIC
				if ((ptotal_T[ipartial]+*phist_T+*(phist_T+1))!=ptotal[ipartial])
				{
					cout<<"\n*****ERROR*****"<<endl;
					cout<<"The total count for the "<<ipartial+1<<". partial is "<<ptotal[ipartial]<<" for the"<<endl;
					cout<<"original and "<<ptotal_T[ipartial]+*phist_T+*(phist_T+1)<<" for the convoluted histogram!"<<endl;
					cout<<"There is an error in the convolution, cannot proceed, exiting..."<<endl;
					CleanExit();
				}
#endif
			}
		}//end of ipartial cycle
	};
#endif

#ifdef _NO_PERIODIC
//loading the central bin index of the atoms
int HistoSet::LoadPosBin()
{ 
	ifstream file;
	int i,idummy;
	SafeOpenTextFile(file,posbinfilename);
	if (CheckFileState(file,"HistoSet::LoadPosBin",posbinfilename)==0)
	{
		cout<<"Histogram and bin indices has to be calculated!"<<endl;
		return(0);//load was not successful
	};
	for (i=0;i<3;i++)
		SkipLine(file);

	for (i=0;i<SimpleCfg::ntotal;i++)
	{
		file>>idummy;
		file>>atom_binind[i];
	}
	file.close();
	return 1;
}
//saving the central bin index of the atoms
void HistoSet::SavePosBin()
{
	int i;
	ofstream file;
	OpenFile(file,posbinfilename,"HistoSet::Save",0);//open file, check, whether it was successfully opened

	file<<"The central indices of the bins, where the atoms can be found\n"<<endl;
	file<<"atom index\t"<<"bin index"<<endl;
	for (i=0;i<SimpleCfg::ntotal;i++)
	{
		file<<i+1<<"\t"<<atom_binind[i]<<endl;
	}
	file.close();

}

//calculating the histogram, which would correspond to periodic boundary conditions
void HistoSet::CalcPerHist(ThreadArg &thread_arg, DataMat &datmat)
{
	int itype, jtype,ipartial,ibin,i;
	int min_ind,max_ind,loc_offset1,*ppos;
	double *phist,*pdV;

	min_ind=thread_arg.min_bin[0];
	max_ind=thread_arg.max_bin[0];

	for (itype=0;itype<ntypes;itype++)
	{
		for (jtype=0;jtype<ntypes;jtype++)
		{
			ipartial=(itype<=jtype ? (itype*ntypes-(itype*(itype+1)/2)+jtype) : (jtype*ntypes-(jtype*(jtype+1)/2)+itype));

			if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
					//in case the main thread already generated the new move
			{
				if (itype<=jtype)//only zero it for the firts occurence of a partial (only for 1-2 but not for 2-1)
				{
					//Calculate the new modified PPCF partials
					phist=periodic_finder[ipartial]+min_ind;
					for (ibin=min_ind;ibin<=max_ind;ibin++)
						*phist++=0.0;
				}

				loc_offset1=SimpleCfg::cumul[itype]*HistoSet::sum_bins_offset+HistoSet::loc_offset[jtype];//to the beginning of the jtype local hist. of the first
					//atom of itype
				ppos=HistoSet::atom_binind+SimpleCfg::cumul[itype];//pointer to the central bin index where each atom can be found 
				for (i=0;i<SimpleCfg::pnatoms[itype];i++) 
				{
					phist=periodic_finder[ipartial]+min_ind;
					//summing the local density
					for (ibin=min_ind;ibin<=max_ind;ibin++)
					{
						*phist+= *(local_count+loc_offset1+ibin) / datmat.dV_in[*ppos *nbins[0]+ibin];
						phist++;
						
					}
					ppos++;
					loc_offset1+=HistoSet::sum_bins_offset;
				}
				
			}//end of if this partial was modified
		}//end of jtype cycle
	}//end of itype cycle

	//multiplying by dV (the volume element in case of periodic bondary  conditions to get 
	//the periodic histogram, and dividing by two, as for ii partials normal hist. only calculated for once for a pair,
	// and in case of ij partials, here this will be the sum os ij and ji
	for (ipartial=0;ipartial<npartials;ipartial++)
	{
		if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
					//in case the main thread already generated the new move
		{

			phist=periodic_finder[ipartial]+min_ind;
			pdV=datmat.dV+min_ind;
			for (ibin=min_ind;ibin<=max_ind;ibin++)
			{
				*phist*=*pdV/2.0;
				pdV++;
				phist++;
			}
		}
	}

};
#endif


