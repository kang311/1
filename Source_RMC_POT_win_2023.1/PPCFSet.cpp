//source PPCFSet.cpp
//Last changed 12.12.2022

#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//classes2.h included through this

//Defining static members
int   PPCFSet::ntypes;//number of atom types
int   PPCFSet::npartials;//number of partials
int   PPCFSet::ndiffbin;//number of different bin sizes
int  *PPCFSet::nbins;//number of points (equals number of hist. bins) for each bin size
int  *PPCFSet::mybins;//the number of used bins
#ifdef _VIBR_AMP
	int	 PPCFSet::nbins_used;//the number of bins where the ppcf will be good,  nbins-3*nbins_conv/5
#endif
#ifdef _NO_PERIODIC
	int PPCFSet::calcmode=1;//whether to calculate from the local hist
#endif
int *PPCFSet::offset;//array elements for each partial :total number of bins  
					//to the beginning of the given partial of the given bin size

double *PPCFSet::deltar;//bins width (in Angstrom) for each bin size dim [ndiffbin]
double *PPCFSet::rmin;//minimum value for the abcissas (in Angstrom) for each bin size dim [ndiffbin]
double *PPCFSet::rmax;//maximum value for the abcissas (in Angstrom) for each bin size dim [ndiffbin]
#ifdef _ATLAS
	double *PPCFSet::gr_0;//ppcf-1, needed for the ATLAS in case of S(Q) and F(Q) anf F(g) sets
#endif 

//----------------constructor with handle initialization----
PPCFSet::PPCFSet(HistoSet &histoset, DataMat &datamat,Move &move_obj)
				:hist(histoset),datmat(datamat),move(move_obj)//constructor
{ 
	int i,ib,sum,imax,index;

	if (ntypes==0)
	{
		cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"PPCFSet constructor"<<endl;
		cout<<"Number of types is 0!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};
	if (mybins[0]<=0)
	{
		cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"PPCFSet constructor"<<endl;
		cout<<"Number of bins was not calculated!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	SetArraysize(&gvalues,npartials*HistoSet::ntot_bins,"gvalues","PPCFSet::PPCFSet");//different ppcf for each different bin size
	SetArraysize(&finder,npartials*ndiffbin,"finder","PPCFSet::PPCFSet");

	//initialising the finder
	sum = 0;
	index = 0;
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (i = 0; i < npartials; i++)//for each partial of each different bin size
		{
			*(finder + index) = gvalues + sum;
			index++;
			sum += mybins[ib];
		}

	}

	//arrays connected to the collecting and averaging of the PPCF-s
	gvalues_sum=0;//default
	finder_sum=0;//default
        
	//ppcf-s has to be averaged
	imax=npartials*HistoSet::ntot_bins;//size of storage array

	if(RunParams::sum_ppcf)
	{
	    ncollect=0;//resetting to zero
	    SetArraysize(&gvalues_sum,imax,"gvalues_sum","PPCFSet::PPCFSet");
	    SetArraysize(&finder_sum,npartials*ndiffbin,"finder_sum","PPCFSet::PPCFSet");
	
	    for (i=0;i<imax;i++)
			gvalues_sum[i]=0.0;//initializing to zero
	    //initialising the finder
		sum = 0;
		index = 0;
		for (ib = 0; ib < ndiffbin; ib++)
		{
			for (i = 0; i < npartials; i++)//for each partial of each different bin size
			{
				*(finder_sum + index) = gvalues_sum + sum;
				sum += mybins[ib];
				index++;
			}

		}
	   
	}

};	

//---- copy constructor with handle initialization------------
PPCFSet::PPCFSet(PPCFSet &source, HistoSet &histoset, DataMat &datamat, Move &move_obj)
	:hist(histoset), datmat(datamat), move(move_obj)
{
	int i, ib, imax, sum, index;

	if (ntypes == 0)
	{
		cout << "\n" << "*****ERROR*****" << endl;
		cout << "PPCFSet copy constructor" << endl;
		cout << "Number of types is 0!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	};
	for (ib = 0; ib < ndiffbin; ib++)
	{
		if (mybins[ib] <= 0)
		{
			cout << "\n" << "*****ERROR*****" << endl;
			cout << "PPCFSet copy constructor" << endl;
			cout << "Number of bins was not calculated!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		};
	}

	SetArraysize(&gvalues, npartials*HistoSet::ntot_bins, "gvalues", "PPCFSet::PPCFSet");//different ppcf for each different bin size
	SetArraysize(&finder, npartials*ndiffbin, "finder", "PPCFSet::PPCFSet");

	
	//copying the PPCF from source to target
	imax=npartials* HistoSet::ntot_bins;//size of storage array
	for (i=0;i<imax;i++)
		gvalues[i]=source.gvalues[i];//PPCF copy

	//initialising the finder
	sum = 0;
	index = 0;
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (i = 0; i < npartials; i++)//for each partial of each different bin size
		{
			*(finder + index) = gvalues + sum;
			index++;
			sum += mybins[ib];

		}
	}
	
	//arrays connected to the collecting and averaging of the PPCF-s
	gvalues_sum=0;//default
	finder_sum=0;//default
	if (RunParams::sum_ppcf)
	{
		ncollect = 0;//resetting to zero
		SetArraysize(&gvalues_sum, imax, "gvalues_sum", "PPCFSet::PPCFSet");
		SetArraysize(&finder_sum, npartials*ndiffbin, "finder_sum", "PPCFSet::PPCFSet");
		ncollect = source.ncollect;
		for (i = 0; i < imax; i++)
			gvalues_sum[i] = source.gvalues_sum[i];//PPCF copy
		//initialising the finder
		sum = 0;
		index = 0;
		for (ib = 0; ib < ndiffbin; ib++)
		{
			for (i = 0; i < npartials; i++)//for each partial of each different bin size
			{
				*(finder_sum + index) = gvalues_sum + sum;
				sum += mybins[ib];
				index++;
			}

		}

	}
	
};	
//-------------Sets the static members------------------------
void PPCFSet::SetPPCFParams(RunParams &rundata)
{
	int i,ib,sum;

	//Checking the parameters
	if (RunParams::ntypes==0)
	{
		cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"PPCFSet::SetPPCFParams: Number of types is 0!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};
	if (rundata.nbins[0]<=0)
	{
		cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"PPCFSet::SetPPCFParams: Number of bins was not calculated!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};
	ntypes=RunParams::ntypes;//number of atom types
	npartials=ntypes*(ntypes+1)/2;
	ndiffbin = RunParams::ndiffbin;
	nbins=rundata.nbins;//number of bins (same for each histogram)
	SetArraysize(&mybins, ndiffbin, "mybins", "PPCFSet::SetPPCFParams");
#ifdef _VIBR_AMP
	nbins_used=HistoSet::nbins_used;
	mybins[0]=nbins_used;
#else
	mybins=nbins;
#endif
	SetArraysize(&rmin, ndiffbin, "rmin", "PPCFSet::SetPPCFParams");
	SetArraysize(&rmax, ndiffbin, "rmax", "PPCFSet::SetPPCFParams");
	
	for (i = 0; i < ndiffbin; i++)
	{
		rmin[i] = RunParams::binshift_diff[i]*RunParams::rspacing_diff[i];//minimum value for the abcissas (in Angstrom)
		rmax[i] = rundata.xmax[i]*RunParams::boxedge;//maximum value for the abcissas (in Angstrom)
		
	}
	deltar = RunParams::rspacing_diff;//bins width (in Angstrom)
	

	SetArraysize(&offset,npartials*ndiffbin,"offset","PPCFSet::SetPPCFParams");
	sum = 0;
	for (ib = 0; ib < ndiffbin; ib++)
	{
		for (i = 0; i < npartials; i++)
		{
			offset[i] = sum;
			sum += mybins[ib];
		}
	}
	
#ifdef _ATLAS
	if ((RunParams::nsq+RunParams::nfq+RunParams::nfg)>0)
		i=npartials*HistoSet::ntot_bins;
	else
		i=1;

	SetArraysize(&gr_0,i,"gr_0","PPCFSet::SetPPCFParams");//only needed if there is S(Q), F(Q) or F(g) data
	
#endif
#ifdef _NO_PERIODIC
	if (RunParams::nek>0)
	 calcmode=0;//do not calculate from the local hist, but from the previously calculated periodic hist
#ifdef _VIBR_AMP
	calcmode=0;
#endif
#endif
};


void PPCFSet::Save(const char *file_name)
{
	int i,ib,ipartial,count;
	ofstream file;

	//save to *.ppcf file
	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE+10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE+10, ppcffilename);

	OpenFile(file,tempfilename,"PPCFSet::Save",0);//open file, check, whether it was successfully opened

	file.setf(ios::fixed, ios::floatfield);
	file.setf(ios::right, ios::adjustfield);
	file.precision(9);
	file<<"This is an object of Class PPCFSet "<<endl;
	file<<ntypes<<"\t number of atom types"<<endl;
	
	
	
	for (ib = 0; ib < ndiffbin; ib++)
	{
		file << "\nFor data set(s) ";
		count = 0;
		for (i = 0; i < RunParams::ntot_datasets; i++)
		{
			if (RunParams::assign_hist[i] == ib)//this binsize is used for this data set
			{
				if (count > 0)
					file << ",";
				file << i + 1;
				count++;
			}
		}
		file << endl;
		file << mybins[ib] << "\t number of discretization bins for the ppcf " << endl;
		file << rmin[ib] << "\t rmin value (Angstroms)" << endl;
		file << rmax[ib] << "\t rmax value (Angstroms)" << endl;
		file << deltar[ib] << "\t discretization step (Angstroms) " << endl;
		file << "\nr value + " << npartials << " partial g(r) at the middle of the histogram bin" << endl;
		for (i = 0; i < mybins[ib]; i++)//for each histogram bin
		{
			file << J10 << rmin[ib] + (i + 0.5)*deltar[ib];
			for (ipartial = 0; ipartial < npartials; ipartial++)
			{
				file << "\t" << J20 << *(finder[ib*npartials+ipartial] + i);
			}
			file << endl;
		}
	}//end of ib cycle 

	if(RunParams::sum_ppcf)//ppcf-s has to be averaged
	{
	    if (ncollect>0)
	    {
			file<<"\npartial g(r)-s has been collected "<<ncollect<<" times, and averaged."<<endl; 
			for (ib = 0; ib < ndiffbin; ib++)
			{
				file << "\nFor data set(s) ";
				count = 0;
				for (i = 0; i < RunParams::ntot_datasets; i++)
				{
					if (RunParams::assign_hist[i] == ib)//this binsize is used for this data set
					{
						if (count > 0)
							file << ",";
						file << i + 1;
						count++;
					}
				}
				file << endl;

				file << "\r value + " << npartials << " averaged partial g(r) at the middle of the histogram bin" << endl;
				for (i = 0; i < mybins[ib]; i++)//for each histogram bin
				{
					file << J10 << rmin[ib] + (i + 0.5)*deltar[ib];
					for (ipartial = 0; ipartial < npartials; ipartial++)
						file << "\t" << J20 << *(finder_sum[ib*npartials+ipartial] + i) / ncollect;
					file << endl;
				}
			}//end of ib cycle
	    }
	    file.precision(6);
		file.unsetf(ios::fixed);
		file.unsetf(ios::right);
	}
	file.close();
	
};


//save the PPCF-s in RMCA format to *.out file
void PPCFSet::SaveOldOut(ofstream &file, const char *file_name) const
{
	//running will not continue, till the saving could be performed
	CheckFileState(file,"PPCFSet::SaveOldOut",file_name);// check, whether it was successfully opened
		
	int i,ib,ipartial;

	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);//formatting
	
	
	file<<"Results from RMC++ ! PARTIALS ARE IN ORDER: 1-1,1-2,1-3...2-2,2-3..,3-3,etc! "<<endl;
	file<<"TITLE "<<endl;
	file<<RunParams::title;
	file<<"\n"<<endl;
	file<<J10<<RunParams::ngr+RunParams::nsq+RunParams::nfq+RunParams::nek<<" experiment; ";
	file<<J10<<npartials<<" partials"<<endl;
	
	//saving the g(r) partials as calculated in the PPCFSet object
	file<<"Partial g(r)'s (from PPCF set)"<<endl;
	file<<J10<<"r"<<" "<<J10<<"g(r)"<<endl;
	file<<"PLOTS"<<endl;
	for (ib = 0; ib < HistoSet::ndiffbin; ib++)
	{
		file << J10 << mybins[ib]<< " " << J10 << npartials << " for bin size " << ib + 1 << endl;
		for (i = 0; i < mybins[ib]; i++)//for each histogram bin
		{
			file << J13 << rmin[ib] + (i + 0.5)*deltar[ib];//r abcissa for i-th point
														//i.e. central value of the ith interval
			for (ipartial = 0; ipartial < npartials; ipartial++)
				file << " " << J13 << *(finder[npartials*ib+ipartial] + i);//i-th point of j-th partial
			file << endl;
		}
	}
	file.unsetf(ios::fixed);//formatting
	file.unsetf(ios::right);
};

//-------calculating the PPCF from the histogram using the normalisation tables------------------------
void PPCFSet::CalcPPCF(ThreadArg &thread_arg)
{

	int ibin,ib;
	
#ifdef _NO_PERIODIC
	int i,itype,jtype,loc_offset1,*ppos;
#ifndef _VIBR_AMP
	double *pd_hist;
	pd_hist=0;
#else
	longint *p_hist;
	p_hist = 0;//to avoid warning
#endif
#else
	longint *p_hist;
	p_hist = 0;//to avoid warning
#endif
	double *p_matrix;

	int ipartial;
	int min_ind,max_ind;
	double *p_ppcf;
#ifdef _VIBR_AMP
	min_ind=thread_arg.min_bin_used;
	max_ind=thread_arg.max_bin_used;
#else
	min_ind=thread_arg.min_bin[0];//will be reset where necessary
	max_ind=thread_arg.max_bin[0];
#endif
	
#ifdef _NO_PERIODIC
	if (calcmode)
	{

		//only calculate it from the local hist count, if no EXAFS and no vibr amp to save the update the otherwise used periodic_hist
		//set it to zero
		for (ipartial=0;ipartial<npartials;ipartial++)
		{
			p_ppcf=finder[ipartial]+min_ind;
			for (ibin=min_ind;ibin<=max_ind;ibin++)
				*p_ppcf++=0;
		}

		for (itype=0;itype<ntypes;itype++)
		{
			for (jtype=0;jtype<ntypes;jtype++)
			{
				ipartial=(itype<=jtype ? (itype*ntypes-(itype*(itype+1)/2)+jtype) : (jtype*ntypes-(jtype*(jtype+1)/2)+itype));
				loc_offset1=(SimpleCfg::cumul[itype])*HistoSet::sum_bins_offset+HistoSet::loc_offset[jtype];//to the beginning of the jtype local hist. of the first
					//atom of itype
				ppos=HistoSet::atom_binind+SimpleCfg::cumul[itype];//pointer to the central bin index where each atom can be found 

				for (i=0;i<SimpleCfg::pnatoms[itype];i++)
				{
					p_ppcf=finder[ipartial]+min_ind;
					//summing the local density
					for (ibin=min_ind;ibin<=max_ind;ibin++)
					{
						*p_ppcf+= *(hist.local_count+loc_offset1+ibin) / datmat.dV_in[*ppos *nbins[0]+ibin];
						p_ppcf++;
						
					}

					ppos++;
					loc_offset1+=HistoSet::sum_bins_offset;

				}
			}//end of jtype cycle
		}//end of itype cycle
		//averaging the locale density and dividing by rho
		for (ipartial=0;ipartial<npartials;ipartial++)
		{
			p_ppcf=finder[ipartial]+min_ind;
			for (ibin=min_ind;ibin<=max_ind;ibin++)
			{
				*p_ppcf*=datmat.ntable[ipartial];
				p_ppcf++;
			}
		}
	}
	else
	{

#endif
		for (ib = 0; ib < ndiffbin; ib++)
		{
#ifdef _VIBR_AMP
			min_ind = thread_arg.min_bin_used;
			max_ind = thread_arg.max_bin_used;
#else
			min_ind = thread_arg.min_bin[ib];
			max_ind = thread_arg.max_bin[ib];
#endif
			for (ipartial = 0; ipartial < npartials; ipartial++)
			{
				p_ppcf = finder[npartials*ib+ipartial] + min_ind;
#ifdef _VIBR_AMP
				p_hist = hist.finder_T[ipartial] + min_ind;//use the convoluted histogram for ppcf calc
#else
#ifdef _NO_PERIODIC
				//only arrives here, if nek>0
				pd_hist = hist.periodic_finder[ipartial] + min_ind;
#else
				p_hist = hist.finder[npartials*ib + ipartial] + min_ind;//use the normal histogram for ppcf calc
#endif
#endif
				p_matrix = datmat.ntable_finder[npartials*ib + ipartial] + min_ind;
				for (ibin = min_ind; ibin <= max_ind; ibin++)
				{
#ifndef _NO_PERIODIC
					*p_ppcf++ = *p_hist++ * *p_matrix++;//normalising the histogram
#else
#ifdef _VIBR_AMP
					*p_ppcf++ = *p_hist++ * *p_matrix++;//normalising the histogram
#else
					*p_ppcf++ = *pd_hist++ * *p_matrix++;//normalising the histogram
#endif
#endif

				}
			}
		}
#ifdef _NO_PERIODIC
	}
#endif
};

//-----------------calculating the PPCF for the modified types from the histogram---------------
void PPCFSet::CalcModPPCF(ThreadArg &thread_arg)
{
	int ipartial,ib,ibin;
	
#ifdef _NO_PERIODIC
	int i,itype,jtype,loc_offset1,*ppos;
#ifndef _VIBR_AMP
	double *pd_hist;
	pd_hist=0;
#else
	longint *p_hist;
	p_hist = 0;//to avoid warning
#endif
#else
	longint *p_hist;
	p_hist = 0;
#endif
	double *p_ntable;
	int min_ind,max_ind;//cycle boundaries
	double *p_ppcf;

#ifdef _VIBR_AMP
	min_ind=thread_arg.min_bin_used;
	max_ind=thread_arg.max_bin_used;
#else
	min_ind=thread_arg.min_bin[0];//reset if necessary
	max_ind=thread_arg.max_bin[0];
#endif

#ifdef _ATLAS
	double *p_ppcf0;
	p_ppcf0=0;//avoid warning
#endif
	
#ifdef _NO_PERIODIC
	if (calcmode)
	{
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
						p_ppcf=finder[ipartial]+min_ind;
						for (ibin=min_ind;ibin<=max_ind;ibin++)
							*p_ppcf++=0;
					}

					loc_offset1=SimpleCfg::cumul[itype]*HistoSet::sum_bins_offset+HistoSet::loc_offset[jtype];//to the beginning of the jtype local hist. of the first
						//atom of itype
					ppos=HistoSet::atom_binind+SimpleCfg::cumul[itype];//pointer to the central bin index where each atom can be found 
					for (i=0;i<SimpleCfg::pnatoms[itype];i++) 
					{
						p_ppcf=finder[ipartial]+min_ind;
						//summing the local density
						for (ibin=min_ind;ibin<=max_ind;ibin++)
						{

							*p_ppcf+= *(hist.local_count+loc_offset1+ibin) / datmat.dV_in[*ppos *nbins[0]+ibin];
							p_ppcf++;
							
						}
						ppos++;
						loc_offset1+=HistoSet::sum_bins_offset;
					}
					
				}//end of if this partial was modified
			}//end of jtype cycle
		}//end of itype cycle

		//averaging the locale density and dividing by rho
		for (ipartial=0;ipartial<npartials;ipartial++)
		{
			if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
						//in case the main thread already generated the new move
			{

				p_ppcf=finder[ipartial]+min_ind;
				for (ibin=min_ind;ibin<=max_ind;ibin++)
				{
					*p_ppcf*=datmat.ntable[ipartial];
					p_ppcf++;
				}
			}
		}
	}
	else
	{

#endif

		for (ib = 0; ib < HistoSet::ndiffbin; ib++)
		{
#ifdef _VIBR_AMP
			min_ind = thread_arg.min_bin_used;
			max_ind = thread_arg.max_bin_used;
#else
			min_ind = thread_arg.min_bin[ib];//reset if necessary
			max_ind = thread_arg.max_bin[ib];
#endif
			for (ipartial = 0; ipartial < npartials; ipartial++)
			{
				if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
							//in case the main thread already generated the new move
				{
					//Calculate the new modified PPCF partials

					//sets the p_ppcf pointer to the beginning of the modified partial in 
					//the PPCFSet gvalues array for this thread
					p_ppcf = *(finder + ib*npartials+ipartial) + min_ind;
#ifdef _ATLAS//ATLAS will be used for Fourier transformation, but ppcf-1 is needed for it
					//sets the p_ppcf0 pointer to the beginning of the modified partial in 
					//the PPCFSet gr_0 array 
					p_ppcf0 = gr_0 + offset[ib*npartials+ipartial] + min_ind;//only used in cas of ATLAS
#endif

#ifdef _VIBR_AMP
					//sets the p_hist pointer to the beginning of the modified partial in 
					//the HistoSet pcounts_T array for this thread
					p_hist = *(hist.finder_T + ipartial) + min_ind;
#else
#ifdef _NO_PERIODIC
					//only arrives here, if nek>0
					pd_hist = hist.periodic_finder[ipartial] + min_ind;
#else
					//sets the p_hist pointer to the beginning of the modified partial in 
					//the HistoSet pcounts array for this thread
					p_hist = *(hist.finder +npartials*ib+ ipartial) + min_ind;
#endif
#endif
					//sets the p_ntable pointer to the beginning of the modified partial in 
					//the DataMat ntable array for this thread
					p_ntable = datmat.ntable_finder[npartials*ib + ipartial] + min_ind;
					for (ibin = min_ind; ibin <= max_ind; ibin++)
					{
#ifndef _NO_PERIODIC
						*p_ppcf = *p_hist++ * *p_ntable++;//normalising the histogram
#else
#ifdef _VIBR_AMP
						*p_ppcf = *p_hist++ * *p_ntable++;//normalising the histogram
#else
						*p_ppcf = *pd_hist++ * *p_ntable++;//normalising the histogram
#endif
#endif

#ifdef _ATLAS//ATLAS will be used for Fourier transformation, but ppcf-1 is needed for it

						//subtract 1 from each ppcf value:
						*p_ppcf0++ = *p_ppcf - 1.0;
#endif
						p_ppcf++;
					}//end of ibin cycle
				}//end of if this partial was modified
			}//end of cycling through the partials ipartial
		}//end of ib cycle
#ifdef _NO_PERIODIC
	}
#endif

}

//Copy the modified parts of the ppcf to target	
void PPCFSet::CopyModified(PPCFSet &target, ThreadArg &thread_arg)
{
	int ipartial,ib,ibin;
	int min_ind,max_ind;//cycle boundaries
	double *p_ppcf,*p_ppcf_target;

	for (ib=0; ib < HistoSet::ndiffbin; ib++)
	{
#ifdef _VIBR_AMP
		min_ind = thread_arg.min_bin_used;
		max_ind = thread_arg.max_bin_used;
#else
		min_ind = thread_arg.min_bin[ib];
		max_ind = thread_arg.max_bin[ib];
#endif
		for (ipartial = 0; ipartial < npartials; ipartial++)
		{
			if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
						//in case the main thread already generated the new move
			{
				//This partial has been modified
				p_ppcf = finder[ib*npartials+ipartial] + min_ind;
				p_ppcf_target = target.finder[ib*npartials+ipartial] + min_ind;
				for (ibin = min_ind; ibin <= max_ind; ibin++)
					*p_ppcf_target++ = *p_ppcf++;
			}//end of if this partial has been modified
		}//end of cycle ipartial
	}//end of ib cycle
}

//Add the present ppcf-s to the sum
void PPCFSet::AddPPCF()
{
	int i;

	//adding the present PPCF to the gvalues_sum array
	int imax=npartials*HistoSet::ntot_bins;//size of storage array
	ncollect++;
	for (i=0;i<imax;i++)
		gvalues_sum[i]+=gvalues[i];//updating the sum

};
