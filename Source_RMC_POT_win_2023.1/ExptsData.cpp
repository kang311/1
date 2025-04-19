//source ExptsData.cpp
//Last changed 28.02.2023

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include"classes2.h"

//defining the static members for ExptsData
int  ExptsData::ngr;//number of g(r) data files
int  ExptsData::nsq;//number of S(Q) data files
int  ExptsData::nfq;//number of F(Q) data files
int  ExptsData::nfg;//number of F(g) data files
int  ExptsData::nek;//number of EXAFS data files E(k)
int  ExptsData::ntot_datasets;//total number of data series
int  ExptsData::ntypes;//number of atom types
int  ExptsData::ndiffbin=1;//number of different bin sizes
int *ExptsData::assign_hist;//pointer to RunParams::assignhist, assigning the different histogram bin sizes to the data sets dim [ntot_datasets]
int	*ExptsData::use_cubic;//array to indicate to use cubic renormalization for the g(r), S(Q) and F(Q) data
double *ExptsData::rspacing;//pointer to Runparams::rspacing

int  ExptsData::rused_tot;//total number of used r points
int *ExptsData::bin_flag;//to indicate, whether the binsize or binshift was reset based on the dr spacing in the experimental gr data set
int *ExptsData::grmin;//indices of first data point to be used (START AT 1!)
int *ExptsData::grmax;//indices of last data point to be used
int *ExptsData::grsize;//number of r data points in g(r) files for each data set
int *ExptsData::grused;//number of used data points in g(r) files for each data set
int *ExptsData::grused_cum;//cumulative number of USED g(r) data points in previous data sets; size:ngr+1
int *ExptsData::grrenorm;//renormalization switches for g(r) files
int *ExptsData::groffset;//offset switches for g(r) files
int *ExptsData::grlinear;//linear background switches for g(r) files
int *ExptsData::grquadratic;//quadratic background switches for g(r) files
int *ExptsData::grcubic;//cubic background switches for g(r) files
int *ExptsData::gruseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default

int  ExptsData::sqused_tot;//total number of used S(Q) points
int *ExptsData::sqmin;//indices of first data point to be used (START AT 1!)
int *ExptsData::sqmax;//indices of last data point to be used
int *ExptsData::sqsize;//number of data points in S(Q) files for each data set
int *ExptsData::sqused;//number of USED data points in S(Q) files for each data set
int *ExptsData::sqused_cum;//cumulative number of USED S(Q) data points in previous data sets; size:nsq+1
int	*ExptsData::sqoffset;//offset switches for S(Q) files
int	*ExptsData::sqrenorm;//renormalization switches for S(Q) files
int *ExptsData::sqlinear;//linear background switches for S(Q) files
int *ExptsData::sqquadratic;//quadratic background switches for S(Q) files
int *ExptsData::sqcubic;//cubic background switches for S(Q) files
int *ExptsData::squseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
int *ExptsData::sqreadcoeffs;//wherther to read the coefficients instead of calculating them
int ExptsData::is_Ncoeff_calc=0;//indicator, whether neutron coeff calculation is used at all 
isotope_count_type *ExptsData::isotope_count;//[nfq*ntypes] number of isotopes for each data set's each component

int  ExptsData::fqused_tot;//total number of used F(Q) points
int *ExptsData::fqmin;//indices of first data point to be used (START AT 1!)
int *ExptsData::fqmax;//indices of last data point to be used
int *ExptsData::fqsize;//number of data points in F(Q) files for each data set
int *ExptsData::fqused;//number of USED data points in F(Q) files for each data set
int *ExptsData::fqused_cum;//cumulative number of USED F(Q) data points in previous data sets; size:nfq+1
int	*ExptsData::fqoffset;//offset switches for F(Q) files
int	*ExptsData::fqrenorm;//renormalization switches for F(Q) files
int *ExptsData::fqlinear;//linear background switches for F(Q) files
int *ExptsData::fqquadratic;//quadratic background switches for F(Q) files
int *ExptsData::fqcubic;//cubic background switches for F(Q) files
int *ExptsData::fqrenalpha;//compton alpha correction switches for F(Q) files
int *ExptsData::fqreadcoeffs;//wherther to read the coefficients instead of calculating them
int *ExptsData::fquseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
int *ExptsData::fqAXS;//0: no AXS, >0 index (STARTING WITH 1) of the anomalous X-ray scattering type/data set
int *ExptsData::fqfprimeindex;//sum_fprime_itype(2**itype) for each set, indicating which types are using f'
int *ExptsData::fqfprimecount;//sum of how many types using f' for the given set
int ExptsData::is_IQ = 0;//indicator, whether I(Q) fitting is used at all,0 no, 1 yes with linear reg, 2 there is non-lin regression
int ExptsData::is_IQmucorr = 0;//indicator, whether there is a mu correction in any of the I(Q) fitting sets
int ExptsData::is_Xcoeff_calc=0;//indicator, whether Xray coeff calculation is used at all 
int ExptsData::is_AXS = 0;//indicator, whether AXS is used
double *ExptsData::fqalpha;//alpha renormalization parameter in I(Q) fit
double *ExptsData::fqa;//a renormalization parameter in I(Q) fit
double *ExptsData::fqb;//b renormalization parameter in I(Q) fit
double *ExptsData::fqmu;//parameter for I(Q) background correction
double *ExptsData::fqdmumax;//max value to vary mu for I(Q) fitting, between fqmu - fqdmumax -> fqmu + fqdmumax
double *ExptsData::fqmuact;//actual value of the mu parameterm between fqmu - fqdmumax -> fqmu + fqdmumax
double *ExptsData::fqfprime;//real part of the f'(E) for each data satom type
double *ExptsData::fqfprimeact;//real part of the actual f'(E) for each atom type
double *ExptsData::fqfprimefactor;//for varying f', values are in the range f'-f'*fqpfactor->f'+f'*fqpfactor
bool *ExptsData::fqfitIQ;//Fit I(Q) instead of F(Q), for each F(Q) data set
bool *ExptsData::fqIQbackgcorr;//use background correction only in case of I(Q) fitting, for each F(Q) data set
bool *ExptsData::fqusecompton;//whether to add the Compton contribution during the I(Q) calculation for the given data set

int  ExptsData::fgused_tot;//total number of used F(g) points
int *ExptsData::fgmin;//indices of first data point to be used (START AT 1!)
int *ExptsData::fgmax;//indices of last data point to be used (START AT 1!)
int *ExptsData::fgsize;//number of data points in F(g) files
int *ExptsData::fgused;//number of used data points in F(g) files
int *ExptsData::fgused_cum;//cumulative number of USED F(g) data points in previous data sets; size:nfg+1
int *ExptsData::fgoffset;//offset switches for F(g) files
int *ExptsData::fgrenorm;//renormalization switch for F(g) files
int *ExptsData::fglinear;//linear background switches for F(g) files
int *ExptsData::fgquadratic;//quadratic background switches for F(g) files
int *ExptsData::fgcubic;//cubic background switches for F(g) files
int *ExptsData::fguseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default

int  ExptsData::ekused_tot;//total number of used E(k) points
int  ExptsData::ek_rused_tot;//total number of r points used in E(k) fitting
int  ExptsData::is_E0shift=0;//indicator whether any EXAFS set uses E0_shift
int *ExptsData::ekmin;//indices of first data point to be used (START AT 1!)
int *ExptsData::ekmax;//indices of last data point to be used
int *ExptsData::ekmin_ori;//indices of first data point to be used originally (needed in case of E0 shift) (START AT 1!)
int *ExptsData::ekmax_ori;//indices of last data point to be used originally (needed in case of E0 shift) (START AT 1!)
int *ExptsData::eksize;//number of data points in E(k) files for each data set
int *ExptsData::ekused;//number of USED data points in E(k) files for each data set
int *ExptsData::ekused_cum;//cumulative number of USED E(k) data points in previous data sets; size:nek+1
int	*ExptsData::ekoffset;//offset switches for E(k) files
int	*ExptsData::ekrenorm;//renormalization switches for E(k) files
int *ExptsData::ek_rused;//number of r points to use
int *ExptsData::ek_rmin;//index of first r point (histogram bin) to use
int *ExptsData::ek_rmax;//index of last r point (histogram bin) to use
int *ExptsData::ek_abstype;//type of the absorbing particle
int *ExptsData::ekchipower;//weighting factor (E(k)*(k power of ekchipower)) for the experimental data points
int *ExptsData::ekuseR;//whether to use R factor instead of chisqure to calculate the squared diff for this set, 0 is default
int *ExptsData::ek_ngrid_in;//number of grid points to read from *.dat including the maximum for k-shift (0-: no shift, 1:only max_shift...)
int *ExptsData::ek_ngrid;//total number of grid points, ek_gridto-ek_gridfrom+1
int *ExptsData::ek_gridind;//the index of the actually chosen grid point during the E0 shift for each E(k) data series 
int *ExptsData::ek_gridstart;//the index of the starting grid point, default is zero (no shift)
int *ExptsData::ek_gridfrom;//index of the first grid point to use in case of assymetric grid (between -NGRID -> +NGRID)
int *ExptsData::ek_gridto;//index of the last grid point to use in case of assymetric grid (between -NGRID -> +NGRID)
int *ExptsData::ek_cf_cum;//cumulative number of coeff matrices for each data set  
int *ExptsData::ek_ngrid_cum;//cumulativ number of gridpoints for each set

double *ExptsData::grsub;//constants to be substracted from data for g(r) files
double *ExptsData::grsigma;//sigma values  for g(r) files for normal chi2
double *ExptsData::gr_coeffs;//partials coefficients of g(r) files

double	*ExptsData::sqsub;//constants to be substracted from data for S(Q) files
double	*ExptsData::sqsigma;//sigma values for S(Q) files for normal chi2
double	*ExptsData::sq_coeffs;//partials coefficients of S(Q) files

double	*ExptsData::fqsub;//constants to be substracted from data for F(Q) files
double	*ExptsData::fqsigma;//sigma values for F(Q) files for normal chi2

double *ExptsData::fgsub;//constants to be subsstracted from data for F(g) files
double *ExptsData::fgsigma;//sigma values for F(g) files


double *ExptsData::eksigma;//sigma values for E(k) files for normal chi2
double *ExptsData::ek_dE0;//maximum E0-shift for both direction
double *ExptsData::ek_max_dkl;//maximum k-shift left for each data point
double *ExptsData::ek_max_dkr;//maximum k-shift right for each data point
double *ExptsData::ek_dE0_grid;//dE0 values of the dE0 grid
double *ExptsData::dE0_actual;//the actual dE0 values belonging to each grid point, this is  not equdistantly set 


char   (*ExptsData::datafilename)[FILE_NAME_SIZE];//names of the experimental data files

//default constructor of the ExptsData class 
//The static members must be first initialised by ExptsData::GetExptsParam before called
ExptsData::ExptsData(RunParams &rundat)
{
	//rundat is needed to check how the rmax for the histogram calculation relates to the rmax
	//of the last experimental data point, and reset them, if necessary
	
	ifstream file;//used for the experimental data file
	int i,j,k,etype,ipartial,itype,nexpts,ntot,npartials,mybins;
#ifdef _NOP_VIBR
	std::streamoff first_r_pos,prev_pos;//to position the get pointer in the file
	first_r_pos = 0;//default
	prev_pos = 0;//default
	double rhist_min, rhist_max;//the central values of the first and last histogram bins
	//the mean values of the first and last histogram bins to be used for this set
	
#endif
	int ecoeff_tot;//number of sofar ued EXAFS coeffs for the finder initialization
	double *pr,*pg,*pq,*ps,*pf, *pf_ori=0,*pIB=0, *pe,*pc,*pk;//pointers used to read data files r, g(r), Q, S(Q), F(Q), E(k) values
	
	char *base_name,*base_name2,*base_name3, *base_name4, *base_name5, *conv_numb=NULL;

	SetArraysize(&base_name, NAME_SIZE, "base_name","ExptsData::ExptsData");
	SetArraysize(&base_name2, NAME_SIZE, "base_name2","ExptsData::ExptsData");
	SetArraysize(&base_name3, NAME_SIZE, "base_name3", "ExptsData::ExptsData");
	SetArraysize(&base_name4, NAME_SIZE, "base_name4", "ExptsData::ExptsData");
	SetArraysize(&base_name5, NAME_SIZE, "base_name5", "ExptsData::ExptsData");
		
	npartials=ntypes*(ntypes+1)/2;//number of partials
	nexpts = ngr + nsq + nfq + nfq + nek;//number of experimental sets
	
	//Checking whether the static members were initialised
	if (nexpts==0 && ::debug)
	{
		cout << "\nWARNING(" << ++warn << "): ExptsData constructor"<<endl;
		cout<<"\tThere is no experimental data!"<<endl;
	}
	
	//Loading the experimental data sets
	cout << "\nNOTE(" << ++note << "): Number of data sets for loading experimental data from: "<<ngr+nsq+nfq+nfg+nek<<endl;

	//g(r) data sets
	if (ngr>0)
	{

	
		//Allocating memory
		SetArraysize(&gr_rvalues,rused_tot,"gr_rvalues","ExptsData::ExptsData");
		SetArraysize(&gr_gvalues,rused_tot,"gr_gvalues","ExptsData::ExptsData");
		SetArraysize(&gr_rfinder,ngr,"gr_rfinder","ExptsData::ExptsData");
		SetArraysize(&gr_gfinder,ngr,"gr_gfinder","ExptsData::ExptsData");
				
		//reset to zero, to calculate it again in case the expt dat file has less data point than in the .dat file
		rused_tot=0;
		//loading the data
		for (i = 0; i < ngr; i++)
		{
			cout << "\nLoading data from g(r) file: " << datafilename[i] << endl;
			SafeOpenTextFile(file, datafilename[i]);
			if (CheckFileState(file, "ExptsData::ExptsData", datafilename[i]) == 0)
			{
				cout << "Cannot load data, exiting..." << endl;
				CleanExit();
			};

			ntot = ReadThisLine(file, 1, 1, "g(r) ntot", "ExptsData::ExptsData", datafilename[i]);//total number of r points in the data file
			

			grused[i] = grmax[i] - grmin[i] + 1;//number of r file points to be used
		

			if (grused[i] < 0)// can happen, if the number of data points are not given in the file
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The number of used data points is smaller, than zero, something must be wrong with the format of the file!" << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}

			SkipLine(file, datafilename[i], 1);//comment line

			for (j = 1; j < grmin[i]; j++)
				SkipLine(file, datafilename[i], 1);//skipping all the unused points

			//Setting the finders for this data set
			gr_rfinder[i] = gr_rvalues + rused_tot;//setting the finder for r-values
			gr_gfinder[i] = gr_gvalues + rused_tot;//setting the finder for g-values

			//  This checks are not needed, if custom histogram will be used, if necessary for each g(r) set to accomodate the range involved in 
			//calculation, only in case of _NO_PERIODIC or _VIBR_AMP, where there are no custom bins
#ifdef _NOP_VIBR
			rhist_min = rundat.xmin[assign_hist[i]] * rundat.boxedge + (0.5 + RunParams::firstbin[i])*rundat.rspacing[i];
			rhist_max = rundat.xmin[assign_hist[i]] * rundat.boxedge + (rundat.nbins[assign_hist[i]] - 0.5)*rundat.rspacing[i];
			pr = gr_rfinder[i];//pointer for r values
			//First just the r values will be loaded to perform some necessary checking, which might
			//mean the resetting of the number of used data points
			int first_point, last_point;//the index of the first and last used data points after checking
			std::streamoff current_pos;
			first_point = -1;//default
			last_point = grused[i] - 1;//default
			first_r_pos = file.tellg();//save the position of the first r (current file position)
			
			for (j = 0; j < grused[i]; j++)//for all used r points
			{
				current_pos = file.tellg();//save the current file position

				*pr = ReadThisLine(file, 1, 1.0, "r for g(r)", "ExptsData::ExptsData", datafilename[i]);//loads the r value
				pr++;

				//The rmin for the data set has to be compared with the rmin for the histogram
				//calculation. 
				//If rmin_histogram is larger, than the rmin_experimental, then 
				//the experimental rmin, and the number of data points has to be reset, to
				//give a positive matrix element for the first expt data point in the convert_bin_dr matrix
				//This has to be done not just to save time, but otherwise the calculated g(r) 
				//will contains 0 at the smallest r values, and this part will be included in the 
				//calculation of chi-square
				if (*(pr - 1) - GRID_TOL > rhist_min)
				{

					//the previous data point (if there is any) should be the first data point
					if (first_point == -1 && j > 0)
					{
						//if this is not the second r data point, the number of data points has to be reset
						first_point = j - 1;//the index of the first data point to use	
						first_r_pos = prev_pos;//save the position of the first used r data point in the file
					}
				}
				prev_pos = current_pos;//stores the position of the r already read 

				//The rmax for the data set has to be compared with the rmax for the histogram
				//calculation. 
				//If rmax_histogram is smaller, than the rmax_experimental, then 
				//the experimental rmax, and the number of data points has to be reset, to
				//give a positive matrix element for the last exp data point in the convert_bin_dr marix
				//This has to be done not just to save time, but otherwise the calculated g(r) 
				//will contains 0 at the large r values, and this part will be included in the 
				//calculation of chi-square
				if (*(pr - 1) + GRID_TOL >= rhist_max)
				{
					//this should be the last r data point
					if (j < grused[i] - 1)
					{
						//this is not the last r data point, the number of data points has to be reset
						last_point = j;
						break;//leave the r reading cycle
					}
				}
			}//end of r checking cycle j

			if (first_point > 0 || last_point < grused[i] - 1)
			{
				cout << "\nWARNING(" << ++warn << "): The number of data point for data set " << i + 1 << " had to be decreased to " << \
					grused[i] - first_point - (grused[i] - 1 - last_point);
				cout << " from " << grused[i] << " data " << endl;
				cout << "\tpoints, by cutting " << first_point << " at the beginning and " << grused[i] - 1 - last_point << " at the end to match the maximal" << endl;
				cout<<"\tcalculated distance of the histogram (" << rundat.rmax << "A) !" << endl;
				grused[i] -= first_point + (grused[i] - 1 - last_point);
				grmin[i] = first_point + 1;
				grmax[i] = grused[i] + grmin[i] - 1;//resetting the index of the last used data point					
			}
			file.clear(ios::goodbit);

			file.seekg(first_r_pos,ios::beg);//resets the position in the file for the reading of the g(r) data
#endif

			pr=gr_rfinder[i];//pointer for r values
			pg=gr_gfinder[i];//position for the g values
		
			for(j=0;j<grused[i];j++)//for all used points
			{
				IntToStr(&conv_numb, j + 1);
				mystrcpy(base_name, NAME_SIZE, "g(r) series used r point ");
				mystrcpy(base_name2, NAME_SIZE, "g(r) used point ");
				mystrcat(base_name, NAME_SIZE, conv_numb);
				mystrcat(base_name2, NAME_SIZE, conv_numb);
				*pr=ReadThisLine(file,3,1.0,base_name, "ExptsData::ExptsData");//reloads the r value

				pr++;
				if (j==grused[i]-1)
					//do not give warning, if end of file was reached during line skipping after data reading, as this can be the last line
					*pg=ReadThisLine(file,1,1.0,base_name2, "ExptsData::ExptsData");//loads the g(r) value for the i-th expt, j-th r 
				else
					//give warning, if end of file was reached during line skipping after data reading
					*pg=ReadThisLine(file,1,1.0,base_name2, "ExptsData::ExptsData",datafilename[i]);//loads the g(r) value for the i-th expt, j-th r 

				(*pg)-=grsub[i];//substracting the constant
				pg++;
			}
		
			if (!CheckReadFileState(file,"ExptsData constructor",datafilename[i]))
				CleanExit();
			next_line_pos = -1;
			file.close();
			rused_tot+=grused[i];//increment the g(r)-points counter
		}//next g(r) constraint
		//setting the cumulative number of used g(r) points in the previous data sets for each set
		grused_cum[0]=0;
		for (i=1;i<=ngr;i++)
			grused_cum[i]=grused_cum[i-1]+grused[i-1];
	}//end of if there are g(r) constraints)

	
	//now S(Q) constraints
	if(nsq>0)//if there are some S(Q) constraints
	{
		//Allocting memory space to arrays
		SetArraysize(&sq_qvalues,sqused_tot,"sq_qvalues","ExptsData::ExptsData");
		SetArraysize(&sq_svalues,sqused_tot,"sq_svalues","ExptsData::ExptsData");
		SetArraysize(&sq_qfinder,nsq,"sq_qfinder","ExptsData::ExptsData");
		SetArraysize(&sq_sfinder,nsq,"sq_sfinder","ExptsData::ExptsData");

		//reset to zero,to calculate it again in case the expt dat file has less point in than the .dat file
		sqused_tot=0;
	
		pq=sq_qvalues;//pointer for Q values
		ps=sq_svalues;//pointer for S values
		//loading the data	
		for(i=0;i<nsq;i++)//for each S(Q) constraints
		{
			cout<<"\nLoading data from S(Q) file: "<<datafilename[ngr+i]<<endl;
			SafeOpenTextFile(file,datafilename[ngr+i]);
			if (CheckFileState(file,"ExptsData::ExptsData",datafilename[ngr+i])==0)
			{
				cout<<"Cannot load data, exiting..."<<endl;
				CleanExit();
			};
					
			ntot=ReadThisLine(file,1,1,"S(Q) ntot", "ExptsData::ExptsData",datafilename[i]);//total number of Q points in the data file
			//if the maximum number of Q points to be used from the .dat file is greater than the number of Q points in the file
			if(ntot<sqmax[i])
			{
				cout << "\nWARNING(" << ++warn << "): The largest index of the Q data points to be used exceeds the total number of Q points in the file! ";
				cout<<"\tMaximal number of Q points to be used is "<<sqmax[i]<<" from the .dat file, "<<ntot<<\
					" from the experimental data file "<<endl;
				cout<<"\tSetting the maximal number of Q points to be used to "<<ntot<<endl;
				sqmax[i]=ntot;
			}
			sqused[i]=sqmax[i]-sqmin[i]+1;//number of file points to be used
			sqsize[i]=ntot;

			SkipLine(file,datafilename[i],1);//comment line
			
			//Setting the finders for this data set
			sq_qfinder[i]=sq_qvalues+sqused_tot;//setting the finder for Q-values
			sq_sfinder[i]=sq_svalues+sqused_tot;//setting the finder for S-values

			for(j=1;j<sqmin[i];j++)
				SkipLine(file,datafilename[i],1);//skipping all the unused Q points
			//The values will be stored the following way in the sq_svalues array:
			//first data set:S(Q1)-S(Qmax),second data set:S(Q1)-S(Qmax)...
			for(j=0;j<sqused[i];j++)//for all used points
			{
				IntToStr(&conv_numb, j + 1);
				mystrcpy(base_name, NAME_SIZE, "S(Q) series used Q point ");
				mystrcpy(base_name2, NAME_SIZE, "S(Q) used point ");
				mystrcat(base_name, NAME_SIZE, conv_numb);
				mystrcat(base_name2, NAME_SIZE, conv_numb);
				*pq=ReadThisLine(file,3,1.0,"base_name", "ExptsData::ExptsData");//loads the Q value
				pq++;
				if (j<sqused[i]-1)
					//give warning, if end of file was reached during line skipping after data reading
					*ps=ReadThisLine(file,1,1.0,"base_name2", "ExptsData::ExptsData",datafilename[i]);//loads the S(Q) value for the i-th expt, j-th Q value
				else
					//do not give warning, if end of file was reached during line skipping after data reading, as this can be the last line
					*ps=ReadThisLine(file,1,1.0,"base_name2", "ExptsData::ExptsData");//loads the g(r) value for the i-th expt, j-th r 
				
				*ps-=sqsub[i];//substracting the constant
				ps++;
			}
		
			
			if (!CheckReadFileState(file,"ExptsData constructor",datafilename[ngr+i]))
				CleanExit();
			next_line_pos = -1;
			file.close();
			sqused_tot+=sqused[i];//increment the data points counter
	
			if (sqreadcoeffs[i] == 0)
				CalcNDCoeffs(i);
		}//next S(Q) constraint 
		
		//setting the cumulative number of used S(Q) points in the previous data sets for each set
		sqused_cum[0]=0;
		for (i=1;i<=nsq;i++)
			sqused_cum[i]=sqused_cum[i-1]+sqused[i-1];
		
		
	}//end of if there are S(Q) constraints)
	
	//now F(Q) constraints
	if(nfq>0)//if there are some F(Q) constraints
	{
		//Allocting memory space to arrays
		SetArraysize(&fq_qvalues,fqused_tot,"fq_qvalues","ExptsData::ExptsData");
		SetArraysize(&fq_fvalues,fqused_tot,"fq_fvalues","ExptsData::ExptsData");
		SetArraysize(&fq_coeffs,fqused_tot*npartials,"fq_coeffs","ExptsData::ExptsData");
		SetArraysize(&fq_qfinder,nfq,"fq_qfinder","ExptsData::ExptsData");
		SetArraysize(&fq_ffinder,nfq,"fq_ffinder","ExptsData::ExptsData");
		SetArraysize(&fq_cfinder,nfq*npartials,"fq_cfinder","ExptsData::ExptsData");
		if (is_Xcoeff_calc || is_IQ > 0)//atomic scattering factors will be needed
		{
			SetArraysize(&fq_scatfactors, fqused_tot *ntypes, "fq_scatfactors", "ExptsData::ExptsData");
			SetArraysize(&fq_sffinder, nfq *ntypes, "fq_scatfactors", "ExptsData::ExptsData");
		}

		if (is_IQ>0)
		{
			//If I(Q) fit is used for any data set, place will be allocated for all the F(Q) sets
			//Anyway it is more likely, that all the F(Q) sets will be handled similarly
			SetArraysize(&fq_A, fqused_tot, "fq_A", "ExptsData::ExptsData");
			SetArraysize(&fq_B, fqused_tot, "fq_B", "ExptsData::ExptsData");
			SetArraysize(&fq_ffromIQrenvalues, fqused_tot, "fq_ffromIQrenvalues", "ExptsData::ExptsData");
			SetArraysize(&fq_ffromIQrenfinder, nfq, "fq_ffromIQrenfinder", "ExptsData::ExptsData");
		}
		if (is_IQmucorr > 0)//only used in case of I(Q) mu background correction
		{
			SetArraysize(&fq_IB, fqused_tot, "fq_IB", "ExptsData::ExptsData");
			SetArraysize(&fq_IBfinder, nfq, "fq_IBfinder", "ExptsData::ExptsData");
			SetArraysize(&fq_IQvalues_ori, fqused_tot, "fq_IQvalues_ori", "ExptsData::ExptsData");
			SetArraysize(&fq_ffinder_ori, nfq, "fq_ffinder_ori", "ExptsData::ExptsData");
		}
		//reset to zero,to calculate it again in case the expt dat file has less point in than the .dat file
		fqused_tot=0;
		fqused_cum[0] = 0;
	
		pq=fq_qvalues;//pointer for Q values
		pf=fq_fvalues;//pointer for F values (even if I(Q) is fitted without background correction, the same arrays used as normally
		if (is_IQmucorr)
			pf_ori = fq_IQvalues_ori;//pointer to I values in case of background correction, the original experimental data is stored here, while the
		     //fq_fvalues will be used for the actual corrected I(Q) used in the chi2 and elsewhere
		//loading the data	
		for(i=0;i<nfq;i++)//for each F(Q) constraints
		{
			cout<<"\nLoading data from F(Q) file: "<<datafilename[ngr+nsq+i]<<endl;
			SafeOpenTextFile(file,datafilename[ngr+nsq+i]);
			if (CheckFileState(file,"ExptsData::ExptsData",datafilename[ngr+nsq+i])==0)
			{
				cout<<"Cannot load data, exiting..."<<endl;
				CleanExit();
			};
					
			ntot=ReadThisLine(file,1,1,"F(Q) ntot", "ExptsData::ExptsData",datafilename[i]);//total number of Q points in the data file
			//if the maximum number of Q points to be used from the .dat file is greater than the number of Q points in the file
			if(ntot<fqmax[i])
			{
				cout << "\nWARNING(" << ++warn << "): The largest index of the Q data points to be used exceeds the total number of Q points in the file! ";
				cout<<"\tMaximal number of Q points to be used is "<<fqmax[i]<<" from the .dat file, "<<ntot<<\
					" from the experimental data file "<<endl;
				cout<<"\tSetting the maximal number of Q points to be used to "<<ntot<<endl;
				fqmax[i]=ntot;
			}
			fqused[i]=fqmax[i]-fqmin[i]+1;//number of file points to be used
			fqsize[i]=ntot;

			SkipLine(file,datafilename[i],1);//comment line
			
			//Setting the finders for this data set
			fq_qfinder[i]=fq_qvalues+fqused_tot;//setting the finder for Q-values
			fq_ffinder[i]=fq_fvalues+fqused_tot;//setting the finder for F-values
			if (is_IQ)
				fq_ffromIQrenfinder[i]= fq_ffromIQrenvalues + fqused_tot;//setting the finder for fq_ffromIQrenvalues, this will contain the F(Q) calculated back from the renormalized ext I(Q)
			if (is_IQmucorr)
			{
				fq_ffinder_ori[i] = fq_IQvalues_ori + fqused_tot;//setting the finder for original I(Q) values
				fq_IBfinder[i] = fq_IB + fqused_tot;//setting the finder for I(Q) background values
			}
			for (ipartial=0;ipartial<npartials;ipartial++)//fqused_tot here only contains the used Q points of the previous data sets
				fq_cfinder[i*npartials+ipartial]=fq_coeffs+fqused_tot*npartials+fqused[i]*ipartial;//setting the finder for the Q dependent coefficients

			if (is_Xcoeff_calc || is_IQ)
			{
				for (itype=0;itype<ntypes;itype++)
					fq_sffinder[i * ntypes + itype] = fq_scatfactors + fqused_tot * ntypes + fqused[i] * itype;//setting the finder for the Q dependent scattering factors

			}
			for(j=1;j<fqmin[i];j++)
				SkipLine(file,datafilename[i],1);//skipping all the unused Q points

			//The values will be stored the following way in the fq_fvalues and if necessary fq_IQvalues_ori array:
			//first data set:F(Q1)-F(Qmax),second data set:F(Q1)-F(Qmax)...
			if (fqfitIQ[i])
				pIB = fq_IB + fqused_cum[i];
			
			for (j = 0; j < fqused[i]; j++)//for all used points
			{
				IntToStr(&conv_numb, j + 1);
				if (fqfitIQ[i])
				{
					mystrcpy(base_name, NAME_SIZE, "I(Q) series used Q point ");
					mystrcpy(base_name2, NAME_SIZE, "I(Q) used point ");
					mystrcpy(base_name3, NAME_SIZE, "IB(Q) used point ");
					mystrcat(base_name3, NAME_SIZE, conv_numb);
					mystrcpy(base_name4, NAME_SIZE, "X-ray scattering factor used point ");
					mystrcat(base_name4, NAME_SIZE, conv_numb);
					mystrcat(base_name4, NAME_SIZE, " type ");

				}
				else
				{
					mystrcpy(base_name, NAME_SIZE, "F(Q) series used Q point ");
					mystrcpy(base_name2, NAME_SIZE, "F(Q) used point ");
					mystrcpy(base_name4, NAME_SIZE, "X-ray coeff used point ");
					mystrcat(base_name4, NAME_SIZE, conv_numb);
					mystrcat(base_name4, NAME_SIZE, " partial ");
				}
				mystrcat(base_name, NAME_SIZE, conv_numb);
				mystrcat(base_name2, NAME_SIZE, conv_numb);
				
				*pq = ReadThisLine(file, 3, 1.0, base_name, "ExptsData::ExptsData");//loads the Q value
				pq++;
				if (!fqIQbackgcorr[i])
				{
					*pf = ReadThisLine(file, 3, 1.0, base_name2, "ExptsData::ExptsData");//loads the F(Q) value for the i-th expt, j-th Q value
					*pf -= fqsub[i];//substracting the constant
				}
				else//constant will not be used, only scalable background correction, which is done later 
					*pf_ori = ReadThisLine(file, 3, 1.0, base_name2, "ExptsData::ExptsData");//loads the I(Q) value for the i-th expt, j-th Q value

				if (fqfitIQ[i] && ((fqIQbackgcorr[i]==0 && *pf < 0) || (fqIQbackgcorr[i] && *pf_ori < 0)))
				{
					cout << "\n*****ERROR*****" << endl;
					
					if (fqIQbackgcorr[i])
						cout << "The " << j + 1 << ". original I(Q) value of the "<<i+1<<". F(Q) data set is negative ( " << *pf_ori << endl;
					else
					{
						cout << "The " << j + 1 << ". I(Q) value of the " << i + 1 << ". F(Q) data set is negative ( " << *pf;
						cout << " ) after the subtraction of the constant!" << endl;
					}
					cout << "The intensity cannot be negative, correct your data file and/or the value of the constant, and try again!" << endl;
					cout << "Exiting..." << endl;
					CleanExit();
				}
				
				if (fqreadcoeffs[i])//read the coeffs or scattering factors from file
				{
					if (fqfitIQ[i] == false)//F(Q) fit, partial coeffs are read
					{
						for (ipartial = 0; ipartial < npartials; ipartial++)
						{
							mystrcpy(base_name5, NAME_SIZE, base_name4);
							IntToStr(&conv_numb, ipartial + 1);
							mystrcat(base_name5, NAME_SIZE, conv_numb);
							if (ipartial < npartials - 1)
								*(fq_cfinder[i * npartials + ipartial] + j) = ReadThisLine(file, 3, 1.0, base_name5, "ExptsData::ExptsData");//reading the coefficients
							else
							{
								if (j < fqused[i] - 1)//give warning, if end of file was reached during line skipping after data reading
									*(fq_cfinder[i * npartials + ipartial] + j) = ReadThisLine(file, 1, 1.0, base_name5, "ExptsData::ExptsData", datafilename[i]);
								else//do not give warning, if end of file was reached during line skipping after data reading, as this can be the last line
									*(fq_cfinder[i * npartials + ipartial] + j) = ReadThisLine(file, 1, 1.0, base_name5, "ExptsData::ExptsData");
							}
						}
					}
					else//I(Q) fitted, f(Q) is read for all types
					{
						for (itype = 0; itype < ntypes; itype++)
						{
							mystrcpy(base_name5, NAME_SIZE, base_name4);
							IntToStr(&conv_numb, itype + 1);
							mystrcat(base_name5, NAME_SIZE, conv_numb);
							if (itype < ntypes - 1)
								*(fq_sffinder[i * ntypes + itype] + j) = ReadThisLine(file, 3, 1.0, base_name5, "ExptsData::ExptsData");//reading the coefficients
							else
							{
								if (j < fqused[i] - 1)//give warning, if end of file was reached during line skipping after data reading
									*(fq_sffinder[i * ntypes + itype] + j) = ReadThisLine(file, 3, 1.0, base_name5, "ExptsData::ExptsData", datafilename[i]);
								else//do not give warning, if end of file was reached during line skipping after data reading, as this can be the last line
									*(fq_sffinder[i * ntypes + itype] + j) = ReadThisLine(file, 3, 1.0, base_name5, "ExptsData::ExptsData");
							}
						}

					}
				}
				else
					if (!fqfitIQ[i])
						SkipLine(file);//go to next line
				if (fqfitIQ[i])
				{
					if (fqIQbackgcorr[i])
					{
						if (j < fqused[i] - 1)
							*pIB = ReadThisLine(file, 1, 1.0, base_name3, "ExptsData::ExptsData", datafilename[i]);//loads the IB(Q) factor for I(Q) fitting
						else
							*pIB = ReadThisLine(file, 1, 1.0, base_name3, "ExptsData::ExptsData");//loads the IB(Q) factor for I(Q) fitting
						*pf = *pf_ori - fqmuact[i] * *pIB;//has to use fqmuact because of continuation this will only hold the actual mu
						
						if (*pf_ori-( fqmu[i] + fqdmumax[i])* *pIB < 0)//the maximum subtraction has to be used for testing!
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The " << j + 1 << ". corrected I(Q) value of the " << i + 1 << ". F(Q) data set after the subtraction of the (mu+dmu_max)*background is negative ( " << *pf_ori - (fqmu[i] + fqdmumax[i])* *pIB << ")" << endl;
							cout << "The intensity cannot be negative, correct your data file and/or the mu and the background, and try again!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
						pIB++;

					}
					else
						SkipLine(file);
				}
				pf++;//increase the pointer regardless whether it was used for this set
				if (is_IQmucorr)
					pf_ori++;//if the array is created increase the pointre to be at the correct place, if needed regardless whether it was used here

			}
			if (!CheckReadFileState(file,"ExptsData constructor",datafilename[ngr+nsq+i]))
				CleanExit();
			next_line_pos = -1;
			file.close();
			if (!fqreadcoeffs[i])
				CalcXrayCoeffs(i);//calc scattering fators, coefficients, A(Q) if I(Q) is fitted
			else
				if (fqfitIQ[i])
					CalcXrayCoeffs(i,1);//only calculate the coefficients and A(Q)
			if (fqfitIQ[i])
				CalcCompton(i);//calculate the Compton B(Q) if necessary, otherwise fill the array with zero
					//no cheking will be performed later, Compton term will be added in case of I(Q) fitting regardless the content

			fqused_tot+=fqused[i];//increment the data points counter
			//setting the cumulative number of used F(Q) points in the previous data sets for each set
			fqused_cum[i+1] = fqused_cum[i] + fqused[i];
		}//next F(Q) constraint 
		
		
	}//end of if there are F(Q) constraints)
	//now F(g) constraints
	if (nfg > 0)//if there are some F(g) constraints
	{
		//Allocting memory space to arrays
		SetArraysize(&fg_gvalues, fgused_tot, "fg_gvalues", "ExptsData::ExptsData");
		SetArraysize(&fg_fvalues, fgused_tot, "fg_fvalues", "ExptsData::ExptsData");
		SetArraysize(&fg_coeffs, fgused_tot*npartials, "fg_coeffs", "ExptsData::ExptsData");
		SetArraysize(&fg_gfinder, nfg, "fg_gfinder", "ExptsData::ExptsData");
		SetArraysize(&fg_ffinder, nfg, "fg_ffinder", "ExptsData::ExptsData");
		SetArraysize(&fg_cfinder, nfg*npartials, "fg_cfinder", "ExptsData::ExptsData");

		//reset to zero,to calculate it again in case the expt dat file has less point in than the .dat file
		fgused_tot = 0;

		pg = fg_gvalues;//pointer for g values
		pf = fg_fvalues;//pointer for F values
		//loading the data	
		for (i = 0; i < nfg; i++)//for each F(g) constraints
		{
			cout << "\nLoading data from F(g) file: " << datafilename[ngr + nsq + nfq + i] << endl;
			SafeOpenTextFile(file, datafilename[ngr + nsq + nfq + i]);
			if (CheckFileState(file, "ExptsData::ExptsData", datafilename[ngr + nsq + nfq + i]) == 0)
			{
				cout << "Cannot load data, exiting..." << endl;
				CleanExit();
			};

			ntot = ReadThisLine(file, 1, 1, "F(g) ntot", "ExptsData::ExptsData", datafilename[i]);//total number of g points in the data file
			//if the maximum number of g points to be used from the .dat file is greater than the number of g points in the file
			if (ntot < fgmax[i])
			{
				cout << "\nWARNING(" << ++warn << "): The largest index of the g data points to be used exceeds " << endl;
				cout<<"\tthe total number of g points in the file!";
				cout << "\tMaximal number of g points to be used is " << fgmax[i] << " from the .dat file, " << ntot << \
					" from the experimental data file " << endl;
				cout << "\tSetting the maximal number of g points to be used to " << ntot << endl;
				fgmax[i] = ntot;
			}
			fgused[i] = fgmax[i] - fgmin[i] + 1;//number of file points to be used
			fgsize[i] = ntot;

			SkipLine(file, datafilename[i], 1);//comment line

			//Setting the finders for this data set
			fg_gfinder[i] = fg_gvalues + fgused_tot;//setting the finder for g-values
			fg_ffinder[i] = fg_fvalues + fgused_tot;//setting the finder for F-values
			for (ipartial = 0; ipartial < npartials; ipartial++)//fgused_tot here only contains the used g points of the previous data sets
				fg_cfinder[i*npartials + ipartial] = fg_coeffs + fgused_tot * npartials + fgused[i] * ipartial;//setting the finder for the g dependent coefficients

			for (j = 1; j < fgmin[i]; j++)
				SkipLine(file, datafilename[i], 1);//skipping all the unused g points

			//The values will be stored the following way in the fg_fvalues array:
			//first data set:F(g1)-F(gmax),second data set:F(g1)-F(gmax)...
			for (j = 0; j < fgused[i]; j++)//for all used points
			{
				IntToStr(&conv_numb, j + 1);
				mystrcpy(base_name, NAME_SIZE, "F(g) series used g point ");
				mystrcpy(base_name2, NAME_SIZE, "F(g) used point ");
				mystrcpy(base_name3, NAME_SIZE, "fgcoeff used point ");
				mystrcat(base_name, NAME_SIZE, conv_numb);
				mystrcat(base_name2, NAME_SIZE, conv_numb);
				mystrcat(base_name3, NAME_SIZE, conv_numb);
				mystrcat(base_name3, NAME_SIZE, " partial ");
				*pg = ReadThisLine(file, 3, 1.0, base_name, "ExptsData::ExptsData");//loads the g value
				pg++;
				*pf = ReadThisLine(file, 3, 1.0, base_name2, "ExptsData::ExptsData");//loads the F(g) value for the i-th expt, j-th g value
				*pf -= fgsub[i];//substracting the constant
				pf++;
				for (ipartial = 0; ipartial < npartials; ipartial++)
				{
					mystrcpy(base_name4, NAME_SIZE,base_name3);
					IntToStr(&conv_numb, ipartial + 1);
					mystrcat(base_name4, NAME_SIZE, conv_numb);

					if (ipartial < npartials - 1)
						*(fg_cfinder[i*npartials + ipartial] + j) = ReadThisLine(file, 3, 1.0, base_name4, "ExptsData::ExptsData");//reading the coefficients
					else
					{
						if (j < fgused[i] - 1)//give warning, if end of file was reached during line skipping after data reading
							*(fg_cfinder[i*npartials + ipartial] + j) = ReadThisLine(file, 1, 1.0, base_name4, "ExptsData::ExptsData", datafilename[i]);
						else//do not give warning, if end of file was reached during line skipping after data reading, as this can be the last line
							*(fg_cfinder[i*npartials + ipartial] + j) = ReadThisLine(file, 1, 1.0, base_name4, "ExptsData::ExptsData");
					}
				}
			}
			if (!CheckReadFileState(file, "ExptsData constructor", datafilename[ngr + nsq + nfq + i]))
				CleanExit();
			next_line_pos = -1;
			file.close();

			fgused_tot += fgused[i];//increment the data points counter

		}//next F(g) constraint 
		//setting the cumulative number of used F(g) points in the previous data sets for each set
		fgused_cum[0] = 0;
		for (i = 1; i <= nfg; i++)
			fgused_cum[i] = fgused_cum[i - 1] + fgused[i - 1];

	}//end of if there are F(g) constraints)

	//now E(k) constraints
	//the number of r points (colums) is not given in the coefficient file, and it will not be cheked!
	//Coeff(k,r) for kmin <= k <= kmax and rmin <= r <= rmax are read from the file
	if(nek>0)//if there are some E(k) constraints
	{
		//Allocting memory space to arrays
		//here only the *_ori arrays will be created, these will be used to read in the data
		//in SetE0ShiftParams the values will be copied to the respective arrays which will be used later on
		//in the calculation. In case of E0 shift the number of data points will change due to the interpolation, and only the data points still
		//available will be copied to the respective arrays for later use. The ekused , ekused_tot and ekused_cum will be reset as well
		//the *_ori arrays will be deleted
		SetArraysize(&ek_kvalues_ori,ekused_tot,"ek_kvalues_ori","ExptsData::ExptsData");
		SetArraysize(&ek_evalues_ori,ekused_tot,"ek_evalues_ori","ExptsData::ExptsData");
		SetArraysize(&ek_coeffs_ori,ekused_tot*ek_rused_tot*ntypes,"ek_coeffs_ori","ExptsData::ExptsData");
		SetArraysize(&ek_kfinder_ori,nek,"ek_kfinder_ori","ExptsData::ExptsData");
		SetArraysize(&ek_efinder_ori,nek,"ek_efinder_ori","ExptsData::ExptsData");
		SetArraysize(&ek_cfinder_ori,nek*ntypes,"ek_cfinder_ori","ExptsData::ExptsData");
		
		//reset to zero,to calculate it again in case the expt dat file has less k point in than the .dat file
		ekused_tot=0;
		ecoeff_tot=0;
	
		pk=ek_kvalues_ori;//pointer for k values
		pe=ek_evalues_ori;//pointer for E values
		pc=ek_coeffs_ori;//pointer for coefficients
		//loading the data	
		for(i=0;i<nek;i++)//for each E(k) constraints
		{
			cout<<"\nLoading data from E(k) file: " << datafilename[ngr + nsq + nfq + nfg + i * 2] << endl;
			SafeOpenTextFile(file, datafilename[ngr + nsq + nfq + nfg + i * 2]);
			if (CheckFileState(file, "ExptsData::ExptsData", datafilename[ngr + nsq + nfq + nfg + i * 2]) == 0)
			{
				cout<<"Cannot load data, exiting..."<<endl;
				CleanExit();
			};
					
			ntot=ReadThisLine(file,1,1,"E(k) ntot", "ExptsData::ExptsData",datafilename[i]);//total number of k points in the data file
			//if the maximum number of k points to be used from the .dat file is greater than the number of k points in the file
			if (ntot<ekmax[i])
			{
				cout << "\nWARNING(" << ++warn << "): The largest index of the k data points to be used exceeds the total number of k points in the file! ";
				cout<<"\tMaximal number of k points to be used is "<<ekmax[i]<<" from the .dat file, "<<ntot<<\
					" from the experimental data file "<<endl;
				cout<<"\tSetting the maximal number of k points to be used to "<<ntot<<endl;
				ekmax[i]=ntot;
			}
			ekused[i]=ekmax[i]-ekmin[i]+1;//number of file points to be used
			eksize[i]=ntot;

			SkipLine(file,datafilename[i],1);//comment line
					
			//Setting the finders for this data set
			ek_kfinder_ori[i]=ek_kvalues_ori+ekused_tot;//setting the finder for k-values
			ek_efinder_ori[i]=ek_evalues_ori+ekused_tot;//setting the finder for E-values
			//setting the finder for the (k,r) dependent coefficients for each atoms type for each EXAFS data series
			j=0;
			for (itype=0;itype<ntypes;itype++)
			{
				ek_cfinder_ori[i*ntypes+itype]=ek_coeffs_ori+ecoeff_tot;//setting the finder for the k dependent coefficients
				ecoeff_tot+=ekused[i]*ek_rused[i*ntypes+itype];//number of the sofar used coeffs
			}

			for(j=1;j<ekmin[i];j++)
				SkipLine(file,datafilename[i],1);//skipping all the unused k points

			//The values will be stored the following way in the ek_evalues array:
			//first data set:E(k1)-E(kmax),second data set:E(k1)-E(kmax)...
			for(j=0;j<ekused[i];j++)//for all used points
			{
				IntToStr(&conv_numb, j + 1);
				mystrcpy(base_name, NAME_SIZE, "E(k)) series used k point ");
				mystrcpy(base_name2, NAME_SIZE, "E(k) used point ");
				mystrcat(base_name, NAME_SIZE, conv_numb);
				mystrcat(base_name2, NAME_SIZE, conv_numb);
				*pk=ReadThisLine(file,3,1.0,base_name, "ExptsData::ExptsData");//loads the k value
			
				if (j<ekused[i]-1)//give warning, if end of file was reached during line skipping after data reading
					*pe=ReadThisLine(file,1,1.0,base_name2, "ExptsData::ExptsData",datafilename[i]);
				else//do not give warning, if end of file was reached during line skipping after data reading, as this can be the last line
					*pe=ReadThisLine(file,1,1.0,base_name2, "ExptsData::ExptsData");
			
				*pe *=pow(*pk,(double)ekchipower[i]);//multiplying with k^ekchipower[i]
				pk++;
				pe++;
			}
		
			if (!CheckReadFileState(file,"ExptsData constructor",datafilename[ngr+nsq+nfq+nfg+i*2]))
				CleanExit();
			next_line_pos = -1;
			file.close();

						
			ekused_tot+=ekused[i];//increment the data points counter
			
			//Checking, whether there are the proper number of k points in the coefficient file
			cout << "\nLoading the coefficients from file: " << datafilename[ngr + nsq + nfq + nfg + i * 2 + 1] << endl;
			SafeOpenTextFile(file, datafilename[ngr + nsq + nfq + nfg + i * 2 + 1]);
			if (CheckFileState(file, "ExptsData::ExptsData", datafilename[ngr + nsq + nfq + nfg + i * 2 + 1]) == 0)
			{
				cout<<"Cannot load data, exiting..."<<endl;
				CleanExit();
			};
			etype=ReadThisLine(file,3,1,"number of atom types in the coeff file", "ExptsData::ExptsData");//number of atom types in the coefficient file		
			file.ignore(50,',');//extracting the comma
			ntot=ReadThisLine(file,1,1,"E(k) ntot coeff file", "ExptsData::ExptsData",datafilename[ngr+nsq+nfq+nfg+i*2+1]);//total number of k points in the coefficient file
						
			//Checking the consistency of the data and coefficient files
			if (etype!=ntypes)
			{
				cout<<"\n*****ERROR****"<<endl;
				cout<<"The number of atom types for the EXAFS edge is "<<etype<<" from the coefficient file, "<<endl;
				cout<<"and "<<ntypes<<" from the configuration file!"<<endl;
				cout<<"Cannot run this way exiting..."<<endl;
				CleanExit();
			}
			if (ntot<ekmax[i])
			{
				cout<<"\n*****ERROR*****"<<endl;
				cout<<"There are not enough k points in the coefficient file!"<<endl;
				cout<<"The last data point to use is "<<ekmax[i]<<" according to the *.dat and E(k) data file, but the"<<endl;
				cout<<"number of k points in the coefficient file is "<<ntot<<endl;
				cout<<"Cannot run this way exiting..."<<endl;
				CleanExit();
			}
			else
			{
				if (ntot!=eksize[i])
				{
					cout << "\nWARNING(" << ++warn << "): The total number of k points in the data file is "<<eksize[i]<<" and in the coefficient file is "<<ntot<<"!"<<endl;
					cout<<"\tTrying to proceed!"<<endl;
				}
			}
			//reading the coefficients
			//the order of the neighbour atom type always follows the order in the .dat and .cfg file 
			//first check, whether there is enough histogram bins in the histogram
#ifdef _VIBR_AMP
			mybins=HistoSet::nbins_used;
#else
			mybins=RunParams::nbins[assign_hist[ngr+nsq+nfq+nfg+i]];//the bin size index for this data set
#endif		
			
			for (itype=0;itype<ntypes;itype++)
			{
				if (mybins<ek_rmax[i*ntypes+itype])
				{
					cout<<"\n*****ERROR*****"<<endl;
					cout<<"The number of available histogram bins is "<<mybins<<", which is smaller, than the ";
					cout<<"index of the largest histogram bin to use ("<< ek_rmax[i*ntypes+itype]<<") given in the ";
					cout<<datfilename<<" for the "<<i+1<<". EXAFS data series!"<<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}
				//Read the number of r points (columns), if it is given
				string l;
				bool rv=static_cast<bool>(getline(file,l));
				stringstream ss(l);
				if(rv) rv=static_cast<bool>(ss>>j);
				if (!rv)//no number could be read, probably the number of colums was not specified
				{
					if (::debug)
					{
						cout << "NOTE(" << ++note << "): The number of r points was not specified for the " << i + 1 << ". EXAFS data set's " << itype + 1 << " partial," << endl;
						cout << "\tso we are attempting to determine it by the number of columns in the E(r,k) datafile!" << endl;
					}
					//file.clear(ios::goodbit);//set it back to goodbit, it is needed, if the file variable is used again				
					std::streampos cpos=file.tellg();
					getline(file,l);
					std::transform(l.begin(),l.end(),l.begin(), [](char ch) {if(ch==',') ch=' '; return ch;});
					file.seekg(cpos);
					stringstream ss2(l);
					j=0;
					double val;
					while(ss2>>val) j++;
					if (::debug)
						cout<<"\tNumber of r columns determined: "<<j<<endl;
				}
				if (j<ek_rmin[i*ntypes+itype])
				{
					cout<<"\n*****ERROR*****"<<endl;
					cout<<"The number of r columns for the "<<i+1<<". EXAFS data set's "<<itype+1<<" partial is smaller"<<endl;
					cout<<"according the the value given in the coefficient file ("<<j<<") than the first r points to be used ("<<ek_rmin[i*ntypes+itype]<<") !"<<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}
				if (j<ek_rmax[i*ntypes+itype])
				{
					cout<<"\n*****ERROR*****"<<endl;
					cout<<"The number of r colums for the "<<i+1<<". EXAFS data set's "<<itype+1<<" partial is smaller"<<endl;
					cout<<"according the the value given in the coefficient file ("<<j<<") than the last r points to be used ("<<ek_rmax[i*ntypes+itype]<<") !"<<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}
				//SkipLine(file,datafilename[ngr+nsq+nfq+nfg+i*2+1],1);//skip the comment belonging to the next type
				for(j=1;j<ekmin[i];j++)
				{
					SkipLine(file,datafilename[ngr+nsq+nfq+nfg+i*2+1],1);//skipping all the unused k points
				}
				for(j=ekmin[i];j<=ekmax[i];j++)//for each USED k data point
				{
					IntToStr(&conv_numb, itype + 1);
					mystrcpy(base_name3, NAME_SIZE, "ekcoeff type ");
					mystrcat(base_name3, NAME_SIZE, conv_numb);
					mystrcat(base_name3, NAME_SIZE, " k point ");
					IntToStr(&conv_numb, j);
					mystrcat(base_name3, NAME_SIZE, conv_numb);
					mystrcat(base_name3, NAME_SIZE, " r column ");
					for(k=1;k<ek_rmin[i*ntypes+itype];k++)
					{
						mystrcpy(base_name4, NAME_SIZE, base_name3);
						IntToStr(&conv_numb, k );
						mystrcat(base_name4, NAME_SIZE, conv_numb);
						*pc=ReadThisLine(file,3,1.0,base_name4, "ExptsData::ExptsData");//reading all the unused r points, but they are discarded
						if (k<ek_rmax[i*ntypes+itype])//to ensure not to try to extract the comma,if this is possibly the last r point, cannot happen theoreticaly
							file.ignore(50,',');//extracting the comma
					}
					for(k=ek_rmin[i*ntypes+itype];k<=ek_rmax[i*ntypes+itype];k++)//for each USED r point
					{
						mystrcpy(base_name4, NAME_SIZE, base_name3);
						IntToStr(&conv_numb, k);
						mystrcat(base_name4, NAME_SIZE, conv_numb);
						*pc=ReadThisLine(file,3,1.0,base_name4, "ExptsData::ExptsData");//reading all the used r points
						pc++;
						if (k<ek_rmax[i*ntypes+itype])//to ensure not to try to extract the comma,if this is possibly the last r point
							file.ignore(50,',');//extracting the comma
					}
					if (!(itype==ntypes-1 && j==ekmax[i]))
						SkipLine(file,datafilename[ngr+nsq+nfq+nfg+i*2+1],1);//goto next k point, if this not the last used k of the lasr partial
				}
				if (itype<ntypes-1)//to ensure not to try to go next line, if this is the last type
				{
					for(j=ekmax[i]+1;j<=ntot;j++)
						SkipLine(file,datafilename[ngr+nsq+nfq+nfg+i*2+1],1);//skipping all the unused k points

				}
			}
			if (!CheckReadFileState(file,"ExptsData constructor",datafilename[ngr+nsq+nfq+ nfg + i*2+1]))
				CleanExit();
			next_line_pos = -1;
			file.close();
				
		}//next E(k) constraint 
		//setting the cumulative number of used E(k) points in the previous data sets for each set
		ekused_cum[0]=0;
		for (i=1;i<=nek;i++)
			ekused_cum[i]=ekused_cum[i-1]+ekused[i-1];
		
	}//end of if there are E(k) constraints)
	
	if (conv_numb!=NULL)
		delete [] conv_numb;
	delete[] base_name;
	delete[] base_name2;
	delete[] base_name3;
	delete[] base_name4;
	delete[] base_name5;
	
//Checking, whether loading was successful
	if (file && ::debug) 
		cout<<"\nExperimental data loading was successful"<<endl;//load was successful
}//end of the constructor	


//Reading the main parameters for the experimental data	
void ExptsData::GetExptParam(ifstream &file)
{
	int i,j,index;
	int npartials;
	char *conv_numb1,*conv_numb2,*conv_numb3,*conv_numb4;
	string myline;
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
	double dumb;
#endif
#endif
	if (CheckFileState(file,"ExptsData::GetExptParam",datfilename)==0)
	{
		cout<<"ExptsData::GetExptParam: "<<datfilename<<" file is not open or non-existing."<<endl;
		CleanExit();
	};

	conv_numb1 = NULL;
	conv_numb2 = NULL;
	conv_numb3=new char[100];
	conv_numb4=new char[100];

	ntypes=RunParams::ntypes;
	npartials=ntypes*(ntypes+1)/2;	
	if (npartials<1)//Check whether it was properly initialised
	{
		cout<<"\n*****ERROR*****"<<endl;
		cout<<"ExptsData::GetExptsParams - The number of partials is smaller than 1! "<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	ngr=RunParams::ngr;
	nsq=RunParams::nsq;
	nfq=RunParams::nfq;
	nfg=RunParams::nfg;
	nek=RunParams::nek;
	ntot_datasets = RunParams::ntot_datasets;
	assign_hist = RunParams::assign_hist;
	rspacing = RunParams::rspacing;

	//first the g(r) then the S(Q),F(Q), F(g) finally the E(k) data files name are stored
	datafilename=new char[ngr+nsq+nfq+nfg+2*nek][FILE_NAME_SIZE];//there data and the coefficients are in separate files for EXAFS
	if (datafilename==NULL)
		NoArray("ExptsData::GetExptsParam","datafilename");

	//Arrays are allocated regardless if there is a constraint of that type
	SetArraysize(&grmin,ngr,"grmin","ExptsData::GetExptsParam");
	SetArraysize(&grmax,ngr,"grmax","ExptsData::GetExptsParam");
	SetArraysize(&grsize,ngr,"grsize","ExptsData::GetExptsParam");
	SetArraysize(&grsigma,ngr,"grsigma","ExptsData::GetExptsParam");
	SetArraysize(&grused,ngr,"grused","ExptsData::GetExptsParam");
	SetArraysize(&grused_cum,ngr+1,"grused_cum","ExptsData::GetExptsParam");
	SetArraysize(&grrenorm,ngr,"grrenorm","ExptsData::GetExptsParam");
	SetArraysize(&groffset,ngr,"groffset","ExptsData::GetExptsParam");
	SetArraysize(&grlinear,ngr,"grlinear","ExptsData::GetExptsParam");
	SetArraysize(&grquadratic,ngr,"grquadratic","ExptsData::GetExptsParam");
	SetArraysize(&grcubic,ngr,"grcubic","ExptsData::GetExptsParam");
	SetArraysize(&grsub,ngr,"grsub","ExptsData::GetExptsParam");
	SetArraysize(&gr_coeffs,ngr*npartials,"gr_coeffs","ExptsData::GetExptsParam");
	SetArraysize(&bin_flag, ngr, "bin_flag", "ExptsData::GetExptsParam");
	SetArraysize(&gruseR, ngr, "gruseR", "ExptsData::GetExptsParam");

	SetArraysize(&sqmin,nsq,"sqmin","ExptsData::GetExptsParam");
	SetArraysize(&sqmax,nsq,"sqmax","ExptsData::GetExptsParam");
	SetArraysize(&sqsize,nsq,"sqsize","ExptsData::GetExptsParam");
	SetArraysize(&sqused,nsq,"sqused","ExptsData::GetExptsParam");
	SetArraysize(&sqused_cum,nsq+1,"sqused_cum","ExptsData::GetExptsParam");
	SetArraysize(&sqrenorm,nsq,"sqrenorm","ExptsData::GetExptsParam");
	SetArraysize(&sqoffset,nsq,"sqoffset","ExptsData::GetExptsParam");
	SetArraysize(&sqlinear,nsq,"sqlinear","ExptsData::GetExptsParam");
	SetArraysize(&sqquadratic,nsq,"sqquadratic","ExptsData::GetExptsParam");
	SetArraysize(&sqcubic,nsq,"sqcubic","ExptsData::GetExptsParam");
	SetArraysize(&sqsub,nsq,"sqsub","ExptsData::GetExptsParam");
	SetArraysize(&sqsigma,nsq,"sqsigma","ExptsData::GetExptsParam");
	SetArraysize(&sq_coeffs,nsq*npartials,"sq_coeffs","ExptsData::GetExptsParam");
	SetArraysize(&squseR, nsq, "squseR", "ExptsData::GetExptsParam");
	SetArraysize(&sqreadcoeffs, nsq, "sqreadcoeffs", "ExptsData::GetExptsParam");
	SetArraysize(&isotope_count, nsq * ntypes, "isotope_count", "ExptsData::GetExptsParamFree");

	SetArraysize(&fqmin, nfq, "fqmin", "ExptsData::GetExptsParam");
	SetArraysize(&fqmax, nfq, "fqmax", "ExptsData::GetExptsParam");
	SetArraysize(&fqsize, nfq, "fqsize", "ExptsData::GetExptsParam");
	SetArraysize(&fqused, nfq, "fqused", "ExptsData::GetExptsParam");
	SetArraysize(&fqused_cum, nfq + 1, "fqused_cum", "ExptsData::GetExptsParam");
	SetArraysize(&fqrenorm, nfq, "fqrenorm", "ExptsData::GetExptsParam");
	SetArraysize(&fqoffset, nfq, "fqoffset", "ExptsData::GetExptsParam");
	SetArraysize(&fqlinear, nfq, "fqlinear", "ExptsData::GetExptsParam");
	SetArraysize(&fqquadratic, nfq, "fqquadratic", "ExptsData::GetExptsParam");
	SetArraysize(&fqcubic, nfq, "fqcubic", "ExptsData::GetExptsParam");
	SetArraysize(&fqsub, nfq, "fqsub", "ExptsData::GetExptsParam");
	SetArraysize(&fqsigma, nfq, "fqsigma", "ExptsData::GetExptsParam");
	SetArraysize(&fqreadcoeffs, nfq, "fqreadcoeffs", "ExptsData::GetExptsParam");
	SetArraysize(&fquseR, nfq, "fquseR", "ExptsData::GetExptsParam");
	SetArraysize(&fqfitIQ, nfq, "fqfitIQ", "ExptsData::GetExptsParam");
	SetArraysize(&fqIQbackgcorr, nfq, "fqIQbackgcorr", "ExptsData::GetExptsParam");
	SetArraysize(&fqfprimecount, nfq, "fqfprimecount", "ExptsData::GetExptsParam");
	SetArraysize(&fqfprime, nfq*ntypes, "fqfprime", "ExptsData::GetExptsParam");
	SetArraysize(&fqfprimeact, nfq*ntypes, "fqfprimeact", "ExptsData::GetExptsParam");
	SetArraysize(&fqfprimeindex, nfq, "fqfprimeindex", "ExptsData::GetExptsParam");
	SetArraysize(&fqfprimefactor, nfq, "fqfprimefactor", "ExptsData::GetExptsParam");
	SetArraysize(&fqAXS, nfq, "fqAXS", "ExptsData::GetExptsParam");
	SetArraysize(&fqmu, nfq, "fqmu", "ExptsData::GetExptsParam");
	SetArraysize(&fqdmumax, nfq, "fqdmumax", "ExptsData::GetExptsParam");
	SetArraysize(&fqmuact, nfq, "fqmuact", "ExptsData::GetExptsParam");
	SetArraysize(&fqrenalpha, nfq, "fqrenalpha", "ExptsData::GetExptsParam");
	SetArraysize(&fqusecompton, nfq, "fqusecompton", "ExptsData::GetExptsParam");
	SetArraysize(&fqalpha, nfq, "fqalpha", "ExptsData::GetExptsParam");
	SetArraysize(&fqa, nfq, "fqa", "ExptsData::GetExptsParam");
	SetArraysize(&fqb, nfq, "fqb", "ExptsData::GetExptsParam");
	
	SetArraysize(&fgmin, nfg, "fgmin", "ExptsData::GetExptsParam");
	SetArraysize(&fgmax, nfg, "fgmax", "ExptsData::GetExptsParam");
	SetArraysize(&fgsize, nfg, "fgsize", "ExptsData::GetExptsParam");
	SetArraysize(&fgused, nfg, "fgused", "ExptsData::GetExptsParam");
	SetArraysize(&fgused_cum, nfg + 1, "fgused_cum", "ExptsData::GetExptsParam");
	SetArraysize(&fgrenorm, nfg, "fgrenorm", "ExptsData::GetExptsParam");
	SetArraysize(&fgoffset, nfg, "fgoffset", "ExptsData::GetExptsParam");
	SetArraysize(&fglinear, nfg, "fglinear", "ExptsData::GetExptsParam");
	SetArraysize(&fgquadratic, nfg, "fgquadratic", "ExptsData::GetExptsParam");
	SetArraysize(&fgcubic, nfg, "fgcubic", "ExptsData::GetExptsParam");
	SetArraysize(&fgsub, nfg, "fgsub", "ExptsData::GetExptsParam");
	SetArraysize(&fgsigma, nfg, "fgsigma", "ExptsData::GetExptsParam");
	SetArraysize(&fguseR, nfg, "fguseR", "ExptsData::GetExptsParam");

	SetArraysize(&ek_rmin,ntypes*nek,"ek_rmin","ExptsData::GetExptsParam");
	SetArraysize(&ek_rmax,ntypes*nek,"ek_rmax","ExptsData::GetExptsParam");
	SetArraysize(&ek_rused,ntypes*nek,"ek_rused","ExptsData::GetExptsParam");
	SetArraysize(&ekmin,nek,"ekmin","ExptsData::GetExptsParam");
	SetArraysize(&ekmax,nek,"ekmax","ExptsData::GetExptsParam");
	SetArraysize(&ekmin_ori, nek, "ekmin_ori", "ExptsData::GetExptsParam");
	SetArraysize(&ekmax_ori, nek, "ekmax_ori", "ExptsData::GetExptsParam");
	SetArraysize(&eksize,nek,"eksize","ExptsData::GetExptsParam");
	SetArraysize(&ekused,nek,"ekused","ExptsData::GetExptsParam");
	SetArraysize(&ekused_cum,nek+1,"ekused_cum","ExptsData::GetExptsParam");
	SetArraysize(&ekrenorm,nek,"ekrenorm","ExptsData::GetExptsParam");
	SetArraysize(&ekoffset,nek,"ekoffset","ExptsData::GetExptsParam");
	SetArraysize(&ek_abstype,nek,"ek_abstype","ExptsData::GetExptsParam");
	SetArraysize(&ekchipower,nek,"ekchipower","ExptsData::GetExptsParam");
	SetArraysize(&eksigma,nek,"eksigma","ExptsData::GetExptsParam");
	SetArraysize(&ekuseR, nek, "ekuseR", "ExptsData::GetExptsParam");
	SetArraysize(&ek_dE0, nek, "ek_dE0", "ExptsData::GetExptsParamFree");//this will not be read by fixed format, only default values are set
	SetArraysize(&ek_ngrid, nek, "ek_ngrid", "ExptsData::GetExptsParam");//this will not be read by fixed format, only default values are set
	SetArraysize(&ek_ngrid_in, nek, "ek_ngrid_in", "ExptsData::GetExptsParam");//this will not be read by fixed format, only default values are set
	SetArraysize(&ek_gridstart, nek, "ek_gridstart", "ExptsData::GetExptsParam");//this will not be read by fixed format, only default values are set
	SetArraysize(&ek_gridfrom, nek, "ek_gridfrom", "ExptsData::GetExptsParam");//this will not be read by fixed format, only default values are set
	SetArraysize(&ek_gridto, nek, "ek_gridto", "ExptsData::GetExptsParam");//this will not be read by fixed format, only default values are set
	SetArraysize(&ek_gridind, nek, "ek_gridind", "ExptsData::GetExptsParam");//index of the presently selected gridpoint
	SetArraysize(&ek_cf_cum, nek+1, "ek_cf_cum", "ExptsData::GetExptsParam");//cumulative number of coeff matrices for each data set, set to default
	SetArraysize(&ek_ngrid_cum, nek + 1, "ek_ngrid_cum", "ExptsData::GetExptsParam");//cumulative number of grid points for each data set, set to default
	//the indicator array for each g(r), F(Q), F(g) and S(Q)  data series whether to use cubic correction, 
	
	SetArraysize(&use_cubic,ngr+nsq+nfq+nfg,"ExptsData::use_cubic","ExptsData::GetExptParams");

	for (int i = 0; i < nsq * ntypes; i++)
	{
		isotope_count[i].count = 0;//default, natural abundance
		isotope_count[i].isotopes = nullptr;
	}
	//There are g(r) constraints
	if (ngr>0)
	{
		for (i=0;i<ngr;i++)
		{
			bin_flag[i] = 0;//default
			IntToStr(&conv_numb1,i+1);
			file>>datafilename[i];//file name
			SkipLine(file,datfilename,1);
			grmin[i]=ReadThisLine(file,3,1,"grmin", "ExptsData::GetExptParams");//first point to use (STARTING at 1)
			grmax[i]=ReadThisLine(file,3,1,"grmax", "ExptsData::GetExptParams",datfilename);//last point to use
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
			dumb= ReadThisLine(file, -8, -1.0, "rspacing", "ExptsData::GetExptParams", datfilename);//last point to use
			
			if (fabs(-1.0 - dumb)>LOAD_TOL)//there was a value
				rspacing[i] = dumb;//there is a custom bin size for this data set
#endif
#endif
			SkipLine(file);
			grsub[i]=ReadThisLine(file,1, SUBTRACT_DEF,"grsub", "ExptsData::GetExptParams",datfilename);//constant to be substracted
			for (j=0;j<npartials;j++)
			{
				IntToStr(&conv_numb2,j+1);
				mystrcpy(conv_numb3,100,"");
				mystrcat(conv_numb3, 100, conv_numb1);
				mystrcat(conv_numb3, 100, ". g(r) data set gr_coeff for partial ");
				mystrcat(conv_numb3, 100, conv_numb2);
				if (j<npartials-1)
					gr_coeffs[i*npartials+j]=ReadThisLine(file,3,1.0,conv_numb3, "ExptsData::GetExptParams");
				else
					gr_coeffs[i*npartials+j]=ReadThisLine(file,1,1.0,conv_numb3, "ExptsData::GetExptParams",datfilename);
			};
			//we have to make sure, that there is data in this line, as if no data and nothing else, it will read faultily from next line
			getline(file, myline);
			stringstream ss(myline);
			ss >> grsigma[i];
			if (ss.fail())
			{
				cout << "\n*****ERROR*****" << endl;
				cout<<"Variable grsigma[" << i + 1 << "] could not be read in ExptsData::GetExptParams, probably end of line was reached!" << endl;
				CleanExit();
			}
			RunParams::ReadLogicalS(ss, gruseR[i], "gruseR", "ExptsData::GetExptParams");//go to next line
			if (grsigma[i]<0)
				ChiSquared::calc_sigma=1;
			//this is the renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,grrenorm[i],"grrenorm", "ExptsData::GetExptParams",3);//not to go to next line
			use_cubic[i] = ReadThisLine(file, 6, USE_CUBIC_DEF, "use_cubic", "ExptsData::GetExptParams", datfilename);//last point to use
			if (use_cubic[i])
			{
				//this is the offset switch (cannot be read directly, as it can be in the old format .false. or .true.,
				//but it can also be 0 for false and 1 for true)
				RunParams::ReadLogical(file,groffset[i],"groffset","ExptsData::GetExptParams");
				//this is the linear renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
				//but it can also be 0 for false and 1 for true)
				RunParams::ReadLogical(file,grlinear[i],"grlinear","ExptsData::GetExptParams");
				//this is the quadratic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
				//but it can also be 0 for false and 1 for true)
				RunParams::ReadLogical(file,grquadratic[i],"grquadratic", "ExptsData::GetExptParams");
				//this is the cubic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
				//but it can also be 0 for false and 1 for true)
				RunParams::ReadLogical(file,grcubic[i],"grcubic", "ExptsData::GetExptParams");
			}
			else
			{
				groffset[i]=RENORM_DEF;
				grlinear[i]= RENORM_DEF;
				grquadratic[i]= RENORM_DEF;
				grcubic[i]= RENORM_DEF;
			}

			//Calculating the number of used g(r) points based on the values given in the .dat file
		
			if (grmax[i]==0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout<<"Maximal number of g(r) points is 0 in the "<<i+1<<". data set!"<<endl;
				cout<<"Cannot run this way, exiting..."<<endl;
				CleanExit();
			}
			rused_tot+=grmax[i]-grmin[i]+1;
		}
	}//end of if there are g(r) constraints

	//If there is S(Q) constraint
	if (nsq>0)
	{
		for (i=0;i<nsq;i++)
		{
			sqreadcoeffs[i] = 1;
			index=ngr+i;
			IntToStr(&conv_numb1,i+1);
			file>>datafilename[ngr+i];//file name 
			SkipLine(file,datfilename,1);
			sqmin[i]=ReadThisLine(file,3,1,"sqmin", "ExptsData::GetExptParams");//first Q point to use (STARTING at 1)
			sqmax[i]=ReadThisLine(file,3,1,"sqmax", "ExptsData::GetExptParams",datfilename);//last Q point to use
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
			dumb = ReadThisLine(file, -8, -1.0, "rspacing", "ExptsData::GetExptParams", datfilename);//last point to use

			if (fabs(-1.0 - dumb) > LOAD_TOL)//there was a value
				rspacing[index] = dumb;//there is a custom bin size for thos data set
#endif
#endif
			SkipLine(file);
			sqsub[i]=ReadThisLine(file,1,SUBTRACT_DEF,"sqsub", "ExptsData::GetExptParams",datfilename);//constant to be substracted
			for (j=0;j<npartials;j++)
			{
				IntToStr(&conv_numb2,j+1);
				mystrcpy(conv_numb3, 100, "");
				mystrcat(conv_numb3, 100, conv_numb1);
				mystrcat(conv_numb3, 100, ". S(Q) data set sq_coeff for partial ");
				mystrcat(conv_numb3, 100, conv_numb2);
				if (j<npartials-1)
					sq_coeffs[i*npartials+j]=ReadThisLine(file,3,1.0,conv_numb3, "ExptsData::GetExptParams");
				else
					sq_coeffs[i*npartials+j]=ReadThisLine(file,1,1.0,conv_numb3, "ExptsData::GetExptParams",datfilename);
			};
				
			//we have to make sure, that there is data in this line, as if no data and nothing else, it will read faultily from next line
			getline(file,myline);
			stringstream ss(myline);
			ss >> sqsigma[i];
			if (ss.fail())
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Variable sqsigma[" <<i+1 << "] could not be read in ExptsData::GetExptParams, probably end of line was reached!" << endl;
				CleanExit();
			}
			RunParams::ReadLogicalS(ss, squseR[i], "squseR", "ExptsData::GetExptParams");//go to next line
		
			if (sqsigma[i]<0)
				ChiSquared::calc_sigma=1;
			//this is the renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,sqrenorm[i],"sqrenorm", "ExptsData::GetExptParams",3);
			use_cubic[index] = ReadThisLine(file, 6, USE_CUBIC_DEF, "use_cubic", "ExptsData::GetExptParams", datfilename);//whether to use cubic
			//this is the offset switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,sqoffset[i],"sqoffset", "ExptsData::GetExptParams");
			//this is the linear renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,sqlinear[i],"sqlinear", "ExptsData::GetExptParams");
			//this is the quadratic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,sqquadratic[i],"sqquadratic","ExptsData::GetExptParams");
			//this is the cubic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			if (use_cubic[index])
				RunParams::ReadLogical(file,sqcubic[i],"sqcubic", "ExptsData::GetExptParams");
			else
				sqcubic[i]=RENORM_DEF;
		
			if (sqmax[i]==0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout<<"Maximal number of Q points is 0 in the "<<i+1<<". S(Q) data set!"<<endl;
				cout<<"Cannot run this way, exiting..."<<endl;
				CleanExit();
			}
			sqused_tot+=(sqmax[i]-sqmin[i]+1);
		}
	}//end of if there are S(Q) constraints

	//If there is F(Q) constraint
	if (nfq>0)
	{
		for (int i = 0; i < nfq * ntypes; i++)
		{
			//fqfprimeact will be added during cieffs calculation regardless it's value, has to be set to 0
			fqfprime[i] = 0.0;
			fqfprimeact[i] = 0.0;
		}
		for (i=0;i<nfq;i++)
		{
			index=ngr+nsq+i;
			fqfitIQ[i] = false;//fixed format cannot read the I(Q) fit realted data
			fqIQbackgcorr[i] = false;
			fqfprimecount[i] = 0;
			fqusecompton[i] = false;
			fqreadcoeffs[i] = 1;
			fqAXS[i] = 0;
			fqfprimeindex[i] = 0;
			fqmuact[i] = 0.0;
			fqrenalpha[i] = 0;
			fqalpha[i] = 0.0;
			fqa[i] = 1.0;
			fqb[i] = 0.0;
			file>>datafilename[ngr+nsq+i];//file name 
			SkipLine(file,datfilename,1);
			fqmin[i]=ReadThisLine(file,3,1,"fqmin", "ExptsData::GetExptParams");//first Q point to use (STARTING at 1)
			fqmax[i]=ReadThisLine(file,3,1,"fqmax", "ExptsData::GetExptParams",datfilename);//last Q point to use
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
			dumb = ReadThisLine(file, -8, -1.0, "rspacing", "ExptsData::GetExptParams", datfilename);//last point to use

			if (fabs(-1.0 - dumb) > LOAD_TOL)//there was a value
				rspacing[index] = dumb;//there is a custom bin size for thos data set
#endif
#endif
			SkipLine(file);
			fqsub[i]=ReadThisLine(file,1,SUBTRACT_DEF,"fqsub", "ExptsData::GetExptParams",datfilename);//constant to be substracted
			//we have to make sure, that there is data in this line, as if no data and nothing else, it will read faultily from next line
			getline(file, myline);
			stringstream ss(myline);
			ss >> fqsigma[i];
			if (ss.fail())
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Variable fqsigma[" << i + 1 << "] could not be read in ExptsData::GetExptParams, probably end of line was reached!" << endl;
				CleanExit();
			}
			RunParams::ReadLogicalS(ss, fquseR[i], "fquseR", "ExptsData::GetExptParams");//go to next line
		
			if (fqsigma[i]<0)
				ChiSquared::calc_sigma=1;
			//this is the renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,fqrenorm[i],"fqrenorm", "ExptsData::GetExptParams",3);
			use_cubic[index] = ReadThisLine(file, 6, USE_CUBIC_DEF, "use_cubic", "ExptsData::GetExptParams", datfilename);//last point to use
			//this is the offset switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,fqoffset[i],"fqoffset", "ExptsData::GetExptParams");
				//this is the linear renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,fqlinear[i],"fqlinear", "ExptsData::GetExptParams");
			//this is the quadratic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,fqquadratic[i],"fqquadratic","ExptsData::GetExptParams");
			//this is the cubic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			if (use_cubic[index])
				RunParams::ReadLogical(file,fqcubic[i],"fqcubic", "ExptsData::GetExptParams");
			else
				fqcubic[i]=RENORM_DEF;

			if (fqmax[i]==0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout<<"Maximal number of Q points is 0 in the "<<i+1<<". F(Q) data set!"<<endl;
				cout<<"Cannot run this way, exiting..."<<endl;
				CleanExit();
			}
			fqused_tot+=(fqmax[i]-fqmin[i]+1);
			
		}
	}//end of if there are F(Q) constraints

	//If there is F(g) constraint
	if (nfg > 0)
	{
		for (i = 0; i < nfg; i++)
		{
			index = ngr + nsq + nfq + i;
			file >> datafilename[ngr + nsq + nfq + i];//file name 
			SkipLine(file, datfilename, 1);
			fgmin[i] = ReadThisLine(file, 3, 1, "fgmin", "ExptsData::GetExptParams");//first g point to use (STARTING at 1)
			fgmax[i] = ReadThisLine(file, 3, 1, "fgmax", "ExptsData::GetExptParams", datfilename);//last g point to use
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
			dumb = ReadThisLine(file, -8, -1.0, "rspacing", "ExptsData::GetExptParams", datfilename);//last point to use

			if (fabs(-1.0 - dumb) > LOAD_TOL)//there was a value
				rspacing[index] = dumb;//there is a custom bin size for thos data set
#endif
#endif
			SkipLine(file);
			fgsub[i] = ReadThisLine(file, 1, SUBTRACT_DEF, "fgsub", "ExptsData::GetExptParams", datfilename);//constant to be substracted
			//we have to make sure, that there is data in this line, as if no data and nothing else, it will read faultily from next line
			getline(file, myline);
			stringstream ss(myline);
			ss >> fgsigma[i];
			if (ss.fail())
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Variable fgsigma[" << i + 1 << "] could not be read in ExptsData::GetExptParams, probably end of line was reached!" << endl;
				CleanExit();
			}
			RunParams::ReadLogicalS(ss, fguseR[i], "fguseR", "ExptsData::GetExptParams");//go to next line
	
			if (fgsigma[i] < 0)
				ChiSquared::calc_sigma = 1;
			//this is the renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file, fgrenorm[i], "fgrenorm", "ExptsData::GetExptParams",3);
			use_cubic[index] = ReadThisLine(file, 6, USE_CUBIC_DEF, "use_cubic", "ExptsData::GetExptParams", datfilename);//last point to use
			//this is the offset switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file, fgoffset[i], "fgoffset", "ExptsData::GetExptParams");
			//this is the linear renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file, fglinear[i], "fglinear", "ExptsData::GetExptParams");
			//this is the quadratic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file, fgquadratic[i], "fgquadratic", "ExptsData::GetExptParams");
			//this is the cubic renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			if (use_cubic[index])
				RunParams::ReadLogical(file, fgcubic[i], "fgcubic", "ExptsData::GetExptParams");
			else
				fgcubic[i] = RENORM_DEF;

			if (fgmax[i] == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Maximal number of g points is 0 in the " << i + 1 << ". F(g) data set!" << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}
			fgused_tot += (fgmax[i] - fgmin[i] + 1);
		}
	}//end of if there are F(g) constraints

	//If there is E(k) constraint
	if (nek>0)
	{
		ek_cf_cum[0] = 0;
		ek_ngrid_cum[0] = 0;
		for (i=0;i<nek;i++)
		{
			index = ngr + nsq + nfq + nfg+i;
			IntToStr(&conv_numb1,i+1);
			//each partial will have its own range, will be stored in the following order: 
			//first exafs series 1, ... partials, second series 1, ... partials, ...
			for (j=0;j<ntypes;j++)
			{

				IntToStr(&conv_numb2,j+1);
				mystrcpy(conv_numb3, 100, "");
				mystrcpy(conv_numb4, 100, "");
				mystrcat(conv_numb3, 100, conv_numb1);
				mystrcat(conv_numb4, 100, conv_numb1);
				mystrcat(conv_numb3, 100, ". E(k) data set ek_rmin for type ");
				mystrcat(conv_numb4, 100, ". E(k) data set ek_rmax for type ");
				mystrcat(conv_numb3, 100, conv_numb2);
				mystrcat(conv_numb4, 100, conv_numb2);
				
				if (j==0)//the values has to be given for at least the first partial
				{
				
					//first r point to (histogram bin) use
					ek_rmin[i*ntypes+j]=ReadThisLine(file,3,1,conv_numb3, "ExptsData::GetExptParams");
					//last r point (histogram bin) to use
					if (j<ntypes-1)
						//last r point (histogram bin) to use
						ek_rmax[i*ntypes+j]=ReadThisLine(file,3,1,conv_numb4, "ExptsData::GetExptParams",datfilename);
					else//last partial, has to go to next line after max reading
						ek_rmax[i*ntypes+j]=ReadThisLine(file,1,1,conv_numb4, "ExptsData::GetExptParams");
				}
				else
				{
					//first r point (histogram bin) to use
					ek_rmin[i*ntypes+j]=ReadThisLine(file,2,ek_rmin[i*ntypes],conv_numb3, "ExptsData::GetExptParams",datfilename);
					if (j<ntypes-1)
						//last r point (histogram bin) to use
						ek_rmax[i*ntypes+j]=ReadThisLine(file,2,ek_rmax[i*ntypes],conv_numb4, "ExptsData::GetExptParams",datfilename);
					else//last partial, has to go to next line after max reading
						ek_rmax[i*ntypes+j]=ReadThisLine(file,0,ek_rmax[i*ntypes],conv_numb4, "ExptsData::GetExptParams",datfilename);
				}
			}
			ekchipower[i]=ReadThisLine(file,1,EXAFS_CHI2_POWER_DEF,"ekchipower", "ExptsData::GetExptParams",datfilename);//(E(k)*(k power of ekchipower)) for the experimental data points
			file>>datafilename[ngr+nsq+nfq+nfg+(i*2)];//file name for E(k) data 
			SkipLine(file,datfilename,1);
			SkipLine(file,datfilename,1);//skipping the number of data points
			ekmin[i]=ReadThisLine(file,3,1,"ekmin", "ExptsData::GetExptParams");//first k point to use (STARTING at 1)
			ekmax[i]=ReadThisLine(file,3,1,"ekmax", "ExptsData::GetExptParams",datfilename);//last k point to use
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
			dumb = ReadThisLine(file, -8, -1.0, "rspacing", "ExptsData::GetExptParams", datfilename);//last point to use

			if (fabs(-1.0 - dumb) > LOAD_TOL)//there was a value
				rspacing[index] = dumb;//there is a custom bin size for thos data set
#endif
#endif
			SkipLine(file);
			ek_abstype[i]=ReadThisLine(file,1,1,"ek_abstype", "ExptsData::GetExptParams",datfilename);//indexing start with 1!
			ek_abstype[i]--;//to start with zero
			file>>datafilename[ngr+nsq+nfq+ nfg + (i*2)+1];//file name for the (k,r)-dependent coefficients 
			SkipLine(file,datfilename,1);
			//we have to make sure, that there is data in this line, as if no data and nothing else, it will read faultily from next line
			getline(file, myline);
			stringstream ss(myline);
			ss >> eksigma[i];
			if (ss.fail())
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Variable eksigma[" << i + 1 << "] could not be read in ExptsData::GetExptParams, probably end of line was reached!" << endl;
				CleanExit();
			}
			RunParams::ReadLogicalS(ss, ekuseR[i], "ekuseR", "ExptsData::GetExptParams");//go to next line
	
			if (eksigma[i]<0)
				ChiSquared::calc_sigma=1;
			//this is the renormalization switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,ekrenorm[i],"ekrenorm", "ExptsData::GetExptParams");
			//this is the offset switch (cannot be read directly, as it can be in the old format .false. or .true.,
			//but it can also be 0 for false and 1 for true)
			RunParams::ReadLogical(file,ekoffset[i],"ekoffset", "ExptsData::GetExptParams");
			if (ekmax[i]==0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout<<"Maximal number of k points is 0 in the "<<i+1<<". E(k) data set!"<<endl;
				cout<<"Cannot run this way, exiting..."<<endl;
				CleanExit();
			}
			ekused_tot+=(ekmax[i]-ekmin[i]+1);
			for (j=0;j<ntypes;j++)
			{
				if (ek_rmax[i*ntypes+1]==0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"Maximal number of r points is 0 in the "<<i+1<<". E(k) data set for the "<<j+1<<". partial!"<<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}
			
				ek_rused[i*ntypes+j]=(ek_rmax[i*ntypes+j]-ek_rmin[i*ntypes+j]+1);
				ek_rused_tot+=ek_rused[i*ntypes+j];
			}
			ek_dE0[i] = MAX_E0_SHIFT_DEF;
			ek_ngrid_in[i] = EXAFS_NGRID_DEF;
			ek_ngrid[i] =2* EXAFS_NGRID_DEF+1;
			ek_gridstart[i] = EXAFS_NGRID_DEF;
			ek_gridfrom[i] = EXAFS_NGRID_DEF;
			ek_gridto[i] = EXAFS_NGRID_DEF;
			ek_cf_cum[i + 1] = ek_cf_cum[i] + ntypes;
			ek_ngrid_cum[i + 1] = ek_ngrid_cum[i] + 1;
			
		}
	}//end of if there are E(k) constraints

	//Checking, whether loading was successful
	if (!CheckReadFileState(file,"AvCoordConst::GetAvCoordConst",datfilename))
			CleanExit();//loading failed
	if (conv_numb1 != NULL)
		delete [] conv_numb1;
	if (conv_numb2 != NULL)
		delete [] conv_numb2;
	if (conv_numb3 != NULL)
		delete [] conv_numb3;
	if (conv_numb4 != NULL)
		delete [] conv_numb4;
};	
		
//==========================save==============
void ExptsData::Save() const
{
	//saves the contents in a file
	
	int i,j,k,ig,index;
	double *px,*py, *pz,*pA,*pB,*pIB=NULL,*pIQ=NULL;
	int npartials,itype;
	ofstream file;
	mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
	mystrcat(tempfilename, FILE_NAME_SIZE + 10, ".expt");
	OpenFile(file,tempfilename,"ExptsData::Save",0);//open file, check, whether it was successfully opened

	npartials=ntypes*(ntypes+1)/2;
	
	file << "This is an object of class ExptsData "<<endl;
	file << "number of g(r) data sets ngr->"<<ngr<<endl;
	file << "number of S(Q) data sets nsq->"<<nsq<<endl;
	file << "number of F(Q) data sets nfq->"<<nfq<<endl;
	file << "number of F(g) data sets nfg->" << nfg << endl;
	file << "number of E(k) data sets nek->"<<nek<<endl;
	
	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::fixed, ios::floatfield);
	file.precision(6);

	if(ngr>0)
		file<<"g(r) data sets details:"<<endl;
	
	px=gr_rvalues;
	py=gr_gvalues;
	pz=gr_coeffs;
	for(i=0;i<ngr;i++)
	{
		file<<"\nset no. "<<i+1<<"/"<<ngr<<endl;
		file<<"original number of points->\t"<<*(grsize+i)<<endl;
		file<<"points used from \t"<<*(grmin+i)<<" to \t"<<*(grmax+i)<<" number of points used->\t"	<<*(grused+i)<<endl;
		file<<"partials coeffs: ";
		for(j=0;j<npartials;j++)
		{
			file<<J10<<*pz++<<"\t";
		}
		file<<"constant to be substracted->\t"<<*(grsub+i)<<endl;
		file<<"weight sigma->\t"<<*(grsigma+i)<<endl;
		file<<"renormalization switch->\t"<<*(grrenorm+i)<<endl;
#ifdef _R_SWITCH_GR_MIN_1
		file << "\n r \t ";
		if (r_switch_power == 1)
			file << "r*";
		else
			file << "(r^" << r_switch_power << ")*";
		file << "[g(r) - 1] is fitted!" << endl;
#else
		file << "\n r \t g(r) " << endl;
#endif
		for(j=0;j<*(grused+i);j++)
		{
			file<<J10<<*px++<<"\t"<<J10<<*py++<<endl;
		}
	}

	if (nsq>0) 
		file<<"\n S(Q) data sets details:"<<endl;

	px=sq_qvalues;
	py=sq_svalues;
	pz=sq_coeffs;
	for(i=0;i<nsq;i++)
	{
		file<<"\nset no. "<<i+1<<"/"<<nsq<<endl;
		file<<"original number of points->\t"<<*(sqsize+i)<<endl;
		file<<"points used from \t"<<*(sqmin+i)<<" to \t"<<*(sqmax+i)<<" number of points used->\t"\
		<<*(sqused+i)<<endl;
		file<<"partials coeffs: ";
		for(j=0;j<npartials;j++)
		{
			file<<J10<<*pz++<<"\t";
		}
		file<<"constant to be substracted->\t"<<*(sqsub+i)<<endl;
		file<<"weight sigma->\t"<<*(sqsigma+i)<<endl;
		file<<"renormalization switch->\t"<<*(sqrenorm+i)<<endl;
		file<<"offset switch->\t"<<*(sqoffset+i)<<endl;
		file<<"linear background switch->\t"<<*(sqlinear+i)<<endl;
		file<<"quadratic background switch->\t"<<*(sqquadratic+i)<<endl;

#ifdef _MUL_SCAT_VECTOR
		file<<"\n r \t Q*S(Q) "<<endl; 
#else
		file<<"\n Q \t S(Q) "<<endl;
#endif
		for(j=0;j<*(sqused+i);j++)
		{
			file<<J10<<*px++<<"\t"<<J10<<*py++<<endl;
		}
	}

	if (nfq>0) file<<"\nF(Q) data sets details:"<<endl;
	file.precision(6);
	file.setf(ios::fixed, ios::floatfield);
	px=fq_qvalues;
	py=fq_fvalues;
	if (is_IQmucorr)
		pIQ = fq_IQvalues_ori;
	pA = fq_A;
	pB = fq_B;
	pIB = fq_IB;
	for(i=0;i<nfq;i++)//for each F(Q) data set
	{
		file<<"\nset no. "<<i+1<<"/"<<nfq<<endl;
		file<<"original number of points->\t"<<*(fqsize+i)<<endl;
		file<<"points used from \t"<<*(fqmin+i)<<" to \t"<<*(fqmax+i)<<" number of points used->\t"	<<*(fqused+i)<<endl;
		file<<"constant to be substracted->\t"<<*(fqsub+i)<<endl;
		file<<"weight sigma->\t"<<*(fqsigma+i)<<endl;
		file<<"renormalization switch->\t"<<*(fqrenorm+i)<<endl;
		file<<"offset switch->\t"<<*(fqoffset+i)<<endl;
		file<<"linear background switch->\t"<<*(fqlinear+i)<<endl;
		file<<"quadratic background switch->\t"<<*(fqquadratic+i)<<endl;
		if (fqIQbackgcorr[i])
			file << "background correction mu actual->\t" << *(fqmuact + i) << endl;
		if (fqAXS[i] > 0)
			file << "\tusing AXS for type " << fqAXS[i] << " with f' actual " << fqfprimeact[i * RunParams::ntypes + fqAXS[i] - 1] << endl;
		if ((fqAXS[i] > 0 && fqfprimecount[i] > 1) || (fqAXS[i] == 0 && fqfprimecount[i] > 0))//there are fixed f' using sets
		{
			file << "\tfixed f' using atom types:" << endl;
			for (int itype = 0; itype < RunParams::ntypes; itype++)
			{
				if ((itype != fqAXS[i] - 1) && fqfprimeindex[i] & (int)pow(2, itype))//this type has fixed f' on it is not AXS
					file << "\t" << itype + 1 << ". type with f': " << fqfprime[i * RunParams::ntypes + itype] << endl;
			}
		}
	

#ifdef _MUL_SCAT_VECTOR
		file<<"\n r \t Q*F(Q) "<<endl; 
#else
		file << "\n Q \t "; 
		if (fqfitIQ[i])
		{
			file << "I(Q) \t ";
			if (fqIQbackgcorr[i])
				file << "I(Q)_corr \t ";
			file<<"partials coeffs \t A(Q) \t B(Q)";
			if (fqIQbackgcorr[i])
				file<< " \t IB(Q)";
		}
		else
			file<<"F(Q) \t partials coeffs ";

		file<< endl;
#endif
		
		for(j=0;j<*(fqused+i);j++)//for each data point
		{
			file.precision(6);
			file.setf(ios::fixed, ios::floatfield);
			file << J10 << *px++ << "\t" << J10;
			if (fqIQbackgcorr[i])
				file << *pIQ++ << "\t" << J10;
			else
				if (is_IQmucorr)
					pIQ++;//to advance the pointer, if needed later
			file<<*py++;//Q value and F(Q) value
			file.setf(ios::scientific, ios::floatfield);
			file.precision(8);
			for(k=0;k<npartials;k++)//for each partial
				file<<"\t"<<J16<<*(fq_cfinder[i*npartials+k]+j);//the partials coeffs for this Q value
			if (fqfitIQ[i])
			{
				file << "\t" << J16 << *pA++ << "\t" << J16 << *pB++;
				if (fqIQbackgcorr[i])
					file << "\t" << J16 << *pIB++;
			}
			else
			{
				if (is_IQ>0)
				{
					pA++;
					pB++;

					if (is_IQmucorr)
						pIB++;
				}
			}
			file<<endl;
		}
	}

	if (nfg > 0) file << "\nF(g) data sets details:" << endl;
	file.precision(6);
	file.setf(ios::fixed, ios::floatfield);
	px = fg_gvalues;
	py = fg_fvalues;
	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		file << "\nset no. " << i + 1 << "/" << nfg << endl;
		file << "original number of points->\t" << *(fgsize + i) << endl;
		file << "points used from \t" << *(fgmin + i) << " to \t" << *(fgmax + i) << " number of points used->\t"\
			<< *(fgused + i) << endl;
		file << "constant to be substracted->\t" << *(fgsub + i) << endl;
		file << "weight sigma->\t" << *(fgsigma + i) << endl;
		file << "renormalization switch->\t" << *(fgrenorm + i) << endl;
		file << "offset switch->\t" << *(fgoffset + i) << endl;
		file << "linear background switch->\t" << *(fglinear + i) << endl;
		file << "guadratic background switch->\t" << *(fgquadratic + i) << endl;
#ifdef _MUL_SCAT_VECTOR
		file << "\n r \t g*F(g) " << endl;
#else
		file << "\n g \t F(g) \t partials coeffs " << endl;
#endif
		for (j = 0; j < *(fgused + i); j++)//for each data point
		{
			file.setf(ios::fixed, ios::floatfield);
			file.precision(6);
			file << J10 << *px++ << "\t" << J10 << *py++;//g value and F(g) value
			file.setf(ios::scientific, ios::floatfield);
			file.precision(8);
			for (k = 0; k < npartials; k++)//for each partial
				file << "\t" << J16 << *(fg_cfinder[i*npartials + k] + j);//the partials coeffs for this g value
			file << endl;
		}
	}

	if (nek>0) file<<"\nE(k) data sets details:"<<endl;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	
	
	py=ek_evalues;
	pz=ek_coeffs;
	index = 0;//the index of the k-set in ek_kfinder
	for(i=0;i<nek;i++)//for each E(k) data set
	{
		px = ek_kfinder_ori[i];//this points to the beginning of the original k values of this data set
		file<<"\nset no. "<<i+1<<"/"<<nek<<endl;
		file << "E0 shift applied with dE0->" << ek_dE0[i] << " eV with number of grid points in one direction->" << ek_ngrid_in[i] << endl;
		file << "Grid points used from->" << ek_gridfrom[i] << " to " << ek_gridto[i] << endl;
		file<<"original number of points->\t"<<*(eksize+i)<<endl;
		file<<"points used from \t"<<*(ekmin+i)<<" to \t"<<*(ekmax+i)<<" number of points used->\t"\
		<<*(ekused+i)<<endl;
		for (j=0;j<ntypes;j++)
			file<<"r points used for "<<j+1<<". partial from \t"<<*(ek_rmin+i*ntypes+j)<<" to \t"<<*(ek_rmax+i*ntypes+j)<<endl;
		file<<"weight sigma->\t"<<*(eksigma+i)<<endl;
		file<<"power for experimental data (chi(k)*k^power->\t"<<*(ekchipower+i)<<endl;
		file<<"renormalization switch->\t"<<*(ekrenorm+i)<<endl;
		file<<"offset switch->\t"<<*(ekoffset+i)<<endl;
		file<<"type of the absorbing particle->\t"<<*(ek_abstype+i)<<endl;
		file<<"\n k \t E(k)=chi(k)*(k/k_max)^"<<*(ekchipower+i)<<endl; 
		
		for(j=0;j<*(ekused+i);j++)//for each data point
			file<<J10<<*px++<<"\t"<<J10<<*py++<<endl;;//k value and E(k) value
		
		if (is_E0shift==false || write_exafs_coeffs)
		{
			file << "\n Normalized coefficients(k,r) for " << i + 1 << ". E(k) set " << endl;
			file << "Number of absorbing particles \t" << SimpleCfg::pnatoms[ek_abstype[i]] << endl;
			
			for (ig = 0; ig < ek_ngrid[i]; ig++)
			{
				for (itype = 0; itype < ntypes; itype++)
				{
					px = ek_kfinder[index];
					if (is_E0shift)
						file << "Grid index "<<ig + ExptsData::ek_gridfrom[i] <<", type " << itype << endl;
					else
						file << "Type " << itype << endl;
					file << "k\tcoeffs" << endl;
					for (j = 0; j < *(ekused + i); j++)//for each k data point
					{
						file << J20<<*px++<<"\t";
						for (k = 0; k < ek_rused[i * ntypes + itype]; k++)//for each r point
							file << J20 << *pz++ << "\t";
						file << endl;
					}
				}
				index++;
			}
		}
		else
			file << "The coefficients will not be written as " << CheckKey("WRITE-EXAFS-COEFFS") << " is false" << endl;
				
	}
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}//end of the save() function

//normalize the EXAFS coefficients with the number of the edge particles, and (k/k_max)^ekchipower to make computation quicker
void ExptsData::NormEKCoeff()
{
	int iexpt, ig,ik, ir, itype, n_edge;
	double *pc, *pk, k_max;


	pc = ek_coeffs;
	for (iexpt = 0; iexpt < nek; iexpt++)//for each E(k) data set
	{
		n_edge = SimpleCfg::pnatoms[ek_abstype[iexpt]];
		k_max = *(ek_kfinder_ori[iexpt] + ekused[iexpt] - 1);
		for (ig = 0; ig < ek_ngrid[iexpt]; ig++)
		{
			for (itype = 0; itype < ntypes; itype++)
			{
				pk = ek_kfinder_ori[iexpt];
				for (ik = 0; ik < ekused[iexpt]; ik++)
				{
					for (ir = 0; ir < ek_rused[iexpt * ntypes + itype]; ir++)
					{
						*pc *= pow(*pk / k_max, double(ekchipower[iexpt])) / n_edge;
						pc++;
					}
					pk++;//next k value
				}
			}
		}
	}
};

//check, whether the data set is consistent, there are enough data point, and correct of possible
//checking the rspacing in the gr data files and resetting the rspacing in RunParams::rspacing to match the rspacing
//in the files, if necessary. It has to be done unfortunately before AssignBinsize is called
void ExptsData::CheckgrSets()
{
	int i, j;
	double r, r_prev=0,  mydr = -1, mydr_min;
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
	double r_first = 0, myrmin = 0, myshift = 0;
#endif
#endif
	ifstream file;
	bool is_equidistant;
	SetArraysize(&RunParams::gr_rspacing_ori, ngr, "RunParams::gr_rspacing_ori", "ExptsData::CheckgrSets");//might be needed, if gr spacing is reset

	for (i = 0; i < ngr; i++)
	{
		mydr_min = 10000;
		is_equidistant = true;
		RunParams::gr_rspacing_ori[i] = rspacing[i];

		cout << "\nLoading the r points from g(r) file: " << datafilename[i] << " to check the spacing" << endl;
		SafeOpenTextFile(file, datafilename[i]);
		if (CheckFileState(file, "ExptsData::ExptsData", datafilename[i]) == 0)
		{
			cout << "Cannot load data, exiting..." << endl;
			CleanExit();
		};

		grsize[i] = ReadThisLine(file, 1, 1, "g(r) ntot", "ExptsData::CheckgrSets", datafilename[i]);//total number of r points in the data file


		if (grsize[i] < 0)// can happen, if the number of data points are not given in the file
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The number of used data points is smaller, than zero, something must be wrong with the format of the file " << datafilename[i] << " !" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}
		if (grsize[i] < grmin[i])// can happen, if the number of data points are not given in the file
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The index of the first used data point (" << grmin[i] << " is smaller, than the number of data points in the file " << grsize[i] << " according to the file header of " << datafilename[i] << " !" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
		}

		SkipLine(file, datafilename[i], 1);//comment line

		for (j = 1; j < grmin[i]; j++)
			SkipLine(file, datafilename[i], 1);//skipping all the unused points

		for (j = 0; j < grmax[i] - grmin[i] + 1; j++)//for all used r points
		{
			//here we should check, whether there is enough data for the datafile
			r = ReadThisLine(file, 4, -1.0, "r for g(r)", "ExptsData::CheckgrSets", datafilename[i]);//loads the r value

			if (fabs(r + 1) < LOAD_TOL)//no data
			{
				cout << "\nWARNING(" << ++warn << "): There is not enough ";
				if (grsize[i] <= grmin[i] + j - 1)
					cout << "used ";//the file is consistent if grsize[i] = grmin[i], only the last used point is outside the range, otherwise there is more points in it as given in he header
				cout << "data point in the " << datafilename[i] << " file, only " << grmin[i] + j - 1 << " data points were found!" << endl;
				cout << "\tResetting the ";
				if (grsize[i] > grmin[i] + j - 1)
					cout << "number of points and ";//the file is consistent, only the last used point is outside the range
				cout << "last point to use to " << grmin[i] + j - 1 << " and trying to continue..." << endl;
				if (grsize[i] <= grmin[i] + j - 1)
					cout << "\tNevertheless correct the " << datfilename << " file, is possible!" << endl;
				else
					cout << "\tNevertheless correct the " << datafilename[i] << " file, is possible!" << endl;

				grmax[i] = grmin[i] + j - 1;
				grsize[i] = grmin[i] + j - 1;
				break;

			}
			SkipLine(file);
			//trying to establish, whether it is equidistant
			
			if (j == 1)
			{
				mydr = r - r_prev;//the difference of the first and second data point
				
			}
			else
			{
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
				if (j == 0)
					r_first = r;
#endif
#endif
				if (j > 1 && is_equidistant)
				{
					if (fabs(mydr - (r - r_prev)) > LOAD_TOL)
						is_equidistant = false;
				}
			}
			if (j>=1 && mydr_min > r - r_prev)
				mydr_min = r - r_prev;//keep the smallest dr
			r_prev = r;

		}//end of r checking cycle j
		//check, whether there is enough point in the file, even if they would not be needed
		for (j = grmax[i]; j < grsize[i]; j++)//for all used r points
		{
			//here we should check, whether there is enough data for the datafile
			r = ReadThisLine(file, 4, -1.0, "r for g(r)", "ExptsData::CheckgrSets", datafilename[i]);//loads the r value
			if (fabs(r + 1) < LOAD_TOL)//no data
			{
				//the file is not consistent, although the missing point are not needed
				cout << "\nWARNING(" << ++warn << "): There is only " << j << " data points in the " << datafilename[i] << " file," << endl;
				cout << "\talthough " << grsize[i] << " data points were given in the header!" << endl;
				cout << "\tAlthough the missing point(s) are not needed for this simulation, correct the file to avoid this warning again!" << endl;
				cout << "\tThe number of points will be reset to " << j << "!" << endl;
				cout << "\tlast point to use to " << grmin[i] + j - 1 << " and trying to continue..." << endl;
				grsize[i] = j;
				break;

			}
			if (j < grsize[i] - 1)
				SkipLine(file);

		}//end of r checking cycle j
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
		if (is_equidistant)
		{

			if (fabs(mydr - rspacing[i]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): The spacing of the " << i + 1 << ". g(r) data set is not the same as the default or custom value for this set!" << endl;
				cout << "\tAs it is equdistant the bin size will be reset to " << mydr << " to avoid smoothing of the data set!\n" << endl;
				rspacing[i] = mydr;
				bin_flag[i] = 1;
			}
		}
		else
		{
			if (mydr_min < rspacing[i])//the chosen rspacing is not adequate, create a new one
			{
				cout << "\nWARNING(" << ++warn << "): The spacing of the non-equidistant " << i + 1 << ". g(r) data set was larger (" << rspacing[i] << ") than the smallest" << endl;
				cout << "\tdifference between consecutive data points (" << mydr_min << ")! Binsize will be reset to " << mydr_min << " to avoid smoothing problems!\n" << endl;
				rspacing[i] = mydr_min;
				bin_flag[i] = 2;
			}
		}
		//now check, whether the default binshift is adequate, or custom has to be applied, to get the proper data r points
		//keep in mind, that the r points are representing the middle of the bins. If this binshift is not equal to the default,
		//set a custom binshift
		myrmin = r_first - rspacing[i] / 2;
		myshift = myrmin / rspacing[i];
		if (myshift < 0)//this is not likely
		{
			cout << "\n*****ERROR******" << endl;
			cout << "The binshift is negative for the " << i + 1 << ". g(r) data set!" << endl;
			cout << "This means that the first bin would not be a whole bin ofr the histogram calculation!" << endl;
			if (is_equidistant)
				cout << "The equidistant spacing ("<<rspacing[i]<<" A) was determined form the data set." << endl;
			else
			{
				cout << "The r-points are not equidistant";
				if (bin_flag[i]==2)
					cout << ", the " << rspacing[i] << " bin size was determined from the smallest dr in the data set." << endl;
				else
					cout << "custom rspacing (" << rspacing[i] << " A) was used." << endl;
				cout << "The first r point is " << r_first << ", which would put the beginning of the first bin at " << myrmin << " A..." << endl;
				cout << "Use smaller custom bin size, or delete the first data point!" << endl;
				cout << "Cannot run this way, exiting..." << endl;
				CleanExit();
			}
		}
		//try to round it

		if (fabs(RunParams::binshift_def - myshift) > LOAD_TOL)
		{
			//custom binshift has to be created

			if (fabs(myshift) < LOAD_TOL)
			{
				myshift = 0;
				RunParams::firstbin[i] = 0;
			}
			else
			{
				//as due to rounding error in case of binshift 1, 2... the (int)binshift is not necessarily correct, an integer firstbin will be calculated here
				if (fabs(myshift - round(myshift)) < LOAD_TOL)//this supposed to be integer
					RunParams::firstbin[i] = (int)round(myshift);
				else
					RunParams::firstbin[i] = (int)myshift;
			}
			RunParams::binshift[i] = myshift;
			RunParams::binshift_flag = true;
		}
	


#endif
#endif
		if (!CheckReadFileState(file, "ExptsData::CheckgrSets", datafilename[i]))
			CleanExit();
		next_line_pos = -1;
		file.close();
		
	}//next g(r) constraint
}


//Reading the main parameters for the experimental data. Returns with false, if incoming datasets are insufficient.
bool ExptsData::GetExptParamFree(std::list<std::list<string> > pool)
{
	int cumul = 0;
	bool retval = true;
	bool Compton_chemsymb = false;//to indicate, whether chemical symbols are needed for Compton calculation 
	string key, params,own;
	char *name=NULL, *conv_numb = NULL;
	
	SetArraysize(&name, NAME_SIZE, "name", "ExptsData::GetExptParamFree");

	ntypes = RunParams::ntypes;
	int npartials = ntypes * (ntypes + 1) / 2;
	ngr = RunParams::ngr;
	nsq = RunParams::nsq;
	nfq = RunParams::nfq;
	nfg = RunParams::nfg;
	nek = RunParams::nek;
	ntot_datasets = RunParams::ntot_datasets;
	assign_hist = RunParams::assign_hist;
	rspacing = RunParams::rspacing;


	//first the g(r) then the S(Q),F(Q), F(g) finally the E(k) data files name are stored
	datafilename = new char[ngr + nsq + nfq + nfg + 2 * nek][FILE_NAME_SIZE];//there data and the coefficients are in separate files for EXAFS
	if (datafilename == NULL)
		NoArray("ExptsData::GetExptsParamFree", "datafilename");

	//Arrays are allocated regardless if there is a constraint of that type
	SetArraysize(&grmin, ngr, "grmin", "ExptsData::GetExptsParamFree");
	SetArraysize(&grmax, ngr, "grmax", "ExptsData::GetExptsParamFree");
	SetArraysize(&grsize, ngr, "grsize", "ExptsData::GetExptsParamFree");
	SetArraysize(&grsigma, ngr, "grsigma", "ExptsData::GetExptsParamFree");
	SetArraysize(&grused, ngr, "grused", "ExptsData::GetExptsParam");
	SetArraysize(&grused_cum, ngr + 1, "grused_cum", "ExptsData::GetExptsParamFree");
	SetArraysize(&grrenorm, ngr, "grrenorm", "ExptsData::GetExptsParamFree");
	SetArraysize(&groffset, ngr, "groffset", "ExptsData::GetExptsParamFree");
	SetArraysize(&grlinear, ngr, "grlinear", "ExptsData::GetExptsParamFree");
	SetArraysize(&grquadratic, ngr, "grquadratic", "ExptsData::GetExptsParamFree");
	SetArraysize(&grcubic, ngr, "grcubic", "ExptsData::GetExptsParamFree");
	SetArraysize(&grsub, ngr, "grsub", "ExptsData::GetExptsParamFree");
	SetArraysize(&gr_coeffs, ngr*npartials, "gr_coeffs", "ExptsData::GetExptsParamFree");
	SetArraysize(&bin_flag, ngr, "bin_flag", "ExptsData::GetExptsParamFree");
	SetArraysize(&gruseR, ngr, "gruseR", "ExptsData::GetExptsParamFree");

	SetArraysize(&sqmin, nsq, "sqmin", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqmax, nsq, "sqmax", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqsize, nsq, "sqsize", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqused, nsq, "sqused", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqused_cum, nsq + 1, "sqused_cum", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqrenorm, nsq, "sqrenorm", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqoffset, nsq, "sqoffset", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqlinear, nsq, "sqlinear", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqquadratic, nsq, "sqquadratic", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqcubic, nsq, "sqcubic", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqsub, nsq, "sqsub", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqsigma, nsq, "sqsigma", "ExptsData::GetExptsParamFree");
	SetArraysize(&sq_coeffs, nsq*npartials, "sq_coeffs", "ExptsData::GetExptsParamFree");
	SetArraysize(&squseR, nsq, "squseR", "ExptsData::GetExptsParamFree");
	SetArraysize(&sqreadcoeffs, nsq, "sqreadcoeffs", "ExptsData::GetExptsParamFree");
	SetArraysize(&isotope_count, nsq*ntypes, "isotope_count", "ExptsData::GetExptsParamFree");

	SetArraysize(&fqmin, nfq, "fqmin", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqmax, nfq, "fqmax", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqsize, nfq, "fqsize", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqused, nfq, "fqused", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqused_cum, nfq + 1, "fqused_cum", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqrenorm, nfq, "fqrenorm", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqoffset, nfq, "fqoffset", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqlinear, nfq, "fqlinear", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqquadratic, nfq, "fqquadratic", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqcubic, nfq, "fqcubic", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqsub, nfq, "fqsub", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqsigma, nfq, "fqsigma", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqreadcoeffs, nfq, "fqreadcoeffs", "ExptsData::GetExptsParamFree");
	SetArraysize(&fquseR, nfq, "fquseR", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqfitIQ, nfq, "fqfitIQ", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqrenalpha, nfq, "fqrenalpha", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqusecompton, nfq, "fqusecompton", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqfprimecount, nfq, "fqfprimecount", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqfprime, nfq * ntypes, "fqfprime", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqfprimeact, nfq * ntypes, "fqfprimeact", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqfprimeindex, nfq, "fqfprimeindex", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqfprimefactor, nfq, "fqfprimefactor", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqAXS, nfq, "fqAXS", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqalpha, nfq, "fqalpha", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqa, nfq, "fqa", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqb, nfq, "fqb", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqIQbackgcorr, nfq, "fqIQbackgcorr", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqmu, nfq, "fqmu", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqdmumax, nfq, "fqdmumax", "ExptsData::GetExptsParamFree");
	SetArraysize(&fqmuact,nfq, "fqmuact", "ExptsData::GetExptsParamFree");
	
	SetArraysize(&fgmin, nfg, "fgmin", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgmax, nfg, "fgmax", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgsize, nfg, "fgsize", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgused, nfg, "fgused", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgused_cum, nfg + 1, "fgused_cum", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgrenorm, nfg, "fgrenorm", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgoffset, nfg, "fgoffset", "ExptsData::GetExptsParamFree");
	SetArraysize(&fglinear, nfg, "fglinear", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgquadratic, nfg, "fgquadratic", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgcubic, nfg, "fgcubic", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgsub, nfg, "fgsub", "ExptsData::GetExptsParamFree");
	SetArraysize(&fgsigma, nfg, "fgsigma", "ExptsData::GetExptsParamFree");
	SetArraysize(&fguseR, nfg, "fguseR", "ExptsData::GetExptsParamFree");

	SetArraysize(&ek_rmin, ntypes*nek, "ek_rmin", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_rmax, ntypes*nek, "ek_rmax", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_rused, ntypes*nek, "ek_rused", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekmin, nek, "ekmin", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekmax, nek, "ekmax", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekmin_ori, nek, "ekmin_ori", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekmax_ori, nek, "ekmax_ori", "ExptsData::GetExptsParamFree");
	SetArraysize(&eksize, nek, "eksize", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekused, nek, "ekused", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekused_cum, nek + 1, "ekused_cum", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekrenorm, nek, "ekrenorm", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekoffset, nek, "ekoffset", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_abstype, nek, "ek_abstype", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekchipower, nek, "ekchipower", "ExptsData::GetExptsParamFree");
	SetArraysize(&eksigma, nek, "eksigma", "ExptsData::GetExptsParamFree");
	SetArraysize(&ekuseR, nek, "ekuseR", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_dE0, nek, "ek_dE0", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_ngrid, nek, "ek_ngrid", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_ngrid_in, nek, "ek_ngrid_in", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_gridind, nek, "ek_gridind", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_gridstart, nek, "ek_gridstart", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_gridfrom, nek, "ek_gridfrom", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_gridto, nek, "ek_gridto", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_cf_cum, nek+1, "ek_cf_cum", "ExptsData::GetExptsParamFree");
	SetArraysize(&ek_ngrid_cum, nek + 1, "ek_ngrid_cum", "ExptsData::GetExptsParamFree");
	  
	for (int i = 0; i < nsq * ntypes; i++)
	{
		isotope_count[i].count = 0;//default, natural abundance
		isotope_count[i].isotopes = nullptr;
	}

	for (int i = 0; i <nfq;i++)
	{
		//set to default
		fqIQbackgcorr[i] = IQBACKGCORR_DEF;
		fqmu[i] = IQ_MU_DEF;
		fqmuact[i] = IQ_MU_DEF;
		fqdmumax[i] = IQ_DMU_MAX_DEF;
		fqfprimeindex[i] = 0;//needed for addition
		fqfprimefactor[i] = 0.0;
		
	}
	for (int i = 0; i < nfq * ntypes; i++)
	{
		fqfprime[i] = 0.0;
		fqfprimeact[i] = 0.0;
	}
	//the indicator array for each g(r), F(Q) and S(Q)  data series whether to use cubic correction, 
	SetArraysize(&use_cubic, ngr + nsq + nfq + nfg, "ExptsData::use_cubic", "ExptsData::GetExptParamsFree");

	std::list<std::list<string>>::iterator it2 = pool.begin();
	int sig;
	for (int i = 0; (i < ngr) && retval; i++, cumul++, it2++)
	{
		sig = 0;
		double liminunits[2] = { -1.0, -1.0 };
		IntToStr(&conv_numb, i + 1);
		grmin[i] = -1;
		grmax[i] = -1;
		bin_flag[i] = 0;
		//Initialize by default values
		grsub[i] = SUBTRACT_DEF;
		grrenorm[i] = RENORM_DEF;
		groffset[i] = grlinear[i] = grquadratic[i] = grcubic[i] = RENORM_DEF;
		use_cubic[cumul] = USE_CUBIC_DEF; 
		gruseR[i] = USE_RFACTOR_DEF;

		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{
			
			WrapFree(*it, key, params);
			stringstream ss(params);
	
			if (key == CheckKey("TYPE"))
			{
				if (!(ss >> own)) return ErrLackValue(*it);
				if (own != "GR")
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Please, check the EXP datasets: only one TYPE entry is allowed for each!" << endl;
					return false;
				}
			}
			else if (key == CheckKey("DATAFILE"))
			{
				params.erase(std::remove(params.begin(), params.end(), ' '), params.end());
				std::copy(params.begin(), params.end(), datafilename[cumul]);
				(datafilename[cumul])[params.length()] = 0;
			}
			else if (key == CheckKey("POINT-RANGE"))
			{
				if (!(GetInt(ss, *it, grmin + i, "ExptsData::GetExptParamFree") && GetInt(ss, *it, grmax + i, "ExptsData::GetExptParamFree")))
				{
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}							
			}
			else if (key == CheckKey("R-RANGE"))
			{
				bool ret = true;
				for (int j = 0; (j < 2) && ret; j++) ret = static_cast<bool>(ss >> liminunits[j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("R-SPACING"))
			{
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
				if (!(ss >> RunParams::rspacing[cumul])) return ErrLackValue(*it);
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
			}
			else if (key == CheckKey("CONST-SUBTRACT"))
			{
				if (!(ss >> grsub[i])) return ErrLackValue(*it);
			}
			else if (key == CheckKey("PARTIAL-COEFFS"))
			{
				bool ret = true;
				for (int j = 0; (j < npartials) && ret; j++) ret = static_cast<bool>(ss >> gr_coeffs[i*npartials + j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide " << npartials << " values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (sig++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple SIGMA(-MASTER)/SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckKeyValue("TYPE","GR") << " constraint!" << endl << "\tThe last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == cumul + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!(ss >> grsigma[i])) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA")||key == CheckKey("SIGMA-MASTER"))
					grsigma[i] = fabs(grsigma[i]);
				if (key == CheckKey("SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1) RunParams::lead_series_ind = cumul + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry was declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						cout << "\tIt will be reset to default " << LEAD_SERIES_IND_DEF << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					grsigma[i] = -fabs(grsigma[i]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					grsigma[i] = fabs(grsigma[i]);
				}
			}
			else if (key == CheckKey("USE-RFACTOR"))
			{
				if (!(ss >> gruseR[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE,"USE-RFACTOR for gr set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(gruseR[i],name,"ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("RENORM"))
			{
				if (!(ss >> grrenorm[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "RENORM for gr set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(grrenorm[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("POLY-BACK-FLAGS"))
			{
				int last = -1;
				bool rv = true;
				int initvals[4];
				for (int j = 0; j < 4; j++)
				{
					if (rv) rv = GetInt(ss, *it, initvals+j, "ExptsData::GetExptParamFree");
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "Too few parameters. You should provide at least one value in the following line: " << endl << *it << endl;
							return false;
						}
						else initvals[j] = 0;
					}
					else last = j;
				}
				if ((last != 4 - 1) && ::debug)
				{
					cout << "\nNOTE(" << ++note << "): polynomial background correction flags from " << last + 1 << ". to 3rd order are/is set to zero at the " << endl;;
					cout<< "\t "<<i + 1 << ". g(r) experimental data." << endl;
				}
				use_cubic[cumul] = (initvals[3] || initvals[2] || initvals[1] || initvals[0]) ? 1 : 0;
				groffset[i] = initvals[0];
				grlinear[i] = initvals[1];
				grquadratic[i] = initvals[2];
				grcubic[i] = initvals[3];
				for (int j = 0; j < 4; j++)
				{
					mystrcpy(name, NAME_SIZE, "POLY-BACK-FLAGS ");
					switch (j)
					{
					case 0:
						mystrcat(name, NAME_SIZE, "0th term for gr set ");
						break;
					case 1:
						mystrcat(name, NAME_SIZE, "1st term for gr set ");
						break;
					case 2:
						mystrcat(name, NAME_SIZE, "2nd term for gr set ");
						break;
					case 3:
						mystrcat(name, NAME_SIZE, "3rd term for gr set ");
						break;
					}
					mystrcat(name, NAME_SIZE, conv_numb);
					Check0_1(initvals[j], name, "ExptsData::GetExptParamFree");
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("EXP"), CheckKeyValue("TYPE", "GR"));//write the appropriate error mesage
				return false;
			}

		}

		//Open experimental data to read to clarify grmin/grmax,etc.
		retval = FindLimits(datafilename[cumul], grmin[i], grmax[i], liminunits[0], liminunits[1]);
		if (retval)
		{
			cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". GR data set point range " << grmin[i] << ".." << grmax[i] << " corresponds to " << liminunits[0] << ".." << liminunits[1] << " in A" << endl;
			rused_tot += grmax[i] - grmin[i] + 1;
		}
		else
			return retval;//no reason to continue
	}

	ExptsData::CheckgrSets();
	isotope_ratio_type *iratio=nullptr;
	for (int i = 0; (i < nsq) && retval; i++, cumul++, it2++)
	{
		sig = 0;
		double liminunits[2] = { -1.0, -1.0 };
		IntToStr(&conv_numb, i + 1);
		sqmin[i] = -1;
		sqmax[i] = -1;
		//Initialize by default values
		sqreadcoeffs[i] = READ_COEFFS_DEF;
		sqsub[i] = SUBTRACT_DEF;
		sqrenorm[i] = RENORM_DEF;
		sqoffset[i] = sqlinear[i] = sqquadratic[i] = sqcubic[i] = RENORM_DEF;
		use_cubic[cumul] = USE_CUBIC_DEF; 
		squseR[i] = USE_RFACTOR_DEF;

		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("TYPE"))
			{
				string own;
				if (!(ss >> own)) return ErrLackValue(*it);
				if (own != "ND")
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Please, check the EXP datasets: only one TYPE entry is allowed for each!" << endl;
					return false;
				}
			}
			else if (key == CheckKey("DATAFILE"))
			{
				params.erase(std::remove(params.begin(), params.end(), ' '), params.end());
				std::copy(params.begin(), params.end(), datafilename[cumul]);
				(datafilename[cumul])[params.length()] = 0;
			}
			else if (key == CheckKey("POINT-RANGE"))
			{
				if (!(GetInt(ss, *it, sqmin + i, "ExptsData::GetExptParamFree") && GetInt(ss, *it, sqmax + i, "ExptsData::GetExptParamFree")))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("Q-RANGE"))
			{
				bool ret = true;
				for (int j = 0; (j < 2) && ret; j++) ret = static_cast<bool>(ss >> liminunits[j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("R-SPACING"))
			{
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
				if (!(ss >> RunParams::rspacing[cumul])) return ErrLackValue(*it);
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
			}
			else if (key == CheckKey("CONST-SUBTRACT"))
			{
				if (!(ss >> sqsub[i])) return ErrLackValue(*it);
			}
			else if (key == CheckKey("PARTIAL-COEFFS"))
			{
				bool ret = true;
				for (int j = 0; (j < npartials) && ret; j++) ret = static_cast<bool>(ss >> sq_coeffs[i*npartials + j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide " << npartials << " values in the following line: " << endl << *it << endl;
					return false;
				}
				sqreadcoeffs[i] = 1;//indicaing, that the coeffs were read
			}
			else if (key == CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS"))
			{
				int mycount;
				if (!(ss >> mycount)) return ErrLackValue(*it);
				if (mycount < 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Number of isotopes to be used should be 0 (natural abundance) or positive number in line!\n" << *it << endl;
					return false;

				}
				else if (mycount > 0)
				{

					SetArraysize(&iratio, mycount, "isotop_ratio", "ExptsData::GetExptParamsFree");
					for (int j = 0; j < mycount; j++)
						if (!(ss >> iratio[j].symbol)) return ErrLackValue(*it);
					double sum = 0;
					for (int j = 0; j < mycount; j++)
					{
						if (!(ss >> iratio[j].ratio)) return ErrLackValue(*it);
						sum += iratio[j].ratio;
					}
					if (fabs(sum - 1.0) > LOAD_TOL)
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The sum of the isotope ratios in line\n" << *it << endl;
						cout << "differs from 1.0 more than the tolerance " << LOAD_TOL << "!" << endl;
						cout << "Correct the " << datfilename << " file, and try again!" << endl;
						return false;
					}
				}
				int npos, cpos;
				string *mysymbols,temp;
				
				//this is for the element symbols
				SetArraysize(&mysymbols, mycount,"mysymbols", "ExptsData::GetExptParamsFree");
				
				for (int j = 0; j < mycount; j++)//check, whether the symbol is in the table 
				{
					mysymbols[j] = "";
					if (N_scat_length.count(iratio[j].symbol) > 0)//the isotope is in the table, extract the element's symbol
					{
						mysymbols[j] = N_scat_length.at(iratio[j].symbol).symbol;
						//check, if this is not a nuclid, but the element symbol only, 
						if (mysymbols[j].compare(iratio[j].symbol) == 0)
						{
							if (mycount > 1)//in case of more than 1 isotopes this would not make sense
							{
								cout << "\n*****ERROR*****" << endl;
								cout << "The " << j + 1 << ". symbol ("<<mysymbols[j]<<") in the line:\n" << *it << endl;
								cout << "is not an isotope, but the element symbol describing natural abundance!" << endl;
								cout << "This does not make sense, if natural abundance should be used, omit this key world!" << endl;
								cout << "Correct the " << datfilename << " file, and try again!" << endl;
								return false;
							}
							else//this key is not necessary, natural abundance sould be used, abandon processing
							{
								mycount = 0;
								isotope_count[i * ntypes + j].count = 0;
								break;
							}

						}
					}
					else//if the symbol was not given in the expected C14 format, convert if possible, if 14C was used
					{
						npos = 0;
						for (unsigned int k = 0; k < iratio[j].symbol.length(); k++)
						{
							if (std::isdigit(iratio[j].symbol[k]))
							{
								npos = k;//position of the first digit
								break;
							}
						}
						if (npos == 0)//most probably 14C format, convert it
						{
							cpos = 0;
							for (unsigned int k = 0; k < iratio[j].symbol.length(); k++)
							{
								if (std::isdigit(iratio[j].symbol[k]) == 0)//first non-digit
								{
									cpos = k;
									break;
								}
							}
							temp = iratio[j].symbol.substr(cpos);
							string temp2 = temp + iratio[j].symbol.substr(0, cpos);
							if (N_scat_length.count(temp2) == 0)//not in the table
							{
								cout << "\n*****ERROR*****" << endl;
								cout<<"Problem in line:\n"<< *it << endl;
								cout << "The isotope symbols should be given in the format like U238 (first the symbol followed by mass number." << endl;
								cout << "The " << iratio[j].symbol << " was converted into " << temp2 << ", but even so it is not found in the neutron coefficient" << endl;
								cout << "calculation table!" << endl;
								GiveNSymbols(0);
								return false;
							}
							else
							{
								//conversion worked, copy back
								iratio[j].symbol = temp2.substr(0);
								mysymbols[j] = N_scat_length.at(iratio[j].symbol).symbol;
							}
						}
					}
					
				}//end of j cycle for isotopes
				if (mycount > 0)
				{
					for (int j = 0; j < mycount; j++)
					{
						if (j > 0)//compare with the first one
						{
							if (mysymbols[j].compare(mysymbols[0]) != 0)
							{
								cout << "\n*****ERROR*****" << endl;
								cout << "Problem in line:\n" << *it << endl;
								cout << iratio[j].symbol << " and " << iratio[0].symbol << " are not isotopes of the same element!" << endl;
								cout << "Correct the " << datfilename << " file, and try again!" << endl;
								return false;
							}
						}
					}
					int found = 0;
					for (int j = 0; j < ntypes; j++)
					{
						//no need to check, if chem_symbols exist, it was checked already in RunParams::ReadFreeExp
						//and if chem_symbols exist, chem_symbols_standard exist too
						if (mysymbols[0].compare(RunParams::chem_symbols_standard[j])==0)//this belongs to the j-th type
						{
							if (isotope_count[i * ntypes + j].count > 0)//it was set already
							{
								cout << "\n*****ERROR*****" << endl;
								cout << "Problem in line:\n" << *it << endl;
								cout << "This key word was already given for element " << mysymbols[0] << " in case of the " << i + 1 << ". ND data set!" << endl;
								cout << "It can be given for each component of the system only once for a data set!" << endl;
								cout << "Correct the " << datfilename << " file, and try again!" << endl;
								return false;
							}
							else
							{
								isotope_count[i * ntypes + j].count = mycount;
								isotope_count[i * ntypes + j].isotopes = iratio;
							}
							found = 1;
						}

					}
					if (found == 0)//no match with chem_symbols
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "Problem in line:\n" << *it << endl;
						cout << "The element " << mysymbols[0] << " cannot be found among the symbols given by " << CheckKey("CHEMICAL-SYMBOLS") << endl;
						cout << "Correct the " << datfilename << " file, and try again!" << endl;
						return false;
					}
				}
				delete[] mysymbols;

			}
			else if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (sig++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple SIGMA(-MASTER)/SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckKeyValue("TYPE", "ND") << " constraint!" << endl << "\tThe last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == cumul + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!(ss >> sqsigma[i])) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA") || key == CheckKey("SIGMA-MASTER"))
					sqsigma[i] = fabs(sqsigma[i]);
				if (key == CheckKey("SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1) RunParams::lead_series_ind = cumul + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry was declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						cout << "\tIt will be reset to default " << LEAD_SERIES_IND_DEF << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					sqsigma[i] = -fabs(sqsigma[i]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					sqsigma[i] = fabs(sqsigma[i]);
				}
			}
			else if (key == CheckKey("USE-RFACTOR"))
			{
				if (!(ss >> squseR[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "USE-RFACTOR for ND set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(squseR[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("RENORM"))
			{
				if (!(ss >> sqrenorm[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "RENORM for ND set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(sqrenorm[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("POLY-BACK-FLAGS"))
			{
				int last = -1;
				bool rv = true;
				int initvals[4];
				for (int j = 0; j < 4; j++)
				{
					if (rv) rv = GetInt(ss, *it, initvals + j, "ExptsData::GetExptParamFree");
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "Too few parameters. You should provide at least one value in the following line: " << endl << *it << endl;
							return false;
						}
						else initvals[j] = 0;
					}
					else last = j;
				}
				if ((last != 4 - 1) && ::debug)
				{
					cout << "\nNOTE(" << ++note << "): polynomial background correction flags from " << last + 1 << ". to 3rd order are/is set to zero at the " << endl;
					cout<<"\t "<< i + 1 << ". S(Q) experimental data." << endl;
				}
				use_cubic[cumul] = initvals[3];
				sqoffset[i] = initvals[0];
				sqlinear[i] = initvals[1];
				sqquadratic[i] = initvals[2];
				sqcubic[i] = initvals[3];
				for (int j = 0; j < 4; j++)
				{
					mystrcpy(name, NAME_SIZE, "POLY-BACK-FLAGS ");
					switch (j)
					{
					case 0:
						mystrcat(name, NAME_SIZE, "0th term for ND set ");
						break;
					case 1:
						mystrcat(name, NAME_SIZE, "1st term for ND set ");
						break;
					case 2:
						mystrcat(name, NAME_SIZE, "2nd term for ND set ");
						break;
					case 3:
						mystrcat(name, NAME_SIZE, "3rd term for ND set ");
						break;
					}
					mystrcat(name, NAME_SIZE, conv_numb);
					Check0_1(initvals[j], name, "ExptsData::GetExptParamFree");
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("EXP"), CheckKeyValue("TYPE", "ND"));//write the appropriate error mesage
				return false;
			}

		}
		//Open experimental data to read to clarify sqmin/sqmax,etc.
		retval = FindLimits(datafilename[cumul], sqmin[i], sqmax[i], liminunits[0], liminunits[1] );
		if (retval)
		{
			cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". ND data set point range " << sqmin[i] << ".." << sqmax[i] << " corresponds to " << liminunits[0] << ".." << liminunits[1] << " in A^-1" << endl;
			sqused_tot += sqmax[i] - sqmin[i] + 1;
		}
		else
			return retval;//no reason to continue
	}

	
	for (int i = 0; (i < nfq) && retval; i++, cumul++, it2++)
	{
		sig = 0;
		double liminunits[2] = { -1.0, -1.0 };
		IntToStr(&conv_numb, i + 1);
		fqmin[i] = -1;
		fqmax[i] = -1;

		//Initialize by default values
		fqsub[i] = SUBTRACT_DEF;
		fqrenorm[i] = RENORM_DEF;
		fqoffset[i] = fqlinear[i] = fqquadratic[i] = fqcubic[i] = RENORM_DEF;
		use_cubic[cumul] = USE_CUBIC_DEF;
		fquseR[i] = USE_RFACTOR_DEF;
		fqrenalpha[i] = RENORM_DEF;
		fqfitIQ[i] = IQFIT_DEF;
		fqIQbackgcorr[i] = IQBACKGCORR_DEF;
		if (RunParams::dat_version == 2)
			fqreadcoeffs[i] = 1;
		else
			fqreadcoeffs[i] = READ_COEFFS_DEF;
		fqAXS[i] = 0;//cannot be changed
		fqusecompton[i] = COMPTON_DEF;
		fqfprimecount[i] = 0;//default no fprime, sum of how many types using fprime for the given set

		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("TYPE"))
			{
				string own;
				if (!(ss >> own)) return ErrLackValue(*it);
				if (own != "XRD")
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Please, check the EXP datasets: only one TYPE entry is allowed for each!" << endl;
					return false;
				}
			}
			else if (key == CheckKey("DATAFILE"))
			{
				params.erase(std::remove(params.begin(), params.end(), ' '), params.end());
				std::copy(params.begin(), params.end(), datafilename[cumul]);
				(datafilename[cumul])[params.length()] = 0;
			}
			else if (key == CheckKey("POINT-RANGE"))
			{
				if (!(GetInt(ss, *it, fqmin + i, "ExptsData::GetExptParamFree") && GetInt(ss, *it, fqmax + i, "ExptsData::GetExptParamFree")))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("Q-RANGE"))
			{
				bool ret = true;
				for (int j = 0; (j < 2) && ret; j++) ret = static_cast<bool>(ss >> liminunits[j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("R-SPACING"))
			{
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
				if (!(ss >> RunParams::rspacing[cumul])) return ErrLackValue(*it);
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
			}
			else if (key == CheckKey("CONST-SUBTRACT"))
			{
				if (!(ss >> fqsub[i])) return ErrLackValue(*it);
			}
			else if (key == CheckKey("I(Q)_A_B_ALPHA"))
			{
#ifdef _MUL_SCAT_VECTOR
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("I(Q)_A_B_ALPHA") << " cannot be used in case of Q*F(Q) fitting" << endl;
				cout << "\t(when the proram was compiled with  _MUL_SCAT_VECTOR compiler option. Ignoring it and continuing...." << endl;
#else
				if (!(ss >> fqfitIQ[i])) return ErrLackValue(*it);
				if (fqfitIQ[i])
				{
					if (is_IQ==0)
						is_IQ = 1;//only set it here, if it was not set before to indicate, there is at least one I(Q) fitting set, it will be later updated, if it is known, what kind of regression
					if (ss >> fqa[i])
					{
						if (!(ss >> fqb[i])) return ErrLackValue(*it);
						if (!(ss >> fqalpha[i])) fqalpha[i] = NO_DVALUE_DEF;//to indicate it was not read
					}
					else
					{
						//no custom values are supplied
						fqa[i] = NO_DVALUE_DEF;
						fqb[i] = NO_DVALUE_DEF;
						fqalpha[i] = NO_DVALUE_DEF;
					}
				}
#endif
			}
			else if (key == CheckKey("I(Q)BACKG_MU_DMU-MAX"))
			{
#ifdef _MUL_SCAT_VECTOR
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("I(Q)BACKG_MU_DMU-MAX") << " cannot be used in case of Q*F(Q) fitting" << endl;
				cout << "\t(when the proram was compiled with  _MUL_SCAT_VECTOR compiler option. Ignoring it and continuing...." << endl;
#else
				if (!(ss >> fqIQbackgcorr[i])) return ErrLackValue(*it);
				if (fqIQbackgcorr[i])
				{
					is_IQmucorr = 1;
					if (ss >> fqmu[i])
					{
						if (fqmu[i]<=0 || fqmu[i]>1)
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The background correction factor mu for the " << i + 1 << ". F(Q) data set is " << fqmu[i] <<","<< endl;
							cout << "which is outside the allowed range 0 < mu <= 1 ! Change it in " <<datfilename<<" file, and try again!"<< endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
						if (ss >> fqdmumax[i])
						{
							if (fqmu[i] - fqdmumax[i] <= 0 || fqmu[i] + fqdmumax[i] > 1)
							{
								cout << "\n*****ERROR*****" << endl;
								if (fqmu[i] - fqdmumax[i] <= 0)
									cout << "The minimum value of the actual background correction factor mu-dmu_max= " << fqmu[i] - fqdmumax[i] << endl;
								else
									cout << "The maximum value of the actual background correction factor mu+dmu_max= " << fqmu[i] + fqdmumax[i] << endl;
								cout << "for the " << i + 1 << ". F(Q) data set is outside the allowed range 0 < mu_actual <= 1 !" << endl;
								cout << "Change mu or dmu_max in " << datfilename << " file, and try again!" << endl;
								cout << "Exiting..." << endl;
								CleanExit();
							}
						}
					}
				}
#endif
			}
			else if (key == CheckKey("FPRIME-INDEX_FPRIME"))
			{
#ifdef _MUL_SCAT_VECTOR
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("FPRIME-INDEX_FPRIME") << " cannot be used in case of Q*F(Q) fitting" << endl;
				cout << "\t(when the proram was compiled with  _MUL_SCAT_VECTOR compiler option. Ignoring it and continuing...." << endl;
#else
				int index;
				if (!GetInt(ss, *it, &index, "ExptsData::GetExptParamFree")) return ErrLackValue(*it);
				if (index > 0)//active key word
				{
					if (fqfprimeindex[i] & (int)pow(2, (index - 1)))//check that it was not yet set for this type
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "For data set" << i + 1 << " the " << index << ". atom type is already involved in using f' during coefficient calculation" << endl;
						cout << "either by a previously used " << CheckKey("FPRIME-INDEX_FPRIME") << " key word or a " << endl;
						cout << CheckKey("AXS-INDEX_FPRIME_FRACTION") << " keyword! An atom type can either be involved in AXS or using fix f' in a data set!" << endl;
						cout << "Correct the " << datfilename << " file and try again!" << endl;
						CleanExit();
					}
					fqfprimecount[i] += 1;
					fqfprimeindex[i] += (int)pow(2, (index - 1));//indices of the sets with f' will be stored as sum_fprime_itype(2**itype) 
					if (!(ss >> fqfprime[i * ntypes + index - 1])) return ErrLackValue(*it);
				}
#endif
			}
			else if (key == CheckKey("COMPTON"))
			{
#ifdef _MUL_SCAT_VECTOR
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("COMPTON") << " cannot be used in case of Q*F(Q) fitting" << endl;
				cout << "\t(when the proram was compiled with  _MUL_SCAT_VECTOR compiler option. Ignoring it and continuing...." << endl;
#else
				if (!(ss >> fqusecompton[i])) return ErrLackValue(*it);
#endif
			}
			else if (key == CheckKey("AXS-INDEX_FPRIME_FRACTION"))
			{
#ifdef _MUL_SCAT_VECTOR
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("AXS-INDEX_FPRIME_FRACTION") << " cannot be used in case of Q*F(Q) fitting" << endl;
				cout << "\t(when the proram was compiled with  _MUL_SCAT_VECTOR compiler option. Ignoring it and continuing...." << endl;
#else
				if (!GetInt(ss, *it, fqAXS + i, "ExptsData::GetExptParamFree")) return ErrLackValue(*it);//0 for no AXS, index of the type starting with 1 for AXS
				if (fqAXS[i] > 0)
				{
					if (fqfprimeindex[i] & (int)pow(2, (fqAXS[i] - 1)))//check that it was not yet set for this type
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "For data set" << i + 1 << " the " << fqAXS[i] << ". atom type is already involved in using f' during coefficient calculation" << endl;
						cout << "either by a previously used " << CheckKey("AXS-INDEX_FPRIME_FRACTION") << " key word or a " << endl;
						cout << CheckKey("FPRIME-INDEX_FPRIME") << " keyword! An atom type can either be involved in AXS or using fix f' in a data set!" << endl;
						cout << "Correct the " << datfilename << " file and try again!" << endl;
						CleanExit();
					}
					is_AXS = 1;
					fqfprimecount[i] += 1;
					fqfprimeindex[i] +=(int)pow(2, (fqAXS[i] - 1));//indices of the sets with f' will be stored as sum_fprime_itype(2**itype) 
					if (!(ss >> fqfprime[i * ntypes + fqAXS[i] - 1])) return ErrLackValue(*it);
					if (!(ss >> fqfprimefactor[i])) return ErrLackValue(*it);
					
					if (fqfprimefactor[i] < 0)
					{
						cout << "\n****ERROR****" << endl;
						cout << "The fraction to vary the f' for data set " << i + 1 << " is " << fqfprimefactor[i] << endl;
						cout << "Has to be  0<=FRACTION<=1! Correct the " << datfilename << " file and try again!" << endl;
						CleanExit();
					}
				}
#endif
			}
			else if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (sig++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple SIGMA(-MASTER)/SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckKeyValue("TYPE", "XRD") << " constraint!" << endl << "\tThe last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == cumul + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!(ss >> fqsigma[i])) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA") || key == CheckKey("SIGMA-MASTER"))
					fqsigma[i] = fabs(fqsigma[i]);
				if (key == CheckKey("SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1) RunParams::lead_series_ind = cumul + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry was declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						cout << "\tIt will be reset to default " << LEAD_SERIES_IND_DEF << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					fqsigma[i] = -fabs(fqsigma[i]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					fqsigma[i] = fabs(fqsigma[i]);
				}
			}
			else if (key == CheckKey("READ-COEFFS"))
			{
				if (!(ss >> fqreadcoeffs[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "READ-COEFFS for XRD set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(fqreadcoeffs[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("USE-RFACTOR"))
			{
				if (!(ss >> fquseR[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "USE-RFACTOR for XRD set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(fquseR[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("RENORM"))
			{
				if (!(ss >> fqrenorm[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "RENORM for XRD set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(fqrenorm[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("POLY-BACK-FLAGS"))
			{
				int last = -1;
				bool rv = true;
				int initvals[4];
				for (int j = 0; j < 4; j++)
				{
					if (rv) rv = GetInt(ss, *it, initvals + j, "ExptsData::GetExptParamFree");
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n****ERROR****" << endl;
							cout<<"Too few parameters. You should provide at least one value in the following line: " << endl << *it << endl;
							return false;
						}
						else initvals[j] = 0;
					}
					else last = j;
				}
				if ((last != 4 - 1) && ::debug)
				{
					cout << "\nNOTE(" << ++note << "): polynomial background correction flags from " << last + 1 << ". to 3rd order are/is set to zero at the " << endl;
					cout<<"\t "<< i + 1 << ". F(Q) experimental data." << endl;
				}
				use_cubic[cumul] = initvals[3];
				fqoffset[i] = initvals[0];
				fqlinear[i] = initvals[1];
				fqquadratic[i] = initvals[2];
				fqcubic[i] = initvals[3];
				for (int j = 0; j < 4; j++)
				{
					mystrcpy(name, NAME_SIZE, "POLY-BACK-FLAGS ");
					switch (j)
					{
					case 0:
						mystrcat(name, NAME_SIZE, "0th term for XRD set ");
						break;
					case 1:
						mystrcat(name, NAME_SIZE, "1st term for XRD set ");
						break;
					case 2:
						mystrcat(name, NAME_SIZE, "2nd term for XRD set ");
						break;
					case 3:
						mystrcat(name, NAME_SIZE, "3rd term for XRD set ");
						break;
					}
					mystrcat(name, NAME_SIZE, conv_numb);
					Check0_1(initvals[j], name, "ExptsData::GetExptParamFree");
				}
			}
			else if (key == CheckKey("COMPTON-SUPR-COEFF"))
			{

				if (!(ss >> fqrenalpha[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "COMPTON-SUPR-COEFF for XRD set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(fqrenalpha[i], name, "ExptsData::GetExptParamFree");
#ifdef _MUL_SCAT_VECTOR
				if (fqrenalpha[i] == 1)
				{
					cout << "\nWARNING(" << ++warn << "): " << CheckKey("COMPTON-SUPR-COEFF") << " cannot be used in case of Q*S(Q) fitting" << endl;
					cout << "\t(when the proram was compiled with  _MUL_SCAT_VECTOR compiler option. Setting it to 0 and continuing...." << endl;
					fqrenalpha[i] = 0;
				}
#endif
			}
			else
			{
				CheckKeyTag(key, CheckTag("EXP"), CheckKeyValue("TYPE", "XRD"));//write the appropriate error mesage
				return false;
			}

		}
		//Open experimental data to read to clarify fqmin/fqmax,etc.
		retval = FindLimits(datafilename[cumul], fqmin[i], fqmax[i], liminunits[0], liminunits[1]);
		if (!retval)
			return retval;//no reason to continue
		//Check. whether all the custom data was supplied for the given renormalization switch combination, if I(Q) is fitted
		if (!fqreadcoeffs[i])
			is_Xcoeff_calc = 1;//there is a set with coeff calculation
		if (fqfitIQ[i])//I(Q) fitting
		{
			int tempflag = 0;
			if (!fqusecompton[i]  && fqrenalpha[i])
				tempflag=1;
			else
			{
				if (RunParams::chem_symbols != nullptr)
				{
					if ((ntypes == 1 && RunParams::chem_symbols[0].compare("dummy") == 0) && fqrenalpha[i])
						tempflag = 1;
				}
			}
			if (tempflag>0)
			{
				cout << "\nWARNING(" << ++warn << "): For the " << i + 1 << ". X-ray data set alpha renormalization was chosen with the " << CheckKey("COMPTON-SUPR-COEFF") << "=1" << endl;
				if (fqusecompton[i] == 0)
					cout << "\tbut no Compton contribution was chosen with the " << CheckKey("COMPTON") << "=0 key word." << endl;
				else
					cout << "\tdummy atom is used in a one component system." << endl;
				cout << "\tNo Compton contribution and therefore no alpha renormalization will be used!" << endl;
				fqrenalpha[i] = 0;
			}
			
			
			if (fqrenorm[i] && !fqoffset[i] && !fqrenalpha[i])//case 1 'a' renorm, custom value for b and alpha needed
			{
				if (fabs(fqb[i] - NO_DVALUE_DEF) < LOAD_TOL)//no b was given
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "I(Q) fitting was chosen for the x-ray data set " << i + 1 << endl;
					cout << "Although 'a' coefficient renormalization was chosen according to the state of the " << CheckKey("RENORM") << "=1, " << endl;
					cout << CheckKey("POLY-BACK-FLAGS") << " constant=0 and " << CheckKey("COMPTON-SUPR-COEFF") << "=0 switches," << endl;
					cout << "no value for b coefficient was supplied at the second place after the " << CheckKey("I(Q)_A_B_ALPHA") << " key!" << endl;
					cout<< "The correct format for this option is "<< CheckKey("I(Q)_A_B_ALPHA")<<" = 1 1.0 b alpha, where the last two are the desired custom parameters!"<<endl; 
					cout << "Correct the *.dat file and try again! Cannot run this way, exiting..." << endl;
					CleanExit();
				}
				if (fabs(fqalpha[i] - NO_DVALUE_DEF) < LOAD_TOL)//no alpha was given
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "I(Q) fitting was chosen for the x-ray data set " << i + 1 << endl;
					cout << "Although 'a' coefficient renormalization was chosen according to the state of the " << CheckKey("RENORM") << "=1, " << endl;
					cout << CheckKey("POLY-BACK-FLAGS") << " constant=0 and " << CheckKey("COMPTON-SUPR-COEFF") << "=0 switches for x-ray data set,"<< endl;
					cout << "no value for alpha coefficient was supplied at the third place after the " << CheckKey("I(Q)_A_B_ALPHA") << " key!" << endl;
					cout << "The correct format for this option is " << CheckKey("I(Q)_A_B_ALPHA") << " = 1 1.0 b alpha, where the last two are the desired custom parameters!" << endl;
					cout << "Correct the *.dat file and try again! Cannot run this way, exiting..." << endl;
					CleanExit();
				}
			}
			if (fqrenorm[i] && fqoffset[i] && !fqrenalpha[i])//case 3 'a' and b renorm, alpha =0, set it
			{
				fqalpha[i] = 0.0;

			}
			if (!fqrenorm[i] && !fqoffset[i] && fqrenalpha[i])//case 32 alpha renorm, custom value for a and b are needed
			{
				if (fabs(fqa[i] - NO_DVALUE_DEF) < LOAD_TOL || fabs(fqa[i] ) < LOAD_TOL)//no a was given or 0 was given
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "I(Q) fitting was chosen for the x-ray data set " << i + 1 << endl;
					cout << "Although 'alpha' coefficient renormalization was chosen according to the state of the " << CheckKey("RENORM") << "=0, " << endl;
					cout << CheckKey("POLY-BACK-FLAGS") << " constant=0 and " << CheckKey("COMPTON-SUPR-COEFF") << "=1 switches for the x-ray data set," << endl;
					cout << "no value or zero for 'a' coefficient was supplied at the first place after the " << CheckKey("I(Q)_A_B_ALPHA") << " key!" << endl;
					cout << "The correct format for this option is " << CheckKey("I(Q)_A_B_ALPHA") << " = 1 a b, where the last two are the desired custom parameters!" << endl;
					cout << "Parameter 'a' should not be zero!" << endl;
					cout << "Correct the *.dat file and try again! Cannot run this way, exiting..." << endl;
					CleanExit();
				}
				if (fabs(fqb[i] - NO_DVALUE_DEF) < LOAD_TOL)//no b was given
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "I(Q) fitting was chosen for the x-ray data set " << i + 1 << endl;
					cout << "Although 'alpha' coefficient renormalization was chosen according to the state of the " << CheckKey("RENORM") << "=0, " << endl;
					cout << CheckKey("POLY-BACK-FLAGS") << " constant=0 and " << CheckKey("COMPTON-SUPR-COEFF") << "=1 switches," << endl;
					cout << "no value for b coefficient was supplied at the second place after the " << CheckKey("I(Q)_A_B_ALPHA") << " key!" << endl;
					cout << "The correct format for this option is " << CheckKey("I(Q)_A_B_ALPHA") << " = 1 a b, where the last two are the desired custom parameters!" << endl;
					cout << "Correct the *.dat file and try again! Cannot run this way, exiting..." << endl;
					CleanExit();
				}
				is_IQ = 2;//nonlin
			}
			if (fqrenorm[i] && !fqoffset[i] && fqrenalpha[i])//case 33 'a' and alpha renorm, b=0
			{
				fqb[i] = 0.0;
				is_IQ = 2;//nonlin
				
			}
			if (fqrenorm[i] && fqoffset[i] && fqrenalpha[i])//case 35 'a', 'b' and alpha renorm
			{
				is_IQ = 2;//nonlin
			}
			fqlinear[i] = 0;
			fqquadratic[i] = 0;
			fqcubic[i] = 0;
			use_cubic[cumul] = 0;
			if (fqIQbackgcorr[i])
				fqmuact[i] = fqmu[i];
		}
		else
		{
			//no I(Q) fitting for this fq series, put the background correction indicator  to def, and the others accordingly
			//regardless what was given
			fqIQbackgcorr[i] = 0;
			fqAXS[i] = 0;
			fqusecompton[i] = 0;
			fqmu[i] = 0;
			fqdmumax[i] = 0;
			fqmuact[i] = 0;
			fqrenalpha[i] = 0;
			if (fqreadcoeffs[i] && fqfprimeindex[i] > 0)//F(q) fitting, coeffs are read, cannot use f'
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "f' should be used for " << i + 1 << ". X-ray data fitting F(Q)!" << endl;
				cout << "This is not possible, either set " << CheckKey("READ-COEFFS") << " to 0" << " indicating that" << endl;
				cout << "coefficients should be calculated and this way f' can be used, or remove " << CheckKey("FPRIME-INDEX_FPRIME") << endl;
				cout << "from the " << datfilename << " file!" << endl;
				CleanExit();
			}
			
		}
		
		//make sure, that AXS or fixed f' is not used for dummy atom
		for (int itype = 0; itype < ntypes; itype++)
		{
			if (RunParams::chem_symbols != nullptr)
			{
				if (RunParams::chem_symbols[itype].compare("dummy") == 0 && (fqfprimeindex[i] & (int)pow(2, (itype))))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "dummy atom was used after the " << CheckKey("CHEMICAL-SYMBOLS") << " key world in the " << endl;
					cout << "\t" << CheckTag("GENERAL") << " tag for the " << itype + 1 << ". type, and for the " << i + 1 << ".XRD data set" << endl;
					if (fqAXS[i] - 1 == itype)
						cout << "\tAXS should be used";
					else
						cout << "\tfixed f' sould be used";
					cout << " according to the " << CheckTag("EXP") << " section of this data set, which has no sense." << endl;
					cout << "\tdummy atoms does not contribute to  X-ray, neutron or Compton scattering, (the X-ray f(Q)=0, f'=0)" << endl;
					if (fqAXS[i] - 1 == itype)
						cout << "\tso f' cannot be varied!" << endl;
					cout << "Correct the *.dat file!" << endl;
					CleanExit();
				}
			}
		}
		if (fqusecompton[i])
			Compton_chemsymb = true;
		if (retval)
		{
			cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". XRD data set point range " << fqmin[i] << ".." << fqmax[i] << " corresponds to " << liminunits[0] << ".." << liminunits[1] << " in A^-1" << endl;
 			fqused_tot += fqmax[i] - fqmin[i] + 1;
		}
		
	}
	//chemical symbols can be needed for S(Q), F(Q) coeff calculation and Compton calculation. Check, if they are present, if needed
	
	if (is_Xcoeff_calc)//chemical symbols are needed, check if they are present
	{
		if (RunParams::chem_symbols==nullptr)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "X-ray coefficients should be calculated for at least one XRD data set, but the " << CheckKey("CHEMICAL-SYMBOLS") << endl;
			cout << "key word is not present in the " << CheckTag("GENERAL") << " section!" << endl;
			cout << "Correct the *.dat file!" << endl;
			CleanExit();
		}
		for (int ii = 0; ii < ntypes; ii++)
		{
			//check the symbol in the table
			if (Xray_coeffs_W.count(RunParams::chem_symbols[ii]) == 0)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The symbol " << RunParams::chem_symbols[ii] << " is not in the Waasmaier-Kirfel table for X-ray coefficient calculation!" << endl;
				GiveWSymbols(0); 
				CleanExit();
			}
		}
	}
	if (Compton_chemsymb)
	{

		if (RunParams::chem_symbols == nullptr)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Compton scattering should be calculated for at least one XRD data set, but the " << CheckKey("CHEMICAL-SYMBOLS") << endl;
			cout << "key word is not present in the " << CheckTag("GENERAL") << " section!" << endl;
			cout << "Correct the *.dat file!" << endl;
			CleanExit();
		}
		for (int ii = 0; ii < ntypes; ii++)
		{
			//check the symbol in the table
			if (Compton.count(RunParams::chem_symbols[ii]) == 0)
			{
				if (Compton.count(RunParams::chem_symbols_standard[ii]) == 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The symbol " << RunParams::chem_symbols[ii];
					if (RunParams::chem_symbols[ii].compare(RunParams::chem_symbols_standard[ii]) != 0)
						cout << " which belongs to element " << RunParams::chem_symbols_standard[ii];
					cout<< "\nis not among the Balyuzi parameters for the Compton scattering calculation!" << endl;
					GiveCSymbols(0);
					CleanExit();
				}
			}
		}
	}
	if (is_Ncoeff_calc)
	{
		for (int ii = 0; ii < ntypes; ii++)
		{
			//check the symbol in the table
			if (N_scat_length.count(RunParams::chem_symbols_standard[ii]) == 0)//there is no natural abundance symbol, but there might be
			{//isotopes given by ISOTOPE-COUNT_SYMBOLS_RATIOS for each data set
				for (int i = 0; i < nsq; i++)
				{
					//if isotope_count>0, then it was already checked, that athe symbol is in the table and metching with chem symbol
					if (isotope_count[i * ntypes + ii].count == 0 && sqreadcoeffs[i] == 0)//should calculate for natural abundance
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The symbol " << RunParams::chem_symbols_standard[ii] << " indicating the natural abudance is not in the Sears table for neutron scattering coefficient calculation," << endl;
						cout << "which should be used at least for the " << ii + 1 << ". atom type of the " << i + 1 << ". neutron data set!" << endl;
						cout << "Either give the partial coefficients for this data set, or give the isotope composition with" << endl;
						cout << CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS") << " using the symbols existing for this isotope!" << endl;
						cout << "Correct the " << datfilename << " file and try again!" << endl;
						GiveNSymbols(1, RunParams::chem_symbols_standard[ii]);
						CleanExit();
					}
				}
			}
		}
	}
	if (is_AXS)
	{
		int countAXS = 0,countIQ=0;
		for (int i = 0; i < nfq; i++)
		{
			if (fqAXS[i] > 0)
				countAXS++;
			if (fqfitIQ[i])
				countIQ++;
		}
		if (2*countAXS > countIQ)//AXS in itself is not meaningfull, there has to be a normal I(Q) set measured with the same setup as well for each
		{
			cout << "\n*****ERROR*****" << endl;
			if (countAXS == 1)
				cout << "There is ";
			else
				cout << "There are ";
			cout << countAXS << " AXS I(Q) set(s), but only "<<countIQ-countAXS<< " normal (not AXS) I(Q) data set(s)! " << endl;
			cout<<"AXS (which is measured close to the absorption edge of the given element) is only meaningful together with" << endl;
			cout << "a normal I(Q) data set measured futher away form the absorption edge at lower energy!" << endl;
			cout << "So there have to be at least twice as much I(Q) fitting data sets as AXS sets!" << endl;
			cout<<"Either remove the AXS set(S), or give an appropriate number of normal I(Q)-fitting data set as well and start again!" << endl;
			cout << "Cannot run this way exiting..." << endl;
			CleanExit();
		}
	}
	if (is_IQ < 2)
		log_nlr_steps = 0;
	for (int i = 0; i < nfq * ntypes; i++)
		fqfprimeact[i] = fqfprime[i];

	
	for (int i = 0; (i < nfg) && retval; i++, cumul++, it2++)
	{
		sig = 0;
		double liminunits[2] = { -1.0, -1.0 };
		IntToStr(&conv_numb, i + 1);
		fgmin[i] = -1;
		fgmax[i] = -1;
		//Initialize by default values
		fgsub[i] = SUBTRACT_DEF;
		fgrenorm[i] = RENORM_DEF;
		fgoffset[i] = fglinear[i] = fgquadratic[i] = fgcubic[i] = RENORM_DEF; 
		use_cubic[cumul] = USE_CUBIC_DEF; 
		fguseR[i] = USE_RFACTOR_DEF;

		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("TYPE"))
			{
				string own;
				if (!(ss >> own)) return ErrLackValue(*it);
				if (own != "EDIFF")
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Please, check the EXP datasets: only one TYPE entry is allowed for each!" << endl;
					return false;
				}
			}
			else if (key == CheckKey("DATAFILE"))
			{
				params.erase(std::remove(params.begin(), params.end(), ' '), params.end());
				std::copy(params.begin(), params.end(), datafilename[cumul]);
				(datafilename[cumul])[params.length()] = 0;
			}
			else if (key == CheckKey("POINT-RANGE"))
			{
				if (!(GetInt(ss, *it, fgmin + i, "ExptsData::GetExptParamFree") && GetInt(ss, *it, fgmax + i, "ExptsData::GetExptParamFree")))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("G-RANGE"))
			{
				bool ret = true;
				for (int j = 0; (j < 2) && ret; j++) ret = static_cast<bool>(ss >> liminunits[j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("R-SPACING"))
			{
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
				if (!(ss >> RunParams::rspacing[cumul])) return ErrLackValue(*it);
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
			}
			else if (key == CheckKey("CONST-SUBTRACT"))
			{
				if (!(ss >> fgsub[i])) return ErrLackValue(*it);
			}
			else if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (sig++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple SIGMA(-MASTER)/SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckKeyValue("TYPE", "EDIFF") << " constraint!" << endl << "\tThe last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == cumul + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!(ss >> fgsigma[i])) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA") || key == CheckKey("SIGMA-MASTER"))
					fgsigma[i] = fabs(fgsigma[i]);
				if (key == CheckKey("SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1) RunParams::lead_series_ind = cumul + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry was declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout<< " section!" << endl;
						cout << "\tIt will be reset to default " << LEAD_SERIES_IND_DEF << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					fgsigma[i] = -fabs(fgsigma[i]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					fgsigma[i] = fabs(fgsigma[i]);
				}
			}
			else if (key == CheckKey("USE-RFACTOR"))
			{
				if (!(ss >> fguseR[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "USE-RFACTOR for EDIFF set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(fguseR[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("RENORM"))
			{
				if (!(ss >> fgrenorm[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "RENORM for EDIFF set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(fgrenorm[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("POLY-BACK-FLAGS"))
			{
				int last = -1;
				bool rv = true;
				int initvals[4];
				for (int j = 0; j < 4; j++)
				{
					if (rv) rv = GetInt(ss, *it, initvals + j, "ExptsData::GetExptParamFree");
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "Too few parameters. You should provide at least one value in the following line: " << endl << *it << endl;
							return false;
						}
						else initvals[j] = 0;
					}
					else last = j;
				}
				if ((last != 4 - 1) && ::debug)
				{
					cout << "\nNOTE(" << ++note << "): polynomial background correction flags from " << last + 1 << ". to 3rd order are/is set to zero at the " << endl;
					cout<<"\t "<<i + 1 << ". F(Q) experimental data." << endl;
				}
				use_cubic[cumul] = initvals[3];
				fgoffset[i] = initvals[0];
				fglinear[i] = initvals[1];
				fgquadratic[i] = initvals[2];
				fgcubic[i] = initvals[3];
				for (int j = 0; j < 4; j++)
				{
					mystrcpy(name, NAME_SIZE, "POLY-BACK-FLAGS ");
					switch (j)
					{
					case 0:
						mystrcat(name, NAME_SIZE, "0th term for EDIFF set ");
						break;
					case 1:
						mystrcat(name, NAME_SIZE, "1st term for EDIFF set ");
						break;
					case 2:
						mystrcat(name, NAME_SIZE, "2nd term for EDIFF set ");
						break;
					case 3:
						mystrcat(name, NAME_SIZE, "3rd term for EDIFF set ");
						break;
					}
					mystrcat(name, NAME_SIZE, conv_numb);
					Check0_1(initvals[j], name, "ExptsData::GetExptParamFree");
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("EXP"),CheckKeyValue("TYPE","EDIFF"));//write the appropriate error mesage
				return false;
			}

		}
		//Open experimental data to read to clarify fgmin/fgmax,etc.
		retval = FindLimits(datafilename[cumul], fgmin[i], fgmax[i], liminunits[0], liminunits[1] );
		if (retval)
		{
			cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". EDIFF data set point range " << fgmin[i] << ".." << fgmax[i] << " corresponds to " << liminunits[0] << ".." << liminunits[1] << " in A^-1" << endl;
			fgused_tot += fgmax[i] - fgmin[i] + 1;
		}
		else
			return retval;//no reason to continue
	}

	
	for (int i = 0; (i < nek) && retval; i++, cumul++, it2++)
	{
		sig = 0;
		double liminunits[2] = { -1.0, -1.0 };
		IntToStr(&conv_numb, i + 1);
		std::vector<double> rmin(ntypes, -1.0);//this is just dummy, as here the binshift was not set, so no rmin rmax can be calculated
		std::vector<double> rmax(ntypes, -1.0);
		for (int j = 0; j < ntypes; j++)
		{
			ek_rmin[i*ntypes + j] = -1;
			ek_rmax[i*ntypes + j] = -1;
		}
		ekmin[i] = -1;
		ekmax[i] = -1;
		//Initialize by default values
		ekchipower[i] = EXAFS_CHI2_POWER_DEF;
		ekrenorm[i] = RENORM_DEF;
		ekoffset[i] = RENORM_DEF;
		ekuseR[i] = USE_RFACTOR_DEF;
		ek_dE0[i] = MAX_E0_SHIFT_DEF;
		ek_ngrid_in[i] = EXAFS_NGRID_DEF; 
		ek_gridstart[i] = EXAFS_NGRID_DEF;
		ek_gridfrom[i] = EXAFS_NGRID_DEF;
		ek_gridto[i] = EXAFS_NGRID_DEF;
		bool gridflag=false;//E0GRID-FROM_TO_START was not read

		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{
			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);
			if (key == CheckKey("TYPE"))
			{
				string own;
				if (!(ss >> own)) return ErrLackValue(*it);
				if (own != "EXAFS")
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Please, check the EXP datasets: only one TYPE entry is allowed for each!" << endl;
					return false;
				}
			}
			else if (key == CheckKey("DATAFILE"))
			{
				//for each EXAFS series first the k-E(k) file name is stored, the the coefficient file
				params.erase(std::remove(params.begin(), params.end(), ' '), params.end());
				std::copy(params.begin(), params.end(), datafilename[cumul + i]);
				(datafilename[cumul + i])[params.length()] = 0;
			}
			else if (key == CheckKey("POINT-RANGE"))
			{
				if (!(GetInt(ss, *it, ekmin + i, "ExptsData::GetExptParamFree") && GetInt(ss, *it, ekmax + i, "ExptsData::GetExptParamFree")))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("K-RANGE"))
			{
				bool ret = true;
				for (int j = 0; (j < 2) && ret; j++) ret = static_cast<bool>(ss >> liminunits[j]);
				if (!ret)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide 2 limiting values in the following line: " << endl << *it << endl;
					return false;
				}
			}
			else if (key == CheckKey("R-SPACING"))
			{
#ifndef _NO_PERIODIC
#ifndef _VIBR_AMP
				if (!(ss >> RunParams::rspacing[cumul])) return ErrLackValue(*it);
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
#else
				cout << "\nWARNING(" << ++warn << "): " << CheckKey("R-SPACING") << " cannot be used in case of no periodic boundary conditions or vibramp, it is ingnored!" << endl;
#endif
			}
			else if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (sig++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple SIGMA(-MASTER)/SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckKeyValue("TYPE", "EXAFS") << " constraint!" << endl << "\tThe last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == cumul + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!(ss >> eksigma[i])) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA") || key == CheckKey("SIGMA-MASTER"))
					eksigma[i] = fabs(eksigma[i]);
				if (key == CheckKey("SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1) RunParams::lead_series_ind = cumul + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry was declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						cout << "\tIt will be reset to default " << LEAD_SERIES_IND_DEF << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					eksigma[i] = -fabs(eksigma[i]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					eksigma[i] = fabs(eksigma[i]);
				}
			}
			else if (key == CheckKey("USE-RFACTOR"))
			{
				if (!(ss >> ekuseR[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "USE-RFACTOR for EXAFS set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(ekuseR[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("DELTAE0_NGRID"))
			{
				if (!(ss >> ek_dE0[i]) || !(ss >> ek_ngrid_in[i]))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide a positive real value for dE0 in eV and the number of grid points in one direction in the following line: " << endl << *it << endl;
					return false;
				}
				if (ek_ngrid_in[i]>0 && fabs(ek_dE0[i])>TOLERANCE2)//only if valid data is given
					is_E0shift = 1;//
			}
			else if (key == CheckKey("E0GRID-FROM_TO_START"))
			{
				if (!(ss >> ek_gridfrom[i]) || !(ss >> ek_gridto[i]))
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "Too few parameters. You should provide two integers between -NGRID -> +NGRID for the grid range to use," << endl;
					cout<<"and optionally a third for the starting grid index in line : " << endl << *it << endl;
					return false;
				}
				if (RunParams::continuation == 0)//if continuation, it will be set to the grid index of the previous run read form state
				{
					if (!(ss >> ek_gridstart[i]))//not needed for exact continuation, dealt with later
					{
						ek_gridstart[i] = ek_gridfrom[i] + (int)(ek_gridto[i] - ek_gridfrom[i]) / 2;
						cout << "\nWARNING(" << ++warn << "): Starting grid index was set to the (left)-middle of the range of FROM: " << ek_gridfrom[i] << " TO: " << ek_gridto[i] << endl;
						cout << "\tto grid index " << ek_gridstart[i] << endl;
					}
				}
				gridflag = true;//parameters were read, goes what it says, if no conflict

			}
			else if (key == CheckKey("RENORM"))
			{
				if (!(ss >> ekrenorm[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "RENORM for EXAFS set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(ekrenorm[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("POLY-BACK-FLAGS"))
			{
				int last = -1;
				bool rv = true;
				int initvals[1];
				for (int j = 0; j < 1; j++)
				{
					if (rv) rv = GetInt(ss, *it, initvals + j, "ExptsData::GetExptParamFree");
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n*****ERROR*****" << endl; 
							cout<<"Too few parameters. You should provide one value in the following line: " << endl << *it << endl;
							return false;
						}
						else initvals[j] = 0;
					}
					else last = j;
				}
				ekoffset[i] = initvals[0];
				mystrcpy(name, NAME_SIZE, "POLY-BACK-FLAGS 0th term for EXAFS set ");
				mystrcat(name, NAME_SIZE, conv_numb);
				Check0_1(ekoffset[i], name, "ExptsData::GetExptParamFree");
			}
			else if (key == CheckKey("ABSORBER-TYPE"))
			{
				if (!(ss >> ek_abstype[i])) return ErrLackValue(*it);
				ek_abstype[i]--;
			}
			else if (key == CheckKey("CHIK-POWER"))
			{
				if (!(ss >> ekchipower[i])) return ErrLackValue(*it);
			}
			else if (key == CheckKey("BACKSCATT-FILE"))
			{
				//for each EXAFS series first the k-E(k) file name is stored, the the coefficient file
				params.erase(std::remove(params.begin(), params.end(), ' '), params.end());
				std::copy(params.begin(), params.end(), datafilename[cumul + i + 1]);
				(datafilename[cumul + i + 1])[params.length()] = 0;
			}
			else if (key == CheckKey("R-POINT-RANGE"))
			{
				int last = -1;
				bool rv = true;
				for (int j = 0; j < ntypes; j++)
				{
					if (rv) rv = (GetInt(ss, *it, ek_rmin + i * ntypes + j, "ExptsData::GetExptParamFree") && GetInt(ss, *it, ek_rmax + i * ntypes + j, "ExptsData::GetExptParamFree"));
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "Too few parameters. You should provide at least two value in the following line: " << endl << *it << endl;
							return false;
						}
						else
						{
							ek_rmin[i*ntypes + j] = ek_rmin[last];
							ek_rmax[i*ntypes + j] = ek_rmax[last];
						}
					}
					else last = j;
				}
				if (last != ntypes - 1 && ::debug)
				{
					cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". EXAFS E(k,r) data file "<<CheckKey("R-POINT-RANGE")<<" for type";
					if (last < ntypes - 2)
						cout << "s ";
					else
						cout << " ";
					for (int ii = 0; ii < ntypes - 1 - last; ii++)
					{
						cout << last + ii + 2;
						if ((last + ii+1) != ntypes - 1)
							cout << ", ";
					}
					cout << " is set to " << ek_rmin[last] << ".." << ek_rmax[last] << endl;
				}
			}
			else if (key == CheckKey("R-RANGE"))
			{
				int last = -1;
				bool rv = true;
				for (int j = 0; j < ntypes; j++)
				{
					if (rv) rv = static_cast<bool>(ss >> rmin[j] >> rmax[j]);
					if (!rv)
					{
						if (last == -1)
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "Too few parameters. You should provide at least two values in the following line: " << endl << *it << endl;
							return false;
						}
						else
						{
							rmin[j] = rmin[last];
							rmax[j] = rmax[last];
						}
					}
					else last = j;
				}
				if (last != ntypes - 1 && debug)
				{
					cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". EXAFS E(k,r) data file " << CheckKey("R-RANGE") << " for type";
					if (last < ntypes - 2)
						cout << "s ";
					else
						cout << " ";
					for (int ii = 0; ii < ntypes - 1 - last; ii++)
					{
						cout << last + ii + 2;
						if ((last + ii + 1) != ntypes - 1)
							cout << ", ";
					}
					cout << " is set to " << rmin[last] << ".." << rmax[last] << " A" << endl;
				}

			}
			else
			{
				CheckKeyTag(key, CheckTag("EXP"), CheckKeyValue("TYPE", "EXAFS"));//write the appropriate error mesage
				return false;
			}

		}
		
		if (gridflag)//grid range was set, check conflict, DELTAE0_NGRID also exist, it was cheked in RunParams::ReadFreeExpt
		{
			if (ek_gridto[i] < ek_gridfrom[i])//swap them
			{
				int temp = ek_gridto[i];
				ek_gridto[i] = ek_gridfrom[i];
				ek_gridfrom[i] = temp;
			}
			if (ek_gridfrom[i]<-ek_ngrid_in[i] || ek_gridfrom[i] > ek_ngrid_in[i])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The first used grid index 'FROM' = " << ek_gridfrom[i] << " in " << CheckKey("E0GRID-FROM_TO_START") << " is outside the acceptabe range" << endl;
				cout << -ek_ngrid_in[i] << " -> " << ek_ngrid_in[i] << " for the "<<i+1<<". EXAFS data set!" << endl;
				return false;
			}
			if (ek_gridto[i]<-ek_ngrid_in[i] || ek_gridto[i] > ek_ngrid_in[i])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The last used grid index 'TO' = " << ek_gridto[i] << " in " << CheckKey("E0GRID-FROM_TO_START") << " is outside the acceptabe range" << endl;
				cout << -ek_ngrid_in[i] << " -> " << ek_ngrid_in[i] << " for the " << i + 1 << ". EXAFS data set!" << endl;
				return false;
			}
			
			if (RunParams::continuation==0 && (ek_gridstart[i]<ek_gridfrom[i] || ek_gridstart[i] > ek_gridto[i]))//in case of continuation start grid will be the finishing grid index
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The starting grid index 'START' = " << ek_gridstart[i] << " in " << CheckKey("E0GRID-FROM_TO_START") << " is outside the acceptabe range" << endl;
				cout << ek_gridfrom[i] << " -> " << ek_gridto[i] << " for the " << i + 1 << ". EXAFS data set!" << endl;
				return false;
			}
		}
		else
		{
			ek_gridfrom[i] = -ek_ngrid_in[i];//only important, if DELTAE0_NGRID was specified
			ek_gridto[i] = ek_ngrid_in[i];
			ek_gridstart[i] = 0;//no shift, it will be reset, if continuation
		}
		//Open experimental data to read to clarify ekmin/ekmax,etc.
		//for each EXAFS series first the k-E(k) file name is stored, the the coefficient file
		retval = FindLimits(datafilename[cumul + i], ekmin[i], ekmax[i], liminunits[0], liminunits[1]);
		if (retval)
		{
			cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". EXAFS data set point range " << ekmin[i] << ".." << ekmax[i] << " corresponds to " << liminunits[0] << ".." << liminunits[1] << " in A^-1" << endl;
			ekused_tot += ekmax[i] - ekmin[i] + 1;
		}
		else
			return retval;//no reason to continue
		//Open backscattering datafile to read to clarify rmin/rmax,etc.
		//for each EXAFS series first the k-E(k) file name is stored, the the coefficient file
		retval = FindLimitsExafsBackscattering(cumul,datafilename[cumul + i + 1], &ek_rmin[i*ntypes], &ek_rmax[i*ntypes], &rmin[0], &rmax[0], ntypes);
		if (retval)
		{
			for (int j = 0; j < ntypes; j++)
			{
				cout << "\nNOTE(" << ++note << "): for the " << i + 1 << ". EXAFS E(k,r) data file for the " << j + 1 << ". partial point range " << ek_rmin[i*ntypes + j] << ".." << ek_rmax[i*ntypes + j]  << endl;
				ek_rused[i*ntypes + j] = (ek_rmax[i*ntypes + j] - ek_rmin[i*ntypes + j] + 1);
				ek_rused_tot += ek_rused[i*ntypes + j];
			}
		}
		else
			return retval;//no reason to continue
	}
	if (conv_numb != NULL)
		delete[] conv_numb;
	delete[] name;
	return retval;
};


bool ExptsData::FindLimitsExafsBackscattering(int expno, const char *filename, int *minlim, int *maxlim, double *minval, double *maxval, const int dimension)
//Looks for the actual limit values given in filename: can be point range(min/maxlim), datarange (min/maxval) or using the whole data. If not mentioned, it should initialize with -1/-1.0
//expno is the index of this data set among all the sets
{
	ifstream dfile(filename);
	for (int i = 0; i < dimension; i++)
	{
		if (minlim[i] > maxlim[i]) { int temp = maxlim[i]; maxlim[i] = minlim[i]; minlim[i] = temp; }
		if (minval[i] > maxval[i]) { double temp = (int)maxval[i]; maxval[i] = minval[i]; minval[i] = temp; }
	}
	if (!dfile)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Occured, when opening of the \"" << filename << "\" file has been attempted." << endl;
		return false;
	}
	int nooflines;
	int dim;
	string l;
	getline(dfile, l);
	std::transform(l.begin(), l.end(), l.begin(), [](char ch) {if (ch == ',') ch = ' '; return ch; });
	stringstream ss2(l);
	if (!(ss2 >> dim >> nooflines))
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Occured at the first line of \"" << filename << "\" file." << endl;
		return false;
	}
	if (dim != dimension)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The expected number of partials (" << dimension << ") is inconsistent with that one provided (" << dim << ") in the \"" << filename << "\"" << endl;
		return false;
	}
	for (int i = 0; i < dim; i++)
	{
		string line;
		getline(dfile, line);//Skip the comment line
		getline(dfile, line);//Get the first real datasets and check, how many values does it have
		std::transform(line.begin(), line.end(), line.begin(), [](char ch) {if (ch == ',') ch = ' '; return ch; });
		stringstream ss(line);
		double value;
		int chkdata = 0;
		while (ss >> value) chkdata++;
		if ((minlim[i] < 0) && (maxlim[i] < 0) && (minval[i] < 0.0) && (maxval[i] < 0.0))
		{
			//Not special, nor absolute units have been given: set the default values
			minlim[i] = 1;
			maxlim[i] = chkdata;
		}
		else if ((minlim[i] < 0) && (maxlim[i] < 0))
		{
			//when called from ReadFreeExp custom binsift was not set yet!
			minlim[i] = (int)(minval[i] / RunParams::rspacing[expno] + 0.999 - RunParams::binshift[expno]);
			maxlim[i] = (int)(maxval[i] / RunParams::rspacing[expno] + 0.999 - RunParams::binshift[expno]);
		}
		minval[i] = RunParams::rspacing[expno]*(RunParams::binshift[expno] + minlim[i] - 0.5);
		maxval[i] = RunParams::rspacing[expno]*(RunParams::binshift[expno] + maxlim[i] - 0.5);
		for (int j = 0; j < nooflines - 1; j++)
		{
			if (!getline(dfile, line))
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "Occured at the " << j + 4 + i * (nooflines + 1) << ". line of the \"" << filename << "\" file." << endl;
				return false;
			}
		}
	}


	dfile.close();
	return true;
}

//Looks for the actual limit values given in the file. Before this the desired limits may have been set, either point range(min/maxlim), or datarange (min/maxval).
//If not set yet, then min/maxlim =-1, min/maxval=-1.0
bool ExptsData::FindLimits(const char *filename, int &minlim, int &maxlim, double &minval, double &maxval)
{
	double actual,first=0;
	bool point_range_set = false;
	ifstream dfile(filename);
	if (minlim > maxlim) { int temp = maxlim; maxlim = minlim; minlim = temp; }
	if (minval > maxval) { double temp = (int)maxval; maxval = minval; minval = temp; }
	if (!dfile)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The \"" << filename << "\" file cannot be opened!" << endl;
		return false;
	}
	

	int nooflines;
	if (!(dfile >> nooflines))
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Occured at the first line of \"" << filename << "\" file." << endl;
		return false;
	}
	if ((minlim < 0) && (maxlim < 0) && (minval < 0.0) && (maxval < 0.0))
	{
		//Neither special, nor absolute units have been given: set the default values
		minlim = 1;
		maxlim = nooflines;
	}
	else
	{
		if (minlim > nooflines)
		{
			cout << "\nWARNING(" << ++warn << "): the first data point to use is larger than the number of data points in the " << filename << " file, which is " << nooflines << "." << endl;
			cout << "\tIt will be reset from " << minlim << " to 1!" << endl;
			minlim = 1;
		}
		if (maxlim>nooflines)
		{
			cout << "WARNING(" << ++warn << "): the last data point to use is larger than the number of data points in the " << filename << " file, which is "<<nooflines <<"."<< endl;
			cout << "\tIt will be reset from " << maxlim << " to " << nooflines << "!"<<endl;
			maxlim = nooflines;
		}
	}
	if ((minlim > 0) && (maxlim > 0))//point range was set
		point_range_set = true;
	double prev = -1.0;
	//j -2 will be the serial number of the data line in the file beginning with 0, but j-1 can be compared to minlim and maxlim, as they begin with 1
	for (int j = 0; j < nooflines + 2; j++)
	{
		string line;
		if (!getline(dfile, line))
		{
			if (j < 2)//there are no data points at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "There are no data points in the "<<filename<<" file"<<endl;
				return false;
			}
			else
			{
				cout << "\nWARNING(" << ++warn << "): The number of data points are less, than given in the first line of the " << filename << " file!" << endl;
				cout << "\tOnly " << j - 2 << " data points were found, although " << nooflines << " points were expected. Setting the number of points to " << j - 2 << "!" << endl;
				cout << "\tProceeding..." << endl;
				nooflines = j - 2;
				dfile.clear(ios::goodbit);
				dfile.close();
				//j-2 is the last existing datapoint beginning with 1
				if (j - 2 <= minlim)
				{
					minlim = 1;
					minval = first;
				}
				if (j - 2 <= maxlim)
				{
					maxlim = j - 2;
					maxval = actual;
				}
				if ((minlim > 0) && (maxlim > 0) && (minval > 0.0) && (maxval > 0.0))//all was set
					return true;
				else
					break;
			}
		}
		stringstream ss(line);
		if (j > 1)
		{
			if (!(ss >> actual))
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "While reading of the " << j + 1 << ". line of the \"" << filename << "\" file." << endl;
				return false;				
			}
			if (j == 2)//first data line
				first = actual;//keep it, in case of minval or maxval out of range
			if (point_range_set)//point range set
			{
				if (j - 1 == minlim)
					minval = actual;//this will be only needed for display
				else
				{
					if (j - 1 == maxlim)
						maxval = actual;//this will be only needed for display
				}
			}
			else
			{
				//value range was set
				if ((minval <= actual) && (j == 2))//value range minimum was smaller than the smallest value in the file, reset
				{
					minval = actual;
					minlim = 1;
				}
				else if ((minval <= actual) && (minval > prev))
				{
					if ((actual - minval) > (minval - prev))
					{
						minval = prev;
						minlim = j - 2;
					}
					else
					{
						minval = actual;
						minlim = j - 1;
					}
				}
					
				if ((maxval >= actual) && (j == nooflines + 1))//last data line in the file
				{
					maxval = actual;
					maxlim = nooflines;
				}
				else if ((maxval <= actual) && (maxval > prev) && j!=2)//not the first line
				{
					if ((actual - maxval) > (maxval - prev))
					{
						maxval = prev;
						maxlim = j - 2;
					}
					else
					{
						maxval = actual;
						maxlim = j - 1;
					}
				}
				prev = actual;
			}
		

		}
	}
	if (minlim < 0)//minval was out of range, minlim not set
	{
		minlim = 1;
		minval = first;
	}
	if (maxlim < 0)//maxval was out of range, maxlim not set
	{
		maxlim = nooflines;
		maxval = actual;
	}
	
	dfile.close();
	return true;
}

//set the params for the E0 shift for EXAFS
void ExptsData::SetE0ShiftParams()
{
	int i,j,jj,ig,itype,ik,ir;
	int my_ig;
	int sum,size,sizef,sizek;
	int *first_k, *last_k;//to determine, which k data points can be first and last to use to make linear interpolation possible 
	double my_k;
	double *pk,*pe,*pc_s,*pc_l,*pdkl, *pdkr,*pc_new;
	double k_max_pow_chiweight;//maximum k value for ran E(k) data set
	bool interpolate, keep_k_ori=false;

	if (nek == 0)
		return;
	SetArraysize(&ek_max_dkl, ekused_tot, "ek_max_dkl", "ExptsData::SetE0ShiftParams");
	SetArraysize(&ek_max_dkr, ekused_tot, "ek_max_dkr", "ExptsData::SetE0ShiftParams");
	
	SetArraysize(&first_k, nek, "first_k", "ExptsData::SetE0ShiftParams");
	SetArraysize(&last_k, nek, "last_k", "ExptsData::SetE0ShiftParams");
	sum = 0;//to calculate the number of used points taking into account the lost points due to interpolation
	size = 0;//the size of the ek_coeffs array, each data set can have different number of grid points!
	sizef = 0;// total number of gridpoints for all the data sets
	sizek = 0;//total number of k points
	ek_cf_cum[0] = 0;
	ek_ngrid_cum[0] = 0; 
	is_E0shift = 0;//decide  again, whether there will be really an E0 shift
	for (i = 0; i < nek; i++)
	{
		ekmin_ori[i] = ekmin[i];
		ekmax_ori[i] = ekmax[i];
		if (ek_ngrid_in[i] < 0)
		{
			cout << "\nWARNING(" << ++warn << "): The number of grid points (" << ek_ngrid_in[i] << ") given by " << CheckKey("DELTAE0_NGRID") << " cannot be negative, " << endl;
			cout << "\tit will be set to 0 and no E0 shift will be performed for the " << i + 1 << ". EXAFS data set!" << endl;
			ek_ngrid_in[i] = 0;
		}
		if (ek_dE0[i] < 0)
		{
			if (fabs(ek_dE0[i] - MAX_E0_SHIFT_DEF) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): The dE0 " << ek_dE0[i] << " given by " << CheckKey("DELTAE0_NGRID") << " cannot be negative, " << endl;
				cout << "\tit will be set to 0 and no E0 shift will be performed for the " << i + 1 << ". EXAFS data set!" << endl;
			}
			ek_dE0[i] = 0.0;
		}
		if (ek_dE0[i] > pow(*ek_kfinder_ori[i],2)/ K_SHIFT_CONST)//dE0 cannot be larger than kmin^2*h_bar^2/2/m_e
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The dE0 (" << ek_dE0[i] << ") given by " << CheckKey("DELTAE0_NGRID") << " for the " << i + 1 << ". EXAFS data set cannot be larger " << endl;
			cout << "than kmin^2*h_bar^2/2/m_e= "<< pow(*ek_kfinder_ori[i], 2) / K_SHIFT_CONST<<" !" << endl;
			cout << "For this dE0 value the first used k value should be larger than " << sqrt(ek_dE0[i] * K_SHIFT_CONST) << " 1/A"<<endl;
			cout << "Either change the dE0 value or the first used k data point, correct the " << datfilename << " file and try again!" << endl;
			cout << "Cannot run this way, exiting..." << endl;
			CleanExit();
			
		}
		if (ek_ngrid_in[i] > 0)
			is_E0shift = 1;
		first_k[i] = 0;
		last_k[i] = ekused[i] - 1;
		if (ek_ngrid_in[i] > 0)//there is E0_shift for this data set
		{
			//Calculte the left and right maximum k-shifts, k1, first itt will be calculated for all the k points, later the number of used
			//points and the connected arrays have to be reset to make sure that there are points on each side to interpolate from
			//all the 2*ngrid_in+1 grid points will be used here regardless how many will actually be used, to have the same number of k points, if
			//the deltaE0 and the number of grid points is the same

			pk = ek_kfinder_ori[i];
			pdkl = ek_max_dkl + ekused_cum[i];
			pdkr = ek_max_dkr + ekused_cum[i];
			for (j = 0; j < ekused[i]; j++)
			{
				*pdkl++ = *pk - sqrt(*pk * *pk - K_SHIFT_CONST * ek_dE0[i]);
				*pdkr++ = -*pk + sqrt(*pk * *pk + K_SHIFT_CONST * ek_dE0[i]);
				pk++;
			}
			//find the first and last point considering the interpolation
			//first k' shoud be larger than ori k_min, last should be smaller than ori k_max
			pk = ek_kfinder_ori[i];
			pdkl = ek_max_dkl + ekused_cum[i];
			pdkr = ek_max_dkr + ekused_cum[i];
			
			jj = 0;
			while (*(pk + jj) - *(pdkl+jj) < *ek_kfinder_ori[i])//compare to first point
				jj++;
			if (jj > first_k[i])
				first_k[i] = jj;
			jj = ekused[i] - 1;
			while (*(pk + jj) + *(pdkr+jj) > *(ek_kfinder_ori[i] + ekused[i] - 1))//compare to last point
				jj--;
			if (jj < last_k[i])
				last_k[i] = jj;
		}
		sum += last_k[i] - first_k[i] + 1;
		//calculate the number of gridpoints
		ek_ngrid[i] = ek_gridto[i] - ek_gridfrom[i] + 1;//ek_ngrid_in[i] * 2 + 1; //this will be the number of used data points
		sizek += (last_k[i] - first_k[i] + 1) * ek_ngrid[i];
		for (j = 0; j < ntypes; j++)
			size += (last_k[i] - first_k[i] + 1) * ek_rused[i * ntypes + j] * ek_ngrid[i];
		sizef += ek_ngrid[i];
		ek_cf_cum[i + 1] = ek_cf_cum[i] + ntypes * ek_ngrid[i];
		ek_ngrid_cum[i+1]=ek_ngrid_cum[i]+ ek_ngrid[i];

	}
	
	
	//calculate the grid coefficient matrix value, this will be used, even if there is no E0 shift
	SetArraysize(&ek_coeffs, size, "ek_coeffs", "ExptsData::SetE0ShiftParams");
	//ek_coeffs contains the data in the order of cycles: data set, grid point, ntypes,k,r
	SetArraysize(&ek_cfinder, sizef*ntypes, "ek_cfinder", "ExptsData::SetE0ShiftParams");
	//will point to each individual matrices
	
	size = 0;
	jj = 0;
	for (i = 0; i < nek; i++)
	{
		for (ig = 0; ig < ek_ngrid[i]; ig++)//only the used points
		{
			for (itype = 0; itype < ntypes; itype++)
			{
				ek_cfinder[jj] = ek_coeffs + size;//setting the finder for the k dependent coefficients
				size += (last_k[i] - first_k[i] + 1) * ek_rused[i * ntypes + itype];
				jj++;
			}
		}
	}

	
	for (i = 0; i < nek; i++)
	{
		if (ek_ngrid_in[i] > 0)
		{
			pk = ek_kfinder_ori[i];
			if (ek_gridfrom[i] * ek_gridto[i] > 0)//the used range does not include the original k-grid
				keep_k_ori = true;//the original ek_kvalues_ori array has to be kept, do not delete
			for (ig = 0; ig < ek_ngrid[i]; ig++)
			{
				my_ig = ig + ek_gridfrom[i];//goes from ek_gridfrom -> ek_gridto
				if (my_ig >= ek_gridfrom[i] && my_ig <= ek_gridto[i])//calculate only for the used grid points specified by the FROM TO range
				{
					for (itype = 0; itype < ntypes; itype++)
					{
						for (ik = first_k[i]; ik <= last_k[i]; ik++)//here go through all the data points remained
						{

							my_k = pk[ik];//the new k value
							if (my_ig != 0)
								my_k = sqrt(pow(pk[ik], 2) + K_SHIFT_CONST * my_ig * ek_dE0[i] / ek_ngrid_in[i]);
							if (my_k < pk[0])//smaller than the first k value, cannot use tis point
								cout << "bad small" << endl;

							j = 0;
							while (j<ekused[i] && my_k > pk[j])//finding the two neighbouring points for iterpolation
								j++;
							if (j >= ekused[i])//the jth point does not exist, cannot interpolate
								cout << "bad" << endl;
							//now my_k is between k[j-1] and k[j]
							interpolate = true;
							if (fabs(my_k - pk[j]) < TOLERANCE)//my_ik is the same as an original, no need to interpolate
								interpolate = false;
							pc_s = ek_cfinder_ori[i * ntypes + itype] + (j - 1) * ek_rused[i * ntypes + itype];//beginning of the coeffs belonging to the j-1. k value 
							pc_l = ek_cfinder_ori[i * ntypes + itype] + j * ek_rused[i * ntypes + itype];//beginning of the coeffs belonging to the j. k value 
							pc_new = ek_cfinder[ek_cf_cum[i] + ig * ntypes + itype] + (ik - first_k[i]) * ek_rused[i * ntypes + itype];//beginning of the k-line to set
							for (ir = 0; ir < ek_rused[i * ntypes + itype]; ir++)
							{
								if (interpolate)
									pc_new[ir] = (my_k - pk[j - 1]) / (pk[j] - pk[j - 1]) * (pc_l[ir] - pc_s[ir]) + pc_s[ir];
								else
									pc_new[ir] = pc_l[ir];

							}

						}
					}
				}
			}
		}
		else
		{
			//this is a set with no E0 shift among sets, where there is a shift, here only copy the coeffs
			for (itype = 0; itype < ntypes; itype++)
			{
				pc_s = ek_cfinder_ori[i * ntypes + itype];//origiginal
				pc_new = ek_cfinder[ek_cf_cum[i] + itype];//copied

				for (ik = 0; ik < ekused[i]; ik++)//for each k data point
				{
					for (ir = 0; ir < ek_rused[i * ntypes + itype]; ir++)//for each r point
						*pc_new++ = *pc_s++;
				}
			}
		}
	}
	
	//sum is the new total number of data points
	SetArraysize(&ek_kvalues, sizek, "ek_kvalues", "ExptsData::SetE0ShiftParams");
	SetArraysize(&ek_evalues, sum, "ek_evalues", "ExptsData::SetE0ShiftParams");
	SetArraysize(&ek_kfinder, sizef, "ek_kfinder", "ExptsData::SetE0ShiftParams");
	SetArraysize(&ek_efinder, nek, "ek_efinder", "ExptsData::SetE0ShiftParams");
	
	//now all the data series were exemined, the used data points and the connected arrays have to be contracted
	pe = ek_evalues;
	for (i = 0; i < nek; i++)
	{

		//divinding the E(k) with k_max^ekchipower[i]
		k_max_pow_chiweight = pow(*(ek_kfinder_ori[i] + last_k[i]), (double)ekchipower[i]);
		for (ik = 0; ik < last_k[i] - first_k[i] + 1; ik++)
			*pe++ = *(ek_efinder_ori[i] + first_k[i] + ik) / k_max_pow_chiweight;

	}
	if (sum < ekused_tot)//there is data point loss due to E0 sift
	{
		sizek = 0;
		jj = 0;
		
		for (i = 0; i < nek; i++)
		{
			for (ig = 0; ig < ek_ngrid[i]; ig++)//used grid points
			{
				my_ig = ig + ek_gridfrom[i];//FROM TO
				ek_kfinder[jj] = ek_kvalues + sizek;
				sizek += (last_k[i] - first_k[i] + 1);
				pk = ek_kfinder[jj];

				if (my_ig < 0)
					pdkl = ek_max_dkl + ekused_cum[i] + first_k[i];
				else
					pdkr = ek_max_dkr + ekused_cum[i] + first_k[i];
				for (ik = 0; ik < last_k[i] - first_k[i] + 1; ik++)
				{
					pk[ik] = *(ek_kfinder_ori[i] + first_k[i] + ik);
					if (my_ig != 0)
						pk[ik] = sqrt(pow(pk[ik], 2) + K_SHIFT_CONST * my_ig * ek_dE0[i] / ek_ngrid_in[i]);
				}
				jj++;
			}
		}
		
		//finally resetting the variables related to the number of used points
		ekused_tot = sum;
		ekused_cum[0] = 0;
		jj = 0;
		double *kori_finder;
		for (i = 0; i < nek; i++)
		{
			if (ek_gridfrom[i] * ek_gridto[i] > 0)//range does not include the original
				kori_finder = ek_kfinder_ori[i] + first_k[i];//the original array remains, only point to the first used point
			else
				kori_finder = ek_kfinder[jj - ek_gridfrom[i]];//this will point to the original k values for each data set
			if (first_k[i]!=0 || ekmax[i] != ekmin[i] + last_k[i])//it was reset
			{
				cout << "\nNOTE(" << ++note << "): The " << i + 1 << ". EXAFS DATAFILE point range " << ekmin[i] << ".." << ekmax[i] << " was reset due to the usage of E0 shift to " << endl;
				cout<< "\t"<<ekmin[i]+first_k[i]<<".."<<ekmin[i]+last_k[i]<<" corresponding to " << *kori_finder << ".." << *(kori_finder+ last_k[i] - first_k[i]) << " in A ^ -1" << endl;
			}
			ekused_tot += ekmax[i] - ekmin[i] + 1;
			ekused[i] = last_k[i] - first_k[i] + 1;
			ekused_cum[i+1] = ekused_cum[i] + ekused[i];
			ekmax[i] = ekmin[i]+last_k[i];
			ekmin[i] += first_k[i];
			ek_efinder[i] = ek_evalues + ekused_cum[i];
			ek_kfinder_ori[i] = kori_finder;
			jj += ek_ngrid[i];
		}
	
	}
	else
	{
		//only copy the values to the new arrays
		for (j = 0; j < ekused_tot; j++)
		{
			
			ek_kvalues[j] = ek_kvalues_ori[j];
			
		}
		for (j = 0; j < size; j++)
			ek_coeffs[j] = ek_coeffs_ori[j];
		for (i = 0; i < nek; i++)
		{
			ek_kfinder[i] = ek_kvalues + ekused_cum[i];
			ek_kfinder_ori[i] = ek_kvalues + ekused_cum[i];//this will point to the original k values for each data set, here there is only originals
			ek_efinder[i] = ek_evalues + ekused_cum[i];
			ek_cfinder[i] = ek_coeffs + (ek_cfinder_ori[i] - ek_cfinder_ori[0]);
		}

	}
	if (!(RunParams::continuation  && is_E0shift))//only not set, if there is E0 shift and continuation, as it was already read
	{
		for (i = 0; i < nek; i++)
			ek_gridind[i] = ek_gridstart[i]-ek_gridfrom[i];//first the coeff matrix will be used, which is indicated by ek_gridstart, if not given
							//the normal coeff matrix will be used, even if there will be E0 shift later 
	}

	if (!keep_k_ori)
		delete[] ek_kvalues_ori;
	delete[] ek_evalues_ori;
	delete[] ek_coeffs_ori;
	delete[] ek_efinder_ori;
	delete[] ek_cfinder_ori;

	
};

//calculate the actual I(Q)_corr in the loop, it should be called only in the steps, where new mu_act was generated
void ExptsData::CalcIQmucorr()
{
	int i, j;
	double *piq, *piqc, *pIB;

	
	for (i = 0; i < nfq; i++)
	{
		if (fqIQbackgcorr[i])
		{
			piq = fq_ffinder_ori[i];
			piqc = fq_ffinder[i];
			pIB = fq_IBfinder[i];

			for (j = 0; j < fqused[i]; j++)
				*piqc++ = *piq++ - fqmuact[i] * *pIB++;
		}
	}
};

//calculate the X-ray coefficients for the given data set
//mode=0: calculating the X-ray scattering factors and the coefficients (F(Q) and I(Q) fit)
//mode=1: calculating the coefficients based on the scattering factors read from file (I(Q) fit only)
//for I(Q) fit calculating A(Q) regardless the mode
void ExptsData::CalcXrayCoeffs(int i, int mode)
{
	int j,k,itype,jtype,ipartial;
	int npartials = ntypes * (ntypes + 1) / 2;
	double f_av=0,f_av2=0,*q,s2,*pA=nullptr,f_i;
	double fp_min;
	
	//it was already checked, that the symbol is in the table if CHEMICAL-SYMBOLS were read
	//check, whether the key word was given at all
	if (RunParams::chem_symbols == nullptr && mode==0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << CheckKey("CHEMICAL-SYMBOLS")<<" key word is needed in the [ " << CheckTag("GENERAL") << " ] section, if Xray atomic scattering " << endl;
		cout << "factors should be calculated in order to calculate the coefficients!" << endl;
		cout << "It is missing, so add it to " << datfilename << " data file and try again!" << endl;
		CleanExit();
	}
	
	q = fq_qfinder[i];
	if (fqfitIQ[i])
		pA = fq_A + fqused_cum[i];
		
	for (j = 0; j < fqused[i]; j++)
	{
		f_av = 0.0;
		if (mode == 0)//calculate the f(Q)
		{
			s2 = *q * *q / 16.0 / PI / PI;
			
			for (itype = 0; itype < ntypes; itype++)
			{
				*(fq_sffinder[i * ntypes + itype] + j) = 0.0;

				for (k = 0; k < 5; k++)
					*(fq_sffinder[i * ntypes + itype] + j) += Xray_coeffs_W.at(RunParams::chem_symbols[itype]).a[k] *
					exp(-Xray_coeffs_W.at(RunParams::chem_symbols[itype]).b[k] * s2);
				*(fq_sffinder[i * ntypes + itype] + j) += Xray_coeffs_W.at(RunParams::chem_symbols[itype]).c;
				fp_min = (fqfprime[i * ntypes + itype] * (1 + fqfprimefactor[i]) < fqfprime[i * ntypes + itype] * (1 - fqfprimefactor[i]) ? 
					fqfprime[i * ntypes + itype] * (1 + fqfprimefactor[i]) : fqfprime[i * ntypes + itype] * (1 - fqfprimefactor[i]));
				if (*(fq_sffinder[i * ntypes + itype] + j) + fp_min < 0)
				{
					cout << "\n*****ERROR*****" << endl;
					{
						cout << "The lowest limit for f(Q)+f' for the " << i + 1 << ". X-ray data set's " << itype+1<<". atom type's "<<j + 1 << ". data point is below zero!" << endl;
						cout << "Q: " << *(fq_qfinder[i] + j) << endl;
						cout << "f(Q) " << *(fq_sffinder[i * ntypes + itype] + j) << endl;
						cout << "f': " << fqfprime[i * ntypes + itype] <<" fraction: "<<fqfprimefactor[i]<< endl;
						cout << "Lowest limit for f(Q)+f': " << *(fq_sffinder[i * ntypes + itype] + j) + fp_min << endl;
						cout << "Give a different f' value or fraction and try again! Exiting..." << endl;
						CleanExit();

					}
				}
				f_av += SimpleCfg::fractions[itype] *( *(fq_sffinder[i * ntypes + itype] + j)+ fqfprimeact[i * ntypes + itype]);//<f+f'>
				
			}//end of itype cycle
			f_av2= f_av * f_av;//it is needed for F(Q) fit coeffcient normalization, in case of I(Q) fit it will calculated at saving only for saving the normalized F(Q)
		}//end of mode==0
				
		ipartial = 0;
		if (fqfitIQ[i])
			*pA = 0.0;//this is <f2>
	
		for (itype = 0; itype < ntypes; itype++)
		{
			f_i = *(fq_sffinder [i * ntypes + itype] + j)+fqfprimeact[i*ntypes+itype];//f' is added regardless if it used or not, it is 0 if not used
			if (fqfitIQ[i])
				*pA += SimpleCfg::fractions[itype] * f_i * f_i;
			
			for (jtype = itype; jtype < ntypes; jtype++)
			{
				*(fq_cfinder[i * npartials + ipartial] + j) = f_i * (*(fq_sffinder[i * ntypes + jtype] + j)+ fqfprimeact[i * ntypes + jtype]) * SimpleCfg::fractions[itype] * SimpleCfg::fractions[jtype];
				if (fqfitIQ[i] == false)//for F(Q) fitting the normalized coefficients are used (divided by <f>**2 to make their sum=1
					*(fq_cfinder[i * npartials + ipartial] + j) /= f_av2;
				if (itype != jtype)
					*(fq_cfinder[i * npartials + ipartial] + j) *= 2;
				ipartial++;
			}//end of jtype cycle
		}//end of itype cycle
		q++;
		
		if (fqfitIQ[i])
			pA++;
	}//end of data point cycle j
};

//calculating the change in the coeffs  due to the change of f' in case of AXS

void ExptsData::CalcXrayCoeffsChange()
{
	int i,j,itype,jtype,ipartial;
	int npartials = ntypes * (ntypes + 1) / 2;
	double f_i,*pA,*q;
		
	for (i = 0; i < nfq; i++)
	{
		if (fqAXS[i] == 0)
			continue;//no change
		q = fq_qfinder[i];
		pA = fq_A + fqused_cum[i];
		itype = fqAXS[i] - 1;//this is the type of the AXS particle
		

		for (j = 0; j < fqused[i]; j++)
		{
			f_i = *(fq_sffinder[i * ntypes + itype] + j) + fqfprimeact[i * ntypes + itype];
			*pA = SimpleCfg::fractions[itype] * f_i * f_i;
						
			for (jtype = 0; jtype < ntypes; jtype++)
			{
				double f_j = *(fq_sffinder[i * ntypes + jtype] + j) + fqfprimeact[i * ntypes + jtype];
				if (itype != jtype)
					*pA += SimpleCfg::fractions[jtype] * f_j * f_j;//this is <f2>
								
				ipartial = (jtype <= itype ? (jtype * ntypes - (jtype * (jtype + 1) / 2) + itype) : (itype * ntypes - (itype * (itype + 1) / 2) + jtype));
				//f' is added regardless if it used or not, it is 0 if not used
				*(fq_cfinder[i * npartials + ipartial] + j) = f_i * f_j	* SimpleCfg::fractions[itype] * SimpleCfg::fractions[jtype];

				if (itype != jtype)
					*(fq_cfinder[i * npartials + ipartial] + j) *= 2;
			}//end of jtype cycle
			q++;
			pA++;
		}//end of data point cycle j
	}//end of expt cycle i
};

//calculate the Compton contribution for the given data se
void ExptsData::CalcCompton(int i)
{
	int j, k, itype;
	double  *pB, *q, s2;
	string *mysymbols;
	pB = fq_B + fqused_cum[i];
	if (!fqusecompton[i])//whether calculate the B(Q)
	{
		for (j = 0; j < fqused[i]; j++)
			*pB++ = 0.0;
		return;
	}
	//it was already checked, that the symbol is in the table if CHEMICAL-SYMBOLS were read
	//check, whether the key word was given at all
	if (RunParams::chem_symbols == nullptr)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << CheckKey("CHEMICAL-SYMBOLS")<<" key word is needed in the [ " << CheckTag("GENERAL") << " ] section, if Xray Compton scattering" << endl;
		cout << "should be calculated! It is missing, so add it to "<<datfilename<<" data file and try again!" << endl;
		CleanExit();
	}
	SetArraysize(&mysymbols, ntypes, "mysymbols", "ExptsData::CalcCompton");//needed as ions should use the atom's values
	for (itype = 0; itype < ntypes; itype++)
	{
		if (Compton.count(RunParams::chem_symbols[itype]) > 0)
			mysymbols[itype] = RunParams::chem_symbols[itype];//atom, copy it
		else
		{
			mysymbols[itype] = RunParams::chem_symbols_standard[itype];
			cout << "\nNOTE(" << ++note << "): The Compton term is calculated for " << RunParams::chem_symbols[itype] << " using the parameters for " << mysymbols[itype] << endl;
		}
	}
	q = fq_qfinder[i];
	
	for (j = 0; j < fqused[i]; j++)
	{
		
		s2 = *q * *q / 16.0 / PI / PI;
		*pB = 0.0;
		for (itype=0;itype<ntypes;itype++)
		{
			double temp = Compton.at(mysymbols[itype]).atomic_number;
			for (k = 0; k < 5; k++)
				temp -= Compton.at(mysymbols[itype]).a[k] * exp(-Compton.at(mysymbols[itype]).b[k] * s2);
			*pB += SimpleCfg::fractions[itype] * temp;
		}//end of itype cycle
		pB++;
		q++;
	}//end of data point cycle j
};

//Determining the neutron coefficients, only call for sets, were it is needed
void ExptsData::CalcNDCoeffs(int i)
{
	//it was already checked, that the symbol is in the table if CHEMICAL-SYMBOLS were read
	//check, whether the key word was given at all
	if (RunParams::chem_symbols == nullptr)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << CheckKey("CHEMICAL-SYMBOLS") << " key word is needed in the [ " << CheckTag("GENERAL") << " ] section, if neutron scattering " << endl;
		cout << " coefficients should be calculated!" << endl;
		cout << "It is missing, so add it to " << datfilename << " data file and try again!" << endl;
		CleanExit();
	}
	int itype,jtype,n;
	int npartials = ntypes * (ntypes + 1) / 2;
	double sum=0;
	complex<double> *bcoh;
	SetArraysize(&bcoh, ntypes, "bcoh", "ExptsData::CalcNDCoeffs");
	//first determining, calculating the bound coherent scattering length for each componenet
	for (itype = 0; itype < ntypes; itype++)
	{
		bcoh[itype] = 0;
		n = isotope_count[i * ntypes + itype].count;
		if ( n> 1)
		{
			for (int k = 0; k < n; k++)
				bcoh[itype] += isotope_count[i * ntypes + itype].isotopes[k].ratio*N_scat_length.at(isotope_count[i * ntypes + itype].isotopes[k].symbol).b_coh;
		}
		else
		{
			//either natural abundance or one isotope, no calculation necessary
			bcoh[itype] = (n == 0 ? N_scat_length.at(RunParams::chem_symbols_standard[itype]).b_coh :
				N_scat_length.at(isotope_count[i * ntypes + itype].isotopes[0].symbol).b_coh);
		}
		//cout << "bcoh " << itype << " " << bcoh[itype] << endl;
	}
	int ipartial = 0;
	for (itype = 0; itype < ntypes; itype++)
	{
		sq_coeffs[npartials * i + ipartial] = pow(SimpleCfg::fractions[itype] * bcoh[itype].real(), 2);
		//cout << itype << " " << itype << " ip " << ipartial << endl;
		sum += sq_coeffs[npartials * i + ipartial];
		ipartial++;
		for (jtype = itype+1; jtype < ntypes; jtype++)
		{
			
			sq_coeffs[npartials * i + ipartial] = 2* SimpleCfg::fractions[itype]* SimpleCfg::fractions[jtype]* bcoh[itype].real() * bcoh[jtype].real();//it was decided to use only the real parts...
			//cout << itype << " " << jtype << " ip " << ipartial << endl;
			sum+=sq_coeffs[npartials * i + ipartial];
			ipartial++;
		}
		
	}
	for (ipartial = 0; ipartial < npartials; ipartial++)
	{
		sq_coeffs[npartials * i + ipartial] /= sum;
		//cout << ipartial << " coeff " << sq_coeffs[npartials * i + ipartial] << endl;
	}
	delete[] bcoh;

}

