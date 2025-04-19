//source CosDistrConst.cpp
//Last changed 24.01.2023

//The number of bins is regulated through the cos(theta) spacing, which is a parameter in the *.dat file.
//A constraint can be "positive", which means, that bond angles are needed in the given region, and "negative" meaning that angles are not wanted in the given region.
//Positive constraint can be set up by specifying either 0 for the calculation method indicating a step function, or 1 meaning Gaussian for the shape of the constraint.
//Step function means having uniform distribution with integral 1 between angle - wcontrol ->angle + wcontrol, (both given in degrees), 0 otherwise. 
//Gaussian indicates a normal distribution with integral 1 having the maximum at angle (given in degree), which is converted to radians and the cosine of it is calculated.
//This cosine(theta) value will give the peak position, and the width is controlled by sigma=wcontrol, where wcontrol is used as it is given in the *.dat file,  so it is 
//a dcos(theta) value. 
//In these cases the theoretical distribution is calculated for all the bins spanning the range -1=<(cosine(theta)=<1, so only one positive constraint (desired angle) can be meaningfully 
//given for a neighbour1-central-neighbour2 triplet, as more than 1 constraint would ruin each others effect.
//Negative constraint, that no angles desired in a given region can be set up using method=2. In this case the constraint is for the given region and not for all the bins, 
//as it could interfere otherwise with a "positive" constraint for the same particle types. A normal distribution curve similarly to method=1 is calculated from that part of 
//-CONF_INT*wcontrol -> + CONF_INT*wcontrol, which is in the -1 -> +1 range, centred around the not desired angle. The calculated curve is stored in theoretical distribution
//but during the chi2 calculation it is used to provide a cos(theta)-dependent weight by multiplying the calculated distribution and calculating the product's squared difference from 0. 
//Make sure, that the interval does not overlap with a positive constraint for the same particle triplet!
//If method=3 is given, then instead of angle and wcontrol a file name is read, and an experinetal distribution is read from this file. In this case the
//number of data points has to be given in the first line of the file. The distribution will be normalized.

#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Move.h"//classes1.h is included through this

int		 CosDistrConst::nconstraints;//the number of constraints
int		 CosDistrConst::ncos_bin;//number of bins 
int		 CosDistrConst::ntot_bins;//total number of data points
int		*CosDistrConst::npoints;//number of data points for each contsraints;
int		*CosDistrConst::cumul_bins;//cumulative number of bins for each constraint
int		*CosDistrConst::max_neigh;//pointer to the maximum number of neighbours for the neighbourlist for an atom
int		*CosDistrConst::method;//array of switches for the calculation method: 0: step function, 1: normal distribution, 2: the specified angle is not desired
int		*CosDistrConst::central;//array of types of central particles
int		*CosDistrConst::neighbour1;//array of types of neighbour1 particles
int		*CosDistrConst::neighbour2;//array of types of neighbour2 particles
longint *CosDistrConst::consttype;//pointer to NeighbourList::consttype, each bit represents a type, containing a 1 for each type that is involved in a constraint
longint CosDistrConst::cosconsttype=0;//each bit represents a type, containing a 1 for each type that is involved in a COS constraint
double	*CosDistrConst::dcos_theta;//width of the cosine theta bin (full range goes -1 -> +1) for each constraint (can differ from the others for method=3 constraints)
double	*CosDistrConst::inv_binwidth;// 1./dcos_theta for each constraint (can differ from the others for method=3 constraints)
double	*CosDistrConst::dmin1;//array of minimal distances for the neighbour1
double	*CosDistrConst::dmin2;//array of minimal distances for the neighbour2
double	*CosDistrConst::dmax1;//array of maximum distances for the neighbour1
double	*CosDistrConst::dmax2;//array of maximum distances for the neighbour2
double	*CosDistrConst::angle;//array of desired bond angle
double	*CosDistrConst::wcontrol;//array of parameters controlling the width of the distribution
double	*CosDistrConst::weights;//array of weights for each CC 
	
int		*CosDistrConst::mod_const;//whteher the constraint was effected by the move
int		*CosDistrConst::neigh_partial;//index of the partial calculated from the two neighbour types
int		*CosDistrConst::first_bin;//index of the first bin, needed for the negative constraints
int		*CosDistrConst::last_bin;//index of the last bin, needed for the negative constraints	
double   CosDistrConst::maxdistsq;//the maximum distance among the max1sq and max2sq arrays, this will be used to determine the limit of the neighbourlist calculation
double	*CosDistrConst::min1sq;//reduced min. dist. squared
double	*CosDistrConst::min2sq;//reduced min. dist. squared
double	*CosDistrConst::max1sq;//reduced max. dist. squared
double	*CosDistrConst::max2sq;//reduced max. dist. squared
double  *CosDistrConst::theordistr;//array of theoretical or experimental distribiutions
char   (*CosDistrConst::datafilename)[FILE_NAME_SIZE];//names of the experimental data files

//----default constructor----------------
CosDistrConst::CosDistrConst(NeighbourList &neigh_list, SimpleCfg &conf)
				:neighlist(neigh_list), config(conf)
{
	if (nconstraints < 1)
	{
		if (::debug)
		{
			cout << "\nWARNING(" << ++warn << "): CosDistrConst constructor" << endl;
			cout << "\tThere are no cosine distribution of the bond angles constraints, or nconstraints was not initialised!" << endl;
		}
	}
	else
	{
		if (ncos_bin<1)
		{
			cout<<"\n****ERROR****"<<endl;
			cout<<"CosDistrConst constructor"<<endl;
			cout<<"Number of bins is "<<ncos_bin<<", cannot run this way, exiting..."<<endl;
			CleanExit();
		}
	}
	//The cos_hist and cosinedistr arrays will be set later in CalcTheoretical

};

//----copy constructor----------------
CosDistrConst::CosDistrConst(CosDistrConst &source, NeighbourList &neigh_list, SimpleCfg &conf)
				:neighlist(neigh_list), config(conf)
{
	int i;
		
	SetArraysize(&cos_hist,ntot_bins,"cos_hist","CosDistrConst::CosDistrConst");//histogram for the calculated distribution
	SetArraysize(&cosinedistr,ntot_bins,"cosinedistr","CosDistrConst::CosDistrConst");//calculated distribution
		
	//copying the arrays
	for (i=0;i<ntot_bins;i++)
	{
		cos_hist[i]=source.cos_hist[i];
		cosinedistr[i]=source.cosinedistr[i];
	}

};


//--------------gets the static parameters------------------------------------
void CosDistrConst::GetCosDistrConst(ifstream &file)
{
	int i;
	double dcos_theta_in,temp;
	char *name,*conv_numb = NULL;
	if (CheckFileState(file,"CosDistrConst::GetCosDistrConst",datfilename)==0)
	{
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	
	nconstraints=RunParams::ncosdistr;
	consttype = &NeighbourList::consttype;
	
	SetArraysize(&name,NAME_SIZE,"name","CosDistrConst::GetCosDistrConst");
	//creating the static arrays
	SetArraysize(&neigh_partial,nconstraints,"neigh_partial","CosDistrConst::GetCosDistrConst");
	SetArraysize(&mod_const,nconstraints,"mod_const","CosDistrConst::GetCosDistrConst");
	SetArraysize(&first_bin,nconstraints,"first_bin","CosDistrConst::GetCosDistrConst");
	SetArraysize(&last_bin,nconstraints,"last_bin","CosDistrConst::GetCosDistrConst");
	SetArraysize(&central,nconstraints,"central","CosDistrConst::GetCosDistrConst");
	SetArraysize(&method,nconstraints,"method","CosDistrConst::GetCosDistrConst");
	SetArraysize(&neighbour1,nconstraints,"neighbour1","CosDistrConst::GetCosDistrConst");
	SetArraysize(&neighbour2,nconstraints,"neighbour2","CosDistrConst::GetCosDistrConst");
	SetArraysize(&dmin1,nconstraints,"dmin1","CosDistrConst::GetCosDistrConst");
	SetArraysize(&dmin2,nconstraints,"dmin2","CosDistrConst::GetCosDistrConst");
	SetArraysize(&dmax1,nconstraints,"dmax1","CosDistrConst::GetCosDistrConst");
	SetArraysize(&dmax2,nconstraints,"dmax2","CosDistrConst::GetCosDistrConst");
	SetArraysize(&angle,nconstraints,"angle","CosDistrConst::GetCosDistrConst");
	SetArraysize(&wcontrol,nconstraints,"wcontrol","CosDistrConst::GetCosDistrConst");
	SetArraysize(&weights,nconstraints,"weights","CosDistrConst::GetCosDistrConst");
	SetArraysize(&min1sq,nconstraints,"min1sq","CosDistrConst::GetCosDistrConst");
	SetArraysize(&min2sq,nconstraints,"min2sq","CosDistrConst::GetCosDistrConst");
	SetArraysize(&max1sq,nconstraints,"max1sq","CosDistrConst::GetCosDistrConst");
	SetArraysize(&max2sq,nconstraints,"max2sq","CosDistrConst::GetCosDistrConst");
	SetArraysize(&npoints,nconstraints,"npoints","CosDistrConst::GetCosDistrConst");
	SetArraysize(&dcos_theta,nconstraints,"dcos_theta","CosDistrConst::GetCosDistrConst");
	SetArraysize(&inv_binwidth,nconstraints,"inv_binwidth","CosDistrConst::GetCosDistrConst");
	SetArraysize(&cumul_bins,nconstraints+1,"cumul_bins","CosDistrConst::GetCosDistrConst");

	datafilename=new char[nconstraints][FILE_NAME_SIZE];
	if (datafilename==NULL)
		NoArray("CosDistrConst::GetCosDistrConst","datafilename");

	maxdistsq=0.0;

	dcos_theta_in=ReadThisLine(file,6,DCOSTHETA_DEF,"cos(theta) spacing", "CosDistrConst::GetCosDistrConst",cosfilename);//spacing in cosine theta space (width of the bin)
	//calculating the number of bins
	temp=2./dcos_theta_in;
	ncos_bin=(int)temp;
	if (temp-ncos_bin>0.5)
		ncos_bin++;//to find the nearest integer
	
	for (i=0;i<nconstraints;i++)
	{
		IntToStr(&conv_numb, i + 1);
		//reset the bin width by default, it will be recalculated for method=3 later
		dcos_theta[i]=2./ncos_bin;
		inv_binwidth[i]=ncos_bin/2.0;//inverse of the bin width

		method[i]=ReadThisLine(file,3,1,"method", "CosDistrConst::GetCosDistrConst");
		if (method[i]<0 || method[i]>3)
		{
			cout<<"CosDistrConst::GetCosDistrConst"<<endl;
			cout<<"The calculation method indicator for the "<<i+1<<". cosine distribution of bond angles constraint is "<<method[i]<<"."<<endl;
			cout<<"It can only be 0, 1, 2 or 3! 1 (positive constraint with normal distribution) will be assumed!"<<endl;
			method[i]=1;	
		}
		mystrcpy(name, NAME_SIZE, "central[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		central[i]=ReadThisLine(file,3,1,name, "CosDistrConst::GetCosDistrConst");
		RunParams::CheckType(central[i],name, "CosDistrConst::GetCosDistrConst");
		central[i]--;
		mystrcpy(name, NAME_SIZE, "neighbour1[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		neighbour1[i]=ReadThisLine(file,3,1,name, "CosDistrConst::GetCosDistrConst");
		RunParams::CheckType(neighbour1[i], name, "CosDistrConst::GetCosDistrConst");
		neighbour1[i]--;
		mystrcpy(name, NAME_SIZE, "neighbour2[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		neighbour2[i]=ReadThisLine(file,3,1,name, "CosDistrConst::GetCosDistrConst");
		RunParams::CheckType(neighbour2[i], name, "CosDistrConst::GetCosDistrConst");
		neighbour2[i]--;
		mystrcpy(name, NAME_SIZE, "dmin1[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmin1[i]=ReadThisLine(file,3,1.0,name, "CosDistrConst::GetCosDistrConst");
		mystrcpy(name, NAME_SIZE, "dmin2[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmin2[i]=ReadThisLine(file,3,1.0,name, "CosDistrConst::GetCosDistrConst");
		if (neighbour1[i] == neighbour2[i]  && fabs(dmin2[i] - dmin1[i]) > LOAD_TOL)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The neighbour types for the " << i + 1 << " cosine distarubution of angles constraint are the same, but the minimum distances are different!" << endl;
			cout << "dmin1 is " << dmin1[i] << " A and dmin2 is " << dmin2[i] << " A." << endl;
			cout << "This does not make sense, and would lead to inconsistency!" << endl;
			cout<<"Check the " << datfilename << " file, and try again! Exiting..." << endl;
			CleanExit();
		}
		mystrcpy(name, NAME_SIZE, "dmax1[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmax1[i]=ReadThisLine(file,3,1.0,name, "CosDistrConst::GetCosDistrConst");
		mystrcpy(name, NAME_SIZE, "dmax2[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmax2[i]=ReadThisLine(file,3,1.0,name, "CosDistrConst::GetCosDistrConst");
		if (neighbour1[i] == neighbour2[i] && fabs(dmin2[i] - dmin1[i]) > LOAD_TOL)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The neighbour types for the " << i + 1 << " cosine distarubution of angles constraint are the same, but the maximum distances are different!" << endl;
			cout << "dmax1 is " << dmax1[i] << " A and dmax2 is " << dmax2[i] << " A." << endl;
			cout << "This does not make sense, and would lead to inconsistency!" << endl;
			cout << "Check the " << datfilename << " file, and try again! Exiting..." << endl;
			CleanExit();
		}
		if (method[i]<3)
		{
			mystrcpy(name, NAME_SIZE, "angle[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, "]");
			angle[i]=ReadThisLine(file,3,1.0,name, "CosDistrConst::GetCosDistrConst");
			mystrcpy(name, NAME_SIZE, "wcontrol[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, "]");
			wcontrol[i]=ReadThisLine(file,3,1.0,name, "CosDistrConst::GetCosDistrConst");
			if (wcontrol[i]==0)
				method[i]=0;//the 0 width case will be dealt with when the method swith is 0! 
		}
		else
		{
			//The data will be read from a file, read the file name
			file>>datafilename[i];//file name
		}
		mystrcpy(name, NAME_SIZE, "sigma[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		weights[i]=ReadThisLine(file,1,1.0,name, "CosDistrConst::GetCosDistrConst",cosfilename);
		if (weights[i]<0)
			ChiSquared::calc_sigma=1;

		//calculating the reduced squared distances
		min1sq[i]=pow(dmin1[i]/RunParams::boxedge,2.);
		min2sq[i]=pow(dmin2[i]/RunParams::boxedge,2.);
		max1sq[i]=pow(dmax1[i]/RunParams::boxedge,2.);
		max2sq[i]=pow(dmax2[i]/RunParams::boxedge,2.);
		
		//finding the largest reduced maximum distance square
		if (max1sq[i]>maxdistsq)
			maxdistsq=max1sq[i];
		if (max2sq[i]>maxdistsq)
			maxdistsq=max2sq[i];

		//The concept is that the neighbourlist is calculated for atoms of each type represented in a cosine distribution
		//constraint or in common or second neighbour constraint. It is required for the neighbours as well, as the list is updated, and not completly recalculated in each step.
		//As each type can be involved in more, than one constraint, all the neighbour lists will be calculated with the
		//same maximum distance, maxdistsq. 
		//Be aware that the central and the neighbour type in cosconsttype types represented by right to left!
		//type 0 is represented by binary 1, type 1 binary 10 and so on  
		//consttype i a pointer to NeighbourList::consttype for all the constarins using the neighbour list
		//cosconsttype only for cos constraints
		if (!(*consttype>>central[i] & 1))//this type is not yet represented
		{
			*consttype+=(longint)pow(2.0,(double)central[i]);
			NeighbourList::nconsttype++;//increase the number of constrained types
		}
		if (!(*consttype>>neighbour1[i] & 1))//this type is not yet represented
		{
			*consttype+=(longint)pow(2.0,(double)neighbour1[i]);
			NeighbourList::nconsttype++;//increase the number of constrained types
		}
		if (!(*consttype>>neighbour2[i] & 1))//this type is not yet represented
		{
			*consttype+=(longint)pow(2.0,(double)neighbour2[i]);
			NeighbourList::nconsttype++;//increase the number of constrained types
		}
		if (!(cosconsttype >> central[i] & 1))//this type is not yet represented
			cosconsttype += (longint)pow(2.0, (double)central[i]);
			
		if (!(cosconsttype >> neighbour1[i] & 1))//this type is not yet represented
			cosconsttype += (longint)pow(2.0, (double)neighbour1[i]);
		
		if (!(cosconsttype >> neighbour2[i] & 1))//this type is not yet represented
			cosconsttype += (longint)pow(2.0, (double)neighbour2[i]);
		
		//Determining the index of the partial calculated from the neighbour types 
		neigh_partial[i]=(neighbour1[i]<=neighbour2[i] ? (neighbour1[i]*RunParams::ntypes-(neighbour1[i]*(neighbour1[i]+1)/2)+neighbour2[i]) : \
			(neighbour2[i]*RunParams::ntypes-(neighbour2[i]*(neighbour2[i]+1)/2)+neighbour1[i]));

		mod_const[i]=1;//to calculate the initial distribution
	}	

	
	ntot_bins=nconstraints*ncos_bin;//by default, this can be increased later, if necessary due to 3-type constraint
	//creating the array for the theoretical distribution 
	SetArraysize(&theordistr,ntot_bins,"theordistr","CosDistrConst::GetCosDistrConst");//theoretical distribution
	if (conv_numb != NULL)
		delete [] conv_numb;
	delete [] name;

	if (!CheckReadFileState(file,"CosDistrConst::GetCosConst",datfilename))
		CleanExit();//loading failed

	
};
//-----------initializing some of the static members-----------------------
void CosDistrConst::SetParams()
{

	max_neigh=&NeighbourList::max_neigh;
	

};
//---------------sets the elements of array cos_hist to 0-----------------
void CosDistrConst::InitBincount()
{
	int i;
	int *p;
	p=cos_hist;
	for (i=0;i<ntot_bins;i++)
		*p++=0;

};

//---------------calculate the desired theoretical distribution-------------------
void CosDistrConst::CalcTheoretical()
{
	int i,ibin;
	int count,first_non0_bin;
	double *p;
	double temp,sigma;
	double step_min,step_max,sum,cos_angle;

	cumul_bins[0]=0;
	for (i=0;i<nconstraints;i++)
	{
		p=theordistr+cumul_bins[i];
					
		switch (method[i])
		{
			case 0: //to have uniform distribution between angle-wcontrol ->angle+wcontrol, 0 otherwise 
			{
				npoints[i]=ncos_bin;
				cumul_bins[i+1]=cumul_bins[i]+ncos_bin;
				first_bin[i]=0;
				last_bin[i]=ncos_bin;
				//angle and wcontrol are given in degree
				step_min=angle[i]-wcontrol[i];//lower limit in degree
				step_min=cos(step_min*PI/180);//lower limit in cosine
				step_max=angle[i]+wcontrol[i];//upper limit in degree
				step_max=cos(step_max*PI/180);//upper limit in cosine
				//swap the boundaries, if necessary 
				if (step_min>step_max)
				{
					temp=step_min;
					step_min=step_max;
					step_max=temp;
				}

				count=0;
				first_non0_bin=-1;

				if (wcontrol[i]==0 || step_max-step_min<dcos_theta[i])//To handle the zero-width and the smaller than binwidth case
				{
					cos_angle=cos(angle[i]*PI/180);//position of the maximum
					//there will be only one non-zero bin
					for (ibin=0;ibin<ncos_bin;ibin++)
					{
						if ((ibin*dcos_theta[i]-1)<=cos_angle && ((ibin+1)*dcos_theta[i]-1)>cos_angle)//cos_angle is in the ibin-th bin
						{
							*p++=1;
							count++;
							if (first_non0_bin==-1)
								first_non0_bin=ibin;//preserve the index of the first non-zero bin
						}
						else
							*p++=0;
					}
				}
				else //the step function is at least one bin wide
				{
					for (ibin=0;ibin<ncos_bin;ibin++)
					{
						if ((ibin+0.5)*dcos_theta[i]-1 >=step_min &&  (ibin+0.5)*dcos_theta[i]-1<=step_max)
						{
							*p++=1;
							count++;
							if (first_non0_bin==-1)
								first_non0_bin=ibin;//preserve the index of the first non-zero bin
						}
						else
							*p++=0;
					}
				}
				p=theordistr+cumul_bins[i]+first_non0_bin;
				//normalize the distribution to have 1 for the integral
				for (ibin=first_non0_bin;ibin<first_non0_bin+count;ibin++)
				{
					*p/=count*dcos_theta[i];
					p++;
				}
				break;
			}
			case 3: //load it from a file
			{
	npoints[i]=ncos_bin;
	cumul_bins[i+1]=cumul_bins[i]+ncos_bin;
	first_bin[i]=0;
	last_bin[i]=ncos_bin;
				LoadCosDistr(i);
	
				break;
			}
			default: //(case 1 and 2: normal distribution, with maximum at angle (given in degree), and sigma=wcontrol (given in cos(theta) unit
			{
				npoints[i]=ncos_bin;
				cumul_bins[i+1]=cumul_bins[i]+ncos_bin;
				//calculte the cosine of angle and sigma 
				cos_angle=cos(angle[i]*PI/180);//position of the maximum
				sigma=wcontrol[i];//sigma
				temp=1/sqrt(2*PI)/sigma;
				if (method[i]==1)
				{
					first_bin[i]=0;
					last_bin[i]=ncos_bin;
				}
				else
				{
					//setting the boundaries of the bin interval according to confidence interval to the CONF_INT level given in iunits.h
					first_bin[i]=int((cos_angle-CONF_INT*sigma+1)/dcos_theta[i]);
					first_bin[i]=(first_bin[i]<0 ? 0 : first_bin[i]);//reset it to zero, if the first would be outside the range
					last_bin[i]=int((cos_angle+CONF_INT*sigma+1)/dcos_theta[i]);
					last_bin[i]=(last_bin[i]>ncos_bin ? ncos_bin : last_bin[i]);//reset it to ncos_bin, if the last would be outside the range
				}
				sum=0;

				for (ibin=0;ibin<first_bin[i];ibin++)
					*p++=0;//just to set it to something, never used
				for (ibin=first_bin[i];ibin<last_bin[i];ibin++)
				{
					*p=temp*exp(-pow(((ibin+0.5)*dcos_theta[i]-1-cos_angle)/sigma,2.)/2);
					sum+=*p;
					p++;
				}
				for (ibin=last_bin[i];ibin<ncos_bin;ibin++)
					*p++=0;//just to set it to something, never used
				//normalize the distribution to have 1 for the integral
				p=theordistr+cumul_bins[i]+first_bin[i];
				for (ibin=first_bin[i];ibin<last_bin[i];ibin++)
				{
					*p/=sum*dcos_theta[i];
					p++;
				}
				p+=ncos_bin-last_bin[i];//to set the pointer to the next data set
				break;
								
			}

		}//end switch
	}//end i cycle
	SetArraysize(&cos_hist,ntot_bins,"cos_hist","CosDistrConst::CosDistrConst");//histogram for the calculated distribution
	SetArraysize(&cosinedistr,ntot_bins,"cosinedistr","CosDistrConst::CosDistrConst");//calculated distribution
	

};

//loading the cosine distribution from a file, if option 3 is given
void CosDistrConst::LoadCosDistr(int i)
{
	int j,mode;
	double n,sum=0;
	double *pcos, cos1=0,cos2=0, dcos1=0,dcos2;
	std::streamoff current_pos;
	ifstream file;//used for the experimental data file
	char *base_name,*conv_numb=NULL;
	SetArraysize(&base_name,NAME_SIZE,"base_name","CosDistrConst::LoadCosDistr");

	cout<<"\nLoading data from cosine distribution of bond angle file: "<<datafilename[i]<<endl;
	SafeOpenTextFile(file,datafilename[i]);
	if (CheckFileState(file,"CosDistrConst::LoadCosDistr",datafilename[i])==0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"Cannot load data, exiting..."<<endl;
		CleanExit();
	};
				
	npoints[i]=ReadThisLine(file,1,1,"cos distr npoints", "CosDistrConst::LoadCosDistr",datafilename[i]);//total number of r points in the data file

	if (npoints[i]<0)// can happen, if the number of data points are not properly given in the file
	{
		cout<<"\n****ERROR****"<<endl;
		cout<<"The number of used data points is smaller, than zero, something must be wrong with the format of the file!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
				
	SkipLine(file,datafilename[i],1);//comment line
	if (npoints[i]>ncos_bin)//it is necessary to resize the array, as the number of data points for this constraint is greater, than ncos_bin
		ResizeArray(&ntot_bins,ntot_bins-ncos_bin+npoints[i],&theordistr,"theordistr","CosDistrConst::LoadCosDistr");
	else
		ntot_bins+=npoints[i]-ncos_bin;

	//first go through the dcos_theta data, to see, if it is equdistant
	current_pos=file.tellg();//save the current file position
	mode=1;//go to next line
	for (j=0;j<npoints[i];j++)
	{
		if (j==npoints[i]-1)
			mode=3;//not to go to next line
		if (j>0)
			cos1=cos2;//cos1 will be the previous data point
		cos2=ReadThisLine(file,mode,1.0,"cos(theta)", "CosDistrConst::LoadCosDistr",datafilename[i]);//loads the r value
		if (!CheckReadFileState(file,"CosDistrConst::LoadCosDistr",datafilename[i]))
			CleanExit();
		if (j==0)
			cos1=cos2;//the first difference will be preserved
		else
		{

			dcos2=cos2-cos1;//there is no previous dcos
			if (j==1)
				dcos1=dcos2;
			else
			{
				if (fabs(dcos1-dcos2)>LOAD_TOL)
				{
					cout<<"\n****ERROR****"<<endl;
					cout<<"The delta cos(theta) between the "<<j<<".  and "<<j+1<<". data points of the "<<i+1<<". cosine disribution"<<endl;;
					cout<<"of bond angles constraint is outside the "<<LOAD_TOL<<" load tolerance compared to the dcos theta value "<<dcos1<<endl;
					cout<<"calculated from the 1. and 2. data point! The cos(theta) points have to be equdistant, and have to denote the"<<endl;
					cout<<"middle value of the bins, the beginning of the first bin starting with -1, and the -1 -> +1 interval had to be "<<endl;
					cout<<"divisible with dcos(theta) without remainder! Correct the "<<datafilename[i]<<" file and try again."<<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

			}
		}
	}
	//now create the binning, the cos(theta) data points has to be
	n=2./dcos1;
	if ( fabs((double)(n-(int)n))>LOAD_TOL || (int)n!=npoints[i])
	{
		cout<<"\n*****ERROR*****"<<endl;
		cout<<"it seems, that the dcos(theta) intervals does not begin with -1 or ends with 1, or the number of data points given  "<<endl;
		cout<<"at the beginning of the "<<datafilename[i]<<" file does not correspond to the dcos(theta) interval "<<endl;
		cout<<"calculated from the cos(theta) values for the "<<i+1<<". cosine disribution"<<endl;;
		cout<<"of bond angles constraint The cos(theta) points have to be equdistant, and have to denote the"<<endl;
		cout<<"middle value of the bins, the beginning of the first bin starting with -1, and the -1 -> +1 interval had to be "<<endl;
		cout<<"divisible with dcos(theta) without remainder! Correct the "<<datafilename[i]<<" file and try again."<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}

	if (!CheckReadFileState(file,"CosDistrConst::LoadCosDistr",datafilename[i]))
		CleanExit();
	dcos_theta[i]=dcos1;	
	inv_binwidth[i]=1./dcos_theta[i];//inverse of the bin width
	file.seekg(current_pos,ios::beg);//resets the position in the file for the reading of the cos data
		
	pcos=theordistr+cumul_bins[i];//pointer for the distribution values

	mode=1;
	for (j=0;j<npoints[i];j++)
	{
		mystrcpy(base_name, NAME_SIZE, "for ");
		IntToStr(&conv_numb, i + 1);
		mystrcat(base_name, NAME_SIZE, conv_numb);
		mystrcat(base_name, NAME_SIZE, ". cosine distribution of bond angle constraint ");
		IntToStr(&conv_numb, j + 1);
		mystrcat(base_name, NAME_SIZE, conv_numb);
		mystrcat(base_name, NAME_SIZE, "distribution value");
	
		if (j==npoints[i]-1)
			mode=3;
		//do not skip line, as this can be the last line
		cos1=ReadThisLine(file,3,0.0,"dcos", "CosDistrConst::LoadCosDistr",datafilename[i]);//loads the cos(theta)
		*pcos=ReadThisLine(file,mode,0.0,base_name, "CosDistrConst::LoadCosDistr",datafilename[i]);//loads the distribution 
		sum+=*pcos;
		pcos++;
	}
	if (conv_numb != NULL)
		delete [] conv_numb;
	if (!CheckReadFileState(file,"CosDistrConst::LoadCosDistr",datafilename[i]))
		CleanExit();
	next_line_pos = -1;
	file.close();

	
	//setting the cumulative number of used g(r) points in the previous data sets for each set
	cumul_bins[i+1]=cumul_bins[i]+npoints[i];
	first_bin[i]=0;
	last_bin[i]=npoints[i];
	pcos=theordistr+cumul_bins[i];//pointer for the distribution values
	for (j=first_bin[i];j<last_bin[i];j++)
	{
		*pcos/=sum*dcos_theta[i];
		pcos++;
	}
	delete [] base_name;

};
//-----------calculating the cosine distribution of the bond angles---------------
void CosDistrConst::CalcDistribution()
{
	int iconst,ibin;
	int *pcos_hist;//pointer to the histogram of the calculated distribution
	double *pdistr;//pointer to the calculated distribution
	double norm;

	for (iconst=0;iconst<nconstraints;iconst++)
	{
		if (mod_const[iconst])
		{

			//the distribution has to be calculated for this constraint
			//calculating the distribution by normalizing the histogram 
			norm=0;
			pcos_hist=cos_hist+cumul_bins[iconst];
			for (ibin=0;ibin<npoints[iconst];ibin++)
				norm+=*pcos_hist++;
			pcos_hist=cos_hist+cumul_bins[iconst];
			pdistr=cosinedistr+cumul_bins[iconst];
			if (norm!=0)
				norm=inv_binwidth[iconst]/norm;//to create a multiplying factor
			for (ibin=0;ibin<npoints[iconst];ibin++)
				*pdistr++=*pcos_hist++ * norm;
		}
	};
};

//-------------------calculate the initial histogram--------------------
void CosDistrConst::CalcHist()
{
	//the pointers and auxiliary variables are used to extract array values only ones to make calculation hopefully quicker
	int iconst,i,inb1,inb2,ibin,icoord;
	int imin, imax;//cycle boundaries
	int partialind;
	int neigh_type1,neigh_type2;//desired type of the neighbours 
	int *pnneigh;//pointer to the number of neighbours 
	int *plist,*pind1,*pind2;//pointers to the indices of the neighbours
	int *ptypelist,*pneigh_type1,*pneigh_type2;//pointers to the type of the neighbours 
	int *pcos_hist;//pointer to the histogram of the calculated distribution

	double costheta;
	double dminsq1,dmaxsq1,dminsq2,dmaxsq2;//minimum and maximum squared distances for the two neighbours
	double *pdistsq,*pdsq1,*pdsq2;//pointers to the distances of the neighbours
	double *pvec_comp,*pvcomp1,*pvcomp2;//pointers to the central-neighbour vector components

	for (iconst=0;iconst<nconstraints;iconst++)
	{
		//setting auxiliary variables and pointers to increase the speed
		dminsq1=min1sq[iconst];
		dminsq2=min2sq[iconst];
		dmaxsq1=max1sq[iconst];
		dmaxsq2=max2sq[iconst];
		neigh_type1=neighbour1[iconst];
		neigh_type2=neighbour2[iconst];
		pcos_hist=cos_hist+cumul_bins[iconst];//beginning of the histogram for this contraint
		pnneigh=neighlist.nneigh_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in nneigh array
		plist=neighlist.neighlist_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in neigh_list array
		ptypelist=neighlist.neightype_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in neigh_type array 
		pdistsq=neighlist.dsq_finder[central[iconst]];//sets the pointer to the beginning ofof the values belonging to the central atom type in dsq array
		pvec_comp=neighlist.vector_finder[central[iconst]];//sets the pointer to the beginning of the vector components of the central atom type
		imin=SimpleCfg::cumul[central[iconst]];
		imax=SimpleCfg::cumul[central[iconst]+1];
		for (i=imin;i<imax;i++)//central atom cycle
		{
			
			//Calculating the histogram
			pneigh_type1=ptypelist;//type of the first neighbour
			pind1=plist;//index of the first neighbour
			pdsq1=pdistsq;//squared distance of the central atom and the first neighbour
			pvcomp1=pvec_comp;//central-neighbour1 vector components
			for (inb1=0;inb1<(*pnneigh-1);inb1++)
			{
				//determine, whether the first neighbour is from the right type, and is in the range
				if ((*pneigh_type1==neigh_type1 && *pdsq1>=dminsq1 && *pdsq1<=dmaxsq1) || \
					(*pneigh_type1==neigh_type2 && *pdsq1>=dminsq2 && *pdsq1<=dmaxsq2))
				{

					pneigh_type2=pneigh_type1+1;//type of the second neighbour
					pind2=pind1+1;//index of the second neighbour
					pdsq2=pdsq1+1;//squared distance of the central atom and the second neighbour 
					pvcomp2=pvcomp1+3;//vector components of the central-neighbour2 vector
					for (inb2=inb1+1;inb2<*pnneigh;inb2++)
					{
						//determining the index of the partial calculated from the neighbours					
						partialind=(*pneigh_type1<=*pneigh_type2 ? (*pneigh_type1*RunParams::ntypes-(*pneigh_type1*(*pneigh_type1+1)/2)+*pneigh_type2) : \
								(*pneigh_type2*RunParams::ntypes-(*pneigh_type2*(*pneigh_type2+1)/2)+*pneigh_type1));
						if (partialind==neigh_partial[iconst])//both neighbours are from the right types
						{
							//determine, whether the second neighbour is in the range
							if ((*pneigh_type2==neigh_type1 && *pdsq2>=dminsq1 && *pdsq2<=dmaxsq1) || \
								(*pneigh_type2==neigh_type2 && *pdsq2>=dminsq2 && *pdsq2<=dmaxsq2))
							{
								//calculating the cosine of the angle of neighbour1-central-neighbour2
								costheta=0;
								for (icoord=0;icoord<3;icoord++)
									costheta+=*(pvcomp1+icoord) * *(pvcomp2+icoord);
								costheta/=sqrt(*pdsq1 * *pdsq2);//dividing by the length of the vectors
								//to avoid rounding errors, put back into the -1 - +1 range
								costheta=(costheta<-1 ? -1 : (costheta>1 ? 1 :costheta));
								ibin=int((costheta+1.)*inv_binwidth[iconst]);
								(*(pcos_hist+ibin))++;//increase the count of his bin
							}//second neighbour is in the range

						}//end if the neighbours are from the right types
						
						pneigh_type2++;
						pind2++;
						pdsq2++;
						pvcomp2+=3;
					}//end of first neighbour cycle inb2
				}//first neighbour is from the right type, and is in the range

				pneigh_type1++;
				pind1++;
				pdsq1++;
				pvcomp1+=3;
			}//end of first neighbour cycle inb1
			pnneigh++;
			plist+= *max_neigh;
			ptypelist+= *max_neigh;
			pdistsq+= *max_neigh;
			pvec_comp+=3* *max_neigh;
		}//end of central atom cycle i
	}//end of iconst cycle
	CalcDistribution();//calculate the distribution
};

//-------------save to a file-----------------
void CosDistrConst::Save(const char *file_name) const
{
	int i,ibin,*phist;
	double *pcalc, *ptheor;
	ofstream file;

	if (strlen(file_name)!=0)
		mystrcpy(tempfilename,FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, cosfilename);
	OpenFile(file,tempfilename,"CosDistrConst::Save",0);//open file, check, whether it was successfully opened


	file.precision(6);
	file.setf(ios::fixed, ios::floatfield);
	file<<"Cosine distribution of bond angles"<<endl; 
	file<<nconstraints<<"\t number of constraints (nconstraints) "<<endl;
	file<<ncos_bin<<"\t number of bins"<<endl;
		
	file<<"\nConstraints properties: "<<endl;
	file<<"index, calculation method, central atom type, neighbour1 types, neighbour2 types "<<endl;
	file<<"inferior distance (1,2), superior distance (1,2), desired bond angle, "<<endl;
	file<<"control parameter for the shape, weighting parameter, width of the cos(theta) bin "<<endl;
	file<<"number of atoms satisfying the constraint "<<endl;
	
	for (i=0;i<nconstraints;i++)
	{
		file<<J4<<i+1<<"\t";
		switch (method[i])
		{
			case 0:
				file<<"Step function"<<"\t";
				break;
			case  1:
				file<<"Normal distr"<<"\t";
				break;
			case 2:
				file<<"No angles"<<"\t";
				break;
			case 3:
				file<<"Experimental distribution"<<"\t";
				break;
		}
			
		file<<J4<<central[i]+1<<"\t";//start with 1 in the file
		file<<J4<<neighbour1[i]+1<<"\t";//start with 1 in the file
		file<<J4<<neighbour2[i]+1<<"\t";//start with 1 in the file
		file<<J10<<dmin1[i]<<"\t";
		file<<J10<<dmin2[i]<<"\t";
		file<<J10<<dmax1[i]<<"\t";
		file<<J10<<dmax2[i]<<"\t";
		if (method[i]==3)
			file<<datafilename[i]<<"\t";
		else
		{
			file<<J10<<angle[i]<<"\t";
			file<<J10<<wcontrol[i]<<"\t";
		}
		file << J10 << weights[i] << "\t";
		file << J10 << dcos_theta[i] << endl;
	}	

	file<<endl;
	file<<"Cosine distribution of bond angles"<<endl;
	pcalc=cosinedistr;
	ptheor=theordistr;
	phist=cos_hist;
	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::fixed, ios::floatfield);
	file.precision(8);
	for (i=0;i<nconstraints;i++)
	{
		file<<"\n"<<i+1<<". constraint"<<endl;
		switch (method[i])
		{
			case 2:
				file<<"cos(theta)\tcalculated\tweight for chi2\tcalc_cos_hist"<<endl;
				break;
			case 3:
				file<<"cos(theta)\tcalculated\texperimental\tcalc_cos_hist"<<endl;
				break;
			default:
				file<<"cos(theta)\tcalculated\ttheoretical\tcalc_cos_hist"<<endl;
		}

		//to set the pointer to the beginning of the next data set
		pcalc+=first_bin[i];
		ptheor+=first_bin[i];
		phist+=first_bin[i];

		for (ibin=first_bin[i];ibin<last_bin[i];ibin++)
			file<<(ibin+0.5)*dcos_theta[i]-1<<"\t"<<*pcalc++<<"\t"<<*ptheor++<<"\t"<<*phist++<<endl;
		
		pcalc+=npoints[i]-last_bin[i];//to set the pointer to the next data set
		ptheor+=npoints[i]-last_bin[i];//to set the pointer to the next data set
		phist+=npoints[i]-last_bin[i];//to set the pointer to the next data set
		file<<endl;
	}
	file.unsetf(ios::right);
	file.unsetf(ios::fixed);
	file.precision(6);
	file.close();

};

//----updating the histogram, if sign switch is -1: removing contribution, if +1 adding contribution--------------------
void CosDistrConst::UpdateHist(Move &move, int sign_switch)
{
	int iconst,imoved,jmoved,inb1,inb2,ibin,icoord,icentral;
	int fixed_neigh,ncentral;
	int partialind;
	int neigh_type1,neigh_type2;//desired type of the neighbours 
	int coffset;
	int skipcycle;//boolean indicator
	int *moved_type,*moved_ind;
	int *pnneigh;//pointer to the number of neighbours 
	int *pindlist,*pind1,*pind2;//pointers to the indices of the neighbours
	int *ptypelist,*pneigh_type1,*pneigh_type2;//pointers to the type of the neighbours 
	int *pcos_hist;//pointer to the histogram of the calculated distribution

	double costheta;
	double dminsq1,dmaxsq1,dminsq2,dmaxsq2;//minimum and maximum squared distances for the two neighbours
	double *pdsq1,*pdsq2;//pointers to the distances of the neighbours
	double *pvcomp1,*pvcomp2;//pointers to the central-neighbour vector components

	if (sign_switch==-1)//to do it only once in each loop
	{
		for (iconst=0;iconst<nconstraints;iconst++)
			mod_const[iconst]=0;//default, not modified
	}
	moved_type=move.types;
	moved_ind=move.indices;
	for (imoved=0; imoved<Move::tot_moved_atoms;imoved++)
	{
		if (!(cosconsttype>>*moved_type & 1))
		{
			moved_type++;
			moved_ind++;
			continue;//this moved atom is not involved in any constraint, go to next
		}

		for (iconst=0;iconst<nconstraints;iconst++)
		{
			if (*moved_type==central[iconst])//the moved atom is central in this constraint
			{
				mod_const[iconst]=1;//this constraint is modified by the move 
				//setting auxiliary variables and pointers to increase the seed
				dminsq1=min1sq[iconst];
				dminsq2=min2sq[iconst];
				dmaxsq1=max1sq[iconst];
				dmaxsq2=max2sq[iconst];
				neigh_type1=neighbour1[iconst];
				neigh_type2=neighbour2[iconst];
				pcos_hist=cos_hist+cumul_bins[iconst];//beginning of the histogram for this contraint
				coffset=*moved_ind-config.cumul[*moved_type];//the offset of the central atom in its own type
				pnneigh=neighlist.nneigh_finder[*moved_type]+coffset;//sets the pointer to the nneigh of the central atom
		
				//setting the pointers for the first neighbour
				coffset*= *max_neigh;//multiplying the offset with the maximum number of neighbours
				pind1=neighlist.neighlist_finder[*moved_type]+coffset;//sets the pointer to the beginning of the neigh_list of the central atom
				pneigh_type1=neighlist.neightype_finder[*moved_type]+coffset;//sets the pointer to the beginning of the neigh_type of the central atom
				pdsq1=neighlist.dsq_finder[*moved_type]+coffset;//sets the pointer to the beginning of the squared distance of the central atom
				pvcomp1=neighlist.vector_finder[*moved_type]+coffset*3;//sets the pointer to the beginning of the central-neighbour1 vector components of the central atom
					
				//going through all the neighbours
				for (inb1=0;inb1<(*pnneigh-1);inb1++)//first neighbour cycle
				{
					//check, whether the neighbour is a moved atom situated later in the move.indices list
					skipcycle=0;
					for (jmoved=imoved+1;jmoved<Move::tot_moved_atoms;jmoved++)
					{
						if (move.indices[jmoved]==*pind1)
						{
							//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
							pneigh_type1++;
							pind1++;
							pdsq1++;
							pvcomp1+=3;
							skipcycle=1;
							break;
						}
					}
					if (skipcycle)
						continue;//go to the next inb1
					//determine, whether the first neighbour is from the right type, and is in the range
					if ((*pneigh_type1==neigh_type1 && *pdsq1>=dminsq1 && *pdsq1<=dmaxsq1) || \
						(*pneigh_type1==neigh_type2 && *pdsq1>=dminsq2 && *pdsq1<=dmaxsq2))
					{

						pneigh_type2=pneigh_type1+1;//type of the second neighbour
						pind2=pind1+1;//index of the second neighbour
						pdsq2=pdsq1+1;//squared distance of the central atom and the second neighbour 
						pvcomp2=pvcomp1+3;//vector components of the central-neighbour2 vector
						for (inb2=inb1+1;inb2<*pnneigh;inb2++)
						{
							//check, whether the neighbour is a moved atom situated later in the move.indices list
							skipcycle=0;
							for (jmoved=imoved+1;jmoved<Move::tot_moved_atoms;jmoved++)
							{
								if (move.indices[jmoved]==*pind2)
								{
									//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
									pneigh_type2++;
									pind2++;
									pdsq2++;
									pvcomp2+=3;
									skipcycle=1;
									break;
								}
							}
							if (skipcycle)
								continue;//go to the next inb2
							//determining the index of the partial calculated from the neighbours					
							partialind=(*pneigh_type1<=*pneigh_type2 ? (*pneigh_type1*RunParams::ntypes-(*pneigh_type1*(*pneigh_type1+1)/2)+*pneigh_type2) : \
								(*pneigh_type2*RunParams::ntypes-(*pneigh_type2*(*pneigh_type2+1)/2)+*pneigh_type1));
							if (partialind==neigh_partial[iconst])//both neighbours are from the right types
							{
								//determine, whether the second neighbour is in the range
								if ((*pneigh_type2==neigh_type1 && *pdsq2>=dminsq1 && *pdsq2<=dmaxsq1) || \
									(*pneigh_type2==neigh_type2 && *pdsq2>=dminsq2 && *pdsq2<=dmaxsq2))
								{
									//calculating the cosine of the angle of neighbour1-central-neighbour2
									costheta=0;
									for (icoord=0;icoord<3;icoord++)
										costheta+=*(pvcomp1+icoord) * *(pvcomp2+icoord);
									costheta/=sqrt(*pdsq1 * *pdsq2);//dividing by the length of the vectors
									//to avoid rounding errors, put back into the -1 - +1 range
									costheta=(costheta<-1 ? -1 : (costheta>1 ? 1 :costheta));
									ibin=int((costheta+1.)*inv_binwidth[iconst]);
									(*(pcos_hist+ibin))+=sign_switch;//update the count of his bin
								}//second neighbour is in the range
							}//both neighbours are from the right type
							
							pneigh_type2++;
							pind2++;
							pdsq2++;
							pvcomp2+=3;
						}//end of second neighbour cycle inb2
					}//first neighbour is from the right type, and is in the range

					pneigh_type1++;
					pind1++;
					pdsq1++;
					pvcomp1+=3;
				}//end of first neighbour cycle inb1
			}//the moved is central
			
			//Decide, whether the moved atom is a neighbour in the constraint 
			fixed_neigh=0;//the moved is not a neighbour in the constraint
			//setting the auxiliary variables the way, that index1 will belong to the fixed neighbour (the moved atom)
			if (*moved_type==neighbour1[iconst])//the moved atom is first neighbour in this constraint
			{
				mod_const[iconst]=1;//this constraint is modified by the move 
				fixed_neigh=1;
				neigh_type1=neighbour1[iconst];//this will be the fixed neighbour
				dminsq1=min1sq[iconst];
				dmaxsq1=max1sq[iconst];
				neigh_type2=neighbour2[iconst];//this will be the "running" neighbour
				dminsq2=min2sq[iconst];
				dmaxsq2=max2sq[iconst];
			}
			else
			{
				if (*moved_type==neighbour2[iconst])//the moved atom is second neighbour in this constraint
				{
					mod_const[iconst]=1;//this constraint is modified by the move 
					fixed_neigh=1;
					neigh_type1=neighbour2[iconst];//this will be the fixed neighbour
					dminsq1=min2sq[iconst];
					dmaxsq1=max2sq[iconst];
					neigh_type2=neighbour1[iconst];//this will be the "running" neighbour
					dminsq2=min1sq[iconst];
					dmaxsq2=max1sq[iconst];
				}
			}
			if (fixed_neigh)
			{
				mod_const[iconst]=1;//this constraint is modified by the move 
				pcos_hist=cos_hist+cumul_bins[iconst];//beginning of the histogram for this contraint
				//have to go through all the neighbours of the moved atom, considering them as central
				coffset=*moved_ind-config.cumul[*moved_type];//the offset of the moved atom in its own type
				ncentral=*(neighlist.nneigh_finder[*moved_type]+coffset);//sets the pointer to the nneigh of the moved atom
				coffset*= *max_neigh;
				pindlist=neighlist.neighlist_finder[*moved_type]+coffset;//sets the pointer to the beginning of the neigh_list of moved atom
				ptypelist=neighlist.neightype_finder[*moved_type]+coffset;//sets the pointer to the beginning of the neigh_type of the moved atom

				for (icentral=0;icentral<ncentral;icentral++)//going through the neighbours of the moved, this will be central
				{
					if (*ptypelist!=central[iconst])
					{
						//this is not the right type to be a central
						pindlist++;
						ptypelist++;
						continue;//go to next central
					}
					//check, whether the central is a moved atom situated later in the move.indices list
					skipcycle=0;
					for (jmoved=imoved+1;jmoved<Move::tot_moved_atoms;jmoved++)
					{
						if (move.indices[jmoved]==*pindlist)
						{
							//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
							ptypelist++;
							pindlist++;
							skipcycle=1;
							break;
						}
					}
					if (skipcycle)
						continue;//go to the next icentral
					coffset=*pindlist-config.cumul[*ptypelist];//the offset of the central in its own type
					pnneigh=neighlist.nneigh_finder[central[iconst]]+coffset;//sets the pointer to the nneigh of the central atom type 
					coffset*= *max_neigh;//multiplying the offset with the maximum number of neighbours

					//the first neighbour will be fixed, the moved atom
					//find the serial index of it among the neighbours of the central
					pind1=neighlist.neighlist_finder[central[iconst]]+coffset;//sets the pointer to the beginning of the neigh_list of the central atom
					for (inb1=0;inb1<(*pnneigh);inb1++)//first neighbour cycle
					{
						if (*pind1==*moved_ind)
							break;
						pind1++;
					}
					//the moved is the inb1-th neighbour for the central 
					pdsq1=neighlist.dsq_finder[central[iconst]]+coffset+inb1;//sets the pointer to the squared distance of the central and moved atom
					if (*pdsq1<dminsq1 || *pdsq1>dmaxsq1) 
					{
						pindlist++;
						ptypelist++;
						continue;//the fixed neighbour and the central is not in the range, go to next central
					}
						
					//setting the pointers for the second neighbour
					pind2=neighlist.neighlist_finder[central[iconst]]+coffset;//sets the pointer to the beginning of the neigh_list of the central atom
					pneigh_type2=neighlist.neightype_finder[central[iconst]]+coffset;//sets the pointer to the beginning of the neigh_type of the central atom
					pdsq2=neighlist.dsq_finder[central[iconst]]+coffset;//sets the pointer to the beginning of the squared distance of the central atom
					pvcomp2=neighlist.vector_finder[central[iconst]]+coffset*3;//sets the pointer to the beginning of the central-neighbour1 vector components of the central atom

					//setting the pointer for the fixed first neighbour, it is the inb1-th neighbour for the central 
					pvcomp1=pvcomp2+inb1*3;//sets the pointer to the beginning of the central-neighbour1 vector components for the moved atom
					for (inb2=0;inb2<*pnneigh;inb2++)//second neighbour cycle
					{
						if (inb2==inb1)
						{
							//this is the first neighbour again, skip
							pneigh_type2++;
							pind2++;
							pdsq2++;
							pvcomp2+=3;
							continue;
						}
						//check, whether the neighbour is a moved atom situated later in the move.indices list
						skipcycle=0;
						for (jmoved=imoved+1;jmoved<Move::tot_moved_atoms;jmoved++)
						{
							if (move.indices[jmoved]==*pind2)
							{
								//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
								pneigh_type2++;
								pind2++;
								pdsq2++;
								pvcomp2+=3;
								skipcycle=1;
								break;
							}
						}
						if (skipcycle)
							continue;//go to the next inb2
						
						//determine, whether the second neighbour is in the range
						if (*pneigh_type2==neigh_type2 && *pdsq2>=dminsq2 && *pdsq2<=dmaxsq2)
						{
							//calculating the cosine of the angle of neighbour1-central-neighbour2
							costheta=0;
							for (icoord=0;icoord<3;icoord++)
								costheta+=*(pvcomp1+icoord) * *(pvcomp2+icoord);
							costheta/=sqrt(*pdsq1 * *pdsq2);//dividing by the length of the vectors
							//to avoid rounding errors, put back into the -1 - +1 range
							costheta=(costheta<-1 ? -1 : (costheta>1 ? 1 :costheta));
							ibin=int((costheta+1.)*inv_binwidth[iconst]);
							(*(pcos_hist+ibin))+=sign_switch;//update the count of his bin
						}//second neighbour is the right type and in the range
						pneigh_type2++;
						pind2++;
						pdsq2++;
						pvcomp2+=3;
					}//end of second neighbour cycle inb2
					pindlist++;
					ptypelist++;
				}//end of icentral cycle	
			}//the moved is neighbour
		}//iconst cycle
		moved_type++;
		moved_ind++;
	}//moved atom cycle imoved

};

//-----------copying the modified histogram parts of the cosine distribution of the bond angles---------------
void CosDistrConst::CopyModified(CosDistrConst &target)
{
	int iconst,ibin;
	int *pi_source,*pi_target;//pointer to the source and target histogram 
	double *pd_source, *pd_target;//pointer to the source and target distribution, had to be copied, as it is 
				//recalculated only if the constraint was effected by the move 
	
	for (iconst=0;iconst<nconstraints;iconst++)
	{
		if (mod_const[iconst])
		{
			//the histogram of this constraint was modified, copy it 
			pi_source=cos_hist+cumul_bins[iconst];
			pi_target=target.cos_hist+cumul_bins[iconst];
			pd_source=cosinedistr+cumul_bins[iconst];
			pd_target=target.cosinedistr+cumul_bins[iconst];
			for (ibin=0;ibin<npoints[iconst];ibin++)
			{
				*pi_target++ = *pi_source++;
				*pd_target++ = *pd_source++;
			}
		}
	};
};
bool CosDistrConst::GetCosDistrConstFree(std::list<std::list<string> > pool )
//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
{
	//int i;
	//double dcos_theta_in,temp;
	nconstraints = RunParams::ncosdistr;
	consttype = &NeighbourList::consttype;

	//creating the static arrays
	SetArraysize(&neigh_partial, nconstraints, "neigh_partial", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&mod_const, nconstraints, "mod_const", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&first_bin, nconstraints, "first_bin", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&last_bin, nconstraints, "last_bin", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&central, nconstraints, "central", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&method, nconstraints, "method", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&neighbour1, nconstraints, "neighbour1", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&neighbour2, nconstraints, "neighbour2", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&dmin1, nconstraints, "dmin1", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&dmin2, nconstraints, "dmin2", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&dmax1, nconstraints, "dmax1", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&dmax2, nconstraints, "dmax2", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&angle, nconstraints, "angle", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&wcontrol, nconstraints, "wcontrol", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&weights, nconstraints, "weights", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&min1sq, nconstraints, "min1sq", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&min2sq, nconstraints, "min2sq", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&max1sq, nconstraints, "max1sq", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&max2sq, nconstraints, "max2sq", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&npoints, nconstraints, "npoints", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&dcos_theta, nconstraints, "dcos_theta", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&inv_binwidth, nconstraints, "inv_binwidth", "CosDistrConst::GetCosDistrConstFree");
	SetArraysize(&cumul_bins, nconstraints + 1, "cumul_bins", "CosDistrConst::GetCosDistrConstFree");

	datafilename = new char[nconstraints][FILE_NAME_SIZE];
	if (datafilename == NULL)
		NoArray("CosDistrConst::GetCosDistrConstFree", "datafilename");

	maxdistsq = 0.0;
	int i = 0;
	bool retval = true;
	bool onlyfile = true;//All datasets should read from file
	double dcostheta_in=-1;//it was not set
	


	for (std::list<std::list<string>>::iterator it2 = pool.begin(); it2 != pool.end() && onlyfile; it2++)
	{
		if ((std::find_if(it2->begin(), it2->end(), [](string c) { return (c.find(CheckKey("DISTRIB-TYPE")) != string::npos) && ((c.find("UNIFORM") != string::npos) || (c.find("GAUSSIAN") != string::npos) || (c.find("ABSENT") != string::npos)); }) != it2->end()) || (std::find_if(it2->begin(), it2->end(), [](string c) { return c.find("DISTRIB-TYPE") != string::npos; }) == it2->end()))
			onlyfile = false;
	}

	for (std::list<std::list<string>>::iterator it2 = pool.begin(); it2 != pool.end() && retval; it2++, i++)//going through all the cos constraints
	{
		
		//Check default values
		if (std::find_if(it2->begin(), it2->end(), [](string c) { return c.find(CheckKey("DISTRIB-TYPE")) != string::npos; }) == it2->end())//no distribution type is given
		{
			method[i] = 0;
			for (int j = 0; j < N_COS_METHOD; j++)
			{
				if (strcmp(cos_method[j], DISTRIB_TYPE_DEF) == 0)
				{
					method[i] = j;
					break;
				}
			}
			cout << "\nWARNING(" << ++warn << ") No DISTRIB-TYPE was found for the " << i + 1 << ". "<<CheckTag("COS")<<" constraint!" << endl << "\tThe DISTRIB-TYPE is set to " << DISTRIB_TYPE_DEF << endl;

		}
		else
			std::iter_swap(it2->begin(), std::find_if(it2->begin(), it2->end(), [](string c) { return c.find(CheckKey("DISTRIB-TYPE")) != string::npos; }));//Insert to the beginning


		int ce = 0, n = 0, si = 0, dt = 0, dd = 0, dw = 0, dc = 0;
		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)//going through the lines of the current constraint
		{

			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("DISTRIB-TYPE"))
			{
				if (dt++ > 0)
				{
					cout << "\nERROR: Multiple " << CheckKey("DISTRIB-TYPE") << " definition in the " << i + 1 << ". "<<CheckTag("COS")<<" constraint!" << endl;
					return false;
				}
				string value;
				if (!(ss >> value)) return ErrLackValue(*it);
				memset(datafilename[i], 0, sizeof(datafilename[i]));
				if ((value == CheckKeyValue("DISTRIB-TYPE", "UNIFORM")) || (value == CheckKeyValue("DISTRIB-TYPE", "GAUSSIAN")) || (value == CheckKeyValue("DISTRIB-TYPE", "ABSENT")))
				{
					for (int j = 0; j < N_COS_METHOD - 1; j++)
					{
						if (strcmp(value.c_str(), cos_method[j]) == 0)
						{
							method[i] = j;
							break;
						}
					}
				}
				else
				{
					method[i] = 3;
					std::copy(value.begin(), value.end(), datafilename[i]);
					
				}
			}
			else if (key == CheckKey("DCOSTH"))
			{
				if (dc++ > 0)
					cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("DCOSTH")<<" definition in the " << i + 1 << ". "<<CheckTag("COS")<<" constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (method[i] == 3)
				{
					if (::debug)
						cout << "\nWARNING(" << ++warn << "): The " << CheckKey("DCOSTH")<<" definition will be omitted in the " << i + 1 << ". "<<CheckTag("COS")<<" constraint, because DISTRIB-TYPE=filename has been set!" << endl;
				}
				else
				{
					double value;
					if (!(ss >> value)) return ErrLackValue(*it);
					if (value <= 0.0)
					{
						cout << "\nERROR: The " << CheckKey("DCOSTH")<< "value should be positive in the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl;
						return false;
					}
					else
						dcostheta_in = value;//will be writte over, if given again
				}
			}
			else if (key == CheckKey("DISTRIB-DEGREES"))
			{
				if (dd++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("DISTRIB-DEGREES")<<" definition in the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (method[i] == 3)
				{
					if (::debug)
						cout << "\nWARNING(" << ++warn << "): The " << CheckKey("DISTRIB-DEGREES")<<" definition will be omitted in the " << i + 1 << ". " << CheckTag("COS")<<" constraint, because DISTRIB-TYPE=filename has been set!" << endl;
				}
				else
				{
					if (!(ss >> angle[i])) return ErrLackValue(*it);
					if ((angle[i] < 0.0) || (angle[i] > 180.0))
					{
						cout << "\nERROR: The " << CheckKey("DISTRIB-DEGREES")<< "value (" << angle[i] << ") belonging to the " << i + 1 << ". " << CheckTag("COS")<<" constraint is out of valid range 0..180 degrees" << endl;
						return false;
					}
				}
			}
			else if (key == CheckKey("DISTRIB-WIDTH"))
			{
				if (dw++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("DISTRIB-WIDTH")<<" definition in the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (method[i] == 3)
				{
					if (::debug)
						cout << "\nWARNING(" << ++warn << "): " << CheckKey("The DISTRIB-WIDTH")<<" definition will be omitted in the " << i + 1 << ". " << CheckTag("COS")<<" constraint, because DISTRIB-TYPE = filename has been set!" << endl;
				}
				else
				{
					if (!(ss >> wcontrol[i])) return ErrLackValue(*it);
					
				}
			}
			else if (key == CheckKey("CENT-TYPE"))
			{
				if (ce++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("CENT-TYPE")<<" definition in the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (!(ss >> central[i])) return ErrLackValue(*it);
				if ((central[i] > RunParams::ntypes) || (central[i] < 1))
				{
					cout << "\nERROR: The " << CheckKey("CENT-TYPE")<<" value (" << central[i] << ") should be in the range: 1.." << RunParams::ntypes << " (number of types)!" << endl << "Please, correct it at the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl;
					return false;
				}
				central[i]--;
			}
			else if (key == CheckKey("NEIGH-TYPE_FROM_TO"))
			{
				if (n++ > 1)
				{
					cout << "\nERROR: "<<CheckKey("NEIGH-TYPE_FROM_TO")<<" is definied more than twice in the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl;
					return false;
				}
				int neighbours; double dmin, dmax;
				if (!GetInt(ss, *it, &neighbours,"CosDistrConst::GetCosDistrConstFree")) return ErrLackValue(*it);
				if (!(ss  >> dmin >> dmax)) return ErrLackValue(*it);
				if ((neighbours > RunParams::ntypes) || (neighbours < 1))
				{
					cout << "\nERROR: " << CheckKey("NEIGH-TYPE_FROM_TO") << " value (" << neighbours << ") should be in the range: 1.." << RunParams::ntypes << " (number of types)!" << endl << "Please, correct it at the " << i + 1 << ". " << CheckTag("COS")<<" constraint!" << endl;
					return false;
				}
				neighbours--;
				if (dmin > dmax)
				{
					cout << "\nWARNING(" << ++warn << "): the second value should be less than the third one in the following line" << endl << *it << endl;
					double temp = dmin;
					dmin = dmax;
					dmax = temp;
				}
				if (n == 1)
				{
					neighbour1[i] = neighbours;
					dmin1[i] = dmin;
					dmax1[i] = dmax;
					//calculating the reduced squared distances
					min1sq[i] = pow(dmin1[i] / RunParams::boxedge, 2.);
					max1sq[i] = pow(dmax1[i] / RunParams::boxedge, 2.);
					//finding the largest reduced maximum distance square
					if (max1sq[i] > maxdistsq) maxdistsq = max1sq[i];
				}
				else
				{
					neighbour2[i] = neighbours;
					dmin2[i] = dmin;
					dmax2[i] = dmax;
					//calculating the reduced squared distances
					min2sq[i] = pow(dmin2[i] / RunParams::boxedge, 2.);
					max2sq[i] = pow(dmax2[i] / RunParams::boxedge, 2.);
					//finding the largest reduced maximum distance square
					if (max2sq[i] > maxdistsq) maxdistsq = max2sq[i];
				}
			}
			else if ((key == CheckKey("SIGMA")) || (key == CheckKey("SIGMA-MASTER")) || (key == CheckKey("SIGMA-SCALABLE")))
			{
				if (si++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple SIGMA(-MASTER)/SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckTag("COS") << " constraint!" << endl << " The last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + i + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!(ss >> weights[i])) return ErrLackValue(*it);
				if (key == CheckKey("SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1)
						RunParams::lead_series_ind = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg+ RunParams::nek + i + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry is declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("SIGMA-SCALABLE"))
				{
					weights[i] = -fabs(weights[i]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					weights[i] = fabs(weights[i]);
				}
			}
			else
			{
				CheckKeyTag(key, CheckTag("COS"));//write the appropriate error mesage
				return false;
			}
		}
		if (retval)
		{
			//make some checks befor proceeding
			if (neighbour1[i] == neighbour2[i]  && fabs(dmin2[i] - dmin1[i]) > LOAD_TOL)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The neighbour types for the " << i + 1 << " cosine distarubution of angles constraint are the same, but the minimum distances are different!" << endl;
				cout << "dmin1 is " << dmin1[i] << " A and dmin2 is " << dmin2[i] << " A." << endl;
				cout << "This does not make sense, and would lead to inconsistency!" << endl;
				cout << "Check the " << datfilename << " file, and try again! Exiting..." << endl;
				CleanExit();
			}
			if (neighbour1[i] == neighbour2[i]  && fabs(dmax2[i] - dmax1[i]) > LOAD_TOL)
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The neighbour types for the " << i + 1 << " cosine distarubution of angles constraint are the same, but the maximum distances are different!" << endl;
				cout << "dmax1 is " << dmax1[i] << " A and dmax2 is " << dmax2[i] << " A." << endl;
				cout << "This does not make sense, and would lead to inconsistency!" << endl;
				cout << "Check the " << datfilename << " file, and try again! Exiting..." << endl;
				CleanExit();
			}
			//The concept is that the neighbourlist is calculated for atoms of each type represented in a cosine distribution
			//constraint. It is required for the neighbours as well, as the list is updated, and not completly recalculated in each step.
			//As each type can be involved in more, than one constraint, all the neighbour lists will be calculated with the
			//same maximum distance, maxdistsq. 
			//marking that the central and the neighbour type in cosconsttype types represented by right to left!
			//type 0 is represented by binary 1, type 1 binary 10 and so on  
			//consttype i a pointer to NeighbourList::consttype for all the constarins using the neighbour list
			//cosconsttype only for cos constraints
			if (!(*consttype >> central[i] & 1))//this type is not yet represented
			{
				*consttype += (longint)pow(2.0, (double)central[i]);
				NeighbourList::nconsttype++;//increase the number of constrained types
			}
			if (!(*consttype >> neighbour1[i] & 1))//this type is not yet represented
			{
				*consttype += (longint)pow(2.0, (double)neighbour1[i]);
				NeighbourList::nconsttype++;//increase the number of constrained types
			}
			if (!(*consttype >> neighbour2[i] & 1))//this type is not yet represented
			{
				*consttype += (longint)pow(2.0, (double)neighbour2[i]);
				NeighbourList::nconsttype++;//increase the number of constrained types
			}
			if (!(cosconsttype >> central[i] & 1))//this type is not yet represented
				cosconsttype += (longint)pow(2.0, (double)central[i]);

			if (!(cosconsttype >> neighbour1[i] & 1))//this type is not yet represented
				cosconsttype += (longint)pow(2.0, (double)neighbour1[i]);

			if (!(cosconsttype >> neighbour2[i] & 1))//this type is not yet represented
				cosconsttype += (longint)pow(2.0, (double)neighbour2[i]);

			//Determining the index of the partial calculated from the neighbour types 
			neigh_partial[i] = (neighbour1[i] <= neighbour2[i] ? (neighbour1[i] * RunParams::ntypes - (neighbour1[i] * (neighbour1[i] + 1) / 2) + neighbour2[i]) : \
				(neighbour2[i] * RunParams::ntypes - (neighbour2[i] * (neighbour2[i] + 1) / 2) + neighbour1[i]));

			mod_const[i] = 1;//to calculate the initial distribution
		}

	}
	double temp;
	if (fabs(dcostheta_in + 1.0) < TOLERANCE2)//no dcostheta was given at all, set the default
		temp = 2. / DCOSTHETA_DEF;//use the default value
	else
		temp = 2. / dcostheta_in;

	ncos_bin = (int)temp;
	if (temp - ncos_bin > 0.5)
		ncos_bin++;//to find the nearest integer
	for (i = 0; i < nconstraints; i++)
	{
		dcos_theta[i] = 2. / ncos_bin;
		inv_binwidth[i] = ncos_bin / 2.0;//inverse of the bin width
	}
	


	if (retval)
	{
		ntot_bins = nconstraints * ncos_bin;//by default, this can be increased later, if necessary due to 3-type constraint
		//creating the array for the theoretical distribution 
		SetArraysize(&theordistr, ntot_bins, "theordistr", "CosDistrConst::GetCosDistrConstFree");//theoretical distribution
	}
	return retval;
};
