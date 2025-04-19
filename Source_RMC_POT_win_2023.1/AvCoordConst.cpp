//source AvCoordConst.cpp
//Last changed 24.01.2023

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC //not to redefine the pointer to the intercation functions
#include "Threads.h" //classes1.h is included through Move.h!

//defining the static members
int  AvCoordConst::nconstraints;//the number of constraints
int  AvCoordConst::nthreads;//the total number of threads to use, set by RunParams::GetParams
int *AvCoordConst::central;//array of types of central particules
int *AvCoordConst::neighbours;//array of types of neighbour particles
int *AvCoordConst::avcctype;//array containing the partial index 

int *AvCoordConst::ncentral;//number of central atoms (in case of periodic only those, which are inside the R0-dmax)
double *AvCoordConst::dmin;//array of minimal distances 
double *AvCoordConst::dmax;//array of maximum distances 
double *AvCoordConst::weights;//array of weights for each av CC
double *AvCoordConst::acnreq;//array of average coordination numbers required
double *AvCoordConst::udminsq;//REDUCED min. dist. squared
double *AvCoordConst::udmaxsq;//REDUCED max. dist. squared	

#ifdef _NO_PERIODIC
double *AvCoordConst::dav_red;//array of average reduced distances 
#endif

//Move *AvCoordConst::move;//this cannot be initialized with the constructor, has to be assigned later		

//---------------------------constructor----------------
AvCoordConst::AvCoordConst()//default constructor
{
	int offset;
	if (nconstraints < 1)
	{
		if (::debug)
		{
			cout << "\nWARNING(" << ++warn << "): AvCoordConst constructor" << endl;
			cout << "\tThere are no average coordination constraints, or nconstraints was not initialised!" << endl;
		}
	}
	//Paddinng will be used, to make sure, that the segments of differet threads will not be cached together.
	offset = (nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*neighbourcount)) / sizeof(*neighbourcount)) : 0);//number of dummy integer elements
	SetArraysize(&neighbourcount,nconstraints*nthreads+(nthreads-1)*offset,"neighboucount","AvCoordConst::AvCoordConst");//array of sum of coord. numbers, for each thread to avoid mutexes locks
			//first for all the constraint for the main thread, then for all the constraint for the auxiliary threads in the order of increasing ID
	
	//setting all the array elements to 'zero' values
	InitNeighbourcount();
	
	
}

//-----------------------copy constructor-----------------------------
AvCoordConst::AvCoordConst(AvCoordConst &source)
{

	int i,offset;

	//Paddinng will be used, to make sure, that the segments of differet threads will not be cached together.
	offset = (nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*neighbourcount)) / sizeof(*neighbourcount)) : 0);//number of dummy integer elements
	SetArraysize(&neighbourcount, nconstraints * nthreads + (nthreads - 1) * offset, "neighboucount", "AvCoordConst::AvCoordConst");//array of sum of coord. numbers, for each thread to avoid mutexes locks
			//first for all the constraint for the main thread, then for all the constraint for the auxiliary threads in the order of increasing ID
	
		//copying the array elements from source for the main thread
	for (i=0;i< nconstraints * nthreads + (nthreads - 1) * offset;i++)
		neighbourcount[i]=source.neighbourcount[i];
}

//--------------Sets the element of neighbourcount array to 0--------------------------
void AvCoordConst::InitNeighbourcount()
{
	int i,offset;
	offset = (nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*neighbourcount)) / sizeof(*neighbourcount)) : 0);//number of dummy integer elements
	//setting all the array elements to 'zero' values
	for (i=0;i<nconstraints*nthreads+(nthreads-1)*offset;i++)
		neighbourcount[i]=0;
}

bool AvCoordConst::GetAvCoordConstFree(std::list<std::list<string> > pool)
//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
{
	nconstraints=RunParams::navcoord;//number of constraints

	//allocating memory
	SetArraysize(&central,nconstraints,"central","AvCoordConst::GetAvCoordConstFree");//array of central atom TYPES
	SetArraysize(&neighbours,nconstraints,"neighbours","AvCoordConst::GetAvCoordConstFree");//array of neighb. at. TYPES
		//all type indices start at 0!
	SetArraysize(&dmin,nconstraints,"dmin","AvCoordConst::GetAvCoordConstFree");//array on minimal dist.
	SetArraysize(&dmax,nconstraints,"dmax","AvCoordConst::GetAvCoordConstFree");//array of maxim. dist.
	SetArraysize(&weights,nconstraints,"weights","AvCoordConst::GetAvCoordConstFree");//array of weights
	SetArraysize(&acnreq,nconstraints,"acnreq","AvCoordConst::GetAvCoordConstFree");//array of required average coord. numb.
	SetArraysize(&udminsq,nconstraints,"udminsq","AvCoordConst::GetAvCoordConstFree");//REDUCED squared min. dist.
	SetArraysize(&udmaxsq,nconstraints,"udmaxsq","AvCoordConst::GetAvCoordConstFree");//REDUCED squared max. dist.
	SetArraysize(&avcctype,nconstraints,"avcctype","AvCoordConst::GetAvCoordConstFree");//array containing the partial index 
	SetArraysize(&ncentral,nconstraints,"ncentral","AvCoordConst::GetAvCoordConstFree");//array of number of central atoms
#ifdef _NO_PERIODIC
	SetArraysize(&dav_red,nconstraints,"dav_red","AvCoordConst::GetAvCoordConstFree");//array of average. red dist.
#endif

	int i=0;
	bool retval=true;
	for(std::list<std::list<string>>::iterator it2=pool.begin(); it2!=pool.end() && retval; it2++, i++)
	{
		int ce=0,n=0,co=0;
		for(std::list<string>::iterator it=it2->begin(); (it!=it2->end()) && retval; it++)
		{
	    
		string key,params;
		WrapFree(*it,key,params);
	    
		stringstream ss(params);
	  
		if(key== CheckKey("CENT-TYPE"))
		{
			if(ce++>0) cout<<"\nWARNING("<<++warn<<"): Multiple "<<CheckKey("CENT-TYPE")<<" definition in the "<<i+1<<". "<<CheckTag("AVCOORD")<<" constraint!"<<endl<<" The last value will be set as the relevant one."<<endl;
			if(!(ss>>central[i])) return ErrLackValue(*it);
			if((central[i]>RunParams::ntypes)||(central[i]<1))
			{
				cout<<"\nERROR: The " << CheckKey("CENT-TYPE")<<" value ("<<central[i]<<") should be in the range: 1.."<<RunParams::ntypes<<" (number of types)!"<<endl<<"Please, correct it at the "<<i+1<<". " << CheckTag("AVCOORD")<<" constraint!"<<endl;
				return false;
			}
				central[i]--;
				ncentral[i]=SimpleCfg::pnatoms[central[i]];//it will be reset in case of no_periodic
		}
		else if(key== CheckKey("NEIGH-TYPE_FROM_TO"))
		{
			if(n++>0) cout<<"\nWARNING("<<++warn<<"): Multiple "<<CheckKey("NEIGH-TYPE_FROM_TO")<<" definition in the "<<i+1<<". " << CheckTag("AVCOORD")<<" constraint!"<<endl<<" The last value will be set as the relevant one."<<endl;
			if (!GetInt(ss, *it, neighbours+i, "AvCoordConst::GetAvCoordConstFree"))  return ErrLackValue(*it);
			if (!(ss>>dmin[i]>>dmax[i])) return ErrLackValue(*it);
			if((neighbours[i]>RunParams::ntypes)||(neighbours[i]<1))
			{
				cout<<"\nERROR: "<<CheckKey("NEIGH-TYPE_FROM_TO")<< "key, the neighbour type value ("<<neighbours[i]<<") should be in the range: 1.."<<RunParams::ntypes<<" (number of types)!"<<endl<<"Please, correct it at the "<<i+1<<". " << CheckTag("AVCOORD")<<" constraint!"<<endl;
				return false;
			}
			neighbours[i]--;
			if(dmin[i]>dmax[i])
			{
				cout<<"\nWARNING("<<++warn<<"): " << CheckKey("NEIGH-TYPE_FROM_TO") << "key, the second value (for rmin) should be less than the third (rmax) in the following line. It will be swapped!"<<endl<<*it<<endl;
				double temp=dmin[i];
				dmin[i]=dmax[i];
				dmax[i]=temp;
			}
			if (RunParams::boxedge>0)
			{
				//Minimum and maximum distances for the constraints squared in reduced units
				udminsq[i]=dmin[i]*dmin[i]/RunParams::boxedge/RunParams::boxedge;
				udmaxsq[i]=dmax[i]*dmax[i]/RunParams::boxedge/RunParams::boxedge;
			}
		}
		else if((key== CheckKey("COORDNUM_SIGMA"))||(key== CheckKey("COORDNUM_SIGMA-MASTER"))||(key== CheckKey("COORDNUM_SIGMA-SCALABLE")))
		{
			if (co++ > 0)
			{
				cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("COORDNUM-SIGMA") << "(-MASTER)/" << CheckKey("COORDNUM_SIGMA-SCALABLE") << " definition in the " << i + 1 << ". " << CheckTag("AVCOORD") << " constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (RunParams::lead_series_ind == RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + CoordNumbConst::tot_subconst + i + 1)
					RunParams::lead_series_ind = -1;
			}
			if(!(ss>>acnreq[i]>>weights[i])) return ErrLackValue(*it);
			if(key== CheckKey("COORDNUM_SIGMA-MASTER"))
			{
				if(RunParams::lead_series_ind==-1) 
					RunParams::lead_series_ind=RunParams::ngr+RunParams::nsq+RunParams::nfq+RunParams::nfg+RunParams::nek+RunParams::ncosdistr+CoordNumbConst::tot_subconst+i+1;
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
			else if(key== CheckKey("COORDNUM_SIGMA-SCALABLE"))
			{
				weights[i]=-fabs(weights[i]);
				ChiSquared::calc_sigma=1;
			}
			else //normal, should not be negative, it does not effect calculation, but would upset *.free...
			{
				weights[i] = fabs(weights[i]);
			}
		}
		else
		{
			CheckKeyTag(key, CheckTag("AVCOORD"));//write the appropriate error mesage
			return false;
		}
		}
		if(retval)
		{
			//Calculating the index of the AvCoordConst::central[iconst],AvCoordConst::neighbours[iconst]-th
			//partial.
			if (central[i]<=neighbours[i])
				avcctype[i]=central[i]*RunParams::ntypes-central[i]*(central[i]+1)/2+neighbours[i];
			else
				avcctype[i]=neighbours[i]*RunParams::ntypes-neighbours[i]*(neighbours[i]+1)/2+central[i];
		}
	}
	
	return retval;
}

//-------Getting the constraint parameters from the .dat file, and sets the static members-----------
void AvCoordConst::GetAvCoordConst(ifstream &file)	
{
	int i;
	char *name, *conv_numb = NULL;

	nconstraints=RunParams::navcoord;//number of constraints
	
	//Checking the number of constraints
	
	    //Checking the main parameters
	if (nconstraints<1)
	{
		cout << "\nNOTE(" << ++note << "): AvCoordConst::GetAvCoordConstThere are no average coordination constraints to read!"<<endl;
		return;
	}

	if (CheckFileState(file,"AvCoordConst::GetAvCoordConst",datfilename)==0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	//allocating memory
	SetArraysize(&name, NAME_SIZE, "name", "CosDistrConst::GetCosDistrConst");
	SetArraysize(&central,nconstraints,"central","AvCoordConst::GetAvCoordConst");//array of central atom TYPES
	SetArraysize(&neighbours,nconstraints,"neighbours","AvCoordConst::GetAvCoordConst");//array of neighb. at. TYPES
		//all type indices start at 0!
	SetArraysize(&dmin,nconstraints,"dmin","AvCoordConst::GetAvCoordConst");//array on minimal dist.
	SetArraysize(&dmax,nconstraints,"dmax","AvCoordConst::GetAvCoordConst");//array of maxim. dist.
	SetArraysize(&weights,nconstraints,"weights","AvCoordConst::GetAvCoordConst");//array of weights
	SetArraysize(&acnreq,nconstraints,"acnreq","AvCoordConst::GetAvCoordConst");//array of required average coord. numb.
	SetArraysize(&udminsq,nconstraints,"udminsq","AvCoordConst::GetAvCoordConst");//REDUCED squared min. dist.
	SetArraysize(&udmaxsq,nconstraints,"udmaxsq","AvCoordConst::GetAvCoordConst");//REDUCED squared max. dist.
	SetArraysize(&avcctype,nconstraints,"avcctype","AvCoordConst::GetAvCoordConst");//array containing the partial index 
	SetArraysize(&ncentral,nconstraints,"ncentral","AvCoordConst::GetAvCoordConst");//array of number of central atoms

#ifdef _NO_PERIODIC
	SetArraysize(&dav_red,nconstraints,"dav_red","AvCoordConst::GetAvCoordConst");//array of average. red dist.

#endif

	//reading the data
	for (i=0;i<nconstraints;i++)
	{
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "central[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		central[i]=ReadThisLine(file,3,1,name, "AvCoordConst::GetAvCoordConst");//central atoms index (Starting with 1)
		RunParams::CheckType(central[i], name, "AvCoordConst::GetAvCoordConst");
		central[i]-=1;//Resetting to start with 0 in program
		ncentral[i]=SimpleCfg::pnatoms[central[i]];//it will be reset in case of no_periodic
		
		mystrcpy(name, NAME_SIZE, "neighbours[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		neighbours[i]=ReadThisLine(file,3,1,name, "AvCoordConst::GetAvCoordConst");//neighbours  atoms index (Starting with 1)
		RunParams::CheckType(neighbours[i],name, "AvCoordConst::GetAvCoordConst");
		neighbours[i]-=1;//Resetting to start with 0 in program
		mystrcpy(name, NAME_SIZE, "dmin[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmin[i]=ReadThisLine(file,3,1.0,name, "AvCoordConst::GetAvCoordConst");//minimal distance
		
		mystrcpy(name, NAME_SIZE, "dmax[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmax[i]=ReadThisLine(file,3,1.0,name, "AvCoordConst::GetAvCoordConst");//maximum distance
		
		mystrcpy(name, NAME_SIZE, "target coord numb[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		acnreq[i]=ReadThisLine(file,3,1.0,name, "AvCoordConst::GetAvCoordConst");//target coordination number
		mystrcpy(name, NAME_SIZE, "sigma");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		weights[i]=ReadThisLine(file,1,1.0,name, "AvCoordConst::GetAvCoordConst",datfilename);//weight
		if (weights[i]<0)
			ChiSquared::calc_sigma=1;

		//Checking, whether loading was successful
		if (!CheckReadFileState(file,"AvCoordConst::GetAvCoordConst",datfilename))
			CleanExit();//loading failed
		
		//Initialising the remainder of the static members
		if (RunParams::boxedge>0)
		{
			//Minimum and maximum distances for the constraints squared in reduced units
			udminsq[i]=dmin[i]*dmin[i]/RunParams::boxedge/RunParams::boxedge;
			udmaxsq[i]=dmax[i]*dmax[i]/RunParams::boxedge/RunParams::boxedge;
		}
		
		//Calculating the index of the AvCoordConst::central[iconst],AvCoordConst::neighbours[iconst]-th
		//partial.
		if (central[i]<=neighbours[i])
			avcctype[i]=central[i]*RunParams::ntypes-central[i]*(central[i]+1)/2+neighbours[i];
		else
			avcctype[i]=neighbours[i]*RunParams::ntypes-neighbours[i]*(neighbours[i]+1)/2+central[i];
			
	}
}

//============load============
int AvCoordConst::Load(double *sigma_percentage)
{//load all the data from the file given in argument
	int i,dumb;
	
	double real_dumb;
	ifstream file;

	SafeOpenTextFile(file,acnfilename);
	if (CheckFileState(file,"AvCoordConst::Load",acnfilename)==0)
	{
		cout << "\nWARNING(" << ++warn << "): There is no open file to load the AvCoordConst object from"<<endl;
		cout << "\toading was not successful, histogram";
		cout << "\tLoading was not successful, histogram";
		if (RunParams::nbvs > 0)
			cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
		else
			cout << " and coordination numbers will be calculated!" << endl;
		return(0);//loading was not successful
	}
	
	SkipLine(file,acnfilename,1);//comment
	SkipLine(file,acnfilename,1);
	cout<<"\nLoading the average coordination constraint object!"<<endl;
	dumb=ReadThisLine(file,1,0,"nconstraints", "AvCoordConst::Load",acnfilename);//number of constraints
	if (dumb!=nconstraints)//checking consistency for number of constraints
	{
		cout << "\nWARNING(" << ++warn << "): in AvCoordConst::Load" << endl;
		cout << "\tThe number of constraint is " << nconstraints << " from the .dat file" << endl;
		cout<<"\tand "<<dumb<<" from the .acn file!"<<endl;
		cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
		return(0);//loading was not successful
	}

	for (i=0;i<4;i++)
		SkipLine(file,acnfilename,1);//skip the comments and blanks
	
	for(i=0;i<nconstraints;i++)
	{
		dumb=ReadThisLine(file,3,1,"index", "AvCoordConst::Load");//index
	
		dumb=ReadThisLine(file,3,1,"central type", "AvCoordConst::Load");//type of central atom (starts with 1 in file!)
		if (central[i]!=(dumb-1))//checking consistency for type of central atom
		{
			cout << "\nWARNING(" << ++warn << "): AvCoordConst::Load"<<endl;
			cout<<"\tThe type of the central atom for "<<i+1<<". constraint is "<<central[i]+1<<" from the .dat file,"<<endl;
			cout<<"\tand "<<dumb<<" from the .acn file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
			return(0);//loading was not successful
		}

		dumb=ReadThisLine(file,3,1,"neighbour type", "AvCoordConst::Load");//type of neighbour atom
		if (neighbours[i]!=(dumb-1))//checking consistency for type of neighbour atom
		{
			cout << "\nWARNING(" << ++warn << "): AvCoordConst::Load"<<endl;
			cout<<"\tThe type of the neighbour atom for the "<<i+1<<". constraint is "<<neighbours[i]+1<<" from the .dat file,"<<endl;
			cout<<"\tand "<<dumb<<" from the .acn file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;	
			return(0);//loading was not successful
		};

		real_dumb=ReadThisLine(file,3,1.0,"minimum distance", "AvCoordConst::Load");//inferior distance
		if (fabs(dmin[i]-real_dumb)>LOAD_TOL)//checking consistency
		{
			cout.precision(12);
			cout.setf(ios::right, ios::adjustfield);
			cout.setf(ios::scientific, ios::floatfield);
			cout << "\nWARNING(" << ++warn << "): AvCoordConst::Load"<<endl;
			cout<<"\tInferior distance for the "<<i+1<<". constraint is "<<dmin[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<real_dumb<<" from the .acn file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
			cout.precision(6);
			cout.unsetf(ios::scientific);
			cout.unsetf(ios::right);
			return(0);//loading was not successful
		};

		real_dumb=ReadThisLine(file,3,1.0,"maximum distance", "AvCoordConst::Load");//superior distance
		if (fabs(dmax[i]-real_dumb)>LOAD_TOL)//checking consistency
		{
			cout.precision(12);
			cout.setf(ios::right, ios::scientific);
			cout << "\nWARNING(" << ++warn << "): AvCoordConst::Load"<<endl;
			cout<<"\tSuperior distance for the "<<i+1<<". constraint is "<<dmax[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<real_dumb<<" from the .acn file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
			cout.precision(6);
			cout.unsetf(ios::scientific);
			cout.unsetf(ios::right);
			return(0);//loading was not successful
		};
		
		real_dumb=ReadThisLine(file,3,1.0,"desired coordination number", "AvCoordConst::Load");//Desired coordination number
		if (acnreq[i]!=real_dumb)//checking consistency
		{
			cout << "\nWARNING(" << ++warn << "): AvCoordConst::Load"<<endl;
			cout<<"\tDesired average coordination number for the "<<i+1<<". constraint is "<<acnreq[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<real_dumb<<" from the .acn file!"<<endl;
			cout<<"\tThe value in the .dat file will be used!"<<endl;
		};

		real_dumb=ReadThisLine(file,1,1.0,"sigma", "AvCoordConst::Load",acnfilename);//weight parameter
		int index = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + ChiSquared::tot_cc_subconst;
		if (weights[i]!=real_dumb && sigma_percentage[index+i]>0)//checking consistency
		{
			cout << "\nWARNING(" << ++warn << "): AvCoordConst::Load"<<endl;
			cout<<"\tWeight parameter for the "<<i+1<<". constraint is "<<weights[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<real_dumb<<" from the .acn file!"<<endl;
			cout<<"\tThe value in the .dat file will be used!"<<endl;
		}

	}//end of cycle through the constraint

	SkipLine(file,acnfilename,1);//skip the comments and blanks
	SkipLine(file,acnfilename,1);

	for(i=0;i<nconstraints;i++)//for each constraint
	{
		SkipLine(file,acnfilename,1);//blank
		SkipLine(file);//constraint's index
		if (i<nconstraints-1)
			*(neighbourcount+i)=ReadThisLine(file,1,1,"neighbourcount", "AvCoordConst::Load",acnfilename);//number of neighbours
		else
			*(neighbourcount+i)=ReadThisLine(file,1,1,"neighbourcount", "AvCoordConst::Load");//number of neighbours
	}

	//Checking, whether loading was successful
	if (CheckReadFileState(file,"AvCoordConst::Load",acnfilename))
		return(1);//load was successful
	else
	{
		cout << "\nWARNING(" << ++warn << "): Loading was not successful, histogram and coordination numbers will be calculated!"<<endl;
		return(0);//loading was not successful
	}
	next_line_pos = -1;
	file.close();
	
	return (1);
}

//============save the average coordintion constraint============
void AvCoordConst::Save(const char *file_name) const
{
	int i;
	ofstream file;

	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10,file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10,acnfilename);
	OpenFile(file,tempfilename,"AvCoordConst::Save",0);//open file, check, whether it was successfully opened

	file<<"Data for Average Coordination Constraints. "\
	<<"file created by AvCoordConst::save \n"<<endl;
	
	file<<nconstraints<<"\t number of constraints (nconstraints)"<<endl;
	file<<"Constraints properties: (RMC-type indices start at 1) "<<endl;
	file<<"index, central atom type, neighbour atom type,"<<endl;
	file<<"inferior distance, superior distance,"<<endl;
	file<<"desired average coordination number, weighting parameter "<<endl;
	
	file.setf(ios::fixed, ios::floatfield);
	file.setf(ios::right, ios::adjustfield);

	//Saving the constraint's parameters
	for(i=0;i<nconstraints;i++)
	{
		file<<J4<<i+1<<"\t";
		file<<J4<<central[i]+1<<"\t";
		file<<J4<<neighbours[i]+1<<"\t";
		file.precision(12);
		file<<J10<<dmin[i]<<"\t";
		file<<J10<<dmax[i]<<"\t";
		file.precision(6);
		file<<J10<<acnreq[i]<<"\t ";
		file<<J10<<weights[i]<<"\t ";
		file<<endl;
	}

	file.unsetf(ios::fixed);
	file.unsetf(ios::right);
	file<<endl;
	file<<"Number of neighbours for the configuration per constraint\n"<<endl;
	
	for(i=0;i<nconstraints;i++)//for each constraint
	{
		file<<i+1<<"-th constraint";
		file<<endl;
		
		file<<J7<<*(neighbourcount+i)<<endl;//number of neighbours
		
		file<<endl;
	}

	file.close();
}
//------Copy the neighbourcount values for the constraints affected--------- 
//--------------by the move from source to target-----------------------
void AvCoordConst::CopyModified(AvCoordConst &target)
{
	int iconst;

	for (iconst=0;iconst<nconstraints;iconst++)
	{
		if (Threads::mod_partial[avcctype[iconst]])
		{
			//Average coordination numbers
			*(target.neighbourcount+iconst)=*(neighbourcount+iconst);
									
		}//end of if this constraint was affected by the move
	}//end of constraint cycle iconst
};
