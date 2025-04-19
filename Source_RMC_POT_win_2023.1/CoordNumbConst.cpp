//source CoordNumbConst.cpp
//Last changed 24.01.2023

//The constraint was generalized. 
//(1)	Each constraint can have more, than one neighbour type. This way it is possible to
//		set up contstraints like the central atom has 2 neighbour types, and should have for example 3 neighbours between 
//		rmin[neigh_type1]-rmax[neigh_type1] for the neighbour type type1 and rmin[neigh_type2]-rmax[neigh_type2] for the 
//		neighbour type type2
//(2)	More than one subconstraint can be given for a constraint. More tha one subconstraint means, that more than one desired coordination
//		number (each with its own desired fraction and sigma) can be specified for a constraint. Subconstraints were introduced 
//		to make the calculation quicker and use less memory and disk space, if constraints should only differ in their desired
//		coordination number. The introduction of subconstraint will not result in additional functionality.


//CHANGE IN THE *.dat file:
//number of constraint; NUMBER OF NEIGHBOUR TYPES/constraint; NUMBER of SUBCONSTRAINTS/constraint
//for each constraint:type of cenral; neighbour type_first;....;neighbour type_last;rmin[first neigh type]...rmin[last neigh type];
//	rmax[first neigh type]...rmax[last neigh type]; target_coord[first subconst]...target_coord[last subconst]; 
//	fraction[first subconst]...fraction[last subconst];weight[first subconst]...weight[last subconst]

//The program is capable of reaing the old format, if no integer values (or not enough for each constraint) is given 
//after the number of constraint for NUMBER OF NEIGHBOUR TYPES, then it assumes 1 neighbour type for each 
//constraint. In this case does not even look for the number of subconstraints, assumes 1 for each constraint!

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h" //classes1.h is included through Move.h!

//Defining the static members
int  CoordNumbConst::nconstraints;//number of constraints
int  CoordNumbConst::nthreads;//the total number of threads to use, set by RunParams::GetParams
int  CoordNumbConst::ncentral_tot;//total number of central atoms for all the constraints
int  CoordNumbConst::tot_neightype;//total number of neighbour types
int  CoordNumbConst::tot_subconst;//total number of subconstraints
int *CoordNumbConst::central;//array of types of central particles
int *CoordNumbConst::n_neightype;//number of neighbour types /constraint
int *CoordNumbConst::neighbours=NULL;//array of types of neighbour particles
int *CoordNumbConst::target_coord;//array of DESIRED coord. numbers
int *CoordNumbConst::ncentral;//number of central atoms in the cfg
int *CoordNumbConst::n_subconst;//number of subconstraint/constraint
int *CoordNumbConst::cum_n_neightype;//cumulative number of neighbour types
int *CoordNumbConst::ncentral_cum;//cumulative number of central atoms before the given constraint for each constraint
int *CoordNumbConst::cum_n_subconst;//number of subconstraint/constraint
longint *CoordNumbConst::cctype;//array containing the partial indices in binary format

double *CoordNumbConst::dmin=NULL;//array of minimal distances 
double *CoordNumbConst::dmax=NULL;//array of maximum distances 
double *CoordNumbConst::weights;//array of weights for each CC 
double *CoordNumbConst::fraction;//array of fraction of atoms desired with the coord. number
double *CoordNumbConst::udminsq=NULL;//array of REDUCED min. dist. squared
double *CoordNumbConst::udmaxsq=NULL;//array of REDUCED max. dist. squared
double CoordNumbConst::maxdistsq=0;//the maximum distance among the udmaxsq arrays,

//if details are written, a separate neighbour list has to be calculated
int CoordNumbConst::max_neigh;
int CoordNumbConst::ncells_nlist;
int CoordNumbConst::overlap_nlist;
bool *CoordNumbConst::write_detail = NULL;//array of whetehr to write detailed cnc information to cncd file
bool CoordNumbConst::is_detail = 0;//whether  there is detail for any constraint
//these are only for the actual central atom's neighbour, not stored, only saved
int *CoordNumbConst::neigh_list;//neighbour indices for each central of each constraint
int *CoordNumbConst::neigh_type;//neighbour types for each central of each constraint
double *CoordNumbConst::dist;//central-neighbour distance in A for each central of each constraint
int **CoordNumbConst::neigh_finder;//finder of neigh_list for each constraint
int **CoordNumbConst::type_finder;//finder of type_list for each constraint
double **CoordNumbConst::dist_finder;//finder of dist_list for each constraint

//Move *CoordNumbConst::move;//this cannot be initialized with the constructor, has to be assigned later		
//constructor
CoordNumbConst::CoordNumbConst(NeighbourList &neigh):neighlist(neigh)//default constructor
{
	int i;
	int offset;
	if (nconstraints < 1)
	{
		if (::debug)
		{
			cout << "\nWARNING(" << ++warn << "): CoordNumbConst constructor" << endl;
			cout << "\tThere are no coordination constraints, or nconstraints was not initialised!" << endl;
		}
	}
	//allocating memory
	//There will be a segment for all the threads, used at histogram change recalculation.
	//Paddinng will be used, to make sure, that the segments of differet threads will not be cached together.
	offset=(tot_subconst>0 ? (int)((CACHE_PADDING-sizeof(*nsatisfy))/sizeof(*nsatisfy)) : 0);//number of dummy integer elements
	SetArraysize(&nsatisfy,tot_subconst*nthreads+(nthreads-1)*offset,"nsatisfy","CoordNumbConst::CoordNumbConst");
	SetArraysize(&coordnumbs,ncentral_tot*nthreads + (nthreads - 1) * offset,"coordnumbs","CoordNumbConst::CoordNumbConst");

	SetArraysize(&finder,nconstraints,"finder","CoordNumbConst::CoordNumbConst");
	SetArraysize(&thread_finder,nthreads,"thread_finder","CoordNumbConst::CoordNumbConst");//for the beginning of each thread's segment
		
	InitCoordnumbs();//Sets the array elements to zero

	//initialising the finder
	thread_finder[nthreads-1]=coordnumbs;//the main has the first segment
	for(i=0;i<nthreads-1;i++)	
		thread_finder[i]=coordnumbs+(i+1)*(ncentral_tot+offset);
	

	for(i=0;i<nconstraints;i++)	
		finder[i]=coordnumbs+ncentral_cum[i];
	
}

//------------------- copy constructor---------------------------
CoordNumbConst::CoordNumbConst(NeighbourList &neigh, CoordNumbConst &source):neighlist(neigh)
{
	int i;
	int offset;
		
	//allocating memory
	//There will be a segment for all the threads, used at histogram change recalculation.
	//Paddinng will be used, to make sure, that the segments of differet threads will not be cached together.
	offset=(tot_subconst>0 ? (int)((CACHE_PADDING-sizeof(*nsatisfy))/sizeof(*nsatisfy)) : 0);//number of dummy integer elements
	SetArraysize(&nsatisfy,tot_subconst*nthreads+(nthreads-1)*offset,"nsatisfy","CoordNumbConst::CoordNumbConst");
	SetArraysize(&coordnumbs,ncentral_tot*nthreads+ (nthreads - 1) * offset,"coordnumbs","CoordNumbConst::CoordNumbConst");

	SetArraysize(&finder,nconstraints,"finder","CoordNumbConst::CoordNumbConst");
	SetArraysize(&thread_finder,nthreads,"thread_finder","CoordNumbConst::CoordNumbConst");//for the beginning of each thread's segment
		
	//Copying the main thread's array elements from source
	for(i=0;i<ncentral_tot;i++)
		coordnumbs[i]=source.coordnumbs[i];

	for(i=0;i<tot_subconst*nthreads+(nthreads-1)*offset;i++)	
		nsatisfy[i]=source.nsatisfy[i];
	
	//initialising the finder
	thread_finder[nthreads-1]=coordnumbs;//the main has the first segment
	for(i=0;i<nthreads-1;i++)	
		thread_finder[i] = coordnumbs + (i + 1) * (ncentral_tot + offset);
	

	for(i=0;i<nconstraints;i++)	
		finder[i]=coordnumbs+ncentral_cum[i];

}

void CoordNumbConst::SetParams()
{
	if (CoordNumbConst::is_detail)//only needed, if details for cnc is calculated, as the normal calculation is done during hist calc
	{
		int safe_add = SAFE_ADD;
		double temp,xmax;
		xmax=CoordNumbConst::maxdistsq;//reduced maximum distance for the neighbourlist 
		temp = int(SimpleCfg::ntotal / 8 * 4 / 3 * pow(xmax, 3.) * PI) + 1;
		//for safety reason increase the number of max_neigh with 10 percent (or at least SAFE_ADD), whichever is more
		max_neigh = (int)(temp + (temp * 0.1 > safe_add ? temp * 0.1 : safe_add));
		ncells_nlist = int(xmax / NeighbourList::cellwidth) + 1;//number of grid cells to test for each direction for the neighbourlist
		overlap_nlist = (ncells_nlist * 2 + 1 - NeighbourList::ngridcells > 0 ? ncells_nlist * 2 + 1 - NeighbourList::ngridcells : 0);

			
		//neigh_list is calculated for one central atom at each time, but not stored after saving
		SetArraysize(&neigh_list, max_neigh, "neigh_list", "CoordNumbConst::SetParams");
		SetArraysize(&neigh_type, max_neigh, "neigh_type", "CoordNumbConst::SetParams");
		SetArraysize(&dist, max_neigh, "dist", "CoordNumbConst::SetParams");

	}

}

//------------Sets the array elements to zero------------------
void CoordNumbConst::InitCoordnumbs()
{
	int i,size1,offset;
	size1=0;
	//Paddinng will be used, to make sure, that the segments of differet threads will not be cached together.
	offset = (tot_subconst > 0 ? (int)((CACHE_PADDING - sizeof(*nsatisfy)) / sizeof(*nsatisfy)) : 0);//number of dummy integer elements

	//Total number of central atoms
	for (i=0;i<nconstraints;i++)
		size1+=ncentral[i];
	size1 *= nthreads;
	size1+=(nthreads - 1) * offset;//each thread has its own segment

	//setting all the array elements to 'zero' values
	for(i=0;i<size1;i++)
		coordnumbs[i]=0;
	
	for (i=0;i<tot_subconst*nthreads+(nthreads-1)*offset;i++)
		nsatisfy[i]=0;
	

};	
//--------------------------------------------------------------------------------
//gets or sets the value for the iconst-th constraints, iatom-th 
//central atom in the coordnumbs array
//iatom index should be in its own type
int &CoordNumbConst::CoordinationNumb(int iconst, int iatom)
{
	return(*(finder[iconst]+iatom));
}


	
//----------------Getting the parameters from the.dat file, and sets the static members-------------------
void CoordNumbConst::GetCoordConst(ifstream &file)
{
	int i,ineightype,isubconst;
	unsigned int npartials=RunParams::ntypes*(RunParams::ntypes+1)/2;
	char *name, *conv_numb = NULL;
	char *conv_numb2=NULL;
	nconstraints=RunParams::nicoord;

	//Checking the main parameters
	if (nconstraints<1)
	{
		cout<<endl;
		cout<<"There are no coordination constraints to read!"<<endl;
		return;
	}

	if (npartials>sizeof (*cctype)*8)
	{
		cout << "\n****ERROR****"<<endl;
		cout<<"CoordNumbConst::GetCoordConst: The program can handle only "<<sizeof (*cctype)*8<<" partials presently!"<<endl;
		cout<<"Change the type of CoordNumbConst::cctype if possible to have more bytes!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if (CheckFileState(file,"CoordNumbConst::GetCoordConst",datfilename)==0)
	{
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	//allocating memory
	SetArraysize(&name, NAME_SIZE, "name", "CosDistrConst::GetCosDistrConst");
	SetArraysize(&n_neightype,nconstraints,"n_neightype","CoordNumbConst::GetCoordConst");
	SetArraysize(&n_subconst,nconstraints,"n_subconst","CoordNumbConst::GetCoordConst");
	SetArraysize(&central,nconstraints,"central","CoordNumbConst::GetCoordConst");
	SetArraysize(&ncentral,nconstraints,"ncentral","CoordNumbConst::GetCoordConst");
	SetArraysize(&cctype,nconstraints,"cctype","CoordNumbConst::GetCoordConst");
	SetArraysize(&write_detail, nconstraints, "write_detail", "CoordNumbConst::GetCoordConst");
	SetArraysize(&cum_n_neightype,nconstraints+1,"cum_n_neightype","CoordNumbConst::GetCoordConst");
	SetArraysize(&cum_n_subconst,nconstraints+1,"cum_n_subconst","CoordNumbConst::GetCoordConst");
	SetArraysize(&ncentral_cum,nconstraints+1,"ncentral_cum","CoordNumbConst::GetCoordConst");
	
	cum_n_subconst[0]=0;//cumulative number of subconstraints before the given constraint
	cum_n_neightype[0]=0;//cumulative number of neighbour type before the given constraint
	maxdistsq = 0;
	//reading the number of neighbour types/constraint
	for (i=0;i<nconstraints;i++)
	{
		write_detail[i] = 0;//only used in free format
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "number of neighbour types[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		n_neightype[i]=ReadThisLine(file,8,1,name, "CoordNumbConst::GetCoordConst");
		
		tot_neightype+=n_neightype[i];//total number of neighbour types
		cum_n_neightype[i+1]=cum_n_neightype[i]+n_neightype[i];//cumulative number of neighbour type before the given constraint
	}
	
	//reading the number of subconstraint/constraint
	for (i=0;i<nconstraints;i++)
	{
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "number of subconstaints/type[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		n_subconst[i]= ReadThisLine(file, 8, 1, name, "CoordNumbConst::GetCoordConst");
			
		tot_subconst+=n_subconst[i];//total number of subconstraints
		cum_n_subconst[i+1]=cum_n_subconst[i]+n_subconst[i];//cumulative number of subconstraint before the given constraint
	}
		
	
	SkipLine(file,datfilename,1);
	//allocating memory
	SetArraysize(&neighbours,tot_neightype,"neighbours","CoordNumbConst::GetCoordConst");
	SetArraysize(&dmin,tot_neightype,"dmin","CoordNumbConst::GetCoordConst");
	SetArraysize(&dmax,tot_neightype,"dmax","CoordNumbConst::GetCoordConst");
	SetArraysize(&udminsq,tot_neightype,"udminsq","CoordNumbConst::GetCoordConst");
	SetArraysize(&udmaxsq,tot_neightype,"udmaxsq","CoordNumbConst::GetCoordConst");
	SetArraysize(&target_coord,tot_subconst,"target_coord","CoordNumbConst::GetCoordConst");
	SetArraysize(&fraction,tot_subconst,"fraction","CoordNumbConst::GetCoordConst");
	SetArraysize(&weights,tot_subconst,"weights","CoordNumbConst::GetCoordConst");

	for (i=0;i<nconstraints;i++)
	{
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "central[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		central[i]=ReadThisLine(file,3,1,name, "CoordNumbConst::GetCoordConst");//central atoms index (Starting with 1)
		RunParams::CheckType(central[i],name, "CoordNumbConst::GetCoordConst");
		central[i]-=1;//Resetting to start with 0 in program
		
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
		{
			IntToStr(&conv_numb2,ineightype + 1);
			mystrcpy(name, NAME_SIZE, "neighbours[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			neighbours[cum_n_neightype[i]+ineightype]=ReadThisLine(file,3,1,name, "CoordNumbConst::GetCoordConst");//neighbours  atoms index (Starting with 1)
			RunParams::CheckType(neighbours[cum_n_neightype[i] + ineightype], name, "CoordNumbConst::GetCoordConst");
			neighbours[cum_n_neightype[i]+ineightype]-=1;//Resetting to start with 0 in program
			
		}
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
		{
			IntToStr(&conv_numb2, ineightype + 1);
			mystrcpy(name, NAME_SIZE, "dmin[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			dmin[cum_n_neightype[i] + ineightype] = ReadThisLine(file, 3, 1.0,name, "CoordNumbConst::GetCoordConst");//minimal distance
		}
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
		{
			IntToStr(&conv_numb2, ineightype + 1);
			mystrcpy(name, NAME_SIZE, "dmax[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			dmax[cum_n_neightype[i] + ineightype] = ReadThisLine(file, 3, 1.0, name, "CoordNumbConst::GetCoordConst");//maximum distance
			for (int jj = 0; jj < ineightype; jj++)//check, whether this type was not already given 
			{
				if (neighbours[cum_n_neightype[i] + ineightype] == neighbours[cum_n_neightype[i] + jj])
				{
					//Check, whether the distances for the duplicate secondary type are the same
					if (dmin[cum_n_neightype[i] + ineightype] != dmin[cum_n_neightype[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". coordination constraint's " << jj + 1 << ". and " << ineightype + 1 << ". neighbour types are the same!" << endl;
						cout << "The minimum distances for the " << jj + 1 << ". neighbour is " << dmin[cum_n_neightype[i] + jj] << " A, while for " << ineightype + 1 << ". neighbour is " << dmin[cum_n_neightype[i] + ineightype] << " A." << endl;
						cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different neighbour atom types for a constraint!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
				}
			}

			for (int jj = 0; jj < ineightype; jj++)//check, whether this type was not already given 
			{
				if (neighbours[cum_n_neightype[i] + ineightype] == neighbours[cum_n_neightype[i] + jj])
				{
					//Check, whether the distances for the duplicate secondary type are the same
					if (dmax[cum_n_neightype[i] + ineightype] != dmax[cum_n_neightype[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". neighbour constraint's " << jj + 1 << ". and " << ineightype + 1 << ". neighbour types are the same!" << endl;
						cout << "The maximum distances for the " << jj + 1 << ". neighbour is " << dmax[cum_n_neightype[i] + jj] << " A, while for " << ineightype + 1 << ". neighbour is " << dmax[cum_n_neightype[i] + ineightype] << " A." << endl;
						cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different neighbour atom types for a constraint!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
					//only the second neighbour type is different, the distance parameters are the same, can continue by reducing the number of second neighbour types 
					cout << "\n*****ERROR*****" << endl;
					cout << "The atom type for the " << i + 1 << ". coordination constraint's " << jj + 1 << ". and " << ineightype + 1 << ". neighbour type is the same!" << endl;
					cout << "Correct the " << datfilename << " file using only different neighbour atom types for a constraint!" << endl;
					cout << "Exiting..." << endl;
					CleanExit();
				}
			}
		}

		for (isubconst = 0; isubconst < n_subconst[i]; isubconst++)
		{
			IntToStr(&conv_numb2, isubconst + 1);
			mystrcpy(name, NAME_SIZE, "target coord numb[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			target_coord[cum_n_subconst[i] + isubconst] = ReadThisLine(file, 3, 1, name, "CoordNumbConst::GetCoordConst");//target coordination number
		}
		for (isubconst = 0; isubconst < n_subconst[i]; isubconst++)
		{
			IntToStr(&conv_numb2, isubconst + 1);
			mystrcpy(name, NAME_SIZE, "fraction[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			fraction[cum_n_subconst[i] + isubconst] = ReadThisLine(file, 3, 1.0, name, "CoordNumbConst::GetCoordConst");//fraction of atoms to have target
		}
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
		{
			IntToStr(&conv_numb2, isubconst + 1);
			mystrcpy(name, NAME_SIZE, "sigma[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			if (isubconst<n_subconst[i]-1)
				weights[cum_n_subconst[i]+isubconst]=ReadThisLine(file,3,1.0,name, "CoordNumbConst::GetCoordConst");//weight
			else//this is the last in the row
				weights[cum_n_subconst[i]+isubconst]=ReadThisLine(file,1,1.0,name, "CoordNumbConst::GetCoordConst",datfilename);//weight
			if (weights[cum_n_subconst[i]+isubconst]<0)
				ChiSquared::calc_sigma=1;
		}
	
		//Checking, whether loading was successful
		if (!CheckReadFileState(file,"CoordNumbConst::GetCoordConst",datfilename))
			CleanExit();//loading failed

		//Initialising the remainder of the static members
		if (RunParams::boxedge>0)
		{
			//Minimum and maximum distances for the constraints squared in reduced units
			for (ineightype=0;ineightype<n_neightype[i];ineightype++)
			{
				udminsq[cum_n_neightype[i]+ineightype]=dmin[cum_n_neightype[i]+ineightype]*dmin[cum_n_neightype[i]+ineightype]/RunParams::boxedge/RunParams::boxedge;
				udmaxsq[cum_n_neightype[i]+ineightype]=dmax[cum_n_neightype[i]+ineightype]*dmax[cum_n_neightype[i]+ineightype]/RunParams::boxedge/RunParams::boxedge;
			}
		}

		ncentral[i]=RunParams::pnatoms[central[i]];//number of central atoms

		//Cumulative number of central atoms before this constraint
		ncentral_cum[i]=ncentral_tot;
		//Total number of central atoms
		ncentral_tot+=ncentral[i];

		//Calculating the indices of the central-neighbours partials for each constraint.
		//the cctype will contain the sum of pow(2,partialind) for each partials involved for the given constraint
		//The partialind has to begin with 1 because of this.
		cctype[i]=0;
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
		{

			if (central[i]<=neighbours[cum_n_neightype[i]+ineightype])
				cctype[i]+=(long)pow(2.0,(double)(central[i]*RunParams::ntypes-central[i]*(central[i]+1)/2+neighbours[cum_n_neightype[i]+ineightype]));
			else
				cctype[i]+=(long)pow(2.0,(double)(neighbours[cum_n_neightype[i]+ineightype]*RunParams::ntypes-neighbours[cum_n_neightype[i]+ineightype]*(neighbours[cum_n_neightype[i]+ineightype]+1)/2+central[i]));
		}
		
	}
	//Cumulative number of central atoms before this constraint
	ncentral_cum[nconstraints]=ncentral_tot;
	if (conv_numb != NULL)
		delete [] conv_numb;
	if (conv_numb2 != NULL)
		delete [] conv_numb2;
	delete [] name;
}

//=======loading the CoordNumbConst object from file============
int CoordNumbConst::Load(double *sigma_percentage)
{
	int i,k,dumb,ineightype,isubconst;
	
	int *calc_nsatisfy;//whether the nsatisfy can be loaded (0) or recalculated (1)
	int *pcoordn;
	ifstream file;
	char dumb_char[52];

	double real_dumb;
	
	SafeOpenTextFile(file,cncfilename);
	if (CheckFileState(file,"CoordNumbConst::Load",cncfilename)==0)
	{
		cout << "\nWARNING(" << ++warn << "): There is no open file to load the CoordNumbConst object from"<<endl;
		cout << "\tLoading was not successful, histogram";
		if (RunParams::nbvs>0)
			cout<<", coordination numbers and bond valence sums will be calculated!"<<endl;
		else
			cout << " and coordination numbers will be calculated!" << endl;
		return(0);//loading was not successful
	}
	
	SetArraysize(&calc_nsatisfy,nconstraints,"calc_nsatisfy","CoordNumbConst::Load");

	cout<<"\nLoading the coordination numbers into the CoordNumbConst object"<<endl;
	//Reading the parameters for the constraint, which were already given in the .dat file
	//The parameters will be compared  to check consistency
	SkipLine(file,cncfilename,1);//comment
	SkipLine(file,cncfilename,1);
	dumb=ReadThisLine(file,1,1,"CoordNumbConst nconstraint", "CoordNumbConst::Load",cncfilename);//number of constraints
	if (dumb!=nconstraints)//checking consistency for number of constraints
	{
		cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
		cout<<"\tNumber of constraint is "<<nconstraints<<" from the .dat file,"<<endl;
		cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
		cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
		return(0);//loading was not successful
	}
	
	for (i=0;i<5;i++)
		SkipLine(file,cncfilename,1);//skipping the comments and blank lines

	for(i=0;i<nconstraints;i++)
	{
		dumb=ReadThisLine(file,3,1,"index", "CoordNumbConst::Load");//index
		dumb=ReadThisLine(file,3,1,"central type", "CoordNumbConst::Load");//type of central atom (starts with 1 in file!)
		if (central[i]!=(dumb-1))//checking consistency for type of central atom
		{
			cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
			cout<<"\tThe type of the central atom for the "<<i+1<<". constraint is "<<central[i]+1<<" from the .dat file,"<<endl;
			cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
			return(0);//loading was not successful
		}
		
		dumb=ReadThisLine(file,3,1,"number of neighbour types", "CoordNumbConst::Load");//number of neighbour atom types
		if (n_neightype[i]!=(dumb))//checking consistency fornumber of neighbour atom types
		{
			cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
			cout<<"\tNumber of neighbour atom types for the "<<i+1<<". constraint is "<<n_neightype[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;	
			return(0);//loading was not successful
		};
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
		{
			dumb=ReadThisLine(file,3,1,"neighbour type", "CoordNumbConst::Load");//type of neighbour atom
			if (neighbours[cum_n_neightype[i]+ineightype]!=(dumb-1))//checking consistency for type of neighbour atom
			{
				cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
				cout<<"\t"<<ineightype+1<<". type of neighbour atom for the "<<i+1<<". constraint is "<<neighbours[cum_n_neightype[i]+ineightype]+1<<" from the .dat file,"<<endl;
				cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
				cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;	
				return(0);//loading was not successful
			}
		};

		
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
		{
			real_dumb=ReadThisLine(file,3,1.0,"minimum distance", "CoordNumbConst::Load");//inferior distance
			if (fabs(dmin[cum_n_neightype[i]+ineightype]-real_dumb)>LOAD_TOL)//checking consistency
			{
				cout.precision(12);
				cout.setf(ios::right, ios::scientific);
				cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
				cout<<"\tInferior distance for the "<<i+1<<". constraint's "<<ineightype+1<<". type of neighbour atom is "<<dmin[cum_n_neightype[i]+ineightype]<<" from the .dat file,"<<endl;
				cout<<"\tand "<<real_dumb<<" from the .cnc file!"<<endl;
				cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
				cout.precision(6);
				cout.unsetf(ios::scientific);
				cout.unsetf(ios::right);
				return(0);//loading was not successful
			}
		};

		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
		{
			real_dumb=ReadThisLine(file,3,1.0,"maximum distance", "CoordNumbConst::Load");//superior distance
			if (fabs(dmax[cum_n_neightype[i]+ineightype]-real_dumb)>LOAD_TOL)//checking consistency
			{
				cout.precision(12);
				cout.setf(ios::right, ios::scientific);
				cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
				cout<<"\tSuperior distance for the "<<i+1<<". constraint's "<<ineightype+1<<". type of neighbour atom is "<<dmax[cum_n_neightype[i]+ineightype]<<" from the .dat file,"<<endl;
				cout<<"\tand "<<real_dumb<<" from the .cnc file!"<<endl;
				cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
				cout.precision(6);
				cout.unsetf(ios::scientific);
				cout.unsetf(ios::right);
				return(0);//loading was not successful
			}
		};
				
		dumb=ReadThisLine(file,3,1,"number of subconstraints", "CoordNumbConst::Load");//number of subconstraints
		if (n_subconst[i]!=(dumb))//checking consistency for number of subconstraints
		{
			cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
			cout<<"\tNumber of subconstraints for the "<<i+1<<". constraint is "<<n_subconst[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;	
			return(0);//loading was not successful
		};

		calc_nsatisfy[i]=0;//default, signal that nsatisfy has to be read
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
		{
			dumb=ReadThisLine(file,3,1,"target coordinatiom number", "CoordNumbConst::Load");//Desired coordination number
			if (target_coord[cum_n_subconst[i]+isubconst]!=dumb)//checking consistency for subconstraints
			{
				cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
				cout<<"\tDesired coordination number for the ";
				cout<<"\t"<<isubconst+1<<". subconstraint of "<<i+1<<". constraint is "<<target_coord[cum_n_subconst[i]+isubconst]<<" from the .dat file,"<<endl;
				cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
				cout<<"\tThe value in the .dat file will be used,"<<endl;
				cout<<"\tand the number of atoms satisfying the constraint will be recalculated!"<<endl;
				calc_nsatisfy[i]=1;//signal that nsatisfy is different, no reading is allowed
			}
		};
	
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
		{
			real_dumb=ReadThisLine(file,3,1.0,"fraction of atoms", "CoordNumbConst::Load");//desired fraction of central atoms satisfying the constraint
			if (fraction[cum_n_subconst[i]+isubconst]!=real_dumb)//checking consistency
			{
				cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
				cout<<"\tDesired fraction of central atoms satisfying the ";
				cout<<"\t"<<isubconst+1<<". subconstraint of "<<i+1<<". constraint is "<<fraction[cum_n_subconst[i]+isubconst]<<" from the .dat file,"<<endl;
				cout<<"\tand "<<real_dumb<<" from the .cnc file!"<<endl;
				cout<<"\tThe value in the .dat file will be used!"<<endl;
			};
		}

		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
		{
			real_dumb=ReadThisLine(file,3,1.0,"sigma", "CoordNumbConst::Load");//weight parameter
			int index = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr;
			if ((weights[cum_n_subconst[i]+isubconst]!=real_dumb) && sigma_percentage[index+ cum_n_subconst[i] + isubconst]>0)//checking consistency
			{
				cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
				cout<<"\tWeight parameter for ";
				cout<<"\t"<<isubconst+1<<". subconstraint of "<<i+1<<". constraint is "<<weights[cum_n_subconst[i]+isubconst]<<" from the .dat file,"<<endl;
				cout<<"\tand "<<real_dumb<<" from the .cnc file!"<<endl;
				cout<<"\tThe value in the .dat file will be used!"<<endl;
			};
		}

		dumb=ReadThisLine(file,1,1,"ncentral", "CoordNumbConst::Load",cncfilename);//number of central atoms
		if (ncentral[i]!=dumb)//checking consistency
		{
			cout << "\nWARNING(" << ++warn << "): CoordNumbConst::Load"<<endl;
			cout<<"\tNumber of central atoms for "<<i+1<<". constraint is "<<ncentral[i]<<" from the .dat file,"<<endl;
			cout<<"\tand "<<dumb<<" from the .cnc file!"<<endl;
			cout<<"\tLoading was not successful, histogram and coordination numbers will be calculated!"<<endl;
			return(0);//loading was not successful
		};
	}

	SkipLine(file,cncfilename,1);//skipping the comments and blank lines
	SkipLine(file,cncfilename,1);//skipping the comments and blank lines

	for(i=0;i<nconstraints;i++)//for each constraint
	{
		SkipLine(file,cncfilename,1);
		SkipLine(file,cncfilename,1);//constraint's index

		if (!calc_nsatisfy[i])
		{
			file>>dumb_char;//nsatisfy
			//Reading the number of atoms satisfying each subconstraint 
			for (isubconst=0;isubconst<n_subconst[i];isubconst++)
			{
				if  (isubconst<n_subconst[i]-1)
					nsatisfy[cum_n_subconst[i]+isubconst]=ReadThisLine(file,3,1,"nsatisfy", "CoordNumbConst::Load");
				else//this is the last in the line
					nsatisfy[cum_n_subconst[i]+isubconst]=ReadThisLine(file,1,1,"nsatisfy", "CoordNumbConst::Load",cncfilename);
			}
		}
		else
			SkipLine(file,cncfilename, 1);
		SkipLine(file,cncfilename,1);

		for (k=0;k<ncentral[i];k++)//for every central atom in this constraint
		{
			dumb=ReadThisLine(file,3,1,"serial index", "CoordNumbConst::Load");
			if (i==nconstraints-1 && k==ncentral[i]-1)//last point of last constraint, no error, if no more line in file
				CoordinationNumb(i,k)=ReadThisLine(file,1,1,"coordination number", "CoordNumbConst::Load");
			else
				CoordinationNumb(i,k)=ReadThisLine(file,1,1,"coordination number", "CoordNumbConst::Load",cncfilename);

				
		}

		if (calc_nsatisfy[i])
		{
			//Calculating the number of atoms satisfying the coordination constraint, if necessary
			pcoordn=finder[i];//sets the pointer to the beginning of the given constraint in the coordnumbs array
			for (k=0;k<ncentral[i];k++)
			{
				for (isubconst=0;isubconst<n_subconst[i];isubconst++)
				{
					if (*pcoordn==target_coord[cum_n_subconst[i]+isubconst])
						nsatisfy[cum_n_subconst[i]+isubconst]++;
				}
				pcoordn++;
			}
		}
	}

	delete [] calc_nsatisfy;

	//Checking, whether loading was successful
	if (CheckReadFileState(file,"CoordNumbConst::Load",cncfilename))
		return(1);//load was successful
	else
	{
		cout << "\nWARNING(" << ++warn << "): Loading was not successful, histogram and coordination numbers will be calculated!"<<endl;
		return(0);//loading was not successful
	}
		
	next_line_pos = -1;
	file.close();
	return(1);//load was successful
}
//============save all the data to the file given in argument============
void CoordNumbConst::Save(const char *file_name) 
{ 
	int i,k,ineightype,isubconst;
	ofstream file;

	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, cncfilename);
	OpenFile(file,tempfilename,"CoordNumbConst::Save",0);//open file, check, whether it was successfully opened

		
	file<<"Data for Coordination Numbers Constraints. "\
	<<"file created by CoordNumbConst::Save \n"<<endl;
	
	file<<nconstraints<<"\t number of constraints (nconstraints) "<<endl;
	file<<"Constraints properties: "<<endl;
	file<<"index, central atom type, number of neighbour types, neighbour atom type(s), inferior distance(s),"<<endl;
	file<<"superior distance(s), number of subconstraints, desired coordination number for each subconstraint, "<<endl;
	file<<"desired fraction of central atoms satisfying the subconstraint for each subconstraint, "<<endl;
	file<<"weighting parameter for each subconstraint, number of central atoms, "<<endl;
	
	for(i=0;i<nconstraints;i++)
	{
		file<<J4<<i+1<<"\t";
		file<<J4<<central[i]+1<<"\t";//start with 1 in the file
		file<<J4<<n_neightype[i]<<"\t";
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
			file<<J4<<neighbours[cum_n_neightype[i]+ineightype]+1<<"\t";//start with 1 in the file
		file.precision(12);
		file.setf(ios::fixed, ios::floatfield);
		file.setf(ios::right, ios::adjustfield);
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
			file<<J10<<dmin[cum_n_neightype[i]+ineightype]<<"\t";
		for (ineightype=0;ineightype<n_neightype[i];ineightype++)
			file<<J10<<dmax[cum_n_neightype[i]+ineightype]<<"\t";
		file.precision(6);
		file.unsetf(ios::fixed);
		file.unsetf(ios::right);
		file<<J4<<n_subconst[i]<<"\t";
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
			file<<J4<<target_coord[cum_n_subconst[i]+isubconst]<<"\t";
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
			file<<J10<<fraction[cum_n_subconst[i]+isubconst]<<"\t";
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
			file<<J10<<weights[cum_n_subconst[i]+isubconst]<<"\t";
		file<<J10<<ncentral[i]<<endl;
	}	

	file<<endl;
	file<<"Number of atoms satisfying each subconstraint of the constraint (nsatisfy)"<<endl;
	file<<"Coordination numbers for central atoms per constraint,"<<endl;

	for(i=0;i<nconstraints;i++)//for each constraint
	{
		file<<i+1<<"-th constraint";

		//Number of atoms in the configuration satisfying the constraint
		file<<"\n"<<J10<<"nsatisfy/subconstraint\t";
		for (isubconst=0;isubconst<n_subconst[i];isubconst++)
			file<<J7<<nsatisfy[cum_n_subconst[i]+isubconst]<<"\t";
		file<<"\n"<<endl;

		for (k=0;k<ncentral[i];k++)//for every central atom in this constraint
			file<<J10<<k+1<<"\t"<<J7<<CoordinationNumb(i,k)<<endl;
		
		file<<endl;
	}
	file.close();
}
//----Copy the coordination numbers and the number of atoms satisfying the ---------
//----the constraints affected by the move from source to target---------
void CoordNumbConst::CopyModified(CoordNumbConst &target, ThreadArg &thread_arg)
{
	int iconst,i,ineightype,isubconst;
	int partialind;
	int min_ind, n_cent;
	int *p_source,*p_target;

	for (iconst=0;iconst<nconstraints;iconst++)
	{
		partialind=0;
		//searching for the first existing partial of this constraint
		for (ineightype=0;ineightype<n_neightype[iconst];ineightype++)
		{
			while (((cctype[iconst] >> partialind) & 1)==0)
				partialind++;
							
			if (Threads::mod_partial[partialind])
			{
				//cc_central_offset_min, cc_central_offset_max contains the segmentation for the central atoms of each constraint
				min_ind=thread_arg.cc_central_offset_min[iconst];//offset of the first central 
				n_cent=thread_arg.cc_central_offset_max[iconst]-min_ind+1;//number of centrals
				
				//Coordination numbers
				p_source=finder[iconst]+min_ind;//sets the pointer to the source's coordnumbs array 
												//to the beginning of the iconst-th constraint for this tread
				p_target=target.finder[iconst]+min_ind;//sets the pointer to the target's coordnumbs array 
											//to the beginning of the iconst-th constraint for this tread
				for (i=0;i<n_cent;i++)
					*p_target++=*p_source++;//copy source coordination numbers to target

				if (thread_arg.thread_index==0)//only the 0-th thread will copy nsatisfy
				{
					//Copy the number of atoms satisfying each subconstraint of the constraint 
					for (isubconst=0;isubconst<n_subconst[iconst];isubconst++)
						target.nsatisfy[cum_n_subconst[iconst]+isubconst]=nsatisfy[cum_n_subconst[iconst]+isubconst];
				}
				break;//at least this partial was modified, the array was copied, no need to search for more modified partial
			}//end of if this constraint was affected by the move
			partialind++;//this partial of the constraint was not modified, goto next
		}//end of neighbour type cycle ineightype 
	}//end of constraint cycle iconst
};

//Updating the coordnumbs array with the contribution of the threads, if there is any, and recalculating the number of atoms satisfying the constraint
int CoordNumbConst::UpdateCoordnumb(ThreadArg &thread_arg)
{
	int iconst,isubconst,ithread,j;
	int neightype;
	int	modified=0;
	int partialind;
	int min_ind, n_cent;
	int *pcoordn,*pcoord_aux;
	longint extract_part;
#ifdef _TEST_MODE
	int index=thread_arg.thread_index;
	std::chrono::duration<double, std::milli> elapsed;
#endif

	//sets the pointer to the beginning of the coordnumbs array
	for (iconst=0;iconst<CoordNumbConst::nconstraints;iconst++)
	{

		extract_part=CoordNumbConst::cctype[iconst];
		partialind=0;
		//searching for the first existing partial of this constraint
		for (neightype=0;neightype<CoordNumbConst::n_neightype[iconst];neightype++)
		{

			while (((extract_part >> partialind) & 1)==0)
				partialind++;
				
			if (Threads::mod_partial[partialind])
			{
#ifdef _TEST_MODE
				Threads::ctime1[index] = std::chrono::high_resolution_clock::now();
#endif
				modified=1;//at least one constraint was modified by the move
				//the constraint was effected by the move
				//first add the contribution of the thread's coordnumbs array to the main
				
				//cc_central_offset_min, cc_central_offset_max contains the segmentation for the central atoms of each constraint
				min_ind=thread_arg.cc_central_offset_min[iconst];//offset of the first central 
				n_cent=thread_arg.cc_central_offset_max[iconst]-min_ind+1;//number of centrals

				for (ithread=0;ithread<nthreads-1;ithread++)
				{
					pcoordn=finder[iconst]+min_ind;//sets the pointer to the first atom of this constraint in the coordnumbs array fot this thread
					pcoord_aux=thread_finder[ithread]+CoordNumbConst::ncentral_cum[iconst]+min_ind;
					for (j=0;j<n_cent;j++)
					{
						*pcoordn+=*pcoord_aux;
						*pcoord_aux++=0;
						pcoordn++;
					}
				}


				//now recalculate the number of atoms satisfying the constraint for the atoms of this tread for all the constraints
				pcoordn=finder[iconst]+min_ind;
				for (isubconst=0;isubconst<n_subconst[iconst];isubconst++)
					*(thread_arg.nsatisfy_finder[iconst]+isubconst)=0;

#ifdef _TEST_MODE
				Threads::ctime2[index] = std::chrono::high_resolution_clock::now();
				elapsed = Threads::ctime2[index]-Threads::ctime1[index];
				Threads::dur_16[index]+=elapsed.count();
#endif

				for (j=0;j<n_cent;j++)
				{
					for (isubconst=0;isubconst<n_subconst[iconst];isubconst++)
					{
						if (*pcoordn==CoordNumbConst::target_coord[cum_n_subconst[iconst]+isubconst])
							*(thread_arg.nsatisfy_finder[iconst]+isubconst)+=1;
			
					}
					pcoordn++;
				}
#ifdef _TEST_MODE
				Threads::ctime1[index] = std::chrono::high_resolution_clock::now();
				elapsed = Threads::ctime1[index]-Threads::ctime2[index];
				Threads::dur_17[index]+=elapsed.count();
#endif
				break;//it is recalculated, no need to search for more modified partials
			}
			partialind++;//this partial of the constraint was not modified, goto next
		}
	}
	return(modified);
};
//----------------Getting the parameters from the.dat file, and sets the static members-------------------
bool CoordNumbConst::GetCoordConstFree(std::list<std::list<string> > pool )
//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
{
	//int i,ineightype,isubconst;
	unsigned int npartials = RunParams::ntypes*(RunParams::ntypes + 1) / 2;
	nconstraints = RunParams::nicoord;
	char *name;
	if (npartials > sizeof(*cctype) * 8)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"The program can handle only " << sizeof(*cctype) * 8 << " partials presently!" << endl;
		cout << "Change the type of CoordNumbConst::cctype if possible to have more bytes!" << endl;
		return false;
	}
	SetArraysize(&name, NAME_SIZE, "name", "CoordNumbConst::GetCoordConstFree");
	//allocating memory
	SetArraysize(&n_neightype, nconstraints, "n_neightype", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&n_subconst, nconstraints, "n_subconst", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&central, nconstraints, "central", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&ncentral, nconstraints, "ncentral", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&cctype, nconstraints, "cctype", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&write_detail, nconstraints, "write_detail", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&cum_n_neightype, nconstraints + 1, "cum_n_neightype", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&cum_n_subconst, nconstraints + 1, "cum_n_subconst", "CoordNumbConst::GetCoordConstFree");
	SetArraysize(&ncentral_cum, nconstraints + 1, "ncentral_cum", "CoordNumbConst::GetCoordConstFree");

	for (int i = 0; i < nconstraints; i++)
		write_detail[i] = WRITE_CNC_DETAIL_DEF;
	maxdistsq = 0;
	cum_n_subconst[0] = 0;//cumulative number of subconstraints before the given constraint
	cum_n_neightype[0] = 0;//cumulative number of neighbour type before the given constraint
	tot_neightype = 0;
	tot_subconst = 0;
	ncentral_tot = 0;

	int i = 0;
	bool retval = true;
	for (std::list<std::list<string>>::iterator it2 = pool.begin(); it2 != pool.end() && retval; it2++, i++)
	{
		//Initialization and allocating memories
		n_neightype[i] =(int)std::count_if(it2->begin(), it2->end(), [](string c) { return c.find(CheckKey("NEIGH-TYPE_FROM_TO")) != string::npos; });
		n_subconst[i] = (int)std::count_if(it2->begin(), it2->end(), [](string c) { return (c.find(CheckKey("SIGMA")) != string::npos); });
		cum_n_neightype[i + 1] = cum_n_neightype[i] + n_neightype[i];//cumulative number of neighbour type before the given constraint
		cum_n_subconst[i + 1] = cum_n_subconst[i] + n_subconst[i];//cumulative number of subconstraint before the given constraint

		//allocating memory-total number of neighbour types
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &neighbours, "neighbours", "CoordNumbConst::GetCoordConstFree"); tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &dmin, "dmin", "CoordNumbConst::GetCoordConstFree"); tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &dmax, "dmax", "CoordNumbConst::GetCoordConstFree"); tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &udminsq, "udminsq", "CoordNumbConst::GetCoordConstFree"); tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &udmaxsq, "udmaxsq", "CoordNumbConst::GetCoordConstFree");
		//total number of subconstraints
		ResizeArray(&tot_subconst, tot_subconst + n_subconst[i], &target_coord, "target_coord", "CoordNumbConst::GetCoordConstFree"); tot_subconst -= n_subconst[i];
		ResizeArray(&tot_subconst, tot_subconst + n_subconst[i], &fraction, "fraction", "CoordNumbConst::GetCoordConstFree"); tot_subconst -= n_subconst[i];
		ResizeArray(&tot_subconst, tot_subconst + n_subconst[i], &weights, "weights", "CoordNumbConst::GetCoordConstFree");

		int ce = 0, n = 0, co = 0;
		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{

			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("CENT-TYPE"))
			{
				if (ce++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple "<< CheckKey("CENT-TYPE")<<" definition in the " << i + 1 << ". "<< CheckTag("COORD")<<" constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (!(ss >> central[i])) return ErrLackValue(*it);
				if ((central[i] > RunParams::ntypes) || (central[i] < 1))
				{
					cout << "\nERROR: The "<< CheckKey("CENT-TYPE")<<" value"<<" (" << central[i] << ") should be in the range: 1.." << RunParams::ntypes << " (number of types)!" << endl << "Please, correct it at the " << i + 1 << ". "<<CheckTag("COORD")<<" constraint!" << endl;
					return false;
				}
				central[i]--;
				ncentral[i] = SimpleCfg::pnatoms[central[i]];//it will be reset in case of no_periodic
				//Cumulative number of central atoms before this constraint
				ncentral_cum[i] = ncentral_tot;
				//Total number of central atoms
				ncentral_tot += ncentral[i];
			}
			else if (key == CheckKey("NEIGH-TYPE_FROM_TO"))
			{
				int aind = cum_n_neightype[i] + n++;//Actual index
				if (!GetInt(ss, *it, neighbours+aind,"CoordNumbConst::GetCoordConstFree")) return ErrLackValue(*it);
				if (!(ss>> dmin[aind] >> dmax[aind])) return ErrLackValue(*it);
				if ((neighbours[aind] > RunParams::ntypes) || (neighbours[aind] < 1))
				{
					cout << "\nERROR: The NEIGH-TYPE(_FROM_TO) value (" << neighbours[aind] << ") should be in the range: 1.." << RunParams::ntypes << " (number of types)!" << endl << "Please, correct it at the " << i + 1 << ". "<<CheckTag("COORD")<<" constraint!" << endl;
					return false;
				}
				neighbours[aind]--;
				if (dmin[aind] > dmax[aind])
				{
					cout << "\nWARNING(" << ++warn << "): the second value should be less than the third one in the following line" << endl << *it << endl;
					double temp = dmin[aind];
					dmin[aind] = dmax[aind];
					dmax[aind] = temp;
				}
				for (int jj = 0; jj < n-1; jj++)//check, whether this type was not already given
				{
					if (neighbours[aind] == neighbours[cum_n_neightype[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmin[aind] != dmin[cum_n_neightype[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". coordination constraint's " << jj + 1 << ". and " << aind - cum_n_neightype[i] + 1 << ". neighbour types are the same!" << endl;
							cout << "The number of neighbour types cannot be reduced, as different minimum distances were given!" << endl;
							cout << "The minimum distances for the " << jj + 1 << ". neighbour is " << dmin[cum_n_neightype[i] + jj] << " A, while for " << aind - cum_n_neightype[i] + 1 << ". neighbour is " << dmin[aind] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different neighbour atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
					}
				}

				for (int jj = 0; jj < n-1; jj++)//check, whether this type was not already given
				{
					if (neighbours[aind] == neighbours[cum_n_neightype[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmax[aind] != dmax[cum_n_neightype[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". neighbour constraint's " << jj + 1 << ". and " << aind - cum_n_neightype[i] + 1 << ". neighbour types are the same!" << endl;
							cout << "The number of neighbour types cannot be reduced, as different maximum neighbour distances were given!" << endl;
							cout << "The maximum distances for the " << jj + 1 << ". neighbour is " << dmax[cum_n_neightype[i] + jj] << " A, while for " << aind - cum_n_neightype[i] + 1 << ". neighbour is " << dmax[aind] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different neighbour atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
						//only the second neighbour type is different, the distance parameters are the same, can continue by reducing the number of second neighbour types 
						cout << "\nWARNING(" << ++warn << "): The atom type for the " << i + 1 << ". coordination constraint's " << jj + 1 << ". and " << aind - cum_n_neightype[i] + 1 << ". neighbour type is the same!" << endl;
						cout << "\tThe number of neighbor types will be reduced!" << endl;
						n_neightype[i]--;
					}
				}
				//Initialising the remainder of the static members
				if (RunParams::boxedge > 0)
				{
					//Minimum and maximum distances for the constraints squared in reduced units
					udminsq[aind] = pow(dmin[aind], 2) / RunParams::boxedge / RunParams::boxedge;
					udmaxsq[aind] = pow(dmax[aind], 2) / RunParams::boxedge / RunParams::boxedge;
					if (udmaxsq[aind] > maxdistsq)
						maxdistsq = udmaxsq[aind];
				}
			}
			else if ((key == CheckKey("COORDNUM_FRACT_SIGMA")) || (key == CheckKey("COORDNUM_FRACT_SIGMA-MASTER")) || (key == CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE")))
			{
				int aind = cum_n_subconst[i] + co++;//Actual index
				if (!GetInt(ss, *it, target_coord+aind,"CoordNumbConst::GetCoordConstFree")) return ErrLackValue(*it);
				if (!(ss >>  fraction[aind])) return ErrLackValue(*it);
				if (!(ss >> weights[aind])) return ErrLackValue(*it);
	
				if (key == CheckKey("COORDNUM_FRACT_SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1)
						RunParams::lead_series_ind = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg+ RunParams::nek + RunParams::ncosdistr + aind + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): Multiple **(..)SIGMA-MASTER** entry is declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
#ifdef _ADVANCED_GEOM_CONST
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#endif
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout<< " section!" << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE"))
				{
					weights[aind] = -fabs(weights[aind]);
					ChiSquared::calc_sigma = 1;
				}
				else //normal, should not be negative, it does not effect calculation, but would upset *.free...
				{
					weights[i] = fabs(weights[i]);
				}
			}
			else if (key == CheckKey("WRITE-CNC-DETAIL"))
			{
				if (!(ss >> write_detail[i])) return ErrLackValue(*it);
				mystrcpy(name, NAME_SIZE, "WRITE-CNC-DETAIL");
				Check0_1(write_detail[i], name, "CoordNumbConst::GetCoordConstFree");
				if (write_detail[i])
					is_detail = 1;
			}
			else
			{
				CheckKeyTag(key, CheckTag("COORD"));//write the appropriate error mesage
				return false;
			}
		}
		if (retval)
		{
			//Calculating the indices of the central-neighbours partials for each constraint.
			//the cctype will contain the sum of pow(2,partialind) for each partials involved for the given constraint
			//The partialind has to begin with 1 because of this.
			cctype[i] = 0;
			for (int ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			{
				if (central[i] <= neighbours[cum_n_neightype[i] + ineightype])
					cctype[i] += (long)pow(2.0, (double)(central[i] * RunParams::ntypes - central[i] * (central[i] + 1) / 2 + neighbours[cum_n_neightype[i] + ineightype]));
				else
					cctype[i] += (long)pow(2.0, (double)(neighbours[cum_n_neightype[i] + ineightype] * RunParams::ntypes - neighbours[cum_n_neightype[i] + ineightype] * (neighbours[cum_n_neightype[i] + ineightype] + 1) / 2 + central[i]));
			}
			
		}
	}
	
	if (name != NULL)
		delete[] name;
	//Cumulative number of central atoms before this constraint
	if (retval) ncentral_cum[nconstraints] = ncentral_tot;
	return retval;
};

void CoordNumbConst::ResizeCoordArrays( int **ind_list, int **type_list, double **dist_list)
{
	ResizeArray(&max_neigh, max_neigh + 5, &neigh_list, "neigh_list", "CoordNumbConst::ResizeCoordArrays");
	max_neigh -= 5;
	ResizeArray(&max_neigh, max_neigh + 5, &neigh_type, "neigh_type", "CoordNumbConst::ResizeCoordArrays");
	max_neigh -= 5;
	ResizeArray(&max_neigh, max_neigh + 5, &dist, "dist", "CoordNumbConst::ResizeCoordArrays");
	*ind_list = neigh_list;
	*type_list = neigh_type;
	*dist_list = dist;

};

//caclculate and save the detials of the neighbours in the *.cncd file
void CoordNumbConst::SaveDetail(const char *file_name)
{
	int i, j,iconst;
	int nneigh;
	ofstream file;

	if (strlen(file_name) != 0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, cncdfilename);
	OpenFile(file, tempfilename, "CoordNumbConst::SaveDetail", 0);//open file, check, whether it was successfully opened
	file.precision(6);
	file.setf(ios::adjustfield, ios::right);
	file.setf(ios::floatfield, ios::fixed);

	file << "Data for Coordination Numbers Constraints details "\
		<< "file created by CoordNumbConst::SaveDetail \n" << endl;

	file << nconstraints << "\t number of constraints (nconstraints) " << endl;
	file << "For each central atom there are four lines:" << endl;
	file << "1. line: central index, number of neighbours" << endl;
	file << "2. line: neighbour indices" << endl;
	file << "3. line: neighbour types" << endl;
	file << "4. line: central-neighbour distances in A" << endl;
	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		if (write_detail[iconst])
		{
			file <<"\n"<< iconst + 1<< ". constraint" << endl;
			for (i = SimpleCfg::cumul[central[iconst]]; i < SimpleCfg::cumul[central[iconst] + 1]; i++)//going through the central atoms of this constraint
			{
				nneigh = neighlist.CalcNeighList(i, neigh_list, neigh_type, nullptr, nullptr, dist, &max_neigh, iconst);//calculating
				file << J10 << i + 1 << "\t" << J4 << nneigh << endl;
				file << "\t" << J10 << " ";
				for (j = 0; j < nneigh; j++)
					file << "\t" << J10 << neigh_list[j] + 1;
				file << endl;
				file << "\t" << J10 << " ";
				for (j = 0; j < nneigh; j++)
					file << "\t" << J10 << neigh_type[j] + 1;
				file << endl;
				file << "\t" << J10 << " ";
				for (j = 0; j < nneigh; j++)
					file << "\t" << J10 << dist[j];
				file << endl;
			}
		}
	}
	file.unsetf(ios::fixed);
	file.unsetf(ios::right);
	file.close();
};