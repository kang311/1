//source CommonNeighConst.cpp
//Last changed 12.12.2022



#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Move.h"//classes1.h is included through this

#ifdef _ADVANCED_GEOM_CONST
int  CommonNeighConst::nconstraints;//the number of constraints
int  CommonNeighConst::tot_sec;//total number of secondary neighbours for all the constraints
int  CommonNeighConst::tot_primary_pairs;//total number possible primary pairs
int  CommonNeighConst::max_hist_size=50;//number of histogram bins(number of common neighbours+1)
int *CommonNeighConst::nsec_neigh=NULL;//number of secondary neighbours for a costraint
int *CommonNeighConst::cumul_sec;//cumulative number of secondary neighbours for each constraint
int *CommonNeighConst::cumul_prim;//cumulative number of primary pairs for each constraint
int *CommonNeighConst::max_neigh;//pointer to the maximum number of neighbours for the neighbourlist for an atom
longint *CommonNeighConst::consttype;//pointer to NeighbourList::consttype, each bit represents a type, containing a 1 for each type that is involved in a constraint
longint CommonNeighConst::comconsttype=0;//each bit represents a type, containing a 1 for each type that is involved in a CONC constraint

//parameters for each constraint

int *CommonNeighConst::primary1;//array of types of primary1 particles [nconstraints]
int *CommonNeighConst::primary2;//array of types of primary2 particles [nconstraints]
int *CommonNeighConst::secondary;//array of types of secondary neighbour particles  [tot_sec]
int *CommonNeighConst::target_coord;//array of DESIRED coord. numbers [nconstraints]
double *CommonNeighConst::dmin;//array of minimal distances between primary1 and primary2 [nconstraints]
double *CommonNeighConst::dmax;//array of maximum distances between primary1 and primary2 [nconstraints]
double *CommonNeighConst::dmin1;//array of minimal distances between primary1 and secondary [tot_sec]
double *CommonNeighConst::dmin2;//array of minimal distances between primary2 and secondary [tot_sec]
double *CommonNeighConst::dmax1;//array of maximum distances between primary1 and secondary [tot_sec]
double *CommonNeighConst::dmax2;//array of maximum distances between primary2 and secondary [tot_sec]
double *CommonNeighConst::weights;//array of weights for each constraint [nconstraints]
double *CommonNeighConst::fraction;//desired fraction for each constraint [nconstraints]

int *CommonNeighConst::mod_const;//whteher the constraint was effected by the move
int *CommonNeighConst::primary_partial;//index of the partial calculated from the two primary types
double *CommonNeighConst::minsq;//reduced min. dist. squared
double *CommonNeighConst::min1sq;//reduced min. dist. squared
double *CommonNeighConst::min2sq;//reduced min. dist. squared
double *CommonNeighConst::maxsq;//reduced min. dist. squared
double *CommonNeighConst::max1sq;//reduced max. dist. squared
double *CommonNeighConst::max2sq;//reduced max. dist. squared
double CommonNeighConst::maxdistsq;//the maximum distance among the maxsq, max1sq and max2sq arrays, this is used determining the limit of the neighbourlist calculation

//----default constructor----------------
CommonNeighConst::CommonNeighConst(NeighbourList &neigh_list, SimpleCfg &conf)
	:neighlist(neigh_list), config(conf)
{
	int i;
	if (nconstraints < 1)
	{
		if (::debug)
		{

			cout << "\nWARNING(" << ++warn << "): CommonNeighConst constructor" << endl;
			cout << "\tThere are no common neighbour constraints, or nconstraints was not initialised!" << endl;
		}
	}
	
	SetArraysize(&nsatisfy, nconstraints, "nsatisfy", "CommonNeighConst::CommonNeighConst");
	SetArraysize(&nprimary, nconstraints, "nprimary", "CommonNeighConst::CommonNeighConst");
	
	//as the number of primary pairs can be large, and we do not need them individually, the histogram will be used instead of coord numbs
	SetArraysize(&hist, max_hist_size*nconstraints, "hist", "CommonNeighConst::CommonNeighConst");
	for (i = 0; i < nconstraints*max_hist_size; i++)
		hist[i] = 0;
};

//copy constructor
CommonNeighConst::CommonNeighConst(CommonNeighConst &source, NeighbourList &neigh_list, SimpleCfg &conf)
	:neighlist(neigh_list), config(conf)
{
	int i;
	
	SetArraysize(&nsatisfy, nconstraints, "nsatisfy", "CommonNeighConst::CommonNeighConst");
	SetArraysize(&nprimary, nconstraints, "nprimary", "CommonNeighConst::CommonNeighConst");
	//as the number of primary pairs can be large, and we do not need them individually, the histogram will be used instead of coord numbs
	SetArraysize(&hist, max_hist_size*nconstraints, "hist", "CommonNeighConst::CommonNeighConst");
	for (i = 0; i < nconstraints*max_hist_size; i++)
		hist[i] = source.hist[i];
	for (i = 0; i < nconstraints; i++)
	{
		nsatisfy[i] = source.nsatisfy[i];
		nprimary[i] = source.nprimary[i];
	}
};


//--------------gets the static parameters------------------------------------
void CommonNeighConst::GetCommonNeighConst(ifstream &file)
{
	int i,j,jj,index;
	char *name, *conv_numb = NULL;
	char *conv_numb2 = NULL;
	if (CheckFileState(file, "CommonNeighConst::GetCommonNeighConst", datfilename) == 0)
	{
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}
	nconstraints = RunParams::ncommonneigh;
	consttype = &NeighbourList::consttype;
	
	SetArraysize(&nsec_neigh, nconstraints, "nsec_neigh", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&cumul_prim, nconstraints+1, "cumul_prim", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&cumul_sec, nconstraints + 1, "cumul_sec", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&name,NAME_SIZE,"name", "CommonNeighConst::GetCommonNeighConst");

	tot_sec = 0;
	cumul_sec[0] = 0;
	for (i = 0; i < nconstraints; i++)
	{
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "nsec_neigh");
		mystrcat(name, NAME_SIZE, conv_numb);
		nsec_neigh[i] = ReadThisLine(file, 8, 1, name, "CommonNeighConst::GetCommonNeighConst");
		tot_sec += nsec_neigh[i];
		cumul_sec[i + 1] = tot_sec;
	}
	SkipLine(file);

	//creating the static arrays
	SetArraysize(&primary_partial, nconstraints, "primary_partial", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&mod_const, nconstraints, "mod_const", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&primary1, nconstraints, "primary1", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&primary2, nconstraints, "primary2", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&dmin, nconstraints, "dmin", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&dmax, nconstraints, "dmax", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&minsq, nconstraints, "minsq", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&maxsq, nconstraints, "maxsq", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&fraction, nconstraints, "fraction", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&target_coord, nconstraints, "target_coord", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&weights, nconstraints, "weights", "CommonNeighConst::GetCommonNeighConst");
	
	SetArraysize(&secondary, tot_sec, "secondary", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&dmin1, tot_sec, "dmin1", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&dmin2, tot_sec, "dmin2", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&dmax1, tot_sec, "dmax1", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&dmax2, tot_sec, "dmax2", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&min1sq, tot_sec, "min1sq", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&min2sq, tot_sec, "min2sq", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&max1sq, tot_sec, "max1sq", "CommonNeighConst::GetCommonNeighConst");
	SetArraysize(&max2sq, tot_sec, "max2sq", "CommonNeighConst::GetCommonNeighConst");
	
	maxdistsq = 0.0;
	tot_primary_pairs = 0;
	cumul_prim[0] = 0;
	for (i = 0; i < nconstraints; i++)
	{
		IntToStr(&conv_numb, i + 1);
		mystrcpy(name, NAME_SIZE, "primary1[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		primary1[i] = ReadThisLine(file, 3, 1,name, "CommonNeighConst::GetCommonNeighConst");
		RunParams::CheckType(primary1[i], name, "CommonNeighConst::GetCommonNeighConst");
		primary1[i]--;
		mystrcpy(name, NAME_SIZE, "primary2[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		primary2[i] = ReadThisLine(file, 3, 1, name, "CommonNeighConst::GetCommonNeighConst");
		RunParams::CheckType(primary2[i], name, "CommonNeighConst::GetCommonNeighConst");
		primary2[i]--;
		//The concept is that the neighbourlist is calculated for atoms of each type represented in a cosine distribution or common neighbour
		//constraint. It is required for the neighbours as well, as the list is updated, and not completly recalculated at each step.
		//As each type can be involved in more, than one constraint, all the neighbourlist will be calculated with the
		//same maximum distance, maxdistsq. 
		//Be aware that the central and the neighbour type in comconsttype types represented by right to left!
		//type 0 is represented by binary 1, type 1 binary 10 and so on  
		//consttype is a pointer to NeighbourList::consttype for all the constraints using the neighbour list
		//comconsttype only for cos constraints
		if (!(*consttype >> primary1[i] & 1))//this type is not yet represented
		{
			*consttype += (longint)pow(2.0, (double)primary1[i]);
			NeighbourList::nconsttype++;//increase the number of constrained types
		}
		if (!(*consttype >> primary2[i] & 1))//this type is not yet represented
		{
			*consttype += (longint)pow(2.0, (double)primary2[i]);
			NeighbourList::nconsttype++;//increase the number of constrained types
		}
		if (!(comconsttype >> primary1[i] & 1))//this type is not yet represented
			comconsttype += (longint)pow(2.0, (double)primary1[i]);
		
		if (!(comconsttype >> primary2[i] & 1))//this type is not yet represented
			comconsttype += (longint)pow(2.0, (double)primary2[i]);
		
		mystrcpy(name, NAME_SIZE, "minimum primary distance, dmin[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmin[i] = ReadThisLine(file, 3, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
		mystrcpy(name, NAME_SIZE, "maximum primary distance, dmax[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		dmax[i] = ReadThisLine(file, 3, 1.0,name, "CommonNeighConst::GetCommonNeighConst");
		//calculating the reduced squared distances
		minsq[i] = pow(dmin[i] / RunParams::boxedge, 2.);
		maxsq[i] = pow(dmax[i] / RunParams::boxedge, 2.);
		if (maxsq[i] > maxdistsq)
			maxdistsq = maxsq[i];
		index = cumul_sec[i];
		for (j = 0; j < nsec_neigh[i]; j++)
		{
			IntToStr(&conv_numb2, j + 1);
			mystrcpy(name, NAME_SIZE, "secondary[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			secondary[index] = ReadThisLine(file, 3, 1, name, "CommonNeighConst::GetCommonNeighConst");
			RunParams::CheckType(secondary[index], name, "CommonNeighConst::GetCommonNeighConst");
			secondary[index]--;
			if (!(*consttype >> secondary[index] & 1))//this type is not yet represented
			{
				*consttype += (longint)pow(2.0, (double)secondary[index]);
				NeighbourList::nconsttype++;//increase the number of constrained types
			}
			if (!(comconsttype >> secondary[index] & 1))//this type is not yet represented
				comconsttype += (longint)pow(2.0, (double)secondary[index]);
			
			index++;
		}
		index = cumul_sec[i];
		for (j = 0; j < nsec_neigh[i]; j++)
		{
			IntToStr(&conv_numb2, j + 1);
			mystrcpy(name, NAME_SIZE, "minimum primary1 - secondary");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, " distance, dmin1[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			dmin1[index] = ReadThisLine(file, 3, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
			for (jj = 0; jj < j; jj++)//check, whether this type was not already given as secondary
			{
				if (secondary[index] == secondary[cumul_sec[i] + jj])
				{
					//Check, whether the distances for the duplicate secondary type are the same
					if (dmin1[index] != dmin1[cumul_sec[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
						cout << "The number of secondary types cannot be reduced, as different minimum primary1-secondary distances were given!" << endl;
						cout << "The distances for the " << jj + 1 << ". secondary and primary1 is " << dmin1[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmin1[index] << " A." << endl;
						cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
				}
			}
			min1sq[index] = pow(dmin1[index] / RunParams::boxedge, 2.);//calculating the reduced squared distances
			index++;
		}
		index = cumul_sec[i];
		for (j = 0; j < nsec_neigh[i]; j++)
		{
			IntToStr(&conv_numb2, j + 1);
			mystrcpy(name, NAME_SIZE, "minimum primary2 - secondary");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, " distance, dmin2[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			dmin2[index] = ReadThisLine(file, 3, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
			for (jj = 0; jj < j; jj++)//check, whether this type was not already given as secondary
			{
				if (secondary[index] == secondary[cumul_sec[i] + jj])
				{
					//Check, whether the distances for the duplicate secondary type are the same
					if (dmin2[index] != dmin2[cumul_sec[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
						cout << "The number of secondary types cannot be reduced, as different minimum primary2-secondary distances were given!" << endl;
						cout << "The distances for the " << jj + 1 << ". secondary and primary2 is " << dmin2[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmin2[index] << " A." << endl;
						cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
				}
			}
			//check, whether in case of the same atom type for both primaries, the dmin1 and dmin2 are the same
			if (primary1[i] == primary2[i] && dmin1[index] != dmin2[index])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The atom types of the two primary neighbours for the "<< i + 1 << ". common neighbour constraint's "<<" are the same, but the secondary-primary1 and secondary-primary2 minimum distance is different." << endl;
				cout << "The secondary-primary1 minimum distance is "<< dmin1[index] << "A, and secondary-primary2 is "<<dmin2[index] << " A!"<<endl;
				cout << "This does not makes sense, give the same value for both, correct the " << datfilename << " file! Exiting..." << endl;
				CleanExit();
			}
			min2sq[index] = pow(dmin2[index] / RunParams::boxedge, 2.);
			index++;
		}
		index = cumul_sec[i];
		for (j = 0; j < nsec_neigh[i]; j++)
		{
			IntToStr(&conv_numb2, j + 1);
			mystrcpy(name, NAME_SIZE, "maximum primary1 - secondary");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, " distance, dmax1[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			dmax1[index] = ReadThisLine(file, 3, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
			for (jj = 0; jj < j; jj++)//check, whether this type was not already given as secondary
			{
				if (secondary[index] == secondary[cumul_sec[i] + jj])
				{
					//Check, whether the distances for the duplicate secondary type are the same
					if (dmax1[index] != dmax1[cumul_sec[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index -cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
						cout << "The number of secondary types cannot be reduced, as different maximum primary1-secondary distances were given!" << endl;
						cout << "The distances for the " << jj + 1 << ". secondary and primary1 is " << dmax1[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmax1[index] << " A." << endl;
						cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
				}
			}
			max1sq[index] = pow(dmax1[index] / RunParams::boxedge, 2.);
			if (max1sq[index] > maxdistsq)//finding the largest reduced maximum distance square
				maxdistsq = max1sq[index];
			index++;
		}
		index = cumul_sec[i];
		for (j = 0; j < nsec_neigh[i]; j++)
		{
			IntToStr(&conv_numb2, j + 1);
			mystrcpy(name, NAME_SIZE, "maximum primary2 - secondary");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, " distance, dmax2[");
			mystrcat(name, NAME_SIZE, conv_numb);
			mystrcat(name, NAME_SIZE, ",");
			mystrcat(name, NAME_SIZE, conv_numb2);
			mystrcat(name, NAME_SIZE, "]");
			dmax2[index] = ReadThisLine(file, 3, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
			//now have to check first, whether the data makes any sense
			for (jj = 0; jj < j; jj++)//check, whether this type was not already given as secondary
			{
				if (secondary[index] == secondary[cumul_sec[i] + jj])
				{
					//Check, whether the distances for the duplicate secondary type are the same
					if (dmax2[index] != dmax2[cumul_sec[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
						cout << "The number of secondary types cannot be reduced, as different minimum primary2-secondary distances were given!" << endl;
						cout << "The distances for the " << jj + 1 << ". secondary and primary2 is " << dmax2[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmax2[index] << " A." << endl;
						cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!"<<endl;
						cout<<"Exiting..." << endl;
						CleanExit();
					}
					//only the secondary type is different, the distance parameters are the same, can continue by reducing the number of secondary types 
					cout << "\nWARNING(" << ++warn << "): The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary type is the same!" << endl;
					cout << "\tThe number of secondary types will be reduced!" << endl;
					nsec_neigh[i]--;
				}
			}
			//check, whether in case of the same atom type for both primaries, the dmax1 and dmax2 are the same
			if (primary1[i] == primary2[i] && dmax1[index] != dmax2[index])
			{
				cout << "\nWARNING(" << ++warn << "): The atom types of the two primary neighbours for the " << i + 1 << ". common neighbour constraint's " <<" are the same, but the secondary-primary1 and secondary-primary2 maximum distance is different." << endl;
				cout << "\tThe secondary-primary1 maximum distance is " << dmax1[index] << "A, and secondary-primary2 is " << dmax2[index] << " A!" << endl;
				cout << "\tThis does not makes sense, give the same value for both, correct the " << datfilename << " file! Exiting..." << endl;
				CleanExit();
			}
			max2sq[index] = pow(dmax2[index] / RunParams::boxedge, 2.);
			if (max2sq[index] > maxdistsq)//finding the largest reduced maximum distance square
				maxdistsq = max2sq[index];
			index++;
		}
					
		mystrcpy(name, NAME_SIZE, "desired coordination number, target_coord[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		target_coord[i] = ReadThisLine(file, 3, 1, name, "CommonNeighConst::GetCommonNeighConst");
		mystrcpy(name, NAME_SIZE, "fraction[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		fraction[i] = ReadThisLine(file, 3, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
		mystrcpy(name, NAME_SIZE, "weights[");
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, "]");
		weights[i] = ReadThisLine(file, 1, 1.0, name, "CommonNeighConst::GetCommonNeighConst");
				
		if (weights[i] < 0)
			ChiSquared::calc_sigma = 1;
	
		
		
		//Determining the index of the partial calculated from the primary types 
		primary_partial[i] = (primary1[i] <=primary2[i] ? (primary1[i] * RunParams::ntypes - (primary1[i] * (primary1[i] + 1) / 2) + primary2[i]) : \
			(primary2[i] * RunParams::ntypes - (primary2[i] * (primary2[i] + 1) / 2) + primary1[i]));

		mod_const[i] = 1;//to calculate the initial distribution

		//calculate the total number of possible primary1 - primary2 pairs
		
		tot_primary_pairs += (primary1[i] == primary2[i] ? (RunParams::pnatoms[primary1[i]] * (RunParams::pnatoms[primary1[i]] - 1) / 2) : RunParams::pnatoms[primary1[i]] * RunParams::pnatoms[primary2[i]]);
		cumul_prim[i + 1] = tot_primary_pairs;
	}//end of i cycle for the constraints
	
	
	if (!CheckReadFileState(file, "CommonNeighConst::GetCosConst", datfilename))
		CleanExit();//loading failed

	if (conv_numb != NULL)
		delete [] conv_numb;
	delete [] name;

};

//--------------gets the static parameters------------------------------------
bool CommonNeighConst::GetCommonNeighConstFree(std::list<std::list<string> > pool )
{
	int i, jj, index;
	char *name, *conv_numb = NULL;
	char *conv_numb2 = NULL;
	bool retval = true;
	
	nconstraints = RunParams::ncommonneigh;
	consttype = &NeighbourList::consttype;
	
	SetArraysize(&nsec_neigh, nconstraints, "nsec_neigh", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&cumul_prim, nconstraints + 1, "cumul_prim", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&cumul_sec, nconstraints + 1, "cumul_sec", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&name, NAME_SIZE, "name", "CommonNeighConst::GetCommonNeighConstFree");

	//creating the static arrays
	SetArraysize(&primary_partial, nconstraints, "primary_partial", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&mod_const, nconstraints, "mod_const", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&primary1, nconstraints, "primary1", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&primary2, nconstraints, "primary2", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&dmin, nconstraints, "dmin", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&dmax, nconstraints, "dmax", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&minsq, nconstraints, "minsq", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&maxsq, nconstraints, "maxsq", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&fraction, nconstraints, "fraction", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&target_coord, nconstraints, "target_coord", "CommonNeighConst::GetCommonNeighConstFree");
	SetArraysize(&weights, nconstraints, "weights", "CommonNeighConst::GetCommonNeighConstFree");

	tot_sec = 0;
	cumul_sec[0] = 0;
	i = 0;
	for (std::list<std::list<string>>::iterator it2 = pool.begin(); it2 != pool.end() && retval; it2++, i++)//going through all the conc constraints
	{
		//Initialization and allocating memories
		nsec_neigh[i] = std::count_if(it2->begin(), it2->end(), [](string c) { return c.find(CheckKey("STYPE_FROM1_TO1_FROM2_TO2")) != string::npos; });
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &secondary, "secondary", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &dmin1, "dmin1", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &dmin2, "dmin2", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &dmax1, "dmax1", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &dmax2, "dmax2", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &min1sq, "min1sq", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &min2sq, "min2sq", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &max1sq, "max1sq", "CommonNeighConst::GetCommonNeighConstFree");
		tot_sec -= nsec_neigh[i];
		ResizeArray(&tot_sec, tot_sec + nsec_neigh[i], &max2sq, "max2sq", "CommonNeighConst::GetCommonNeighConstFree");

		int si = 0, pt = 0,  nsec=0;
		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)//going through the lines of the current constraint
		{

			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("PTYPE1_PTYPE2_FROM_TO"))
			{
				if (pt++ > 0)
				{
					cout << "\nERROR: Multiple " << CheckKey("PTYPE1_PTYPE2_FROM_TO") << " definition in the " << i + 1 << ". "<<CheckTag("CONC")<<" constraint!" << endl;
					return false;
				}
				string value;
				if (!GetInt(ss, *it, primary1+i,"CommonNeighConst::GetCommonNeighConstFree"))  return ErrLackValue(*it);
				if (!GetInt(ss, *it, primary2+i,"CommonNeighConst::GetCommonNeighConstFree"))  return ErrLackValue(*it);
				if (!(ss >>  dmin[i] >> dmax[i])) return ErrLackValue(*it);
				IntToStr(&conv_numb, i + 1);
				mystrcpy(name, NAME_SIZE, "primary1[");
				mystrcat(name, NAME_SIZE, conv_numb);
				mystrcat(name, NAME_SIZE, "]");
				RunParams::CheckType(primary1[i], name, "CommonNeighConst::GetCommonNeighConstFree");
				primary1[i]--;
				mystrcpy(name, NAME_SIZE, "primary2[");
				mystrcat(name, NAME_SIZE, conv_numb);
				mystrcat(name, NAME_SIZE, "]");
				RunParams::CheckType(primary2[i], name, "CommonNeighConst::GetCommonNeighConstFree");
				primary2[i]--;
				//The concept is that the neighbourlist is calculated for atoms of each type represented in a cosine distribution or common neighbour
				//constraint. It is required for the neighbours as well, as the list is updated, and not completly recalculated at each step.
				//As each type can be involved in more, than one constraint, all the neighbourlist will be calculated with the
				//same maximum distance, maxdistsq. 
				//Be aware that the central and the neighbour type in comconsttype types represented by right to left!
				//type 0 is represented by binary 1, type 1 binary 10 and so on  
				//consttype i a pointer to NeighbourList::consttype for all the constarins using the neighbour list
				//comconsttype only for cos constraints
				if (!(*consttype >> primary1[i] & 1))//this type is not yet represented
				{
					*consttype += (longint)pow(2.0, (double)primary1[i]);
					NeighbourList::nconsttype++;//increase the number of constrained types
				}
				if (!(*consttype >> primary2[i] & 1))//this type is not yet represented
				{
					*consttype += (longint)pow(2.0, (double)primary2[i]);
					NeighbourList::nconsttype++;//increase the number of constrained types
				}
				if (!(comconsttype >> primary1[i] & 1))//this type is not yet represented
					comconsttype += (longint)pow(2.0, (double)primary1[i]);

				if (!(comconsttype >> primary2[i] & 1))//this type is not yet represented
					comconsttype += (longint)pow(2.0, (double)primary2[i]);

				if (dmin[i] > dmax[i])
				{
					cout << "\nWARNING(" << ++warn << "): the minimum distance value should be less than the maximum distance in the following line, they will be exchanged!" << endl << *it << endl;
					double temp = dmin[i];
					dmin[i] = dmax[i];
					dmax[i] = temp;
				}

				//calculating the reduced squared distances
				minsq[i] = pow(dmin[i] / RunParams::boxedge, 2.);
				maxsq[i] = pow(dmax[i] / RunParams::boxedge, 2.);
				//finding the largest reduced maximum distance square
				if (maxsq[i] > maxdistsq) maxdistsq = maxsq[i];
			}
			else if (key == CheckKey("STYPE_FROM1_TO1_FROM2_TO2"))
			{

				index = cumul_sec[i]+nsec;
				if (!GetInt(ss, *it, secondary+index,"CommonNeighConst::GetCommonNeighConstFree")) return ErrLackValue(*it);
				if (!(ss >> dmin1[index] >> dmax1[index])) return ErrLackValue(*it);
				if (!(ss >> dmin2[index] >> dmax2[index]))
				{
					if (::debug)
					{
						cout << "\nWARNING(" << ++warn << "): No minimum and maximum values were given for the primary2-secondary distance in the " << i + 1 << ". " << CheckTag("COS") << " constraint!" << endl;
						cout << "\tThe values for the primary1-secondary distance will be used!" << endl;
					}
					dmin2[index] = dmin1[index];
					dmax2[index] = dmax1[index];

				}
				IntToStr(&conv_numb2, index + 1);
				mystrcpy(name, NAME_SIZE, "secondary[");
				mystrcat(name, NAME_SIZE, conv_numb);
				mystrcat(name, NAME_SIZE, ",");
				mystrcat(name, NAME_SIZE, conv_numb2);
				mystrcat(name, NAME_SIZE, "]");
				RunParams::CheckType(secondary[index], name, "CommonNeighConst::GetCommonNeighConstFree");
				secondary[index]--;
				if (!(*consttype >> secondary[index] & 1))//this type is not yet represented
				{
					*consttype += (longint)pow(2.0, (double)secondary[index]);
					NeighbourList::nconsttype++;//increase the number of constrained types
				}
				if (!(comconsttype >> secondary[index] & 1))//this type is not yet represented
					comconsttype += (longint)pow(2.0, (double)secondary[index]);

				if (dmin1[index] > dmax1[index])
				{
					cout << "\nWARNING(" << ++warn << "): the 2. (dmin1) value should be less than the 3. (dmax1) one in the following line, they will be exchanged! " << endl << *it << endl;
					double temp = dmin1[index];
					dmin1[index] = dmax1[index];
					dmax1[index] = temp;
				}
				if (dmin2[index] > dmax2[index])
				{
					cout << "\nWARNING(" << ++warn << "): the 4. (dmin2) value should be less than the 5. (dmax2) one in the following line, they will be exchanged! " << endl << *it << endl;
					double temp = dmin2[index];
					dmin2[index] = dmax2[index];
					dmax2[index] = temp;
				}
				for (int jj = 0; jj < index - cumul_sec[i]; jj++)//check, whether this type was not already given as secondary
				{
					if (secondary[index] == secondary[cumul_sec[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmin1[index] != dmin1[cumul_sec[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
							cout << "The number of secondary types cannot be reduced, as different minimum primary1-secondary distances were given!" << endl;
							cout << "The distances for the " << jj + 1 << ". secondary and primary1 is " << dmin1[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmin1[index] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
					}
				}
				for (jj = 0; jj < index - cumul_sec[i]; jj++)//check, whether this type was not already given as secondary
				{
					if (secondary[index] == secondary[cumul_sec[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmin2[index] != dmin2[cumul_sec[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
							cout << "The number of secondary types cannot be reduced, as different minimum primary2-secondary distances were given!" << endl;
							cout << "The distances for the " << jj + 1 << ". secondary and primary2 is " << dmin2[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmin2[index] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
					}
				}
				//check, whether in case of the same atom type for both primaries, the dmin1 and dmin2 are the same
				if (primary1[i] == primary2[i] && dmin1[index] != dmin2[index])
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The atom types of the two primary are the same, but the secondary-primary1 and secondary-primary2 minimum distance is different." << endl;
					cout << "The secondary-primary1 minimum distance is " << dmin1[index] << "A, and secondary-primary2 is " << dmin2[index] << " A!" << endl;
					cout << "This does not makes sense, give the same value for both, correct the " << datfilename << " file! Exiting..." << endl;
					CleanExit();
				}
				for (jj = 0; jj < index - cumul_sec[i]; jj++)//check, whether this type was not already given as secondary
				{
					if (secondary[index] == secondary[cumul_sec[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmax1[index] != dmax1[cumul_sec[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
							cout << "The number of secondary types cannot be reduced, as different maximum primary1-secondary distances were given!" << endl;
							cout << "The distances for the " << jj + 1 << ". secondary and primary1 is " << dmax1[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmax1[index] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
					}
				}
				for (jj = 0; jj < index - cumul_sec[i]; jj++)//check, whether this type was not already given as secondary
				{
					if (secondary[index] == secondary[cumul_sec[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmax2[index] != dmax2[cumul_sec[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary types are the same!" << endl;
							cout << "The number of secondary types cannot be reduced, as different maximum primary2-secondary distances were given!" << endl;
							cout << "The distances for the " << jj + 1 << ". secondary and primary2 is " << dmax2[cumul_sec[i] + jj] << " A, while for " << index - cumul_sec[i] + 1 << ". secondary and primary1 is " << dmax2[index] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different secondary atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
						//only the secondary type is different, the distance parameters are the same, can continue by reducing the number of secondary types 
						cout << "\nWARNING(" << ++warn << "): The atom type for the " << i + 1 << ". common neighbour constraint's " << jj + 1 << ". and " << index - cumul_sec[i] + 1 << ". secondary type is the same!" << endl;
						cout << "\tThe number of secondary types will be reduced!" << endl;
						nsec_neigh[i]--;
					}
				}
				//check, whether in case of the same atom type for both primaries, the dmax1 and dmax2 are the same
				if (primary1[i] == primary2[i] && dmax1[index] != dmax2[index])
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The atom types of the two primaryy are the same, but the secondary-primary1 and secondary-primary2 maximum distance is different." << endl;
					cout << "The secondary-primary1 maximum distance is " << dmax1[index] << "A, and secondary-primary2 is " << dmax2[index] << " A!" << endl;
					cout << "This does not makes sense, give the same value for both, correct the " << datfilename << " file! Exiting..." << endl;
					CleanExit();
				}
				//Seems to be a valid secondary neighbour
				min1sq[index] = pow(dmin1[index] / RunParams::boxedge, 2.);//calculating the reduced squared distances
				min2sq[index] = pow(dmin2[index] / RunParams::boxedge, 2.);
				max1sq[index] = pow(dmax1[index] / RunParams::boxedge, 2.);
				max2sq[index] = pow(dmax2[index] / RunParams::boxedge, 2.);
				if (max1sq[index] > maxdistsq)//finding the largest reduced maximum distance square
					maxdistsq = max1sq[index];
				if (max2sq[index] > maxdistsq)//finding the largest reduced maximum distance square
					maxdistsq = max2sq[index];
				nsec++;
			}
			else if ((key == CheckKey("COORDNUM_FRACT_SIGMA")) || (key == CheckKey("COORDNUM_FRACT_SIGMA-MASTER")) || (key == CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE")))
			{
				if (si++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple ...SIGMA(-MASTER)/...SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckTag("CONC") << " constraint!" << endl << "\tThe last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + CoordNumbConst::tot_subconst + RunParams::navcoord + i + 1)
						RunParams::lead_series_ind = -1;
				}
				if (!GetInt(ss, *it, target_coord+i, "CommonNeighConst::GetCommonNeighConstFree")) return ErrLackValue(*it);
				if (!(ss  >> fraction[i] >> weights[i])) return ErrLackValue(*it);
				if (key == CheckKey("COORDNUM_FRACT_SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1)
						RunParams::lead_series_ind = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + CoordNumbConst::tot_subconst + RunParams::navcoord + i + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry is declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE"))
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
				CheckKeyTag(key, CheckTag("CONC"));//write the appropriate error mesage
				return false;
			}
			cumul_sec[i + 1] = cumul_sec[i] + nsec_neigh[i];
		}
		if (retval)
		{
			//Determining the index of the partial calculated from the primary types 
			primary_partial[i] = (primary1[i] <= primary2[i] ? (primary1[i] * RunParams::ntypes - (primary1[i] * (primary1[i] + 1) / 2) + primary2[i]) : \
				(primary2[i] * RunParams::ntypes - (primary2[i] * (primary2[i] + 1) / 2) + primary1[i]));

			mod_const[i] = 1;//to calculate the initial distribution

			//calculate the total number of possible primary1 - primary2 pairs

			tot_primary_pairs += (primary1[i] == primary2[i] ? (RunParams::pnatoms[primary1[i]] * (RunParams::pnatoms[primary1[i]] - 1) / 2) : RunParams::pnatoms[primary1[i]] * RunParams::pnatoms[primary2[i]]);
			cumul_prim[i + 1] = tot_primary_pairs;
		}
	}

	if (conv_numb != NULL)
		delete [] conv_numb;
	delete [] name;
	return retval;

};
//-----------initializing some of the static members-----------------------
void CommonNeighConst::SetParams()
{

	max_neigh = &NeighbourList::max_neigh;
	

};
//-------------------calculate the initial histogram--------------------
void CommonNeighConst::CalcHist()
{
	//the pointers and auxiliary variables are used to extract array values only ones to make calculation hopefully quicker
	int iconst, isec, i, inb,inb1, inb2;
	int imin, imax;//cycle boundaries
	int *pnneigh,*pnneigh2;//pointer to the number of neighbours 
	int *plist, *pind, *pind1, *pind2;//pointers to the indices of the neighbours
	int *ptypelist, *pneigh_type, *pneigh_type1, *pneigh_type2;//pointers to the type of the neighbours 
	int count;//common neighbour count for a primary pair
	int offset;

	double dminsq, dmaxsq;//minimum and maximum squared distances for the two neighbours
	double *pdistsq, *pdsq, *pdsq1, *pdsq2;//pointers to the distances of the neighbours
	
	//First find the primary1 -primary2 pairs taking primary1 as central. If a pair found, go through the neighbours of them, and find the common one(s).
	//A histogam count is added for each primary pair, to the zero bin, if there is no common, to the ith bin, with i common neighbours
	//nprimary is the total number of primary pairs for the given constraint, which is needed to calculate the fraction satisfying the constraint. 
	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		//setting auxiliary variables and pointers to increase the speed
		dminsq = minsq[iconst];
		dmaxsq = maxsq[iconst];
		nprimary[iconst] =0;

		pnneigh = neighlist.nneigh_finder[primary1[iconst]];//sets the pointer to the beginning of the values belonging to the primary1 atom type in nneigh array
		plist = neighlist.neighlist_finder[primary1[iconst]];//sets the pointer to the beginning of the values belonging to the primary1 atom type in neigh_list array
		ptypelist = neighlist.neightype_finder[primary1[iconst]];//sets the pointer to the beginning of the values belonging to the primary1 atom type in neigh_type array 
		pdistsq = neighlist.dsq_finder[primary1[iconst]];//sets the pointer to the beginning of the values belonging to the primary1 atom type in dsq array
		imin = SimpleCfg::cumul[primary1[iconst]];
		imax = SimpleCfg::cumul[primary1[iconst] + 1];
		for (i = imin; i < imax; i++)//primary1 atom cycle
		{
			pneigh_type = ptypelist;//type of the first neighbour, possible primary2
			pind = plist;//index of the neighbour
			pdsq = pdistsq;//squared distance of the neighbour
			for (inb = 0; inb < *pnneigh; inb++)//neighbour cycle to find primary2
			{
				//determine, whether the neighbour is from the right type, and is in the range
				if (*pneigh_type == primary2[iconst] && *pdsq >= dminsq && *pdsq <= dmaxsq)
				{
					if ((primary1[iconst] != primary2[iconst]) || (primary1[iconst] == primary2[iconst] && i < *pind))
					{//to make sure, that in case of the same types, a pair is included only once, when the primary1's index is smaller
						nprimary[iconst]++;
						pneigh_type1 = ptypelist;//type of the first neighbour, possible secondary
						pind1 = plist;//index of the neighbour
						pdsq1 = pdistsq;//squared distance of the neighbour
						count = 0;
						for (inb1 = 0; inb1 < *pnneigh; inb1++)//going through again the neighbourlist of the primary1
						{
							if (*pind != *pind1)//this is not primary2
							{
								offset = *pind - SimpleCfg::cumul[primary2[iconst]];
								pnneigh2 = neighlist.nneigh_finder[primary2[iconst]]+offset;//number of neighbours of the of primary2
								offset *= *max_neigh;
								pneigh_type2 = neighlist.neightype_finder[primary2[iconst]]+offset;//type of the neighbours of primary2
								pind2 = neighlist.neighlist_finder[primary2[iconst]]+offset;//index of the neighbours of primary2
								pdsq2 = neighlist.dsq_finder[primary2[iconst]]+offset;//squared distanceof of the neighbours of primary2
								for (inb2 = 0; inb2 < *pnneigh2; inb2++)//going through the neighbourlist of the primary2
								{
									if (*pind2 != i)//this is not primary1
									{
										for (isec = 0; isec < nsec_neigh[iconst]; isec++)
										{
											if (*pind1 == *pind2 && *pneigh_type1 == secondary[cumul_sec[iconst] + isec] && \
												*pdsq1 >= min1sq[cumul_sec[iconst] + isec] && *pdsq1 <= max1sq[cumul_sec[iconst] + isec] && \
												*pdsq2 >= min2sq[cumul_sec[iconst] + isec] && *pdsq2 <= max2sq[cumul_sec[iconst] + isec])//this is a common neighbour
											{
												count++;
											}
										}

									}
									pneigh_type2++;
									pind2++;
									pdsq2++;
								}

							}
							pneigh_type1++;
							pind1++;
							pdsq1++;
						}
						hist[iconst*max_hist_size + count]++;//update the histogram
					}//this is a new primary pair
				}
				pneigh_type++;
				pind++;
				pdsq++;
			}//end of neighbour cycle inb1
			pnneigh++;
			plist += *max_neigh;
			ptypelist += *max_neigh;
			pdistsq += *max_neigh;
		}//end of central atom cycle i for primary1
		nsatisfy[iconst] = hist[iconst*max_hist_size + target_coord[iconst]];
	}//end of iconst cycle
};


//----updating the common neighbour histogram, if sign switch is -1: removing contribution, if +1 adding contribution--------------------
void CommonNeighConst::UpdateHist(Move& move, int sign_switch)
{
	int iconst, isec, isec2, imoved, jmoved, inb, inb1, inb2, is1, is2, icoord;
	int index;//index for the secondary types
	int fixed_neigh;
	int partialind;
	int mytype2;//actual type of the 
	int myind1 = 0, myind2 = 0;//actual index of the primaries in its own type 
	int coffset;
	int skipcycle;//boolean indicator
	int* moved_type, * moved_ind;
	int* pnneigh, * pnneigh2, * psnneigh1, * psnneigh2;//pointer to the number of neighbours 
	int* pind, * pind1, * pind2, * psind1, * psind2;//pointers to the indices of the neighbours
	int* pneigh_type, * pneigh_type1, * pneigh_type2, * ps_type1, * ps_type2;//pointers to the type of the neighbours 
	int count;

	double distsq, dminsq, dmaxsq;//minimum and maximum squared distances for the two neighbours
	double dminsq1, dminsq2, dmaxsq1, dmaxsq2;
	double* pdsq, * pdsq1, * pdsq2, * psdsq1, * psdsq2;//pointers to the distances of the neighbours
	double* pvcomp1, * pvcomp2;//poiters to the central-neibour vector components
	double* pmin1sq, * pmax1sq, * pmin2sq, * pmax2sq;//these are pointers for the appropriate min1sq, max1sq, min2sq and max2sq


	if (sign_switch == -1)//to do it only once in each loop
	{
		for (iconst = 0; iconst < nconstraints; iconst++)
			mod_const[iconst] = 0;//default, not modified
	}


	//first check, whether the moved atom is a primary
	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		//setting auxiliary variables and pointers to increase the speed
		dminsq = minsq[iconst];
		dmaxsq = maxsq[iconst];
		moved_type = move.types;
		moved_ind = move.indices;
		for (imoved = 0; imoved < Move::tot_moved_atoms; imoved++)
		{
			if (!(comconsttype >> *moved_type & 1))
			{
				moved_type++;
				moved_ind++;
				continue;//this moved atom is not involved in this constraint, go to next
			}
			fixed_neigh = 0;//the moved is not a primary in the constraint
			//setting the auxiliary variables the way, that index1 will belong to the fixed neighbour (the moved atom)
			if (*moved_type == primary1[iconst])//the moved atom can be a first primary in this constraint
			{
				mod_const[iconst] = 1;//this constraint is modified by the move 
				fixed_neigh = 1;//this is the first primary, this will be the fixed neighbour
				pmin1sq = min1sq;
				pmax1sq = max1sq;
				mytype2 = primary2[iconst];//this is for the other primary
				pmin2sq = min2sq;
				pmax2sq = max2sq;

			}
			else
			{
				if (*moved_type == primary2[iconst])//the moved atom is second primary in this constraint
				{
					mod_const[iconst] = 1;//this constraint is modified by the move 
					fixed_neigh = 2;//this is the second primary
					pmin1sq = min2sq;
					pmax1sq = max2sq;
					mytype2 = primary1[iconst];//this is for the other primary
					pmin2sq = min1sq;
					pmax2sq = max1sq;

				}
			}
			if (fixed_neigh > 0)//the moved is a possible primary
			{
				coffset = *moved_ind - config.cumul[*moved_type];//the offset of the central atom in its own type
				pnneigh = neighlist.nneigh_finder[*moved_type] + coffset;//sets the pointer to the nneigh of the central atom

				//setting the pointers for the first neighbour of the moved atom, possible primary
				coffset *= *max_neigh;//multiplying the offset with the maximum number of neighbours
				pind = neighlist.neighlist_finder[*moved_type] + coffset;//sets the pointer to the beginning of the neigh_list of the central atom
				pneigh_type = neighlist.neightype_finder[*moved_type] + coffset;//sets the pointer to the beginning of the neigh_type of the central atom
				pdsq = neighlist.dsq_finder[*moved_type] + coffset;//sets the pointer to the beginning of the squared distance of the central atom


				for (inb = 0; inb < *pnneigh; inb++)//neighbour cycle to find the other primary
				{
					//check, whether this is a moved atom situated later in the move.indices list
					skipcycle = 0;
					for (jmoved = imoved + 1;jmoved < Move::tot_moved_atoms;jmoved++)
					{
						if (move.indices[jmoved] == *pind)
						{
							//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
							pneigh_type++;
							pind++;
							pdsq++;
							skipcycle = 1;
							break;
						}
					}
					if (skipcycle)
						continue;//go to the next inb1
					//determine, whether the neighbour is from the right type, and is in the range
					if (*pneigh_type == mytype2 && *pdsq >= dminsq && *pdsq <= dmaxsq)
					{
						nprimary[iconst] += sign_switch;//this is a primary pair

						//going through again the neighbourlist of the primary1 to find secondary
						coffset = *moved_ind - config.cumul[*moved_type];//the offset of the central atom in its own type
						coffset *= *max_neigh;//multiplying the offset with the maximum number of neighbours
						pind1 = neighlist.neighlist_finder[*moved_type] + coffset;//sets the pointer to the beginning of the neigh_list of the central atom
						pneigh_type1 = neighlist.neightype_finder[*moved_type] + coffset;//sets the pointer to the beginning of the neigh_type of the central atom
						pdsq1 = neighlist.dsq_finder[*moved_type] + coffset;//sets the pointer to the beginning of the squared distance of the central atom

						count = 0;
						for (inb1 = 0; inb1 < *pnneigh; inb1++)
						{
							//check, whether the neighbour is a moved atom situated later in the move.indices list
							skipcycle = 0;
							for (jmoved = imoved + 1;jmoved < Move::tot_moved_atoms;jmoved++)
							{
								if (move.indices[jmoved] == *pind1)
								{
									//this will be evaluated for the jmoved-th moved atom to prevent the same pair calculated multiple times 
									pneigh_type1++;
									pind1++;
									pdsq1++;
									skipcycle = 1;
									break;
								}
							}
							if (skipcycle)
								continue;//go to the next inb1
							if (*pind != *pind1)//this is not the other primary
							{
								coffset = *pind - SimpleCfg::cumul[mytype2];
								pnneigh2 = neighlist.nneigh_finder[mytype2] + coffset;//number of neighbours of the of primary2
								coffset *= *max_neigh;
								pneigh_type2 = neighlist.neightype_finder[mytype2] + coffset;//type of the neighbours of primary2
								pind2 = neighlist.neighlist_finder[mytype2] + coffset;//index of the neighbours of primary2
								pdsq2 = neighlist.dsq_finder[mytype2] + coffset;//squared distanceof of the neighbours of primary2
								for (inb2 = 0; inb2 < *pnneigh2; inb2++)//going through the neighbourlist of the primary2
								{
									if (*pind2 != *moved_ind)//this is not the moved
									{
										for (isec = 0; isec < nsec_neigh[iconst]; isec++)
										{
											if (*pind1 == *pind2 && *pneigh_type1 == secondary[cumul_sec[iconst] + isec] && \
												* pdsq1 >= pmin1sq[cumul_sec[iconst] + isec] && *pdsq1 <= pmax1sq[cumul_sec[iconst] + isec] && \
												* pdsq2 >= pmin2sq[cumul_sec[iconst] + isec] && *pdsq2 <= pmax2sq[cumul_sec[iconst] + isec])//this is a common neighbour
											{
												count++;

											}
										}

									}
									pneigh_type2++;
									pind2++;
									pdsq2++;
								}//end of inb2 cycle 

							}
							pneigh_type1++;
							pind1++;
							pdsq1++;
						}//end of inb1 cycle
						hist[iconst * max_hist_size + count] += sign_switch;//update the histogram
					}
					pneigh_type++;
					pind++;
					pdsq++;
				}//end of neighbour cycle inb
			}//end of the moved can be a primary

			//now see, if the moved can be a secondary
			for (isec = 0; isec < nsec_neigh[iconst]; isec++)
			{
				index = cumul_sec[iconst] + isec;
				if (*moved_type == secondary[index])
				{
					mod_const[iconst] = 1;//this constraint is modified by the move 
						//setting auxiliary variables and pointers to increase the speed
					dminsq1 = min1sq[index];
					dminsq2 = min2sq[index];
					dmaxsq1 = max1sq[index];
					dmaxsq2 = max2sq[index];
					coffset = *moved_ind - config.cumul[*moved_type];//the offset of the central (now secondary atom) in its own type
					pnneigh = neighlist.nneigh_finder[*moved_type] + coffset;//sets the pointer to the nneigh of the central atom

					//setting the pointers for the first neighbour, possible primary
					coffset *= *max_neigh;//multiplying the offset with the maximum number of neighbours
					pind1 = neighlist.neighlist_finder[*moved_type] + coffset;//sets the pointer to the beginning of the neigh_list of the central atom
					pneigh_type1 = neighlist.neightype_finder[*moved_type] + coffset;//sets the pointer to the beginning of the neigh_type of the central atom
					pdsq1 = neighlist.dsq_finder[*moved_type] + coffset;//sets the pointer to the beginning of the squared distance of the central atom
					pvcomp1 = neighlist.vector_finder[*moved_type] + coffset * 3;//sets the pointer to the beginning of the central-neighbour1 vector components of the central atom
					//going through all the neighbours
					for (inb1 = 0; inb1 < (*pnneigh - 1); inb1++)//first neighbour cycle
					{
						//check, whether the neighbour is a moved atom situated later in the move.indices list
						skipcycle = 0;
						for (jmoved = imoved + 1; jmoved < Move::tot_moved_atoms; jmoved++)
						{
							if (move.indices[jmoved] == *pind1)
							{
								//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
								pneigh_type1++;
								pind1++;
								pdsq1++;
								pvcomp1 += 3;
								skipcycle = 1;
								break;
							}
						}
						if (skipcycle)
							continue;//go to the next inb1
						//determine, whether the first neighbour is from the right type, and is in the range, so it is a primary
						if ((*pneigh_type1 == primary1[iconst] && *pdsq1 >= dminsq1 && *pdsq1 <= dmaxsq1) || \
							(*pneigh_type1 == primary2[iconst] && *pdsq1 >= dminsq2 && *pdsq1 <= dmaxsq2))
						{
							//							mytype1 = (*pneigh_type1 == primary1[iconst] ? primary1[iconst] : primary2[iconst]);
							if (*pneigh_type1 == primary1[iconst])
								myind1 = *pind1 - SimpleCfg::cumul[primary1[iconst]];
							else
								myind2 = *pind1 - SimpleCfg::cumul[primary2[iconst]];
							pneigh_type2 = pneigh_type1 + 1;//type of the second neighbour
							pind2 = pind1 + 1;//index of the second neighbour
							pdsq2 = pdsq1 + 1;//squared distance of the central atom and the second neighbour 
							pvcomp2 = pvcomp1 + 3;//vector components of the central atom and the second neighbour 
							for (inb2 = inb1 + 1; inb2 < *pnneigh; inb2++)
							{
								//check, whether the neighbour is a moved atom situated later in the move.indices list
								skipcycle = 0;
								for (jmoved = imoved + 1; jmoved < Move::tot_moved_atoms; jmoved++)
								{
									if (move.indices[jmoved] == *pind2)
									{
										//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
										pneigh_type2++;
										pind2++;
										pdsq2++;
										pvcomp2 += 3;
										skipcycle = 1;
										break;
									}
								}
								if (skipcycle)
									continue;//go to the next inb2
								//determining the index of the partial calculated from the neighbours					
								partialind = (*pneigh_type1 <= *pneigh_type2 ? (*pneigh_type1 * RunParams::ntypes - (*pneigh_type1 * (*pneigh_type1 + 1) / 2) + *pneigh_type2) : \
									(*pneigh_type2 * RunParams::ntypes - (*pneigh_type2 * (*pneigh_type2 + 1) / 2) + *pneigh_type1));
								if (partialind == primary_partial[iconst])//both neighbours are from the right types, these can be a primary pair
								{
									//determine, whether the second neighbour is in the range
									if ((*pneigh_type2 == primary1[iconst] && *pdsq2 >= dminsq1 && *pdsq2 <= dmaxsq1) || \
										(*pneigh_type2 == primary2[iconst] && *pdsq2 >= dminsq2 && *pdsq2 <= dmaxsq2))
									{
										//mytype2 = (*pneigh_type2 == primary1[iconst] ? primary1[iconst] : primary2[iconst]);
										if (*pneigh_type2 == primary2[iconst])
											myind2 = *pind2 - SimpleCfg::cumul[primary2[iconst]];
										else
											myind1 = *pind2 - SimpleCfg::cumul[primary1[iconst]];
										//calculating the distance between neighours, as the length of pvcomp1-pvcomp2 vector
										distsq = 0;
										for (icoord = 0; icoord < 3; icoord++)
											distsq += pow(*(pvcomp1 + icoord) - *(pvcomp2 + icoord), 2);
										if (distsq >= dminsq && distsq <= dmaxsq) //these are primaries, this is a valid triplett, update the count
										{
											count = 1;
											//must find all the common neighbours of this primary pair
											//first go through all the neighbours of the first primary
											psnneigh1 = neighlist.nneigh_finder[primary1[iconst]] + myind1;//sets the pointer to the nneigh of the central atom
											coffset = myind1 * *max_neigh;//multiplying the offset with the maximum number of neighbours
											psind1 = neighlist.neighlist_finder[primary1[iconst]] + coffset;//sets the pointer to the beginning of the neigh_list of the central atom
											ps_type1 = neighlist.neightype_finder[primary1[iconst]] + coffset;//sets the pointer to the beginning of the neigh_type of the central atom
											psdsq1 = neighlist.dsq_finder[primary1[iconst]] + coffset;//sets the pointer to the beginning of the squared distance of the central atom

											for (is1 = 0; is1 < *psnneigh1; is1++)//going through again the neighbourlist of the primary1
											{
												//check, whether the this a moved atom situated later in the move.indices list
												skipcycle = 0;
												for (jmoved = imoved + 1; jmoved < Move::tot_moved_atoms; jmoved++)
												{
													if (move.indices[jmoved] == *pind2)
													{
														//this will be evaluated for the jmoved-th moved atom to prevent the same triplet to be calculated multiple times 
														ps_type1++;
														psind1++;
														psdsq1++;
														skipcycle = 1;
														break;
													}
												}
												if (*psind1 != *pind1 && *psind1 != *pind2 && *psind1 != *moved_ind)//this is not primary or the moved atom
												{
													psnneigh2 = neighlist.nneigh_finder[primary2[iconst]] + myind2;//number of neighbours of the of primary2
													coffset = myind2 * *max_neigh;
													ps_type2 = neighlist.neightype_finder[primary2[iconst]] + coffset;//type of the neighbours of primary2
													psind2 = neighlist.neighlist_finder[primary2[iconst]] + coffset;//index of the neighbours of primary2
													psdsq2 = neighlist.dsq_finder[primary2[iconst]] + coffset;//squared distanceof of the neighbours of primary2
													for (is2 = 0; is2 < *psnneigh2; is2++)//going through the neighbourlist of the primary2
													{
														//have to check all the possible secondary types
														for (isec2 = 0; isec2 < nsec_neigh[iconst]; isec2++)
														{
															if (*psind2 == *psind1 && *ps_type1 == secondary[cumul_sec[iconst] + isec2] && \
																* psdsq1 >= min1sq[cumul_sec[iconst] + isec2] && *psdsq1 <= max1sq[cumul_sec[iconst] + isec2] && \
																* psdsq2 >= min2sq[cumul_sec[iconst] + isec2] && *psdsq2 <= max2sq[cumul_sec[iconst] + isec2])//this is a common neighbour
															{
																count++;
															}
														}
														ps_type2++;
														psind2++;
														psdsq2++;
													}//end of is2 cycle
												}

												ps_type1++;
												psind1++;
												psdsq1++;
											}//end of is1 cycle

											hist[iconst * max_hist_size + count] += sign_switch;//this triplett has to be removed for sign=-1, added is sign=+1
											//the count has to be at least 1 (the mmoved as secondary, the primary pair is not effected
											hist[iconst * max_hist_size + count - 1] -= sign_switch;//the other common neighbours for this primary are not effected
										}//end of this is a triplett
									}//second neighbour is in the range
								}//both neighbours are from the right type

								pneigh_type2++;
								pind2++;
								pdsq2++;
								pvcomp2 += 3;
							}//end of second neighbour cycle inb2
						}//first neighbour is from the right type, and is in the range

						pneigh_type1++;
						pind1++;
						pdsq1++;
						pvcomp1 += 3;
					}//end of first neighbour cycle inb1
				}//end of is moved can be a secondary
			}//end of secondary cycle isec
			moved_type++;
			moved_ind++;
		}//end of moved atom cycle imoved 
		nsatisfy[iconst] = hist[iconst * max_hist_size + target_coord[iconst]];
	}//end of iconst cycle

};


//-----------copying the modified histogram parts of the cosine distribution of the bond angles---------------
void CommonNeighConst::CopyModified(CommonNeighConst &target)
{
	int iconst, i;
	int *pi_source, *pi_target;//pointer to the source and target histogram 
	
	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		if (mod_const[iconst])
		{
			//the histogram of this constraint was modified, copy it 
			pi_source = hist+iconst*max_hist_size;
			pi_target = target.hist+iconst*max_hist_size;
			for (i = 0; i < max_hist_size; i++)
				*pi_target++ = *pi_source++;
			
			target.nsatisfy[iconst] = nsatisfy[iconst];
			target.nprimary[iconst] = nprimary[iconst];
		}
	};
};
void CommonNeighConst::SaveHist(const char *file_name) const
{
	int i, j;
	int *pcount;
	ofstream file;

	if (strlen(file_name) != 0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, concfilename);
	OpenFile(file, tempfilename, "CommonNeighConst::Save", 0);//open file, check, whether it was successfully opened


	file.precision(6);
	file.setf(ios::fixed, ios::floatfield);
	file << "Common neighbour constraint" << endl;
	file << nconstraints << "\t number of constraints (nconstraints) " << endl;

	file << "\nConstraints properties: " << endl;

	file << "index, primary type1, primary type2, primary dmin, primary dmax, number of secondary types, secondary type for each constraint, " << endl;
	file << "minimal primary1-secondary distance for each secondary type, minimal primary2 - secondary distance for each secondary type, " << endl;
	file << "maximal primary1-secondary distance for each secondary type, maximal primary2 - secondary distance for each secondary type, " << endl;
	file << "desired coordination number, fraction, weight " << endl;

	for (i = 0; i < nconstraints; i++)
	{
		file << J4 << i + 1 << "\t";


		file << J4 << primary1[i] + 1 << "\t";//start with 1 in the file
		file << J4 << primary2[i] + 1 << "\t";//start with 1 in the file
		file << J10 << dmin[i] << "\t";
		file << J10 << dmax[i] << "\t";
		file << J4 << nsec_neigh[i] << "\t";
		for (j = 0; j < nsec_neigh[i]; j++)
			file << J4 << secondary[cumul_sec[i]+j] + 1 << "\t";//start with 1 in the file
		for (j = 0; j < nsec_neigh[i]; j++)
			file << J10 << dmin1[cumul_sec[i] + j] << "\t";
		for (j = 0; j < nsec_neigh[i]; j++)
			file << J10 << dmin2[cumul_sec[i] + j] << "\t";
		for (j = 0; j < nsec_neigh[i]; j++)
			file << J10 << dmax1[cumul_sec[i] + j] << "\t";
		for (j = 0; j < nsec_neigh[i]; j++)
			file << J10 << dmax2[cumul_sec[i] + j] << "\t";
		file << J4 << target_coord[i] << "\t";
		file << J10 << fraction[i] << "\t";
		file << J10 << weights[i] << "\t";
		file << endl;
	}

	file << endl;
	file << "Common neighbour histogram for the actual primary pairs" << endl;
	pcount = hist;
	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::fixed, ios::floatfield);
	file.precision(8);
	for (i = 0; i < nconstraints; i++)
	{
		file << "\n" << i + 1 << ". constraint" << endl;
		file << "Number of primary  pairs: " << nprimary[i] << endl;
		file << "Number of primary  pairs with "<<target_coord[i]<<" common neighbours: " << hist[i*max_hist_size+target_coord[i]] << endl;
		file << "\nbin index \tcount" << endl;

		for (j = 0; j < max_hist_size; j++)
		{
			file << j  << "\t" << *pcount++ << endl;
		}
	}
	file.unsetf(ios::right);
	file.unsetf(ios::fixed);
	file.precision(6);
	file.close();
};

#endif // _ADVANCED_GEOM_CONST