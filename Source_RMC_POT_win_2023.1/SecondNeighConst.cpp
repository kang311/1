//source SecondNeighConst.cpp
//Last changed 12.12.2022



#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Move.h"//classes1.h is included through this

#ifdef _ADVANCED_GEOM_CONST
int  SecondNeighConst::nconstraints;//the number of constraints
int  SecondNeighConst::ncentral_tot;//total number of central atoms for all the constraints
int  SecondNeighConst::nthreads;//the total number of threads to use = RunParams::nthreads
int  SecondNeighConst::tot_sec;//total number of second neighbours for all the constraints
int *SecondNeighConst::nsec_type = NULL;//number of second neighbour types for each costraint
int *SecondNeighConst::nsec_type_cum;//cumulative number of second neighbour types for each constraint
int *SecondNeighConst::max_neigh;//pointer to the maximum number of neighbours for the neighbourlist for an atom

int *SecondNeighConst::ncentral;//the number of central atoms for each constraints
int *SecondNeighConst::ncentral_cum;//cumulative number of central atoms before the given constraint for each constraint
int *SecondNeighConst::central;//array of types of central particles [nconstraints]
int *SecondNeighConst::fneighbour;//array of types of first neighbour particles  [nconstraints]
int *SecondNeighConst::sneighbours;//array of types of second neighbour particles  [tot_sec]
int *SecondNeighConst::target_coord;//array of DESIRED coord. numbers [nconstraints]

int *SecondNeighConst::mod_const;//whteher the constraint was effected by the move
int *SecondNeighConst::cp_partial;//index of the partial calculated from the central and first neighbour types

int SecondNeighConst::involved_count;//number of involved indices
int *SecondNeighConst::involved_indices=NULL;//indices of the atoms to recalculate for (the moved atoms and their first and second neighbours)

longint *SecondNeighConst::consttype;//pointer to NeighbourList::consttype, each bit represents a type, containing a 1 for each type that is involved in a constraint
longint SecondNeighConst::secconsttype=0;//each bit represents a type, containing a 1 for each type that is involved in a second neigh constraint

double *SecondNeighConst::dmin1;//array of minimal distances between central and first neighbour [nconstraints]
double *SecondNeighConst::dmax1;//array of maximum distances between central and first neighbour [nconstraints]
double *SecondNeighConst::dmin2;//array of minimal distances between first and second neighbours [tot_sec]
double *SecondNeighConst::dmax2;//array of maximum distances between first and second neighbours [tot_sec]
double *SecondNeighConst::weights;//array of weights for each constraint [nconstraints]
double *SecondNeighConst::fraction;//desired fraction for each constraint [nconstraints]

double *SecondNeighConst::min1sq;//reduced min. dist. squared
double *SecondNeighConst::min2sq;//reduced min. dist. squared
double *SecondNeighConst::max1sq;//reduced max. dist. squared
double *SecondNeighConst::max2sq;//reduced max. dist. squared
double  SecondNeighConst::maxdistsq;//the maximum distance among the maxsq, max1sq and max2sq arrays, this is used determining the limit of the neighbourlist calculation


//----default constructor----------------
SecondNeighConst::SecondNeighConst(NeighbourList &neigh_list, SimpleCfg &conf)
	:neighlist(neigh_list), config(conf)
{
	int i;
	if (nconstraints < 1)
	{
		if (::debug)
		{
			cout << "\nWARNING(" << ++warn << "): SecondNeighConst constructor" << endl;
			cout << "\tThere are no second neighbour constraints, or nconstraints was not initialised!" << endl;
		}
	}
	//allocating memory
	SetArraysize(&nsatisfy, nconstraints, "nsatisfy", "SecondNeighConst::SecondNeighConst");

	SetArraysize(&coordnumbs, ncentral_tot, "coordnumbs", "SecondNeighConst::SecondNeighConst");
	SetArraysize(&finder, nconstraints, "finder", "SecondNeighConst::SecondNeighConst");
	
	InitCoordnumbs();//Sets the array elements to zero

	//initialising the finder
	for (i = 0; i < nconstraints; i++)
		finder[i] = coordnumbs + ncentral_cum[i];
	SetArraysize(&nsatisfy, nconstraints, "nsatisfy", "SecondNeighConst::SecondNeighConst");
	for (i = 0; i < nconstraints; i++)
		nsatisfy[i] = 0;
	
};

//copy constructor
SecondNeighConst::SecondNeighConst(SecondNeighConst &source, NeighbourList &neigh_list, SimpleCfg &conf)
	:neighlist(neigh_list), config(conf)
{
	int i;
	
	//allocating memory
	SetArraysize(&nsatisfy, nconstraints, "nsatisfy", "SecondNeighConst::SecondNeighConst");

	SetArraysize(&coordnumbs, ncentral_tot, "coordnumbs", "SecondNeighConst::SecondNeighConst");
	SetArraysize(&finder, nconstraints, "finder", "SecondNeighConst::SecondNeighConst");

	for (i = 0; i < ncentral_tot; i++)
		coordnumbs[i] = source.coordnumbs[i];

	//initialising the finder
	for (i = 0; i < nconstraints; i++)
		finder[i] = coordnumbs + ncentral_cum[i];
	SetArraysize(&nsatisfy, nconstraints, "nsatisfy", "SecondNeighConst::SecondNeighConst");
		
	for (i = 0; i < nconstraints; i++)
		nsatisfy[i] = source.nsatisfy[i];
		
};

//--------------------------------------------------------------------------------
//gets or sets the value for the iconst-th constraints, iatom-th 
//central atom in the coordnumbs array
//iatom index should be in its own type
int &SecondNeighConst::CoordinationNumb(int iconst, int iatom)
{
	return(*(finder[iconst] + iatom));
}

//------------Sets the array elements to zero------------------
void SecondNeighConst::InitCoordnumbs()
{
	int i;
		
	//setting all the array elements to 'zero' values
	for (i = 0; i < ncentral_tot; i++)
		coordnumbs[i] = 0;
	
	for (i = 0; i < nconstraints; i++)
		nsatisfy[i] = 0;
};
//-----------initializing some of the static members-----------------------
void SecondNeighConst::SetParams()
{

	max_neigh = &NeighbourList::max_neigh;
};

//--------------gets the static parameters------------------------------------
bool SecondNeighConst::GetSecondNeighConstFree(std::list<std::list<string> > pool )
{
	int i, index;
	char *name, *conv_numb = NULL;
	char *conv_numb2 = NULL;
	bool retval = true;

	nconstraints = RunParams::nsecondneigh;
	consttype = &NeighbourList::consttype;
	
	//creating the static arrays
	SetArraysize(&central, nconstraints, "central", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&fneighbour, nconstraints, "fneighbour", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&dmin1, nconstraints, "dmin1", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&dmax1, nconstraints, "dmax1", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&min1sq, nconstraints, "min1sq", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&max1sq, nconstraints, "max1sq", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&fraction, nconstraints, "fraction", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&target_coord, nconstraints, "target_coord", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&weights, nconstraints, "weights", "SecondNeighConst::GetSecondNeighConstFree");
	

	SetArraysize(&cp_partial, nconstraints, "cp_partial", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&mod_const, nconstraints, "mod_const", "SecondNeighConst::GetSecondNeighConst");
	SetArraysize(&ncentral, nconstraints, "ncentral", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&ncentral_cum, nconstraints + 1, "ncentral", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&nsec_type, nconstraints, "nsec_type", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&nsec_type_cum, nconstraints + 1, "nsec_type", "SecondNeighConst::GetSecondNeighConstFree");
	SetArraysize(&name, NAME_SIZE, "name", "SecondNeighConst::GetSecondNeighConstFree");

	tot_sec = 0;
	nsec_type_cum[0] = 0;
	i = 0;
	for (std::list<std::list<string>>::iterator it2 = pool.begin(); it2 != pool.end() && retval; it2++, i++)//going through all the snc constraints
	{
		//Initialization and allocating memories
		nsec_type[i] = std::count_if(it2->begin(), it2->end(), [](string c) { return c.find(CheckKey("STYPE_FROM_TO")) != string::npos; });
		ResizeArray(&tot_sec, tot_sec + nsec_type[i], &sneighbours, "sneighbours", "SecondNeighConst::GetSecondNeighConstFree");
		tot_sec -= nsec_type[i];
		ResizeArray(&tot_sec, tot_sec + nsec_type[i], &dmin2, "dmin2", "SecondNeighConst::GetSecondNeighConstFree");
		tot_sec -= nsec_type[i];
		ResizeArray(&tot_sec, tot_sec + nsec_type[i], &dmax2, "dmax2", "SecondNeighConst::GetSecondNeighConstFree");
		tot_sec -= nsec_type[i];
		ResizeArray(&tot_sec, tot_sec + nsec_type[i], &min2sq, "min2sq", "SecondNeighConst::GetSecondNeighConstFree");
		tot_sec -= nsec_type[i];
		ResizeArray(&tot_sec, tot_sec + nsec_type[i], &max2sq, "max2sq", "SecondNeighConst::GetSecondNeighConstFree");

		int si = 0, pt = 0, nsec = 0;
		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)//going through the lines of the current constraint
		{

			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("CTYPE_FTYPE_FROM_TO"))
			{
				if (pt++ > 0)
				{
					cout << "\nERROR: Multiple " << CheckKey("CTYPE_FTYPE_FROM_TO") << " definition in the " << i + 1 << ". " << CheckTag("SNC") << " constraint!" << endl;
					return false;
				}
				string value;
				if (!GetInt(ss, *it, central+i,"SecondNeighConst::GetSecondNeighConstFree")) return ErrLackValue(*it);
				if (!GetInt(ss, *it, fneighbour + i, "SecondNeighConst::GetSecondNeighConstFree")) return ErrLackValue(*it);
				if (!(ss >> dmin1[i] >> dmax1[i])) return ErrLackValue(*it);
				IntToStr(&conv_numb, i + 1);
				mystrcpy(name, NAME_SIZE, "central[");
				mystrcat(name, NAME_SIZE, conv_numb);
				mystrcat(name, NAME_SIZE, "]");
				RunParams::CheckType(central[i], name, "SecondNeighConst::GetSecondNeighConstFree");
				central[i]--;
				ncentral[i] = SimpleCfg::pnatoms[central[i]];//it will be reset in case of no_periodic
				//Cumulative number of central atoms before this constraint
				ncentral_cum[i] = ncentral_tot;
				//Total number of central atoms
				ncentral_tot += ncentral[i];

				mystrcpy(name, NAME_SIZE, "fneighbour[");
				mystrcat(name, NAME_SIZE, conv_numb);
				mystrcat(name, NAME_SIZE, "]");
				RunParams::CheckType(fneighbour[i], name, "SecondNeighConst::GetSecondNeighConstFree");
				fneighbour[i]--;
				//The concept is that the neighbourlist is calculated for atoms of each type represented in a cosine distribution, second or common neighbour
				//constraint. It is required for the neighbours as well, as the list is updated, and not completly recalculated at each step.
				//As each type can be involved in more than one constraint, all the neighbourlist will be calculated with the
				//same maximum distance, maxdistsq. 
				//Be aware that the central and the neighbour type in secconsttype types represented by right to left!
				//type 0 is represented by binary 1, type 1 binary 10 and so on  
				//consttype pointer is to NeighbourList::consttype, contains all the types involved in Cos, CONC or SNC
				//secconstype contains the types involved only in SNC
				if (!(*consttype >> central[i] & 1))//this type is not yet represented
				{
					*consttype += (longint)pow(2.0, (double)central[i]);
					NeighbourList::nconsttype++;//increase the number of constrained types
				}
				if (!(*consttype >> fneighbour[i] & 1))//this type is not yet represented
				{
					*consttype += (longint)pow(2.0, (double)fneighbour[i]);
					NeighbourList::nconsttype++;//increase the number of constrained types
				}
				if (!(secconsttype >> central[i] & 1))//this type is not yet represented
					secconsttype += (longint)pow(2.0, (double)central[i]);
				
				if (!(secconsttype >> fneighbour[i] & 1))//this type is not yet represented
					secconsttype += (longint)pow(2.0, (double)fneighbour[i]);
				
				if (dmin1[i] > dmax1[i])
				{
					cout << "\nWARNING(" << ++warn << "): the minimum distance value should be less than the minimum distance in the following line, they will be exchanged!" << endl << *it << endl;
					double temp = dmin1[i];
					dmin1[i] = dmax1[i];
					dmax1[i] = temp;
				}

				//calculating the reduced squared distances
				min1sq[i] = pow(dmin1[i] / RunParams::boxedge, 2.);
				max1sq[i] = pow(dmax1[i] / RunParams::boxedge, 2.);
				//finding the largest reduced maximum distance square
				if (max1sq[i] > maxdistsq) maxdistsq = max1sq[i];
			}
			else if (key == CheckKey("STYPE_FROM_TO"))
			{

				index = nsec_type_cum[i] + nsec;
				if (!GetInt(ss, *it, sneighbours + index, "SecondNeighConst::GetSecondNeighConstFree")) return ErrLackValue(*it);
				if (!(ss >> dmin2[index] >> dmax2[index])) return ErrLackValue(*it);
				IntToStr(&conv_numb2, index + 1);
				mystrcpy(name, NAME_SIZE, "sneighbours[");
				mystrcat(name, NAME_SIZE, conv_numb);
				mystrcat(name, NAME_SIZE, ",");
				mystrcat(name, NAME_SIZE, conv_numb2);
				mystrcat(name, NAME_SIZE, "]");
				RunParams::CheckType(sneighbours[index], name, "SecondNeighConst::GetSecondNeighConstFree");
				sneighbours[index]--;
				if (!(*consttype >> sneighbours[index] & 1))//this type is not yet represented
				{
					*consttype += (longint)pow(2.0, (double)sneighbours[index]);
					NeighbourList::nconsttype++;//increase the number of constrained types
				}
				if (!(secconsttype >> sneighbours[index] & 1))//this type is not yet represented
					secconsttype += (longint)pow(2.0, (double)sneighbours[index]);
				
				if (dmin2[index] > dmax2[index])
				{
					cout << "\nWARNING(" << ++warn << "): the 3. (minimum distance) value should be less than the 4. (maximum distance) one in the following line, they will be exchanged!" << endl << *it << endl;
					double temp = dmin1[index];
					dmin1[index] = dmax1[index];
					dmax1[index] = temp;
				}
				
				for (int jj = 0; jj < index - nsec_type_cum[i]; jj++)//check, whether this type was not already given as secondary
				{
					if (sneighbours[index] == sneighbours[nsec_type_cum[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmin2[index] != dmin2[nsec_type_cum[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". second neighbour constraint's " << jj + 1 << ". and " << index - nsec_type_cum[i] + 1 << ". second neighbour types are the same!" << endl;
							cout << "The number of second neighbour types cannot be reduced, as different minimum first-second neighbour distances were given!" << endl;
							cout << "The distances for the " << jj + 1 << ". first and second neighbour is " << dmin2[nsec_type_cum[i] + jj] << " A, while for " << index - nsec_type_cum[i] + 1 << ". first and second neighbour is " << dmin2[index] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different second neighbour atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
					}
				}
				
				for (int jj = 0; jj < index - nsec_type_cum[i]; jj++)//check, whether this type was not already given as secondary
				{
					if (sneighbours[index] == sneighbours[nsec_type_cum[i] + jj])
					{
						//Check, whether the distances for the duplicate secondary type are the same
						if (dmax2[index] != dmax2[nsec_type_cum[i] + jj])
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "The atom type for the " << i + 1 << ". second neighbour constraint's " << jj + 1 << ". and " << index - nsec_type_cum[i] + 1 << ". second neighbour types are the same!" << endl;
							cout << "The number of second neighbour types cannot be reduced, as different maximum first-second neighbour distances were given!" << endl;
							cout << "The distances for the " << jj + 1 << ". first and second neighbour is " << dmax2[nsec_type_cum[i] + jj] << " A, while for " << index - nsec_type_cum[i] + 1 << ". first and second neighbour is " << dmax2[index] << " A." << endl;
							cout << "Cannot decide, which one to use, correct the " << datfilename << " file using only different second neighbour atom types for a constraint!" << endl;
							cout << "Exiting..." << endl;
							CleanExit();
						}
						//only the second neighbour type is different, the distance parameters are the same, can continue by reducing the number of second neighbour types 
						cout << "\nWARNING(" << ++warn << "): The atom type for the " << i + 1 << ". second neighbour constraint's " << jj + 1 << ". and " << index - nsec_type_cum[i] + 1 << ". second neighbour type is the same!" << endl;
						cout << "\tThe number of second neighbor types will be reduced!" << endl;
						nsec_type[i]--;
					}
				}
								
				//Seems to be a valid secondary neighbour
				min2sq[index] = pow(dmin2[index] / RunParams::boxedge, 2.);
				max2sq[index] = pow(dmax2[index] / RunParams::boxedge, 2.);
				if (max2sq[index] > maxdistsq)//finding the largest reduced maximum distance square
					maxdistsq = max2sq[index];
				nsec++;
			}
			else if ((key == CheckKey("COORDNUM_FRACT_SIGMA")) || (key == CheckKey("COORDNUM_FRACT_SIGMA-MASTER")) || (key == CheckKey("COORDNUM_FRACT_SIGMA-SCALABLE")))
			{
				if (si++ > 0)
				{
					cout << "\nWARNING(" << ++warn << "): Multiple ...SIGMA(-MASTER)/...SIGMA-SCALABLE definition in the " << i + 1 << ". " << CheckTag("SNC") << " constraint!" << endl << " The last value will be set as the relevant one." << endl;
					if (RunParams::lead_series_ind == RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + CoordNumbConst::tot_subconst + RunParams::navcoord + RunParams::ncommonneigh + i + 1)
						RunParams::lead_series_ind = -1;
				}

				if (!GetInt(ss, *it, target_coord + i, "SecondNeighConst::GetSecondNeighConstFree")) return ErrLackValue(*it);
				if (!(ss >> fraction[i] >> weights[i])) return ErrLackValue(*it);
				if (key == CheckKey("COORDNUM_FRACT_SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1)
						RunParams::lead_series_ind = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + CoordNumbConst::tot_subconst + RunParams::navcoord + RunParams::ncommonneigh+i + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): More than one **(..)SIGMA-MASTER** entry is declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS") << ", " << CheckTag("CONC");
						cout << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
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
				CheckKeyTag(key, CheckTag("SNC"));//write the appropriate error mesage
				return false;
			}
			nsec_type_cum[i + 1] = nsec_type_cum[i] + nsec_type[i];
		}
		if (retval)
		{
			//Determining the index of the partial calculated from the central first neighbour types 
			cp_partial[i] = (central[i] <= fneighbour[i] ? (central[i] * RunParams::ntypes - (central[i] * (central[i] + 1) / 2) + fneighbour[i]) : \
				(fneighbour[i] * RunParams::ntypes - (fneighbour[i] * (fneighbour[i] + 1) / 2) + central[i]));

			mod_const[i] = 1;//to calculate the initial distribution

		}
	}

	if (conv_numb != NULL)
		delete[] conv_numb;
	delete[] name;
	return retval;

};
//calculate the first and second neighbours
void SecondNeighConst::CalcNeigh()
{
	//the pointers and auxiliary variables are used to extract array values only ones to make calculation hopefully quicker
	int iconst, isec, i, inb, inb2;
	int imin, imax;//cycle boundaries
	int* pnneigh, * pnneigh2;//pointer to the number of neighbours 
	int* plist, * pind, * pind2;//pointers to the indices of the neighbours
	int* ptypelist, * pneigh_type, * pneigh_type2;//pointers to the type of the neighbours 
	int offset;

	double dminsq, dmaxsq;//minimum and maximum squared distances for the neighbours
	double* pdistsq, * pdsq, * pdsq2;//pointers to the distances of the neighbours
	bool quit_search;

	//First find the central-firts neighbour pairs. If a pair found, go through the neighbours of the first neigbours and find the appropriate second neighbours.
	//Calculate the number of centrals, which satisfy the consttraints. 
	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		//setting auxiliary variables and pointers to increase the speed
		dminsq = min1sq[iconst];
		dmaxsq = max1sq[iconst];
		
		pnneigh = neighlist.nneigh_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in nneigh array
		plist = neighlist.neighlist_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in neigh_list array
		ptypelist = neighlist.neightype_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in neigh_type array 
		pdistsq = neighlist.dsq_finder[central[iconst]];//sets the pointer to the beginning of the values belonging to the central atom type in dsq array
		imin = SimpleCfg::cumul[central[iconst]];
		imax = SimpleCfg::cumul[central[iconst] + 1];
		for (i = imin; i < imax; i++)//central atom cycle
		{
			quit_search = false;//serach for neighbours
			pneigh_type = ptypelist;//type of the first neighbour
			pind = plist;//index of the neighbour
			pdsq = pdistsq;//squared distance of the neighbour
			for (inb = 0; inb < *pnneigh; inb++)//neighbour cycle to find first neighbours
			{
				//determine, whether the neighbour is from the right type, and is in the range
				if (*pneigh_type == fneighbour[iconst] && *pdsq >= dminsq && *pdsq <= dmaxsq)
				{
					
					//set the parameters to search for the neighbours of the first neighbour (second neighbours of central)
					offset = *pind - SimpleCfg::cumul[fneighbour[iconst]];
					pnneigh2 = neighlist.nneigh_finder[fneighbour[iconst]] + offset;//number of neighbours of first neighbour
					offset *= *max_neigh;
					pneigh_type2 = neighlist.neightype_finder[fneighbour[iconst]] + offset;//type of the second neighbour
					pind2 = neighlist.neighlist_finder[fneighbour[iconst]] + offset;//index of second neighbour
					pdsq2 = neighlist.dsq_finder[fneighbour[iconst]] + offset;//squared distanceof of second neighbour
					for (inb2 = 0; inb2 < *pnneigh2; inb2++)//going through the neighbourlist of the first neighbour
					{
						if (*pind2 != i)//this is not central 
						{
							for (isec = 0; isec < nsec_type[iconst]; isec++)
							{
								if (*pneigh_type2 == sneighbours[nsec_type_cum[iconst] + isec] && \
									*pdsq2 >= min2sq[nsec_type_cum[iconst] + isec] && *pdsq2 <= max2sq[nsec_type_cum[iconst] + isec])//this is an appropriate second neighbour
								{
									//the atom index should be in its own type
									if (i - SimpleCfg::cumul[central[iconst]] >= ncentral_tot)
										cout << i- SimpleCfg::cumul[central[iconst]] << endl;
									CoordinationNumb(iconst,i- SimpleCfg::cumul[central[iconst]])++;//increase the coordination number for the central atom
									quit_search = true;//if a second neighbour found from the given types, do not have to look for another one
									break;
								}
							}

						}
						pneigh_type2++;
						pind2++;
						pdsq2++;
						if (quit_search)
						{
							quit_search = false;
							break;
						}
					}
				}
				pneigh_type++;
				pind++;
				pdsq++;
			}//end of neighbour cycle inb1
			pnneigh++;
			plist += *max_neigh;
			ptypelist += *max_neigh;
			pdistsq += *max_neigh;
			if (CoordinationNumb(iconst, i - SimpleCfg::cumul[central[iconst]]) == target_coord[iconst])
				nsatisfy[iconst]++;
		}//end of central atom cycle i for primary1
		
	}//end of iconst cycle
};

void SecondNeighConst::UpdateNeigh(Move& move, int sign_switch)
{
	//the pointers and auxiliary variables are used to extract array values only ones to make calculation hopefully quicker
	int iconst, imoved, isec, i, inb1, inb2;
	int involved_type;//ype of the central atom
	int *pnneigh, * pnneigh1;//pointer to the number of neighbours 
	int *pind1, * pind2;//pointers to the indices of the neighbours
	int *pneigh_type1, * pneigh_type2;//pointers to the type of the neighbours 
	int offset;
	int* moved_type, *moved_ind;

	double* pdsq1, * pdsq2;
	double dminsq, dmaxsq;//minimum and maximum squared distances for the neighbours
	
	bool quit_search;

	moved_type = move.types;
	moved_ind = move.indices;
	if (sign_switch == -1)
		involved_count = 0;
	for (imoved = 0; imoved < Move::tot_moved_atoms; imoved++)
	{
		if (!(secconsttype >> *moved_type & 1))
		{
			moved_type++;
			moved_ind++;
			continue;//this moved atom is not involved in this constraint, go to next
		}


		
		if (sign_switch == -1)
			AddToArray(*moved_ind, &involved_count, &involved_indices, "involved_indices", "SecondNeighConst::UpdateNeigh");
				
		//first collect the indices of the moved atoms and their first and second neighbours, as the constraints has to be recalclated for them 
		//all the neighbours on the neighbour list involved in any Second Neigh constraint will be collected, regardless the ranges of the constraits
		offset = *moved_ind - config.cumul[*moved_type];//the offset of the central atom in its own type
		pnneigh = neighlist.nneigh_finder[*moved_type] + offset;//sets the pointer to the nneigh of the central atom

		//setting the pointers to find the first neighbour of the moved atom
		offset *= *max_neigh;//multiplying the offset with the maximum number of neighbours
		pind1 = neighlist.neighlist_finder[*moved_type] + offset;//sets the pointer to the beginning of the neigh_list of the central atom, possible first neighbou
		pneigh_type1 = neighlist.neightype_finder[*moved_type] + offset;//type of the first neighbour
		for (inb1 = 0; inb1 < *pnneigh; inb1++)//neighbour cycle to find  first neighbours
		{
			if (!(secconsttype >> *pneigh_type1 & 1))//not involved in constraints
			{
				pind1++;
				pneigh_type1++;
				continue;
			}
			//add the first neighbour
			AddToArray(*pind1, &involved_count, &involved_indices, "involved_indices", "SecondNeighConst::UpdateNeigh");
			offset = *pind1 - SimpleCfg::cumul[*pneigh_type1];
			pnneigh1 = neighlist.nneigh_finder[*pneigh_type1] + offset;//number of neighbours for the first neighbour
			offset *= *max_neigh;
			pind2 = neighlist.neighlist_finder[*pneigh_type1] + offset;//index of the neighbours for first neighbour
			pneigh_type2 = neighlist.neightype_finder[*pneigh_type1] + offset;//type of the first neighbour

			for (inb2 = 0; inb2 < *pnneigh1; inb2++)//going through the neighbourlist of first neighbour
			{
				if (*pind2 != *moved_ind)//this is not the moved
				{
					if (!(secconsttype >> *pneigh_type2 & 1))//not involved in constraints
					{
						pind2++;
						pneigh_type2++;
						continue;
					}
					AddToArray(*pind2, &involved_count, &involved_indices, "involved_indices", "SecondNeighConst::UpdateNeigh");
				}
				pneigh_type2++;
				pind2++;
			}//end of inb2 cyle
		
			pind1++;
			pneigh_type1++;
		}//end of neighbour cycle inb1
		moved_type++;
		moved_ind++;
	}//end of the imoved cycle

	if (sign_switch == 1)
	{
		
		if (involved_count == 0)
			return;//the moved atoms are not involved in any constarint, no need to update

		//Update the number of centrals, which satisfy the consttraints. 
		for (iconst = 0; iconst < nconstraints; iconst++)
		{
			//setting auxiliary variables and pointers to increase the speed
			dminsq = min1sq[iconst];
			dmaxsq = max1sq[iconst];

			for (i = 0; i < involved_count; i++)//central atom cycle
			{
				involved_type = 0;
				while (involved_indices[i] >= config.cumul[(involved_type)+1])
					involved_type++;
				if (involved_type == central[iconst])//this can be a central in this constraint
				{
					if (CoordinationNumb(iconst, involved_indices[i] - SimpleCfg::cumul[central[iconst]]) == target_coord[iconst])
						nsatisfy[iconst]--;
					CoordinationNumb(iconst, involved_indices[i] - SimpleCfg::cumul[central[iconst]]) = 0;
					//involved_type now points to the type of the moved atom
					quit_search = false;//serach for neighbours
					offset = involved_indices[i] - SimpleCfg::cumul[involved_type];
					pnneigh = neighlist.nneigh_finder[involved_type] + offset;//number of neighbours of central
					offset *= *max_neigh;
					pneigh_type1 = neighlist.neightype_finder[involved_type] + offset;//type of first neighbour
					pind1 = neighlist.neighlist_finder[involved_type] + offset;//index of firstneighbour
					pdsq1 = neighlist.dsq_finder[involved_type] + offset;//squared distance of of central and first neighbour

					for (inb1 = 0; inb1 < *pnneigh; inb1++)//neighbour cycle to find first neighbours
					{
						//determine, whether the neighbour is from the right type, and is in the range
						if (*pneigh_type1 == fneighbour[iconst] && *pdsq1 >= dminsq && *pdsq1 <= dmaxsq)
						{

							//set the parameters to search for the neighbours of the first neighbour (second neighbours of central)
							offset = *pind1 - SimpleCfg::cumul[fneighbour[iconst]];
							pnneigh1 = neighlist.nneigh_finder[fneighbour[iconst]] + offset;//number of neighbours of first neighbour
							offset *= *max_neigh;
							pneigh_type2 = neighlist.neightype_finder[fneighbour[iconst]] + offset;//type of the second neighbour
							pind2 = neighlist.neighlist_finder[fneighbour[iconst]] + offset;//index of second neighbour
							pdsq2 = neighlist.dsq_finder[fneighbour[iconst]] + offset;//squared distanceof of second neighbour
							for (inb2 = 0; inb2 < *pnneigh1; inb2++)//going through the neighbourlist of the first neighbour
							{
								if (*pind2 != involved_indices[i])//this is not central 
								{
									for (isec = 0; isec < nsec_type[iconst]; isec++)
									{
										if (*pneigh_type2 == sneighbours[nsec_type_cum[iconst] + isec] && \
											* pdsq2 >= min2sq[nsec_type_cum[iconst] + isec] && *pdsq2 <= max2sq[nsec_type_cum[iconst] + isec])//this is an appropriate second neighbour
										{
											//the atom index should be in its own type
											if (involved_indices[i] - SimpleCfg::cumul[central[iconst]] >= ncentral_tot)
												cout << i - SimpleCfg::cumul[central[iconst]] << endl;
											CoordinationNumb(iconst, involved_indices[i] - SimpleCfg::cumul[central[iconst]])++;//increase the coordination number for the central atom
											
											quit_search = true;//if a second neighbour found from the given types, do not have to look for another one
											break;
										}
									}

								}
								pneigh_type2++;
								pind2++;
								pdsq2++;
								if (quit_search)
								{
									quit_search = false;
									break;
								}
							}
						}
						pneigh_type1++;
						pind1++;
						pdsq1++;
					}//end of neighbour cycle inb1
					
					
					if (CoordinationNumb(iconst, involved_indices[i] - SimpleCfg::cumul[central[iconst]]) == target_coord[iconst])
						nsatisfy[iconst]++;
				}//this can be central
			}//end of central atom cycle i for primary1

		}//end of iconst cycle
		if (involved_indices != NULL)
		{
			delete[] involved_indices;
			involved_indices = NULL;
			involved_count = 0;


		}
	}
}
//save to a file
void SecondNeighConst::Save(const char *file_name )
{
	int i, k, ineightype;
	ofstream file;

	if (strlen(file_name) != 0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, sncfilename);
	OpenFile(file, tempfilename, "SecondNeighConst::Save", 0);//open file, check, whether it was successfully opened


	file << "Data for Second Neighbour Constraints. "\
		<< "file created by SecondNeighConst::Save \n" << endl;

	file << nconstraints << "\t number of constraints (nconstraints) " << endl;
	file << "Constraints properties: " << endl;
	file << "index, central atom type, first neighbour atom type, inferior and superior central-first neigh. distance," << endl;
	file << "number of second neighbour types, second neighbour type(s), inferior and superior first and second neigh. distance," << endl;
	file << "desired coordination number, desired fraction of central atoms,  weighting parameter, number of central atoms, " << endl;
	file << "number of cenral atoms satisfying the constraint " << endl;

	for (i = 0; i < nconstraints; i++)
	{
		file << J4 << i + 1 << "\t";
		file << J4 << central[i] + 1 << "\t";//start with 1 in the file
		file << J4 << fneighbour[i] + 1 << "\t";//start with 1 in the file
		file.precision(12);
		file.setf(ios::fixed, ios::floatfield);
		file.setf(ios::right, ios::adjustfield);
		file << J10 << dmin1[i] << "\t";//minimum distance between central and first neighbour
		file << J10 << dmax1[i] << "\t";//maximum distance between central and first neighbour
		file << J4 << nsec_type[i] << "\t";//number of second neighbour types
		for (ineightype = 0; ineightype < nsec_type[i]; ineightype++)
			file << J4 << sneighbours[nsec_type_cum[i] + ineightype] + 1 << "\t";//start with 1 in the file
		for (ineightype = 0; ineightype < nsec_type[i]; ineightype++)
			file << J10 << dmin2[nsec_type_cum[i] + ineightype] << "\t";
		for (ineightype = 0; ineightype < nsec_type[i]; ineightype++)
			file << J10 << dmax2[nsec_type_cum[i] + ineightype] << "\t";
		file.precision(6);
		file.unsetf(ios::fixed);
		file.unsetf(ios::right);
		file << J4  << target_coord[i] << "\t";
		file << J10 << fraction[i] << "\t";
		file << J10 << weights[i] << "\t";
		file << J10 << ncentral[i] << endl;
	}

	file << endl;
	file << "Number of atoms satisfying each constraint (nsatisfy)" << endl;
	file << "Coordination numbers for central atoms per constraint," << endl;

	for (i = 0; i < nconstraints; i++)//for each constraint
	{
		file << i + 1 << "-th constraint";

		//Number of atoms in the configuration satisfying the constraint
		file << "\n" << J10 << "nsatisfy/constraint\t";
		file << J7 << nsatisfy[i] << "\t";
		file << "\n" << endl;

		for (k = 0; k < ncentral[i]; k++)//for every central atom in this constraint
		{
			file << J10 << k + 1 << "\t";
			file << J7 << CoordinationNumb(i,k);
			file << endl;
		}
		file << endl;
	}
	file.close();
}

//----Copy the coordination numbers and the number of atoms satisfying the ---------
//----constraints affected by the move from source to target---------
void SecondNeighConst::CopyModified(SecondNeighConst &target)
{
	int iconst, i;
	int *p_source, *p_target;

	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		if (mod_const[iconst])
		{
			//the coordnumb array of this constraint was modified, copy it 
			p_source = coordnumbs + ncentral_cum[iconst];
			p_target = target.coordnumbs + ncentral_cum[iconst];
			for (i = 0; i < ncentral[iconst]; i++)
				*p_target++ = *p_source++;
			
			target.nsatisfy[iconst] = nsatisfy[iconst];
		}

	}//end of constraint cycle iconst
};







#endif // _ADVANCED_GEOM_CONST