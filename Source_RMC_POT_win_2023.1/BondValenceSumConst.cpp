//source BondValenceSumConst.cpp
//Last changed 24.01.2023

//(1)	Each constraint can have more, than one neighbour type with their own rmax, R_0 and b parameters
//(2)	If a central atom has more than one valence state and therefore desired valence, than each has to be given as a separate constraint
//		with the target fraction and target valence. Alternatively the difference valence states can have their own RMC type.


#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//classes1.h is included through this

#ifdef _ADVANCED_GEOM_CONST
int  BondValenceSumConst::nconstraints;//the number of constraints
int  BondValenceSumConst::nthreads;//the total number of threads to use, set by RunParams::GetParams
int  BondValenceSumConst::tot_neightype;//total number of neighbour types for all the constraints

longint *BondValenceSumConst::bvsctype;//array containing the partial indices in binary format

//parameters for each constraint

int *BondValenceSumConst::central;//array of types of central particles [nconstraints]
int *BondValenceSumConst::neighbours;//array of types of neighbour particles  [nconstraints]
int *BondValenceSumConst::n_neightype;//number of bvs neighbours types for a costraint
int *BondValenceSumConst::ncentral;//number of central atoms for each constarint
int BondValenceSumConst::ncentral_tot;//total number of central atoms for all the constarint
int *BondValenceSumConst::ncentral_cum;//cumulative number of central atoms before the given constraint for each constraint
int *BondValenceSumConst::cum_n_neightype;//cumulative number of neighbour types

int *BondValenceSumConst::central_charge;//assumed oxidation state of the central atom  [nconstraints]
int *BondValenceSumConst::neigh_charge;//assumed oxidation state of the neighbour atom [tot_bvs_neigh] 

 
double *BondValenceSumConst::dmax;//array of maximum distances between central and neighbour [tot_bvs_neigh]
double *BondValenceSumConst::R0;//ideal bond length for each neighbour of each constraint [tot_bvs_neigh]
double *BondValenceSumConst::b;//b parameter for each neighbour of each constraint [tot_bvs_neigh]
double *BondValenceSumConst::R0_in;//ideal bond length for each neighbour of each constraint [tot_bvs_neigh] to keep it for *.free 
double *BondValenceSumConst::b_in;//b parameter for each neighbour of each constraint [tot_bvs_neigh] to keep it for *.free 

double *BondValenceSumConst::target_valence;//array of DESIRED valence  [nconstraints]
double *BondValenceSumConst::weights;//array of weights for each subconstraint  [nconstraints]

int BondValenceSumConst::is_set_def=0;//whether to set default R0 or b for any constraint
int *BondValenceSumConst::set_def;//whether to set default for a subconstraint[tot_bvs_neigh]
int *BondValenceSumConst::mod_const;//whteher the constraint was effected by the move
double *BondValenceSumConst::maxsq;//reduced min. dist. squared
double BondValenceSumConst::maxdistsq;//the maximum distance in the maxsq, this is used determining the limit of the neighbourlist calculation
string *BondValenceSumConst::pair_name;//the names with valence states for the atom pairs, if it is given [tot_bvs_neigh] 
									  //----default constructor----------------
BondValenceSumConst::BondValenceSumConst()
{
	int i,offset;
	if (nconstraints < 1)
	{
		if (::debug)
		{
			cout << "\nWARNING(" << ++warn << "): BondValenceSumConst constructor" << endl;
			cout << "\tThere are no bond valence sum constraints, or nconstraints was not initialised!" << endl;
		}
	}
	//allocating memory
	//CACHE padding should be applied
	offset = (nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*valence)) / sizeof(*valence)) : 0);//number of dummy double elements
	SetArraysize(&valence, ncentral_tot * nthreads + (nthreads - 1) * offset, "valence", "BondValenceSumConst::BondValenceSumConst");

	SetArraysize(&finder, nconstraints, "finder", "BondValenceSumConst::BondValenceSumConst");
	SetArraysize(&thread_finder, nthreads, "thread_finder", "BondValenceSumConst::BondValenceSumConst");//for the beginning of each thread's segment

	InitValence();//Sets the array elements to zero

	//initialising the finder
	thread_finder[nthreads - 1] = valence;//the main has the first segment
	for (i = 0; i < nthreads - 1; i++)
		thread_finder[i] = valence + (i + 1) * (ncentral_tot+offset);

	//initialising the finder
	for (i = 0; i < nconstraints; i++)
		finder[i] = valence + ncentral_cum[i];
};

//copy constructor
BondValenceSumConst::BondValenceSumConst(BondValenceSumConst &source)
{
	int i,offset;

	//allocating memory
	//CACHE padding should be applied
	offset = (nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*valence)) / sizeof(*valence)) : 0);//number of dummy double elements
	SetArraysize(&valence, ncentral_tot * nthreads+(nthreads-1)*offset, "valence", "BondValenceSumConst::BondValenceSumConst");

	SetArraysize(&finder, nconstraints, "finder", "BondValenceSumConst::BondValenceSumConst");
	SetArraysize(&thread_finder, nthreads, "thread_finder", "BondValenceSumConst::BondValenceSumConst");//for the beginning of each thread's segment
	
	for (i = 0; i < ncentral_tot; i++)
		valence[i] = source.valence[i];

	//initialising the finder
	thread_finder[nthreads - 1] = valence;//the main has the first segment
	for (i = 0; i < nthreads - 1; i++)
		thread_finder[i] = valence + (i + 1) * (ncentral_tot+offset); 

	//initialising the finder for the constraints
	for (i = 0; i < nconstraints; i++)
		finder[i] = valence + ncentral_cum[i];
	
};
void BondValenceSumConst::InitValence()
{
	int i, size1, offset;
	size1 = 0;

	offset = (nconstraints > 0 ? (int)((CACHE_PADDING - sizeof(*valence)) / sizeof(*valence)) : 0);//number of dummy double elements
	//Total number of central atoms
	for (i = 0; i < nconstraints; i++)
		size1 += ncentral[i];
	size1 *= nthreads;
	size1+=(nthreads - 1) * offset;//each thread has its own segment

	//setting all the array elements to 'zero' values
	for (i = 0; i < size1; i++)
		valence[i] = 0.0;


};
//----------------Getting the parameters from the.dat file, and sets the static members-------------------
bool BondValenceSumConst::GetBondValenceSumConstFree(std::list<std::list<string> > pool)
//Process free format strings related to setup. Returns with false, if incoming datasets are insufficient.
{
	int i;
	unsigned int npartials = RunParams::ntypes * (RunParams::ntypes + 1) / 2;
	nconstraints =  RunParams::nbvs;
	char *name;
	if (npartials > sizeof(*bvsctype) * 8)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The program can handle only " << sizeof(*bvsctype) * 8 << " partials presently!" << endl;
		cout << "Change the type of BondValenceSumConst::cctype if possible to have more bytes!" << endl;
		return false;
	}
	SetArraysize(&name, NAME_SIZE, "name", "BondValenceSumConst::GetBondValenceSumConstFree");
	//allocating memory
	SetArraysize(&central, nconstraints, "central", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&central_charge, nconstraints, "central_charge", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&target_valence, nconstraints, "target_valence", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&weights, nconstraints, "weights", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&n_neightype, nconstraints, "n_neightype", "BondValenceSumConst::GetBondValenceSumConstFree");
	
	SetArraysize(&ncentral, nconstraints, "ncentral", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&bvsctype, nconstraints, "bvsctype", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&cum_n_neightype, nconstraints + 1, "cum_n_neightype", "BondValenceSumConst::GetBondValenceSumConstFree");
	SetArraysize(&ncentral_cum, nconstraints + 1, "ncentral_cum", "BondValenceSumConst::GetBondValenceSumConstFree");
		
	

	
	maxdistsq = 0;
	cum_n_neightype[0] = 0;//cumulative number of neighbour type before the given constraint
	tot_neightype = 0;
	ncentral_tot = 0;

	i = 0;
	bool retval = true;
	for (std::list<std::list<string>>::iterator it2 = pool.begin(); it2 != pool.end() && retval; it2++, i++)
	{
		//Initialization and allocating memories
		n_neightype[i] = (int)std::count_if(it2->begin(), it2->end(), [](string c) { return c.find(CheckKey("NEIGH-TYPE_RMAX_R0_B_CHARGE")) != string::npos; });
		cum_n_neightype[i + 1] = cum_n_neightype[i] + n_neightype[i];//cumulative number of neighbour type before the given constraint
		
		//allocating memory-total number of neighbour types
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &neighbours, "neighbours", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &dmax, "dmax", "BondValenceSumConst::GetBondValenceSumConstFree"); 
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &maxsq, "maxsq", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &R0, "R0", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &b, "b", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &R0_in, "R0_in", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &b_in, "b_in", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &neigh_charge, "neigh_charge", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &set_def, "set_def", "BondValenceSumConst::GetBondValenceSumConstFree");
		tot_neightype -= n_neightype[i];
		ResizeArray(&tot_neightype, tot_neightype + n_neightype[i], &pair_name, "pair_name", "BondValenceSumConst::GetBondValenceSumConstFree");
		for (int j = cum_n_neightype[i]; j < cum_n_neightype[i + 1]; j++)
		{
			set_def[j] = 0;//do not set the default value
			pair_name[j].append("");
		}

		int ce = 0, n = 0, co = 0;
		for (std::list<string>::iterator it = it2->begin(); (it != it2->end()) && retval; it++)
		{

			string key, params;
			WrapFree(*it, key, params);

			stringstream ss(params);

			if (key == CheckKey("CENT-TYPE_CHARGE"))
			{
				if (ce++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("CENT-TYPE_CHARGE") << " definition in the " << i + 1 << ". " << CheckTag("BVS") << " constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (!GetInt(ss, *it, central+i, "BondValenceSumConst::GetBondValenceSumConstFree")) return ErrLackValue(*it);
				if ((central[i] > RunParams::ntypes) || (central[i] < 1))
				{
					cout << "\nERROR: The " << CheckKey("CENT-TYPE_CHARGE") << " value" << " (" << central[i] << ") should be in the range: 1.." << RunParams::ntypes << " (number of types)!" << endl << "Please, correct it at the " << i + 1 << ". " << CheckTag("COORD") << " constraint!" << endl;
					return false;
				}
				if (!(ss >> central_charge[i])) central_charge[i]=NO_IVALUE_DEF;//no charge was given, set_def is not set here
				central[i]--;
				ncentral[i] = SimpleCfg::pnatoms[central[i]];//it will be reset in case of no_periodic
				//Cumulative number of central atoms before this constraint
				ncentral_cum[i] = ncentral_tot;
				//Total number of central atoms
				ncentral_tot += ncentral[i];
			}
			else if (key == CheckKey("NEIGH-TYPE_RMAX_R0_B_CHARGE"))
			{
				int aind = cum_n_neightype[i] + n++;//Actual index
				if (!GetInt(ss, *it, neighbours + aind, "BondValenceSumConst::GetBondValenceSumConstFree")) return ErrLackValue(*it);
				if ((neighbours[aind] > RunParams::ntypes) || (neighbours[aind] < 1))
				{
					cout << "\nERROR: The NEIGH-TYPE_RMAX_R0_B_CHARGE) value (" << neighbours[aind] << ") should be in the range: 1.." << RunParams::ntypes << " (number of types)!" << endl << "Please, correct it at the " << i + 1 << ". " << CheckTag("BVS") << " constraint!" << endl;
					return false;
				}
				neighbours[aind]--;
				for (int jj = 0; jj < n-1; jj++)//check, whether this type was not already given as neighbour
				{
					if (neighbours[aind] == neighbours[cum_n_neightype[i] + jj])
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "The atom type for the " << i + 1 << ". bond valence sum constraint's " << jj + 1 << ". and " << aind - cum_n_neightype[i] + 1 << ". neighbour types are the same!" << endl;
						cout << "Cannot give the same neighbour type twice, correct the " << datfilename << " file using only different neighbour atom types for a constraint!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
				}
				if (!(ss >> dmax[aind])) return ErrLackValue(*it);
				if (!(ss >> R0[aind])) return ErrLackValue(*it);
				R0_in[aind] = R0[aind];
				if (R0[aind]< 0)
				{
					is_set_def = 1;//default should be set for some value
					set_def[aind] += 1;//adding 1 means set default R0
				}
				if (!(ss >> b[aind]))
				{
					if (set_def[aind] > 0)//R0 is not given, default should be extracted from the table
						return ErrLackValue(*it);
					else//R0 is given, simple b deafult is set, no need here for table lookup, charge not needed
						b[aind] = BVS_B_DEF;
				}
				else//something was read for b
				{
					if (b[aind] < 0) 
					{
						if (set_def[aind] == 0)//R0 is given
							b[aind] = BVS_B_DEF;//no need for default value looktable
						else
							set_def[aind] += 2;//adding 2 means set default b form table
					}
				}
				b_in[aind] = b[aind];
							
				if (set_def[aind] > 0)//table lookup needed, read charge
				{
					if (!(ss >> neigh_charge[aind])) return ErrLackValue(*it);
				}
				else
					neigh_charge[aind] = NO_IVALUE_DEF;
	
				//Initialising the remainder of the static members
				if (RunParams::boxedge > 0)
				{
					//Minimum and maximum distances for the constraints squared in reduced units
					maxsq[aind] = pow(dmax[aind], 2) / RunParams::boxedge / RunParams::boxedge;
					if (maxsq[aind] > maxdistsq)
						maxdistsq = maxsq[aind];
				}
			}
			else if ((key == CheckKey("VALENCE_SIGMA")) || (key == CheckKey("VALENCE_SIGMA-MASTER")) || (key == CheckKey("VALENCE_SIGMA-SCALABLE")))
			{
				if (co++ > 0) cout << "\nWARNING(" << ++warn << "): Multiple " << CheckKey("VALENCE_SIGMA...") << " definition in the " << i + 1 << ". " << CheckTag("BVS") << " constraint!" << endl << " The last value will be set as the relevant one." << endl;
				if (!(ss >> target_valence[i])) return ErrLackValue(*it);
						if (!(ss >> weights[i])) return ErrLackValue(*it);
								
				if (key == CheckKey("VALENCE_SIGMA-MASTER"))
				{
					if (RunParams::lead_series_ind == -1)
						RunParams::lead_series_ind = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr + CoordNumbConst::tot_subconst + RunParams::navcoord + RunParams::ncommonneigh + RunParams::nsecondneigh+ i + 1;
					else
					{
						cout << "\nWARNING(" << ++warn << "): Multiple **(..)SIGMA-MASTER** entry is declared in " << CheckTag("EXP") << ", " << CheckTag("COORD") << ", " << CheckTag("AVCOORD") << ", " << CheckTag("COS");
						cout << ", " << CheckTag("CONC") << ", " << CheckTag("SNC") << ", " << CheckTag("BVS");
#ifdef LOCAL_INV
						cout << ", " << CheckTag("LOCINV");
#endif
						cout << " section!" << endl;
						RunParams::lead_series_ind = -2;//not set correctly
					}
				}
				else if (key == CheckKey("VALENCE_SIGMA-SCALABLE"))
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
				CheckKeyTag(key, CheckTag("BVS"));//write the appropriate error mesage
				return false;
			}
		}
		if (retval)
		{
			//Calculating the indices of the central-neighbours partials for each constraint.
			//the bvcsctype will contain the sum of pow(2,partialind) for each partials involved for the given constraint
			//The partialind has to begin with 1 because of this.
			
			bvsctype[i] = 0;
			for (int ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			{
				if (central[i] <= neighbours[cum_n_neightype[i] + ineightype])
					bvsctype[i] += (long)pow(2.0, (double)(central[i] * RunParams::ntypes - central[i] * (central[i] + 1) / 2 + neighbours[cum_n_neightype[i] + ineightype]));
				else
					bvsctype[i] += (long)pow(2.0, (double)(neighbours[cum_n_neightype[i] + ineightype] * RunParams::ntypes - neighbours[cum_n_neightype[i] + ineightype] * (neighbours[cum_n_neightype[i] + ineightype] + 1) / 2 + central[i]));
				
			}
			//Checking, whether the central charge is given, if default value lookup is needed
			if (central_charge[i] == NO_DVALUE_DEF)
			{
				for (int ineightype = cum_n_neightype[i]; ineightype < cum_n_neightype[i + 1]; ineightype++)
				{
					if (set_def[ineightype] > 0)
					{
						cout << "\nERROR: In the " << i + 1 << ". Bond valence sum constraint no charge was given for the central atom, although" << endl;
						cout << "default value for R0 (and for b) is needed for the " << ineightype - cum_n_neightype[i] + 1 << ". neighbour type!" << endl;
						cout << "Amend the " << datfilename << " file and try again! Cannot run this way, exiting..." << endl;
						CleanExit();
						return false;
					}
				}
			}

		}
	}
	if (is_set_def>0)
		GetDefParams();
	CheckR0();//check, whether the R0 is inside the Rmax range 
	if (name != NULL)
		delete[] name;
	//Cumulative number of central atoms before this constraint
	if (retval) ncentral_cum[nconstraints] = ncentral_tot;
	return retval;
};

//get the default R0 and b from the BV_param table
void BondValenceSumConst::GetDefParams()
{
	int  iconst, ineigh;
	int abs_charge;
	char *name1, *name2;
	string full_name;
	bool found = false;
	if (RunParams::chem_symbols_standard == nullptr)
	{
	

		cout << "\n*****ERROR*****" << endl;
		cout << CheckKey("CHEMICAL-SYMBOLS") << " key word is needed in the [ " << CheckTag("GENERAL") << " ] section, if bond valence sum constraint " << endl;
		cout << " is used!" << endl;
		cout << "It is missing, so add it to " << datfilename << " data file and try again!" << endl;
		CleanExit();
	}
	SetArraysize(&name1, NAME_SIZE, "name1", "BondValenceSumConst::GetDefParams");
	SetArraysize(&name2, NAME_SIZE, "name2", "BondValenceSumConst::GetDefParams");
	//going through the constraints
	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		string cname(RunParams::chem_symbols_standard[central[iconst]]);
		abs_charge=abs(central_charge[iconst]);
		IntToStr(&name1,abs_charge);
		cname.append(name1);
		if (central_charge[iconst] < 0)
			cname.append("-");
		else
			if (central_charge[iconst] > 0)
				cname.append("+");
		for (ineigh = cum_n_neightype[iconst]; ineigh < cum_n_neightype[iconst + 1]; ineigh++)
		{
			found = false;
			if (set_def[ineigh] > 0)//there is a value to look up in the table
			{
				string neiname(RunParams::chem_symbols_standard[neighbours[ineigh]]);
				abs_charge = abs(neigh_charge[ineigh]);
				IntToStr(&name2, abs_charge);
				neiname.append(name2);
				if (neigh_charge[ineigh] < 0)
					neiname.append("-");
				else
					if (neigh_charge[ineigh] > 0)
						neiname.append("+");
				//most of the first compounds in the BVS_param table are cations, and the second anion, concatenate this way first
				bool first_positive=false;
				if (central_charge[iconst] >= 0)
					first_positive = true;
				full_name.clear();
				if (first_positive)
				{
					full_name.append(cname);
					full_name.append(" ");
					full_name.append(neiname);
				}
				else
				{
					full_name.append(neiname);
					full_name.append(" ");
					full_name.append(cname);
				}

				if (BV_param.count(full_name) > 0)//it is in the table
				{
					R0[ineigh] = BV_param.at(full_name).R0;
					if (set_def[ineigh] == 3)//b is needed as well
						b[ineigh] = BV_param.at(full_name).b;
					found = true;
					pair_name[ineigh].append(full_name);
				}
				else//central - neigh not in the table, try neigh - central
				{
					full_name.clear();
					if (first_positive)
					{
						full_name.append(neiname);
						full_name.append(" ");
						full_name.append(cname);
					}
					else
					{
						full_name.append(cname);
						full_name.append(" ");
						full_name.append(neiname);
					}

					if (BV_param.count(full_name) > 0)//it is in the table
					{
						R0[ineigh] = BV_param.at(full_name).R0;
						if (set_def[ineigh] == 3)//b is needed as well
							b[ineigh] = BV_param.at(full_name).b;
						found = true;
						pair_name[ineigh].append(full_name);
					}
					//if no exact oxidation state is found, can search among the nonspecified as well, but not implemented now //GO
					if (!found)
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "There are no default parameters in the Bond valence parameters table for " << endl;
						cout << full_name << " given for the " << iconst + 1 << ". BVS constraint's " << ineigh + 1-cum_n_neightype[iconst] << ". neighbour." << endl;
						cout << "Entries containing the central neighbour pair " << RunParams::chem_symbols_standard[central[iconst]] << " ";
						cout << RunParams::chem_symbols_standard[neighbours[ineigh]] <<" are the following:"<< endl;
						GiveBVSymbols(3, RunParams::chem_symbols_standard[central[iconst]],RunParams::chem_symbols_standard[neighbours[ineigh]]);
						cout << "Either correct the oxidation states or give the R0 and b values in the " << datfilename << " file!" << endl;
						cout << "Cannot run this way, exiting!" << endl;
						CleanExit();
					}
				}
			}
		}
	}
};

//check whether the R0 is inside rmax
void BondValenceSumConst::CheckR0()
{
	int i,j;
	for (i = 0;i<nconstraints  ;i++)
		for (j=cum_n_neightype[i];j< cum_n_neightype[i+1];j++)
			if (R0[j] > dmax[j])
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "The R0 value (" << R0[j] << " A) for the " << i + 1 << ". BVS constraint's " << j- cum_n_neightype[i] + 1 << ". neighbour is outside the cutoff ";
				cout << dmax[j] << " A !" << endl;
				cout<<"Specify adequate values in the "<<datfilename<<" file!" << endl;
				cout << "Cannot run this way, exiting!" << endl;
				CleanExit();
			}

};

//============save all the data to the file given in argument============
void BondValenceSumConst::Save(const char *file_name)
{
	int i, k, ineightype;
	ofstream file;

	if (strlen(file_name) != 0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, bvsfilename);
	OpenFile(file, tempfilename, "BondValenceSumConst::Save", 0);//open file, check, whether it was successfully opened


	file << "Data for Bond Valnce Sum Constraints. "\
		<< "file created by BondValenceSumConst::Save \n" << endl;

	file << nconstraints << "\t number of constraints (nconstraints) " << endl;
	file << "Constraints properties: " << endl;
	file << "index, central atom type, number of neighbour types, neighbour atom type(s), cutoff distance(s), R0(s), b(s)," << endl;
	file << "desired valence, weighting parameter, number of central atoms" << endl;

	for (i = 0; i < nconstraints; i++)
	{
		file << J4 << i + 1 << "\t";
		file << J4 << central[i] + 1 << "\t";//start with 1 in the file
		file << J4 << n_neightype[i] << "\t";
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			file << J4 << neighbours[cum_n_neightype[i] + ineightype] + 1 << "\t";//start with 1 in the file
		file.precision(12);
		file.setf(ios::fixed, ios::floatfield);
		file.setf(ios::right, ios::adjustfield);
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			file << J10 << dmax[cum_n_neightype[i] + ineightype] << "\t";
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			file << J10 << R0[cum_n_neightype[i] + ineightype] << "\t";
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			file << J10 << b[cum_n_neightype[i] + ineightype] << "\t";

		file.precision(6);
		file.unsetf(ios::fixed);
		file.unsetf(ios::right);
		
		file << J4 << target_valence[i] << "\t";
		file << J10 << weights[i] << "\t";
		file << J10 << ncentral[i] << endl;
	}

	file << endl;
	file << "Bond valence sum for central atoms per constraint," << endl;

	for (i = 0; i < nconstraints; i++)//for each constraint
	{
		file << i + 1 << "-th constraint";
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
			if (pair_name[cum_n_neightype[i] + ineightype].compare("") != 0)
				file << ", " << J10 << pair_name[cum_n_neightype[i] + ineightype];
			else
				file << ", unknown";
		file << endl;
	
		for (k = 0; k < ncentral[i]; k++)//for every central atom in this constraint
			file << J10 << k + 1 << "\t" << J7 << BondValence(i, k) << endl;

		file << endl;
	}

	file.close();
}

//--------------------------------------------------------------------------------
//gets or sets the value for the iconst-th constraints, iatom-th 
//central atom in the valence array
//iatom index should be in its own type
double &BondValenceSumConst::BondValence(int iconst, int iatom)
{
	return(*(finder[iconst] + iatom));
}

//----Copy the valence for the constraints affected by the move from source to target---------
void BondValenceSumConst::CopyModified(BondValenceSumConst &target, ThreadArg &thread_arg)
{
	int iconst, i, ineightype;
	int partialind;
	int min_ind, n_cent;
	double *p_source, *p_target;

	for (iconst = 0; iconst < nconstraints; iconst++)
	{
		partialind = 0;
		//searching for the first existing partial of this constraint
		for (ineightype = 0; ineightype < n_neightype[iconst]; ineightype++)
		{
			while (((bvsctype[iconst] >> partialind) & 1) == 0)
				partialind++;

			if (Threads::mod_partial[partialind])
			{
				//bvs_central_offset_min, bvs_central_offset_max contains the segmentation for the central atoms of each constraint
				min_ind = thread_arg.bvs_central_offset_min[iconst];//offset of the first central 
				n_cent = thread_arg.bvs_central_offset_max[iconst] - min_ind + 1;//number of centrals

				//Valences
				p_source = finder[iconst] + min_ind;//sets the pointer to the source's valence array 
												//to the beginning of the iconst-th constraint for this tread
				p_target = target.finder[iconst] + min_ind;//sets the pointer to the target's valence array 
											//to the beginning of the iconst-th constraint for this tread
				for (i = 0; i < n_cent; i++)
					*p_target++ = *p_source++;//copy source valence to target

				break;//at least this partial was modified, the array was copied, no need to search for more modified partial
			}//end of if this constraint was affected by the move
			partialind++;//this partial of the constraint was not modified, goto next
		}//end of neighbour type cycle ineightype 
	}//end of constraint cycle iconst
};

//Updating the valence array with the contribution of the threads, if there is any
int BondValenceSumConst::UpdateValence(ThreadArg &thread_arg)
{
	int iconst, ithread, j;
	int neightype;
	int	modified = 0;
	int partialind;
	int min_ind, n_cent;
	double *pval, *pval_aux;
	longint extract_part;
#ifdef _TEST_MODE
	int index = thread_arg.thread_index;
	std::chrono::duration<double, std::milli> elapsed;
#endif

	//sets the pointer to the beginning of the valence array
	for (iconst = 0; iconst < BondValenceSumConst::nconstraints; iconst++)
	{

		extract_part = BondValenceSumConst::bvsctype[iconst];
		partialind = 0;
		//searching for the first existing partial of this constraint
		for (neightype = 0; neightype < BondValenceSumConst::n_neightype[iconst]; neightype++)
		{

			while (((extract_part >> partialind) & 1) == 0)
				partialind++;

			if (Threads::mod_partial[partialind])
			{
#ifdef _TEST_MODE
				Threads::btime1[index] = std::chrono::high_resolution_clock::now();
#endif
				modified = 1;//at least one constraint was modified by the move
				//the constraint was effected by the move
				//first add the contribution of the thread's valence array to the main

				//bvs_central_offset_min, bvs_central_offset_max contains the segmentation for the central atoms of each constraint
				min_ind = thread_arg.bvs_central_offset_min[iconst];//offset of the first central 
				n_cent = thread_arg.bvs_central_offset_max[iconst] - min_ind + 1;//number of centrals

				for (ithread = 0; ithread < nthreads - 1; ithread++)
				{
					pval = finder[iconst] + min_ind;//sets the pointer to the first atom of this constraint in the valence array fot this thread
					pval_aux = thread_finder[ithread] + BondValenceSumConst::ncentral_cum[iconst] + min_ind;
					for (j = 0; j < n_cent; j++)
					{
						*pval += *pval_aux;
						*pval_aux++ = 0;
						pval++;
					}
				}


#ifdef _TEST_MODE
				Threads::btime2[index] = std::chrono::high_resolution_clock::now();
				elapsed = Threads::btime2[index] - Threads::btime1[index];
				Threads::dur_18[index] += elapsed.count();
#endif
				break;//it is recalculated, no need to search for more modified partials
			}
			partialind++;//this partial of the constraint was not modified, goto next
		}
	}
	return(modified);
};

//loading the bond valence sum const from the *.bbv file
int BondValenceSumConst::LoadBinary(double *sigma_percentage) const
{
	int i, j, temp, *temp_i;
	double *temp_d;
	string warning;
	if (RunParams::potential>0)
		warning.append("\tHistogram, potential and related constraints has to be recalculated!");
	else
		warning.append("\tHistogram and related constraints has to be recalculated!");
	
	ifstream file;

	CleanOpen(file, bbvfilename, 1);//open file
	if (CheckFileState(file, " BondValenceSumConst::LoadBinary", bbvfilename) == 0)//check, whether it was successfully opened
	{
		cout << "\nWARNING(" << ++warn << "): There is no open " << bbvfilename << " file to load the valence from" << endl;
		cout << warning << endl;
		return(0);//loading was not successful
	}


	//reading the number of constraints
	if (ReadBin(file, &temp, 1, "nconstraints", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp != nconstraints)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of bond valence sum constraints is " << temp << endl;
		cout << "\tfrom the " << bbvfilename << " file, and " << nconstraints << " based on the "<<datfilename<<" file!" << endl;
		cout << warning << endl;
		return(0);
	}
	if (nconstraints < 1)
		return(0);
	SetArraysize(&temp_i, tot_neightype, "temp_i", "BondValenceSumConst::LoadPotBinary");
	SetArraysize(&temp_d, tot_neightype, "temp_d", "BondValenceSumConst::LoadPotBinary");

	if (ReadBin(file, temp_i, nconstraints, "central", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		if (temp_i[i] != central[i])
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The central type for the " << i + 1 << ". bond valence sum constraint is " << temp_i[i] << endl;
			cout << "\tfrom the " << bbvfilename << " file, and " << central[i] << " based on the " << datfilename << " file!" << endl;
			cout << warning << endl;
			return(0);
		}
	}
	
	if (ReadBin(file, temp_d, nconstraints, "target_valence", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		if (fabs(temp_d[i] -target_valence[i])>LOAD_TOL)
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The target valence for the "<<i+1<<". bond valence sum constraints is " << temp_d[i] << endl;
			cout << "\tfrom the " << bbvfilename << " file, and " << target_valence[i] << " based on the " << datfilename << " file!" << endl;
			cout << warning << endl;
			return(0);
		}
	}
	
	if (ReadBin(file, temp_d, nconstraints, "weights", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	
	int index = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr+CoordNumbConst::tot_subconst+RunParams::navcoord;
	index += RunParams::ncommonneigh + RunParams::nsecondneigh;
	for (i = 0; i < nconstraints; i++)
	{
		if (fabs(temp_d[i] - weights[i]) > LOAD_TOL && sigma_percentage[index + i] > 0)
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The weight for the " << i + 1 << ". bond valence sum constraints is " << temp_d[i] << endl;
			cout << "\tfrom the " << bbvfilename << " file, and " << weights[i] << " based on the " << datfilename << " file!" << endl;
			cout << warning << endl;
			return(0);
		}
	}
	if (ReadBin(file, temp_i, nconstraints, "ncentral", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		if (temp_i[i] != ncentral[i])
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of central atoms for the " << i + 1 << ". bond valence sum constraint is " << temp_i[i] << endl;
			cout << "\tfrom the " << bbvfilename << " file, and " << ncentral[i] << " based on the " << datfilename << " file!" << endl;
			cout << warning << endl;
			return(0);
		}
	}
	if (ReadBin(file, temp_i, nconstraints, "n_neightype", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		if (temp_i[i] != n_neightype[i])
		{
			cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of neighbour types for the " << i + 1 << ". bond valence sum constraint is " << temp_i[i] << endl;
			cout << "\tfrom the " << bbvfilename << " file, and " << n_neightype[i] << " based on the " << datfilename << " file!" << endl;
			cout << warning << endl;
			return(0);
		}
	}
	if (ReadBin(file, temp_i, tot_neightype, "neighbours", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		for (j = cum_n_neightype[i]; i < cum_n_neightype[i + 1]; i++)
		{
			if (temp_i[j] != neighbours[j])
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The neighbour type for the " << i + 1 << ". bond valence sum constraint's ";
				cout<<j- cum_n_neightype[i]+1<<". neighbour is " << temp_i[j] << endl;
				cout << "\tfrom the " << bbvfilename << " file, and " << neighbours[j] << " based on the " << datfilename << " file!" << endl;
				cout << warning << endl;
				return(0);
			}
		}
	}
	if (ReadBin(file, temp_d, tot_neightype, "dmax", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		for (j = cum_n_neightype[i]; i < cum_n_neightype[i + 1]; i++)
		{
			if (fabs(temp_d[j]- dmax[j])>LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The maximum distance for the " << i + 1 << ". bond valence sum constraint's ";
				cout << j - cum_n_neightype[i] + 1 << ". neighbour is " << temp_d[j] << " A"<<endl;
				cout << "\tfrom the " << bbvfilename << " file, and " << dmax[j] << " A based on the " << datfilename << " file!" << endl;
				cout << warning << endl;
				return(0);
			}
		}
	}
	if (ReadBin(file, temp_d, tot_neightype, "R0", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		for (j = cum_n_neightype[i]; i < cum_n_neightype[i + 1]; i++)
		{
			if (fabs(temp_d[j] - R0[j]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The R0 for the " << i + 1 << ". bond valence sum constraint's ";
				cout << j - cum_n_neightype[i] + 1 << ". neighbour is " << temp_d[j] << " A" << endl;
				cout << "\tfrom the " << bbvfilename << " file, and " << R0[j] << " A based on the " << datfilename << " file!" << endl;
				cout << warning << endl;
				return(0);
			}
		}
	}
	if (ReadBin(file, temp_d, tot_neightype, "b", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	for (i = 0; i < nconstraints; i++)
	{
		for (j = cum_n_neightype[i]; i < cum_n_neightype[i + 1]; i++)
		{
			if (fabs(temp_d[j] - b[j]) > LOAD_TOL)
			{
				cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The b parameter for the " << i + 1 << ". bond valence sum constraint's ";
				cout << j - cum_n_neightype[i] + 1 << ". neighbour is " << temp_d[j] << " A" << endl;
				cout << "\tfrom the " << bbvfilename << " file, and " << b[j] << " A based on the " << datfilename << " file!" << endl;
				cout << warning << endl;
				return(0);
			}
		}
	}
	if (ReadBin(file, valence, ncentral_tot, "valence", "BondValenceSumConst::LoadBinary", bbvfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	
	delete[] temp_i;
	delete[] temp_d;
	
	if (::debug)
		cout << "Loading of the bond valence sum constraint was successful" << endl;
	return (1);//loading was successful
};

//saving the valence in binary format
void BondValenceSumConst::SaveBinary() const
{
	ofstream file;
		OpenFile(file, bbvfilename, "BondValenceSumConst::SaveBinary", 1);//open file, check, whether it was successfully opened

	file.write(reinterpret_cast<const char *>(&nconstraints), sizeof(nconstraints));
	if (nconstraints < 1)
	{
		file.close();
		return;
	}
	file.write(reinterpret_cast<const char *>(central), nconstraints * sizeof(*central));
	file.write(reinterpret_cast<const char *>(target_valence), nconstraints * sizeof(*target_valence));
	file.write(reinterpret_cast<const char *>(weights), nconstraints * sizeof(*weights));
	file.write(reinterpret_cast<const char *>(ncentral), nconstraints * sizeof(*ncentral));

	file.write(reinterpret_cast<const char *>(n_neightype),  nconstraints* sizeof(*n_neightype));

	file.write(reinterpret_cast<const char *>(neighbours), tot_neightype* sizeof(*neighbours));
	file.write(reinterpret_cast<const char *>(dmax), tot_neightype * sizeof(*dmax));
	file.write(reinterpret_cast<const char *>(R0), tot_neightype * sizeof(*R0));
	file.write(reinterpret_cast<const char *>(b), tot_neightype * sizeof(*b));

	file.write(reinterpret_cast<const char *>(valence), ncentral_tot * sizeof(*valence));

	file.close();
};


//=======loading the BondValenceSumConst object from file============
//this is not used, only the binary Load, but was not removed
/*int BondValenceSumConst::Load(double *sigma_percentage)
{
	int i, k, dumb, ineightype;
	ifstream file;

	double real_dumb;

	SafeOpenTextFile(file, bvsfilename);
	if (CheckFileState(file, "BondValenceSumConst::Load", bvsfilename) == 0)
	{
		cout << "\nWARNING(" << ++warn << "): There is no open file to load the BondValenceSumConst object from" << endl;
		if (RunParams::nicoord > 0)
			cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
		else
			cout << " and bond valences will be calculated!" << endl;
		return(0);//loading was not successful
	}

	cout << "\nLoading the valences into the BondValenceSumConst object" << endl;
	//Reading the parameters for the constraint, which were already given in the .dat file
	//The parameters will be compared  to check consistency
	SkipLine(file, bvsfilename, 1);//comment
	SkipLine(file, bvsfilename, 1);
	dumb = ReadThisLine(file, 1, 1, "BondValenceSumConst nconstraint", "BondValenceSumConst::Load", bvsfilename);//number of constraints
	if (dumb != nconstraints)//checking consistency for number of constraints
	{
		cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
		cout << "\tNumber of constraint is " << nconstraints << " from the .dat file," << endl;
		cout << "\tand " << dumb << " from the .bvs file!" << endl;
		cout << "\tLoading was not successful, histogram";
		if (RunParams::nicoord > 0)
			cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
		else
			cout << " and bond valences will be calculated!" << endl; 
		return(0);//loading was not successful
	}

	for (i = 0; i < 3; i++)
		SkipLine(file, bvsfilename, 1);//skipping the comments and blank lines

	for (i = 0; i < nconstraints; i++)
	{
		dumb = ReadThisLine(file, 3, 1, "index", "BondValenceSumConst::Load");//index
		dumb = ReadThisLine(file, 3, 1, "central type", "BondValenceSumConst::Load");//type of central atom (starts with 1 in file!)
		if (central[i] != (dumb - 1))//checking consistency for type of central atom
		{
			cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
			cout << "\tThe type of the central atom for the " << i + 1 << ". constraint is " << central[i] + 1 << " from the .dat file," << endl;
			cout << "\tand " << dumb << " from the .bvs file!" << endl;
			cout << "\tLoading was not successful, histogram";
			if (RunParams::nicoord > 0)
				cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
			else
				cout << " and bond valences will be calculated!" << endl;
			return(0);//loading was not successful
		}

		dumb = ReadThisLine(file, 3, 1, "number of neighbour types", "BondValenceSumConst::Load");//number of neighbour atom types
		if (n_neightype[i] != (dumb))//checking consistency fornumber of neighbour atom types
		{
			cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
			cout << "\tNumber of neighbour atom types for the " << i + 1 << ". constraint is " << n_neightype[i] << " from the .dat file," << endl;
			cout << "\tand " << dumb << " from the .bvs file!" << endl;
			cout << "\tLoading was not successful, histogram";
			if (RunParams::nicoord > 0)
				cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
			else
				cout << " and bond valences will be calculated!" << endl;
			return(0);//loading was not successful
		};
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
		{
			dumb = ReadThisLine(file, 3, 1, "neighbour type", "BondValenceSumConst::Load");//type of neighbour atom
			if (neighbours[cum_n_neightype[i] + ineightype] != (dumb - 1))//checking consistency for type of neighbour atom
			{
				cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
				cout << "\t" << ineightype + 1 << ". type of neighbour atom for the " << i + 1 << ". constraint is " << neighbours[cum_n_neightype[i] + ineightype] + 1 << " from the .dat file," << endl;
				cout << "\tand " << dumb << " from the .bvs file!" << endl;
				cout << "\tLoading was not successful, histogram";
				if (RunParams::nicoord > 0)
					cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
				else
					cout << " and bond valences will be calculated!" << endl;
				return(0);//loading was not successful
			}
		};


		
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
		{
			real_dumb = ReadThisLine(file, 3, 1.0, "maximum distance", "BondValenceSumConst::Load");//superior distance
			if (fabs(dmax[cum_n_neightype[i] + ineightype] - real_dumb) > LOAD_TOL)//checking consistency
			{
				cout.precision(12);
				cout.setf(ios::right, ios::scientific);
				cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
				cout << "\tSuperior distance for the " << i + 1 << ". constraint's " << ineightype + 1 << ". type of neighbour atom is " << dmax[cum_n_neightype[i] + ineightype] << " from the .dat file," << endl;
				cout << "\tand " << real_dumb << " from the .bvs file!" << endl;
				cout << "\tLoading was not successful, histogram";
				if (RunParams::nicoord > 0)
					cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
				else
					cout << " and bond valences will be calculated!" << endl;
				cout.precision(6);
				cout.unsetf(ios::scientific);
				cout.unsetf(ios::right);
				return(0);//loading was not successful
			}
		};
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
		{
			real_dumb = ReadThisLine(file, 3, 1.0, "R0", "BondValenceSumConst::Load");//superior distance
			if (fabs(R0[cum_n_neightype[i] + ineightype] - real_dumb) > LOAD_TOL)//checking consistency
			{
				cout.precision(12);
				cout.setf(ios::right, ios::scientific);
				cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
				cout << "\tR0 for the " << i + 1 << ". constraint's " << ineightype + 1 << ". type of neighbour atom is " << dmax[cum_n_neightype[i] + ineightype] << " from the .dat file," << endl;
				cout << "\tand " << real_dumb << " from the .bvs file!" << endl;
				cout << "\tLoading was not successful, histogram";
				if (RunParams::nicoord > 0)
					cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
				else
					cout << " and bond valences will be calculated!" << endl;
				cout.precision(6);
				cout.unsetf(ios::scientific);
				cout.unsetf(ios::right);
				return(0);//loading was not successful
			}
		};
		for (ineightype = 0; ineightype < n_neightype[i]; ineightype++)
		{
			real_dumb = ReadThisLine(file, 3, 1.0, "b parameter", "BondValenceSumConst::Load");//superior distance
			if (fabs(b[cum_n_neightype[i] + ineightype] - real_dumb) > LOAD_TOL)//checking consistency
			{
				cout.precision(12);
				cout.setf(ios::right, ios::scientific);
				cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
				cout << "\tb parameter for the " << i + 1 << ". constraint's " << ineightype + 1 << ". type of neighbour atom is " << dmax[cum_n_neightype[i] + ineightype] << " from the .dat file," << endl;
				cout << "\tand " << real_dumb << " from the .bvs file!" << endl;
				cout << "\tLoading was not successful, histogram";
				if (RunParams::nicoord > 0)
					cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
				else
					cout << " and bond valences will be calculated!" << endl;
				cout.precision(6);
				cout.unsetf(ios::scientific);
				cout.unsetf(ios::right);
				return(0);//loading was not successful
			}
		};
		
		real_dumb = ReadThisLine(file, 3, 1, "target valence", "BondValenceSumConst::Load");//Desired valence
		if (fabs(target_valence[i]-real_dumb)>LOAD_TOL)//checking consistency for subconstraints
		{
			cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
			cout << "\tDesired valence for the " << i+ 1 << ". constraint is " << target_valence[i] << " from the .dat file," << endl;
			cout << "\tand " << real_dumb << " from the .bvs file!" << endl;
			cout << "\tThe value in the .dat file will be used!" << endl;
		};

		real_dumb = ReadThisLine(file, 3, 1.0, "sigma", "BondValenceSumConst::Load");//weight parameter
		int index = RunParams::ngr + RunParams::nsq + RunParams::nfq + RunParams::nfg + RunParams::nek + RunParams::ncosdistr+ CoordNumbConst::tot_subconst + RunParams::navcoord ;
		index += RunParams::ncommonneigh + RunParams::nsecondneigh;
		if (fabs(weights[i] - real_dumb)>LOAD_TOL && sigma_percentage[index + i] > 0)//checking consistency
		{
			cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
			cout << "\tWeight parameter for the " << i + 1 << " constraint is " << weights[i] << " from the .dat file," << endl;
			cout << "\tand " << real_dumb << " from the .bvs file!" << endl;
			cout << "\tThe value in the .dat file will be used!" << endl;
		};
		
		dumb = ReadThisLine(file, 1, 1, "ncentral", "BondValenceSumConst::Load", bvsfilename);//number of central atoms
		if (ncentral[i] != dumb)//checking consistency
		{
			cout << "\nWARNING(" << ++warn << "): BondValenceSumConst::Load" << endl;
			cout << "\tNumber of central atoms for the " << i + 1 << " constraint is " << ncentral[i] << " from the .dat file," << endl;
			cout << "\tand " << dumb << " from the .bvs file!" << endl;
			cout << "\tLoading was not successful, histogram";
			if (RunParams::nicoord > 0)
				cout << ", coordination numbers and bond valence sums will be calculated!" << endl;
			else
				cout << " and bond valences will be calculated!" << endl;
			return(0);//loading was not successful
		};
	}

	SkipLine(file, bvsfilename, 1);//skipping the comments and blank lines

	for (i = 0; i < nconstraints; i++)//for each constraint
	{
		SkipLine(file, bvsfilename, 1);
		SkipLine(file, bvsfilename, 1);//constraint's index

		for (k = 0; k < ncentral[i]; k++)//for every central atom in this constraint
		{
			dumb = ReadThisLine(file, 3, 1, "serial index", "BondValenceSumConst::Load");
			if (i == nconstraints - 1 && k == ncentral[i] - 1)//last point of last constraint, no error, if no more line in file
				BondValence(i, k) = ReadThisLine(file, 1, 1.0, "valence", "BondValenceSumConst::Load");
			else
				BondValence(i, k) = ReadThisLine(file, 1, 1.0, "valence", "BondValenceSumConst::Load", bvsfilename);
		}
		
	}

	//Checking, whether loading was successful
	if (CheckReadFileState(file, "BondValenceSumConst::Load", bvsfilename))
		return(1);//load was successful
	else
	{
		cout << "\nWARNING(" << ++warn << "): Loading was not successful, histogram, coordination numbers and bond valence sums will be calculated!" << endl;
		return(0);//loading was not successful
	}
	next_line_pos = -1;
	file.close();
	return(1);//load was successful
}*/

#endif //_ADVANCED_GEOM_CONST