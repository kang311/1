//source ChiSquared.cpp
//Last changed 03.03.2023

//chi2 is defined, that in case of renormalization the experimental data is renormalized

//If the potential option is on to use  non-bonded potential as well,
//total2 is used for all the potential related total chi2 contributions, 
//as the non-bonded pot can be negative, and its contribution would decrease, or even cancel out the chi2 contribution of the data sets
//A separate acceptance process decide, whether the move is acceptable based on the chenge in total2 similarly to the orignal total 
//If only bonded potential is used for keeping the molecules together, then their contribution will be added to the normal total
//as those contribution cannot be negative

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//classes1.h is included through this

//Defining the static members
int ChiSquared::chisize;//number of elements in the chicomp array
int ChiSquared::ngr;//number of g(r) data sets
int ChiSquared::nsq;//number of S(Q) data sets
int ChiSquared::nfq;//number of F(Q) data sets
int ChiSquared::nfg;//number of F(g) data sets
int ChiSquared::nek;//number of E(k) data sets
int ChiSquared::nacc;//number of average coordination constraints
int ChiSquared::nicc;//number of individual coordination constraints
int ChiSquared::ncos;//number of individual cosine distribution of bond angle constraints
int ChiSquared::ncommonneigh;//number of individual common neighbour constraints
int ChiSquared::nsecondneigh;//number of individual second neighbour constraints
int ChiSquared::nbvs;//number of bond valence sum constraints
int ChiSquared::tot_cc_subconst;//total number of coordination subconstraints
int ChiSquared::calc_sigma;//indicator, whether the sigma should be calculated
int ChiSquared::failed_nlr=0;//number of failed nonlinear regression since the last successful 
int *ChiSquared::tot_failed_nlr;//total number for failed non-lin regressions for an I(Q) set
double ChiSquared::total;//sum of the chi squared components


int ChiSquared::nbonds=0;//number of bond types to use in chisquare calculation (if there are bonds, Topology::bond_weight_mode==0 ->1, or nbond_types, 0 otherwise)
int ChiSquared::nangles=0;//number of angle types to use in chisquare calculation(if there are angles,  Topology::angle_weight_mode==0 ->1, or nangle_types, 0 otherwise)
int ChiSquared::nperdihs=0;//number of periodic dihedral types to use in chisquare calculation(if there are perdihs,  Topology::perdih_weight_mode==0 ->1, or nperdihedral_types, 0 otherwise)
int ChiSquared::nharmdihs=0;//number of harmonic dihedral types to use in chisquare calculation(if there are harmdihs,  Topology::harmdih_weight_mode==0 ->1, or nharmdihedral_types, 0 otherwise)
int ChiSquared::nRBdihs=0;//number of RB dihedral types to use in chisquare calculation(if there are RBdihs,  Topology::RBdih_weight_mode==0 ->1, or nRBdihedral_types, 0 otherwise)
int ChiSquared::nNB=0;//number of non-bonded interaction types or 1 for Aenet
int ChiSquared::npair_types=0;//number of 1-4 non-bonded interaction types
int ChiSquared::chi_potstart=0;//beginning of the first potential related member in the chicomp array
double *ChiSquared::p_pot_chitotal=0;//pointer to the total chi2 value to be usd for the potential related components
double ChiSquared::total2=0.0;//sum of the chi squared potential components if non-bonded potential option is on
double ChiSquared::pot_chi2_low_limit=0.0;//if the potential does not have to be 


#ifdef _LOCAL_INV
	int ChiSquared::nlocint;//number of local inv intervals
	int ChiSquared::sum_loc_natoms_tot;//sum of loc_natoms, needed for some arrays
	int *ChiSquared::sum_loc_natoms;//number of atoms for each interval
	int *ChiSquared::cumalpart_loc_atoms;//cumulative number of loc atoms for each partials, where 1-2<> 2-1 (ntypes^2 partials)
	int *ChiSquared::nused_loc;//number of used threads for local inv.
	int *ChiSquared::loc_natoms_min;//first atom to be used from each type in case of loc_chi2_mode=1
	int *ChiSquared::first_loc_ind;//the first histogram bin to use for a central atom/per central atom for each thread
	int *ChiSquared::starting_count;//the starting count for the first used bin/per central atom for each thread
	int *ChiSquared::loc_natoms;//X atom to be used from each type in case of loc_chi2_mode=1
	int *ChiSquared::av_loc_hist_offset;//offset to the begininnig of each type in av_loc_hist for a thread from the beginning of this thread's segment
	double *ChiSquared::av_loc_hist;//the averaged histogram bins according to the distance of the i-th atom of a type
	
	ChiSquared::TCalcLocChi2Thread ChiSquared::CalcLocChi2Thread;//pointer to the actual local invariance chi2 function
#endif


#ifdef _AENET
	ChiSquared::ChiSquared(ExptsData &edat,CosDistrConst &cosc, HistoSet &histog, Threads &thread_object, Aenet &aen)
		:edata(edat), cosconst(cosc), hist(histog), thread_obj(thread_object), aenet(aen)
#else
	ChiSquared::ChiSquared(ExptsData &edat,CosDistrConst &cosc, HistoSet &histog, Threads &thread_object)
		:edata(edat), cosconst(cosc), hist(histog), thread_obj(thread_object)
#endif
{ 
	//constructor
	int i,index,j,size;
	double *pexpt,*pq,*pr,*pg;
	
	SetArraysize(&name,NAME_SIZE,"name","ChiSquared::ChiSquared");
	SetArraysize(&chicomp,chisize,"chicomp","ChiSquared::ChiSquared");
	
	//initialisation of chisq. values to 0.0
	for (i = 0; i < chisize; i++)
		*(chicomp + i) = 0.0;

	
	//copying the original sigma values in case it should be reset
	SetArraysize(&sigma_percentage,chisize,"sigma_percentage","ChiSquared::ChiSquared");
	
#ifdef _LOCAL_INV
	int padding_offset;
	loc_bin_width = (HistoSet::max_loc_r - HistoSet::min_loc_r) / HistoSet::max_loc_nbins * SimpleCfg::boxedge;//in A
	//values for the local invariance chi2 for each thread, padding has to be used
	padding_offset = (RunParams::nthreads > 1 ? (int)((CACHE_PADDING - sizeof(*loc_chi)) / sizeof(*loc_chi)) : 0);//number of dummy double elements
	//although not all the threads write to each local inv intervals chi2, all will have an array element
	SetArraysize(&loc_chi, RunParams::nthreads * nlocint + (RunParams::nthreads - 1) * padding_offset, "loc_chi", "ChiSquared::ChiSquared");
	SetArraysize(&loc_chi_finder, RunParams::nthreads, "loc_chi_finder", "ChiSquared::ChiSquared");
	//loc_chi_finder will show the beginning of this thread of the loc_chi array

	for (i = 0; i < RunParams::nthreads; i++)
		loc_chi_finder[i] = loc_chi + i * (nlocint + padding_offset);



	if (RunParams::loc_chi2_mode == 1)
	{
		for (i = 0; i < nlocint; i++)
			RunParams::loc_npoints[i] = sum_loc_natoms[i] * SimpleCfg::ntotal;
		SetArraysize(&min_loc_r_now, (int)pow(SimpleCfg::ntypes, 2) * nlocint, "min_loc_r_now", "ChiSquared::ChiSquared");
		SetArraysize(&max_loc_r_now, (int)pow(SimpleCfg::ntypes, 2) * nlocint, "max_loc_r_now", "ChiSquared::ChiSquared");
		for (i = 0; i < pow(SimpleCfg::ntypes, 2) * nlocint; i++)
		{
			min_loc_r_now[i] = 0;
			max_loc_r_now[i] = 0;
		}
	}
	else
		RunParams::loc_npoints[0] = HistoSet::sum_bins * SimpleCfg::ntotal;
	for (i = 0; i < nlocint; i++)
		RunParams::ntotal_points += RunParams::loc_npoints[i];


#endif
	
	index=0;
	for( i=0;i<ngr;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=ExptsData::grsigma[i];
			if (RunParams::lead_series_ind==i)
			{
				if (ExptsData::grsigma[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th g(r) series) cannot be negative (no rescaling for this!)"<<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=ExptsData::grsigma[i];
			}
			if (ExptsData::grsigma[i]<0)
				ExptsData::grsigma[i]=1.0;
			
			
		}
		index++;
	}
	for( i=0;i<nsq;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=ExptsData::sqsigma[i];
			if (RunParams::lead_series_ind==index)
			{
				if (ExptsData::sqsigma[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th S(Q) series) cannot be scalable!" <<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=ExptsData::sqsigma[i];
			}
			if (ExptsData::sqsigma[i]<0)
				ExptsData::sqsigma[i]=1.0;

		}
		index++;
	}
	for( i=0;i<nfq;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=ExptsData::fqsigma[i];
			if (RunParams::lead_series_ind==index)
			{
				if (ExptsData::fqsigma[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th F(Q) series) cannot be scalable!" <<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=ExptsData::fqsigma[i];
			}
			if (ExptsData::fqsigma[i]<0)
				ExptsData::fqsigma[i]=1.0;
		}
		index++;
	}
	for (i = 0; i < nfg; i++)
	{
		if (calc_sigma == 1)
		{
			sigma_percentage[index] = ExptsData::fgsigma[i];
			if (RunParams::lead_series_ind == index)
			{
				if (ExptsData::fgsigma[i] < 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The sigma value for the leading series (" << i+1 << "-th F(g) series) cannot be scalable!"  << endl;
					cout << "Cannot run this way, exiting..." << endl;
					CleanExit();
				}

				lead_sigma = ExptsData::fgsigma[i];
			}
			if (ExptsData::fgsigma[i] < 0)
				ExptsData::fgsigma[i] = 1.0;
		}
		index++;
	}
	for( i=0;i<nek;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=ExptsData::eksigma[i];
			if (RunParams::lead_series_ind==index)
			{
				if (ExptsData::eksigma[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th E(k) series) cannot be scalable!" <<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=ExptsData::eksigma[i];
			}
			if (ExptsData::eksigma[i]<0)
				ExptsData::eksigma[i]=1.0;
		}
		index++;
	}
	for( i=0;i<ncos;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=CosDistrConst::weights[i];
			if (RunParams::lead_series_ind==index)
			{
				if (CosDistrConst::weights[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th cosine distribution of angles constraint) cannot be scalable!" <<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=CosDistrConst::weights[i];
			}
			if (CosDistrConst::weights[i]<0)
				CosDistrConst::weights[i]=1.0;
		}
		index++;
	}


	for( i=0;i<tot_cc_subconst;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=CoordNumbConst::weights[i];
			if (RunParams::lead_series_ind==index)
			{
				if (CoordNumbConst::weights[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th coordination number constraint) cannot be scalable!" <<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=CoordNumbConst::weights[i];
			}
			if (CoordNumbConst::weights[i]<0)
				CoordNumbConst::weights[i]=1.0;
		}
		index++;
	}
	
	for( i=0;i<nacc;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=AvCoordConst::weights[i];
			if (RunParams::lead_series_ind==index)
			{
				if (AvCoordConst::weights[i]<0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout<<"The sigma value for the leading series ("<<i+1<<"-th average coordination constraint) cannot be scalable!" <<endl;
					cout<<"Cannot run this way, exiting..."<<endl;
					CleanExit();
				}

				lead_sigma=AvCoordConst::weights[i];
			}
			if (AvCoordConst::weights[i]<0)
				AvCoordConst::weights[i]=1.0;
		}
		index++;
	}
#ifdef _ADVANCED_GEOM_CONST
	for (i = 0; i < ncommonneigh; i++)
	{
		if (calc_sigma == 1)
		{
			sigma_percentage[index] = CommonNeighConst::weights[i];
			if (RunParams::lead_series_ind == index)
			{
				if (CommonNeighConst::weights[i] < 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The sigma value for the leading series (" << i+1 << "-th common neighbour constraint) cannot be scalable!" << endl;
					cout << "Cannot run this way, exiting..." << endl;
					CleanExit();
				}

				lead_sigma = CommonNeighConst::weights[i];
			}
			if (CommonNeighConst::weights[i] < 0)
				CommonNeighConst::weights[i] = 1.0;
		}
		index++;
	}
	for (i = 0; i < nsecondneigh; i++)
	{
		if (calc_sigma == 1)
		{
			sigma_percentage[index] = SecondNeighConst::weights[i];
			if (RunParams::lead_series_ind == index)
			{
				if (SecondNeighConst::weights[i] < 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The sigma value for the leading series (" << i+1 << "-th second neighbour constraint) cannot be scalable!" << endl;
					cout << "Cannot run this way, exiting..." << endl;
					CleanExit();
				}

				lead_sigma = SecondNeighConst::weights[i];
			}
			if (SecondNeighConst::weights[i] < 0)
				SecondNeighConst::weights[i] = 1.0;
		}
		index++;
	}
	for (i = 0; i < nbvs; i++)
	{
		if (calc_sigma == 1)
		{
			sigma_percentage[index] = BondValenceSumConst::weights[i];
			if (RunParams::lead_series_ind == index)
			{
				if (BondValenceSumConst::weights[i] < 0)
				{
					cout << "\n*****ERROR*****" << endl;
					cout << "The sigma value for the leading series (" << i + 1 << "-th bond valence sum const constraint) cannot be scalable!"  << endl;
					cout << "Cannot run this way, exiting..." << endl;
					CleanExit();
				}

				lead_sigma = BondValenceSumConst::weights[i];
			}
			if (BondValenceSumConst::weights[i] < 0)
				BondValenceSumConst::weights[i] = 1.0;
		}
		index++;
	}

#endif	
#ifdef _LOCAL_INV
	for (i=0;i<nlocint;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=RunParams::loc_inv_sigma[i];
			if (RunParams::lead_series_ind==index)
				{
					if (RunParams::loc_inv_sigma[i]<0)
					{
						cout << "\n*****ERROR*****" << endl;
						cout<<"The sigma value for the leading series ("<<i+1<<"-th local inveraianvce) cannot be scalable!" <<endl;
						cout<<"Cannot run this way, exiting..."<<endl;
						CleanExit();
					}

					lead_sigma=RunParams::loc_inv_sigma[i];
				}
			if (RunParams::loc_inv_sigma[i]<0)
				RunParams::loc_inv_sigma[i]=1.0;
		}
		index++;
	}
#endif
	if (RunParams::potential>0)
		lead_ser_pot=RunParams::lead_series_ind2;
	else
		lead_ser_pot=RunParams::lead_series_ind;
#ifdef _AENET
	if  (RunParams::potential==20)
	{
		if (calc_sigma==1)
			{
				sigma_percentage[index]=RunParams::aenet_weight;
				if (RunParams::aenet_weight<0)
					RunParams::aenet_weight=1.0;
			}
			index++;
	}
#endif	
	if (RunParams::potential==1 || RunParams::potential==10)
	{
		for( i=0;i<nNB;i++)
		{
			if (calc_sigma==1)
			{
				sigma_percentage[index]=RunParams::vdW_weight[i];
				if (RunParams::vdW_weight[i]<0)
					RunParams::vdW_weight[i]=1.0;
			}
			index++;
		}
	}
	if (RunParams::potential == 1)
	{
		for (i = 0; i < nNB; i++)
		{
			if (calc_sigma == 1)
			{
				sigma_percentage[index] = RunParams::Coulomb_weight[i];
				if (RunParams::Coulomb_weight[i] < 0 || fabs(RunParams::Coulomb_weight[i] - 0.0) < TOLERANCE)
					RunParams::Coulomb_weight[i] = 1.0;
			}
			index++;
		}
		if (Topology::npair_types > 0)
		{
			for (i = 0; i < nNB; i++)
			{
				if (calc_sigma == 1)//1-4 vdW
					sigma_percentage[index] = sigma_percentage[index - 2 * nNB];//the same is used as for the normal vdW

				index++;
			}

			for (i = 0; i < nNB; i++)
			{
				if (calc_sigma == 1)//1-4 Coulomb
					sigma_percentage[index] = sigma_percentage[index - 2 * nNB];//the same is used as for the normal Coulomb

				index++;
			}
		}
	}
	for( i=0;i<nbonds;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=Topology::bond_sigma[i];
			if (Topology::bond_sigma[i]<0 || fabs(Topology::bond_sigma[i]-0.0)<TOLERANCE)
				Topology::bond_sigma[i]=1.0;
		}
		index++;
	}

	for( i=0;i<nangles;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=Topology::angle_sigma[i];
			if (Topology::angle_sigma[i]<0 || fabs(Topology::angle_sigma[i]-0.0)<TOLERANCE)
				Topology::angle_sigma[i]=1.0;
		}
		index++;
	}

	for( i=0;i<nperdihs;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=Topology::perdihedral_sigma[i];
			if (Topology::perdihedral_sigma[i]<0 || fabs(Topology::perdihedral_sigma[i]-0.0)<TOLERANCE)
				Topology::perdihedral_sigma[i]=1.0;
		}
		index++;
	}
	for( i=0;i<nharmdihs;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=Topology::harmdihedral_sigma[i];
			if (Topology::harmdihedral_sigma[i]<0 || fabs(Topology::harmdihedral_sigma[i]-0.0)<TOLERANCE)
				Topology::harmdihedral_sigma[i]=1.0;
		}
		index++;
	}
	for( i=0;i<nRBdihs;i++)
	{
		if (calc_sigma==1)
		{
			sigma_percentage[index]=Topology::RBdihedral_sigma[i];
			if (Topology::RBdihedral_sigma[i]<0|| fabs(Topology::RBdihedral_sigma[i]-0.0)<TOLERANCE)
				Topology::RBdihedral_sigma[i]=1.0;
		}
		index++;
	}

	
	size = ngr + nsq + nfq + nfg + nek + ncos;//size of the sums arrays having an element for the CosDistrConst

	//CosDistrConst has to be the last type of constraint, as it has no element in all the other arrays!!!
	SetArraysize(&ff,size,"ff","ChiSquared::ChiSquared");//sums of experimental values squared
	SetArraysize(&ss,size,"ss","ChiSquared::ChiSquared");//sums of calculated values squared
	SetArraysize(&fs,size,"fs","ChiSquared::ChiSquared");//sums of cross products
	
	for(i=0;i<size;i++)
	{
		*(ff+i)=0.0;
		*(fs+i)=0.0;
		*(ss+i)=0.0;
	}

	size=ngr+nsq+nfq+nfg+nek;//size of the sums arrays (CosDistrConst cannot be renormalized, so these arrays are not needed for it)
	//presently EXAFS is not normalized, but there is a place for it
							 //X means r or Q depending on the data set
	SetArraysize(&f,size,"f","ChiSquared::ChiSquared");//sums of experimental values
	SetArraysize(&fx,size,"fx","ChiSquared::ChiSquared");//sums of experimental values multiplied by X
	SetArraysize(&fx2,size,"fx2","ChiSquared::ChiSquared");//sums of experimental values multiplied by X square
	SetArraysize(&fx3,size,"fx3","ChiSquared::ChiSquared");//sums of experimental values multiplied by X ^3
	SetArraysize(&s,size,"s","ChiSquared::ChiSquared");//sums of calculated values
	SetArraysize(&sx,size,"sx","ChiSquared::ChiSquared");//sums of calculated values multiplied by X
	SetArraysize(&sx2,size,"sx2","ChiSquared::ChiSquared");//sums of calculated values multiplied by X square
	SetArraysize(&sx3,size,"sx3","ChiSquared::ChiSquared");//sums of calculated values multiplied by X ^3
	SetArraysize(&x,size,"x","ChiSquared::ChiSquared");//sums of X values 
	SetArraysize(&x2,size,"x2","ChiSquared::ChiSquared");//sums of X^2 values 
	SetArraysize(&x3,size,"x3","ChiSquared::ChiSquared");//sums of X^3 values 
	SetArraysize(&x4,size,"x4","ChiSquared::ChiSquared");//sums of X^4 values 
	SetArraysize(&x5,size,"x5","ChiSquared::ChiSquared");//sums of X^5values 
	SetArraysize(&x6,size,"x6","ChiSquared::ChiSquared");//sums of X^6values 
	SetArraysize(&a,size,"a","ChiSquared::ChiSquared");//multiplicative factors
	SetArraysize(&b,size,"b","ChiSquared::ChiSquared");//additive factors
	SetArraysize(&c,size,"c","ChiSquared::ChiSquared");//linear coefficient
	SetArraysize(&d,size,"d","ChiSquared::ChiSquared");//quadratic coefficient
	SetArraysize(&e,size,"e","ChiSquared::ChiSquared");//cubic coefficient
	if (ExptsData::is_IQ>0)
	{
		SetArraysize(&alpha, nfq, "alpha", "ChiSquared::ChiSquared");//alpha coeff for Xray I(Q) fit
		SetArraysize(&C, nfq, "C", "ChiSquared::ChiSquared");//for sum of  B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j)
		SetArraysize(&CC, nfq, "CC", "ChiSquared::ChiSquared");//for sum of [B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j)]^2
		SetArraysize(&fC, nfq, "fC", "ChiSquared::ChiSquared");//for sum of  Iexp_j*B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j)
		SetArraysize(&sC, nfq, "sC", "ChiSquared::ChiSquared");//for sum of  Icalc_j*B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j)
		
		SetArraysize(&last_nlr_chi2fr, nfq, "last_nlr_chi2fr", "ChiSquared::ChiSquared");//for the chi2_end/chi2_start of the last accepted move during non-lin reg
		SetArraysize(&maxdec_nlr_chi2fr, nfq, "maxdec_nlr_chi2fr", "ChiSquared::ChiSquared");//for the chi2_end/chi2_start with the largest decrease during non-lin reg during the whole simualtion
		
	}
	else
	{
		alpha = NULL;
		C = NULL;
		CC = NULL;
		fC = NULL;
		sC = NULL;
		last_nlr_chi2fr = NULL;
		maxdec_nlr_chi2fr = NULL;

	}
	SetArraysize(&fit_index, size, "fit_index", "ChiSquared::ChiSquared");//fit_index to decide about the renormalization formula coefficient
	SetArraysize(&Rw,size,"Rw","ChiSquared::ChiSquared");//measurement of the goddness of the fit
	//n.b. the g(r) data are not allowed to be changed by an additive factor
	//hence the first ngr values of this array will always be 0.0
	
	
	//initialisation of these arrays
	for(i=0;i<size;i++)
	{
		*(f+i)=0.0;
		*(fx+i)=0.0;
		*(fx2+i)=0.0;
		*(fx3+i)=0.0;
		*(s+i)=0.0;
		*(sx+i)=0.0;
		*(sx2+i)=0.0;
		*(sx3+i)=0.0;
		*(x+i)=0.0;
		*(x2+i)=0.0;
		*(x3+i)=0.0;
		*(x4+i)=0.0;
		*(x5+i)=0.0;
		*(x6+i)=0.0;
		*(a+i)=1.0;
		*(b+i)=0.0;//this can be reset to custom later
		*(c+i)=0.0;
		*(d+i)=0.0;
		*(e+i)=0.0;
	}
	if (ExptsData::is_IQ>0)
	{
		for (i = 0; i < nfq; i++)//initialization to defaultm custom values will be set later, if there is any
		{
			*(alpha + i) = IQ_ALPHA_DEF;
			maxdec_nlr_chi2fr[i] = 1e20;
			last_nlr_chi2fr[i] = 0.0;
		}
		
	}
	
	//Calculation of the arrays depending upon only the experimental data
	
	//g(r), only if cubic is used
	pexpt=edata.gr_gvalues;//sets the pointer to the beginning of g(r) data
	pr=edata.gr_rvalues;//sets the pointer to the beginning of r data for g(r) sets
	for (i=0;i<ngr;i++)//for each data set
	{
		index=i;
		//set fit_index
		fit_index[index] = *(ExptsData::grrenorm + i) + *(ExptsData::groffset + i) * 2 + *(ExptsData::grlinear + i) * 4 + *(ExptsData::grquadratic + i) * 8 + *(ExptsData::grcubic + i) * 16;
		
		//for every data point 
		for (j=0;j<edata.grused[i];j++)
		{
			ff[index]+=*pexpt* (*pexpt);
			if (ExptsData::use_cubic[index])
			{
				f[index]+=*pexpt;
				fx[index]+=(*pexpt * *pr);
				fx2[index]+=(*pexpt * *pr * *pr);
				fx3[index]+=(*pexpt * *pr * *pr * *pr);
				x[index]+=*pr;
				x2[index]+=(*pr * *pr);
				x3[index]+=(*pr * *pr * *pr);
				x4[index]+=(*pr * *pr * *pr * *pr);
				x5[index]+=(*pr * *pr * *pr * *pr * *pr);
				x6[index]+=(*pr * *pr * *pr * *pr * *pr * *pr);
			
			}
			pexpt++;
			pr++;
		}
	}

	//S(Q)
	pexpt=edata.sq_svalues;//sets the pointer to the beginning of S(Q) data
	pq=edata.sq_qvalues;//sets the pointer to the beginning of Q data for S(Q) sets
	for (i=0;i<nsq;i++)//for each data set
	{
		index=i+ngr;
		//for every data point 
		fit_index[index] = *(ExptsData::sqrenorm + i) + *(ExptsData::sqoffset + i) * 2 + *(ExptsData::sqlinear + i) * 4 + *(ExptsData::sqquadratic + i) * 8 + *(ExptsData::sqcubic + i) * 16;
		
		for (j=0;j<edata.sqused[i];j++)
		{
			f[index]+=*pexpt;
			ff[index]+=*pexpt* (*pexpt);
			fx[index]+=(*pexpt * *pq);
			fx2[index]+=(*pexpt * *pq * *pq);
			
			x[index]+=*pq;
			x2[index]+=(*pq * *pq);
			x3[index]+=(*pq * *pq * *pq);
			x4[index]+=(*pq * *pq * *pq * *pq);

			if (ExptsData::use_cubic[index])
			{
				fx3[index]+=(*pexpt * *pq * *pq * *pq);
				x5[index]+=(*pq * *pq * *pq * *pq * *pq);
				x6[index]+=(*pq * *pq * *pq * *pq * *pq * *pq);
			}
			pexpt++;
			pq++;
		}
	}

	//F(Q)
	pexpt=edata.fq_fvalues;//sets the pointer to the beginning of F(Q) data
	pq=edata.fq_qvalues;//sets the pointer to the beginning of Q data for F(Q) sets
	for (i=0;i<nfq;i++)//for each data set
	{
		index=i+ngr+nsq;
		fit_index[index] = *(ExptsData::fqrenorm + i) + *(ExptsData::fqoffset + i) * 2 + *(ExptsData::fqlinear + i) * 4 + *(ExptsData::fqquadratic + i) * 8 + *(ExptsData::fqcubic + i) * 16 + ExptsData::fqrenalpha[i] * 32;
		
		if (ExptsData::fqfitIQ[i])
		{//set the custom values if there is any
			switch (fit_index[index])
			{
				case 1://'a' is fitted only with linear regression, set b and alpha
					b[index] = ExptsData::fqb[i];
				case 3://'a' and 'b' are fitted with linear regression
					alpha[i] = ExptsData::fqalpha[i];//presently only alpha=0 is handled, but this way it could be custom
					break;
				case 32://alpha is fitted with nonlin regression, set custom 'a' and 'b'
					a[index] = ExptsData::fqa[i];
					b[index] = ExptsData::fqb[i];
					if (abs(RunParams::continuation))
						alpha[i] = ExptsData::fqalpha[i];//this will be read from state file, set it here to use as initial guess
					break;
				case 33://a, alpha is fitted with nonlin regression, b=0
				case 35://a,b, alpha is fitted
					if (abs(RunParams::continuation))
					{
						a[index] = ExptsData::fqa[i];
						b[index] = ExptsData::fqb[i];
						alpha[i] = ExptsData::fqalpha[i];//this will be read from state file, set it here to use as initial guess
					}
					break;
			}
			
		}
		
		//for every data point, they are calculated, even if they are not needed
		for (j = 0; j < edata.fqused[i]; j++)//in case of I(Q) mu corr it will be recalculated, if necessary, and the values 
			//for the last accepted move will be stored as for the fmuact
		{
			
			f[index] += *pexpt;
			ff[index] += *pexpt * (*pexpt);
			fx[index] += (*pexpt * *pq);
			fx2[index] += (*pexpt * *pq * *pq);
			x[index] += *pq;
			x2[index] += (*pq * *pq);
			x3[index] += (*pq * *pq * *pq);
			x4[index] += (*pq * *pq * *pq * *pq);

			if (ExptsData::use_cubic[index])
			{
				fx3[index] += (*pexpt * *pq * *pq * *pq);
				x5[index] += (*pq * *pq * *pq * *pq * *pq);
				x6[index] += (*pq * *pq * *pq * *pq * *pq * *pq);
			}
			pexpt++;
			pq++;
		}
		
	}
	//F(g)
	pexpt = edata.fg_fvalues;//sets the pointer to the beginning of F(g) data
	pg = edata.fg_gvalues;//sets the pointer to the beginning of g data for F(g) sets
	for (i = 0; i < nfg; i++)//for each data set
	{
		index = i + ngr + nsq + nfq;
		fit_index[index] = *(ExptsData::fgrenorm + i) + *(ExptsData::fgoffset + i) * 2 + *(ExptsData::fglinear + i) * 4 + *(ExptsData::fgquadratic + i) * 8 + *(ExptsData::fgcubic + i) * 16;
		//for every data point 
		for (j = 0; j < edata.fgused[i]; j++)
		{
			f[index] += *pexpt;
			ff[index] += *pexpt* (*pexpt);
			fx[index] += (*pexpt * *pg);
			fx2[index] += (*pexpt * *pg * *pg);
			x[index] += *pg;
			x2[index] += (*pg * *pg);
			x3[index] += (*pg * *pg * *pg);
			x4[index] += (*pg * *pg * *pg * *pg);

			if (ExptsData::use_cubic[index])
			{
				fx3[index] += (*pexpt * *pg * *pg * *pg);
				x5[index] += (*pg * *pg * *pg * *pg * *pg);
				x6[index] += (*pg * *pg * *pg * *pg * *pg * *pg);
			}
			pexpt++;
			pg++;
		}
	}
	//E(k)
	pexpt = edata.ek_evalues;//sets the pointer to the beginning of E(k) data
	for (i = 0; i < nek; i++)//for each data set
	{
		index = i + ngr + nsq + nfq + nfg;

		//for every data point 
		for (j=0;j<edata.ekused[i];j++)
		{
			f[index]+=*pexpt;
			ff[index]+=*pexpt* (*pexpt);
			pexpt++;
		}
	}

	//CosDistrConst
	pexpt=cosconst.theordistr;
	for (i=0;i<ncos;i++)//for each data set
	{
		index=i+ngr+nsq+nfq+nfg+nek;
		if (CosDistrConst::method[i]==2)
		{
			//The desired theoretical values are zero for each point, no summing is necessary
			pexpt+=CosDistrConst::ncos_bin;//advance the pointer
		}
		else
		{
			pexpt+=CosDistrConst::first_bin[i];//set it to the first meaningful point
			for(j=CosDistrConst::first_bin[i];j<CosDistrConst::last_bin[i];j++)//for all data points
			{
				//array f is not needed, as no renormalisation is allowed
				ff[index]+=*pexpt* (*pexpt);
				pexpt++;
			}
			pexpt+=CosDistrConst::npoints[i]-CosDistrConst::last_bin[i];//to skip the unused points of the this constraint
		}
	}
};	
		
//---------------sets the static members---------------
void ChiSquared::SetChiSquaredParams(RunParams &rundat)
{
	ngr=rundat.ngr;//number of g(r) data sets
	nsq=rundat.nsq;//number of S(Q) data sets
	nfq=rundat.nfq;//number of F(Q) data sets
	nfg = rundat.nfg;//number of F(g) data sets
	nek=rundat.nek;//number of E(k) data sets
	nacc=rundat.navcoord;//number of average coordination constraints
	nicc=rundat.nicoord;//number of individual coordination constraints
	ncos=rundat.ncosdistr;//number of cosine distribution of bond angles constraint
	ncommonneigh = rundat.ncommonneigh;//number of common neighbour constraints
	nsecondneigh = rundat.nsecondneigh;//number of second neighbour constraints
	nbvs = rundat.nbvs;//number of bond valence sum constraints
	tot_cc_subconst=CoordNumbConst::tot_subconst;//total number of subconstraints (effectively the number of coordination constraint
		
	if (ExptsData::is_IQ>0)
	{
		SetArraysize(&tot_failed_nlr, nfq, "tot_failed_nlr", "ChiSquared::SetChiSquaredParams");
		for (int i = 0; i < nfq; i++)
			tot_failed_nlr[i] = 0;
	}
	else
		tot_failed_nlr = NULL;

	//in case of non-bonded potential the leading series will be potential related, most probably the vDW to keep the potential components
	//proportional to each other, but this potential lead series itself can be scaled to the leading series of the data sets
	//if there is no non-bonded interaction, then the bonded interactions will be scaled to the lead series of the data sets 
	if (RunParams::fnc==4)
	{
		nbonds=(Topology::nbond_types>0 ? (Topology::bond_weight_mode ? Topology::nbond_types : 1) : 0);//number of bond types
		nangles=(Topology::nangle_types>0 ? (Topology::angle_weight_mode ? Topology::nangle_types: 1) : 0);//number of angle types
		nRBdihs=(Topology::nRBdihedral_types>0 ? (Topology::RBdih_weight_mode ? Topology::nRBdihedral_types :1) : 0);//number of RBdihedral types
		nperdihs=(Topology::nperdihedral_types >0 ? (Topology::perdih_weight_mode ? Topology::nperdihedral_types : 1) : 0);//number of periodic dihedral types
		nharmdihs=(Topology::nharmdihedral_types>0 ?( Topology::harmdih_weight_mode ? Topology::nharmdihedral_types :1) : 0);//number of harmonic dihedral types
		npair_types=Topology::npair_types;//number of 1-4 interaction types
	}
	nNB = 0;
	switch (RunParams::potential)
	{
	case 1:
	{
		nNB = (RunParams::NB_weight_mode == 1 ? SimpleCfg::npartials : 1);//number of nonbonded interactions /type(either 1, or npartials)
		break;
	}
	case 10:
	{
		nNB = RunParams::nused_potpartials;
		break;
	}
	case 20:
		nNB=1;
		break;
	}
	
	if (RunParams::potential>0)
		p_pot_chitotal=&total2;//all bonded and non-bonded potential related chi2 components (aenet is handled this way too) wil be summed separately
	else if (RunParams::fnc==4)
		p_pot_chitotal=&total;//the bonded contribution will be summed together with the chi2 of the data sets
	//the order of storage in these arrays is the order the data are given in
	//the .dat file, which is reflected by the order in the next line 
	chisize = ngr + nsq + nfq + nfg + nek + ncos + nacc + tot_cc_subconst;//total number of chi squared parts
#ifdef _ADVANCED_GEOM_CONST
	chisize += ncommonneigh + nsecondneigh+nbvs;
#endif
#ifdef _LOCAL_INV
	int i, j,k, padding_offset;
	nlocint=RunParams::nlocint;
	chisize+=nlocint;//for the locale invariance contribution

	if (RunParams::loc_chi2_mode==1)
		CalcLocChi2Thread=&ChiSquared::CalcLocChiSquaredThread_dist;
	else
		CalcLocChi2Thread=&ChiSquared::CalcLocChiSquaredThread_bin;

	SetArraysize(&sum_loc_natoms,nlocint,"sum_loc_natoms","ChiSquared::SetChiSquaredParams");

	if (RunParams::loc_chi2_mode)
	{
		SetArraysize(&loc_natoms_min,SimpleCfg::ntypes*(nlocint+1),"loc_natoms_min","ChiSquared::SetChiSquaredParams");
		SetArraysize(&loc_natoms,SimpleCfg::ntypes,"loc_natoms","ChiSquared::SetChiSquaredParams");
		SetArraysize(&nused_loc,SimpleCfg::ntypes,"nused_loc","ChiSquared::SetChiSquaredParams");
		for (i=0;i<SimpleCfg::ntypes;i++)
		{
			//loc_natoms neighbour atoms beginning with loc_natoms_min will be used from each type

			for (j=0;j<=nlocint;j++)
			{
				//(first all interval of type 1, then type2...)
				loc_natoms_min[i*(nlocint+1)+j]=int(SimpleCfg::pnatoms[i]*RunParams::loc_at_ratio[j]);
				if (loc_natoms_min[i*(nlocint+1)+j]>0)
					loc_natoms_min[i*(nlocint+1)+j]--;//indexing start with 0
				if (j==0)
					loc_natoms[i]=int(SimpleCfg::pnatoms[i]*RunParams::loc_at_ratio[nlocint])-1-loc_natoms_min[i*(nlocint+1)]+1;
			}
		}

		padding_offset=0;//;(RunParams::nthreads>1 ? (int)((CACHE_PADDING-sizeof(*first_loc_ind))/sizeof(*first_loc_ind)) : 0);//number of dummy double elements
		SetArraysize(&first_loc_ind,SimpleCfg::ntotal*SimpleCfg::ntypes*RunParams::nthreads+(RunParams::nthreads-1)*padding_offset,"first_loc_ind","ChiSquared::SetChiSquaredParams");
		SetArraysize(&starting_count,SimpleCfg::ntotal*SimpleCfg::ntypes*RunParams::nthreads+(RunParams::nthreads-1)*padding_offset,"starting_count","ChiSquared::SetChiSquaredParams");
		SetArraysize(&cumalpart_loc_atoms,(int)pow((double)SimpleCfg::ntypes,2)+1,"cumalpart_loc_atoms","ChiSquared::SetChiSquaredParams");

		//Here there has to be separate av_hist for all the different partials (1-2 and 2-1 are not the same)
		for (i=0;i<nlocint;i++)
		{
			sum_loc_natoms[i]=0;
			for (j=0;j<SimpleCfg::ntypes;j++)
			{
				sum_loc_natoms[i]+=loc_natoms_min[j*(nlocint+1)+i+1]-loc_natoms_min[j*(nlocint+1)+i];
				if (i==nlocint-1)
					sum_loc_natoms[i]++;//the last atom will be in the last interval
				
				if (i==0)
					sum_loc_natoms_tot+=loc_natoms[j];
			}
		}

		SetArraysize(&av_loc_hist,sum_loc_natoms_tot*SimpleCfg::ntypes*RunParams::nthreads,"av_loc_hist","ChiSquared::SetChiSquaredParams");
		SetArraysize(&av_loc_hist_offset,(int)pow((double)SimpleCfg::ntypes,2)+1,"av_loc_hist_offset","ChiSquared::SetChiSquaredParams");
		av_loc_hist_offset[0]=0;
		cumalpart_loc_atoms[0]=0;
		k=0;
		for (i=0;i<SimpleCfg::ntypes;i++)
		{
			for (j=0;j<SimpleCfg::ntypes;j++)
			{
				cumalpart_loc_atoms[k+1]=cumalpart_loc_atoms[k]+SimpleCfg::pnatoms[i];
				av_loc_hist_offset[k+1]=av_loc_hist_offset[k]+loc_natoms[j];
				k++;
			}
		}
	}


#endif

	chi_potstart=chisize;//beginning of the first potential related member in the chicomp array
	if (RunParams::potential==1)//there is LJ potential
	{
		chisize+=2*nNB;//vdW+Coulomb (NB is 1, if RunParams::NB_weight_mode=0,2, npartials if 1)
		if (Topology::npairs>0)//there is 1-4 interactions
			chisize+=2*nNB;//vdW14+Coulomb14 contribution
	}
	else if (RunParams::potential == 10)//there is tabulated potential, for which only vdW is used
	{
		chisize += nNB;//vdW (NB is nused_potpartials, RunParams::NB_weight_mode=1)
	}
	else if (RunParams::potential == 20)//there is aenet
	{
		chisize +=1;
	}
	//here the order is first the non-bonded interactions, then the bonded interactions in the order they are given here
	chisize+=nbonds+nangles+nperdihs+nharmdihs+nRBdihs;//for aenet does not matter, these will be 0
	
	if (RunParams::potential>0)
	{
		switch (RunParams::potential)
		{
			case 1:
			{
				RunParams::ntotal_points2 = nbonds + nangles + nperdihs + nharmdihs + nRBdihs;//all the bonded and non-bonded interactions will be counted separately from the data sets
				RunParams::ntotal_points2 += 2 * nNB + 2 * RunParams::tot_used14part;
				break;
			}
			case 10:
			{
				RunParams::ntotal_points2 = nbonds + nangles + nperdihs + nharmdihs + nRBdihs;//all the bonded and non-bonded interactions will be counted separately from the data sets
				RunParams::ntotal_points2+=nNB;
				break;
			}
			case 20:
			{
				RunParams::ntotal_points2=1;
			}
			
		}
		

	}
	else if (RunParams::fnc==4) RunParams::ntotal_points+=nbonds+nangles+nperdihs+nharmdihs+nRBdihs;
};

void ChiSquared::Save(ofstream &file) const//save to a file
{
	double *browser;
	int i;
	
		
	file<<"This is an object of Class ChiSquared "<<endl;
	file << "Number of constraints: g(r), S(Q), F(Q), F(g), E(k), cosine distr. c., individual c.c., average c.c.";
#ifdef _ADVANCED_GEOM_CONST
	file << " common neighbour c., second neighbour c., bond valence sum c.";
#endif
	
#ifdef _LOCAl_INV
	file<<"Local invariance is calulated"<<endl;
#endif
	if (RunParams::potential>0)
		file<<endl;
	else if (RunParams::fnc==4)
		file<<", bond, angle, dihedral: (periodic, harmonic, RB)"<<endl;
	else file<<endl;

	file << ngr << "  " << nsq << "  " << nfq << "  " << nfg << "  " << nek << "  " << ncos << "  " << tot_cc_subconst << "  " << nacc;
#ifdef _ADVANCED_GEOM_CONST
	file << "  " << ncommonneigh << "  " << nsecondneigh << "  " << nbvs;
#endif
	if (RunParams::potential>0)
	{
		file<<endl;
		file<<J10<<total<<" !total chi squared value"<<endl;
		if (RunParams::potential == 1)
		{
			file << "Total number of bonded and non-bonded interactions: vdW, Coulomb, 1-4 vDW, 1-4 Coulomb, bond, angle, dihedral: (periodic, harmonic, RB)" << endl;
			file << nNB << "  " << nNB << "  " << RunParams::tot_used14part << "  " << RunParams::tot_used14part << "  " << nbonds << "  " << nangles << "  " << nperdihs << "  " << nharmdihs << "  " << nRBdihs << endl;
			
		}
		else if (RunParams::potential ==10)
		{
			file << "Total number of tabulated potential interactions";
			if (RunParams::fnc == 4)
				file << " and bonded  interactions: bond, angle, dihedral: (periodic, harmonic, RB)";
			file<< endl;
			file << nNB;
			if (RunParams::fnc == 4)
				file << "  " << nbonds << "  " << nangles << "  " << nperdihs << "  " << nharmdihs << "  " << nRBdihs;
			file<< endl;

		}
		else if (RunParams::potential ==20)
		{
			file<<"Total number of ANN potential interactions"<<endl;
			file<< 1<<endl;
		}
		file << J10 << total2 << " !total bonded and non-bonded, aenet-related chi squared value" << endl;

	}
	else 
	{
		if (RunParams::fnc==4) file<<"  "<<nbonds<<"  "<<nangles<<"  "<<nperdihs<<"  "<<nharmdihs<<"  "<<nRBdihs<<endl;
		else file<<endl;
		file<<J10<<total<<" !total chi squared value"<<endl;
	}
	
	file<<"Chi squared details "<<endl;
	
	browser=chicomp;
	file.precision(4);
	file.setf(ios::floatfield, ios::fixed);
	file.setf(ios::adjustfield, ios::right);
	if(ngr>0)
		file<<"\ng(r) sets"<<endl;
	for(i=0;i<ngr;i++)
	{
		file << *browser;
		if (ExptsData::gruseR[i])
			file << "\t\tRw was used";
		file<< endl;
		browser++;
	}		
	if(nsq>0)
		file<<"\nS(Q) sets"<<endl;
	for(i=0;i<nsq;i++)
	{
		file<<*browser;
		if (ExptsData::squseR[i])
			file << "\t\tRw was used";
		file << endl;
		browser++;
	}		
	if(nfq>0)
		file<<"\nF(Q) sets"<<endl;
	for(i=0;i<nfq;i++)
	{
		file<<*browser;
		if (ExptsData::fquseR[i])
			file << "\t\tRw was used";
		file << endl;
		browser++;
	}
	if (nfg > 0)
		file << "\nF(g) sets" << endl;
	for (i = 0; i < nfg; i++)
	{
		file << *browser;
		if (ExptsData::fguseR[i])
			file << "\t\tRw was used";
		file << endl;
		browser++;
	}
	if(nek>0)
		file<<"\nE(k) sets"<<endl;
	
	for(i=0;i<nek;i++)
	{
		file<<*browser;
		if (ExptsData::ekuseR[i])
			file << "\t\tRw was used";
		file << endl;
		browser++;
	}
	if(ncos>0)
		file<<"\nCosine distr. c."<<endl;
	for(i=1;i<=ncos;i++)
	{
		file<<*browser<<endl;
		browser++;
	}
	if(tot_cc_subconst>0)
		file<<"\nIndividual c.c."<<endl;
	for(i=1;i<=tot_cc_subconst;i++)
	{
		file<<*browser<<endl;
		browser++;
	}	
	if(nacc>0)
		file<<"\nAverage c.c."<<endl;
	for(i=1;i<=nacc;i++)
	{
		file<<*browser<<endl;
		browser++;
	}
#ifdef _ADVANCED_GEOM_CONST
	if (ncommonneigh > 0)
		file << "\nCommon neighbour c." << endl;
	for (i = 1; i <= ncommonneigh; i++)
	{
		file << *browser << endl;
		browser++;
	}
	if (nsecondneigh > 0)
		file << "\nSecond neighbour c." << endl;
	for (i = 1; i <= nsecondneigh; i++)
	{
		file << *browser << endl;
		browser++;
	}
	if (nbvs > 0)
		file << "\nIndividual Bond valence sum c." << endl;
	for (i = 1; i <= nbvs; i++)
	{
		file << *browser << endl;
		browser++;
	}
#endif
#ifdef _LOCAL_INV
	file<<"\nLocal invariance"<<endl;
	file<<*browser<<endl;
	browser++;
#endif

	if (RunParams::potential==1)
	{
		file<<"\nvdW potential"<<endl;
	
		for (i=0;i<nNB;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
		file<<"\nCoulomb potential"<<endl;
		for (i=0;i<nNB;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
		if (Topology::npair_types>0)
		{
			file<<"\nvdW 1-4 potential"<<endl;
			for (i=0;i<nNB;i++)
			{
				file<<*browser<<endl;
				browser++;
			}
			file<<"\nCoulomb 1-4 potential"<<endl;
			for (i=0;i<nNB;i++)
			{
				file<<*browser<<endl;
				browser++;
			}
		}
	}
	else if (RunParams::potential == 10)
	{
		file << "\nTabulated potential" << endl;

		for (i = 0; i < nNB; i++)
		{
			file << *browser << endl;
			browser++;
		}
	}
	else if (RunParams::potential == 20)
	{
		file << "\nANN potential" << endl;
		file << *browser << endl;
			browser++;
	}
	
	if (nbonds>0)
	{
		file<<"\nBond potential"<<endl;
		for (i=0;i<nbonds;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
	}
	if (nangles>0)
	{
		file<<"\nAngle potential"<<endl;
		for (i=0;i<nangles;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
	}
	if (nperdihs>0)
	{
		file<<"\nPeriodic dihedral potential"<<endl;
		for (i=0;i<nperdihs;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
	}
	if (nharmdihs>0)
	{
		file<<"\nHarmonic dihedral potential"<<endl;
		for (i=0;i<nharmdihs;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
	}
	if (nRBdihs>0)
	{
		file<<"\nRB dihedral potential"<<endl;
		for (i=0;i<nRBdihs;i++)
		{
			file<<*browser<<endl;
			browser++;
		}
	}

	//might be additional saving to the file, not closed here
}

//Due to the introduction of background correction anf AXS for I(Q) data sets, which is done in a step without move
//for F(Q) sets it has to be checked here, whether to calculate or not.
//For other data sets if no calculaton is necessary, this routine will not be called
void ChiSquared::Coefficient(int nset, int first_set, int *renorm, int *offset, int *linear, int *quadratic, int *cubic, int *nused, bool checkIQ)
{//Calcule each coefficients

	int i,j;
	bool calc;

	for (i=0;i<nset;i++)//for each data set of this type
	{
		calc = 1;
		if (checkIQ)
		{
			calc = 0;
			if (!Move::nomove || (Move::makemucorr && edata.fqIQbackgcorr[i]) || (Move::makefprimeshift && edata.fqAXS[i]>0))
				calc = 1;
		}
		if (calc == false)
			continue;//go to next set
		j= first_set+i;//used for pointer positioning, j is the index of the set among all the data sets

        switch (fit_index[j])
        {
	       	case 15:
			//Solve the		x4  x3  x2  fx2  sx2
			//				x3  x2  x    fx   sx
			//				x2  x   x0  f    s
			//				fx2  fx  f   ff   fs
				four(*(x4+j), *(x3+j), *(x2+j), *(fx2+j), *(sx2+j), *(x2+j), *(x+j), *(fx+j), *(sx+j), (double) *(nused+i),\
					 *(f+j), *(s+j), *(ff+j), *(fs+j), *(d+j), *(c+j), *(b+j), *(a+j));//renorm+constant+linear+quadratic
				*( e+j)=0.0;
   				break;

			case  0:
				*( a+j)=1.0;//nothing
    			*( b+j)=0.0;
        		*( c+j)=0.0;
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

			case  1:
				*( a+j)= ( *( fs+j) / *( ff+j) );//renorm
				if (!(checkIQ && ExptsData::fqfitIQ[i]))//in case of of F(Q) sets, in case of I(Q) fitting b can have custon defined value, should not be overwritten
					*( b+j)=0.0;
        		*( c+j)=0.0;
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;
        
        	case  2:
				*( b+j)= ( *( s+j) - *( f+j) ) / ((double) *(nused+i));//constant
    			*( a+j)=1.0;
        		*( c+j)=0.0;
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

			case  3:
				two((double) *(nused+i),*( f+j),*( s+j),*( ff+j),*( fs+j),*( b+j),*( a+j));//renorm+constant
   				*( c+j)=0.0;
      			*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case  4:
				*( c+j)= ( *( sx+j) - *( fx+j) ) / ( *( x2+j) );//linear
    			*( a+j)=1.0;
        		*( b+j)=0.0;
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case  5:
				two(*( x2+j),*( fx+j),*( sx+j),*( ff+j),*( fs+j),*( c+j),*( a+j));//renorm+linear
    			*( b+j)=0.0;
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case  6:
				two(*( x2+j),*( x+j),*( sx+j) - *( fx+j),(double) *(nused+i),*( s+j) - *( f+j),\
    				*( c+j),*( b+j));//constant+linear
        		*( a+j)=1.0;
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case  7:
				three(*( x2+j),*( x+j),*( fx+j),*( sx+j),(double) *(nused+i),*( f+j),*( s+j),\
    				*( ff+j),*( fs+j),*( c+j),*( b+j),*( a+j));//renorm+constant+linear
        		*( d+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case  8:
				*( d+j)=( *( sx2+j) - *( fx2+j) ) / ( *( x4+j) );//quadratic
        		*( a+j)=1.0;
        		*( b+j)=0.0;
        		*( c+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case  9:
				two(*( x4+j),*( fx2+j),*( sx2+j),*( ff+j),*( fs+j),*( d+j),*( a+j));//renorm+quadratic
        		*( b+j)=0.0;
        		*( c+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case 10:
				two(*( x4+j),*( x2+j),*( sx2+j) - *( fx2+j),(double) *(nused+i),*( s+j) - *( f+j),\
        			*( d+j),*( b+j));//constant+quadratic
				*( a+j)=1.0;
        		*( c+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case 11:
				three(*( x4+j),*( x2+j),*( fx2+j),*( sx2+j),(double) *(nused+i),*( f+j),*( s+j),\
        			*( ff+j),*( fs+j),*( d+j),*( b+j),*( a+j));//renorm+constant+quadratic
        		*( c+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case 12:
				two(*( x4+j),*( x3+j),*( sx2+j) - *( fx2+j),*( x2+j),*( sx+j) - *( fx+j),\
        		*( d+j),*( c+j));//linear+quadratic
        		*( a+j)=1.0;
        		*( b+j)=0.0;
				*( e+j)=0.0;
        		break;

        	case 13:
				three(*( x4+j),*( x3+j),*( fx2+j),*( sx2+j),*( x2+j),*( fx+j),*( sx+j),\
        			*( ff+j),*( fs+j),*( d+j),*( c+j),*( a+j));//renorm+linear+quadratic
        		*( b+j)=0.0;
				*( e+j)=0.0;
        		break;

			case 14:
				three(*( x4+j),*( x3+j),*( x2+j),*( sx2+j)-*( fx2+j),*( x2+j),*( x+j),\
					*( sx+j)-*( fx+j),(double) *(nused+i),*( s+j) - *( f+j),\
				*( d+j),*( c+j),*( b+j));//constant+linear+quadratic
        		*( a+j)=1.0;
				*( e+j)=0.0;
        		break;

			case 16:
				*( a+j)=1.0;
				*( b+j)=0.0;
        		*( c+j)=0.0;
        		*( d+j)=0.0;
				*( e+j)=( *( sx3+j) - *( fx3+j))/ *(x6+j); 
				break;

			case 17:
				two(*( x6+j),*( fx3+j),*( sx3+j),*( ff+j),*( fs+j),	*( e+j),*( a+j));//renorm, cubic
				*( b+j)=0.0;
        		*( c+j)=0.0;
        		*( d+j)=0.0;
				break;
			

			case 18:
				two(*( x6+j),*( x3+j),*( sx3+j) - *( fx3+j),(double)*( nused+i), *(s+j) - *(f+j),*( e+j),*( b+j));//constant, cubic
				*( a+j)=1.0;
        		*( c+j)=0.0;
        		*( d+j)=0.0;
				break;

			case 19:
				three(*( x6+j),*( x3+j),*( fx3+j),*( sx3+j),(double)*( nused+i),*( f+j),*( s+j),\
					*( ff+j),*( fs+j),*( e+j),*( b+j),*( a+j));//renorm+constant+cubic
        		*( c+j)=0.0;
				*( d+j)=0.0;
        		break;
			
			case 20:
				two(*( x6+j),*( x4+j),*( sx3+j) - *( fx3+j),*( x2+j), *(sx+j) - *(fx+j),*( e+j),*( c+j));//linear, cubic
				*( a+j)=1.0;
        		*( b+j)=0.0;
        		*( d+j)=0.0;
				break;

			case 21:
				three(*( x6+j),*( x4+j),*( fx3+j),*( sx3+j),*( x2+j),*( fx+j),*( sx+j),\
					*( ff+j),*( fs+j),*( e+j),*( c+j),*( a+j));//renorm+linear+cubic
        		*( b+j)=0.0;
				*( d+j)=0.0;
        		break;

			case 22:
				three(*( x6+j),*( x4+j),*( x3+j),*( sx3+j)- *(fx3+j),*( x2+j),*( x+j),*( sx+j)- *(fx+j),\
					(double)*( nused+i),*( s+j)- *( f+j),*( e+j),*( c+j),*( b+j));//constant+linear+cubic
        		*( a+j)=1.0;
				*( d+j)=0.0;
        		break;

			case 23:
				four(*(x6+j), *(x4+j), *(x3+j), *(fx3+j), *(sx3+j), *(x2+j), *(x+j), *(fx+j), *(sx+j), (double) *(nused+i),\
					 *(f+j), *(s+j), *(ff+j), *(fs+j), *(e+j), *(c+j), *(b+j), *(a+j));//renorm+constant+linear+cubic
				*( d+j)=0.0;
   				break;
			case 24:
				two(*( x6+j),*( x5+j),*( sx3+j) - *( fx3+j),*( x4+j), *(sx2+j) - *(fx2+j),*( e+j),*( d+j));//quadratic, cubic
				*( a+j)=1.0;
        		*( b+j)=0.0;
        		*( c+j)=0.0;
				break;

			case 25:
				three(*( x6+j),*( x5+j),*( fx3+j),*( sx3+j), *(x4+j),*( fx2+j),*( sx2+j),*( ff+j),\
					*( fs+j),*( e+j),*( d+j),*( a+j));//renorm+quadratic+cubic
        		*( b+j)=0.0;
				*( c+j)=0.0;
        		break;

			case 26:
				three(*( x6+j),*( x5+j),*( x3+j),*( sx3+j) - * (fx3+j), *(x4+j),*(x2+j),*(sx2+j) - *(fx2+j),(double)*(nused+i),\
					*( s+j) - *(f+j),*( e+j),*( d+j),*( b+j));//constant+quadratic+cubic
        		*( a+j)=1.0;
				*( c+j)=0.0;
        		break;

			case 27:
				four(*(x6+j), *(x5+j), *(x3+j), *(fx3+j), *(sx3+j), *(x4+j), *(x2+j), *(fx2+j), *(sx2+j), (double) *(nused+i),\
					 *(f+j), *(s+j), *(ff+j), *(fs+j), *(e+j), *(d+j), *(b+j), *(a+j));//renorm+constant+quadratic+cubic
				*( c+j)=0.0;
   				break;

			case 28:
				three(*( x6+j),*( x5+j),*( x4+j),*( sx3+j) - * (fx3+j), *(x4+j),*(x3+j),*(sx2+j) - *(fx2+j),*(x2+j),\
					*( sx+j) - *(fx+j),*( e+j),*( d+j),*( c+j));//linear+quadratic+cubic
        		*( a+j)=1.0;
				*( b+j)=0.0;
        		break;


			case 29:
				four(*(x6+j), *(x5+j), *(x4+j), *(fx3+j), *(sx3+j), *(x4+j), *(x3+j), *(fx2+j), *(sx2+j), *(x2+i),\
					 *(fx+j), *(sx+j), *(ff+j), *(fs+j), *(e+j), *(d+j), *(c+j), *(a+j));//renorm+linear+quadratic+cubic
				*( b+j)=0.0;
   				break;

			case 30:
				four(*(x6+j), *(x5+j), *(x4+j), *(x3+j), *(sx3+j) - *(fx3+j), *(x4+j), *(x3+j), *(x2+j), *(sx2+j)- *(fx2+j),\
					*(x2+j), *(x+j), *(sx+j) - *(fx+j), (double) *(nused+i), *(s+j) - *(f+j), *(e+j), *(d+j), *(c+j), *(b+j));//constant+linear+quadratic+cubic
				*( a+j)=1.0;
   				break;

			case 31:
				five(*(x6+j), *(x5+j), *(x4+j), *(x3+j), *(fx3+j), *(sx3+j), *(x4+j), *(x3+j), *(x2+j), *(fx2+j), *(sx2+j),\
					*(x2+j), *(x+j), *(fx+j), *(sx+j), (double) *(nused+i), *(f+j), *(s+j), *(ff+j), *(fs+j),\
					*(e+j), *(d+j), *(c+j), *(b+j), *(a+j));//renorm+constant+linear+quadratic+cubic
				break;
			case 32://nonlin, will be calculated later by NonLinReg
				//a, b custom, alpha will be determined
				
			case 33://nonlin, will be calculated later by NonLinReg
				//b=0, a, alpha will be determined
				
			case 35://nonlin, will be calculated later by NonLinReg
				//a, b, alpha will be determined
				*(c + j) = 0.0;
				*(d + j) = 0.0;
				*(e + j) = 0.0;
				break;
			default:
				//this can only happen in case of Xray-I(Q) fit, where not all the combinations implemented
				MissingRenorCombError(i);

		}
	}
};

// The out1, out2, out3... variables will be in the order of a3, a2, a1, a0, alpha (e, d, c, b, a) or any subset of them
void ChiSquared::two(double A, double B, double C, double D, double E, double &out1, double &out2)
{
	//Calculate two coefficients
	//Solve the   A B C
	//	          B D E type matrix
	
	out2=( A * E - B * C ) / ( A * D - B * B ); //y
	out1=( C - out2 * B ) / A; //x
};

void ChiSquared::three(double A, double B, double C, double D, double E, double F, double G, double H, double I, double &out1, double &out2, double &out3)
{
	//Calculate three coefficients
	//Solve the   A B C D
	//		      B E F G
	//		      C F H I type matrix

	double a = A * E - B * B;
	double b = A * F - B * C;
	double c = A * G - B * D;
	double d = A * H - C * C;
	double e = A * I - C * D;

	two(a,b,c,d,e,out2,out3);

	out1 = ( D - C * out3 - B * out2 ) / A; //x
};

void ChiSquared::four(double j, double k, double l, double m, double n, double o, double p, double q, double r, double s, double t,\
					  double u, double v, double w, double &out1, double &out2, double &out3, double &out4)
{
	//Calculate four coefficients 
	//Solve the j	k	l	m	n
	//			k	o	p	q	r
	//			l	p	s	t	u
	//			m	q	t	v	w matrix

	double A = j * o - k * k;
	double B = j * p - k * l;
	double C = j * q - k * m;
	double D = j * r - k * n;
	double E = j * s - l * l;
	double F = j * t - l * m;
	double G = j * u - l * n;
	double H = j * v - m * m;
	double I = j * w - m * n;

	three(A,B,C,D,E,F,G,H,I,out2,out3,out4);
	out1 = ( n - out4 * m - out3 * l - out2 * k )/ ( j );

}; 

void ChiSquared::five(double A5, double B5, double C5, double D5, double E5, double F5, double G5, double H5,\
					  double I5, double J5, double K5, double L5, double M5, double N5, double O5, double P5,\
					  double Q5, double R5, double S5, double T5, double &out1, double &out2, double &out3, \
					  double &out4, double &out5)
{
	

	//Calculate all coefficients in case of five coeff
	//Solve the		A5	B5	C5	D5	E5	F5
	//				B5	G5	H5	I5	J5	K5
	//				C5	H5	L5	M5	N5	O5
	//				D5	I5	M5	P5	Q5	R5
	//				E5	J5	N5	Q5	S5	T5


	double j = A5 * G5 - B5 * B5;
	double k = A5 * H5 - B5 * C5;
	double l = A5 * I5 - B5 * D5;
	double m = A5 * J5 - B5 * E5;
	double n = A5 * K5 - B5 * F5;
	double o = A5 * L5 - C5 * C5;
	double p = A5 * M5 - C5 * D5;
	double q = A5 * N5 - C5 * E5;
	double r = A5 * O5 - C5 * F5;
	double s = A5 * P5 - D5 * D5;
	double t = A5 * Q5 - D5 * E5;
	double u = A5 * R5 - D5 * F5;
	double v = A5 * S5 - E5 * E5;
	double w = A5 * T5 - E5 * F5;


	four(j,k,l,m,n,o,p,q,r,s,t,u,v,w,out2,out3,out4,out5);

	out1 = ( F5 - out5 * E5 - out4 * D5 - out3 * C5 - out2 * B5)/ (A5) ;
};

//calculate the chi squared
#ifdef _ADVANCED_GEOM_CONST
	longint ChiSquared::CalcChiSquared(CalcData &calc, CoordNumbConst &icc, AvCoordConst &avcc, CommonNeighConst &conc, SecondNeighConst &snc, BondValenceSumConst &bvs)
#else
	longint ChiSquared::CalcChiSquared(CalcData &calc, CoordNumbConst &icc, AvCoordConst &avcc)
#endif
{
	
		int i,j,k,index,isubconst;
		longint calc_method;
		int non_lin_failure;
		int *pi1, *pi2; //pointers used in C.C. add-on
		double diff,temp,background, exprenorm;
		double *pexp, *pcalc, *pq, *pg, *pr;//pointers to data (expt and calculated, q,r,compton,coeffs) values
		double *browser;//pointer to the chicomp array
		double *pd1,*pd2;//used in C.C. add-on
		bool exprenorm_calculated;//	whether the background was already calculated (only needed if R factor is used instead of the normal chisquare
#ifdef _TEST_MODE	
#ifdef _LOCAL_INV	
		std::chrono::duration<double, std::milli> elapsed;//for time diff calc
#endif
#endif
		calc_method=0;//default
	
		if (!Move::nomove)//has to calculated at least once, even if there is no atomic movement at all
		{

			//first compute the sums for each data set
			//the g(r) components, array s is not needed!
			pexp = edata.gr_gvalues;//pointer initialisation
			pcalc = calc.grvalues;
			for (i = 0; i < ngr; i++)//for each g(r) data set
			{
				ss[i] = 0.0;//sum of calculated*calculated data
				fs[i] = 0.0;//sum of calculated*experimental data
				for (j = 0; j < *(edata.grused + i); j++)//for all data points
				{

					ss[i] += (*pcalc * *pcalc);
					fs[i] += (*pcalc * *pexp);
					pexp++;
					pcalc++;
				}
			}

			//the r-dependent terms are only calculated, if cubic correction is used, as normally, the g(r) can only be normalized by a multiplier
			for (i = 0; i < ngr; i++)//for each g(r) data set
			{

				if (ExptsData::use_cubic[i])
				{
					pcalc = calc.grfinder[i];//pointer initialisation
					pr = edata.gr_rfinder[i];
					s[i] = 0.0;//sum of calculated data
					sx[i] = 0.0;//sums of calculated*r
					sx2[i] = 0.0;//sums of calculated*r*r
					sx3[i] = 0.0;//sums of calculated*r*r*r
					for (j = 0; j < *(edata.grused + i); j++)//for all data points
					{
						s[i] += *pcalc;//sum of calculated data
						sx[i] += (*pcalc * *pr);
						sx2[i] += (*pcalc * *pr * *pr);
						sx3[i] += (*pcalc * *pr * *pr * *pr);
						pcalc++;
						pr++;
					}
				}
			}

			//the S(Q) components
			pexp = edata.sq_svalues;//pointer initialisation
			pcalc = calc.sqvalues;
			pq = edata.sq_qvalues;
			for (i = 0; i < nsq; i++)//for each s(q) data set
			{
				index = ngr + i;
				s[index] = 0.0;//sum of calculated data
				ss[index] = 0.0;//sum of calculated*calculated data
				fs[index] = 0.0;//sum of calculated*experimental data
				sx[index] = 0.0;//sums of calculated*Q
				sx2[index] = 0.0;//sums of calculated*Q*Q
				for (j = 0; j < *(edata.sqused + i); j++)//for all data points
				{
					s[index] += *pcalc;//sum of calculated data
					ss[index] += (*pcalc * *pcalc);
					fs[index] += (*pcalc * *pexp);
					sx[index] += (*pcalc * *pq);
					sx2[index] += (*pcalc * *pq * *pq);
					pexp++;
					pcalc++;
					pq++;
				}
			}
			//the cubic terms will only be calculated, if they are needed
			for (i = 0; i < nsq; i++)//for each s(q) data set
			{
				index = ngr + i;
				if (ExptsData::use_cubic[index])
				{
					pexp = edata.sq_sfinder[i];//pointer initialisation
					pcalc = calc.sqfinder[i];
					pq = edata.sq_qfinder[i];
					sx3[index] = 0.0;
					for (j = 0; j < *(edata.sqused + i); j++)//for all data points
					{
						sx3[index] += (*pcalc * *pq * *pq * *pq);
						pcalc++;
						pq++;
					}
				}
			}

		}

	
		//the F(Q) components
		for (i = 0; i < nfq; i++)//for each F(Q) data set
		{
			index = ngr + nsq + i;
			//only calculate, if nothing is prevented the atomic move, if there is no move only calculate for sets with background corr or fprimeshift
			if (!Move::nomove || (Move::makemucorr && edata.fqIQbackgcorr[i]) || (Move::makefprimeshift && edata.fqAXS[i] > 0))
			{ 
				if (ExptsData::fqfitIQ[i])//I(Q) fitting
				{
					pcalc = calc.iqfinder[i];
					pexp = edata.fq_ffinder[i];//pointer initialisation
					//pq = edata.fq_qfinder[i];
		
					//if mucorr is used with no alpha renorm, then only the experimental is changing the calc not
					if (!(Move::makemucorr && ExptsData::fqIQbackgcorr[i]))
					{
						switch (fit_index[index])
						{
						case 1://linear:a optimized with custom b and alpha, s=I(Q)_calc-b
						{
							s[index] = 0.0;//sum of calculated*calculated data, needed, as b can be non-zero
							ss[index] = 0.0;//sum of calculated*calculated data
							fs[index] = 0.0;//sum of calculated*experimental data
							for (j = 0; j < *(edata.fqused + i); j++)//for all data points
							{//s[index] is not needed, will not be calculted
								s[index] += *pcalc - b[index];
								ss[index] += (pow(*pcalc - b[index], 2));
								fs[index] += ((*pcalc - b[index]) * *pexp);
								pexp++;
								pcalc++;
							}
							break;
						}
						case 3://lnear:here 'a' and 'b' is fitted, but I(Q) should be used instead of F(Q)
							s[index] = 0.0;//sum of calculated data
							ss[index] = 0.0;//sum of calculated*calculated data
							fs[index] = 0.0;//sum of calculated*experimental data
							for (j = 0; j < *(edata.fqused + i); j++)//for all data points
							{
								s[index] += *pcalc;
								ss[index] += (*pcalc * *pcalc);
								fs[index] += (*pcalc * *pexp);
								pexp++;
								pcalc++;
							}
							break;
						//if alpha renorm is used, calulated data changes during iteration, those arrays containing s have to be calculated during
							//iteration, here only f containing but no alpha containing arrays are calculated
						case 32://nonlinear alpha fit, a,b custom; the I(Q)-containing sums have to be recalculated during the iteration, not calculated here
						case 33://nonlinear a, alpha fit, b=0
						case 35://nonlinear a,b, alpha fit, 
							break;
						default:
							MissingRenorCombError(i);
						}
					}
					else//mu correction, only expt-containing arrays change if alpha is constant
					{
						//the actual experimental data is changed, recalculate the arrays depending on it
						switch (fit_index[index])
						{
						case 1://linear:a is optimized with custom b and alpha, s=I(Q)_calc-b
							pexp = edata.fq_ffinder[i];
							pcalc = calc.iqfinder[i];
							f[index] = 0;
							ff[index] = 0;
							fs[index] = 0.0;//sum of calculated*experimental data
							//fx[index] = 0;
							//fx2[index] = 0;
							for (j = 0; j < *(edata.fqused + i); j++)//for all data points
							{
								f[index] += *pexp;
								ff[index] += *pexp * (*pexp);
								fs[index] += (*pcalc-b[index]) * *pexp;
								pcalc++;
								pexp++;
							}
							break;
						case 3://linear:here 'a' and 'b' is fitted, but I(Q) should be used instead of F(Q)
							//all is needed for chi2 calc
							pexp = edata.fq_ffinder[i];
							pcalc = calc.iqfinder[i];
							f[index] = 0;
							ff[index] = 0;
							fs[index] = 0.0;//sum of calculated*experimental data
							//fx[index] = 0;
							//fx2[index] = 0;
							for (j = 0; j < *(edata.fqused + i); j++)//for all data points
							{
								f[index] += *pexp;
								ff[index] += *pexp * (*pexp);
								fs[index] += (*pcalc * *pexp);
								//	fx[index] += (*pexp * *pq);these are presently not used, as no lin or quadratic or cubic
								//	fx2[index] += (*pexp * *pq * *pq);
								pcalc++;
								pexp++;
								//pq++;
							}
							break;
						//the I(Q)-containing sums and fC have to be recalculated during the iteration, not calculated here
						case 35://nonlinear a,b, alpha fit
						case 33://nonlinear a, alpha fit, b=0
							pexp = edata.fq_ffinder[i];
							f[index] = 0;//the only exp containing sums do not change during lin reg, should be calculated here
							ff[index] = 0;
							for (j = 0; j < *(edata.fqused + i); j++)//for all data points
							{
								ff[index] += *pexp * (*pexp);
								f[index] += *pexp++;
							}
							break;
						case 32://nonlinear alpha fit, a,b custom;
							break;
						default:
							MissingRenorCombError(i);

						}//end switch
					}//end mucorr
				}
				else
				{
					s[index] = 0.0;//sum of calculated data
					ss[index] = 0.0;//sum of calculated*calculated data
					fs[index] = 0.0;//sum of calculated*experimental data
					sx[index] = 0.0;//sums of calculated*Q
					sx2[index] = 0.0;//sums of calculated*Q*Q
					//normal fit, pointer initialisation
					pcalc = calc.fqfinder[i];
					pexp = edata.fq_ffinder[i];
					pq = edata.fq_qfinder[i];
					for (j = 0; j < *(edata.fqused + i); j++)//for all data points
					{
						s[index] += *pcalc;
						ss[index] += (*pcalc * *pcalc);
						fs[index] += (*pcalc * *pexp);
						sx[index] += (*pcalc * *pq);
						sx2[index] += (*pcalc * *pq * *pq);
						pexp++;
						pcalc++;
						pq++;
					}
				}
			}//end of only calculate if there is change in the data
		
		}
	
		if (!Move::nomove)//only calculate, if nothing is prevented the atomic move
		{
			//the cubic terms will only be calculated, if they are needed, no cubic term is allowed for I(Q) data
			for (i = 0; i < nfq; i++)//for each F(q) data set
			{
				index = ngr + nsq + i;

				if (ExptsData::use_cubic[index])
				{
					pexp = edata.fq_ffinder[i];//pointer initialisation
					pcalc = calc.fqfinder[i];
					pq = edata.fq_qfinder[i];
					sx3[index] = 0.0;
					for (j = 0; j < *(edata.fqused + i); j++)//for all data points
					{
						sx3[index] += (*pcalc * *pq * *pq * *pq);
						pcalc++;
						pq++;
					}
				}
			}
	
		
			//the F(g) components
			pexp = edata.fg_fvalues;//pointer initialisation
			pcalc = calc.fgvalues;
			pg = edata.fg_gvalues;
			for (i = 0; i < nfg; i++)//for each F(g) data set
			{
				index = ngr + nsq + nfq + i;
				s[index] = 0.0;//sum of calculated data
				ss[index] = 0.0;//sum of calculated*calculated data
				fs[index] = 0.0;//sum of calculated*experimental data
				sx[index] = 0.0;//sums of calculated*g
				sx2[index] = 0.0;//sums of calculated*g*g
				for (j = 0; j < *(edata.fgused + i); j++)//for all data points
				{
					s[index] += *pcalc;
					ss[index] += (*pcalc * *pcalc);
					fs[index] += (*pcalc * *pexp);
					sx[index] += (*pcalc * *pg);
					sx2[index] += (*pcalc * *pg * *pg);
					pexp++;
					pcalc++;
					pg++;
				}
			}

			//the cubic terms will only be calculated, if they are needed
			for (i = 0; i < nfg; i++)//for each F(g) data set
			{
				index = ngr + nsq + nfq + i;
				if (ExptsData::use_cubic[index])
				{
					pexp = edata.fg_ffinder[i];//pointer initialisation
					pcalc = calc.fgfinder[i];
					pq = edata.fg_gfinder[i];
					sx3[index] = 0.0;
					for (j = 0; j < *(edata.fgused + i); j++)//for all data points
					{
						sx3[index] += (*pcalc * *pg * *pg * *pg);
						pcalc++;
						pg++;
					}
				}
			}
		}
		
		//the E(k) components
		pexp = edata.ek_evalues;//pointer initialisation
		pcalc = calc.ekvalues;
		for (i = 0; i < nek; i++)//for each E(k) data set
		{
			if (!Move::nomove || (Move::makeE0shift && edata.ek_ngrid_in[i] > 0))
			{//calculate if there is no E0 shift in this step, or if there is only for those, where there is
			//resets the value to zero(needed in the main loop)
				index = ngr + nsq + nfq + nfg + i;
				s[index] = 0.0;//sum of calculated data
				ss[index] = 0.0;//sum of calculated*calculated data
				fs[index] = 0.0;//sum of calculated*experimental data
				for (j = 0; j < *(edata.ekused + i); j++)//for all data points
				{
					s[index] += *pcalc;
					ss[index] += (*pcalc * *pcalc);
					fs[index] += (*pcalc * *pexp);
					pexp++;
					pcalc++;
				}
			}
		}

		if (!Move::nomove)
		{
			//the CosDistrConst components, array s is not needed!
			pexp = CosDistrConst::theordistr;//pointer initialisation for theoretical
			pcalc = cosconst.cosinedistr;//pointer initialisation for calculated
			for (i = 0; i < ncos; i++)//for each CosDistrConst data set
			{
				index = ngr + nsq + nfq + nfg + nek + i;
				ss[index] = 0.0;//sum of calculated data
				fs[index] = 0.0;//sum of calculated*theoretical
				pexp += CosDistrConst::first_bin[i];//set it to the first meaningful point
				pcalc += CosDistrConst::first_bin[i];//set it to the first meaningful point
				if (CosDistrConst::method[i] == 2)//negative constraint
				{
					fs[index] = 0.0;//The desired theoretical values are zero for each point 
					for (j = CosDistrConst::first_bin[i]; j < CosDistrConst::last_bin[i]; j++)//for all used data points
					{
						ss[index] += (*pcalc * *pcalc) * (*pexp * *pexp);//the calculated data is weigthed by the cos_theta-dependent weighting factor
						pexp++;
						pcalc++;
					}
				}
				else
				{
					for (j = CosDistrConst::first_bin[i]; j < CosDistrConst::last_bin[i]; j++)//for all used data points
					{

						ss[index] += (*pcalc * *pcalc);
						fs[index] += (*pcalc * *pexp);
						pexp++;
						pcalc++;
					}
				}
				pexp += CosDistrConst::npoints[i] - CosDistrConst::last_bin[i];//to skip the unused points of the this constraint
				pcalc += CosDistrConst::npoints[i] - CosDistrConst::last_bin[i];//to skip the unused points of the this constraint
			}

			//now calculating the renormalization coeffs
			//for the g(r) data
			Coefficient(ngr, 0, edata.grrenorm, edata.groffset, edata.grlinear, edata.grquadratic, edata.grcubic, edata.grused, false);//Calculate the g(r) coeffs

			//the S(Q) components
			Coefficient(nsq, ngr, edata.sqrenorm, edata.sqoffset, edata.sqlinear, edata.sqquadratic, edata.sqcubic, edata.sqused, false);//Calculate the S(Q) coeffs

		}
	
	
		//the F(Q) components, it will be decided inside Coefficient, whether calculate or not for a given data set
		Coefficient(nfq, ngr + nsq, edata.fqrenorm, edata.fqoffset, edata.fqlinear, edata.fqquadratic, edata.fqcubic, edata.fqused, true);//Calculate the F(Q) coeffs
		if (ExptsData::is_IQ == 2)
		{
			non_lin_failure = NonLinReg(calc);
			if (non_lin_failure)
			{
				failed_nlr++;//increase the number of failed non-lin regressions since the last good one
				calc_method = -1;
				return calc_method;//no need to calculate further, this move cannot be accepeted 
			}
			else
				failed_nlr = 0;//reset
		
		}
	
		if (RunParams::continuation == -1)
				RunParams::continuation = 1;//set it back to show it is not the initial calculation
	

		//the F(g) components
		if (!Move::nomove)
			Coefficient(nfg, ngr + nsq + nfq, edata.fgrenorm, edata.fgoffset, edata.fglinear, edata.fgquadratic, edata.fgcubic, edata.fgused, false);//Calculate the F(g) coeffs
	
		//the E(k) component
		//as E(k) is calculated only for a subset of the histogram points, it can happen, that there are no counts at all
		//in the subset, if the configuration is far from the desired. (Theroretically it can happen for the other type of data sets
		//if the histogram is calculated only for a small range, but not likely, and will not be dealt with.
		//care has to be taken, that no renormalization should be attempted in this case
		for(i=0;i<nek;i++)//for each E(k) data
		{
			if (!Move::nomove || (Move::makeE0shift && edata.ek_ngrid_in[i] > 0))
			{//calculate if there is no E0 shift in this step, or if there is only for those, where there is
			//resets the value to zero(needed in the main loop)
				index = ngr + nsq + nfq + nfg + i;//used for pointer positioning
				if (*(edata.ekrenorm + i) == 1)//renormalization allowed
				{
					if (*(edata.ekoffset + i) == 1)//offset allowed
					{
						*(a + index) = (*(edata.ekused + i) * *(fs + index) - *(s + index) * *(f + index)) / \
							(*(edata.ekused + i) * *(ff + index) - *(f + index) * *(f + index));
						*(b + index) = (*(s + index) - *(a + index) * *(f + index)) / (*(edata.ekused + i));
					}
					else//offset not allowed
					{
						*(a + index) = *(fs + index) / (*(ff + index));
						*(b + index) = 0.0;
					}
					if (*(ff + index) == 0)//all the calculated data points are zero
						*(a + index) = 1.0;//no renormalization
				}
				else//renormalization not allowed
				{
					if (*(edata.ekoffset + i) == 1)//offset allowed
					{
						*(a + index) = 1.0;
						*(b + index) = (*(s + index) - *(f + index)) / (*(edata.ekused + i));
					}
					else//data must  not be adjusted
					{
						*(a + index) = 1.0;
						*(b + index) = 0.0;
					}
				}
			}
		}	

		//values of the chi squared components (one for each data set)
		browser=chicomp;//pointer initialisation
		if (!Move::nomove)
		{
			//the g(r) components
			index = 0;
			for (i = 0; i < ngr; i++)//for each g(r) data set
			{
				exprenorm_calculated = 0;//only needed in case of Rw used instead of normal chisquare
				if (ExptsData::use_cubic[i])
				{
					//the developed sum
					*browser = *(ss + i) + *(a + i) * *(a + i) * *(ff + i) + \
						* (edata.grused + i) * *(b + i) * *(b + i) + \
						* (x2 + i) * (*(c + i) * *(c + i) + 2. * *(b + i) * *(d + i)) + \
						* (d + i) * *(d + i) * *(x4 + i) + \
						2. * (*(a + i) * (*(b + i) * *(f + i) + *(c + i) * *(fx + i) + \
							* (d + i) * *(fx2 + i) - *(fs + i)) + \
							* (b + i) * (*(c + i) * *(x + i) - *(s + i)) + \
							* (c + i) * (*(d + i) * *(x3 + i) - *(sx + i)) - \
							* (d + i) * *(sx2 + i) + \
							//here begins the cubic correction part
							* (e + i) * (*(a + i) * *(fx3 + i) + *(b + i) * *(x3 + i) + *(c + i) * *(x4 + i) + \
								* (d + i) * *(x5 + i) - *(sx3 + i))) + *(e + i) * *(e + i) * *(x6 + i);
				}
				else
					//there is no cubic term, no need for the lengthy calculation of the developed sum, only renormalization is used
					*browser = *(ss + i) + *(a + i) * *(a + i) * *(ff + i) - 2 * *(a + i) * *(fs + i);
				if (*browser <= 0)
				{
					//The chi_square is close to zero, and it cannot be calculated according to the usual
					//way because of the rounding errors, it has to be recalculated

					//for every data point
					pexp = edata.gr_gfinder[i];//sets the pointer to the beginning of g(r) data for i-th expt data series
					pcalc = calc.grfinder[i];//sets the pointer to the beginning of the calculated g(r) data for i-th data series
					pr = edata.gr_rfinder[i];//sets the pointer to the beginning of the r values of the i-th exp g(r) data set
					diff = 0;
					exprenorm = 0;
					if (ExptsData::use_cubic[i])
					{
						for (j = 0; j < edata.grused[i]; j++)
						{
							background = (b[i] + c[i] * *pr + d[i] * *pr * *pr + e[i] * *pr * *pr * *pr);
							exprenorm += pow(a[i] * *pexp + background, 2);
							diff += pow(*pcalc - (a[i] * *pexp + background), 2);
							pcalc++;
							pexp++;
							pr++;
						}
						exprenorm_calculated = 1;
					}
					else
					{
						//only renormalization coeff
						for (j = 0; j < edata.grused[i]; j++)
						{
							diff += pow(*pcalc - a[i] * *pexp, 2);
							pcalc++;
							pexp++;

						}
					}

					*browser = diff;
					//Calc_method is an integer indicator,if the iexpt-th bit is 1, indicates, that for
					//this data series the chi2 has to be calculated using the square of the diff method 
					//if calc_method=0, chi2 could be calculated normally for each data series, if it is
					// greater than zero, than then the data series(es) denoted by non-zero bits were calculated
					//as the square of the difference, and the value of the renormalization coeff a, if 
					//it was calculated at all, is questionable
					//in case of nonlinear fit for a data set 0 means convergence, 1 means no convergence
					//
					if (i == 0)
						calc_method++;
					else
					{
						k = 1;
						for (j = 0; j < i; j++)
							k *= 2;
						calc_method += k;
					}
				}
				if (ExptsData::gruseR[i])
				{
					if (!exprenorm_calculated)//it was not yet calculated
					{
						//for every data point
						pexp = edata.gr_gfinder[i];//sets the pointer to the beginning of g(r) data for i-th expt data series
						pr = edata.gr_rfinder[i];//sets the pointer to the beginning of the r values of the i-th exp g(r) data set
						exprenorm = 0;

						if (ExptsData::use_cubic[i])
						{
							for (j = 0; j < edata.grused[i]; j++)
							{
								background = (b[i] + c[i] * *pr + d[i] * *pr * *pr + e[i] * *pr * *pr * *pr);
								exprenorm += pow(a[i] * *pexp + background, 2);
								pexp++;
								pr++;
							}
						}
						else
						{
							//only renormalization coeff
							for (j = 0; j < edata.grused[i]; j++)
							{
								exprenorm += pow(a[i] * *pexp, 2);
								pcalc++;
								pexp++;

							}
						}
					}
					if (fabs(exprenorm) > TOLERANCE)
					{
						Rw[i] = sqrt(*browser / exprenorm);
						*browser = Rw[i] / *(edata.grsigma + i);//calculating Rw and dividing by sigma 
						Rw[i] *= 100;//this is what was outputed sofar
					}
					else
						Rw[i] = -1000;
				}
				else
				{
					//dividing by the standard dev.squared
					*browser /= *(edata.grsigma + i) * *(edata.grsigma + i);//calculating Rw and dividing by sigma 
				}
				total += *browser;//updating the global result
				if (calc_sigma && RunParams::lead_series_ind == i)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "g(r) ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}


				browser++;
			}
			index = ngr;

			//the S(Q) components
			for (i = 0; i < nsq; i++)//for each S(Q) data set
			{
				//the developed sum
				exprenorm_calculated = 0;
				*browser = *(ss + index) + *(a + index) * *(a + index) * *(ff + index) + \
					* (edata.sqused + i) * *(b + index) * *(b + index) + \
					* (x2 + index) * (*(c + index) * *(c + index) + 2. * *(b + index) * *(d + index)) + \
					* (d + index) * *(d + index) * *(x4 + index) + \
					2. * (*(a + index) * (*(b + index) * *(f + index) + *(c + index) * *(fx + index) + \
						* (d + index) * *(fx2 + index) - *(fs + index)) + \
						* (b + index) * (*(c + index) * *(x + index) - *(s + index)) + \
						* (c + index) * (*(d + index) * *(x3 + index) - *(sx + index)) - \
						* (d + index) * *(sx2 + index));
				if (ExptsData::use_cubic[index])
				{
					//here begins the cubic correction part
					*browser += 2 * *(e + index) * (*(a + index) * *(fx3 + index) + *(b + index) * *(x3 + index) + *(c + index) * *(x4 + index) + \
						* (d + index) * *(x5 + index) - *(sx3 + index)) + *(e + index) * *(e + index) * *(x6 + index);
				}
				if (*browser <= 0)
				{
					//The chi_square is close to zero, and it cannot be calculated according to the usual
					//way because of the rounding errors, it has to be recalculated

					//for every data point
					pq = edata.sq_qfinder[i];
					pexp = edata.sq_sfinder[i];//sets the pointer to the beginning of S(Q) data for i-th expt data series
					pcalc = calc.sqfinder[i];//sets the pointer to the beginning of the calculated S(Q) data for i-th data series
					diff = 0;
					exprenorm = 0;
					for (j = 0; j < edata.sqused[i]; j++)
					{
						background = (b[index] + c[index] * *pq + d[index] * *pq * *pq + e[index] * *pq * *pq * *pq);
						exprenorm += pow(a[index] * *pexp + background, 2);
						diff += pow(*pcalc - (a[index] * *pexp + background), 2);
						pcalc++;
						pexp++;
						pq++;

					}
					exprenorm_calculated = 1;
					chicomp[index] = diff;
					//Calc_method is an integer indicator,if the iexpt-th bit is 1, indicates, that for
					//this data series the chi2 has to be calculated using the square of the diff method 
					//if calc_method=0, chi2 could be calculated normally for each data series, if it is
					// greated than zero, than then the data series(es) denoted by non-zero bits were calculated
					//as the square of the difference, and the value of the renormalization coeff a, if 
					//it was calculated at all, is questionable
					//
					if (i == 0 && ngr == 0)
						calc_method++;//this will be the first experimental series
					else
					{
						k = 1;
						for (j = 0; j < index; j++)
							k *= 2;
						calc_method += k;
					}
				}
				if (ExptsData::squseR[i])
				{
					if (!exprenorm_calculated)//it was not yet calculated
					{
						//for every data point
						pq = edata.sq_qfinder[i];
						pexp = edata.sq_sfinder[i];//sets the pointer to the beginning of S(Q) data for i-th expt data series
						exprenorm = 0;
						for (j = 0; j < edata.sqused[i]; j++)
						{
							background = (b[index] + c[index] * *pq + d[index] * *pq * *pq + e[index] * *pq * *pq * *pq);
							exprenorm += pow(a[index] * *pexp + background, 2);
							pexp++;
							pq++;
						}
					}
					if (fabs(exprenorm) > TOLERANCE)
					{
						Rw[index] = sqrt(*browser / exprenorm);
						*browser = Rw[index] / *(edata.sqsigma + i);//calculating Rw and dividing by sigma 
						Rw[index] *= 100;//this is what was outputed sofar
					}
					else
						Rw[index] = -1000;
				}
				else
				{
					//dividing by the standard dev.squared
					*browser /= *(edata.sqsigma + i) * *(edata.sqsigma + i);//dividing by sigma squared
				}
				total += *browser;//updating the global result
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "S(Q) ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				browser++;
				index++;
			}
		}
		browser = chicomp + ngr + nsq;//reset, as it might not have been increased
		index = ngr + nsq;
		//the F(Q) components
		for (i = 0; i < nfq; i++)//for each F(Q) data set
		{
			if (!Move::nomove || (Move::makemucorr && edata.fqIQbackgcorr[i])|| (Move::makefprimeshift && edata.fqAXS[i] > 0))
			{
				exprenorm_calculated = 0;
				if (fit_index[index] < 32)//chi2 was already calculated during nonlin regression, no need to calculate it here
				{
					//the developed sum
					*browser = *(ss + index) + *(a + index) * *(a + index) * *(ff + index) + \
						* (edata.fqused + i) * *(b + index) * *(b + index) + \
						* (x2 + index) * (*(c + index) * *(c + index) + 2. * *(b + index) * *(d + index)) + \
						* (d + index) * *(d + index) * *(x4 + index) + \
						2. * (*(a + index) * (*(b + index) * *(f + index) + *(c + index) * *(fx + index) + \
							* (d + index) * *(fx2 + index) - *(fs + index)) + \
							* (b + index) * (*(c + index) * *(x + index) - *(s + index)) + \
							* (c + index) * (*(d + index) * *(x3 + index) - *(sx + index)) - \
							* (d + index) * *(sx2 + index));
					
					//here begins the cubic correction part
					if (ExptsData::use_cubic[index])
					{
						*browser += 2 * *(e + index) * (*(a + index) * *(fx3 + index) + *(b + index) * *(x3 + index) + *(c + index) * *(x4 + index) + \
							* (d + index) * *(x5 + index) - *(sx3 + index)) + *(e + index) * *(e + index) * *(x6 + index);
					}

					if (*browser <= 0)
					{
						//The chi_square is close to zero, and it cannot be calculated according to the usual
						//way because of the rounding errors, it has to be recalculated

						//for every data point
						pq = edata.fq_qfinder[i];
						pexp = edata.fq_ffinder[i];//sets the pointer to the beginning of F(Q) data for i-th expt data series
						pcalc = calc.fqfinder[i];//sets the pointer to the beginning of the calculated F(Q) data for i-th data series
						diff = 0;
						exprenorm = 0;
						for (j = 0; j < edata.fqused[i]; j++)
						{
							background = (b[index] + c[index] * *pq + d[index] * *pq * *pq + e[index] * *pq * *pq * *pq);
							exprenorm += pow(a[index] * *pexp + background, 2);
							diff += pow(*pcalc - (a[index] * *pexp + background), 2);
							pcalc++;
							pexp++;
							pq++;
						}
						exprenorm_calculated = 1;
						chicomp[index] = diff;
						//Calc_method is an integer indicator,if the iexpt-th bit is 1, indicates, that for
						//this data series the chi2 has to be calculated using the square of the diff method 
						//if calc_method=0, chi2 could be calculated normally for each data series, if it is
						// greated than zero, than then the data series(es) denoted by non-zero bits were calculated
						//as the square of the difference, and the value of the renormalization coeff a, if 
						//it was calculated at all, is questionable
						//
						if (i == 0 && ngr + nsq == 0)
							calc_method++;//this will be the first experimental series
						else
						{
							k = 1;
							for (j = 0; j < index; j++)
								k *= 2;
							calc_method += k;
						}
					}
				}
				if (ExptsData::fquseR[i])
				{
					if (!exprenorm_calculated)//it was not yet calculated
					{
						//for every data point
						pq = edata.fq_qfinder[i];
						pexp = edata.fq_ffinder[i];//sets the pointer to the beginning of F(Q) data for i-th expt data series
						exprenorm = 0;
						for (j = 0; j < edata.fqused[i]; j++)
						{
							background = (b[index] + c[index] * *pq + d[index] * *pq * *pq + e[index] * *pq * *pq * *pq);
							exprenorm += pow(a[index] * *pexp + background, 2);
							pexp++;
							pq++;
						}
					}
					if (fabs(exprenorm) > TOLERANCE)
					{
						Rw[index] = sqrt(*browser / exprenorm);
						*browser = Rw[index] / *(edata.fqsigma + i);//calculating Rw and dividing by sigma 
						Rw[index] *= 100;//this is what was outputed sofar
					}
					else
						Rw[index] = -1000;

				}
				else
				{
					//dividing by the standard dev.squared
					*browser /= *(edata.fqsigma + i) * *(edata.fqsigma + i);//dividing by sigma squared
				}
				total += *browser;//updating the global result
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "F(Q) ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}
			}
			browser++;
			index++;
		}

		if (!Move::nomove)
		{
			//the F(g) components
			for (i = 0; i < nfg; i++)//for each F(g) data set
			{
				//the developed sum
				exprenorm_calculated = 0;
				*browser = *(ss + index) + *(a + index) * *(a + index) * *(ff + index) + \
					* (edata.fgused + i) * *(b + index) * *(b + index) + \
					* (x2 + index) * (*(c + index) * *(c + index) + 2. * *(b + index) * *(d + index)) + \
					* (d + index) * *(d + index) * *(x4 + index) + \
					2. * (*(a + index) * (*(b + index) * *(f + index) + *(c + index) * *(fx + index) + \
						* (d + index) * *(fx2 + index) - *(fs + index)) + \
						* (b + index) * (*(c + index) * *(x + index) - *(s + index)) + \
						* (c + index) * (*(d + index) * *(x3 + index) - *(sx + index)) - \
						* (d + index) * *(sx2 + index));
				//here begins the cubic correction part
				if (ExptsData::use_cubic[index])
				{
					*browser += 2 * *(e + index) * (*(a + index) * *(fx3 + index) + *(b + index) * *(x3 + index) + *(c + index) * *(x4 + index) + \
						* (d + index) * *(x5 + index) - *(sx3 + index)) + *(e + index) * *(e + index) * *(x6 + index);
				}

				if (*browser <= 0)
				{
					//The chi_sguare is close to zero, and it cannot be calculated according to the usual
					//way because of the rounding errors, it has to be recalculated

					//for every data point
					pg = edata.fg_gfinder[i];
					pexp = edata.fg_ffinder[i];//sets the pointer to the beginning of F(g) data for i-th expt data series
					pcalc = calc.fgfinder[i];//sets the pointer to the beginning of the calculated F(g) data for i-th data series
					diff = 0;
					exprenorm = 0;
					for (j = 0; j < edata.fgused[i]; j++)
					{
						background = (b[index] + c[index] * *pg + d[index] * *pg * *pg + e[index] * *pg * *pg * *pg);
						exprenorm += pow(a[index] * *pexp + background, 2);
						diff += pow(*pcalc - (a[index] * *pexp + background), 2);
						pcalc++;
						pexp++;
						pg++;

					}
					exprenorm_calculated = 1;
					chicomp[index] = diff;
					//Calc_method is an integer indicator,if the iexpt-th bit is 1, indicates, that for
					//this data series the chi2 has to be calculated using the square of the diff method 
					//if calc_method=0, chi2 could be calculated normally for each data series, if it is
					// greated than zero, than then the data series(es) denoted by non-zero bits were calculated
					//as the square of the difference, and the value of the renormalization coeff a, if 
					//it was calculated at all, is questionable
					//
					if (i == 0 && ngr + nsq + nfq == 0)
						calc_method++;//this will be the first experimental series
					else
					{
						k = 1;
						for (j = 0; j < index; j++)
							k *= 2;
						calc_method += k;
					}
				}
				if (ExptsData::fguseR[i])
				{
					if (!exprenorm_calculated)//it was not yet calculated
					{
						//for every data point
						pg = edata.fg_gfinder[i];
						pexp = edata.fg_ffinder[i];//sets the pointer to the beginning of F(g) data for i-th expt data series
						exprenorm = 0;
						for (j = 0; j < edata.fgused[i]; j++)
						{
							background = (b[index] + c[index] * *pg + d[index] * *pg * *pg + e[index] * *pg * *pg * *pg);
							exprenorm += pow(a[index] * *pexp + background, 2);
							pexp++;
							pg++;
						}
					}
					if (fabs(exprenorm) > TOLERANCE)
					{
						Rw[index] = sqrt(*browser / exprenorm);
						*browser = Rw[index] / *(edata.fgsigma + i);//calculating Rw and dividing by sigma 
						Rw[index] *= 100;//this is what was outputed sofar
					}
					else
						Rw[index] = -1000;

				}
				else
				{
					//dividing by the standard dev.squared
					*browser /= *(edata.fgsigma + i) * *(edata.fgsigma + i);//dividing by sigma squared
				}
				total += *browser;//updating the global result
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "F(g) ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				browser++;
				index++;
			}//end of i cycle for F(g) series
		}

		browser = chicomp + ngr + nsq + nfq + nfg;//reset, as it might not have been increased
		index = ngr + nsq + nfq + nfg;
		//the E(k) components
		for (i = 0; i < nek; i++)//for each E(k) data set
		{
			if (!Move::nomove || (Move::makeE0shift && edata.ek_ngrid_in[i] > 0))
			{//calculate if there is no E0 shift in this step, or if there is only for those, where there is
			//resets the value to zero(needed in the main loop)

				//the developed sum
				exprenorm_calculated = 0;
				*browser = *(ss + index) + *(a + index) * *(a + index) * *(ff + index) + \
					* (edata.ekused + i) * *(b + index) * *(b + index) + \
					2 * *(a + index) * *(b + index) * *(f + index) - 2 * *(a + index) * *(fs + index)\
					- 2 * *(b + index) * *(s + index);
				if (*browser <= 0)
				{
					//The chi_square is close to zero, and it cannot be calculated according to the usual
					//way because of the rounding errors, it has to be recalculated

					//for every data point
					pexp = edata.ek_efinder[i];//sets the pointer to the beginning of E(k) data for i-th expt data series
					pcalc = calc.ekfinder[i];//sets the pointer to the beginning of the calculated E(k) data for i-th data series
					diff = 0;
					exprenorm = 0;
					for (j = 0; j < edata.ekused[i]; j++)
					{
						exprenorm += pow(a[index] * *pexp + b[index], 2);
						diff += pow(*pcalc - (a[index] * *pexp + b[index]), 2);
						pcalc++;
						pexp++;

					}
					exprenorm_calculated = 1;
					chicomp[index] = diff;
					//Calc_method is an integer indicator,if the iexpt-th bit is 1, indicates, that for
					//this data series the chi2 has to be calculated using the square of the diff method 
					//if calc_method=0, chi2 could be calculated normally for each data series, if it is
					// greater than zero, than the data series denoted by non-zero bits were calculated
					//as the square of the difference, and the value of the renormalization coeff a, if 
					//it was calculated at all, is questionable
					//
					if (i == 0 && ngr + nsq + nfq + nfg == 0)
						calc_method++;//this will be the first experimental series
					else
					{
						k = 1;
						for (j = 0; j < index; j++)
							k *= 2;
						calc_method += k;
					}
				}
				if (ExptsData::ekuseR[i])
				{
					if (!exprenorm_calculated)//it was not yet calculated
					{
						//for every data point
						pexp = edata.ek_efinder[i];//sets the pointer to the beginning of S(Q) data for i-th expt data series
						exprenorm = 0;
						for (j = 0; j < edata.ekused[i]; j++)
						{
							exprenorm += pow(a[index] * *pexp + b[index], 2);
							pexp++;
						}
					}
					if (fabs(exprenorm) > TOLERANCE)
					{
						Rw[index] = sqrt(*browser / exprenorm);
						*browser = Rw[index] / *(edata.eksigma + i);//calculating Rw and dividing by sigma 
						Rw[index] *= 100;//this is what was outputed sofar
					}
					else
						Rw[index] = -1000;
				}
				else
				{
					//dividing by the standard dev.squared
					*browser /= pow(*(edata.eksigma + i), 2);//dividing by sigma squared 
				}
				total += *browser;//updating the global result
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "E(k) ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}
			}
	
			browser++;
			index++;
		}

		if (!Move::nomove)
		{
			for (i = 0; i < ncos; i++)//for each CosDistrConst data set
			{
				//the developed sum, there is no renoemalization!
				*browser = *(ff + index) + *(ss + index) - 2 * *(fs + index);
				if (*browser <= 0)
				{
					//The chi_square is close to zero, and it cannot be calculated according to the usual
					//way because of the rounding errors, it has to be recalculated

					//for every data point
					pexp = cosconst.theordistr + CosDistrConst::cumul_bins[i] + CosDistrConst::first_bin[i];//sets the pointer to the beginning of theoretical distribution for i-th data series
					pcalc = cosconst.cosinedistr + CosDistrConst::cumul_bins[i] + CosDistrConst::first_bin[i];//sets the pointer to the beginning of the calculated distribution for the i-th data series
					diff = 0;
					if (CosDistrConst::method[i] == 2)
					{
						for (j = CosDistrConst::first_bin[i]; j < CosDistrConst::last_bin[i]; j++)//for all used data points
						{
							diff += *pcalc * *pcalc * *pexp * *pexp;//the squared difference of the calculated series weighted by cos_theta-dependent weigth from 0 is calculated  
							pcalc++;
							pexp++;
						}

					}
					else
					{
						//The squared difference between the calculated and theoretical is determined 
						for (j = 0; j < CosDistrConst::npoints[i]; j++)
						{
							//Renormalization is not yet included GO!!!!
							diff += (*pcalc - *pexp) * (*pcalc - *pexp);
							pcalc++;
							pexp++;

						}
					}
					*browser = diff;
					//Calc_method is an integer indicator,if the iexpt-th bit is 1, indicates, that for
					//this data series the chi2 has to be calculated using the square of the diff method 
					//if calc_method=0, chi2 could be calculated normally for each data series, if it is
					// greater than zero, than then the data series(es) denoted by non-zero bits were calculated
					//as the square of the difference, and the value of the renormalization coeff a, if 
					//it was calculated at all, is questionable
					//
					if (i == 0 && ngr + nsq + nfq + nfg + nek == 0)
						calc_method++;//this will be the first experimental series
					else
					{
						k = 1;
						for (j = 0; j < index; j++)
							k *= 2;
						calc_method += k;
					}
				}
				//dividing by the standard dev.squared
				*browser /= *(CosDistrConst::weights + i) * *(CosDistrConst::weights + i);//dividing by sigma squared
				total += *browser;//updating the global result
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "cosine disribution of angle constraint ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				browser++;
				index++;
			}

			//add the coordinate constraints contributions
			pi1 = icc.nsatisfy;
			pi2 = icc.ncentral;
			pd1 = icc.fraction;
			pd2 = icc.weights;
			for (i = 0; i < nicc; i++)
			{
				for (isubconst = 0; isubconst < icc.n_subconst[i]; isubconst++)
				{ 
					if (*pi2 == 0)//here this cannot happen...
						temp = 0;
					else
						temp = double(*pi1) / double(*pi2);//fraction of atoms satisfying the CC
					temp -= *pd1;
					*browser = temp * temp / (*pd2 * *pd2);
					if (calc_sigma && RunParams::lead_series_ind == index)
					{
						if (fabs(*browser - 0.0) < TOLERANCE)
						{
							mystrcpy(name, NAME_SIZE, "coordination number constraint ");
							IntToStr(&name_ext, i + 1);
							mystrcat(name, NAME_SIZE, name_ext);
							ZeroChi2Error(name);
						}
						lead_chi2 = *browser;
					}

					total += *browser++;
					pd1++;
					pd2++;
					pi1++;
					index++;
				}
				pi2++;

			}

			//add the average coordination constraints contribution
			pi1 = avcc.neighbourcount;//points to the total coordination number
						//for the config.
			pi2 = avcc.ncentral;//points to the type of the CC central atom
			pd1 = avcc.acnreq;//required average coordination number
			pd2 = avcc.weights;//weights for the av. CC.
			for (i = 0; i < nacc; i++)
			{
				if (*pi2 == 0)//here this cannot happen...
					temp = 0;
				else
					temp = double(*pi1) / double(*pi2);
				temp -= *pd1;
				*browser = temp * temp / (*pd2 * *pd2);
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "average coordination constraint ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				total += *browser++;
				pi1++;
				pi2++;
				pd1++;
				pd2++;
				index++;
			}
		
		

#ifdef _ADVANCED_GEOM_CONST
			//add the common neighbour constraints contributions
			pi1 = conc.nsatisfy;
			pi2 = conc.nprimary;
			pd1 = conc.fraction;
			pd2 = conc.weights;
			for (i = 0; i < ncommonneigh; i++)
			{
				if (*pi2 == 0)
					temp = 0;
				else
					temp = double(*pi1) / double(*pi2);//fraction of atoms satisfying the CONC
				temp -= *pd1;
				*browser = temp * temp / (*pd2 * *pd2);
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "common neighbour constraint ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				total += *browser++;
				pd1++;
				pd2++;
				pi1++;
				pi2++;
				index++;
			}
			//add the second neighbour constraints contributions
			pi1 = snc.nsatisfy;
			pi2 = snc.ncentral;
			pd1 = snc.fraction;
			pd2 = snc.weights;
			for (i = 0; i < nsecondneigh; i++)
			{
				if (*pi2 == 0)//here this cannot happen...
					temp = 0;
				else
					temp = double(*pi1) / double(*pi2);//fraction of atoms satisfying the CONC
				temp -= *pd1;
				*browser = temp * temp / (*pd2 * *pd2);
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "second neighbour constraint ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				total += *browser++;
				pd1++;
				pd2++;
				pi1++;
				pi2++;
				index++;
			}
			//add the bond valence sum neighbour constraints contributions
			pi2 = bvs.ncentral;
			logfile.precision(15);
			pd1 = bvs.target_valence;
			pd2 = bvs.weights;
			double *pd3 = bvs.valence;
			for (i = 0; i < nbvs; i++)
			{
				temp = 0;
				for (int ic = 0; ic < bvs.ncentral[i]; ic++)
					temp += pow(*pd3++ - *pd1, 2);
					
				*browser = temp / (*pd2 * *pd2);
				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "bond valence sum constraint ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}
					lead_chi2 = *browser;
				}

				total += *browser++;
				pd1++;
				pd2++;
				pi2++;
				index++;
			}
#endif
#ifdef _LOCAL_INV
			double *pmain = 0, *paux = 0;

			//the contribution to chi2 will be calulated parallel, as it takes quite a lot of time
				//------------------MULTITHREADING-------------------

			thread_obj.LOC_count = 0;//counter indicating how many thread finished with the calculation, reset it to 0
#ifdef _TEST_MODE
			Threads::lctime1[RunParams::nthreads - 1] = std::chrono::high_resolution_clock::now();
#endif
			if (RunParams::loc_chi2_mode == 1 && RunParams::nthreads > 1)//signal to threads to start calculation 
			{
				std::unique_lock<std::mutex> guard(thread_obj.inv_mutex);//lock mutex
				thread_obj.start_flag = 1;//time to proceed
				thread_obj.check_inv_start.notify_all();
			}

			//calculate the chisqure contribution coming from the main thread, has to be done in two step in case of loc_chi_mode=1
			//as different segmentation has to be used 

			if (RunParams::loc_chi2_mode == 1)
				CalcLocAvThread_dist(thread_obj.thread_arg[RunParams::nthreads - 1]);

			//Threads have to wait, till all of them finished
			if (RunParams::loc_chi2_mode == 1 && RunParams::nthreads > 1)
			{
				std::unique_lock<std::mutex> guard(thread_obj.LOC_count_mutex);//lock mutex
				thread_obj.LOC_count++;
				if (thread_obj.LOC_count == RunParams::nthreads)//time to proceed
				{
					thread_obj.start_flag = 0;//reset it to 0
					thread_obj.LOC_count = 0;//reset it to 0
					thread_obj.LOC_count2 = 0;//reset it to 0
					thread_obj.LOC_check_count.notify_all();
				}
				else thread_obj.LOC_check_count.wait(guard);
			}

#ifdef _TEST_MODE
			Threads::lctime2[RunParams::nthreads - 1] = std::chrono::high_resolution_clock::now();
			elapsed = Threads::lctime2[RunParams::nthreads - 1]- Threads::lctime1[RunParams::nthreads - 1];
			Threads::dur_chiloc1[RunParams::nthreads - 1]+=elapsed.count();
#endif

			//Update and average the av_loc_hist
			if (RunParams::loc_chi2_mode == 1)
			{

				for (i = 1; i < RunParams::nthreads; i++)
				{
					pmain = av_loc_hist;//this is the 0.th thread segment, this will hold the total
					paux = thread_obj.thread_arg[i].av_loc_hist_finder;
					for (j = 0; j < sum_loc_natoms_tot * SimpleCfg::ntypes; j++)
					{
						*pmain += *paux++;
						pmain++;
					}

				}
				pmain = av_loc_hist;
				for (i = 0; i < SimpleCfg::ntypes; i++)
				{
					for (j = 0; j < SimpleCfg::ntypes; j++)
					{
						for (k = 0; k < loc_natoms[j]; k++)
						{
							*pmain = RoundNearInt(*pmain / SimpleCfg::pnatoms[i]);//calculate the average
							pmain++;
						}
					}
				}
			}


#ifdef _TEST_MODE
			Threads::lctime1[RunParams::nthreads - 1] = std::chrono::high_resolution_clock::now();
			elapsed = Threads::lctime1[RunParams::nthreads - 1]- Threads::lctime2[RunParams::nthreads - 1];
			Threads::dur_chiloc2[RunParams::nthreads - 1]+=elapsed.count();
#endif

			thread_obj.LOC_count2 = 0;//counter indicating how many thread finished with the calculation, reset it to 0
			if (RunParams::nthreads > 1)//signal to threads to start calculation 
			{
				std::unique_lock<std::mutex> guard(thread_obj.inv_mutex);//lock mutex
				thread_obj.start_flag = 1;//time to proceed
				thread_obj.check_inv_start.notify_all();
			}

			//calculate the chisquare contribution coming from the main thread
			(this->*CalcLocChi2Thread)(thread_obj.thread_arg[RunParams::nthreads - 1]);


#ifdef _TEST_MODE
			Threads::lctime2[RunParams::nthreads - 1] = std::chrono::high_resolution_clock::now();
			elapsed = Threads::lctime2[RunParams::nthreads - 1]- Threads::lctime1[RunParams::nthreads - 1];
			Threads::dur_chiloc3[RunParams::nthreads - 1]+=elapsed.count();
#endif

			//Threads have to wait, till all of them finished, as each of them will use the av_loc_hist values
			if (RunParams::nthreads > 1)
			{
				std::unique_lock<std::mutex> guard(thread_obj.LOC_count_mutex2);//lock mutex
				thread_obj.LOC_count2++;
				if (thread_obj.LOC_count2 == RunParams::nthreads)//time to proceed
				{
					thread_obj.start_flag = 0;//reset it to 0
					thread_obj.LOC_count = 0;//reset it to 0
					thread_obj.LOC_count2 = 0;//reset it to 0
					thread_obj.LOC_check_count2.notify_all();
				}
				else    thread_obj.LOC_check_count2.wait(guard);
			};
#ifdef _TEST_MODE
			Threads::lctime1[RunParams::nthreads - 1] = std::chrono::high_resolution_clock::now();
			elapsed = Threads::lctime1[RunParams::nthreads - 1]- Threads::lctime2[RunParams::nthreads - 1];
			Threads::dur_chiloc4[RunParams::nthreads - 1]+=elapsed.count();
#endif
			//Summing the contributions

			for (i = 0; i < nlocint; i++)
			{
				*browser = 0.0;
				for (j = 0; j < RunParams::nthreads; j++)
					*browser += *(loc_chi_finder[j] + i);

				if (calc_sigma && RunParams::lead_series_ind == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "local invariance ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
						ZeroChi2Error(name);
					}

					lead_chi2 = *browser;
				}

				*browser /= pow(RunParams::loc_inv_sigma[i], 2);
				total += *browser++;
				index++;
			}

#endif

		}
		else
		{

			total = 0.0;
			for (int ii = 0; ii < chi_potstart; ii++)
				total += chicomp[ii];
		}
		return(calc_method);
	}

#ifdef _LOCAL_INV
	//calculating the chisquared contribution by the threads 
	void ChiSquared::CalcLocChiSquaredThread_bin(ThreadArg &thread_arg)
	{
		int i,ibin,ipartial,itype,jtype,thread_id;
		int *phloc;
		int min_ind,max_ind;//cycle limits
		longint *phsum;
		double n,*sum;

		//calculate the chi2 coming from the local invariance
		thread_id=thread_arg.thread_index;//the ID of the thread
		sum=loc_chi_finder[thread_id];//pointer to the chi2 contribution of the local invariance of this thread
		*sum=0.0;
		for (itype=0;itype<HistoSet::ntypes;itype++)
		{
			min_ind=thread_arg.min_atom_index[itype];//the index of the first atom in its own type
			max_ind=thread_arg.max_atom_index[itype];//the index of the last atom in its own type
			
			phloc=hist.local_count+SimpleCfg::cumul[itype]*HistoSet::sum_bins_offset;//local histogram for the first atom of this type 
			for (jtype=0;jtype<HistoSet::ntypes;jtype++)
			{
				n = SimpleCfg::pnatoms[itype];//for calculating the average: number of centre atoms*1 for pure and *2 for mixed (local sum is counted for both atoms as well, so sum should be divided by two for mixed partials
				if (itype != jtype)
					n *= 2.;
				ipartial=(itype<=jtype ? (itype*HistoSet::ntypes-(itype*(itype+1)/2)+ jtype) : (jtype*HistoSet::ntypes-(jtype*(jtype+1)/2)+ itype));
				
				phsum=hist.local_sum_count+hist.loc_offset[ipartial];//first bin to use from this partial in the average histogram
				for (ibin=0;ibin<HistoSet::max_loc_nbins;ibin++)
				{
					for (i=min_ind;i<=max_ind;i++)//going through all the atoms of this type
						*sum+=pow(*(phloc+i*HistoSet::sum_bins_offset)-*phsum/n,2)/DataMat::dV_loc[ibin];	
					phloc++;
					phsum++;
				}
			}
		}
	};

	//calculating the distance_based averega histogram, multi-thraeding cannot be introduced sufficiently
	void ChiSquared::CalcLocAvThread_dist(ThreadArg &thread_arg)
	{
		int i,ibin,itype,jtype,ipartial,bin_count,sum_count;
		int ithread;
		int *phloc,next;
		int *my_first_bin,*my_start_count;
		int min_ind,max_ind,min_neigh_ind,max_neigh_ind;//cycle limits
		double *my_av,*av_start;

#ifdef _TEST_MODE
		int thread_id;
		thread_id=thread_arg.thread_index;
		double tloc[2];
		std::chrono::duration<double, std::milli> elapsed;
		std::chrono::time_point<std::chrono::high_resolution_clock> t1,t2;

		for (i=0;i<2;i++)
			tloc[i]=0;
#endif				
		//setting to average to zero for all the partials of this thread
		av_start=thread_arg.av_loc_hist_finder;//starting point for this thread in av_loc_hist array
		my_av=av_start;
		for (i=0;i<sum_loc_natoms_tot*SimpleCfg::ntypes;i++)
			*my_av++=0;

		for (itype=0;itype<HistoSet::ntypes;itype++)
		{
			for (jtype=0;jtype<HistoSet::ntypes;jtype++)
			{
				ipartial=itype*HistoSet::ntypes+jtype;//1-2 and 2-1 are not the same!
								

				min_ind=thread_arg.min_atom_index[itype];//the index of the first central atom in its own type to calculate for
				max_ind=thread_arg.max_atom_index[itype];//the index of the last central atom in its own type
			
				if (min_ind<0 &&max_ind<0)
					continue;
			
				
				//as the segmentation for the local histogram in the av_loc_hist calculation is vertical (according to ncentral)
				//and later in case of loc chi2 calculation is horizontal (neighbour based), starting points has to be establised for 
				//the loc chi2 calculation, so for a central we will go through the neighbours according to the loc chi2 calculation segments
				//to set the first_loc_ind and starting_count arrays for each later threads 
				//That is the reason of the ithread cycle
				//first calculate the average
				for (i=min_ind;i<=max_ind;i++)//going through the central atoms of this type for this thread to handle
				{
					my_av=av_start+av_loc_hist_offset[ipartial];//starting point for this thread in av_loc_hist array
					//starting point for the first thread in the loc chi2 calculation, the offset is thread_obj.thread_arg[ithread].loc_hist_offset for the other threads
					my_first_bin=first_loc_ind+cumalpart_loc_atoms[ipartial]+i;//starting point in first_loc_ind array for this partial
					my_start_count=starting_count+cumalpart_loc_atoms[ipartial]+i;//starting point for this thread in start_count array for this partia
					//local histogram for the i-th atom of this type 
					phloc=hist.local_count+SimpleCfg::cumul[itype]*HistoSet::sum_bins_offset+HistoSet::loc_offset[jtype]+i*HistoSet::sum_bins_offset;
					//fist skipping the neighbours used by other threads
					sum_count=0;//min ind will start with 0
					ibin=0;
					bin_count=*phloc;//set the first in case it does not go into the while loop
					for (ithread=0;ithread<RunParams::nthreads;ithread++)
					{
						min_neigh_ind=thread_obj.thread_arg[ithread].min_loc_atom_index[jtype];	
						max_neigh_ind=thread_obj.thread_arg[ithread].max_loc_atom_index[jtype];	
						if (min_neigh_ind<0 && max_neigh_ind<0)
							continue;
						if (itype==jtype && max_neigh_ind==SimpleCfg::pnatoms[itype]-1)//this is the thread handling the largest indexed neighbours 
						{
							max_neigh_ind--;
							if (max_neigh_ind<min_neigh_ind)//can happen, if there are only few atoms
								continue;//go to next jtype
						}
#ifdef _TEST_MODE
						t1 = std::chrono::high_resolution_clock::now();
#endif
						
						next=1;//starting value

						while (sum_count<min_neigh_ind)//going through the neighbour-containing bins
						{
							while (bin_count>0 && sum_count<min_neigh_ind)//while there is count in this bin
							{
								next=1;
								bin_count--;
								sum_count++;
								if (bin_count>0)
									next=0;
							}
							if (next)
							{
								phloc++;//go to the next bin for this atom in local histogram
								bin_count=*phloc;
								ibin++;//index of next bin
							}
							
						}
						while (bin_count==0)//searching for the first bin with count
						{
							phloc++;//go to the next bin for this atom in local histogram
							bin_count=*phloc;
							ibin++;//index of next bin
						}
						
						*(my_first_bin+thread_obj.thread_arg[ithread].loc_hist_offset)=ibin;
						*(my_start_count+thread_obj.thread_arg[ithread].loc_hist_offset)=bin_count;
#ifdef _TEST_MODE
						t2 = std::chrono::high_resolution_clock::now();
						elapsed = t2 - t1;
						tloc[0]+=elapsed.count();
												
#endif
						//now calculate the average, which means the average bin for the i-th neighbour of each central atom
						while (sum_count<=max_neigh_ind)//going through the neighbour-containing bins
						{
							while (bin_count>0 && sum_count<=max_neigh_ind)//while there is count in this bin
							{
								next=1;
								
								*my_av+=ibin;//adding the index of this bin to the average
								my_av++;
								bin_count--;
								sum_count++;
								if (bin_count>0)
									next=0;
							}
							if (next)
							{
								phloc++;//go to the next bin for this atom in local histogram
								bin_count=*phloc;
								ibin++;//index of next bin
							}
						}
#ifdef _TEST_MODE
						t1 = std::chrono::high_resolution_clock::now();
						elapsed = t1 - t2;
						tloc[1]+=elapsed.count();
#endif
					}//end of ithread cycle
					my_first_bin++;
					my_start_count++;
				}//end of i cycle for central atoms
			}//end jtype
		}//end itype
#ifdef _TEST_MODE
		Threads::dur_loc1[thread_id]+=tloc[0];
		Threads::dur_loc2[thread_id]+=tloc[1];
#endif

	}

	//calculating the chisquared contribution by the threads 
	//the segmentaion will be done the way, that each thread will handle a consecutive segment of neighbour atoms (sorted by the distance)
	//which caused a count in the local histogram. Which is the first and last histogram bin to handle cannot be said beforhand, so the threads
	//have to run through the local histogram bins of an atom, to find the first one to start calculation for. 
	void ChiSquared::CalcLocChiSquaredThread_dist(ThreadArg &thread_arg)
	{
		int i,j,k,ibin,itype,jtype,thread_id,bin_count,ipartial, nused_thread,ilocint,start_int;
		int *phloc;
		int *my_first_bin,*my_start_count;
		int min_ind,max_ind,neigh_offset;//cycle limits
		int *my_first_loc_ind,*my_starting_count;
		double *my_av,*av_start;
		double *sum;

		thread_id=thread_arg.thread_index;//the ID of the thread
		my_first_loc_ind=first_loc_ind+thread_arg.loc_hist_offset;
		my_starting_count=starting_count+thread_arg.loc_hist_offset;
		
		sum=loc_chi_finder[thread_id];

		for (ilocint=0;ilocint<nlocint;ilocint++)
			sum[ilocint]=0.0;

		
#ifdef _TEST_MODE
		double tloc[2];
		std::chrono::duration<double, std::milli> elapsed;
		std::chrono::time_point<std::chrono::high_resolution_clock> t1,t2;
		for (i=0;i<2;i++)
			tloc[i]=0;
		t1 = std::chrono::high_resolution_clock::now();
#endif
		for (itype=0;itype<HistoSet::ntypes;itype++)
		{
			for (jtype=0;jtype<HistoSet::ntypes;jtype++)
			{
				min_ind=thread_arg.min_loc_atom_index[jtype];//the index of the first neighbour atom in its own type to calculate for
				max_ind=thread_arg.max_loc_atom_index[jtype];//the index of the last neighbour atom in its own type
				if (min_ind<0 &&max_ind<0)
					continue;
				nused_thread=abs(nused_loc[jtype]);//the index of the last used thread for this type
				if (itype==jtype && max_ind==SimpleCfg::pnatoms[itype]-1)//this is the thread handling the largest indexed neighbours 
				{
					if (nused_loc[jtype]<0 && thread_id==nused_thread-1)
						nused_thread--;//this will be the last thread to use
					max_ind--;
					if (max_ind<min_ind)//can happen, if there are only few atoms
						continue;//go to next jtype
				}
				
				neigh_offset=min_ind-thread_obj.thread_arg[0].min_loc_atom_index[jtype];//the offset for this thread from the beginnig of this partial in av_loc_hist
				
				ipartial=itype*HistoSet::ntypes+jtype;//1-2 and 2-1 are not the same!
				av_start=av_loc_hist+av_loc_hist_offset[ipartial]+neigh_offset;//starting point for this thread in av_loc_hist array											

				//now calculating the chi2
				//starting point for this thread 
				my_first_bin=my_first_loc_ind+cumalpart_loc_atoms[ipartial];//starting point in first_loc_ind array for this partial
				my_start_count=my_starting_count+cumalpart_loc_atoms[ipartial];//starting point for this thread in start_count array for this partia
				
				//set the first interval to calculate for for this thread
				start_int=0;
				while (min_ind>=loc_natoms_min[jtype*(nlocint+1)+start_int] && start_int<nlocint)
					start_int++;
				start_int--;

				for (i=0;i<SimpleCfg::pnatoms[itype];i++)//going through all the atoms of this type
				{
					//fist skipping the neighbours used by other threads
					ibin=*my_first_bin;
					bin_count=*my_start_count;
					//local histogram for the i-th atom of this type for the bij to start the calculation with
					phloc=hist.local_count+SimpleCfg::cumul[itype]*HistoSet::sum_bins_offset+HistoSet::loc_offset[jtype]+i*HistoSet::sum_bins_offset+ibin;
					//now calculate the chi2
					my_av=av_start;//starting point for this thread in av_loc_hist array

					ilocint=start_int;

					//going through the neighbours to calculate for 
					for (j=min_ind;j<=max_ind;j++)
					{
						//find the next bin with count
						while (bin_count==0)//while there is no count in this bin
						{
							phloc++;//go to the next bin for this atom in local histogram
							bin_count=*phloc;
							ibin++;//index of next bin
						}
						if (j==loc_natoms_min[jtype*(nlocint+1)+ilocint+1] && ilocint<nlocint-1)
							ilocint++;//this is the beginning of the new interval, care was to be taken, that do not increase, if it is the end of the last interval

						sum[ilocint]+=pow((ibin-*my_av)*loc_bin_width,2)/DataMat::dV_loc[(int)*my_av];	
						my_av++;
						bin_count--;
					}
					my_first_bin++;
					my_start_count++;
				}//end of i cycle

			
	/*			if (thread_id==0 )//&& itype==jtype)
				{
					for (i=0;i<nlocint;i++)
					{
						j=(i==nlocint-1 ? 0 : 1);
						min_loc_r_now[itype*ntypes*nlocint+jtypes*nlocint+i]=av_loc_hist[av_loc_hist_offset*loc_binwidth;
	*/
							//*(av_start+loc_natoms_min[jtype*(nlocint+1)+i]-loc_natoms_min[jtype*(nlocint+1)]) *loc_bin_width;//the starting average distance for this type in local invariance calc for this interval
	/*					max_loc_r_now[itype*ntypes*nlocint+jtypes*nlocint+i]=*(av_start+loc_natoms_min[jtype*(nlocint+1)+i+1]-j-loc_natoms_min[jtype*(nlocint+1)]) *loc_bin_width;//the starting average distance for this type in local invariance calc for this interval
cout<<"min "<<itype<<" "<<jtype<<"  "<<i<<"  "<<min_loc_r_now[itype*ntypes*nlocint+jtypes*nlocint+i]<<" "<<loc_natoms_min[jtype*(nlocint+1)+i]<<" "<<loc_natoms_min[jtype*(nlocint+1)]<<endl;
cout<<"max "<<itype<<" "<<jtype<<"  "<<i<<"  "<<max_loc_r_now[itype*ntypes*nlocint+jtypes*nlocint+i]<<" "<<loc_natoms_min[jtype*(nlocint+1)+i+1]<<" "<<loc_natoms_min[jtype*(nlocint+1)]<<endl;
					}
				}*/
				


			}//end of jtype cycle
		}//end of itype cycle

		k=0;
		for (itype=0;itype<HistoSet::ntypes;itype++)
		{
			
			for (jtype=0;jtype<HistoSet::ntypes;jtype++)
			{
				for (i=0;i<nlocint;i++)
				{
					j=(i==nlocint-1 ? 0 : 1);
					min_loc_r_now[itype*HistoSet::ntypes*nlocint+jtype*nlocint+i]=av_loc_hist[k+av_loc_hist_offset[jtype]+loc_natoms_min[jtype*(nlocint+1)+i]]*loc_bin_width;;//the starting average distance for this type in local invariance calc for this interval
					max_loc_r_now[itype*HistoSet::ntypes*nlocint+jtype*nlocint+i]=av_loc_hist[k+av_loc_hist_offset[jtype]+loc_natoms_min[jtype*(nlocint+1)+i+1]-j]*loc_bin_width;;//the finishing average distance for this type in local invariance calc for this interval
//cout<<"min "<<itype<<" "<<jtype<<"  "<<i<<"  "<<min_loc_r_now[itype*HistoSet::ntypes*nlocint+jtype*nlocint+i]<<" "<<k<<" "<<av_loc_hist_offset[jtype]<<" "<<loc_natoms_min[jtype*(nlocint+1)+i]<<endl;
//cout<<"max "<<itype<<" "<<jtype<<"  "<<i<<"  "<<max_loc_r_now[itype*HistoSet::ntypes*nlocint+jtype*nlocint+i]<<" "<<k<<" "<<av_loc_hist_offset[jtype]<<" "<<loc_natoms_min[jtype*(nlocint+1)+i+1-j]<<endl;
					
				}
			}
			k+=av_loc_hist_offset[HistoSet::ntypes];
		}
#ifdef _TEST_MODE
		t2 = std::chrono::high_resolution_clock::now();
		elapsed = t2 - t1;
		tloc[0]+=elapsed.count();
		
		//it is better to do here at the end, and use local variables in the routine to prevent false sharing
			
		Threads::dur_loc3[thread_id]+=tloc[0];		
		Threads::dur_loc4[thread_id]+=tloc[1];		
	
#endif
	
	};
#endif

void ChiSquared::ZeroChi2Error(const char *name, int mode)
{
	cout << "\n*****ERROR*****" << endl;
	cout << "The chi2 of the leading series (" << name << ") is around zero, so it cannot be used for scaling the other series'" << endl;
	cout<<"chi2 to it. Choose an other ";
	if (mode)
		cout<<"potential-related ";
	cout<<"leading series, or set all the other ";
	
	if (mode)
		cout<<"potential-related ";
	cout<<"sigma parameters to a positive value to avoid scaling! ";
	cout<<"Cannot run this way, exiting...";
	CleanExit();

}


//calculating the chi2 for the potential related components	
void ChiSquared::CalcPotChiSquared()
{
	int i,index=0;
	double *browser;//pointer to the chicomp array
	double *ppot,*pweight;//pointer to the potential and weigth
	double *my_lead_chi2;
	
	total2=0.0;//the potential related chi2 total contribution, if there is non-bonded/aenet potential
	if (RunParams::potential>0)
		my_lead_chi2=&lead_chi2_pot;//it is scaled  to lead_series_ind2
	else
	{
		my_lead_chi2=&lead_chi2;//it is scaled  to lead_series_ind
		index=chi_potstart;
	}

	browser=chicomp+chi_potstart;
	//calc the chi2 for the vdW interaction, in case of potential=10 (tabulated potential) for the tabulated potential
	//or in case of potential =20 for ANN
	#ifdef _AENET
	if (RunParams::potential==20)
	{
			pweight=&RunParams::aenet_weight;
			ppot=&aenet.E_tot;
	}
	else
	{
	#endif
		pweight=RunParams::vdW_weight;//weights for the vdW interaction
		ppot=(fabs(RunParams::NB_weight_mode)==1 ? FNC_POT::vdW_pot : &FNC_POT::vdW_tot_pot);
	#ifdef _AENET
	}
	#endif
	for (i=0;i<nNB;i++)
	{
		*browser=*ppot/(*pweight * *pweight);
		if (calc_sigma && RunParams::lead_series_ind2==index)
		{
			if (fabs(*browser-0.0)<TOLERANCE)
			{
				switch (RunParams::potential)
				{
				case 1:
				{
					mystrcpy(name, NAME_SIZE, "van der Waals potential ");
					if (nNB > 1)
					{
						mystrcat(name, NAME_SIZE, " partial ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
					}
					break;
				}
				case 10:
				{
					mystrcpy(name, NAME_SIZE, "tabulated potential ");
					mystrcat(name, NAME_SIZE, " partial ");
					IntToStr(&name_ext, FNC_POT::tabpot_index[i] + 1);
					mystrcat(name, NAME_SIZE, name_ext);
					break;
				}
				case 20:
				{
					mystrcpy(name, NAME_SIZE, "ANN potential");
					break;
				}
				}
				if (RunParams::potential!=20)//there cannot be other potential components
					ZeroChi2Error(name,1);
			}	
			lead_chi2_pot=*browser;
		}
		total2+=*browser++;
		ppot++;
		pweight++;
		index++;
	}
	if (RunParams::potential == 1)
	{
		//calc the chi2 for the Coulomb interaction
		ppot = (fabs(RunParams::NB_weight_mode) == 1 ? FNC_POT::Coulomb_pot : &FNC_POT::Coulomb_tot_pot);
		pweight = RunParams::Coulomb_weight;//weights for the Coulomb interaction
		for (i = 0; i < nNB; i++)
		{
			*browser = *ppot / (*pweight * *pweight);
			if (calc_sigma && RunParams::lead_series_ind2 == index)
			{
				if (fabs(*browser - 0.0) < TOLERANCE)
				{
					mystrcpy(name, NAME_SIZE, "Coulomb potential ");
					if (nNB > 1)
					{
						mystrcat(name, NAME_SIZE, " partial ");
						IntToStr(&name_ext, i + 1);
						mystrcat(name, NAME_SIZE, name_ext);
					}
					ZeroChi2Error(name, 1);
				}
				lead_chi2_pot = *browser;
			}
			total2 += *browser++;
			ppot++;
			pweight++;
			index++;
		}

		//calc the chi2 for the vdW14 interaction
		if (Topology::npair_types > 0)
		{
			ppot = (fabs(RunParams::NB_weight_mode) == 1 ? FNC_POT::vdW14_pot : &FNC_POT::vdW14_tot_pot);
			pweight = RunParams::vdW_weight;//weights for the vdW14 interaction, same as for vdW
			for (i = 0; i < nNB; i++)
			{
				*browser = *ppot / (*pweight * *pweight);
				if (calc_sigma && RunParams::lead_series_ind2 == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "van der Waals 1-4 potential ");
						if (nNB > 1)
						{
							mystrcat(name, NAME_SIZE, " partial ");
							IntToStr(&name_ext, i + 1);
							mystrcat(name, NAME_SIZE, name_ext);
						}
						ZeroChi2Error(name, 1);
					}
					lead_chi2_pot = *browser;
				}
				total2 += *browser++;
				ppot++;
				pweight++;
				index++;
			}

			//calc the chi2 for the Coulomb14 interaction
			ppot = (fabs(RunParams::NB_weight_mode) == 1 ? FNC_POT::Coulomb14_pot : &FNC_POT::Coulomb14_tot_pot);
			pweight = RunParams::Coulomb_weight;//weights for the Coulomb14 interaction, same as for Coulomb
			for (i = 0; i < nNB; i++)
			{
				*browser = *ppot / (*pweight * *pweight);
				if (calc_sigma && RunParams::lead_series_ind2 == index)
				{
					if (fabs(*browser - 0.0) < TOLERANCE)
					{
						mystrcpy(name, NAME_SIZE, "Coulomb 1-4 potential ");
						if (nNB > 1)
						{
							mystrcat(name, NAME_SIZE, " partial ");
							IntToStr(&name_ext, i + 1);
							mystrcat(name, NAME_SIZE, name_ext);
						}
						ZeroChi2Error(name, 1);
					}
					lead_chi2_pot = *browser;
				}
				total2 += *browser++;
				ppot++;
				pweight++;
				index++;
			}
		}
	}
	double bonded_chitotal = 0;//to calc the total bonded chi2, in case of only flexible molecules and RB dihedrals needed to check whether it is positive or not
	//calc the chi2 for the BONDS interaction
	ppot = (Topology::bond_weight_mode == 0 ? &FNC_POT::bond_tot_pot : FNC_POT::bond_pot);
	pweight = Topology::bond_sigma;//weights for the BONDS
	for (i = 0; i < nbonds; i++)
	{
		*browser = *ppot / (*pweight * *pweight);
		if (calc_sigma && lead_ser_pot == index)
		{
			if (fabs(*browser - 0.0) < TOLERANCE)
			{
				mystrcpy(name, NAME_SIZE, "bond streching potential ");
				if (nNB > 1)
				{
					mystrcat(name, NAME_SIZE, " partial ");
					IntToStr(&name_ext, i + 1);
					mystrcat(name, NAME_SIZE, name_ext);
				}
				ZeroChi2Error(name, 1);
			}
			*my_lead_chi2 = *browser;
		}
		bonded_chitotal += *browser;
		*p_pot_chitotal += *browser++;
		ppot++;
		pweight++;
		index++;
	}

	//calc the chi2 for the ANGLES interaction
	ppot = (Topology::angle_weight_mode == 0 ? &FNC_POT::angle_tot_pot : FNC_POT::angle_pot);
	pweight = Topology::angle_sigma;//weights for the ANGLES
	for (i = 0; i < nangles; i++)
	{
		*browser = *ppot / (*pweight * *pweight);
		if (calc_sigma && lead_ser_pot == index)
		{
			if (fabs(*browser - 0.0) < TOLERANCE)
			{
				mystrcpy(name, NAME_SIZE, "angle bending potential ");
				if (nNB > 1)
				{
					mystrcat(name, NAME_SIZE, " partial ");
					IntToStr(&name_ext, i + 1);
					mystrcat(name, NAME_SIZE, name_ext);
				}
				ZeroChi2Error(name, 1);
			}
			*my_lead_chi2 = *browser;
		}
		bonded_chitotal += *browser;
		*p_pot_chitotal += *browser++;
		ppot++;
		pweight++;
		index++;
	}

	//calc the chi2 for the periodic DIHEDRALS interaction
	ppot = (Topology::perdih_weight_mode == 0 ? &FNC_POT::perdihedral_tot_pot : FNC_POT::perdihedral_pot);
	pweight = Topology::perdihedral_sigma;//weights for the perodic DIHEDRALs
	for (i = 0; i < nperdihs; i++)
	{
		*browser = *ppot / (*pweight * *pweight);
		if (calc_sigma && lead_ser_pot == index)
		{
			if (fabs(*browser - 0.0) < TOLERANCE)
			{
				mystrcpy(name, NAME_SIZE, "periodic dihedral potential ");
				if (nNB > 1)
				{
					mystrcat(name, NAME_SIZE, " partial ");
					IntToStr(&name_ext, i + 1);
					mystrcat(name, NAME_SIZE, name_ext);
				}
				ZeroChi2Error(name, 1);
			}
			*my_lead_chi2 = *browser;
		}
		bonded_chitotal += *browser;
		*p_pot_chitotal += *browser++;
		ppot++;
		pweight++;
		index++;
	}

	//calc the chi2 for the harmonic DIHEDRALS interaction
	ppot = (Topology::harmdih_weight_mode == 0 ? &FNC_POT::harmdihedral_tot_pot : FNC_POT::harmdihedral_pot);
	pweight = Topology::harmdihedral_sigma;//weights for the harmonic DIHEDRALs
	for (i = 0; i < nharmdihs; i++)
	{
		*browser = *ppot / (*pweight * *pweight);
		if (calc_sigma && lead_ser_pot == index)
		{
			if (fabs(*browser - 0.0) < TOLERANCE)
			{
				mystrcpy(name, NAME_SIZE, "harmonic dihedral potential ");
				if (nNB > 1)
				{
					mystrcat(name, NAME_SIZE, " partial ");
					IntToStr(&name_ext, i + 1);
					mystrcat(name, NAME_SIZE, name_ext);
				}
				ZeroChi2Error(name, 1);
			}
			*my_lead_chi2 = *browser;
		}
		bonded_chitotal += *browser;
		*p_pot_chitotal += *browser++;
		ppot++;
		pweight++;
		index++;
	}

	//calc the chi2 for the RB DIHEDRALS interaction
	ppot = (Topology::RBdih_weight_mode == 0 ? &FNC_POT::RBdihedral_tot_pot : FNC_POT::RBdihedral_pot);
	pweight = Topology::RBdihedral_sigma;//weights for the RB DIHEDRALs
	for (i = 0; i < nRBdihs; i++)
	{
		*browser = *ppot / (*pweight * *pweight);
		if (calc_sigma && lead_ser_pot == index)
		{
			if (fabs(*browser - 0.0) < TOLERANCE)
			{
				mystrcpy(name, NAME_SIZE, "RB dihedral potential ");
				if (nNB > 1)
				{
					mystrcat(name, NAME_SIZE, " partial ");
					IntToStr(&name_ext, i + 1);
					mystrcat(name, NAME_SIZE, name_ext);
				}
				ZeroChi2Error(name, 1);
			}
			*my_lead_chi2 = *browser;
		}
		bonded_chitotal += *browser;
		*p_pot_chitotal += *browser++;
		ppot++;
		pweight++;
		index++;
	}
	if (RunParams::potential == 0 && Topology::RB_used && bonded_chitotal < 0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "The total chi2 contribution of the bonded interaction is negative probably due to the RB-potential!" << endl;
		cout << "Flexible molecules can only be used without non-bonding potential, if their contribution is positive." << endl;
		cout << "Use non-bonding potetial togethet with flexible molecules in this case, even if its contribution is scaled down" << endl;
		cout << "not to influence the simulation, as this will result in separate potential related chi2 term, which can be negative!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();

	}
	
}		


	


//Calculate the sigma values, if they were given as a percentage 
void ChiSquared::CalcSigma()
{
	int i,index=0;

	double *browser;

	browser=chicomp;//pointer initialisation
	total=0;
	total2=0.0;//the potential related chi2 total contribution, if there is non-bonded potential

	//the g(r) components
	index=0;
	for(i=0;i<ngr;i++)//for each g(r) data set
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				edata.grsigma[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				edata.grsigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=edata.grsigma[i] *edata.grsigma[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
	//the S(Q) components
	for(i=0;i<nsq;i++)//for each S(Q) data set
	{
		//dividing by the standard dev.squared
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				edata.sqsigma[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				edata.sqsigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=edata.sqsigma[i] *edata.sqsigma[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
	//the F(Q) components
	for(i=0;i<nfq;i++)//for each F(Q) data set
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				edata.fqsigma[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				edata.fqsigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=edata.fqsigma[i] *edata.fqsigma[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
	//the F(g) components
	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		if (sigma_percentage[index] < 0)
		{
			if (fabs(*browser - 0.0) > TOLERANCE)//to prevent calculation, if chisquare is 0
				edata.fgsigma[i] = sqrt(*browser / lead_chi2 / fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				edata.fgsigma[i] = lead_sigma * fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser /= edata.fgsigma[i] * edata.fgsigma[i];
		}
		total += *browser;//updating the global result
		browser++;
		index++;
	}
	//the E(k) components
	for(i=0;i<nek;i++)//for each E(k) data set
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				edata.eksigma[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				edata.eksigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=edata.eksigma[i] *edata.eksigma[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
	
	//the Cosine distribution of bond angles constraint
	for(i=0;i<ncos;i++)//for each CosDistr data set
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				CosDistrConst::weights[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				CosDistrConst::weights[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=CosDistrConst::weights[i] * CosDistrConst::weights[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
	//the Coordination number constraint
	for(i=0;i<tot_cc_subconst;i++)//for each Coordination number constraint
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				CoordNumbConst::weights[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				CoordNumbConst::weights[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=CoordNumbConst::weights[i] * CoordNumbConst::weights[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
	//the Average coordination constraint
	for(i=0;i<nacc;i++)//for each average coordination number constraint
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				AvCoordConst::weights[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				AvCoordConst::weights[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=AvCoordConst::weights[i] * AvCoordConst::weights[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
#ifdef _ADVANCED_GEOM_CONST
	//Common neighbour constraint
	for (i = 0; i < ncommonneigh; i++)//for each CommonNeigh constraint
	{
		if (sigma_percentage[index] < 0)
		{
			if (fabs(*browser - 0.0) > TOLERANCE)//to prevent calculation, if chisquare is 0
				CommonNeighConst::weights[i] = sqrt(*browser / lead_chi2 / fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				CommonNeighConst::weights[i] = lead_sigma * fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser /= CommonNeighConst::weights[i] * CommonNeighConst::weights[i];
	}
		total += *browser;//updating the global result
		browser++;
		index++;
	}
	//Second neighbour constraint
	for (i = 0; i < nsecondneigh; i++)//for each SecondNeigh constraint
	{
		if (sigma_percentage[index] < 0)
		{
			if (fabs(*browser - 0.0) > TOLERANCE)//to prevent calculation, if chisquare is 0
				SecondNeighConst::weights[i] = sqrt(*browser / lead_chi2 / fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				SecondNeighConst::weights[i] = lead_sigma * fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser /= SecondNeighConst::weights[i] * SecondNeighConst::weights[i];
		}
		total += *browser;//updating the global result
		browser++;
		index++;
	}
	//BVS constraint
	for (i = 0; i < nbvs; i++)//for each BVS constraint
	{
		if (sigma_percentage[index] < 0)
		{
			if (fabs(*browser - 0.0) > TOLERANCE)//to prevent calculation, if chisquare is 0
				BondValenceSumConst::weights[i] = sqrt(*browser / lead_chi2 / fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				BondValenceSumConst::weights[i] = lead_sigma * fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser /= BondValenceSumConst::weights[i] * BondValenceSumConst::weights[i];
		}
		total += *browser;//updating the global result
		browser++;
		index++;
	}
#endif
#ifdef _LOCAL_INV
	for(i=0;i<nlocint;i++)//for each local inv interval
	{
		if (sigma_percentage[index]<0)
		{
			if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
				RunParams::loc_inv_sigma[i]=sqrt(*browser/lead_chi2/ fabs(sigma_percentage[index]));//dividing by the standard dev.squared
			else
				RunParams::loc_inv_sigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
			*browser/=RunParams::loc_inv_sigma[i] * RunParams::loc_inv_sigma[i];
		}
		total+=*browser;//updating the global result
		browser++;
		index++;
	}
#endif

	//In case of non-bonded/aenet potential total2 is used for all the potential related contributions 
	//as the non-bonded pot can be negative, and its contribution would decrease, or even cancel out the chi2 contribution of the data sets
	//A separate acceptance process decide, whether the move is acceptable based on the change in total2 similarly to the orignal total 
	//If only bonded potential is used for keeping the molecules together, then their contribution will be added to the normal total
	//as those contribution cannot be negative

	int index2=0, redo, start,imoltype,iinttype;
	double my_lead_chi2;
	
	//it is possible, that the lead series of the potential interaction has to be rescaled to the normal data set leading series
	//this has to be done first in case this is not the first potential interaction

	if (RunParams::lead_series_ind2<1 || fabs(RunParams::NB_weight_mode==2))
		redo=0;//it can be done in one go, if there is no non-bonded pot, or the vdW_weight[0] will be used for all the pot contribution's sigma
	else
		redo=2;//first just rescale the leading series then scale the others according to this
		
	//to see, if there is at all zero sigma for any interaction, so the first non-zero sigma for the bond types should be collected
	
	start=index;//to keep the value
	do
	{
		browser=chicomp+start;
		index=start;
		index2=0;
		double *myweight;
		#ifdef _AENET
		if (RunParams::potential==20)
			myweight=&RunParams::aenet_weight;
		else
		#endif
			myweight=RunParams::vdW_weight;
		//the vdW interaction
		for(i=0;i<nNB;i++)//for each vdW-interaction 
		{
			if ((RunParams::lead_series_ind2==index2 && redo==2) || (RunParams::lead_series_ind2!=index2 && redo==1)|| redo==0)
			{
				if (index2 == RunParams::lead_series_ind2)
				{
					my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
					lead_sigma_pot = *myweight;
				}
				else
					my_lead_chi2=lead_chi2_pot;//scaled to the leading pot chi2

				if (sigma_percentage[index]<0)
				{
					if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
						*myweight=sqrt(fabs(*browser/my_lead_chi2/ fabs(sigma_percentage[index])));//dividing by the standard dev.squared
					else
					{
						if (RunParams::lead_series_ind2==index2 && redo==2)//this is the leadseries for the pot
							*myweight=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
						else
							*myweight=lead_sigma_pot*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
					}
					*browser/=*myweight * *myweight;
					if (index2==RunParams::lead_series_ind2)
					{
						lead_chi2_pot=*browser;//reset the lead pot series chi2 to the rescaled value
						lead_sigma_pot=*myweight;
					}
				}
				total2+=*browser;//updating the global result
			}
			browser++;
			if (i<nNB-1)
				myweight++;
			index++;
			index2++;
		}
		if (RunParams::potential == 1)
		{
			//the Coulomb interaction
			for (i = 0; i < nNB; i++)//for each Coulomb-interaction type with different sigma
			{
				if ((RunParams::lead_series_ind2 == index2 && redo == 2) || (RunParams::lead_series_ind2 != index2 && redo == 1) || redo == 0)
				{
					if (index2 == RunParams::lead_series_ind2)
					{
						my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
						lead_sigma_pot = RunParams::Coulomb_weight[i];
					}
					else
						my_lead_chi2 = lead_chi2_pot;//scaled to the leading pot chi2
					if (sigma_percentage[index] < 0 || fabs(sigma_percentage[index] - 0.0) < TOLERANCE)
					{
						if (fabs(RunParams::NB_weight_mode)!= 2 && sigma_percentage[index] < 0)
						{
							if (fabs(*browser - 0.0) > TOLERANCE)//to prevent calculation, if chisquare is 0
								RunParams::Coulomb_weight[i] = sqrt(fabs(*browser / my_lead_chi2 / fabs(sigma_percentage[index])));//dividing by the standard dev.squared
							else
							{
								if (RunParams::lead_series_ind2 == index2 && redo == 2)//this is the lead series for the pot
									RunParams::Coulomb_weight[i] = lead_sigma * fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
								else
									RunParams::Coulomb_weight[i] = lead_sigma_pot * fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
							}
						}
						else
						{
							RunParams::Coulomb_weight[0] = RunParams::vdW_weight[0];
						}
						*browser /= RunParams::Coulomb_weight[i] * RunParams::Coulomb_weight[i];
						if (index2 == RunParams::lead_series_ind2)
						{
							lead_chi2_pot = *browser;//reset the lead pot series chi2 to the rescaled value
							lead_sigma_pot = RunParams::Coulomb_weight[i];
						}

					}
					total2 += *browser;//updating the global result
				}
				browser++;
				index++;
				index2++;
			}
			if (Topology::npair_types > 0)
			{
				//the sigma is already recalculated, if it was necessary, only calculate the chi2, if necessary and the total
				for (i = 0; i < nNB; i++)//for each 1-4 vdW-interaction 
				{
					if (redo < 2)
					{
						if (sigma_percentage[index] < 0)
							*browser /= RunParams::vdW_weight[i] * RunParams::vdW_weight[i];

						total2 += *browser;//updating the global result
					}
					browser++;
					index++;
					index2++;
				}
				//the Coulomb interaction
				for (i = 0; i < nNB; i++)//for each Coulomb-interaction type with different sigma
				{
					if (redo < 2)
					{
						if (sigma_percentage[index] < 0 || fabs(sigma_percentage[index] - 0.0) < TOLERANCE)
							*browser /= RunParams::Coulomb_weight[i] * RunParams::Coulomb_weight[i];

						total2 += *browser;//updating the global result
					}
					browser++;
					index++;
					index2++;
				}
			}
		}
	
		//the BONDs
		for(i=0;i<nbonds;i++)//for each different BOND type
		{
			if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
			{
				index++;
				index2++;
				browser++;
				continue; //setting the sigma for the zero sigma series will be done later, when all the other sigma is available
			}
			if ((RunParams::lead_series_ind2==index2 && redo==2) || (RunParams::lead_series_ind2!=index2 && redo==1)|| redo==0)
			{
				if (RunParams::potential == 0 || index2 == RunParams::lead_series_ind2)
				{
					my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
					lead_sigma_pot = Topology::bond_sigma[i];
				}
				else
					my_lead_chi2=lead_chi2_pot;//scaled to the leading pot chi2
				if (sigma_percentage[index]<0)
				{
					if (fabs(RunParams::NB_weight_mode)!=2)
					{
						if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
							Topology::bond_sigma[i]=sqrt(fabs(*browser/my_lead_chi2/ fabs(sigma_percentage[index])));//dividing by the standard dev.squared
						else
						{
							if (RunParams::lead_series_ind2==index2 && redo==2)//this is the lead series for the pot
								Topology::bond_sigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
							else
								Topology::bond_sigma[i]=lead_sigma_pot*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
						}
					}
					else
						Topology::bond_sigma[i]=RunParams::vdW_weight[0];

					*browser/=Topology::bond_sigma[i] * Topology::bond_sigma[i];
					if (index2==RunParams::lead_series_ind2)
					{
						lead_chi2_pot=*browser;//reset the lead pot series chi2 to the rescaled value
						lead_sigma_pot=Topology::bond_sigma[i];
					}
				}
				*p_pot_chitotal+=*browser;//updating the global result
			}
			
			browser++;
			index++;
			index2++;
		}
	
		//the ANGLEs
		for(i=0;i<nangles;i++)//for each different ANGLE type
		{
			if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
			{
				index++;
				index2++;
				browser++;
				continue; //setting the sigma for the zero sigma series will be done later, when all the other sigma is available
			}
			if ((RunParams::lead_series_ind2 == index2 && redo == 2) || (RunParams::lead_series_ind2 != index2 && redo == 1) || redo == 0)
			{
				if (RunParams::potential == 0 || index2 == RunParams::lead_series_ind2)
				{
					my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
					lead_sigma_pot = Topology::angle_sigma[i];
				}
				else
					my_lead_chi2=lead_chi2_pot;//scaled to the leading pot chi2
				if (sigma_percentage[index]<0)
				{
					if (fabs(RunParams::NB_weight_mode) != 2)
					{
						if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
							Topology::angle_sigma[i]=sqrt(fabs(*browser/my_lead_chi2/ fabs(sigma_percentage[index])));//dividing by the standard dev.squared
						else
						{
							if (RunParams::lead_series_ind2==index2 && redo==2)//this is the lead series for the pot
								Topology::angle_sigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
							else
								Topology::angle_sigma[i]=lead_sigma_pot*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
						}
					}
					else
						Topology::angle_sigma[i]=RunParams::vdW_weight[0];
					*browser/=Topology::angle_sigma[i] * Topology::angle_sigma[i];
					if (index2==RunParams::lead_series_ind2)
					{	lead_chi2_pot=*browser;//reset the lead pot series chi2 to the rescaled value
						lead_sigma_pot=Topology::angle_sigma[i];
					}
				}
				*p_pot_chitotal+=*browser;//updating the global result
			}
			browser++;
			index++;
			index2++;
		}

		//the periodic DIHEDRALs
		for(i=0;i<nperdihs;i++)//for each different periodic DIHEDRAL type
		{
			if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
			{
				index++;
				index2++;
				browser++;
				continue; //setting the sigma for the zero sigma series will be done later, whne all the other sigma is available
			}
			if ((RunParams::lead_series_ind2==index2 && redo==2) || (RunParams::lead_series_ind2!=index2 && redo==1)|| redo==0)
			{
				if (RunParams::potential == 0 || index2 == RunParams::lead_series_ind2)
				{
					my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
					lead_sigma_pot = Topology::perdihedral_sigma[i];
				}
				else
					my_lead_chi2=lead_chi2_pot;//scaled to the leading pot chi2
				if (sigma_percentage[index]<0)
				{
					if (fabs(RunParams::NB_weight_mode) != 2)
					{
						if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
							Topology::perdihedral_sigma[i]=sqrt(fabs(*browser/my_lead_chi2/ fabs(sigma_percentage[index])));//dividing by the standard dev.squared
						else
						{
							if (RunParams::lead_series_ind2==index2 && redo==2)//this is the leadbseries for the pot
								Topology::perdihedral_sigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
							else
								Topology::perdihedral_sigma[i]=lead_sigma_pot*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
						}
					}
					else
						Topology::perdihedral_sigma[i]=RunParams::vdW_weight[0];
					*browser/=Topology::perdihedral_sigma[i] * Topology::perdihedral_sigma[i];
					if (index2==RunParams::lead_series_ind2)
					{
						lead_chi2_pot=*browser;//reset the lead pot series chi2 to the rescaled value
						lead_sigma_pot=Topology::perdihedral_sigma[i];
					}
				}
				*p_pot_chitotal+=*browser;//updating the global result
			}
			browser++;
			index++;
			index2++;
		}

		//the harmonic DIHEDRALs
		for(i=0;i<nharmdihs;i++)//for each different periodic DIHEDRAL type
		{
			if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
			{
				index++;
				index2++;
				browser++;
				continue; //setting the sigma for the zero sigma series will be done later, whne all the other sigma is available
			}
			if ((RunParams::lead_series_ind2==index2 && redo==2) || (RunParams::lead_series_ind2!=index2 && redo==1)|| redo==0)
			{
				if (RunParams::potential == 0 || index2 == RunParams::lead_series_ind2)
				{
					my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
					lead_sigma_pot = Topology::harmdihedral_sigma[i];
				}
				else
					my_lead_chi2=lead_chi2_pot;//scaled to the leading pot chi2
				if (sigma_percentage[index]<0)
				{
					if (fabs(RunParams::NB_weight_mode) != 2)
					{
						if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
							Topology::harmdihedral_sigma[i]=sqrt(fabs(*browser/my_lead_chi2/ fabs(sigma_percentage[index])));//dividing by the standard dev.squared
						else
						{
							if (RunParams::lead_series_ind2==index2 && redo==2)//this is the leadbseries for the pot
								Topology::harmdihedral_sigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
							else
								Topology::harmdihedral_sigma[i]=lead_sigma_pot*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
						}
					}
					else
						Topology::harmdihedral_sigma[i]=RunParams::vdW_weight[0];
					*browser/=Topology::harmdihedral_sigma[i] * Topology::harmdihedral_sigma[i];
					if (index2==RunParams::lead_series_ind2)
					{
						lead_chi2_pot=*browser;//reset the lead pot series chi2 to the rescaled value
						lead_sigma_pot=Topology::harmdihedral_sigma[i];
					}
				}
				*p_pot_chitotal+=*browser;//updating the global result
			}
			browser++;
			index++;
			index2++;
		}

		//the RB DIHEDRALs
		for(i=0;i<nRBdihs;i++)//for each different RB DIHEDRAL type
		{
			if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
			{
				index++;
				index2++;
				browser++;
				continue; //setting the sigma for the zero sigma series will be done later, whne all the other sigma is available
			}
			if ((RunParams::lead_series_ind2==index2 && redo==2) || (RunParams::lead_series_ind2!=index2 && redo==1)|| redo==0)
			{
				if (RunParams::potential == 0 || index2 == RunParams::lead_series_ind2)
				{
					my_lead_chi2 = lead_chi2;//this is the leading chi2 for the potentials, but this will be scaled, if it is scaled to the normal data sets
					lead_sigma_pot = Topology::RBdihedral_sigma[i];
				}
				else
					my_lead_chi2=lead_chi2_pot;//scaled to the leading pot chi2
				if (sigma_percentage[index]<0)
				{
					if (fabs(RunParams::NB_weight_mode) != 2)
					{
						if (fabs(*browser-0.0)>TOLERANCE)//to prevent calculation, if chisquare is 0
							Topology::RBdihedral_sigma[i]=sqrt(fabs(*browser/my_lead_chi2/ fabs(sigma_percentage[index])));//dividing by the standard dev.squared
						else
						{
							if (RunParams::lead_series_ind2==index2 && redo==2)//this is the leadbseries for the pot
								Topology::RBdihedral_sigma[i]=lead_sigma*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
							else
								Topology::RBdihedral_sigma[i]=lead_sigma_pot*fabs(sigma_percentage[index]);//sets the sigma as the percentage of the leading series' sigma
						}
					}
					else
						Topology::RBdihedral_sigma[i]=RunParams::vdW_weight[0];
					*browser/=Topology::RBdihedral_sigma[i] * Topology::RBdihedral_sigma[i];
					if (index2==RunParams::lead_series_ind2)
					{
						lead_chi2_pot=*browser;//reset the lead pot series chi2 to the rescaled value
						lead_sigma_pot=Topology::RBdihedral_sigma[i];
					}
				}
				*p_pot_chitotal+=*browser;//updating the global result
			}
			browser++;
			index++;
			index2++;
		}

		redo--;
	}
	while (redo>0);

	//to start with the bonds
	switch (RunParams::potential)
	{
	case 0:
	{
		index = start;
	}
	case 1:
	{
		index = start + (Topology::npair_types > 0 ? 4 * nNB : 2 * nNB);
		break;
	}
	case 10:
	{
		index = start + nNB;
		break;
	}
	//not needed for aenet
	}
	
	browser=chicomp+index;
	//only calculate, if there is zero sigma, and weight mode is 1, (separate sigmas are used for each bond type)
	if (Topology::zero_bond_sigma!=nullptr)//not created for aenet
	{
		if (Topology::zero_bond_sigma[0] && Topology::bond_weight_mode)// array is created if topology is read, even if there are no bonds
		{
			i=0;
			for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
			{
				for (iinttype=0;iinttype<Topology::nbondtypes_per_type[imoltype];iinttype++)
				{
					if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
					{
						Topology::bond_sigma[Topology::cumul_bondtypes[imoltype]+iinttype]=Topology::bond_sigma[Topology::zero_bond_sigma[imoltype+1]];
						*browser/=Topology::bond_sigma[i] * Topology::bond_sigma[i];
						*p_pot_chitotal+=*browser;//updating the global result	
					}
					index++;
					i++;
					browser++;
				}
			}
		}
		else
		{
			index+=nbonds;
			browser+=nbonds;
		}
	}
	
	if (Topology::zero_angle_sigma)
	{
		if (Topology::angle_weight_mode)
		{
			i=0;
			for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
			{
				for (iinttype=0;iinttype<Topology::nangletypes_per_type[imoltype];iinttype++)
				{
					if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
					{
						Topology::angle_sigma[Topology::cumul_angletypes[imoltype]+iinttype]=Topology::bond_sigma[Topology::zero_bond_sigma[imoltype+1]];
						*browser/=Topology::angle_sigma[i] * Topology::angle_sigma[i];
						*p_pot_chitotal+=*browser;//updating the global result	
					}
					index++;
					i++;
					browser++;
				}
			}	
		}
		else
		{
			//only could be collapsed, if the bond sigmas for each molecule type are the same, according to it this has to be set
			Topology::angle_sigma[0]=Topology::bond_sigma[Topology::zero_bond_sigma[1]];//only the first sigma is set here, this will be copied to the others later
			*browser/=Topology::angle_sigma[0] * Topology::angle_sigma[0];
			*p_pot_chitotal+=*browser;//updating the global result	
			index++;
			browser++;
		}

	}
	else
	{
		index+=nangles;
		browser+=nangles;
	}

	if (Topology::zero_perdih_sigma)
	{
		if (Topology::perdih_weight_mode)
		{
			i=0;
			for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
			{
				for (iinttype=0;iinttype<Topology::nperdihtypes_per_type[imoltype];iinttype++)
				{
					if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
					{
						Topology::perdihedral_sigma[Topology::cumul_perdihtypes[imoltype]+iinttype]=Topology::bond_sigma[Topology::zero_bond_sigma[imoltype+1]];
						*browser/=Topology::perdihedral_sigma[i] * Topology::perdihedral_sigma[i];
						*p_pot_chitotal+=*browser;//updating the global result	
					}
					index++;
					i++;
					browser++;
				}
			}
		}
		else
		{
			//only could be collapsed, if the bond sigmas for each molecule type are the same, according to it this has to be set
			Topology::perdihedral_sigma[0]=Topology::bond_sigma[Topology::zero_bond_sigma[1]];//only the first sigma is set here, this will be copied to the othe later
			*browser/=Topology::perdihedral_sigma[0] * Topology::perdihedral_sigma[0];
			*p_pot_chitotal+=*browser;//updating the global result	
			index++;
			browser++;
		}
	}
	else
	{
		index+=nperdihs;
		browser+=nperdihs;
	}

	if (Topology::zero_harmdih_sigma)
	{
		if (Topology::harmdih_weight_mode)
		{
			i=0;
			for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
			{
				for (iinttype=0;iinttype<Topology::nharmdihtypes_per_type[imoltype];iinttype++)
				{
					if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
					{
						Topology::harmdihedral_sigma[Topology::cumul_harmdihtypes[imoltype]+iinttype]=Topology::bond_sigma[Topology::zero_bond_sigma[imoltype+1]];
						*browser/=Topology::harmdihedral_sigma[i] * Topology::harmdihedral_sigma[i];
						*p_pot_chitotal+=*browser;//updating the global result	
					}
					index++;
					i++;
					browser++;
				}
			}
		}
		else
		{
			//only could be collapsed, if the bond sigmas for each molecule type are the same, according to it this has to be set
			Topology::harmdihedral_sigma[0]=Topology::bond_sigma[Topology::zero_bond_sigma[1]];//only the first sigma is set here, this will be copied to the othe later
			*browser/=Topology::harmdihedral_sigma[0] * Topology::harmdihedral_sigma[0];
			*p_pot_chitotal+=*browser;//updating the global result	
			index++;
			browser++;
		}
	}
	else
	{
		index+=nharmdihs;
		browser+=nharmdihs;
	}
	if (Topology::zero_RBdih_sigma)
	{
		if (Topology::RBdih_weight_mode)
		{
			i=0;
			for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
			{
				for (iinttype=0;iinttype<Topology::nRBdihtypes_per_type[imoltype];iinttype++)
				{
					if (fabs(sigma_percentage[index]-0.0)<TOLERANCE)
					{
						Topology::RBdihedral_sigma[Topology::cumul_RBdihtypes[imoltype]+iinttype]=Topology::bond_sigma[Topology::zero_bond_sigma[imoltype+1]];
						*browser/=Topology::RBdihedral_sigma[i] * Topology::RBdihedral_sigma[i];
						*p_pot_chitotal+=*browser;//updating the global result	
					}
					index++;
					i++;
					browser++;
				}
			}
		}
		else
		{
			//only could be collapsed, if the bond sigmas for each molecule type are the same, according to it this has to be set
			Topology::RBdihedral_sigma[0]=Topology::bond_sigma[Topology::zero_bond_sigma[1]];//only the first sigma is set here, this will be copied to the othe later
			*browser/=Topology::RBdihedral_sigma[0] * Topology::RBdihedral_sigma[0];
			*p_pot_chitotal+=*browser;//updating the global result	
			index++;
			browser++;
		}
	}
	else
	{
		index+=nRBdihs;
		browser+=nRBdihs;
	}

	//setting the sigma for the other bond types, only needed for History Topology parameters save
	if (Topology::bond_weight_mode==0)
	{
		for (i=1;i<Topology::nbond_types;i++)
			Topology::bond_sigma[i]=Topology::bond_sigma[0];
	}
	//setting the sigma for the other angle types, only needed for History Topology parameters save
	if (Topology::angle_weight_mode==0)
	{
		for (i=1;i<Topology::nangle_types;i++)
			Topology::angle_sigma[i]=Topology::angle_sigma[0];
	}
	//setting the sigma for the other periodic dihedral types, only needed for History Topology parameters save
	if (Topology::perdih_weight_mode==0)
	{
		for (i=1;i<Topology::nperdihedral_types;i++)
			Topology::perdihedral_sigma[i]=Topology::perdihedral_sigma[0];
	}
	//setting the sigma for the other harmonic dihedral types, only needed for History Topology parameters save
	if (Topology::harmdih_weight_mode==0)
	{
		for (i=1;i<Topology::nharmdihedral_types;i++)
			Topology::harmdihedral_sigma[i]=Topology::harmdihedral_sigma[0];
	}
	//setting the sigma for the other RB dihedral types, only needed for History Topology parameters save
	if (Topology::RBdih_weight_mode==0)
	{
		for (i=1;i<Topology::nRBdihedral_types;i++)
			Topology::RBdihedral_sigma[i]=Topology::RBdihedral_sigma[0];
	}

};

//calculate the sums for chi2 and the derivates for non-lin regression for the ith F(Q) data set
//the parameters are for the nonline reg.
//1: a
//2: alpha
//3: b
//FD-s are the first derivates, SD_s are the second derivates
void ChiSquared::CalcDerivates(CalcData &calc,int i)
{
	int j,index;
	double *pcalc, *pexp, *pq, *pcompt;
	double comptQQ;// compt *Q *Q, (a*Iexp_j+b-Scalc_j-A_j)/B_j
	index = ngr + nsq + i;
	//unfortunately all the sums contain alpha, so have to be recalculated at each step
	if (log_nlr_steps)
		logfile << "\nCalculating derivates for F(Q) data set " << i + 1 << endl;
	
	s[index] = 0.0;
	fs[index] = 0.0; 
	C[i] = 0.0;//sum of B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j)
	CC[i] = 0.0;//sum of (B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j))^2
	sC[i] = 0.0;// sum of(I(Q)_calc * (B_j * Q_j * Q_j * exp(-alpha * Q_j * Q_j))
	fC[i] = 0.0;//sum of (experimental*( B_j*Q_j*Q_j*exp(-alpha*Q_j*Q_j))
	
	switch (fit_index[index])
	{
		case 35://nonlinear a,b, alpha fit,  case 33 and 35 are connected, calculated in optimal order
			pcalc = calc.iqfinder[i];
			for (j = 0; j < *(edata.fqused + i); j++)//for all data points
				s[index] += *pcalc++;
		case 33://nonlinear a, alpha fit, b=0
			pcalc = calc.iqfinder[i];
			pexp = edata.fq_ffinder[i];//pointer initialisation
			for (j = 0; j < *(edata.fqused + i); j++)//for all data points
				fs[index] += *pcalc++ * *pexp++;
		case 32://nonlinear alpha fit, a,b custom
			pcalc = calc.iqfinder[i];
			pexp = edata.fq_ffinder[i];//pointer initialisation
			pq = edata.fq_qfinder[i];
			pcompt = calc.comptfinder[i];
			for (j = 0; j < *(edata.fqused + i); j++)//for all data points
			{
				comptQQ = *pcompt++ * *pq * *pq;
				C[i] += comptQQ;
				CC[i] += comptQQ * comptQQ;
				sC[i] += *pcalc++ * comptQQ;
				fC[i] += *pexp++ * comptQQ++;
				pq++;
			}
			break;
	}

	switch (fit_index[index])
	{
		case 35://derivates containing db (and the others)
			FD3 = (a[index] * f[index] + b[index] * edata.fqused[i] - s[index]);//dchi2/db
			SD13 = f[index];//dchi2/da*dchi2/db
			SD23 = C[i];//dchi2/dalpha/db
			SD33 = edata.fqused[i];//dchi2/db*dchi2/db
		case 33://derivates containing a (and alpha)
			FD1 = (a[index] * ff[index] + b[index] * f[index] - fs[index]);//dchi2/da
			SD11 = ff[index];//dchi2/da*dchi2/da
			SD12 = fC[i];//dchi2/da*dchi2/dalpha
		case 32:
			FD2 = (a[index] * fC[i] + b[index] * C[i] - sC[i]);//dchi2/dalpha
			SD22 = CC[i];//dchi2/dalpha*dchi2/dalpha
			break;
	}
};

//calculating the parameters with non-linear regression
//The formula used without sigma: chi2=summa_k[a*I_expt(Q_k)+b-I_calc(Q_k)]^2, where
//I_calc(Q_k)=F_tot(Q_k)+A(Q_k)+B(Q_k)exp(-alpha*Q_k^2)
//this has to be solved with iterative non-linear regression using Levenberg-Marquardt method
//https://e-maxx.ru/bookz/files/numerical_recipes.pdf  15.5.2
//gives back 0 if all the data sets converged, 1 otherwise
//it is checked, whether the denominators in the d_X calculations are non-zero, and try to recover from it
//skip 0: no problem
//skip 1: derivates has to be recalculated after changing lambda
//skip 2: only lambda changes
int ChiSquared::NonLinReg(CalcData &calc)
{
	int i,j,index,it;
	int fail_flag;
	int result = 0, skip;
	double delta_a, delta_b, delta_alpha=0;
	double *pcalc, *pexp;
	double D,E,F,G;
	double chi2=0,chi2_old,chi2_start;
	double old_a, old_b, old_alpha;
	double best_a=1, best_b=0, best_alpha=0, chi2_best;
	double lambda;
	bool cont;

	for (i = 0; i < nfq; i++)
	{
		
		index = ngr + nsq + i;
		if (ExptsData::fqfitIQ[i] && fit_index[index] > 31)//non-lin
		{
			if (!(Move::nomove) || (Move::makefprimeshift && edata.fqAXS[i] > 0) || (Move::makemucorr && edata.fqIQbackgcorr[i]))
			{
				cont = true;
				chi2_best = 1e10;
				lambda = RunParams::lambda_nonlin;
				delta_a = 0;//might be needed for log output
				delta_b = 0;

				//unfortunately all the sums contain alpha, so have to be recalculated at each step
				CalcDerivates(calc, i);//calculates sums and the derivates

				//calculate chi2

				//chi2_old = pow(a[index], 2) * ff[index] + ss[index] + pow(b[index], 2) * edata.fqused[i] + 2 * (a[index] * b[index] * f[index] - a[index] * fs[index] - b[index] * s[index]);
				//if (chi2_old <= 0.0)//close to zero, recalculate
				//{
				//as later the new chi2 is calculated with the square of the difference method for the sake of continuation runs this will be used here as well 
				chi2_old = 0.0;
				pcalc = calc.iqfinder[i];
				pexp = edata.fq_ffinder[i];//pointer initialisation
				for (j = 0; j < edata.fqused[i]; j++)
					chi2_old += pow(a[index] * *pexp++ + b[index] - *pcalc++, 2);
				//}
				if (RunParams::continuation == -1)//do not calculate in case of the initial chi2 calculation of a continuation run, as the renorm parameters are read from state file
				{
					chicomp[index] = chi2_old;
					continue;
				}
				chi2_start = chi2_old;
				if (chi2_start < chi2_best)
				{
					chi2_best = chi2_start;
					best_a = a[index];
					best_b = b[index];
					best_alpha = alpha[i];

				}

				if (log_nlr_steps)
					logfile << "start chi2 " << chi2_old << " a " << a[index] << " alpha " << alpha[i] << " b " << b[index] << endl;
				old_a = a[index];//preserving the old values
				old_b = b[index];
				old_alpha = alpha[i];
				//now iterate 
				it = 1;
				chi2 = chi2_start;
				do
				{
					fail_flag = 0;
					skip = 0;//this is for preventing division by zero
					switch (fit_index[index])
					{
					case 32:
						if (fabs(SD22) < TOLERANCE2)//SD22 is zero, not depending on lambda, do not calculate, leave previous value
						{
							skip = 1;
							fail_flag += 4;
						}
						else
							delta_alpha = FD2 / SD22 / (1 + lambda);
						break;
					case 33://fit a,alpha
						if (fabs(SD12 * SD12 - SD11 * (1 + lambda) * SD22 * (1 + lambda)) < TOLERANCE2)// zero, depending on lambda
						{
							if (SD12 * SD12 < TOLERANCE2)//it is close-zero in both term, lambda change itself cannot help, the previous delta_alpha will be used
							{
								skip = 1; //old delta_alpha will be used
								fail_flag += 4;
							}
							else
							{
								skip = 2;
								break;
							}
						}
						else
							delta_alpha = (FD1 * SD12 - FD2 * SD11 * (1 + lambda)) / (SD12 * SD12 - SD11 * (1 + lambda) * SD22 * (1 + lambda));
						if (fabs(SD11) < TOLERANCE2)//SD11 is zero, not depending on lambda
						{
							skip = 1;//old delta_a will be used
							fail_flag += 1;
						}
						else
							delta_a = (FD1 - SD12 * delta_alpha) / SD11 / (1 + lambda);
						a[index] -= delta_a;
						break;
					case 35://fit a,b,alpha
						G = SD23 * SD11 * (1 + lambda) - SD13 * SD12;
						F = FD3 * SD11 * (1 + lambda) - FD1 * SD13;
						if (fabs(G)<TOLERANCE2)//G is zero, depending on lambda, new D cannot be calculated, it is needed for all the three params
						{
							if (fabs(SD13 * SD12) < TOLERANCE2) // it is close - zero in both term, lambda change itself cannot help
							{
								skip = 1;
								fail_flag = 7;
								a[index] -= delta_a;//use the old delta values again, alpha is set later
								b[index] -= delta_b;
							}
							else
								skip = 2;
							break;
						}
						D = (SD33 * (1 + lambda) * SD11 * (1 + lambda) - SD13 * SD13) / G;
						E = SD11 * (1 + lambda) * SD22 * (1 + lambda) - SD12 * SD12;
						if (fabs(G - D * E) < TOLERANCE2)// zero, depending on lambda
						{
							skip = 2;
							break;
						}
						delta_b = (SD11 * (1 + lambda) * FD2 - FD1 * SD12 - F / G * E) / (G - D * E);
						delta_alpha = F / G - delta_b * D;
						delta_a = (FD1 - delta_alpha * SD12 - delta_b * SD13) / SD11 / (1 + lambda);
						a[index] -= delta_a;
						b[index] -= delta_b;
						break;
					}
					if (skip<2)//paramaters change, although delta value only in case of skip=0
					{
						alpha[i] -= delta_alpha;
						//as alpha has changed, recalculate the I(Q)
						calc.RecalcIQ(i);

						//calculate new chi2, the sums will not be recalculated here!
						chi2 = 0.0;
						pcalc = calc.iqfinder[i];
						pexp = edata.fq_ffinder[i];//pointer initialisation
						for (j = 0; j < edata.fqused[i]; j++)
							chi2 += pow(a[index] * *pexp++ + b[index] - *pcalc++, 2);
						if (log_nlr_steps)
						{
							logfile << "n_iteration " << it;
							if (skip)
								logfile << " failed with division by zero using old delta value for";
							if (fail_flag & 1)
								logfile << " a";
							if (fail_flag & 2)
								logfile << " b";
							if (fail_flag & 4)
								logfile << " alpha";
							logfile << " chi2 " << chi2 << " chi2_old " << chi2_old << " chi2_old - chi2 " << chi2_old - chi2 << " d_a " << delta_a << " d_alpha " << delta_alpha << " delta_b " << delta_b << " a " << a[index] << " alpha " << alpha[i] << " b " << b[index] << " lambda " << lambda << endl;
						}
						if (chi2 < chi2_best)//preserv the lowest
						{
							chi2_best = chi2;
							best_a = a[index];
							best_b = b[index];
							best_alpha = alpha[i];

						}
					}
					else
					{
						if (log_nlr_steps)
							logfile << "n_iteration "<< it <<" failed with division by zero changing only labda, chi2 " << chi2 << " chi2_old " << chi2_old << " chi2_old-chi2 " << chi2_old - chi2 << " d_a " << delta_a << " d_alpha " << delta_alpha << " delta_b " << delta_b << " a " << a[index] << " alpha " << alpha[i] << " b " << b[index] << " lambda " << lambda << endl;
					}
					if (skip<2 && (fabs(chi2_old - chi2) < RunParams::epsilon_nonlin || fabs((chi2_old - chi2) / chi2) < RunParams::epsilon_nonlin))
					{//for skip 2 there was no chi2 calc
						//do not continue any further, the chi2 change is below the limit 
						cont = false;

					}
					else
					{
						if (chi2 >= chi2_old  || skip==2)//increased, or skipped due to zero division,set back the old parameters, only lambda will change, no need to recalc derivates
						{
							lambda *= RunParams::factor_nonlin;
							alpha[i] = old_alpha;//set it back
							if (fit_index[index] > 32)
								a[index] = old_a;//set it back
							if (fit_index[index] == 35)
								b[index] = old_b;
						}
						else//decreased, copy params and chi2 to old
						{
							lambda /= RunParams::factor_nonlin;
							old_alpha = alpha[i];
							if (fit_index[index] > 32)
								old_a = a[index];
							if (fit_index[index] == 35)
								old_b = b[index];
							chi2_old = chi2;

							CalcDerivates(calc, i);//recalculate the new sums and the derivates
						}
					}
					it++;
				} while (cont && it < RunParams::niter_nonlin);
				//calcdat will contain the last I(Q) and chi2 is the last chi2
				if ((chi2 > chi2_best && chi2 - chi2_best > RunParams::epsilon_nonlin) || skip>0)//only significantly smaller chi2_best should be used
				{
					chi2 = chi2_best;
					a[index] = best_a;
					b[index] = best_b;
					alpha[i] = best_alpha;
					calc.RecalcIQ(i);
					if (log_nlr_steps)
						logfile << "Smallest chi2 during iteration replaced the last one" << endl;

				}
				if (cont)//terminated without convergence
				{
					tot_failed_nlr[i]++;
					result = 1;
					return result;//there is no point to try with the other data sets, if there is any, this step will not be accepted anyhow
				}
				else
				{
					//regression was successful
					chicomp[index] = chi2;
					if (chi2 / chi2_start < maxdec_nlr_chi2fr[i])
						maxdec_nlr_chi2fr[i] = chi2 / chi2_start;//for statistic
					last_nlr_chi2fr[i] = chi2 / chi2_start;
					if (log_nlr_steps)
						logfile << "end " << " chi2 " << chicomp[index] << " a " << a[index] << " alpha " << alpha[i] << " b " << b[index] << endl;
				}
			}//there was change
		}//fit I(Q)
	}//data set cycle
	return result;
};

//in case of Xray I(Q) fitting only some of the possible combinations is meaningful and therefore implemented, otherwise this message is given
void ChiSquared::MissingRenorCombError(int i)
{
	cout << "\n*****ERROR*****" << endl;
	cout << "Unimplemented renormalization combination for " << i + 1 << ". X-ray data set I(Q) fitting:" << endl;
	cout.setf(std::ios::boolalpha);
	cout << "\tFit parameter a    : " << (bool)ExptsData::fqrenorm[i] << endl;
	cout << "\tFit parameter b    : " << (bool)ExptsData::fqlinear[i] << endl;
	cout << "\tFit parameter alpha: " << (bool)ExptsData::fqrenalpha[i] << endl;
	cout.unsetf(std::ios::boolalpha);
	cout << "The implemented combinations are:" << endl;
	cout << "\tFitted parameter a, with constant custom parameters b and alpha" << endl;
	cout << "\tFitted parameter alpha, with constant custom parameters a and b" << endl;
	cout << "\tFitted parameters a and b , with alpha=0" << endl;
	cout << "\tFitted parameters a and alpha, with b=0" << endl;
	cout << "\tFitted parameters a, b and alpha" << endl;
	cout << "Choose form the above possibilities, and try again!" << endl;
	cout << "Cannot run this way, exiting..." << endl;
	CleanExit();
};