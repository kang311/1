//source History.cpp
//Last changed 27.02.2023

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "classes2.h"

int			 History::configtype;//1,3 if the coordinates come from the binary .bcf file, 0 if from the text .cfg file 
int			 History::swap;//boolean, whether there are swaps
int			 History::buffsize;//buffer size
int			 History::chisize;//number of chi_squared parameters to save
int			 History::stepratio;//how often should values be kept
int			 History::fill;//current number of buffered steps
longint		*History::potaccepted=0;//number of accepted moves
int			 History::nNB=0;//number of non-bonded interactions with different sigma
int			 History::chi_offset1=0;//last element before 1-4 poptential params in chicomp array
int			 History::chi_offset2=0;//first element after 1-4 poptential params in chicomp array
int			 History::chisize2=0;//number of standard, (not potential related) chi squared parameters to save
int			 History::datapoint_tot2=0;//total number of data points in case of non-bonded pot present for the bonded and non-bonded interacions
#ifdef _LOCAL_INV
int			 History::nlocint;//number of local invariance intervals
int			*History::loc_npoints;//number of data points for local invariance (depend on the calculation mode)
#endif
int 		 History::datapoint_tot;//total number of data points
longint		*History::generated;//number of generated moves
longint		*History::swapgen;//number of generated swaps
longint		*History::E0shiftgen;//number of generated shifts
longint		*History::fprimeshiftgen;//number of generated f' shifts
longint		*History::mucorrgen;//number of generated I(Q) mu corrections
longint		*History::tried;//number of tried moves
longint		*History::goodnonlinreg;//number of moves ended in good nonlin reg
longint		*History::accepted;//number of accepted moves
longint		*History::swapacc;//number of accepted swaps
longint		*History::E0shiftacc;//number of accepted shifts
longint		*History::fprimeshiftacc;//number of accepted f' shifts
longint		*History::mucorracc;//number of accepted I(Q) mu corrections
longint		History::last_gen;//number of generated moves since the last buffering, cannot use rundat's as displaying and history buffering happens at differnt times
longint		History::last_tried;//number of tried moves since the last buffering
longint		History::last_acc;//number of accepted moves since the last buffering

longint	 History::calc_method;//indicating the calculation method of the normal chi2

double		*History::time;//time values
double		*History::ir1;//instant ratio of number of tried/generated moves
						//over the last nstep moves (defined by the run parameters)
double		*History::ir2;//instant ratio of number of accepted/tried moves
						//over the last nstep moves (defined by the run parameters)
double		*History::chicomp;//chi squared components (one per data set + 1 per CC+
														//the total)
														//one per velocity dist. const+the total)
//constructor
History::History()
{	
	//this class defines a `run history'
	//which allows to save the evolution of the RMC++ run
	

	//number of chi squared parameters to save
	//the total + the values for the expts data sets +
	//the values for the coordination constraints 
	
	//declaration of arrays
	SetArraysize(&generated,buffsize,"generated","History::History");//array of number of generated moves
	SetArraysize(&tried,buffsize,"tried","History::History");//array of number of tried moves
	SetArraysize(&accepted,buffsize,"accepted","History::History");//array of number of accepted moves
	SetArraysize(&potaccepted,buffsize,"potaccepted","History::History");//array of number of accepted moves
	SetArraysize(&swapgen,buffsize,"swapgen","History::History");//array of number of generated swaps
	SetArraysize(&swapacc,buffsize,"swapacc","History::History");//array of number of accepted swaps
	SetArraysize(&E0shiftgen, buffsize, "E0shiftgen", "History::History");//array of number of generated E0 shifts
	SetArraysize(&E0shiftacc, buffsize, "E0shiftacc", "History::History");//array of number of accepted E0 shifts
	SetArraysize(&fprimeshiftgen, buffsize, "fprimeshiftgen", "History::History");//array of number of generated f' shifts
	SetArraysize(&fprimeshiftacc, buffsize, "fprimeshiftacc", "History::History");//array of number of accepted f' shifts
	SetArraysize(&mucorrgen, buffsize, "mucorrgen", "History::History");//array of number of generated mu corrections
	SetArraysize(&mucorracc, buffsize, "mucorracc", "History::History");//array of number of accepted mu corrections
	SetArraysize(&goodnonlinreg, buffsize, "goodnonlinreg", "History::History");//array of number of good nonlinear regression steps
	SetArraysize(&time,buffsize,"time","History::History");//array of time values
	SetArraysize(&ir1,buffsize,"ir1","History::History");//array of instant ratio tried/generated
	SetArraysize(&ir2,buffsize,"ir2","History::History");//array of instant rations accepted/tried
	SetArraysize(&chicomp,buffsize*chisize,"chicomp","History::History");// array of chi squared values
	
};		

//--------Initialising the static members-------------
void History::SetHistoryParams(RunParams &rundat)
{
	int i=0;
	buffsize=rundat.histbuffsize;//buffer size
	stepratio=rundat.histstepratio;//stepratio
	
	if (RunParams::swap_fraction>0)
		swap=1;
	else
		swap=0;

	//number of chi_squared parameters to save the total + the values for the 
	//expts data sets + the values for the coordination constraints
	chisize = rundat.ngr + rundat.nsq + rundat.nfq + rundat.nfg + rundat.nek + rundat.ncosdistr + CoordNumbConst::tot_subconst + rundat.navcoord;
#ifdef _ADVANCED_GEOM_CONST
	chisize += rundat.ncommonneigh + rundat.nsecondneigh+rundat.nbvs;
#endif
	if (chisize>0)
		chisize++;
#ifdef _LOCAL_INV
	nlocint=RunParams::nlocint;
	chisize+=nlocint;//for the locale invariance contribution
#endif
	if (RunParams::potential>0)
		chisize++;//there will be one more entry for the separately handled total bonded and non-bonded/aenet potential contributions
		//total number of data points
	datapoint_tot=RunParams::ntotal_points;//total number of data points
	if (RunParams::potential>0 || RunParams::fnc==4)
	{
		nNB = 0;
		switch (RunParams::potential)
		{
			case 1:
			{
				nNB = (fabs(RunParams::NB_weight_mode) == 1 ? SimpleCfg::npartials : 1);//number of nonbonded interactions /type(either 1, or npartials)
				chisize += 2 * nNB;//vdW+Coulomb (NB is 1, if RunParams::NB_weight_mode=0 or 2, npartials if 1)
				chi_offset1 = chisize - 1;//last element before 1-4
				chi_offset2 = chi_offset1 + 1;//next element after 1-4, default
				if (Topology::npairs > 0)//there is 1-4 interactions
				{
					if (nNB == 1)
						chi_offset2 = chi_offset1;//there is only 1 value, no need for cheking for used partials
					else
						chi_offset2 = chisize + 2 * nNB;//next element after 1-4 
					chisize += 2 * nNB;//vdW14+Coulomb14 contribution
				}
				datapoint_tot2 = RunParams::ntotal_points2;//total number of data points
				break;
			}
			case 10:
			{
				nNB = RunParams::nused_potpartials;
				chisize += nNB;//tabulated (NB is nused_potpartials,  RunParams::NB_weight_mode 1)
				chi_offset1 = chisize - 1;//last element before 1-4
				chi_offset2 = chi_offset1 + 1;//next element after 1-4, default
				datapoint_tot2 = RunParams::ntotal_points2;//total number of data points, 
				break;
			}
			case 20:
			{
				nNB = 1;
				chisize += 1;
				chi_offset1 = chisize - 1;//last element before 1-4
				chi_offset2 = chi_offset1 + 1;//next element after 1-4, default
				datapoint_tot2=1;
				break;
			}
		}
		
		//here the order is first the non-bonded interactions, then the bonded interactions in the order they are given here
		i=0;
		i+=(Topology::nbond_types>0 ? (Topology::bond_weight_mode ? Topology::nbond_types : 1) : 0);
		i+=(Topology::nangle_types>0 ? (Topology::angle_weight_mode ? Topology::nangle_types: 1) : 0);
		i+=(Topology::nRBdihedral_types>0 ? (Topology::RBdih_weight_mode ? Topology::nRBdihedral_types :1) : 0);
		i+=(Topology::nperdihedral_types >0 ? (Topology::perdih_weight_mode ? Topology::nperdihedral_types : 1) : 0);
		i+=(Topology::nharmdihedral_types>0 ?( Topology::harmdih_weight_mode ? Topology::nharmdihedral_types :1) : 0);
		chisize+=i;
	}
#ifdef _LOCAL_INV
	int j;
	SetArraysize(&loc_npoints,nlocint,"loc_npoints","SetHistoryParams");
	for (j=0;j<nlocint;j++)
		loc_npoints[j]=RunParams::loc_npoints[j];//number of points for local invariance
#endif



	
};
//------------------------function inifile-----------------------

void History::Inifile(ofstream &file, RunParams &rundat, ExptsData &edata)
{
	int i,isubconst;
	char *name,*name2,*name_ext=NULL;

	SetArraysize(&name,NAME_SIZE,"name","History::Inifile");
	SetArraysize(&name2,NAME_SIZE,"name2","History::Inifile");

	//running will not continue, till the saving could be performed
	CheckFileState(file,"History::EndNotes",hstfilename);
	
	file<<"This file contains History data produced by History::save "<<endl;
	if (rundat.custmove==0)
	{
		file<<RunParams::nthreads<<"-threaded RMC++ with "<<rundat.nmoved<<" moved atom(s)"<<endl;
#ifdef _USE_INT32
	file<<"\n_USE_INT32 option is ON: "<<sizeof(longint)<<" byte integers will be used for longint variables!"<<endl;
#else
	file<<"\n_USE_INT32 option is OFF: "<<sizeof(longint)<<" byte integers will be used for longint variables!"<<endl;
#endif
	
#ifdef _TEST_MODE 
	if (RunParams::runmode==1)
 		file<<"\nRunning in TEST MODE to ngenerated="<<-(int)RunParams::runlimit<<endl;
	else
		file << "\nRunning in TEST MODE to naccepted=" << (int)RunParams::runlimit << endl;
#else
	switch (RunParams::runmode)
	{
	case 1:
		file << "\nRunning to ngenerated=" << -(int)RunParams::runlimit << endl;
		break;
	case 2:
		file << "\nRunning to naccepted=" << (int)RunParams::runlimit << endl;
		break;
	default:
		file << " Running to " << RunParams::runlimit << " minutes" << endl;
	}
#endif
#ifdef _ADVANCED_GEOM_CONST
	file << "_ADVANCED_GEOM_CONST is ON, common neighbour, second neighbour and bond valence sum constraits are available!" << endl;
#endif
#ifdef _NO_PERIODIC
	file<<"\n_NO_PERIODIC option is ON, the sample will be in the middle of the simulation box!"<<endl;
#endif
#ifdef _VIBR_AMP
	int j,ipartial;
	file<<"\n_VIBR_AMP option is ON. The thermal atomic vibrations will be calculated by the "<<endl;
	file<<"convolution of the histogram with a Gaussian distribution with sigma parameter"<<endl;
	file<<"for the different partials:"<<endl;
	ipartial=0;
	for (i=0; i<RunParams::ntypes; i++)//for each first atom type
	{
		for(j=i;j<RunParams::ntypes;j++)//for each second atom type
		{
			file<<i<<" - "<<j<<" partial: "<<HistoSet::T_corr_sigma[ipartial]<<" Angstrom "<<endl;
			ipartial++;
		}
	}
#endif
#ifdef _LOCAL_INV 
	file<<"\n_LOCAL_INV option is used to calculate local invariance!"<<endl;
#endif
#ifdef _ATLAS 
	file<<"\nATLAS is used for matrix-vector multiplication!"<<endl;
#endif 
#ifdef _MUL_SCAT_VECTOR 
	file << "\n_MUL_SCAT_VECTOR option is ON: For all the S(Q), F(Q) and F(g) data sets the normally calculated S(Q) and F(Q) will be multiplied by Q and F(g) multiplied by g, and fitted!" << endl;
#endif
#ifdef _R_SWITCH_GR_MIN_1 
	file << "\n_R_SWITCH_GR_MIN_1 option is ON: For all the g(r) data set instead of the normally calculated g(r), ";
	if (r_switch_power == 1)
		file << "r*";
	else
		file << "(r^" << r_switch_power << ")*";
	file << "[g(r) - 1] is fitted!" << endl;
#endif
#ifdef _NEI 
	file<<"\n_NEI option is ON: The neighbourlist";
#ifdef _NEIE
	file<<" the sqaured distance and vector components";
#endif
	file<<" will be saved into the *.nei file !"<<endl;
#endif
#ifdef _AV_MOVE
	file<<"\n_AV_MOVE option is ON, the average moved distance of the atoms will be written at each screen displayand saved into the *.avm file!"<<endl;
#endif

#ifdef _CODE_WARRIOR_MAC
	file<<"\n_CODE_WARRIOR_MAC option is ON!"<<endl;
#endif

		if (rundat.fnc>0)
			file<<"FNC was used with option "<<rundat.fnc<<endl;
		else
		{
			if (RunParams::ntypes>1)
			{
				file<<"Fraction of swaps: "<<RunParams::swap_fraction<<endl;
				file<<"Number of swap pair types: "<<RunParams::nswap_pairs<<endl;
				file<<"Type of atoms in allowed swap pairs: "<<endl;
				for (i=0;i<RunParams::nswap_pairs;i++)
					file<<"\t"<<i+1<<". swap pair type "<<RunParams::swap_type1[i]+1<<" - "<<RunParams::swap_type2[i]+1<<endl;
			}
		}
	}
	else
		file<<RunParams::nthreads<<"-threaded molecular RMC++ with "<<rundat.nmoved<<" moved atoms"<<endl;
	file<<"Old format *.out file will "<<(RunParams::old_out ? "" : "not ")<<"be created."<<endl;
	file<<(RunParams::sum_ppcf ? "Summed " : "Oridinary ")<<"ppcf's will be saved into .ppcf file."<<endl;
#ifdef _ATLAS 
	file<<"ATLAS is used for matrix-vector multiplication!"<<endl;
#endif
	if (RunParams::chem_symbols != nullptr)
	{
		file << "Composition of the system: " << endl;
		for (i = 0; i < RunParams::ntypes; i++)
			file << RunParams::chem_symbols[i] << " fraction: " << (double)(SimpleCfg::pnatoms[i]) / SimpleCfg::ntotal << endl;
	}
	switch (configtype)
	{
		case 1:
			file<<"The .cfg and .bcf files were consistent!"<<endl;
			file<<"The cooordinates from the binary type .bcf file were used!"<<endl;
			break;
		case 2:
			file<<"The .cfg and .bcf files were not consistent!"<<endl;
			file<<"The cooordinates from the text type .cfg file were used!"<<endl;
			break;
		case 3:
			file<<"There was no text type .cfg file!"<<endl;
			file<<"The cooordinates from the binary type .bcf file were used!"<<endl;
			break;
		default: 
			file<<"There was no binary .bcf file!"<<endl;
			file<<"The cooordinates from the text type .cfg file were used!"<<endl;
			break;
	}
	switch (RunParams::potential)
	{
		case 0:
		{
			file << "\nPotential will not be used." << endl;
			break;
		}
		case 1:
		{
			file << "\n\t 6-" << RunParams::LJ_rep_N << " Lennard-Jones potential will be used" << endl;
			file << "\t\t Force field name: " << Topology::force_field_name << endl;
			file << "\t\t LJ cutoff:  " << RunParams::vdW_cutoff << " A" << endl;
			file << "\t\t Coulomb cutoff:  " << RunParams::Coulomb_cutoff << " A" << endl;;
			file << "\t\t Scaling factor for 1-4 LJ interactions: " << RunParams::vdW14_fudge << endl;
			file << "\t\t Scaling factor for 1-4 Coulomb interactions: " << RunParams::Coulomb14_fudge << endl;
			file << "\t\t For each different GROMACS partial:" << endl;
			file << "\t\t " << J13 << "Atom type1" << J13 << "Atom type2" << J16 << "C6 (kj/mol*A^6)" << J4 << "C" << RunParams::LJ_rep_N << " (kj/mol*A^" << RunParams::LJ_rep_N << endl;
			file.setf(ios::fixed, ios::floatfield);
			file.setf(ios::right, ios::adjustfield);

			int j, k = 0;
			for (i = 0; i < RunParams::nGRtypes; i++)
			{
				for (j = i; j < RunParams::nGRtypes; j++)
				{
					if (RunParams::vdW_comb_rule == 0|| RunParams::vdW_comb_rule == 1)
						file << "\t\t " << J13 << Topology::atom_type[i] << J13 << Topology::atom_type[j] << J16 << RunParams::vdW_pot1[k] << J20 << RunParams::vdW_pot2[k] << endl;
					else
						file << "\t\t " << J13 << Topology::atom_type[i] << J13 << Topology::atom_type[j] << J16 << 4 * RunParams::vdW_pot2[k] * pow(RunParams::vdW_pot1[k], 6.0) << J20 << 4 * RunParams::vdW_pot2[k] * pow(RunParams::vdW_pot1[k], RunParams::LJ_rep_N) << endl;
					k++;
				}
			}
			file << "\t\t " << J16 << "   " << "Weights for the LJ " << J16 << " Coulomb interaction:" << endl;
			if (fabs(RunParams::NB_weight_mode) != 1)
				file << "\t\t " << J16 << "   " << J16 << RunParams::vdW_weight[0] << J16 << RunParams::Coulomb_weight[0] << endl;
			else
			{
				for (i = 0; i < SimpleCfg::npartials; i++)
					file << "\t\t " << J2 << i + 1 << J16 << ". RMC partial   " << J16 << RunParams::vdW_weight[i] << J16 << RunParams::Coulomb_weight[i] << endl;
			}

			file.unsetf(ios::fixed);
			file.unsetf(ios::right);
			break;
		}
		case 10:
		{
			file << "\nTabulated potential will be used for " << RunParams::nused_potpartials << " partial" << endl;

			file << "\t\t Teperature: " << RunParams::temperature << " K " << endl;
			file << "\t\t For partials:" << endl;
			for (i = 0; i < RunParams::nused_potpartials; i++)
			{
				file << FNC_POT::tabpot_index[i] + 1 << ". partial" << endl;
				file << "\tCutoff: " << FNC_POT::tab_cutoff[i] * SimpleCfg::boxedge << " weight parameter : " << RunParams::vdW_weight[i] << endl;
			}
			break;
		}
#ifdef _AENET
		case 20:
			file<<"\nANN potential is used"<<endl;
			file << "\tCutoff: " << Aenet::cutoff<<endl;
			file<<"\tANN potential is calculated at each "<<RunParams::aenet_step<<". step"<<endl;
			break;
#endif
	}
	if (RunParams::fnc==4)
		Topology::PrintParams(file);

	file<<"Number of used data points:"<<endl;
	for (i=0;i<rundat.ngr;i++)
	{

		mystrcpy(name, NAME_SIZE, "g(r) set #");
		IntToStr(&name_ext,i+1);
		mystrcat(name, NAME_SIZE, name_ext);
		file<<J35<<name<<"\t"<<edata.grused[i]<<endl;
	}
	for (i=0;i<rundat.nsq;i++)
	{
		mystrcpy(name, NAME_SIZE, "S(Q) set #");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		file<<J35<<name<<"\t"<<edata.sqused[i]<<endl;
	}
	for (i=0;i<rundat.nfq;i++)
	{
		mystrcpy(name, NAME_SIZE, "F(Q) set #");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		file<<J35<<name<<"\t"<<edata.fqused[i]<<endl;
	}
	for (i = 0; i < rundat.nfg; i++)
	{
		mystrcpy(name, NAME_SIZE, "F(g) set #");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		file << J35 << name << "\t" << edata.fgused[i] << endl;
	}
	for (i=0;i<rundat.nek;i++)
	{
		mystrcpy(name, NAME_SIZE, "E(k) set #");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		file<<J35<<name<<"\t"<<edata.ekused[i]<<endl;
	}
	if (rundat.ncosdistr>0)
	{
		file<<J20<<"Cosine distribution constraint \t"<<CosDistrConst::ncos_bin<<" for each constraint"<<endl;
	}
#ifdef _LOCAL_INV
	if (RunParams::loc_chi2_mode==0)
		file<<J35<<"Bin-based local invariance points \t"<<loc_npoints[0]<<", number of atoms: "<<SimpleCfg::ntotal<<", number of bins: "<<HistoSet::sum_bins<<endl;
	else
	{
		for (i=0;i<nlocint;i++)
			file<<J35<<"Distance-based local invariance points for interval "<<i+1<<":\t"<<loc_npoints[i]<<", number of central atoms: "<<SimpleCfg::ntotal<<", total number of neighbours: "<<ChiSquared::sum_loc_natoms[i]<<endl;
	}
#endif
	file<<"Total number of used data points: (Total number of g(r), S(Q), F(Q), F(g), E(k), Cosine distribution constraint, Coord. subconstraints, Average coord. const.";
#ifdef _ADVANCED_GEOM_CONST
	file << ", Common neigh. const, Second neigh const. Bond valence sum const.";
#endif
#ifdef _LOCAL_INV
	file<<", Local Invariance";
#endif
	if (RunParams::potential>0)
	{
		file<<endl;
		file<<datapoint_tot<<endl;
		switch (RunParams::potential)
		{
			case 1:
			{
				file << "Total number of non-bonded and bonded interaction related data points (vDW+Coulomb potential";
				break;
			}
			case 10:
			{
				file << "Total number of tabulated and bonded potential interaction related data points (";
				break;
			}
			case 20:
			{
				file<<"Total number of ANN potential related data points";
				break;
			}
		}
		if (fabs(RunParams::NB_weight_mode)==1 && RunParams::potential==1)
			file<<"for each partial";
		if (Topology::npair_types > 0)
			file << ", vDW+Coulomb 1-4 potential";
		if (fabs(RunParams::NB_weight_mode == 1))
			file << "for each used partial";
		if (RunParams::potential!=20)
			file << ")";
		if (RunParams::fnc == 4)
			file << ", bond types, angle types, dihedral types";
		
		file<<"): \t";
		file<<"\n"<<datapoint_tot2<<endl;

				
	}
	else
	{
		if (RunParams::fnc==4)
			file<<", bond types, angle types, dihedral types";
		file<<"): \t";
		file<<datapoint_tot<<endl;
	}
	file << "Leading series for the simulation: " << RunParams::lead_series_name << endl;
	if (rundat.ngr>0)
	{
		file<<"Sigma values for g(r) data sets:"<<endl;
		for (i=0;i<rundat.ngr;i++)
			file<<"\t"<<i+1<<". data set\t"<<ExptsData::grsigma[i]<<endl;
	}
	if (rundat.nsq>0)
	{
		file<<"Sigma values for S(Q) data sets:"<<endl;
		for (i=0;i<rundat.nsq;i++)
			file<<"\t"<<i+1<<". data set\t"<<ExptsData::sqsigma[i]<<endl;
	}
	if (rundat.nfq>0)
	{
		file<<"Sigma values for F(Q) data sets:"<<endl;
		for (i=0;i<rundat.nfq;i++)
			file<<"\t"<<i+1<<". data set\t"<<ExptsData::fqsigma[i]<<endl;
	}
	if (rundat.nfg > 0)
	{
		file << "Sigma values for F(g) data sets:" << endl;
		for (i = 0; i < rundat.nfg; i++)
			file << "\t" << i + 1 << ". data set\t" << ExptsData::fgsigma[i] << endl;
	}
	if (rundat.nek>0)
	{
		file<<"Sigma values for E(k) data sets:"<<endl;
		for (i=0;i<rundat.nek;i++)
			file<<"\t"<<i+1<<". data set\t"<<ExptsData::eksigma[i]<<endl;
	}
	if (rundat.ncosdistr>0)
	{
		file<<"Sigma values for Cosine distribution constraints:"<<endl;
		for (i=0;i<rundat.ncosdistr;i++)
			file<<"\t"<<i+1<<". constraint\t"<<CosDistrConst::weights[i]<<endl;
	}
	if (rundat.nicoord>0)
	{
		file<<"Sigma values for the coordination constraints:"<<endl;
		for (i=0;i<rundat.nicoord;i++)
			for (isubconst=0;isubconst<CoordNumbConst::n_subconst[i];isubconst++)
				file<<"\t"<<isubconst+1<<". subconstraint of "<<i+1<<". constraint\t"<<CoordNumbConst::weights[CoordNumbConst::cum_n_subconst[i]+isubconst]<<endl;
	}
	if (rundat.navcoord>0)
	{
		file<<"Sigma values for the average coordination constraints:"<<endl;
		for (i=0;i<rundat.navcoord;i++)
			file<<"\t"<<i+1<<". constraint\t"<<AvCoordConst::weights[i]<<endl;
	}
#ifdef _ADVANCED_GEOM_CONST
	if (rundat.ncommonneigh > 0)
	{
		file << "Sigma values for the common neighbour constraints:" << endl;
		for (i = 0; i < rundat.ncommonneigh; i++)
			file << "\t" << i + 1 << ". constraint\t" << CommonNeighConst::weights[i] << endl;
	}
	if (rundat.nsecondneigh > 0)
	{
		file << "Sigma values for the second neighbour constraints:" << endl;
		for (i = 0; i < rundat.nsecondneigh; i++)
			file << "\t" << i + 1 << ". constraint\t" << SecondNeighConst::weights[i] << endl;
	}
	if (rundat.nbvs > 0)
	{
		file << "Sigma values for the bond valence sum constraints:" << endl;
		for (i = 0; i < rundat.nbvs; i++)
			file << "\t" << i + 1 << ". constraint\t" << BondValenceSumConst::weights[i] << endl;
	}
#endif
#ifdef _LOCAL_INV
		file<<"Sigma values for the local invariance:"<<endl;
		for (i=0;i<nlocint;i++)
			file<<"\t"<<i+1<<". interval:\t"<<RunParams::loc_inv_sigma[i]<<endl;
#endif
	if (RunParams::potential > 0)
	{
		file << "Leading potential series for the simulation: " << RunParams::lead_series_name2 << endl;
		if (fabs(RunParams::NB_weight_mode != 1))
		{
			if (RunParams::potential==1)
			{
				file << "Sigma value for the " << RunParams::vdW_name << " potential (1-4 interaction as well, if there is any):" << endl;
				file << "\t" << RunParams::vdW_weight[0] << endl;
			}
#ifdef _AENET
			else if (RunParams::potential==20)
			{
				file << "Sigma value for the ANN potential potential: "<< endl;
				file << "\t" << RunParams::aenet_weight << endl;
			}
#endif
		}
		else
		{//this is for tabulated as well
			file << "Sigma values for the " << RunParams::vdW_name << " potential:" << endl;
			for (i = 0; i < FNC_POT::pot_dim; i++)
				file << "\t" << i + 1 << ". partial\t" << RunParams::vdW_weight[i] << endl;
		}

		if (RunParams::potential == 1)
		{
			if (fabs(RunParams::NB_weight_mode != 1))
			{
				file << "Sigma value for the Coulomb potential (1-4 interaction as well, if there is any):" << endl;
				file << "\t" << RunParams::Coulomb_weight[0] << endl;
			}
			else
			{
				file << "Sigma values for the Coulomb potential partials:" << endl;
				for (i = 0; i < FNC_POT::pot_dim; i++)
					file << "\t" << i + 1 << ". partial\t" << RunParams::Coulomb_weight[i] << endl;
			}
		}
	}
	if (Topology::nbond_types>0)
	{
		if (Topology::bond_weight_mode==0)
		{
			file<<"Sigma value for the bond interaction:"<<endl;
			file<<"\t"<<Topology::bond_sigma[0]<<endl;
		}
		else
		{
			file<<"Sigmas values for the bond types"<<endl;
			for (i=0; i<Topology::nbond_types;i++)
				file<<"\t"<<i+1<<". type\t"<<Topology::bond_sigma[i]<<endl;
		}
		

	}
	if (Topology::nangle_types>0)
	{
		if (Topology::angle_weight_mode==0)
		{
			file<<"Sigma value for the angle interaction:"<<endl;
			file<<"\t"<<Topology::angle_sigma[0]<<endl;
		}
		else
		{
			file<<"Sigmas values for the angle types"<<endl;
			for (i=0; i<Topology::nangle_types;i++)
				file<<"\t"<<i+1<<". type\t"<<Topology::angle_sigma[i]<<endl;
		}
	}
	if (Topology::nperdihedral_types>0)
	{
		if (Topology::force_field_type==OPLSAA)
			mystrcpy(name, NAME_SIZE, "improper dihedral");
		else
			mystrcpy(name, NAME_SIZE, "proper dihedral");
		if (Topology::perdih_weight_mode==0)
		{
			file<<"Sigma value for the "<<name<<" interaction:"<<endl;
			file<<"\t"<<Topology::perdihedral_sigma[0]<<endl;
		}
		else
		{
			file<<"Sigmas values for the "<<name<<" types"<<endl;
			for (i=0; i<Topology::nperdihedral_types;i++)
				file<<"\t"<<i+1<<". type\t"<<Topology::perdihedral_sigma[i]<<endl;
		}
	}
	if (Topology::nharmdihedral_types>0)
	{
		if (Topology::harmdih_weight_mode==0)
		{
			file<<"Sigma value for the improper dihedral interaction:"<<endl;
			file<<"\t"<<Topology::harmdihedral_sigma[0]<<endl;
		}
		else
		{
			file<<"Sigmas values for the improper dihedral types"<<endl;
			for (i=0; i<Topology::nharmdihedral_types;i++)
				file<<"\t"<<i+1<<". type\t"<<Topology::harmdihedral_sigma[i]<<endl;
		}
	}
	if (Topology::nRBdihedral_types>0)
	{
		if (Topology::force_field_type==OPLSAA)
			mystrcpy(name, NAME_SIZE, "proper dihedral");
		else
			mystrcpy(name, NAME_SIZE, "RB dihedral");
		if (Topology::RBdih_weight_mode==0)
		{
			file<<"Sigma value for the "<<name<<" interaction:"<<endl;
			file<<"\t"<<Topology::RBdihedral_sigma[0]<<endl;
		}
		else
		{
			file<<"Sigmas values for the "<<name<<" types"<<endl;
			for (i=0; i<Topology::nRBdihedral_types;i++)
				file<<"\t"<<i+1<<". type\t"<<Topology::RBdihedral_sigma[i]<<endl;
		}
	}
	if (ExptsData::is_E0shift)
	{
		for (i = 0; i < ExptsData::nek; i++)
		{
			if (edata.ek_ngrid_in[i] > 0)
			{
				file << "EXAFS data set " << i + 1 << " maximum E0 shift: " << edata.ek_dE0[i] << " eV, for number of grid points in one direction: " << edata.ek_ngrid_in[i] << endl;
				file << "\tgrid points used from " << ExptsData::ek_gridfrom[i] << " to " << ExptsData::ek_gridto[i] << ", with k points used from " << ExptsData::ekmin[i] << " to " << ExptsData::ekmax[i] << endl;
			}
		}
	}
	
	for (i = 0; i < ExptsData::nfq; i++)
	{
		file << "For F(Q) data set " << i + 1;
		if (ExptsData::fqfitIQ[i])
		{
			file<<"\tI(Q) fitting was used"<<endl;
			if (ExptsData::fqusecompton[i])
				file << "\tCompton term will be used" << endl;
			if (ExptsData::fqIQbackgcorr[i])
				file << "\tusing scalable mu correction with mu: " << ExptsData::fqmu[i] << " and dmu_max " << ExptsData::fqdmumax[i]<<endl;

			if (ExptsData::fqAXS[i] > 0)
				file << "\tusing AXS for type " << ExptsData::fqAXS[i] << " with f' " << ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1] << endl;

		}
		else
			file << "\tF(Q) fitting was used" << endl;
		if ((ExptsData::fqAXS[i] > 0 && ExptsData::fqfprimecount[i] > 1) || (ExptsData::fqAXS[i] == 0 && ExptsData::fqfprimecount[i] > 0))//there are fixed f' using sets
		{
			file <<"\tfixed f' using atom types:" << endl;
			for (int itype = 0; itype < RunParams::ntypes; itype++)
			{
				if ((itype != ExptsData::fqAXS[i] - 1) && ExptsData::fqfprimeindex[i] & (int)pow(2, itype))//this type has fixed f' on it is not AXS
					file << "\t" << itype + 1 << ". type with f': " << ExptsData::fqfprime[i * RunParams::ntypes + itype] << endl;
			}
		}
		
	}
	file.unsetf(ios::adjustfield);
	file.setf(ios::left, ios::adjustfield);
	file<<"\n"<<J16<<"generated"<<"\t"<<J16<<"tried"<<"\t"<<J16<<"accepted";
	if (swap)
		file<<"\t"<<J16<<"generated_swap"<< "\t" << J16<<"accepted_swap";
	if (RunParams::potential>0)
		file<< "\t" << J16<<"potaccpeted";
	if (ExptsData::is_E0shift)
		file<< "\t" << J20<<"generated_E0shifts" << "\t" << J20 << "accepted_E0shifts" ;
	if (ExptsData::is_IQ == 2)
		file << "\t" << J20 << "good_nonlin._reg.";
	if (ExptsData::is_IQmucorr)
		file << "\t" << J16 << "generated_mucorr" << "\t" << J16 << "accepted_mucorr" ;
	if (ExptsData::is_AXS)
		file << "\t" << J20 << "generated_f'shifts" << "\t" << J20 << "accepted_f'shifts" ;
	file<< "\t" << J10<<"time"<< "\t" << J16<<"tried/gen"<< "\t" << J16<<"acc/tried";

	if (datapoint_tot>0)
		file<< "\t" << J27<<"chi2_tot/data_point;";
	if (datapoint_tot2>0)
		file<< "\t" << J30<<"chi2_pot_tot/data_point2;";
	for (i=0;i<rundat.ngr;i++)
	{
		if (ExptsData::gruseR[i])
			mystrcpy(name, NAME_SIZE, "Rw/g(r)_set#");
		else
			mystrcpy(name, NAME_SIZE, "chi2/g(r)_set#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27<<name;
	}
	for (i=0;i<rundat.nsq;i++)
	{
		if (ExptsData::squseR[i])
			mystrcpy(name, NAME_SIZE, "Rw/S(Q)_set#");
		else
			mystrcpy(name, NAME_SIZE, "chi2/S(Q)_set#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i=0;i<rundat.nfq;i++)
	{
		if (ExptsData::fquseR[i])
			mystrcpy(name, NAME_SIZE, "Rw/F(Q)_set#");
		else		
			mystrcpy(name, NAME_SIZE, "chi2/F(Q)_set#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i = 0; i < rundat.nfg; i++)
	{
		if (ExptsData::fguseR[i])
			mystrcpy(name, NAME_SIZE, "Rw/F(g)_set#");
		else
			mystrcpy(name, NAME_SIZE, "chi2/F(g)_set#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i=0;i<rundat.nek;i++)
	{
		if (ExptsData::ekuseR[i])
			mystrcpy(name, NAME_SIZE, "Rw/E(k)_set#");
		else
			mystrcpy(name, NAME_SIZE, "chi2/E(k)_set#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i=0;i<rundat.ncosdistr;i++)
	{
		mystrcpy(name, NAME_SIZE, "chi2/cosine_distr_C#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i=0;i<rundat.nicoord;i++)
	{
		for (isubconst=0;isubconst<CoordNumbConst::n_subconst[i];isubconst++)
		{
			mystrcpy(name, NAME_SIZE, "chi2/CoordConst#");
			IntToStr(&name_ext, i + 1);
			mystrcat(name, NAME_SIZE, name_ext);
			mystrcat(name, NAME_SIZE, "/");
			IntToStr(&name_ext, isubconst + 1);
			mystrcat(name, NAME_SIZE, name_ext);
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
	}
	for (i=0;i<rundat.navcoord;i++)
	{
		mystrcpy(name, NAME_SIZE, "chi2/AvCoordConst#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
#ifdef _ADVANCED_GEOM_CONST
	for (i = 0; i < rundat.ncommonneigh; i++)
	{
		mystrcpy(name, NAME_SIZE, "chi2/CommonNeighConst#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i = 0; i < rundat.nsecondneigh; i++)
	{
		mystrcpy(name, NAME_SIZE, "chi2/SecondNeighConst#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
	for (i = 0; i < rundat.nbvs; i++)
	{
		
		mystrcpy(name, NAME_SIZE, "chi2/BVSConst#");
		IntToStr(&name_ext, i + 1);
		mystrcat(name, NAME_SIZE, name_ext);
		mystrcat(name, NAME_SIZE, ";");
		file << "\t" << J27 << name;
	}
#endif
#ifdef _LOCAL_INV
	if (RunParams::loc_chi2_mode)
	{
		for (i=0;i<nlocint;i++)
		{
			mystrcpy(name, NAME_SIZE, "chi2/Local inv#");
			IntToStr(&name_ext, i + 1);
			mystrcat(name, NAME_SIZE, name_ext);
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
	}
	else
		file<<"\t"<<J27<<"chi2/Local inv.;";
#endif
	//Potential
	if (RunParams::potential>0)
	{
		//vdW, tabulated and aenet
		mystrcpy(name,NAME_SIZE, RunParams::vdW_name);
		mystrcat(name, NAME_SIZE, "_pot");
		if (fabs(RunParams::NB_weight_mode != 1))
		{
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
		else
		{
			
			for (i=0;i<FNC_POT::pot_dim;i++)
			{
				mystrcpy(name2, NAME_SIZE, name);
				IntToStr(&name_ext, i + 1);
				mystrcat(name2, NAME_SIZE, name_ext);
				mystrcat(name2, NAME_SIZE, ";");
				file << "\t" << J27 << name2;
			}
		}
		if (RunParams::potential==1)
		{

			if (fabs(RunParams::NB_weight_mode != 1))
				file<<"\t"<<J27<<"Coulomb_pot;";
			else
			{
				mystrcpy(name, NAME_SIZE, "Coulomb_pot");
				for (i=0;i<FNC_POT::pot_dim;i++)
				{
					mystrcpy(name2, NAME_SIZE, name);
					IntToStr(&name_ext, i + 1);
					mystrcat(name2, NAME_SIZE, name_ext);
					mystrcat(name2, NAME_SIZE, ";");
					file << "\t" << J27 << name2;
				
				}
			}

			if (Topology::npair_types > 0)
			{
				mystrcpy(name, NAME_SIZE, "1-4_");
				mystrcat(name, NAME_SIZE, RunParams::vdW_name);
				mystrcat(name, NAME_SIZE, "_pot");
				if (fabs(RunParams::NB_weight_mode != 1))
				{
					mystrcat(name, NAME_SIZE, ";");
					file << "\t" << J27 << name;
				}
				else
				{
					for (i = 0; i < FNC_POT::pot_dim; i++)
					{
						if (RunParams::used_14partials[i])//only display, if it was used
						{
							mystrcpy(name2, NAME_SIZE, name);
							IntToStr(&name_ext, i + 1);
							mystrcat(name2, NAME_SIZE, name_ext);
							mystrcat(name2, NAME_SIZE, ";");
							file << "\t" << J27 << name2 ;
						}
					}
				}

				mystrcpy(name, NAME_SIZE, "1-4_");
				mystrcat(name, NAME_SIZE, "Coulomb_pot");
				if (fabs(RunParams::NB_weight_mode != 1))
				{
					mystrcat(name, NAME_SIZE, ";");
					file << "\t" << J27 << name;
				}
				else
				{
					for (i = 0; i < FNC_POT::pot_dim; i++)
					{
						if (RunParams::used_14partials[i])//only display, if it was used
						{
							mystrcpy(name2, NAME_SIZE, name);
							IntToStr(&name_ext, i + 1);
							mystrcat(name2, NAME_SIZE, name_ext);
							mystrcat(name2, NAME_SIZE, ";");
							file << "\t" << J27 << name2 ;
						}
					}
				}
			}
		}
		file.unsetf(ios::scientific);
	}
	
	if (Topology::nbond_types>0)
	{
		mystrcpy(name, NAME_SIZE, "Bond");
		if (Topology::bond_weight_mode==0)
		{
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
		else
		{
			mystrcat(name, NAME_SIZE, "_type");
			for (i=0; i<Topology::nbond_types;i++)
			{
				mystrcpy(name2, NAME_SIZE, name);
				IntToStr(&name_ext, i + 1);
				mystrcat(name2, NAME_SIZE, name_ext);
				mystrcat(name2, NAME_SIZE, ";");
				file << "\t" << J27<<name2;
			}
		}
	}

	if (Topology::nangle_types>0)
	{
		mystrcpy(name, NAME_SIZE, "Angle");
		if (Topology::angle_weight_mode==0)
		{
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
		else
		{
			mystrcat(name, NAME_SIZE, "_type");
			for (i=0; i<Topology::nangle_types;i++)
			{
				mystrcpy(name2, NAME_SIZE, name);
				IntToStr(&name_ext, i + 1);
				mystrcat(name2, NAME_SIZE, name_ext);
				mystrcat(name2, NAME_SIZE, ";");
				file << "\t" << J27 << name2;
			}
		}
	}

	if (Topology::nperdihedral_types>0)
	{
		if (Topology::force_field_type==OPLSAA)
			mystrcpy(name, NAME_SIZE, "Improper_dihedral");
		else
			mystrcpy(name, NAME_SIZE, "Proper_dihedral");
		if (Topology::perdih_weight_mode==0)
		{
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
		else
		{
			mystrcat(name, NAME_SIZE, "_type");
			for (i=0; i<Topology::nperdihedral_types;i++)
			{
				mystrcpy(name2, NAME_SIZE, name);
				IntToStr(&name_ext, i + 1);
				mystrcat(name2, NAME_SIZE, name_ext);
				mystrcat(name2, NAME_SIZE, ";");
				file << "\t" << J27 << name2;
				
			}
		}
	}

	if (Topology::nharmdihedral_types>0)
	{
		mystrcpy(name, NAME_SIZE, "Improper_dihedral");
		if (Topology::harmdih_weight_mode==0)
		{
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
		else
		{
			mystrcat(name, NAME_SIZE, "_type");
			for (i=0; i<Topology::nharmdihedral_types;i++)
			{
				mystrcpy(name2, NAME_SIZE, name);
				IntToStr(&name_ext, i + 1);
				mystrcat(name2, NAME_SIZE, name_ext);
				mystrcat(name2, NAME_SIZE, ";");
				file << "\t" << J27 << name2;
			}
		}
	}

	if (Topology::nRBdihedral_types>0)
	{
		if (Topology::force_field_type==OPLSAA)
			mystrcpy(name, NAME_SIZE, "Poper_dihedral");
		else
			mystrcpy(name, NAME_SIZE, "RB_dihedral");
		if (Topology::RBdih_weight_mode==0)
		{
			mystrcat(name, NAME_SIZE, ";");
			file << "\t" << J27 << name;
		}
		else
		{
			mystrcat(name, NAME_SIZE, "_type");
			for (i=0; i<Topology::nRBdihedral_types;i++)
			{
				mystrcpy(name2, NAME_SIZE, name);
				IntToStr(&name_ext, i + 1);
				mystrcat(name2, NAME_SIZE, name_ext);
				mystrcat(name2, NAME_SIZE, ";");
				file << "\t" << J27 << name2;
			}
		}
	}

	delete [] name;
	delete [] name2;
	if (name_ext != NULL)
		delete [] name_ext;
	file<<endl;
	file.unsetf(ios::left);

	//might be additional saving to the file, not closed here
};

//------------------------function save---------------------------
void History::Save(ofstream &file)
{
	//running will not continue, till the saving could be performed
	CheckFileState(file,"History::Save",hstfilename);

	
	cout<<"**********************************************"<<endl;
	cout<<"*flushing the history record to the .hst file*"<<endl;
	cout<<"**********************************************"<<endl;
	
	int i,j;
	longint *pl1, *pl2, *pl3,*pl4,*pl5,*pl6,*pl7,*pl8,*pl9, *pl10, *pl11, *pl12, *pl13;
	double *pd1, *pd2, *pd3, *pd4;
	int k=0;
	
	//pointers positioning
	pl1=generated;
	pl2=tried;
	pl3=accepted;
	pl4=swapgen;
	pl5=swapacc;
	pl6 = potaccepted;
	pl7 = E0shiftgen;
	pl8 = E0shiftacc;
	pl9 = goodnonlinreg;
	pl10 = mucorrgen;
	pl11 = mucorracc;
	pl12 = fprimeshiftgen;
	pl13 = fprimeshiftacc;
	
	pd1=time;
	pd2=ir1;
	pd3=ir2;
	pd4=chicomp;
	
	for(i=0;i<fill;i++)//for all buffered steps
	{
		file.setf(ios::right, ios::adjustfield);
		file<< J16<<*pl1++;//generated moves
		file<< "\t" << J16<<*pl2++;//tried moves
		file<< "\t" << J16<<*pl3++;//accepted moves
		if (swap)
		{
			file<< "\t" << J16<<*pl4++;//generated swaps
			file<< "\t" << J16<<*pl5++;//accepted swaps
		}
		if (RunParams::potential>0)
			file<< "\t" << J16<<*pl6++;//potaccepted moves
		if (ExptsData::is_E0shift)
		{
			file << "\t" << J20 << *pl7++ ;//generated E0 shifts
			file << "\t" << J20 << *pl8++ ;//accepted E0 shifts
		}
		if (ExptsData::is_IQ==2)
			file << "\t" << J20 << *pl9++;//good nonlin reg
		if (ExptsData::is_IQmucorr)
		{
			file << "\t" << J16 << *pl10++ ;//generated mucorr
			file << "\t" << J16 << *pl11++;//accepted mucorr
		}
		if (ExptsData::is_AXS)
		{
			file << "\t" << J20 << *pl12++;//generated f' shifts
			file << "\t" << J20 << *pl13++;//accepted f' shifts
		}
		file.precision(0);
	
		file.setf(ios::fixed, ios::floatfield);
		file<< "\t" << J10<<*pd1++;//time
		file.precision(6);
		file<< "\t" << J16<<*pd2++;//instant ratio tried/generated
		file<< "\t" << J16<<*pd3++;//instant ratio accepted/tried
		//now the buffered chi squared values
		file.precision(15);
		file.setf(ios::scientific, ios::floatfield);
		file.unsetf(ios::right);
		file.setf(ios::left, ios::adjustfield);
		if (chisize>1)
		{
			for(j=0;j<chisize;j++)
			{
				if (j>chi_offset1 && j<chi_offset2)
				{
					if (RunParams::used_14partials[k])
						file<< "\t" << J27<<*pd4++;//this is a used partial, write the value
					else
						pd4++;//this partial is not used, skip it
					if (k==SimpleCfg::npartials-1)//end of vdW interaction, set back k to 0 for Coulomb
						k=0;
					else
						k++;

				}
				else
					if(j==1 && RunParams::potential > 0)
						file << "\t" << J30 << *pd4++;
					else
						file << "\t" << J27 << *pd4++;
			}
		}
		file<<endl;
	}
	file.precision(6);
	file.unsetf(ios::fixed);
	file.unsetf(ios::left);
	file.unsetf(ios::scientific);
	fill=0;//reset the number of buffered lines
};

//-----------------function buffline-----------------------------------------

//function used to buffer one line of the .hst file
//subtract is only used presently for ANN pot calculation, where if step>1, it is possible, that the last accepted move 
//of the simulation did not calculate the ANN potential. Therefore at the end if necessary the ANN potential is recalculated
//and if the line was already buffered, subtract=-1, and the buffer line will be rewritten 
void History::Buffline(double &chron, RunParams &rundat, ChiSquared &chi, ExptsData &edata)
{
	if (fill>=buffsize)//safety check -this should never happened-
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"History::buffline attempting to buffer a line beyond the size "<<endl;
		cout<<"allocated for the history object:"<<endl;
		cout <<"buffer size: "<<buffsize<<"; filling position: "<<fill<<endl;
		cout <<"exiting..."<<endl;
		CleanExit();
	}
	
	int j,k,start=1;
	double temp,*pchi;

	

	*(generated+fill)=rundat.n_generated;//generated moves
	*(tried+fill)=rundat.n_tried;//tried moves
	*(accepted+fill)=rundat.n_accepted;//accepted moves
	if (RunParams::potential>0)
		*(potaccepted+fill)=rundat.n_potaccepted;//potential-accepted moves
	if (swap)
	{
		*(swapgen+fill)=rundat.n_swapgen;//generated swaps
		*(swapacc+fill)=rundat.n_swapacc;//accepted swaps
	}
	if (ExptsData::is_E0shift)
	{
		*(E0shiftgen + fill) = rundat.n_E0shiftgen;//generated E0 shifts
		*(E0shiftacc + fill) = rundat.n_E0shiftacc;//accepted E0 shifts
	}
	if (ExptsData::is_IQ==2)
		*(goodnonlinreg + fill) = rundat.n_goodnonlinreg;//tried steps ended with good nonlin reg
	if (ExptsData::is_IQmucorr)
	{
		*(mucorrgen + fill) = rundat.n_mucorrgen;//generated mu correction
		*(mucorracc + fill) = rundat.n_mucorracc;//accepted mu correction
	}
	if (ExptsData::is_AXS)
	{
		*(fprimeshiftgen + fill) = rundat.n_fprimeshiftgen;//generated f' shifts
		*(fprimeshiftacc + fill) = rundat.n_fprimeshiftacc;//accepted f' shifts
	}
	*(time+fill)=chron;//time
	if(rundat.n_generated>last_gen)//safety to avoid division by 0 at start, the last_.. values refer to since the last buffering
		temp=double(rundat.n_tried-last_tried)/double(rundat.n_generated-last_gen);
	else 
		temp=0;
	*(ir1+fill)=temp;//instant ratio tried/generated
	if(rundat.n_tried>last_tried)//safety to avoid division by 0 at start
		temp=double(rundat.n_accepted-last_acc)/double(rundat.n_tried-last_tried);//instant ratio accepted/tried
	else 
		temp=0;
	*(ir2+fill)=temp;

	//now the chi squared components:
	j=fill*chisize;//number of elements already loaded in the chicommp array
				//used for positioning the pointer
	pchi=chicomp+j;
	if (datapoint_tot>0)
		*pchi++=chi.total/datapoint_tot;//total chi squared value;
	if (RunParams::potential>0)
	{
		*pchi++=chi.total2/RunParams::ntotal_points2;
		start++;
	}
	//individual contributions of experimental sets and CC's
	//this loop copies the chicomp array of the chi squared object in the
	//chicomp array of the history object -note that the chicomp array of
	//the history object has one more `column' for the chi.total
	// hence the range of the index, and the -1
	for(k=start;k<chisize;k++)
		*pchi++=*(chi.chicomp+k-start);

	fill+=1;//increment the number of buffered lines
};

void History::EndNotes(ofstream &file, ChiSquared &chi) const
{
	int i;

	//running will not continue, till the saving could be performed
	CheckFileState(file,"History::EndNotes",hstfilename);
		
	file.precision(6);
	file.setf(ios::scientific, ios::floatfield);
	file<<"\nThe Rw values for the data sets are:"<<endl;
	for (i=0;i<chi.ngr;i++)
		file<<J20<<"g(r) set #"<<i+1<<"\t"<<J13<<chi.Rw[i]<<endl;
	for (i=0;i<chi.nsq;i++)
		file<<J20<<"S(Q) set #"<<i+1<<"\t"<<J13<<chi.Rw[i+chi.ngr]<<endl;
	for (i=0;i<chi.nfq;i++)
		file<<J20<<"F(Q) set #"<<i+1<<"\t"<<J13<<chi.Rw[i+chi.ngr+chi.nsq]<<endl;
	for (i = 0; i < chi.nfg; i++)
		file << J20 << "F(g) set #" << i + 1 << "\t" << J13 << chi.Rw[i + chi.ngr + chi.nsq + chi.nfq] << endl;
	for (i = 0; i < chi.nek; i++)
		file << J20 << "E(k) set #" << i + 1 << "\t" << J13 << chi.Rw[i + chi.ngr + chi.nsq + chi.nfq + chi.nfg] << endl;
	file.unsetf(ios::scientific);

	if (RunParams::potential == 10)
	{
		file << "\nFinal E/kT for the partial in case of tabulated potential" << endl;
		for (i = 0; i < nNB; i++)
		{
			file << "\tTabulated potential partial  " << FNC_POT::tabpot_index[i]+1 << ".\t" << FNC_POT::vdW_pot[i] * 1000 / k_BOLTZMANN / RunParams::temperature << endl;//using kJ for U in the program
		}
	}
	for (i = 0; i < chi.nfq; i++)
	{
		if (chi.fit_index[chi.ngr+chi.nsq+i] > 31)//nonlin regression
		{
			file << "\tNon-linear regression last accepted chi2_end/chi2_start for F(Q) set "<<i+1<<":\t" << chi.last_nlr_chi2fr[i] << endl;
			file << "\tNon-linear regression largest decrease in chi2_end/chi2_start for F(Q) set " << i + 1 << " :" << ":\t" << chi.maxdec_nlr_chi2fr[i] << endl;

		}
	}
	if (chi.failed_nlr==RunParams::terminate_nonlin)//forced termination
	{
		file << "\nThe program is terminating, because the number of successive failed non-linear regressions reached " << RunParams::terminate_nonlin << "!" << endl;
		for (i = 0; i < ExptsData::nfq; i++)
			if (chi.fit_index[chi.ngr + chi.nsq + i]> 31)//non-lin reg
				cout << "Total number of failed non-linear regression attempts for x-ray data set  " << i + 1 << " is: " << ChiSquared::tot_failed_nlr[i] << endl;
		file << "\nIt has to be mentioned that after non-linear regression failed for a data set, the further chi2 calculation" << endl;
		file << "and parameter optimization for the subsequent data sets is not performed!" << endl;
		file << "You can increase the failure limit by including the TERMINATE-NONLIN = XXX into the *.dat file!" << endl;
	}
	else
	{
		for (i = 0; i < ExptsData::nfq; i++)
		{
			if (chi.fit_index[chi.ngr + chi.nsq + i] > 31)//non-lin reg
			{
				if (ChiSquared::tot_failed_nlr[i] > 0)
					file << "Total number of failed non-linear regression attempts for x-ray data set  " << i + 1 << " is: " << ChiSquared::tot_failed_nlr[i] << endl;
			}
		}
	
		file << "\nThe calculation method for the chi_square/data sets \n" << endl;
		for (i = 0; i < chi.ngr; i++)
		{
			file << J4 << i + 1 << " g(r) data series :   ";
			if (1 & calc_method)
				file << "square of the difference" << endl;
			else
				file << "normal" << endl;

			calc_method = calc_method >> 1;
		}
		for (i = 0; i < chi.nsq; i++)
		{
			file << J4 << i + 1 << " S(Q) data series :   ";
			if (1 & calc_method)
				file << "square of the difference" << endl;
			else
				file << "normal" << endl;

			calc_method = calc_method >> 1;
		}
		for (i = 0; i < chi.nfq; i++)
		{
			file << J4 << i + 1 << " F(Q) data series :   ";
			if (1 & calc_method)
				file << "square of the difference" << endl;
			else
				file << "normal" << endl;

			calc_method = calc_method >> 1;
		}
		for (i = 0; i < chi.nfg; i++)
		{
			file << J4 << i + 1 << " F(g) data series :   ";
			if (1 & calc_method)
				file << "square of the difference" << endl;
			else
				file << "normal" << endl;

			calc_method = calc_method >> 1;
		}
		for (i = 0; i < chi.nek; i++)
		{
			file << J4 << i + 1 << " E(k) data series :   ";
			if (1 & calc_method)
				file << "square of the difference" << endl;
			else
				file << "normal" << endl;

			calc_method = calc_method >> 1;
		}
		for (i = 0; i < chi.ncos; i++)
		{
			file << J4 << i + 1 << " Cosine distribution constraint :   ";
			if (1 & calc_method)
				file << "square of the difference" << endl;
			else
				file << "normal" << endl;

			calc_method = calc_method >> 1;
		}
	}
};

//extremely large potential values occured, give warning
void History::PotWarning(ostream &file, int pind, const char *name)
{
	file << "\nWARNING(" << ++warn << "): Larger than " << POT_WARNING << " "<<name<<" potential value occured last for the " << pind << ". partial!" << endl;
	file << "\tThis means that excluded atoms or virtual sites were unreasonably close to each other!" << endl;
	file << "\tOnly the index of the partial this problem last occured is displayed, might have " << endl;
	file << "\toccured for other partials as well... The closeness caused very high potential contribution," << endl;
	file << "\twhich means that due to the precision of the floating point representaion there can be rounding" << endl;
	file << "\terrors in the potential update! In double precision ~15 decimal digits can be represented!" << endl;
	file << "\tYou can check the severity of the problem by recalculating the potential from the last configuration " << endl;
	file << "\tof the simulation not using the load histogram (RELOAD) option!" << endl;
	return;
};
