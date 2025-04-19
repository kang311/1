//source RMC_POT.cpp
//Last changed 01.03.2023

//United to handle both atomic and molecular RMC++ 
//Multi-threaded version with long threads
//possibility of ATLAS is included
//possibility of flexibly handling molecules kept together with MD-like bonds, angles and dihedrals
//virtual sites included
#define RMC_VERSION "2023.1"

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"

int acceptable;//boolean indicator, whether the move is acceptable

int main(int argc, char* argv[])
{
#ifdef _TEST_MODE
	//for speed measurement
	double dur0, dur1, dur2, dur3, dur3a, dur3b, dur4, dur5, dur6, dur7;//GO
	std::chrono::duration<double, std::milli> elapsed;
	std::chrono::time_point<std::chrono::high_resolution_clock> time1;
	std::chrono::time_point<std::chrono::high_resolution_clock> time2;
#endif

	//---------------------------------------------------------------
	//declarations of parameters needed to read RMC files

	int i, k;//cycle variables
	int count = 0;//counting th range error for E(k) sets
	int ithread = 1;
	int current_index;//current index of the toclose atom in the move.tclindices array
	int *tcltemp;//temporary array to store the 'bad' atoms below cutoffs
	int silent_quit = 0;//waiting for response at the end of the program, default; if 1, than quit without asking

	char key;

	//related to the reading of thh data
	int ntooclose;//Number of atoms below cutoff
	int load_status;//whether the loading of the histogram or the coordination numbers were successful,
					//and further calculation is necessary
					//0: recalculate the whole histogram,
					//1: loading was successful, no further calculation is necessary
					//it is also used for loading status of the neighbourlist


	//===========declaration of variables used in main loop=============

	int history_buffered = 0;//Whether History was buffered for the given step in the loop

	int nsaves;//number of saves since last history buffering
	int mode=0;//	0: default, timelimit in minutes was given for the run, or running to runlimit accepted number of steps (both in normal or test mode), last is an accepted step anyway
			   //	1: runnig to -runlimit generated steps, either in _TEST_MODE or in normal mode (-runlimit is used for compatibility with fixed format input) 
	int cfg_collect;//the number of configurations yet to be collected 
	int cfg_frequency;//collect the configuration at each coll_frequency save, if 0, not collected
	int save_count = 1;//counter to count the number of saves between collecting the configurations
	
	longint &naccepted = RunParams::n_accepted;//number of accepted moves
	longint &npotaccepted = RunParams::n_potaccepted;//number of potaccepted moves
	longint &ntried = RunParams::n_tried;//number of tried moves
	longint &ngenerated = RunParams::n_generated;//number of generated moves		

	longint seed;//seed for the random number generator

	double old_chitotal;//the chi square total of the previous accepted move
	double chidiff;//difference of the new and old chi values and its exponential
	double duration;//the duration of the main program cycle
	double run_limit;//duration of the run, either in seconds or number of steps
	double time_save;//time period for saving in seconds
	double savelapse;//the elapsed time since last save


	//temp* Only used in _TEST_MODE or if timelimit is negative (the program runs till maximum ngenerated or accepted steps)
	double *tempchi=nullptr;//for the chisquared of the last acceptable move
	double *tempa;//for the renorm coeff 'a' of the last acceptable move
	double *tempb;//for the renorm coeff 'b' of the last acceptable move
	double *tempc;//for the renorm coeff 'c' of the last acceptable move
	double *tempd;//for the renorm coeff 'd' of the last acceptable move
	double *tempe;//for the renorm coeff 'e' of the last acceptable move
	double *tempalpha;//for the renorm coeff 'alpha' of the last acceptable move
	double *templast_nlr_chi2fr;//for nonln regression last chi2_end/chi2_start of the iteration for all F(Q) data sets
	double *tempmaxdec_nlr_chi2fr;//for nonln regression smallest chi2_end/chi2_start of the iteration for all F(Q) data sets
	
	//to preserve the values of the last accepted move in case of E0 shift for EXAFS 
	int *old_ek_gridind=NULL;

	//to preserve the value of the last acceppted muact in case of I(Q) mu correction
	double *old_muact=NULL;//for the actual mu of the last acceptable move
	//for updating some  arrays of ChiSquared, as in case of I(Q) mucorr exp data without calc data arrays will be only calculated 
	//at mucorr steps, so has to be restored if the move is rejected
	//as only s-containing arrays are not calculated at only mucorr steps, if the previous atomic move was rejected, than the last accepted
	//s is needed in case of fit_index=3
	double *old_f = NULL;//f
	double *old_ff = NULL;//ff
	//double *old_fx = NULL;//fx, these are not used presently
	//double *old_fx2 = NULL;//fx2
	double *old_s = NULL;//s
	double *old_ss = NULL;//ss
	double *old_fs = NULL;//fs

	//to preserve the value of the last acceppted f' shift in case of I(Q) AXS correction
	double *old_fprimeact = NULL;//for the actual f' of the last acceptable move
	double old_chitotal2 = 0;//the chi square total2 of the previous accepted move in case of non-bonded potential
	double nsteps_d;//double precision ngenerated or naccepted depending on the runmode
	double *present;//pointer to the duraton of the loop or the limiting steps

	bool calc_neighlist = false;//only calculate if there is cos distr, common or second neighbour constraint
	bool forced_termination = false;

//used in initialization as well
	int n_fnc = 1;//the number of fnc instances, 1 by default
	int use_old_pot=0;//use the old potential loaded, even if the histogram (and pot) had to be recalculated, only if continuation and bcf exist 
	int extended_fnc = 0;//1 if NB or tabulated potential or fnc==4 is used, 0 if neither
	int is_ANN=0;//whether ANN potential is used, it is set even if not compiled with AENET
	int is_pot_loaded = 0;//whether the ANN pot is loaded even if histogram or other load failed
#ifdef _AENET
	
	is_ANN=1;//if copiled with _AENET, cannot be used without it, it will terminate at parameter reading
#endif
	bool bcf_exist=1;//whether the bcf file exists
	T_interaction_type ti;
#ifdef _USE_LOCAL_INV
	int load_status2=0;
#endif
#ifdef _NO_PERIODIC
	int load_status3;
#endif

	ofstream outfile;

	time_t t_start, t_current, t_lastsave;//time values
	
	//files streams 
	ofstream hstfile, chifile, newcfgfile,log2file;
	/*OpenFile(log2file, "log2","main",0);
	log2file.precision(16);
	log2file.setf(ios::scientific, ios::floatfield);*/
	//to avoid warning in case of consecutive compilation
	//set the temporary array to mark the atoms below cutoff for zero, it will be created later, if needed
	tcltemp = NULL;
	ntooclose = 0;
	cout << "=====================================================================================================================" << endl;
	cout << "==================================================                     ==============================================" << endl;
	cout << "==================================================  RMC_POT++ program  ==============================================" << endl;
	cout << "==================================================   VERSION " << J7<< RMC_VERSION << "   ==============================================" << endl;
	cout << "==================================================                     ==============================================" << endl;
	cout << "=====================================================================================================================" << endl;
	
#ifdef _USE_INT32
	cout << "\n_USE_INT32 option is ON:   " << sizeof(longint) << " byte integers will be used for longint variables!" << endl;
#else
	cout << "\n_USE_INT32 option is OFF:  " << sizeof(longint) << " byte integers will be used for longint variables!" << endl;
#endif	

#ifdef _TEST_MODE 
	cout << "\nRunning in TEST MODE"<< endl;
	
#endif
#ifdef _ADVANCED_GEOM_CONST
	cout << "\n_ADVANCED_GEOM_CONST option is ON, common neighbour, second neighbour and bond valence sum constraints can be used!" << endl;
#endif
#ifdef _AENET
	cout<<"\nAENET option is ON, ANN potentials can be used!"<<endl;

#endif
#ifdef _NO_PERIODIC
	cout << "\n_NO_PERIODIC option is ON, the sample will be in the middle of the simulation box!" << endl;
#ifdef _WRITE_MESSAGE
	cout << "\nLOCAL INVARIACE calculation cannot be used with non-periodic boundary conditions, " << endl;
	cout << "_LOCAL_INV option will be switched off!!!" << endl;
#endif
#endif
#ifdef _VIBR_AMP
	cout << "\n_VIBR_AMP option is ON, the histogram will be convoluted with a Gaussian " << endl;
	cout << "\tfunction to account for the atomic vibrations." << endl;
#endif
#ifdef _LOCAL_INV 
	cout << "\n_LOCAL_INV option is used to calculate local invariance!" << endl;
#endif
#ifdef _ATLAS 
	cout << "\nATLAS is used for matrix-vector multiplication!" << endl;
#endif 
#ifdef _MUL_SCAT_VECTOR 
	cout << "\n_MUL_SCAT_VECTOR option is ON: For all the S(Q) F(Q) and F(g) data sets the normally calculated S(Q), F(Q) and F(g) will be multiplied by Q or g!" << endl;
#endif
#ifdef _R_SWITCH_GR_MIN_1 
	if (r_switch_power < 1)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "r_switch_power has to be equal or greater than 1!" << endl;
		cout << "Cannor run this way, exiting..." << endl;
		CleanExit();
	}
	cout << "\n_R_SWITCH_GR_MIN_1 option is ON: For all the g(r) data set instead of the" << endl;
	cout << "\tnormally calculated g(r), ";
	cout << "(r^N * [g(r) - 1] will be expected, calculated and outputted!" << endl;
	cout << "\tThe value of N can be modified in the free format *.dat by the R-SWITCH-POWER key word, presently N=2!" << endl;
#endif
#ifdef _NEI 
	cout << "\n_NEI option is ON: The neighbourlist";
#ifdef _NEIE
	cout << " the sqaured distance and vector components";
#endif
	cout << " will be saved into the *.nei file !" << endl;
#endif
#ifdef _AV_MOVE
	cout << "\n_AV_MOVE option is ON, the average moved distance of the atoms will be written at each screen display and saved into the *.avm file!" << endl;
#endif

#ifdef _CODE_WARRIOR_MAC
	cout << "\n_CODE_WARRIOR_MAC option is ON!" << endl;
#endif


	//========getting information phase=========================
	//the command line arguments can be found in data.cpp, in the command_arg table!
	//the program can be started the following ways:
	//- EXECUTABLE: the file name will be asked
	//- EXECUTABLE filename 
	//- EXECUTABLE filename cont
	//- EXECUTABLE filename y
	//- EXECUTABLE filename cont y
	//where cont means to continue the run, y stand for any single character to show, that the program can exit silently

	//general help, how the program can be started
	//- EXECUTABLE -help

	//help related to the free format data input
	//- EXECUTABLE -helptags
	//- EXECUTABLE -helpkeys
	//- EXECUTABLE -helpkeysman
	//- EXECUTABLE -helpkeysop
	//- EXECUTABLE -helpkeys tag
	//- EXECUTABLE -helpkey key

	//help related to the X-ray scattering factor calculation
	//- EXECUTABLE -helpXcoeffs
	//- EXECUTABLE -helpXcoeffsat
	//- EXECUTABLE -helpXcoeffsion
	//- EXECUTABLE -helpXcoeff symbol
	//- EXECUTABLE -helpXcoeffsS symbol
	//- EXECUTABLE -helpXcoeffsZ Z

	//help related to the X-ray Compton scattering calculation
	//- EXECUTABLE -helpCompton
	//- EXECUTABLE -helpComptonS symbol
	//-	EXECUTABLE -helpComptonZ Z

	//help related to neutron scattering length table
	//- EXECUTABLE -helpNscatlength
	//- EXECUTABLE -helpNscatlengthS symbol
	//-	EXECUTABLE -helpNscatlengthZ Z

	//help related to bond valence parameters table
	//- EXECUTABLE -helpBVparam
	//- EXECUTABLE -helpBVparamS symbol
	//-	EXECUTABLE -helpBVparamZ Z
	if (argc >= 2)
	{
		std::string s(argv[1]);
		CheckCommandlineParams(s, argc - 2);
	}
	switch (argc)
	{
	case 1://no command line arguments were specified
	{
		cout << "\nEnter the file name (WITHOUT EXTENSION!) for the simulation:" << endl;
		cin >> filename;
		break;
	}
	case 2:
	{
		if (strcmp(argv[1], "-helptags")==0)
		{
			GiveTagTable();//gives the whole tag table
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpkeys") == 0)
		{
			GiveKeyTable(0);//gives the whole key table
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpkeysop") == 0)
		{
			GiveKeyTable(2);//gives the optional key table
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpkeysman") == 0)
		{
			GiveKeyTable(1);//gives the mandatory key table
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpXcoeffs") == 0)
		{
			GiveWSymbols(0);//gives the symbols implemented for the calculation of the Xray coeffs
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpXcoeffsat") == 0)
		{
			GiveWSymbols(1);//gives the symbols implemented for the calculation of the Xray coeffs
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpXcoeffsion") == 0)
		{
			GiveWSymbols(2);//gives the symbols implemented for the calculation of the Xray coeffs
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpCompton") == 0)
		{
			GiveCSymbols(0);//gives the symbols implemented for the calculation of the Xray coeffs
			CleanExit();
		}
		else if (strcmp(argv[1], "-helpNscatlength") == 0)
		{
			GiveNSymbols(0);//gives the symbols implemented for the calculation of the neutron coeffs
			CleanExit();
		}
#ifdef _ADVANCED_GEOM_CONST
		else if (strcmp(argv[1], "-helpBVparam") == 0)
		{
			GiveBVSymbols(0);//gives the symbols implemented for the default R0 and b for bond valence sum constraint
			CleanExit();
		}
#endif
		else if (strcmp(argv[1], "-help") == 0)
		{
			Help();//general help
			CleanExit();
		}
		mystrcpy(filename, FILE_NAME_SIZE, argv[1]);
		silent_quit = 1;
		break;
	}
	case 3://the second can be either the "cont" option for the continuation of the run, 
			//or a single character to 	collapse the program window without asking
	{
		if (strcmp(argv[2], "cont") == 0)
			RunParams::continuation = 1;//try to continue the run
		else
		{
			 if (strcmp(argv[1], "-helpkey") == 0)
			 {
				 GiveKeyHelp(argv[2]);//gives info about the specified key
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpkeys") == 0)
			 {
				 GiveKeyTable(0,argv[2]);//gives the keys for the specified tag
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpXcoeff") == 0)
			 {
				 GiveWaasmaierPar(0, argv[2]);//gives the data for a symbol from the Xray coeffs table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpXcoeffsS") == 0)
			 {
				 GiveWSymbols(3, argv[2]);//gives the data for the entries beginning with symbol from the Xray coeffs table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpXcoeffsZ") == 0)
			 {
				 GiveWSymbols(4, argv[2]);//gives the data for the given atomic number from the Xray coeffs table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpXcoeff") == 0)
			 {
				 GiveComptonPar(0, argv[2]);//gives the data for a symbol from the Compton table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpComptonS") == 0)
			 {
				 GiveCSymbols(1, argv[2]);//gives the data for the entries beginning with symbol from the Compton table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpComptonZ") == 0)
			 {
				 GiveCSymbols(2, argv[2]);//gives the data for the given atomic number from the Compton table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpNscatlengthS") == 0)
			 {
				 GiveNSymbols(1,argv[2]);//gives the data for the entries beginning with symbol from the neutron scattering lenght table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpNscatlengthZ") == 0)
			 {
				 GiveNSymbols(2, argv[2]);//gives the data for the given atomic number from the neutron scattering lenght table
				 CleanExit();
			 }
#ifdef _ADVANCED_GEOM_CONST
			 else if (strcmp(argv[1], "-helpBVparamS") == 0)
			 {
				 GiveBVSymbols(1,argv[2]);//gives the data for the entries biginning with symbol from bond valence parameters table
				 CleanExit();
			 }
			 else if (strcmp(argv[1], "-helpBVparamZ") == 0)
			 {
				 GiveBVSymbols(2, argv[2]);//gives the data for the given atomic number from bond valence parameters table
				 CleanExit();
			 }
#endif
			 else
			 {
				 if (strlen(argv[2]) == 1)
					 silent_quit = 1;
				 else
					 Help();
			 }
		}
		mystrcpy(filename, FILE_NAME_SIZE, argv[1]);
		break;
	}
	case 4://the second command line should be "cont" the third a single character 
	{
		if (strcmp(argv[2], "cont") == 0)
			RunParams::continuation = 1;//try to continue the run
		else
			Help();
		if (strlen(argv[3]) == 1)
			silent_quit = 1;
		else
			Help();

		mystrcpy(filename, FILE_NAME_SIZE, argv[1]);
		break;
	}
	default:
	{
		Help();
	}
	}


	cout << "\n****************************************************************************" << endl;
	cout << "The program will run with file name: " << filename << endl;
	if (RunParams::continuation)
		cout << "The program will continue the previous run with the same file name." << endl;
	cout << "****************************************************************************\n" << endl;

	//make file names
	mystrcpy(datfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(datfilename, FILE_NAME_SIZE + 5, datext);
	mystrcpy(freerfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(freerfilename, FILE_NAME_SIZE + 5, freerext);
	mystrcpy(freersfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(freersfilename, FILE_NAME_SIZE + 5, freersext);
	mystrcpy(freefilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(freefilename, FILE_NAME_SIZE + 5, freeext);
	mystrcpy(cfgfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(cfgfilename, FILE_NAME_SIZE + 5, cfgext);
	mystrcpy(bincfgfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(bincfgfilename, FILE_NAME_SIZE + 5, bcfext);
	mystrcpy(fncfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(fncfilename, FILE_NAME_SIZE + 5, fncext);
	mystrcpy(hgmfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(hgmfilename, FILE_NAME_SIZE + 5, hgmext);
	mystrcpy(hstfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(hstfilename, FILE_NAME_SIZE + 5, hstext);
	mystrcpy(pgrfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(pgrfilename, FILE_NAME_SIZE + 5, pgrext);
	mystrcpy(psqfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(psqfilename, FILE_NAME_SIZE + 5, psqext);
	mystrcpy(pfqfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(pfqfilename, FILE_NAME_SIZE + 5, pfqext);
	mystrcpy(pfgfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(pfgfilename, FILE_NAME_SIZE + 5, pfgext);
	mystrcpy(pekfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(pekfilename, FILE_NAME_SIZE + 5, pekext);
	mystrcpy(movefilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(movefilename, FILE_NAME_SIZE + 5, movext);
	mystrcpy(cncfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(cncfilename, FILE_NAME_SIZE + 5, cncext);
	mystrcpy(cncdfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(cncdfilename, FILE_NAME_SIZE + 5, cncdext);
	mystrcpy(acnfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(acnfilename, FILE_NAME_SIZE + 5, acnext);
	mystrcpy(tcafilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(tcafilename, FILE_NAME_SIZE + 5, tcaext);
	mystrcpy(ppcffilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(ppcffilename, FILE_NAME_SIZE + 5, ppcfext);
	mystrcpy(gridfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(gridfilename, FILE_NAME_SIZE + 5, gridext);
#ifdef _NEI
	mystrcpy(neighfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(neighfilename, FILE_NAME_SIZE + 5, ".nei");
#endif
	mystrcpy(cosfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(cosfilename, FILE_NAME_SIZE + 5, cosext);
	mystrcpy(fitfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(fitfilename, FILE_NAME_SIZE + 5, fitext);
	mystrcpy(chifilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(chifilename, FILE_NAME_SIZE + 5, chiext);
	mystrcpy(statefilename, FILE_NAME_SIZE + 10, filename);
	mystrcat(statefilename, FILE_NAME_SIZE + 10, stateext);
	mystrcpy(logfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(logfilename, FILE_NAME_SIZE + 5, logext);
	mystrcpy(logfilename1, FILE_NAME_SIZE + 5, filename);
	mystrcat(logfilename1, FILE_NAME_SIZE + 5, logext);
	mystrcat(logfilename1, FILE_NAME_SIZE + 5, "1");
	mystrcpy(cfgcollectfilename, FILE_NAME_SIZE + 10, filename);
	mystrcat(cfgcollectfilename, FILE_NAME_SIZE + 10, "_coll");
#ifdef _ADVANCED_GEOM_CONST
	mystrcpy(concfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(concfilename, FILE_NAME_SIZE + 5, concext);
	mystrcpy(sncfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(sncfilename, FILE_NAME_SIZE + 5, sncext);
	mystrcpy(bvsfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(bvsfilename, FILE_NAME_SIZE + 5, bvsext);
	mystrcpy(bbvfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(bbvfilename, FILE_NAME_SIZE + 5, bbvext);
#endif
#ifdef _USE_LOCAL_INV
	mystrcpy(lhgmfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(lhgmfilename, FILE_NAME_SIZE + 5, lhgmext);
#endif
	mystrcpy(topfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(topfilename, FILE_NAME_SIZE + 5, topext);
	mystrcpy(exclfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(exclfilename, FILE_NAME_SIZE + 5, exclext);
	mystrcpy(potfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(potfilename, FILE_NAME_SIZE + 5, potext);
	mystrcpy(tabpotfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(tabpotfilename, FILE_NAME_SIZE + 5, tabpotext);
#ifdef _NO_PERIODIC
	mystrcpy(posbinfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(posbinfilename, FILE_NAME_SIZE + 5, binext);
#endif
#ifdef _AV_MOVE
	mystrcpy(avmfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(avmfilename, FILE_NAME_SIZE + 5, avmext);
#endif
#ifdef _AENET
	mystrcpy(energyfilename, FILE_NAME_SIZE + 5, filename);
	mystrcat(energyfilename, FILE_NAME_SIZE + 5, enerext);
#endif

	//determine the text file ending
	CheckEndLine();

	//Read the parameters related to the configurations from the text type .cfg file
	SimpleCfg::cfg_exist = SimpleCfg::GetParamsCfg();

	if (!SimpleCfg::cfg_exist)
		SimpleCfg::GetParamsBcf();



	RunParams::SetParams(SimpleCfg::ntotal, SimpleCfg::ntypes, SimpleCfg::boxedge, SimpleCfg::pnatoms);
	RunParams::GetParams();

#ifdef _TEST_MODE 
	cout << "NOTE(" << ++note << "): Running in TEST MODE up to ";
	if (RunParams::runmode==1)
		cout<< -RunParams::runlimit << " generated steps" << endl;
	else
		cout << RunParams::runlimit << " accepted steps" << endl;
#endif
	//Creating the logfie, if write_log is true
	if (write_log)
	{
		OpenFile(logfile, logfilename, "RMC_POT::main", 0);//open file, check, whether it was successfully opened
		/*OpenFile(logfile0, "log0", "RMC_POT::main", 0);//open file, check, whether it was successfully opened
		OpenFile(logfile1, "log1", "RMC_POT::main", 0);//open file, check, whether it was successfully opened
		OpenFile(logfile2, "log2", "RMC_POT::main", 0);//open file, check, whether it was successfully opened*/
	}

	if (RunParams::continuation && ExptsData::is_IQ)
		RunParams::continuation = -1;//this will be set back to 1 after the initial chi2 calculation, it is needed to differentiate 
			//the initial chi2 calc from the others, as there non-lin regression will not be performed, as the a,b,alpha is read from the state file

	
	if (RunParams::potential == 10)
		FNC_POT::SaveTabulatedPotential(tabpotfilename);

	if (RunParams::old_out)
	{
		mystrcpy(outfilename, FILE_NAME_SIZE + 5, filename);
		mystrcat(outfilename, FILE_NAME_SIZE + 5, outext);
	}

	if (RunParams::fnc == 4 && RunParams::potential == 0)
		RunParams::NB_weight_mode = 0;//readthe sigmas for the bonded interactions from topology
	//Cannot be done before the swap fraction is known
	//all the possible gen and acc types will be represented at all time
	SetArraysize(&SimpleCfg::nmoves, SimpleCfg::nmoves_dim, "SimpleCfg::nmoves", "main.cpp");

	*SimpleCfg::nmoves = &ngenerated;
	*(SimpleCfg::nmoves + 1) = &ntried;
	*(SimpleCfg::nmoves + 2) = &naccepted;
	*(SimpleCfg::nmoves + 3) = &npotaccepted;
	*(SimpleCfg::nmoves + 4) = &RunParams::n_swapgen;
	*(SimpleCfg::nmoves + 5) = &RunParams::n_swapacc;
	*(SimpleCfg::nmoves + 6) = &RunParams::n_E0shiftgen;
	*(SimpleCfg::nmoves + 7) = &RunParams::n_E0shiftacc;
	*(SimpleCfg::nmoves + 8) = &RunParams::n_mucorrgen;
	*(SimpleCfg::nmoves + 9) = &RunParams::n_mucorracc;
	*(SimpleCfg::nmoves + 10) = &RunParams::n_fprimeshiftgen;
	*(SimpleCfg::nmoves + 11) = &RunParams::n_fprimeshiftacc;
	*(SimpleCfg::nmoves + 12) = &RunParams::n_goodnonlinreg;
		
	SimpleCfg config;

	if (SimpleCfg::cfg_exist)
		config.Load();

	//GO to write out the coordinates in binary format
	//config.Save_binary();

	//Reading the parameters related to the configurations from the binary type
	//.bcf file, if it exists
	SimpleCfg config2;
	longint *moves=nullptr;//temporary array for the number of moves in the bcf, might be needed in case of mismatch
	bcf_exist=config2.LoadBinary(&moves);
	if (bcf_exist)
	{
		if (SimpleCfg::cfg_exist)//if cfg file exist
		{
			//The loading of the binary file was successful, the data has to be
			//compared with the data coming from the text .cfg file for consistency
			if (config2.Comp(config))//Comparing the coordinates)
			{
				//the coordinates are the same within  the given tolerance
				config.Copy(config2);//the "binary" coordinates will be used
				if (!RunParams::continuation)
					cout << "\nNOTE(" << ++note << "): The content of the .cfg and .bcf match, the binary file will be used!" << endl;
				else
				{
					if (moves != nullptr)//there was a mismatch in the number of moves, copy the binary values
						for (i = 0; i < SimpleCfg::nmoves_dim; i++)
							*(SimpleCfg::nmoves[i]) = moves[i];
				}
				History::configtype = 1;//the coordinates come from the binary .bcf file
			}
			else
			{
				//not consistent,
				if (!RunParams::continuation)
					History::configtype = 2;//the coordinates come from the text .cfg file
				else
				{
					History::configtype = 1;//the coordinates come from the binary .bcf file
					config.Copy(config2);//the "binary" coordinates will be used
					if (moves != nullptr)//there was a mismatch in the number of moves, copy the binary values
						for (i = 0; i < SimpleCfg::nmoves_dim; i++)
							*(SimpleCfg::nmoves[i]) = moves[i];
				}

			}
		}
		else
		{
			//there is no .cfg file
			config.Copy(config2);//the "binary" coordinates will be used
			if (!RunParams::continuation)
				cout << "\nNOTE(" << ++note << "): Text type .cfg file does nor exist, the values from the binary file will be used!" << endl;
			History::configtype = 3;//the coordinates come from the binary .bcf file
		}

	}
	config2.SimpleCfg::~SimpleCfg();//the config2 will not be needed any more
	config.PutToBox();

	extended_fnc = (((RunParams::potential > 0 && RunParams::potential<20)|| RunParams::fnc == 4) ? 1 : 0);//1 if calssical potential or fnc=4 is used
	//not used for AENET!
	if ((RunParams::potential == 1 || RunParams::fnc == 4))//not needed for only tabulated potential
	{
		//molecules are defined here, kept together by bond, angle and dihedral potential, read the topology
		//topology file is needed, even in case of atomic systems, if non-bonded, not tabulated potential is used to render the potential parameters to the RMC atoms
		//In case of tabulated ptential, no topology is needed!
		Topology::SetParams();
		Topology::GetParams();
		Topology::CheckRMCType();//check, whether an atom has the same RMC_type in all the instances of a molecule type
		Topology::CheckVirtualSites();//check, whether all the virtual sites of the [ atoms ] section are defined
		Topology::CheckExclusions();//check the exclusions given in the exclusion lines for inconsistencies
		if (Topology::n_gr_virtuals > 0)
		{
			SimpleCfg::SetVirtual();//set the params connected to the virtual sites
			config.ResizePosition();//resize the position array to make space for the virtual coordinates
			config.CalcVirtualCoords();//calculating the initial coordinates of the virtual sites
			
		}
		
		SimpleCfg::SetGRtype();//setting the GROMACS-RMC type arrays
		

	}
	if (RunParams::potential == 1 && RunParams::fnc == 4)
		Topology::CreatevdW14Params();//creating the 1-4 parameter list per GROMACS types 

	
	RunParams rundata;
		
	rundata.AssignBinsize();
	rundata.SetXmax();//not needed, if periodic boundary condition is not used
#ifdef _NO_PERIODIC
	//check and rescale the sample if necessary, as the reduced coordinates have to be between -0.5 and 0.5
	config.CheckSample();

#endif	
	
	Topology topology;//creating the topology object
	RunParams::SetBondedParams();

	if (RunParams::continuation)
	{
		cout << "\nLoading the " << statefilename << " file for exact continuation!" << endl;
		RunParams::LoadState();//loading the sigma values
	}

	HistoSet::SetHistParams(rundata);//set the static members
#ifdef _VIBR_AMP
	HistoSet::CalcConvFunction();//calculate the convolution function
#endif
	
	ExptsData edata(rundata);
	edata.SetE0ShiftParams();
	edata.NormEKCoeff();
	edata.Save();
	
	RunParams::SetPoints();//set the total number of data points
	//Get FNC parameters, if it is needed
	if (RunParams::custmove > 0 && RunParams::fnc != 1)
	{
		cout << "\nWARNING(" << ++warn << "): The custom move indicator in the .dat file is " << rundata.custmove << endl;
		cout << "\tbut the fnc indicator has not been switched on in the .dat file!" << endl;
		cout << "\tSwitching fnc indicator to 1, trying to run molecular RMC++" << endl;
		RunParams::fnc = 1;
	}
	if (RunParams::fnc > 0 && RunParams::fnc != 4)
		FNC_POT::GetFNCParams();

	

	if (extended_fnc)
		FNC_POT::SetPotParams();//setting the static arrays 

	//Create FNC object(s)
	FNC_POT *pfnc;
	if (RunParams::fnc == 4)
	{
		FNC_POT::natoms = SimpleCfg::ntotal;
		FNC_POT::nvirtuals = SimpleCfg::ntotvirtual;
		n_fnc = N_POT_TYPE;
	}

	pfnc = new FNC_POT[n_fnc];//to make it possible to use it with normal fnc
	
	if ((RunParams::fnc == 4) || (RunParams::potential > 0))//it is needed for tabulated as well, as the infrastructure for updating is done by pfnc[0]
	{
		for (i = 0; i < n_fnc; i++)//an instance will be created for each force type
		{
			ti = (T_interaction_type)i;
			pfnc[i].Init(ti);
		}
	}
	FNC_POT::config = &config;
	FNC_POT::rundat = &rundata;
	FNC_POT &fnc = pfnc[0];

	if (RunParams::fnc > 0 && RunParams::fnc != 4)
	{	//Load fnc constraints
		fnc.RawLoad();

#ifdef _TEST_MODE
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, fncext);
		fnc.Save(tempfilename);
#endif
		//Checking, whether the configurations satisfy the FNC constraints,
		//and whether there the FNC constraint is below the cutoff
		rundata.fnc_conflict = fnc.CheckFNC(rundata.pcutsq);
	}
	else
	{
		if (extended_fnc)
		{
			//Creating the index list for the different forces, setting other parameters
			for (i = 0; i < n_fnc; i++)
			{
				topology.CreateList(i, pfnc[i]);
				mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
				mystrcat(tempfilename, FILE_NAME_SIZE + 10, interaction_types[i]);
				mystrcat(tempfilename, FILE_NAME_SIZE + 10, fncext);
				pfnc[i].Save(tempfilename);
			}
		}
		//Create the exclusion index list and charge group centres for the non-bonded topology-based interactions, if there are molecules
		if (RunParams::fnc == 4 && RunParams::potential == 1)
		{
			config.SetChargeGroupCentre(pfnc[0]);//setting up the charge centres
			pfnc[0].GenerateExclusion(topology);
#ifdef _TEST_MODE
			pfnc[0].SaveExclusions(exclfilename);
#endif
		}
	}
	if (RunParams::auto_cutoff > 0)
	{
		for (i = 0; i < SimpleCfg::npartials; i++)
		{
			RunParams::pcutoff[i] = 4 * SimpleCfg::boxedge;//just to put it to a large value, for determining the minimum later 
			rundata.pcutsq[i] = 16;
		}
	}
	else if (RunParams::auto_cutoff <0)//read from state
	{
		for (i = 0; i < SimpleCfg::npartials; i++)
			rundata.pcutsq[i] = RunParams::pcutoff[i] * RunParams::pcutoff[i] / SimpleCfg::boxedge / SimpleCfg::boxedge;
	}

#ifdef _AENET
	Aenet::AenetInit();//initialize Aenet
	Aenet::SetCutoff();
#endif

	//Creating the gridding for the neighbourlist
	NeighbourList::SetParams();
	NeighbourList neighlist(config, pfnc, rundata);	

#ifdef _AENET
	Aenet aenetnew(config, rundata, neighlist);
	if (RunParams::potential!=20)
		aenetnew.Aenet::~Aenet();
	aenetnew.SetNeighList();
#endif
	load_status = neighlist.LoadGrid();//trying to load the neighbourlist
	if (load_status == 0)
	{
		//setting up the grid, if load failed
		neighlist.MakeGrid();
		neighlist.SaveGrid();

		//saving with a different name to keep the initial values
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, gridext);
		neighlist.SaveGrid(tempfilename);
	}

	config.Save(rundata.title);//saving the configuration in text format  
	CoordNumbConst iccnew(neighlist);
	AvCoordConst avccnew;
	CosDistrConst::SetParams();//set some static members
	CosDistrConst cosnew(neighlist, config);
	if (RunParams::ncosdistr)
		calc_neighlist = true;
#ifdef _ADVANCED_GEOM_CONST
	CommonNeighConst::SetParams();//set some static members
	CommonNeighConst concnew(neighlist, config);
	if (RunParams::ncommonneigh)
		calc_neighlist = true;
	SecondNeighConst::SetParams();//set some static members
	SecondNeighConst sncnew(neighlist, config);
	if (RunParams::nsecondneigh)
		calc_neighlist = true;
	BondValenceSumConst bvsnew;
#endif

	if (RunParams::ncosdistr)
	{
		cosnew.CalcTheoretical();//calculating the theoretical distribution
		cosnew.InitBincount();//initilaizing the histogram to zero 
	}
	if (calc_neighlist)
	{
		neighlist.InitNeighList();//calculating the neighbourlist, if there is cosine distribution of bond angles, conc or snc constraint
#ifdef _NEI
		neighlist.SaveNeilist();
#endif
	}
	if (RunParams::ncosdistr)
	{
		cosnew.CalcHist();//calculating the initial cosine distribution of bond angles
		cosnew.Save();

		//saving with a different name to keep the initial values
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, cosext);
		cosnew.Save(tempfilename);
	}
	CosDistrConst cosold(cosnew, neighlist, config);//creating the second instance by copying cosnew
#ifdef _ADVANCED_GEOM_CONST
	if (RunParams::ncommonneigh)
	{
		
		concnew.CalcHist();//calculating the initial common neighbour constraint
		concnew.SaveHist();

		//saving with a different name to keep the initial values
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, concext);
		concnew.SaveHist(tempfilename);
	}
	CommonNeighConst concold(concnew, neighlist, config);//creating the second instance by copying concnew
	if (RunParams::nsecondneigh)
	{

		sncnew.CalcNeigh();//calculating the initial second neighbour constraint
		sncnew.Save();

		//saving with a different name to keep the initial values
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, sncext);
		sncnew.Save(tempfilename);
	}
	SecondNeighConst sncold(sncnew, neighlist, config);//creating the second instance by copying sncnew
#endif


	//declaration of references to make the using of some parameters easier and quicker
	longint &nswapacc = rundata.n_swapacc;//number of accepted swaps
	longint &nswapgen = rundata.n_swapgen;//number of generated swaps
	longint &nE0shiftacc = rundata.n_E0shiftacc;//number of accepted shifts
	longint &nE0shiftgen = rundata.n_E0shiftgen;//number of generated shifts
	longint &nmucorracc = rundata.n_mucorracc;//number of accepted mu corrections
	longint &nmucorrgen = rundata.n_mucorrgen;//number of generated mu corrections
	longint &nfprimeshiftacc = rundata.n_fprimeshiftacc;//number of accepted f' shift steps
	longint &nfprimeshiftgen = rundata.n_fprimeshiftgen;//number of generated f' shift steps
	longint &ngoodnlr = rundata.n_goodnonlinreg;//number of tried moves ended with good non-lin regression steps
	longint &last_gen = rundata.last_gen;//number of generated moves at last display
	longint &last_tried = rundata.last_tried;//number of tried moves at last display
	longint &last_acc = rundata.last_accepted;//number of accepted moves at last display
	
	int &ngr = rundata.ngr;//the number of g(r) data sets
	int &nsq = rundata.nsq;//the number of S(Q) data sets
	int &nfq = rundata.nfq;//the number of F(Q) data sets
	int &nfg = rundata.nfg;//the number of F(g) data sets
	int &nek = rundata.nek;//the number of E(k) data sets
	int &nthreads = rundata.nthreads;//total number of threads

	//Saving the coordinates in text format, if there were no .cfg file
	if (History::configtype == 3)
	{
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, cfgext);
		config.Save("Starting configuration", tempfilename);//saving the configuration
	}
	
	//Create the Threads object for multithreading
	Threads::SetThreadsParams();
	Threads thread_obj;
	FNC_POT::thread_object = &thread_obj;
	//setting the handles for the object already created
	thread_obj.thread_p = &thread_obj;//this is needed for the static function ThreadLoop to see thread_obj
	thread_obj.edata_p = &edata;
	thread_obj.iccnew_p = &iccnew;
	thread_obj.avccnew_p = &avccnew;
#ifdef _ADVANCED_GEOM_CONST
	thread_obj.bvsnew_p = &bvsnew;
#endif
#ifdef _AENET
	thread_obj.aenetnew_p = &aenetnew;
	NeighbourList::thread = &thread_obj;
#endif

	//Create and initialise the histogram	
#ifdef _AENET
	HistoSet histnew(config, rundata, iccnew, avccnew, aenetnew, pfnc, thread_obj);
#ifdef _ADVANCED_GEOM_CONST
	histnew.bns = &bvsnew;
#endif
#else
	HistoSet histnew(config, rundata, iccnew, avccnew, pfnc, thread_obj);
#ifdef _ADVANCED_GEOM_CONST
	histnew.bvs = &bvsnew;
#endif
#endif
	thread_obj.histnew_p = &histnew;//setting the pointer

		//Create ChiSquared object
	ChiSquared::SetChiSquaredParams(rundata);
	RunParams::CheckLeadSerInd();//checking the lead series indices

	//Creating the complete (.free) and reduced (.freer) free format control files ,lead series index had to be reset for this to start with 1
	RunParams::lead_series_ind++;//to start with 1
	RunParams::lead_series_ind2++;//to start with 1
	RunParams::CreateFreeFormat(0);//*.free
	RunParams::CreateFreeFormat(1);//*.freer
	RunParams::lead_series_ind--;//to start with 0
	RunParams::lead_series_ind2--;//to start with 0
	if (ExptsData::is_Ncoeff_calc == 1)
		ExptsData::is_Ncoeff_calc = 2;//the coeffs were calculated, set it for *.freers saving
#ifdef _AENET
	ChiSquared chisquare(edata, cosnew, histnew, thread_obj,aenetnew);//needed for CalcData constructor
#else
	ChiSquared chisquare(edata, cosnew, histnew, thread_obj);//needed for CalcData constructor
#endif
#ifdef _ADVANCED_GEOM_CONST
	if (BondValenceSumConst::is_set_def == 1)
		BondValenceSumConst::is_set_def = 2;//for freers output
#endif


#ifdef _LOCAL_INV
	thread_obj.chi_p = &chisquare;//setting the pointer
#endif

	//splitting the work among the threads
	thread_obj.SplitToThreads();//the thread object's pointers to the histogram, iccmew, and avccnew has to be initialized already 

	Move::SetMoveParams();//Sets the static members, neeeded in case of the too close option
	//Loading the histogram 
	//Loading also the coordination numbers and the too close atom (if there is any)
	load_status = 0;//set the load_status to false
	
	if (RunParams::reload)
	{
		load_status = 1;//set the load_status to 1
		
		//loading the potential related parameters and the potential
		//even if the others cannot be loaded, for exact continuation the loaded potential can be used
		if (extended_fnc)
		{
			load_status = pfnc[0].LoadPotBinary();
			is_pot_loaded=load_status;
		}
#ifdef _AENET
		else if (is_ANN)
		{
			is_pot_loaded=aenetnew.LoadPotBinary();//this is not calculated during histogram calc
			if (RunParams::continuation && bcf_exist && is_pot_loaded)
				use_old_pot = 1;//Does not matter, whether the others were loaded, if bcf exist, the recalculated
				//hist and the others will be the same as if it were loaded (maybe the recalculated coeffs not)
		}
#endif
		if ((extended_fnc || is_ANN ) && is_pot_loaded == 0)//if the load was unsuccessful
		{
			if (RunParams::continuation)
				cout << "\nWARNING(" << ++warn << "): Exact continuation is not possible without loading successfully the binary potential file due to the rounding errors." << endl;
			if (extended_fnc)
				FNC_POT::InitPot();//setting potential back to zero
			//ANN does not need resetting
			
		}
		//if the moveout option is chosen the too close atoms are loaded
		if ( load_status > 0 && RunParams::moveout)
			load_status = Move::LoadTooClose();

		if (CoordNumbConst::nconstraints > 0 && load_status > 0)
		{
			//loading the coordination numbers from file if there is any  and the previous load was successful
			load_status = iccnew.Load(chisquare.sigma_percentage);//if load was unsuccessful, load_status is false
			if (load_status == 0)//if the load was unsuccessful
			{
				iccnew.InitCoordnumbs();//sets the element of the coordnumbs array back to zero
				Move::ResetTooClose();//deleting the "too close"-related arrays, if they were created
			}
		}
		//if there is average coordination constraint and the previous load was successful
		if (AvCoordConst::nconstraints > 0 && load_status > 0)
		{
			//loading the average coordination numbers from file if there is any
			load_status = avccnew.Load(chisquare.sigma_percentage);//if load was unsuccessful, load_status is false
			if (load_status == 0)//if the load was unsuccessful
			{
				avccnew.InitNeighbourcount();//sets the element of the neighbourcount array back to zero	
				if (CoordNumbConst::nconstraints > 0)
					iccnew.InitCoordnumbs();//sets the element of the coordnumbs array back to zero	
				Move::ResetTooClose();//deleting the "too close"-related arrays, if they were created
			}
		}//end of loading the average coordination constraint

#ifdef _ADVANCED_GEOM_CONST
		//if there is bond valence sum constraint and the previous load was successful
		if (BondValenceSumConst::nconstraints > 0 && load_status > 0)
		{
			//loading the bond valence sums from file if there is any
			load_status = bvsnew.LoadBinary(chisquare.sigma_percentage);//if load was unsuccessful, load_status is false
			if (load_status == 0)//if the load was unsuccessful
			{
				avccnew.InitNeighbourcount();//sets the element of the neighbourcount array back to zero	
				if (CoordNumbConst::nconstraints > 0)
					iccnew.InitCoordnumbs();//sets the element of the coordnumbs array back to zero	
				Move::ResetTooClose();//deleting the "too close"-related arrays, if they were created
			}
		}//end of loading the average coordination constraint
#endif

		//if the previous loads were successful
		if (load_status > 0)
		{
			load_status = histnew.Load();//if load was unsuccessful, load_status is false

#ifdef _USE_LOCAL_INV
			load_status2 = histnew.Load(1);//if load was unsuccessful, load_status is false
#ifdef _NO_PERIODIC
			load_status3 = histnew.LoadPosBin();//loading the central bin index for each atom
			if (load_status3 == 0)//if the load was unsuccessful
			{
				load_status = 0;
				load_status2 = 0;
			}
#endif
			if (load_status == 0 || load_status2 == 0)//if the load was unsuccessful
			{
				load_status = 0;//set it to zero, even if only the local histogram was missing
				histnew.InitLocHist();//sets the element of the histograms and total counts arrays back to zero	
#else

			if (load_status == 0)//if the load was unsuccessful
			{
#endif//end of _USE_LOCAL_INV
				if (CoordNumbConst::nconstraints > 0)
					iccnew.InitCoordnumbs();//sets the element of the coordnumbs array back to zero	
				if (AvCoordConst::nconstraints > 0)
					avccnew.InitNeighbourcount();//sets the element of the neighbourcount array back to zero
				histnew.InitHist();//sets the element of the histograms and total counts arrays back to zero	
				Move::ResetTooClose();//deleting the "too close"-related arrays, if they were created
				if (extended_fnc)//if the load was unsuccessful
				{
					if (RunParams::continuation && bcf_exist)//pot was loaded, otherwise woulld not reach this point
					{
						use_old_pot = 1;//although the potentials will be calculated, the pot read from the file will be used
									//if there is not a big difference, and the bcf file exist, so the other 
									//things which has to be recalculated are calculated from the bcf, so consistent
									//only the exotic things then the calculated coeffs might be a bit different
						for (i = 0; i < n_fnc; i++)
							pfnc[i].UpdatePot(0);//Copying the potential components to the 'old' arrays	to preserve them for comparing
					}
					FNC_POT::InitPot();//setting potential back to zero
				}
			}
		}//end of loading the histogram
	};//end of if the histogram and coordination constraint should be loaded
	if (::debug)
	{
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "_load.hgm");
		histnew.Save(tempfilename);
	}
	//create and initialise a temporary array to mark the atoms below cutoff if histogram will be calculated
	//Dim: number of configurations x number of atoms in a configuration
	if (load_status!=1)
		k=SimpleCfg::ntotal;
	else
		k=0;
	SetArraysize(&tcltemp,k,"tcltemp","main");
	
	for (i=0;i<k;i++)
		tcltemp[i]=0;

	//Creating the threads for the whole duration of the program
	//when paralell calculation is necessary, work will be assigned to these threads
	if (load_status==0)
		thread_obj.calc_hist_flag=1;//the load was unsuccessful, calculate the histogram
	
	for (ithread=0;ithread<nthreads-1;ithread++)
	{
		//assigning the parameters for the new threads
		thread_obj.thread_arg[ithread].thread_index=ithread;//the index of the thread
		try
		{
		  thread_obj.threads[ithread]=std::thread(&Threads::ThreadEntry,std::ref(thread_obj.thread_arg[ithread]));
		}
		catch(const std::system_error& e)
		{
			//there was an error at the thread creation, terminate the run
			cout << "\n*****ERROR*****" << endl;
			cout<<"The creation of the "<<ithread<<". thread failed with error code "<<e.code()<<endl;
			cout<<"Cannot run this way, exiting. Press any key!"<<endl;
			cin>>key;
			CleanExit();
		}
		/*if (status>0)
		{
			//there was an error at the thread creation, terminate the run
			cout<<"The creation of the "<<ithread<<". thread failed with error code "<<status<<endl;
			cout<<"Cannot run this way, exiting. Press any key!"<<endl;
			cin>>key;
			CleanExit();
		}*/
	}
	
	thread_obj.thread_arg[nthreads-1].thread_index=nthreads-1;//the index of the main thread
#ifdef _AENET
	if (is_ANN)
	{
		//THE POT WAS LOADED and there is bcf and continuation
		if (use_old_pot)//check the newly calculated pot against the one loaded earlier
		{
			if (!(aenetnew.CheckLoadedPot(thread_obj.thread_arg[RunParams::nthreads-1])))//see if it is inside tolerable error, probably this pot belongs to the configuration
			{
				is_pot_loaded=0;
			}
			Aenet::calc_ANN=0;//main thread's segment was calculated
			
		}
		if (is_pot_loaded)
		{
			thread_obj.calc_ANN=0;
			if (use_old_pot==0)
				Aenet::calc_ANN=0;
		}
	}
#endif
	if (load_status!=1)// if the load was unsuccessful or not attempted
	{
		histnew.tclarray=tcltemp;//to make it available to the histnew object
		histnew.HistCalc(ntooclose);//Calculating the histogram, multi-threading is managed by this routine
	}

	//Create DataMat object
	DataMat::SetDataMatParams(rundata);//set the static members
 	DataMat datmat(edata);
#ifdef _TEST_MODE
	datmat.Save();
#endif

#ifdef _NO_PERIODIC
	//---------------MULTITHREADING for calculating the periodic histogram-------------------
	thread_obj.NOP_count=0;//counter indicating how many thread finished with the calculation, reset it to 0
	if (nthreads>1)//signal to threads to start calculation 
	{
		std::unique_lock<std::mutex> guard(thread_obj.start2_mutex);//lock mutex
		thread_obj.start_flag=1;//time to proceed
		thread_obj.check_cd_start.notify_all();
	}

	//this has to preceed the thermal correction, if there is any
	if (datmat.calcmode==0)
		histnew.CalcPerHist(thread_obj.thread_arg[nthreads-1],datmat);	//calculate the peiodic histogram, if necessary

	//Threads have to wait, till all of them finished, as each of them will use all the bins for Fourie-transformation
	//in case of Q-space fitting, or for cin->r conversion in case of g(r) fitting. (If there is more than one data set, 
	//it would not always be possible to assign a bin unambigously to a thread during the smoothing!)
	if (nthreads>1)
	{
		std::unique_lock<std::mutex> guard(thread_obj.NOP_count_mutex);//lock mutex
		thread_obj.NOP_count++;
		if (thread_obj.NOP_count==nthreads)//time to proceed
		{
			thread_obj.start_flag=0;//reset it to 0
			thread_obj.NOP_count=0;//reset it to 0
			thread_obj.NOP_check_count.notify_all();
		}
		else thread_obj.NOP_check_count.wait(guard);
	}
#endif //_NO_PERIODIC

#ifdef _VIBR_AMP
#ifdef _USE_LOCAL_INV // it has to be calculated, even if the histogram was loaded, because hist_T is not loaded in this case
	histnew.ThermalCorr();//calculates the convoluted histogram to include atomic vibrations
#else
	if (load_status != 1)// if no local inv is used and hist_T was loaded, then no calculation necessary 
		histnew.ThermalCorr();//calculates the convoluted histogram to include atomic vibrations
#endif
#endif
	
	//Had to wait with saving for the periodic histogram calculation, if _NO_PERIODIC is present
	if (load_status!=1 || HistoSet::write_hist)// if the load was unsuccessful or not attempted, or it was successful, 
		//									but the data set histogram bin assignment changed, while the ndiffbin and binsizes are the same
	{
		//Save the histogram
		if (::debug)
			cout<<"\nSaving the histogram"<<endl;
		histnew.Save();
		 //Save the histogram with a different name to keep the starting histogram
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, hgmext);
		histnew.Save(tempfilename);

#ifdef _USE_LOCAL_INV
		//localc histogram is only saved after chi2 calculation, as in cae of distance based method average hist is only calculated during chi2 calc
		//average local histogram is not loaded! 
#ifdef _NO_PERIODIC
		if (::debug)
			cout<<"\nSaving the central indices for the bins, where the atoms can be found"<<endl;
		histnew.SavePosBin();
#endif
#endif
		//the too close atoms will be saved after the Move object is created

		//Save the result of the coordination number constraint calculation
		if (CoordNumbConst::nconstraints>0)
		{
			if (::debug)
				cout<<"\nSaving the result of the  coordination number constraint calculation"<<endl;
			iccnew.Save();
			if (CoordNumbConst::is_detail)
			{
				CoordNumbConst::SetParams();
				iccnew.SaveDetail();
			}
		}


		//Save the result of the average coordination number calculation
		if (AvCoordConst::nconstraints>0)
		{
			if (::debug)
				cout<<"\nSaving the result of the average coordination number constraint calculation"<<endl;
			avccnew.Save();
		}
#ifdef _ADVANCED_GEOM_CONST
		//Save the result of the bond valence sum calculation
		if (BondValenceSumConst::nconstraints > 0)
		{
			if (::debug)
				cout << "\nSaving the result of the bond valence sum constraint calculation" << endl;
			bvsnew.Save();
			bvsnew.SaveBinary();
		}
#endif
	}

	//Save the result of the coordination number constraint calculation with a different name to keep the starting values
	if (CoordNumbConst::nconstraints>0)
	{
		if (::debug)
			cout<<"\nSaving the starting values of the coordination number constraint calculation"<<endl;
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, cncext);
		iccnew.Save(tempfilename);
		if (CoordNumbConst::is_detail)
		{
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
			mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
			mystrcat(tempfilename, FILE_NAME_SIZE + 10, cncdext);
			iccnew.SaveDetail(tempfilename);
		}
	}
	//Save the result of the average coordination number calculation with a different name to keep the starting values
	if (AvCoordConst::nconstraints>0)
	{
		if (::debug)
			cout<<"\nSaving the starting values of the average coordination number constraint calculation"<<endl;
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, acnext);
		avccnew.Save(tempfilename);
	}
#ifdef _ADVANCED_GEOM_CONST
	//Save the result of the bond valence sum calculation with a different name to keep the starting values
	if (BondValenceSumConst::nconstraints > 0)
	{
		if (::debug)
			cout << "\nSaving the result of the bond valence sum constraint calculation" << endl;
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, bvsext);
		bvsnew.Save(tempfilename);
		bvsnew.SaveBinary();//no need to kepp the start
	}
#endif
	//Create and initialise the second HistoSet object based on the first one
#ifdef _AENET
	HistoSet histold(config,rundata,iccnew,avccnew,aenetnew,pfnc,thread_obj,histnew);//only the histogram is used from this object
#else
	HistoSet histold(config,rundata,iccnew,avccnew,pfnc,thread_obj,histnew);
#endif
	thread_obj.histold_p=&histold;//sets the pointer for the thread object
#ifdef _ADVANCED_GEOM_CONST
	histold.bvs = &bvsnew; //GO
	//Create and initialise the second BondValenceSumConst object based on the first one
	BondValenceSumConst bvsold(bvsnew);
	thread_obj.bvsold_p = &bvsold;
#endif
	//Create and initialise the second CoordNumbConst object based on the first one
	CoordNumbConst iccold(neighlist,iccnew);
	thread_obj.iccold_p=&iccold;

	//Create and initialise the second AvCoordConst object based on the first one
	AvCoordConst avccold(avccnew);
	thread_obj.avccold_p=&avccold;
	
	//Copy the values of the potential to the 'old' arrays to preserve them
	if (extended_fnc)
	{
		if (load_status!=1)// if the load was unsuccessful or not attempted
		{
			//calculating the bonded interactions
			if (RunParams::fnc==4)
			{
				for (i=0;i<N_POT_TYPE;i++)
					(pfnc[i].*InteractionFunct[i])();
			}
		
			if (use_old_pot)//check the newly calculated pot against the one loaded earlier, only in case of continuation
			{
				if (FNC_POT::CheckPot())//inside tolerable error, probably this pot belongs to the configuration
				{
					for (i=0;i<n_fnc;i++)
						pfnc[i].UpdatePot(1);//Copying the loaded potential components back to the 'new' arrays	to use the loaded ones
				}
				else
				{
					//warniing is given by CheckPot
					use_old_pot=0;
			
				}
			}
			if (use_old_pot==0)
				pfnc[0].SavePotBinary();//saving the potential related parameters and the potential
				
		}
		if (use_old_pot==0)//	No need to copy, it was already done if loaded is used
			for (i=0;i<n_fnc;i++)//not necessary in case of continuation with hgm recalulation
				pfnc[i].UpdatePot(0);//Copying the potential components to the 'old' arrays	to preserve them in case of the move is not accepted
		
	}

	//Create the Move object
	Move move(config);
	if (RunParams::custmove==0) //standard atomic RMC++
	{
		Move move1(ntooclose, tcltemp,config);
		move1.Copy(move);
	}
	else //Molecular RMC++
	{
		cout << "\nNOTE(" << ++note << "): Custom move: inspecting the .cus file and creating the Move"<<endl;
		Move move1(fnc,ntooclose, tcltemp,config);
		move1.Copy(move);
		cout<<"\tCustom move will be applied with "<<RunParams::nmoved<<" atoms moved simultaneously "<<endl;
		cout<<"\tmove custom parameters:"<<endl;
		cout.setf(ios::scientific, ios::floatfield);//formatting
		for (i=0;i<move.nparams;i++)
			cout<<i<<"->"<<*(move.custparams+i)<<endl;
		cout.unsetf(ios::scientific);//formatting
	}
	delete [] tcltemp;//this big array is not needed anymore
	
	//setting the pointers
	HistoSet::move=&move;
	//CoordNumbConst::move=&move;
	//AvCoordConst::move=&move;
	thread_obj.move_p=&move;//setting the Threads::move_p pointer
	
	if (load_status==0 && RunParams::moveout)
	{
		//Save the atoms below cutoff 
		if (::debug)
			cout<<"\nSaving the atoms below cutoff"<<endl;
		Move::SaveTooClose();
		
		//Save the atoms below cutoff with a different name to keep the starting values
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, tcaext);
		Move::SaveTooClose(tempfilename);
#ifdef _TEST_MODE
		if (histnew.ntcp>0)
		{
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
			mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start.tcp");
			histnew.SaveTooClosePairs(tempfilename);
		}	
#endif
	}

	if (RunParams::auto_cutoff > 0)
	{
		//calculate the cutoffs
		for (i = 0; i < SimpleCfg::npartials; i++)
		{
			rundata.pcutsq[i] -= TOLERANCE;
			RunParams::pcutoff[i] = sqrt(rundata.pcutsq[i]) * SimpleCfg::boxedge;
		}
		RunParams::auto_cutoff = 2;//indicating, that the new cutoffs were set
	}
	else if (RunParams::auto_cutoff <0) RunParams::auto_cutoff = 2;//indicating, that the new cutoffs were already loaded from state
	NeighbourList::SetCutoff();//set the cutoff dependent params

	//Calculate the PPCF
	PPCFSet::SetPPCFParams(rundata);
	PPCFSet ppcfnew(histnew,datmat,move);
	thread_obj.ppcfnew_p=&ppcfnew;//setting the Threads::ppcfnew_p pointer

	//Create CalcPart object
	CalcPart::SetCalcPartParams();
	CalcPart calcpnew(edata,datmat,ppcfnew,histnew,move);
	thread_obj.calcpnew_p=&calcpnew;//setting the Threads::calcpnew_p pointer



	//Create CalcData object
	CalcData::SetCalcDataParams();
	CalcData calcdat(edata,calcpnew,chisquare);
	thread_obj.calcdat_p=&calcdat;//setting the Threads::calcdat pointer

	//---------------MULTITHREADING for calculating the ANN energy if present, PPCF, partials and totals-------------------
	thread_obj.complete_count2=0;//counter indicating how many thread finished with the calculation, reset it to 0
	if (nthreads>1)//signal to threads to start calculation 
	{
		std::unique_lock<std::mutex> guard(thread_obj.start2_mutex);//lock mutex
		thread_obj.start_flag=1;//time to proceed
		thread_obj.check_cd_start.notify_all();
	}
#ifdef _AENET
#ifdef _TEST_MODE
	time1 = std::chrono::high_resolution_clock::now();
#endif
	//calculating the initial ANN potential
	if (Aenet::calc_ANN)
		aenetnew.CalcANN(thread_obj.thread_arg[nthreads-1]);
#ifdef _TEST_MODE
	time2 = std::chrono::high_resolution_clock::now();
	elapsed = time2 - time1;
	thread_obj.dur_ann1[nthreads-1]=elapsed.count();
#endif
#endif
	ppcfnew.CalcPPCF(thread_obj.thread_arg[nthreads-1]);

	//Threads have to wait, till all of them finished, as each of them will use all the bins for Fourie-transformation
	//in case of Q or g-space fitting, or for cin->r conversion in case of g(r) fitting. (If there is more than one data set, 
	//it would not always be possible to assign a bin unambigously to a thread during the smoothing!)
	if (nthreads>1)
	{
		std::unique_lock<std::mutex> guard(thread_obj.count_mutex);//lock mutex
		thread_obj.complete_count2++;
		if (thread_obj.complete_count2==nthreads)//time to proceed
		{
			thread_obj.start_flag=0;//reset it to 0
			thread_obj.complete_count2=0;//reset it to 0
			thread_obj.check_count.notify_all();
		}
		else thread_obj.check_count.wait(guard);
	}
#ifdef _AENET
	//copy the threads E_i to the array, only the changed ones are copied
	if (Aenet::calc_ANN)
	{
		aenetnew.E_tot=0;
		for (ithread=0;ithread<nthreads;ithread++)
		{
			{
				for (i=0;i<thread_obj.thread_arg[ithread].aenet_at_count;i++)
				{
					aenetnew.E_i[thread_obj.thread_arg[ithread].aenet_at_index[i]]=thread_obj.thread_arg[ithread].aenet_E_i[i];
					aenetnew.E_tot+=aenetnew.E_i[thread_obj.thread_arg[ithread].aenet_at_index[i]];
				}
			}		
		}
	}
	/*logfile.precision(12);
	logfile.setf(ios::scientific, ios::floatfield);
	logfile<<"E_tot "<<aenetnew.E_tot<<endl;*/
	if (Aenet::write_energy)
	{
		aenetnew.SaveEnergy();//save the atomic energies
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, enerext);//preserving the starting energy with a different name
		aenetnew.SaveEnergy(tempfilename);
	}
	Aenet aenetold(aenetnew,config,rundata,neighlist);
	thread_obj.aenetold_p=&aenetold;
	if (is_pot_loaded==0)
		aenetnew.SavePotBinary();//saving the potential related parameters and the potential

#endif
	
	
	//calculating the partials by the main thread
	if (ngr>0)//g(r) fitting
		calcpnew.CalcPartialgr(thread_obj.thread_arg[nthreads-1]);//converts the ppcfs from bins to r data point
	
	if (nsq>0)
		calcpnew.CalcPartialSq(thread_obj.thread_arg[nthreads-1]);//calculating the partial S(Q)-s by Fourier transformation
			
	if (nfq>0)
		calcpnew.CalcPartialFq(thread_obj.thread_arg[nthreads-1]);//calculating the partial F(Q)-s by Fourier transformation

	if (nfg > 0)
		calcpnew.CalcPartialFg(thread_obj.thread_arg[nthreads - 1]);//calculating the partial F(g)-s by Fourier transformation

	if (nek > 0)
		calcpnew.CalcPartialEk(thread_obj.thread_arg[nthreads - 1]);//calculating the partial E(k)-s by Fourier transformation
	
	calcdat.CalcTotal(thread_obj.thread_arg[nthreads-1]);//calculate the total from the partials by the main thread
	
	//waiting for all threads to finish
	if (nthreads>1)
	{
		std::unique_lock<std::mutex> guard(thread_obj.T_count_mutex);//lock mutex
		thread_obj.T_count++;
		if (thread_obj.T_count==nthreads)//time to proceed
		{
			thread_obj.start_flag=0;//reset it to 0
			thread_obj.T_count=0;//reset it to 0
			thread_obj.T_check_count.notify_all();
		}
		else thread_obj.T_check_count.wait(guard);
	};
	//Checking, whether range error (no non-zero data ponts) occued for E(k) sets
	k = 0;
	for (i = 0; i < nek; i++)
	{
		count = 0;
		for (ithread = 0; ithread < nthreads; ithread++)
		{

			if (thread_obj.thread_arg[ithread].ek_range_error[i] == 1)
				count++;
		}
		if (count==nthreads)//all the segments are all zero
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "All the E(k) data points are zero for " << i + 1 << ". EXAFS data set!" << endl;
			cout << "Probably the non-zero coefficient range used does not overlap with non-zero histogram range." << endl;
			cout << "The " << RunParams::assign_hist[ngr + nsq + nfq + nfg + i] + 1 << ". histogram width " << RunParams::rspacing[ngr + nsq + nfq + nfg + i] << " A bin size is used for this data set!" << endl;
			cout << "Check the " << hgmfilename << " file and the coefficients, and try again..." << endl;
			k = 1;
			
		}
		
	}
	if (k == 1)
	{
		cout << "No point running the simulation this way, exiting..." << endl;
		CleanExit();
	}
//---------------------------END OF MULTITHREADING----------------------------------------------------
	
	//Saving the ppcf
	ppcfnew.ncollect=0;//to set it to zero, before collecting the ppcf-s begin, if the option is chosen
	ppcfnew.Save();

#ifdef _TEST_MODE
	mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
	mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
	mystrcat(tempfilename, FILE_NAME_SIZE + 10, ppcfext);
	ppcfnew.Save(tempfilename);
#endif
	
	//Create the second PPCFSet object based on the first one
	PPCFSet ppcfold(ppcfnew,histnew,datmat,move);
	thread_obj.ppcfold_p=&ppcfold;//setting the Threads::ppcfold_p pointer

	//Create the second CalcPart object object based on the first one
	CalcPart calcpold(calcpnew,edata,datmat,ppcfnew,histnew,move);
	thread_obj.calcpold_p=&calcpold;//setting the Threads::calcpold_p pointer

	//saving the partials
	if (ngr>0)//g(r) fitting
		calcpnew.Savegr();
			
	if (nsq>0)
		calcpnew.SaveSQ();
			
	if (nfq>0)
		calcpnew.SaveFQ();
	
	if (nfg > 0)
		calcpnew.SaveFg();

	if (nek>0)
		calcpnew.SaveEK();
	
		
	//calculate the chi square value
#ifdef _ADVANCED_GEOM_CONST
	History::calc_method=chisquare.CalcChiSquared(calcdat,iccnew,avccnew,concnew,sncnew,bvsnew);
#else
	History::calc_method = chisquare.CalcChiSquared(calcdat, iccnew, avccnew);
#endif

#ifdef  _USE_LOCAL_INV
#ifndef _LOCAL_INV_NS//only save, if this was not given, may need a lot of time, it has to be saved here, as in case of distance-based calculation average histogram is only calculated during chi2 calculation
	if (load_status2 == 0) //local histogram loading was not successful
	{
		if (::debug)
			cout << "\nSaving the local histograms" << endl;
		histnew.SaveLoc(lhgmfilename);
		//Save the local histogram with a different name to keep the starting histogram
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
		mystrcat(tempfilename, FILE_NAME_SIZE + 10, lhgmext);
		histnew.SaveLoc(tempfilename);
	}
#endif
#endif
#ifdef _LOCAL_INV
	if (RunParams::loc_chi2_mode==1 && nthreads>1)
		thread_obj.ResetSegment();//try to improve load balancing for the segmentation of the neighbour atoms
#endif
	
	//saving the totals
	calcdat.Save();//might need the optimized I(Q), has to be after chi2 calc

	if (extended_fnc || is_ANN)//potential or fnc==4 is used
		chisquare.CalcPotChiSquared();//calculated the chi2 components fo rthe potentials, if non-bonded potential is present
	
	//create the *.freers file, to output calculated parameters, as scalable sigma, calculated cutoffs, calculated Ncoeffs, looked up 
	//lead series_ind and lead_series_ind2 is starting with 0 for the further calculations
	if (ChiSquared::calc_sigma && RunParams::continuation == 0) 
		chisquare.CalcSigma();//The initial sigma values has to be set for some of the serieses, lead_series_ind and lead_series_ind2 should start with 0
		
	//CreateFreeFormat and its routines were written for lead_series_ind and lead_series_ind2 starting with 1
	//here auto_cutoff should be already set to 2, if there is autocutoff det
#ifdef _ADVANCED_GEOM_CONST
	if ((ChiSquared::calc_sigma && RunParams::continuation == 0) || RunParams::auto_cutoff > 0 || ExptsData::is_Ncoeff_calc > 0 || BondValenceSumConst::is_set_def>0)
#else
	if ((ChiSquared::calc_sigma && RunParams::continuation == 0) || RunParams::auto_cutoff > 0 || ExptsData::is_Ncoeff_calc > 0)
#endif
	{
		RunParams::lead_series_ind++;//to start with 1
		RunParams::lead_series_ind2++;//to start with 1
		RunParams::CreateFreeFormat(2);//.freers reduced file with the scaled sigma parameters
		RunParams::lead_series_ind--;//to start with 0
		RunParams::lead_series_ind2--;//to start with 0

	}
	if (ChiSquared::calc_sigma && RunParams::continuation == 0)
		ChiSquared::calc_sigma = 0;//not to trigger error message, if the ledaing series' chi2 decreases to zero during the run
	
									   //has to be set after the last CreateFreeFormat
	cout << "\n******************************************" << endl;
	cout << "**** Number of warning messages: " << J4 << warn << " ****" << endl;
	cout << "**** Number of notes:            " << J4 << note << " ****" << endl;
	cout << "******************************************\n" << endl;
	rundata.PrintRunParams(edata);

	if (RunParams::potential==1 || RunParams::fnc==4)
		Topology::PrintParams(cout);

	CalcData calcdatold(calcdat,edata,calcpnew,chisquare);//needed for the saving of the totals, if the last move(s) before saving was/were
									//unacceptable, this will contain the last acceptable
#ifdef _TEST_MODE
	mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
	mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start");
	mystrcat(tempfilename, FILE_NAME_SIZE + 10, fitext);
	calcdat.SaveResult(tempfilename);//saving the totals
#endif
	if (RunParams::runmode == 1)
		mode = 1;//needed for copying back the values of the last accepted move
				//runmode=0 o 2, runlimit or to naccepted steps, last is accepetd, can run normally
	
	//needed for copying back the values of the last accepted move
	SetArraysize(&tempchi, chisquare.chisize, "tempchi", "main");
	SetArraysize(&tempa, ngr + nsq + nfq + nfg + nek, "tempa", "main");
	SetArraysize(&tempb, ngr + nsq + nfq + nfg + nek, "tempb", "main");
	SetArraysize(&tempc, ngr + nsq + nfq + nfg + nek, "tempc", "main");
	SetArraysize(&tempd, ngr + nsq + nfq + nfg + nek, "tempd", "main");
	SetArraysize(&tempe, ngr + nsq + nfq + nfg + nek, "tempe", "main");

	if (ExptsData::is_IQ > 0)
	{
		SetArraysize(&tempalpha, nfq, "tempalpha", "main");
		SetArraysize(&tempmaxdec_nlr_chi2fr, nfq, "tempmaxdec_nlr_chi2fr", "main");
		SetArraysize(&templast_nlr_chi2fr, nfq, "templast_nlr_chi2fr", "main");
		if (ExptsData::is_IQmucorr)
		{
			SetArraysize(&old_muact, nfq, "old_muact", "main");
			SetArraysize(&old_f, nfq, "old_f", "main");
			SetArraysize(&old_ff, nfq, "old_ff", "main");
			//SetArraysize(&old_fx, nfq, "old_fx", "main");
			//SetArraysize(&old_fx2, nfq, "old_fx2", "main");
			SetArraysize(&old_s, nfq, "old_s", "main");
			SetArraysize(&old_fs, nfq, "old_fs", "main");
			SetArraysize(&old_ss, nfq, "old_ss", "main");


		}
		if (ExptsData::is_AXS>0)
			SetArraysize(&old_fprimeact, nfq, "old_fprimeact", "main");
	}
	else
	{
		tempalpha = NULL;
		templast_nlr_chi2fr = NULL;
		tempmaxdec_nlr_chi2fr = NULL;
	}

	if (ExptsData::is_E0shift)
	{
		SetArraysize(&old_ek_gridind, nek, "old_ek_gridind", "main");
		for (i = 0; i < nek; i++)
			old_ek_gridind[i] = ExptsData::ek_gridind[i];//default the original 
	}

	if (mode || (ExptsData::is_AXS || ExptsData::is_E0shift || ExptsData::is_IQmucorr))//program will run to given number of ngenerated steps, so the last step might be not acceptable, so the chi2
		//of the last acceptable steps had to be copied for displaying
		//or it is a hybrid run with possible no atomic moves
	{
		for (i = 0; i < ChiSquared::chisize; i++)
		{
			tempchi[i] = chisquare.chicomp[i];
		
		}
	
		for (i = 0; i < ngr + nsq + nfq + nfg + nek; i++)
		{
			tempa[i] = chisquare.a[i];
			tempb[i] = chisquare.b[i];
			tempc[i] = chisquare.c[i];
			tempd[i] = chisquare.d[i];
			tempe[i] = chisquare.e[i];
		}
	}
	else
	{
		//in case of nonlin regression they are needed anyhow
		if (ExptsData::is_IQ > 0)//
		{
			for (i = 0; i < nfq; i++)
			{
				tempa[ngr + nsq + i] = chisquare.a[ngr + nsq + i];
				tempb[ngr + nsq + i] = chisquare.b[ngr + nsq + i];
			}
		}
	}
	if (ExptsData::is_IQ > 0)//
	{
		for (i = 0; i < nfq; i++)
		{
			tempalpha[i] = chisquare.alpha[i];
			templast_nlr_chi2fr[i] = chisquare.last_nlr_chi2fr[i];
			tempmaxdec_nlr_chi2fr[i]=chisquare.maxdec_nlr_chi2fr[i];
			if (ExptsData::is_IQmucorr)
			{
				old_muact[i] = ExptsData::fqmuact[i];//has to be muact in case of continuation
				old_f[i] = chisquare.f[ngr + nsq + i];
				old_ff[i] = chisquare.ff[ngr + nsq + i];
				//old_fx[i] = chisquare.fx[ngr + nsq + i];//not used
				//old_fx2[i] = chisquare.fx2[ngr + nsq + i];
				old_s[i] = chisquare.s[ngr + nsq + i];
				old_fs[i] = chisquare.fs[ngr + nsq + i];
				old_ss[i] = chisquare.ss[ngr + nsq + i];
			}
			if (ExptsData::fqAXS[i]>0)
				old_fprimeact[i] = ExptsData::fqfprimeact[i*RunParams::ntypes+ExptsData::fqAXS[i]-1];
		}
	}
	
	
	if (History::calc_method>0)
	{
		cout << "\nNOTE(" << ++note << "): The chi2 calculation is done according to square-of-the-difference method for "<<endl;
		cout<<"\tsome of the data sets:    "<<endl;
		for (i=0;i<ngr;i++)
		{
			cout<<J4<<i+1<<"\t\tg(r) data series :   ";
			if (1&History::calc_method)
				cout<<"square of the difference"<<endl;
			else
				cout<<"normal"<<endl;

			History::calc_method=History::calc_method>>1;
		}
		for (i=0;i<nsq;i++)
		{
			cout<<J4<<i+1<<"\t\tS(Q) data series :   ";
			if (1&History::calc_method)
				cout<<"square of the difference"<<endl;
			else
				cout<<"normal"<<endl;

			History::calc_method=History::calc_method>>1;
		}
		for (i=0;i<nfq;i++)
		{
			cout<<J4<<i+1<<"\t\tF(Q) data series :   ";
			if (1&History::calc_method)
				cout<<"square of the difference"<<endl;
			else
				cout<<"normal"<<endl;

			History::calc_method=History::calc_method>>1;
		}
		for (i = 0; i < nfg; i++)
		{
			cout << J4 << i + 1 << "\t\tF(g) data series :   ";
			if (1 & History::calc_method)
				cout << "square of the difference" << endl;
			else
				cout << "normal" << endl;

			History::calc_method = History::calc_method >> 1;
		}
		for (i=0;i<nek;i++)
		{
			cout<<J4<<i+1<<"\t\tE(k) data series :   ";
			if (1&History::calc_method)
				cout<<"square of the difference"<<endl;
			else
				cout<<"normal"<<endl;

			History::calc_method=History::calc_method>>1;
		}
	}
	old_chitotal=chisquare.total;//the value of the total chi square
	
	if (RunParams::potential>0)
		old_chitotal2=chisquare.total2;
	if (chisquare.total2>0)
		chisquare.pot_chi2_low_limit=chisquare.total2*RunParams::pot_chi2_low_lim_fraction;
	else
		chisquare.pot_chi2_low_limit=chisquare.total2+(chisquare.total2*(1.0-RunParams::pot_chi2_low_lim_fraction));
	
	OpenFile(chifile,chifilename,"main",0);//open file, check, whether it was successfully opened
	chisquare.Save(chifile);
	
#ifdef _WRITE_CHI2_DETAIL
		chifile.precision(6);
		chifile.setf(ios::scientific, ios::floatfield);
		chifile.setf(ios::right, ios::adjustfield);
		chifile << "\nChi2 details" << endl;
		chifile << "         \tngen\tntried\tnacc\tnpotacc\tnswapgen\tswapacc\tnE0shiftgen\tnE0shiftacc\tmucorrgen\tmucorracc\tf'shiftgen\tf'shiftacc\tchi component as in *.hst" << endl;
		chifile << "first    \t " << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;
		for (i=0;i<chisquare.chisize;i++)
			chifile<<"\t"<<chisquare.chicomp[i];
		chifile<<endl;
		
#else		
	chifile.close();//no detailed chi2 will be saved in case of rejection, file is closed
#endif

	
	//Create History object
	History::SetHistoryParams(rundata);//sets the static members
	History history;
	if (history.buffsize>0)//if history is recorded
	{
		if (RunParams::continuation)
			CleanOpen(hstfile,hstfilename,2);//open in append mode
		else
			CleanOpen(hstfile,hstfilename);//open in overwrite mode
	
		history.Inifile(hstfile,rundata,edata);//initialise the history file
	
#ifdef _TEST_MODE
#ifdef	WRITE_TH_DURATION
		hstfile.setf(ios::scientific, ios::floatfield);
		hstfile.setf(ios::right, ios::adjustfield);

		for (ithread = 0; ithread < nthreads - 1; ithread++)
		{
			hstfile << "dur_" << ithread << "0 \t" << Threads::dur_10[ithread] / 1000 << endl;
			hstfile << "dur_" << ithread << "1 \t" << Threads::dur_11[ithread] / 1000 << endl;
			hstfile << "dur_" << ithread << "2 \t" << Threads::dur_12[ithread] / 1000 << endl;
#ifdef _AENET
			hstfile << "dur_ann1[" << nthreads-1 << "] \t" << Threads::dur_ann1[nthreads - 1] / 1000 << endl;
#endif
		}
		hstfile << "dur_" << nthreads - 1 << "0 \t" << Threads::dur_00 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "1 \t" << Threads::dur_01 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "2 \t" << Threads::dur_02 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "3 \t" << Threads::dur_03 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "4 \t" << Threads::dur_04 / 1000 << endl;
#ifdef _AENET
			hstfile << "dur_ann1[" << nthreads-1 << "] \t" << Threads::dur_ann1[nthreads - 1] / 1000 << endl;
#endif
		hstfile.unsetf(ios::scientific);
		hstfile.unsetf(ios::right);
		
#endif
#endif
		duration=0;
		i=0;
		history.Buffline(duration,rundata,chisquare,edata);//put the initial chi values into the buffer
		if (history.fill>=history.buffsize)
			//it is time to flush the history buffer to the file
			history.Save(hstfile);//fill is reset to zero
	}
	
	//-------------------------------------------------------------------
	//initialisation of loop structures
	//-------------------------------------------------------------------

	logfile.precision(16);
	logfile.setf(ios::scientific, ios::floatfield);
	

	duration=0;//duration of the main cycle
	nsaves=0;//number of saves since last history buffering
	
	if (!RunParams::continuation)
	{
		ngenerated=0;
		ntried=0;
		naccepted=0;
		npotaccepted=0;
	}


	
	
	switch (rundata.runmode)
	{
	case 1:
		run_limit = -rundata.runlimit;//negative value means number of generated steps, run_limit always positive
		present = &nsteps_d;//here this will be ngenerated
		nsteps_d = 0;
		break;
	case 2:
		run_limit = rundata.runlimit;//to naccepted
		present = &nsteps_d;//here it will be naccepeted
		break;
	default:
		run_limit = rundata.runlimit * 60;//limiting duration (in seconds)
		present = &duration;

	}
	time_save=rundata.timesave*60;//saving period in seconds
	cfg_collect=rundata.cfgnumb;//number of configurations to collect
	cfg_frequency=rundata.coll_frequency;//collect the configuration at each coll_frequency save, if 0, not collecting
	
	time(&t_start);//initialise the start time value
	time(&t_lastsave);//initialise the last save time value
#ifdef _ADVANCED_GEOM_CONST
	rundata.PrintStatus(edata,chisquare,iccnew,avccnew, concnew, sncnew, bvsnew, config,0);
#else
	rundata.PrintStatus(edata, chisquare, iccnew, avccnew, config, 0);
#endif

	rundata.SaveState(chisquare);//save the data for continuation
	calcdat.SaveResult();//saving the total
	if (nthreads>1)//signal to threads to call ThreadLoop, (all the necessary pointer was initialited)
	{
		std::unique_lock<std::mutex> guard(thread_obj.start3_mutex);//lock mutex
		thread_obj.start_flag=1;//time to proceed
		thread_obj.check_loop_start.notify_all();
	}

#ifdef _TEST_MODE
	//for recording the speed of certain parts of the program
	dur0=0;
	dur1=0;
	dur2=0;
	dur3=0;
	dur3a=0;
	dur3b=0;
	dur4=0;
	dur5=0;
	dur6=0;
	dur7=0;

	Threads::dur_00=0;
	Threads::dur_00a=0;
	Threads::dur_01=0;
	Threads::dur_02=0;
	Threads::dur_03=0;
	Threads::dur_04=0;
	Threads::dur_05=0;
	Threads::dur_06=0;

	for (ithread=0;ithread<nthreads;ithread++)
	{
			Threads::dur_10[ithread]=0;
			Threads::dur_11[ithread]=0;
			Threads::dur_12[ithread]=0;
			Threads::dur_12a[ithread]=0;
			Threads::dur_12b[ithread]=0;
			Threads::dur_13[ithread]=0;
			Threads::dur_13a[ithread]=0;
			Threads::dur_13b[ithread]=0;
			Threads::dur_14[ithread]=0;
			Threads::dur_15[ithread]=0;
			Threads::dur_15b[ithread]=0;
			Threads::dur_16[ithread]=0;
			Threads::dur_17[ithread]=0;
#ifdef _ADVANCED_GEOM_CONST
			Threads::dur_18[ithread] = 0;
#endif
#ifdef _USE_LOCAL_INV
			Threads::dur_loc1[ithread]=0;
			Threads::dur_loc2[ithread]=0;
			Threads::dur_loc3[ithread]=0;
			Threads::dur_loc4[ithread]=0;
#endif
#ifdef _LOCAL_INV
			Threads::dur_chiloc1[ithread]=0;
			Threads::dur_chiloc2[ithread]=0;
			Threads::dur_chiloc3[ithread]=0;
			Threads::dur_chiloc4[ithread]=0;
#endif
		
	}
#endif
	/*
	logfile<<chisquare.total<<" old "<<old_chitotal<<   endl;
	for (int ic = 0; ic < chisquare.chisize; ic++)
		logfile << "chi" << ic << " " << chisquare.chicomp[ic] << endl;
		*/
  	cout<<"\n====================================================================================================================="<<endl;
	cout << "===========================================                               ===========================================" << endl;
	cout<<  "===========================================  METROPOLIS LOOP HAS STARTED  ==========================================="<<endl;
#ifdef _ATLAS 
	cout<<  "=================================  ATLAS is used for matrix-vector multiplication!  ================================="<<endl;
#endif
	cout << "===========================================                               ===========================================" << endl;
	cout<<"=====================================================================================================================\n"<<endl;
	
	HistoSet::ngen = &ngenerated;
#ifdef _TEST_MODE //program runs in test mode
	seed=-100000;//initialising the random number generator to start the same way at each run
	cout << "RUNNING IN TEST MODE till ";
	if (RunParams::runmode==1)
		cout<<"ngenerated = "<<run_limit<<" !!!"<<endl;
	else
		cout << "naccepted = " << run_limit << " !!!" << endl;
	

	/*logfile.precision(15);//GO
	logfile.setf(ios::scientific, ios::floatfield);
	logfile << "vdWbeg ";
	for (i = 0; i < Move::npartials; i++)
		logfile << " " << FNC_POT::vdW_pot[i];
	logfile << endl;
	logfile << "Coulbeg";
	for (i = 0; i < Move::npartials; i++)
		logfile << " " << FNC_POT::Coulomb_pot[i];
	logfile << endl;
	logfile.precision(6);//GO
	logfile.unsetf(ios::scientific);
	logfile << "potacc " << acceptable << endl;*/

	
	while (*present<run_limit && !forced_termination)//the run will end reaching ngenerated or naccepeted = run_limit depending on the mode, configurations are not collected
#else //normal running mode
	seed=(longint)-time(nullptr);//initialising the random number generator
	//if there are configurations to collect, collection starts when duration reached time limit
	//and the run will end, when the required number of configurations are collected
	while ((*present < run_limit || cfg_collect) && !forced_termination)//present is either duration or ngenerated_d
#endif	

	{

#ifdef _TEST_MODE
		time1 = std::chrono::high_resolution_clock::now();
#endif

		ngenerated++;//increase the number of generated moves
		//logfile << "\n"<<ngenerated << endl;

		if (log_nlr_steps)
			logfile << "\nngenerated " << ngenerated << endl;

		if (RunParams::runmode==1)
			nsteps_d = (double)ngenerated;//convert it to double in case it used in the loop length determination
		
		
		history_buffered = 0;//History was not buffered for this step by default

		Move::nomove = false;//there is atomic move by default
		if (ExptsData::is_E0shift && (ngenerated % RunParams::exafs_shiftstep == 0))
		{
			Move::makeE0shift = 1;
			Move::nomove = true;//to make the if statements shorter
			nE0shiftgen++;
		}
		else
			Move::makeE0shift = 0;
		
		if (ExptsData::is_IQmucorr && (RunParams::IQ_backg_corr_step>0) && (ngenerated % RunParams::IQ_backg_corr_step == 0))
		{
			Move::makemucorr = 1;
			Move::nomove = true;//to make the if statements shorter
			nmucorrgen++;

		}
		else
			Move::makemucorr = 0;

		if (ExptsData::is_AXS && (ngenerated % RunParams::AXS_shiftstep == 0))
		{
			Move::makefprimeshift = 1;
			Move::nomove = true;//to make the if statements shorter
			nfprimeshiftgen++;
		}
		else
			Move::makefprimeshift = 0;
#ifdef _AENET
		if (is_ANN && (ngenerated % RunParams::aenet_step == 0))
		{
			if (RunParams::aenet_step==1)
			{
				if (Move::nomove)
					Aenet::calc_ANN=0;
				else
					Aenet::calc_ANN=1;
			}
			else
				Aenet::calc_ANN=2;//not calculated in each step, needed to differentiate which array is passed to SplitAenetAtoms
		}
		else
			Aenet::calc_ANN=0;
			//logfile<<"aeant calc "<<Aenet::calc_ANN<<endl;
#endif
		//create the random move (for the nmoved atoms involved)
		if (RunParams::custmove == 0)
			move.MakeMove(config, seed, ((RunParams::potential == 1 && RunParams::fnc == 4) ? &fnc : 0));//standard RMC++
		else
			move.MakeMove(config, fnc, seed);//custom move
		/*logfile << "moved";
		for (int mi = 0; mi < move.tot_moved_atoms; mi++)
			logfile << " " << move.indices[mi];
		logfile << endl;*/
		if (Move::makemucorr)
		{
			//logfile << "mucorr " << edata.fqmuact[0] << endl;
			edata.CalcIQmucorr();
		}
		if (Move::makefprimeshift)
			edata.CalcXrayCoeffsChange();//calculate the change in the coeffs
		
		if (Move::swap)
			nswapgen++;//Increasing the number of generated swaps

		/*logfile << "\nngen "<<ngenerated<<" moved ";//GO

		if (!Move::nomove)
		{
			for (i = 0; i < Move::tot_moved_atoms + Move::tot_moved_virtuals; i++)
				logfile << " " << move.indices[i];
			logfile << endl;
		}
		else
			logfile<<" AXS "<<Move::makefprimeshift<<endl;*/
		/*logfile.precision(15);//GO
		logfile.setf(ios::scientific, ios::floatfield);
		logfile << "vdWbefcalc ";
		for (i = 0; i < Move::npartials; i++)
			logfile << " " << FNC_POT::vdW_pot[i];
		logfile << endl;
		logfile << "Coulbefcalc ";
		for (i = 0; i < Move::npartials; i++)
			logfile << " " << FNC_POT::Coulomb_pot[i];
		logfile << endl;
		logfile.precision(6);//GO
		logfile.unsetf(ios::scientific);*/

		//check, whether the FNC is satisfied, if necessary
		if (RunParams::fnc > 0 && RunParams::fnc != 4 && !Move::nomove)
		{
			if (!fnc.CheckFNCChange(move))//Check the FNC among the moved atoms and its neighbours
			{
				move.ResetCoord(extended_fnc ? &fnc : 0);//copy back the old coordinates
				continue;//the move is not acceptable, starts the main loop again, and makes a new move
			}
		}

		if (!Move::nomove)
		{
#ifdef _AENET
			//collect the old neighbours of the moved atoms, even if E is not calculated in this setp 
			Aenet::nchanged=0;
			aenetold.Update(move, 0,thread_obj.thread_arg[nthreads-1]);//this is serial, only the arrays are used		
#endif
			//check, whether the cutoffs are satisfied among the moved atoms, and the other atoms in their vicinity
			neighlist.UpdateGrid(move);//update the grid to reflect the new positions created by the move

#ifdef _TEST_MODE
			time2 = std::chrono::high_resolution_clock::now();
			elapsed = time2 - time1;
			dur0+=elapsed.count();
#endif

			//checking, whether the move is acceptable regarding the cutoff distances
			if (!neighlist.CheckCutoff(move))
			{
				move.ResetCoord(extended_fnc ? &fnc : 0);//copy back the old coordinates
				neighlist.ResetGrid(move.indices);//restore the grid to the original state

#ifdef _TEST_MODE
				time1 = std::chrono::high_resolution_clock::now();
				elapsed = time1 - time2;
				dur1+=elapsed.count();
#endif
				continue;//the move is not acceptable on the ground of the cutoff distance, make a new move
			}

#ifdef _AENET
			//collect the new neighbours of the moved atoms, and calculate the new energy for the moved atoms, is necessary
			aenetnew.Update(move, 1,thread_obj.thread_arg[nthreads-1]);//this is serial, only the main threads's array segments are used
			if (Aenet::calc_ANN==2)
				Aenet::UpdateStored();//add the new neighbours to the stored array
#endif

#ifdef _TEST_MODE
			time1 = std::chrono::high_resolution_clock::now();
			elapsed = time1 - time2;
			dur1+=elapsed.count();
#endif
			//Calculate the change in the histogram
			histnew.HistCalcChange();//here no cancellation can happen anymore!
			//HistCalcChange will call the aenetnew.CalcANN for the main thread
#ifdef _TEST_MODE
			time2 = std::chrono::high_resolution_clock::now();
			elapsed = time2 - time1;
			dur2+=elapsed.count();
#endif
#ifdef _AENET
			//now the ANN energy is calculated, calculate the total E 
			if (Aenet::calc_ANN)
			{
				aenetnew.E_tot=0;
						
				for (i=0;i<Move::tot_moved_atoms;i++)
					aenetnew.E_i[move.indices[i]]=Aenet::E_moved[i];
				for (ithread=0;ithread<nthreads;ithread++)
				{
					for (i=0;i<thread_obj.thread_arg[ithread].aenet_at_count;i++)
						aenetnew.E_i[thread_obj.thread_arg[ithread].aenet_at_index[i]]=thread_obj.thread_arg[ithread].aenet_E_i[i];
				
				}
				//calculate the total
				for (i=0;i<Aenet::nAtoms;i++)
					aenetnew.E_tot+=aenetnew.E_i[i];
			}
#endif
		}
		
		ntried++;
		//logfile<<"Ntried "<<ntried<<endl;
		if (!Move::nomove)
		{

			//it can be used by chisquare.CalcPotChiSquared
			chisquare.total = 0.0;//the normal chi2 total, it is needed, if there is no non-bonded interaction

			if (RunParams::potential == 1)
			{
				if (Topology::nexclusion_max > 0)
				{
					/*logfile.precision(15);//GO
					logfile.setf(ios::scientific, ios::floatfield);
					logfile << "bef excl "<<FNC_POT::vdW_pot[0] <<endl;
					logfile.precision(6);//GO
					logfile.unsetf(ios::scientific);*/

					pfnc[0].CalcNBExclusion(move);//subtract the exclusions, and update the 1-4 interaction, if there is any

					/*logfile.precision(15);//GO
					logfile.setf(ios::scientific, ios::floatfield);
					logfile << "aft excl " << FNC_POT::vdW_pot[0]<<endl;
					logfile.precision(6);//GO
					logfile.unsetf(ios::scientific);*/
				}

				pfnc[0].UpdateTotPot();//calculating the total 
			}
			if (RunParams::fnc == 4)
			{
				//calculating the bonded interactions
				for (i = 0; i < N_POT_TYPE; i++)
					(pfnc[i].*InteractionFunctChange[i])(move);

			}


			//If there is non-bonded potential or ANN, the chi2 components coming from all the bonded and non-bonded interactions
			//will be summed separately from the data sets, as the non-bonded potental can be negative
			//Decide, whether the move is acceptable, based on this total2
			//if the moveout option is on, and there is no moved atom, than the Metropolis test is performed, if moved atom
			//is involved, than the move is acceptable, regardless the chi2, so no check is performed
			//If there is not non-bonded-potential only bonded, than the chi2 contributions of them is added to the normal chi2

			if (extended_fnc || is_ANN)
				chisquare.CalcPotChiSquared();
				

			if (RunParams::potential > 0 && move.ntcl_selected < 1)//if this was not a 'moved out' atom, then  
			{			//the Metropolis test must be performed
#ifdef _NO_POT_CHI2	
			//the floating point errors of the potential calculation can be slightly different in case of different thread numbers
			//and can cause, that in one case a move is accepted, in the other is not, so the results cannot be compared
			//therefor in test mode every 10th step is accepted regardless the potential

				if (ngenerated / 10.0 != (double)int(ngenerated / 10.0))
					acceptable = 0;

#else

				chidiff = 0.5 * (old_chitotal2 - chisquare.total2);
				if (chidiff < -20)
					acceptable = 0;
				else
				{
					if (chidiff < 0)
					{
						if (exp(chidiff) < Ran1(seed))
							acceptable = 0;
					}
					else
					{
						if (RunParams::pot_chi2_low_lim_fraction < 1 && chisquare.total2 < ChiSquared::pot_chi2_low_limit)
							acceptable = 0;//do not let the chi2_pot go under the limit
					}
				}
#endif
			}


			/*logfile.precision(15);//GO
			logfile.setf(ios::scientific, ios::floatfield);
			logfile << "vdW ";
			for (i = 0; i < Move::npartials; i++)
				logfile << " " << FNC_POT::vdW_pot[i];
			logfile << endl;
			logfile << "Coul";
			for (i = 0; i < Move::npartials; i++)
				logfile << " " << FNC_POT::Coulomb_pot[i];
			logfile << endl;
			logfile.precision(6);//GO
			logfile.unsetf(ios::scientific);
			logfile << "potacc " << acceptable << endl;*/


			if (!acceptable)
			{
				//logfile<<"pot reject"<<endl;
				// if the move has been rejected then the modified parts during this try 
				//must be restored to their original state
				move.ResetCoord(&fnc);//copy back the old coordinates
#ifdef _NO_PERIODIC
				move.ResetBinInd();//reset the central bin indices to the original state
#endif

				neighlist.ResetGrid(move.indices);//restore the grid to the original state

				avccold.CopyModified(avccnew);//Copy the modified av. coordination numbers 

				//------------------MULTI-THREADING-------------------
				iccold.CopyModified(iccnew, thread_obj.thread_arg[nthreads - 1]);//Copy the modified coordination numbers 
#ifdef _ADVANCED_GEOM_CONST
				bvsold.CopyModified(bvsnew, thread_obj.thread_arg[nthreads - 1]);//Copy the modified valences
#endif
				histold.CopyModified(histnew, thread_obj.thread_arg[nthreads - 1]);//Copying the modified histogram parts

#ifdef _AENET
				if (Aenet::calc_ANN)
					aenetold.Copy(aenetnew,move.indices,thread_obj.thread_arg[nthreads - 1]);//Copy the modified energies back

#endif		
				//------------------END OF MULTI-THREADING-------------------


				for (i = 0; i < n_fnc; i++)
					pfnc[i].UpdatePot(1);//Copying the potential components from the 'old' arrays to restore them to the values of the 
										 //last accepted move
#ifdef _WRITE_CHI2_DETAIL
				chifile << "pot_based\t " << ngenerated << "\t" << ntried << "\t" << naccepted << "\t" << npotaccepted << "\t" << nswapgen << "\t" <<
					nswapacc << "\t" << nE0shiftgen << "\t" << nE0shiftacc << "\t" << nmucorrgen << "\t" << nmucorracc << "\t" << nfprimeshiftgen << "\t" << nfprimeshiftacc ;

				for (i = 0; i < chisquare.chisize; i++)
					chifile << "\t"<< chisquare.chicomp[i] ;
				chifile << endl;

#endif			
				//------------------MULTITHREADING-------------------
				//thread_obj.complete_count2 = 0;//counter indicating how many thread finished with the calculation, reset it to 0
#ifdef _TH_ORDER
				cout << "r" << endl;
#endif
				if (nthreads > 1)//signal to threads to start calculation 
				{
					std::unique_lock<std::mutex> guard(thread_obj.start2_mutex);//lock mutex
					thread_obj.start_flag = 1;//time to proceed
					thread_obj.check_cd_start.notify_all();
#ifdef _TH_ORDER
					cout << "-" << endl;
#endif
				}

				continue;//the move is not acceptable, let start the main loop again, and make a new move
			}//end of if the move was rejected

			if (RunParams::potential > 0)
				npotaccepted++;

#ifdef _NO_PERIODIC
			//------------------MULTITHREADING-------------------
			thread_obj.NOP_count = 0;//counter indicating how many thread finished with the calculation, reset it to 0

			if (nthreads > 1)//signal to threads to start calculation 
			{
				std::unique_lock<std::mutex> guard(thread_obj.start2_mutex);//lock mutex
				thread_obj.start_flag = 1;//time to proceed
				thread_obj.check_cd_start.notify_all();
			}


			if (datmat.calcmode == 0)
				histnew.CalcPerHist(thread_obj.thread_arg[nthreads - 1], datmat);//calculate the peiodic histogram, if necessary	

#ifdef _TEST_MODE
			time1 = std::chrono::high_resolution_clock::now();
			elapsed = time1 - time2;
			dur3+=elapsed.count();
#endif
			//Threads have to wait, till all of them finished, as each of them will use all the PPCF values
			if (nthreads > 1)
			{
				std::unique_lock<std::mutex> guard(thread_obj.NOP_count_mutex);//lock mutex
				thread_obj.NOP_count++;
				if (thread_obj.NOP_count == nthreads)//time to proceed
				{
					thread_obj.start_flag = 0;//reset it to 0
					thread_obj.NOP_count = 0;//reset it to 0
					thread_obj.NOP_check_count.notify_all();
				}
				else thread_obj.NOP_check_count.wait(guard);
			}
#endif //_NO_PERIODIC


#ifdef _VIBR_AMP		
			histnew.ThermalCorr();//calculate the convoluted histogram to include atomic vibrations
#endif
		
		
			//------------------MULTITHREADING-------------------
			thread_obj.CMP_count=0;//counter indicating how many thread finished with the calculation, reset it to 0
#ifdef _TH_ORDER
			cout<<"z"<<endl;
#endif
			if (nthreads>1)//signal to threads to start calculation 
			{
				std::unique_lock<std::mutex> guard(thread_obj.start2_mutex);//lock mutex
				thread_obj.start_flag=1;//time to proceed
				thread_obj.check_cd_start.notify_all();
#ifdef _TH_ORDER
				cout << "/" << endl;
#endif
			}

			//calculate the modidified parts of the ppcf
			ppcfnew.CalcModPPCF(thread_obj.thread_arg[nthreads-1]);

#ifdef _TEST_MODE
			time1 = std::chrono::high_resolution_clock::now();
			elapsed = time1 - time2;
			dur3+=elapsed.count();
#endif
			//Threads have to wait, till all of them finished, as each of them will use all the PPCF values
			if (nthreads>1)
			{
				std::unique_lock<std::mutex> guard(thread_obj.CMP_count_mutex);//lock mutex
				thread_obj.CMP_count++;
#ifdef _TH_ORDER
				cout << "t" << endl;
				logfile << "CMP " << thread_obj.CMP_count << endl;
#endif
				if (thread_obj.CMP_count==nthreads)//time to proceed
				{
					thread_obj.start_flag=0;//reset it to 0
					thread_obj.CMP_count=0;//reset it to 0
					thread_obj.CMP_check_count.notify_all();
#ifdef _TH_ORDER
					cout << "!" << endl;
#endif
				}
				else
				{
					thread_obj.CMP_check_count.wait(guard);
#ifdef _TH_ORDER
					cout << "+"<<endl;
#endif
				}
#ifdef _TH_ORDER
				logfile << "ac " << thread_obj.CMP_count<<endl;
#endif
			}
		}
		else//no atomic move, no histogram calculation, set here some variables 
		 //which are otherwise set in HistCalcChange and handle threading to remain synchronized with ThreadLoop
		{
			if (nthreads > 1)//signal to threads to start calculation 
			{
				//wait, if not all the new threads finished with the previous loop
				std::unique_lock<std::mutex> guard(thread_obj.count_mutex);//lock mutex
				thread_obj.start_count++;
#ifdef _TH_ORDER
				cout << "u" << endl;
				logfile << "IHCC " << thread_obj.start_count << endl;
#endif
				if (thread_obj.start_count == nthreads)//time to proceed, this was the last thread to reach this point
				{
					acceptable = 1;//the move is acceptable, if it does not prove otherwise
					if (Move::makeE0shift)
					{
						Threads::is_E0_shift = 1;//this is needed that the actual value should be preserved for the duration of the thread loop
					}
					else
						Threads::is_E0_shift = 0;
					if (Move::makemucorr)
					{
						Threads::is_IQ_mucorr = 1;
					}
					else
						Threads::is_IQ_mucorr = 0;

					if (Move::makefprimeshift)
					{
						Threads::is_fprime_shift = 1;
					}
					else
						Threads::is_fprime_shift = 0;
					for (i = 0; i < nek; i++)
						Threads::EXAFS_gridind[i] = ExptsData::ek_gridind[i];

					for (i = 0; i < ExptsData::nfq; i++)
						Threads::IQ_muact[i] = ExptsData::fqmuact[i];

					for (i = 0; i < ExptsData::nfq; i++)
					{
						if (ExptsData::fqAXS[i] > 0)
							Threads::IQ_fprimeact[i] = ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1];
					}
#ifdef _AENET
					//if ANN is not calculated in each step, then if there were accepted moves without ANN
					//calculation since the last ANN calculation, then has to calculate ANN
					if (Aenet::last_gen_accepted != Aenet::last_gen_aenet_accepted && Aenet::calc_ANN == 2)
					{
						Threads::calc_ANN = Aenet::calc_ANN;
						//split has to be done here after all the threads finished with updates
						thread_obj.SplitAenetAtoms(Aenet::nstored, Aenet::aenet_stored_ind);//use the stored indices, al the atoms involved in change since last calc
					}
					else
						Threads::calc_ANN = 0;//no need to calculate if ANN is calculated in each step
#endif
					thread_obj.start_flag = 0;//reset it for later use
					thread_obj.start_count = 0;//reset it to 0
					if (!RunParams::custmove && SimpleCfg::ntypes > 1)//for custmove or ntypes=1 no update is necessary, as all the partials and types are modified
					{
						for (i = 0; i < Move::npartials; i++)
							Threads::mod_partial[i] = 0;

						if (Move::makeE0shift)//technically only those partialy modified, which contains the edge particle, but the others are not used
						{
							for (i = 0; i < Move::npartials; i++)
								Threads::mod_partial[i] = 1;
						}
						//mucorr and f' shift does not alter the partial F(Q)-s
					}

					thread_obj.check_count.notify_all();
#ifdef _TH_ORDER
					cout << "$" << endl;
#endif
				}
				else thread_obj.check_count.wait(guard);
			}
			else
			{
				//Compiled as MULTI, but used only with nthreads==1
				acceptable = 1;//the move is acceptable, if it does not prove otherwise

				if (Move::makeE0shift)
				{
					Threads::is_E0_shift = 1;//this is needed that the actual value should be preserved for the duration of the thread loop
				}
				else
					Threads::is_E0_shift = 0;
				if (Move::makemucorr)
				{
					Threads::is_IQ_mucorr = 1;
				}
				else
					Threads::is_IQ_mucorr = 0;

				if (Move::makefprimeshift)
				{
					Threads::is_fprime_shift = 1;
				}
				else
					Threads::is_fprime_shift = 0;
				for (i = 0; i < nek; i++)
					Threads::EXAFS_gridind[i] = ExptsData::ek_gridind[i];

				for (i = 0; i < ExptsData::nfq; i++)
					Threads::IQ_muact[i] = ExptsData::fqmuact[i];

				for (i = 0; i < ExptsData::nfq; i++)
				{
					if (ExptsData::fqAXS[i] > 0)
						Threads::IQ_fprimeact[i] = ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1];
				}
#ifdef _AENET
				//if ANN is not calculated in each step, then if there were accepted moves without ANN
				//calculation since the last ANN calculation, then has to calculate ANN
				if (Aenet::last_gen_accepted != Aenet::last_gen_aenet_accepted && Aenet::calc_ANN == 2)
				{
					Threads::calc_ANN = Aenet::calc_ANN;
					//split has to be done here after all the threads finished with updates
					thread_obj.SplitAenetAtoms(Aenet::nstored, Aenet::aenet_stored_ind);//use the stored indices, al the atoms involved in change since last calc
				}
				else
					Threads::calc_ANN = 0;//no need to calculate if ANN is calculated in each step

#endif
				if (!RunParams::custmove && SimpleCfg::ntypes > 1)//for custmove or ntypes=1 no update is necessary, as all the partials and types are modified
				{
					for (i = 0; i < SimpleCfg::npartials; i++)
						Threads::mod_partial[i] = 0;//decide which should be calculated at all
					if (Move::makeE0shift)//technically only those partials are modified, which contain the edge particle, but the others are not used
						for (i = 0; i < Move::npartials; i++)
							Threads::mod_partial[i] = 1;
					//mucorr and f' shift does not alter the partial F(Q)-s
				}

			}
#ifdef _AENET
			if (Aenet::calc_ANN)
				aenetnew.CalcANN(thread_obj.thread_arg[nthreads - 1]);
#endif

			//this would be the end of histogram change
			if (nthreads > 1)
			{
				std::unique_lock<std::mutex> guard(thread_obj.CC_count_mutex);//lock mutex
				thread_obj.CC_count++;
#ifdef _TH_ORDER
				cout << "y" << endl;
				logfile << "ICC " << thread_obj.CC_count << endl;
#endif
				if (thread_obj.CC_count == nthreads)//time to proceed
				{
					thread_obj.start_flag = 0;//reset it to 0
					thread_obj.CC_count = 0;//reset it to 0
					thread_obj.CC_check_count.notify_all();
#ifdef _TH_ORDER
					cout << "," << endl;
#endif
				}
				else thread_obj.CC_check_count.wait(guard);
			};
			//------END OF MULTI-THREADING------
#ifdef _AENET
			if (Aenet::last_gen_accepted != Aenet::last_gen_aenet_accepted && Aenet::calc_ANN == 2)
			{
				//now the ANN energy is calculated, calculate the total E 
				aenetnew.E_tot = 0;


				for (ithread = 0; ithread < nthreads; ithread++)
				{
					for (i = 0; i < thread_obj.thread_arg[ithread].aenet_at_count; i++)
						aenetnew.E_i[thread_obj.thread_arg[ithread].aenet_at_index[i]] = thread_obj.thread_arg[ithread].aenet_E_i[i];
				}
				//calculate the total
				for (i = 0; i < Aenet::nAtoms; i++)
					aenetnew.E_tot += aenetnew.E_i[i];
				//logfile<<ngenerated<<" ANN Etot old "<<aenetold.E_tot<<" new "<<aenetnew.E_tot<<endl;
				chisquare.CalcPotChiSquared();
			}
#endif
			
		}//end of no atomic move
#ifdef _TEST_MODE
		time2 = std::chrono::high_resolution_clock::now();
		elapsed = time2 - time1;
		dur3a+=elapsed.count();
#endif
																
		//calculate the modidified partial g(r) in case of g(r) fitting and S(Q) for S(Q) fitting and  F(g) for F(g) fitting and E(k) for E(k) fitting
		//in case of no atomic move calculating for E0 shift sets, Mucorr and f' shift does not alter the partials
		calcpnew.CalcModPartial(thread_obj.thread_arg[nthreads-1]);

#ifdef _TEST_MODE
			time1 = std::chrono::high_resolution_clock::now();
			elapsed = time1 - time2;
			dur3b+=elapsed.count();
#endif

		//Calculate the new calculated data CalcData for g(r), S(Q) F(Q) F(g) and E(k) fitting
		//in case of no atomic move calculating for E0 shift sets, f' shift sets, for mucorr I(Q)  only will change later during nonlinreg if apha renorm is used 
		calcdat.CalcTotal(thread_obj.thread_arg[nthreads-1]);//calculate the total from the partials

#ifdef _TEST_MODE
		time2 = std::chrono::high_resolution_clock::now();
		elapsed = time2 - time1;
		dur4+=elapsed.count();
#endif
		//waiting for all threads to finish calculation 
		if (nthreads>1)
		{
			std::unique_lock<std::mutex> guard(thread_obj.T_count_mutex);//lock mutex
			thread_obj.T_count++;
#ifdef _TH_ORDER
			cout << "v" << endl;
			logfile << "T " << thread_obj.T_count << endl;
#endif
			if (thread_obj.T_count==nthreads)//time to proceed
			{
				thread_obj.start_flag=0;//reset it to 0
				thread_obj.T_count=0;//reset it to 0
				thread_obj.T_check_count.notify_all();
#ifdef _TH_ORDER
				cout << "~" << endl;
#endif
			}
			else 
				thread_obj.T_check_count.wait(guard);
		}
		//------------------END OF MULTITHREADING-------------------
		if (!Move::nomove)
		{
			//update the CosDistrConst, if there is a constraint
			if (CosDistrConst::nconstraints > 0)
				cosnew.UpdateHist(move, -1);//removing the contribution of the atoms effected by the move using the original neighbourlist
#ifdef _ADVANCED_GEOM_CONST
			if (CommonNeighConst::nconstraints > 0)
				concnew.UpdateHist(move, -1);
			if (SecondNeighConst::nconstraints > 0)
				sncnew.UpdateNeigh(move, -1);
#endif
			if (calc_neighlist)
				neighlist.UpdateNeighList(move);//update the neighbourlist
			if (CosDistrConst::nconstraints > 0)
			{
				cosnew.UpdateHist(move, +1);//adding the contribution of the atoms effected by the move using the updated neighbourlist
				cosnew.CalcDistribution();//calculating the distribution
			}
#ifdef _ADVANCED_GEOM_CONST
			if (CommonNeighConst::nconstraints > 0)
				concnew.UpdateHist(move, +1);
			if (SecondNeighConst::nconstraints > 0)
				sncnew.UpdateNeigh(move, +1);
#endif

		}
		//Calculate the chi square value
#ifdef _ADVANCED_GEOM_CONST
		History::calc_method=chisquare.CalcChiSquared(calcdat,iccnew,avccnew,concnew,sncnew,bvsnew);
#else
		History::calc_method = chisquare.CalcChiSquared(calcdat, iccnew, avccnew);
#endif
		
		if (History::calc_method != -1)//it was normally calculated, and if non-lin regression was used, it has convereged
		{
			//------------------------------------------------
			//decide whether to accept move
			//------------------------------------------------
			if (move.ntcl_selected < 1)//if this was not a 'moved out' atom, then  
			{			//the Metropolis test must be performed
				chidiff = 0.5 * (old_chitotal - chisquare.total);

				if (chidiff < -20)
					acceptable = 0;
				else
				{
					if (chidiff < 0)
					{
						if (exp(chidiff) < Ran1(seed))
							acceptable = 0;
					}
				}
			}
			else
			{
				for (i = 0; i < move.ntcl_selected; i++)
				{
					//this was a "too close" atom, the move is accepted, regardless the chi square
					//to make sure, that this tooclose atom was not yet removed by RemoveTooClose as the tooclose pair 
					//of the (0->i-1)-th moved tooclose atom in this step
					if (move.tooclose_index[i] == move.tclindices[move.itooclose[i]])//the atom is still at its origonal place on the tooclose list
						current_index = move.itooclose[i];//index of the tooclose atom on the tooclose list
					else
					{
						//This tooclose atom was either removed from the list, or it was the last one and was swapped with an
						//other removed tooclose atom. The list has to be searched for this tooclose atom.
						current_index = -1;//not on the list by default
						for (k = 0; k < move.ntooclose; k++)//cycling through the atoms of the too close list
						{
							if (move.tclindices[k] == move.tooclose_index[i])
							{
								current_index = k;//index of the tooclose atom on the tooclose list
								break;
							}
						}
					}
					if (current_index > -1)//still on the list
					{
						if (move.tooclose_remove[i])//all the too close pairs disappeared for the given atom
						{
							//the atom now satisfies the cutoff, it can be removed from the list
							if (current_index == move.ntooclose - 1)//this is the last on the list
							{
								move.tclindices[move.ntooclose - 1] = -1;//erase the index, to prevent to cause any problems
								move.ntooclose--;//decrease the number of too close atoms
							}
							else
							{
								//replace the index of this atom with the last one on the list
								move.tclindices[current_index] = move.tclindices[move.ntooclose - 1];
								//replace the too close pair count of this atom with the last one on the list
								move.tclpaircount[current_index] = move.tclpaircount[move.ntooclose - 1];
								move.tclindices[move.ntooclose - 1] = -1;//erase the index, to prevent to cause any problems

								move.ntooclose--;//decrease the number of too close atoms
							}
						}//end of removing the moved atom from the too close list
						else
							move.tclpaircount[move.itooclose[i]] -= move.ntcl_pair[i];//decrease the number of bad pairs for the moved atom, but
											//the atom still has to remain on the list
					}//end of this too close atom was still on the list 
					//remove the other atom(s) of the too close pairs now above cutoff, if there is any

					if (move.ntcl_pair[i] > 0)
						move.RemoveTooClosePairs(i);//remove the other atom of the too close pairs now above cutoff

					if (move.ntooclose == 0)
					{
						cout << "\nNOTE(" << ++note << "): All atoms satisfy the cutoffs" << endl;
						move.moveout = 0;//switch off the moveout option
					}//end of else (i.e. if moved_out>=0)
				}
			}//end of if too close
		}
		else
		{
			acceptable = 0;//cannot be accepeted, non-lin regression failure
			if (chisquare.failed_nlr == rundata.terminate_nonlin)//the number of successive non-lin reg failures reached the limit, terminate
				forced_termination = true;
		}

#ifdef _TEST_MODE
		time1 = std::chrono::high_resolution_clock::now();
		elapsed = time1 - time2;
		dur5+=elapsed.count();
#endif

		//Updating the different arrays depending on whether the move was accepeted or rejected
		//signal to threads, that the copying can start
		if (nthreads>1) 
		{
			std::unique_lock<std::mutex> guard(thread_obj.start3_mutex);//lock mutex
			thread_obj.start_flag=1;//time to proceed
			thread_obj.check_cp_start.notify_all();
#ifdef _TH_ORDER
			cout << "^" << endl;
#endif
		}
		
		if (log_nlr_steps && !Move::nomove)
			logfile << "acceptable " << acceptable << endl;

		//logfile << "acceptable " << acceptable << endl;//GO
/*
#ifdef _AENET
		logfile << "b chi " << chisquare.total2 << " o " << old_chitotal2 << " etot " << aenetnew.E_tot << " olde" << aenetold.E_tot << endl;
#else
		logfile<<"b chi "<<chisquare.total2<< " o "<<old_chitotal2<< " etot "<<FNC_POT::vdW_pot[0]<<" olde"<< FNC_POT::vdW_pot[0] <<endl;
#endif*/
		/*logfile << " chitot " << chisquare.total << " old chi " << old_chitotal << endl;
		for (int ic = 0; ic < chisquare.chisize; ic++)
			logfile << "chi" << ic << " " << chisquare.chicomp[ic] << endl;*/
		if (acceptable)//if the move has been accepted
		{
			naccepted++;
			//logfile << ngenerated<<" acceptable " << naccepted << endl;
		
			if (RunParams::runmode == 2)
				nsteps_d = (double)naccepted;//convert it to double in case it used in the loop length determination
								
			if (ExptsData::is_IQ == 2)//there was non-linear regression, it converged
				ngoodnlr++;
			if (!Move::nomove)//only excute in case of atomic move
			{
				/*logfile.precision(15);//GO
				logfile.setf(ios::scientific, ios::floatfield);
				logfile << "accvdW ";
				for (i = 0; i < Move::npartials; i++)
					logfile << " " << FNC_POT::vdW_pot[i];
				logfile << endl;
				logfile << "accCoul ";
				for (i = 0; i < Move::npartials; i++)
					logfile << " " << FNC_POT::Coulomb_pot[i];
				logfile << endl;
				logfile.precision(6);//GO
				logfile.unsetf(ios::scientific);*/

				for (i = 0; i < FNC_POT::nconstraints; i++)
				{
					FNC_POT::nout_of_range[i] -= FNC_POT::dec_out_of_range[i];
					FNC_POT::sum_out_of_range -= FNC_POT::dec_out_of_range[i];
				}

				if (Move::swap)
					nswapacc++;//Increasing the number of accepted swaps
			}
			//The move was finally accepted, copying new object to the old
			old_chitotal=chisquare.total;//copying the new chi total to the old
			
			if (!Move::nomove)
			{
				avccnew.CopyModified(avccold);//Copy the modified av. coordination numbers

				if (CosDistrConst::nconstraints > 0)
					cosnew.CopyModified(cosold);//Copying the modified cos_hist histogram parts
#ifdef _ADVANCED_GEOM_CONST
				if (CommonNeighConst::nconstraints > 0)
					concnew.CopyModified(concold);//Copy the modified common neighbour histogram part
				if (SecondNeighConst::nconstraints > 0)
					sncnew.CopyModified(sncold);//Copy the modified second neighbour coordnumbs part
#endif

			//------------------MULTITHREADING-------------------
				iccnew.CopyModified(iccold, thread_obj.thread_arg[nthreads - 1]);//Copy the modified coordination numbers
#ifdef _ADVANCED_GEOM_CONST
				bvsnew.CopyModified(bvsold, thread_obj.thread_arg[nthreads - 1]);//Copy the modified valences
#endif
				histnew.CopyModified(histold, thread_obj.thread_arg[nthreads - 1]);//Copying the modified histogram parts
				ppcfnew.CopyModified(ppcfold, thread_obj.thread_arg[nthreads - 1]);//Copy the modified ppcf parts 
			}
#ifdef _AENET
			Aenet::last_gen_accepted=ngenerated;
			if (Aenet::calc_ANN)
			{
				Aenet::last_gen_aenet_accepted=ngenerated;
				Aenet::nstored=0;//reset the list
				aenetnew.Copy(aenetold,move.indices,thread_obj.thread_arg[nthreads - 1]);//Copy the modified energy to aenetold
				old_chitotal2 = chisquare.total2;//copying the new chi total2 to the old
			}
			else
			{
				//no ANN pot calculation in this step, store the indices
				if (!Move::nomove)
				{
					Aenet::UpdateStored();
				/*	logfile<<"aenet_stored_ind acc "<<Aenet::nstored<<endl;
					for (i=0;i<Aenet::nstored;i++)
						logfile<<Aenet::aenet_stored_ind[i]<<endl;*/ //GO
				}
			}
#endif
			calcpnew.CopyModified(calcpold,thread_obj.thread_arg[nthreads-1]);//Copy the modified partial parts 
			//------------------END OF MULTI-THREADING-------------------
			if (!Move::nomove)
			{
				if (extended_fnc)
				{
					for (i = 0; i < n_fnc; i++)//n_fnc ids the number of used fnc instances
					{
						pfnc[i].UpdatePot(0);//Copying the potential components to the 'old' arrays	
					}
				if (RunParams::potential > 0)
					old_chitotal2 = chisquare.total2;//copying the new chi total2 to the old
				}
			}
			
			
			if (Move::makeE0shift)
			{
				//preserve the value of the ek_gridind
				for (i = 0; i < nek; i++)
				{
					old_ek_gridind[i] = edata.ek_gridind[i];
					//logfile << " new E0grid  " << edata.ek_gridind[i] << "  ";
				}
				nE0shiftacc++;
				//logfile << "E0 acc " << nE0shiftacc << " chi "<< chisquare.chicomp[2];
			}
			
			if (Move::makemucorr)
			{
				for (i = 0; i < nfq; i++)
				{						
					if (ExptsData::fqIQbackgcorr[i])
					{
						old_f[i] = chisquare.f[ngr + nsq + i];
						old_ff[i] = chisquare.ff[ngr + nsq + i];
						old_fs[i] = chisquare.fs[ngr + nsq + i];
						//old_fx[i] = chisquare.fx[ngr + nsq + i];//not used presently
						//old_fx2[i] = chisquare.fx2[ngr + nsq + i];
						old_muact[i] = ExptsData::fqmuact[i];
					}
				}
				nmucorracc++;
				
			}
			else if (ExptsData::is_IQmucorr)//if there is move, preserve the last acc s, as it needed only in case for fit_index 3
			{							//and if not accepted move is followed by an accepted mucorr, the s would not be recalculated 
				for (i = 0; i < nfq; i++)
				{
					if (ExptsData::fqIQbackgcorr[i])
					{
						old_s[i] = chisquare.s[ngr + nsq + i];
						old_fs[i] = chisquare.fs[ngr + nsq + i];
						old_ss[i] = chisquare.ss[ngr + nsq + i];
					}
				}
			}
			if (Move::makefprimeshift)
			{
				nfprimeshiftacc++;
				for (i = 0; i < nfq; i++)
				{
					if (ExptsData::fqAXS[i] > 0)
					{
						old_fprimeact[i] = ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1];
						
					}
				}
				
			}
			
			
			//if mode==1 the program will run to given number of ngenerated steps, so the last step might not be acceptable, so the chi2
			//of the last acceptable steps had to be copied for displaying
			//preserving the accepted total values, might be needed for the final save in TEST MODE, if the last move(s) was/were rejected
			//or in case of rejected E0, f' shift or mucorr(here the total only changes in case of I(Q) nonlinreg alpha corr) steps to copy back the last accepted
			if (mode || (ExptsData::is_AXS || ExptsData::is_E0shift || ExptsData::is_IQmucorr))
			{
				calcdat.Copy(calcdatold);
							
				//preserving the chi2 and normalization values too
				for (i = 0; i < ChiSquared::chisize; i++)
					tempchi[i] = chisquare.chicomp[i];
				for (i = 0; i < ngr + nsq + nfq + nfg + nek; i++)
				{
					tempa[i] = chisquare.a[i];
					tempb[i] = chisquare.b[i];
					tempc[i] = chisquare.c[i];
					tempd[i] = chisquare.d[i];
					tempe[i] = chisquare.e[i];
				}
			}
			else
			{
				if (ExptsData::is_IQ > 0)
				{
					for (i = 0; i < nfq; i++)
					{
						tempa[ngr + nsq + i] = chisquare.a[ngr + nsq + i];
						tempb[ngr + nsq + i] = chisquare.b[ngr + nsq + i];
					}
				}
			}
			if (ExptsData::is_IQ > 0)
			{
				for (i = 0; i < nfq; i++)
				{
					tempalpha[i] = chisquare.alpha[i];
					templast_nlr_chi2fr[i] = chisquare.last_nlr_chi2fr[i];
					tempmaxdec_nlr_chi2fr[i] = chisquare.maxdec_nlr_chi2fr[i];
				}
			}
			
#ifdef _LOCAL_INV

			if (RunParams::loc_chi2_mode==1 && nthreads>1 && naccepted % LOC_LOAD_BALANCE_STEP==0)
				thread_obj.ResetSegment();//try to improve load balancing for the segmentation of the neighbour atoms
#endif
	
		}//end of if the move has been accepted
		else//the move is rejected
		{
			
#ifdef _WRITE_CHI2_DETAIL
			chifile << "data     \t " << ngenerated << "\t" << ntried << "\t" << naccepted << "\t" << npotaccepted << "\t" << nswapgen << "\t" <<
				nswapacc << "\t" << nE0shiftgen << "\t" << nE0shiftacc << "\t" << nmucorrgen << "\t" << nmucorracc << "\t" << nfprimeshiftgen << "\t" << nfprimeshiftacc;

			for (i = 0; i < chisquare.chisize; i++)
				chifile << "\t" << chisquare.chicomp[i];
			chifile << endl;
#endif	
			if (!Move::nomove)
			{
				// if the move has been rejected then the modified parts during this try 
				//must be restored to their original state
				move.ResetCoord(extended_fnc ? &fnc : 0);//copy back the old coordinates

				neighlist.ResetGrid(move.indices);//restore the grid to the original state
#ifdef _NO_PERIODIC
				move.ResetBinInd();//reset the central bin indices to the original state
#endif
				if (calc_neighlist)
					neighlist.ResetList();//copy back the neighbourlist values for the effected atoms
				if (CosDistrConst::nconstraints > 0)
					cosold.CopyModified(cosnew);//Copying the modified cos_hist histogram parts
#ifdef _ADVANCED_GEOM_CONST
				if (CommonNeighConst::nconstraints > 0)
					concold.CopyModified(concnew);//Copying the modified histogram parts for common neighbour constraint
				if (SecondNeighConst::nconstraints > 0)
					sncold.CopyModified(sncnew);//Copying the modified coordnumbs parts for second neighbour constraint
#endif
				avccold.CopyModified(avccnew);//Copy the modified av. coordination numbers 

				//------------------MULTI-THREADING-------------------
				iccold.CopyModified(iccnew, thread_obj.thread_arg[nthreads - 1]);//Copy the modified coordination numbers 
#ifdef _ADVANCED_GEOM_CONST
				bvsold.CopyModified(bvsnew, thread_obj.thread_arg[nthreads - 1]);//Copy the modified valences
#endif
				histold.CopyModified(histnew, thread_obj.thread_arg[nthreads - 1]);//Copying the modified histogram parts
				ppcfold.CopyModified(ppcfnew, thread_obj.thread_arg[nthreads - 1]);//Copy the modified ppcf parts 
			}
#ifdef _AENET
			if (Aenet::calc_ANN)
			{
				aenetold.Copy(aenetnew,move.indices,thread_obj.thread_arg[nthreads - 1]);//Copy the modified energy 
			}
#endif
			calcpold.CopyModified(calcpnew,thread_obj.thread_arg[nthreads-1]);//Copy the modified partial parts 


			//------------------END OF MULTI-THREADING-------------------
			if (!Move::nomove)
			{
				if (extended_fnc)
				{
					for (i = 0; i < n_fnc; i++)
						pfnc[i].UpdatePot(1);//Copying the potential components from the 'old' arrays to restore them to the values of the 
											 //last accepted move
				}
			}
			
			if (ExptsData::is_AXS || ExptsData::is_E0shift || ExptsData::is_IQmucorr)//there can be 'no atomic moves', where only certain data sets
				//are updated,copy back the last accepted moves's values, which have a algorithmic role, or can be needed in case of saving 
			{
				calcdatold.Copy(calcdat);//copy back the last accepted in case there will be a saving for the next accepted nonatomic move
				if (RunParams::potential > 0)
					chisquare.total2 = old_chitotal2;//copying the old back, as it is not recalculated
				for (i = 0; i < chisquare.chisize; i++)
					chisquare.chicomp[i] = tempchi[i];
				
				for (i = 0; i < ngr + nsq + nfq + nfg + nek; i++)
				{
					chisquare.a[i] = tempa[i];
					chisquare.b[i] = tempb[i];
					chisquare.c[i] = tempc[i];
					chisquare.d[i] = tempd[i];
					chisquare.e[i] = tempe[i];
				}
				
				if (ExptsData::is_IQ > 0)
				{
					for (i = 0; i < nfq; i++)
					{
						chisquare.alpha[i] = tempalpha[i];
						chisquare.last_nlr_chi2fr[i] = templast_nlr_chi2fr[i];
						chisquare.maxdec_nlr_chi2fr[i] = tempmaxdec_nlr_chi2fr[i];
					}
					if (ExptsData::is_IQmucorr)
					{
						for (i = 0; i < nfq; i++)//the last acceptable move's mmuact will be restored
						{
							ExptsData::fqmuact[i] = old_muact[i];
							chisquare.f[ngr + nsq + i] = old_f[i];
							chisquare.ff[ngr + nsq + i] = old_ff[i];
							//chisquare.fx[ngr + nsq + i] = old_fx[i];//not used presently for I(Q)
							//chisquare.fx2[ngr + nsq + i] = old_fx2[i];
							chisquare.s[ngr + nsq + i]=old_s[i];
							chisquare.fs[ngr + nsq + i]=old_fs[i];
							chisquare.ss[ngr + nsq + i] = old_ss[i];
						}
						edata.CalcIQmucorr();
					}
					if (Move::makefprimeshift)
					{
						for (i = 0; i < nfq; i++)
						{
							if (ExptsData::fqAXS[i] > 0)
							{
								ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1] = old_fprimeact[i];
							}
						}
						edata.CalcXrayCoeffsChange();
					}
				}
			}
			else
			{
				//normal run, only atomic moves
				if (ExptsData::is_IQ > 0)//if there is change in the I(Q)
				{
					for (i = 0; i < nfq; i++)//the last acceptable move's values will be the initial values for the non-lin regression
					{
						chisquare.a[ngr + nsq + i] = tempa[ngr + nsq + i];
						chisquare.b[ngr + nsq + i] = tempb[ngr + nsq + i];
						chisquare.alpha[i] = tempalpha[i];
						chisquare.last_nlr_chi2fr[i] = templast_nlr_chi2fr[i];
						chisquare.maxdec_nlr_chi2fr[i] = tempmaxdec_nlr_chi2fr[i];
						
					}
				}
			}
			
			if (Move::makeE0shift)
			{
				//restore the value of the ek_gridind
				for (i = 0; i < nek; i++)
					edata.ek_gridind[i] = old_ek_gridind[i];
			}
			
			if (Move::makemucorr)
			{
				for (i = 0; i < nfq; i++)//the last acceptable move's mmuact will be restored
				{
					ExptsData::fqmuact[i] = old_muact[i];
					chisquare.f[ngr + nsq + i]= old_f[i];
					chisquare.ff[ngr + nsq + i] = old_ff[i];
					//chisquare.fx[ngr + nsq + i] = old_fx[i];//not used presently for I(Q)
					//chisquare.fx2[ngr + nsq + i] = old_fx2[i];
					chisquare.s[ngr + nsq + i] = old_s[i];
					chisquare.fs[ngr + nsq + i] = old_fs[i];
					chisquare.ss[ngr + nsq + i] = old_ss[i];
				}
				edata.CalcIQmucorr();
			}
			if (Move::makefprimeshift)
			{
				for (i = 0; i < nfq; i++)
				{
					if (ExptsData::fqAXS[i] > 0)
					{
						ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1] = old_fprimeact[i];
					}
				}
				edata.CalcXrayCoeffsChange();
			}

			/*logfile.precision(15);//GO
			logfile.setf(ios::scientific, ios::floatfield);
			logfile << "backvdW ";
			for (i = 0; i < Move::npartials; i++)
				logfile << " " << FNC_POT::vdW_pot[i];
			logfile << endl;
			logfile << "backCoul ";
			for (i = 0; i < Move::npartials; i++)
				logfile << " " << FNC_POT::Coulomb_pot[i];
			logfile << endl;
			logfile.precision(6);//GO
			logfile.unsetf(ios::scientific);*/
		}//end of if the move was rejected
		cout.precision(15);
		//GO
/*#ifdef _AENET
		if (fabs(chisquare.total2 - aenetnew.E_tot) > 1e-6)
			cout << ngenerated << " " << acceptable << " e chi " << chisquare.total2 << " o " << old_chitotal2 << " etot " << aenetnew.E_tot << " olde" << aenetold.E_tot << endl;
		logfile << "e chi " << chisquare.total2 << " o " << old_chitotal2 << " etot " << aenetnew.E_tot << " olde" << aenetold.E_tot << endl;
#else
		if (fabs(chisquare.total2 - FNC_POT::vdW_pot[0]) > 1e-6)
			cout<<ngenerated<<" "<<acceptable << " e chi " << chisquare.total2 << " o " << old_chitotal2 << " etot " << FNC_POT::vdW_pot[0] << " olde" << FNC_POT::vdW_pot[0] << endl;
		logfile << "e chi " << chisquare.total2 << " o " << old_chitotal2 << " etot " << FNC_POT::vdW_pot[0] << " olde" << FNC_POT::vdW_pot[0] << endl;
#endif*/


		
#ifdef _TEST_MODE
		time2 = std::chrono::high_resolution_clock::now();
		elapsed = time2 - time1;
		dur6+=elapsed.count();
#endif
		//finally, check for time, display, save etc.....
		if (acceptable)//if the move has been accepted
		{
			//update time variables
			time(&t_current);
			duration=difftime(t_current,t_start);//that is the elapsed time since
									//the entry in the main loop
			savelapse=difftime(t_current,t_lastsave);
		
		
			//Write results
			if ((ngenerated-last_gen)>=rundata.printstep)
			{
				//time for displaying results
#ifdef _ADVANCED_GEOM_CONST
				rundata.PrintStatus(edata,chisquare,iccnew,avccnew,concnew,sncnew,bvsnew, config,1);
#else
				rundata.PrintStatus(edata, chisquare, iccnew, avccnew, config, 1);
#endif
				//update the number of moves since last display
				last_gen=ngenerated;
				last_tried = ntried;
				last_acc = naccepted;
				
			}//end of time for display
						
			

			if (savelapse>=time_save)
			{
				cout<<"\nTotal time in the loop: "<<duration<<" seconds"<<endl;
				cout <<"Time/generated move: "<<1000*duration/ngenerated<<" ms"<<endl;
				t_lastsave=t_current;//resets the time of last save
				//Save the configuration
				cout<<"**************************************************************"<<endl;
				cout<<"****** Saving the configuration and the results to disc ******"<<endl;
				cout<<"******     DO NOT INTERRUPT UNTIL SAVE IS COMPLETED     ******"<<endl;
				cout<<"**************************************************************"<<endl;
				
				if (FNC_POT::sum_out_of_range>0)
				{
					for (i=0;i<FNC_POT::nconstraints;i++)
						cout<<"The number of FNC pairs still out of FNC range for constraint "<<i+1<<" : "<<FNC_POT::nout_of_range[i]<<endl;
				}
				for (i = 0; i < ExptsData::nfq; i++)
					if (chisquare.fit_index[ngr + nsq + i] >31)//non-lin reg
						cout << "Total number of failed non-linear regression attempts for x-ray data set  " << i + 1 << " is: " << ChiSquared::tot_failed_nlr[i] << endl;
				
				config.Save(rundata.title);//saving the configuration in text format
				config.Save_binary();//saving the configurations in binary format
				rundata.SaveState(chisquare);//save the data for continuation
				histnew.Save();//saving the histogram
#ifdef _USE_LOCAL_INV
#ifndef _LOCAL_INV_NS//only save, if this was not given, may need a lot of time
				histnew.SaveLoc(lhgmfilename);
#endif
#endif
#ifdef _NO_PERIODIC
				histnew.SavePosBin();//save the central indices of the atoms
#endif
				//Save the atoms below cutoff 
				if (RunParams::moveout)//saved, if there were tooclose atoms at the beginning, and the moveout was on,
					//even if they all moved above cutoff by now
				{
					if (::debug)
						cout<<"\nSaving the atoms below cutoff"<<endl;
					Move::SaveTooClose();
				}
				//Save the grid
				neighlist.SaveGrid();

#ifdef _NEI
				neighlist.SaveNeilist();//save the neighbourlist
#endif

				if (CoordNumbConst::nconstraints > 0)
				{
					iccnew.Save();//saving the coordination constraint
					if (CoordNumbConst::is_detail)
						iccnew.SaveDetail();
				}

				if (AvCoordConst::nconstraints > 0)
					avccnew.Save();//saving the average coordination constraint

				if (CosDistrConst::nconstraints > 0)
					cosnew.Save();

#ifdef _ADVANCED_GEOM_CONST	
				if (CommonNeighConst::nconstraints > 0)
					concnew.SaveHist();
				if (SecondNeighConst::nconstraints > 0)
					sncnew.Save();
				if (BondValenceSumConst::nconstraints > 0)
				{
					bvsnew.Save();
					bvsnew.SaveBinary();
				}
#endif
				
				calcdat.SaveResult();//saving the totals
			
				if(RunParams::old_out)//RMCA format
				{
				  OpenFile(outfile,outfilename,"RMC_combined main",0);
				  ppcfnew.SaveOldOut(outfile,outfilename);//saving the ppcf-s
				  calcpnew.SaveOldOut(outfile,outfilename);//saving the partial S(Q)-s, F(Q)-s, F(g)-s and E(k)-s
				  calcdat.SaveOldOut(outfile,outfilename);//saving the totals
				  outfile.close();
				}
				
				if(RunParams::sum_ppcf) ppcfnew.AddPPCF();//updating the sum
				ppcfnew.Save();//saving the ppcf-s
				
				if (ngr>0)//saving the g(r) partials if there is g(r) data
					calcpnew.Savegr();
				
				if (nsq>0)//saving the S(Q) partials, if there is S(Q) data
					calcpnew.SaveSQ();
					
				if (nfq>0)//saving the F(Q) partials, if there is F(Q) data
					calcpnew.SaveFQ();
				
				if (nfg > 0)//saving the F(g) partials, if there is F(g) data
					calcpnew.SaveFg();

				if (nek>0)//saving the E(k) partials, if there is E(k) data
					calcpnew.SaveEK();
#ifdef _AV_MOVE
				config.SaveMovedDist();//saving the moved distance of the atoms
#endif

				//write history, if it is required
				if (history.buffsize>0)
				{
					//history is recorded
					if (nsaves==history.stepratio)
					{
						//it is time to write history
						nsaves=0;//reset the number of saves since last history recording
 						history.Buffline(duration,rundata,chisquare,edata);
						History::last_gen = ngenerated;
						History::last_tried = ntried;
						History::last_acc = naccepted;
						history_buffered=1;//History was buffered for this step
						if (history.fill>=history.buffsize)
							//it is time to flush the history buffer to the file
							history.Save(hstfile);//fill is reset to zero
					}
					else
						nsaves++;//increment the number of saves since last history buffering
				}//end of if history is recorded
				
				if (duration>run_limit && cfg_collect>0)
				{
					if (save_count%cfg_frequency==0)
					{
						save_count=1;//reset
						//collecting configurations
						cout<<"\nCollecting configurations "<<rundata.cfgnumb-cfg_collect+1<<"/"<<rundata.cfgnumb<<endl;
						string cfgcollactfn(cfgcollectfilename);
						cfgcollactfn+="_"+std::to_string(rundata.cfgnumb-cfg_collect+1)+cfgext;
						config.Save(RunParams::title,cfgcollactfn.c_str(),2);
						cfg_collect--;//decrease the number of configurations remained to be collected
					}
					else
						save_count++;//increase the counter

				}
				if (extended_fnc)
					pfnc[0].SavePotBinary();//saving the potential related parameters and the potential
#ifdef _AENET
				if (is_ANN)
					aenetnew.SavePotBinary();

				if (Aenet::write_energy)
					aenetnew.SaveEnergy();
#endif				
				cout<<"Saving completed!"<<endl;
			}//end of it is time to save results to disc
		}//end if acceptable
		
#ifdef _TEST_MODE
		time1 = std::chrono::high_resolution_clock::now();
		elapsed = time1 - time2;
		dur7+=elapsed.count();
		
#endif
	}//end of main loop
	
	//signal to threads to abandon ThreadLoop
	Threads::loop_flag=0;//set the loop indicator to zero
	//As the threads most probably already in the conditional waiting state for start_count to reach nthreads, signal them
	{
		std::unique_lock<std::mutex> guard(thread_obj.count_mutex);//lock mutex
		thread_obj.start_count++;
#ifdef _TH_ORDER
		cout << "w" << endl;
		logfile << "EM " << thread_obj.start_count << endl;
#endif
		if (thread_obj.start_count == nthreads) thread_obj.check_count.notify_all();
		else thread_obj.check_count.wait(guard);

	}

	if (mode)//program will run to given number of ngenerated steps, so the last step might be not acceptable, so the chi2
	//of the last acceptable steps had to be copied for displaying
	{
		Move::nomove = false;//set this as it is atomic move, to copy everything
		calcdatold.Copy(calcdat);//copying back the preserved total values of the last accepted move for the final save
			//needed, if the last move was rejected 
		chisquare.total=old_chitotal;//copying the old chi total to the new
		if (RunParams::potential>0)
			chisquare.total2=old_chitotal2;//copying the old chi total to the new
		for (i=0;i<chisquare.chisize;i++)
			chisquare.chicomp[i]=tempchi[i];
		for (i=0;i<ngr+nsq+nfq+ nfg+nek;i++)
		{
			chisquare.a[i]=tempa[i];
			chisquare.b[i]=tempb[i];
			chisquare.c[i]=tempc[i];
			chisquare.d[i]=tempd[i];
			chisquare.e[i] = tempe[i];
		}
		if (ExptsData::is_IQ > 0)
		{
			for (i = 0; i < nfq; i++)
			{
				chisquare.alpha[i] = tempalpha[i];
				chisquare.last_nlr_chi2fr[i] = templast_nlr_chi2fr[i];
				chisquare.maxdec_nlr_chi2fr[i] = tempmaxdec_nlr_chi2fr[i];
			}
			if (ExptsData::is_IQmucorr)
			{
				for (i = 0; i < nfq; i++)//the last acceptable move's mmuact will be restored
				{
					ExptsData::fqmuact[i] = old_muact[i];
					chisquare.f[ngr + nsq + i] = old_f[i];
					chisquare.ff[ngr + nsq + i] = old_ff[i];
					//chisquare.fx[ngr + nsq + i] = old_fx[i];//not used presently for I(Q)
					//chisquare.fx2[ngr + nsq + i] = old_fx2[i];
					chisquare.s[ngr + nsq + i] = old_s[i];
					chisquare.fs[ngr + nsq + i] = old_fs[i];
					chisquare.ss[ngr + nsq + i] = old_ss[i];
				}
				edata.CalcIQmucorr();
			}
			if (Move::makefprimeshift)
			{
				for (i = 0; i < nfq; i++)
				{
					if (ExptsData::fqAXS[i] > 0)
					{
						ExptsData::fqfprimeact[i * RunParams::ntypes + ExptsData::fqAXS[i] - 1] = old_fprimeact[i];
					}
				}
				edata.CalcXrayCoeffsChange();
			}
		}
				
		if (CosDistrConst::nconstraints)
		{
			//recalculate the cosine distribution of bond angles, in case the last move was rejected
			for (i=0;i<CosDistrConst::nconstraints;i++)
				cosnew.mod_const[i]=1;//to calculate the initial distribution
			cosnew.CalcDistribution();
		}

	}
#ifdef _AENET
	if (RunParams::aenet_step>1)
	{
		if (Aenet::last_gen_accepted==Aenet::last_gen_aenet_accepted)//there were not accepted steps after the last energy calculation
			Aenet::nstored=0;

		//because of threading the ThreadAenetFinal has to be called, if no calculation necessary, there will not be atoms
		thread_obj.SplitAenetAtoms(Aenet::nstored,Aenet::aenet_stored_ind);
		if (nthreads>1)//signal to threads to call ThreadAenetFinal
		{
			std::unique_lock<std::mutex> guard(thread_obj.aenet_mutex);//lock mutex
			thread_obj.aenetstart_flag=1;//time to proceed
			thread_obj.check_aenet_start.notify_all();
		}
		aenetnew.CalcANN(thread_obj.thread_arg[nthreads-1]);

		if (nthreads>1)
		{
			std::unique_lock<std::mutex> guard(thread_obj.count_mutex);//lock mutex
			thread_obj.complete_count2++;

			if (thread_obj.complete_count2==nthreads)//time to proceed
			{
				thread_obj.aenetstart_flag=0;//reset it to 0
				thread_obj.complete_count2=0;//reset it to 0
				thread_obj.check_count.notify_all();
			}
			else thread_obj.check_count.wait(guard);
			
		}
		if (Aenet::last_gen_accepted!=Aenet::last_gen_aenet_accepted)
		{
			chisquare.total2-=chisquare.chicomp[ChiSquared::chi_potstart];//subtract the old value
			aenetnew.E_tot=0;
			for (ithread=0;ithread<nthreads;ithread++)
			{
//				cout<<"ith "<<ithread<<endl;
				for (i=0;i<thread_obj.thread_arg[ithread].aenet_at_count;i++)
				{
//					cout<<"\ti "<<i<<" ind "<<thread_obj.thread_arg[ithread].aenet_at_index[i]<<" old "<<aenetnew.E_i[thread_obj.thread_arg[ithread].aenet_at_index[i]]<<" new "<<thread_obj.thread_arg[ithread].aenet_E_i[i]<<endl;
					aenetnew.E_i[thread_obj.thread_arg[ithread].aenet_at_index[i]]=thread_obj.thread_arg[ithread].aenet_E_i[i];
				}
			}
			//calculate the total
			for (i=0;i<Aenet::nAtoms;i++)
				aenetnew.E_tot+=aenetnew.E_i[i];
			//the chicomp should be recalculated for ANN
			
			chisquare.chicomp[ChiSquared::chi_potstart]=aenetnew.E_tot/pow(RunParams::aenet_weight,2);
			chisquare.total2+=chisquare.chicomp[ChiSquared::chi_potstart];
			if (history.buffsize>0)//there is history
				if (history_buffered)
				{
					if (History::fill>0)
						History::fill--;//already buffered, still in buffer, rewrite it
					else
						hstfile<<"	ANN potential was calculated and saved below for the last accepted configuration, which was not an ANN calculation step" <<endl;
					history_buffered=0;//add to buffer
				}
		}
	}
#endif
#ifdef _ADVANCED_GEOM_CONST
	rundata.PrintStatus(edata, chisquare, iccnew, avccnew, concnew, sncnew, bvsnew,config, 1);
#else
	rundata.PrintStatus(edata, chisquare, iccnew, avccnew, config, 1);
#endif
	cout << "\n\n=====================================================================================================================" << endl;
	cout << "==========================================                               ============================================" << endl;
	cout << "==========================================  THE SIMULATION IS COMPLETED  ============================================" << endl;
	cout << "==========================================                               ============================================" << endl;
	cout << "=====================================================================================================================" << endl;
	cout<<"\nTOTAL RUNNING TIME: "<<duration<<" seconds"<<endl;
	if (ngenerated > 0)
		cout << "TIME/GENERATED MOVES: " << 1000 * duration / ngenerated << " ms" << endl;
	cout << "TOTAL NUMBER OF MOVES:" << endl;
	cout << J20 << "GENERATED" << J13 << "TRIED" << J13 << "ACCEPTED";
	if (RunParams::swap_fraction > 0)
		cout << J20 << "SWAP GENERATED" << J20 << "SWAP ACCEPTED";
	if (RunParams::potential > 0)
		cout << J20 << "POTENTIAL ACCEPETD";
	cout << endl;
	cout << J20 << ngenerated << J13 << ntried << J13 << naccepted;
	if (RunParams::swap_fraction > 0)
		cout << J20 << nswapgen << J20 << nswapacc;
	if (RunParams::potential > 0)
		cout << J20 << npotaccepted;
	cout << endl;

	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << "I(Q) RELATED:  ";
	if (ExptsData::is_IQ == 2)
		cout << J23 << "COVERGED NON-LIN REG.";
	if (ExptsData::is_IQmucorr)
		cout << J20 << "MU CORR. GENERATED" << J20 << "MU CORR. ACCEPTED";
	if (ExptsData::is_AXS > 0)
		cout << J20 << "f'-SHIFT GENERATED" << J20 << "f'-SHIFT ACCEPTED";
	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << endl;

	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << "               ";
	if (ExptsData::is_IQ == 2)
		cout << J23 << ngoodnlr;
	if (ExptsData::is_IQmucorr)
		cout << J20 << nmucorrgen << J20 << nmucorracc;
	if (ExptsData::is_AXS > 0)
		cout << J20 << nfprimeshiftgen << J20 << nfprimeshiftacc;
	if (ExptsData::is_IQ == 2 || ExptsData::is_AXS > 0 || ExptsData::is_IQmucorr)
		cout << endl;

	if (ExptsData::is_E0shift)
	{
		cout << "EXAfS RELATED: " << J23 << "E0-SHIFT GENERATED" << J20 << "E0-SHIFT ACCEPTED" << endl;
		cout << "               " << J23 << RunParams::n_E0shiftgen << J20 << RunParams::n_E0shiftacc << endl;
	}

	//Save the configuration
	cout<<"\n****************************************************************"<<endl;
	cout<<"**** Saving the final configuration and the results to disc ****"<<endl;
	cout<<"****        DO NOT INTERRUPT UNTIL SAVE IS COMPLETED        ****"<<endl;
	cout<<"****************************************************************"<<endl;
	if (FNC_POT::sum_out_of_range>0)
	{
		for (i=0;i<FNC_POT::nconstraints;i++)
			cout<<"\nThe number of FNC pairs still out of FNC range for constraint "<<i+1<<" : "<<FNC_POT::nout_of_range[i]<<endl;
	}
	if (Move::ntooclose > 0)
		cout << "\nThe number of atoms below cutoff: " << Move::ntooclose << endl;
	
	config.Save(rundata.title);//saving the configuration
	config.Save_binary();//saving the configurations in binary format

	rundata.SaveState(chisquare);//save the data for continuation
	histnew.Save();//saving the histogram
#ifdef _NO_PERIODIC
	histnew.SavePosBin();//save the central indices of the atoms
#endif
#ifdef _USE_LOCAL_INV
#ifndef _LOCAL_INV_NS//only save, if this was not given, may need a lot of time
	histnew.SaveLoc(lhgmfilename);
#endif
#endif
	
	//Save the atoms below cutoff 
	if (RunParams::moveout)//saved, if there were tooclose atoms at the beginning, and the moveout was on,
					//even if they all moved above cutoff by now
	{
		if (::debug)
			cout<<"\nSaving the atoms below cutoff"<<endl;
		Move::SaveTooClose();
	}
	//Save the grid
	neighlist.SaveGrid();
#ifdef _NEI
	neighlist.SaveNeilist();//save the neighbourlist
#endif

	if (CoordNumbConst::nconstraints > 0)
	{
		iccnew.Save();//saving the coordination constraint
		if (CoordNumbConst::is_detail)
			iccnew.SaveDetail();
	}

	if (AvCoordConst::nconstraints>0)
		avccnew.Save();//saving the average coordination constraint
	
	if (CosDistrConst::nconstraints>0)
		cosnew.Save();
	
#ifdef _ADVANCED_GEOM_CONST	
	if (CommonNeighConst::nconstraints > 0)
		concnew.SaveHist();
	if (SecondNeighConst::nconstraints > 0)
		sncnew.Save();
	if (BondValenceSumConst::nconstraints > 0)
	{
		bvsnew.Save();
		bvsnew.SaveBinary();
	}
#endif
	calcdat.SaveResult();//saving the total

	if(RunParams::old_out)//RMCA format
	{
	  OpenFile(outfile,outfilename,"RMC_combined main",0);
	  ppcfnew.SaveOldOut(outfile,outfilename);//saving the ppcf-s
	  calcpnew.SaveOldOut(outfile,outfilename);//saving the partial S(Q)-s, F(Q)-s F(g)-s and E(k)-s
	  calcdat.SaveOldOut(outfile,outfilename);//saving the totals
	  outfile.close();
	}

	if(RunParams::sum_ppcf) ppcfnew.AddPPCF();//updating the sum
	/*logfile.precision(15);//GO
	logfile.setf(ios::scientific, ios::floatfield);
	logfile << "endvdW ";
	for (i = 0; i < Move::npartials; i++)
		logfile << " " << FNC_POT::vdW_pot[i];
	logfile << endl;
	logfile << "endCoul ";
	for (i = 0; i < Move::npartials; i++)
		logfile << " " << FNC_POT::Coulomb_pot[i];
	logfile << endl;
	logfile.precision(6);//GO
	logfile.unsetf(ios::scientific);*/

	ppcfnew.Save();//saving the ppcf-s
	if (ngr>0)//saving the g(r) partials if there is g(r) data
		calcpnew.Savegr();
	
	if (nsq>0)//saving the S(Q) partials, if there is S(Q) data
		calcpnew.SaveSQ();
		
	if (nfq>0)//saving the F(Q) partials, if there is F(Q) data
		calcpnew.SaveFQ();
	
	if (nfg > 0)//saving the F(g) partials, if there is F(g) data
		calcpnew.SaveFg();

	if (nek>0)//saving the E(k) partials, if there is E(k) data
		calcpnew.SaveEK();
	
#ifdef _AV_MOVE
	config.SaveMovedDist();//saving the moved distance of the atoms
#endif
	//Flush the history, if it is necessary
	if (history.buffsize>0)
	{
		if(!history_buffered)
			history.Buffline(duration,rundata,chisquare,edata);//Put the results of the last step to the buffer
					
		history.Save(hstfile);
		history.EndNotes(hstfile,chisquare);
#ifdef _TEST_MODE
#ifdef WRITE_TH_DURATION
		hstfile.setf(ios::right, ios::adjustfield);
		hstfile.setf(ios::scientific, ios::floatfield);
		hstfile << "Clock_per_second" << CLOCKS_PER_SEC << endl;
		hstfile << "dur0 \t" << dur0 / 1000 << endl;
		hstfile << "dur1 \t" << dur1 / 1000 << endl;
		hstfile << "dur2 \t" << dur2 / 1000 << endl;
		hstfile << "dur3 \t" << dur3 / 1000 << endl;
		hstfile << "dur3a \t" << dur3a / 1000 << endl;
		hstfile << "dur3b \t" << dur3b / 1000 << endl;
		hstfile << "dur4 \t" << dur4 / 1000 << endl;
		hstfile << "dur5 \t" << dur5 / 1000 << endl;
		hstfile << "dur6 \t" << dur6 / 1000 << endl;
		hstfile << "dur7 \t" << dur7 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "a \t" << Threads::dur_00a / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "0 \t" << Threads::dur_00 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "1 \t" << Threads::dur_01 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "2 \t" << Threads::dur_02 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "3 \t" << Threads::dur_03 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "4 \t" << Threads::dur_04 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "5 \t" << Threads::dur_05 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "6 \t" << Threads::dur_06 / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "7 \t" << Threads::dur_16[nthreads - 1] / 1000 << endl;
		hstfile << "dur_" << nthreads - 1 << "8 \t" << Threads::dur_17[nthreads - 1] / 1000 << endl;
#ifdef _ADVANCED_GEOM_CONST
		hstfile << "dur_" << nthreads - 1 << "9 \t" << Threads::dur_18[nthreads - 1] / 1000 << endl;
#endif
#ifdef _USE_LOCAL_INV
		hstfile << "dur_loc1[" << nthreads - 1 << "] \t" << Threads::dur_loc1[nthreads - 1] / 1000 << endl;
		hstfile << "dur_loc2[" << nthreads - 1 << "] \t" << Threads::dur_loc2[nthreads - 1] / 1000 << endl;
		hstfile << "dur_loc3[" << nthreads - 1 << "] \t" << Threads::dur_loc3[nthreads - 1] / 1000 << endl;
		hstfile << "dur_loc4[" << nthreads - 1 << "] \t" << Threads::dur_loc4[nthreads - 1] / 1000 << endl;
#endif
#ifdef _LOCAL_INV
		hstfile << "dur_chiloc1[" << nthreads - 1 << "] \t" << Threads::dur_chiloc1[nthreads - 1] / 1000 << endl;
		hstfile << "dur_chiloc2[" << nthreads - 1 << "] \t" << Threads::dur_chiloc2[nthreads - 1] / 1000 << endl;
		hstfile << "dur_chiloc3[" << nthreads - 1 << "] \t" << Threads::dur_chiloc3[nthreads - 1] / 1000 << endl;
		hstfile << "dur_chiloc4[" << nthreads - 1 << "] \t" << Threads::dur_chiloc4[nthreads - 1] / 1000 << endl;
#endif
#ifdef _AENET
		hstfile << "dur_ann2[" << nthreads-1 << "] \t" << Threads::dur_ann2[nthreads - 1] / 1000 << endl;
#endif
		for (ithread = 0; ithread < nthreads - 1; ithread++)
		{
			hstfile << "dur_0[" << ithread << "] \t" << Threads::dur_10[ithread] / 1000 << endl;
			hstfile << "dur_1[" << ithread << "] \t" << Threads::dur_11[ithread] / 1000 << endl;
			hstfile << "dur_2[" << ithread << "] \t" << Threads::dur_12[ithread] / 1000 << endl;
			hstfile << "dur_2a[" << ithread << "] \t" << Threads::dur_12a[ithread] / 1000 << endl;
			hstfile << "dur_2b[" << ithread << "] \t" << Threads::dur_12b[ithread] / 1000 << endl;
			hstfile << "dur_3[" << ithread << "] \t" << Threads::dur_13[ithread] / 1000 << endl;
			hstfile << "dur_3a[" << ithread << "] \t" << Threads::dur_13a[ithread] / 1000 << endl;
			hstfile << "dur_3b[" << ithread << "] \t" << Threads::dur_13b[ithread] / 1000 << endl;
			hstfile << "dur_4[" << ithread << "] \t" << Threads::dur_14[ithread] / 1000 << endl;
			hstfile << "dur_5[" << ithread << "] \t" << Threads::dur_15[ithread] / 1000 << endl;
			hstfile << "dur_5b[" << ithread << "] \t" << Threads::dur_15b[ithread] / 1000 << endl;
			hstfile << "dur_6[" << ithread << "] \t" << Threads::dur_16[ithread] / 1000 << endl;
			hstfile << "dur_7[" << ithread << "] \t" << Threads::dur_17[ithread] / 1000 << endl;
#ifdef _ADVANCED_GEOM_CONST
			hstfile << "dur_8[" << ithread << "] \t" << Threads::dur_18[ithread] / 1000 << endl;
#endif
#ifdef _USE_LOCAL_INV
			hstfile << "dur_loc1[" << ithread << "] \t" << Threads::dur_loc1[ithread] / 1000 << endl;
			hstfile << "dur_loc2[" << ithread << "] \t" << Threads::dur_loc2[ithread] / 1000 << endl;
			hstfile << "dur_loc3[" << ithread << "] \t" << Threads::dur_loc3[ithread] / 1000 << endl;
			hstfile << "dur_loc4[" << ithread << "] \t" << Threads::dur_loc4[ithread] / 1000 << endl;
#endif
#ifdef _LOCAL_INV
			hstfile << "dur_chiloc1[" << ithread << "] \t" << Threads::dur_chiloc1[ithread] / 1000 << endl;
			hstfile << "dur_chiloc2[" << ithread << "] \t" << Threads::dur_chiloc2[ithread] / 1000 << endl;
			hstfile << "dur_chiloc3[" << ithread << "] \t" << Threads::dur_chiloc3[ithread] / 1000 << endl;
			hstfile << "dur_chiloc4[" << ithread << "] \t" << Threads::dur_chiloc4[ithread] / 1000 << endl;
#endif
#ifdef _AENET
			hstfile << "dur_ann1[" << ithread << "] \t" << Threads::dur_ann1[ithread] / 1000 << endl;
			hstfile << "dur_ann2[" << ithread << "] \t" << Threads::dur_ann2[ithread] / 1000 << endl;
#endif
		}
		hstfile.unsetf(ios::scientific);
		hstfile.unsetf(ios::right);
#endif
#endif
		if (FNC_POT::Coul_warning > 0)
		{
			History::PotWarning(hstfile, FNC_POT::Coul_warning, "Coulomb");
			History::PotWarning(cout, FNC_POT::Coul_warning, "Coulomb");

		}
		if (FNC_POT::vdW_warning > 0)
		{
			History::PotWarning(hstfile, FNC_POT::vdW_warning, "vdW");
			History::PotWarning(cout, FNC_POT::vdW_warning, "vdW");

		}
		hstfile.close();
	}
	if (extended_fnc)
		pfnc[0].SavePotBinary();//saving the potential related parameters and the potential
#ifdef _AENET
	if (is_ANN)
		aenetnew.SavePotBinary();
#endif
#ifdef _AENET
	if (Aenet::write_energy)
		aenetnew.SaveEnergy();
#endif	
#ifdef _AENET
	aenet_final(&i);
#endif
	delete [] tempchi;
	delete [] tempa;
	delete [] tempb;
	delete [] tempc;
	delete [] tempd;
	delete [] tempe;
	if (tempalpha != NULL)
		delete[] tempalpha;
	if (templast_nlr_chi2fr != NULL)
		delete[] templast_nlr_chi2fr;
	if (tempmaxdec_nlr_chi2fr != NULL)
		delete[] tempmaxdec_nlr_chi2fr;
	if (old_ek_gridind != NULL)
		delete[] old_ek_gridind;
	if (old_muact != NULL)
		delete[] old_muact;
	if (old_f != NULL)
		delete[] old_f;
	if (old_ff != NULL)
		delete[] old_ff;
	/*if (old_fx != NULL)
		delete[] old_fx;
	if (old_fx2 != NULL)
		delete[] old_fx2;*/
	if (old_s != NULL)
		delete[] old_s;
	if (old_fs != NULL)
		delete[] old_fs;
	if (old_ss != NULL)
		delete[] old_ss;


#ifdef _WRITE_CHI2_DETAIL
	chifile.close();
#endif

	if (write_log)
	{
		logfile.close();
		/*logfile0.close();//GO
		logfile1.close();
		logfile2.close();*/
	}
		
	
	if (ExptsData::is_IQ == 2)
	{
		if (forced_termination)
		{
			cout << "\n*****ERROR*****" << endl;
			cout<<"The program is terminating, because the number of successive failed non - linear regressions reached " << rundata.terminate_nonlin << "!" << endl;
			for (i = 0; i < ExptsData::nfq; i++)
				if (chisquare.fit_index[ngr + nsq + i] > 31)//non-lin reg
					cout << "Total number of failed non-linear regression attempts for x-ray data set  " << i + 1 << " is: " << ChiSquared::tot_failed_nlr[i] << endl;
			cout << "\nIt has to be mentioned that after non-linear regression failed for a data set, the further chi2 calculation" << endl;
			cout << "and parameter optimization for the sunsequent data sets is not performed!" << endl;
			cout << "You can increase the failure limit by including the TERMINATE-NONLIN = XXX into the *.dat file!" << endl;
		}
		else
		{
			for (i = 0; i < ExptsData::nfq; i++)
			{
				if (chisquare.fit_index[ngr + nsq + i] > 31)//non-lin reg
				{
					if (ChiSquared::tot_failed_nlr[i] > 0)
						cout << "\nNOTE(" << ++note << "): Total number of failed non-linear regression attempts for x-ray data set  " << i + 1 << " is: " << ChiSquared::tot_failed_nlr[i] << endl;
				}
			}
		}
	}

	cout << "\n******************************************" << endl;
	cout << "**** Number of warning messages: " <<J4<< warn <<" ****"<< endl;
	cout << "**** Number of notes:            " << J4<<note <<" ****"<< endl;
	cout << "******************************************\n" << endl;

	if (RunParams::custmove == 0)
		cout << "\nEND OF ATOMIC RMC++ PROGRAM! " << endl;
	else
		cout<< "\nEND OF MOLECULAR RMC++ PROGRAM! "<<endl;
	
	if (!silent_quit)//if silent quit exiting without question, useful, if it is running in the background 
	{
		//waiting for a character to quit
		cout<<"\tQUIT (y)? ->";
		cin>>key;
	}

	return(0);

}
