//source SimpleCfg.cpp
//Last changed 27.02.2023

#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include"classes1.h"

int  SimpleCfg::ntotal;//total number of atoms
int  SimpleCfg::ntypes;// number of atom types
int  SimpleCfg::npartials;//number of partials
int  SimpleCfg::ntotvirtual;//total number of virtual sites in the configuration
int  SimpleCfg::nvirtualtypes;//number of virtual types in the configuration
int  SimpleCfg::cfg_exist;//whether the *.cfg file exist, yes by default
int  SimpleCfg::nmoves_dim=13;//dim of the nmoves array
int  SimpleCfg::max_diff_virtual;//maximum number of virtual site an atom is involved in
int *SimpleCfg::ndiff_virtual;//number of virtual sites a given atom is involved in [ntotal]
int *SimpleCfg::virtual_ind;//the RMC atom indices of the virtual sites for this atom [max_diff_virtual*ntotal]
int *SimpleCfg::virtual_source;//the RMC indices ofthe atoms building the virtual site for each virtual site [ntotvirtual*4]
int *SimpleCfg::virtual_host_type;//RMC type of the host atom for the virtual sites, which decide, which partials the potentials contribute to
int *SimpleCfg::pnatoms=0;
int *SimpleCfg::cumul=0;

int SimpleCfg::n_sepGRtypes=0;//n_sepGRtypes>=nGRtypes
int *SimpleCfg::pGRatoms=0;//pointer to number of atoms (per n_sepGRtype)
int *SimpleCfg::cumul_GR=0;//array of cumulative number of atoms split to GROMACS types (dimension n_sepGRtypes)
int *SimpleCfg::GRsep_GR_type=0;//array for each to show the GROMACS type for each separate GROMACS type segment (dim: n_sepGRtypes)

longint **SimpleCfg::nmoves;//pointer to number of generated, tried and accepeted moves

double SimpleCfg::boxedge;
double *SimpleCfg::fractions;//molar fractions
#ifdef _NO_PERIODIC
double SimpleCfg::boxedge_ori;//keep the original boxedge, it maybe needed for output
#endif

//-----------constructor----------
//BEWARE: objects created by the constructor are meaningless 
//and they should be initialised by loading data from a file	
SimpleCfg::SimpleCfg()//default constructor
{
	int i;
	SetArraysize(&positions,3*ntotal,"positions","SimpleCfg::SimpleCfg");
	SetArraysize(&finder,ntypes,"finder","SimpleCfg::SimpleCfg");//pointers to the coord of first atom of a
								// given type in the coordinates array
		
	if (pnatoms==NULL)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Number of atoms is zero. Sorry, cannot create SimpleCfg. Exiting..."<<endl;
		CleanExit();
	}
	//initialising the finder
	for(i=0;i<ntypes;i++)
		*(finder+i)=positions+(3* cumul[i]);
	
#ifdef _AV_MOVE
	SetArraysize(&start_pos,3*SimpleCfg::ntotal,"start_pos","SimpleCfg::SimpleCfg");
	SetArraysize(&moved_dist,SimpleCfg::ntotal,"moved_dist","SimpleCfg::SimpleCfg");
		
	for (i=0; i<3*SimpleCfg::ntotal; i++)
		positions[i]=0;
	for (i=0; i<SimpleCfg::ntotal; i++)
		moved_dist[i]=0;
	
		
#endif


	charge_gr_centre=NULL;//set later
}

SimpleCfg::SimpleCfg(SimpleCfg &source)//copy constructor, source to the new object
{
	int i;
	double *pos1, *pos2;

	SetArraysize(&positions,3*ntotal,"positions","SimpleCfg::SimpleCfg");
	SetArraysize(&finder,ntypes,"finder","SimpleCfg::SimpleCfg");//pointers to the coord of first atom of a
								// given type in the coordinates array
	
	pos1=source.positions;
	pos2=positions;
	for(i=0;i<ntotal;i++)
	{
		*pos2++=*pos1++;
		*pos2++=*pos1++;
		*pos2++=*pos1++;
	}

	//initialising the finder
	for(i=0;i<ntypes;i++)
		*(finder+i)=positions+(3* cumul[i]);
	
#ifdef _AV_MOVE
	SetArraysize(&start_pos,3*SimpleCfg::ntotal,"start_pos","SimpleCfg::SimpleCfg");
	SetArraysize(&moved_dist,SimpleCfg::ntotal,"moved_dist","SimpleCfg::SimpleCfg");

	
	for (i=0; i<3*SimpleCfg::ntotal; i++)
		start_pos[i]=source.start_pos[i];
	for (i=0; i<SimpleCfg::ntotal; i++)
		moved_dist[i]=source.moved_dist[i];
	
		
#endif
}

//- - - - - - - - - - -  
//Reading the static parameters (number of atoms, types) from the cfg file
//and initializing the static arrays pnatoms, cumul
int SimpleCfg::GetParamsCfg()
{
	int *typ;
	int i;
	ifstream file;

	SafeOpenTextFile(file,cfgfilename);
	if (RunParams::continuation)
	{
		if (CheckFileState(file, "SimpleCfg::GetParamsCfg") == 0)
			return (0);
	}
	else
	{
		if (CheckFileState(file, "SimpleCfg::GetParamsCfg", cfgfilename) == 0)
			return (0);
	}
	
	SkipLine(file,cfgfilename,1);//comment
	SkipLine(file,cfgfilename,1);//component
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);//moves generated, tried, accepted...
	SkipLine(file,cfgfilename,1);//configurations saved
	SkipLine(file,cfgfilename,1);
	ntotal = ReadThisLine(file, 1, 1, "ntotal", "SimpleCfg::GetParamsCfg", cfgfilename);//total number of atoms
	ntypes=ReadThisLine(file,1,1,"ntypes", "SimpleCfg::GetParamsCfg",cfgfilename);//number of atom types
	SkipLine(file,cfgfilename,1);//largest number of atoms in a molecule
	SkipLine(file,cfgfilename,1);//Euler angles
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);//box shape
	SkipLine(file,cfgfilename,1);//defining vectors are:
	boxedge=ReadThisLine(file,1,1.0,"boxedge", "SimpleCfg::GetParamsCfg",cfgfilename);//half-length of the box
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);

	//Setting the size of the arrays pnatos,cumul
	SetArraysize(&pnatoms,ntypes,"pnatoms","SimpleCfg::GetParamsCfg");
	SetArraysize(&fractions, ntypes, "fractions", "SimpleCfg::GetParamsCfg");
	SetArraysize(&cumul,ntypes+1,"cumul","SimpleCfg::GetParamsCfg");
	
	npartials=ntypes*(ntypes+1)/2;

	//loading the number of atoms for each type
	typ=pnatoms;
	for(i=0;i<ntypes;i++)
	{
		*typ=ReadThisLine(file,1,1,"pnatoms", "SimpleCfg::GetParamsCfg",cfgfilename);
		fractions[i] = (double)*typ / ntotal;
		typ++;
		SkipLine(file,cfgfilename,1);
		SkipLine(file,cfgfilename,1);
		SkipLine(file,cfgfilename,1);
	}
	//initialising the cumul array
	*cumul=0;//first element
	for(i=1;i<=ntypes;i++)
	{
		*(cumul+i)=*(cumul+i-1)+*(pnatoms+i-1);
	}

#ifdef _NO_PERIODIC
	boxedge_ori=boxedge;//keep the original boxedge, it maybe needed for output
#endif
	//Checking, whether loading was successful
	if (!CheckReadFileState(file,"SimpleCfg::GetParamsCfg",cfgfilename))
		CleanExit();//loading failed
	next_line_pos = -1;
	file.close();
	return (1);//it was succesfull
	
}

//---------------------------------------------------------------------
//Reading the static parameters (number of atoms, types) from the bcf file, if the cfg file is not present
//if there is a problem, cannot continue
//and initializing the static arrays pnatoms, cumul
void SimpleCfg::GetParamsBcf()
{
	int i,flag=0;
	int mode = 0;
	longint temp_l;//data type has to be the same, as nmoves
	ifstream file;

	CleanOpen(file,bincfgfilename,1);
	if (CheckFileState(file,"SimpleCfg::GetParamsBcf",bincfgfilename)==0)
	{
		//the *.bcf could not be opened, probably dos not exist
		cout<<"Cannot run this way, exiting"<<endl;
		CleanExit();
	}
	if (RunParams::continuation)
		mode = 2;
	temp_l=0;

	//reading the length of the longint variable
	//file.read(reinterpret_cast<char *>(&i), sizeof(i));
	
	ReadBin(file, &i, 1,"size_longint", "SimpleCfg::GetParamsBcf", bincfgfilename, mode);
	if( i!=sizeof(temp_l))
	{
		cout << "\n****ERROR****"<<endl;
		cout<<"INCONSISTENCY! The size of the typedef longint variable is "<<i<<" from the "<<bincfgfilename<<" file,"<<endl;
		cout<<"and "<<sizeof(temp_l)<<" according to this compilation of the code!"<<endl;
	
		if (i==4)
		{
			cout<<"Use a code compiled witht the I32=0 passed to make for Linux, or with _USE_INT32 switched ON in altern.h for Windows!"<<endl;
			flag=1;
		}
		else
		{
			if (i==8)
			{
				cout<<"Use a code compiled without the I32=0 passed to make for Linux, or with _USE_INT32 switched OFF in altern.h for Windows!"<<endl;
				flag=1;
			}
			else
			{
				cout<<"The length of the typedef longint is not what it is expected by the programme (4 or 8). "<<endl;
				if (RunParams::potential>0 || RunParams::fnc==4 )
					flag=1;//cannot continue
#ifdef _USE_LOCAL_INV

				flag=1;
#endif
				flag+=2;//this shows, that the value will be treated as ntotal
				ntotal=i;
				cout<<"This may be an older *.bcf file, where the first item is the number of atoms."<<endl;
			}				
		}
		
		if (flag==1 || (RunParams::continuation && flag&1))//cannot continue, if the wrong compilation is used, or if it is an old file and 
			//continuation should be used with new features as potential or local inv.
		{
			cout<<"Make sure, that the same compilation is used for the continuation of the run, as was used"<<endl;
			cout<<"for the original run!"<<endl;
			if (RunParams::continuation)
				cout<<"Cannot do exact continuation this way!"<<endl;
			CleanExit();
		
		}
		else
			cout << "\nWARNING(" << ++warn << "): There is no way to check the size of the longint, trying to continue..."<<endl;
	}
	

	//number of atoms
	if (flag<2)//only read, if longint size was just read before
		ReadBin(file, &ntotal, 1,"ntotal", "SimpleCfg::GetParamsBcf", bincfgfilename, mode);
			
	//number of atom types
	ReadBin(file, &ntypes, 1, "ntypes", "SimpleCfg::GetParamsBcf", bincfgfilename, mode);
	
	npartials=ntypes*(ntypes+1)/2;

	//Setting the size of the arrays pnatos,cumul
	SetArraysize(&pnatoms,ntypes,"pnatoms","SimpleCfg::GetParamsBcf");
	SetArraysize(&fractions, ntypes, "fractions", "SimpleCfg::GetParamsBcf");
	SetArraysize(&cumul,ntypes+1,"cumul","SimpleCfg::GetParamsBcf");
	
	//loading the number of atoms for each type

	//number of atoms/type
	for (i=0;i<ntypes;i++)
	{
		ReadBin(file, &pnatoms[i], 1, "pnatoms", "SimpleCfg::GetParamsBcf", bincfgfilename, mode);
		fractions[i] = (double)pnatoms[i] / ntotal;
		
	}
	
	//half size of the simulation cell
	ReadBin(file, &boxedge, 1, "boxedge", "SimpleCfg::GetParamsBcf", bincfgfilename, mode);

#ifdef _NO_PERIODIC
	boxedge_ori=boxedge;//keep the original boxedge, it maybe needed for output
#endif
	//initialising the cumul array
	*cumul=0;//first element
	for(i=1;i<=ntypes;i++)
	{
		*(cumul+i)=*(cumul+i-1)+*(pnatoms+i-1);
	}
		
	file.close();

	
}

//---------------Copy source to target---------------	
void SimpleCfg::Copy(SimpleCfg &source) 
{
	int i;
	double *pos1, *pos2;
	
	pos1=source.positions;
	pos2=positions;
	for(i=0;i<ntotal;i++)
	{
		*pos2++=*pos1++;
		*pos2++=*pos1++;
		*pos2++=*pos1++;
	}
}
	
//---------------save---------------	

void SimpleCfg::Save(const string title,const char *file_name,int flag) const
{
	//flag to determine how to open the file, or whether the original box edge shold be used
	//0: overwrite mode
	//2: append mode
	//3: overwrite mode, boxedge_ori should be used
	int i;
	double *pos,myboxedge;
	int *numb;
	int dumbint=0;
	ofstream file;

	myboxedge=boxedge;
#ifdef _NO_PERIODIC
	if (flag==3)
		myboxedge=boxedge_ori;//the original boxedge will be saved
#endif

	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, cfgfilename);
	
	OpenFile(file,tempfilename,"SimpleCfg::Save",flag);//open file, check, whether it was successfully opened

	file.setf(ios::fixed, ios::floatfield);
	file.setf(ios::right, ios::adjustfield);
#ifdef _TEST_MODE
	if (ntotvirtual>0)
		file<<" (Version 4 format configuration file) !file created by SimpleCfg::save !"<<endl;
	else
#endif
		file << " (Version 3 format configuration file) !file created by SimpleCfg::save !" << endl;
	file << title << endl;
	file<<endl;
	file<<endl;
	for (i=0;i<nmoves_dim;i++)
		file<<J16<<*nmoves[i];
	file << " moves gen, tried, acc, pot-acc, swap-gen, swap-acc, E0shift-gen, E0shift-acc, I(Q)mu-gen, I(Q)mu-acc, f'shift-gen, f'shift-acc, goodnonlinreg" << endl; 	
	file<<" "<<J10<<
	dumbint<<"              configurations saved\n "<<endl;
	file << " " << J10 << ntotal;
#ifdef _TEST_MODE
	if (ntotvirtual > 0)
		file <<" "<< J10 << ntotvirtual;
#endif
	file << " molecules";
#ifdef _TEST_MODE
	if (ntotvirtual > 0)
		file << ", virtual sites";
#endif
	file<<" of all types" << endl;
	file << " " << J10 << ntypes;
#ifdef _TEST_MODE
	if (ntotvirtual > 0)
		file << " " << J10 << nvirtualtypes;
#endif
	file << " types of molecules ";
#ifdef _TEST_MODE
	if (ntotvirtual > 0)
		file << "and virtual sites";
#endif
	file<< endl;
	file<<" "<<J10<<1<<" is the largest number of atoms in a molecule "<<endl;
	file<<" "<<J10<<0<<" Euler angles are provided \n "<<endl;
	file<<" "<<J10<<"F"<<" (box is cubic)	"<<endl;
	file<<"            Defining vectors are: "<<endl;

	file.precision(6);

	file<<" "<<J10<<" "<<" "<<J10<<myboxedge<<" "<<J10<<0.0<<" "<<J10<<0.0<<endl;
	file<<" "<<J10<<" "<<" "<<J10<<0.0<<" "<<J10<<myboxedge<<" "<<J10<<0.0<<endl;
	file<<" "<<J10<<" "<<" "<<J10<<0.0<<" "<<J10<<0.0<<" "<<J10<<myboxedge<<endl;
	file<<endl;

	numb=pnatoms;
	for(i=1;i<=ntypes;i++)
	{
		file<<" "<<J10<<*numb++<<" molecules of type"<<"  "<<i<<endl;
		file<<" "<<J10<<1<<" atomic sites "<<endl;
		file<<" "<<J10<<" "<<" "<<J10<<0.00<<" "<<J10<<0.00<<" "<<J10<<0.00<<"\n "<<endl;
	}
#ifdef _TEST_MODE
	for (i = 1; i <= nvirtualtypes; i++)
	{
		file << " " << J10 << *numb++ << " virtual sites of type" << "  " << i << endl;
		file << endl;
		file << endl;
		file << endl;
	}
#endif	
	file.precision(15);
	pos=positions;
	for(i=0;i<ntotal;i++)
	{
		file<<" "<<J10<<*pos++<<" ";
		file<<" "<<J10<<*pos++<<" ";
		file<<" "<<J10<<*pos++<<endl;
	};
#ifdef _TEST_MODE
	for (i = 0; i < ntotvirtual; i++)
	{
		file << " " << J10 << *pos++ << " ";
		file << " " << J10 << *pos++ << " ";
		file << " " << J10 << *pos++ << endl;
	};
#endif
	cout.unsetf(ios::fixed);
	cout.unsetf(ios::right);
	if (flag==2)
		file<<endl;//add a new line append mode
	file.close();
};

//----------------Save the SimpleCfg object to a binary file-----------------------
void SimpleCfg::Save_binary() const
{
	int i;
	longint configsize;

	ofstream file;
	OpenFile(file,bincfgfilename,"SimpleCfg::Save_binary",1);//open file, check, whether it was successfully opened

	i=sizeof(longint);
	file.write(reinterpret_cast<const char *>(&i), sizeof(i));
	file.write(reinterpret_cast<const char *>(&ntotal), sizeof(ntotal));
	file.write(reinterpret_cast<const char *>(&ntypes), sizeof(ntypes));
	file.write(reinterpret_cast<const char *>(pnatoms), (ntypes*sizeof(*pnatoms)));
	file.write(reinterpret_cast<const char *>(&boxedge), sizeof(boxedge));
	
	configsize=3*ntotal*sizeof(*positions);
	file.write(reinterpret_cast<const char *>(positions),configsize);
	//They have to be written separately, as they are references to variables in the main loop, not continuos in memory!!!
	file.write(reinterpret_cast<const char *>(&nmoves_dim), sizeof(nmoves_dim));
	for (i=0;i<nmoves_dim;i++)
			file.write(reinterpret_cast<const char *>(&(*nmoves[i])),sizeof(**nmoves));
	file.close();

}

void SimpleCfg::PutToBox()
{
	//put back the coordinates to the box if they are out

	int i;
	for (i=0;i<3*ntotal;i++)
	{
		if (positions[i]>1)
			positions[i]-=2;
		if (positions[i]<-1)
			positions[i]+=2;

	}
}

//---------------Load---------------	
//it checks, whether virtual sites are present, but ignores them at the reading
//all info about virtual sites comes from the topology
void SimpleCfg::Load()
{//loads from a file

	int *typ;
	int i;//,dumb;
	double *ppos;
	char *conv_numb=NULL,*conv_tx,*conv_ty,*conv_tz;
	ifstream file;

	//the file exist and is convereted if necessary, as it was checked before by GetParams, just in case check it again
	CleanOpen(file,cfgfilename);//no need for SaveOpenTetxtFile, i
	
	if (CheckFileState(file,"SimpleCfg::Load",cfgfilename)==0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"The "<<cfgfilename<<" file was found earlier, but now it is not available, cannot continue, exiting..."<<endl;
		CleanExit();//the *.cfg could not be opened, probably dos not exist
	}
	
	SetArraysize(&conv_tx,100,"conv_tx","SimpleCfg::Load");
	SetArraysize(&conv_ty,100,"conv_ty","SimpleCfg::Load");
	SetArraysize(&conv_tz,100,"conv_tz","SimpleCfg::Load");

	SkipLine(file,cfgfilename,1);//comment
	SkipLine(file,cfgfilename,1);//component
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);
	//as the type of the template function is determined for some reason by the default value and not the return value's type (although both InT in the template),
	//the default value has to be promoted to longint
	if (RunParams::continuation)
	{
		*nmoves[0] = ReadThisLine(file, 3,(longint)0, "ngenerated needed for continuation", "SimpleCfg::Load",cfgfilename);//ngenerated
		*nmoves[1] = ReadThisLine(file, 3,(longint)0, "ntried needed for continuation", "SimpleCfg::Load",cfgfilename);//ntried
		*nmoves[2] = ReadThisLine(file, 3, (longint)0, "naccepted needed for continuation", "SimpleCfg::Load",cfgfilename);//naccepted
		*nmoves[3] = ReadThisLine(file, 3, (longint)0, "npotaccepted needed for continuation", "SimpleCfg::Load", cfgfilename);//npotaccepted
		*nmoves[4] = ReadThisLine(file, 3, (longint)0, "nswap-generated needed for continuation", "SimpleCfg::Load", cfgfilename);//nswap-gen
		*nmoves[5] = ReadThisLine(file, 3, (longint)0, "nswap-accepted needed for continuation", "SimpleCfg::Load", cfgfilename);//nswap-acc
		*nmoves[6] = ReadThisLine(file, 3, (longint)0, "nE0shift-gen needed for continuation", "SimpleCfg::Load",cfgfilename);//nE0shift-gen
		*nmoves[7] = ReadThisLine(file, 3, (longint)0, "nE0shift-acc needed for continuation", "SimpleCfg::Load",cfgfilename);//nE0shift-acc
		*nmoves[8] = ReadThisLine(file, 3, (longint)0, "nmucorr-gen needed for continuation", "SimpleCfg::Load", cfgfilename);//nmucorr-gen
		*nmoves[9] = ReadThisLine(file, 3, (longint)0, "nmucorr-acc needed for continuation", "SimpleCfg::Load", cfgfilename);//nmucorr-acc
		*nmoves[10] = ReadThisLine(file, 3, (longint)0, "nfprimeshift-gen needed for continuation", "SimpleCfg::Load", cfgfilename);//nfrprimeshift-gen
		*nmoves[11] = ReadThisLine(file, 3, (longint)0, "nfprimeshift-acc needed for continuation", "SimpleCfg::Load", cfgfilename);//nfprimeshift-acc
		*nmoves[12] = ReadThisLine(file, 1, (longint)0, "ngoodnonlinreg needed for continuation", "SimpleCfg::Load", cfgfilename);//ngoodnonlinreg needed for nonlin regression
	}
	else
		SkipLine(file,cfgfilename,1);//number of gen, tried, acc moves...
	SkipLine(file,cfgfilename,1);//configurations saved
	SkipLine(file,cfgfilename,1);
	ntotal=ReadThisLine(file,3,0,"number of atoms", "SimpleCfg::Load",cfgfilename);//total number of atoms
	ntotvirtual = ReadThisLine(file, 4, 0, "number of virtual sites", "SimpleCfg::Load", cfgfilename);//total number of virtual sites
	SkipLine(file, cfgfilename, 1);
	ntypes=ReadThisLine(file,3,0,"ntypes", "SimpleCfg::Load",cfgfilename);//number of atom types
	nvirtualtypes = ReadThisLine(file, 4, 0, "nvirtualtypes", "SimpleCfg::Load", cfgfilename);//number of virtual types
	SkipLine(file, cfgfilename, 1);//skip the reamining of the types line
	SkipLine(file,cfgfilename,1);//largest number of atoms in a molecule
	SkipLine(file,cfgfilename,1);//Euler angles
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);//box shape
	SkipLine(file,cfgfilename,1);//defining vectors are:
	//dumb=
	ReadThisLine(file,1,1.0,"boxedge", "SimpleCfg::Load",cfgfilename);//half-length of the box is not read in case it was reset
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);
	SkipLine(file,cfgfilename,1);

	//loading the number of atom types for each type
	typ=pnatoms;
	for(i=1;i<=ntypes;i++)
	{
		IntToStr(&conv_numb,i);
		mystrcpy(conv_tx,100,"number of atoms for type ");
		*typ++=ReadThisLine(file,1,1,mystrcat(conv_tx, 100, conv_numb), "SimpleCfg::Load",cfgfilename);
		SkipLine(file,cfgfilename,1);
		SkipLine(file,cfgfilename,1);
		SkipLine(file,cfgfilename,1);
	}
	//skip the virtual types, if there is any
	for (i = 1; i <= nvirtualtypes; i++)
	{
		SkipLine(file, cfgfilename, 1);
		SkipLine(file, cfgfilename, 1);
		SkipLine(file, cfgfilename, 1);
		SkipLine(file, cfgfilename, 1);
	}
	//set it back
	ntotvirtual = 0;
	nvirtualtypes = 0;
	//now loading the positions	
	ppos=positions;
	for(i=0;i<ntotal;i++)
	{
		mystrcpy(conv_tx, 100, "x coordinate of atom ");
		mystrcpy(conv_ty, 100, "y coordinate of atom ");
		mystrcpy(conv_tz, 100, "z coordinate of atom ");
		
		IntToStr(&conv_numb,i+1);
		*ppos=ReadThisLine(file,3,1.0,mystrcat(conv_tx, 100, conv_numb), "SimpleCfg::Load",cfgfilename);
		if (*ppos < -3 || *ppos>3)
			ConfigError(*ppos,cfgfilename, conv_tx);
		ppos++;
		*ppos=ReadThisLine(file,3,1.0,mystrcat(conv_ty, 100, conv_numb), "SimpleCfg::Load",cfgfilename);
		if (*ppos < -3 || *ppos>3)
			ConfigError(*ppos,cfgfilename, conv_ty);
		ppos++;
		if (i<ntotal-1)
			*ppos=ReadThisLine(file,1,1.0,mystrcat(conv_tz, 100, conv_numb), "SimpleCfg::Load",cfgfilename);
		else
			*ppos=ReadThisLine(file,1,1.0,mystrcat(conv_tz, 100, conv_numb), "SimpleCfg::Load");
		if (*ppos < -3 || *ppos>3)
			ConfigError(*ppos,cfgfilename, conv_tz);
		ppos++;
	}
#ifdef _AV_MOVE
	for (i=0; i<3*SimpleCfg::ntotal; i++)
		start_pos[i]=positions[i];
#endif
	if (conv_numb != NULL)
		delete [] conv_numb;
	delete[] conv_tx;
	delete[] conv_ty;
	delete[] conv_tz;
	next_line_pos = -1;
	file.close();

};

//----Loads the configurations written in rmc version1 binary format----
//------------gives back 1, if it was succesfull, 0 if failed---------------
int SimpleCfg::LoadBinary(longint **moves)			
{
	ifstream file;
	char *conv_numb = nullptr,conv_text[100];
	
	CleanOpen(file,bincfgfilename,1);
	if (CheckFileState(file,"SimpleCfg::LoadBinary",bincfgfilename)==0)
	{
		if (RunParams::continuation)
		{
			cout << "\n*****ERROR*****" << endl;
			cout<<"Cannot perform run continuation without a binary coordinate file! Exiting..."<<endl;
			CleanExit();
		}
		else
		{
			cout << "\nNOTE(" << ++note << "): The configurations from the text type .cfg file will be used!"<<endl;
			return(0);
		}
	};
	
	int i,flag=0,term_status=0;
	int temp_i;//data type has to be the same, as the parameters, that is read into it!
	longint temp_l;//data type has to be the same, as nmoves
	double temp_d;//data type has to be the same as boxedge's
	char *text,*text2;
	string warning = "\tThe binary file will be ignored, the data from the .cfg file will be used!";

	SetArraysize(&text,NAME_SIZE,"text","SimpleCfg::LoadBinary");
	SetArraysize(&text2,NAME_SIZE,"text2","SimpleCfg::LoadBinary");
	if (RunParams::continuation)//the bcf file is needed, terminate,if there is problem with reading data, different message is given depending on term_status
		term_status = 2;
	if( cfg_exist == 0)//the bcf file is needed, terminate,if there is problem with reading data
		term_status = 1;
	//reading the length of the longint variable
	if (ReadBin(file, &i, 1, "size_longint", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0) { cout << warning << endl; file.close(); return(0); }

	if (i!=sizeof(temp_l))
	{
		if (RunParams::continuation)
		{
			cout << "\n****ERROR****"<<endl;
			cout<<"INCONSISTENCY! The size of the typedef longint variable is "<<i<<" from the "<<bincfgfilename<<" file,"<<endl;
			cout<<"and "<<sizeof(temp_l)<<" according to this compilation of the code!"<<endl;
		}
		if (i==4)
		{
			if (RunParams::continuation)
			{
				cout<<"Use a code compiled with the I32=0 passed to make for Linux, or with _USE_INT32 switched ON in altern.h for Windows!"<<endl;
				flag=1;
			}
		}
		else
		{
			if (i==8)
			{
				if (RunParams::continuation)
				{
					cout<<"Use a code compiled without the I32=0 passed to make for Linux, or with _USE_INT32 switched OFF in altern.h for Windows!"<<endl;
					flag=1;
				}
			}
			else
			{
#ifndef _USE_LOCAL_INV
				if (!RunParams::continuation && !(RunParams::potential > 0 || RunParams::fnc == 4))
					cout << "\nWARNING(" << ++warn << "): " << endl;
#endif
				cout<<"The length of the typedef longint is not what it is expected by the programme (4 or 8). "<<endl;
				if (RunParams::potential>0 || RunParams::fnc==4 )
					flag=1;//cannot continue
#ifdef _USE_LOCAL_INV

				flag=1;
#endif
				if (i==ntotal)
				{
					flag+=2;//this shows, that the value will be treated as ntotal
					cout<<"This must be an older *.bcf file, where the first item is the number of atoms."<<endl;
				}
			}				
		}
		
		if (flag==1 || (RunParams::continuation && flag&1))//cannot continue, if the wrong compilation is used, or if it is an old file and con-
			//tinuation should be used with new features as potential or local inv.
		{
			cout<<"Make sure, that the same compilation is used for the continuation of the run, as was used"<<endl;
			cout<<"for the original run!"<<endl;
			cout<<"Cannot do exact continuation this way!"<<endl;
			CleanExit();
		
		}
		else
			cout << "\nWARNING(" << ++warn << "): There is no way to check the size of the longint, trying to continue..."<<endl;
	}
	//Reading and checking the static parameters
	//number of atoms
	if (flag<2)//only read, if longint size was just read before
	{
		if (ReadBin(file, &temp_i, 1, "ntotal", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0) { cout << warning << endl; file.close(); return(0); }
		
		if (temp_i!=ntotal)
		{
			cout << "\n****ERROR****" << endl;
			cout<<"INCONSISTENCY! The number of atoms is "<<ntotal<<" from the .cfg file, and "<<endl;
			cout<<temp_i<<" from the binary .bcf file!"<<endl;
			cout<<"The binary file will be ignored, the data from the .cfg file will be used!"<<endl;
			return(0);
		}
	}
	//number of atom types
	if (ReadBin(file, &temp_i, 1, "ntypes", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0) { cout << warning << endl; file.close(); return(0); }
	if (temp_i!=ntypes)
	{
		cout << "\n****ERROR****" << endl;
		cout<<"INCONSISTENCY! The number of atom types is "<<ntypes<<" from the .cfg file, and "<<endl;
		cout<<temp_i<<" from the binary .bcf file!"<<endl;
		cout<<"The binary file will be ignored, the data from the .cfg file will be used!"<<endl;
		return(0);
	}
	//number of atoms/type
	for (i=0;i<ntypes;i++)
	{
		if (ReadBin(file, &temp_i, 1, "pnatoms", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0) { cout << warning << endl; file.close(); return(0); }
		if (temp_i!=pnatoms[i])
		{
			cout << "\n****ERROR****" << endl;
			cout<<"INCONSISTENCY! The number of atoms for type "<<i+1<<" is "<<pnatoms[i]<<" from the .cfg file, and "<<endl;
			cout<<temp_i<<" from the binary .bcf file! (Starting with 1)"<<endl;
			cout<<"The binary file will be ignored, the data from the .cfg file will be used!"<<endl;
			return(0);
		}
	}

	//half size of the simulation cell
	if (ReadBin(file, &temp_d, 1, "boxedge", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0) { cout << warning << endl; file.close(); return(0); }
	
	//Loading of the coordinates
	if (ReadBin(file, positions, 3 * ntotal, "positions", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0) { cout << warning << endl; file.close();	return(0); }
	
	double *ppos = positions;
	for (i = 0; i < ntotal;i++)
	{
		IntToStr(&conv_numb, i + 1);
		for (int j = 0; j < 3; j++)
		{
			if (*ppos <= -3 || *ppos > 3)
			{
				switch (j)
				{
					case 0: 
						mystrcpy(conv_text, 100, "x coordinate of atom ");
						mystrcat(conv_text, 100, conv_numb);
						break;
					case 1: 
						mystrcpy(conv_text, 100, "y coordinate of atom ");
						mystrcat(conv_text, 100, conv_numb);
						break;
					default: 
						mystrcpy(conv_text, 100, "z coordinate of atom ");
						mystrcat(conv_text, 100, conv_numb);
				}
				ConfigError(*ppos, bincfgfilename, conv_text);

				ppos++;
			}
		}
	}
	if (conv_numb!=nullptr)
		delete[] conv_numb;

#ifdef _AV_MOVE
	for (i=0; i<3*SimpleCfg::ntotal; i++)
		start_pos[i]=positions[i];
#endif
	if (RunParams::continuation)//the run should continue, number of steps has to be read
	{
		SetArraysize(moves, nmoves_dim, "moves", "SimpleCfg::LoadBinary");
		flag = 0;//this will be used, if nmoves_dim is not a match
		if (ReadBin(file, &temp_i, 1, "nmoves_dim", "SimpleCfg::LoadBinary", bincfgfilename, term_status) == 0)
		{
			cout << warning << endl; 
			file.close();
			return(0);
		}
		if (temp_i != nmoves_dim)
		{
			cout << "\n****ERROR****" << endl;
			cout << "INCONSISTENCY! The dimension of the nmoves array is different in the program (" << nmoves_dim;
			cout<< ") and in the .bcf file ("<<temp_i<<")!" << endl;
			cout << "The run cannot continue, start it again without the continuation option! Exiting..." << endl;
			CleanExit();
		}
		for (i=0;i<nmoves_dim;i++)
		{
			file.read(reinterpret_cast<char *>((*moves)+i),sizeof(**nmoves));//this may not been written for old bcfs
			//check, whether it was given
			mystrcpy(text,NAME_SIZE,"SimpleCfg::LoadBinary number of ");
			switch (i)
			{
				case 0:
				{
					mystrcpy(text2,NAME_SIZE,"generated moves");
					break;
				}
				case 1:
				{
					mystrcpy(text2, NAME_SIZE, "tried moves");
					break;
				}
				case 2:
				{
					mystrcpy(text2, NAME_SIZE, "accepted moves");
					break;
				}
				case 3:
				{
					mystrcpy(text2, NAME_SIZE, "potential-accepted moves");
					break;
				}
				case 4:
				{
					mystrcpy(text2, NAME_SIZE, "swap-generated moves");
					break;
				}
				case 5:
				{
					mystrcpy(text2, NAME_SIZE, "swap-accepted moves");
					break;
				}
				case 6:
				{
					mystrcpy(text2, NAME_SIZE, "E0shift-generated moves");
				}
				case 7:
				{
					mystrcpy(text2, NAME_SIZE, "E0shift-accepted moves");
				}
				case 8:
				{
					mystrcpy(text2, NAME_SIZE, "I(Q)mucorr-generated moves");
				}
				case 9:
				{
					mystrcpy(text2, NAME_SIZE, "I(Q)mucorr-accepted moves");
				}
				case 10:
				{
					mystrcpy(text2, NAME_SIZE, "f'shift-generated moves");
				}
				case 11:
				{
					mystrcpy(text2, NAME_SIZE, "f'shift-accepted moves");
				}
				case 12:
				{
					mystrcpy(text2, NAME_SIZE, "good linear-regression moves");
				}
			}
			mystrcat(text,100,text2);
			
			if (!CheckReadFileState(file,text,bincfgfilename))
			{
				cout << "\n*****ERROR*****" << endl;
				cout<<"The  run cannot continue, start it again without the continuation option! Exiting..."<<endl;
				CleanExit();
			}
		
			if (cfg_exist)
			{//check consistency with the cfg file
				if ((*moves)[i] != *nmoves[i])
				{
					cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number " << text2 << " is " << *nmoves[i] << " from the .cfg file," << endl;
					cout << "\tand " << (*moves)[i] << " from the binary .bcf file!" << endl;

					//cout<<"The  run cannot continue, start it again without the continuation option! Exiting..."<<endl;
					//CleanExit();

					flag = 1;
				}
			}
			else
				*nmoves[i]=(*moves)[i];//could not be read in GetParams_binary, as the arrays could not be initialised in that time
		}
		if (!flag)//no mismatch in the numer of moves
		{
			delete[] * moves;
			*moves = nullptr;
		}
	}

	file.close();
	
	return(1);//the load succesfull, there is a mismatch in the number of moves, if the moves array is not NULL
}
//----------//compare the coordinates of the caller object with object conf------
int SimpleCfg::Comp(const SimpleCfg &totconf) const
{
	int j;
	double *coord1,*coord2,diff;

	
	coord1=positions;
	coord2=totconf.positions;
	for (j=0;j<3*ntotal;j++)
	{
		diff=*coord1++ - *coord2++;
		if (fabs(diff)>TOLERANCE)
		{
			cout << "\nWARNING(" << ++warn << "): The "<<(j % 3)+1<<" coordinate of the "<<j/3+1<<" atom differs "<<diff<<" in the text and binary files!"<<endl;
			cout<<"\tThat is more, than the tolerance "<<TOLERANCE<<endl; 
			cout<<"\tThe two files might not come from the same run!"<<endl;
			cout<<"\tThe binary file will be ignored!"<<endl;
			cout.precision(15);
			cout<<"\t"<<*(coord1-1)<<"  "<<*(coord2-1)<<endl;
			cout.precision(6);
			return(0);
		}
	}
	return(1);
};

#ifdef _AV_MOVE
//calculate the moved distance and the average of it for the atoms
double SimpleCfg::CalcAvMove()
{
	int i,j;
	double sum,d,*p1,*p2, *moved,av_d=0;
	p1=start_pos;
	p2=positions;
	moved=moved_dist;

	for (i=0; i<ntotal; i++)
	{
		sum=0;
		for (j=0; j<3; j++)
		{
			d=*p1 - *p2;
			//taking into account the minimum image
			if (d>=1)
				d-=2;
			else
			{
				if (d<=-1)
					d+=2;
			}
			
			sum+=d * d;
			p1++;
			p2++;
		}
		*moved=sqrt(sum);
		if (*moved>TOLERANCE)	
			av_d+=*moved++;
		else
			*moved++=0;//only rounding error
	}
	return av_d/ntotal;
}
#endif



//setting up the arrays connected to the GROMACS types
void SimpleCfg::SetGRtype()
{
	int i,imoltype,iat,imol,ind1,ind2,iv ;
	int *temp_GRtype;

	SetArraysize(&temp_GRtype,ntotal+ntotvirtual,"temp_GRtype","SimpleCfg::SetGRtype");
	for (i=0;i<ntotal+ntotvirtual;i++)
		temp_GRtype[i]=-1;//not assigned
	iv = 0;
	for (imoltype=0;imoltype<Topology::nmoltype;imoltype++)
	{
		for (iat=0;iat<Topology::natoms_per_type[imoltype];iat++)
		{
			//iat is the GROMACS index of the atom in the molecule, it start with 0
			//first_index[cumul_atom[imoltype]+iat] is the RMC_index of the first occurence of this atomtype, if this molecule
			//	type would be the first in the system
			//RMCindex_offset[imoltype] is the offset of the RMC index to take into account that might not be the first moleculetype
			//imol*delta_index[cumul_atom[imoltype]+iat] is the offset of this occurence's index relative to the first
			
			//RMC_index of the atom in the first instance
			if (Topology::first_index[Topology::cumul_atoms[imoltype] + iat] == 0)//virtual site
			{
				ind1 = cumul[ntypes+ iv];
				iv++;
			}
			else
				ind1=Topology::first_index[Topology::cumul_atoms[imoltype]+iat]-1+Topology::RMCindex_offset[imoltype];
			for (imol=0;imol<Topology::nmol_per_type[imoltype];imol++)
			{
				ind2=ind1+imol*Topology::delta_index[Topology::cumul_atoms[imoltype]+iat];
				temp_GRtype[ind2]=Topology::GROMACS_type[Topology::cumul_atoms[imoltype]+iat];
			}
		}
	}
	//Check,if each atom is accounted for
	for (i=0;i<ntotal+ ntotvirtual;i++)
	{
		if (temp_GRtype[i]==-1)
		{
			//This atom was not assigned
			cout<<"RMC atom "<<i+1<<" was not assigned a GROMACS type, something is wrong with the topology-configuration"<<endl;
			cout<<"consistency! Reset the topology RMC_indices and try again! Exiting..."<<endl;
			CleanExit();
		}
	}
	//Finding the consecutive GROMACS type segments
	n_sepGRtypes=1;
	SetArraysize(&pGRatoms,n_sepGRtypes,"pGRatoms","SimpleCfg::SetGRtype");//number of atoms for each separate GROMACS type
	SetArraysize(&GRsep_GR_type,n_sepGRtypes,"GRsep_GR_type","SimpleCfg::SetGRtype");//showing the real GROMACS type for each separate GROMACS type segment
	SetArraysize(&cumul_GR,n_sepGRtypes+1,"cumul_GR","SimpleCfg::SetGRtype");//cumulative number of atoms for each separate GROMACS type
	pGRatoms[0]=1;
	cumul_GR[0]=0;
	GRsep_GR_type[0]=temp_GRtype[0];
	for (i=1;i<ntotal+ntotvirtual;i++)
	{
		if (temp_GRtype[i]==temp_GRtype[i-1])//this has the same GROMACS type, as the previous
			pGRatoms[n_sepGRtypes-1]++;//increase the numer of atoms for this type
		else
		{
			cumul_GR[n_sepGRtypes]=cumul_GR[n_sepGRtypes-1]+pGRatoms[n_sepGRtypes-1];
			ResizeArray(&n_sepGRtypes,n_sepGRtypes+1,&pGRatoms,"pGRatoms","SimpleCfg::SetGRtype");
			n_sepGRtypes--;//put it back, as cumul is 1 longer
			ResizeArray(&n_sepGRtypes,n_sepGRtypes+1,&GRsep_GR_type,"GRsep_GR_type","SimpleCfg::SetGRtype");
			ResizeArray(&n_sepGRtypes,n_sepGRtypes+1,&cumul_GR,"cumul_GR","SimpleCfg::SetGRtype");
			n_sepGRtypes--;//put it back, as cumul is 1 longer
			pGRatoms[n_sepGRtypes-1]=1;
			GRsep_GR_type[n_sepGRtypes-1]=temp_GRtype[i];

		}
	}
	cumul_GR[n_sepGRtypes]=cumul_GR[n_sepGRtypes-1]+pGRatoms[n_sepGRtypes-1];

	delete [] temp_GRtype;
};

//setting up the arrays connected to the virtual sites
//one atom can be involved in building more than one virtual sites
void SimpleCfg::SetVirtual()
{
	int i, ii,j,k, imoltype, vi, iat, imol;
	int ind1, ind2,delta;
	int found;
	
	if (Topology::n_gr_virtuals == 0)
		return;
	
	//set the SimpleCfg values connected to the virtual sites
	ntotvirtual = 0;//total number of virtual atoms
	nvirtualtypes = Topology::n_gr_virtuals;//each virtual site given in topology [virtual_stesX ] forms a virtual type
	i = ntypes;
	ResizeArray(&i, ntypes + nvirtualtypes, &pnatoms, "pnatoms", "SimpleCfg::SetVirtual");
	i = ntypes + 1;
	ResizeArray(&i, ntypes + nvirtualtypes+1, &cumul, "cumul", "SimpleCfg::SetVirtual");
	SetArraysize(&ndiff_virtual, ntotal, "ndiff_virtual", "SimpleCfg::SetVirtual");//number of virtual sites an atom is involved in
	SetArraysize(&virtual_host_type, nvirtualtypes, "virtual_RMC_type", "SimpleCfg::SetVirtual");
	i = ntypes;//RMC index of the virtual type in the configuration
	k = 0;//consecutive virtual type index
	for (imoltype = 0; imoltype < Topology::nmoltype; imoltype++)
	{
		ntotvirtual += Topology::nvirtuals_per_type[imoltype] * Topology::nmol_per_type[imoltype];
		for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)
		{
			pnatoms[i] = Topology::nmol_per_type[imoltype];
			cumul[i + 1] = cumul[i] + Topology::nmol_per_type[imoltype];
			//determine the RMC type for each virtual site to be able to add the potentials to the partials determined by it
			ind2 = Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].host_index;//gromacs index of the host in its own mol
			ind1 = Topology::first_index[Topology::cumul_atoms[imoltype]+ind2-1] - 1 + Topology::RMCindex_offset[imoltype];
			virtual_host_type[k] = 0;
			while (ind1 >= cumul[virtual_host_type[k] + 1])
				virtual_host_type[k]++;
			k++;
			i++;
		}
	}
		
	//this will hold the RMC indices of the atoms building the virtual site for each virtual site
	SetArraysize(&virtual_source, ntotvirtual*4, "virtual_source", "SimpleCfg::SetVirtual");
	for (i = 0; i < ntotvirtual * 4; i++)
		virtual_source[i] = -1;
	//first only detemnine the maximum how many virtual site an atom can be involved in
	max_diff_virtual = 0;
	i = 0;
	for (imoltype = 0; imoltype < Topology::nmoltype; imoltype++)
	{
		for (iat = 0; iat < Topology::natoms_per_type[imoltype]; iat++)
		{
			if (Topology::first_index[i] != 0)
			{
				found = 0;
				for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
				{
					for (j = 0; j < Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].tnumb; j++)//going through the indices of the atoms making the site 
					{
						if (iat + 1 == Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].indices[j])//the iat atom is among the ones boulding this virtual site
						{
							found++;
							break;//j cycle
						}
					}
				}//end of vi cycle
				if (found > max_diff_virtual)
					max_diff_virtual = found;
				//set the number of virtual sites for all the instances of this GROMACS atom
				//RMC_index of the atom building the site in the first instance, and the increment for the following ones
				/*ind1 = Topology::first_index[i] - 1 + Topology::RMCindex_offset[imoltype];
				delta = Topology::delta_index[i];
				for (imol = 0; imol < Topology::nmol_per_type[imoltype]; imol++)
					ndiff_virtual[ind1 + imol * delta] = found;*/
			}//not a virtual site
			i++;
		}//end of iat cycle
	}//end of imoltype cycle
		
		
	//now set the arrays and get the virtual_site's RMC indices for the atoms, and the indices of the atoms building the virtual sites for each virtual site
	SetArraysize(&virtual_ind, ntotal * max_diff_virtual, "virtual_ind", "SimpleCfg::SetVirtual");
	for (i = 0; i < ntotal * max_diff_virtual; i++)
		virtual_ind[i] = -1;
	for (i = 0; i < ntotal; i++)
		ndiff_virtual[i] = 0;
	
	for (imoltype = 0; imoltype < Topology::nmoltype; imoltype++)
	{
		ii = Topology::cumul_atoms[imoltype];
		found = 0;
		for (vi = 0; vi < Topology::nvirtuals_per_type[imoltype]; vi++)//go through the virtual sites of this molecule type
		{
			for (j = 0; j < Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].tnumb; j++)//going through the indices of the atoms making the site 
			{
				i = ii + Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].indices[j] - 1;//consecutive GROMACS index of the atom involved in the virtual site 
				ind1 = Topology::first_index[i] - 1 + Topology::RMCindex_offset[imoltype]; //RMC index
				delta = Topology::delta_index[i];
				for (imol = 0; imol < Topology::nmol_per_type[imoltype]; imol++)
				{
					ind2 = ind1 + imol * delta;
					if (ind2 >=SimpleCfg::ntotal)
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "In SimpleCfg::SetVirtual the RMC index of the " << imoltype + 1 << ". molecule type's " << Topology::virtual_site[Topology::cumul_virtuals[imoltype] + vi].indices[j] << ". atoms is " << ind2+1 << " in the ";
						cout<< imol + 1 << ". molecule, which is larger than the number of atoms in the molecule (" << ntotal << ")" << endl;
						cout << "Something is wrong with the RMC index assignment in the topology, correct it!" << endl;
						cout << "Cannot continue, exiting..." << endl;
						CleanExit();
					}
					virtual_ind[ind2 * max_diff_virtual + ndiff_virtual[ind2]] = SimpleCfg::cumul[SimpleCfg::ntypes + Topology::cumul_virtuals[imoltype] + vi] + imol;
					ndiff_virtual[ind2]++;
				}//end of imol cycle
			}//end of j cycle
		}//end of vi cycle			
	}//end of imoltype cycle
	
	//setting the virtual source array with the RMC indices of the atoms building each virtual site
	for (imoltype = 0; imoltype < Topology::nmoltype; imoltype++)
	{
		for (imol = 0; imol < Topology::nmol_per_type[imoltype]; imol++)
		{
			for (j = 0; j < Topology::nvirtuals_per_type[imoltype]; j++)
			{
				vi = Topology::cumul_virtuals[imoltype]+j;//index of the virtual site in the virtual_site array
				i = ntypes + vi;//this will be the RMC type of the virtual site
				ind1 = cumul[i]-cumul[ntypes]+imol;//index of the virtual site in the virtual_source array
				
				for (k = 0; k < Topology::virtual_site[vi].tnumb; k++)
				{	
					//the GROMACS index in its own molecule of the atom involved in building this site is Topology::virtual_site[vi].indices[k]
					//the GRomas index of it in the configuration is Topology::cumul_atoms[imoltype] + Topology::virtual_site[vi].indices[k]
					//ind2 is the RMC index of the building atom in the first instance
					ind2=Topology::first_index[Topology::cumul_atoms[imoltype] + Topology::virtual_site[vi].indices[k]-1] - 1 + Topology::RMCindex_offset[imoltype];
					delta= Topology::delta_index[Topology::cumul_atoms[imoltype] + Topology::virtual_site[vi].indices[k] - 1];
					virtual_source[4*ind1 + k] = ind2+delta*imol;
					
				}
				
				
			}
		}
	}
		
};

//the virtual sites has to be handled as the atoms if there is any
void SimpleCfg::SetChargeGroupCentre(FNC_POT &pfnc)
{
	int i,j,icoord;
	int nat;
	int first_ind_group;
	int *numb,*poffset,*pind;
	double pos;

	//Assigning the atoms to the atom_for_charge_gr array and setting the pointer
	SetArraysize(&numb,Topology::nRMC_charge_centre,"numb","SimpleCfg::SetChargeGroupCentre");
	for (i=0;i<Topology::nRMC_charge_centre;i++)
		numb[i]=0;
	for (i=0;i<ntotal+ntotvirtual;i++)
		numb[pfnc.charge_centre[i]]++;//counting the atoms for each charge_centre

	pfnc.atom_for_charge_gr_finder[0]=pfnc.atom_for_charge_gr;
	for (i=1;i<Topology::nRMC_charge_centre;i++)//setting the finder
		pfnc.atom_for_charge_gr_finder[i]=pfnc.atom_for_charge_gr_finder[i-1]+numb[i-1];
	for (i=0;i<ntotal + ntotvirtual;i++)
	{
		poffset=pfnc.atom_for_charge_gr_finder[pfnc.charge_centre[i]];//position for the first atom for this charge centre
		while (*poffset>-1)
			poffset++;//skipping the position already filled for this group
		*poffset=i;
	}
	
	//some quick checking, whether each was assigned, though probably that not all were assigned cannot happen
	for (i=0;i<ntotal + ntotvirtual;i++)
	{
		if (pfnc.atom_for_charge_gr[i]==-1)
		{
			cout << "\n****ERROR****"<<endl;
			cout << "The assigment of the atoms ";
			if (ntotvirtual > 0)
				cout << "and charge groups ";
			cout<<"to the charge groups (array atom_for_charge_gr of the FNC_POT class)"<<endl;
			cout<<"has failed, no atoms was assigned for array element "<<i<<" (index start with 0)!"<<endl;
			cout<<"Maybe something is wrong with the charge groups in the Topology, try again!"<<endl;
			CleanExit();
		}
	}

	
	//this will be the position of the charge group centres
	//there is no clear indications, whether in case of virtual sites, the virtual site coordinates should be used to contribute
	//when the position of the charge group centre is calculated. Making a lot of trials with grompp, as it gives the maximum van der Waals and
	//Coulomb largest charge centre radii it was conluded, that the virtual sites are included in the calculations, and the 
	//same charge centre is used for both van der Waals and Coulomb, but the radii can be different, as those atoms, which have zero 
	//LJ parameters, does not considered during the LJ radii calculation, and maybe those with no charge during the Coulomb radii
	//calculation. The radii does not matter in RMC... 
	SetArraysize(&charge_gr_centre,Topology::nRMC_charge_centre*3,"charge_gr_centre","SimpleCfg::SetChargeGroupCentre");

	for (i=0;i<Topology::nRMC_charge_centre*3;i++)
		charge_gr_centre[i]=0;
	
	
	for (i=0;i<Topology::nRMC_charge_centre;i++)
	{
		pind=pfnc.atom_for_charge_gr_finder[i];//pointer to the first atom of this charge centre
		//the charge group will be centred on the first atom of the charge group, which will be shifted to 0,0,0	
		first_ind_group=*pind;//first atom of this charge centre
		nat = Topology::natoms_per_charge_group[pfnc.charge_group[first_ind_group]];//number of atoms in the charge group
		pind++;//points to the second atom of the charge group centre
		for (j=1;j<nat;j++)//going through the atoms of this RMC charge group
		{
			//calculating the centre of the charge group
			for (icoord=0;icoord<3;icoord++)//sum
			{
				pos=positions[3* *pind +icoord]-positions[3*first_ind_group+icoord];
				if (pos<-1)
					pos+=2;
				else
				{
					if (pos>1)
						pos-=2;
				}

				charge_gr_centre[3*i+icoord]+=pos;
			}
			pind++;
		}
		for (icoord=0;icoord<3;icoord++)
		{
			charge_gr_centre[3*i+icoord]/=nat;//average
			charge_gr_centre[3*i+icoord]+=positions[3*first_ind_group+icoord];//translate back
		}
		GetMinImage(charge_gr_centre + 3 * i);
		
				
	}
};

//creating the coordinates of the virtual sites with RMC index ind
void SimpleCfg::CalcVirtualCoord(int ind)//ind is the index of the  virtual site in the configuration
{
	int ivtype,icoord;
	double a, b,c, d,theta,dist, s_ij_jk;
	double pos_i[3], pos_j[3], pos_k[3], pos_l[3];
	double r_ja[3],r_jb[3],r_m[3],d_r_m,r_comp[3];

	ivtype= 0;
	while (ind >= cumul[(ivtype) + 1])
		(ivtype)++;
	//now ivtype is the RMC type of the virtual atom
	ivtype -= ntypes;//this is the position in the virtual sites array
	//broken molecules has to be made whole befoe calculation, the molecules will be translated that the first source is at zero
	for (icoord = 0; icoord < 3; icoord++)
	{
		pos_i[icoord] = *(positions + 3 * virtual_source[(ind - ntotal) * 4] + icoord);//positions of  atom_i
		pos_j[icoord] = *(positions + 3 * virtual_source[(ind - ntotal) * 4 + 1] + icoord) - pos_i[icoord];//position of atom_j -pos_i
		if (pos_j[icoord] < -1)
			pos_j[icoord] += 2;
		else
		{
			if (pos_j[icoord] > 1)
				pos_j[icoord] -= 2;
		}
	}
	//during the calulations pos_i=0
	switch (Topology::virtual_site[ivtype].tnumb)
	{
	case 2: //"2"

		//calculating the coordinates of the virtula site and translating back to its palce
		for (icoord = 0; icoord < 3; icoord++)
			positions[3 * ind + icoord] = *(pos_i+icoord)+Topology::virtual_site[ivtype].params[0] * (*(pos_j + icoord));
		break;
	case 3:
		for (icoord = 0; icoord < 3; icoord++)
		{
			pos_k[icoord] = *(positions + 3 * virtual_source[(ind - ntotal) * 4 + 2] + icoord) - pos_i[icoord];//position of atom_k -pos_i
			if (pos_k[icoord] < -1)
				pos_k[icoord] += 2;
			else
			{
				if (pos_k[icoord] > 1)
					pos_k[icoord] -= 2;
			}
		}
		switch (Topology::virtual_site[ivtype].type)
		{
		case 1:   //"3"
			a = Topology::virtual_site[ivtype].params[0];//parameter a
			b = Topology::virtual_site[ivtype].params[1];//parameter b
			for (icoord = 0; icoord < 3; icoord++)
				positions[3 * ind + icoord] = *(pos_i + icoord) + a * *(pos_j + icoord) + b * *(pos_k + icoord);
			break;
		case 2:  //"3fd"
			a = Topology::virtual_site[ivtype].params[0];//parameter a
			b = Topology::virtual_site[ivtype].params[1];//parameter b
			dist = 0;
			for (icoord = 0; icoord < 3; icoord++)
				dist += pow(pos_j[icoord] + a * (pos_k[icoord] - pos_j[icoord]), 2);
			dist = sqrt(dist);//|r_ij+a*r_jk|
			for (icoord = 0; icoord < 3; icoord++)
				positions[3 * ind+ icoord]= *(pos_i + icoord) + b / dist * (*(pos_j + icoord) + a * (*(pos_k + icoord) - *(pos_j + icoord)));
			break;
		case 3:  //"3fad"
			theta = Topology::virtual_site[ivtype].params[0] * PI / 180;//parameter theta in radian
			d = Topology::virtual_site[ivtype].params[1];//parameter d
			a = pow(*pos_j, 2) + pow(*(pos_j + 1), 2) + pow(*(pos_j + 2), 2);//|r_ij|2
			s_ij_jk = (*pos_j ) * (*pos_k - *pos_j) + (*(pos_j + 1)) * (*(pos_k + 1) - *(pos_j + 1)) + (*(pos_j + 2) ) * (*(pos_k + 2) - *(pos_j + 2));//r_ij*r_jk scalar
			b = 0;
			for (icoord = 0; icoord < 3; icoord++)//components of right_angle
			{
				r_comp[icoord] = (*(pos_k + icoord) - *(pos_j + icoord)) - s_ij_jk / a * *(pos_j + icoord);
				b += pow(r_comp[icoord], 2);
			}
			b = sqrt(b);//|right_angle|
			for (icoord = 0; icoord < 3; icoord++)
				positions[3 * ind+icoord] = *(pos_i + icoord) + d * cos(theta) / sqrt(a) * (*(pos_j + icoord)) + d * sin(theta) / b * r_comp[icoord];
			break;
		case 4:  //"3out"
			a = Topology::virtual_site[ivtype].params[0];//parameter a
			b = Topology::virtual_site[ivtype].params[1];//parameter b
			c = Topology::virtual_site[ivtype].params[2];//parameter c (in nm-1) but this is really just a scaling factor, the nm-1 is to make nm from nm2
			//vector product r_ij x r_ik
			r_comp[0] = *(pos_j + 1) * *(pos_k + 2) - *(pos_j + 2) * *(pos_k + 1);//x component
			r_comp[1] = *(pos_j + 2) * *(pos_k)     - *(pos_j)     * *(pos_k + 2);//y component
			r_comp[2] = *(pos_j)     * *(pos_k + 1) - *(pos_j + 1) * *(pos_k);//z component
			for (icoord = 0; icoord < 3; icoord++)
				positions[3 * ind+icoord]=*(pos_i + icoord) + a * (*(pos_j + icoord)) + b * (*(pos_k + icoord)) + c * r_comp[icoord];
			break;

		}
		break;
	case 4:  //"4fdn"
		for (icoord = 0; icoord < 3; icoord++)
		{
			pos_k[icoord] = *(positions + 3 * virtual_source[(ind - ntotal) * 4 + 2] + icoord) - pos_i[icoord];//position of atom_k -pos_i
			if (pos_k[icoord] < -1)
				pos_k[icoord] += 2;
			else
			{
				if (pos_k[icoord] > 1)
					pos_k[icoord] -= 2;
			}
			pos_l[icoord] = *(positions + 3 * virtual_source[(ind - ntotal) * 4 + 3] + icoord) - pos_i[icoord];//position of atom_l -pos_i
			if (pos_l[icoord] < -1)
				pos_l[icoord] += 2;
			else
			{
				if (pos_l[icoord] > 1)
					pos_l[icoord] -= 2;
			}
		}
		
		a = Topology::virtual_site[ivtype].params[0];//parameter a
		b = Topology::virtual_site[ivtype].params[1];//parameter b
		c = Topology::virtual_site[ivtype].params[2];//parameter c
		for (icoord = 0; icoord < 3; icoord++)
		{
			r_ja[icoord] = a * *(pos_k + icoord) - *(pos_j + icoord);
			r_jb[icoord] = b * *(pos_l + icoord) - *(pos_j + icoord);
		}
		//cross  product
		r_m[0] = r_ja[1] * r_jb[2] - r_ja[2] * r_jb[1];
		r_m[1] = r_ja[2] * r_jb[0] - r_ja[0] * r_jb[2];
		r_m[2] = r_ja[0] * r_jb[1] - r_ja[1] * r_jb[0];
		d_r_m = sqrt(r_m[0] * r_m[0] + r_m[1] * r_m[1] + r_m[2] * r_m[2]);
		for (icoord = 0; icoord < 3; icoord++)
			positions[3 * ind+icoord] = *(pos_i+icoord) + c * r_m[icoord] / d_r_m;
	}
	GetMinImage(positions+3 * ind);//calculating minimum image


};

//calculate the coordinates of the virtual sites
void SimpleCfg::CalcVirtualCoords()
{
	int i;
	//calculates the initial coordinates of the virtual sites
	if (::debug)
		cout << "\nNOTE(" << ++note << "): Calculating the virtual site coordinates" << endl;
	for (i = 0; i < ntotvirtual; i++)
		CalcVirtualCoord(i + ntotal);
		
};

// in case of virtual sites resize the positions array and the finders
void SimpleCfg::ResizePosition()
{
	int i = 3*ntotal;
	ResizeArray(&i, 3*(ntotal + ntotvirtual), &positions, "positions", "SimpleCfg::ResizePosition");
	i = ntypes;
	ResizeArray(&i, ntypes + nvirtualtypes, &finder, "finder", "SimpleCfg::ResizePosition");
	//initialising the finder
	for (i = 0; i < ntypes+nvirtualtypes; i++)
		*(finder + i) = positions + (3 * cumul[i]);


};

void SimpleCfg::ConfigError(double pos, const char *file_name, const char *text)
{
	cout << "\n*****ERROR*****" << endl;
	cout << text << " is "<<pos<<" which is outside the tolarable range -3 -> +3, so even after trying to put" << endl;
	cout << "it back to the simulation box, it is outside the range of -1 -> +1!" << endl;
	cout << "Either the configuration file does not follow the RMC standard, or there is something wrong with" << endl;
	cout << "the loading of the file due to inappropriate line endings, which the program cannot resolve!" << endl;
	cout << "Check the " << file_name << " file, and try again!" << endl;
	cout << "Cannot continue, exiting..." << endl;
	CleanExit();
};

#ifdef _NO_PERIODIC
	//check and rescale the sample if necessary, as the reduced coordinates have to be inside the sample
	//basically between -0.5 and 0.5.
	void SimpleCfg::CheckSample()
	{
		int i,j,count=0,flag=0;
		double max=0,dist;
		double *pos,av[3]={0,0,0};
		double factor;

		if (RunParams::recentre_flag)
		{
			pos=positions;
			//first check the coordinates
			for (i=0; i<ntotal;i++)
			{
				for (j=0;j<3;j++)
				{
					av[j]+=*pos;//determining the centre of the sample
					pos++;
				}
			}
			
			//centering the simulation box, if necessary
			for (j=0;j<3;j++)
			{
				av[j]/=ntotal;//calculate the average
				if (fabs(av[j])>RunParams::rspacing[0]/boxedge/4)
				{
					cout << "\nWARNING(" << ++warn << "): The average reduced coordinates along the";
					switch (j)
					{
						case 0:
						{
							cout<<" x";
							break;
						}
						case 1:
						{
							cout<<" y";
							break;
						}
						case 2:
						{
							cout<<" z";
							break;
						}
					}
					cout << " axis are " << av[j] * boxedge << "," << endl;;
					cout<<"\twhich is outside the quarter of the binsize ";
					cout<<RunParams::rspacing[0]/4<<" A."<<endl;
					cout<<"\tThe coordinates will be centered."<<endl; 
					if (flag==0)//it was not saved yet
					{
						mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
						mystrcat(tempfilename, FILE_NAME_SIZE + 10, "_original_start.cfg");
						Save("\tOriginal coordinates before rescaling",tempfilename,3);
						flag=1;
					}
					pos=positions+j;
					for (i=0; i<ntotal;i++)
					{
						*pos-=av[j];
						pos+=3;
					}

				}
			}
		}
		//now see, if the sample is inside the sphere
		pos=positions;
		for (i=0; i<ntotal;i++)
		{
			dist=0;
			for (j=0;j<3;j++)
			{
				dist+=pow(*pos,2);
				pos++;
			}
			dist=sqrt(dist);
			if (dist*boxedge>RunParams::R0)
			{
				count++;
			}
			if (dist>max)
				max=dist;//find the largest
		}
		//The atom farthest from the centre should be 20th of the binsize form the surface for safety reasons
		factor=(RunParams::R0-(RunParams::rspacing[0]/20))/boxedge/max;
		if (count>0 || (factor*0.5-0.5)>TOLERANCE)
		{
			cout << "\nWARNING(" << ++warn << "): ";
			if (count>0)
			{
				cout<<"The coordinates of "<<count<<" atoms are outside the spherical sample of radius "<<RunParams::R0<<endl;
				cout<<"\tIt is assumed, that the sample is spherical, and the coordinates will be rescaled to fit inside "<<endl;
				cout<<"\tthe "<<RunParams::R0<<" radius sphere!"<<endl;
			}
			else
			{
				cout<<"The maximum atomic distance from the centre is "<<max*boxedge<<" A, which is smaller, than "; 
				cout<<"\tthe (radius of the spherical sample - rspacing[0]/20) "<<RunParams::R0-(RunParams::rspacing[0]/20)<<" A!"<<endl;
				cout<<"\tIt is assumed, that the sample is spherical, and the coordinates will be rescaled!"<<endl;
			}
						
			if (flag==0)//it was not saved yet
			{
				mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
				mystrcat(tempfilename, FILE_NAME_SIZE + 10, "_original_start.cfg");
				Save("\tOriginal coordinates before rescaling",tempfilename,3);
				flag=1;
			}
					
			//rescale the coordinates, keeping the R0 and boxedge fixed
			pos=positions;
		
			for (i=0;i<3*ntotal;i++)
			{

				*pos*=factor;
				pos++;
			}
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
			mystrcat(tempfilename, FILE_NAME_SIZE + 10, "start.cfg");
			Save("\tStarting configuration after rescaling",tempfilename);
			Save(RunParams::title);//write out the rescaled coordinates
		}
		
	};
#endif

#ifdef _AV_MOVE
	void SimpleCfg::SaveMovedDist()
	{
		int i;
		double *pos;
		ofstream file;
		CalcAvMove();//calculate the move for each atom 
		OpenFile(file,avmfilename,"SimpleCfg::SaveMovedDist",0);//open file, check, whether it was successfully opened

		file.setf(ios::scientific, ios::floatfield);
		file.setf(ios::right, ios::adjustfield);
		
		file<<" Moved distance for each atom"<<endl;
		
		file.precision(8);
		pos=moved_dist;
		for(i=0;i<ntotal;i++)
		{
			file<<i+1<<"\t"<<*pos++<<endl;
		};
		file.unsetf(ios::scientific);
		file.unsetf(ios::right);
		file.close();
	};
#endif

	