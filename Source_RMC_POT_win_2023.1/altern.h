//header altern.h
//Last changed 01.03.2023

//Might need some alteration in case of other platforms! 
//FOR LINUX platform the supplied Makefile will define _GNU_LINUX, so this part will be skipped,
//does not matter, which option is here switched on or off. The command line arguments passed to the Makefile
//will define, which option is switched on.

//_AENET can be compiled only on LINUX, so it is not given below 
#ifndef _ALTERN_H
#define _ALTERN_H

//---FOR WINDOWS-MICROSOFT VC-----

#ifndef _GNU_LINUX //not LINUX platform
#define _MICROSOFT_VC//If switched on, compiling using MS Visual C++ with WINDOWS is assumed. 
//THE FOLLOWING PARAMETERS WILL REGULATE THE WAY THE EXE FILE IS CREATED
//IF AN OPTION (A #define _parameter) IS LEFT HERE IN THE CODE, THAN THE
//CODE SEGMENTS BELONGING TO THIS PARAMETER WILL BE BUILT IN TO THE CODE 
//IF AN OPTION IS NOT NEEDED IT HAS TO BE COMMENTED OUT HERE, AND THE RELATED
//CODE SEGMENTS WILL NOT BE ADDED TO THE EXE FILE TO MAKE EXECUTION FASTER.

#ifndef _CL_COMP//will be compiled with CL online Microsoft compiler, disregard the option switches
//- - - - - - -
//the recent RMC++ code was never tested on MAC, there is no guarantee, that it can be compiled without
//any change, this was only kept from historical reasons!!!
//#define _CODE_WARRIOR_MAC //DEFAULT is the PC or UNIX version, this option has to be commented out normally
//***********THESE ARE FOR DEVELOPERS, not all can be switched in linux Makefile****************
//in _TEST_MODE:
//	-	the random number generator will be initialized with the same value, so the 
//		random numbers will follow the same pattern at each run to make testing easier
//	-	the run will stop at ngenerated or naccepted steps, can be set by run time (-/+ value) with fixed format,
//		with free format RUN-TIME not processed, LAST-GEN or LAST-ACC should be used (defaults in units.h) 
//#define _TEST_MODE //this has to be disabled, if the program is used in normal running mode
//
//#define _NO_CHI2_POT_ACC//The value of chi2_pot will not determine the acceptance, every 10th step will be accpeted based on the chi2_pot, 
//						    so the rounding error does not distroy testing by causing the simalation path to fork, thi sis not in Makefile or manual
//- - - - - - - - -
//#define _WRITE_CHI2_DETAIL//if it is on, the chi2 components will be written in to the *.chi file for every rejected move, not in manual or Makefile
//------
//#define WRITE_TH_DURATION 1//Only has effect, if _TEST_MODE is ON. If 1, writes the times of different tasks for 
//                             each thread in the *.hst file. If zero, it is nor written.
//-------
//#define _TH_ORDER //for thread order print out in case of threading problems
//****************END OF DEVELOPERS SECTION*********************************
//****************options concerning the working of the program*************
//- - - - - - - 
//#define _ADVANCED_GEOM_CONST//use advanced geometric constraints as CommonNeighConst, Second Neighbour Constraint and Bond Valence Sum Constraint
//- - - - - - -
//#define _LOCAL_INV//use local invariance added to chi2
//- - - - - - -
//#define _LOCAL_INV_NS//do not save local invariance histogram
//- - - - - - -
//#define _MUL_SCAT_VECTOR//if it is on, then for all the S(Q), F(Q) and F(g) data sets the normally calculated S(Q), F(Q) and F(g) will be multiplied by Q and F(g) by g
//- - - - - - -
//#define _NO_PERIODIC//if it is on, periodic boundary conditions not used
//- - - - - - - 
//#define _READ_EDIFF//if it is on, the program will read the number of ediff data sets among Xray and EXAFS in case of fixed format input, otherwise not, so the old format can be used for the fixed format *.dat files. Free format can read ediff regerdless this option
//- - - - - - - 
//#define _R_SWITCH_GR_MIN_1// If this is option is on, then r^N*g(r) is assumed at reading the experimental data for g(r) sets, and this will be calculated and fitted. 
//                             N is specified by R_SWITCH_POWER_DEF which is located in units.h. 
//                             Default value can be also changed with key R-SWITCH-POWER in case of free format parameter file.
//- - - - - - - 
//#define _VIBR_AMP//if it is on, than the histogram will be convoluted with a Gauss-distribution
//*************options connected to outputs*********************************
//- - - - - - -
//#define _AV_MOVE//if it is on, the move in each direction separately accumulated for every atom and displayed 
//- - - - - - -
//_NEI, _NEIE only has effect, if neighbour list is calculated 
//#define _NEI //if it is on, saving the neighbourlist into the *.nei file
//#define _NEIE //if it is on, saving the neighbourlist and also the squared neighbour distances and vector components into the *.nei file
//**********options, which were important for early versions****************
//- - - - - - - - - - 
//#define _ATLAS//use the ATLAS libraries for matrix operations
//- - - - - - - - 
//#define _OLD_HEADER //if the old style header (name.h) are used, THE NEW STYLE headers are used by default
//                      this option has to be commented out normally
//- - - - - - -
//Changed in version 2.5 form _USE_INT64 to eliminate the need to switch this on for int64 compilation, as nowdays this is the default.
//Now include this if you DO NOT WANT to use 64 bit integers
//#define _USE_INT32//NOT to use integer 64 (8byte) for longint (4 byte will be used), otherwise int64 will be used(8 byte) 


#endif //end of _CL, disregard option switches
#endif//end of not LINUX platform


//______________________________________________________________________________________________________________________
//DO NOT TOUCH THE DEFINE STATEMENTS BELOW THIS!!!!
//to set up the local histogram calculations, some arrays are shared with _NO_PERIODIC 
#ifdef _LOCAL_INV
#define _USE_LOCAL_INV
#endif

//To make sure, that local invariance in not used
#ifdef _NO_PERIODIC
#define _NOP_VIBR //noperiodic or vibr amp or both is defined
#define _USE_LOCAL_INV//to be able to use the local histograms
#ifdef _LOCAL_INV
	//cannot use local invariance with non-poriodic boundary conditions, LOCAL INV will be switched off
#define _WRITE_MESSAGE
#undef _LOCAL_INV
#endif
#endif

#ifdef _VIBR_AMP
#define _NOP_VIBR //noperiodic or vibr amp or both is defined
#endif

#ifdef _NO_CHI2_POT_ACC
#ifdef _TEST_MODE
#define _NO_POT_CHI2 //needed to make sure, that only in case of both switched on this is activated
#endif
#endif
//- - - - - - -

#include "global.h"

//----------------------------------------------------------------

//to make sure, that _NEI is defined, if _NEIE is defined
#ifdef _NEIE
#ifndef _NEI
#define _NEI
#endif
#endif

#ifdef _OLD_HEADER 
	#include<iostream.h>
	#include<fstream.h>
	#include<iomanip.h>
#else
	#include<iostream>
	#include<fstream>
	#include<iomanip>
	#include<complex>
	//These are needed for the new style header file, std::list cannot be used, as it can be mixed with using declaration
	using std::ifstream;
	using std::ofstream;
	using std::istream;
	using std::ostream;
	using std::stringstream;
	using std::cin;
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::ios;
	using std::string;
	using namespace std::complex_literals;

#endif

#ifdef _USE_INT32
	typedef int longint;//simple integer will be used for longint as well
#else
	//64bit inetegr will be used, names are operation system dependent before int64_t. Does not want to change to int64_t yet to preserve downward compatibility
	#ifdef _GNU_LINUX//GNU-LINUX systems: only included in the code, if _GNU_LINUX switch is defined in altern.h!!!
		typedef long int longint;
	#else//--End of GNU/LINUX specific--
	
		#ifdef _MICROSOFT_VC //Microsoft Visula C++: only included in the code, if _MICROSOFT_VC switch is defined in altern.h!!!
			typedef  __int64 longint;//---MICROSOFT Visual C++ Specific-----
 
			//There is a compiler bug for Microsoft Visual C++ for earlier versions, << and >> operators were not defined for __int64!
			//So they are declared here  For Visual Studio 2015 seemingly it has been fixed, and if these lines are remaining in the code, 
			//this generates error for ambigouos operators PUT them back, if needed
	//		ostream& operator<<(ostream& os, __int64 i );
	//		istream &operator>>(istream& is, __int64 &i );
		#endif//--End of Microsoft Visual C++ specific--
	#endif
#endif

#ifdef _GNU_LINUX//GNU-LINUX systems: only included in the code, if _GNU_LINUX switch is defined in altern.h!!!
  #include <pthread.h>
#endif





#ifdef _MICROSOFT_VC //Microsoft Visula C++: only included in the code, if _MICROSOFT_VC switch is defined in altern.h!!!

inline char * mystrcpy(char *target, int length, const char *source)
{
	strcpy_s(target, length, source);//no checking
	return target;
}
inline char * mystrncpy(char *target, int length, const char *source, size_t count)
{
	strncpy_s(target, length, source, count);//no checking
	return target;
}
inline char * mystrcat(char *target, int length, const  char *source)
{
	strcat_s(target, length, source);//no checking
	return target;
}
inline char *mystrncat(char *target, int length, const  char *source,size_t count)
{
	strncat_s(target, length, source,count);//no checking
	return target;
}
inline int mysprintf(char *buffer, size_t sizeOfBuffer, const char *format, longint data)
{
	return sprintf_s(buffer, sizeOfBuffer, format, data);
}
#else		
inline char * mystrcpy(char *target, int length, const char *source)
{
	strcpy(target, source);//no checking
	return target;

}
inline char * mystrncpy(char *target, int length, const char *source, size_t count)
{
	strncpy(target, source, count);//no checking
	return target;

}
inline char * mystrcat(char *target, int length, const char *source)
{
	strcat(target, source);//no checking
	return target;
}
inline char *mystrncat(char *target, int length, const char *source,size_t count)
{
	strncat(target, source,count);//no checking
	return target;
}
inline int mysprintf(char *buffer, size_t sizeOfBuffer, const char *format, longint data)
{
	return sprintf(buffer, format, data);//no checking

}
#endif




//functions with different source files depending on the compiler

//to open files for writing, the way it is done depending on flag, replaces file.open(filename....)
//COMPILER DEPENDENT:NEW or OLD style headers
//ios::nocreate is not recognised by codewarrior and is necessary for VisualC++, if it is using the old .h header files, 
//but it is not part of the new standard C++ library files
//VisualC++ only if the old libraries are used

template < class InFile >
void CleanOpen(InFile &file, const char *file_name, int flag=0)
{
	//flag: 0:text mode, 1: binary mode, 2: text mode append
	//to replace file.open(filename)
#ifdef _OLD_HEADER //only in case of OLD STYLE VisualC++ and UNIX/LINUX header files (iomanip.h, iostream.h fstream.h ...)
	switch (flag)
	{
		case 1://open binary files
		{
			file.open(file_name,ios::binary,ios::nocreate);
			break;
		}
		case 2://open text files for appending
		{
			file.open(file_name,ios::app,ios::nocreate);
			break;
		}
		default:
		{
			file.open(file_name,ios::nocreate);
		}
	}
#else //DEFAULT: NEW SYTLE HEADERS: iomanip, iostream, fstream ...),nocreate is not a member of the standard C++ library!!!
	switch (flag)
	{
		case 1://open binary files
		{
			file.open(file_name,ios::binary);
			break;
		}
		case 2://open text files for appending
		{
			file.open(file_name,ios::app);
			break;
		}
		default:
		{
			file.open(file_name);
		}
	}
#endif
}



void CleanExit();//to replace exit(EXIT_FAILURE);
//needed because on windows, the console windows collapses at exit and therefore
//prevents the reading of any error message
#endif


