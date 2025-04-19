//source: altern.cpp
//Last changed 27.02.2023

//functions with different source files depending on the compiler

#include"altern.h"


	
//================================================================================

#ifdef _CODE_WARRIOR_MAC //these are the CodeWarrior for Mac versions
	void CleanExit()
	{
		//to replace exit(EXIT_FAILURE);
		//needed because on windows, the console windows collapses at exit and therefore
		//prevents the reading of any error message
		cout<<"\n bye...\n";
		exit(EXIT_FAILURE);
	}
#else
	void CleanExit()//THIS IS THE DEFAULT (PC, UNIX/LINUX)
	{
		//to replace exit(EXIT_FAILURE);
		//needed because on windows, the console windows collapses at exit and therefore
		//prevents the reading of any error message
		char key;
		cerr<<"Quit (y)?->\n";
		cin>>key;
		exit(EXIT_FAILURE);
	};	
#endif
	
//==================================================================================



//There is a compiler bug for Microsoft Visual C++, << and >> operators were not defined 
//for __int64!
//So they are defined here
#ifdef _MICROSOFT_VC
#ifndef _USE_INT32
	ostream& operator<<(ostream& os, __int64 i )
	{
	   char buf[21];
	   mysprintf(buf,21,"%I64d", i );
	   os << buf;
	   return os;
	}

	//this can handle the line break, the space and the tabs
	istream &operator>>(istream& is, __int64 &i )
	{
		char separator;
		char instring[21];
		do 
		{
			is.get(separator);
		}
		while (separator==' ' || separator=='\t' || separator =='\n');
		is.putback(separator);
		for (i=0;i<20;i++)
		{
			is.get(separator);
			if (separator=='\t' || separator=='\n' || separator==' ')
			{
				is.putback(separator);
				instring[i]='\0';
				break;
			}
			else
				instring[i]=separator;
		}
		i=_atoi64(instring);
		return is;
	}
#endif
#endif//--end of Microsoft Visula C++ specific
