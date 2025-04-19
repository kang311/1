//source: utilities.cpp
//Last changed 24.01.2023

// functions used by various classes

#include<math.h>
#include <map>
#define _DEF_FILES //not redefine the file names included through global.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "utilities.h"

//if filename is given, then warning is displayed, if end of file was reached during skipping
int SkipLine(ifstream &file, const char *file_name, int termstat)
{
	int code;
	int lineend_code = 10;
	if (end_flag == 2)
		lineend_code = 13;
	if(!file) 
		return(0);
	else
	{
		do
		{
			code=file.get();
			if (code==EOF)
			{
				if (strlen(file_name)!=0 && ::debug)//only gives warning, if filename is given
					cout << "\nWARNING(" << ++warn << "): End of file "<<file_name<<" was reached during line skipping!"<<endl;
				if (termstat)//this is a serious error, data should follow
				{
					cout<<"\n*****ERROR*****"<<endl;
					cout<<"There are not enough data, terminating program..."<<endl;
					CleanExit();
				}
				else
					return(0);
			}
		}
		while (code!=lineend_code);
				
	}
	
	return(1);//no error
};


//----------------J2 formatting utility to write numbers cleanly----
//this utility creates a printing space of 2 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in test of normalisation

ostream &J2(ostream &os) {
	os.width(2);
	os.fill(' ');
 	return os;
}
//----------------J4 formatting utility to write numbers cleanly----
//this utility creates a printing space of 4 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in test of normalisation

ostream &J4(ostream &os) {
	os.width(4);
	os.fill(' ');
 	return os;
}
//----------------J10 formatting utility to write numbers cleanly----
//this utility creates a printing space of 10 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J10(ostream &os) {
	os.width(10);
	os.fill(' ');
 	return os;
}




//----------------J7 formatting utility to write numbers cleanly----
//this utility creates a printing space of 7 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in test of normalisation

ostream &J7(ostream &os) {
	os.width(7);
	os.fill(' ');
 	return os;
}

//----------------J12 formatting utility to write numbers cleanly----
//this utility creates a printing space of 12 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J12(ostream &os) {
	os.width(12);
	os.fill(' ');
 	return os;
}

//----------------J13 formatting utility to write numbers cleanly----
//this utility creates a printing space of 13 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J13(ostream &os) {
	os.width(13);
	os.fill(' ');
	return os;
}
//----------------J16 formatting utility to write numbers cleanly----
//this utility creates a printing space of 16 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J16(ostream &os) {
	os.width(16);
	os.fill(' ');
 	return os;
}


//----------------J20 formatting utility to write numbers cleanly----
//this utility creates a printing space of 20 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J20(ostream &os) {
	os.width(20);
	os.fill(' ');
 	return os;
}

//----------------J23 formatting utility to write numbers cleanly----
//this utility creates a printing space of 23 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J23(ostream &os) {
	os.width(23);
	os.fill(' ');
	return os;
}

//----------------J27 formatting utility to write numbers cleanly----
//this utility creates a printing space of 27 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J27(ostream &os) {
	os.width(27);
	os.fill(' ');
	return os;
}

//----------------J30 formatting utility to write numbers cleanly----
//this utility creates a printing space of 30 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J30(ostream &os) {
	os.width(30);
	os.fill(' ');
	return os;
}
//----------------J35 formatting utility to write numbers cleanly----
//this utility creates a printing space of 35 characters filled with spaces
//for writing the next output ( with << )
//justifying and adjusting the precision of real numbers is done externally
//see use in SimpleCfg::save(ifstream &file)

ostream &J35(ostream &os) {
	os.width(35);
	os.fill(' ');
 	return os;
}
//------------rounds number to the nearest integer----------------
int RoundNearInt(double number)
{
	int result;
	result=(int)number;
	if ((number-(int)number)>=0.5)
		result++;
	return(result);

};

void GetMinImage(double *coord)
{
	int i;
	for (i=0;i<3;i++)
	{
		if (coord[i]>=1) 
			coord[i]-=2;
		else 
			if (coord[i]<=-1) 
				coord[i]+=2;
	}
}
//----------------random number generator from numerical recipies--------
//returns a number flatly distributed between 0 and 1
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Ran1(longint &idum)
{
        longint j;
        longint k;
        static longint iy=0;
        static longint iv[NTAB];
        double temp;
        
        if(idum <=0 || !iy) {
                if (-idum < 1) idum=1;
                else idum = -idum;
                for (j=NTAB+7;j>=0;j--) {
                        k=idum/IQ;
                        idum=IA*(idum-k*IQ)-IR*k;
                        if (idum<0) idum += IM;
                        if (j<NTAB) iv[j]=idum;
                }
                iy=iv[0];
        }
        k=idum/IQ;
        idum=IA*(idum-k*IQ)-IR*k;
        if (idum <0) idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j]=idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
};

//Error message, when no new array could be created
void NoArray(const char *routine, const char *array_name)
{
	cout << "\n*****ERROR*****" << endl;
	cout<<array_name<<" array in "<<routine<<" routine cannot be created most probably because of the shortage of memory."<<endl; 
	CleanExit();
};



//-------checking the file status after reading--------
int CheckReadFileState(ifstream &file, const char *routine_name, const char* file_name)
{
	if (file.rdstate() & file.failbit)//there was some error during reading
	{
		cout << "\nWARNING(" << ++warn << "): In routine "<<routine_name<<" data could not be read from "<<file_name<<" file!"<<endl;
		if (file.rdstate() & file.eofbit)
			cout<<"\tEnd of file was reached during loading!"<<endl;
		else
			cout<<"\tLoading failed, possibly due to wrong file format!"<<endl;
		
		return(0);
	}
	else
	{
		if (file.rdstate() & file.eofbit)//end of file was reached after the last required data, clear rdstate
			file.clear(ios::goodbit);//set it back to goodbit, it is needed, if the file variable is used again
		
	}
	return(1);

};


//------------------checking the file status, if the file is already open for saving--------------
//checking the file status after reading
int CheckFileState(ifstream &file, const char *routine_name, const char* file_name)
{
	if (file.is_open()==0 || !file)   //verification that file is open
	{
		if (strcmp(file_name,"")!=0)
			cout << "\nWARNING(" << ++warn << "): In routine "<<routine_name<<" file "<<file_name<<" is not open for reading!"<<endl;
		file.clear(ios::goodbit);//set it back to goodbit, it is needed, if the file variable is used again
		return(0);//the file is not open 
		
	}
	else
		return(1);//the file is open

}

//------------------checking the file status, if the file is already open for saving--------------
void CheckFileState(ofstream &file, const char *routine_name, const char* file_name)
{
	int repeat_status;
	char key;
	char *fname;
	SetArraysize(&fname,FILE_NAME_SIZE,"fname","CheckFiletate"),
	mystrcpy(fname,FILE_NAME_SIZE,file_name);

	
	do 
	{
		repeat_status=0;//do not repeat opening by default
		if (file.is_open()==0 || !file)   //verification that file is open
		{
			cout << "\nWARNING(" << ++warn << "): In routine "<<routine_name<<" file "<<file_name<<" could not be opened for saving!"<<endl;
			cout<<"\tIf the file is used by another application, close it first."<<endl;
			cout<<"\tTry again? (Press any key!)"<<endl;
			cin>>key;
			file.clear(ios::goodbit);//set it back to goodbit, it is needed, if the file variable is used again
			repeat_status=1;
			CleanOpen(file,fname);
		}
	}
	while (repeat_status);
};
//------------------checking the file status after opening it for saving----------------
void OpenFile(ofstream &file, const char *file_name, const char *routine_name,int flag)
{
	//flag: 0:text mode, 1: binary mode, 2: text mode append
	
	int repeat_status;
	char key;

	
	do 
	{
		repeat_status=0;//do not repeat opening by default
		if (flag==1)
			CleanOpen(file,file_name,1);//open it in binary mode
		else
		{
			if (flag==0 || flag==3)
				CleanOpen(file,file_name);//open it in text mode 
			else
				CleanOpen(file,file_name,2);//open it in text mode to append
		}

		if (file.is_open()==0 || !file)   //verification that file is open
		{
			repeat_status=1;
			cout << "\nWARNING(" << ++warn << "): In "<<routine_name<<" file "<<file_name<<" could not be opened for saving!"<<endl;
			cout<<"\tIf the file is used by another application, close it first."<<endl;
			cout<<"\tTry again? (Press any key!)"<<endl;
			cin>>key;
			file.clear(ios::goodbit);//set it back to goodbit, it is needed, if the file variable is used again
			
		};
	}
	while (repeat_status);

};

//check the text file ending
//-2: this is not an expected ending
//-1: could not be determined due to file reading error
// 0: linux
// 1: windows
// 2: mac 
void CheckEndLine()
{
	ofstream file;
	ifstream infile;
	char c1, c2;

	end_flag = -2;//default
	CleanOpen(file, "t.txt");
	//check the text file ending
	if (file.is_open() == 0 || !file)   //verification that file is open
	{
		end_flag = -1;//could not be determined
		return;
	}
	file << endl;
	file.close();
	CleanOpen(infile, "t.txt",1);
	if (infile.is_open() == 0 || !infile)   //verification that file is open
	{
		end_flag = -1;//could not be determined
		return;
	}
	c1 = infile.get();
	c2 = infile.get();
	if (c2 == EOF)
	{
		if (c1 == 10)
			end_flag = 0;
		else
		{
			if (c1 == 13)
				end_flag = 2;
			else
				end_flag = -1;
		}
	}
	else
	{
		if (c1 == 13 && c2 == 10)
			end_flag = 1;
	}
	infile.close();
	remove("t.txt");
	return;
};

//converting the text file to have endings corrsponding to the platform and opening it for reading
//THIS HAS TO BE USED for opening a text file the first time
//write_flag = 0 only reading will be performed later with this file,
//			   1 writing also occurs, file renaming has to be performed, temporary name cannot be used later (default)
void SafeOpenTextFile(ifstream &file, const char *filename,  int write_flag)
{
	int repeat_status;
	int count_r = 0, count_n = 0, count_rn = 0;
	string text;
	char c1,c2=-1;
	char key;
	char newname[FILE_NAME_SIZE + 15];
	char tempname[FILE_NAME_SIZE + 15];
	bool convert = false;
	bool file_error = false;

	ofstream outfile;
	
	CleanOpen(file, filename, 1);//open as binary, to be extract the exact line endings
	if (CheckFileState(file, filename) == 0)
		return;//let the calling function handle it
	
	
	c1 = file.get();
	if (c1 == EOF)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << filename << " file is empty, cannot proceed, exiting..." << endl;
		CleanExit();
	}
	//First check the line endings
	switch (c1)
	{
	case 10:
		count_n++;
		break;
	case 13:
		count_r++;
	}

	do
	{
		c2 = file.get();
		switch (c2)
		{
		case 10:
			if (c1 == 13)
			{
				count_rn++;
				count_r--;
			}
			else
				count_n++;
			break;
		case 13:
			count_r++;
		}
		c1 = c2;
	} while (c2 != EOF);
	file.close();
	switch (end_flag)
	{
	case 0: //linux
		if (count_rn > 0 || count_r > 0)
			convert = true;
		else
			if (count_n == 0)
				file_error = true;
		break;
	case 1: //windows
		if (count_n > 0 || count_r > 0)
			convert = true;
		else
			if (count_rn == 0)
				file_error = true;
		break;
	case 2: //mac
		if (count_rn > 0 || count_n > 0)
			convert = true;
		else
			if (count_r == 0)
				file_error = true;
		break;
	default: //the right ending could not be determined, try to convert anyhow
		convert = true;
	}

	if (file_error)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << filename << " text file does not contain any recognized line ending (neither linux \n, windows \r\n nor mac \r)!" << endl;
		cout << "Replace the file with a suitable one, end try again! Exiting..." << endl;
	}

	if (!convert)
	{
		CleanOpen(file, filename);
		return;
	}

	//The file has to be converted
	CleanOpen(file, filename, 1);
	mystrcpy(tempname, FILE_NAME_SIZE + 15, filename);
	mystrcat(tempname, FILE_NAME_SIZE + 15, "#temp");
	CleanOpen(outfile, tempname);
	if (outfile.is_open() == 0 || !outfile)   //verification that file is open
	{
		cout << "\nWARNING(" << ++warn << "): Input text file cannot be converted to acceptable format!" << endl;
		return;
	}
	

	c1 = file.get();
	if (c1 != 10 && c1 != 13)
		outfile<<c1;
	do
	{

		c2 = file.get();
		switch (c2)
		{
		case 10:
			outfile << endl;
			break;
		case 13:
			if (c1 == 13)
				outfile << endl;
			break;//wait to see, whether the next is 10
		default:
			if (c1 == 13)
				outfile << endl;
			if (c2 != EOF)
				outfile << c2;
		}
		c1 = c2;
	} while (c2 != EOF);
	file.clear(ios::goodbit);
	file.close();
	outfile.close();
	mystrcpy(newname, FILE_NAME_SIZE + 15, filename);
	mystrcat(newname, FILE_NAME_SIZE + 15,"#ori");
	ifstream infile;
	CleanOpen(infile, newname);
	if (CheckFileState(infile, "SafeOpenTextFile") == 1)//file exist
	{
		infile.close();
		repeat_status = remove(newname);//do not repeat opening by default
		do
		{
			if ( repeat_status!= 0)//could not remove it
			{
				cout << "\nWARNING(" << ++warn << "): In routine SafeModeTextFile file " << newname << " could not be overwritten" << endl;
				cout << "\twith a new one during text file conversion for file " <<filename<<". The file most probably is open, close it!"<< endl;
				cout << "\tTry again? (Press any key!)" << endl;
				cin >> key;
				repeat_status = remove(newname);
			}
		} while (repeat_status);
	}
	
	repeat_status = rename(filename, newname);
	do
	{
		if (repeat_status != 0)
		{
			cout << "\nWARNING(" << ++warn << "): Input text file " << filename << " could not be renamed to " << newname << endl;
			cout << "\tSo it cannot be replaced with the converted one containing proper line endings!" << endl;
			cout << "\tMost probably the " << filename <<" file is in use, so close it." << endl;
			cout << "\tTry again? (Press any key!)" << endl;
			cin >> key;
			repeat_status = rename(filename, newname);
		}
	} while (repeat_status);
	
	repeat_status = rename(tempname,filename);
	do
	{
		if (repeat_status != 0)//this not likely to happen, but check anyway
		{
			cout << "\nWARNING(" << ++warn << "): Input text file " << filename << " could not be replaced by the converted one containing " << endl;
			cout << "\tthe proper line endings! Close the "<<tempfilename<<" file and press any key to try again!" << endl;
			cin >> key;
			repeat_status = rename(tempname, filename);

		}
	} while (repeat_status);

	if (::debug || strcmp(filename, cfgfilename) == 0 || strcmp(filename, datfilename) == 0)//debug is not read for cfg and dat
	{
		cout << "\nWARNING(" << ++warn << "): Input text file " << filename << " contained unsuitable line endings:" << endl;
		switch (end_flag)
		{
		case 0:
			if (count_rn > 0)
				cout << "\t\tNumber of windows-style \\r\\n line endings: " << count_rn << endl;
			if (count_r > 0)
				cout << "\t\tNumber of mac-style \\r line endings: " << count_r << endl;
			break;
		case 1:
			if (count_n > 0)
				cout << "\t\tNumber of linux-style \\n line endings: " << count_n << endl;
			if (count_r > 0)
				cout << "\t\tNumber of mac-style \\r line endings: " << count_r << endl;
			break;
		case 2:
			if (count_rn > 0)
				cout << "\t\tNumber of windows-style \\r\\n line endings: " << count_rn << endl;
			if (count_n > 0)
				cout << "\t\tNumber of linux-style \\n line endings: " << count_n << endl;
			break;
		default://this should not happen
			if (count_n > 0)
				cout << "\t\tNumber of linux-style \\n line endings: " << count_n << endl;
			if (count_rn > 0)
				cout << "\t\tNumber of windows-style  \\r\\n line endings: " << count_rn << endl;
			if (count_r > 0)
				cout << "\t\tNumber of mac-style  \\r line endings: " << count_r << endl;
		}
	}
		
	CleanOpen(file, filename);
	if (::debug || strcmp(filename, cfgfilename) == 0 || strcmp(filename, datfilename) == 0)
	{
		cout << "\tThe file was converted into the native file format of the operating system," << endl;
		cout << "\tand the original file was renamed to " << newname << endl;
	
	}
};

//converts a number's digit into their character representation, and put it into conv_numb array
//for example: integer=4012 to *char="4012"; integer=-521 to *char"-521" 
void IntToStr(char **conv_number, int number)
{
	int i, digits, length, remainder;
	int pos_number;
	
	if (*conv_number != NULL)
		delete [] *conv_number;
	
	if (number == 0)
	{
		SetArraysize(conv_number, 2, "conv_number", "IntToStr");
		(*conv_number)[0] = '0';
		(*conv_number)[1] = '\0';
		return;
	}
	
	pos_number = abs(number);
	digits = (int)log10(pos_number) + 1;
	if (number < 0)
		length = digits + 2;//to have a place for the minus 
	else
		length = digits + 1;
	SetArraysize(conv_number, length, "conv_number", "IntToStr");//the last place is for the terminating character
	if (number < 0)
		(*conv_number) [0]= '-';
	(*conv_number)[length - 1] = '\0';//just to put the terminating character to its place,

	//this cycle is cutting the last digit, and places it into its proper place in the character array
	for (i = length - 2; i >= length - digits - 1; i--)
	{
		remainder = (pos_number % 10);
		switch (remainder)
		{
		case 0:
		{
			(*conv_number)[i] = '0';
			break;
		}
		case 1:
		{
			(*conv_number)[i] = '1';
			break;
		}
		case 2:
		{
			(*conv_number)[i] = '2';
			break;
		}
		case 3:
		{
			(*conv_number)[i] = '3';
			break;
		}
		case 4:
		{
			(*conv_number)[i] = '4';
			break;
		}
		case 5:
		{
			(*conv_number)[i] = '5';
			break;
		}
		case 6:
		{
			(*conv_number)[i] = '6';
			break;
		}
		case 7:
		{
			(*conv_number)[i] = '7';
			break;
		}
		case 8:
		{
			(*conv_number)[i] = '8';
			break;
		}
		case 9:
		{
			(*conv_number)[i] = '9';
			break;
		}


		}//end switch
		pos_number -= remainder;
		pos_number /= 10;//cut the last digit
	}
	return;
};



bool ErrLackValue(std::string &inputline)
//This function complains about lack of parameter in inputline and stops the run
{
  cout<<"\nERROR: lack of value in the following line:"<<endl<<inputline<<endl;
  return false;
}

void WrapFree(const std::string &input, std::string &keyword, std::string &parameters)
//Splits input format into keyword and parameters
{
  std::size_t le=input.find_first_of('=');
  keyword=input.substr(0,le);
  keyword.erase(std::remove(keyword.begin(),keyword.end(),' '),keyword.end());
  std::transform(keyword.begin(), keyword.end(), keyword.begin(), ::toupper);
  parameters=input.substr(le+1,input.length()-le-1);
  std::transform(parameters.begin(), parameters.end(), parameters.begin(), [](char ch) {if((ch=='\r')||(ch=='\t')) ch=' '; return ch;});
  Ltrim(parameters);
  Rtrim(parameters);
  return;
}

string ExtractChemSymbol(string &symbol)
{
	unsigned int i;
	bool digit2 = false, valence2 = false;
	string st_chem;
	if (symbol.compare("dummy") == 0)
	{
		return symbol;
	}
	//check if it is a valid chemical symbol, first character should be upper case, if there is second character for the symbol should be lower case,
	//can be followed by some specifier : nujber and +/- for ions or  V for valence state 
	for (i = 0; i < symbol.length(); i++)
	{
		switch (i)
		{
		case 0:
			if (std::isupper(symbol[i]) == 0)//not an upper case
				ErrBadSymbol(symbol);
			else
			{
				if (symbol.compare("D") == 0)
				{
					st_chem.push_back('H');
					return (st_chem);
				}
				st_chem.push_back(symbol[i]);
				if (symbol.length() == i + 1)//there is nothing else
					return (st_chem);
			}
			break;
		case 1:
			if (std::islower(symbol[i]) > 0)//lower case
			{
				st_chem.push_back(symbol[i]);//can be the second char of the symbol
				if (symbol.length() == i + 1)//there is nothing else
					return (st_chem);
			}
			else
			{
				if (std::isdigit(symbol[i]) != 0)
					digit2 = true;//digit
				else
				{
					if (symbol.compare(0, 2, "CV")==0)
					{
						valence2 = true;//V for the valence
						if (symbol.length() != i + 1)//there cannot be anything else
							ErrBadSymbol(symbol);
						else
							return (st_chem);
					}
					else
					{
						if (symbol.compare(0, 2, "HS")==0)
						{
							if (symbol.length() != i + 1)//there cannot be anything else
								ErrBadSymbol(symbol);
							else
								return (st_chem);
						}
						else
						{
							if (symbol.compare(i, 1, "+") == 0 || symbol.compare(i, 1, "-") == 0)//no digit was given before +/- sign, mend it
							{
								if (symbol.length() != i + 1)//there cannot be anything else
									ErrBadSymbol(symbol);
								else
								{
									symbol.insert(i, 1, '1');
									return (st_chem);
								}
							}
							else
								ErrBadSymbol(symbol);
						}
					}
				}
			}
			break;
		case 2:
			if (digit2 && (symbol.compare(i, 1, "+")==0 || symbol.compare(i, 1, "-")==0)) //F1-
			{
				if (symbol.length() != i + 1)//there cannot be anything else
					ErrBadSymbol(symbol);
				else
					return (st_chem);
			}
			else if ((digit2 == false && valence2 == false) && symbol.compare(i, 1, "V")==0)//SiV
			{
				if (symbol.length() != i + 1)//there cannot be anything else
					ErrBadSymbol(symbol);
				else
					return (st_chem);
			}
			else if (symbol.compare(i, 1, "+") == 0 || symbol.compare(i, 1, "-") == 0)//Li+ , no digit was given before +/- sign, mend it
			{
				if (symbol.length() != i + 1)//there cannot be anything else
					ErrBadSymbol(symbol);
				else
				{
					symbol.insert(i, 1, '1');
					return (st_chem);
				}
			}
			else if (!((digit2 == false && valence2 == false) && std::isdigit(symbol[i]) != 0))//Li1+
			{
				ErrBadSymbol(symbol);

			}
			break;
		case 3:
			if (symbol.compare(i, 1, "+")==0 || symbol.compare(i, 1, "-")==0)
				return (st_chem);
			else
				ErrBadSymbol(symbol);
			break;
		default:
			ErrBadSymbol(symbol);
		}

	}
	return st_chem;
}






