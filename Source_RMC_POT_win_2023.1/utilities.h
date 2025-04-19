//header utilities.h
//Last changed 26.02.2023

// prototypes of  functions used in the definition of classes

#include "data.h"
#include <algorithm>
#include <string>

enum T_interaction_type {BOND, ANGLE, DIHEDRAL};
enum T_comp_dir {d_include, d_ifdef, d_else, d_endif};//these are directives for the preprocessor in the topology file, some might be neede
//the force field files, nor used, just skipped
enum T_force_field {OPLSAA, Encads, Encadv, CHARMM,AMBER,UNKNOWN_GENERATE_PAIRS,UNKNOWN_READ_PAIRS};


//skip a line in reading a file
//if filename is given,than warning is displayed, if end of file was reached during skipping
int SkipLine(ifstream &file, const char *file_name="", int termstat=0);

//converting the text file to have endings corrsponding to the platform, and opening it
void SafeOpenTextFile(ifstream &file,const char *filename,int write_flag=1);

//check the text file ending
void CheckEndLine();

//calculating the minimum image vector
void GetMinImage(double *coord);

//template function, has to be defined here, as it will not compile, if it is in utilities.cpp
//only can read from the present line, gives warning, if not enough data, and do not try to read from the next line!
//The work of the fuction depend on the status parameter.
//IMPORTANT: IT SEEMS, that the actual type for formal type InT in the parenthesis (in other words, the type of the def_value
//parameter) has precedence in the determination of the actual type over the return type of the function. Although the value 
//of the def_value parameter is only important in case of the status=0, 2, 4 (optional parameters),its type DETERMINES the
//type of InT, so always numerical values matching the type of the variable we want 
//to read in (the return type of the function) has to be given, otherwise data truncation can happen during reading!!!
//Other possibility to give the same variable in the function call for def_type as the variable the return value of the
//function is stored (the variable, we want to read).
//The >> operator is not sensitive to slight type mismatch, so for example double values can be read truncated into integer variables 
//without any warning!!!
template < class InT >
InT ReadThisLine(ifstream &file, int status, InT def_value, const char *var_name, const char *routine_name, const char* file_name="")
{
	//status can be
	//	 0: for optional data, if not present, only warning is generated and default value is set; file pos will be set to the next line
	//   1: for mandatory data, if not present, error message is given, and the program terminates, if present go to next line
	//	 2: for optional data, if not present, only warning is generated and default value is set, file pos set as before reading, if present file pos stays as after reading
	//   3: for mandatory data, if not present, error message is given, and the program terminates, if present do not go to next line
	//	 4: for optional data, if not present, file pos set as before reading, if present file pos stays as after reading, no warning
	//(-)6: for optional data, this handles optional data read after a mandatory from *.dat, trying to avoid problems coming from not
	//      having comment or anything else on the line. Goes to next line ater reading, gives warning for 6, not for -6
	//(-)8: for optional data, same as 6, but do not go to next line
	int result;
	int term_stat=0;//whether to terminate, if there is no more line in the file
	std::streamoff first_pos,current_pos,after_pos;//to position the get pointer in the file
	InT input;
	bool nodata = false;//indicating there was no data at all after the previous reading in the line

	switch (status)
	{
		case 0:
		case 1:
			if (strlen(file_name)!=0)
				term_stat=1;
			break;
	}

	first_pos=file.tellg();//storing the file position
	//has to make sure, that in case of optional data reading after mandatory data reading in a line, if nothing 
	//followed the mandatory data, then the file pointer can go to next line instead of staying in the same
	//it is not clear, when and why it goes to next line, when it stays after data reading
	//even if it stays and nothing follows, that would cause a problem (reading from the next line), which is handled later
	if (first_pos==next_line_pos && (fabs((double)status==6) || (fabs((double)status == 8) )))
		nodata=true;//do not try to read
	else
	{
		result = SkipLine(file, file_name, term_stat);//more data should follow, if filename was given, program will terminate, if there is not enough line

			if (result == 0)
				file.clear(ios::goodbit);

			next_line_pos = file.tellg();//storing the file position for the next line
			file.seekg(first_pos, ios::beg);//resets the position in the file for the reading of the data
			file >> input;

		after_pos = file.tellg();//storing the file position after reading, -1, if could not be read
		if (after_pos == -1)
		{
			file.seekg(first_pos, ios::beg);//setting back to before reading, as could not read
			file.clear(ios::goodbit);
			nodata = true;
		}
		else
		{
			if (after_pos > next_line_pos)//read something, but from the next line, as no comment was given on the targat line
			{
				if (status==8)
					file.seekg(first_pos, ios::beg);//setting to next line
				else
					file.seekg(next_line_pos, ios::beg);//setting to next line
				file.clear(ios::goodbit);
				nodata = true;
			}
		}

		result = SkipLine(file);
		if (result == 0)
			file.clear(ios::goodbit);

		current_pos = file.tellg();//getting the file position after reading to see if match with the beginning of the next line

	}
	if ((file.rdstate() & file.failbit) || nodata)
	{
		if (!((status==4) || (status==-6) || (status == -8)))
		{
			if (status & 1)//status odd
				cout << "\n****ERROR***" << endl;
			else
				cout << "\nWARNING(" << ++warn << "):";
				
			cout<<"Variable "<<var_name<<" could not be read in "<<routine_name<<"!"<<endl;
			cout<<"\tProbably end of line was reached!"<<endl;
		}
		if (status&1)
		{
			cout<<"\tCannot run this way, exiting "<<endl;//" rd "<<file.rdstate()<<endl;
			CleanExit();
		}
		else
		{
			if (!((status == 4) || (status != -6) || (status != -8)))
				cout<<"\tIt will be set to default value "<<def_value<<"!"<<endl;
			input=def_value;
			file.clear(ios::goodbit);
			file.seekg(first_pos,ios::beg);//resets the position in the file to the place before reading, as after faulty read the position can be lost
			
			if (status==0 || fabs((double)status)==6)
				SkipLine(file);//go to next line 
			
		}
		
	}
	else
	{
		if (current_pos != next_line_pos && status < 2)//in case of 0 and 1 data is expected in the next line
			cout << "\nWARNING(" << ++warn << "): Could not go to next line after reading " << var_name << " in routine " << routine_name << "!" << endl;
		if (fabs((double)status) >= 2 && (fabs((double)status!=6)))
			file.seekg(after_pos, ios::beg);//resets the position in the file to the place after reading, as after faulty read the position can be lost
	}
	return(input);
};

//read the binary file and notice, if there is no data
//flag=-1 no warning, return for optional data
//flag=0 only warning and return
//flag=1,2 ERROR, quit
template < class InT >
int ReadBin(ifstream &file, InT *read_value, longint size, const char *var_name, const char *routine_name, const char *file_name, int flag=0 )
{
	file.read(reinterpret_cast<char *>(read_value),size* sizeof(*read_value));
	if ((file.rdstate() & file.failbit))
	{
		file.clear(ios::goodbit);
		if (flag>0)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "The variable " << var_name << " could not be read from file " << file_name << " in routine " << routine_name << endl;
			cout << "as end of file was reached!" << endl;
			if (flag == 2)
				cout << "Exact continuation is not possible! ";
			cout<<"Cannot run this way, exiting..." << endl;
			CleanExit();
		}
		else
		{
			if ( flag==0)
			{
			cout << "\nWARNING(" << ++warn << "): Variable " << var_name << " could not be read in " << routine_name << "!" << endl;
			cout << "\tThe binary file " << file_name << " is either damaged. or was created with different parameters." << endl;
			}
			return 0;
		}
	}
	return 1;
}

template < class InT >
void ResizeArray(int *max, int new_max, InT **array, const char *array_name, const char *routine_name, int multiply = 1)
{
	int i;
	bool isarray = true;
	InT *temp_array, *newp, *oldp;

	if (*array ==  NULL)//was not created
		isarray = false;
	if (isarray)
	{
		temp_array = new InT[*max * multiply];//the original array will be copied here
		if (temp_array == NULL)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Resizing of " << array_name << " array in ResizeArray called from " << routine_name << " routine was not successful," << endl;
			cout << "as temporary array cannot be created, most probably because of the shortage of memory." << endl;
			CleanExit();
		}

		newp = temp_array;
		oldp = *array;
		for (i = 0; i < *max * multiply; i++)
			*newp++ = *oldp++;//copy the original to temp

		delete[] * array;
	}

	*array=new InT[new_max * multiply];//creating the larger new array
	if (*array==NULL)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"Resizing of "<<array_name<<" array in ResizeArray called from "<<routine_name<<" routine was not successful,"<<endl;
		cout<<"as the new, larger array cannot be created, most probably because of the shortage of memory."<<endl; 
		CleanExit();
	}
	if (isarray)
	{
		newp = *array;
		oldp = temp_array;
		for (i = 0; i < *max *multiply; i++)
		{
			*newp++ = *oldp++;//copy the values to the new array
		}

		delete [] temp_array;
	}
	*max=new_max;
}

void NoArray(const char *routine, const char *array_name);//Error message, when no new array could be created

//-------------Allocating memory for static integer arrays----------------------
template < class InT >
void SetArraysize(InT **apointer, int size, const char *array_name, const char *routine_name)
{
	
	*apointer=new InT[size];
	if (*apointer==NULL)
		NoArray(routine_name,array_name);
	return;
}

ostream &J2(ostream &os);//formatting utility
ostream &J4(ostream &os);//formatting utility
ostream &J7(ostream &os);//formatting utility
ostream &J10(ostream &os);//formatting utility
ostream &J12(ostream &os);//formatting utility
ostream &J13(ostream &os);//formatting utility
ostream &J16(ostream &os);//formatting utility
ostream &J20(ostream &os);//formatting utility
ostream &J23(ostream &os);//formatting utility
ostream &J27(ostream &os);//formatting utility
ostream &J30(ostream &os);//formatting utility
ostream &J35(ostream &os);//formatting utility

int RoundNearInt(double number);//rounds number to the nearest integer
double Ran1(longint  &idum);//returns a random number flatly distributed between 0 and 1
//the call is usually Ran1(SEED)


//file handling functions
int CheckReadFileState(ifstream &file, const char *routine_name, const char *file_name="");//checking the status of the file after reading

int CheckFileState(ifstream &file, const char *routine_name, const char *file_name="");//checking the file status, if the file is already open for reading
void CheckFileState(ofstream &file, const char *routine_name, const char *file_name="");//checking the file status, if the file is already open for saving
void OpenFile(ofstream &file, const char *file_name, const char *routine_name, int flag);//checking the file status after opening it for saving
const char *ConvertTextFile(ifstream &file);//converting text file to have the proper line end

//Although itoa works for Visual C++, but it does not recognized by gcc, it was necessary to create a universally applicable
//routine, *conv number array has to be NULL or already created
void IntToStr(char **conv_number, int number);//converts a number's digit into their character representation

//Keyword related stuffs
void WrapFree(const string &input, string &keyword, string &parameters);//Splits input format into keyword and parameters
bool ErrLackValue(string &inputline);//This function complains about lack of parameter in inputline and stops the run

inline string& Ltrim(string& str, const string& chars = "\t\n\v\f\r ")
{
	str.erase(0, str.find_first_not_of(chars));
	return str;
}

inline string& Rtrim(string& str, const string& chars = "\t\n\v\f\r ")
{
	str.erase(str.find_last_not_of(chars) + 1);
	return str;
}

//To handle if there is integer in normal form, float form, +/-NeNNNN or +/-NENNNN form in the stream, and there would be other reading after that
//the >> operator would stop at the point or e, E, and the remainder part would be read into the next variable
//This will decide, whether it is an integer represented as integer, float or proper scientific form, and continue if possible. 
//return with true if reading was ok, false if it failed as there was no data, terminates if unapropriate data was present
inline bool GetInt(stringstream &ss, string &line, int *intvalue, const char *routine_name)
{
	int power;
	char nextchar;
	string nextitem,remainder;
	int i,posE=-100;
	bool good=true,is_real=false;
	
	ss >> nextitem;//extract the next item from the stream
	stringstream mystream(nextitem);
	if (!(mystream >> *intvalue))
		return false;
	
	if (mystream.rdstate() == ios::eofbit)
		return true;
	nextchar = mystream.peek();
	switch (nextchar)
	{
	case '.':
		mystream >> nextchar;
		mystream >> remainder;
		for (i = 0; i < (int)remainder.length(); i++)
		{
			if (remainder[i]!='0')//it is not all 0
			{
				if (remainder[i]=='e' || remainder[i] == 'E')
				{
					if (posE != -100)
					{
						good = false;
						break;//this is not  valid number e appears twice
					}
					posE = i;
				}
				else
				{
					if (remainder[i] > '0' && remainder[i] <= '9')
						is_real = true;//there are non-zero digits after the decimal point
					else
					{
						good = false;
						break;//this cannot be a number at all
					}
				}
			}
		}
		
		//all is 0, can proceed
		break;
	case 'e':
	case 'E':
		posE = -1;//there is no decimal point
		mystream >> nextchar;
		mystream >> remainder;
		if (remainder[0] == '-')
		{
			good = false;
			break;//c
		}
		for (i = 0; i < (int)remainder.length(); i++)
		{
			if (!(remainder[i] >= '0' && remainder[i] <= '9'))
			{	
				good = false;
				break;//this is not  valid number + appears twice
			}
		}
		break;
	default:
			good = false;
	}
	if (is_real && posE == -100)//not scientific only non-zero fraction
		good = false;
	else
	{
		if (posE >= -1)//this might be scientific format
		{
			if (posE == -1) //NNNeNNN or NNNe-NNN format where N is a digit, only integers in format  +/-NeNNNN  accepted, the number of trailing digits can be more
			{
				if (*intvalue > 9 || *intvalue < -9)
					good = false;
				else
				{
					power = atoi(remainder.c_str());
					*intvalue *= (int)pow(10, power);
				}


			}
			else
			{
			//there is a decimal point as well
			int fraction = 0;
			mystream.clear();
			mystream.str(remainder);
			if (posE > 0)
				mystream >> fraction;//this is the fraction part
			string fractionstring = std::to_string(fraction);
			mystream >> nextchar;//extract the e or E
			mystream >> power;
			string finalstring;
			finalstring = std::to_string(*intvalue) + std::to_string(fraction);
			finalstring.append(power - fractionstring.length(), '0');
			mystream.clear();
			mystream.str(finalstring);
			mystream >> *intvalue;
			}
		}
	}
	if (!good)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "In " << routine_name << " during processing the line " << endl;
		cout << line << endl;
		cout << "in the " << datfilename << " file " << nextitem << " was encountered instead of an integer value!" << endl;
		cout << "Only format 1234, 1.234e3, 1.234E3 or 1234.00 (any number of trailing zeros) can be extracted!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	}

	return true;
}


//extracting the satndard chemical symbol from what was given in *.dat
string ExtractChemSymbol(string &symbol);

//Check the tags in the tag_table, to prevent the usage of not listed tags, for developers
//Do not use explicit tag comparision,, always use tags with this function!
inline string CheckTag(string name)
{
	if (tag_table.count(name) > 0)
		return name;
	
	else
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Invalid tag \"" << name << "\" was used in the code!" << endl;
		cout << "Update the tag_table in global.cpp!. Exiting..." << endl;
		CleanExit();
		return "";
	}
	
}

//Check the keys in the key_table to prevent the usage of not listed keys, for developers
//Do not use explicit key comparision, always use keys with this function!
inline string CheckKey(string name)
{
	if (key_table.count(name) > 0)
		return name;

	else
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Invalid key \"" << name << "\" was used in the code!" << endl;
		cout << "Update the key_table in global.cpp!. Exiting..." << endl;
		CleanExit();
		return "";
	}

}

//Check acceptable values in the key_table to prevent the usage of not listed values, for developers
//Do not use explicit values, always use values with this function!
inline string CheckKeyValue(string key, string value)
{
	if (key_table.count(key) == 0)//invalid key
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Invalid key \"" << key << "\" was used in the code!" << endl;
		cout << "Update the key_table in global.cpp!. Exiting..." << endl;
		CleanExit();
		return "";
	}
	else
	{
		if (key_table.count(key) > 0)
		{
			if (key_table.at(key).values.count(value) > 0)
				return value;
		}
		else
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Invalid value \"" << value << "\" was used for KEY: \""<<key<<"\" in the code!" << endl;
			cout << "Check the key_table in global.cpp!. Exiting..." << endl;
			CleanExit();
			return "";
		}
		return "";
	}
	return "";

}

//Check, whether the key is implemented in the key_table at all, and for what tag it is acceptable to indicate the usage of not listed keys
//mode is only needed, if for a tag there are keys, which cannot be used in all cases
inline bool CheckKeyTag(string key, string tag,string mode="")
{
	if (key_table.count(key) == 0)//invalid key
	{
		cout << "\nERROR: Key /"<<key<<"/ was used in the "<<datfilename<<" ***"<<tag<<"*** section!" << endl;
		return false;
	}
	else
	{
		//it is in the key table
		
			
		int i = 0;
		string valid_tags;
		bool is_intag = false;

		for (auto&& it : key_table.at(key).tag_names)
		{
			if (i > 0)
				valid_tags+=", ";
			valid_tags += "***";
			valid_tags += it;
			valid_tags += "***";
			if (it == tag)
				is_intag = true;//the key is in the NB section
			++i;
			
		}
		if (is_intag)
		{
			cout << "\nERROR: Key /" << key << "/ cannot be used in the [ " << tag << " ] section";
			if (!mode.empty())
				cout << " in case of " << mode << "!";
			else
				cout << endl;
			
		}
		else
		{
			cout << "\nERROR: Key /" << key << "/ can only be used for tag(s) " << valid_tags;
			cout << " and not for tag ***" << tag << "*** !" << endl;
		}
		return true;
	}
	return "";

}


//give back the acceptable values for a key in the key_table. This is to prevent the usage of not listed values, for developers
//Do not use explicit values, always use values with this function!
inline string GiveKeyValues(string key)
{
	string mystring;
	if (key_table.count(key) == 0)//invalid key
	{
		cout << "\nERROR: Invalid key \"" << key << "\" was used in the code!" << endl;
		cout << "Update the key_table in global.cpp!. Exiting..." << endl;
		CleanExit();
		return "";
	}
	else
	{
		if (key_table.count(key) > 0)
		{
			if (key_table.at(key).values.size() > 0)
			{
				int i = 0;

				//std::set<unsigned long>::iterator it;
				for (auto&& it : key_table.at(key).values)
				{
					if (i > 0)
						mystring += ", ";
					mystring += (it);
					++i;
				}

				return mystring;
			}
		}
		return "";
	}
}
//free format input help related functions
inline void GiveTagTable()
{
	cout << "Free format data input\nImplemented tags:" << endl;
	for (auto&& it : tag_table)
	{
		cout << "Tag name: " << J12<<it.first <<J16<< "  Description:  " << it.second.description <<endl;
		cout<< "\t\t\tMultiplicity: " << it.second.multiplicity << endl;
		
	}
		
}

// gives information about the specific key
inline void GiveKeyHelp(string key="")
{
	string mystring;
	if (key != "")
	{
		//std::transform(key.begin(), key.end(), key.begin(), ::toupper);
		if (key_table.count(key) == 0)//invalid key
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "No key \"" << key << "\" is implemented in the program!" << endl;
			CleanExit();
	
		}
		else
		{
			if (key_table.count(key) > 0)
			{
				cout << key << endl;
				cout << "\tDescription: " << key_table.at(key).description << endl;
				if ((key_table.at(key).alternative_key.size()) > 0 )
				{
					
					int i = 0;
					bool isvalue = false;
					//std::set<unsigned long>::iterator it;
					for (auto&& it : key_table.at(key).alternative_key)
					{
						if (it != "")
						{
							isvalue = true;
							if (i > 0)
								cout << ", ";
							else
								cout << "\tAlternatives: ";
							cout << (it);
						}
						++i;
					}
					if (isvalue)
						cout << endl;
				}
				cout << "\tStatus: " << key_table.at(key).status << endl;
				if (key_table.at(key).tag_names.size() > 0)
				{
					cout << "\tTag name(s) where it can be used: ";
					int i = 0;
					//std::set<unsigned long>::iterator it;
					for (auto&& it : key_table.at(key).tag_names)
					{
						if (i > 0)
							cout << ", ";
						cout << (it);
						++i;
					}
					cout << endl;
				}
				if (key_table.at(key).values.size() > 0)
				{
					cout << "\tAcceptable values: ";
					int i = 0;
					//std::set<unsigned long>::iterator it;
					for (auto&& it : key_table.at(key).values)
					{
						if (i > 0)
							cout << ", ";
						cout << (it);
						++i;
					}
					cout << endl;
				}

			//	std::transform(key_table.at(key).default_value.begin(), key_table.at(key).default_value.end(), mystring, ::toupper);
				if (key_table.at(key).default_value != "")
					cout << "\tDefault value: " << key_table.at(key).default_value << endl;
			}
		
		}
	}
}
//writing the key_table
//status=0 whole
//status =1 mandatory keys
//status=2 optional keys
//if tag given, all the keywords for that tag regardless the status 
inline void GiveKeyTable(int status, string tag="")
{
	if (tag != "")
	{
		status = 0;
		//std::transform(tag.begin(), tag.end(), tag.begin(), ::toupper);
		if (tag_table.count(tag) < 1)
		{
			cout << "No tag " << tag << " is implemented in the program!" << endl;
			CleanExit();
		}
	}
	cout << "Free format data input" << endl;
	switch (status)
	{
		case 1: cout << "Implemented mandatory keywords : " << endl; break;
		case 2: cout << "Implemented optional keywords : " << endl; break;
		default: 
			cout << "Implemented keywords";
			if (tag != "")
				cout<<" for tag " << tag;
			cout<<": " << endl;
	}
	string mystatus;
	
	for (auto&& it : key_table)
	{
		if (status == 0)
		{
			if (tag != "")
			{
				for (auto&& it2 : key_table.at(it.first).tag_names)
				{
					if (it2 == tag)
						GiveKeyHelp(it.first);
				}
			}
			else
				GiveKeyHelp(it.first);
		}
		else
		{
			mystatus = it.second.status;
			//std::transform(mystatus.begin(), mystatus.end(), mystatus.begin(), ::toupper);
			if (status == 1)
			{
				if (mystatus.find("MANDATORY") != string::npos)
					GiveKeyHelp(it.first);
						
			}
			else
			{
				if (status == 2)
				{
					if (mystatus.find("OPTIONAL") != string::npos)
						GiveKeyHelp(it.first);
				}
				else
					GiveKeyHelp(it.first);
			}
		}
		
	}

}

//gives the informatiom about a specific symbol in the Waasmaier xray coeff calculation table
//mode=0 no header
//mode=1 with header
inline void GiveWaasmaierPar(int mode, string symbol = "")
{
	int i;
	string mystring;
	if (symbol != "")
	{
		if (Xray_coeffs_W.count(symbol) == 0)//invalid symbol
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "No symbol \"" << symbol << "\" is found in the extended Waasmaier-Kirfel Xray coeff calculation table!" << endl;
			cout << "Cannot calculate the X-ray coefficients, exiting..." << endl;
			CleanExit();

		}
		else
		{
			if (Xray_coeffs_W.count(symbol) > 0)
			{
				if (mode == 1)
				{
					cout << "The  according to the extended Waasmaier-Kirfel  for the atomic scattering factor calculations" << endl;
					cout << "Acta Cryst. (1995). A51,416-431" << endl;
				}
				cout <<"\nAtomic number: "<< Xray_coeffs_W.at(symbol).atomic_number<<" "<<symbol << "  "<< Xray_coeffs_W.at(symbol).description<<endl;
				for (i = 0; i < 5; i++)
					cout << "\tParameters a[" << i + 1 << "]: " << Xray_coeffs_W.at(symbol).a[i] << endl;
				for (i = 0; i < 5; i++)
					cout << "\tParameters b[" << i + 1 << "]: " << Xray_coeffs_W.at(symbol).b[i] << endl;
				cout << "\tParameters c: " << Xray_coeffs_W.at(symbol).c << endl;
			}
		}
	}
}

//give the list of chemical symbols and other info for Xray scattering factor calculation based on the Waasmaier-Kirfel table
//mode=0 all records (no par)
//mode=1 only neutral (atoms, valence state, radical) (no par)
//mode=2 ions (no par)
//mode=3 all record beginning with par
//mode=4 all record for atomic number represented by par
inline void GiveWSymbols(int mode, string par = "")
{
	int i,N;
	stringstream ss(par);
	if (par != "")
	{
		switch (mode)
		{
		case 3://gives the entries beginning with par
			i = 1;//give header first time
			for (auto &&it : (Xray_coeffs_W))
			{
				if (it.first.rfind(par, 0) == 0)
				{
					GiveWaasmaierPar(i, it.first);
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic symbol " << par << endl;
				CleanExit();
			}
			break;
		case 4://find all entries with N atomic number
			if (!(ss >> N))
			{
				cout << par << " cannot be an atomic number, exiting..." << endl;
				CleanExit();
			}
			i = 1;//give header first time
			for (auto &&it : (Xray_coeffs_W))
			{
				if (Xray_coeffs_W.at(it.first).atomic_number == N)
				{
					GiveWaasmaierPar(i, it.first);
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic number " << N << endl;
				CleanExit();
			}
			break;
		default:
			cout << "\n*****ERROR*****" << endl;
			cout << "This combination is not implemented" << endl;
			CleanExit();
		}

	}
	else
	{
		cout << "The program can calculate the X-ray scattering factors according to the extended Waasmaier-Kirfel tables" << endl;
		cout << "(Acta Cryst. (1995). A51,416-431), and from these the X-ray coefficients as well." << endl;
		switch (mode)
		{
		case 0:
			cout << "The following symbols are implemented:" << endl;
			break;
		case 1:
			cout << "The following atomic symbols (radicals, valence states) are implemented:" << endl;
			break;
		case 2:
			cout << "The following ions  are implemented:" << endl;
			break;

		}

		for (auto &&it : (Xray_coeffs_W))
		{
			switch (mode)
			{
			case 0:
				cout << it.first << endl;
				break;
			case 1:
				if (Xray_coeffs_W.at(it.first).description != "ion")
					cout << it.first << endl;
				break;
			case 2:
				if (Xray_coeffs_W.at(it.first).description == "ion")
					cout << it.first << endl;
				break;
			case 3:
				break;
			}

		}
	}
}


//gives the information about a specific symbol in Compton scattering calculation table
//mode=0 no header
//mode=1 with header
inline void GiveComptonPar(int mode, string symbol = "")
{
	int i;
	string mystring;
	if (symbol != "")
	{
		if (Compton.count(symbol) == 0)//invalid symbol
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "No symbol \"" << symbol << "\" is found in the Compton scattering calculation table!" << endl;
			cout << "Cannot calculate it, exiting..." << endl;
			CleanExit();

		}
		else
		{
			if (Compton.count(symbol) > 0)
			{
				if (mode == 1)
				{
					cout << "The Balyuzi parameters for the Compton scattering calculations" << endl;
					cout << "Acta Cryst. (1975). A31, 600 " << endl;
				}
				cout << "\nAtomic number: " << Compton.at(symbol).atomic_number << " " << symbol << endl;
				for (i = 0; i < 5; i++)
					cout << "\tParameters a[" << i + 1 << "]: " << Compton.at(symbol).a[i] << endl;
				for (i = 0; i < 5; i++)
					cout << "\tParameters b[" << i + 1 << "]: " << Compton.at(symbol).b[i] << endl;
				
			}
		}
	}
}

//give the list of chemical symbols and other info for Compton scattering based on the Balyuzi table
//mode=0 all records (no par)
//mode=1 all record beginning with par
//mode=2 all record for atomic number represented by par
inline void GiveCSymbols(int mode, string par = "")
{
	int i, N;
	stringstream ss(par);
	if (par != "")
	{
		switch (mode)
		{
		case 1://gives the entries beginning with par
			i = 1;//give header first time
			for (auto &&it : (Compton))
			{
				if (it.first.rfind(par, 0) == 0)
				{
					GiveComptonPar(i, it.first);
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic symbol " << par << endl;
				CleanExit();
			}
			break;
		case 2://find all entries with N atomic number
			if (!(ss >> N))
			{
				cout << "\n*****ERROR*****" << endl;
				cout << par << " cannot be an atomic number, exiting..." << endl;
				CleanExit();
			}
			i = 1;//give header first time
			for (auto &&it : (Compton))
			{
				if (Compton.at(it.first).atomic_number == N)
				{
					GiveComptonPar(i, it.first);
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic number " << N << endl;
				CleanExit();
			}
			break;
		default:
			cout << "\n*****ERROR*****" << endl;
			cout << "This combination is not inplemented" << endl;
			CleanExit();
		}

	}
	else
	{
		cout << "The program can calculate the Compton scattering term according to Balyuzi" << endl;
		cout << "Acta Cryst. (1975). A31, 600" << endl;
		cout << "The following symbols are implemented:" << endl;
				
		for (auto &&it : (Compton))
			cout << it.first << endl;
			
	}
}

//gives the information about a specific symbol in coherent neutron scattering length table
//mode=0 no header
//mode=1 with header
inline void GiveNscatlengthPar(int mode, string symbol = "")
{
	string mystring;
	if (symbol != "")
	{
		if (N_scat_length.count(symbol) == 0)//invalid symbol
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "No symbol \"" << symbol << "\" is found in the coherent neutron scattatering length table!" << endl;
			cout << "Cannot proceed, exiting..." << endl;
			CleanExit();

		}
		else
		{
			if (N_scat_length.count(symbol) > 0)
			{
				if (mode == 1)
				{
					cout << "The coherent neutron scattering length parameters is based on Sears," << endl;
					cout << "Neutron News (1992) 3, Nr. 3, 26-37 " << endl;
				}
				cout << "\nAtomic number: " << N_scat_length.at(symbol).atomic_number << " " << symbol << endl;
				cout << "\tb_coh real part: " << N_scat_length.at(symbol).b_coh.real() << endl;
				

			}
		}
	}
}

//give the list of chemical symbols and other info for coherent neutron scattering length based on the Sears table
//mode=0 all records (no par) only key
//mode=1 all record beginning with par
//mode=2 all record for atomic number represented by par
//mode=3 all record (no par) Symbol and key
inline void GiveNSymbols(int mode, string par = "")
{
	int i, N;
	stringstream ss(par);
	if (par != "")
	{
		switch (mode)
		{
		case 1://gives the entries beginning with par
			i = 1;//give header first time
			for (auto &&it : (N_scat_length))
			{
				if (it.first.rfind(par, 0) == 0)
				{
					GiveNscatlengthPar(i, it.first);
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic symbol " << par << endl;
				CleanExit();
			}
			break;
		case 2://find all entries with N atomic number
			if (!(ss >> N))
			{
				cout << "\n*****ERROR*****" << endl;
				cout << par << " cannot be an atomic number, exiting..." << endl;
				CleanExit();
			}
			i = 1;//give header first time
			for (auto &&it : (N_scat_length))
			{
				if (N_scat_length.at(it.first).atomic_number == N)
				{
					GiveNscatlengthPar(i, it.first);
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic number " << N << endl;
				CleanExit();
			}
			break;
		default:
			cout << "\n*****ERROR*****" << endl;
			cout << "This combination is not implemented" << endl;
			CleanExit();
		}

	}
	else
	{
		cout << "The program can calculate the neutron scattering coefficients based on the" << endl;
		cout << "coherent neutron scattering lengthes published by Sears" << endl;
		cout << "Neutron News (1992) 3, Nr. 3, 26-37 " << endl;
		cout << "The following symbols are implemented:" << endl;

		if (mode==3)
			for (auto &&it : (N_scat_length))
				cout << it.second.symbol<<" for isotope "<<it.first << endl;
		else
			for (auto &&it : (N_scat_length))
				cout <<it.first << endl;

	}
}

#ifdef _ADVANCED_GEOM_CONST
//gives the information about a specific symbol (key) in BV param table
//mode=0 no header
//mode=1 with header
inline void GiveBVPar(int mode, string symbol = "", string symbol_rev="")
{
	string mystring;
	if (symbol != "" && symbol_rev!="")
	{
		if (BV_param.count(symbol) == 0 && BV_param.count(symbol_rev)==0)//invalid symbol
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "No parameters for either \"" << symbol << "\" or \"" << symbol_rev << "\" is found in the bond valence parameters table!" << endl;
			cout << "Cannot proceed, exiting..." << endl;
			CleanExit();

		}
		else
		{
			if (BV_param.count(symbol) > 0)
				mystring.append(symbol);
			else
				mystring.append(symbol_rev);
			if (mode == 1)
			{
				cout << "The bond valence R0 and b parameters are based on I. David Brown" << endl;
				cout << "Brockhouse Institute for Materials Research, McMaster University, Hamilton, Ontario Canada." << endl;
				cout<<"idbrown@mcmaster.ca" << endl;
			}
			cout << "\n"<<mystring<<" Atomic numbers: ";
			if (BV_param.count(symbol) > 0) 
				cout<< BV_param.at(mystring).central_atomic_number << " " << BV_param.at(mystring).neigh_atomic_number << endl;
			else
				cout << BV_param.at(mystring).neigh_atomic_number << " " << BV_param.at(mystring).central_atomic_number << endl;
			cout << "\tR0: " << BV_param.at(mystring).R0<<" A b: "<< BV_param.at(symbol).b << " A"<<endl;
		}
	}
}

//give the list of chemical symbols and other info for bond valence R0 and b parameters based on the I. David Brown table
//mode=0 all records (no par) only oxidation state names
//mode=1 all records beginning with par
//mode=2 all record for atomic number represented by par
//mode=3 all record for par par2 pairs oxidation state names 
inline void GiveBVSymbols(int mode, string par = "",string par2 ="")
{
	int i, N;
	stringstream ss(par);
	if (par != "")
	{
		switch (mode)
		{
		case 1://gives the entries beginning with par
			i = 1;//give header first time
			for (auto &&it : (BV_param))
			{
				unsigned int pos = it.first.rfind(" ");
				if (it.first.rfind(par, 0) == 0)
				{
					GiveBVPar(i, it.first,"no");//no need for reverse, it was found
					i = 0;//not to give header again
				}
				else
				{
					if (it.first.rfind(par, pos+1) == pos+1)
					{
						GiveBVPar(i, it.first, "no");//no need for reverse, it was found
						i = 0;//not to give header again
					}
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found for " << par << endl;
				CleanExit();
			}
			break;
		case 2://find all entries with N atomic number
			if (!(ss >> N))
			{
				cout << "\n*****ERROR*****" << endl;
				cout << par << " cannot be an atomic number, exiting..." << endl;
				CleanExit();
			}
			i = 1;//give header first time
			for (auto &&it : (BV_param))
			{
				if (BV_param.at(it.first).central_atomic_number == N)
				{
					GiveBVPar(i, it.first,"no");
					i = 0;//not to give header again
				}
				else
				{
					if (BV_param.at(it.first).neigh_atomic_number == N)
					{
						GiveBVPar(i, it.first,"no");
						i = 0;//not to give header again
					}
				}

			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found with atomic number " << N << endl;
				CleanExit();
			}
			break;
		case 3://gives the containing the par par2 pairs
			if (par2 == "")
			{
				cout << "No second symbol was given for GiveBVSymbols with mode=3" << endl;
				cout << "Exiting..." << endl;
				CleanExit();
			}
			i = 1;//give header first time
			for (auto &&it : (BV_param))
			{
				if ((it.second.central_standard_symbol.compare(par) == 0 && it.second.neigh_standard_symbol.compare(par2)==0) || (it.second.central_standard_symbol.compare(par2) == 0 && it.second.neigh_standard_symbol.compare(par) == 0))
				{
					GiveBVPar(i, it.first, "no");//no need for reverse, it was found
					i = 0;//not to give header again
				}
			}
			if (i == 1)//no record found at all
			{
				cout << "\n*****ERROR*****" << endl;
				cout << "No entry was found for " << par <<" and "<<par2<< endl;
				CleanExit();
			}
			break;
		default:
			cout << "\n*****ERROR*****" << endl;
			cout << "This combination is not implemented" << endl;
			CleanExit();
		}

	}
	else
	{
		cout << "The program uses parameters from the table of I. David Brown" << endl;
		cout << "Brockhouse Institute for Materials Research, McMaster University, Hamilton, Ontario Canada." << endl;
		cout << "idbrown@mcmaster.ca" << endl; 
		cout << "for the default values of the R0 and b parameters for bond valence sum calculation" << endl;
		cout << "The parameters for the following pairs can be found:" << endl;

		/*if (mode == 3)
			for (auto &&it : (BV_param))
				cout << it.second.symbol << " for isotope " << it.first << endl;
		else*/
			for (auto &&it : (BV_param))
				cout << it.first << endl;

	}
}
#endif
inline void Help()
{
	cout << "The the program can be started the following ways:" << endl;
	cout << "\tEXECUTABLE\t\t\t(the file name will be asked)" << endl;
	cout << "\tEXECUTABLE filename" << endl;
	cout << "\tEXECUTABLE filename cont" << endl;
	cout << "\tEXECUTABLE filename y" << endl;
	cout << "\tEXECUTABLE filename cont y" << endl;
	cout << "\twhere cont means to continue the run, y stand for any single character to show, that the program can exit silently" << endl;

	cout << "\nGeneral help, how the program can be started" << endl;
	cout << "\tEXECUTABLE -help" << endl;

	cout << "\nHelp related to the free format data input" << endl;
	cout << "\tEXECUTABLE -helptags\n\t\tInformation about the tags" << endl;
	cout << "\tEXECUTABLE -helpkeys\n\t\tInformation about the keywords" << endl;
	cout << "\tEXECUTABLE -helpkeys tag\n\t\tInformation about the keywords of the specified tag" << endl;
	cout << "\tEXECUTABLE -helpkey key\n\t\tInformation about the specified key" << endl;
	cout << "\tEXECUTABLE -helpkeyman\n\t\tInformation about the mandatory keywords" << endl;
	cout << "\tEXECUTABLE -helpkeyop\n\t\tInformation about the optional keywords" << endl;

	cout << "\nHelp related to the Xray coefficient calculation" << endl;
	cout << "\tEXECUTABLE -helpXcoeffs\n\t\tListing all the available symbols implemented in the extended Waasmaier-Kirfel table" << endl;
	cout << "\tEXECUTABLE -helpXcoeffsat\n\t\tListing the available neutral entries implemented in the extended Waasmaier-Kirfel table" << endl;
	cout << "\tEXECUTABLE -helpXcoeffsion\n\t\tListing the available ions implemented in the extended Waasmaier-Kirfel table" << endl;
	cout << "\tEXECUTABLE -helpXcoeff symbol\n\t\tListing the parameters for the atomic symbol 'symbol' in the extended Waasmaier-Kirfel table" << endl;
	cout << "\tEXECUTABLE -helpXcoeffsS symbol\n\t\tListing the parameters for the entries beginning with 'symbol' in the extended Waasmaier-Kirfel table" << endl;
	cout << "\tEXECUTABLE -helpXcoeffsZ Z\n\t\tListing the parameters for the entries with atomic number 'Z' in the extended Waasmaier-Kirfel table" << endl;

	cout << "\nHelp related to the Xray Compton scattering calculation" << endl;
	cout << "\tEXECUTABLE -helpCompton\n\t\tListing the parameters in the extended Balyuzi table" << endl;
	cout << "\tEXECUTABLE -helpComptonS symbol\n\t\tListing the parameters for the entries with atomic symbol 'S' in the extended Balyuzi table" << endl;
	cout << "\tEXECUTABLE -helpComptonZ Z\n\t\tListing the parameters for the entries with atomic number 'Z' in the extended Balyuzi table" << endl;

	cout << "\nHelp related to the neutron scattering length table" << endl;
	cout << "\tEXECUTABLE -helpNcatlenth\n\t\tListing the b_coh in the Sears neutron scattering parameters table" << endl;
	cout << "\tEXECUTABLE -helpNcatlenthS symbol\n\t\tListing the b_coh for the entries with atomic symbol 'S' in the Sears neutron scattering parameters table" << endl;
	cout << "\tEXECUTABLE -helpNcatlenthZ Z\n\t\tListing the b_coh for the entries with atomic number 'Z' in the Sears neutron scattering parameters table" << endl;
#ifdef _ADVANCED_GEOM_CONST
	cout << "\nHelp related to the bond valence parameters table" << endl;
	cout << "\tEXECUTABLE -helpBVparam\n\t\tListing the R0 and b in the I. David Brown bond valence parameters table" << endl;
	cout << "\tEXECUTABLE -helpBVparamS symbol\n\t\tListing the R0 and b for the entries with atomic symbol 'S' in the I. David Brown bond valence parameters table" << endl;
	cout << "\tEXECUTABLE -helpBVparamZ Z\n\t\tListing the R0 and b for the entries with atomic number 'Z' in the I. David Brown bond valence parameters table" << endl;
#endif
	CleanExit();
}
//add an unique value to the array and resize it
inline void AddToArray(int newind, int* count, int** array, const char* array_name, const char* routine_name)
{
	int i;
	bool isarray = true;
	int* temp_array, * newp, * oldp;
	newp = *array;
	for ( i = 0;i < *count;i++)
	{
		if (*newp == newind )
			return;//no double record
		newp++;

	}

	if (*array == NULL)//was not created
		isarray = false;
	if (isarray)
	{
		temp_array = new int [*count];//the original array will be copied here
		if (temp_array == NULL)
		{
			cout << "\n*****ERROR*****" << endl;
			cout << "Resizing of " << array_name << " array in ResizeArray called from " << routine_name << " routine was not successful," << endl;
			cout << "as temporary array cannot be created, most probably because of the shortage of memory." << endl;
			CleanExit();
		}

		newp = temp_array;
		oldp = *array;
		for (i = 0;i < *count; i++)//copy the original to temp
		{
			*newp = *oldp;
			newp++;
			oldp++;
		}

		delete[] *array;
	}

	*array = new int[*count + 1];//creating the larger new array

	if (*array == NULL)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "Resizing of " << array_name << " array in ResizeArray called from " << routine_name << " routine was not successful," << endl;
		cout << "as the new, larger array cannot be created, most probably because of the shortage of memory." << endl;
		CleanExit();
	}
	newp = *array;
	if (isarray)
	{
		oldp = temp_array;
		for (i = 0; i < *count; i++)//copy the values to the new array
		{
			*newp = *oldp;
			newp++;
			oldp++;
		}

		delete[] temp_array;
	}
	*count += 1;
	*newp =  newind;
};

inline void Check0_1(int value, const char *name, const char *routine_name)
{
	if (!(value == 0 || value == 1))
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "For " << name <<  " input value " << value << " was read in " << routine_name << " from file "<<datfilename<<"!" << endl;
		cout << "Only 0 and 1 can be given, cannot run this way, exiting..." << endl;
		CleanExit();
	}
};

inline void CheckCommandlineParams(string name, int npar)
{
	command_arg_type it;
	if (command_arg.count(name) > 0)//in the table, check, if the number of parameters are sufficient
	{
		it = command_arg.at(name);
		for (unsigned int j = 0; j < it.myparams.size(); j++)
		{
			
			if ((int)it.myparams[j].param_name.size() == npar)//there is enough parameters for this argument
			{
				return;
			}
		}
		//there is not enough command line parameters
		cout << "\n*****ERROR*****" << endl;
		cout << "The number of parameters for command line argument " << name << " is not adequate!"<<endl;
		cout << "The right usage of this command line argument: ";
		cout << name << " ";
		if (it.myparams.size() == 1)
		{
			for (unsigned int k = 0; k < it.myparams[0].param_name.size(); k++)
				cout << " " << it.myparams[0].param_name[k];
			cout << endl;
			cout << "\tDescription: " << it.myparams[0].description << endl;
		}
		else
		{
			cout << endl;
			cout << "There are " << it.myparams.size() << " variations:" << endl;
			for (unsigned int j = 0; j < it.myparams.size(); j++)
			{
				cout << "  " << name;
				cout << "\t";
				for (unsigned int k = 0; k < it.myparams[j].param_name.size(); k++)
					cout << " " << it.myparams[j].param_name[k];
				cout << endl;
				cout << "\tDescription:  " << it.myparams[j].description << endl;
				
			}
			
		}
		CleanExit();
		/*for (auto &&it : command_arg)
		{
			cout << it.first << " ";
			if (it.second.myparams.size() == 1)
			{
				for (unsigned int k = 0; k < it.second.myparams[0].param_name.size(); k++)
					cout << " " << it.second.myparams[0].param_name[k];
				cout << endl;
				cout << "\tDescription: " << it.second.myparams[0].description << endl;
			}
			else
			{
				cout << endl;
				cout << "  There are " << it.second.myparams.size() << " variations:" << endl;
				for (unsigned int j = 0; j < it.second.myparams.size(); j++)
				{
					cout << "\t" << it.first;
					cout << "\t";
					for (unsigned int k = 0; k < it.second.myparams[j].param_name.size(); k++)
						cout << " " << it.second.myparams[j].param_name[k];
					cout << endl;
					cout << "\tDescription:  " << it.second.myparams[j].description << endl;
				}
			}

		}*/
	}
	
}

inline void ErrBadSymbol(string symbol)
{
	cout << "\n*****ERROR*****" << endl;
	cout << "The " << symbol << " given for the " << CheckKey("CHEMICAL-SYMBOLS") << " in " << datfilename << " file is not a proper chemical symbol!" << endl;
	cout << "Use proper chemical symbol in the format first upper case, if there is second character lower case," << endl;
	cout<<"for ions give the number, even if it 1!" << endl;
	cout << "The symbols given in " << CheckKey("CHEMICAL-SYMBOLS") << " are mainly for the X-ray coefficient calculation," << endl;
	cout << "so they should match the proper state of the element, which should be used for this system! But they are used" << endl;
	cout<<"as well for Compton calculation, if it is needed, where in case of ions the values for the corresponding atoms are used."<<endl;
	cout << "In case of neutron coefficient calculation the corresponding element symbols are extracted and used for natural abundance," << endl;
	cout<<"other isotope compositions should be specified with the " << CheckKey("ISOTOPE-COUNT_SYMBOLS_RATIOS") << endl;
	cout << "key word in the " << CheckTag("EXP") << " section of the given data set!" << endl;
	cout << "Examples: He, Li1+, (special cases: dummy; CV or SiV for valence state carbon or silicon)" << endl;
	GiveWSymbols(0);
	CleanExit();

}