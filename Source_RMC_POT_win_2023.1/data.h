//header data.h
//Last changed 15.12.2022

//for tabulated data concerning the X-ray coefficient and Compton scattering calculation and sfactorcube
#include "altern.h"
#include <set>
#include <unordered_set>
#include <map>
#include <string>
#include <complex>
#include <vector>


using std::ios;
using std::string; 
using std::complex;
using namespace std::complex_literals;

struct Xray_coeffs_W_type {
	double a[5];
	double b[5];
	double c;
	string description;
	int atomic_number;
	string standard_symbol;
};
struct Compton_type {
	double a[5];
	double b[5];
	int atomic_number;
};
struct N_scat_length_type
{
	int atomic_number;
	string symbol;
	complex<double> b_coh;

};
struct isotope_ratio_type {
	string symbol;
	double ratio;
};
struct isotope_count_type {
	int count;
	isotope_ratio_type *isotopes;
};
#ifdef _ADVANCED_GEOM_CONST
struct BV_param_type{
	int central_atomic_number;
	int neigh_atomic_number;
	double R0;//in Angstrom
	double b;//in Angstrom
	string central_standard_symbol;
	string neigh_standard_symbol;
	string citation;
	string comment;
	string multiple_flag;
};
#endif
//for checking the commad line arguments
//exename command_argument param1 param2...
//variation means same command line argument, different nunber of param, for example variation 2:
//- EXECUTABLE -helpkeys
//- EXECUTABLE -helpkeys tag
// or variation 4
//- EXECUTABLE filename 
//- EXECUTABLE filename cont
//- EXECUTABLE filename y
//- EXECUTABLE filename cont y
struct param_type {  
	std::vector<string> param_name;//if there is two, here goes cont y for example, number of element vector.size()
	string description;
};

struct command_arg_type {
	string group;
	std::vector<param_type> myparams;
};



//not to rewrite the handling of the surface factor, the x and y values should be kept in separate arrays
extern double surfacefact_x[SFACTOR_SIZE_DEF];
extern double surfacefact_y[SFACTOR_SIZE_DEF];
extern const  std::map<std::string, Xray_coeffs_W_type> Xray_coeffs_W;
extern const  std::map<std::string, Compton_type> Compton;
extern const  std::map<std::string, N_scat_length_type> N_scat_length; 
#ifdef _ADVANCED_GEOM_CONST
extern const  std::map<std::string, BV_param_type> BV_param;
#endif
extern const  std::map<std::string, command_arg_type> command_arg;









