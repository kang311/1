//source DataMat.cpp
//Last changed 24.01.2023

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include"classes2.h"

//Defining the static members
int   DataMat::ngr;//number of g(r) data sets
int   DataMat::nsq;//number of S(Q) data sets
int   DataMat::nfq;//number of F(Q) data sets
int   DataMat::nfg;//number of F(g) data sets
int   DataMat::npartials;//number of partials
int   DataMat::ntypes;//number of atom types
int   DataMat::ndiffbin;//number of different bin sizes
int   DataMat::calcmode = 0;//whether to calculate from the local hist, 0 is default
int   DataMat::sfactor_size;//number of points in the surface factor
int  *DataMat::nbins;//number of bins in the g(r) histograms
int  *DataMat::mybins;//the number of used bins
int	 *DataMat::assign_hist;//assigning the different histogram bin sizes to the data sets dim [ntot_datasets]

#ifdef _VIBR_AMP
	int DataMat::nbins_used;//number of bins used in the g(r) histograms
#endif

int *DataMat::pnatoms;//pointer to number of atoms per type

double DataMat::boxedge;//size of the simulation cell in Angstrom
double DataMat::dr;//r spacing in Angstrom, value can change  for the different sets 
double DataMat::volume;//sample volume
double DataMat::rho;//total number density
double *DataMat::xmin;//pointer to RunParams::xmin, minimum value for the histogram in reduced units, can differ for aech different bin size dim [ndiffbin]
double *DataMat::xmax;//pointer to RunParams::xmax, maximum value for the histogram in reduced units for each different bin size dim [ndiffbin]

#ifdef _LOCAL_INV
	double *DataMat::dV_loc;//the histogram bin volumes has to be stored in case of locale invariance calculation
#endif

#ifdef _NO_PERIODIC
	double  DataMat::R0;//the radius of the speric sample inside the box
	double *DataMat::dV_in;//the actual volume of the ribbon inside sphere A for all the possible origins of spfere B and C
	double *DataMat::dV;//the volume of the normal volume element, needed for exafs
#endif

DataMat::DataMat(ExptsData &edata)//constructor
{
	int i,j, m, n, size, inc, k, nmax, ntot_bins,myused_bins;
	int is_hole;//boolean indicator, that there is a hole is the bin->dr assignment
	int message_flag;//boolean indicator to write warning
	int left_bin, right_bin, last_assigned, first_flag;//for the bin->r conversion
	int *gr_flag;//for the bin->r conversion
	double r, rmin, r1, r2, fourpi, tol1, tol2, max;
	double dx=0;
	double *ptemp, *pmatrix, *pqvalue, *pgvalue, *prvalue;
	double rinf, rsup;
	double *rtable = NULL, *sftable = NULL;


	ifstream file;

	//to avoid compiler warning
	//dx=0;
	ptemp = 0;

	if (ntypes == 0)
	{
		cout << "\n" << "*****ERROR*****" << endl;
		cout << "DataMat constructor: Number of types is 0!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	};
	if (nbins[0] <= 0)
	{
		cout << "\n" << "*****ERROR*****" << endl;
		cout << "DataMat constructor: Number of bins was not calculated!" << endl;
		cout << "Cannot run this way, exiting..." << endl;
		CleanExit();
	};
	grsize = edata.grused;//array of data sets sizes i.e. number of points
	sqsize = edata.sqused;//used in each data set
	fqsize = edata.fqused;
	fgsize = edata.fgused;

	//Constants
	fourpi = 4 * PI;//4*pi

	//decide, whether to load the sfactorcube
	max = 0;
	for (i = 0; i < RunParams::ndiffbin; i++)
	{
		if (xmax[i] > max)
			max = xmax[i];
	}
	i = 0;
	if (max > 1.41421356)
		i = 1;
#ifdef _USE_LOCAL_INV
	if (HistoSet::max_loc_r > 1.41421356)
		i = 1;
#endif
	if (i == 1)
	{
		//data for the tabulated values of the surface factor
		if (RunParams::use_custom_sftable)//should be read from file, custom values, should be equdistant
		{
			SafeOpenTextFile(file, sffilename);
			if (CheckFileState(file, "DataMat::DataMat constructor", sffilename) == 0)
			{
				cout << "\nWARNING(" << ++warn << "): Built-in surface factors will be used!" << endl;
				RunParams::use_custom_sftable = false;
			}
		}
		if (RunParams::use_custom_sftable)//the file already open
		{
			sfactor_size=ReadThisLine(file, 1, 1, "surface factor size", "DataMat::DataMat");
			SetArraysize(&rtable, sfactor_size, "rtable", "DataMat::DataMat");// r values for the discretized s factor values
			SetArraysize(&sftable, sfactor_size, "sftable", "DataMat::DataMat");//sfactor values
			
			prvalue = rtable;
			ptemp = sftable;
			SkipLine(file, sffilename, 1);//skipping the first lines used for comment
			
			for (i = 0; i < sfactor_size; i++)//filling up the arrays
			{
				*prvalue = ReadThisLine(file, 3, 1.0, "x for surface factor", "DataMat::DataMat", sffilename);
				if (i == 1)//see, whether it is equidistant
					dx = *prvalue - *(prvalue - 1);//first spacing
				else
				{
					if (fabs(*prvalue - *(prvalue - 1) - dx) > TOLERANCE2)
					{
						cout << "\nWARNING(" << ++warn << "): The x spacing in the custom surface factor file " << sffilename << " is not equidistant!" << endl;
						cout << "\tThe built-in surface factor table will be used instead!" << endl;
						RunParams::use_custom_sftable = false;
						delete[] rtable;
						delete[] sftable;
						rtable = nullptr;
						sftable = nullptr;
						break;
					}
				}
				if (i == sfactor_size-1)//last point, no problem, if there is not more line in the file 
					*ptemp = ReadThisLine(file, 1, 1.0, "surface factor", "DataMat::DataMat");
				else
					*ptemp = ReadThisLine(file, 1, 1.0, "surface factor", "DataMat::DataMat", sffilename);
				prvalue++;
				ptemp++;
			}
			next_line_pos = -1;
			file.close();
		}

		if (!RunParams::use_custom_sftable)//use the built in surface factor table
		{
			sfactor_size = SFACTOR_SIZE_DEF;
			rtable = surfacefact_x;
			sftable= surfacefact_y;

		}
	}



	//memory allocation for the finders
	SetArraysize(&grfinder, ngr, "grfinder", "DataMat::DataMat");
	SetArraysize(&fqfinder, nfq, "fqfinder", "DataMat::DataMat");
	SetArraysize(&fgfinder, nfg, "fgfinder", "DataMat::DataMat");
	SetArraysize(&sqfinder, nsq, "sqfinder", "DataMat::DataMat");
	//implementing the matrix grmij that allows to compute partial
	//pair correlation functions from the histograms contents
	//this is not a very efficient way to compute these matrices
	//but they are not too big and it has to be done only once

	//to avoid compiler warning
	left_bin = 0;
	right_bin = 0;
	last_assigned = 0;
	size = 0;
	for (i = 0; i < ngr; i++)
		size += grsize[i] * (mybins[assign_hist[i]]-RunParams::firstbin[i]);	//calculating the total size, as the number of mins can differ from set to set
	SetArraysize(&grmij, size, "grmij", "DataMat::DataMat");
	//positioning the finder
	size = 0;
	nmax = 0;
	for (i = 0; i < ngr; i++)
	{

		*(grfinder + i) = grmij + size;
		size += grsize[i] * (mybins[assign_hist[i]]-RunParams::firstbin[i]);
		if (*(grsize + i) > nmax)
			nmax = *(grsize + i);//determining the size of the longest g(r) series
	}

	SetArraysize(&gr_flag, nmax, "gr_flag", "DataMat::DataMat");//temporary array to show the holes in the conversion (those r values, which did not got aany 
			//bin assigned to
	pmatrix = grmij;//positioning the pointer to the matrices elements
	prvalue = edata.gr_rvalues;//pointer to the r values of the g(r) data
	for (i = 0; i < ngr; i++)//for each g(r) data set
	{
		dr = RunParams::rspacing[i];
		rmin = boxedge * xmin[assign_hist[i]]+RunParams::firstbin[i]*dr;//minimum r value for this data set (in Angstroms)
		myused_bins = mybins[assign_hist[i]] - RunParams::firstbin[i];//as the histogram is always calculated from the beginning, the matrix will be calculated only for the
		//range described by the binshift
		is_hole = 0;//no holes by default
		for (m = 0; m < nmax; m++)
			gr_flag[m] = 0;//default, no holes in the assignment

		for (m = 0; m < *(grsize + i); m++)//line (r in data set) index
		{
			if (m == 0)//first point in data set
			{
				r1 = *prvalue;
				r2 = (*(prvalue + 1) + r1)*0.5;
				tol1 = GRID_TOL;//tolerance for the lower limit, to prevent miscalculation due to rounding error
				tol2 = 0;
			}
			else
			{
				if (m == *(grsize + i) - 1)//last r value
				{
					r2 = *prvalue;
					r1 = 0.5*(*prvalue + *(prvalue - 1));
					tol1 = 0;
					tol2 = GRID_TOL;//tolerance for the higher limit, to prevent miscalculation due to rounding error
				}
				else
				{
					r1 = 0.5*(*prvalue + *(prvalue - 1));
					r2 = 0.5*(*prvalue + *(prvalue + 1));
					tol1 = 0;
					tol2 = 0;
				}
			}

			prvalue++;//go to next data point
			inc = 0;//reset the number of bins belonging to this r interval
			for (n = 0; n < myused_bins; n++)//column index (i.e. bins in histograms)
			{
				r = rmin + (n + 0.5)*dr;//central r value of the n-th bin

				if (r<(r1 - tol1) || r>(r2 + tol2))
					*pmatrix = 0;
				else
				{
					if (inc == 0)
						ptemp = pmatrix;//ptemp will be the address of the first bin in this r interval
					inc++;//that is the total number of bins included in this interval

				}
				
				pmatrix++;//this will goto the next line and next matrix automatically 
			
			}//next column (n cycle)
			
			//set matrix components and shift n
			for (k = 0; k < inc; k++)//if there are more than 1 bin included
				*ptemp++ = 1.0 / (inc);//set correponding matrix components

			if (inc == 0)//there would be a 'hole' (no bin assigned to this g(r) value due to rounding
			{
				gr_flag[m] = 1;//indicate, that this r value has no bin assigned to it
				is_hole = 1;//there is a hole in the assignment of this data set
			}
			else
				last_assigned = m;//last _assigned will contain the index of the last assigned r value at the end

		}//next line

		if (is_hole)//if there is a hole in the assignment
		{
			message_flag = 0;//no warning by default
			inc = 1;
			for (m = 0; m < *(grsize + i); m += inc)//going through the r values, except the first and the last, as this would not be considered a hole
			{
				inc = 1;//counter for adjacent holes
				first_flag = 1;//to show, that no right bin is yet assigned
				if (gr_flag[m])//this is a hole
				{
					do
					{
						if (m + inc < *(grsize + i) && gr_flag[m + inc])//this is a hole too
							inc++;
					} while (m + inc < *(grsize + i) && gr_flag[m + inc]);
					//inc will contain the number of adjacent holes
					if (inc > 1)//the bin size is larger than dr
						message_flag = 1;//write warning message at the end 

					//determine the largest index (left) bin assigned before the hole and the smallest index (right) bin after the hole
					for (n = 0; n < myused_bins; n++)//column index (i.e. bins in histograms)
					{
						//first dealing with the inner holes
						//there is at least one r value to the left to this one with assigned bin 
						if (m > 0 && *(grfinder[i] + (m - 1)*myused_bins + n) > 0)
							left_bin = n;//this will be last left assigned bin

						//there is at least one r value with bin assignement following this 
						if (m < last_assigned && *(grfinder[i] + (m + inc)*myused_bins + n)>0 && first_flag)
						{
							first_flag = 0;
							right_bin = n;//this will be the first right assigned bin
						}
						//dealing with the empty first bins, if there is any
						if (m == 0)
						{
							left_bin = (right_bin > 0 ? right_bin - 1 : 0);//the left will be the one before the right bin, if possible
						}
						if (m > last_assigned)
						{
							right_bin = (left_bin < (myused_bins - 1) ? left_bin + 1 : myused_bins - 1);//the left will be the one before the right bin, if possible
						}
					}

					if (right_bin - left_bin > 1)
					{
						cout << "\n*****ERROR*****" << endl ;
						cout << "DataMat constructor: Error during the calculation of the bin->r conversion table of the " << i + 1 << ". g(r) set!" << endl;
						cout << "Use different bin size, or a different bin_shift value!" << endl;
						cout << "Cannot run this way, exiting..." << endl;
						CleanExit();
					}
					if (inc == 1 && right_bin == left_bin)//can happen for the first and last r value
						*(grfinder[i] + (m)*myused_bins + left_bin) = 1;
					else
					{
						switch (inc)
						{
						case 1://only one unassigned r
							*(grfinder[i] + (m)*myused_bins + left_bin) = 0.5;
							*(grfinder[i] + (m)*myused_bins + right_bin) = 0.5;
							break;

						case 2://two unassigned r, can happen, if dr is smaller than bin width, but should be avoided!
							*(grfinder[i] + (m)*myused_bins + left_bin) = 1;
							*(grfinder[i] + (m + 1)*myused_bins + right_bin) = 1;
							break;

						default://more than two unassigned r, can happen, if dr is smaller than bin width, but should be avoided!
							for (n = 0; n < inc; n++)//go through al the holes
							{
								if (n < inc / 2)
									*(grfinder[i] + (m + n)*myused_bins + left_bin) = 1;
								else
									*(grfinder[i] + (m + n)*myused_bins + right_bin) = 1;
							}
							break;

						}//end of switch
					}//end of else
				}//end of this is a hole
			}//r point cycle m
			if (message_flag)
			{
				cout << "\nWARNING(" << ++warn << "): DataMat constructor: bin->r conversion table of the " << i + 1 << ". g(r) set!" << endl;
				cout << "\tThe bin size is larger than dr! Consider to use smaller bins!" << endl;
				cout << "\tContinuing, but g(r) will be less detailed, than the experimental data would require!" << endl;
			}
		}//end of if hole

	}//next data set
	delete [] gr_flag;//deleting the temporary array
	//implementing the Sine Fourier Transform matrices for S(Q) data
	//matrix elements will contain number density!!!
	size = 0;
	for (i = 0; i < nsq; i++)
		size += sqsize[i] *( mybins[assign_hist[i + ngr]]-RunParams::firstbin[ngr+i]);
	SetArraysize(&sqmij, size, "sqmij", "DataMat::DataMat");

	//positioning the finder
	size = 0;
	for (i = 0; i < nsq; i++)
	{
		*(sqfinder + i) = sqmij + size;
		size += sqsize[i] *( mybins[assign_hist[i + ngr]] - RunParams::firstbin[ngr + i]);
	}
	//fourpi=12.5663704;//rmca value
	ptemp = sqmij;//initialisation of the pointer to the matrices elements
	pqvalue = edata.sq_qvalues;//pointer to the Q values of the S(Q) data	
	for (i = 0; i < nsq; i++)//for each S(Q) data set
	{
		dr = RunParams::rspacing[ngr + i];
		rmin = boxedge * xmin[assign_hist[ngr+i]]+ RunParams::firstbin[ngr + i]*dr;//minimum r value (in Angstroms)
		for (m = 0; m < *(sqsize + i); m++)//line (q) index
		{
			for (n = 0; n < mybins[assign_hist[i + ngr]]-RunParams::firstbin[ngr+i]; n++)//column (r) index
			{
				//r=rmin+n*dr;//lower limit

				//r=rmin+(n+0.5)*dr;//central approx r value...
				//*ptemp=fourpi*r*sin(*pqvalue*r)/(*pqvalue)*dr;//...and atrix element


				//integral approx
				rinf = rmin + n * dr;
				rsup = rinf + dr;
				
				*ptemp = fourpi * rho / *pqvalue *(sin(*pqvalue * rsup) / (*pqvalue) - rsup * cos(*pqvalue * rsup) - \
					sin(*pqvalue * rinf) / (*pqvalue) + rinf * cos(*pqvalue * rinf));//it is Q*S(Q) 

#if !defined _MUL_SCAT_VECTOR
				*ptemp /= *pqvalue;//dividing it with Q for normal calculation mode (if _MUL_SCAT_VECTOR is switched off)
#endif

			ptemp++;
			}
			pqvalue++;//next Q value/line (goes to next data set automatically)
		}
	}//next data set

	//implementing the Sine Fourier Transform matrices for F(Q) data
	//matrix elements will contain number density!!!
	size = 0;
	for (i = 0; i < nfq; i++)
		size += fqsize[i] * (mybins[assign_hist[i + ngr + nsq]]- RunParams::firstbin[ngr+nsq+i]);
	SetArraysize(&fqmij, size, "fqmij", "DataMat::DataMat");

	//positioning the finder
	size = 0;
	for (i = 0; i < nfq; i++)
	{
		*(fqfinder + i) = fqmij + size;
		size += fqsize[i] * (mybins[assign_hist[i + ngr + nsq]]- RunParams::firstbin[ngr + nsq + i]);
	}

	ptemp = fqmij;//initialisation of the pointer to the matrices elements
	pqvalue = edata.fq_qvalues;//pointer to the Q values of the F(Q) data

	for (i = 0; i < nfq; i++)//for each F(Q) data set
	{
		dr = RunParams::rspacing[ngr + nsq + i];
		rmin = boxedge * xmin[assign_hist[ngr + nsq+i]]+ RunParams::firstbin[ngr +nsq+ i]*dr;//minimum r value (in Angstroms)
		for (m = 0; m < *(fqsize + i); m++)//line (q) index
		{
			for (n = 0; n < mybins[assign_hist[i + ngr + nsq]]- RunParams::firstbin[ngr + nsq + i]; n++)//column (r) index
			{
				//r=rmin+n*dr;//lower limit

				//r=rmin+(n+0.5)*dr;//central approx r value...
				//*ptemp=fourpi*r*sin(*pqvalue*r)/(*pqvalue)*dr;//...and matrix element


				//integral approx
				rinf = rmin + n * dr;
				rsup = rinf + dr;
				*ptemp = fourpi * rho / *pqvalue *(sin(*pqvalue * rsup) / (*pqvalue) - rsup * cos(*pqvalue * rsup) - \
					sin(*pqvalue * rinf) / (*pqvalue) + rinf * cos(*pqvalue * rinf));//it is Q*S(Q) 

#if !defined _MUL_SCAT_VECTOR
				*ptemp /= *pqvalue;//dividing it with Q for normal calculation mode (if _MUL_SCAT_VECTOR is switched off)
#endif
				ptemp++;
			}
			pqvalue++;//next Q value/line (goes to next data set automatically)
		}
	}//next data set

	//implementing the Sine Fourier Transform matrices for F(g) data
	//matrix elements will contain number density!!!

	size = 0;
	for (i = 0; i < nfg; i++)
		size += fgsize[i] *( mybins[assign_hist[i + ngr + nsq + nfq]]- RunParams::firstbin[ngr + nsq + nfq+i]);
	SetArraysize(&fgmij, size, "fgmij", "DataMat::DataMat");
	//positioning the finder
	size = 0;
	for (i = 0; i < nfg; i++)
	{
		*(fgfinder + i) = fgmij + size;
		size += fgsize[i] *( mybins[assign_hist[i + ngr + nsq + nfq]]- RunParams::firstbin[ngr + nsq + nfq+i]);
	}

	ptemp = fgmij;//initialisation of the pointer to the matrices elements
	pgvalue = edata.fg_gvalues;//pointer to the g values of the F(g) data

	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		dr = RunParams::rspacing[ngr + nsq + nfq + i];
		rmin = boxedge * xmin[assign_hist[ngr +nsq+nfq+ i]]+ RunParams::firstbin[ngr + nsq+nfq+i]*dr;//minimum r value (in Angstroms)
		for (m = 0; m < *(fgsize + i); m++)//line (g) index
		{
			for (n = 0; n < mybins[assign_hist[i + ngr + nsq + nfq]]- RunParams::firstbin[ngr + nsq +nfq+ i]; n++)//column (r) index
			{
				//r=rmin+n*dr;//lower limit

				//r=rmin+(n+0.5)*dr;//central approx r value...
				//*ptemp=fourpi*r*sin(*pgvalue*r)/(*pgvalue)*dr;//...and matrix element


				//integral approx
				rinf = rmin + n * dr;
				rsup = rinf + dr;
				*ptemp = fourpi * rho / *pgvalue *(sin(*pgvalue * rsup) / (*pgvalue) - rsup * cos(*pgvalue * rsup) - \
					sin(*pgvalue * rinf) / (*pgvalue) + rinf * cos(*pgvalue * rinf));//it is g*S(g) 

#if !defined _MUL_SCAT_VECTOR
				*ptemp /= *pgvalue;//dividing it with g for normal calculation mode (if _MUL_SCAT_VECTOR is switched off)
#endif
				ptemp++;
			}
			pgvalue++;//next g value/line (goes to next data set automatically)
		}
	}//next data set
	//implementing the tables used to convert histogram counts to PPCF values
	SetArraysize(&nfactor, npartials, "nfactor", "DataMat::DataMat");

	if (calcmode)//only in case of _NO_PERIODIC
		n = npartials;
	else
	{
		ntot_bins = 0;
		for (i = 0; i < ndiffbin; i++)
			ntot_bins += nbins[i];

		n = npartials * ntot_bins;//always set using nbins, also in case of _VIBR_AMP only nbins_used is 
		//used, but always starts using it from the beginning of the partial+min_ind for the thread
	}

	SetArraysize(&ntable, n, "ntable", "DataMat::DataMat");//there will be an ntable part for each different bin size,
	//following each other in{ for binsize {for partial{for bins}} } order
	n=npartials * ndiffbin;
	
	SetArraysize(&ntable_finder, n, "ntable_finder", "DataMat::DataMat");//finders for  each partial ntable
	n = 0;
	size = 0;
	for (i = 0; i < ndiffbin; i++)//initializing the finder
	{
		for (j = 0; j < npartials; j++)
		{
			ntable_finder[n] = ntable + size;
			size+=nbins[i];
			n++;
		}
	}
	
	CalcdV(0, rtable, sftable);//calculating the ntable (normalization table for the histograms)
#ifdef _LOCAL_INV
	SetArraysize(&dV_loc,HistoSet::max_loc_nbins,"dV_loc","DataMat::DataMat");
	CalcdV(1,rtable, sftable);//calculating the dV_loc array for the local invariance histogram bins
#endif

	if (RunParams::use_custom_sftable)
	{
		if (rtable != NULL)
			delete[] rtable;
		if (sftable != NULL)
			delete[] sftable;
	}
};//end of constructor
	
//calculating the renormalization table in case of mode=0, or storing the dV for mode=1 
//if no extra features like _VIBR_AMP ot _NO_PERIODIC is used, then there can be more than one bin sizes,
//Therefor there will be a normalization table for every one. In case of _VIBR_AMP and/or _NO_PERIODIC there is only one binsize, and the
//ib_size cycle is executed only once
void DataMat::CalcdV(int mode, double *rtable, double *sftable)
{
	int i,j,k,ib_size,n1,n2,my_nbins=0;
	double my_xmin=0,*pf=0,deltax=0,dv=0;
	double fourpi;
		
#ifdef _NO_PERIODIC
	int calc_mode;//flag to indicate what type of correction is needed
	int dV_flag=0;//whether to store the dV
	double *pn = 0;
	double x1;//the height of the sphere cap A with base radius y1
	double x2;//the height of the sphere cap A with base radius y2
	double m1;//the height of the sphere cap B
	double m2;//the height of the sphere cap C
	double rBj;//the radius of sphere B
	double rCj;//the radius of sphere C
	double z;//rBj-m1+x1 or rCj-m2+x2
	double dvj;//the volume of the rBj sphere

	double	V_cap_Ay1;//the volume of the flat based cap of sphere A with base radius y1
	double	V_cap_Ay2;//the volume of the flat based cap of sphere A with base radius y2
	double	V_cap_By1;//the volume of the flat based cap of sphere B with base radius y1
	double	V_cap_Cy2;//the volume of the flat based cap of sphere C with base radius y2
	double	V_cap_B;//the volume of V_cap_By1-V_cap_Ay1
	double	V_cap_C;//the volume of V_cap_Cy2-V_cap_Ay2
	
#else
	int surface_ind1,surface_ind2,m;
	double sfactor,temp;
	double xinf=0, xsup=0;
	double boxedge3=boxedge*boxedge*boxedge;

#endif
#ifndef _NO_PERIODIC
	double dx = 0;
	double sq2;
	if (rtable != NULL)
		dx = *(rtable + 1) - *rtable;// discretizing interval in the sfactor tabulation
	sq2 = 1.4142135623730951;//squareroot of 2
#endif
	fourpi=4*PI;
	
	if (mode==0)//calculating the normalization table for the histograms
	{
	
#ifdef _NO_PERIODIC
		pn = ntable;
		if (calcmode==0)
				dV_flag=1;

#endif
		pf=nfactor;

		for(i=0;i<ntypes;i++)//type of central atom
		{
			n1=*(pnatoms+i);//number of central atoms
			for(j=i;j<ntypes;j++)//type of neighbour atom
			{
				n2=*(pnatoms+j);//number of neighbour atoms
				if(i==j)//central type and neighbour type are the same
				{
					*pf=volume/n1/(n1-1);//part of the table element independent
#ifdef _NO_PERIODIC
					if (calcmode)
					{
						*pn=*pf;
						pn++;
					}
#endif
					*pf*=2;//as normally the i-i type histogram is calculated only for one of the pair being central
					pf++;
				}
				else
				{
					//of  the r  value
					*pf=volume/n1/n2;
#ifdef _NO_PERIODIC
					if (calcmode)
					{
						*pn=*pf/2;//has to be divided by two, as during the ppcf calculation the sum of i-j and j-i local densities are calculated
						pn++;
					}
#endif
					pf++;
				}
			}
		}
	//	pn=ntable;
	}
#ifdef _LOCAL_INV
	else//calculating the dV for the local invariance bins
	{
		deltax=(HistoSet::max_loc_r-HistoSet::min_loc_r)/HistoSet::max_loc_nbins;///bins width (in reduced units)
		my_nbins=HistoSet::max_loc_nbins;
		my_xmin=HistoSet::min_loc_r;
	}
#endif
	for (ib_size = 0; ib_size < ndiffbin; ib_size++)//if _NO_PERIODIC and/or _VIBR_AMP then ndiffbin is 1, and only executed once as before different bin size was introduced
	{
		if (mode == 0)
		{
			deltax = (xmax[ib_size] - xmin[ib_size]) / nbins[ib_size];///bins width (in reduced units)
			my_nbins = mybins[ib_size];
			my_xmin = xmin[ib_size];
		}
#ifdef _NO_PERIODIC

		for (j = 0; j < nbins[0]; j++)
		{
			rBj = (my_xmin + (j + 1)*deltax)*boxedge;//radius of sphere B in Angstrom
			rCj = (my_xmin + j * deltax)*boxedge;//radius of sphere C in Angstrom
			dvj = fourpi / 3.0*(pow(rBj, 3) - pow(rCj, 3));//in Ang3
			if (dV_flag)
			{
				dV[j] = dvj;
				dv = dvj;//needed for normal ntable calculation
			}
#endif
			for (k = 0; k < my_nbins; k++)//for each table element
			{
#ifndef _NO_PERIODIC
				//k is running through the bins containing the centre of the B and C sphere
				xinf = my_xmin + k * deltax;
				xsup = xinf + deltax;
				//find the value of the surface factor
				if (xsup <= 1)
				{

					dv = fourpi / 3.0*boxedge3*(xsup*xsup*xsup - xinf * xinf*xinf);//integral
#else
				//It assumes, that the R0 is divisible by dr
				z = R0 - (k + 0.5)*deltax*boxedge;//the distance of the middle of the sphere B and C to the surface of sphere A in Angstrom 
				calc_mode = 2;//no cerrection necessary
				if (rBj - z > deltax*boxedge)
					calc_mode = 0;//both sphere B and C protruding
				else
				{
					if (rBj - z > 0)
					{
						calc_mode = 1;//only sphere B is protruding
						V_cap_C = 0;
					}
					else
						dV_in[k*my_nbins + j] = dvj;
				}

				switch (calc_mode)
				{
					case 0://calculate the protruding part of sphere C
					{
						x2 = (z*z - rCj * rCj) / (2 * z - 2 * R0);//in Angstrom
						m2 = rCj - z + x2;//height of protruding sphere cap of sphere C in an Angstrom
						V_cap_Cy2 = PI / 3 * m2*m2*(3 * rCj - m2);//protruding cap of sphere C with base radius y2 in Ang3
						V_cap_Ay2 = PI / 3 * x2*x2*(3 * R0 - x2);//protruding cap of sphere A with base radius y2 in Ang3
						V_cap_C = V_cap_Cy2 - V_cap_Ay2;//protruding part of sphere C from sphere A in Ang3
					}
					case 1://caculate the protruding part of sphere B
					{
						x1 = (z*z - rBj * rBj) / (2 * z - 2 * R0);//in Angstrom
						m1 = rBj - z + x1;//height of protruding sphere cap of sphere B in an Angstrom
						V_cap_By1 = PI / 3 * m1*m1*(3 * rBj - m1);//protruding cap of sphere B with base radius y1 in Ang3
						V_cap_Ay1 = PI / 3 * x1*x1*(3 * R0 - x1);//protruding cap of sphere A with base radius y1 in Ang3
						V_cap_B = V_cap_By1 - V_cap_Ay1;//protruding part of sphere B from sphere A in Ang3
						dV_in[k*my_nbins + j] = dvj - (V_cap_B - V_cap_C);//protruding part of sphere ring B-C from sphere A in Ang3
						break;
					}
				}
#endif
#ifndef _NO_PERIODIC
				}
			else
			{
				//xsup>1
				if (xsup <= sq2)
				{
					if (xinf >= 1)
						//dv=fourpi/3.0*boxedge*boxedge*boxedge*(xsup*xsup*xsup-xinf*xinf*xinf);//RMCA mimick
						dv = -2.0*fourpi / 3.0*(xsup*xsup*xsup - xinf * xinf*xinf)*boxedge3 + \
						1.5*fourpi*(xsup*xsup - xinf * xinf)*boxedge3;
					else
					{
						//dv=fourpi/3.0*boxedge*boxedge*boxedge*(xsup*xsup*xsup-xinf*xinf*xinf);//RMCA mimick
						dv = fourpi / 3.0*boxedge3*(1 - xinf * xinf*xinf);
						dv += -2.0*fourpi / 3.0*(xsup*xsup*xsup - 1)*boxedge3 + \
							1.5*fourpi*(xsup*xsup - 1)*boxedge3;
					}
				}
				else
				{
					//xsup<=sq3)

					if (xinf <= sq2)//The calculation has to be devided into two parts
					{
						//x is below sq2
						//multiplication with boxedge raised to the 3-d power will come at the end
						dv = -2.0*fourpi / 3.0*(sq2*sq2*sq2 - xinf * xinf*xinf) + \
							1.5*fourpi*(sq2*sq2 - xinf * xinf);
						//second part:x is above sq2
						//sets the lower boundary for the integration of the surface
						//factor to the first element of the rtable array
						surface_ind1 = 0;
					}
					else
					{
						//The xinf-->xsup intervallum will be divided into three parts:
						//1.) xinf-->rtable[surface_ind1] (tabulated element just above xinf)
						//2.) rtable[surface_ind1]-->rtable[surface_ind2] (tabulated element just below xsup)
						//3.) rtable[surface_ind2]-->xsup
						//The surface factors belonging two xinf and xsup has to be calculated
						//by linear interpolation
						//The volume contributions will be determined by the trapezoidal rule

						//index of the array element just above the xinf
						surface_ind1 = int((xinf - rtable[0]) / dx) + 1;

						//This is a precaution, because the rtable intervals
						//are not exactly the same
						//Check, whether the index is right for the element just below the xinf
						while (xinf < rtable[surface_ind1 - 1])//set it right
							surface_ind1--;
						//Check, whether the index is right for the element just above the xinf
						while (xinf > rtable[surface_ind1])//set it right
							surface_ind1++;

						if (surface_ind1 > sfactor_size - 1)//cannot really happen
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "DataMat constaructor: surface_ind1, r exceeds the tabulated values in sfactorcube!" << endl;
							cout << "Use smaller xmax in *.dat file! Cannot run this way, exiting..." << endl;
							CleanExit();
						}

						//linear interpolation between surface_ind1-1 and surface_ind1
						// in the sfactor array to get the sfactor belonging to xinf
						sfactor = sftable[surface_ind1 - 1] + (xinf - rtable[surface_ind1 - 1]) / (rtable[surface_ind1] - rtable[surface_ind1 - 1])*\
							(sftable[surface_ind1] - sftable[surface_ind1 - 1]);
						//calculate the volume contribution from the xinf -->rtable[surface_ind1]
						dv = (rtable[surface_ind1] - xinf)*(sftable[surface_ind1] + sfactor) / 2;

					}

					//index of the array element just below the xsup
					surface_ind2 = int((xsup - rtable[0]) / dx);

					//This is a precaution, because the rtable intervals
					//are not exactly the same
					//Check, whether the index is right for the element just below the xsup
					if (surface_ind2 > sfactor_size - 1 && k == my_nbins - 1)//this is te last bin, for local inv it can be over sqrt3, set the index back
					{
						surface_ind2 = sfactor_size - 1;//the last element
						if (surface_ind2 < surface_ind1)//this cannot really happen
						{
							cout << "\n*****ERROR*****" << endl;
							cout << "DataMat constaructor: surface_ind2 is larger than surface_in1!" << endl;
							cout << "Cannot run this way, exiting.." << endl;
							CleanExit();
						}
					}

					while (xsup < rtable[surface_ind2])//set it right
						surface_ind2--;
					//Check, whether the index is right for the element just above the xsup
					while (xsup > rtable[surface_ind2 + 1] && (k < my_nbins - 1 || surface_ind2 < sfactor_size - 1))//set it right, it can happen , that the xsup for the last bin 
						//exceeds rtable[surface_ind2], in this case do not increase, as there is no more element in the sftable, and the contribution would be zero anyhow
						surface_ind2++;

					if (surface_ind2 > sfactor_size - 1)//cannot really happen
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "DataMat constaructor: surface_ind2, r exceeds the tabulated values in sfactorcube!" << endl;
						cout << "Use smaller xmax in *.dat file! Cannot run this way, exiting..." << endl;
						CleanExit();
					}

					//linear interpolation between surface_ind2 and surface_ind2+1
					// in the sfactor array to get the sfactor belonging to xsup
					if (surface_ind2 + 1 < sfactor_size)//do not do it, if xsup is over the sftable
					{
						sfactor = sftable[surface_ind2] + (xsup - rtable[surface_ind2]) / (rtable[surface_ind2 + 1] - rtable[surface_ind2])*\
							(sftable[surface_ind2 + 1] - sftable[surface_ind2]);
						//calculate the volume contribution from the xsup -->rtable[surface_ind2]
						dv += (xsup - rtable[surface_ind2])*(sftable[surface_ind2] + sfactor) / 2;
					}

					//integration over the elements between surface_ind1 - surface_ind2 
					//to get the volume contribution using the trapezoidal rule
					temp = 0;
					for (m = surface_ind1 + 1; m < surface_ind2; m++)
						temp += sftable[m];
					temp += (sftable[surface_ind1] + sftable[surface_ind2]) / 2;
					dv += temp * dx;
					dv *= boxedge3;//to get Angstrom
				}//end of xsup<=3
			}//end of if (xsup>1)
			if (dv <= 0)
				dv = 0; //this is just a precaution to deal with the lower limit dv
#endif
#ifdef _NO_PERIODIC

			if (calcmode == 0 && k == 0)//only go in, if VIBR_AMP and/or EXAFS is also used
			{
#else
			j = k;
#endif				

			if (mode == 0)
			{
				for (i = 0; i < npartials; i++)//setting all the partials for this bin
				{

					*(ntable_finder[ib_size*npartials + i] + j) = nfactor[i] / dv;//divide by the s factor)
				}
			}
#ifdef _NO_PERIODIC
			}

#endif
#ifdef _LOCAL_INV
		if (mode)
			dV_loc[k] = dv;
#endif

			}//end of bin cycle k
#ifdef _NO_PERIODIC
		}//end of bin cycle j
#endif
	}//end of ib_size cycle

};

void DataMat::Save() const
{
	

	int i,j,m,n,ib,usebins;
	double *ptemp;

	ofstream file;
	//save to *.datmat file
	mystrcpy(tempfilename, FILE_NAME_SIZE+10, filename);
	mystrcat(tempfilename, FILE_NAME_SIZE+10, ".datmat");
	OpenFile(file,tempfilename,"DataMat::Save",0);//	open file, check, whether it was successfully opened
	file<<"This is an object of class DataMat with dimensions:"<<endl;
	file<<"Number of atom types: "<<ntypes<<endl;
	file<<"number of partials: "<<npartials<<endl;
#ifdef _NO_PERIODIC
	file<<"sample volume:"<<volume<<endl;
#else
	file<<"cell volume:"<<volume<<endl;
#endif
	file<<"Total number density: "<<rho<<endl;
	file<<"Number of g(r) data sets->"<<ngr<<endl;

	for(i=0;i<ngr;i++)//for each g(r) data set
	{
		file<<"sets size (number of data points): ";
		file<<*(grsize+i)<<"  "<<endl;
	}
	file<<endl;

	file<<"Number of S(Q) data sets->"<<nsq<<endl;
	for(i=0;i<nsq;i++)//for each S(Q) data set
	{
		file<<"sets size (number of data points): ";
		file<<*(sqsize+i)<<"  "<<endl;
	}
	file<<endl;

	file<<"Number of F(Q) data sets->"<<nfq<<endl;

	for(i=0;i<nfq;i++)//for each F(Q) data set
	{
		file<<"sets size (number of data points): ";
		file<<*(fqsize+i)<<"  "<<endl;
	}
	file<<endl;

	file << "Number of F(g) data sets->" << nfg << endl;

	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		file << "sets size (number of data points): ";
		file << *(fgsize + i) << "  " << endl;
	}
	file << endl;

	//now the matrices contents

	file<<"smoothing/sampling table to g(r) data:"<<endl;
	ptemp=grmij;//initialisation of the pointer to the matrices elements
	for(i=0;i<ngr;i++)//for each g(r) data set
	{
		file<<"g(r) set #"<<i+1<<"/"<<ngr<<endl;
		for(m=0;m<*(grsize+i);m++)//line index
		{
#ifdef _VIBR_AMP
			usebins = nbins_used;//number of bins used in the g(r) histograms
#else
			usebins = nbins[assign_hist[i]];
#endif
			usebins -= RunParams::firstbin[i];//to account for the binshift
			for(n=0;n<usebins;n++)//column index
			{
				file<<*ptemp<<" ";
				ptemp++;
			}
			file<<endl;
		}
	}//next data set

	file<<"\n Normalizing tables for the histograms "<<endl;
#ifdef _NO_PERIODIC
	ptemp=dV_in;
	//the histogram is calculated from the beginnnig regardless the binshift, since custom binsize is introduced, even if in this case there can be only one binsize and binshift
	usebins = mybins[0];
	file<<"\t";
	for (i=0;i<usebins;i++)
		file<<i<<"\t";
	file<<endl;

	for (j=0;j< usebins;j++)
	{
		file<<j<<"\t";
		for (i=0;i< usebins;i++)
		{
			file<<*ptemp++<<"\t";
		}
		file<<endl;

	}
	if (calcmode==0)//this happens, if ek is used or _VIBR_AMP is set
	{
		file<<"\n ntable tables for the histograms "<<endl;
#endif
	
	
		ptemp=ntable;
		file.precision(16);
		for (ib=0;ib<ndiffbin;ib++)
		{
#ifdef _VIBR_AMP
			usebins = nbins_used;//number of bins used in the g(r) histograms
#else
			usebins = nbins[ib];
#endif
			for (i = 0; i < npartials; i++)
			{
				file << "\n binsize #"<<ib+1<<" binwidth  "<<RunParams::rspacing_diff[ib]<<" Ang, partial #" << i + 1 << "/" << npartials << "  bin index, normalising value" << endl;

				for (j = 0; j < usebins; j++)
				{
					file << J7 << j + 1 << "\t" << J10 << *ptemp << endl;
					ptemp++;
				}
			}
		}
#ifdef _NO_PERIODIC
	}
#endif
	
#ifdef _LOCAL_INV
	file<<"\n dV for the local histograms "<<endl;
	ptemp=dV_loc;
	file<<"\n local bin index, dV"<<endl;
	for(j=0;j<HistoSet::max_loc_nbins;j++)
	{
		file<<J7<<j+1<<"\t"<<J10<<*ptemp<<endl;
		ptemp++;
	}
#endif
		
	file<<"\n Sine Fourier tables: "<<endl;
	ptemp=sqmij;//initialisation of the pointer to the matrices elements
	for(i=0;i<nsq;i++)//for each S(Q) data set
	{
		file<<"S(Q) set #"<<i+1<<"/"<<nsq<<endl;
		for(m=0;m<*(sqsize+i);m++)//line index
		{
			for(n=0;n<mybins[assign_hist[ngr+i]]-RunParams::firstbin[ngr+i];n++)//column index
			{
				file<<J10<<*ptemp<<"\t";
				ptemp++;
			}
			file<<endl;
		}
	}//next data set
		
	//now the F(Q) data
	ptemp=fqmij;//initialisation of the pointer to the matrices elements
	for(i=0;i<nfq;i++)//for each F(Q) data set
	{
		file<<"F(Q) set #"<<i+1<<"/"<<nfq<<endl;
		for(m=0;m<*(fqsize+i);m++)//line index
		{
			for(n=0;n<mybins[assign_hist[ngr +nsq+ i]]-RunParams::firstbin[ngr +nsq+ i];n++)//column index
			{
				file<<J10<<*ptemp<<"\t";
				ptemp++;
			}
			file<<endl;
		}
	}//next data set	
	//now the F(g) data
	ptemp = fgmij;//initialisation of the pointer to the matrices elements
	for (i = 0; i < nfg; i++)//for each F(g) data set
	{
		file << "F(g) set #" << i + 1 << "/" << nfg << endl;
		for (m = 0; m < *(fgsize + i); m++)//line index
		{
			for (n = 0; n < mybins[assign_hist[ngr + nsq+nfq+i]]-RunParams::firstbin[ngr +nsq+nfq+ i]; n++)//column index
			{
				file << J10 << *ptemp << "\t";
				ptemp++;
			}
			file << endl;
		}
	}//next data set	
	file.precision(8);
	file.close();
};		
//----------------------sets the static members---------------
void DataMat::SetDataMatParams(RunParams &rundat)
{
	ndiffbin = RunParams::ndiffbin;
	nbins=rundat.nbins;//number of bins in the g(r) histograms, can differ from each other normally
	assign_hist = rundat.assign_hist;
	SetArraysize(&mybins, ndiffbin, "mybins", "DataMat::SetDataMatParams");
	//if _VIBR_AMP or _NO_PERIODIC is used, then only one bin size can be used, as originally
#ifdef _VIBR_AMP
	nbins_used=HistoSet::nbins_used;
#ifdef _NO_PERIODIC
	int i;
	mybins[0]=nbins[0];//dV_in has to have dimensions nbins*nbins
#else
	mybins[0]=nbins_used;//number of bins used in the g(r) histograms
#endif
#else
	int i;
	for (i=0;i<ndiffbin;i++)
		mybins[i]=nbins[i];
#endif
	
	ntypes=RunParams::ntypes;//number of atom types
	npartials=ntypes*(ntypes+1)/2;;//number of partials
	pnatoms=RunParams::pnatoms;//pointer to number of atoms per type
	xmin=rundat.xmin;//minimum value for the histogram in reduced units
	xmax=RunParams::xmax;//maximum value for the histogram in reduced units
	boxedge=RunParams::boxedge;//size of the simulation cell in Angstrom
#ifdef _NO_PERIODIC
	R0=RunParams::R0;//the radius of the speric sample inside the box
	volume=4.0/3.0*PI*R0*R0*R0;//volume of the spherical sample
	SetArraysize(&dV_in,nbins[0]*nbins[0],"dV_in","DataMat::SetDataMatParams");
	calcmode=1;
	if (RunParams::nek>0)
	{
		calcmode=0;//do not calculate from the local hist, but from the previously calculated periodic hist
		i=nbins[0];
	}
#ifdef _VIBR_AMP
	calcmode=0;
	i=nbins[0];
#endif

	//dV is needed to calculate a periodic-like histogram, if _VIBR_AMP and or EXAFS is used
	//the dimension is always nbins[0]
	SetArraysize(&dV,i,"dV","DataMat::SetDataMatParams");

#else
	volume=8.0*boxedge*boxedge*boxedge;//cell volume
#endif
	rho=RunParams::rho;//total number density
	ngr=RunParams::ngr;//number of g(r) data sets
	nsq=RunParams::nsq;//number of S(Q) data sets
	nfq=RunParams::nfq;//number of F(Q) data sets
	nfg = RunParams::nfg;//number of F(g) data sets
}



