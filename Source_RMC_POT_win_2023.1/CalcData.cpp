//source CalcData.cpp
//Last changed 24.02.2022

#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//classes2.h included through this

int CalcData::ntypes;//number of types
int CalcData::npartials;//number of partials
int CalcData::ngr;//number of g(r) data sets
int CalcData::nsq;//number of S(Q) data sets
int CalcData::nfq;//number of F(Q) data set
int CalcData::nfg;//number of F(g) data set
int CalcData::nek;//number of F(Q) data set

//------------constructor-------
CalcData::CalcData(ExptsData &edat, CalcPart &calcpartial, ChiSquared &chisq)
					:edata(edat), calcp(calcpartial), chi(chisq)
{
	//constructor of the CalcData class 
	//the objects are for handle initialization
	int i;
	
	//memory allocation for data arrays
	SetArraysize(&grvalues,edata.rused_tot,"grvalues","CalcData::CalcData");//g(r) values
	SetArraysize(&sqvalues,edata.sqused_tot,"sqvalues","CalcData::CalcData");//S(Q) values
	SetArraysize(&fqvalues,edata.fqused_tot,"fqvalues","CalcData::CalcData");//F(Q) values
	SetArraysize(&fgvalues, edata.fgused_tot, "fgvalues", "CalcData::CalcData");//F(g) values
	SetArraysize(&ekvalues,edata.ekused_tot,"ekvalues","CalcData::CalcData");//E(k) values

	//finder pointers
	SetArraysize(&grfinder,ngr,"grfinder","CalcData::CalcData");//finder of g(r) data
	SetArraysize(&sqfinder,nsq,"sqfinder","CalcData::CalcData");//finder of S(Q) data
	SetArraysize(&fqfinder,nfq,"fqfinder","CalcData::CalcData");//finder of F(Q) data	
	SetArraysize(&fgfinder, nfg, "fgfinder", "CalcData::CalcData");//finder of F(g) data	
	SetArraysize(&ekfinder,nek,"ekfinder","CalcData::CalcData");//finder of E(k) data
	if (ExptsData::is_IQ>0)
	{
		SetArraysize(&iqvalues, edata.fqused_tot, "iqvalues", "CalcData::CalcData");//I(Q) values
		SetArraysize(&iqfinder, nfq, "iqfinder", "CalcData::CalcData");//finder of I(Q) data	
		SetArraysize(&compton, edata.fqused_tot, "compton correction", "CalcData::CalcData");//compton correction
		SetArraysize(&comptfinder, nfq, "comptfinder", "CalcData::CalcData");//finder of compton correction	
	}
	else
	{
		iqvalues = NULL;
		iqfinder = NULL;
		compton = NULL;
		comptfinder = NULL;
	}

	//initialize the finder arrays
	for(i=0;i<ngr;i++)
		grfinder[i]=grvalues+edata.grused_cum[i];
		
	for(i=0;i<nsq;i++)
		sqfinder[i]=sqvalues+edata.sqused_cum[i];
	
	for(i=0;i<nfq;i++)
		fqfinder[i]=fqvalues+edata.fqused_cum[i];

	for (i = 0; i < nfg; i++)
		fgfinder[i] = fgvalues + edata.fgused_cum[i];

	for(i=0;i<nek;i++)
		ekfinder[i]=ekvalues+edata.ekused_cum[i];
	if (ExptsData::is_IQ>0)
	{
		for (i = 0; i < nfq; i++)
		{
			iqfinder[i] = iqvalues + edata.fqused_cum[i];//I(Q) is calculated form F(Q)
			comptfinder[i]=compton+ edata.fqused_cum[i];//compton
		}

	}
				
}//end of the constructor	

//----copy constructor
CalcData::CalcData(CalcData &source, ExptsData &edat,CalcPart &calcpartial, ChiSquared &chisq)
			:edata(edat), calcp(calcpartial), chi(chisq) 
{
	//copy constructor of the CalcData class 
	//The source is for providing he data, the other objects are for handle initialization
	int i;

	//memory allocation for data arrays
	SetArraysize(&grvalues,edata.rused_tot,"grvalues","CalcData::CalcData");//g(r) values
	for(i=0;i<edata.rused_tot;i++)
		*(grvalues+i)=*(source.grvalues+i);//copying from source

	SetArraysize(&sqvalues,edata.sqused_tot,"sqvalues","CalcData::CalcData");//S(Q) values
	for(i=0;i<edata.sqused_tot;i++)
		*(sqvalues+i)=*(source.sqvalues+i);//copying from source

	SetArraysize(&fqvalues,edata.fqused_tot,"fqvalues","CalcData::CalcData");//F(Q) values
	for(i=0;i<edata.fqused_tot;i++) 
		*(fqvalues+i)=*(source.fqvalues+i);//copying from source
	
	SetArraysize(&fgvalues, edata.fgused_tot, "fgvalues", "CalcData::CalcData");//F(g) values
	for (i = 0; i < edata.fgused_tot; i++)
		*(fgvalues + i) = *(source.fgvalues + i);//copying from source

	SetArraysize(&ekvalues,edata.ekused_tot,"ekvalues","CalcData::CalcData");//E(k) values
	for(i=0;i<edata.ekused_tot;i++) 
		*(ekvalues+i)=*(source.ekvalues+i);//setting array contents to 0
	
	//finder pointers
	SetArraysize(&grfinder,ngr,"grfinder","CalcData::CalcData");//finder of g(r) data
	SetArraysize(&sqfinder,nsq,"sqfinder","CalcData::CalcData");//finder of S(Q) data
	SetArraysize(&fqfinder,nfq,"fqfinder","CalcData::CalcData");//finder of F(Q) data	
	SetArraysize(&fgfinder, nfg, "fgfinder", "CalcData::CalcData");//finder of F(g) data
	SetArraysize(&ekfinder,nek,"ekfinder","CalcData::CalcData");//finder of E(k) data	
			
	//initialize the finder arrays
	for(i=0;i<ngr;i++)
		grfinder[i]=grvalues+edata.grused_cum[i];
		
	for(i=0;i<nsq;i++)
		sqfinder[i]=sqvalues+edata.sqused_cum[i];
	
	for(i=0;i<nfq;i++)
		fqfinder[i]=fqvalues+edata.fqused_cum[i];

	for (i = 0; i < nfg; i++)
		fgfinder[i] = fgvalues + edata.fgused_cum[i];

	for(i=0;i<nek;i++)
		ekfinder[i]=ekvalues+edata.ekused_cum[i];

	if (ExptsData::is_IQ>0)
	{
		SetArraysize(&iqvalues, edata.fqused_tot, "iqvalues", "CalcData::CalcData");//I(Q) values
		for (i = 0; i < edata.fqused_tot; i++)
			*(iqvalues + i) = *(source.iqvalues + i);//copying from source
		SetArraysize(&compton, edata.fqused_tot, "compton correction", "CalcData::CalcData");//compton correction
		for (i = 0; i < edata.fqused_tot; i++)
			*(compton + i) = *(source.compton + i);//copying from source
		SetArraysize(&iqfinder, nfq, "iqfinder", "CalcData::CalcData");//finder of I(Q) data	
		
		SetArraysize(&comptfinder, nfq, "comptfinder", "CalcData::CalcData");//finder of compton correction	
		for (i = 0; i < nfq; i++)
		{
			iqfinder[i] = iqvalues + edata.fqused_cum[i];
			comptfinder[i] = compton + edata.fqused_cum[i];//compton
		}
	}
	else
	{
		iqvalues = NULL;
		iqfinder = NULL;
		compton = NULL;
		comptfinder = NULL;
	}
	
}//end of the constructor	

//----------sets the static parameters--------------
void CalcData::SetCalcDataParams()
{
	ntypes=RunParams::ntypes;
	npartials=ntypes*(ntypes+1)/2;
	ngr=ExptsData::ngr;//number of g(r) data sets
	nsq=ExptsData::nsq;//number of S(Q) data sets
	nfq=ExptsData::nfq;//number of F(Q) data sets
	nfg = ExptsData::nfg;//number of F(g) data sets
	nek=ExptsData::nek;//number of E(k) data sets
};

//==========================save==============
void CalcData::Save() const
{
	int i,j;
	double *px,*py,*pi=NULL;
	ofstream file;
	mystrcpy(tempfilename, FILE_NAME_SIZE + 10, filename);
	mystrcat(tempfilename, FILE_NAME_SIZE + 10, ".calcdat");
	OpenFile(file,tempfilename,"CalcData::Save",0);//open file, check, whether it was successfully opened
	//saves contents in a file
	
	file<<"This is an object of class CalcData "<<endl;
	file<<ngr<<"\t !number of g(r) data sets ngr"<<endl;
	file<<nsq<<"\t !number of S(Q) data sets nsq"<<endl;
	file<<nfq<<"\t !number of F(Q) data sets nfq"<<endl;
	file << nfg << "\t !number of F(g) data sets nfg" << endl;
	file<<nek<<"\t !number of E(k) data sets nek"<<endl;
		
	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	if(ngr>0)
		file<<"\ng(r) data sets details:"<<endl;
		
	for(i=0;i<ngr;i++)
	{
		px=*(edata.gr_rfinder+i);
		py=*(grfinder+i);
		file << "\nset no. " << i + 1 << "/" << ngr << endl;
		switch (ExptsData::bin_flag[i])
		{
		case 1:
			file << "***WARNING***: The bin size was reset to " << RunParams::rspacing_diff[RunParams::assign_hist[i]] << " Angstrom from " << RunParams::gr_rspacing_ori[i] << " Angstrom to be equal to the r points spacing in the " << ExptsData::datafilename[i] << " file." << endl;
			break;
		case 2:
			file << "***WARNING***: The bin size was reset to " << RunParams::rspacing_diff[RunParams::assign_hist[i]] << " Angstrom from " << RunParams::gr_rspacing_ori[i] << " Angstrom to be equal to the smallest dr in the non-equdistant spacing in the " << ExptsData::datafilename[i] << " file." << endl;
		}
		file<<"Original number of points->"<<*(edata.grsize+i)<<endl;
		file<<"Points used from "<<*(edata.grmin+i)<<" to "<<*(edata.grmax+i)<<\
		", number of points used->"<<*(edata.grused+i)<<endl;
#ifdef _R_SWITCH_GR_MIN_1
		file << "\n r \t ";
		if (r_switch_power == 1)
			file << "r*";
		else
			file << "(r^" << r_switch_power << ")*";
		file << "[g(r) - 1]" << endl;
#else
		file << "\n r \t g(r) " << endl;
#endif
		for(j=0;j<*(edata.grused+i);j++)
		{
			file.precision(6);
			file<<J12<<*px++<<"\t";
			file.precision(12);
			file<<J20<<*py++<<endl;
		}
	}
	
	if(nsq>0)
		file<<"\nS(Q) data sets details:"<<endl;
	
	for(i=0;i<nsq;i++)
	{
		px=*(edata.sq_qfinder+i);
		py=*(sqfinder+i);
		file<<"\nset no. "<<i+1<<"/"<<nsq<<endl;
		file<<"original number of points->"<<*(edata.sqsize+i)<<endl;
		file<<"points used from "<<*(edata.sqmin+i)<<" to "\
		<<*(edata.sqmax+i)<<", number of points used->"<<*(edata.sqused+i)<<endl;
		
#ifdef _MUL_SCAT_VECTOR
		file << "\n Q \t Q*S(Q) " << endl;
#else
		file << "\n Q \t S(Q) " << endl;
#endif
		for(j=0;j<*(edata.sqused+i);j++)
		{
			file.precision(6);
			file<<J12<<*px++<<"\t";
			file.precision(12);
			file<<J20<<*py++<<endl;
		}
	}
	 
	if(nfq>0)	 
		file<<"\nF(Q) data sets details:"<<endl;
	
	for(i=0;i<nfq;i++)
	{
		px=*(edata.fq_qfinder+i);
		py=*(fqfinder+i);
				
		file << "\nset no. " << i + 1 << "/" << nfq;
		if (ExptsData::fqfitIQ[i])
		{
			pi = *(iqfinder + i);
			file << " I(Q) is used, ";
			switch (chi.fit_index[ngr + nsq + i])
			{//hopefully if fault combination was given, it was already dealt with
			case 1:
				file << "fitted parameter a, with constant custom parameters b and alpha";
				break;
			case 32:
				file << "fitted parameter alpha, with constant custom parameters a and b";
				break;
			case 3:
				file << "fitted parameters a and b, with alpha=0";
				break;
			case 33:
				file << "fitted parameters a and alpha, with b=0";
				break;
			case 35:
				file << "fitted parameters a, b and alpha";
			}
		}
		file<< endl;
		file<<"original number of points->"<<*(edata.fqsize+i)<<endl;
		file<<"points used from "<<*(edata.fqmin+i)<<" to "<<*(edata.fqmax+i)<<", number of points used->"<<*(edata.fqused+i)<<endl;

#ifdef _MUL_SCAT_VECTOR
		file<<"\n Q \t Q*F(Q) "<<endl;
#else
		file << "\n Q";
		if (ExptsData::fqfitIQ[i])
			file<<" \t I(Q)"; 
		file<< " \t F(Q)";
		file<< endl;
#endif
		for(j=0;j<*(edata.fqused+i);j++)
		{
			file.precision(6);
			file<<J12<<*px++<<"\t";// Q value
			file.precision(12);
			if (ExptsData::fqfitIQ[i])
			{
				double f_av = 0;
				for (int itype = 0; itype < ntypes; itype++)
				{
					f_av += SimpleCfg::fractions[itype] * (*(edata.fq_sffinder[i * ntypes + itype] + j) + edata.fqfprimeact[i * ntypes + itype]);//<f+f'>
				}
				double f_av2 = f_av * f_av;
				file << J20 << *pi++<<J20 << *py++/f_av2; // I(Q),F(Q) (here the coefficients are not normalized, that is why it has to be divided by f_av2
			}
			else
				file << J20 << *py++; // F(Q) value
			file<< endl;
		}
	}
	if (nfg > 0)
		file << "\nF(g) data sets details:" << endl;

	for (i = 0; i < nfg; i++)
	{
		px = *(edata.fg_gfinder + i);
		py = *(fgfinder + i);
		file << "\nset no. " << i + 1 << "/" << nfg << endl;
		file << "original number of points->" << *(edata.fgsize + i) << endl;
		file << "points used from " << *(edata.fgmin + i) << " to " << *(edata.fgmax + i)\
			<< ", number of points used->"<< *(edata.fgused + i) << endl;
#ifdef _MUL_SCAT_VECTOR
		file << "\n g \t g*F(g) " << endl;
#else
		file << "\n g \t F(g)" << endl;
#endif
		for (j = 0; j < *(edata.fgused + i); j++)
		{
			file.precision(6);
			file << J12 << *px++ << "\t";
			file.precision(12);
			file << J20 << *py++ << endl;//g value and F(g) value
		}
	}

	if(nek>0)	 
		file<<"\nE(k) data sets details:"<<endl;
	
	for(i=0;i<nek;i++)
	{
		px=*(edata.ek_kfinder_ori+i);
		py=*(ekfinder+i);
		file << "\nset no. " << i + 1 << "/" << nek;
		if (ExptsData::ek_ngrid_in[i] > 0)
			file << " calculated with coefficient grid point index " << ExptsData::ek_gridind[i] + ExptsData::ek_gridfrom[i];
		file << endl;
		file<<"original number of points->"<<*(edata.eksize+i)<<endl;
		file<<"points used from "<<*(edata.ekmin+i)<<" to "<<*(edata.ekmax+i)\
		<<", number of points used->"		<<*(edata.ekused+i)<<endl;
		
		file<<"\n k \t E(k)"<<endl; 
		for(j=0;j<*(edata.ekused+i);j++)
		{
			file.precision(6);
			file<<J12<<*px++<<"\t";
			file.precision(12);
			file<<J20<<*py++<<endl;//k value and E(k) value
		}
	}
	file.precision(6);
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}//end of the save function

//==========================save==============
void CalcData::SaveResult(const char *file_name) const
{
	int i,j,index;
	double background;//background correction for S(Q)
	double *px=0;//to avoid warning
	double *py=0;
	double *pexpt=0;
	double *pi = NULL, *pcompt = NULL, *pi_ori = NULL, *pIB = NULL, *pA=NULL;
	double *a,*b,*c,*d,*e,*alpha=NULL;
	char *name, *name2,*conv_numb = NULL;
	ofstream file;
	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10,file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10,fitfilename);
	OpenFile(file,tempfilename,"CalcData::SaveResult",0);//open file, check, whether it was successfully opened
	SetArraysize(&name, NAME_SIZE, "name", "CalcData::SaveResult");
	SetArraysize(&name2, NAME_SIZE, "name2", "CalcData::SaveResult");
	//saves contents in a file
			
	file<<"Calculated and experimental total data for each data sets"<<endl;
	file<<ngr<<"\t g(r) data sets"<<endl;
	file<<nsq<<"\t S(Q) data sets"<<endl;
	file<<nfq<<"\t F(Q) data sets"<<endl;
	file << nfg << "\t F(g) data sets" << endl;
	file<<nek<<"\t E(k) data sets"<<endl;
		
	a=chi.a;//multiplicative factor
	b=chi.b;//constant
	c=chi.c;//linear
	d=chi.d;//quadratic
	e=chi.e;//cubic
	

	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	if(ngr>0)
	{
		px=*edata.gr_rfinder;
		py=*grfinder;
		pexpt=*edata.gr_gfinder;
	}
#ifdef _R_SWITCH_GR_MIN_1
	if (r_switch_power == 1)
		mystrcpy(name, NAME_SIZE, "r*[g(r) - 1]");
	else
	{
		mystrcpy(name, NAME_SIZE, "(r^");
		IntToStr(&conv_numb, r_switch_power);
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, ")*[g(r) - 1]");
	}
#else
	mystrcpy(name, NAME_SIZE, "g(r)");
#endif
	mystrcpy(name2, NAME_SIZE, name);
	mystrcat(name2, NAME_SIZE, " RMC++");

	for(i=0;i<ngr;i++)
	{
		file.precision(6);
		file << "\nTotal "<<name<<" for data set #" << i + 1 << "/" << ngr << endl;
		switch (ExptsData::bin_flag[i])
		{
		case 1:
			file << "***WARNING***: The bin size was reset to " << RunParams::rspacing_diff[RunParams::assign_hist[i]] << " Angstrom from " << RunParams::gr_rspacing_ori[i] << " Angstrom to be equal to the r points spacing in the " << ExptsData::datafilename[i] << " file." << endl;
			break;
		case 2:
			file << "***WARNING***: The bin size was reset to " << RunParams::rspacing_diff[RunParams::assign_hist[i]] << " Angstrom from " << RunParams::gr_rspacing_ori[i] << " Angstrom to be equal to the smallest dr in the non-equdistant spacing in the " << ExptsData::datafilename[i] << " file." << endl;
		}

		
		file<<"Renormalization factor a=\t"<<*a<<endl;
		if (ExptsData::use_cubic[i])
		{
			file<<"Renormalization factor b=\t"<<*b<<endl;
			file<<"Renormalization factor c=\t"<<*c<<endl;
			file<<"Renormalization factor d=\t"<<*d<<endl;
			file<<"Renormalization factor e=\t"<<*e<<endl;
		}
		
		file << J12<<"r"<<"\t"<<J20<<name2<<"\t"<<J20<<"renorm exp. g(r)"<<"\t"<<J20<<"original exp. g(r)";
		if (ExptsData::use_cubic[i])
			file<<"\t"<<J20<<"background";
		file<<endl;
		for(j=0;j<*(edata.grused+i);j++)
		{
			file.precision(6);
			background=(*b + *c * (*px)+ *d * (*px) * (*px) + *e * (*px) * (*px) * (*px));
			file<<J12<<*px++<<"\t";
			file.precision(12);
			if (ExptsData::use_cubic[i])		
				file<<J20<<*py++<<"\t"<<J20<<(*a * *pexpt+ background)<<"\t"<<J20<<*pexpt<<"\t"<<J20<<background<<endl;
			else
				file<<J20<<*py++<<"\t"<<J20<<*a * *pexpt <<"\t"<<J20<<*pexpt<<endl;
			pexpt++;
		}
		
		a++;//go to next
		b++;//go to next
		c++;//go to next
		d++;//go to next
		e++;//go to next
	}
	
	if(nsq>0)
	{
		px=*edata.sq_qfinder;
		py=*sqfinder;
		pexpt=*edata.sq_sfinder;
	}
	
	for(i=0;i<nsq;i++)
	{
		index=ngr+i;
#ifdef _MUL_SCAT_VECTOR
		file<<"\nTotal Q*S(Q) for data set #"<<i+1<<"/"<<nsq<<endl;
#else
		file<<"\nTotal S(Q) for data set #"<<i+1<<"/"<<nsq<<endl;
#endif
		file<<"Renormalization factor a=\t"<<*a<<endl;
		file<<"Renormalization factor b=\t"<<*b<<endl;
		file<<"Renormalization factor c=\t"<<*c<<endl;
		file<<"Renormalization factor d=\t"<<*d<<endl;
		if (ExptsData::use_cubic[index])
			file<<"Renormalization factor e=\t"<<*e<<endl;

#ifdef _MUL_SCAT_VECTOR
		file << J12 << "Q" <<"\t"<<J20<<"Q*S(Q) RMC++"<<"\t"<<J20<<"renorm exp. S(Q)"<<"\t"<<J20<<"original exp. S(Q)"<<"\t"<<J20<<"background" << endl;
#else
		file << J12 << "Q" <<"\t"<<J20<<"S(Q) RMC++"<<"\t"<<J20<<"renorma exp. S(Q)"<<"\t"<<J20<<"original exp. S(Q)"<<"\t"<<J20<<"background" << endl;
#endif
		for(j=0;j<*(edata.sqused+i);j++)
		{
			file.precision(6);
			background=(*b + *c * (*px)+ *d * (*px) * (*px) + *e * (*px) * (*px) * (*px));
			file<<J12<<*px++<<"\t";
			file.precision(12);
			file<<J20<<*py++<<"\t"<<J20<<(*a * *pexpt+ background)<<"\t"<<J20<<*pexpt<<"\t"<<J20<<background<<endl;
			pexpt++;
		}
		
		a++;//go to next 
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next

	}

	alpha = chi.alpha;//alpha parameter in case of I(Q) fitting 
	if(nfq>0)	 
	{
		px=*edata.fq_qfinder;
		py=*fqfinder;
		pexpt=*edata.fq_ffinder;
	}
	
	for(i=0;i<nfq;i++)
	{
		index=ngr+nsq+i;

#ifdef _MUL_SCAT_VECTOR
		file<<"\nTotal Q*F(Q) for data set #"<<i+1<<"/"<<nfq<<endl;
#else
		if (ExptsData::fqfitIQ[i])
			file << "\nTotal I(Q) and F(Q)";
		else
			file << "\nTotal F(Q)";
		
		file<< " for data set #" << i + 1 << "/" << nfq << endl;
		if (ExptsData::fqfitIQ[i])
		{
			pi = *(iqfinder + i);//RMC calc I(Q)
			pcompt = *(comptfinder + i);
			pA = edata.fq_A + edata.fqused_cum[i];
			if (ExptsData::fqIQbackgcorr[i])
			{
				pi_ori = *(edata.fq_ffinder_ori + i);//ori expt I(Q)
				pIB = *(edata.fq_IBfinder + i);//original background
			}
			if (ExptsData::fqAXS[i] > 0)
				file << "AXS ";
			file << " I(Q) is used for fitting, ";
			switch (chi.fit_index[ngr + nsq + i])
			{//hopefully if faulty combination was given, it was already dealt with
			case 1:
				file << "fitted parameter a, with constant custom parameters b and alpha";
				break;
			case 32:
				file << "fitted parameter alpha, with constant custom parameters a and b";
				break;
			case 3:
				file << "fitted parameters a and b, with alpha=0";
				break;
			case 33:
				file << "fitted parameters a and alpha, with b=0";
				break;
			case 35:
				file << "fitted parameters a, b and alpha";
			}
			if (ExptsData::fqIQbackgcorr[i])
				file << ", scalable background correction with actual mu=\t" << edata.fqmuact[i];
			file << endl;
			if (ExptsData::fqAXS[i] > 0)
				file << "For the " << ExptsData::fqAXS[i] << ". type actual f'=\t" << ExptsData::fqfprimeact[i * ntypes + ExptsData::fqAXS[i] - 1] << endl;;
			
			file << "Compton refers to B(Q)*exp(-alpha*Q*Q)";
			if (ExptsData::fqIQbackgcorr[i])
				file << " IB(Q) is the original background supplied in the expt data files";
			file<< endl;
		}
#endif

		file<<"Renormalization factor a=\t"<<*a<<endl;
		file<<"Renormalization factor b=\t"<<*b<<endl;
		file<<"Renormalization factor c=\t"<<*c<<endl;
		file<<"Renormalization factor d=\t"<<*d<<endl;
		if (ExptsData::use_cubic[index])
			file<<"Renormalization factor e=\t"<<*e<<endl;
		if (ExptsData::fqfitIQ[i])
		{
			file << "Renormalization factor alpha=\t" << *alpha << endl;

		}

#ifdef _MUL_SCAT_VECTOR
		file << J12<<"Q"<<"\t"<<J20<<"Q*F(Q) RMC++"<<"\t"<<J20<<"renorm exp. F(Q)"<<"\t"<<J20<<"original exp. F(Q)"<<"\t"<<J20<<"background" << endl;
#else
		file << J12<<"Q"<<"\t";
		if (ExptsData::fqfitIQ[i])
		{
			if (ExptsData::fqIQbackgcorr[i])
				file << J20<<"I(Q) RMC++"<<"\t"<<J23<<"renorm exp. I(Q)_corr"<<"\t"<<J20<<"original exp. I(Q)"<<"\t"<<J20<<"mu_act*IB(Q)"<<"\t"<<J20<<"Compton/a"<<"\t" <<J20<<"<f(Q)^2>/a"<<"\t"<< J20<<"<f(Q)>^2/a"<<"\t"<<J20<<"F(Q) RMC++"<<"\t"<<J20<<"renorm exp. F(Q)";
			else
				file << J20<<"I(Q) RMC++"<<"\t"<<J23<<"renorm exp. I(Q)_corr"<<"\t"<<J20<<"original exp. I(Q)"<<"\t"<<J20<<"Compton/a"<<"\t"<<J20<<"<f(Q)^2>/a"<<"\t"<<J20<<"<f(Q)>^2/a"<<"\t"<<J20<<"F(Q) RMC++"<<"\t"<<J20<<"renorm exp. F(Q)";
		}
		else
			file << J20<<"F(Q) RMC++"<<"\t"<<J23<<"renorm exp. F(Q)"<<"\t"<<J20<<"original exp. F(Q)"<<"\t"<<J20<<"background";
		
		file<< endl;
#endif


		for(j=0;j<*(edata.fqused+i);j++)
		{
			file.precision(6);
			background=(*b + *c * (*px)+ *d * (*px) * (*px)+ *e * (*px) * (*px) * (*px));
			file<<J12<<*px++<<"\t";
			file.precision(12);
			if (ExptsData::fqfitIQ[i])
			{
				//calculate f_av2, only needed for display, as in case of AXS the coefficients can change, it is recalculated 
				double f_av = 0;
				for (int itype = 0; itype < ntypes; itype++)
				{
					f_av += SimpleCfg::fractions[itype] * (*(edata.fq_sffinder[i * ntypes + itype] + j) + edata.fqfprimeact[i * ntypes + itype]);//<f+f'>
				}
				double f_av2 = f_av * f_av;

				if (ExptsData::fqIQbackgcorr[i])
					file << J20 << *pi++ << "\t" << J23 << (*a * *pexpt + background) << "\t" << J20 << *pi_ori++ << "\t" << J20 << ExptsData::fqmuact[i] * *pIB++ << "\t" << J20 << *pcompt/ *a 
						 <<	"\t" << J20 << *pA / *a << "\t" << J20 << f_av2 / *a << "\t" << J20 << *py / f_av2 << "\t" << J20 << ((*a * *pexpt + background) - *pA - *pcompt)/f_av2;
				else
					file << J20 << *pi++ << "\t" << J23 << (*a * *pexpt + background) << "\t" << J20 << *pexpt << "\t" << J20 << *pcompt / *a 
						 << "\t" << J20 << *pA / *a << "\t" << J20 << f_av2 / *a << "\t" << J20 << *py / f_av2 << "\t" << J20 << ((*a * *pexpt + background) - *pA - *pcompt) / f_av2;
					
				pcompt++;
				py++;
				pA++;
			}
			else
				file << J20 << *py++ << "\t" << J23 << (*a * *pexpt + background) << "\t" << J20 << *pexpt << "\t" << J20 << background;
			
			file<< endl;
			pexpt++;
		}
	
		a++;//go to next 
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next
		alpha++;
	}

	if (nfg > 0)
	{
		px = *edata.fg_gfinder;
		py = *fgfinder;
		pexpt = *edata.fg_ffinder;
	}

	for (i = 0; i < nfg; i++)
	{
		index = ngr + nsq + nfq + i;
#ifdef _MUL_SCAT_VECTOR
		file << "\nTotal g*F(g) for data set #" << i + 1 << "/" << nfg << endl;
#else
		file << "\nTotal F(g) for data set #" << i + 1 << "/" << nfg << endl;
#endif

		file << "Renormalization factor a=\t" << *a << endl;
		file << "Renormalization factor b=\t" << *b << endl;
		file << "Renormalization factor c=\t" << *c << endl;
		file << "Renormalization factor d=\t" << *d << endl;
		if (ExptsData::use_cubic[index])
			file << "Renormalization factor e=\t" << *e << endl;

#ifdef _MUL_SCAT_VECTOR
		file <<J12<< "g"<<"\t"<<J20<<"g*F(g) RMC++"<<"\t"<<J20<<"renorm exp. F(g)"<<"\t"<<J20<<"original exp. F(g)"<<"\t"<<J20<<"background" << endl;
#else
		file << J12<<"g"<<"\t"<<J20<<"F(g) RMC++"<<"\t"<<J20<<"renorm exp. F(g)"<<"\t"<<J20<<"original exp. F(g)"<<"\t"<<J20<<"background" << endl;
#endif
		for (j = 0; j < *(edata.fgused + i); j++)
		{
			file.precision(6);
			background = (*b + *c * (*px) + *d * (*px) * (*px) + *e * (*px) * (*px) * (*px));
			file << J12 << *px++ << "\t";
			file.precision(12);
			file << J20 << *py++ << "\t" << J20 << (*a * *pexpt + background) << "\t" << J20 << *pexpt << "\t" << J20 << background << endl;
			pexpt++;
		}

		a++;//go to next 
		b++;//go to next 
		c++;//go to next 
		d++;//go to next 
		e++;//go to next
	}

	if(nek>0)	 
	{
		py=*ekfinder;
		pexpt=*edata.ek_efinder;
	}
	double k_k_max_pow_chiweight,k_max_pow_chiweight;
	
	for (i = 0; i < nek; i++)
	{
		mystrcpy(name, NAME_SIZE, "(k^");
		IntToStr(&conv_numb, edata.ekchipower[i]);
		mystrcat(name, NAME_SIZE, conv_numb);
		mystrcat(name, NAME_SIZE, ")*chi(k) RMC++");
		mystrcpy(name2, NAME_SIZE, "renorm (k^");
		mystrcat(name2, NAME_SIZE, conv_numb);
		mystrcat(name2, NAME_SIZE, ")*chi(k) exp.");
		px = edata.ek_kfinder_ori[i];
		file << "\nTotal E(k) for data set #" << i + 1 << "/" << nek;
		if (ExptsData::ek_ngrid_in[i] > 0)//the k-axis and E0 axis are going in opposite directions
			file << " calculated with coefficient matrix with grid point index " << ExptsData::ek_gridind[i] + ExptsData::ek_gridfrom[i] <<" corresponding to dE0 actual: "<<-(ExptsData::ek_gridind[i] - ExptsData::ek_ngrid_in[i] )* ExptsData::ek_dE0[i] / ExptsData::ek_ngrid_in[i]<<" eV";
		file << endl;
		file<<"Renormalization factor a=\t"<<*a<<endl;
		file<<"Renormalization factor (calculated with k_max^eq_chiweight normalization) b=\t"<<*b<<endl;
		file<<"Applied k_max durinig the chisquare calculation: \t"<<*(edata.ek_kfinder_ori[i]+edata.ekused[i]-1)<<endl;
		file << "BEWARE, that E(k) is fitted, and the chi2 is calculated from E(k) RMC++ and renormalized experimental E(k), 6-7 column" << endl;
		
		file << J12<<"k"<<"\t"<<J20<<name<<"\t"<<J27<<name2 << "\t"<<J20<<"chi(k) RMC++"<<"\t"<<J20<<"original chi(k) exp.";
		file << "\t"<<J20<<"E(k) RMC++"<<"\t"<<J20<<"renorm exp. E(k)"<<"\t"<<J20<<"original exp. E(k)"<<endl;
	
		//for exafs the origial chi(k)=E(k)*(k_max/k)^ekchipower[i] and will be also given for the calculated and original experimental data as well
		// k^CHIK-POWER * chi(k) RMC++ and  k^CHIK-POWER * chi(k) expt 

		k_max_pow_chiweight= pow(*(edata.ek_kfinder_ori[i] + edata.ekused[i] - 1), (double)edata.ekchipower[i]);
		for(j=0;j<*(edata.ekused+i);j++)
		{
			k_k_max_pow_chiweight = pow(*(edata.ek_kfinder_ori[i] + edata.ekused[i] - 1)/ *px, (double)edata.ekchipower[i]);
			file.precision(6);
			file<<J12<<*px;
			file.precision(12);
			file << "\t" << J20 << *py * k_max_pow_chiweight << "\t" << J27 << (*a * *pexpt + *b) * k_max_pow_chiweight;
			file << "\t" << J20 << *py * k_k_max_pow_chiweight << "\t" << J20 << *pexpt * k_k_max_pow_chiweight;
			file << "\t" << J20 << *py <<"\t" << J20 << (*a * *pexpt + *b) << "\t" << J20 << *pexpt<< endl;
			px++;
			pexpt++;
			py++;
		}
		a++;//go to next 
		b++;//go to next 
	}
	delete[] name;
	delete[] name2;
	file.precision(6);
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}//end of the save function

//-----save the totals in RMCA format to *.out file-------------------------
void CalcData::SaveOldOut(ofstream &file, const char *file_name) const
{
	//running will not continue, till the saving could be performed
	CheckFileState(file,"CalcData::SaveOldOut",file_name);// check, whether it was successfully opened
	
	int i,k;
	double background;//background correction for S(Q)
	double *p1,*p2,*p3;//pointers
	double *a, *b, *c, *d, *e;

	file.setf(ios::scientific, ios::floatfield);//formatting
	file.setf(ios::right, ios::adjustfield);

	a=chi.a;//renormalization multiplicative factor
	b=chi.b;//renormalization additive factor, 
	c=chi.c;//linear 
	d=chi.d;//quadratic
	e=chi.e;//cubic

	if(ngr>0)	 
	{
		p1=*(edata.gr_rfinder);//pointer to r values
		p2=*(edata.gr_gfinder);//pointer to g values
		p3=grvalues;//pointer to the RMC calculated g(r)
	}

	//now the comparison between total calculated and experimental g(r) data if any
	for(k=0;k<ngr;k++)//for each g(r) data set
	{
		file<<"ENDGROUP "<<endl;
		file<<"Total g(r) for data set #"<<k+1<<"/"<<ngr<<endl;
		
		file<<"r, g(r) RMC++, renormalised g(r) expt, original g(r) exp ";
		if (ExptsData::use_cubic[k])
			file<<"background";
		file<<endl;
		file<<"CURVES "<<endl;
	
		file<<J12<<*(edata.grused+k);
		file<<" "<<J12<<2<<endl;
							
		for(i=0;i<*(edata.grused+k);i++)//for each data point
		{
			if (ExptsData::use_cubic[k])
			{
				background=(*b + *c * (*p1)+ *d * (*p1) * (*p1) + *e * (*p1) * (*p1) * (*p1));
				file<<J13<<*p1++<<" "<<J13<<*p3++<<" "<<J13<<(*a * *p2 + background) <<" "<<J13<<*p2<<" "<<J13<<background<<endl;
			}
			else
				file<<J13<<*p1++<<" "<<J13<<*p3++<<" "<<J13<<*p2 * *a<<" "<<J13<<*p2<<endl;
			p2++;
		}
		
		a++;//goto next
		b++;//goto next
		c++;//goto next
		d++;//goto next	
		e++;//goto next
			
	}
	
	
	//now the comparison between total calculated and experimental S(Q) data if any
	if(nsq>0)
	{
		p1=*(edata.sq_qfinder);//pointer to Q values
		p2=*(edata.sq_sfinder);//pointer to S(Q) values 
		p3=*(sqfinder);//pointer to the RMC calculated S(Q) 
	}

	for(k=0;k<nsq;k++)//for each S(Q) data set
	{
		file<<"ENDGROUP "<<endl;
		file<<"Total S(Q) for data set #"<<k+1<<"/"<<nsq<<endl;
		
		file<<"Q S(Q) RMC++, renormalised S(Q) expt , original S(Q) expt, background"<<endl;
		file<<"CURVES "<<endl;
	
		file<<J12<<*(edata.sqused+k);
		file<<" "<<J12<<4<<endl;
	
		for(i=0;i<*(edata.sqused+k);i++)//for each data point
		{
			background=(*b + *c * (*p1)+ *d * (*p1) * (*p1) + *e * (*p1) * (*p1) * (*p1));
			file<<J13<<*p1++<<" "<<J13<<*p3++<<" "<<J13<<(*a * *p2 + background) <<" "<<J13<<*p2<<" "<<J13<<background<<endl;
			p2++;

		}
		a++;//goto next
		b++;//goto next
		c++;//goto next
		d++;//goto next	
		e++;//goto next
	} 
	
	//now the comparison between total calculated and experimental F(Q) data if any
	if(nfq>0)
	{
		p1=*(edata.fq_qfinder);//pointer to Q values 
		p2=*(edata.fq_ffinder);//pointer to F(Q) data set
		p3=*(fqfinder);//pointer to the RMC calculated F(Q) 
	}		
	for(k=0;k<nfq;k++)//for each F(Q) data set
	{
		file<<"ENDGROUP "<<endl;
		file<<"Total F(Q) for data set #"<<k+1<<"/"<<nfq<<endl;
		

		file<<"Q, F(Q) RMC++, renormalised F(Q) expt, original F(Q) expt, background"<<endl;
		file<<"CURVES "<<endl;
	
		file<<J12<<*(edata.fqused+k);
		file<<" "<<J12<<2<<endl;
		
		for(i=0;i<*(edata.fqused+k);i++)//for each data point
		{
			background=(*b + *c * (*p1)+ *d * (*p1) * (*p1) + *e * (*p1) * (*p1) * (*p1));
			file<<J13<<*p1++<<" "<<J13<<*p3++<<" "<<J13<<(*a * *p2 + background) <<" "<<J13<<*p2<<" "<<J13<<background<<endl;
			p2++;
		}//next data point
		
		a++;//goto next
		b++;//goto next
		c++;//goto next
		d++;//goto next	
		e++;//goto next
	}//next data set
	
	//now the comparison between total calculated and experimental F(g) data if any
	if (nfg > 0)
	{
		p1 = *(edata.fg_gfinder);//pointer to g values 
		p2 = *(edata.fg_ffinder);//pointer to F(g) data set
		p3 = *(fgfinder);//pointer to the RMC calculated F(g) 
	}
	for (k = 0; k < nfg; k++)//for each F(g) data set
	{
		file << "ENDGROUP " << endl;
		file << "Total F(g) for data set #" << k + 1 << "/" << nfg << endl;


		file << "g, F(g) RMC++, renormalised F(g) expt, original F(g) expt, background" << endl;
		file << "CURVES " << endl;

		file << J12 << *(edata.fgused + k);
		file << " " << J12 << 2 << endl;

		for (i = 0; i < *(edata.fgused + k); i++)//for each data point
		{
			background = (*b + *c * (*p1) + *d * (*p1) * (*p1) + *e * (*p1) * (*p1) * (*p1));
			file << J13 << *p1++ << " " << J13 << *p3++ << " " << J13 << (*a * *p2 + background) << " " << J13 << *p2 << " " << J13 << background << endl;
			p2++;
		}//next data point

		a++;//goto next
		b++;//goto next
		c++;//goto next
		d++;//goto next	
		e++;//goto next
	}//next data set

	//now the comparison between total calculated and experimental E(k) data if any
	if(nek>0)
	{
		p2=*(edata.ek_efinder);//pointer to E(k) data set
		p3=*(ekfinder);//pointer to the RMC calculated E(k) 
	}		
	for(k=0;k<nek;k++)//for each E(k) data set
	{
		p1 = edata.ek_kfinder_ori[k];//pointer to k values 
		file<<"ENDGROUP "<<endl;
		file<<"Total E(k) for data set #"<<k+1<<"/"<<nek<<endl;
		
		file<<"k, E(k) RMC++, renormalised E(k) expt, original  E(k) expt"<<endl;
		file<<"CURVES "<<endl;
		
		file<<J12<<*(edata.ekused+k);
		file<<" "<<J12<<2<<endl;
			
		for(i=0;i<*(edata.ekused+k);i++)//for each data point
		{
			file<<J13<<*p1++<<" "<<J13<<*p3++<<" "<<J13<<(*a * *p2 + *b)<<" "<<J13<<*p2<<endl;
			p2++;
		}//next data point
		
		a++;//goto next
		b++;//goto next
	}//next data set
	
	file.unsetf(ios::scientific);//formatting
	file.unsetf(ios::right);
};

//-------Calculating the totals from the partials---------------	
void CalcData::CalcTotal(ThreadArg &thread_arg)
{
	int iexpt,ir,iq,ik,ig,ipartial,itype;
	int *min_ind, *max_ind;
	double *pg,*ps,*pf,*pe,*pcoeff,*p_part;
	double *pq,*pi,*pA, *pB, *pcompt;
	//here the Threads parameters for no move should be used
	if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift))//only calculate the simulated data in normal steps
	{
		//g(r) constraint
		//min_r_index, max_r_index contains the segmentation 
		min_ind = thread_arg.min_r_index;
		max_ind = thread_arg.max_r_index;
		if (ngr > 0)
		{
			//sets the pointer to the beginning of the coefficients
			pcoeff = edata.gr_coeffs;//coefficients
			for (iexpt = 0; iexpt < ngr; iexpt++)//for each g(r) data set
			{
				//resets the value to zero(needed in the main loop)
				pg = grfinder[iexpt] + *min_ind;
				for (ir = *min_ind; ir <= *max_ind; ir++)
					*pg++ = 0;
				for (ipartial = 0; ipartial < npartials; ipartial++)
				{
					//sets the pointers to the beginning of this thread in the CalcPart grvalues arrays
					p_part = *(calcp.grfinder + iexpt * npartials + ipartial) + *min_ind;
					pg = *(grfinder + iexpt) + *min_ind;//pointer to the g(r) values to be computed
					for (ir = *min_ind; ir <= *max_ind; ir++)//for each data point for this this thread
					{
						*(pg++) += *pcoeff * (*p_part++);

					}//end of ir cycle
					pcoeff++;//sets the coherent coefficient pointer to the next partial
				}//end of ipartial cycle
				min_ind++;
				max_ind++;
			}//iexpt cycle
		}//if there is g(r) data

		//min_Q_g_index, max_Q_g_index contains the segmentation for all Q and g-based data sets!
		min_ind = thread_arg.min_Q_g_index;
		max_ind = thread_arg.max_Q_g_index;

		//S(Q) constraint
		if (nsq > 0)
		{

			//sets the pointer to the beginning of the coefficients
			pcoeff = edata.sq_coeffs;//coefficients

			for (iexpt = 0; iexpt < nsq; iexpt++)
			{
				//resets the value to zero(needed in the main loop)
				ps = sqfinder[iexpt] + *min_ind;
				for (iq = *min_ind; iq <= *max_ind; iq++)
					*ps++ = 0;

				for (ipartial = 0; ipartial < npartials; ipartial++)
				{
					//sets the pointers to the beginning of Q-values for this thread for this partial and data set in the CalcPart sqvalues arrays
					p_part = *(calcp.sqfinder + iexpt * npartials + ipartial) + *min_ind;

					//resets the pointer to the beginning of Q-values for this thread of the iexpt-th data set in the sqvalues array
					ps = sqfinder[iexpt] + *min_ind;

					for (iq = *min_ind; iq <= *max_ind; iq++)
					{
						*(ps++) += *pcoeff * *p_part++;
					}//end of data point cycle iq
					pcoeff++;//sets the coherent coefficient pointer to the next partial
				}//end of partial cycle ipartial
				min_ind++;
				max_ind++;
			}//end of data set cycle iexpt	
		}
	}//if there is atomic move

	//F(Q) constraint
	if (nfq > 0)
	{
		min_ind = thread_arg.min_Q_g_index + nsq;//if no atomic move, it was not set
		max_ind = thread_arg.max_Q_g_index + nsq;
		for (iexpt = 0; iexpt < nfq; iexpt++)
		{
			//only calculate, if nothing is prevented the atomic move, if there is no move only calculate for sets with background corr or AXS f' change
			if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift) || (Threads::is_fprime_shift && edata.fqAXS[iexpt] > 0))
			{
				
				//resets the value to zero(needed in the main loop)
				pf = fqfinder[iexpt] + *min_ind;
				for (iq = *min_ind; iq <= *max_ind; iq++)
					*pf++ = 0;
				
				for (ipartial = 0; ipartial < npartials; ipartial++)
				{
					//sets the pointers to the beginning of Q-values for this thread for this partial and data set in the CalcPart fqvalues arrays
					p_part = *(calcp.fqfinder + iexpt * npartials + ipartial) + *min_ind;
					//sets the pointer to the beginning of Q-dependent coefficients for this thread for this partial and data set in the coefficients array	
					pcoeff = *(edata.fq_cfinder + iexpt * npartials + ipartial) + *min_ind;
					//resets the pointer to the beginning of this thread the iexpt-th data set in the fqvalues array
					pf = fqfinder[iexpt] + *min_ind;

					for (iq = *min_ind; iq <= *max_ind; iq++)
					{
						*(pf++) += *pcoeff++ * *p_part++;
						
					}//end of data point cycle iq
				}//end of partial cycle ipartial
			}//if calculate for this set
			min_ind++;
			max_ind++;
		}//end of data set cycle iexpt	
		if (ExptsData::is_IQ > 0)
		{
			min_ind -= nfq;//set it back to the first F(Q) set
			max_ind -= nfq;
			//calculate the I(Q), where necessary
			for (iexpt = 0; iexpt < nfq; iexpt++)
			{
				//only calculate, if nothing is prevented the atomic move, if there is no move only calculate for sets with background corr
				if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift) || (Threads::is_fprime_shift && edata.fqAXS[iexpt] > 0))
				{
					if (ExptsData::fqfitIQ[iexpt])
					{
						pf = fqfinder[iexpt] + *min_ind;
						pi = iqfinder[iexpt] + *min_ind;
						pcompt = comptfinder[iexpt] + *min_ind;
						pq = edata.fq_qfinder[iexpt] + *min_ind;
						pA = edata.fq_A + edata.fqused_cum[iexpt] + *min_ind;
						pB = edata.fq_B + edata.fqused_cum[iexpt] + *min_ind;

						for (iq = *min_ind; iq <= *max_ind; iq++)
						{
							*pcompt = *pB++ * exp(-chi.alpha[iexpt] * *pq * *pq);
							*(pi++) = *pf++ + *pA++ + *pcompt++;
							pq++;
						}//end of data point cycle iq
					}
				}
				min_ind++;
				max_ind++;
			}//end of data set cycle iexpt	
		}
	}

	if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift))//only calculate the simulated data in normal steps
	{
		//F(g) constraint
		if (nfg > 0)
		{
			min_ind = thread_arg.min_Q_g_index + nsq+nfq;//if no atomic move, it was not set
			max_ind = thread_arg.max_Q_g_index + nsq+nfq;
			for (iexpt = 0; iexpt < nfg; iexpt++)
			{
				//resets the value to zero(needed in the main loop)
				pf = fgfinder[iexpt] + *min_ind;
				for (ig = *min_ind; ig <= *max_ind; ig++)
					*pf++ = 0;

				for (ipartial = 0; ipartial < npartials; ipartial++)
				{
					//sets the pointers to the beginning of g-values for this thread for this partial and data set in the CalcPart fgvalues arrays
					p_part = *(calcp.fgfinder + iexpt * npartials + ipartial) + *min_ind;
					//sets the pointer to the beginning of g-dependent coefficients for this thread for this partial and data set in the coefficients array	
					pcoeff = *(edata.fg_cfinder + iexpt * npartials + ipartial) + *min_ind;
					//resets the pointer to the beginning of this thread the iexpt-th data set in the fgvalues array
					pf = fgfinder[iexpt] + *min_ind;

					for (ig = *min_ind; ig <= *max_ind; ig++)
					{
						*(pf++) += *pcoeff++ * *p_part++;
					}//end of data point cycle ig
				}//end of partial cycle ipartial
				min_ind++;
				max_ind++;
			}//end of data set cycle iexpt	
		}
	}//end of if there is atomic move

	
	//min_k_index, max_k_index contains the segmentation for all k-based data sets!
	min_ind = thread_arg.min_k_index;
	max_ind = thread_arg.max_k_index;

	//E(k) constraint
	if (nek > 0)
	{
		//the partials are already containing the (k,r) dependent coefficients, only the E_i=sum(E_ij) for each j type 
		//has to be calculated, where i is the type of the absorbing particle

		for (iexpt = 0; iexpt < nek; iexpt++)
		{	//only calculate, if nothing is prevented the atomic move, if there is no move only calcuate for the shifted sets
			//here the Threads parameters for no move should be used
			if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift) || (Threads::is_E0_shift && edata.ek_ngrid_in[iexpt] > 0))
			{//calculate if there is no E0 shift in this step, or if there is only for those, where there is
			//resets the value to zero(needed in the main loop)
				pe = ekfinder[iexpt] + *min_ind;//sets the pointers to the beginning of this thread in the ekvalues arrays
				for (ik = *min_ind; ik <= *max_ind; ik++)
					*pe++ = 0;

				for (itype = 0; itype < ntypes; itype++)
				{
					p_part = *(calcp.ekfinder + iexpt * ntypes + itype) + *min_ind;//beginning of k-values for this thread for this partial and data set
					//resets the pointer to the beginning of this thread the iexpt-th data set in the ekvalues array
					pe = ekfinder[iexpt] + *min_ind;

					for (ik = *min_ind; ik <= *max_ind; ik++)
					{
						*(pe++) += *p_part++;
					}//end of data point cycle ik
				}//end of type cycle itype
			}//end of calculate for this set
			min_ind++;
			max_ind++;
		}//end of iext loop
	}//end if no ek set

};

//recalc I(Q) data sets during non-linear regression iteration steps
//threading would have been difficult to handle
void CalcData::RecalcIQ(int i)
{
	int iq;
	double *pf, *pi, *pq, *pA, *pB, *pcompt;
	//F(Q) constraint
	if (nfq > 0)
	{
		//calculate the I(Q), where necessary
		if (chi.fit_index[ngr+nsq+i]>31)//non_lin
		{
			pf = fqfinder[i];
			pi = iqfinder[i];
			pcompt = comptfinder[i];
			pq = edata.fq_qfinder[i];
			pA = edata.fq_A + edata.fqused_cum[i];
			pB = edata.fq_B + edata.fqused_cum[i];
			//logfile << "Calciq " << i + 1 << " alpha " << chi.alpha[i] << " F(Q) "<<*pf<<" IQ_old "<<*pi ;//GO
			for (iq = 0; iq < edata.fqused[i]; iq++)
			{
				//logfile << *pq << " B " << *pB <<" A " << *pA << " f " << *pf ;
				*pcompt = *pB++ * exp(-chi.alpha[i] * *pq * *pq);
				//logfile << " pcomp " << *pcompt;
				*pi = *pf++ + *pA++ + *pcompt++;
				//logfile<< " i " << *pi << endl;
				pi++;
				pq++;
				
			}//end of data point cycle iq
			//logfile << " IQ_new " << *iqfinder[i] << endl;//GO
		}
			
	}//end of data set cycle iexpt	
	
}
//----------copying the totals, needed only for TEST MODE or in case of runs with possible no atomic moves----------------
//as this is not threading, what to copy is decided by the main threads Move:nomove, Move::makefprimeshift, Move::makemucorr, Move::makeE0shift
//values, in case of TEST MODE, for the final saving the Move:nomove is set to false for normal atomic move, indicating the everything has to be copied
//everything has to be saved in case of runs with possible no atomic move steps, when normal atomic move is rejected, as if this is followed by an accepted
//no atomic move saving can occur, and has to contain the last accepted move's data
void CalcData::Copy(CalcData &target)
{
	int i;
	double *p_source,*p_target;

	if (!Move::nomove)//atomic move
	{
		p_source = grvalues;
		p_target = target.grvalues;
		for (i = 0; i < edata.rused_tot; i++)
			*p_target++ = *p_source++;

		p_source = sqvalues;
		p_target = target.sqvalues;
		for (i = 0; i < edata.sqused_tot; i++)
			*p_target++ = *p_source++;
	}
	if (!Move::nomove || Move::makefprimeshift)
	{
		p_source = fqvalues;
		p_target = target.fqvalues;
		for (i = 0; i < edata.fqused_tot; i++)
			*p_target++ = *p_source++;
	}

	if (!Move::nomove || Move::makemucorr || Move::makefprimeshift)
	{
		if (ExptsData::is_IQ > 0)
		{
			p_source = iqvalues;
			p_target = target.iqvalues;
			for (i = 0; i < edata.fqused_tot; i++)
				*p_target++ = *p_source++;
			p_source = compton;
			p_target = target.compton;
			for (i = 0; i < edata.fqused_tot; i++)
				*p_target++ = *p_source++;
		}
	}
	if (!Move::nomove)
	{
		p_source = fgvalues;
		p_target = target.fgvalues;
		for (i = 0; i < edata.fgused_tot; i++)
			*p_target++ = *p_source++;
	}
	if (!Move::nomove || Move::makeE0shift)
	{
		p_source = ekvalues;
		p_target = target.ekvalues;
		for (i = 0; i < edata.ekused_tot; i++)
			*p_target++ = *p_source++;
	}
};
