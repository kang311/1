//source CalcPart.cpp
//Last changed 04.03.2021

#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//classes2.h included through this

int  CalcPart::ntypes;//number of types
int  CalcPart::npartials;//number of partials
int *CalcPart::nbins;////pointer to RunParams::nbins number of bins
#ifdef _VIBR_AMP
	int CalcPart::nbins_used;
#endif
int  CalcPart::ngr;//number of g(r) data sets
int  CalcPart::nsq;//number of S(Q) data sets
int  CalcPart::nfq;//number of F(Q) data set
int  CalcPart::nfg;//number of F(g) data set
int  CalcPart::nek;//number of E(k) data set
int *CalcPart::grused;//array of number of g(r) points/set
int *CalcPart::sqused;//array of number of S(Q) points/set
int *CalcPart::fqused;//array of number of F(Q) points/set
int *CalcPart::fgused;//array of number of F(g) points/set
int *CalcPart::ekused;//array of number of E(k) points/set
	

//----------------constructor with handle initialization----
CalcPart::CalcPart(ExptsData &edat, DataMat &datamat, PPCFSet &ppcfset,  HistoSet &histog, Move &move_obj)
					:edata(edat), datmat(datamat), ppcf(ppcfset), hist(histog), move(move_obj)
{

	//Checking the parameters
	if (npartials==0)
	{
		cout<<"\n"<<"****ERROR****"<<endl;
		cout<<"CalcPart constructor: Number of types is 0!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	if (edata.ngr+edata.nsq+edata.nfq+edata.nfg+edata.nek<=0)
	{
		cout << "\nWARNING(" << ++warn << "): CalcPart constructor: there is no experimental data!"<<endl;
	}
	double *ptemp;
	int i,j;
	
	//computation of the size and memory allocation
	SetArraysize(&grvalues,npartials*edata.rused_tot,"grvalues","CalcPart::CalcPart");
	SetArraysize(&sqvalues,npartials*edata.sqused_tot,"sqvalues","CalcPart::CalcPart");
	SetArraysize(&fqvalues,npartials*edata.fqused_tot,"fqvalues","CalcPart::CalcPart");
	SetArraysize(&fgvalues, npartials*edata.fgused_tot, "fgvalues", "CalcPart::CalcPart");
	SetArraysize(&ekvalues,ntypes*edata.ekused_tot,"ekvalues","CalcPart::CalcPart");//this is different from the others, for each data set only those partials
					//exist, which includes the absorbing particle type
	//finder pointers
	SetArraysize(&grfinder,ngr*npartials,"grfinder","CalcPart::CalcPart");//finder of g(r) partials
	SetArraysize(&sqfinder,nsq*npartials,"sqfinder","CalcPart::CalcPart");//finder of S(Q) partials
	SetArraysize(&fqfinder,nfq*npartials,"fqfinder","CalcPart::CalcPart");//finder of F(Q) partials
	SetArraysize(&fgfinder, nfg*npartials, "fgfinder", "CalcPart::CalcPart");//finder of F(g) partials
	SetArraysize(&ekfinder,nek*ntypes,"ekfinder","CalcPart::CalcPart");//finder of E(k) partials, only ntypes for each data set!
	
	//finder initialisation for each partials 
	for(i=0;i<ngr;i++)//for each g(r) set
	{
		for(j=0;j<npartials;j++)//for each partial
			*(grfinder+i*npartials+j)=grvalues+edata.grused_cum[i]*npartials+*(grused+i)*j;
			//pointer to the g(r) values shifted by the current cumulative size
			//(i.e. to the first g(r) value of the current set) and
			//furthermore shifted to the j-th partial
	}
		
	for(i=0;i<nsq;i++)//for each S(Q) set
	{
		for(j=0;j<npartials;j++)//for each partial
			*(sqfinder+i*npartials+j)=sqvalues+edata.sqused_cum[i]*npartials+*(sqused+i)*j;
			//pointer to the S(Q) values shifted by the current cumulative size
			//(i.e. to the first S(Q) value of the current set) and
			//furthermore shifted to the j-th partial
	}

	for(i=0;i<nfq;i++)//for each F(Q) set
	{
		for(j=0;j<npartials;j++)//for each partial
			*(fqfinder+i*npartials+j)=fqvalues+edata.fqused_cum[i]*npartials+*(fqused+i)*j;
			//pointer to the F(Q) values shifted by the current cumulative size
			//(i.e. to the first F(Q) value of the current set) and
			//furthermore shifted to the j-th partial
	}

	for (i = 0; i < nfg; i++)//for each F(g) set
	{
		for (j = 0; j < npartials; j++)//for each partial
			*(fgfinder + i * npartials + j) = fgvalues + edata.fgused_cum[i] * npartials + *(fgused + i)*j;
		//pointer to the F(g) values shifted by the current cumulative size
		//(i.e. to the first F(g) value of the current set) and
		//furthermore shifted to the j-th partial
	}

	for(i=0;i<nek;i++)//for each E(k) set
	{
		for(j=0;j<ntypes;j++)//for each existing partial
			*(ekfinder+i*ntypes+j)=ekvalues+edata.ekused_cum[i]*ntypes+*(ekused+i)*j;
			//pointer to the E(k) values shifted by the current cumulative size
			//(i.e. to the first E(k) value of the current set) and
			//furthermore shifted to the j-th partial
	}
		
	//filling arrays with zeros
	ptemp=grvalues;
	for(i=0;i<npartials*edata.rused_tot;i++)
		*ptemp++=0.0;
		
	
	ptemp=sqvalues;
	for(i=0;i<npartials*edata.sqused_tot;i++)
		*ptemp++=0.0;
		
	ptemp=fqvalues;
	for(i=0;i<npartials*edata.fqused_tot;i++)
		*ptemp++=0.0;

	ptemp = fgvalues;
	for (i = 0; i < npartials*edata.fgused_tot; i++)
		*ptemp++ = 0.0;

	ptemp=ekvalues;
	for(i=0;i<ntypes*edata.ekused_tot;i++)
		*ptemp++=0.0;
	
};//end of constructor

//----------------copy constructor with handle initialization----
CalcPart::CalcPart(CalcPart &source, ExptsData &edat, DataMat &datamat, PPCFSet &ppcfset,  HistoSet &histog, Move &move_obj)
					:edata(edat), datmat(datamat), ppcf(ppcfset), hist(histog), move(move_obj)
{
	//Checking the parameters
	if (npartials==0)
	{
		cout<<"\n*****ERROR*****"<<endl;
		cout<<"CalcPart constructor: Number of types is 0!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	if (edata.ngr+edata.nsq+edata.nfq+edata.nfg+edata.nek<=0)
	{
		cout << "\nWARNING(" << ++warn << "): CalcPart copy constructor: there is no experimental data!"<<endl;
	}
	double *ptemp,*ptemp2;
	int i,j;
		
	//computation of the size and memory allocation
	SetArraysize(&grvalues,npartials*edata.rused_tot,"grvalues","CalcPart::CalcPart");
	SetArraysize(&sqvalues,npartials*edata.sqused_tot,"sqvalues","CalcPart::CalcPart");
	SetArraysize(&fqvalues,npartials*edata.fqused_tot,"fqvalues","CalcPart::CalcPart");
	SetArraysize(&fgvalues,npartials*edata.fgused_tot,"fgvalues","CalcPart::CalcPart");
	SetArraysize(&ekvalues,ntypes*edata.ekused_tot,"ekvalues","CalcPart::CalcPart");//this is different from the others, for each data set only those partials
	
	//finder pointers
	SetArraysize(&grfinder,ngr*npartials,"grfinder","CalcPart::CalcPart");//finder of g(r) partials
	SetArraysize(&sqfinder,nsq*npartials,"sqfinder","CalcPart::CalcPart");//finder of S(Q) partials
	SetArraysize(&fqfinder,nfq*npartials,"fqfinder","CalcPart::CalcPart");//finder of F(Q) partials
	SetArraysize(&fgfinder,nfg*npartials,"fgfinder","CalcPart::CalcPart");//finder of F(g) partials
	SetArraysize(&ekfinder,nek*ntypes,"ekfinder","CalcPart::CalcPart");//finder of E(k) partials, only ntypes for each data set!
	
	//finder initialisation for each partials 
	for(i=0;i<ngr;i++)//for each g(r) set
	{
		for(j=0;j<npartials;j++)//for each partial
			*(grfinder+i*npartials+j)=grvalues+edata.grused_cum[i]*npartials+*(grused+i)*j;
			//pointer to the g(r) values shifted by the current cumulative size
			//(i.e. to the first g(r) value of the current set) and
			//furthermore shifted to the j-th partial
	}
		
	for(i=0;i<nsq;i++)//for each S(Q) set
	{
		for(j=0;j<npartials;j++)//for each partial
			*(sqfinder+i*npartials+j)=sqvalues+edata.sqused_cum[i]*npartials+*(sqused+i)*j;
			//pointer to the S(Q) values shifted by the current cumulative size
			//(i.e. to the first S(Q) value of the current set) and
			//furthermore shifted to the j-th partial
	}

	for(i=0;i<nfq;i++)//for each F(Q) set
	{
		for(j=0;j<npartials;j++)//for each partial
			*(fqfinder+i*npartials+j)=fqvalues+edata.fqused_cum[i]*npartials+*(fqused+i)*j;
			//pointer to the F(Q) values shifted by the current cumulative size
			//(i.e. to the first F(Q) value of the current set) and
			//furthermore shifted to the j-th partial
	}
	
	for (i = 0; i < nfg; i++)//for each F(g) set
	{
		for (j = 0; j < npartials; j++)//for each partial
			*(fgfinder + i * npartials + j) = fgvalues + edata.fgused_cum[i] * npartials + *(fgused + i)*j;
		//pointer to the F(g) values shifted by the current cumulative size
		//(i.e. to the first F(g) value of the current set) and
		//furthermore shifted to the j-th partial
	}

	for(i=0;i<nek;i++)//for each E(k) set
	{
		for(j=0;j<ntypes;j++)//for each existing partial
			*(ekfinder+i*ntypes+j)=ekvalues+edata.ekused_cum[i]*ntypes+*(ekused+i)*j;
			//pointer to the E(k) values shifted by the current cumulative size
			//(i.e. to the first E(k) value of the current set) and
			//furthermore shifted to the j-th partial
	}
		
	//copying arrays
	ptemp=grvalues;
	ptemp2=source.grvalues;
	for(i=0;i<npartials*edata.rused_tot;i++)
		*ptemp++=*ptemp2++;

	ptemp=sqvalues;
	ptemp2=source.sqvalues;
	for(i=0;i<npartials*edata.sqused_tot;i++)
		*ptemp++=*ptemp2++;

	ptemp=fqvalues;
	ptemp2=source.fqvalues;
	for(i=0;i<npartials*edata.fqused_tot;i++)
		*ptemp++=*ptemp2++;

	ptemp = fgvalues;
	ptemp2 = source.fgvalues;
	for (i = 0; i < npartials*edata.fgused_tot; i++)
		*ptemp++ = *ptemp2++;

	ptemp=ekvalues;
	ptemp2=source.ekvalues;
	for(i=0;i<ntypes*edata.ekused_tot;i++)
		*ptemp++=*ptemp2++;
};//end of constructor

//-----------Set static parameters----------------------
void CalcPart::SetCalcPartParams()
{
	ntypes=RunParams::ntypes;
	npartials=ntypes*(ntypes+1)/2;//number of partials
	nbins=RunParams::nbins;//pointer to number of bins
#ifdef _VIBR_AMP
	nbins_used=HistoSet::nbins_used;
#endif
	ngr=ExptsData::ngr;//number of g(r) data sets
	nsq=ExptsData::nsq;//number of S(Q) data sets
	nfq=ExptsData::nfq;//number of F(Q) data sets
	nfg = ExptsData::nfg;//number of F(g) data sets
	nek=ExptsData::nek;//number of E(k) data sets
	//number of used data points /data set
	grused=ExptsData::grused;
	sqused=ExptsData::sqused;
	fqused=ExptsData::fqused;
	fgused = ExptsData::fgused;
	ekused=ExptsData::ekused;

}


//-----save the g(r) partials to a file-------------------------
void CalcPart::Savegr() const
{
	double *p1,*p2,*p3;
	int i,j,k;
	ofstream file;
	
	OpenFile(file,pgrfilename,"CalcPart::Savegr",0);//open file, check, whether it was successfully opened
	
	file<<J10<<ngr<<" g(r) data set(s)"<<endl;
	file<<J10<<npartials<<" partial g(r)-s"<<endl;
	file<<endl;	
		
	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	p1=edata.gr_rvalues;//initialisation of pointer to r values
	for(i=0;i<ngr;i++)//for each g(r) set
	{
#ifdef _R_SWITCH_GR_MIN_1
		file << "\ng(r) set " << i + 1 << "/" << ngr << " r value + " << npartials << " partial ";
		if (r_switch_power == 1)
			file << "r*";
		else
			file << "(r^" << r_switch_power << ")*";
		file << "[g(r) - 1]" << endl;
#else
		file << "\ng(r) set " << i + 1 << "/" << ngr << " r value + " << npartials << " partial g(r)" << endl;
		
#endif
		switch (ExptsData::bin_flag[i])
		{
		case 1:
			file << "***WARNING***: The bin size was reset to " << RunParams::rspacing_diff[RunParams::assign_hist[i]] << " Angstrom from " << RunParams::gr_rspacing_ori[i] << " Angstrom to be equal to the r points spacing in the " << ExptsData::datafilename[i] << " file." << endl;
			break;
		case 2:
			file << "***WARNING***: The bin size was reset to " << RunParams::rspacing_diff[RunParams::assign_hist[i]] << " Angstrom from " << RunParams::gr_rspacing_ori[i] << " Angstrom to be equal to the smallest dr in the non-equdistant spacing in the " << ExptsData::datafilename[i] << " file." << endl;
		}
		for(k=0;k<*(grused+i);k++)//for each data point
		{
			file.precision(6);
			file<<J10<<*p1<<"\t";//that's the r position
			p1++;
			file.precision(12);
			for(j=0;j<npartials;j++)//for each partial
			{
				//pointer to g(r) values
				p2=*(grfinder+i*npartials+j);//that's the position of the
				//first g(r) for the j-th partial of the i-th set
				p3=p2+k;//that's the position of th k-th g(r) value 
				file<<J10<<*p3<<"\t";
			}//next partial
			file<<endl;
		}//next r value
	}//next g(r) data set
	file.precision(6);
	file.unsetf(ios::right);
	file.unsetf(ios::scientific);
	file.close();
}

//-----save the S(Q) partials to a file-------------------------
void CalcPart::SaveSQ() const
{
	double *p1,*p2,*p3;
	int i,j,k;
	ofstream file;
	
	OpenFile(file,psqfilename,"CalcPart::SaveSQ",0);//open file, check, whether it was successfully opened

	file<<J10<<nsq<<" S(Q) data set(s)"<<endl;
	file<<J10<<npartials<<" partial S(Q)-s"<<endl;
	file<<endl;	

	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	p1=edata.sq_qvalues;//initialisation of pointer to Q values
	for(i=0;i<nsq;i++)//for each S(Q) set
	{
#ifdef _MUL_SCAT_VECTOR
		file<<"\nS(Q) set "<<i+1<<"/"<<nsq<<" Q value + "<<npartials<<" partial Q*S(Q)"<<endl;
#else
		file<<"\nS(Q) set "<<i+1<<"/"<<nsq<<" Q value + "<<npartials<<" partial S(Q)"<<endl;
#endif
		
		
		for(k=0;k<*(sqused+i);k++)//for each data point
		{
			file.precision(6);
			file<<J10<<*p1<<"\t";//that's the Q position
			p1++;
			file.precision(12);
			for(j=0;j<npartials;j++)//for each partial
			{
				//pointer to S(Q) values
				p2=*(sqfinder+i*npartials+j);//that's the position of the
				//first S(Q) for the j-th partial of the i-th set
				p3=p2+k;//that's the position of th k-th S(Q) value 
				file<<J10<<*p3<<"\t";
			}//next partial
			file<<endl;
		}//next Q value
	}//next S(Q) data set
	file.precision(6);
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}
	
//-----save the F(Q) partials to a file-------------------------
void CalcPart::SaveFQ() const
{
	double *p1,*p2,*p3;
	int i,j,k;
	ofstream file;
	
	OpenFile(file,pfqfilename,"CalcPart::SaveFQ",0);//open file, check, whether it was successfully opened

	file<<J10<<nfq<<" F(Q) data set(s)"<<endl;
	file<<J10<<npartials<<" partial F(Q)-s"<<endl;
	file<<endl;	

	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	p1=edata.fq_qvalues;//initialisation of pointer to Q values
	for(i=0;i<nfq;i++)//for each F(Q) set
	{
#ifdef _MUL_SCAT_VECTOR
		file<<"\nF(Q) set "<<i+1<<"/"<<nfq<<" Q value + "<<npartials<<" partial Q*F(Q)"<<endl;
#else
		file<<"\nF(Q) set "<<i+1<<"/"<<nfq<<" Q value + "<<npartials<<" partial F(Q)"<<endl;
#endif
				
		for(k=0;k<*(fqused+i);k++)//for each data point
		{
			file.precision(6);
			file<<J10<<*p1<<"\t";//that's the Q position
			p1++;
			file.precision(12);
			for(j=0;j<npartials;j++)//for each partial
			{
				//pointer to F(Q) values
				p2=*(fqfinder+i*npartials+j);//that's the position of the
				//first F(Q) for the j-th partial of the i-th set
				p3=p2+k;//that's the position of the k-th F(Q) value 
				file<<J10<<*p3<<"\t";
			}//next partial
			file<<endl;
		}//next Q value
	}//next F(Q) data set
	file.precision(6);
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}//end of save function

//-----save the F(g) partials to a file-------------------------
void CalcPart::SaveFg() const
{
	double *p1, *p2, *p3;
	int i, j, k;
	ofstream file;

	OpenFile(file, pfgfilename, "CalcPart::SaveFg", 0);//open file, check, whether it was successfully opened

	file << J10 << nfg << " F(g) data set(s)" << endl;
	file << J10 << npartials << " partial F(g)-s" << endl;
	file << endl;

	file.setf(ios::right, ios::adjustfield);
	file.setf(ios::scientific, ios::floatfield);
	p1 = edata.fg_gvalues;//initialisation of pointer to g values
	for (i = 0; i < nfg; i++)//for each F(g) set
	{
#ifdef _MUL_SCAT_VECTOR
		file << "\nF(g) set " << i + 1 << "/" << nfg << " g value + " << npartials << " partial g*F(g)" << endl;
#else
		file << "\nF(g) set " << i + 1 << "/" << nfg << " g value + " << npartials << " partial F(g)" << endl;
#endif

		for (k = 0; k < *(fgused + i); k++)//for each data point
		{
			file.precision(6);
			file << J10 << *p1 << "\t";//that's the g position
			p1++;
			file.precision(12);
			for (j = 0; j < npartials; j++)//for each partial
			{
				//pointer to F(g) values
				p2 = *(fgfinder + i * npartials + j);//that's the position of the
				//first F(g) for the j-th partial of the i-th set
				p3 = p2 + k;//that's the position of the k-th F(g) value 
				file << J10 << *p3 << "\t";
			}//next partial
			file << endl;
		}//next g value
	}//next F(g) data set
	file.precision(6);
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}//end of save function

//-----save the E(k) partials to a file-------------------------
void CalcPart::SaveEK() const
{
	double *p1,*p2,*p3;
	int i,j,k;
	ofstream file;
	
	OpenFile(file,pekfilename,"CalcPart::SaveEK",0);//open file, check, whether it was successfully opened

	file<<J10<<nek<<" E(k) data set(s)"<<endl;
	file<<J10<<ntypes<<" partial E(k)-s containing the absorbing particle"<<endl;
	file<<endl;	

	file.setf(ios::right,ios::adjustfield);
	file.setf(ios::scientific,ios::floatfield);
	
	for(i=0;i<nek;i++)//for each E(k) set
	{
		p1 = edata.ek_kfinder_ori[i];//initialisation of pointer to k values
		file<<"\nE(k) set "<<i+1<<"/"<<nek<<\
		" k value + "<<ntypes<<" partial E(k) containing the absorbing particle"<<endl;
		
		for(k=0;k<*(ekused+i);k++)//for each data point
		{
			file.precision(6);
			file<<J10<<*p1<<"\t";//that's the k position
			p1++;
			file.precision(12);
			for(j=0;j<ntypes;j++)//for each partial
			{
				//pointer to E(k) values
				p2=*(ekfinder+i*ntypes+j);//that's the position of the
				//first E(k) for the j-th type of the i-th set
				p3=p2+k;//that's the position of the k-th E(k) value 
				file<<J10<<*p3<<"\t";
			}//next type
			file<<endl;
		}//next k value
	}//next E(k) data set
	file.precision(6);
	file.unsetf(ios::scientific);
	file.unsetf(ios::right);
	file.close();
}//end of save function

//-----save the partials in RMCA format to *.fit file-------------------------
void CalcPart::SaveOldOut(ofstream &file, const char *file_name) const
{
	double *p1,*p2,*p3;
	int i,j,k;

	//running will not continue, till the saving could be performed
	CheckFileState(file,"CalcPart::SaveOldOut",file_name);// check, whether it was successfully opened
		
	file.setf(ios::scientific,ios::floatfield);//formatting
	file.setf(ios::right,ios::adjustfield);

	//now the partial S(Q)-s for each data set
	p1=edata.sq_qvalues;//initialisation of pointer to Q values
	for(i=0;i<nsq;i++)//for each S(Q) set
	{
		file<<"ENDGROUP "<<endl;
		file<<"S(Q) partials for data set #"<<i+1<<"/"<<nsq<<endl;
		file<<"Q, S(Q) partials IN RMCA ORDER! "<<endl;
		file<<"PLOTS "<<endl;
		file<<J10<<*(sqused+i);
		file<<" "<<J10<<npartials<<endl;
			
		for(k=0;k<*(sqused+i);k++)//for each data point
		{
			file<<J13<<*p1++;//that's the q-value
			for (j=0;j<npartials;j++)//for each partial
			{
				//pointer to S(Q) values
				p2=*(sqfinder+i*npartials+j);//that's the position of the
				//first S(Q) for the j-th partial of the i-th set
				p3=p2+k;//that's the position of th k-th S(Q) value 
				file<<" "<<J13<<*p3;
			}//next partial
			file<<endl;
		}
	} 

	//now the partial F(Q)-s for each data set
	p1=edata.fq_qvalues;//initialisation of pointer to Q values
	for(i=0;i<nfq;i++)//for each F(Q) set
	{
		file<<"ENDGROUP "<<endl;
		file<<"F(Q) partials for data set #"<<i+1<<"/"<<nfq<<endl;
		file<<"Q, F(Q) partials IN RMCA ORDER! "<<endl;
		file<<"PLOTS "<<endl;
		file<<J10<<*(fqused+i);
		file<<" "<<J10<<npartials<<endl;
		
		for(k=0;k<*(fqused+i);k++)//for each data point
		{
			file<<J13<<*p1++;//that's the q-value
			for (j=0;j<npartials;j++)//for each partial
			{
				//pointer to F(Q) values
				p2=*(fqfinder+i*npartials+j);//that's the position of the
				//first F(Q) for the j-th partial of the i-th set
				p3=p2+k;//that's the position of th k-th F(Q) value 
				file<<" "<<J13<<*p3;
			}//next partial
			file<<endl;
		}
	} 

	//now the partial F(g)-s for each data set
	p1 = edata.fg_gvalues;//initialisation of pointer to g values
	for (i = 0; i < nfg; i++)//for each F(g) set
	{
		file << "ENDGROUP " << endl;
		file << "F(g) partials for data set #" << i + 1 << "/" << nfg << endl;
		file << "g, F(g) partials IN RMCA ORDER! " << endl;
		file << "PLOTS " << endl;
		file << J10 << *(fgused + i);
		file << " " << J10 << npartials << endl;

		for (k = 0; k < *(fgused + i); k++)//for each data point
		{
			file << J13 << *p1++;//that's the g-value
			for (j = 0; j < npartials; j++)//for each partial
			{
				//pointer to F(g) values
				p2 = *(fgfinder + i * npartials + j);//that's the position of the
				//first F(g) for the j-th partial of the i-th set
				p3 = p2 + k;//that's the position of th k-th F(g) value 
				file << " " << J13 << *p3;
			}//next partial
			file << endl;
		}
	}

	//now the partial E(k)-s containing the absorbing particle for each data set
	p1=edata.ek_kvalues;//initialisation of pointer to k values
	for(i=0;i<nek;i++)//for each E(k) set
	{
		file<<"ENDGROUP "<<endl;
		file<<"E(k) partials containing the absorbing particle for data set #"<<i+1<<"/"<<nek<<endl;
		file<<"k, E(k) partials IN RMCA ORDER! "<<endl;
		file<<"PLOTS "<<endl;
		file<<J10<<*(ekused+i);
		file<<" "<<J10<<ntypes<<endl;
		
		for(k=0;k<*(ekused+i);k++)//for each data point
		{
			file<<J13<<*p1++;//that's the k-value
			for (j=0;j<ntypes;j++)//for each partial
			{
				//pointer to E(k) values
				p2=*(ekfinder+i*ntypes+j);//that's the position of the
				//first E(k) for the j-th type of the i-th set
				p3=p2+k;//that's the position of the k-th E(k) value 
				file<<" "<<J13<<*p3;
			}//next type
			file<<endl;
		}
	} 

	file.unsetf(ios::scientific);//formatting
	file.unsetf(ios::right);
};

//-----Calculate the partial pair correlation functions for the same data points,-------
//--------------------------- as the experimental, in case of g(r) fitting----------------
void CalcPart::CalcPartialgr(ThreadArg &thread_arg)
{
	int iexpt,ipartial,ir,ibin,mybins,offset;//cycle variables
	int *min_ind,*max_ind;//cycle boundaries
	double *p_matrix,*p_ppcf,*p_calcpart;//pointer to the conversion matrix, and  g(r) arrays
#ifdef _R_SWITCH_GR_MIN_1
	double *pr;//pointer to r values
#endif


	//min_r_index, max_r_index contains the segmentation 
	min_ind=thread_arg.min_r_index;
	max_ind=thread_arg.max_r_index;
	for (iexpt=0; iexpt<edata.ngr; iexpt++)//for each data set
	{
#ifdef _VIBR_AMP
		mybins = nbins_used;
		offset = 0;
#else
		mybins = nbins[RunParams::assign_hist[iexpt]];//set the histogram bins for this data set, if binshift greater 1, the first few bin will not be used, but it will be handled later
		offset = RunParams::assign_hist[iexpt] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppsf
#endif
		mybins -= RunParams::firstbin[iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
		for (ipartial=0; ipartial<npartials; ipartial++)//for each partial
		{
			//positions the pointer to the beginning of the iext-th experimet,
			//ipartial-th partial in the CalcPart grvalues, *min_ind data point
			p_calcpart=*(grfinder+iexpt*npartials+ipartial)+ *min_ind;			
			//positions the pointer to the beginning for this thread of the iexpt-th bin->dr 
			//conversion matrix of the DataMat object, has to be set for each r point, it was calculated only for the bins used by this set
			p_matrix = datmat.grfinder[iexpt] + *min_ind * mybins;
			
#ifdef _R_SWITCH_GR_MIN_1
			pr = edata.gr_rvalues + *min_ind;//initialisation of pointer to r values
#endif
		
			for (ir= *min_ind; ir<= *max_ind; ir++)//for each r value
			{
				//positions the pointer to the beginning of the PPCFSet's g(r) values,
				//for the ipartial-th partial
				p_ppcf=ppcf.finder[offset+ipartial]+ RunParams::firstbin[iexpt];
				for (ibin= 0; ibin<mybins; ibin++)//for each used bin for this dataset
				{
					//cout<<"ir "<<ir<<" ibin "<<ibin<<" p matrix "<<p_matrix<<endl;
					*p_calcpart+=(*p_matrix++)*(*p_ppcf++);//convert bins to r data point

				}//end of bin cycle ibin
#ifdef _R_SWITCH_GR_MIN_1
				*p_calcpart -= 1.0;//to have (r^r_switch_power)*[g(r)-1]
				*p_calcpart *= pow(*pr, r_switch_power);
				pr++;

#endif
				p_calcpart++;//goes to the next r data point in CalcPart g(r) array
			}//end of r data point cycle ir
		}//end of partial cycle ipartial
		min_ind++;
		max_ind++;
	}//end of data set cycle iexpt 
};

	
//calculate the partial S(Q) functions by Fourier-transformation in case of S(Q) fitting
void CalcPart::CalcPartialSq(ThreadArg &thread_arg)
{
	int iexpt,ipartial,ir,iq,mybins,offset;//cycle variables
	int *min_ind,*max_ind;//cycle boundaries
	double *p_ppcf;//pointer to the partial g(r) arrays
	double *psqvalue;//pointer to the values of the CalcData::sqvalues array    
	double *pfmatrix;//pointer to the Fourier  matrices in DataMat object
	//min_Q_g_index, max_Q_g_index contains the segmentation for all Q-based data sets!
	min_ind=thread_arg.min_Q_g_index;
	max_ind=thread_arg.max_Q_g_index;		
	for (iexpt=0;iexpt<nsq;iexpt++)
	{
#ifdef _VIBR_AMP
		mybins = nbins_used;
		offset = 0;
#else
		mybins = nbins[RunParams::assign_hist[iexpt+ngr]];//set the binsize for this data set
		offset = RunParams::assign_hist[iexpt+ngr] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppcf
#endif
		
		mybins -= RunParams::firstbin[ngr+iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
		for (ipartial=0; ipartial<npartials; ipartial++)//for each partial
		{
			psqvalue=*(sqfinder+iexpt*npartials+ipartial)+ *min_ind;//beginning of Q-values for this thread	for this partial and data set
			pfmatrix=*(datmat.sqfinder+iexpt)+ *min_ind * mybins;//sets the pointer to the first Q value for this thread of sqfinder of the iexpts Fourier matrix
			for (iq= *min_ind;iq<= *max_ind;iq++)
			{
				//sets the pointer to the beginning of the given partial 
				//in the PPCFSet gvalues array
				p_ppcf=ppcf.finder[offset+ipartial]+ RunParams::firstbin[ngr + iexpt];
				
				for (ir=0;ir<mybins;ir++)
				{
					*psqvalue+=(*p_ppcf-1) * *pfmatrix++;//summing 
					p_ppcf++;
		
				}//end of data point cycle ir	
				psqvalue++;//next S(Q) data point
			}//end of Q points cycle iq
		}//end of partial cycle ipartial
		min_ind++;
		max_ind++;
	}//end of data set cycle iexpt
};

//calculate the partial F(Q) functions by Fourier-transformation in case of F(Q) fitting
void CalcPart::CalcPartialFq(ThreadArg &thread_arg)
{
	int iexpt,ipartial,ir,iq,mybins,offset;//cycle variables
	int *min_ind,*max_ind;//cycle boundaries
	double *p_ppcf;//pointer to the partial g(r) arrays
	double *pfqvalue;//pointer to the values of the CalcData::fqvalues array    
	double *pfmatrix;//pointer to the Fourier  matrices in DataMat object

	//min_Q_g_index, max_Q_g_index contains the segmentation for all Q-based data sets!
	min_ind=thread_arg.min_Q_g_index+nsq;
	max_ind=thread_arg.max_Q_g_index+nsq;
		
	for (iexpt=0;iexpt<nfq;iexpt++)
	{
#ifdef _VIBR_AMP
		mybins = nbins_used;
		offset = 0;
#else
		mybins = nbins[RunParams::assign_hist[iexpt + ngr+nsq]];//set the binsize for this data set
		offset = RunParams::assign_hist[iexpt + ngr+nsq] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppcf
#endif
		mybins -= RunParams::firstbin[ngr+nsq+iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
		for (ipartial=0; ipartial<npartials; ipartial++)//for each partial
		{
			pfqvalue=*(fqfinder+iexpt*npartials+ipartial)+ *min_ind;//beginning of Q-values for this thread	for this partial and data set
			pfmatrix=*(datmat.fqfinder+iexpt)+ *min_ind * mybins;//sets the pointer to the first Q value for this thread of fqfinder of the iexpts Fourier matrix
			
			for (iq= *min_ind;iq<= *max_ind;iq++)
			{
				//sets the pointer to the beginning of the given partial 
				//in the PPCFSet gvalues array
				p_ppcf=ppcf.finder[offset+ipartial]+ RunParams::firstbin[ngr + nsq+iexpt];
				
				for (ir=0;ir<mybins;ir++)
				{
					*pfqvalue+=(*p_ppcf-1) * *pfmatrix++;//summing 
					p_ppcf++;
		
				}//end of data point cycle ir	
				pfqvalue++;//next F(Q) data point
			}//end of Q points cycle iq
		}//end of partial cycle ipartial
		min_ind++;
		max_ind++;
	}//end of data set cycle iexpt
};

//calculate the partial F(g) functions by Fourier-transformation in case of F(g) fitting
void CalcPart::CalcPartialFg(ThreadArg &thread_arg)
{
	int iexpt, ipartial, ir, ig, mybins,offset;//cycle variables
	int *min_ind, *max_ind;//cycle boundaries
	double *p_ppcf;//pointer to the partial g(r) arrays
	double *pfgvalue;//pointer to the values of the CalcData::fgvalues array    
	double *pfmatrix;//pointer to the Fourier  matrices in DataMat object

	//min_g_index, max_g_index contains the segmentation for all g-based data sets!
	min_ind = thread_arg.min_Q_g_index + nsq + nfq;
	max_ind = thread_arg.max_Q_g_index + nsq + nfq;

	for (iexpt = 0; iexpt < nfg; iexpt++)
	{
#ifdef _VIBR_AMP
		mybins = nbins_used;
		offset = 0;
#else
		mybins = nbins[RunParams::assign_hist[iexpt + ngr+nsq+nfq]];//set the binsize for this data set
		offset = RunParams::assign_hist[iexpt + ngr + nsq+nfq] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppcf
#endif
		mybins -= RunParams::firstbin[ngr+nsq+nfq+iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
		for (ipartial = 0; ipartial < npartials; ipartial++)//for each partial
		{
			pfgvalue = *(fgfinder + iexpt * npartials + ipartial) + *min_ind;//beginning of g-values for this thread	for this partial and data set
			pfmatrix = *(datmat.fgfinder + iexpt) + *min_ind * mybins;//sets the pointer to the first g value for this thread of fgfinder of the iexpts Fourier matrix

			for (ig = *min_ind; ig <= *max_ind; ig++)
			{
				//sets the pointer to the beginning of the given partial 
				//in the PPCFSet gvalues array
				p_ppcf = ppcf.finder[offset + ipartial]+ RunParams::firstbin[ngr + nsq+nfq+iexpt];

				for (ir = 0; ir < mybins; ir++)
				{
					*pfgvalue += (*p_ppcf - 1) * *pfmatrix++;//summing 
					p_ppcf++;

				}//end of data point cycle ir	
				pfgvalue++;//next F(g) data point
			}//end of g points cycle ig
		}//end of partial cycle ipartial
		min_ind++;
		max_ind++;
	}//end of data set cycle iexpt
};


//calculate the partial E(k) functions by multiplying the histogram by ek_coeffs for E(k) fitting
//the coefficients are given for the bins used by the set, not for bins which are left out according to binshift
//if _NO_PERIODIC option is used and EXAFS has to be calculated, then periodic_histogram is used
void CalcPart::CalcPartialEk(ThreadArg &thread_arg)
{
	int iexpt,itype,ir,ik;//cycle variables
	int *min_ind,*max_ind;//cycle boundaries
	int partialind;//index of the partial histogram  
	int myfactor;//to take into account that the histogram contains only the unique pairs, and for the pure partials should be multiplied by 2 to cnsider all atoms as central
	int abs_type;//type of the absorbing particle
	double *pekvalue;//pointer to the values of the CalcData::ekvalues array    
	double *pcoeff;//pointer to the (k,r) dependent coefficient matrices in ExptsData object
	bool all_zero;//to check, if there is a non-zero data point, as it can happen, that if the coeff range and the histogram range does not coincide

#ifdef _NO_PERIODIC
#ifndef _VIBR_AMP
	double *pd_hist;
	pd_hist=0;
#else
	longint *phist;//pointer to the partial histogram arrays
	phist = 0;
#endif
#else
#ifndef _VIBR_AMP
	int offset;//to the beginning of the given data set's beginning
#endif
	longint *phist;//pointer to the partial histogram arrays
	phist=0;
#endif		
	//min_k_index, max_k_index contains the segmentation for all k-based data sets!
	min_ind=thread_arg.min_k_index;
	max_ind=thread_arg.max_k_index;
		

	for (iexpt=0;iexpt<nek;iexpt++)
	{
		all_zero = true;
		thread_arg.ek_range_error[iexpt] = 0;//no error by default
#ifndef _VIBR_AMP
#ifndef _NO_PERIODIC
		offset = npartials*RunParams::assign_hist[iexpt + ngr + nsq + nfq + nfg];//This gives the index to the first partial of this data set's histtogram of the hist finder
#endif
#endif
		abs_type=edata.ek_abstype[iexpt];//type of the absorbing particle
		for (itype=0; itype<ntypes; itype++)//for each partial containing the absorbing particle type
		{
			//Determining the index of the partial histogram  
			partialind=(itype<=abs_type ? (itype*ntypes-(itype*(itype+1)/2)+ abs_type) : \
						(abs_type*ntypes-(abs_type*(abs_type+1)/2)+ itype));

			pekvalue=*(ekfinder+iexpt*ntypes+itype)+ *min_ind;//beginning of E(k)-values for this thread	for this partial and data set
			pcoeff=*(edata.ek_cfinder+ edata.ek_cf_cum[iexpt]+edata.ek_gridind[iexpt]*ntypes+itype) + *min_ind * edata.ek_rused[iexpt*ntypes+itype];//sets the pointer to the first coefficient of 
			myfactor = 1;
			if (abs_type == itype)
				myfactor = 2;
			for (ik= *min_ind;ik<= *max_ind;ik++)
			{
				
				//sets the pointer to the first used bin of the given partial 
				//in the histogram
#ifdef _VIBR_AMP
				phist=+*(hist.finder_T+partialind)+edata.ek_rmin[iexpt*ntypes+itype]-1+RunParams::firstbin[ngr + nsq + nfq + nfg + iexpt];
#else
				
#ifdef _NO_PERIODIC
				pd_hist=*(hist.periodic_finder+partialind)+edata.ek_rmin[iexpt*ntypes+itype]-1 + RunParams::firstbin[ngr + nsq + nfq + nfg + iexpt];
				
#else
				phist=*(hist.finder+offset+partialind)+edata.ek_rmin[iexpt*ntypes+itype]-1 + RunParams::firstbin[ngr + nsq + nfq + nfg + iexpt];
				
#endif
#endif				
	
				for (ir=0;ir<edata.ek_rused[iexpt*ntypes+itype];ir++)
				{
#ifndef _NO_PERIODIC
					*pekvalue+=(*phist++) * *pcoeff++;//summing 
#else
#ifdef _VIBR_AMP
					*pekvalue+=(*phist++) * *pcoeff++;//summing
#else
					*pekvalue+= (*pd_hist++)  * *pcoeff++;//summing    			
#endif
#endif
				}//end of data point cycle ir	
				if (fabs(*pekvalue) > TOLERANCE)
					all_zero = false;//there is at least one non-zero value
				*pekvalue *= myfactor;
				pekvalue++;//next E(k) data point
			}//end of k points cycle ik
		}//end of partial cycle ipartial
		min_ind++;
		max_ind++;
		if (all_zero)
		{
			thread_arg.ek_range_error[iexpt] = 1;//all_zero occured for this data sets k segment handled by this thread
			
		}
	}//end of data set cycle iexpt

};

//-----calculate the modified part of the partial pair correlation functions for the same data points, as the experimental
void CalcPart::CalcModPartial(ThreadArg &thread_arg)
{
	int iexpt,itype,ipartial,ir,ik,mybins;
	int offset;
	int myfactor;//to take into account that the histogram contains only the unique pairs, and for the pure partials should be multiplied by 2 to cnsider all atoms as central
	int *min_ind,*max_ind;//cycle boundaries
	int abs_type;
	
	double *p_calcpart,*p_matrix,*p_ppcf,*psqvalue,*pfqvalue,*pfgvalue,*pekvalue;
#ifdef _R_SWITCH_GR_MIN_1
	double *pr;//pointer to r values
#endif

#ifdef _NO_PERIODIC
	offset = 0;
#ifndef _VIBR_AMP
	double *pd_hist;
	pd_hist=0;
#else
	longint *phist;//pointer to the partial histogram arrays
	phist = 0;
#endif
#else
	longint *phist;
	phist = 0;
#endif

#ifndef _ATLAS
	int ibin, iq, ig;
	ibin=0;//to avoid warning
	iq=0;
	ig = 0;
#endif	
	
	for (ipartial = 0; ipartial < npartials; ipartial++)
	{
		if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
					//in case the main thread already generated the new move
		{
			//here the Threads parameters for no atomic move have to be used
			if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift))//only calculate for, if its is a normal atomic move
			{

				//g(r) fitting: Calculating the new  CalcPart modified partials
				//Each thread handle a different segment of the r points, but can use all the bins
				min_ind = thread_arg.min_r_index;
				max_ind = thread_arg.max_r_index;
				for (iexpt = 0; iexpt < ngr; iexpt++)//for each data set
				{
#ifdef _VIBR_AMP
					mybins = nbins_used;
					offset = 0;
#else
					mybins = nbins[RunParams::assign_hist[iexpt]];//set the binsize for the histogram of this data set, not necessarily the binsize of this data set
					offset = RunParams::assign_hist[iexpt] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppsf
#endif
					mybins -= RunParams::firstbin[iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
					//positions the pointer to the beginning of the iext-th experimet,
					//ipartial-th partial in the CalcPart grvalues, *min_ind data point
					p_calcpart = *(grfinder + iexpt * npartials + ipartial) + *min_ind;
					//positions the pointer to the beginning for this thread of the iexpt-th bin->dr 
					//conversion matrix of the DataMat object, it was only calculated for the bins used by this set
					p_matrix = datmat.grfinder[iexpt] + *min_ind * mybins;
#ifdef _R_SWITCH_GR_MIN_1
					pr = edata.gr_rvalues + *min_ind;//initialisation of pointer to r values
#endif
#ifdef _ATLAS//ATLAS is used, it will make the matrix * vector muliplication

					//sets the pointer to the beginning of the given partial 
					//in the PPCFSet grvalues array
					p_ppcf = *(ppcf.finder + offset + ipartial) + RunParams::firstbin[iexpt];

					//BLAS routine call: g(r)=bin->dr_converter*ppcf
					cblas_dgemv(CblasRowMajor, CblasNoTrans, (*max_ind - *min_ind + 1), mybins, 1.0, p_matrix, mybins, p_ppcf, 1, 0, p_calcpart, 1);

#else//ATLAS is not used
					for (ir = *min_ind; ir <= *max_ind; ir++)//for each r value
					{
						//positions the pointer to the beginning of the PPCFSet's g(r) values,
						//for the ipartial-th partial
						p_ppcf = *(ppcf.finder + offset + ipartial) + RunParams::firstbin[iexpt];
						*p_calcpart = 0;//resets the value to zero before summing

						for (ibin = 0; ibin < mybins; ibin++)//for each bin
						{
							//cout<<"ir "<<ir<<" ibin "<<ibin<<" p matrix "<<p_matrix<<endl;
							*p_calcpart += (*p_matrix++) * (*p_ppcf++);//convert bins to r data point
						}//end of bin cycle ibin
#ifdef _R_SWITCH_GR_MIN_1
						*p_calcpart -= 1.0;//to have (r^r_switch_power)*[g(r)-1]
						*p_calcpart *= pow(*pr, r_switch_power);
						pr++;
#endif
						p_calcpart++;//goes to the next r data point in CalcPart g(r) array
					}//end of r data point cycle ir
#endif
					min_ind++;
					max_ind++;
				}//end of data set cycle iexpt 

				//Each thread handle a different segment of the Q points, but uses all the r-points
				//min_Q_g_index, max_Q_g_index contains the segmentation for all Q-based data sets!
				min_ind = thread_arg.min_Q_g_index;
				max_ind = thread_arg.max_Q_g_index;

				//S(Q) fitting
				for (iexpt = 0; iexpt < nsq; iexpt++)//for each data set
				{
#ifdef _VIBR_AMP
					mybins = nbins_used;
					offset = 0;
#else
					mybins = nbins[RunParams::assign_hist[iexpt + ngr]];//set the binsize for this data set
					offset = RunParams::assign_hist[iexpt + ngr] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppsf
#endif
					mybins -= RunParams::firstbin[ngr + iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
					psqvalue = *(sqfinder + iexpt * npartials + ipartial) + *min_ind;//beginning of Q-values for this thread	for this partial and data set
					p_matrix = *(datmat.sqfinder + iexpt) + *min_ind * mybins;//sets the pointer to the first Q value for this thread of sqfinder of the iexpts Fourier matrix
#ifdef _ATLAS//ATLAS is used, it will make the matrix * vector muliplication

				//sets the pointer to the beginning of the given partial 
				//in the PPCFSet gr_0 array
					p_ppcf = PPCFSet::gr_0 + PPCFSet::offset[offset + ipartial] + RunParams::firstbin[ngr + iexpt];
					//BLAS routine call: A(Q)=Fouriermatrix*(ppcf-1)
					cblas_dgemv(CblasRowMajor, CblasNoTrans, (*max_ind - *min_ind + 1), mybins, 1.0, p_matrix, mybins, p_ppcf, 1, 0, psqvalue, 1);

#else//ATLAS is not used
					for (iq = *min_ind; iq <= *max_ind; iq++)
					{
						//sets the pointer to the beginning of the given partial 
						//in the PPCFSet grvalues array
						p_ppcf = *(ppcf.finder + offset + ipartial) + RunParams::firstbin[ngr + iexpt];

						*psqvalue = 0;//resets to zero before summing
						for (ir = 0; ir < mybins; ir++)
						{
							*psqvalue += (*p_ppcf - 1) * *p_matrix++;//summing 
							p_ppcf++;

						}//end of data point cycle ir	
						psqvalue++;//next S(Q) data point
					}//end of Q points cycle iq
#endif
					min_ind++;
					max_ind++;
				}//end of data set cycle iexpt 

			

				//reset it, as it might not have been set
				min_ind = thread_arg.min_Q_g_index+nsq;
				max_ind = thread_arg.max_Q_g_index+nsq;
				//F(Q) fitting
				for (iexpt = 0; iexpt < nfq; iexpt++)//for each data set
				{
				
#ifdef _VIBR_AMP
					mybins = nbins_used;
					offset = 0;
#else
					mybins = nbins[RunParams::assign_hist[iexpt + ngr + nsq]];//set the binsize for this data set
					offset = RunParams::assign_hist[iexpt + ngr + nsq] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppsf
#endif
					mybins -= RunParams::firstbin[ngr + nsq + iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
					pfqvalue = *(fqfinder + iexpt * npartials + ipartial) + *min_ind;//beginning of Q-values for this thread	for this partial and data set
					p_matrix = *(datmat.fqfinder + iexpt) + *min_ind * mybins;//sets the pointer to the first Q value for this thread of fqfinder of the iexpts Fourier matrix
#ifdef _ATLAS//ATLAS is used, it will make the matrix * vector muliplication

					//sets the pointer to the beginning of the given partial 
					//in the PPCFSet gr_0 array
					p_ppcf = PPCFSet::gr_0 + PPCFSet::offset[offset + ipartial] + RunParams::firstbin[ngr + nsq + iexpt];
					//BLAS routine call: A(Q)=Fouriermatrix*(ppcf-1)
					cblas_dgemv(CblasRowMajor, CblasNoTrans, (*max_ind - *min_ind + 1), mybins, 1.0, p_matrix, mybins, p_ppcf, 1, 0, pfqvalue, 1);

#else//ATLAS is not used
					for (iq = *min_ind; iq <= *max_ind; iq++)
					{
						//sets the pointer to the beginning of the given partial 
						//in the PPCFSet grvalues array
						p_ppcf = *(ppcf.finder + offset + ipartial) + RunParams::firstbin[ngr + nsq + iexpt];

						*pfqvalue = 0;//resets to zero befora summing
						for (ir = 0; ir < mybins; ir++)
						{
							*pfqvalue += (*p_ppcf - 1) * *p_matrix++;//summing 
							p_ppcf++;

						}//end of data point cycle ir	
						pfqvalue++;//next F(Q) data point
					}//end of Q points cycle iq
#endif			
				
					min_ind++;
					max_ind++;
				}//end of data set cycle iexpt 

				//F(g) fitting
				for (iexpt = 0; iexpt < nfg; iexpt++)//for each data set
				{
#ifdef _VIBR_AMP
					mybins = nbins_used;
					offset = 0;
#else
					mybins = nbins[RunParams::assign_hist[iexpt + ngr + nsq + nfq]];//set the binsize for this data set
					offset = RunParams::assign_hist[iexpt + ngr + nsq + nfq] * npartials;//for the ppcf finder to point to the firts partial of this data set's ppsf
#endif
					mybins -= RunParams::firstbin[ngr + nsq + nfq + iexpt];//start with the first bin of this set, as the histogram is calculated from the beginning regardless binshift
					pfgvalue = *(fgfinder + iexpt * npartials + ipartial) + *min_ind;//beginning of g-values for this thread	for this partial and data set
					p_matrix = *(datmat.fgfinder + iexpt) + *min_ind * mybins;//sets the pointer to the first g value for this thread of fgfinder of the iexpts Fourier matrix
#ifdef _ATLAS//ATLAS is used, it will make the matrix * vector muliplication

		//sets the pointer to the beginning of the given partial 
		//in the PPCFSet gr_0 array
					p_ppcf = PPCFSet::gr_0 + PPCFSet::offset[offset + ipartial] + RunParams::firstbin[ngr + nsq + nfq + iexpt];
					//BLAS routine call: A(g)=Fouriermatrix*(ppcf-1)
					cblas_dgemv(CblasRowMajor, CblasNoTrans, (*max_ind - *min_ind + 1), mybins, 1.0, p_matrix, mybins, p_ppcf, 1, 0, pfgvalue, 1);

#else//ATLAS is not used
					for (ig = *min_ind; ig <= *max_ind; ig++)
					{
						//sets the pointer to the beginning of the given partial 
						//in the PPCFSet grvalues array
						p_ppcf = *(ppcf.finder + offset + ipartial) + RunParams::firstbin[ngr + nsq + nfq + iexpt];

						*pfgvalue = 0;//resets to zero befora summing
						for (ir = 0; ir < mybins; ir++)
						{
							*pfgvalue += (*p_ppcf - 1) * *p_matrix++;//summing 
							p_ppcf++;

						}//end of data point cycle ir	
						pfgvalue++;//next F(g) data point
					}//end of g points cycle ig
#endif
					min_ind++;
					max_ind++;
				}//end of data set cycle iexpt 
			}//end of whether to calculate 
		}//end of if this partial was modified
	}//end of cycling through the partials ipartial
	
	//Each thread handle a different segment of the k points, but uses all the r-points
	//min_k_index, max_k_index contains the segmentation for all k-based data sets!
	min_ind=thread_arg.min_k_index;
	max_ind=thread_arg.max_k_index;
	//E(k) fitting, done separately, as not all the partials exist
	//sets the pointer to the beginning of the finder for the r-->k matrix
	for (iexpt=0; iexpt<nek; iexpt++)//for each data set
	{
		//here the Threads parameters for the no atomic moves have to be used!
		if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift) || (Threads::is_E0_shift && edata.ek_ngrid_in[iexpt] > 0))
		{//calculate if there is atomic move, or for those sets where there is change due to E0 shift 
			offset = npartials * RunParams::assign_hist[iexpt + ngr + nsq + nfq + nfg];//This gives the index to the first partial of this data set's histtogram of the hist finder
			abs_type = edata.ek_abstype[iexpt];
			//Determining the index of the partial histogram  
			for (itype = 0; itype < ntypes; itype++)//for each partial containing the absorbing particle type
			{
				//index of the partial		
				ipartial = (itype <= abs_type ? (itype * ntypes - (itype * (itype + 1) / 2) + abs_type) : \
					(abs_type * ntypes - (abs_type * (abs_type + 1) / 2) + itype));
				if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
						//in case the main thread already generated the new move
				{
					pekvalue = *(ekfinder + iexpt * ntypes + itype) + *min_ind;//beginning of k-values for this thread	for this partial and data set
					p_matrix = *(edata.ek_cfinder + edata.ek_cf_cum[iexpt] + Threads::EXAFS_gridind[iexpt] * ntypes + itype) + *min_ind * edata.ek_rused[iexpt * ntypes + itype];//sets the pointer to the first coefficient of 
					myfactor = 1;
					if (abs_type == itype)
						myfactor = 2;
					//the iexpts coefficient matrix for itype neighbour atom for this thread
					for (ik = *min_ind; ik <= *max_ind; ik++)
					{
						//sets the pointer to the first used bin of the given partial 
						//in the histogram
#ifdef _VIBR_AMP
						phist = *(hist.finder_T + ipartial) + edata.ek_rmin[iexpt * ntypes + itype] - 1 + RunParams::firstbin[ngr + nsq + nfq + nfg + iexpt];
#else
#ifdef _NO_PERIODIC
						pd_hist = *(hist.periodic_finder + ipartial) + edata.ek_rmin[iexpt * ntypes + itype] - 1 + RunParams::firstbin[ngr + nsq + nfq + nfg + iexpt];

#else
						phist = *(hist.finder + offset + ipartial) + edata.ek_rmin[iexpt * ntypes + itype] - 1 + RunParams::firstbin[ngr + nsq + nfq + nfg + iexpt];

#endif
#endif

						*pekvalue = 0;//resets to zero befora summing
						for (ir = 0; ir < edata.ek_rused[iexpt * ntypes + itype]; ir++)
						{
#ifndef _NO_PERIODIC							
							*pekvalue += (*phist++) * *p_matrix++;//summing 
#else
#ifdef _VIBR_AMP
							*pekvalue += (*phist++) * *p_matrix++;//summing 
#else
							*pekvalue += (*pd_hist++) * *p_matrix++;//summing  
#endif
#endif

						}//end of data point cycle ir	
						*pekvalue *= myfactor;
						pekvalue++;//next E(k) data point
					}//end of k points cycle ik
				}//end if this partial has been modified
			}//end of type cycle itype
		}//end of calculate for this set
		min_ind++;
		max_ind++;
	}//end of data set cycle iexpt 
}

//Copy the modified parts of the partials to target	
void CalcPart::CopyModified(CalcPart &target, ThreadArg &thread_arg)
{
	int ipartial,itype,iexpt,ir,iq,ig,ik;
	int *min_ind,*max_ind;//cycle boundaries
	int abs_type;
	double *p_calcpart,*p_calcpart_target;

	if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift))
	{
		for (ipartial = 0; ipartial < npartials; ipartial++)
		{
			if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
						//in case the main thread already generated the new move
			{
				//This partial has been modified
				//Each thread handle a different segment of the r points
				min_ind = thread_arg.min_r_index;
				max_ind = thread_arg.max_r_index;
				for (iexpt = 0; iexpt < ngr; iexpt++)//for each g(r) data set
				{
					//positions the pointer to the beginning of the iexpt-th experimet,
					//ipartial-th partial in the CalcPart grvalues, *min_ind data point
					p_calcpart = *(grfinder + iexpt * npartials + ipartial) + *min_ind;
					p_calcpart_target = *(target.grfinder + iexpt * npartials + ipartial) + *min_ind;
					for (ir = *min_ind; ir <= *max_ind; ir++)//for each r value
						*p_calcpart_target++ = *p_calcpart++;
					min_ind++;
					max_ind++;
				}//end of cycle iexpt

				//Each thread handle a different segment of the Q points
				min_ind = thread_arg.min_Q_g_index;
				max_ind = thread_arg.max_Q_g_index;
				for (iexpt = 0; iexpt < nsq; iexpt++)//for each S(Q) data set
				{
					//positions the pointer to the beginning of the iext-th experiment,
					//ipartial-th partial in the CalcPart sqvalues, *min_ind data point
					p_calcpart = *(sqfinder + iexpt * npartials + ipartial) + *min_ind;
					p_calcpart_target = *(target.sqfinder + iexpt * npartials + ipartial) + *min_ind;
					for (iq = *min_ind; iq <= *max_ind; iq++)//for each Q value
						*p_calcpart_target++ = *p_calcpart++;
					min_ind++;
					max_ind++;
				}//end of cycle iexpt

				for (iexpt = 0; iexpt < nfq; iexpt++)//for each F(Q) data set
				{
					//positions the pointer to the beginning of the iext-th experiment,
					//ipartial-th partial in the CalcPart fqvalues, *min_ind data point
					p_calcpart = *(fqfinder + iexpt * npartials + ipartial) + *min_ind;
					p_calcpart_target = *(target.fqfinder + iexpt * npartials + ipartial) + *min_ind;
					for (iq = *min_ind; iq <= *max_ind; iq++)//for each Q value
						*p_calcpart_target++ = *p_calcpart++;
					min_ind++;
					max_ind++;
				}//end of cycle iexpt
				for (iexpt = 0; iexpt < nfg; iexpt++)//for each F(g) data set
				{
					//positions the pointer to the beginning of the iext-th experiment,
					//ipartial-th partial in the CalcPart fgvalues, *min_ind data point
					p_calcpart = *(fgfinder + iexpt * npartials + ipartial) + *min_ind;
					p_calcpart_target = *(target.fgfinder + iexpt * npartials + ipartial) + *min_ind;
					for (ig = *min_ind; ig <= *max_ind; ig++)//for each g value
						*p_calcpart_target++ = *p_calcpart++;
					min_ind++;
					max_ind++;
				}//end of cycle iexpt

			}//end of if this partial has been modified
		}//end of cycle ipartial
	}

	//E(k) fitting, done separately, as not all the partials exist
	//Each thread handle a different segment of the k points
	min_ind=thread_arg.min_k_index;
	max_ind=thread_arg.max_k_index;
	for (iexpt=0; iexpt<nek; iexpt++)//for each E(k) data set
	{
		abs_type=edata.ek_abstype[iexpt];
		//Determining the index of the partial histogram  
		if (!(Threads::is_E0_shift || Threads::is_IQ_mucorr || Threads::is_fprime_shift) || (Threads::is_E0_shift && edata.ek_ngrid_in[iexpt] > 0))
		for (itype=0; itype<ntypes; itype++)//for each partial containing the absorbing particle type
		{
			//index of the partial		
			ipartial=(itype<=abs_type ? (itype*ntypes-(itype*(itype+1)/2)+ abs_type) : \
						(abs_type*ntypes-(abs_type*(abs_type+1)/2)+ itype));
			if (Threads::mod_partial[ipartial])//this has to be used instead  of the original move.modpart
					//in case the main thread already generated the new move
			{
				//positions the pointer to the beginning of the iext-th experiment,
				//itype-th partial in the CalcPart ekvalues, *min_ind data point
				p_calcpart=*(ekfinder+iexpt*ntypes+itype)+ *min_ind;																		
				p_calcpart_target=*(target.ekfinder+iexpt*ntypes+itype)+ *min_ind;																		
				for (ik= *min_ind; ik<= *max_ind; ik++)//for each k value
					*p_calcpart_target++=*p_calcpart++;
			
			}//end of if this partial has been modified
		}//end of cycle itype
		min_ind++;
		max_ind++;
	}//end of iexpt cycle
};
