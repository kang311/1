//source Aenet.cpp
//Last changed 12.12.2022

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#ifndef _DEF_THREADS
    #include "Threads.h"
    #define _DEF_THREADS
#endif

//xtern int aenet_Rc_max;
int Aenet::nAtoms=-1;
int Aenet::ntypes;
int Aenet::max_neigh;
int Aenet::ncells_nlist;//for the neighbour list, use 2*ncells_nlist+1 cells for neighbour search
int Aenet::overlap_nlist;//for the neighbousearch
int Aenet::max_tot_moved;
int Aenet::nchanged=0;//number of elemet in the *aenet_changed_neigh_ind array
int Aenet::nstored=0;//the number of stored atoms in the *aenet_stored_ind array
int Aenet::last_gen_aenet_accepted;//ngenerated for the last accepted ANN pot calculation step
int Aenet::last_gen_accepted;//ngenerated for the last accepted calculation step
int *Aenet::aenet_changed_neigh_ind=nullptr; 
int *Aenet::aenet_stored_ind;//if E is nor calculated in every step, store the indices of the moved atoms and their neighbours
double Aenet::cutoff=-1;//cutoff to pass to Aenet library
double Aenet::red_cutoff_sq=-1;//reduced cutoff squared
double *Aenet::E_moved;//temporary storage for the new enery of the moved atoms
bool Aenet::write_energy=false;//whether to save the atomic energies
bool Aenet::is_def_cutoff=1;//whether the original cutoff of the potential is used
int Aenet::calc_ANN=1;//indicator, that ANN should be calculated, used for initial as well, in loop 0: do not calculate in this step; 1: calculate in this step, E is calculated in every step, 2: calculate in this step, but E is not calculated in every step, 
char (*Aenet::aenet_type_names)[3];//name of the Aenet types, not necessarily the chem_symbols
//char (**Aenet::aenet_type_names_finder)[2];//finder to pass to aenet_init
char **Aenet::aenet_type_names_finder;//finder to pass to aenet_init
char(*Aenet::ANN_filenames)[FILE_NAME_SIZE];//ANN potential file names for each atom type
 

Aenet::Aenet(SimpleCfg &config, RunParams &rundata, NeighbourList &neigh) : conf(config),rundat(rundata),neighlist(neigh)//default constructor
{
    if (SimpleCfg::ntotal<=0)
    {
        cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"Aenet constructor"<<endl;
		cout<<"Number of atoms is zero!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
    }
    SetArraysize(&E_i,nAtoms,"Aenet::nAtoms","Aenet::Aenet");
    E_tot=0;
    if (::debug)
        cout<<"Aenet object is created"<<endl;
}
Aenet::Aenet(Aenet &source, SimpleCfg &config, RunParams &rundata, NeighbourList &neigh) : conf(config),rundat(rundata),neighlist(neigh)//default constructor
{
    if (SimpleCfg::ntotal<=0)
    {
        cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"Aenet constructor"<<endl;
		cout<<"Number of atoms is zero!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
    }
    SetArraysize(&E_i,nAtoms,"Aenet::nAtoms","Aenet::Aenet");
    for (int i=0;i<nAtoms;i++)
        E_i[i]=source.E_i[i];
    E_tot=source.E_tot;
    if (::debug)
        cout<<"Aenet object is created"<<endl;
}

void Aenet::Copy(Aenet &target, int *moved_indices, ThreadArg &thread_arg)//copy the energy
{
    int i;
    if (thread_arg.thread_index==RunParams::nthreads-1)
    {
        //this is done by the main thread it is safe to use moved indices here
        target.E_tot=E_tot;
        if (!Move::nomove)//only if there was atomic move
        {
            for (i=0;i<Move::tot_moved_atoms;i++)
                target.E_i[moved_indices[i]]=E_i[moved_indices[i]];
        }
    }
    for (i=0;i<thread_arg.aenet_at_count;i++)//only copy the changed ones
        target.E_i[thread_arg.aenet_at_index[i]]= E_i[thread_arg.aenet_at_index[i]];
    
}
void Aenet::SetParams()
{
    nAtoms =SimpleCfg::ntotal;
    ntypes=SimpleCfg::ntypes;
    ANN_filenames=new char[ntypes][FILE_NAME_SIZE];
    if (ANN_filenames==NULL)
        NoArray("Aenet::SetParams","ANN_filenames");
    for (int i =0;i<ntypes;i++)
        ANN_filenames[i][0]='\0';
    
    aenet_type_names=new char[ntypes][3];
    if (aenet_type_names==NULL)
        NoArray("Aenet::SetParams","aenet_type_names");
    //aenet_type_names_finder=new char[ntypes][2];
    SetArraysize(&aenet_type_names_finder, ntypes, "aenet_type_name_finder", "Aenet::SetParams");;
    if (aenet_type_names_finder==NULL)
        NoArray("Aenet::SetParams","aenet_type_names_finder");
    for (int i=0;i<ntypes;i++)
        aenet_type_names_finder[i]=(char*)aenet_type_names+i*3;
   
}

//initialize aenet, load the potentials
void Aenet::AenetInit()
{
    int stat;
    write_energy=RunParams::write_E;
   	aenet_init(ntypes,Aenet::aenet_type_names_finder,&stat);
	if (stat!=0)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"libaenet could not be initialized!"<<endl;
		cout<<"Exiting..."<<endl;
        CleanExit();
	}
    for (int i=0;i<ntypes;i++)
    {
        aenet_load_potential(i+1,(char*)(ANN_filenames+i),&stat);//fortran start array indexing with 1!
        if (stat!=0)
        {
            cout << "\n*****ERROR*****" << endl;
            cout<<"ANN potential "<<ANN_filenames[i]<<" could not be loaded!"<<endl;
            cout<<"Exiting..."<<endl;
            aenet_final(&stat);
            CleanExit();
        }
    }
};


//create the array for the neigh list
void Aenet::SetNeighList()
{
    max_tot_moved=RunParams::nmoved;
    if (RunParams::swap_fraction>0)
        max_tot_moved++;
     //calculate the average number of neighbours inside the cutoff, multiplying with 1.1 to accont for fluctuantions
    max_neigh=(int)(nAtoms/8.*4./3.*pow(cutoff/SimpleCfg::boxedge,3)*PI*1.1);
	ncells_nlist = int(cutoff / SimpleCfg::boxedge / NeighbourList::cellwidth) + 1;//number of grid cells to test for each direction for the neighbourlist
	overlap_nlist = (ncells_nlist * 2 + 1 - NeighbourList::ngridcells > 0 ? ncells_nlist * 2 + 1 - NeighbourList::ngridcells : 0);
    SetArraysize(&aenet_changed_neigh_ind, max_neigh*max_tot_moved,"aenet_changed_neigh_ind", "Aenet::SetNeighList");
    if (RunParams::aenet_step>1)//only needed, if E is not calcuilated in every move to store the atoms to recalculate old energy in the next E-calc step 
         SetArraysize(&aenet_stored_ind, nAtoms,"aenet_stored_ind", "Aenet::SetNeighList");
    SetArraysize(&E_moved, max_tot_moved,"E_moved", "Aenet::SetNeighList");

};



//calculate the potential
//this is parallelized for the threads
//as it is used both for the initial calculation and in the loop, not segments, but atom index lists are used 
void Aenet::CalcANN(ThreadArg &thread_arg)
{
    int k,iat;
    //int *neigh_list,*type_list;
    int c_type,nneigh;
    int stat;
    int count;
    //double *neigh_coord;
    double central_coord[3],*pE;
           
    pE=thread_arg.aenet_E_i;
    count=thread_arg.aenet_at_count;

    /*std::ofstream *plog;
    plog=&logfile;
    if (RunParams::nthreads>1)   
    {
        switch (thread_arg.thread_index)
        {
        case 0:
            plog=&logfile0;
            break;
        case 1:
            plog=&logfile1;
            break;
        case 2:
            plog=&logfile2;
            break;
        }
    }
    (*plog).precision(16);*/
    
    //*plog<<"b aenet count "<<count<<" first "<<thread_arg.aenet_at_index[0]<< " last "<<thread_arg.aenet_at_index[count-1]<<endl;
    for (k=0;k<count;k++)
    {
        iat=thread_arg.aenet_at_index[k];
        c_type=0;
		while (iat>=SimpleCfg::cumul[c_type +1])
		    c_type++;
        c_type++;//to have the type in the 1->ntypes range for the fortran
        for (int j=0;j<3;j++)
        {
            central_coord[j]=(*(conf.positions+iat*3+j)+1)*SimpleCfg::boxedge;//fill the central atom's coordinates in A
           /* central_coord[j]*=1e6;
            central_coord[j]=(longint)central_coord[j];
            central_coord[j]/=1e6;        */
        }
        
        nneigh=neighlist.CalcNeighList(iat,thread_arg.aenet_neighlist,thread_arg.aenet_neightype,nullptr,nullptr,nullptr,&thread_arg.aenet_max_nneigh,thread_arg.thread_index,thread_arg.aenet_neighcoord);
    
        /*if (iat==46)
        {
            *plog<<"atom index "<<iat<<" nneigh "<<nneigh<<endl; 
            for (int i=0;i<nneigh;i++)
                *plog<<"\tneighind "<<thread_arg.aenet_neighlist[i]<<" neightype "<<thread_arg.aenet_neightype[i]<<" "<<thread_arg.aenet_neighcoord[i*3]<<" "<<thread_arg.aenet_neighcoord[i*3+1]<<" "<<thread_arg.aenet_neighcoord[i*3+2]<<endl;
        }*/

        aenet_atomic_energy(central_coord,c_type,nneigh,thread_arg.aenet_neighcoord,thread_arg.aenet_neightype,pE,&stat);
        switch (stat)
        {
            case 1://cannot happen
                cout << "\n*****ERROR*****" << endl;
                cout<<"Aenet was not initialized! Exiting..."<<endl;
                CleanExit();
                break;
            case 2:
                cout << "\n*****ERROR*****" << endl;
                cout<<"Aenet could not allocate the required arrays! Exiting..."<<endl;
                CleanExit();
            break;
        }

        /*if (iat==445)
        *plog<<"CalcANN iat "<<iat<<" E_i["<<iat<<"] "<<*pE<<endl;*/
        pE++;
        
    }
    //*plog<<"a aenet count "<<count<<" first "<<thread_arg.aenet_at_index[0]<< " last "<<thread_arg.aenet_at_index[count-1]<<endl;
    
};


//set the cutoff and the reduced square
void Aenet::SetCutoff()
{
    if (RunParams::aenet_cutoff<0 )//not set
       cutoff=aenet_Rc_max;//set the cutoff of the potential
    else if  (RunParams::aenet_cutoff>aenet_Rc_max)
    {
        cout << "\nWARNING(" << ++warn << "): The ANN potential cutoff from the "<<datfilename <<" file is"<<endl;
        cout<<RunParams::aenet_cutoff<<" A, which is larger than the cutoff used during the potential training ("<<aenet_Rc_max<<" A)"<<endl;
        cout<<"The value used during the potential training wil be used!"<<endl;
        cutoff=aenet_Rc_max;
    }
    else
    {
        cutoff=RunParams::aenet_cutoff;
        is_def_cutoff=false;
    }
    red_cutoff_sq = pow(cutoff/SimpleCfg::boxedge, 2.0);//reduced, squared
};

void Aenet::SaveEnergy(const char *file_name)
{
    int i;
    ofstream file;

    //save to *.en file
	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE+10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE+10, energyfilename);
    OpenFile(file,tempfilename,"Aenet::Save",0);//open file, check, whether it was successfully opened
    file.precision(8);
    file.setf(ios::right, ios::adjustfield);
	file.setf(ios::fixed, ios::floatfield);
    file<<"This is an object of Class Aenet containing the atomic energies"<<endl;
	file<<J20<<nAtoms<<"\t number of atoms"<<endl;
    file<<J20<<E_tot<<"\t eV total energy\n"<<endl;
    
    file<<J10<<"index"<<J23<<"E (eV)"<<endl;
    for (i=0;i<nAtoms;i++)
        file<<J10<<i+1<<J23<<E_i[i]<<endl;

    file.unsetf(ios::fixed);
	file.unsetf(ios::right);

}
//saving the potential related parameters to the same *.pot file, which is used by FNC_POT for the other potential
void Aenet::SavePotBinary() const
{
    int i;
    ofstream file;
	OpenFile(file,potfilename,"Aenet::SavePotBinary",1);//open file, check, whether it was successfully opened

    file.write(reinterpret_cast<const char *>(&RunParams::potential), sizeof(RunParams::potential));//to be able to distiguish from other potential files
    i=1;//there is only one potential component presently, but just in case have the possibility for more
    file.write(reinterpret_cast<const char *>(&i), sizeof(i));
    file.write(reinterpret_cast<const char *>(&E_tot), sizeof(E_tot));
    file.write(reinterpret_cast<const char *>(&nAtoms), sizeof(nAtoms));
    file.write(reinterpret_cast<const char *>(E_i), nAtoms*sizeof(*E_i));
    file.write(reinterpret_cast<const char *>(&cutoff), sizeof(cutoff));

}

//loading *.pot
 int Aenet::LoadPotBinary()
 {
     int itemp;
     double temp;
     string warning = "\tHistogram and potentials has to be recalculated!";

	ifstream file;

	CleanOpen(file,potfilename,1);//open file
	if (CheckFileState(file,"Aenet::LoadPotBinary",potfilename)==0)//check, whether it was successfully opened
	{
		cout << "\nWARNING(" << ++warn << "): There is no open "<<potfilename<<" file to load the potentials from"<<endl;
		cout<<"\tPotential has to be recalculated!"<<endl;
		return(0);//loading was not successful
	}

     //reading the RunParams::potential
	if (ReadBin(file, &itemp, 1, "RunParams::potential", "Aenet::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (itemp!=RunParams::potential)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The potential switch is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<RunParams::potential<<" (ANN) based on the "<<datfilename<<" file!"<<endl;
		cout<<"\tPotential has to be recalculated!"<<endl;
		return(0);
	}
    if (ReadBin(file, &itemp, 1, "pot dimension", "Aenet::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (itemp!=1)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The ANN potential's dimension is "<<temp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and it should be 1 (only E_tot)"<<endl;
		cout<<"\tPotential has to be recalculated!"<<endl;
		return(0);
	}
    if (ReadBin(file, &E_tot, 1, "E_tot", "Aenet::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
    if (ReadBin(file, &itemp, 1, "nAtoms", "Aenet::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
	if (itemp!=nAtoms)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The number of atoms is "<<itemp<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<nAtoms<<" from the "<<cfgfilename<< " file!"<<endl;
		cout<<"\tPotential has to be recalculated!"<<endl;
		return(0);
	}
    if (ReadBin(file, E_i, nAtoms, "E_i", "Aenet::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
    if (ReadBin(file, &temp, 1, "cutoff", "Aenet::LoadPotBinary", potfilename) == 0) { cout << warning << endl; file.close(); return(0); }
    if (fabs(cutoff-Aenet::cutoff)>LOAD_TOL)
	{
		cout << "\nWARNING(" << ++warn << "): INCONSISTENCY! The ANN potential cutoff is "<<temp<<" A"<<endl;
		cout<<"\tfrom the "<<potfilename<<" file, and "<<Aenet::cutoff<<" A in this simulation!"<<endl;
		cout<<"\tPotential has to be recalculated!"<<endl;
		return(0);
	}
    return 1;
 };
 
 //this is not parallel, only it will use the main thread's parameters
 //first it is called with the old grid to determine the old neighbours which are stored in Aenet::aenet_changed_neigh_ind
 //then it is called for the new grid, to calculate the new energy for the moved atoms, and add the neighbours, which are not already in the 
 //aenet_changed_indices_array to it.
 //(Then the indices in the Aenet::aenet_changed_neigh_ind array are split among the threads
 //and the new E is calculated parallel by different routine)
 //mode=0: old grid, only to determine the neighbours, 
 //mode=1: new grid, calculate E for moved atoms and determine neighbors
 void Aenet::Update(Move &move, int mode, ThreadArg &thread_arg)
 {
    int imoved, i,j,iat,c_type,nneigh,stat,store;
    double central_coord[3],*pmoved_pos;
    
    logfile.precision(16); 

    if (mode)
        pmoved_pos=move.newpos;
    else
        pmoved_pos=move.oldpos;
    for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
    {
        iat=move.indices[imoved];
        c_type=0;
		while (iat>=SimpleCfg::cumul[c_type +1])
		    c_type++;
        c_type++;//to have the type in the 1->ntypes range for the fortran
        for (int j=0;j<3;j++)
        {
            central_coord[j]=(pmoved_pos[imoved*3+j]+1)*SimpleCfg::boxedge;//fill the central atom's coordinates in A
            /*central_coord[j]*=1e6;
            central_coord[j]=(longint)central_coord[j];
            central_coord[j]/=1e6;*/
        }

        //it will use the main thread's segment
        if (mode)
            nneigh=neighlist.CalcNeighList(iat,thread_arg.aenet_neighlist,thread_arg.aenet_neightype,nullptr,nullptr,nullptr,&thread_arg.aenet_max_nneigh,thread_arg.thread_index,thread_arg.aenet_neighcoord);
        else
            nneigh=neighlist.CalcNeighList(iat,thread_arg.aenet_neighlist,thread_arg.aenet_neightype,nullptr,nullptr,nullptr,&thread_arg.aenet_max_nneigh,thread_arg.thread_index,thread_arg.aenet_neighcoord,&move,imoved);
        //logfile<<"update mode "<<mode<<" atom index "<<iat<<" nneigh "<<nneigh<<endl; 
        
        for (i=0;i<nneigh;i++)
        {
            store=1;
            for (j=0;j<nchanged;j++)
            {
                if (aenet_changed_neigh_ind[j]==thread_arg.aenet_neighlist[i])
                {
                    store=0;
                    break;
                }
            }
            if (store==1)
            {
                if (nchanged>=max_tot_moved*max_neigh)
                    ResizeArray(&max_neigh,max_neigh+5,&aenet_changed_neigh_ind,"aenet_changed_neigh_ind","Aenet::Update",max_tot_moved);
                aenet_changed_neigh_ind[nchanged]=thread_arg.aenet_neighlist[i];
                nchanged++;
            }
        }
       
        if (mode && Aenet::calc_ANN)
        {
            //cout<<"nneigh "<<nneigh<<endl; 
            aenet_atomic_energy(central_coord,c_type,nneigh,thread_arg.aenet_neighcoord,thread_arg.aenet_neightype,E_moved+imoved,&stat);
            
            switch (stat)
            {
                case 1://cannot happen
                    cout << "\n*****ERROR*****" << endl;
                    cout<<"Aenet was not initialized! Exiting..."<<endl;
                    CleanExit();
                    break;
                case 2:
                    cout << "\n*****ERROR*****" << endl;
                    cout<<"Aenet could not allocate the required arrays! Exiting..."<<endl;
                    CleanExit();
                break;
            }
        }
        if (mode && Aenet::calc_ANN==0)//no E calculation in this step
        {
            store=1;
            // add the moved atom to the list as well
            for (j=0;j<nchanged;j++)
            {
                if (aenet_changed_neigh_ind[j]==iat)
                {
                    store=0;
                    break;
                }
            }
            if (store==1)
            {
                if (nchanged>=max_tot_moved*max_neigh)
                    ResizeArray(&max_neigh,max_neigh+5,&aenet_changed_neigh_ind,"aenet_changed_neigh_ind","Aenet::Update",max_tot_moved);
                aenet_changed_neigh_ind[nchanged]=iat;
                nchanged++;
            }
        }
       /*   if (mode)
            {
                logfile<<"\nstored indices for moved "<<move.indices[imoved]<<" nchanged "<<nchanged<<endl;
                for (j=0;j<nchanged;j++)
                    logfile<<"\t"<<aenet_changed_neigh_ind[j]<<endl;
            }
        //}*/  //GO
       
    }
 };

//update the list of stored indices to calculate the old energy for in the next ANN pot calc step, if not not calculated in every step 
void Aenet::UpdateStored()
{
    int i,j,store;
    for (i=0;i<nchanged;i++)
    {
        store=1;
        for (j=0;j<nstored;j++)
        {
            if (aenet_changed_neigh_ind[i]==aenet_stored_ind[j])
            {
                store=0;
                break;
            }
        }
        if (store==1)
        {
            aenet_stored_ind[nstored]=aenet_changed_neigh_ind[i];
            nstored++;
        }
    }
};

//recalculate the pot for the first atom of each atom type and check against the loaded 
int Aenet::CheckLoadedPot(ThreadArg &thread_arg)
{
    double norm;
	//the allocation was already made for the main thread, those atoms will be checked
    CalcANN(thread_arg);//the calculated energies will be in the thread_arg.E_i
    //check the energies
    for (int i=0;i<thread_arg.aenet_at_count;i++)
    {
    	if (fabs(thread_arg.aenet_E_i[i]-0)>TOLERANCE)
			norm=thread_arg.aenet_E_i[i];
		else
			norm=1.0;
        if (fabs((thread_arg.aenet_E_i[i]-E_i[thread_arg.aenet_at_index[i]])/norm)>LOAD_TOL)
		{
			cout.precision(15);
            cout.setf(ios::right, ios::adjustfield);
            cout.setf(ios::scientific, ios::floatfield);
            cout << "\nWARNING(" << ++warn << "): The loaded ("<<E_i[thread_arg.aenet_at_index[i]]<<") and recalculated ("<<thread_arg.aenet_E_i[i]<<") ANN potential is outside the tolerable margin of errors"<<endl;
           	cout << "\tmost probably due to a mismatch between the coordination and *.pot file!"<<endl;
			cout << "\tPotential will be recalculated!"<<endl;
            cout.precision(6);
            cout.unsetf(ios::scientific);
            cout.unsetf(ios::right);
			return(0);
		}
    }
    return 1;
};

