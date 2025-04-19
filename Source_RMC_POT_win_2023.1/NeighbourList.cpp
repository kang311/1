//source NeighbouList.cpp
//Last changed 24.01.2023

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include "Threads.h"//Move.h"//NeighbourList.h is included through this

int    NeighbourList::ngridcells;//number of gird cells in one dimension
int	   NeighbourList::ngridcells2;//square of number of gird cells
int	   NeighbourList::max_gridatom;//maximum number of atoms in a grid cell
int    NeighbourList::natoms;//total number of atoms
int    NeighbourList::ntypes;// number of atom types
int	   NeighbourList::nconsttype;//number of constrained types
int    NeighbourList::ncells_cutoff;//number of grid cells to test for each direction for cutoff check
int	   NeighbourList::overlap;//number of gridcells to leave out from checking of the cutoff in each positiv direction due to overlapping
int	   NeighbourList::overlap_nlist;//number of gridcells to leave out from calculation of the neighbourlist in each positiv direction due to overlapping
int    NeighbourList::ncells_nlist;//number of grid cells to test for each direction for neighbourlist
int    NeighbourList::max_neigh;//maximum number of neighbopurs for the neighbourlist for an atom
int	   NeighbourList::copy_changed;//indicator, whether in case of resizing the arrays connected to max_neigh the changed_* arrays should be copied
int	   NeighbourList::usedatom;//number of atoms involved in neighbourlist calculation
longint NeighbourList::consttype=0;//each bit represents an atom type, containing a 1 for each type that is involved in a constraint, which require the neighbour list
double NeighbourList::xmax;//the maximum distance for the neighbourlist in reduced unit
double NeighbourList::xmax2;//square of the maximum distance for the neighbourlist in reduced unit
double NeighbourList::max_red_cutoff;//maximum cutoff distance in reduced unit
double NeighbourList::cellwidth;//the width of the grid cell in reduced units
double NeighbourList::half_ngridcells;//the half of ngridcells, to increase speed
#ifdef _AENET
	Threads *NeighbourList::thread=nullptr;
#endif

//------------constructor---------------
NeighbourList::NeighbourList(SimpleCfg &conf, FNC_POT* &fncneigh, RunParams &rundata)
				:config(conf), fnc(fncneigh), rundat(rundata)
{
	int i,j,offset;
	int tot_moved=RunParams::nmoved;

	if (ngridcells==0)
	{	cout<<"\n"<<"*****ERROR*****"<<endl;
		cout<<"NeighbourList constructor: Number of grid cells is 0!"<<endl;
		cout<<"Cannot run this way, exiting..."<<endl;
		CleanExit();
	};

	if (RunParams::swap_fraction>0)
		tot_moved++;

	SetArraysize(&grid_index,3*natoms,"grid_index","NeighbourList::NeighbourList");//stored in order atom_0:x_cell,y_cell,z_cell;atom_1:x_cell,y_cell,z_cell...

	SetArraysize(&grid_count,ngridcells*ngridcells*ngridcells,"grid_count","NeighbourList::NeighbourList");//count stored in the following order, cell indices denoted as x_i,y_j,z_k for the i-j-k cell
					//x0-y0-z0,x1-y0-z0,x2-y0-z0....x_last-y0-z0,x0-y1-z0,x1-y1-z0....x_last-y_last-z0,x0-y0-z1....
	SetArraysize(&atom_indices,ngridcells*ngridcells*ngridcells*max_gridatom,"atom_indices","NeighbourList::NeighbourList");//cells stored in the same order, as above, for each cell max_gridatom element
				//x0-y0-z0:atom(0)->atom(max_gridatom-1),x1-y0-z0:atom(0)->atom(max_gridatom-1),x2-y0-z0:atom(0)->atom(max_gridatom-1)...
	SetArraysize(&old_grid,3*tot_moved,"old_grid","NeighbourList::NeighbourList");//indices of the grid cell for each moved atom before the move

	//calculating the necessary dimension of the neighbourlist and the squared distance list array
	if (RunParams::ncosdistr ||RunParams::ncommonneigh||RunParams::nsecondneigh)
	{
		usedatom = -1;
		if (nconsttype==ntypes)//all the types are involved, neighbourlist is needed for all the atoms
			usedatom=natoms*max_neigh;

		if (usedatom==-1)
		{
			//determining, how many different types are involved in the cosine distribution of bond angle, common or second  neighbour constraint, as
			//the neighbourlist has to be kept for each of their atoms
			//all the neighbours will be included into the list, regardless their type!
			
			j=0;
			usedatom = 0;
			for (j = 0; j < SimpleCfg::ntypes; j++)
			{
				if (consttype>>j & 1)//this j-th type is involved
						usedatom += SimpleCfg::pnatoms[j];//number of atoms of this type				
			}
		}
		
	}


	SetArraysize(&nneigh,usedatom,"nneigh","NeighbourList::NeighbourList");
	//there will be maxneigh element for each atom
	SetArraysize(&neigh_list,usedatom*max_neigh,"neigh_list","NeighbourList::NeighbourList");
	SetArraysize(&neigh_type,usedatom*max_neigh,"neigh_type","NeighbourList::NeighbourList");
	SetArraysize(&dsq,usedatom*max_neigh,"dsq","NeighbourList::NeighbourList");
	SetArraysize(&vector_comp,usedatom*max_neigh*3,"vector_comp","NeighbourList::NeighbourList");

	//creating the finder for the number of neighbours, the neighbourlist, the squared distance list, and the 
	//Central->Neighbour vector component list. There will be a finder for each type, even if this type is not
	//represented on the neighbourlist, but in this case the finder value will be set to NULL
	SetArraysize(&nneigh_finder,ntypes,"nneigh_finder","NeighbourList::NeighbourList");
	SetArraysize(&neighlist_finder,ntypes,"neighlist_finder","NeighbourList::NeighbourList");
	SetArraysize(&neightype_finder,ntypes,"neightype_finder","NeighbourList::NeighbourList");
	SetArraysize(&dsq_finder,ntypes,"dsq_finder","NeighbourList::NeighbourList");
	SetArraysize(&vector_finder,ntypes,"vector_finder","NeighbourList::NeighbourList");

	//arrays for storing the old values of the neigh lists, while it is decided, that the move is acceptable or not
	SetArraysize(&changed_indices,tot_moved*(1+2*max_neigh),"changed_indices","NeighbourList::NeighbourList");//array of atoms (moved, and their new and old neighbours), which are effected
						//by the move
	SetArraysize(&changed_types,tot_moved*(1+2*max_neigh),"changed_types","NeighbourList::NeighbourList");//array of types of the atoms (moved, and their new and old neighbours), which are effected
						//by the move
	//for each central atom represented by the changed_indices list
	SetArraysize(&changed_nneigh,tot_moved*(1+2*max_neigh),"changed_nneigh","NeighbourList::NeighbourList");//number of neighbours 
	SetArraysize(&changed_neigh_list,tot_moved*(1+2*max_neigh)*max_neigh,"changed_neigh_list","NeighbourList::NeighbourList");//indices of the neighbours
	SetArraysize(&changed_neigh_type,tot_moved*(1+2*max_neigh)*max_neigh,"changed_neigh_type","NeighbourList::NeighbourList");//types of the neighbours
	SetArraysize(&changed_dsq,tot_moved*(1+2*max_neigh)*max_neigh,"changed_dsq","NeighbourList::NeighbourList");//squared distance for the neighbours
	SetArraysize(&changed_vector_comp,tot_moved*(1+2*max_neigh)*max_neigh*3,"changed_vector_comp","NeighbourList::NeighbourList");//components of the Central-Neighbour vector

	offset=0;
	for (i=0;i<ntypes;i++)
	{
		if (consttype>>i & 1)//this type is represented
		{
			nneigh_finder[i]=nneigh+offset;
			neighlist_finder[i]=neigh_list+max_neigh*offset;
			neightype_finder[i]=neigh_type+max_neigh*offset;
			dsq_finder[i]=dsq+max_neigh*offset;
			vector_finder[i]=vector_comp+max_neigh*offset*3;
			offset+=SimpleCfg::pnatoms[i];//increase offset with the number of atoms of this type * max_neigh
		}
		else
		{
			//this type is not represented, the finder will not be called for it normally
			nneigh_finder[i]=NULL;
			neighlist_finder[i]=NULL;
			neightype_finder[i]=NULL;
			dsq_finder[i]=NULL;
			vector_finder[i]=NULL;
		}



	}
	for (i=0;i<ngridcells*ngridcells*ngridcells;i++)
		grid_count[i]=0;//setting to zero
};

//initializing the static members
void NeighbourList::SetParams()
{
	int temp;
	int safe_add = SAFE_ADD;
	double gridcount;

	natoms = SimpleCfg::ntotal;
	ntypes = SimpleCfg::ntypes;
	max_gridatom = (natoms < RunParams::max_gridatom ? natoms : RunParams::max_gridatom);//maximum number of atoms in a grid cell

	//calculating the number of grid cells
	gridcount = natoms / max_gridatom;
	ngridcells = int(pow(gridcount, 1. / 3.));//determining the nearest cubic number
	if (gridcount - pow((double)ngridcells, 3.) > pow((double)ngridcells + 1, 3.) - gridcount)
		ngridcells++;

	temp = (int)(natoms / pow((double)ngridcells, 3.)) + 1;//resetting the maximum number of atoms in a gridcell
	//for safety reason increase the number of max_gridatoms with 20 percent (or at least SAFE_ADD), whichever is more
	max_gridatom = (int)(temp + (temp * 0.2 > safe_add ? temp * 0.2 : safe_add));
	if (::debug)
	{
		cout << "\nNOTE(" << ++note << "): For safety reasons " << max_gridatom << " atom/grid will be used in array dimensions!" << endl;
		cout << "\t" << ngridcells << " will be the number of grid cells in each dimension." << endl;
	}

	//setting the additional parameters
	ngridcells2 = ngridcells * ngridcells;//to make calculation quicker
	half_ngridcells = ngridcells / 2.0;
	cellwidth = 2. / ngridcells;//the width of a grid cell

	//PRESENTLY THE NEIGHBOUR LIST IS APPLIED FOR CosDistrConst, CommonNeighConst, SecondNeighConst, CoordNumbConst detail and AENET
	//IF IT WILL BE APPLIED FOR OTHER TASKS, THIS HAS TO BE INCLUDED HERE FOR THE CALCULATION OF THE LARGEST DISTANCE
	//setting the maximum number of neighbours for the neighbourlist for an atom
	//as the neighbourlist can be used with different maximum calculation distances, the largest distance 
	//applied as the maximum for the calculation of the CosDistrConst, CommonNeighConst or SecondNeighConst will be used for calculating the array dimension of neighlist
	//calculating the average number of neighbours inside the sphere 
	//the calcultion uses the reduced number density, so all the neighbours regardless the type is included in it 
	xmax2 = CosDistrConst::maxdistsq;

#ifdef _ADVANCED_GEOM_CONST
	xmax2 = (CommonNeighConst::maxdistsq > xmax2 ? CommonNeighConst::maxdistsq : xmax2);//reduced squared maximum distance for the neighbourlist 
	xmax2 = (SecondNeighConst::maxdistsq > xmax2 ? SecondNeighConst::maxdistsq : xmax2);//reduced squared maximum distance for the neighbourlist 
#endif
	
//#ifdef _AENET//GO
//	xmax2 = (Aenet::red_cutoff_sq > xmax2 ? Aenet::red_cutoff_sq : xmax2);//reduced squared maximum distance for the neighbourlist 
//#endif
	xmax = sqrt(xmax2);//reduced maximum distance for the neighbourlist 

	temp = int(natoms / 8 * 4 / 3 * pow(xmax, 3.) * PI) + 1;
	//for safety reason increase the number of max_neigh with 10 percent (or at least SAFE_ADD), whichever is more
	max_neigh = (int)(temp + (temp * 0.1 > safe_add ? temp * 0.1 : safe_add));
	ncells_nlist = int(xmax / cellwidth) + 1;//number of grid cells to test for each direction for the neighbourlist
	overlap_nlist = (ncells_nlist * 2 + 1 - ngridcells > 0 ? ncells_nlist * 2 + 1 - ngridcells : 0);
}

void NeighbourList::SetCutoff()
{
	int i;
	int ntypes = SimpleCfg::ntypes;
	//calculating the maximum of the reduced cutoff distances for checking, whteher the move is acceptable on the 
	//basis of cutoffs
	for (i=0;i<ntypes*(ntypes+1)/2;i++)
	{
		if (RunParams::pcutoff[i]>max_red_cutoff)
			max_red_cutoff=RunParams::pcutoff[i];
	}
	max_red_cutoff/=RunParams::boxedge;
	ncells_cutoff=int(max_red_cutoff/cellwidth)+1;//number of grid cells to test for each direction to check cutoff
	//number of gridcells to leave out from checking in each positiv direction due to overlapping
	overlap=(ncells_cutoff*2+1 -ngridcells>0?ncells_cutoff*2+1 -ngridcells:0);

};


//--------setting or reding grid_count, the number of atoms in the (i,j,k) grid cell---------
int &NeighbourList::GridCount(int i,int j,int k)
{
	return(*(grid_count+k*ngridcells2+j*ngridcells+i));
};

//--------spointer to the number of atoms in the (i,j,k) grid cell---------
int *NeighbourList::PtrGridCount(int i,int j,int k)
{
	return(grid_count+k*ngridcells2+j*ngridcells+i);
};

//----the index of the atom, this is the atom_ind_grid-th atom for this grid
int &NeighbourList::AtomIndex(int i,int j,int k, int atom_ind_grid)
{
	return(*(atom_indices+(k*ngridcells2+j*ngridcells+i)*max_gridatom+atom_ind_grid));
};

//----pointer to the index of the atom, this is the atom_ind_grid-th atom for this grid
int *NeighbourList::PtrAtomIndex(int i,int j,int k, int atom_ind_grid)
{
	return(atom_indices+(k*ngridcells2+j*ngridcells+i)*max_gridatom+atom_ind_grid);
};
//------------Creating the grid, determining the gridcell for each atom----------------
void NeighbourList::MakeGrid()
{
	int i,ic, jc, kc;
	int old_max_gridatom;
	int atom_count;
	int *pgrid_ind;
	int safe_add=SAFE_ADD;
	double *pcoord;


	pcoord=config.positions;//position of the first atom
	pgrid_ind=grid_index;
	for (i=0;i<natoms;i++)
	{
		ic=int((*pcoord +1)*half_ngridcells);
		pcoord++;
		jc=int((*pcoord +1)*half_ngridcells);
		pcoord++;
		kc=int((*pcoord +1)*half_ngridcells);
		pcoord++;
		//if the coordinate component is +1, it is the same, as -1, should go to the first grid cell in this dimension (with index=0)
		if (ic>ngridcells-1)
			ic-=ngridcells;
		if (jc>ngridcells-1)
			jc-=ngridcells;
		if (kc>ngridcells-1)
			kc-=ngridcells;
		atom_count=GridCount(ic,jc,kc);//to save the time of calculating it several times

		if (atom_count==max_gridatom)
		{
			if (::debug)
			{
				cout << "\nWARNING(" << ++warn << "): NeighbourList::MakeGrid: Gridding problem, too many particles in a grid cell!" << endl;
				cout << "\tResize array and continue..." << endl;
			}
			old_max_gridatom=max_gridatom;
			max_gridatom=(int)(max_gridatom+ (max_gridatom*0.2> safe_add ? max_gridatom*0.2 : safe_add));//new value 
			ResizeIArray(old_max_gridatom,max_gridatom,(int)pow((double)ngridcells,3),(int)pow((double)ngridcells,3),&atom_indices,NULL);
		}

		//the index of the atom, this is the GridCount(ic,jc,kc)-th atom for this grid
		AtomIndex(ic,jc,kc,atom_count)=i;
		
		//indices of this grid cell for this atom 
		*pgrid_ind++=ic;
		*pgrid_ind++=jc;
		*pgrid_ind++=kc;
		//increasing the number of atoms belonginig to this grid

		GridCount(ic,jc,kc)++;
	}
};
//----------updating the gridding for the moved atom(s) ----------------
void NeighbourList::UpdateGrid(Move &move)
{
	int imoved;
	int i,ic, jc, kc;
	int old_max_gridatom;
	int safe_add=SAFE_ADD;
	int *atom_count,*old_count;
	int *pgrid_ind, *poldgrid_ind,*moved_index,*patom_ind;
	double *pcoord;

	pcoord=move.newpos;//position of the moved atoms
	moved_index=move.indices;//index of the moved atoms
	poldgrid_ind=old_grid;//to store the indices of the old grid cells for each moved atom in case the move is rejected
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		
		pgrid_ind=grid_index+3* (*moved_index);//position at the first gridcell coordinate of this atom
		*poldgrid_ind++=*pgrid_ind;//store the old grid cell indices
		*poldgrid_ind++=*(pgrid_ind+1);
		*poldgrid_ind++=*(pgrid_ind+2);

		ic=int((*pcoord +1)*half_ngridcells);
		pcoord++;
		jc=int((*pcoord +1)*half_ngridcells);
		pcoord++;
		kc=int((*pcoord +1)*half_ngridcells);
		pcoord++;
		//if the coordinate component is +1, it is the same, as -1, should go to the first grid cell in this dimension (with index=0)
		if (ic==ngridcells)
			ic=0;
		if (jc==ngridcells)
			jc=0;
		if (kc==ngridcells)
			kc=0;

		//to see, if the atom changed grid cell
		if (ic== *pgrid_ind && jc== *(pgrid_ind+1) && kc== *(pgrid_ind+2))
		{
			moved_index++;//remained in the same grid, go to next moved atom
			continue;
		}

		atom_count=PtrGridCount(ic,jc,kc);//to save the time of calculating it several times

		if (*atom_count==max_gridatom)
		{
			cout << "\nWARNING(" << ++warn << "): NeighbourList::UpdateGrid: Gridding problem, too many particles in a grid cell!"<<endl;
			cout<<"\tResize array and continue..."<<endl;
			old_max_gridatom=max_gridatom;
			max_gridatom=(int)(max_gridatom+ (max_gridatom*0.2> safe_add ? max_gridatom*0.2 : safe_add));//new value 
			ResizeIArray(old_max_gridatom,max_gridatom,(int)pow((double)ngridcells,3),(int)pow((double)ngridcells,3),&atom_indices,NULL);
		}

		//remove the atom's contribution from the old grid cell
		old_count=PtrGridCount(*pgrid_ind,*(pgrid_ind+1),*(pgrid_ind+2));//set the pointer to the count of the old grid cell
		(*old_count)--;//decrease the atom count for the old grid cell
		
		//remove the index of the atom from the grid cell belonging to the original position
		patom_ind=PtrAtomIndex(*pgrid_ind,*(pgrid_ind+1),*(pgrid_ind+2),0);//pointer to the index of first atom for the  grid
		for (i=0;i<*old_count;i++)
		{
			if (*patom_ind==*moved_index)
			{
				//change the index of the moved atom wih the last on the list
				*patom_ind=AtomIndex(*pgrid_ind,*(pgrid_ind+1),*(pgrid_ind+2),*old_count);
				break;
			}
			patom_ind++;
		}

		//set the new values
		AtomIndex(ic,jc,kc,*atom_count)=*moved_index;
		
		//indices of this grid cell for this atom 
		*pgrid_ind++=ic;
		*pgrid_ind++=jc;
		*pgrid_ind++=kc;
		//increasing the number of atoms belonginig to this grid
		(*atom_count)++;
		moved_index++;//go to next moved atom

	}

};

//----------updating the gridding for the moved atom(s) ----------------
void NeighbourList::ResetGrid(int *moved_indices)
{
	int imoved;
	int i;
	int old_max_gridatom;
	int safe_add=SAFE_ADD;
	int *atom_count,*old_count;
	int *pgrid_ind, *poldgrid_ind,*moved_index,*patom_ind;

	moved_index=moved_indices;//index of the moved atoms
	poldgrid_ind=old_grid;//to the indices of the grid cells before the move for each moved atom 
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		
		pgrid_ind=grid_index+3* (*moved_index);//position at the first gridcell coordinate of this atom

		//to see, if the atom changed grid cell
		if (*poldgrid_ind == *pgrid_ind && *(poldgrid_ind+1) == *(pgrid_ind+1) && *(poldgrid_ind+2)== *(pgrid_ind+2))
		{
			moved_index++;//remained in the same grid, go to next moved atom
			poldgrid_ind+=3;
			continue;
		}

		//pointer to the number of atoms in the grid cell the atom should go back to
		atom_count=PtrGridCount(*poldgrid_ind,*(poldgrid_ind+1),*(poldgrid_ind+2));//to save the time of calculating it several times

		if (*atom_count==max_gridatom)//this theoretically cannot happen, as this was checked, when the atom originally got into this cell
		{
			cout << "\nWARNING(" << ++warn << "): NeighbourList::ResetGrid: Gridding problem, too many particles in a grid cell!"<<endl;
			cout<<"\tResize array and continue..."<<endl;
			old_max_gridatom=max_gridatom;
			max_gridatom=(int)(max_gridatom+ (max_gridatom*0.2> safe_add ? max_gridatom*0.2 : safe_add));//new value 
			ResizeIArray(old_max_gridatom,max_gridatom,(int)pow((double)ngridcells,3),(int)pow((double)ngridcells,3),&atom_indices,NULL);
		}

		//remove the atom's contribution from the old grid cell
		old_count=PtrGridCount(*pgrid_ind,*(pgrid_ind+1),*(pgrid_ind+2));
		(*old_count)--;//decrease the atom count for the old grid cell
		
		//remove the index of the atom from the old grid cell
		patom_ind=PtrAtomIndex(*pgrid_ind,*(pgrid_ind+1),*(pgrid_ind+2),0);//index of first atom for the old grid
		for (i=0;i<*old_count;i++)
		{
			if (*patom_ind==*moved_index)
			{
				//change the index of the moved atom wih the last on the list
				*patom_ind=AtomIndex(*pgrid_ind,*(pgrid_ind+1),*(pgrid_ind+2),*old_count);
				break;
			}
			patom_ind++;
		}

		//set back the original values before the move
		AtomIndex(*poldgrid_ind,*(poldgrid_ind+1),*(poldgrid_ind+2),*atom_count)=*moved_index;
		
		//indices of this grid cell for this atom 
		*pgrid_ind++=*poldgrid_ind++;
		*pgrid_ind++=*poldgrid_ind++;
		*pgrid_ind=*poldgrid_ind++;
		//increasing the number of atoms belonginig to this grid
		(*atom_count)++;
		moved_index++;//go to next moved atom

	}

};
//------------------saving the grid to a file------------------
void NeighbourList::SaveGrid(const char *file_name)
{
	int i,ix,iy,iz,count;
	int *pint;
	ofstream file;
	if (strlen(file_name)!=0)
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
	else
		mystrcpy(tempfilename, FILE_NAME_SIZE + 10, gridfilename);
		
	OpenFile(file,tempfilename,"Grid::Save",0);//open file, check, whether it was successfully opened
		
	file<<"This is an object of class NeighbourList "<<endl;
	file<<natoms<<"\t number of atoms"<<endl;
	file<<ntypes<<"\t number of atom types"<<endl;
	file<<max_gridatom<<"\t maximum number of atoms in a grid cell"<<endl;
	file<<ngridcells<<"\t number of grid cells in each dimension"<<endl;

	//saving the indices of the gridcells for each atoms
	file<<"\nSaving the indices of the gridcells for each atoms"<<endl;
	pint=grid_index;
	for (i=0;i<natoms;i++)
	{
		file<<*pint++<<"\t";//grid index in x direction 
		file<<*pint++<<"\t";//grid index in y direction
		file<<*pint++<<endl;;//grid index in z direction
	}
		
	//Saving the number and the indices of the atoms belonging to each grid cells
	//cell order:x0-y0-z0:atom(0)->atom(cellcount-1),x1-y0-z0:atom(0)->atom(cellcount-1),x2-y0-z0:atom(0)->atom(cellcount-1)...
	file<<"\nSaving the number and the indices of the atoms belonging to each grid cells"<<endl;
	pint=atom_indices;
	for (iz=0;iz<ngridcells;iz++)//grid index in the direction of z-axis
	{
		for (iy=0;iy<ngridcells;iy++)//grid index in the direction of y-axis
		{
			for (ix=0;ix<ngridcells;ix++)//grid index in the direction of x-axis
			{
				count=GridCount(ix,iy,iz);
				file<<count;
				for (i=0;i<count;i++)//atom indices belonging to this grid
					file<<"\t"<<*pint++;
				pint+=max_gridatom-count;//setting to the beginning of the next grid cell
				file<<endl;
			}
		}

	}
	file.close();
}

#ifdef _NEI
	//------------------saving the neighbourlist to a file------------------
	void NeighbourList::SaveNeilist(const char *file_name )
	{
		int i,j,count,itype;
		int *pnneigh,*pind,*ptype;
		double *pdsq,*pvec;
		ofstream file;
		if (strlen(file_name)!=0)
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, file_name);
		else
			mystrcpy(tempfilename, FILE_NAME_SIZE + 10, neighfilename);
		
		OpenFile(file,tempfilename,"Grid::Save",0);//open file, check, whether it was successfully opened
		
		
		file<<"This is an object of class NeighbourList calculated to "<<xmax*RunParams::boxedge<<" A"<<endl;
		//saving the neighbourlist
		pnneigh=nneigh;//number of neighbours
		pind=neigh_list;//indices of the neighbours
		ptype=neigh_type;//type of the neighbours
		pdsq=dsq;
		pvec=vector_comp;
		file.precision(15);
		file.setf(ios::fixed, ios::floatfield);
		file.setf(ios::right, ios::adjustfield);
		for (itype=0;itype<ntypes;itype++)
		{
			if (consttype>>itype & 1)//this type is involved
			{
				file<<"\n"<<itype+1<<" type central index then neighbours"<<endl;
				for (i=SimpleCfg::cumul[itype];i<SimpleCfg::cumul[itype+1];i++)//going through the atoms of this type
				{
					count=*pnneigh++;
					file<<i;
					for (j=0;j<count;j++)
					{
						if ((consttype>>*(ptype+j) & 1))
							file<<"\t"<<*pind++; //this is a constrained type so for this the list is saved for it
						else
							pind++;//unconstrained type particles's index is not saved, as the list is not updated for them
					}
					pind+=max_neigh-count;//go to the beginning of the next atoms's first neighbour
					file<<endl;
#ifdef _NEIE
					//writing squared distances
					file<<i;
					for (j=0;j<count;j++)
					{
						if ((consttype>>*(ptype+j) & 1))
							file<<"\t"<<*pdsq++; //this is a constrained type so for this the list is updated
						else
							pdsq++;//unconstrained type particles's index is not saved, as the list is not updated for them
					}
					pdsq+=max_neigh-count;//go to the beginning of the next atoms's first first neighbours's squard distance
					file<<endl;

					//writing vector components
					file<<i;
					for (j=0;j<count;j++)
					{
						if ((consttype>>*(ptype+j) & 1))
						{
							//this is a constrained type so for this the list is updated
							file<<"\t"<<*pvec++;//x comp
							file<<"\t"<<*pvec++;//y comp
							file<<"\t"<<*pvec++;//z comp
						}
						else
							pvec+=3;//unconstrained type particles's vector comps are not saved, as the list is not updated for them
					}
					
					pvec+=3*(max_neigh-count);//go to the beginning of the next atoms's first neighbour's vector component
					file<<endl;
#endif				
				ptype+=max_neigh;//go to the beginning of the next atoms's first neighbours's type 	
				}
			}//end if this type is involved
		}//end of itype cycle
		file.precision(6);
		file.unsetf(ios::fixed);
		file.unsetf(ios::right);
		file.close();
	};
#endif

//------------loading from file---------------
int NeighbourList::LoadGrid()
{
	int dumb,i,ix,iy,iz,count;
	int *pint;
	char *conv_numb1,*conv_numb2,*conv_numb3,*conv_numb4;
	char  *conv_tx,*conv_ty,*conv_tz, *conv_t;
	ifstream file;

	SafeOpenTextFile(file,gridfilename);
	if (CheckFileState(file,"NeighbourList::LoadGrid",gridfilename)==0)
	{
		cout << "\tNeighbourList has to be calculated!"<<endl;
		return(0);//load was not successful
	};
	conv_numb1 = NULL;
	conv_numb2 = NULL;
	conv_numb3 = NULL;
	conv_numb4 = NULL;
	SetArraysize(&conv_t,100,"conv_t","NeighbourList::LoadGrid");
	SetArraysize(&conv_tx,100,"conv_tx","NeighbourList::LoadGrid");
	SetArraysize(&conv_ty,100,"conv_ty","NeighbourList::LoadGrid");
	SetArraysize(&conv_tz,100,"conv_tz","NeighbourList::LoadGrid");

	
	SkipLine(file,gridfilename,1);
	cout<<"\nLoading the NeighbourList object!"<<endl;
	dumb=ReadThisLine(file,1,1,"total number of atoms","NeighbouList::LoadGrid",gridfilename);
	if (dumb!=natoms)
	{
		cout << "\nWARNING(" << ++warn << "): The number of atoms is "<<natoms<<" in .cfg file and "<<dumb<<\
			" in the .nei file are not the same!"<<endl;
		cout<<"\tNeighbourList has to be calculated!"<<endl;
		return(0);//load was not successful
	}
	
	dumb=ReadThisLine(file,1,1,"ntypes","NeighbouList::LoadGrid",gridfilename);
	if (dumb!=ntypes)
	{
		cout << "\nWARNING(" << ++warn << "): The number of atom types is "<<ntypes<<" in .cfg file and "<<dumb<<\
			" in the .nei file are not the same!"<<endl;
		cout<<"\tNeighbourList has to be calculated!"<<endl;
		return(0);//load was not successful
	}

	dumb=ReadThisLine(file,1,1,"max_gridatom","NeighbouList::LoadGrid",gridfilename);
	if (dumb!=max_gridatom)
	{
		cout << "\nWARNING(" << ++warn << "): The maximum number of atoms in a grid cell was " << RunParams::max_gridatom << endl;
		cout << "\tin the .dat file, which was reset to " << max_gridatom << " during the number of grid cell calculation," << endl;
		cout<<"\tand "<<dumb<<" in the .nei file. The last two are not the same! NeighbourList has to be calculated!"<<endl;
		return(0);//load was not successful
	}

	dumb=ReadThisLine(file,1,1,"ngridcells","NeighbouList::LoadGrid",gridfilename);
	if (dumb!=ngridcells)
	{
		cout << "\nWARNING(" << ++warn << "):  The number of grid cell is "<<ngridcells<<" calculated by the program, and ";
		cout<<"\t"<<dumb<<" in the .nei file are not the same! NeighbourList has to be calculated!"<<endl;
		return(0);//load was not successful
	}
	
	SkipLine(file,gridfilename,1);
	SkipLine(file,gridfilename,1);

	//loading the indices of the gridcells for eech atoms
	pint=grid_index;
	for (i=0;i<natoms;i++)
	{
		IntToStr(&conv_numb1,i+1);
		mystrcpy(conv_tx,100,"grid index x-axis for atom ");
		mystrcpy(conv_ty, 100, "grid index y-axis for atom ");
		mystrcpy(conv_tz, 100, "grid index z-axis for atom ");
		mystrcat(conv_tx, 100, conv_numb1);
		mystrcat(conv_ty, 100, conv_numb1);
		mystrcat(conv_tz, 100, conv_numb1);
		*pint++=ReadThisLine(file,3,1,conv_tx,"NeighbouList::LoadGrid");//grid index in x direction 
		*pint++=ReadThisLine(file,3,1,conv_ty,"NeighbouList::LoadGrid");//grid index in y direction
		*pint++=ReadThisLine(file,1,1,conv_tz,"NeighbouList::LoadGrid",gridfilename);//grid index in z direction
	}
	SkipLine(file,gridfilename,1);
	SkipLine(file,gridfilename,1);

	//Loading the number and the indices of the atoms belonging to each grid cells
	//order:x0-y0-z0:atom(0)->atom(cellcount-1),x1-y0-z0:atom(0)->atom(cellcount-1),x2-y0-z0:atom(0)->atom(cellcount-1)...
	pint=atom_indices;
	for (iz=0;iz<ngridcells;iz++)//grid index in the direction of z-axis
	{
		IntToStr(&conv_numb3,iz+1);
		mystrcpy(conv_tz, 100, ",");
		mystrcat(conv_tz, 100, conv_numb3);
		for (iy=0;iy<ngridcells;iy++)//grid index in the direction of y-axis
		{
			IntToStr(&conv_numb2,iy+1);
			mystrcpy(conv_ty, 100, ",");
			mystrcat(conv_ty, 100, conv_numb2);
			mystrcat(conv_ty, 100, conv_tz);
			for (ix=0;ix<ngridcells;ix++)//grid index in the direction of x-axis
			{
				IntToStr(&conv_numb1,ix+1);
				mystrcpy(conv_tx, 100, conv_numb1);
				mystrcat(conv_tx, 100, conv_ty);
				mystrcpy(conv_t, 100, "grid count for cell ");
				mystrcat(conv_t, 100, conv_tx);
								
				GridCount(ix,iy,iz)=ReadThisLine(file,3,1,conv_t,"NeighbouList::LoadGrid");
				count=GridCount(ix,iy,iz);
				if (count==0)
					SkipLine(file,gridfilename);//there is no more data in this line
				for (i=0;i<count;i++)//atom indices belonging to this grid
				{

					IntToStr(&conv_numb4,i+1);
					mystrcpy(conv_t, 100, conv_numb4);
					mystrcat(conv_t, 100, ". atom index for gridcell ");
					mystrcat(conv_t, 100, conv_tx);
					if (i<count-1)
						*pint++=ReadThisLine(file,3,1,conv_t,"NeighbouList::LoadGrid");
					else
					{
						//last number of the line
						if (ix==ngridcells-1 && iy==ngridcells-1 && iz==ngridcells-1)//last line of the file
							*pint++=ReadThisLine(file,1,1,conv_t,"NeighbouList::LoadGrid");
						else
							*pint++=ReadThisLine(file,1,1,conv_t,"NeighbouList::LoadGrid",gridfilename);
					}
				}
				pint+=max_gridatom-count;//setting to the beginning of the next grid cell
			}
		}

	}
	if (conv_numb1 != NULL)
		delete[] conv_numb1;
	if (conv_numb2 != NULL)
		delete[] conv_numb2;
	if (conv_numb3 != NULL)
		delete[] conv_numb3;
	if (conv_numb4 != NULL)
		delete[] conv_numb4;
	next_line_pos = -1;
	if (CheckReadFileState(file, "NeighbourList::LoadGrid", gridfilename))
	{
		file.close();
		if (::debug)
			cout << "Loading of the grid was successful" << endl;
		return(1);//load was successful
	}
	else
	{
		cout<<"\tLoading of the grid was not successful, neighbour list will be calculated!"<<endl;
		file.clear();
		file.close();
		return(0);//loading was not successful
	}
};

//------determining, whether the move is acceptable on the basis of cutoffs-------------
//--if below cutoff, it is checked, if it an FNC pair, as FNC has priority-------
//--in case of flexible molecules (fnc=4) if bonded distance or 1-3 distance of angle is below cutoff, still accepted 
//----if too close is still below cutoff, increased distance move is accepted----------
//---returning 1 if acceptable, zero otherwise-----------------
//--The concept is, that the grid is updated, but the neighbourlist is not yet, using the grid directly-- 
int NeighbourList::CheckCutoff(Move &move)
{
	//parameters
	//- move: move object 

	int imoved,jmoved,i,inb,gx,gy,gz,ix,iy,iz;
	int acc;
	int fnc_ind;//indicator, that the too close pair is an FNC pair
	int skip_cycle;//skip the cycle indicator
	int moved_tooclose;//flag to indicate, whther the imoved moved atom is a too close
	int neigh_moved_tooclose;//flag to indicate, whether the neighbour atom is a moved  tooclose atom
	int *atom_index,*moved_type;
	int ipartial,neightype;
	int n_atom;//number of atoms in the given grid cell
	int gridx,gridy,gridz;//index of the grid cell with the central atom
	int offset;//offset of the given atom
	int neigh_ind,neigh_offset,*p_fnc_neighind;//index and offset of the neighbour, pointer for the index of the FNC neighbour
	int	*p_fnc_neighind2;
	double x,y,z,*oldpos;//coordinates of the central atom
	double dx, dy, dz, d_square,dold_square;//distance components, squared distance
	double *x_max2;//pointer to reduced cutoff squared
	double *neigh_old_coord;//pointer to the old positions of the neighbour atom
	
	atom_index=move.indices;
	moved_type=move.types;
	x_max2=rundat.pcutsq;

	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		offset=3* *atom_index;
	
					
		//atoms can be below cutoff in the starting configuration, even if the moveout option is not used
		//in this case only those moves of this atoms can be accepeted, where they move above cutoff!
		moved_tooclose=0;//the moved not too close atom by default
		if (move.tcl_ind[imoved]>-1)//if the imoved-th atom is a too close atom
			moved_tooclose=1;//it can happen, only if the moveout option is on
			
		//determine the grid cell of the central atom 
		gridx=*(grid_index+offset);
		gridy=*(grid_index+offset+1);
		gridz=*(grid_index+offset+2);
		//collecting the coordinates of the moved atom
		x=*(config.positions+offset);
		y=*(config.positions+offset+1);
		z=*(config.positions+offset+2);

		//Normally the cells to check should go from gz=gridz-ncells_cutoff to gz<=gridz+ncells_cutoff, if there is a lot of gridcells, and the max_cutoff
		//would imply the checking only a few gridcells in each direction.
		//If the number of grid cells in any direction is small, or the max_cutoff and consequently the number of cells to check in a direction is large,
		//it can happen, that for example the gridz-ncells_cutoff would overlap with gz<=gridz+ncells_cutoff after taking into account the periodicity.
		//This would cause, that the overlapping cells would be checked twice, which is a waste of time, and in case of the moveout option it can cause
		//error during the update of the tclpair_ind array.
		//Therefore overlapping has to be prevented!
		for (gz=gridz-ncells_cutoff;gz<=gridz+ncells_cutoff-overlap;gz++)
		{
			//to take into account the periodicity
			iz=gz;
			if (gz<0)
				iz=gz+ngridcells;
			if (gz>ngridcells-1)
				iz=gz-ngridcells;

			for (gy=gridy-ncells_cutoff;gy<=gridy+ncells_cutoff-overlap;gy++)
			{
				//to take into account the periodicity
				iy=gy;
				if (gy<0)
					iy=gy+ngridcells;
				if (gy>ngridcells-1)
					iy=gy-ngridcells;

				for (gx=gridx-ncells_cutoff;gx<=gridx+ncells_cutoff-overlap;gx++)
				{
					//to take into account the periodicity
					ix=gx;
					if (gx<0)
						ix=gx+ngridcells;
					if (gx>ngridcells-1)
						ix=gx-ngridcells;

					n_atom=GridCount(ix,iy,iz);//number of atoms in this grid

					for (i=0;i<n_atom;i++)
					{
						neigh_ind=AtomIndex(ix,iy,iz,i);//index of the neighbour
						neigh_offset=3*neigh_ind;//offset of the neighbour
						neigh_old_coord=config.positions+neigh_offset;//beginning of the old coordinates of the neighbour in the 
											//config.positions array

						//check, whether the neighbour is not the imoved moved atom itself
						//check too, if the neighbour is not the previously checked moved atoms, for which all 
						//the necessary neighbours has been already checked to make sure, that each atom pair is checked
						//only once, important for the update of the tooclose list
						skip_cycle=0;
						for (jmoved=0; jmoved<=imoved; jmoved++)
						{
							if (neigh_ind==move.indices[jmoved])
								skip_cycle=1;//set, that not to calculate for this pair
						}
						if (skip_cycle)
							continue;//skip, this is the moved atom, or this pair was already calculated

						//the moved-moved pairs with neighbour having greater imoved index, than the present will be calculated,
						//set the flag
						neigh_moved_tooclose=0;//the neighbour is not a moved tooclose atom by default
						for (jmoved=imoved+1; jmoved<Move::tot_moved_atoms; jmoved++)
						{
							if (neigh_ind==move.indices[jmoved])
							{
								if (move.tcl_ind[jmoved]>-1)//if the jmoved-th atom is a too close atom
									neigh_moved_tooclose=1;//it can happen, only if the moveout option is on
								
								neigh_old_coord=move.oldpos+3*jmoved;//reset it, use the old coordinates from the oldpos array,
												//as this is a moved atom too, and the config.positions array contains the
												//new positions for the moved atoms
								break;//break the jmoved cycle, to preserve the value of jmoved
							}
						}

						
						//type of the neighbour atom
						neightype=0;
						while (neigh_ind>=SimpleCfg::cumul[(neightype) +1])
							neightype++;
						//neightype now points to the type of the neighbour atom
						//determining the partial index
						ipartial=(*moved_type<=neightype ? (*moved_type*ntypes-(*moved_type*(*moved_type+1)/2)+ neightype) : \
								(neightype*ntypes-(neightype*(neightype+1)/2)+ *moved_type));
								
						dx=x-*(config.positions+neigh_offset);
						dy=y-*(config.positions+neigh_offset+1);
						dz=z-*(config.positions+neigh_offset+2);
						//taking into account the periodical boundary conditions 
						if (dx>=1)
							dx-=2;
						else
							if (dx<=-1)
								dx+=2;
						if (dy>=1) 
							dy-=2;
						else
							if (dy<=-1) 
								dy+=2;
						if (dz>=1)
							dz-=2;
						else 
							if (dz<=-1)
								dz+=2;
						d_square=dx*dx+dy*dy+dz*dz;

						if (d_square<x_max2[ipartial])//distance is below cutoff 
						{
							acc=0;
							if (rundat.fnc==4 || rundat.fnc_conflict)//There is a conflict between the cutoffs and FNC
										//or in flexible molecules
							{
								//check, if it is an FNC neighbour, or in case of a BOND 
								p_fnc_neighind=fnc[0].neighbours+*(fnc[0].converter+ *atom_index);
								for (inb=0;inb<*(fnc[0].numbneigh+*atom_index);inb++)//for each neighbour of this central atom
								{		
									if (*p_fnc_neighind+FNC_POT::index_offset==neigh_ind)
									{
										acc=1;//this is FNC neighbour
										break;//(inb cycle)
									}
									p_fnc_neighind++;//next neighbour index
								}//next FNC neighbour
							}

							//check, if they are 1-3 atoms (two neighbour for the middle atom of the angle) in case of flexible ANGLE 
							//the 1-2 or 2-3 pairs are not checked, as it was checked during BOND checking
							if (rundat.fnc==4 && acc==0)
							{
								p_fnc_neighind=fnc[1].neighbours+*(fnc[1].converter+ *atom_index);
								p_fnc_neighind2=fnc[1].neighbours2+*(fnc[1].converter+ *atom_index);
								for (inb=0;inb<*(fnc[1].numbneigh+*atom_index);inb++)//for each neighbour of this central atom
								{		
									if (*p_fnc_neighind<0 && *p_fnc_neighind2 + FNC_POT::index_offset ==neigh_ind)//- sign in neighbour array shows that that atom is the middle in the angle
									{
										acc=1;//this is 1-3 atom in the angle, the potential will take care of it keeping them at a proper distance
										break;//(inb cycle)
									}
									p_fnc_neighind++;//next neighbour index
									p_fnc_neighind2++;//next neighbour index
								}//next FNC neighbour
							}
							//no FNC, or it was not an FNC neighbour, check if it is a tooclose pair
							//for this the imoved moved atom has to be a tooclose atom
							if (acc==0 && moved_tooclose)
							{
								//old distance can be below cutoff only for too close pairs, no need to check for that explicitly
								//(and for conflicting FNC pairs, but that was already excluded)
								oldpos=move.oldpos+3*imoved;//old x coordinate of the moved atom 
								//check, if it a too close,, and the distance increased
								dx=*oldpos++ -*neigh_old_coord;
								dy=*oldpos++ -*(neigh_old_coord+1);
								dz=*oldpos-*(neigh_old_coord+2);
								//taking into account the periodical boundary conditions 
								if (dx>=1)
									dx-=2;
								else
									if (dx<=-1)
										dx+=2;
								if (dy>=1) 
									dy-=2;
								else
									if (dy<=-1) 
										dy+=2;
								if (dz>=1)
									dz-=2;
								else 
									if (dz<=-1)
										dz+=2;
								dold_square=dx*dx+dy*dy+dz*dz;

								if (d_square<dold_square)//the distance decreased
								{
									return(0);//the move is not acceptable
								}

								acc=1;//the distance increased or not changed for the too close atom
								if (moved_tooclose)//the imoved atom is a tooclose
									move.tooclose_remove[move.tcl_ind[imoved]] = 0;//cannot be removed from the "tooclose" list
								
								//If the neighbour atom is a moved atom too, then its tooclose_remove has to be updated 
								if (neigh_moved_tooclose)//the neighbour atom is the jmoved moved atom, and it is a tooclose atom
									move.tooclose_remove[move.tcl_ind[jmoved]] = 0;//cannot be removed from the "tooclose" list
							}
							if (!acc)
							{
								
								return(0);//the move is not acceptable, as it is not a conflicting FNC or a increased_dist tooclose with moveout option on
							}
						}//end if new distance is below cutoff
						else
						{

							//new distance is acceptable
							//check if it is a tooclose pair, for this the imoved moved atom has to be a tooclose atom
							if (moved_tooclose)
							{
								//see, if it is a too close pair, (old distance is below cutoff)
								oldpos=move.oldpos+3*imoved;//old x coordinate of the moved atom 
								//check, if it a too close,, and the distance increased
								dx=*oldpos++ -*(neigh_old_coord);
								dy=*oldpos++ -*(neigh_old_coord+1);
								dz=*oldpos-*(neigh_old_coord+2);
								//taking into account the periodical boundary conditions 
								if (dx>=1)
									dx-=2;
								else
									if (dx<=-1)
										dx+=2;
								if (dy>=1) 
									dy-=2;
								else
									if (dy<=-1) 
										dy+=2;
								if (dz>=1)
									dz-=2;
								else 
									if (dz<=-1)
										dz+=2;
								dold_square=dx*dx+dy*dy+dz*dz;
								if (dold_square<x_max2[ipartial])//old distance is below cutoff, the new is above
								{
									fnc_ind=0;//default, this is not an fnc pair
									if (rundat.fnc==4 || rundat.fnc_conflict)//There is a conflict between the cutoffs and FNC
										//or in flexible molecules
									{
										//check whether this is an FNC pair or a BOND, as these are not counted among the too close
										p_fnc_neighind=fnc[0].neighbours+*(fnc[0].converter+*atom_index);
										for (inb=0;inb<*(fnc[0].numbneigh+*atom_index);inb++)
										{							
											//for each neighbour of this central atom
											if (*p_fnc_neighind+FNC_POT::index_offset==neigh_ind)
											{
												fnc_ind=1;//this is a FNC neighbour, not included in the neighbourlist
												break;//(inb cycle)
											}
											p_fnc_neighind++;//next neighbour index
										}//next fnc neighbour
									}//end of if conflict
									//check, if they are 1-3 atoms (two neighbour for the middle atom of the angle) in case of flexible ANGLE 
									//the 1-2 or 2-3 pairs are not checked, as it was checked during BOND checking
									if (rundat.fnc==4)
									{
										p_fnc_neighind=fnc[1].neighbours+*(fnc[1].converter+ *atom_index);
										p_fnc_neighind2=fnc[1].neighbours2+*(fnc[1].converter+ *atom_index);
										for (inb=0;inb<*(fnc[1].numbneigh+*atom_index);inb++)//for each neighbour of this central atom
										{		
											if (*p_fnc_neighind<0 && *p_fnc_neighind2 + FNC_POT::index_offset ==neigh_ind)//- sign in neighbour array shows that that atom is the middle in the angle
											{
												fnc_ind=1;//this is 1-3 atom in the angle, the potential will take care of it keeping them at a proper distance
												break;//(inb cycle)
											}
											p_fnc_neighind++;//next neighbour index
											p_fnc_neighind2++;//next neighbour index
										}//next FNC neighbour
									}
									if (fnc_ind==0)//this is not an fnc_neighbour, too close list has to be updated
									{
										//too close list has to be updated
										move.tclpair_ind[move.tcl_ind[imoved]*move.max_pairs+move.ntcl_pair[move.tcl_ind[imoved]]]=neigh_ind;//index of the neighbour atom (between 0 and ntotal)
										move.ntcl_pair[move.tcl_ind[imoved]]++;//increase the number of pairs moved above cutoff 
									}
								}
							}//end of one of the pair is a too close atom
						}//end if new distance is above cutoff
					}//cycling through the atoms of this grid
				}//going throgh the grid cells in x-axis direction
			}//going throgh the grid cells in y-axis direction
		}//going through the grid cells in z-axis direction
		atom_index++;
		moved_type++;
	};//going through the moved atoms, imoved cycle
	
	return(1);//there are no atoms below cutoff, which are not FNC neighbours or too close atoms with increased distance
};

//--------------determining the neighbours of atom atom_index inside xmax regardless their type------------------
//----------------------------placing them into the array specified by plist
//normally working with the config positions array, which in the loop contains the new coordinates
//used for COS, CONC, SNC with parameters dist2_list, v_comp and no coord_list,
//	 in this case called by serial calculation

//used for ANN potential neighbour list with  coord_list and no dist2_list and v_comp
//	parallel, called by threads
//	as it is possible that in case of larger inhomogenities, the neighbour list has to be resized, max neighbours have to be passed,
//	as in case of ANN pot the max neighbours can differ for the different threads, as it would be unfeasible to resize all the separate arrays
//	in case of ANN used by the different threads, therefore, no global variable for max_neigh can be used
//	as the neighboutlist is not stored, for the calculation of the old neighboulist the old coordinates of the moved atoms have to be used!!!
int NeighbourList::CalcNeighList(int atom_index, int *ind_list, int *type_list, double *dist2_list, double *v_comp, double *dist_list, int *my_max_neigh, int id, double *coord_list, Move *move, int imoved)
{
	//parameters
	//- atom_index:  	index of the central atom for which the neighbourlist is calculated
	//- *ind_list:   	pointer to the starting element of the array, where the neighbours list is stored 
	//- *type_list:	 	pointer to the starting element of the array, where the neighbour type list stored
	//- *dist2_list: 	USED FOR COS, CONC, SNC: pointer to the starting point of the squared reduced distance list for the neighbours
	//- *v_comp: 	 	USED FOR COS, CONC, SNC: pointer to the starting point of the reduced vector components of the Central->Neighbour vector
	//- *dist_list		default: nullptr, USED FOR CNC detail: pointer to the starting point of the central-neighbour distance in A
	//- *my_max_neigh:	default: NeighbourList::max_neigh, FOR CNC: maximum number of cnc neighbours, FOR ANN pot, the maximum number of neighbours for this thread 
	//- id:				default: 0,       USED FOR CNC detail: constarint index, USED FOR ANN pot: thread id 
	//- *coord_list: 	default: nullptr, ONLY USED FOR ANN pot: pointer to the starting point of the aenet-format coordinates of the neighbours
	//- *move:			default: nullptr, ONLY USED FOR ANN pot: pointet to the move object, only needed in the loop to collect the old coordinates
	//- int imoved:		default: -1,      ONLY USED FOR ANN pot: index of the atom_index atom in the moved arrays of the move object, only needed in the loop to collect the old coordinates
	
	// the neighbourlist is calculated for the same xmax distance for each atom 

	int i,j,gx,gy,gz,ix,iy,iz;
	int mode=0;//0:COS, CONC, SNC; 1: ANN, 2: CNC detail
	int n_atom;//number of atoms in the given grid cell
	int gridx,gridy,gridz;//index of the grid cell with the central atom
	int neigh=0;//number of neighbours
	int offset;//offset of the given atom
	int neigh_ind;//index of the neighbour
	int *pneigh_type;//pointer to the neighbour's type
	int my_ncells_nlist=ncells_nlist;//default, mode=0
	int my_overlap_nlist=overlap_nlist;//default, mode=0
	double central_coord[3];//coordinates of the central atom
	#ifdef _AENET
	double aenet_central[3];
	#endif
	double dr[3], d_square;//distance components, squared distance
	double my_xmax2=xmax2;//default for COS, CONC, SNC
	double *pc;//v_comp for norma, coord_list for aenet neighbours
	double *pneigh_coord;//pointer to teh neighbour coordinates
	bool good_type = 1;//normally neighbour list is calculated regardless the type, but in case of CNC, type matters
#ifdef _AENET
	int set_old=0;//set the old coordinates for the moved atoms
	int im;
#endif
	pc=v_comp;//by default
//the validity of my_max_neigh is not checked!
#ifndef _AENET
	if ((dist2_list==nullptr || v_comp==nullptr) && dist_list==nullptr)
	{
		cout << "\n*****ERROR*****" << endl;
		cout << "NeighbourList::CalcNeighList was called without dist2_list or v_comp, and also without dist_list!" << endl;
		cout << "Either use dist2_list and v_comp and dist_lis=nullptr for COS, CONC and SNC, "  << endl;
		cout<<"or call with dist2_list=nullptr, v_comp=nullptr and dist_lis for CNC detail!"<<"\nCannot run this way, exiting..."<<endl;
		CleanExit();
	}
	

#else
	if (((dist2_list==nullptr || v_comp==nullptr) && dist_list==nullptr) && coord_list==nullptr)
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"NeighbourList::CalcNeighList was called without dist2_list or v_comp and coord_list as well, cannot run this way, exiting..."<<endl;
		CleanExit();
	}
	if ((move!=nullptr && imoved==-1) || (move==nullptr && imoved!=-1))
	{
		cout << "\n*****ERROR*****" << endl;
		cout<<"NeighbourList::CalcNeighList was called without initializing both move and moved_ind, exiting..."<<endl;
		CleanExit();
	}
	if ((dist2_list==nullptr && v_comp==nullptr) && coord_list!=nullptr)
	{
		mode=1;//ANN neighbour calc
		my_xmax2=Aenet::red_cutoff_sq;
		my_ncells_nlist=Aenet::ncells_nlist;
		my_overlap_nlist=Aenet::overlap_nlist;
	}
#endif
	if ((dist2_list == nullptr && v_comp == nullptr) && dist_list != nullptr)
	{
		mode = 2;//CNC detail
		my_xmax2 = CoordNumbConst::maxdistsq;//this is the largest, it will checed again, when type is known
		good_type = 0;//it is set for each neighbour
		my_ncells_nlist=CoordNumbConst::ncells_nlist;
		my_overlap_nlist=CoordNumbConst::overlap_nlist;
	}
	offset=3*atom_index;
	//determine the grid cell of the central atom 
	gridx=*(grid_index+offset);
	gridy=*(grid_index+offset+1);
	gridz=*(grid_index+offset+2);
	//collecting the coordinates of the moved atom
	for (j=0;j<3;j++)
	{
		central_coord[j]=*(config.positions+offset+j);
	#ifdef _AENET
		if (imoved!=-1)//calulate for the old neighbours of the moved atom in loop 
		{
			central_coord[j]=move->oldpos[imoved*3+j];
			if (move->tot_moved_atoms>1)//the neighbour can be a moved atom as well
				set_old=1;
		}
		aenet_central[j]=(central_coord[j]+1)*SimpleCfg::boxedge;
	#endif
	}


	//Normally the cells to check should go from gz=gridz-ncells_nlist to gz<=gridz+ncells_nlist, if there is a lot of gridcells, and the xmax distance
	//determining the radius of the sphere to check in would imply the checking only a few gridcells in each direction.
	//If the number of grid cells in any direction is small, or the xmax and consequently the number of cells to check in a direction is large,
	//it can happen, that for example the gridz-ncells_nlist would overlap with gz<=gridz+ncells_nlist after taking into account the periodicity.
	//This would cause, that the overlapping cells would be checked twice, adding the same neighbours twice to the neighbourlist
	//overfilling most probably the list and the program would stop with the message of too many neighbours.
	//Therefore overlapping has to be prevented!
	for (gz=gridz-my_ncells_nlist;gz<=gridz+my_ncells_nlist-my_overlap_nlist;gz++)
	{
		//to take into account the periodicity
		iz=gz;
		if (gz<0)
			iz=gz+ngridcells;
		if (gz>ngridcells-1)
			iz=gz-ngridcells;

		for (gy=gridy-my_ncells_nlist;gy<=gridy+my_ncells_nlist-my_overlap_nlist;gy++)
		{
			//to take into account the periodicity
			iy=gy;
			if (gy<0)
				iy=gy+ngridcells;
			if (gy>ngridcells-1)
				iy=gy-ngridcells;
			for (gx=gridx-my_ncells_nlist;gx<=gridx+my_ncells_nlist-my_overlap_nlist;gx++)
			{
				//to take into account the periodicity
				ix=gx;
				if (gx<0)
					ix=gx+ngridcells;
				if (gx>ngridcells-1)
					ix=gx-ngridcells;

				n_atom=GridCount(ix,iy,iz);//number of atoms in this grid

				for (i=0;i<n_atom;i++)
				{
					neigh_ind=AtomIndex(ix,iy,iz,i);//index of the neighbour
					if (neigh_ind==atom_index)//this is the same atom
						continue;
					
					pneigh_coord=config.positions+3*neigh_ind;
#ifdef _AENET
					if (set_old)//calculating for the old neighbours
					{
						for (im=0;im<move->tot_moved_atoms;im++)
						{//see if the neighbour is a moved atom as well, as in this case the old coordinates has to be uesd
							if (neigh_ind==move->indices[im])//the neighbour is a moved atom as well, collectthe old coordinates
					 		{
								pneigh_coord=move->oldpos+3*im;
								break;
							}
						}
					}
#endif
					if (mode == 2)
						good_type = 0;
					d_square=0;
					for (j=0;j<3;j++)
					{
						dr[j]=pneigh_coord[j]-central_coord[j];
						//taking into account the periodical boundary conditions 
						if (dr[j]>=1)
							dr[j]-=2;
						else
							if (dr[j]<=-1)
								dr[j]+=2;
						d_square+=pow(dr[j],2);
					}
					
					if (d_square<=my_xmax2)//to make the least change for CNC calculation, the neighbour is added first to the least, and if not the right
						//type, the counter is not increased
					{
						if (neigh==*my_max_neigh)
						{
							if (::debug)
							{
								cout << "NeighbourList::CalcNeighList: Too many neighbours (" << neigh << ") for the " << atom_index << ". atom!" << endl;
								cout << "Resize arrays and continue..." << endl;
							}
							if (mode==0)
							{
								ResizeArrays(atom_index, &ind_list, &type_list, &dist2_list, &v_comp);//resizing the arrays having max_neigh in their dimension
								pc=v_comp+3*neigh;
							}
							else if (mode == 2)
							{
								CoordNumbConst::ResizeCoordArrays(&ind_list, &type_list, &dist_list);//this might be unnecessary, but does not matter
							}
#ifdef _AENET
							if (mode==1)
								thread->ResizeAenetArrays(id, &ind_list, &type_list,&coord_list);
#endif

						}

						//Determining the type of the neighbour atoms
						pneigh_type=type_list+neigh;//has to be set as it could have been resized
						*pneigh_type=0;
						while (neigh_ind>=SimpleCfg::cumul[(*pneigh_type) +1])
							(*pneigh_type)++;
							//*pneigh_type now points to the type of the neighbour atom
						if (mode == 2)//see, if this is a valid neighbour type for the cnc constraint
						{
							for (int ineightype = 0; ineightype < CoordNumbConst::n_neightype[id]; ineightype++)
								if (CoordNumbConst::neighbours[CoordNumbConst::cum_n_neightype[id] + ineightype] == *pneigh_type)
								{
									if (d_square<CoordNumbConst::udmaxsq[CoordNumbConst::cum_n_neightype[id] + ineightype] && d_square > CoordNumbConst::udminsq[CoordNumbConst::cum_n_neightype[id] + ineightype])
										good_type = 1;//this is a good type, and in range
									break;
								}
						}
#ifdef _AENET
						if (mode==1)
							(*pneigh_type)++;//to have the type in the 1->ntypes range for the fortran
#endif
						*(ind_list+neigh)=neigh_ind;//index of the neighbour 
						//setting the other list elements for this neighbour
						if (mode==0)
						{
							*(dist2_list+neigh)=d_square;//squared distance
							for (j=0;j<3;j++)
								*pc++=dr[j];//vector components
						
						}
						else if (mode == 2)
						{
							*(dist_list + neigh) = sqrt(d_square)* SimpleCfg::boxedge;//distance in A
						}
#ifdef _AENET
						else
						{
							pc=coord_list+3*neigh;//has to be set as it could have been resized
							for (j=0;j<3;j++)
							{
								*pc=(pneigh_coord[j]+1)*SimpleCfg::boxedge;//fill the neigh coordinates in A
								/**pc*=1e6;
                				*pc=(longint)(*pc);
                				*pc/=1e6;*/
                				//have to use the neighbour image inside cutoff
								double diff=aenet_central[j]-*pc;
								if (diff>SimpleCfg::boxedge)
									*pc+=2*SimpleCfg::boxedge;
								else if (diff<-SimpleCfg::boxedge)
									*pc-=2*SimpleCfg::boxedge;
								pc++;
							}
						}
#endif
						if (!(mode==2 && good_type==0))
							neigh++;//increase the counter, only not, if CNC calculation and this is a bad type
					}
				}//cycling through the atoms of this grid
			}//going throgh the grid cells in x-axis direction
		}//going throgh the grid cells in y-axis direction
	}//going throgh the grid cells in z-axis direction
	return(neigh);//number of neighbours
};

//calculating the neighbourlist for each atom belonging to the types involved in CosDistrCons, CommonNeighConst, SecondNeighConst
//all the neighbours inside the specified range are renturned regardless of their type, as atom 'i' may be
//involved in more, than one constraint, as central or neighbour
void NeighbourList::InitNeighList()
{
	int itype,i;
	int offset,rel_index;
	int *pnneigh;

	pnneigh=nneigh;//number of neighbours for each atom
	for (itype=0;itype<ntypes;itype++)
	{
		if (consttype>>itype & 1)
		{
			offset=0;//the beginning of the neighbourlist for this atom relative to the beginning of this type
			rel_index=0;
			for (i=SimpleCfg::cumul[itype];i<SimpleCfg::cumul[itype+1];i++)//going through the atoms of this type
			{
				offset=rel_index*max_neigh;//set offset, max_neigh might cjange, if the arrays are resized in CalcNeighList
				
				*pnneigh=CalcNeighList(i,neighlist_finder[itype]+offset,neightype_finder[itype]+offset,dsq_finder[itype]+offset,vector_finder[itype]+offset*3);
				rel_index++;
				pnneigh++;
			}
		}//end if this type is involved
	}//end of itype cycle
	
};

//----------------updating the neighbourlist to reflect the move--------------------
//-------------------it assumes, that the grid was already updated------------------- 
void NeighbourList::UpdateNeighList(Move &move)
{
	int imoved,ichanged,inb1,i;
	int nmoved_to_update=0;//number of moved atoms, which are involved in cosine or common neighbour constraint as has to be updated
	int index,offset;
	int *moved_index,*moved_type;
	int *pchanged_ind,*pchanged_type;
	int *pnneigh1, *plist1, *ptypelist1;
	int *pchanged_nneigh, *pchanged_neigh_list, *pchanged_neigh_type;
	double  *pdistsq1, *pvec_comp1;
	double  *pchanged_dsq, *pchanged_vector_comp;

	nchanged=0;//number of atoms effected by the move
	moved_index=move.indices;
	moved_type=move.types;
	
	pchanged_ind=changed_indices;
	pchanged_type=changed_types;
	pchanged_nneigh=changed_nneigh;

	//First store the old values for the moved atoms
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		if (!(consttype>>*moved_type & 1))
		{
			moved_type++;
			moved_index++;
			continue;//this moved atom is not involved in any constraint, go to next
		}
		nmoved_to_update++;//this moved atom has to be updated

		//storing the indices and types of the moved atoms, this is done to keep the parameters of all the changed atoms together
		*pchanged_ind++=*moved_index;
		*pchanged_type++=*moved_type;
		//set the pointers for the values to be copied
		index=*moved_index-config.cumul[*moved_type];//the index of the moved atom in its own type
		offset=index*max_neigh;//offset considering the maximum number of neighbours
		
		pnneigh1=nneigh_finder[*moved_type]+index;//sets the pointer to the number of neighbours of the moved atom
		plist1=neighlist_finder[*moved_type]+offset;//sets the pointer to the beginning of the neigh_list of the moved atom
		ptypelist1=neightype_finder[*moved_type]+offset;//sets the pointer to the beginning of the neigh_type of the moved atom
		pdistsq1=dsq_finder[*moved_type]+offset;//sets the pointer to the beginning of the dsq of the moved atom
		pvec_comp1=vector_finder[*moved_type]+offset*3;//sets the pointer to the beginning of the vector components of the moved atom

		//set the pointers for the arrays storing the original values
		offset=nchanged*max_neigh;//offset for the changed lists
		*pchanged_nneigh++=*pnneigh1;//sets the number of neighbours
		pchanged_neigh_list=changed_neigh_list+offset;//pointer to the indices of the neighbours
		pchanged_neigh_type=changed_neigh_type+offset;//pointer to the types of the neighbours
		pchanged_dsq=changed_dsq+offset;//pointer to the squared distance for the neighbours
		pchanged_vector_comp=changed_vector_comp+offset*3;//pointer to the vector components

		//copying the values
		for (inb1=0;inb1<*pnneigh1;inb1++)
		{

			*pchanged_neigh_list++=*plist1++;
			*pchanged_neigh_type++=*ptypelist1++;
			*pchanged_dsq++=*pdistsq1++;
			for (i=0;i<3;i++)
				*pchanged_vector_comp++=*pvec_comp1++;
		}//end of neighbour atom cycle for imoved
		nchanged++;
		
		moved_index++;
		moved_type++;
	}//end of imoved cycle
	
	// store the old values for the old neighbours of the moved atoms
	StoreList(move);

	//Now recalculate the neighbourlist for the moved atoms
	moved_index=move.indices;
	moved_type=move.types;
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		if (!(consttype>>*moved_type & 1))
		{
			moved_type++;
			moved_index++;
			continue;//this moved atom is not involved in any constraint, go to next
		}
		index=*moved_index-config.cumul[*moved_type];//the index of the moved atom in its own type
		offset=index*max_neigh;//offset considering the maximum number of neighbours
		*(nneigh_finder[*moved_type]+index)=CalcNeighList(*moved_index,neighlist_finder[*moved_type]+offset,neightype_finder[*moved_type]+\
											offset,dsq_finder[*moved_type]+offset,vector_finder[*moved_type]+offset*3);
		moved_index++;
		moved_type++;
	}//end of imoved cycle
	
	// store the old values for the new neighbours
	StoreList(move);

	//Now recalculate the list for the old and new neighbours
	pchanged_ind=changed_indices+nmoved_to_update;
	pchanged_type=changed_types+nmoved_to_update;
	for (ichanged=nmoved_to_update;ichanged<nchanged;ichanged++)//the first nmoved_to_update, it was recalculated for them already, if they had to be updated
	{
		index=*pchanged_ind-config.cumul[*pchanged_type];//the index of the neighbour of the moved atom in its own type
		offset=index*max_neigh;//offset considering the maximum number of neighbours
		*(nneigh_finder[*pchanged_type]+index)=CalcNeighList(*pchanged_ind,neighlist_finder[*pchanged_type]+offset,neightype_finder[*pchanged_type]+\
											offset,dsq_finder[*pchanged_type]+offset,vector_finder[*pchanged_type]+offset*3);
		pchanged_ind++;
		pchanged_type++;
	}


};

//store the values of the neighbourlist connected arrays for the neighbours of the moved atoms
void NeighbourList::StoreList(Move &move)
{
	int imoved,ichanged,inb1,inb2,i;
	int store;//boolean indicator
	int index,offset;
	int *moved_index,*moved_type;
	int *pnneigh1, *plist1, *ptypelist1;
	int *pnneigh2, *plist2, *ptypelist2;
	int *pchanged_nneigh, *pchanged_neigh_list, *pchanged_neigh_type;
	double  *pdistsq2, *pvec_comp2;
	double  *pchanged_dsq, *pchanged_vector_comp;

	moved_index=move.indices;
	moved_type=move.types;
	pchanged_nneigh=changed_nneigh+nchanged;
	for (imoved=0;imoved<Move::tot_moved_atoms;imoved++)
	{
		if (!(consttype>>*moved_type & 1))
		{
			moved_type++;
			moved_index++;
			continue;//this moved atom is not involved in any constraint, go to next
		}
		//store the old values for the neighbours
		index=*moved_index-config.cumul[*moved_type];//the index of the moved atom in its own type
		offset=index*max_neigh;//offset considering the maximum number of neighbours
		//set the pointers for the lists of the moved atom	
		pnneigh1=nneigh_finder[*moved_type]+index;//sets the pointer to the number of neighbours of the moved atom
		plist1=neighlist_finder[*moved_type]+offset;//sets the pointer to the beginning of the neigh_list of the moved atom
		ptypelist1=neightype_finder[*moved_type]+offset;//sets the pointer to the beginning of the neigh_type of the moved atom

		for (inb1=0;inb1<*pnneigh1;inb1++)//go through the neighbours of the moved atom
		{
			if (!(consttype>>*ptypelist1 & 1))//there is no neighbourlist for this neighbour, skip it
			{
				ptypelist1++;
				plist1++;
				continue;//this moved atom is not involved in any constraint, go to next
			}
			//as it is possible, that an atom is neighbour to more than one moved atoms, or especially in case of molecular move
			//the moved atoms can be each others neighbour, care is taken, that each atom is represented only once on the list
			store=1;//store this atom by default
			//check, whether it is already on the list
			for (ichanged=0;ichanged<nchanged;ichanged++)
			{
				if (changed_indices[ichanged]==*plist1)
				{	
					store=0;
					break;
				}
			}
			if (store)
			{
				//this atom should be stored, set the pointers for the values to be copied
				changed_indices[nchanged]=*plist1;//the index of this atom, to show, that it was stored
				changed_types[nchanged]=*ptypelist1;//the type of this atom
				index=*plist1-config.cumul[*ptypelist1];//the index of the neighbour of the moved atom in its own type
				offset=index*max_neigh;//offset considering the maximum number of neighbours
				
				pnneigh2=nneigh_finder[*ptypelist1]+index;//sets the pointer to the number of neighbours of the atom
				plist2=neighlist_finder[*ptypelist1]+offset;//sets the pointer to the beginning of the neigh_list of the atom
				ptypelist2=neightype_finder[*ptypelist1]+offset;//sets the pointer to the beginning of the neigh_type of the atom
				pdistsq2=dsq_finder[*ptypelist1]+offset;//sets the pointer to the beginning of the dsq of the atom
				pvec_comp2=vector_finder[*ptypelist1]+offset*3;//sets the pointer to the beginning of the vector components of the atom

				//set the pointers for the arrays storing the original values
				offset=nchanged*max_neigh;//offset for the changed lists
				*pchanged_nneigh++=*pnneigh2;//sets the number of neighbours
				pchanged_neigh_list=changed_neigh_list+offset;//pointer to the indices of the neighbours
				pchanged_neigh_type=changed_neigh_type+offset;//pointer to the types of the neighbours
				pchanged_dsq=changed_dsq+offset;//pointer to the squared distance for the neighbours
				pchanged_vector_comp=changed_vector_comp+offset*3;//pointer to the vector components

				//copying the values
				for (inb2=0;inb2<*pnneigh2;inb2++)
				{
					*pchanged_neigh_list++=*plist2++;
					*pchanged_neigh_type++=*ptypelist2++;
					*pchanged_dsq++=*pdistsq2++;
					for (i=0;i<3;i++)
						*pchanged_vector_comp++=*pvec_comp2++;
				}//end of neighbour atom cycle for imoved
				nchanged++;
			}//end of if this atom's values should be stored
			plist1++;
			ptypelist1++;
		}//end of neighbour cycle inb1
		moved_index++;
		moved_type++;
	}//end of imoved cycle
};

//------------if the move was rejected, reset the neighbour lists using the stored old values--------------
void NeighbourList::ResetList()
{
	int ichanged,inb1,i;
	int index,offset;
	int *pchanged_ind,*pchanged_type;
	int *pnneigh1, *plist1, *ptypelist1;
	int *pchanged_nneigh, *pchanged_neigh_list, *pchanged_neigh_type;
	double  *pdistsq1, *pvec_comp1;
	double  *pchanged_dsq, *pchanged_vector_comp;

	//Copy back the stored old values of the lists for the atoms effected by the move
	pchanged_ind=changed_indices;
	pchanged_type=changed_types;
	pchanged_nneigh=changed_nneigh;
	for (ichanged=0;ichanged<nchanged;ichanged++)
	{
		index=*pchanged_ind-config.cumul[*pchanged_type];//the index of the neighbour of the moved atom in its own type
		offset=index*max_neigh;//offset considering the maximum number of neighbours
	
		//restore the list, set the pointers for the lists of the effected atoms	
		pnneigh1=nneigh_finder[*pchanged_type]+index;//sets the pointer to the number of neighbours of the effected atom
		plist1=neighlist_finder[*pchanged_type]+offset;//sets the pointer to the beginning of the neigh_list of the effected atom
		ptypelist1=neightype_finder[*pchanged_type]+offset;//sets the pointer to the beginning of the neigh_type of the effected atom
		pdistsq1=dsq_finder[*pchanged_type]+offset;//sets the pointer to the beginning of the dsq of the effected atom
		pvec_comp1=vector_finder[*pchanged_type]+offset*3;//sets the pointer to the beginning of the vector components of the effected atom

		//set the pointers for the arrays storing the old values
		offset=ichanged*max_neigh;//offset for the changed lists
		*pnneigh1=*pchanged_nneigh++;//sets the number of neighbours
		pchanged_neigh_list=changed_neigh_list+offset;//pointer to the indices of the neighbours
		pchanged_neigh_type=changed_neigh_type+offset;//pointer to the types of the neighbours
		pchanged_dsq=changed_dsq+offset;//pointer to the squared distance for the neighbours
		pchanged_vector_comp=changed_vector_comp+offset*3;//pointer to the vector components

		for (inb1=0;inb1<*pnneigh1;inb1++)//go through the neighbours of the effected atom
		{
			*plist1++=*pchanged_neigh_list++;
			*ptypelist1++=*pchanged_neigh_type++;
			*pdistsq1++=*pchanged_dsq++;
			for (i=0;i<3;i++)
				*pvec_comp1++=*pchanged_vector_comp++;
		}//end of neighbour atom cycle for atom ichanged
		
		pchanged_ind++;
		pchanged_type++;
	}

};

//---------if due to large inhomogeneity in the sample the atom count of one gridcell or the size of the neighbourlist -------
//--------------would exceed the maximum, resize the array to be able to continue
void NeighbourList::ResizeIArray(int max, int new_max, int oldsize, int newsize, int **array, int **finder)
{
	int i,j;
	int offset;
	int *temp_array, *newp, *oldp;
	
	SetArraysize(&temp_array,oldsize*max,"temp_array","NeighbourList::ResizeIArray");//the original array will be copied here
	
	newp=temp_array;
	oldp=*array;
	for (i=0;i<oldsize*max;i++)
		*newp++=*oldp++;//copy the original to temp

	delete [] *array;
	SetArraysize(array,newsize*new_max,"array","NeighbourList::ResizeIArray");
	
	newp=*array;
	oldp=temp_array;
	for (i=0;i<oldsize;i++)
	{
		for (j=0;j<max;j++)
			*newp++=*oldp++;//copy the values to the new array
		newp+=new_max-max;//skip the new empty element
	}
	
	delete [] temp_array;

	offset=0;
	if (finder!=NULL)//set the finder, if there is any
	{
		for (i=0;i<ntypes;i++)
		{
			if (consttype>>i & 1)//this type is represented
			{
				finder[i]=*array+new_max*offset;
				offset+=SimpleCfg::pnatoms[i];//increase offset with the number of atoms of this type * max_neigh
			}
			else
				//this type is not represented, the finder will not be called for it normally
				finder[i]=NULL;

		}
	}
};

//---------if due to large inhomogeneity in the sample the atom count of one gridcell or the size of the neighbourlist -------
//--------------would exceed the maximum, resize the array to be able to continue
void NeighbourList::ResizeDArray(int max, int new_max, int oldsize, int newsize, int multiply, double **array, double **finder)
{
	int i,j;
	int offset;
	double *temp_array, *newp, *oldp;
	
	SetArraysize(&temp_array,oldsize*max*multiply,"temp_array","NeighbourList::ResizeDArray");//the original array will be copied here

	newp=temp_array;
	oldp=*array;
	for (i=0;i<oldsize*max*multiply;i++)
		*newp++=*oldp++;//copy the original to temp

	delete [] *array;

	SetArraysize(array,newsize*new_max*multiply,"array","NeighbourList::ResizeDArray");//creating the larger new array
	
	newp=*array;
	oldp=temp_array;
	for (i=0;i<oldsize;i++)
	{
		for (j=0;j<max*multiply;j++)
			*newp++=*oldp++;//copy the values to the new array
		newp+=multiply*(new_max-max);//skip the new empty element
	}
	
	delete [] temp_array;

	offset=0;
	if (finder!=NULL)//set the finder, if there is any
	{
		for (i=0;i<ntypes;i++)
		{
			if (consttype>>i & 1)//this type is represented
			{
				finder[i]=*array+new_max*offset*multiply;
				offset+=SimpleCfg::pnatoms[i];//increase offset with the number of atoms of this type * max_neigh *multiply
			}
			else
				//this type is not represented, the finder will not be called for it normally
				finder[i]=NULL;

		}
	}

};

//---------if due to large inhomogeneity in the sample the the size of the neighbourlist -------
//--------------would exceed the maximum, resize the arrays to be able to continue
void NeighbourList::ResizeArrays(int atom_index, int **ind_list, int **type_list, double **dist2_list, double **v_comp)
{
	int i, itype;
	int offset;
	int old_neigh=max_neigh;
	int safe_add=SAFE_ADD;
	int tot_moved=RunParams::nmoved;
	int *temp_array;

	max_neigh=(int)(max_neigh+ (max_neigh*0.1> safe_add ? max_neigh*0.1 : safe_add));//new value 

	ResizeIArray(old_neigh,max_neigh,usedatom,usedatom,&neigh_list,neighlist_finder);
	ResizeIArray(old_neigh,max_neigh,usedatom,usedatom,&neigh_type,neightype_finder);
	ResizeDArray(old_neigh, max_neigh, usedatom, usedatom, 1, &dsq, dsq_finder);
	ResizeDArray(old_neigh, max_neigh, usedatom, usedatom, 3, &vector_comp, vector_finder);
		
	//reset the values of the pointers belonging to the atom where the error occured
	//find the type (0 to ntypes-1) of the moved atom
	itype=0;
	while (atom_index>=SimpleCfg::cumul[(itype) +1])
			(itype)++;
	offset=max_neigh*(atom_index-SimpleCfg::cumul[itype]);

	*ind_list=neighlist_finder[itype]+offset;
	*type_list=neightype_finder[itype]+offset;
	*dist2_list=dsq_finder[itype]+offset;
	*v_comp=vector_finder[itype]+3*offset;

	//resizing the changed_* arrays
	SetArraysize(&temp_array,tot_moved*(1+2*old_neigh),"temp_array","NeighbourList::ResizeArrays");
	
	for (i=0;i<tot_moved*(1+2*old_neigh);i++)
		temp_array[i]=changed_indices[i];
	delete [] changed_indices;
	SetArraysize(&changed_indices,tot_moved*(1+2*max_neigh),"changed_indices","NeighbourList::ResizeArrays");
	
	if (copy_changed)
		for (i=0;i<tot_moved*(1+2*old_neigh);i++)
			changed_indices[i]=temp_array[i];

	for (i=0;i<tot_moved*(1+2*old_neigh);i++)
		temp_array[i]=changed_types[i];
	delete [] changed_types;
	SetArraysize(&changed_types,tot_moved*(1+2*max_neigh),"changed_types","NeighbourList::ResizeArrays");
	
	if (copy_changed)
		for (i=0;i<tot_moved*(1+2*old_neigh);i++)
			changed_types[i]=temp_array[i];

	for (i=0;i<tot_moved*(1+2*old_neigh);i++)
		temp_array[i]=changed_nneigh[i];
	delete [] changed_nneigh;
	SetArraysize(&changed_nneigh,tot_moved*(1+2*max_neigh),"changed_nneigh","NeighbourList::ResizeArrays");
	
	if (copy_changed)
		for (i=0;i<tot_moved*(1+2*old_neigh);i++)
			changed_nneigh[i]=temp_array[i];

	ResizeIArray(old_neigh,max_neigh,tot_moved*(1+2*old_neigh),tot_moved*(1+2*max_neigh),&changed_neigh_list,NULL);
	ResizeIArray(old_neigh,max_neigh,tot_moved*(1+2*old_neigh),tot_moved*(1+2*max_neigh),&changed_neigh_type,NULL);
	ResizeDArray(old_neigh,max_neigh,tot_moved*(1+2*old_neigh),tot_moved*(1+2*max_neigh),1,&changed_dsq,NULL);
	ResizeDArray(old_neigh,max_neigh,tot_moved*(1+2*old_neigh),tot_moved*(1+2*max_neigh),3,&changed_vector_comp,NULL);

};


