//header NeighbourList.h
//Last changed 16.02.2022

#include "classes1.h"

//gridding the simulation box by divinding it to ngridcell cell in each dimension,
//determining, which atom belongs to which cell, calculating the neighbourlist
class NeighbourList
{
	public:
	NeighbourList(SimpleCfg &conf, FNC_POT* &fncneigh, RunParams &rundata);//constructor

	~NeighbourList();//destructor

	static int ngridcells;//number of gird cells in one dimension, 1: no gridding is applied
	static int ngridcells2;//square of number of gird cells
	static int max_gridatom;//maximum number of atoms in a grid cell
	static int natoms;//the number of atoms in the cell
	static int ntypes;//the number of atom type (NECESSARY for memory allocation)
	static int nconsttype;//number of constrained types
	static int ncells_cutoff;//number of grid cells to test for each direction for cutoff check
	static int overlap;//number of gridcells to leave out from checking the cutoff in each positiv direction due to overlapping
	static int overlap_nlist;//number of gridcells to leave out from calculation of the neighbourlist in each positiv direction due to overlapping
	static int ncells_nlist;//number of grid cells to test for each direction for neighbourlist
	static int max_neigh;//maximum number of neighbopurs for the neighbourlist for an atom
	static int copy_changed;//indicator, whether in case of resizing the arrays connected to max_neigh the changed_* arrays should be copied
	static int usedatom;//number of atoms involved in neighbourlist calculation
	static longint consttype;//each bit represents an atom type, containing a 1 for each type that is involved in a constraint, which require the neighbour list
	static double xmax;//the maximum distance for the neighbourlist in reduced unit
	static double xmax2;//square of the maximum distance for the neighbourlist in reduced unit
	static double max_red_cutoff;//maximum cutoff distance in reduced unit
	static double cellwidth;//the width of the grid cell in reduced units
	static double half_ngridcells;//the half of ngridcells, to increase speed

	//arrays for the gridding
	int *grid_index;//array of grid index in each direction for each atoms
	int *atom_indices;//array of indices of atoms for each grid
	int *grid_count;//number of atoms in a grid cell
	int *old_grid;//indices of the old grid cell (the cell before the move) for each moved atom

	//arrays connected to the neighbourlist
	int *nneigh;//number of neighbours for each central atom represented by the neigh_list
	int *neigh_list;//indices of the neighbours
	int *neigh_type;//types of the neighbours
	double *dsq;//squared distance for the neighbours
	double *vector_comp;//components of the Central-Neighbour vector
	
	//arrays connected to storing the changed array parts of the neighbour lists, while the acceptance
	//of the move is decided
	//not the whole lists are copied, as that would be time and memory consuming
	int nchanged;//number of atoms effected by the move
	int *changed_indices;//array of atoms (moved, and their new and old neighbours), which are effected
						//by the move
	int *changed_types;//array of the type of the atoms (moved, and their new and old neighbours), which are effected
						//by the move
	int *changed_nneigh;//number of neighbours for each central atom represented by the changed_indices list
	int *changed_neigh_list;//indices of the neighbours
	int *changed_neigh_type;//types of the neighbours
	double *changed_dsq;//squared distance for the neighbours
	double *changed_vector_comp;//components of the Central-Neighbour vector


	int **nneigh_finder;//finder of the beginning of the neighbourlist of this central type/type
	int **neighlist_finder;//finder of the beginning of the neighbourlist of this central type/type
	int **neightype_finder;//finder of the beginning of the neighbour types of this central type/type
	double **dsq_finder;//finder of the beginning of the squared distane list of this central type/type
	double **vector_finder;//finder of the beginning of the vector components of this central type/type
	

	//handles to objects
	SimpleCfg &config;
	FNC_POT* &fnc;
	RunParams &rundat;
#ifdef _AENET
	static Threads *thread;
#endif

	static void SetParams();//initializing the static members
	static void SetCutoff();//set the cutoff related params, had to be removed from SetParams, to call it after hist calc in case of automatic cutoff determination
	void MakeGrid();//Creating the grid, determining the gridcell for each atom
	void UpdateGrid(Move &moved);//updating the gridding for the moved atom(s) 
	void ResetGrid(int *moved_indices);//resetting the gridding for the moved atom(s) in case of the move was rejected 
	int &GridCount(int i,int j,int k);//setting or reding grid_count, the number of atoms in the (i,j,k) grid cell
	int *PtrGridCount(int i,int j,int k);//pointer to the number of atoms in the (i,j,k) grid cell
	//the index of the atom, this is the atom_ind_grid-th atom for this grid cell
	int &AtomIndex(int i,int j,int k, int atom_ind_grid);
	//pointer to the index of the atom, this is the atom_ind_grid-th atom for this grid cell
	int *PtrAtomIndex(int i,int j,int k, int atom_ind_grid);
	void SaveGrid(const char *file_name="");//saving the object to a file
	int LoadGrid();//loading from file
	void ResizeIArray(int max, int new_max, int oldsize, int newsize, int **array, int **finder);//if due to large inhomogeneity in the sample the atom count of one gridcell or the size of the neighbourlist 
				//would exceed the maximum, resize integer array to be able to continue
	void ResizeDArray(int max, int new_max, int oldsize, int newsize, int multiply, double **array, double **finder);//if due to large inhomogeneity in the sample the atom count of one gridcell or the size of the neighbourlist 
				//would exceed the maximum, resize double array to be able to continue
	void ResizeArrays(int atom_index, int **ind_list, int **type_list, double **dist_list, double **v_comp);
	
	void InitNeighList();//calculating the neighbourlist for each atom belonging to the types involved in CosDistrConst
	//determining the neighbours of atom atom_index inside xmax
	int CalcNeighList(int atom_index, int *ind_list, int *type_list, double *dist2_list, double *v_comp, double *dist_list=nullptr, int *my_max_neigh=&NeighbourList::max_neigh, int id=0, double *coord_list=nullptr, Move *move=nullptr, int imoved=-1);
	void UpdateNeighList(Move &move);//updating the lists to reflect the move, store the old values for the effected atoms
	void StoreList(Move &move);//store the values of the neighbourlist connected arrays for the neighbours of the moved atom(s) 	
	void ResetList();//if the move was rejected, reset the neighbour lists using the stored old values
	void SaveNeilist(const char *file_name="");//saving the object to a file

	int CheckCutoff(Move &move);
				

};

inline NeighbourList::~NeighbourList()
{
	delete [] grid_index;
	delete [] atom_indices;
	delete [] grid_count;
	delete [] old_grid;
	delete [] nneigh;
	delete [] neigh_list;
	delete [] dsq;
	delete [] vector_comp;
	delete [] nneigh_finder;
	delete [] neighlist_finder;
	delete [] dsq_finder;
	delete [] vector_finder;
	if(::debug) cout<<"NeighbourList destructor"<<endl;
};
