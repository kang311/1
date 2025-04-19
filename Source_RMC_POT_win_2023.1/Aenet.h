//header Aenet.h
//Last changed 12.12.2022

#include "Topology.h"



//extern function declaration for Aenet 
extern double aenet_Rc_max;
extern double aenet_Rc_min;

extern "C" 
{	 
    void aenet_init(int ntypes, char **atom_types,int *stat);
    void aenet_final(int *stat);
    void aenet_load_potential(int index, char *ANNfilename, int *stat);
    void aenet_atomic_energy(double coo_i[3], int type_i, int n_j, double *coo_j, int *type_j, double *E_i, int *stat);
}


//forward declaration 
class RunParams;
class SimpleCfg;
class NeighbourList;
class Move;
struct ThreadArg;

class Aenet
{
    public:
	Aenet(SimpleCfg &config, RunParams &rundata, NeighbourList &neigh);
    Aenet(Aenet &source, SimpleCfg &config, RunParams &rundata, NeighbourList &neigh);//copy constructor
	
	~Aenet(); //destructor
    static int nAtoms;
    static int ntypes;
    static int max_neigh;
    static int ncells_nlist;//for the neighbour list, use 2*ncells_nlist+1 cells for neighbour search
    static int overlap_nlist;//for the neighbousearch
    static int max_tot_moved;
    static int last_gen_aenet_accepted;//ngenerated for the last accepted ANN pot calculation step
    static int last_gen_accepted;//ngenerated for the last accepted calculation step
    static int nchanged;//number of elemet in the *aenet_changed_neigh_ind array
    static int *aenet_changed_neigh_ind;//atoms for which the energy should be recalculated: old and new neighbours of moved atoms 
    static int nstored;//the number of stored atoms in the *aenet_stored_ind array
    static int *aenet_stored_ind;//if E is not calculated in every step, store the indices of the moved atoms and their neighbours for the accepted moves
            //needed to recalc old energy in the next E calc step

    static double cutoff;//cutoff to pass to Aenet library
    static double red_cutoff_sq;//reduced cutoff squared

    static double *E_moved;//temporary storage for the new enery of the moved atoms
    
    static char (*aenet_type_names)[3];//name of the Aenet types, not necessarily the chem_symbols
    //static char (**aenet_type_names_finder)[2];//finder to pass to aenet_init
    static char **aenet_type_names_finder;//finder to pass to aenet_init
    static char(*ANN_filenames)[FILE_NAME_SIZE];//ANN potential file names for each atom type
    
    static bool write_energy;//whether to save the atomic energies
    static bool is_def_cutoff;//whether the original cutoff of the potential is used, needed for output
    static int calc_ANN;//indicator, that ANN should be calculated 
    
   
    double E_tot;//total energy
    double *E_i=nullptr;//energy/atom
    


    //Handles to objects
	SimpleCfg &conf;
	RunParams &rundat;
    NeighbourList &neighlist;

    static void SetParams();//Setting the static members
    static void AenetInit();//initialize aenet, load the potentials
    static void SetCutoff();//set the cutoff and the reduced square
    static void SetNeighList();//create the array for the neigh list
    static void UpdateStored();//update the list of stored indices to calculate the old energy for in the next ANN pot calc step, if not not calculated in every step 
    void Copy(Aenet &target, int *moved_indices, ThreadArg &thread_arg);//copy the energy
    void ResizeIArray(int max, int new_max, int oldsize, int newsize, int **array);
    void CalcANN(ThreadArg &thread_arg);//calculate the initial potential
    void Update(Move &move, int mode, ThreadArg &thread_arg);//determine the atoms to calculate for, and calculate the energy
    void SaveEnergy(const char *file_name="");//save to a file
    int LoadPotBinary();//loading *.pot
    int CheckLoadedPot(ThreadArg &thread_arg);//recalculate the pot for the first atom of each atom type and check against the loaded 
    void SavePotBinary() const;//saving the ptential related parameters to the same *.pot file, which is used by FNC_POT for the other potential
    
    //void my_init(char **atom_types, int stat);
};


inline Aenet::~Aenet()
{
   
    if (E_i!=nullptr)
        delete[] E_i;
    if(::debug) cout<<"Aenet destructor"<<endl;
}