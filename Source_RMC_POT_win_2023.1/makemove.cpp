//source MakeMove.cpp
//Last changed 06.02.2022

//this file defines the Move which is to be applied at each passage in
//the main RMC loop. It can be adapted for special needs

#define _DEF_FILES //not redefine the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include"Move.h"

void Move::MakeMove(SimpleCfg &config, longint &seed, FNC_POT *fnc)//if potential is used, FNC is needed
{
	//this routine defines the move for the Move class;
	//e.g. it defines which atoms are moved and how.
	//It will store their old and new positions in the
	// corresponding Move arrays
	//and update the array of modified partials (*modpart)
	//accordingly.
	//This routine is to be adapted for different needs
	//-----this is the 'standard' RMCA-type move--------
	//arguments are:
	//config: the configuration
	//seed: the seed for the random number generator

	//Swap_fraction controls the fraction of calls to the MakeMove function, where atoms are not just simply moved, but 
	//the last atom of the nmoved moved atoms (index nmoved-1) is swapped with a chosen atom of a different type, as the (nmoved-1)-th. 
	//Swaps are only allowed between type pairs specified in the *.dat file. In case of
	//swap there are nmoved+1 atom changing position! Swap cannot be applied with FNC, (and because of that in case of molecular
	//move), in case of monoatomic system (as here swapping atoms would not change anything), but it can be applied for 
	//multi-atomic moves, and together with the moveout option.


	//In case of moveout option, if the atom was chosen on the basis of its non complience with the cutoffs 
	//keeps the index of the moved tooclose atom in the Move.tclindices array (0 to ntooclose-1)
	//in the itooclose array (mind you: this is NOT the SimpleCfg index!)
	//There can be MORE THAN one tooclose atom among the nmoved moved atom in case of moveout option! 
	
	//Care is taken to ensure that different atoms are chosen in case of atomic RMC with nmoved>1

	//local variables
	int *ptype,*pind;
	int i,m,icoord;
	int found,new_virt;
	int max_cycle;
	int newatom=0; //To ensure, that in case of nmoved>1 the same atom cannot be chosen twice
	double *pold,*pnew, *ppos1,*ppos2,dmove;
	
	int nat,offset,offset2,*piat_ch=0, first_ind_group,j,*pcharge_c_indices=0,jmoved;
	double *pold_charge_c_pos=0,pos,*pold_ch_pos=0;

#ifdef _NO_PERIODIC
//	int	make_new;
	int *p_atom_binind_old;
	double dist,factor;
#endif
	//there can be only one swap regardless the number of moved atoms
	swap = 0;//default: no swap
	tot_moved_atoms = nmoved;
	max_cycle = nmoved;//by default
	tot_moved_virtuals = 0;

	ntcl_selected = 0;//sets the the default, no too close atom is chosen

	for (i = 0; i < tot_moved_atoms; i++)
	{
		tooclose_remove[i] = 1;//the moved atom can be removed from the "too close" list by default,
		tooclose_index[i] = -1;//sets the default for the index of the moved atom among all the atoms of a configuration
						  //if the move is acceptable
	}

	if (makeE0shift)//no atomic move only E0 shift will be in this step
	{
		for (i = 0; i < ExptsData::nek; i++)
		{
			double tr = Ran1(seed);
			ExptsData::ek_gridind[i] = int(ExptsData::ek_ngrid[i] * tr);
		}
	}

	if (makemucorr)//no atomic move, only dmu change, between mu-dmumax->mu+dmumax
	{
		//logfile << "mu generated " << endl;
		for (i = 0; i < ExptsData::nfq; i++)
		{
			if (ExptsData::fqIQbackgcorr[i])
				ExptsData::fqmuact[i] = ExptsData::fqmu[i] + 2 * ExptsData::fqdmumax[i] * Ran1(seed) - ExptsData::fqdmumax[i];
		}
	}
	if (makefprimeshift)//no atomic move, only f' change, between f'-fprimefactor*f' -> f'+fprimefactor*f' 
	{
		for (i = 0; i < ExptsData::nfq; i++)
		{
			if (ExptsData::fqAXS[i]>0)
				ExptsData::fqfprimeact[i * ntypes + ExptsData::fqAXS[i] - 1] = ExptsData::fqfprime[i * ntypes + ExptsData::fqAXS[i] - 1] *
					(1+ (2 * ExptsData::fqfprimefactor[i] * Ran1(seed) - ExptsData::fqfprimefactor[i]));
		}
	}
	if (makeE0shift || makemucorr || makefprimeshift)
	{
		nomove = true;//no atomic move
		return;
	}
	nomove = false;
	
	ptype=types;//pointer to the types array of the Move
	pind=indices;//pointer to the indices of the moved atoms
	pold=oldpos;//pointer to the old positions of the moved atoms
	pnew=newpos;//pointer to the new positions of the moved atoms
	if(fnc)
	{
	  pcharge_c_indices=charge_c_indices;//pointer to the charge centre indices of the moved atoms
	  pold_charge_c_pos=old_charge_c_pos;//pointer to the old positions of the charge group centres of the moved atoms
	}

#ifdef _NO_PERIODIC
	p_atom_binind_old=atom_binind_old;
#endif	
	

	if (Ran1(seed)<swap_fraction)
	{
		//The (nmoved-1)-th atom will be swapped with an extra atom
		swap=1;
		tot_moved_atoms++;
		max_cycle--;//the atoms involved in the swap will be handled separately
	}

	for (m=0;m<max_cycle;m++)//for all atoms to be moved in a single Move, except the swap pair, if there is one
	{	
		do		//To ensure, that in case of nmoved>1 the same atom cannot be chosen twice
		{
			if (newatom==1 && tcl_ind[m]>-1)
			{
				//a new atom has to be selected for the m-th moved atom, as it was already chosen
				ntcl_selected--;//set back the number of selected tooclose atoms
			}
			newatom=0;
			//choose the atom index in the configuration
			if(moveout==1)//if the moveout option is applied
			{
				if(Ran1(seed)<RunParams::too_close_fraction)//choose preferably the atoms too close
				{
					itooclose[ntcl_selected]=int(ntooclose*Ran1(seed));
					tooclose_index[ntcl_selected]=tclindices[itooclose[ntcl_selected]];//sets the index for the index of the moved atom among all the atoms of a configuration
					*pind=*(tclindices+itooclose[ntcl_selected]);
					ntcl_pair[ntcl_selected]=0;//reset the number of too close pairs moved above cutoff for the moved atom
					tcl_ind[m]=ntcl_selected;//indicator showing, that this moved atom is the ntcl_selected-th tooclose moved atom
					ntcl_selected++;//increase the number of too close atoms among the moved
				}
				else
				{
					*pind=int(config.ntotal*Ran1(seed));//(index between 0 and ntotal-1), can be too close!!!
					tcl_ind[m]=-1;//indicator showing, that this moved atom is not a too close by default
					
					//has to check, if it is a too close atom
					for (i=0;i<ntooclose;i++)
					{
						if (*pind==tclindices[i])//this is a tooclose
						{
							itooclose[ntcl_selected]=i;
							tooclose_index[ntcl_selected]=tclindices[itooclose[ntcl_selected]];//sets the index for the index of the moved atom among all the atoms of a configuration
							ntcl_pair[ntcl_selected]=0;//reset the number of too close pairs moved above cutoff for the moved atom
							tcl_ind[m]=ntcl_selected;//indicator showing, that this moved atom is the ntcl_selected-th tooclose moved atom
							ntcl_selected++;//increase the number of too close atoms among the moved
						}
					}
				}

			}//end if moveout option is applied
			else 
			{
				*pind=int(config.ntotal*Ran1(seed));//(index between 0 and ntotal-1)
				tcl_ind[m]=-1;//indicator showing, that this moved atom is not a too close
			}

			for (i=0;i<m;i++)	//To check, that this atom was not chosen before
			{
				if (*(indices+i)==*pind) 
				{
					newatom=1;
					break;
				}
			}
		}
		while (newatom==1);
		ppos1=config.positions+3*(*pind);//pointer to the old 
									//position of the moved atom
								
		//find the type (0 to ntypes-1) of the moved atom
		*ptype=0;
		while (*pind>=config.cumul[(*ptype) +1])
			(*ptype)++;
		//*ptype now points to the type of the moved atom
		//collect the old position of the atom and
		//define the new position of the moved atom
		//\nWARNING! In this version the NEW POSITIONS of the moved atoms are copied to the config.positions
		//array to make the creation of the neighbourlist quicker. If the move is rejected, the old positions are copied
		//back from the oldpos array!!!


		//first the concept was to keep the atoms inside the spherical sample by generating moves only moving the atom inside the sample
		//This meant, that atoms close to the surface had less chance making valid moves
		//the remnants of this approach was kept in the code
		//Therefore this concept was changed the way, that if the atom moved out of the sample, then it would be put back at the other side 
		//of the sample along the line going through the middle of the simulation box and the new position of the atom, and its 
		//distance from the surface will be the same as the distance it went out of the sample
//#ifdef _NO_PERIODIC
		//		make_new=0;
//		do
//		{
//#endif
			for (icoord=0;icoord<3;icoord++)
			{
/*#ifdef _NO_PERIODIC
				if (!make_new)//do not collect old positions from the positions array, if it is not the first move for this atom
					//as the array was already modified to contain the new positions, and the old positions are in the oldpos array already
				{
#endif*/
					*pold=*ppos1;//collect old position
//#ifdef _NO_PERIODIC

//				}
//#endif
				dmove=(2*Ran1(seed)-1)*(*(unscamp+*ptype));//compute new position
				*pnew=*pold +dmove;//add the move to the total move of this atom
				
#ifndef _NO_PERIODIC
				//check for the toroidal box structure
				if(*pnew>=1)
					*pnew-=2;
				else
					if(*pnew<-1)
						*pnew+=2;
#endif
				*ppos1=*pnew;//copy the new positions to the config positions array
				ppos1++;
				pold++;
				pnew++;
			}//next direction

#ifdef _NO_PERIODIC
//			make_new=0;
			//have to make sure, that the atoms stay inside the sample
			dist=sqrt(pow(*(pnew-3),2)+pow(*(pnew-2),2)+pow(*(pnew-1),2));
			if (dist-R0_red>0)
			{
/*				//this move cannot be excepted, as it would move the atom out of the sample
				make_new=1;
				pnew-=3;//set the pointer back
				pold-=3;
				ppos1-=3;*/
				// put it back in the other side of the box
				factor=(dist-2*R0_red)/dist;
				for (icoord=-3;icoord<0;icoord++)
				{
					*(pnew+icoord)=*(pnew+icoord)*factor;
					*(ppos1+icoord)=*(pnew+icoord);
				}

			}
//		}
//		while (make_new);
		//store the old bin index for the imoved moved atom
		*p_atom_binind_old++=HistoSet::atom_binind[*pind];
#endif
		//if there are virtual sites, update their position, if necessary
		new_virt = 0;
		if (SimpleCfg::ntotvirtual > 0)
		{
			
			for (j = 0; j < SimpleCfg::ndiff_virtual[*pind]; j++)
			{
				found = -1;
				for (i = 0; i < tot_moved_virtuals; i++)//check, whether the same virtual site was not added already to the list
				{
					if (SimpleCfg::virtual_ind[*pind * SimpleCfg::max_diff_virtual + j] == indices[tot_moved_atoms + i])
					{
						found = tot_moved_atoms + i;//index of the virtual already in the arrays
						break;
					}
				}
				if (found == -1)//not added
				{
					new_virt++;//new virtual site is added for this moved atom
					if (tot_moved_atoms + tot_moved_virtuals >= tot_moved_atoms * (SimpleCfg::max_diff_virtual + 1))
					{
						cout << "\n*****ERROR*****" << endl;
						cout << "More virtula sites should have been used in the move than the dimension of Move::indices in makemove!" << endl;
						cout << "Exiting..." << endl;
						CleanExit();
					}
					indices[tot_moved_atoms + tot_moved_virtuals] = SimpleCfg::virtual_ind[*pind * SimpleCfg::max_diff_virtual + j];
					if (fnc && update_charge_c)
						charge_c_indices[tot_moved_atoms + tot_moved_virtuals] = fnc->charge_centre[indices[tot_moved_atoms + tot_moved_virtuals]];//index of the charge group centre of the moved atom
					types[tot_moved_atoms + tot_moved_virtuals] = 0;
					while (indices[tot_moved_atoms + tot_moved_virtuals] >= config.cumul[(types[tot_moved_atoms + tot_moved_virtuals]) + 1])
						types[tot_moved_atoms + tot_moved_virtuals]++;//now points to the type of the moved virtual site
					for (icoord = 0; icoord < 3; icoord++)
						oldpos[(tot_moved_atoms + tot_moved_virtuals)*3 + icoord] = config.positions[indices[tot_moved_atoms + tot_moved_virtuals] * 3 + icoord];
					config.CalcVirtualCoord(indices[tot_moved_atoms + tot_moved_virtuals]);//calculating the new coordinates and putting them into the config.positions arra
					for (icoord = 0; icoord < 3; icoord++)
						newpos[(tot_moved_atoms + tot_moved_virtuals)*3 + icoord] = config.positions[indices[tot_moved_atoms + tot_moved_virtuals] * 3 + icoord];
					tot_moved_virtuals++;
				}
				else
					//only the new coordinates has to be recalculated
				{
					config.CalcVirtualCoord(indices[found]);//calculating the new coordinates and putting them into the config.positions arra
					for (icoord = 0; icoord < 3; icoord++)
						newpos[found*3 + icoord] = config.positions[indices[found] * 3 + icoord];
				}
				
			}
		}
		
		//Updating the charge centre positions, if necessary
		if (fnc && update_charge_c)
		{
			*pcharge_c_indices=fnc->charge_centre[*pind];//index of the charge group centre of the moved atom
			offset=3* *pcharge_c_indices;//beginning of the coordinates of the charge group centre
			pold_ch_pos=config.charge_gr_centre+offset;//the beginnig of the charge group centre coordinate 
			//check, whether this charge group centre was not already modified for a previously moved atom
			for (jmoved=0;jmoved<m;jmoved++)
			{
				if (*pcharge_c_indices==charge_c_indices[jmoved])//this charge centre was already modified
				{
					pold_ch_pos=old_charge_c_pos+3 * jmoved;//the original coordinates of this charge centre has to be used
					break;
				}
			}
			
			for (icoord=0;icoord<3;icoord++)
			{
				pold_charge_c_pos[icoord] = *pold_ch_pos++;//copying the original old coordinates of the moved atom's charge centre
				config.charge_gr_centre[offset+icoord]=0;

			}
			for (j = 0; j < new_virt; j++)//fill the old_charge_c_pos array for the newly added virtual sites (they belong to the same cahrge group as the moved atom by definition)
			{
				for (icoord = 0; icoord < 3; icoord++)
					old_charge_c_pos[(tot_moved_atoms+tot_moved_virtuals-new_virt+j)*3+icoord]=pold_charge_c_pos[icoord];

			}
			pold_charge_c_pos += 3;
			piat_ch=fnc->atom_for_charge_gr_finder[fnc->charge_centre[*pind]];//pointer to the first atom of this charge group
			first_ind_group=*piat_ch;//index of the first atom of this charge group
			offset2 = 0;
			//calculating the new centre for the charge group of the moved atom	
			nat = Topology::natoms_per_charge_group[fnc->charge_group[*pind]];
			piat_ch++;//set the pointer to the second atom of this charge group, if there is any
			
			for (j=1+offset2;j<nat;j++)//going through the atoms of this RMC charge group
			{
				//calculating the centre of the charge group, making, sure, that the atoms are inside minimum image range
				for (icoord=0;icoord<3;icoord++)//sum
				{
					pos=config.positions[3* *piat_ch+icoord]-config.positions[3*first_ind_group+icoord];
					if (pos<-1)
						pos+=2;
					else
					{
						if (pos>1)
							pos-=2;
					}
				
					config.charge_gr_centre[offset+icoord]+=pos;
				}
				piat_ch++;
				
			}
		
			for (icoord=0;icoord<3;icoord++)
			{
				config.charge_gr_centre[offset+icoord]/=nat;//average
				config.charge_gr_centre[offset+icoord]+=config.positions[3*first_ind_group+icoord];//translate back
				if (config.charge_gr_centre[offset+icoord]<-1)
					config.charge_gr_centre[offset+icoord]+=2;
				else
				{
					if (config.charge_gr_centre[offset+icoord]>1)
						config.charge_gr_centre[offset+icoord]-=2;
				}
				
				
			}
			pcharge_c_indices++;
		}//end of if update charge centre				
		
		ptype++;
		pind++;
	}//next moved atom of the Move
	
	//to determine the swap pair of the last moved atom, in case of swaps
	if (swap)
	{
		//this was notset earlier
		tooclose_remove[tot_moved_atoms-1] = 1;//the moved atom can be removed from the "too close" list by default,
		tooclose_index[tot_moved_atoms-1] = -1;//sets the default for the index of the moved atom among all the atoms of a configuration
		//choose the type of the swap pair
		m=int(nswap_pairs*Ran1(seed));//index of the swap pair type in the swap_typeX array, chosen from the allowed swap pair types
		
	
		*ptype=swap_type1[m];//type of the first atom of the swap pair
		*(ptype+1)=swap_type2[m];//type of the second atom of the swap pair

		for (m=max_cycle;m<max_cycle+2;m++) //choose the two atoms
		{
			do//to make sure, that this swap atom was not yet selected
			{
				newatom=0;
	
				//choose the index of the first atom of the swap pair
				*pind=config.cumul[*ptype]+int(config.pnatoms[*ptype]*Ran1(seed));//(index between config.cumul[*ptype] and config.cumul[*ptype+1]-1)
				tcl_ind[m]=-1;//indicator showing, that this moved atom is not a too close by default
				
				for (i=0;i<nmoved-1;i++)	//To check, that this atom was not chosen before
				{
					if (*(indices+i)==*pind) 
					{
						newatom=1;//choose another
						break;
					}
				}
				if (!newatom && moveout)//check, whether it is a tooclose atom in case of moveout
				{
					for (i=0;i<ntooclose;i++)
					{
						if (*pind==tclindices[i])//this is a tooclose
						{
							itooclose[ntcl_selected]=i;
							tooclose_index[ntcl_selected]=tclindices[itooclose[ntcl_selected]];//sets the index for the index of the moved atom among all the atoms of a configuration
							ntcl_pair[ntcl_selected]=0;//reset the number of too close pairs moved above cutoff for the moved atom
							tcl_ind[m]=ntcl_selected;//indicator showing, that this moved atom is the ntcl_selected-th too close
							ntcl_selected++;//increase the number of seletcted tooclose atoms among the moved
						}
					}
				}//end of if not newatom && moveout
			}
			while (newatom);
			ptype++;
			pind++;
		};//end of m cycle

		//Swap the coordinates
		ppos1=config.positions+3*(indices[nmoved-1]);//pointer to the old position of the first swap atom
		ppos2=config.positions+3*(indices[nmoved]);//pointer to the old position of the second swap atom
				//pold,pnew points to the oldpos,newpos of the (nmoved-1)-th atom

		for(icoord=0;icoord<3;icoord++)
		{
			
			*pold=*ppos1;//collect old position of the (nmoved-1)-th atom
			*pnew=*ppos2;//set new position of the (nmoved-1)-th atom
			*(pold+3)=*ppos2;//collect old position of the swap pair
			*(pnew+3)=*ppos1;//set new position of the swap pair
			
			//swap the coordinates in the config.positions array to contain the newpositions
			*ppos1=*ppos2;
			*ppos2=*pold;
			ppos1++;
			ppos2++;
			pold++;
			pnew++;
		}
#ifdef _NO_PERIODIC
		//store the old bin index for the swapped atom
		*p_atom_binind_old++=HistoSet::atom_binind[indices[nmoved-1]];
		*p_atom_binind_old=HistoSet::atom_binind[indices[nmoved]];
		//swap the bin indices of the atoms
		HistoSet::atom_binind[indices[nmoved]]=HistoSet::atom_binind[indices[nmoved-1]];
		HistoSet::atom_binind[indices[nmoved-1]]=*p_atom_binind_old;
		
#endif
	}//if swap
	
	return;
	
}//end of the Move::MakeMove function 	
