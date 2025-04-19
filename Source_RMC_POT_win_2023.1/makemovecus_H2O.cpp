//source makemove.cpp for the H20 molecule
//Last changed 06.03.2021

//this file defines the Move which is to be applied at each passage in
//the main RMC loop.
//This Move is designed for the H2O molecules
//i.e. nmoved=3, no preferred rotation axis
//move= rotation+translation+shaking of H atoms
//WARNING! In this version the NEW POSITIONS of the moved atoms are copied to the config.positions
//array to make the creation of the neighbour list quicker. If the move is rejected, the old positions are copied
//back from the oldpos array!!!

#define _DEF_FILES //not redifen the file names included through files.h
#define _DEF_INTERACTION_FUNC//not to redefine the pointer to the intercation functions
#include"Move.h"


//================================================================
//
//--------custom molecular Move constructor-----------------------
//
//================================================================


//overloaded constructor for custom 'molecular' move
Move::Move(FNC_POT &fnc, int &ntoocl, int *tclarray, SimpleCfg &conf)
	:config(conf)
{
	//this class defines the Move applied at each passage in the
	//main RMC loop 
	//some of the moves parameters appear in the .dat file
	//(e.g. move amplitude for each atom type)
	//if a non-standard (i.e non-RMCA move) is to be defined
	//then the possible additional parameters can be added
	//here 
	
	//arguments are:
	//1) a .mov file which provide the necessary parameters for the
	//		custom move
	//2) FNC constarint
	//3) ntoocl: the number of atoms not satisfying the cutoffs calculated by HistCalc (loaded too 
	//close atoms not included, they are represented by ntooclose
	//(this is initially needed for memory allocation for the tclindices array
	//but it is later updated as these atoms are moved and their number decreases)
	
	//4) tclarray is an array for all atoms containing how many too close pairs they are inviolved calculated by HistCalc
	//5) for initializing the handle to the SimpleCfg object
	//check if the .mov file is open
	
	//check if the .mov file is open
	ifstream movefile;
	if(RunParams::cus_entry.size()==0)
	{
	  cleanopen(movefile,movefilename);
	  if (CheckFileState(movefile,"Move::Move custom constructor",movefilename)==0)
	  {
		  cout << "\n*****ERROR*****" << endl;
		  cout<<"Move:: custom constructor does not find the "<<movefilename<<" file"<<endl;
		  cout<<"Cannot run this way, exiting..."<<endl;
		  cleanexit();
	  }
	}

	if (nmoved!=3)
	{
		cout << "WARNING(" << ++warn << "): Number of moved atoms is "<<nmoved<<" from the .dat file, but it has to be 3 for H2O molecular move!"<<endl;
		cout<<"\tNumber of moved atoms is reset to 3!"<<endl;
		nmoved=3;//reset it for H2O
	}
	
	
	//constructor
	int i,k;
	//double temp;

	max_pairs=0;
	//Setting the parameters for the moveout option	
	if (ntooclose+ntoocl==0 && moveout==1)//if all atoms satisfy the cutoffs
	{
		cout << "NOTE(" << ++note << "): All atoms satisfy the cutoffs. Moveout option switched off. "<<endl;
		moveout=0;	//switch off the moveout option
		RunParams::moveout=0;// in the run object as well
	}

	//creation of the tclindices and tclpaircount arrays, if the atoms below cutoff were not loaded
	if (tooclose_flag==0 && moveout)
	{
		//the too close atoms (if there is any) were not loaded, initialize the connected parameters
		ntooclose=ntoocl;//number of atoms not satisfying the cutoffs
		
		if (ntooclose!=0)
		{
			tclindices=new int[ntooclose];//array of indices of 'bad' atoms
			if (tclindices==NULL)
				no_array("Move::Move","tclindices");
			tclpaircount=new int[ntooclose];//array of too close pairs the given atom is involved
			if (tclpaircount==NULL)
				no_array("Move::Move","tclpaircount");
			tcl_ind=new int[nmoved];//array, the serial number of this moved too close atoms among the moved too 
									//close atoms beginning with 0, -1 if the moved atom is not too close
			if (tcl_ind==NULL)
				no_array("Move::Move","tcl_ind");
		}
	

		i=0;//there were no loaded atoms below cutoff
		//Finding the atoms not satisfying the cutoffs to build the tclindices and tclpaircount arrays
		for (k=0;k<SimpleCfg::ntotal;k++)//for all atoms
		{
			if (tclarray[k]>0)
			{ 
				tclindices[i]=k;//put its index in the tclindices array
				tclpaircount[i]=tclarray[k];//put its atom pairs count in the tclpaircount array
				i++;
			}
		}	
	
	}
	if (moveout)
	{
		for (i=0;i<ntooclose;i++)
			if (tclpaircount[i]>max_pairs)
				max_pairs=tclpaircount[i];//determining the maximum number of the pairs for an atom
	}
	tclpair_ind=new int[max_pairs*tot_moved];//creating the array for the indices of the removable pairs of the moved atom
	if (tclpair_ind==NULL)
		no_array("Move::Move","tclpair_ind");

	//read the .mov file
	if(RunParams::cus_entry.size()==0)
	{
	  skipline(movefile, movefilename,1);
	  nparams=readthisline(movefile,1,1,"nparams","Move::Move",movefilename);//move parameters
	}
	else nparams=RunParams::cus_entry.size()-std::count_if(RunParams::cus_entry.begin(),RunParams::cus_entry.end(),[](std::string s){return s.find("NMOVED-ATOMS")!=std::string::npos;});//Get the number of real parameters
	//array of custom move parameters
	custparams=new double[nparams];
	if (custparams==NULL)
		no_array("Move::Move","custparams");

	//arrays of move amplitudes - 1 per atom type

	unscamp=new double[nparams];//array of reduced move amplitudes
	if (unscamp==NULL)
		no_array("Move::Move","unscamp");
	
	itooclose=new int[nmoved];//indices of the moved atom(s) among the too close atoms, if it was a "too close" one, -1 otherwise
	if (itooclose==NULL)
		no_array("Move::Move","itooclose");
	tooclose_index=new int[nmoved];//indices of the moved too close atom(s) in the configuration (starting with 0)
	if (tooclose_index==NULL)
		no_array("Move::Move","tooclose_index");
	tooclose_remove=new int[nmoved];//indicator, whether the moved too close atom(s) can be removed from the list
	if (tooclose_remove==NULL)
		no_array("Move::Move","tooclose_remove");
	ntcl_pair=new int[nmoved];//number of too close pairs the given moved tooclose atom(s) is/are involved in, which can be removed from the list 
	if (ntcl_pair==NULL)
		no_array("Move::Move","ntcl_pair");
	indices=new int[nmoved];//indices of moved atoms in the cfg (start at 0)
	if (indices==NULL)
		no_array("Move::Move","indices");
	types=new int[nmoved];//types of atoms moved (start at 0)
	if (types==NULL)
		no_array("Move::Move","types");
	tcl_ind=new int[nmoved];//array, the serial number of this moved too close atoms among the moved too close atoms beginning 
							//with 0, -1 if the moved atom is not too close
	if (tcl_ind==NULL)
		no_array("Move::Move","tcl_ind");
	oldpos=new double[nmoved*3];//old coordinates of moved atoms
	if (oldpos==NULL)
		no_array("Move::Move","oldpos");
	newpos=new double[nmoved*3];//new coordinates of the moved atoms
	if (newpos==NULL)
		no_array("Move::Move","newpos");
	//the order of the coordinates (in box coordinates system) is
	//{x0,y0,z0,x1,y1,z1,x2,...,x(nmoved-1),y(nmoved-1),z(nmoved-1)}
	//initialisation of the arrays modified at each move
	for(i=0;i<nmoved;i++)
		*(indices+i)=0;

	//array of types of moved atoms (ALWAYS THE SAME)
	//in program format: first type start at 0
	types[0]=0;//Oxygene
	types[1]=1;//H 1
	types[2]=1;//H 2

	for(i=0;i<3*nmoved;i++)
	{
		*(oldpos+i)=0.0;
		*(newpos+i)=0.0;
	}
	
	//initialisation of permanent arrays from the RunParams in argument
	for(i=0;i<ntypes;i++)
		*(unscamp+i)=*(moveamp+i)/RunParams::boxedge;//REDUCED move amp		
		
	//Initializing potential related features
	charge_c_indices=0;
	old_charge_c_pos=0;
	update_charge_c=0;
	
//******************************************************
//this part is to be adapted to the molecule under study
//( for reading the .cus file)
//******************************************************
	if(RunParams::cus_entry.size()==0)
	{
	  custparams[0]=readthisline(movefile,1,1.0,"translational amplitude","Move::Move",movefilename)/RunParams::boxedge;//translation amplitude read in Angstrom than reduced
	  custparams[1]=readthisline(movefile,1,1.0,"rotational amplitude","Move::Move",movefilename)*0.0174533;//rotation amplitude read in degrees, transformed to radians
	  custparams[2]=readthisline(movefile,1,1.0,"move amplitude for hydrogen","Move::Move",movefilename)/RunParams::boxedge;//individual move amplitude for hydrogen read in angstrom than reduced
	  custparams[3]=readthisline(movefile,1,1.0,"safety distance","Move::Move",movefilename)/RunParams::boxedge;//safety distance read in angstrom than reduced
	  movefile.close();
	}
	else
	{
	  int counter=0;
	  for(std::list<std::string>::iterator it=RunParams::cus_entry.begin(); it!=RunParams::cus_entry.end(); it++)
	  {
	    std::string key,params;
	    WrapFree(*it,key,params);
	    
	    std::stringstream ss(params);
	    if(key=="CUSTOM-LINE")
	    {
	      double temp;
	      if(!(ss>>temp)&& !ErrLackValue(*it)) cleanexit();
	      else
	      {
		switch(counter)
		{
		  case 0:							//translation amplitude read in Angstrom than reduced
		  case 2: 							//individual move amplitude for hydrogen read in angstrom than reduced
		  case 3:custparams[counter]=temp/RunParams::boxedge; break;  	//safety distance read in angstrom than reduced
		  case 1: custparams[counter]=temp*0.0174533; break;		//rotation amplitude read in degrees, transformed to radians
		  default:
		    std::cout<<"ERROR: more parameter in makemovecus_H2O.cpp than available slots."<<std::endl;
		    cleanexit();
		}
	      }
	      counter++;
	    }
	  }
	}
};	


//================================================================
//
//           The makemove function
//
//================================================================


void Move::makemove(SimpleCfg &config, FNC_POT &fnc, long int &seed)
{
	//local variables
	int i,j,it;
	//int count;//number of too close atoms in the molecule
	int *pind, *pi3;
	double *pd1, *pd2, *pp1;
	int across;//indicator if the molecule is across the limit of the cell
	double temp;
	double safedist;//safety distance to the cell border beyond which the 
				//molecule is within the box 
	double xc,yc,zc;
	double ampl;//move random amplitude
	double amplx, amply,amplz;//move amplitude along directions
	double theta, phi;//angles defing the raotation axis
	double ctheta,stheta,cphi,sphi, omega, comega, somega;//angular parameters
	double intpos1[9];//intermediate positions array
	double intpos2[9];//intermediate positions array

	ntcl_selected=0;//sets the the default for the number of the too close atoms among the moved atoms 
	for (i=0;i<nmoved;i++)
	{
		tooclose_remove[i]=1;//the moved atom can be removed from the "too close" list by default,
		tooclose_index[i]=-1;//sets the default for the index of the moved atom among all the atoms of a configuration
						  //if the move is acceptable
		tcl_ind[i]=-1;//indicator showing, that this moved atom is not a too close by default
	}

	//defining the safety distance (largest possible O-H distance)
	safedist=*(custparams+3);//last element of the custom parameters array
							//(read from the .cus file by the constructor)

	
	//array of indices of moved atoms (defined by the FNC)
	pind=indices;//pointer to the indices of the moved atoms
	//choose the index of the molecule
	if(moveout==1)//if the moveout option is applied
	{
		//count=0;//default, number of tooclose atoms in the molecule
		if (ran1(seed)<0.1)
		{
			//choose from the too close atoms
			itooclose[0]=int(ntooclose*ran1(seed));//might be reset later, if the tooclose is not a carbon
			tooclose_index[0]=tclindices[itooclose[0]];//might be reset later, if the tooclose is not a carbon
			if (tooclose_index[0]<(config.ntotal/nmoved))//between 0 and the number of molecules 
			{
				//the too close atom is the oxygene  of the molecules with all the fnc neighbours given in the *.fnc file
				*pind=tooclose_index[0];//oxygene of moved molecule
				i=*(fnc.converter+tooclose_index[0]);//offset of the first neighbour in the neighbours array 

			}
			else
			{
				//the too close atom is hydrogene
				i=*(fnc.converter+tooclose_index[0]);//offset of the first neighbour (first H) in the neighbours array
				*pind=fnc.neighbours[i];//index of the first H atom, as it is always the first neighbour
				i=*(fnc.converter+ indices[0]);//offset of the first neighbour in the neighbours array 
			}
		}//chosen through too close atom
		else
		{
			//any molecule is chosen, can be too close atom among the atoms!
			*pind=int(config.ntotal/nmoved*ran1(seed));//[index between 0 and (number of molecules-1) ]
	
			i=*(fnc.converter+*pind);//offset of the first neighbour in the neighbours array
		}//chosen randomly

		//has to set the variables connected with the tooclose atoms, if there is/are too close atom(s) in the molecule
		for (it=0;it<ntooclose;it++)
		{
			if (*pind==tclindices[it])//the O atom is a too close
			{
				itooclose[ntcl_selected]=it;//index of the too close atom in the tcl arrays
				tooclose_index[ntcl_selected]=tclindices[it];//sets the index for the index of the moved atom among all the atoms of a configuration
				ntcl_pair[ntcl_selected]=0;//reset the number of too close pairs moved above cutoff for the moved atom
				tcl_ind[0]=ntcl_selected;//indicator showing, that this moved atom is the ntcl_selected-th tooclose moved atom
				ntcl_selected++;//increase the number of tooclose atoms among the moved
				continue;//go to next too close
			}
			for (j=0;j<nmoved-1;j++)//check the two H-s
			{
		
				if (*(fnc.neighbours+i+j)==tclindices[it])//this is a tooclose
				{
					itooclose[ntcl_selected]=it;//index of the too close atom in the tcl arrays
					tooclose_index[ntcl_selected]=tclindices[it];//sets the index for the index of the moved atom among all the atoms of a configuration
					ntcl_pair[ntcl_selected]=0;//reset the number of too close pairs moved above cutoff for the moved atom
					tcl_ind[j+1]=ntcl_selected;//indicator showing, that this moved atom is the ntcl_selected-th tooclose moved atom
					ntcl_selected++;//increase the number of tooclose atoms among the moved
				}
			}

		}//end of it cycle for tooclose atoms
	}//end of if moveout is applied
	else
	{
		*pind=int(config.ntotal/nmoved*ran1(seed));//[index between 0 and (numb. of molecules)-1 ]
	
		i=*(fnc.converter+*pind);//offset of the first neighbour in the neighbours array
	}


	
	pind++;//move to the second element of the indices array
	pi3=fnc.neighbours+i;//positioning the pointer in the neighbours array
	*pind++=*pi3++;//first neighbour=first HYDROGEN
	*pind=*pi3;//second neighbour=second HYDROGEN

	//collect the old positions of the moved atoms
	pd1=oldpos;//pointer to the old positions array
	pd2=newpos;//pointer to the new positions array
	pind=indices;//pointer to moved atoms indices
	for(i=0;i<nmoved;i++)//for each atom moved
	{
		pp1=config.positions+3*(*pind++);//pointer to the old 
								//position of the moved atom
		for(j=0;j<3;j++)//for each coordinate
		{
			*pd1++=*pp1;//collect old position
			*pd2++=*pp1++;//collect old position
		}
	
	}
	//check for molecule crossing the cell limit
	across=0;
	if (1-*newpos<safedist)
		across=1;
	else
	{
		if (*newpos+1<safedist)
			across=1;
		else
		{
			if (1-*(newpos+1)<safedist)
				across=1;
			else
			{
				if (*(newpos+1)+1<safedist)
					across=1;
				else
				{
					if (1-*(newpos+2)<safedist)
						across=1;
					else
						if (*(newpos+2)+1<safedist)
							across=1;
				}
			}
		}
	}
	//correct positions if necessary
	if (across>0)	
	{
		pd2=newpos+3;//points to the position of the first hydrogen	
		for(i=1; i<nmoved;i++)//for all neighbours of the oxygen
		{
			pd1=newpos;//points to the position of the first carbon
			for(j=0;j<3;j++)//for all coordinates
			{
				temp=*pd2-*pd1;
				if(temp>1) 
					*pd2-=2;
				else
					if(temp<-1) 
						*pd2+=2;
				pd1++;
				pd2++;
			}//next coordinate
		}//next neighbour
	}//end of if across>0
		
	//translations
	//amplitudes
	amplx=(2*ran1(seed)-1)*(*(custparams));
	amply=(2*ran1(seed)-1)*(*(custparams));
	amplz=(2*ran1(seed)-1)*(*custparams);

	//translation along the x-axis
	for(i=0;i<nmoved;i++)
		*(newpos+3*i)+= amplx;
	//translation along the y-axis
	for(i=0;i<nmoved;i++)
		*(newpos+3*i+1)+= amply;
	//translation along the z-axis
	for(i=0;i<nmoved;i++)
		*(newpos+3*i+2)+= amplz;

	//individual atomic movements (2 hydrogen only)
	ampl=*(custparams+2);//amplitude of individual moves
	pd1=newpos+3;//point to the new hydrogen positions
	for(i=0;i<2;i++)// for each hydrogen
		*pd1++ +=(2*ran1(seed)-1)*ampl;

	//compute the  'rotation center' position (oxygen atom)
	pd1=newpos;//first carbon position
	xc=*pd1++;
	yc=*pd1++;
	zc=*pd1;

	//rotation
	//first translate  the molecule to the origin
	pd1=intpos1;//points to intermediate positions array
	pd2=newpos;
	for(i=0;i<nmoved;i++)
	{
		*pd1++=*pd2++ -xc;
		*pd1++=*pd2++ -yc;
		*pd1++=*pd2++ -zc;
	}

	//compute the angular parameters for the rotation
	//no preferential rotation axis
	theta=(2*ran1(seed)-1)*PI;
	phi=2*ran1(seed)*PI;

	ctheta=cos(theta);//cos theta
	stheta=sin(theta);//sin theta
	if(stheta==0)
	{
		cphi=1;
		sphi=0;//subrotation matrix =identity
	}
	else 
	{
		cphi=cos(phi);//cphi
		sphi=sin(phi);//sphi
	}

	//rotation angle
	omega=(2*ran1(seed)-1)*(*(custparams+1));
	comega=cos(omega);
	somega=sin(omega);	
		
	//apply subrotations (TO HYDROGEN ATOMS ONLY)
	for(i=1;i<nmoved;i++)
	{
		pd1=intpos1+3*i;//x coordinate of the i-th atom
		pd2=intpos2+3*i;
		//Rz(-phi)
		*pd2=*pd1*cphi+*(pd1+1)*(sphi);
		*(pd2+1)=*pd1*(-sphi)+*(pd1+1)*cphi;
		*(pd2+2)=*(pd1+2);
		//Ry(-theta)
		*pd1=ctheta*(*pd2)-stheta*(*(pd2+2));
		*(pd1+1)=*(pd2+1);
		*(pd1+2)=stheta*(*pd2)+ctheta*(*(pd2+2));
		//Rz(omega)
		*pd2=*pd1*comega+*(pd1+1)*-somega;
		*(pd2+1)=*pd1*somega+*(pd1+1)*comega;
		*(pd2+2)=*(pd1+2);
		//Ry(theta)
		*pd1=ctheta*(*pd2)+stheta*(*(pd2+2));
		*(pd1+1)=*(pd2+1);
		*(pd1+2)=-stheta*(*pd2)+ctheta*(*(pd2+2));
		//Rz(phi)
		*pd2=*pd1*cphi+*(pd1+1)*(-sphi);
		*(pd2+1)=*pd1*sphi+*(pd1+1)*cphi;
		*(pd2+2)=*(pd1+2);
	}//next atom

	//translate Hydrogen atoms back to original position
	//(the new positions are in intpos2)	
	pd2=intpos2+3;//points to the 1st Hydrogen in intermediate positions array
	pd1=newpos+3;//points to the first H in new positions array
	for(i=1;i<nmoved;i++)
	{
		*pd1++=*pd2++ +xc;
		*pd1++=*pd2++ +yc;
		*pd1++=*pd2++ +zc;
	}

	//apply toroidal cell structure
	pd1=newpos;
	for(i=0;i<9;i++)//3 times 3 atoms
	{
		if(*pd1>1) 
			*pd1-=2;
		else 
			if (*pd1<-1) 
				*pd1+=2;
		pd1++;
	}
	
	//update the config.positions array with the new positions!
	pd1=newpos;
	pind=indices;//pointer to moved atoms indices
	for(i=0;i<nmoved;i++)//for each atom moved
	{
		pp1=config.positions+3*(*pind);//pointer to the old 
								//position of the moved atom
		pind++;
		for(j=0;j<3;j++)//for each coordinate
			*pp1++=*pd1++;//copy new positions

	}
	
}//end of the Move::makemove function 	
