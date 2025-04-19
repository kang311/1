//source makemovecus.cpp for CCl4 molecule
//Last changed 04.03.2021

//this file defines the Move which is to be applied at each passage in
//the main RMC loop.
//This Move is designed for the CCl4 molecules
//i.e. nmoved=5, preferred rotation axis
//move= rotation+translation+shaking 
//\nWARNING! In this version the NEW POSITIONS of the moved atoms are copied to the config.positions
//array to make the creation of the neighbourlist quicker. If the move is rejected, the old positions are copied
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
	//2) FNC constraint
	//3) ntoocl: the number of atoms not satisfying the cutoffs calculated by HistCalc (loaded too 
	//close atoms not included, they are represented by ntooclose
	//(this is initially needed for memory allocation for the tclindices array
	//but it is later updated as these atoms are moved and their number decreases)
	
	//4) tclarray is an array for all atoms containing how many too close pairs they are inviolved calculated by HistCalc
	//5) for initializing the handle to the SimpleCfg object
	//check if the .mov file is open
	ifstream movefile;
	if(RunParams::cus_entry.size()==0)
	{
	  CleanOpen(movefile,movefilename);
	  if (CheckFileState(movefile,"Move::Move custom constructor",movefilename)==0)
	  {
		  cout << "\n*****ERROR*****" << endl;
		  cout<<"Move:: custom constructor does not find the "<<movefilename<<" file"<<endl;
		  cout<<"Cannot run this way, exiting..."<<endl;
		  CleanExit();
	  }
	}
	if (nmoved!=5)
	{
		cout << "\nWARNING(" << ++warn << "): Number of moved atoms is "<<nmoved<<" from the .dat file, but it has to be 5 for CCl4 molecular move!"<<endl;
		cout<<"\tNumber of moved atoms is reset to 5!"<<endl;
		nmoved=5;//reset it for CCl4
	}
	
	//constructor
	int i,k;
	
	max_pairs=0;

	//Setting the parameters for the moveout option	
	if (ntooclose+ntoocl==0 && moveout==1)//if all atoms satisfy the cutoffs
	{
		cout << "\nNOTE(" << ++note << "): All atoms satisfy the cutoffs. Moveout option switched off. "<<endl;
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
			SetArraysize(&tclindices,ntooclose,"tclindices","Move::Move");//array of indices of 'bad' atoms
			
			SetArraysize(&tclpaircount,ntooclose,"tclpaircount","Move::Move");//array of too close pairs the given atom is involved
			SetArraysize(&tcl_ind,nmoved,"tcl_ind","Move::Move");//array, the serial number of this moved too close atoms among the moved too 
									//close atoms beginning with 0, -1 if the moved atom is not too close
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
	SetArraysize(&tclpair_ind,max_pairs*tot_moved_atoms,"tclpair_ind","Move::Move");//creating the array for the indices of the removable pairs of the moved atom
			
	//read the .mov file
	if(RunParams::cus_entry.size()==0)
	{
	  SkipLine(movefile, movefilename,1);
	  nparams=ReadThisLine(movefile,1,1,"nparams","Move::Move",movefilename);//move parameters
	}
	else 
		nparams= (int)(RunParams::cus_entry.size()-std::count_if(RunParams::cus_entry.begin(),RunParams::cus_entry.end(),[](string s){return s.find("NMOVED-ATOMS")!=string::npos;}));//Get the number of real parameters
	//array of custom move parameters
	SetArraysize(&custparams,nparams,"custparams","Move::Move");

	//arrays of move amplitudes - 1 per atom type
	SetArraysize(&unscamp,nparams,"unscamp","Move::Move");//array of reduced move amplitudes
	//parameters updated at each application of the move
	SetArraysize(&itooclose,nmoved,"itooclose","Move::Move");//indices of the moved atom(s) among the too close atoms, if it was a "too close" one, -1 otherwise
	
	SetArraysize(&tooclose_index,nmoved,"tooclose_index","Move::Move");//indices of the moved too close atom(s) in the configuration (starting with 0)
	SetArraysize(&tooclose_remove,nmoved,"tooclose_remove","Move::Move");//indicator, whether the moved too close atom(s) can be removed from the list
	SetArraysize(&ntcl_pair,nmoved,"ntcl_pair","Move::Move");//number of too close pairs the given moved tooclose atom(s) is/are involved in, which can be removed from the list 
	SetArraysize(&indices,nmoved,"indices","Move::Move");//indices of moved atoms in the cfg (start at 0)
	SetArraysize(&types,nmoved,"types","Move::Move");//types of atoms moved (start at 0)
	SetArraysize(&tcl_ind,nmoved,"tcl_ind","Move::Move");//array, the serial number of this moved too close atoms among the moved too close atoms beginning 
							//with 0, -1 if the moved atom is not too close
	SetArraysize(&oldpos,nmoved*3,"oldpos","Move::Move");//old coordinates of moved atoms
	SetArraysize(&newpos,nmoved*3,"newpos","Move::Move");//new coordinates of the moved atoms
	
	//the order of the coordinates (in box coordinates system) is
	//{x0,y0,z0,x1,y1,z1,x2,...,x(nmoved-1),y(nmoved-1),z(nmoved-1)}
	//initialisation of the arrays modified at each move
	for(i=0;i<nmoved;i++)
		*(indices+i)=0;
	
	//array of types of moved atoms (ALWAYS THE SAME)
	//in program format: first type start at 0
	types[0]=0;//carbon 1
	for (i=1;i<5;i++)
		types[i]=1;//chlorine 1-4

	for(i=0;i<3*nmoved;i++)
	{
		*(oldpos+i)=0.0;
		*(newpos+i)=0.0;
	}

	//initialisation of permanent arrays from the RunParams in argument
	for(i=0;i<ntypes;i++)
		*(unscamp+i)=*(moveamp+i)/RunParams::boxedge;//REDUCED move amplitude
		
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
	  custparams[0]=ReadThisLine(movefile,1,1.0,"translational amplitude","Move::Move",movefilename)/RunParams::boxedge;//translation amplitude read in Angstrom than reduced
	  custparams[1]=ReadThisLine(movefile,1,1.0,"rotational amplitude","Move::Move",movefilename)*0.0174533;//rotation amplitude read in degrees, transformed to radians
	  custparams[2]=ReadThisLine(movefile,1,1.0,"move amplitude for carbon","Move::Move",movefilename)/RunParams::boxedge;//individual move amplitude for carbon read in angstrom than reduced
	  custparams[3]=ReadThisLine(movefile,1,1.0,"move amplitude for chlorine","Move::Move",movefilename)/RunParams::boxedge;//individual move amplitude for chlorine read in angstrom than reduced
	  custparams[4]=ReadThisLine(movefile,1,1.0,"safety distance","Move::Move",movefilename)/RunParams::boxedge;//safety distance read in angstrom than reduced
	  movefile.close();
	}
	else
	{
	  int counter=0;
	  for(std::list<string>::iterator it=RunParams::cus_entry.begin(); it!=RunParams::cus_entry.end(); it++)
	  {
	    string key,params;
	    WrapFree(*it,key,params);
	    
	    stringstream ss(params);
	    if(key=="CUSTOM-LINE")
	    {
	      double temp;
	      if(!(ss>>temp)&& !ErrLackValue(*it)) CleanExit();
	      else
	      {
		switch(counter)
		{
		  case 0:							//translation amplitude read in Angstrom than reduced
		  case 2: 							//individual move amplitude for carbon read in angstrom than reduced
		  case 3:							//individual move amplitude for chlorine read in angstrom than reduced
		  case 4:custparams[counter]=temp/RunParams::boxedge; break;  	//safety distance read in angstrom than reduced
		  case 1: custparams[counter]=temp*0.0174533; break;		//rotation amplitude read in degrees, transformed to radians
		  default:
			  cout << "\n*****ERROR*****" << endl;
			  cout<<"\tMore parameter in makemovecus_CCl4.cpp than available slots."<<endl;
			  CleanExit();
		}
	      }
	      counter++;
	    }
	  }
	}
};	


//================================================================
//
//           The MakeMove function
//
//================================================================

//The moveout option is implemented, there can be more than one too close atom among the moved atoms.
void Move::MakeMove(SimpleCfg &config, FNC_POT &fnc, longint &seed)
{
	//local variables
	int i,j,it;
	//int count;//number of too close atom sin the molecule
	int across;//indicator if the molecule is across the limit of the cell
	int *pind, *pi3;
	double *pd1, *pd2, *pp1;
	double temp;
	double safedist;//safety distance to the cell border beyond which the 
				//molecule is within the box 
	double xc,yc,zc;
	double ampl;//move random amplitude
	double amplx, amply,amplz;//move amplitude along directions
	double theta, phi;//angles defing the raotation axis
	double ctheta,stheta,cphi,sphi, omega, comega, somega;//angular parameters
	double intpos1[15];//intermediate positions array
	double intpos2[15];//intermediate positions array

	

	ntcl_selected=0;//sets the the default for the number of the too close atoms among the moved atoms 
	for (i=0;i<5;i++)
	{
		tooclose_remove[i]=1;//the moved atom can be removed from the "too close" list by default,
		tooclose_index[i]=-1;//sets the default for the index of the moved atom among all the atoms of a configuration
						  //if the move is acceptable
		tcl_ind[i]=-1;//indicator showing, that this moved atom is not a too close by default
	}

	//defining the safety distance (largest possible C-Cl distance)
	safedist=*(custparams+4);//last element of the custom parameters array
							//(read from the .cus file by the constructor)


	//array of indices of moved atoms (defined by the FNC)
	pind=indices;//pointer to the indices of the moved atoms
	//choose the index of the molecule
	if(moveout==1)//if the moveout option is applied
	{
		//count=0;//default, number of tooclose atoms in the molecule
		if (Ran1(seed)<0.1)
		{
			//choose from the too close atoms
			itooclose[0]=int(ntooclose*Ran1(seed));//might be reset later, if the tooclose is not a carbon
			tooclose_index[0]=tclindices[itooclose[0]];//might be reset later, if the tooclose is not a carbon
			if (tooclose_index[0]<*config.pnatoms)
			{
				//the too close atom is carbon
				*pind=tooclose_index[0];//carbon of moved molecule
				i=*(fnc.converter+tooclose_index[0]);//offset of the first neighbour in the neighbours array 

			}
			else
			{
				//the too close atom is chlorine
				i=*(fnc.converter+tooclose_index[0]);//offset of the first neighbour (C atom) in the neighbours array
				*pind=fnc.neighbours[i];//index of the C atom, as it is always the first neighbour
				i=*(fnc.converter+ indices[0]);//offset of the first neighbour in the neighbours array 
			}
		}//chosen through too close atom
		else
		{
			//any molecule is chosen, can be too close atom among the atoms!
			*pind=int(*config.pnatoms*Ran1(seed));//[index between 0 and (numb. of C)-1 ]
	
			i=*(fnc.converter+*pind);//offset of the first neighbour in the neighbours array
		}//chosen randomly

		//has to set the variables connected with the tooclose atoms, if there is/are too close atom(s) in the molecule
		for (it=0;it<ntooclose;it++)
		{
			if (*pind==tclindices[it])//the central atom is a too close
			{
				itooclose[ntcl_selected]=it;//index of the too close atom in the tcl arrays
				tooclose_index[ntcl_selected]=tclindices[it];//sets the index for the index of the moved atom among all the atoms of a configuration
				ntcl_pair[ntcl_selected]=0;//reset the number of too close pairs moved above cutoff for the moved atom
				tcl_ind[0]=ntcl_selected;//indicator showing, that this moved atom is the ntcl_selected-th tooclose moved atom
				ntcl_selected++;//increase the number of tooclose atoms among the moved
				continue;//go to next too close
			}
			for (j=0;j<4;j++)//check the chlorines
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

		}//end of it cycle fortooclose atoms
	}//end of if moveout is applied
	else
	{
		*pind=int(*config.pnatoms*Ran1(seed));//[index between 0 and (numb. of C)-1 ]
	
		i=*(fnc.converter+*pind);//offset of the first neighbour in the neighbours array
	}


	pind++;//move to the second element of the indices array
	pi3=fnc.neighbours+i;//positioning the pointer in the neighbours array
	*pind++=*pi3++;//first neighbour=first chlorine
	*pind++=*pi3++;//second neighbour=second chlorine
	*pind++=*pi3++;//third neighbour=third chlorine
	*pind=*pi3;//fourth neighbour =fourth chlorine

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
		pd2=newpos+3;//points to the position of the first neighbour	
		for(i=1; i<nmoved;i++)//for all neighbours of the carbon
		{
			pd1=newpos;//points to the position of the first neighbour
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
		
	//translations (first because they don't change the orientation vectors)
	//amplitudes
	amplx=(2*Ran1(seed)-1)*(*(custparams));
	amply=(2*Ran1(seed)-1)*(*(custparams));
	amplz=(2*Ran1(seed)-1)*(*custparams);

	//translation along the x-axis
	for(i=0;i<nmoved;i++)
		*(newpos+3*i)+= amplx;
	//translation along the y-axis
	for(i=0;i<nmoved;i++)
		*(newpos+3*i+1)+= amply;
	//translation along the z-axis
	for(i=0;i<nmoved;i++)
		*(newpos+3*i+2)+= amplz;

	//individual atomic movements
	//carbon
	ampl=*(custparams+2);//amplitude of individual moves
	pd1=newpos;
	for(i=0;i<3;i++)// 1 carbon atoms
		*pd1++ +=(2*Ran1(seed)-1)*ampl;
	//chlorine atoms
	ampl=*(custparams+3);
	for(i=0;i<12;i++)// 3 times 4 chlorine atoms
		*pd1++ +=(2*Ran1(seed)-1)*ampl;

	//compute the  'center' position (carbon atom)
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
	theta=(2*Ran1(seed)-1)*PI;
	phi=2*Ran1(seed)*PI;

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

	omega=(2*Ran1(seed)-1)*(*(custparams+1));
	comega=cos(omega);
	somega=sin(omega);	
		
	//apply subrotations
	for(i=0;i<nmoved;i++)
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


	//translate back to original position
	//(the new positions are in intpos2)	
	pd2=intpos2;//points to intermediate positions array
	pd1=newpos;
	for(i=0;i<nmoved;i++)
	{
		*pd1++=*pd2++ +xc;
		*pd1++=*pd2++ +yc;
		*pd1++=*pd2++ +zc;
	}

	//apply toroidal cell structure
	pd1=newpos;
	for(i=0;i<15;i++)//3 times 5 atoms
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
		pp1=config.positions+3*(*pind++);//pointer to the old 
								//position of the moved atom
		for(j=0;j<3;j++)//for each coordinate
			*pp1++=*pd1++;//copy new positions
	}
	/*
	cout<<"final atoms positions:"<<endl;
	for(i=0;i<nmoved;i++)
	{
		for(j=0;j<3;j++)
		{
			cout<<J10<<*(newpos+3*i+j)<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	*/	

}//end of the Move::MakeMove function 	
