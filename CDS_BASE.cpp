#include "CDS_BASE.h"
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <limits>
//#include "mpi.h"
#define MAX_Partices      100000
#define Random_min       -0.05
#define Random_max        0.05 
#define MAX_LINE_LENGTH   200
#define gamma(i,j)        gamma[i][j]
#define PHI_old(i,j)      PHI_old[i][j]
#define PHI(i,j)          PHI[i][j]
#define Laplacian2(i,j)   Laplacian2[i][j]
#define POLYMER(i,j)      POLYMER[i*Ny+j]
#define PP3(i,j)          PP3[i][j]
#define PPP(i,j)          PPP[i][j]
#define P3(i,j)           P3[i][j]
#define P2(i,j)           P2[i][j]
#define CP4(i,j)          CP4[i][j]
#define P2_local(i,j)     P2_local[i*nlocaly +j]
#define index(ic,nc)      (ic[0]*nc[1] + ic[1])


BASE::BASE()
{
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);      
	double rationals[26];
	int integers[8];
	if(rank==0)
	{
		Read_input_parameters(integers, rationals);
	}
   MPI_Bcast(integers,8,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(rationals,26,MPI_DOUBLE,0,MPI_COMM_WORLD);
   Nx   = integers[0];     Ny    = integers[1]; //Nz=integers[2];
   Procx=integers[2];      Procy = integers[3]; //Procz =integers[5];
   Nbig =integers[4];      Nsmall= integers[5];
   Nparticles=integers[6];
   M = rationals[0]; u = rationals[1]; B = rationals[2];
   D = rationals[3]; A = rationals[4]; v = rationals[5];
   tau = rationals[6]; F = rationals[7]; r = rationals[8];
   Max_time_iteration=rationals[9];
   U1 = rationals[10]; CSI = rationals[11];
   ENNE = rationals[12];      ALPHA = rationals[13];
   beta = rationals[14];     GAMMs = rationals[15];
   GAMMb = rationals[16];    temp1 = rationals[17];
   R1small = rationals[18];   R1big = rationals[19];
   ALP1 = rationals[20];     dxx   = rationals[21];
   EMME = rationals[22];     sigma = rationals[23];
   P00s = rationals[24];     P00b  = rationals[25];
   //B2 =rationals[26]; 
	delta_t=0.01; 
   //Nbig = integers[4];
   //Nparticles=integers[5];
   RR1small = R1small * pow((1.0 + (1.0 / log( 2.0 ) ) ) ,(1.0 / ALP1));
   RR1big   = R1big * pow((1.0 + (1.0 / log(2.0))), (1.0/ALP1));
   RCUT = RR1big; //<---This is arbitrary
   /*--Compute the number of cells for linked cell lists---*/
   nc1 = new int[2];
   L   = new double[2];
   Num_cell_x = (int)floor(Nx / RCUT);
   Num_cell_y = (int)floor(Ny / RCUT);
	int coords[2], dims[2], periods[2];

      MPI_Comm new_comm = CreatCartesianTopology();
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Cart_get(new_comm, 2, dims, periods, coords);
      assert(Procx*Procy == size);
      printf("dimensions of topology = (%d, %d)\n",dims[0], dims[1]);
      
      if(Nx % Procx != 0)
	{
         if(coords[0] < (dims[0] - 1 ))
           {
             nlocalx=(int)Nx / Procx;
	   }
	 else
          {
            nlocalx = (int)Nx / Procx + (int)Nx%Procx;

          }
	}
      else
	{
		
          nlocalx = (int)Nx/Procx;
        }

      if(Ny % Procy != 0)
	{
	  if(coords[1] < (dims[1] - 1))
	    {
	      nlocaly = (int)Ny / Procy;
	    }
	 else
	    {
	      nlocaly = (int)Ny / Procy + (int)Ny % Procy;
	    }
	 }
      else
	 {
             nlocaly = (int)Ny / Procy;

          }
   
//----------------------------Calculate the local grid on a process------------------------------------//
      locala_x = 0.0 + coords[0] * nlocalx * (1.0 / Nx); //hx=1.0/Nx, [0.0, 1.0]x[0.0,1.0], global a =0.0, global b =1.0
      //localb_x= 


/*--Compute Number of linked cells in x and y direction -------*/
    /*  LinkedCells.resize(2);
      CellSize.resize(2);
      for (int direction=0; direction<2; direction++) 
          {
	    LinkedCells[direction] = (int) (al[direction]/RCUT) ; 
	    CellSize[direction]    = (double)(al[direction]/LinkedCells[direction]);
	    printf("Cellsize[%d]=%lf\n",direction, LinkedCells[direction]);
          }
      //---Add ghost cells to number of linked cells---/
      for (int direction=0; direction<2; direction++)
	   {
              LinkedCells[direction]=LinkedCells[direction] + 2;
	   }

      head.resize(LinkedCells[0]*LinkedCells[1]);
      LinkedCell_list.resize(Nparticles);
     */
    //-----------------------------------------------------------------
    //Establish the local number of particles on each MPI process
    //---------------------------------------------------------------
           /*
	   nbig_local= (int)(Nbig/size);
	   nbig_local=(my2drank<Nbig%size)?nbig_local +1 : nbig_local ; //<-----check this
	   nsmall_local= (int)(Nsmall/size);
	   nsmall_local=(my2drank<Nsmall%size)?nsmall_local +1 : nsmall_local ;

           nlocal_particles =Nparticles/size;
           N_start = (int)my2drank*nlocal_particles;
	   if(Nparticles!=(Nbig + Nsmall))
	   {
            printf("particles don't match; Nbig=%d, Nsmall=%d, Nparticles=%d\n", Nbig, Nsmall, Nparticles);
	    printf("Condition Nparticles = Nbig + Nsmall, is not satisfied");
	    assert(Nparticles==(Nbig + Nsmall));
	   }
	   */
       blocklengths[0] = 1;
       offsets[0] = 0;
       types[0] = MPI_INT;
       MPI_Type_extent(MPI_INT, &extent);
       offsets[1] = 1 * extent;
       types[1] = MPI_DOUBLE;
       blocklengths[1] = 9;
       MPI_Type_struct(2, blocklengths, offsets, types, &particletype);
       MPI_Type_commit(&particletype);
      /*     nlocal_particles =Nparticles/size;
           
	   N_start = (int)my2drank*nlocal_particles;
      if(Nparticles%size!=0)
	{
	  if(coords[0]<(dims[0]-1))
	    {
	      nlocal_particles = (int)Nparticles/size;
	      N_start = (int)my2drank*nlocal_particles;

	    }
	 else
	    {
	      nlocal_particles = (int)Nparticles/size + (int)Nparticles%size;
	      N_start = (int)my2drank*nlocal_particles;
	      
	    }
	 }
      else
	 {
             nlocal_particles = (int)Nparticles/size;
             N_start = (int)my2drank*nlocal_particles;

          }
         */
  

	   if(Nparticles % size)
	    {
             nlocal_particles = (my2drank < Nparticles % size) ? nlocal_particles + 1 : nlocal_particles ;
	     	N_start = N_start - ((my2drank < ((Nparticles) % size) ) ? my2drank : ((Nparticles) % size));
	    }
	   //printf("Rank=%d, N_start=%d, nlocalParticles=%d\n",my2drank, N_start,nlocal_particles);
      //     printf("Rank=%d,Nlocalx=%d, Nlocaly=%d, total=%d\n", my2drank, nlocalx,nlocaly,nlocalx*nlocaly);	    
            

  /*
  recvcounts=new int[size];
  displacement=new int[size];
  recvcounts_polymer=new int[size];
  disp_polymer      =new int[size];
 
  disp_polymer[0]=0;

  for(int i=0;i<size;i++)
     {
     recvcounts_polymer[i]=((i<(size-1))?(int)Nx/Procx:(int)Nx/Procx + (int)Nx%Procx)*((i<0)?(int)Ny/Procy:(int)Ny/Procy + (int)Ny%Procy) ;
  
  
  
    disp_polymer[i]= i*((int)Nx/Procx)*((int)Ny/Procy);//recvcounts_polymer[i-1] + disp_polymer[i-1];    
   
   printf("Rank=%d::disp_poly[%d]=%d, count[%d]=%d\n",my2drank,i, disp_polymer[i],i,recvcounts_polymer[i]); 

   }
    //printf("Rank=%d,countx=%d, county=%d, total=%d\n", my2drank, count_localx1,count_localy1,count_localx1*count_localy1);

*/

//--------------------------------------------------------------------------------------
    /*  MPI_Status status;
      partitionx =new int*[size];
      partitiony =new int*[size];
      for(int i=0;i<size;i++)
	 {
            partitionx[i]=new int[2];
	    partitiony[i]=new int[2];
	 }*/
     
      Nglobal_startx = 3 + (int) Nx / Procx*coords[0]; 
      Nglobal_endx = Nglobal_startx + nlocalx - 1;
   
      Nglobal_starty = 3 + (int) Ny / Procy * coords[1];
      Nglobal_endy = Nglobal_starty + nlocaly - 1;
      
/*      partitionx[my2drank][0]=Nglobal_startx;
      partitionx[my2drank][1]=Nglobal_endx;
      partitiony[my2drank][0]=Nglobal_starty;
      partitiony[my2drank][1]=Nglobal_endy;

      for(int i=0;i<size;i++)
	 {
           if(my2drank!=i)
	     {
                MPI_Sendrecv(partitionx[my2drank],2,MPI_INT,i,100,partitionx[i],2,MPI_INT,i,100,new_comm,&status);
		MPI_Sendrecv(partitiony[my2drank],2,MPI_INT,i,100,partitiony[i],2,MPI_INT,i,100,new_comm,            &status);    

	     }
	 }*/
  
      nc1[0] = (int)ceil(nlocalx / RCUT);
      nc1[1] = (int)ceil(nlocaly / RCUT);
      L[0] = (double)nlocalx;
      L[1] = (double)nlocaly;
      Grid.resize((nc1[0] + 2) * (nc1[1] + 2 ) ); //+2 ghost cells
      Grid.resize(nc1[0] + 2);
      for(int i = 0; i < (nc1[0] + 2); i++)
	 {
           Grid[i].resize(nc1[1] + 2);
	 }
      //printf("Rank=%d, (Nglobal_startx =%d, Nglobal_endx=%d)\n", my2drank, Nglobal_startx, Nglobal_endx);
      printf("RCUT=%lf, Number_cells_in_X=%d, Number_cells_in_Y=%d\n",RCUT,nc1[0], nc1[1]);
   
//------------------------------------------------------------------------------------

 

       PHI        = new double*[nlocalx+6];
       PHI_old    = new double*[nlocalx+6];
       gamma      = new double*[nlocalx+6];
       Laplacian2 = new double*[nlocalx+6];
       PP3        = new double*[nlocalx+6];
       P2         = new double*[nlocalx+6];//[Nx*Ny];

       //P2         = new double*[nlocalx+6];
       CP4        = new double*[nlocalx+6];
       //P2_local   = new double*[nlocalx+6];
       PPP        = new double*[nlocalx+6];
       P3         = new double*[nlocalx+6];
       PHI_local_result = new double[nlocalx*nlocaly];
       PPP_local        = new double[nlocalx*nlocaly];
       PPP_local2       = new double[nlocalx*nlocaly];
       P2_local         = new double[nlocalx*nlocaly];
       POLYMER = new double [Nx*Ny]; //<----used for storing polymer points for coupling.
       matrix_lower = new double[3*nlocaly];
       matrix_upper = new double[3*nlocaly];
       matrix_left  = new double[3*(nlocalx+6)];
       matrix_right = new double[3*(nlocalx+6)];
      // std::shared_ptr<Node> nd( new Node);//malloc(sizeof(Node));

       for(int i=0;i<nlocalx+6; i++)
	  {
	   PHI[i]            = new double[nlocaly+6];
	   PHI_old[i]        = new double[nlocaly+6]; 
	   gamma[i]          = new double[nlocaly+6];
	   Laplacian2[i]     = new double [nlocaly+6];
	   PP3[i]            = new double [nlocaly+6];
	  // PP3[i]  =new double [nlocaly+6];
	   P2[i]             = new double [nlocaly+6];
	   CP4[i]            = new double [nlocaly+6];
	   PPP[i]            = new double [nlocaly+6];
	   P3[i]             = new double [nlocaly+6];


          }
      

       for(int i = 0;i < nlocalx + 6; i++)
	   {
            for(int j = 0;j < nlocaly + 6; j++)
			{
				PHI[i][j] = 0.0;
				PHI_old[i][j] = 0.0;
	       		CP4[i][j] = 0.0;
	       		PP3[i][j] = 0.0;
	       		P2[i][j]  = 0.0;
	       		gamma[i][j] = 0.0;
	       		Laplacian2[i][j] = 0.0;
	       		PPP[i][j] = 0.0;
	       		P3[i][j] = 0.0;

	    }

	   }
       //std::memset(P2_local, 0, sizeof(double)*nlocalx*nlocaly);
       std::memset(PPP_local, 0, sizeof(double)*nlocalx*nlocaly);
       std::memset(PPP_local2, 0, sizeof(double)*nlocalx*nlocaly);
       //std::memset(PPP, 0, sizeof(double)*Nx*Ny);
       //std::memset(PP3, 0, sizeof(double)*Nx*Ny);

    // std::vector<int> counter;
    number_particles.resize(size);



//   XI.resize(nlocal_particles);
   bd1.resize(Nparticles);
   //bd.resize(Nparticles); 
   //bd = new BODY[Nparticles];
   //bd1= new BODY[Nparticles];
//   send_list  =new BODY*[size];
//   recv_list  =new BODY*[size];
   send_buffer=new BODY*[size];
   for(int i=0;i<size;i++)
     {
  //     send_list[i]  =new BODY[Nparticles];
   //    recv_list[i]  =new BODY[Nparticles];
       send_buffer[i]=new BODY[Nparticles];
     }
   Posx.resize(Nparticles);
   Posy.resize(Nparticles);
   Fx.resize(Nparticles);
   Fy.resize(Nparticles);
   Vx.resize(Nparticles);
   Vy.resize(Nparticles);
   Global_Index_X.resize(Nx);
   Global_Index_Y.resize(Ny);
   UPX.resize(Nx);
   DOWNX.resize(Nx);
   //Recv_bodies= new BODY[nlocal_particles];
   //Send_bodies= new BODY[nlocal_particles];
   //down_nbr_list.resize(nlocal_particles);
   //right_nbr_list.resize(nlocal_particles);
   //up_nbr_list.resize(nlocal_particles);
   //left_nbr_list.resize(nlocal_particles);
   //Virtual_list.resize(nlocal_particles);  
   //Recv_bodies.resize(nlocal_particles);
   //send_bodies.resize(nlocal_particles);
   R12.resize(Nparticles);
   Dx.resize(Nparticles);
   Dy.resize(Nparticles);
   RI.resize(Nparticles);
   Dx1.resize(Nx);
   Dy1.resize(Nx);
   //PPP.resize(nlocalx);
   POS.resize(Nx);//<--- rethink about the size of this, it is bether to use a linked list
   VerletList_local.resize(Nparticles);
   VerletTable_local.resize(Nparticles);
   TOTX.resize(10*Nparticles);
   TOTY.resize(10*Nparticles);
   for(int i=0;i<Nx;i++)
      {
	POS[i].resize(Ny);
	Dx1[i].resize(Ny);
	Dy1[i].resize(Ny);
      }

   for(int i=0; i<Nparticles;i++)
      {
        Dx[i].resize(Nparticles);
        Dy[i].resize(Nparticles);
        R12[i].resize(Nparticles);
        RI[i].resize(Nparticles);
        VerletTable_local[i].resize(Nparticles);

        }
     


  up1x.resize(nlocaly+6);
  down1x.resize(nlocaly+6);


     printf("Base Construction made. \n");
  FreeCommunicator(new_comm); 

 }

BASE::~BASE()
     {
        for(int i=0;i<nlocalx+6;i++)
           {
             delete [] PHI[i];
	     delete [] PHI_old[i];
	     delete [] gamma[i];   
	     delete [] Laplacian2[i];
	     delete [] PP3[i];
	     delete [] P2[i];
	     delete [] CP4[i];
	     delete [] PPP[i];
	     delete [] P3[i];
	     //delete [] P2_local[i];

	   }
	for(int i=0;i<size;i++)
	   {
           //  delete [] send_list[i];
             //delete [] recv_list[i];
	     delete [] send_buffer[i];
           }
	//  delete [] send_list;
	  delete [] send_buffer;
        //  delete [] recv_list;
	  delete [] PHI;
	  delete [] PHI_old;     
	  delete [] Laplacian2;
	  delete [] gamma;	    
	  delete [] P3;
	 // delete [] bd;
	 // delete [] bd1;
	 // delete [] Recv_bodies;
	 // delete [] Send_bodies;
	 // delete [] recvcounts;
	 // delete [] displacement;
	//  delete [] recvcounts_polymer;
	//  delete [] disp_polymer;
	  delete [] PP3;
	  delete [] P2;
	  delete [] CP4;
	  delete [] PPP;
	  delete [] PHI_local_result;
	  delete [] POLYMER;
	  delete [] PPP_local;
          delete [] P2_local;
	  delete [] PPP_local2;

	  delete [] matrix_lower;
	  delete [] matrix_upper;
	  delete [] matrix_left;
	  delete [] matrix_right;
	  delete [] nc1;
	  delete [] L;

	 // delete  nd;
	  std::cout<<"Destructor in base class called"<<std::endl;
	//  std::cout<<"Destructor in base class called: colloids destroyed"<<std::endl; 

          //MPI_Comm_free(&new_comm);
	//  FreeCommunicator(new_comm);


     }
void BASE::InitializePolymer(MPI_Comm new_comm)
    {

      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);

     
      srand(time(NULL) + my2drank);
       for(int i=3; i<=nlocalx+2;i++)
            {
             for(int j=3;j<=nlocaly+2;j++)
	        {
	         double range =Random_max-Random_min;
	         double div =RAND_MAX / range;
	         PHI_old[i][j]=Random_min + (rand()/div);
		// PHI_local_result[(i-3)*nlocaly + j-3]= P00b;//PHI_old[i][j];
                 //PPP_local[(i-3)*nlocaly + j-3]=PPP_local[(i-3)*nlocaly + j-3] + B2*(PHI_old[i][j]-P00b);
		 //printf("Rank=%d, PHI_old[%d][%d]=%lf\n", my2drank, i,j, PHI_old[i][j]);
	        }
	    }
   
     }

void BASE::setIndex(MPI_Comm new_comm)
{
   int coords[2], dims[2], periods[2];
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_get(new_comm,2,dims,periods, coords );
   if(dims[1]==1)
   {
   for(int s=3; s<nlocaly+3; s++)
      {    
        up1x[s]=s+1;
        down1x[s]=s-1;
      }
    down1x[3]=nlocaly+2;
    up1x[nlocaly+2]=3;
   }
   else
   {
    for(int s=3; s<nlocaly+3; s++)
       {    
         up1x[s]=s+1;
         down1x[s]=s-1;
       }
   }
 }

void BASE::setGlobalIndex(MPI_Comm new_comm)
{
    MPI_Comm_size(new_comm, &size);

    for(auto i=0;i<Nx;i++)
       {
	UPX[i] =i+1;
	DOWNX[i]=i-1;
       }
    UPX[Nx-1]=0;
    DOWNX[0]=Nx-1;


}

void BASE::setLaplacianBase(MPI_Comm new_comm)
    {
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
         
         setIndex(new_comm);
         //cout<<"gamma(0,0)="<<gamma(0,0)<<"\t";                                                                       
         double AP=0.0, BP=0.0,  ATP=0.0;//gtemp=0.0;                                                           
         for(int i=3;i<=nlocalx+2;i++)
            {
            for(int j=3;j<=nlocaly+2;j++)
	       {
	        AP=(1.0/6.0)*(PHI_old(i+1,j) + PHI_old(i-1,j)	+ PHI_old(i,j+1) + PHI_old(i,j-1));
                BP=(1.0/12.0)*(PHI_old(i-1,j+1) + PHI_old(i-1,j-1) + PHI_old(i+1,j+1) + PHI_old(i+1,j-1));
                ATP = AP + BP;
	        gamma(i,j)=g(PHI_old(i,j))-PHI_old(i,j)+ dxx*D*(ATP-PHI_old(i,j));
               }
            }
    }
void BASE::SetCP4(MPI_Comm new_comm)
{
	 MPI_Comm_rank(new_comm, &my2drank);
	 MPI_Comm_size(new_comm, &size);
         //polymer colloid coupling                                                                                
         double AP3(0.0), BP3(0.0);
         for(int i=3;i<=nlocalx+2;i++)
            {
             for(int j=3;j<=nlocaly+2;j++)
	       {
                 CP4(i,j)=0.0;
	        AP3=(1.0/6.0)*(P2(i+1,j)+P2(i-1,j) + P2(i,j+1) + P2(i,j-1));
	        BP3=(1.0/12.0)*(P2(i-1,j-1) + P2(i-1,j+1) + P2(i+1,j-1) + P2(i+1,j+1));
                CP4(i,j)=AP3 + BP3-P2(i,j);
               //  printf("CP4(%d,%d)=%lf\n",i,j,CP4(i,j));                                                         
	       }
            }
        
 }

double BASE::g(double phi)
{
    double q(0.0);
    q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
    return q;
}

void BASE::SetSecondLaplacian2(MPI_Comm new_comm)
{
      int coords[2], dims[2], periods[2];
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
      MPI_Cart_get(new_comm,2,dims,periods, coords );
              setIndex(new_comm);
          double AP=0.0, BP=0.0;
          for(int i=3;i<=nlocalx+2;i++)
             {
	      for(int j=3;j<=nlocaly+2;j++)
	         {
	           AP=(1.0/6.0)*(gamma(i+1,j)    + gamma(i-1,j) + gamma(i,j+1)    + gamma(i,j-1));
	           BP=(1.0/12.0)*(gamma(i-1,j+1) + gamma(i-1,j-1) +gamma(i+1,j+1)  + gamma(i+1,j-1));
	           Laplacian2(i,j) = AP + BP;
	         }
              }
}
void BASE::UpdateSolution(MPI_Comm new_comm)
{
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
      
          for(int i=3;i<nlocalx+3;i++)
	     {
	      for(int j=3;j<nlocaly+3;j++)
	         {
	          PHI_old[i][j]=PHI[i][j];
	          P3[i][j]=PHI_old[i][j];
		//  PPP_local[(i-3)*nlocaly + j-3]=PPP_local[(i-3)*nlocaly + j-3] + B2*(PHI_old[i][j]-P00b);
	         }
	     } 
   

  }

void BASE::FiniteDifferenceScheme(MPI_Comm new_comm)
{
      int coords[2];
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
      MPI_Cart_coords(new_comm,my2drank,2,coords);


      std::cout<<"in finite difference scheme parallel"<<std::endl;
     
	
         for(int i=3;i<nlocalx+3;i++)
            {
             for(int j=3;j<nlocaly+3;j++)
                {
		// PHI[i][j]= PHI_old(i,j) - B*(PHI_old(i,j)-1.0 + 2*r)+ gamma(i,j)-Laplacian2(i,j);
	         // int X=i-3+coords[0]*nlocalx; //<--global x index
		 // int Y=j-3+coords[1]*nlocaly;
		  PHI[i][j]= PHI_old(i,j)+EMME*delta_t*(-B*(1.0-PPP(i,j))*PHI_old(i,j)+ (dxx)*(gamma(i,j)-Laplacian2(i,j)) + (dxx)*CP4[i][j]);
	          PHI_local_result[(i-3)*nlocaly + j-3]= PHI[i][j];
	//	  printf("B=%lf\n",B);
                // printf("PHI[%d][%d]=%lf\n",i,j,PHI[i][j]); 
		}
             }
        
}

void BASE::ExchangeData(MPI_Comm new_comm, double **array)
{
     MPI_Status status;
     MPI_Comm_rank(new_comm,&my2drank);
     int coords[2], upper_nbr(0), down_nbr(0), tag, left_nbr(0), right_nbr(0) ;
     tag=201;
     int dims[2];
     int periods[2];
     MPI_Cart_coords(new_comm,my2drank,2,coords);
     MPI_Cart_get(new_comm,2,dims,periods, coords );
     //copy boundary points for exchange among mpi processes
    if((dims[0]>1)&&(dims[1]>1))
      {
       if(dims[0]>1)
        {
	   for(int i=3;i<6;i++)
	      {
                for(int j=3;j<=nlocaly+2; j++)
                   {
                    matrix_upper[(i-3)*nlocaly + j-3]=array[i][j];
		    //matrix_upper[j-1]=array[1][j];

                   }
	      } 
	   for(int i=nlocalx; i<nlocalx+3 ; i++)
	      {
               for(int j=3;j<=nlocaly+2; j++)
   	          {
                   matrix_lower[(i-nlocalx)*nlocaly + j-3]=array[i][j];
		   // matrix_lower[j-1]=array[nlocalx][j];
                  }
	      }
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr ); //move along x direction
        MPI_Sendrecv_replace(matrix_upper, 3*nlocaly,  MPI_DOUBLE, upper_nbr, tag, down_nbr,tag, new_comm,&status);
        MPI_Sendrecv_replace(matrix_lower, 3*nlocaly,  MPI_DOUBLE, down_nbr, tag, upper_nbr,tag, new_comm,&status);
        //copy received ghost points into matrix, array
	for(int i=0;i<3;i++)
	   {
            for(int j=3; j<=nlocaly+2; j++)
               {
                array[i][j]=matrix_lower[(i)*nlocaly+j-3];
	       }
	   }
	for(int i=nlocalx+3; i<nlocalx+6; i++)
	   {
	     for(int j=3;j<nlocaly+3;j++)
		 {
                  array[i][j]=matrix_upper[(i-(nlocalx+3))*nlocaly + j-3];
		 }
	   }
     	   
      }
     if(dims[1]>0)
       {
         // for(int j=2;j<4;j++)
           //  {
	       for(int i=0;i<=nlocalx+5;i++)
                  {
	             matrix_left[i]=array[i][5];
                     matrix_left[(nlocalx+6)+i]=array[i][4];
		     matrix_left[2*(nlocalx+6)+i]=array[i][3];

		    }
             //   }
	  

        //for(int j=nlocaly;j<nlocaly+2;j++)
	  // {
             for(int i=0;i<=nlocalx+5;i++)
  	        {
		 matrix_right[ i]=array[i][nlocaly];
                 matrix_right[ (nlocalx+6)+i]=array[i][nlocaly+1];
		 matrix_right[2*(nlocalx+6) + i]=array[i][nlocaly+2];

		}
          // }
	 MPI_Comm_rank(new_comm,&my2drank);
	 MPI_Cart_coords(new_comm,my2drank,2,coords);
	 MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
	 MPI_Sendrecv_replace(matrix_right, 3*(nlocalx+6),  MPI_DOUBLE, right_nbr, tag, left_nbr,tag, new_comm,&status);
         MPI_Sendrecv_replace(matrix_left, 3*(nlocalx+6),  MPI_DOUBLE, left_nbr, tag, right_nbr,tag, new_comm,&status);
        //copy received ghost points into matrix
         
	    //for(int j=0;j<2;j++)
              // {
	        for(int i=0;i<=nlocalx+5; i++)
     	           {
         	     array[i][0]=matrix_right[i];
		     array[i][1]=matrix_right[(nlocalx+6)+i];
                     array[i][2]=matrix_right[2*(nlocalx+6)+i];
	           }
	      // }
	    //for(int j=nlocaly+3;j>(nlocaly+1);j--)
	      // {
	         for(int i=0;i<nlocalx+6;i++)
	            {
	             array[i][nlocaly+5]=matrix_left[i];
		     array[i][nlocaly+4]=matrix_left[(nlocalx+6)+i];
		     array[i][nlocaly+3]=matrix_left[2*(nlocalx+6)+i];

    	            }
	   
           //}
   }
 }
 else if((dims[0]>1)&&(dims[1]==1))
        {
         for(int j=3;j<=nlocaly+2; j++)
            {
              matrix_upper[ j-3]         =array[5][j];
	      matrix_upper[nlocaly + j-3]=array[4][j];
	      matrix_upper[2*nlocaly+ j-3]=array[3][j];
            }
	      //} 
	  // for(int i=nlocalx; i<nlocalx+2 ; i++)
	   //   {
               for(int j=3;j<=nlocaly+2; j++)
   	          {
                    matrix_lower[ j-3]         = array[nlocalx][j];
		    matrix_lower[nlocaly + j-3]= array[nlocalx+1][j];
                    matrix_lower[2*nlocaly+j-3]= array[nlocalx+2][j];
                  }
	     // }
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr ); //move along x direction
        MPI_Sendrecv_replace(matrix_upper, 3*nlocaly,  MPI_DOUBLE, upper_nbr, tag, down_nbr,tag, new_comm,&status);
        MPI_Sendrecv_replace(matrix_lower, 3*nlocaly,  MPI_DOUBLE, down_nbr, tag, upper_nbr,tag, new_comm,&status);
        //copy received ghost points into matrix, array
	//for(int i=0;i<2;i++)
	  // {
            for(int j=3; j<=nlocaly+2; j++)
               {
                array[0][j]=matrix_lower[j-3];
		array[1][j]=matrix_lower[nlocaly+j-3];
		array[2][j]=matrix_lower[2*nlocaly+j-3];
	       }
	  // }
      //	for(int i=nlocalx+2; i<nlocalx+4; i++)
	//   {
	     for(int j=3;j<nlocaly+3;j++)
		 {
                   array[nlocalx+5][j] = matrix_upper[ j-3];
	           array[nlocalx+4][j] = matrix_upper[nlocaly + j-3];
		   array[nlocalx+3][j] = matrix_upper[2*nlocaly + j-3];
		   
		 }
	  // }
	
        }
 else if((dims[0]==1)&&(dims[1]>1))
       {
        std:: cout<<"Not yet implemented. Realigned the processors. That is, let Procy=1 and Procx=size"<<std::endl;
       }

    }



 

//-----------------END of POLYMER PART--------------------------------------------------------
void BASE::ArrangeParticles_Cells(std::vector<std::vector<std::list<BODY> >  > &grid, int *nc,std::list<BODY> &colloid_list, double *l,MPI_Comm new_comm)
{
   int ic[2], kc[2], count=0;
   std::list<BODY>::iterator it;
   std::list<BODY> temp;

   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   int coords[2], dims[2],periods[2];
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   MPI_Cart_get(new_comm, 2,dims, periods,coords);
   //int OFFSET[2];
   const int OFFSETX= coords[0]*(Nx/Procx);
   const int OFFSETY= coords[1]*(Ny/Procy);
 //  double local_coordsx= bd[i].r[0]-floor(coords[0]*(OFFSETX));

 //reset grid
  for(ic[0]=0;ic[0]<nc[0];ic[0]++)
     {
      for(ic[1]=0;ic[1]<nc[1];ic[1]++)
         {
           grid[ic[0]][ic[1]].clear(); //<----reset the list on every grid point to zero(ie no element)       
         }
     }
  
  for(it=colloid_list.begin();it!=colloid_list.end();++it)
     {
      //for(int d=0;d<2;d++)
        // {
	 // if((it->r[0]!=0.0)&&(it->r[1]!=0.0))
	   // {
              //std::cout<<"it-r[x]="<<it->r[0]-(double)OFFSETX<<std::endl;
	      //std::cout<<"it-r[y]="<<it->r[1]-(double)OFFSETY<<std::endl;
          //kc[d]=(int)(floor((it->r[d]-OFFSET[d])*nc[d]/l[d] ));
          // printf("kd[%d]=%d\n",d,kc[d]);
         //}
             kc[0]=(int)(floor((it->r[0]-(double)OFFSETX)*(nc[0]/l[0]) ));
             kc[1]=(int)(floor((it->r[1]-(double)OFFSETY)*(nc[1]/l[1] )));


             //int Gridsize=(int)grid.size();
             //printf("Rank=%d, kd[0]=%d, kd[1]=%d, cell index=%d, GridSize=%d\n",my2drank, kc[0],kc[1], index(kc,nc),Gridsize);
	  //   if((index(kc,nc))==Gridsize)
	    //   {
              //  printf("Rank=%d, kd[0]=%d, kd[1]=%d, cell index=%d\n",my2drank, kc[0],kc[1], index(kc,nc));
	    //   }
             grid[kc[0]][kc[1]].push_back(*it);
            // count++;
         // }
     }

}






void BASE::ComputeForce(MPI_Comm new_comm)
{
   int coords[2],dims[2],periods[2];
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   //MPI_Status status;
   //MPI_Request request;

   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_get(new_comm,2,dims,periods, coords);
   double G4(0.0), AA(0.0), BB(0.0), NEW1(0.0), NEW2(0.0), G2(0.0), G3(0.0), G3A(0.0),C(0.0);
   double DSUMX(0.0),DSUMY(0.0), DSUMXY(0.0);

for(int i=N_start;i<N_start+ nlocal_particles-1;i++)
     {
     for(int j=i+1;j<Nparticles;j++)
        {
               
	  double dx=bd[i].r[0]-bd[j].r[0];
	  double dy=bd[i].r[1]-bd[j].r[1];
	  //printf("dx=%lf, dy=%lf\n", dx,dy);
	 
	  double R2=sqrt(pow(dx,2)+pow(dy,2));
	  if(R2<2.0*R1big)
	    {
	     //std::cout<<"I am here"<<std::endl;
	     Dx[i][j]=dx;
	     Dy[i][j]=dy;
	     RI[i][j]=R2;
	     RI[j][i]=R2;
	     //C = U1/R12[i][j];
	     C = U1/(2.0*R1big);
	     //AA=RI[i][j]/R12[i][j];
	     AA=R2/(2.0*R1big);
	     G2 = pow(AA + beta,2.0);
	     G3 = 1.0 + ALPHA*(AA + beta);
	     BB = (AA-1.0)*(-ALPHA);
	     G3A = exp(BB);
	     NEW1=pow((1.0 + (ENNE/ALPHA)+ beta),2.0);
	     NEW2=(exp(-ENNE))*((1.0 + ALPHA*(beta + (ENNE/ALPHA)+1.0)))/(NEW1);
	     G4 = ((C*(G3A)*G3)/G2)-(C*NEW2);
	     bd[i].Force[0]=bd[i].Force[0] + G4*(dx/R2);//(Dx[i][j]/RI[i][j]);
	     bd[i].Force[1]=bd[i].Force[1] + G4*(dy/R2);//(Dy[i][j]/RI[i][j]);
	     //printf("Colloid Force (fx=%lf, fy=%lf) \n", bd[i].Force[0], bd[i].Force[1]);
	     bd[j].Force[0]=bd[j].Force[0] + G4*(Dx[j][i])/RI[i][j];
	     bd[j].Force[1]=bd[j].Force[1] + G4*(Dy[j][i])/RI[i][j];

	   }
	 }
   }
   //Compute random force
    double csi1(0.0), csi2(0.0);
    srand(time(NULL) + my2drank);
    double GAMM(0.0);
    for(int i=N_start;i<N_start + nlocal_particles; i++)
       {
        double range =Random_max-Random_min;
        double div =RAND_MAX / range;
        csi1 = Random_min + (rand()/div);
        csi2 = Random_min + (rand()/div);
    //  printf("Random number (csi1=%lf, csi2=%lf) \n", csi1, csi2);
        if(bd[i].index==0)
          {
       	   GAMM= GAMMs;
           bd[i].v[0]= (1.0/GAMM)*(bd[i].Force[0]+ bd[i].TOT[0])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
	   bd[i].v[1]= (1.0/GAMM)*(bd[i].Force[1]+ bd[i].TOT[1])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
	 // printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);
	  }
	if(bd[i].index==1)
	  {
	   GAMM= GAMMb;
	   bd[i].v[0]= (1.0/GAMM)*(bd[i].Force[0]+ bd[i].TOT[0])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
	   bd[i].v[1]= (1.0/GAMM)*(bd[i].Force[1]+ bd[i].TOT[1])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
	   //printf("Colloid (TOTX =%lf, TOTY=%lf) \n", TOTX[i], TOTY[i]);
	   //printf("Colloid force (fx=%lf, fy=%lf) \n", bd[i].Force[0], bd[i].Force[1]);

           }
	   //position update
	  bd[i].r[0] =bd[i].r[0] + bd[i].v[0]; //<---bd[i].vx is already multiplied by t. 
	  bd[i].r[1] =bd[i].r[1] + bd[i].v[1];

          //bc condition
          
	  
	   bd[i].r[0]=bd[i].r[0]-Nx*( (int)(bd[i].r[0]/Nx + 1)-1);
	 
	   bd[i].r[1]=bd[i].r[1]-Ny*( (int)(bd[i].r[1]/Nx + 1)-1);
           //if(my2drank==0)	  
 	   //printf("Force_x[%d]=%lf, Force_y[%d]=%lf\n",i, bd[i].Force[0],i,bd[i].Force[1]);

	
        //computing mean displacement
        Posx[i]=Posx[i] + bd[i].v[0];
        Posy[i]=Posy[i] + bd[i].v[1];
        DSUMX  = DSUMX + pow(Posx[i],2.0);
        DSUMY  = DSUMY + pow(Posy[i],2.0);
        DSUMXY = DSUMXY + Posx[i]*Posy[i];
      //printf("Colloid positions (bd[100].x=%lf, bd[100].y=%lf) \n", bd[100].r[0], bd[100].r[1]);
    }
    // DSUMX=DSUMX/
}





void BASE::WriteToFile_MPI(MPI_Comm new_comm)
{
   MPI_Status status;
   MPI_File     fh;
   MPI_Datatype filetype;
   int dims[2],coords[2], periods[2], start_indices[2], localsizes[2];
   int globalsizes[2];
   globalsizes[0]=Nx;
   globalsizes[1]=Ny;
   localsizes[0]=(Nx/Procx); localsizes[1]=(Ny/Procy) ;
   int local_array_size =nlocalx*nlocaly;
   MPI_Cart_get (new_comm,2,dims,periods, coords );
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   start_indices[0]=coords[0]*localsizes[0];
   start_indices[1]=coords[1]*localsizes[1];
   printf("rank=%d,r=%d,start_indices=(%d,%d)\n",my2drank,Nx%Procx,start_indices[0],start_indices[1]);
   MPI_Type_create_subarray(2, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);
   MPI_Type_commit(&filetype);
   char outputfilename[]="datafile_mpi";
   char filememview[]="native";
   MPI_File_open(new_comm, outputfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
   MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, filememview,MPI_INFO_NULL);
   MPI_File_write_all(fh, PHI_local_result,local_array_size , MPI_DOUBLE, &status);
   MPI_File_close(&fh);
   MPI_Type_free(&filetype);

}


void BASE::Read_input_parameters(int *integers, double *rationals)
{
       //int value(0);
       FILE* file;
       char Data[MAX_LINE_LENGTH], *string;
	    if((file=fopen("ParameterFile.dat","r"))==NULL)
              {
	         printf("Error opening ParameterFile.dat\n");
	         return;
	      }
        
	string=fgets(Data,MAX_LINE_LENGTH,file);
          		
	int value=fscanf(file,"%d\n",&integers[0]);

	 string=fgets(Data,MAX_LINE_LENGTH,file);
         
	 value=fscanf(file,"%d\n",&integers[1]);
	 
         string=fgets(Data,MAX_LINE_LENGTH,file);
	 
	 value=fscanf(file,"%d\n",&integers[2]);
	 
	 string=fgets(Data,MAX_LINE_LENGTH,file);
         
	 value=fscanf(file,"%d\n",&integers[3]);
	 
	 string=fgets(Data,MAX_LINE_LENGTH,file);
	 
	 value=fscanf(file,"%d\n",&integers[4]);
	 
	string=fgets(Data,MAX_LINE_LENGTH,file);
         
	value=fscanf(file,"%d\n",&integers[5]);
	 
	string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%d\n",&integers[6]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%d\n",&integers[7]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
 	value=fscanf(file,"%lf\n",&rationals[0]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[1]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[2]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[3]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[4]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[5]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[6]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[7]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[8]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%lf\n",&rationals[9]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[10]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%lf\n",&rationals[11]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[12]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[13]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[14]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[15]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[16]);

	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[17]);

	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[18]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[19]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[20]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[21]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[22]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[23]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[24]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[25]);
//	string=fgets(Data,MAX_LINE_LENGTH,file);
//        value=fscanf(file,"%lf\n",&rationals[26]);
	
	
       // fgets(Data,MAX_LINE_LENGTH,file);
	//printf("value=%d,character=%s\n",value,string[0]);
	if((value!=0)||(string!=NULL))
	  {
           printf("failed to read in\n");
	  }
        fclose(file);
 
}



MPI_Comm BASE::CreatCartesianTopology()
      {
         MPI_Comm   new_comm;
         if(Procx*Procy!=size)
           {
             printf("Number of processors is not factorizable. That is Px*Py=%d, is not equal to the number of processors=%d\n ",Procx*Procy,size);
             printf("Hint: Check Parameter.dat file to change the number of cores\n");

		
	    exit(0);
	  }
        const int ROWS=0;const int COL=1;	     
        int periods[2];
        periods[0]=1; periods[1]=1;
        int coords[2];
        int dims[2];
        dims[ROWS]=Procx;
        dims[COL] =Procy;
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Cart_create( MPI_COMM_WORLD, 2,dims, periods, 0, &new_comm );
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Comm_size( new_comm, &size );

         return new_comm;

 }

void BASE::GenerateColloids(std::list<BODY>& bodies_list,MPI_Comm new_comm)
{
   
    std::cout<<"----------------------- Generating Colloids--------------------"<<std::endl;
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);
    MPI_Status status;

    int coords[2], periods[2];
    int dims[2];
    MPI_Cart_coords(new_comm,my2drank,2,coords);
    MPI_Cart_get(new_comm,2,dims,periods, coords );

    int tag=101;
    //MPI_Status status;
    
    // srand(time(NULL) + my2drank);
     char fname[200];
     sprintf(fname, "solution%d.dat",my2drank );
     FILE *fp;
     fp=fopen(fname, "w");
     if(my2drank==0)
       {
           std::random_device rd; // obtain a random number from hardware
           std::mt19937 eng(rd() ); // seed the generator
           std::uniform_real_distribution<> distr(3.0, (double)(Nx-3)); //distr(3.0, (double)(Nx-3)) 
  
          srand(time(NULL) + my2drank);
          
               
      	
           bd1[0].r[0]=distr(eng);//(((double)rand()/RAND_MAX))*(Nx-1);//*nlocalx + (double)coords[0]*(nlocalx); 
           bd1[0].r[1]=distr(eng); //(((double)rand()/RAND_MAX))*(Ny-1);//*nlocaly+ (double)coords[1]*(nlocaly);
		   //(double)coords[0]*nlocalx*(1.0/Nx) is the offset in x direction
           bd1[0].v[0]=0.0;
           bd1[0].v[1]=0.0;
           bd1[0].index=1;
	   bd1[0].TOT[0]=0.0;
           bd1[0].TOT[1]=0.0;
	   bd1[0].Force[0]=0.0;
	   bd1[0].Force[1]=0.0;
	    
	    // printf("Colloid (rank=%d, dx=%lf, dy=%lf) \n",my2drank, bd[i].r[0], bd[i].r[1]);
        int i=0;       
       for(i=1; i<Nparticles; i++)
          {
label1:     bd1[i].r[0]= distr(eng);//(((double)rand()/RAND_MAX))*(Nx-1);//(nlocalx)+ (double)coords[0]*(nlocalx) ; 
	    bd1[i].r[1]= distr(eng);//(((double)rand()/RAND_MAX))*(Ny-1);//nlocaly + (double)coords[1]*(nlocaly) ;
	    bd1[i].v[0]=0.0;
            bd1[i].v[1]=0.0;
            bd1[i].index=1;
	    bd1[i].TOT[0]=0.0;
	    bd1[i].TOT[1]=0.0;
            bd1[i].Force[0]=0.0;
	    bd1[i].Force[1]=0.0;

           for(int j=0;j<i;j++)
	      {
               double RX=bd1[i].r[0]-bd1[j].r[0];
	       double RY=bd1[i].r[1]-bd1[j].r[1];
	       double R2=sqrt(pow(RX,2)+pow(RY,2));
	       if(R2<RR1big*2.0)//<---they should not occupy the same position
	         {
		   //std::cout<<"Going to label 1"<<std::endl;
		   goto label1;
		  }
	      }
	      // fprintf(fp, "%d  %lf   %lf\n",i, bd[i].r[0], bd[i].r[1]);
	   }
         
              
          for(int i=0;i<Nparticles;i++)
	     {
               int processID=Get_Process_ID(i,bd1);
	   //    printf("ProcessID=%d\n",processID);
	       send_buffer[processID][number_particles[processID]]=bd1[i];
	       number_particles[processID]++;
	       //printf("counter[%d]=%d\n",processID,counter[processID]);

	     }
          for(int cores=0;cores<size;cores++)
	     {
               MPI_Send(&number_particles[cores],1,MPI_INT,cores,tag,new_comm);
	       MPI_Send(send_buffer[cores],number_particles[cores],particletype,cores,tag,new_comm);
	     }
      }

  //else
 // {
       MPI_Recv(&number_particles[my2drank],1,MPI_INT,0,tag,new_comm,&status);
       std::vector<BODY> recv_buffer;
       recv_buffer.resize(number_particles[my2drank]);
       //bd.resize(counter[my2drank]);
    //   printf("counter=%d\n",counter[my2drank]);
       MPI_Recv(&recv_buffer[0],number_particles[my2drank],particletype,0,tag,new_comm,&status);
      //MPI_Recv(&bd[0],number_particles[my2drank],particletype,0,tag,new_comm,&status);
       //bodies_list.resize(number_particles[my2drank]);
       //std::copy(recv_buffer.begin(), recv_buffer.end(), std::back_inserter(bodies_list));
       for(int i=0;i<number_particles[my2drank];i++)
	  {
            bodies_list.push_back(recv_buffer[i]);

	  }
      printf("Rank=%d, counter[%d]=%d\n",my2drank, my2drank,number_particles[my2drank]);
       int count_particles=0;
       std::list<BODY>::iterator iter; 
       for(iter= bodies_list.begin();iter!=bodies_list.end();++iter)               //(int i=0; i<number_particles[my2drank];i++)
          {
	    
            fprintf(fp, "%d   %lf   %lf\n",count_particles, iter->r[0], iter->r[1]);
	    count_particles++;
          }
      

  //}
//       for(int i=0; i<counter[my2drank];i++)
//	   {
  //            fprintf(fp, "%d  %lf   %lf\n",i, recv_buffer[i].r[0], recv_buffer[i].r[1]);
//	   }

        
 	 fclose(fp);
   
	 std::cout<<"---------------------End of Colloid Generation-----------------"<<std::endl;


}
int BASE::Get_Process_ID(int index,const std::vector<BODY> &bodies) const
{
 
  //int processID=0;	
/*  double xcoord=0.0, ycoord=0.0;
  int grid_location=0;
  xcoord=bd[index].r[0];
  ycoord=bd[index].r[1];

  for(int cores=0;cores<size;cores++)
     {
       if((xcoord<partitionx[cores][1])&&(xcoord>partition[cores][0]))
	 {
           if((ycoord<partitiony[cores][1])&&(ycoord>partitiony[cores][0]))
	     {
               grid_location=cores;
	     }
	  
	 }

     }
 */
   
  int  processID= ((int)(bodies[index].r[0]/(Nx/Procx)))*Procy + (int)(bodies[index].r[1]/(Ny/Procy));

return processID;

}


 void  BASE::Solver()
  {
    double t=0.0;
    int count =0;
    MPI_Comm new_comm;
    new_comm = CreatCartesianTopology();
    InitializePolymer(new_comm);
    GenerateColloids(particles_list, new_comm);
    
    while(t<Max_time_iteration)
      {
        std::cout<<"in while loop solver"<<std::endl;
	
              
	 ExchangeData(new_comm,PHI_old);
         //ParticleManager(new_comm); 
	 //ExchangeData(new_comm,P2);

	 Coupling(particles_list,new_comm);
	 ExchangeData(new_comm,P2);
         SetCP4( new_comm);
	 
	 if((int)t%10==0)
           {
            

           // ComputeForce( new_comm);
	   //to compute force
	    //--prepare the grid ->use ArrangeParticles_Cells
	     ArrangeParticles_Cells(Grid, nc1,particles_list,L,new_comm);
	    //--Exchange boundary cells-> Use SendGridPoints_LC
	     SendGridPoints_LC(Grid,nc1, new_comm);

	    //--Compute Force in Lc
	    ComputeForce_LC(Grid,nc1);

	  // After remove colloids that donot belong to the process->use SendColloids
	     SendColloids(particles_list, new_comm);
	   }
	 
	  //ExchangeData(new_comm,P2_local);
	  setLaplacianBase( new_comm);
	  ExchangeData(new_comm,gamma);
	  SetSecondLaplacian2(new_comm);
	 
	  //ExchangeData(new_comm,Laplacian2);
	  FiniteDifferenceScheme(new_comm);
	  UpdateSolution(new_comm);
	 
         
      t+=1.0;
      count++;
      printf("time=%lf, count=%d\n",t,count);
      }
   WriteToFile_MPI( new_comm);
 
 FreeCommunicator(new_comm);

 }




/*
void BASE::SetEffectiveRadius(MPI_Comm new_comm)
{
  // printf("IN EFFECTIVE RAD FUNCTION");
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   

  for( auto i=N_start;i<N_start + nlocal_particles-1;i++)//(int)body->size()-1;i++)
      {
	for(size_t j=i+1; j<bd.size();j++) //Nparticles becos of verlet list
	   {
	   // if(((*body)[i].index==0)&&((*body)[j].index==0))
	     // {
	       //R12[i][j]=RR1small*2.0;

	      //}
	    if((bd[i].index==1)&&((bd[j].index==1)))
	      {
	       R12[i][j]=RR1big*2.0;
	      }
	   // if(((*body)[i].index==1)&&((*body)[j].index==0))
	    //{
	      //R12[i][j]=RR1big + RR1small;
	    //}
	    //if(((*body)[i].index==0)&&((*body)[j].index==1))
	    //{
	      //R12[i][j]=RR1big + RR1small;

	    //}
	 }
     }

}
*/
void BASE:: InsertList( ParticleList **root_list, ParticleList *i )
{

   i->next = *root_list;
   *root_list=i;
}


void BASE::ParticleManager(MPI_Comm new_comm)
{
   int coords[2];// send_offset, recv_offset;	
   //int dims[2];, pred, succ;
   //MPI_Status status;
   //MPI_Request send_request, recv_request;
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   
 
   
    
   for (int i = 0; i < size; i++ )
     { 
        recvcounts[i] = (( i < (Nparticles%size) )? nlocal_particles + 1: nlocal_particles); 
	
        displacement[i] = ((int)(i*nlocal_particles) + ((i < Nparticles%size ) ? i : (Nparticles%size)));       
     }

   /*succ =(my2drank+1)%size;
   pred =(my2drank-1+size)%size;
   send_offset=N_start;//my2drank*nlocal_particles;
   recv_offset=((my2drank-1+size)%size)*recvcounts[pred];
   for(int i=0;i<size-1; i++)
      {
       MPI_Isend(&bd[send_offset], nlocal_particles, particletype, succ, 0, MPI_COMM_WORLD, &send_request);
       MPI_Irecv(&bd[displacement[pred]], recvcounts[pred], particletype, pred, 0, MPI_COMM_WORLD, &recv_request);
       send_offset=((my2drank-i-1+size)%size)*nlocal_particles;
       recv_offset=((my2drank-i-2+size)%size)*nlocal_particles;
       MPI_Wait(&send_request, &status);
       MPI_Wait(&recv_request, &status);
 
      }
   */
   
   
        //MPI_Comm_rank(new_comm,&my2drank);
       // MPI_Cart_coords(new_comm,my2drank,2,coords);
        //MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr ); //move along x direction

        //MPI_Isend(&bd[N_start], nlocal_particles,particletype, right_nbr, 0, new_comm, &send_request);
        //MPI_Irecv(&bd[displacement[left_nbr]],recvcounts[left_nbr],particletype, left_nbr, 0,new_comm, &recv_request);
       
       // MPI_Wait(&send_request, &status);
	//MPI_Wait(&recv_request, &status);
       //compute force here
	//MPI_Comm_rank(new_comm,&my2drank);
        //MPI_Cart_coords(new_comm,my2drank,2,coords);
        //MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr );
       
	//MPI_Isend(&bd[N_start], (nlocal_particles+recvcounts[left_nbr]),particletype, down_nbr, 0, new_comm, &send_request);
        //MPI_Irecv(&bd[displacement[left_nbr]],recvcounts[left_nbr],particletype, left_nbr, 0,new_comm, &recv_request);
     



   




   //for(int i=0;i<Virt)
  // if (size>1)
   //{
    // MPI_Allgatherv(bd+N_start, nlocal_particles,particletype, bd, recvcounts, displacement,particletype, new_comm );
 

  // }
   
   //Virtual_list.resize(1);
   
   //std::cout<<"Size On Each Core="<<<<std::endl;
//delete ParList;


 
}



void BASE::GenerateVerlet(MPI_Comm new_comm)
{
    //check this implementation. Check if the total number of particles per process is Nparticles. 
//Nparticles becos, after the call of MPI_AllgatherV, all processes should have Nparticles and generate only a local verlet list.
printf("IN VERLET FUNCTION");
 MPI_Comm_rank(new_comm, &my2drank);

 MPI_Comm_size(new_comm, &size);
    for(int i=N_start; i<N_start+nlocal_particles;i++)
       {
         VerletList_local[i]=0;
         for(int j=i+1;j<Nparticles;j++)
            {
             double dx=bd[i].r[0]-bd[j].r[0];
             double dy=bd[i].r[1]-bd[j].r[1];
             dx=dx-nlocalx*(int((dx/(nlocalx/2)+3)/2)-1);
             dy=dy-nlocaly*(int((dy/(nlocaly/2)+3)/2)-1);
             double R2=sqrt(pow(dx,2)+pow(dy,2));
             if(R2<R12[i][j])
               {
                VerletList_local[i]=VerletList_local[i] +1;
                VerletTable_local[i][VerletList_local[i]]=j;
	        Dx[i][j]=dx; //distance between body i and body j
	 //       Dx[j][i]=-dx;
	        Dy[i][j]=dy;
	   //     Dy[j][i]=-dy;
	        RI[i][j]=R2;
	     //   RI[j][i]=R2;
	      }
	    }
      }
	

}
void BASE:: FreeCommunicator(MPI_Comm new_comm)
{
         MPI_Comm_free(&new_comm);
		    

}

void BASE::SendColloids(std::list<BODY> &bodies1, MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   //int processID=0;
   std::vector<int> recv_counter;
   recv_counter.resize(size);
   std::vector<int> counter;
   counter.resize(size);
   //int number_colloids= (int) bodies.size();
   //std::memset(counter, 0, sizeof(int)*size);
 //  counter = new int[size];
   for(int i=0;i<size;i++)
      {
	//counter[i]=0;
	//printf("Rank=%d, counter[%d]=%d\n", my2drank,counter[i]);
	std::cout<<"Counter="<<counter[i]<<std::endl;
      }
   //std::memset(recv_counter, 0, sizeof(int)*size);
   std::vector<BODY> Bodies;
   int length= (int) bodies1.size();
   printf("Initial list size in SendColloid=%d,  Rank=%d\n",length, my2drank);
   Bodies.resize(length);
  // std::copy(bodies1.begin(),bodies1.end(), back_inserter(bodies));
    for(std::list<BODY>::iterator it=bodies1.begin();it!=bodies1.end();++it)
       {
         Bodies.push_back(*it);
       }

   
   int tag=100;
   
   //int number;
   MPI_Request send_request, recv_request;
   MPI_Status  status;
   std::vector<std::list<BODY> > SendList;
   SendList.resize(size);
   // work with list instead of vectors 
   std::list<BODY>::iterator itr = bodies1.begin();
   //for(std::list<BODY>::iterator itr = bodies1.begin();itr != bodies1.end();)
     while(itr != bodies1.end())
      {
         
        //int processID=Get_Process_ID(i,Bodies);
	int processID = ((int)(itr->r[0]/(Nx/Procx)))*Procy + (int)(itr->r[1]/(Ny/ Procy));
	printf("Rank=%d, ProcessID=%d\n", my2drank, processID);
	 if(processID==my2drank)
           {
	     ++itr;	      
	   } 
	 else
	 {
             //printf("Rank=%d, ProcessID=%d\n", my2drank, processID);
             // send_list[processID][counter[processID]]=bodies[i];  //<---Use a  list for this
              SendList[processID].push_back(*itr);
	     
	      counter[processID]++;
	      /*--delete colloid from current process---*/
              itr = bodies1.erase(itr);  //<--returns the next element
	      //Bodies.erase(Bodies.begin() + i);
	      //ResizeArray(bodies,i,  new_comm);

	 }
     }
   //copy the elements of the list into an array or vector, for it to be used as a parameter in MPI call
 // std::cout<<"I am Here"<<std::endl;
  BODY **send_list=new BODY*[size];
   for(int i=0;i<size;i++)
     {
       //printf("Rank=%d, counter[%d]=%d\n", my2drank,counter[i]);
       int amount=counter[i];
       printf("Amount=%d\n",amount);
       if(amount>=1)
       {     
	
	//printf("Amount=%d\n",amount);
        send_list[i]  =new BODY[amount];
       }
       else
	  continue;
     }
   // copy the elements from the list to the array
   
   for(int i=0;i<size;i++)
      {
       int j=0;
       for(std::list<BODY>::iterator it=SendList[i].begin();it!=SendList[i].end();++it)
	  {
	   
              send_list[i][j] = *it;
	      j++;
	  }
       printf("j=%d\n",j);
      }

    
   for(int cores=0;cores<size;cores++)
      {
	if(cores!=my2drank)
	  {
           MPI_Isend(&counter[cores], 1,MPI_INT,cores,tag, new_comm, &send_request);
	  // MPI_Isend(send_list[cores],counter[cores],particletype,cores,tag, new_comm, &send_request);
	   MPI_Recv(&recv_counter[cores], 1,MPI_INT,cores,tag,new_comm,&status);
	   printf("Rank=%d, Recv_counter[%d]=%d\n",my2drank,  cores, recv_counter[cores]);

           MPI_Wait(&send_request, &status);
           //MPI_Wait(&recv_request, &status);

         }
      }

   BODY **recv_list=new BODY*[size];
   for(int i=0;i<size;i++)
     {
       int amount=recv_counter[i];
       //printf("Amount=%d\n",amount);
       if(amount>=1)
        {
          recv_list[i]  = new BODY[amount];
	}
       else
	  continue;
     }
   
     
     for(int cores=0;cores<size;cores++)
      {
	if(cores!=my2drank)
	  {

           MPI_Isend(send_list[cores],counter[cores],particletype,cores,tag, new_comm, &send_request);
           MPI_Recv(recv_list[cores],recv_counter[cores],particletype,cores,tag, new_comm, &status);
	   //printf("Rank=%d, Recv_counter[%d]=%d\n",my2drank,  cores, recv_counter[cores]);
           MPI_Wait(&send_request, &status);
           //MPI_Wait(&recv_request, &status);
	  }
      }

 for(int cores=0;cores<size;cores++)
      {

	   if(counter[cores]!=0)
             {
               for(int recvelts=0;recvelts<recv_counter[cores];recvelts++)
		  {
	            Bodies.push_back(recv_list[cores][recvelts]);		  
		   }
	     }
     }
      

   
//copy the elements from the vector back to the list. The list is now modified
// First we clear the list to avoid a duplicate of an element since we are using back_inserter

//bodies1.clear();
//std::copy( Bodies.begin(), Bodies.end(), std::back_inserter( bodies1 ) );
//int length2= (int)bodies1.size();
//printf("final  list size in SendC=%d\n",length2);

for(int i=0;i<size;i++)
{
       //printf("Rank=%d, counter[%d]=%d\n", my2drank,counter[i]);
       int amount=counter[i];
       //printf("Amount=%d\n",amount);
       if(amount >=1)
       {      
         delete [] send_list[i];
       }
       int amount1=recv_counter[i];
       if(amount1>=1)
	{
          delete [] recv_list[i];
	}
//   delete [] send_buffer[i];
   }
 delete [] send_list;
 delete [] recv_list;
 
//delete [] counter;
}

void BASE::ComputeForce_LC(std::vector<std::vector<std::list<BODY> >  > &grid, int *nc)
{
  int ic[2], nb_c[2], my_cell_index;
  double G4(0.0), AA(0.0), BB(0.0), NEW1(0.0), NEW2(0.0), G2(0.0), G3(0.0), G3A(0.0),C(0.0);

  //scan inner cells
  for(ic[0]=1;ic[0]<nc[0]-1; ic[0]++)
     {
      for(ic[1]=1;ic[1]<nc[1]-1;ic[1]++)
         {
          //compute scalar cell index
          //my_cell_index= ic[0]*nc[1] + ic[1];
          for(std::list<BODY>::iterator List_it=grid[ic[0]][ic[1]].begin();List_it!=grid[ic[0]][ic[1]].end();List_it++)
	     {
                // scan neighbors
	       //scan the neighbor cells
	        for(nb_c[0]=ic[0]-1;nb_c[0]<ic[0]+1;nb_c[0]++)
	           {
		    for(nb_c[1]=ic[1]-1;nb_c[1]<ic[1]+1;nb_c[1]++)
		       {
		        for(std::list<BODY>::iterator neighbor_List=grid[nb_c[0]][nb_c[1]].begin();neighbor_List!=grid[nb_c[0]][nb_c[1]].end(); neighbor_List++)
		           {
		             if(List_it!=neighbor_List)
			       {
		                 double dx=List_it->r[0]-neighbor_List->r[0];
		    		 double dy=List_it->r[1]-neighbor_List->r[1];
				 double R2=sqrt(pow(dx,2)+pow(dy,2));
				 if(R2<2.0*R1big)
				   {
				    //Dx[i][j]=dx;
				   // Dy[i][j]=dy;
				    //RI[i][j]=R2;
				    //RI[j][i]=R2;
				    C = U1/(2.0*R1big);
				    AA=R2/(2.0*R1big);
				    G2 = pow(AA + beta,2.0);
				    G3 = 1.0 + ALPHA*(AA + beta);
				    BB = (AA-1.0)*(-ALPHA);
				    G3A = exp(BB);
				    NEW1=pow((1.0 + (ENNE/ALPHA)+ beta),2.0);
				    NEW2=(exp(-ENNE))*((1.0 + ALPHA*(beta + (ENNE/ALPHA)+1.0)))/(NEW1);
				    G4 = ((C*(G3A)*G3)/G2)-(C*NEW2);
				    List_it->Force[0]=List_it->Force[0] + G4*(dx/R2);
				    List_it->Force[1]=List_it->Force[1] + G4*(dx/R2);
                                }
			  }
		  }
	    }
	 }
      }
    }
  }

//Compute random force

 double csi1(0.0), csi2(0.0);
 srand(time(NULL) + my2drank);
 double GAMM(0.0);
 for(ic[0]=1;ic[0]<nc[0]-1; ic[0]++)
    {
     for(ic[1]=1;ic[1]<nc[1]-1;ic[1]++)
        {
          for(std::list<BODY>::iterator List_it=grid[ic[0]][ic[1]].begin();List_it!=grid[ic[0]][ic[1]].end();List_it++)
  	     {
     	       double range =Random_max-Random_min;
               double div =RAND_MAX / range;
               csi1 = Random_min + (rand()/div);
               csi2 = Random_min + (rand()/div);
               //  printf("Random number (csi1=%lf, csi2=%lf) \n", csi1, csi2);
               if(List_it->index==0)
                 {
         	  GAMM= GAMMs;
		  List_it->v[0]= (1.0/GAMM)*(List_it->Force[0]+ List_it->TOT[0])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
		  List_it->v[1]= (1.0/GAMM)*(List_it->Force[1]+ List_it->TOT[1])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
		 // printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);
		 }
	       if(List_it->index==1)
	         {
	          GAMM= GAMMb;
	          List_it->v[0]= (1.0/GAMM)*(List_it->Force[0]+ List_it->TOT[0])*delta_t + sqrt(2.0*temp1/GAMM)*csi1* sqrt(delta_t);
		  List_it->v[1]= (1.0/GAMM)*(List_it->Force[1]+ List_it->TOT[1])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*  sqrt(delta_t);
	 //printf("Colloid force (fx=%lf, fy=%lf) \n", bd[i].Force[0], bd[i].Force[1]);
	 }
	 //position update
         List_it->r[0] =List_it->r[0] + List_it->v[0]; //<---bd[i].vx is already multiplied by t. 
	 List_it->r[1] =List_it->r[1] + List_it->v[1];
	 /*---------------------------------------------
	 #############################################################################################################
         The boundary conditions on the colloids is implicitly applied here. That is, 
	 in the SendColloid function, the function scans the list of colloids
	 belonging to each process. If a colloid does not belong that that process (due to position update),
	 it sends it to the appropriate process.

         ######################################################################################################
	 */

	 //bc condition
	// bd[i].r[0]=bd[i].r[0]-Nx*( (int)(bd[i].r[0]/Nx + 1)-1);
	 // bd[i].r[1]=bd[i].r[1]-Ny*( (int)(bd[i].r[1]/Nx + 1)-1);
        //if(my2drank==0)	  
       //printf("Force_x[%d]=%lf, Force_y[%d]=%lf\n",i, bd[i].Force[0],i,bd[i].Force[1]);
       //computing mean displacement
     //  Posx[i]=Posx[i] + List_it->v[0];
     /*  Posy[i]=Posy[i] + List_it->v[1];
       DSUMX  = DSUMX + pow(Posx[i],2.0);
       DSUMY  = DSUMY + pow(Posy[i],2.0);
       DSUMXY = DSUMXY + Posx[i]*Posy[i];*/ 
       //printf("Colloid positions (bd[100].x=%lf, bd[100].y=%lf) \n", bd[100].r[0], bd[100].r[1]);
  }
										      // DSUMX=DSUMX/


	}
    }

}

void BASE::SendGridPoints_LC(std::vector<std::vector<std::list<BODY> > > &mygrid,const int *nc,MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   int coords[2], left_nbr(0), right_nbr(0), up_nbr(0), down_nbr(0), left_count(0),right_count(0);
   std::list<BODY>::iterator it;
   int tag=100, k;
   int recv_count=0;

   
      MPI_Request send_request, recv_request,recv_request1;
      MPI_Status  status;
  //-----Beginning of sending to left neighbour and receiving from right neighbour------//
   
   //----Count number of colloids to send to neighbours--------//
   for(auto i=1;i<nc[0]-1;i++)
      {
       int j=1;
       
       left_count +=(int)mygrid[i][j].size();
          
      }
      int ToSend=left_count;
     for(auto i=1;i<nc[0]-1;i++)
      {
       int j=nc[1]-2; //last col
       
       right_count +=(int)mygrid[i][j].size();
          
      }
     int ToSend_right=right_count;
  //send data to  left and receive from right

  MPI_Cart_shift(new_comm, 1,1, &right_nbr, &left_nbr);  
  MPI_Isend(&ToSend, 1, MPI_INT, left_nbr,  tag,new_comm,  &send_request);
  MPI_Recv(&recv_count, 1, MPI_INT, right_nbr,tag, new_comm, &status);
  MPI_Wait(&send_request,&status);
  //MPI_Wait(&recv_request,&status);
 
  //copy the grid portion to be sent
  std::vector<BODY> Send_particles_left;
  Send_particles_left.resize(left_count*nc[0]);
  std::vector<BODY> Receive_particles_right;
  Receive_particles_right.resize(nc[0]*recv_count);
   k=0;
 
  for(int i=1;i<nc[0]-1;i++)
     {
        int j=1;
        
        for(std::list<BODY>::iterator it=mygrid[i][j].begin();it!=mygrid[i][j].end();it++)
           {
             Send_particles_left[ k]=*it;
             k++;
           }
         

      }
 MPI_Isend(&Send_particles_left[0], ToSend,particletype,left_nbr, tag, new_comm, &send_request);
 MPI_Recv(&Receive_particles_right[0],recv_count, particletype,right_nbr,tag, new_comm, &status);
 MPI_Wait(&send_request,&status);                  
 //MPI_Wait(&recv_request1,&status); 
 //int INDEX=0;
 int partition_cells_right=(int)recv_count/nc[0];
 for(int i=1;i<nc[0]-1;i++)
    {
      int j= nc[1]-1;
      
      
      for(int k=0;k<partition_cells_right;k++)
         {
           mygrid[i][j].push_back(Receive_particles_right[i*partition_cells_right + k]);
       
         }


    } 
 if(recv_count%nc[0]!=0) //if the number of particles in the recv array is not divisible by the grid number nc[1]  
   {
      //if the remainder is more than one particles, add the rest of the particles to the last grid list
	
	  auto  amount=recv_count%nc[0];
	  for(auto i=0;i<amount;i++ )
	     {
               mygrid[(nc[0]-2)][nc[1]-1].push_back(Receive_particles_right[(nc[0])*(partition_cells_right)+i]);
             }
     
   }

   
 //---------send to right neighbour and receive from left neighbour---------
  tag=200;
  int recv_count_left=0;
  MPI_Cart_shift(new_comm, 1,1, &left_nbr, &right_nbr);
  MPI_Isend(&ToSend_right, 1, MPI_INT, right_nbr, tag,new_comm,  &send_request);
  MPI_Recv(&recv_count_left, 1, MPI_INT, left_nbr,tag, new_comm, &status);
  MPI_Wait(&send_request,&status);
  //MPI_Wait(&recv_request,&status);
  
  std::vector<BODY> Receive_particles_left;
  Receive_particles_left.resize(nc[0]*recv_count_left);
  std::vector<BODY> Send_particles_right;
  Send_particles_right.resize(right_count*nc[0]);
  k=0;
   for(int i=1;i<nc[0]-1;i++)
     {
      
        int j=nc[0]-2;
        
        for(std::list<BODY>::iterator it=mygrid[i][j].begin();it!=mygrid[i][j].end();++it)
           {
              Send_particles_right[ k]=*it;
              k++;
           }
         

      }
 MPI_Isend(&Send_particles_right[0], ToSend_right,particletype,right_nbr, tag+1, new_comm, &send_request);
 MPI_Recv(&Receive_particles_left[0],recv_count_left, particletype,left_nbr,tag+1, new_comm, &status);
 MPI_Wait(&send_request,&status);
 //MPI_Wait(&recv_request1,&status);

 int partition_cells_left=(int)recv_count_left/nc[0];
 for(int i=1;i<nc[0]-1;i++)
    {
      int j= 0;
      
      
      for(int k=0;k<partition_cells_left;k++)
         {
           mygrid[i][j].push_back(Receive_particles_left[i*partition_cells_left + k]);
       
         }

    }

 if(recv_count_left%nc[0]!=0) //if the number of particles in the recv array is not divisible by the grid number nc[1]  
   {
      //if the remainder is more than one particles, add the rest of the particles to the last grid list
	
	  auto  amount=recv_count_left%nc[0];
	  for(auto i=0;i<amount;i++ )
	     {
               mygrid[(nc[0]-2)][0].push_back(Receive_particles_left[(nc[0])*(partition_cells_left)+i]);
             }
     
   }





 /*
 //---Add particles to ghost grid cells cells before communicating up and down
// use ArrangeParticles_cells function. But we need a list
int size_k = (int) (recv_count_left);
INDEX =0;
 for(int i=1;i<nc[0]-1;i++)
    {
      int j= 0;
      int I=i*nc[1] + j;
      for(int k=0;k<size_k;k++)
         {
           mygrid[I].push_back(Receive_particles_left[(INDEX)*size_k + k]);
       
         }
       if(INDEX<nc[0])
         {
           INDEX++;
         }

    } 
*/

//------END OF LEFT and RIGHT COMMUNICATION OF PARTICLES-----------------//
// ######################################################################
// -------Begin Up and down communication--------------------------------//
//
int up_count(0), down_count(0);
for(auto j=0;j<nc[1];j++)
   {
       int i=1;
       
       up_count +=(int)mygrid[i][j].size();
          
   }
      
     for(auto j=0;j<nc[1];j++)
      {
       int i=nc[0]-2; //last col
       
       down_count +=(int)mygrid[i][j].size();
          
      }
     
  tag=301;
  int recv_count_up(0);
//--------Sending to down and receiving from up--------------------------------//
  MPI_Cart_shift(new_comm, 0,1, &up_nbr, &down_nbr);  
  MPI_Isend(&down_count, 1, MPI_INT, down_nbr,  tag,new_comm,  &send_request);
  MPI_Recv(&recv_count_up, 1, MPI_INT, up_nbr,tag, new_comm, &status);
  MPI_Wait(&send_request,&status);
  //MPI_Wait(&recv_request,&status);
// printf("Rank=%d, send=%d, recv=%d: to=%d, From=%d\n", my2drank, down_count,recv_count_up, down_nbr,up_nbr);
  std::vector<BODY> Send_down;
  Send_down.resize(down_count*nc[1]);
  std::vector<BODY> Receive_particles_up;
  Receive_particles_up.resize(recv_count_up*nc[1]);
   k=0;  
   for(int j=0;j<nc[1];j++)
     {
       // int k=0;
        int i=nc[0]-2;
        
        
        for(std::list<BODY>::iterator it=mygrid[i][j].begin();it!=mygrid[i][j].end();++it)
           {
              Send_down[ k]=*it;
              k++;
           }
         

      }
 MPI_Isend(&Send_down[0], down_count,particletype,down_nbr, tag+1, new_comm, &send_request);
 MPI_Recv(&Receive_particles_up[0],recv_count_up, particletype,up_nbr,tag+1, new_comm, &status);
 MPI_Wait(&send_request,&status);
 //MPI_Wait(&recv_request1,&status);

//------Insert particles into upper ghost cells of LC grid--------------------
 int partition_cells=(int)(recv_count_up/nc[1]); //<-- ie divide the particles into cells
 for(int j=0;j<nc[1];j++)
    {
      int i= 0;
      
      for(int k=0;k<partition_cells;k++)
         {
          mygrid[i][j].push_back(Receive_particles_up[(j)*partition_cells + k]);
         }
    } 

 if(recv_count_up%nc[1]!=0) //if the number of particles in the recv array is not divisible by the grid number nc[1]  
   {
      //if the remainder is more than one particles, add the rest of the particles to the last grid list
	
	  auto  amount=recv_count_up%nc[1];
	  for(auto i=0;i<amount;i++ )
	     {
               mygrid[0][nc[1]-1].push_back(Receive_particles_up[(nc[1])*(partition_cells)+i]);
             }
     
   }
//-------sending to up neighbour and receiving from down-------------------
//
// ###### check here -----------
  tag=tag+19;
  int recv_count_down1(0);
  MPI_Cart_shift(new_comm, 0,1, &down_nbr, &up_nbr);  
  MPI_Isend(&up_count, 1, MPI_INT, down_nbr,  tag,new_comm,  &send_request);
  MPI_Recv(&recv_count_down1, 1, MPI_INT, down_nbr,tag, new_comm, &status);
  MPI_Wait(&send_request,&status);
  //MPI_Wait(&recv_request,&status);
 

  std::vector<BODY> Send_up;
  Send_up.resize(up_count*nc[1]);
  std::vector<BODY> Receive_particles_down;
  Receive_particles_down.resize(recv_count_down1*nc[1]);
   k=0; 
   for(int j=0;j<nc[1];j++)
     {
        
        int i=1;
        
        for(std::list<BODY>::iterator it=mygrid[i][j].begin();it!=mygrid[i][j].end();++it)
           {
              Send_up[ k]=*it;
              k++;
           }
     }

 MPI_Isend(&Send_up[0], up_count,particletype,up_nbr, tag+12, new_comm, &send_request);
 MPI_Recv(&Receive_particles_down[0],recv_count_down1, particletype,down_nbr,tag+12, new_comm, &status);
 MPI_Wait(&send_request,&status);
 //MPI_Wait(&recv_request1,&status);
//------Insert particles into Lower  ghost cells of LC grid--------------------
 int partition_cells1=(int)(recv_count_down1/nc[1]); //<-- ie divide the particles into cells
 for(int j=0;j<nc[1];j++)
    {
      int i= nc[0]-1;
      
      for(int k=0;k<partition_cells1;k++)
         {
           mygrid[i][j].push_back(Receive_particles_down[(j)*partition_cells1 + k]);
         }
    }
 if(recv_count_down1%nc[1]!=0) //if the number of particles in the recv array is not divisible by the grid number nc[1]  
   {
      //if the remainder is more than one particles, add the rest of the particles to the last grid list
	
	  auto  amount=recv_count_down1%nc[1];
	  for(auto i=0;i<amount;i++ )
	     {
               mygrid[nc[0]-1][nc[1]-1].push_back(Receive_particles_down[(nc[1])*(partition_cells1)+i]);
             }
     
   }


}



void BASE::ResizeArray(std::vector<BODY> &bodies,int position,  MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   if( position >=nlocal_particles )
     printf("Deletion not possible.\n");
   else
      {
       for(int c = position-1; c<nlocal_particles-1;c++)
	  {
	    bodies[c] = bodies[c+1];
	  }
     }
       
} 

void BASE::InsertColloid(BODY *bodies,int postion, MPI_Comm new_comm)
{
  MPI_Comm_rank(new_comm, &my2drank);
  MPI_Comm_size(new_comm, &size);






}



void BASE::Coupling(std::list<BODY>& bodies_list, MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   int coords[2], dims[2],periods[2];
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   MPI_Cart_get(new_comm, 2,dims, periods,coords);
  // std::vector<auto> integration_index; //<--stores indices for integration
  char fname[200],fname2[200];
  sprintf(fname, "coupling%d.dat",my2drank );
  sprintf(fname2, "gridpoints%d.dat",my2drank );

  FILE *fp,*fp2;
  fp=fopen(fname, "w");
  fp2=fopen(fname2, "w");

  //Extract the particles from the list into bd vector

   bd.resize(bodies_list.size());
   //std::copy(std::begin(bodies_list),std::end(bodies_list), std::back_inserter(bd));
   for(std::list<BODY>::iterator it=bodies_list.begin();it!=bodies_list.end();++it)
      {
        bd.push_back(*it);
      }



   int N1(0), X(0), Y(0);
   double RX1(0.0),RY1(0.0), R2(0.0), B1(0.0),B2(0.0), D1(0.0), D2(0.0),GT1(0.0),GT2(0.0);
   //int counter=0;
   for(int i=0;i<nlocalx+6;i++)
     {
      for(int j=0;j<nlocaly+6;j++)
         {
            PP3(i,j)=0.0;
	  }
    }
   
  //for(int i=0;i<number_particles[my2drank];i++)
    for(int i=0; i<(int)bd.size();i++)
     {
       //int i =N_start + index_particle; 
      if (bd[i].index==1 )
       {

	 /*here, we should consider only particles that are in the current process. We start by checking
	  weather the particles are within the MPI process.  */
	  double dx=0.5; //<--grid spacing
          P00=P00b;
          //N1=floor((RR1big/dx)+2);
	  N1= (int)(RR1big)+1;
	  bd[i].TOT[0]=0.0;
	  bd[i].TOT[1]=0.0;
	  const int OFFSETX = coords[0]*(nlocalx);
	  const int OFFSETY = coords[1]*(nlocaly);
	  double local_coordsx= bd[i].r[0]-(OFFSETX);
	  //int N_grid_x= (int)local_coordsx/dx;
	  double local_coordsy= bd[i].r[1]-OFFSETY;
	  //int N_grid_y=(int)local_coordsy/dx;
	  //double min_x = local_coordsx - N1;
	  //double max_x = local_coordsx + N1;
	  //int min_y = (int)local_coordsy-N1;
	  //int max_y = (int)local_coordsy+N1;
	 /* if(local_coordsx<RR1big)     //<---check this, maybe R1big instead of RR1big
	    {
             min_x=-N1;
	    }
	 else
	    {
	      min_x=local_coordsx-N1;	    
	    }
	 if(local_coordsx+RR1big>nlocalx)
	   {
	     max_x=N1;

	   }
	 else
	   {
             max_x= local_coordsx+N1;
	   }
	  int MIN_r =(int)min_x;
	  int MAX_r =(int)max_x;
    //REDO THIS LOOP
    */     //std::cout<<"offset="<<coords[0]<<std::endl;
           for(int j=-N1;j<N1;j++)
	      {  
		 
		 // std::cout<<"offset="<<OFFSETX<<std::endl;
		 // std::cout<<"offsety="<<OFFSETY<<std::endl;

		  GT1=bd[i].r[0]-OFFSETX+j;
           if(GT1>(double)(nlocalx+2))
	           {
	           // printf("Rank=%d, GT>nlocalx:: GT=(%lf,%lf ), J_index=%d\n ",my2drank,bodies[i].r[0]-OFFSETX , GT1,j);
                    //GT1=GT1-((Nx))*((int)(GT1/((Nx)) +1)-1);
                 GT1=GT1-3.0;
	           }
	        if(GT1<3.0)
	          {
                    GT1= GT1+3.0;
	          }
             
	         X = (int )(GT1);  //<--- global index
		 //X =X-OFFSETX+3;   //<--  local index
               /*  if(X<0)
                   {
                     //int indices=X-(nlocalx+2);
                     X=X+3;//indices+nlocalx; //<---how many indices into the halo region
                   }
                 if(X>nlocalx+2)
                   {
                    X=X-3;
                   }
		*/



     for(int k=-N1;k<N1;k++)
        {
                    
            GT2= bd[i].r[1]-OFFSETY +k;
		   //printf("Rank=%d, GT=(%lf,%lf)\n ",my2drank, bodies[i].r[0]-OFFSETX,bodies[i].r[1]-OFFSETY);

				    
       
           if(GT2>(double)(nlocaly+2))
	     {
           //GT2=GT2-(Ny)*((int)(GT2/((Ny)) +1)-1);
               GT2=GT2-3.0;
	     }
           if(GT2<3.0)
            {
              GT2= GT2+3.0;
            }
               
	   Y=(int)GT2;
		    //Y = Y-OFFSETY+3; // +3 becos our grid starts from 3 (0,1,2 are halo regions)
                /*    if(Y<0)
                      {
                         Y = Y+3;
                      }
                    if(Y>nlocaly+2)
                      {
                       Y=Y-3;
                      }
                */

                    
       // printf("nlocalx=%d, nlocaly=%d, X=%d, Y=%d\n", nlocalx, nlocaly,X,Y); 
		     
                  // increase the halo regions and check again		    
		  
          RX1=bd[i].r[0]-OFFSETX-GT1;
		   // printf(" Before :: RX1=%lf\n", RX1);

		    //RX1=RX1-((Nx-1)*((int)((RX1/((Nx-1)/2)+3)/2)-1));
          RY1=bd[i].r[1]-OFFSETY -GT2;
		    //printf(" Before :: RY1=%lf\n", RY1);

		   // RY1=RY1-(Ny-1)*(floor((RY1/((Ny-1)/2)+3)/2)-1);
                    //printf("After :: RX1=%lf, RY1=%lf\n", RX1, RY1);
                    		  
	 R2=sqrt(pow(RX1,2)+pow(RY1,2));
		    //R2=(double)sqrt(pow(j-local_coordsx,2)+pow(k-local_coordsy,2));
		   // printf("R2=%lf,RR1big=%lf\n",R2,RR1big);
                    //printf("Rank=%d, GT1=%lf, GT2=%lf, bodies.r[0]=%lf, bodies.r[1]=%lf\n", my2drank, GT1,GT2, bodies[i].r[0], bodies[i].r[1]);


	 if(R2<RR1big)
	   {
	               //Dx1[j][k]=RX1;
	               //Dy1[j][k]=RY1;
	       B1 = 1.0-pow((R2/RR1big), ALP1);
	       B2 =exp(1.0-(1.0/B1));  //<----- Phi_c,i equation 16
	       D1 =ALP1/pow(RR1big,ALP1);
	       D2 =-(1.0)*(D1*pow(R2,ALP1-2.0)/(pow(B1,2)))*B2;//<--check this
		       //POS[j][k]=i;
	  //     printf("Rank=%d:: X=%d, Y=%d\n", my2drank,X,Y);
	       fprintf(fp2, "%d  %lf   %lf  %d  %d  %d %d\n",my2drank,local_coordsy,local_coordsx, X, Y, j ,k );

                //printf("Rank=%d, X=%d,Y=%d, B2=%lf) \n",my2drank,X, Y, B2*(PHI_old(X,Y)-P00));
               // double PP3_aux=B2*(PHI_old(j,k)-P00);
	      PP3(X,Y)=  PP3(X,Y) + B2*(P3(X,Y)-P00);

		       
		       
	      
	      if(R2<R1big)
	       {
	         PPP(X,Y)=1.0;
		
               }
	               
		       
	     bd[i].TOT[0]+= -sigma*RX1*D2*pow((P3(X,Y)-P00),2);
	     bd[i].TOT[1]+= -sigma*RY1*D2*pow((P3(X,Y)-P00),2);
	     fprintf(fp, "%d  %lf   %lf\n",i, bd[i].TOT[0], bd[i].TOT[1]);
             //printf("Rank=%d,  counter[%d]=%d) \n",my2drank, my2drank,counter[my2drank]);
	//     printf("Rank=%d,R1big=%lf,  (X =%d, Y=%d,PP3[%d][%d]=%lf), Phi_old=%lf \n",my2drank,RR1big,X, Y,X,Y, PP3(X,Y),PHI_old(X,Y) );
		      // printf("Colloid (X =%d, Y=%d,TOTX=%lf,TOTY=%lf) \n",X, Y, bd[i].TOT[0],bd[i].TOT[1]);

	   }
		   // printf("Rank=%d,  counter[%d]=%d) \n",my2drank, my2drank,counterx[my2drank]);
         }
	     		
      }
	     
             

    }
 } //end of particle loop
       fclose(fp);
       fclose(fp2);
   //computation equation 29
        
   for(int i=3;i<nlocalx+3;i++)
      {
	
       for(int j=3;j<nlocaly+3;j++)
	  {
           if(POS[i][j]!=0)
            {
             int z= POS[i][j];
	      if(bd[z].index==0)
		 {
                  P00=P00s;
		 }
	      else if(bd[z].index==1)
		 {
                   P00=P00b;
		 }
	      else
		 {
	          continue;
		 }
	    }

	   //CP4(i,j)=0.0;
	   P2(i,j)=0.0;
           P2(i,j)=sigma*2.0*PP3(i ,j );
	  // printf("Rank=%d, PP3(%d,%d)=%lf\n",my2drank, i,j,PP3(i,j));
	   //printf("PP3[%d][%d]=%lf\n",Global_Index_X[countx] ,Global_Index_Y[county], PP3(Global_Index_X[countx],Global_Index_Y[county]));
	  }
      
   
      }
/*
 if(my2drank==size-1)
    {	  
     for(int i=0;i<nlocalx;i++)
        {
	 int x = i + coords[0]*nlocalx;

    	 for(int j=0;j<nlocaly;j++)
	    {
	       int y = j + coords[1]*nlocaly;

	      printf("PP3(%d,%d)=%lf\n",i,j,PP3(x,y) );
	  }
      }
    }*/
 
 
   //calculate laplacian for use with polymer

   //FIX THIS PLACE. THE PROBLEM IS ABOUT THE BOUNDARY CONDITION ON P2

 /*  int glob_index_x=0;
   int glob_index_y=0;
   double AP3(0.0), BP3(0.0);

   for(int i=3;i<=nlocalx+2;i++)
      {
       for(int j=3;j<=nlocaly+2;j++)
	  {

                AP3=(1.0/6.0)*(P2(i+1,glob_index_y )+P2(DOWNX[glob_index_x],glob_index_y ) + P2(glob_index_x ,UPX[glob_index_y]) + P2(glob_index_x ,DOWNX[glob_index_y]));
		BP3=(1.0/12.0)*(P2(DOWNX[glob_index_x],DOWNX[glob_index_y]) + P2(DOWNX[glob_index_x],UPX[glob_index_y]) + P2(UPX[glob_index_x], DOWNX[glob_index_y]) + P2(UPX[glob_index_x],UPX[glob_index_y]));
	//	CP4(i,j)=AP3 + BP3-P2(glob_index_x,glob_index_y);
	      //}
             //CP4(i,j)=0.0;//AP3 + BP3-P2(glob_index_x,glob_index_y);
	    //printf("UPX[%d]=%d, DOWNX[%d]=%d\n",glob_index_x, glob_index_x, UPX[glob_index_x],DOWNX[glob_index_x]);
	  }

      }*/

//bodies_list.clear();
//std::copy(bd.begin(),bd.end(), std::back_inserter(bodies_list));
bd.clear();
}
/*
void BASE::WriteColliodsToFile_MPI(MPI_Comm new_comm)
{
   MPI_Status status;
   MPI_File     fh;
   MPI_Datatype filetype;
   int dims[2],coords[2], periods[2], start_indices[2], localsizes[2];
   int globalsizes[2];
   globalsizes[0]=Nx;
   globalsizes[1]=Ny;
   localsizes[0]=(Nx/Procx); localsizes[1]=(Ny/Procy) ;
   int local_array_size =nlocalx*nlocaly;
   MPI_Cart_get (new_comm,2,dims,periods, coords );
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   start_indices[0]=coords[0]*localsizes[0];
   start_indices[1]=coords[1]*localsizes[1];
   printf("rank=%d,r=%d,start_indices=(%d,%d)\n",my2drank,Nx%Procx,start_indices[0],start_indices[1]);
   MPI_Type_create_subarray(2, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);
   MPI_Type_commit(&filetype);
   char outputfilename[]="datafile_mpi";
   char filememview[]="native";
   MPI_File_open(new_comm, outputfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
   MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, filememview,MPI_INFO_NULL);
   MPI_File_write_all(fh, PHI_local_result,local_array_size , MPI_DOUBLE, &status);
   MPI_File_close(&fh);
   MPI_Type_free(&filetype);

}

*/



