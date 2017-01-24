#ifndef _BASECLASS_
#define _BASECLASS_

//#include "CDS_Polymer.h"
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <assert.h>
#include "mpi.h"
#include <memory>
#include <cstddef>
#include <algorithm>
#include <cmath>

#define EMPTY -1
#define NBMAX 100000  //Maximum # of copied boundary atoms per neighbor.

typedef struct{                  
 int index; //determines if body is small or larg ie index=1 corresponds to large colloid and index=0 otherwise 
 double r[2];    
 //double y;     
 double v[2];
 //double vy; 
 double mass;
 double TOT[2]; //<------Coupling with polymer
 double Force[2];
 }BODY;



typedef struct ParticleList{
   BODY bd;
   struct ParticleList *next;
}ParticleList;


typedef struct Node{
   double cx;
   double cy;
   double cmx;
   double cmy;
   double nmass;
   double xmin, xmax;
   double ymin, ymax;
   struct Node* nw;
   struct Node* ne;
   struct Node* sw;
   struct Node* se;
   BODY* bd;
   int has_childs;
 } Node;


//U1 energy scale, kbT
class BASE
{

  protected:
	    int Nx, Ny, Procx, Procy,nlocalx,nlocaly,my2drank,rank,size, Nbig,Nsmall, Nparticles, nlocal_particles, N_start;
	    int nbig_local,nsmall_local,  *recvcounts, *displacement, *recvcounts_polymer, *disp_polymer;;
	    double **PHI, **PHI_old, **Laplacian2, **gamma,**PP3,**P2,**CP4,**P3,  *PHI_local_result, **PPP, *POLYMER,*PPP_local, *P2_local, *PPP_local2;
            int InitUcell[2]={16,16}; // 
	    double al[2]={1.0,1.0};  //<---box length in x and y direction
	    double  RCUT;
	    int Num_cell_x, Num_cell_y;
            double locala_x, localb_x, locala_y, localb_y; //<--- parameters to determine process grid portion [a,b].
	    int lsb[6][NBMAX];
	    int Nglobal_endx, Nglobal_startx, Nglobal_endy, Nglobal_starty; //<-----parameters to determine the global indices of the left and right boundaries of a process. Useful in Coupling member function. 
	    std::vector<int> LinkedCells;  //<--number of linked cells in x and y direction
	    std::vector<double> CellSize;
	    std::vector<int> head;
	    std::vector<int> LinkedCell_list;
	    int **partitionx, **partitiony;
	    int * nc1; // -----Parameter for the linked list-----
	    double * L;  //-------Parameter for the linked list----
	   // BODY *bd;
	   // BODY *bd1;
//	    BODY **send_list;
    //        BODY **recv_list;
	    BODY **send_buffer;
	    std::list<BODY> particles_list;
	    std::vector<BODY> bd;
	    std::vector<BODY> bd1;
	    std::vector<BODY> down_nbr_list;
	    std::vector<BODY> right_nbr_list;
	    std::vector<BODY> up_nbr_list;
	    std::vector<BODY> left_nbr_list;
	    std::vector<BODY> Virtual_list;
	    // use for creating a grid wherein each grid cell is a list
	    //std::vector<std::vector<std::list<BODY> > > Grid;
	    //std::vector<std::list<BODY>  > Grid;
            std::vector<std::vector<std::list< BODY> > > Grid; 
	    Node *nd;
	   // std::vector<BODY> bd;
	    double ALP1, RR1small, R1small, R1big, RR1big, P00, P00s,sigma, P00b;
	    double ENNE, ALPHA, GAMMs, GAMMb, U1, beta, delta_t, temp1, dxx, EMME, CSI, B2 ;

            int MPI_return;

	    double D, M, u, A , v, r, F, tau, Max_time_iteration, B;
	    std::vector<std::vector<double> > R12; //effective radius 
	    std::vector<int >VerletList_local;
	    std::vector<std::vector<int> >VerletTable_local;
	    std::vector<std::vector<double> > Dx;
	    std::vector<std::vector<double> > Dx1;
	    std::vector<std::vector<double> > Dy;
            std::vector<std::vector<double> > Dy1;
	    std::vector<std::vector<double> > RI;
	    //std::vector<std::vector<double> > r;
	  //  std::vector<std::vector<double> > PPP;
	    std::vector<std::vector<int>  > POS;
	    std::vector<double> TOTX;
	    std::vector<double> TOTY;
	    double *matrix_lower;
	    double *matrix_upper;
	    double *matrix_left;
	    double *matrix_right;
            std::vector<int> number_particles;

	    std::vector<double> Posx;
	    std::vector<double> Posy;
	    std::vector<double> Vx;
	    std::vector<double> Vy;
	    std::vector<double> Fx;
	    std::vector<double> Fy;
	    std::vector<int> up1x;
	    std::vector<int> UPX;
	    std::vector<int> DOWNX;
	    std::vector<int> down1x;
	    std::vector<int> Global_Index_X;
	    std::vector<int> Global_Index_Y;
	    const int nitems=5;
	    int  blocklengths[2] ;//= {1,1,1,1,1,1};
	     
	    MPI_Datatype types[2] ;//= {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
            MPI_Datatype particletype;
            MPI_Aint     offsets[2], extent;//={offsetof(BODY, index),offsetof(BODY, x),offsetof(BODY, y),offsetof(BODY, vx),offsetof(BODY, vy),offsetof(BODY, mass)};
	     //MPI_DOUBLE   
	    	     
             //MPI_Type_create_struct(6,blocklengths,offsets,types,&mpi_body_type);
	    // MPI_Type_create_stuct(6,blocklengths,offsets,types,&mpi_body_type);
	    // MPI_Type_commit(&mpi_body_type);
            

  public:
        explicit   BASE();
        virtual   ~BASE();
 
        void Read_input_parameters(int *integers, double *rationals);

        MPI_Comm CreatCartesianTopology();

        void GenerateColloids(std::list<BODY>& bodies_list, MPI_Comm new_comm);

        void DistributeColloids(MPI_Comm new_comm);

        void ParticleManager(MPI_Comm new_comm);

        void GenerateVerlet(MPI_Comm new_comm);

        void SetEffectiveRadius(MPI_Comm new_comm);

        void FreeCommunicator(MPI_Comm new_comm);

        void Coupling(std::list<BODY>& bodies_list, MPI_Comm new_comm);

        int  OnBoundary( int ku,int i);

        void ComputeForce(MPI_Comm new_comm);

        void ComputeForce_LC(std::vector<std::vector< std::list<BODY> >  > &grid, int *nc);

        void InsertList( ParticleList **root_list, ParticleList *i );

        int  CheckBoundary(std::vector<BODY> body, int i);

        void WriteColliodsToFile_MPI(MPI_Comm new_comm);

        void SendColloids(std::list<BODY> &bodies1, MPI_Comm new_comm);

        void SendGridPoints_LC(std::vector<std::vector<std::list<BODY> > > &mygrid, const int *nc, MPI_Comm new_comm);

        void ResizeArray(std::vector<BODY> &bodies,int position, MPI_Comm new_comm);

        void InsertColloid(BODY *bodies,int postion, MPI_Comm new_comm);

        void ArrangeParticles_Cells(std::vector<std::vector<std::list<BODY> > > &grid, int *nc,std::list<BODY> &colloid_list, double *l,MPI_Comm new_comm);

       //	void DistributeColloids(MPI_Comm new_comm);

       void CreateVectorGrid(std::shared_ptr<BODY> grid);	 
	  /*--------------
       ####### polymer part  #######################
	  ---------------*/
       void UpdateSolution(MPI_Comm new_comm);

       void InitializePolymer(MPI_Comm new_comm);

       void setLaplacianBase(MPI_Comm new_comm);

       double g(double phi);

       void SetSecondLaplacian2(MPI_Comm new_comm);

       void FiniteDifferenceScheme(MPI_Comm new_comm);

       void ExchangeData(MPI_Comm new_comm, double **array);

       void WriteToFile_MPI(MPI_Comm new_comm);

       void Solver(); 

       void setIndex(MPI_Comm new_comm);

       void setGlobalIndex(MPI_Comm new_comm);

       int  Get_Process_ID(int index,const std::vector<BODY> &bodies) const;/* takes index of a paricles and returns its partition*/

       void SetCP4(MPI_Comm new_comm);

       // Node *CreateNode(double xim, double xmax,double ymin, double ymax);
       //void AddBodyToTree(Node* root, BODY* bd);
//       friend void Polymer2D::ExchangeData(MPI_Comm new_comm, double **array);		 






 };
#endif
