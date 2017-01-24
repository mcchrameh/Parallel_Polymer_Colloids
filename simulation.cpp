#include "CDS_BASE.h"
//#include "CDS_Polymer.h"
//#include "CDS_Colloids.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h> 
#include "mpi.h"
 


int main(int argc, char **argv)
{
 
 
 MPI_Init(&argc, &argv);
 BASE code;
 code.Solver();

 MPI_Finalize();


 return 0;

}
