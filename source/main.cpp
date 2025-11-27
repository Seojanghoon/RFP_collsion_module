#include <iostream>
#include <cstdlib>
#include <mpi.h>

#include "gyrokinetic.h"

using namespace std;

int main(int argc, char **argv) 
{
    /* init MPI parallelization */
    MPI_Init(&argc, &argv);
    int rank, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	//Gyrokinetic class initialization
	Gyrokinetic gyrokinetic(sokuri, MPI_COMM_WORLD);

	//Gyrokinetic main loop
	gyrokinetic.procs_mult();

	MPI_Finalize();

	return EXIT_SUCCESS;
}

