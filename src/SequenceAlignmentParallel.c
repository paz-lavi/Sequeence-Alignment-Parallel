/*
 ============================================================================
 Name        : SequeenceAlignmentParallel.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello MPI World in C 
 ============================================================================
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "helpers.h"
#include "main_func.h"
#include "comparition.h"

#define MAX_PROC 2
int main(int argc, char *argv[])
{
	int my_rank; /* rank of process */
	int numOfProcs; /* number of processes */

	/* start up MPI */

	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);

	if (numOfProcs != MAX_PROC)
	{
		printf("Program requires %d processes\n", MAX_PROC);
		MPI_Finalize();
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (my_rank == MASTER)
	{
		master(argc, argv); // start master work
	}
	else if (my_rank == SLAVE)
	{
		slave(); // start slave work
	}

	/* shut down MPI */
	MPI_Finalize();



	return 0;
}

