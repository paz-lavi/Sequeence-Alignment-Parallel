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
#include <time.h>

#define MAX_PROC 2
int main(int argc, char *argv[]) {
	clock_t begin = clock();
	int my_rank; /* rank of process */
	int numOfProcs; /* number of processes */
	//int source;   /* rank of sender */
	//int dest;     /* rank of receiver */
//	int tag=0;    /* tag for messages */
	//char message[100];        /* storage for message */
	//MPI_Status status; /* return status for receive */

	/* start up MPI */

	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);

	if (numOfProcs != MAX_PROC) {
		printf("Program requires %d processes\n", MAX_PROC);
		MPI_Finalize();
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (my_rank == MASTER) {
		master(argc, argv);
	} else if (my_rank == SLAVE) {
		slave();
	}

	/* shut down MPI */
	MPI_Finalize();

	if (my_rank == MASTER) {
		clock_t end = clock();
		double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
		printf("runtime = %f", time_spent);
	}
	return 0;
}

