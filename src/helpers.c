/*
 * helpers.c
 *
 *  Created on: 13 Jul 2020
 *      Author: paz
 */
#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/*calloc char array*/
char* myCharArrCalloc(int size)
{
	/* space allocation*/
	char *arr = (char*) calloc(size, sizeof(char));
	/*check if allocation succeeded. if not the program abort*/
	if (arr == NULL)
	{
		fprintf(stderr, "Could not allocate array");

	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return arr;
}

/*set all char to '\0'*/
void resetCharsArray(int size, char **arr)
{
	int i;
	for (i = 0; i < size; i++)
		(*arr)[i] = '\0';

}

/*calloc arrray of Strings*/
char** myStringArrCalloc(int size)
{
	/*space allocation*/
	char **arr = (char**) calloc(size, sizeof(char*));
	/*check if allocation succeeded. if not the program abort*/
	if (arr == NULL)
	{
		fprintf(stderr, "Could not allocate array");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return arr;
}
/*malloc  ms_results*/
struct ms_results* resultMalloc()
{
	/*space allocation*/
	struct ms_results *res = (struct ms_results*) malloc(sizeof(struct ms_results) * 1);
	/*check if allocation succeeded. if not the program abort*/
	if (res == NULL)
	{
		fprintf(stderr, "Could not allocate ms_results");
			MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return res;
}

/*malloc  ms_results array */

struct ms_results** resultArrayMalloc(int size)
{
	/*space allocation*/
	struct ms_results **res = (struct ms_results**) malloc(
			sizeof(struct ms_results*) * size);
	/*check if allocation succeeded. if not the program abort*/
	if (res == NULL)
	{
		fprintf(stderr, "Could not allocate ms_results");
			MPI_Abort(MPI_COMM_WORLD, 1);
	}
	int i;
	for (i = 0; i < size; i++)
	{
		/*space allocation for each ms_results*/
		res[i] = resultMalloc();
	}
	return res;
}

