/*
 * helpers.c
 *
 *  Created on: 13 Jul 2020
 *      Author: paz
 */
#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>

//calloc char array
char* myCharArrCalloc(int size) {
	char *arr = (char*) calloc(size, sizeof(char)); // space allocation
	if (arr == NULL) { //check if allocation succeeded. if not the program abort
		fprintf(stderr, "Could not allocate array");
		exit(1);
		//	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return arr;
}

//set all char to '\0'
void resetArray(int size, char **arr) {
	//printf("reset\n");
	int i;
	for (i = 0; i < size; i++)
		(*arr)[i] = '\0';

}

//calloc arrray of Strings
char** myStringArrCalloc(int size) {
	char **arr = (char**) calloc(size, sizeof(char*)); // space allocation
	if (arr == NULL) { //check if allocation succeeded. if not the program abort
		fprintf(stderr, "Could not allocate array");
		exit(1);
		//	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return arr;
}
struct ms_results* resultMalloc() {
	struct ms_results *res = (struct ms_results*) malloc(
			sizeof(struct ms_results) * 1);
	if (res == NULL) { //check if allocation succeeded. if not the program abort
		fprintf(stderr, "Could not allocate ms_results");
		exit(1);
		//	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return res;
}

struct ms_results** resultArrayMalloc(int size) {
	struct ms_results **res = (struct ms_results**) malloc(
			sizeof(struct ms_results*) * size);
	if (res == NULL) { //check if allocation succeeded. if not the program abort
		fprintf(stderr, "Could not allocate ms_results");
		exit(1);
		//	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	int i;
	for (i = 0; i < size; i++) {
		res[i] = resultMalloc();
	}
	return res;
}


