/*
 * main_func.c
 *
 *  Created on: 12 Jul 2020
 *      Author: paz
 */
#include <stdio.h>
#include <stdlib.h>
#include "helpers.h"
#include "comparition.h"
#include "main_func.h"
#include "cuda.h"
#include "mpi.h"
#include <string.h>

/* master program*/
void master(int argc, char *argv[])
{
	double start = MPI_Wtime();

	if (argc != 2)
	{
		printf(argc < 2 ? "Data file is required!" :    //handle input error
				"Only data file is required!");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	struct file_data fd;

	/*read data from file*/
	initDataFromFile(argv[1], &fd);

	int i;
	/* allocate char array for creating MS(i)*/
	char *mk = myCharArrCalloc(SEQ2_MAX_LENGTH + 1);

	/*allocate array of result. each index represent result of one sequence*/
	struct ms_results **res = resultArrayMalloc(fd.sizeOfSeq2Arr);

	/*calculate the best score for all seq2 in the file */
	calculateParallel(res, fd);

	/* write results to output.txt*/
	writeResultsToFile(res, fd.sizeOfSeq2Arr);

	/*free all allocated vars*/

	free(mk);
	free(fd.seq1);
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{
		free(fd.arrOfSeq2[i]);
		free(res[i]);
	}
	free(res);
	free(fd.arrOfSeq2);
	double end = MPI_Wtime();
	printf("running time = %lf", end - start);
}
/* Slave program */
void slave()
{

	struct file_data fd; /*get data from master into fd for using in slave program*/
	struct ms_results **res; /*array of results of slave calculations */
	int i;
	double weight[4]; /* array for getting weight from master*/
	MPI_Status status;

	/*receive numbers of sequences  to slave*/
	MPI_Recv(&fd.sizeOfSeq2Arr, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);

	/*space allocation for seq1*/
	fd.seq1 = myCharArrCalloc(SEQ1_MAX_LENGTH + 1);

	/*receive seq1 from master*/
	MPI_Recv(fd.seq1, SEQ1_MAX_LENGTH + 1, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD,
			&status);

	fd.arrOfSeq2 = myStringArrCalloc(fd.sizeOfSeq2Arr); /*space allocation for seq2 array*/
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{
		/*space allocation for each seq2 */
		fd.arrOfSeq2[i] = myCharArrCalloc(SEQ2_MAX_LENGTH + 1);

		/*receive each seq2 from master*/
		MPI_Recv(fd.arrOfSeq2[i], SEQ2_MAX_LENGTH + 1, MPI_CHAR, MASTER, 0,
		MPI_COMM_WORLD, &status);
	}
	/*receive weights from master*/
	MPI_Recv(&weight, 4, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);

	fd.w1 = weight[0];
	fd.w2 = weight[1];
	fd.w3 = weight[2];
	fd.w4 = weight[3];

	/*space allocation for results array*/
	res = resultArrayMalloc(fd.sizeOfSeq2Arr);

	/*perform all slave's calculations with cuda*/
	cudaCalculeteAll(fd, res, SLAVE);

	/* for each seq2 send the result to master*/
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{

		MPI_Send(&res[i]->score, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(&res[i]->offset, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(&res[i]->k, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

	}

	/*Release memory */
	free(fd.seq1);
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
		free(fd.arrOfSeq2[i]);
	free(fd.arrOfSeq2);
}

void calculateParallel(struct ms_results **res, struct file_data fd)
{
	MPI_Status status;
	int i;

	/*create array of weight that will be send to slave */
	double weight[] =
	{ fd.w1, fd.w2, fd.w3, fd.w4 };

	/*space allocation for slave's results array */
	struct ms_results **temp_res = resultArrayMalloc(fd.sizeOfSeq2Arr);
	;
	/*send num of sequences of  to slave*/
	MPI_Send(&fd.sizeOfSeq2Arr, 1, MPI_INT, SLAVE, 0, MPI_COMM_WORLD);

	/*send seq1 to slave*/
	MPI_Send(fd.seq1, SEQ1_MAX_LENGTH + 1, MPI_CHAR, SLAVE, 0, MPI_COMM_WORLD);

	/*send each seq2 to slave*/
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
		MPI_Send(fd.arrOfSeq2[i], SEQ2_MAX_LENGTH + 1, MPI_CHAR, SLAVE, 0,
		MPI_COMM_WORLD);

	/*send weights to slave*/
	MPI_Send(&weight, 4, MPI_DOUBLE, SLAVE, 0, MPI_COMM_WORLD);

	/*perform all master's calculations with openmp */
	startCalculateOpenMP(res, fd, MASTER);

	/*receive all results from slave*/
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{
		MPI_Recv(&temp_res[i]->score, 1, MPI_DOUBLE, SLAVE, 0, MPI_COMM_WORLD,
				&status);
		MPI_Recv(&temp_res[i]->offset, 1, MPI_INT, SLAVE, 0, MPI_COMM_WORLD,
				&status);
		MPI_Recv(&temp_res[i]->k, 1, MPI_INT, SLAVE, 0, MPI_COMM_WORLD,
				&status);

	}

	/*compare master (openmp) and slave (cuda) results of each seq2. set res[i] to the best score/ */
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{
		if (temp_res[i]->score > res[i]->score)
		{
			res[i]->score = temp_res[i]->score;
			res[i]->offset = temp_res[i]->offset;
			res[i]->k = temp_res[i]->k;
		}
	}

	/*Release memory */
	free(temp_res);
}

void startCalculateOpenMP(struct ms_results **res, struct file_data fd,
		int proc)
{
	/*space allocation Mutant Sequence */
	char *mk = myCharArrCalloc(SEQ2_MAX_LENGTH + 1);
	int i, start, end;

	//char **s1 = myStringArrCalloc(fd.sizeOfSeq2Arr);

	/*find best score for each sequence with openmp*/
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{
		/*set start and end for each seq2*/
		start = proc == MASTER ? 1 : strlen((fd.arrOfSeq2[i])) / 2;
		end = proc == MASTER ?
				strlen(fd.arrOfSeq2[i]) / 2 + 1 : strlen(fd.arrOfSeq2[i]);

		/*reset res[i]*/
		res[i]->k = -1;
		res[i]->offset = -1;
		res[i]->score = LENGTH_ERROR;  //reset values
		findBestOpenMP(&fd.seq1, &fd.arrOfSeq2[i], fd, res[i], &mk, start, end); /*find best score for seq2 in proc range*/

	}

	/*Release memory */
	free(mk);
}

/*reading the file. the the return to caller through fd*/
void initDataFromFile(char *filename, struct file_data *fd)
{
	FILE *file = fopen(filename, "r"); /* open input file*/
	if (!file)
	{
		fprintf(stderr, "Could not open file: %s\n", filename);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int i, numOfSeq;
	double w1, w2, w3, w4; /* define weight vars*/
	if (fscanf(file, "%lf%lf%lf%lf\n", &w1, &w2, &w3, &w4) != 4) /* read first line into weight vars */
	{
		fprintf(stderr, "row %d is not represented correctly!",
		__LINE__);
		MPI_Abort(MPI_COMM_WORLD, 1);

	}
	char *seq1 = myCharArrCalloc(SEQ1_MAX_LENGTH + 1); /*space allocation  for seq1*/
	;
	readLine(file, &seq1, SEQ1_MAX_LENGTH); /* read seq1 from file */

	if (fscanf(file, "%d\n", &numOfSeq) != 1) /*read number of NS2 of sequences seq2 to check against seq1*/
	{
		fprintf(stderr, "row %d is not represented correctly!",
		__LINE__);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	char **arrOfSeq2 = myStringArrCalloc(numOfSeq); /* space allocation  for arr of seq2 lines*/

	for (i = 0; i < numOfSeq; i++)
	{
		arrOfSeq2[i] = myCharArrCalloc(SEQ2_MAX_LENGTH); /*space allocation  for  seq2*/
		readLine(file, &(arrOfSeq2[i]), SEQ2_MAX_LENGTH); /*read seq2*/
	}

	fclose(file);

	/*insert the result into fd*/
	fd->w1 = w1;
	fd->w2 = w2;
	fd->w3 = w3;
	fd->w4 = w4;
	fd->seq1 = seq1;
	fd->arrOfSeq2 = arrOfSeq2;
	fd->sizeOfSeq2Arr = numOfSeq;

}

/*read line of strings (seq)*/
void readLine(FILE *file, char **str, int max_size)
{
	char c;
	int i = 0;

	if (fscanf(file, "%c", &c) != 1)/*read first char*/
	{
		fprintf(stderr, "char %d is not represented correctly!",
		__LINE__);

		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	/* while not end of line*/
	while (c != '\n' && i < max_size)
	{
		(*str)[i++] = c; /*append c to string*/

		if (fscanf(file, "%c", &c) != 1) /*read next char*/
		{
			fprintf(stderr, "char %d is not represented correctly!",
			__LINE__);

			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	(*str)[i] = '\0'; /*set end  of string*/

}

/* write result to output.txt file */
void writeResultsToFile(struct ms_results **res, int size)
{
	FILE *file = fopen(OUTPUT, "w"); /*open output.txt : override if file exist else create it  */
	if (!file)
	{
		fprintf(stderr, "Could generate file: %s\n", OUTPUT);
		MPI_Abort(MPI_COMM_WORLD, 1);

	}
	int i;

	fprintf(file, "Score For Each Sequence Alignment:\n\n");

	/*write result of each sequence*/
	for (i = 0; i < size; i++)
	{
		fprintf(file, "Sequence Alignment #%d:\n"
				"Best Score: %f\n"
				"offset = %d\n"
				"MS(k): k = %d\n\n", i, res[i]->score, res[i]->offset,
				res[i]->k);

	}
	fclose(file);
}
