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
#include "mpi.h"
#include <string.h>
void master(int argc, char *argv[]) {

	if (argc != 2) {
		printf(argc < 2 ? "Data file is required!" :    //handle input error
				"Only data file is required!");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	struct file_data fd;

//read data from file
	initDataFromFile(argv[1], &fd);

	int i;
// allocate char array for creating MS(i)
	char *mk = myCharArrCalloc(SEQ2_MAX_LENGTH + 1);

	//allocate array of result. each index represent result of one sequence
	struct ms_results **res = resultArrayMalloc(fd.sizeOfSeq2Arr);

	//find best score for each sequence
//	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
//		res[i]->k = -1;
//		res[i]->offset = -1;
//		res[i]->score = LENGTH_ERROR;  //reset values
//		findBest(&fd.seq1, &fd.arrOfSeq2[i], fd, res[i], &mk); // calculate
//
//		printf("seq2 # %d :  score = %f , offset = %d , k = %d \n", i,
//				res[i]->score, res[i]->offset, res[i]->k); // print result to console
//
//	}
//
	calculateParallel(res, fd);

	// write results to output.txt
	writeResultsToFile(res, fd.sizeOfSeq2Arr);

	//free all allocated vars
	free(mk);
	free(fd.seq1);
	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
		free(fd.arrOfSeq2[i]);
		free(res[i]);
	}
	free(res);
	free(fd.arrOfSeq2);

}

void slave(MPI_Status status) {
	//receive num of sequences of  to slave
	struct file_data fd;
	struct ms_results **res;
	int i;
	float weigth[4];

	MPI_Recv(&fd.sizeOfSeq2Arr, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);

	fd.seq1 = myCharArrCalloc(SEQ1_MAX_LENGTH + 1);
	MPI_Recv(fd.seq1, SEQ1_MAX_LENGTH + 1, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD,
			&status);

	fd.arrOfSeq2 = myStringArrCalloc(fd.sizeOfSeq2Arr);
	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
		fd.arrOfSeq2[i] = myCharArrCalloc(SEQ2_MAX_LENGTH + 1);
		MPI_Recv(fd.arrOfSeq2[i], SEQ2_MAX_LENGTH + 1, MPI_CHAR, MASTER, 0,
				MPI_COMM_WORLD, &status);
	}
	MPI_Recv(&weigth, 4, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, &status);

	fd.w1 = weigth[0];
	fd.w2 = weigth[1];
	fd.w3 = weigth[2];
	fd.w4 = weigth[3];

	res = resultArrayMalloc(fd.sizeOfSeq2Arr);
	startCalculate(res, fd, SLAVE);

	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
		MPI_Send(&res[i]->score, 1, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(&res[i]->offset, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		MPI_Send(&res[i]->k, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

	}
	printf("all send back to master\n");
	free(fd.seq1);
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
		free(fd.arrOfSeq2[i]);
	free(fd.arrOfSeq2);
}

void calculateParallel(struct ms_results **res, struct file_data fd) {
	int i;
	float weigth[] = { fd.w1, fd.w2, fd.w3, fd.w4 };
	MPI_Status status;
	struct ms_results **temp_res = resultArrayMalloc(fd.sizeOfSeq2Arr);
	;
	//send num of sequences of  to slave
	MPI_Send(&fd.sizeOfSeq2Arr, 1, MPI_INT, SLAVE, 0, MPI_COMM_WORLD);
	MPI_Send(fd.seq1, SEQ1_MAX_LENGTH + 1, MPI_CHAR, SLAVE, 0, MPI_COMM_WORLD);

	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
		MPI_Send(fd.arrOfSeq2[i], SEQ2_MAX_LENGTH + 1, MPI_CHAR, SLAVE, 0,
				MPI_COMM_WORLD);
	MPI_Send(&weigth, 4, MPI_FLOAT, SLAVE, 0, MPI_COMM_WORLD);

	startCalculate(res, fd, MASTER);

	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
		MPI_Recv(&temp_res[i]->score, 1, MPI_FLOAT, SLAVE, 0, MPI_COMM_WORLD,
				&status);
		MPI_Recv(&temp_res[i]->offset, 1, MPI_INT, SLAVE, 0, MPI_COMM_WORLD,
				&status);
		MPI_Recv(&temp_res[i]->k, 1, MPI_INT, SLAVE, 0, MPI_COMM_WORLD,
				&status);

	}
	printf("master recived all\n");
	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
		if (temp_res[i]->score > res[i]->score) {
			res[i]->score = temp_res[i]->score;
			res[i]->offset = temp_res[i]->offset;
			res[i]->k = temp_res[i]->k;
		}
	}
	free(temp_res);
}

void startCalculate(struct ms_results **res, struct file_data fd, int proc) {
	char *mk = myCharArrCalloc(SEQ2_MAX_LENGTH + 1);
	int i, start, end;
	//find best score for each sequence
	for (i = 0; i < fd.sizeOfSeq2Arr; i++) {
		start = proc == MASTER ? 1 : strlen((fd.arrOfSeq2[i])) / 2;
		end = proc == MASTER ?
				strlen(fd.arrOfSeq2[i]) / 2 + 1 : strlen(fd.arrOfSeq2[i]);
		res[i]->k = -1;
		res[i]->offset = -1;
		res[i]->score = LENGTH_ERROR;  //reset values
		//findBest(&fd.seq1, &fd.arrOfSeq2[i], fd, res[i], &mk, start, end); // calculate
		findBestOpenMP(&fd.seq1, &fd.arrOfSeq2[i], fd, res[i], &mk, start, end); // calculate
		printf("seq2 # %d :  score = %f , offset = %d , k = %d \n", i,
				res[i]->score, res[i]->offset, res[i]->k); // print result to console

	}
	free(mk);
}

//reading the file. the the return to caller through fd
void initDataFromFile(char *filename, struct file_data *fd) {
	FILE *file = fopen(filename, "r");  // open input file
	if (!file) {
		fprintf(stderr, "Could not open file: %s\n", filename);
		exit(1);
//	MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int i, numOfSeq;
	float w1, w2, w3, w4; // define weight vars
	if (fscanf(file, "%f%f%f%f\n", &w1, &w2, &w3, &w4) != 4) { // read first line into weight vars
		fprintf(stderr, "row %d is not represented correctly!",
		__LINE__);
		exit(1);
		//MPI_Abort(MPI_COMM_WORLD, 1);

	}
	char *seq1 = myCharArrCalloc(SEQ1_MAX_LENGTH + 1); //space allocation  for seq1;
	readLine(file, &seq1, SEQ1_MAX_LENGTH); // read seq1 from file

	if (fscanf(file, "%d\n", &numOfSeq) != 1) { //read number of NS2 of sequences seq2 to check against seq1
		fprintf(stderr, "rowwww %d is not represented correctly!",
		__LINE__);
		exit(1);
		//MPI_Abort(MPI_COMM_WORLD, 1);
	}
	char **arrOfSeq2 = myStringArrCalloc(numOfSeq); // space allocation  for arr of seq2 lines

	for (i = 0; i < numOfSeq; i++) {
		arrOfSeq2[i] = myCharArrCalloc(SEQ2_MAX_LENGTH); // space allocation  for  seq2
		readLine(file, &(arrOfSeq2[i]), SEQ2_MAX_LENGTH); //read seq2
	}

	fclose(file);
	//insert the result into fd
	fd->w1 = w1;
	fd->w2 = w2;
	fd->w3 = w3;
	fd->w4 = w4;
	fd->seq1 = seq1;
	fd->arrOfSeq2 = arrOfSeq2;
	fd->sizeOfSeq2Arr = numOfSeq;

}

//read line of strings (seq)
void readLine(FILE *file, char **str, int max_size) {
	char c;
	int i = 0;

	if (fscanf(file, "%c", &c) != 1) {   //read first char
		fprintf(stderr, "char %d is not represented correctly!",
		__LINE__);
		exit(1);
		//	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	while (c != '\n' && i < max_size) { // while not end of line
		(*str)[i++] = c;
		//strcpy((*str + i++) , c); // insert line into char array

		if (fscanf(file, "%c", &c) != 1) { // read next char
			fprintf(stderr, "char %d is not represented correctly!",
			__LINE__);
			exit(1);
			//	MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	(*str)[i] = '\0'; //set end  of string;

}

// write result to output.txt file
void writeResultsToFile(struct ms_results **res, int size) {
	FILE *file = fopen(OUTPUT, "w");
	if (!file) {
		fprintf(stderr, "Could generate file: %s\n", OUTPUT);
		//MPI_Abort(MPI_COMM_WORLD, 1);
		exit(1);
	}
	int i;

	fprintf(file, "Score For Each Sequence Alignment:\n\n");

	//write result of each sequence
	for (i = 0; i < size; i++) {
		fprintf(file, "Sequence Alignment #%d:\n"
				"Best Score: %f\n"
				"offset = %d\n"
				"MS(k): k = %d\n\n", i, res[i]->score, res[i]->offset,
				res[i]->k);

	}
	fclose(file);
}
