/*
 * main_func.h
 *
 *  Created on: 12 Jul 2020
 *      Author: paz
 */

#ifndef MAIN_FUNC_H_
#define MAIN_FUNC_H_
#include <stdio.h>
#include "helpers.h"
#define OUTPUT "output.txt"


void master(int argc, char* argv[]);
void slave();
void initDataFromFile(char *filename, struct file_data *fd);
void readLine(FILE *file, char **str, int max_size);
void writeResultsToFile(struct ms_results **res, int size);
void startCalculate(struct ms_results **res,struct file_data fd, int proc);
void calculateParallel(struct ms_results **res,struct file_data fd);

#endif /* MAIN_FUNC_H_ */
