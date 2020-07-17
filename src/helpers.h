/*
 * helpers.h
 *
 *  Created on: 13 Jul 2020
 *      Author: paz
 */

#ifndef HELPERS_H_
#define HELPERS_H_
#define MASTER 0
#define SLAVE 1
#define SEQ1_MAX_LENGTH 3000
#define SEQ2_MAX_LENGTH 2000
struct ms_results
{
	double score;
	int k;
	int offset;
};

struct file_data
{
	double w1, w2, w3, w4;
	char *seq1;
	char **arrOfSeq2;
	int sizeOfSeq2Arr;
};

char*
myCharArrCalloc(int size);
char**
myStringArrCalloc(int size);
void
resetCharsArray(int size, char **arr);
struct ms_results*
resultMalloc();
struct ms_results**
resultArrayMalloc(int size);

#endif /* HELPERS_H_ */
