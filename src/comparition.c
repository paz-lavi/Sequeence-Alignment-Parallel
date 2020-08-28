/*
 * comparition.c
 *
 *  Created on: 12 Jul 2020
 *      Author: paz
 */

#include "comparition.h"
#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* define CONSERVATIVE_GROUPS */
const char *CONSERVATIVE_GROUPS[] =
{ "NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF" };

/* define SEMI CONSERVATIVE_GROUPS */
const char *SEMI_CONSERVATIVE_GROUPS[] =
{ "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK",
		"HFY", "FVLIM" };
/* define HYPHEN ('-') */

const char HYPHEN = '-';

/*find best score of seq2 with ms(k) and offset = n*/
void findBestOpenMP(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk, int start, int end)
{

	/*get length of seq1 and seq2*/
	int l1 = strlen(*seq1);
	int l2 = strlen(*seq2);

	/*calculate the maximum offset for this 2 sequences */
	int max_offset = l1 - l2;

	/*array of seq2 after mk(k)  for each openmp thread */
	char **t_mk;
	/*array of result for each openmp thread*/
	struct ms_results **temp_res;

#pragma omp parallel
	{
		double d;
		const int tid = omp_get_thread_num();
		const int nthreads = omp_get_num_threads();

#pragma omp single
		{
			/*space allocation - performing only once*/
			t_mk = myStringArrCalloc(nthreads);
			temp_res = resultArrayMalloc((nthreads));


		}

#pragma omp for
		/*space allocation. each thread has is own string of seq2 after mk*/
		for (int i = 0; i < nthreads; i++)
			t_mk[i] = myCharArrCalloc(SEQ2_MAX_LENGTH);

#pragma omp for
		/*find best score for current seq2 in the proc range*/
		for (int j = start; j < end; j++)
		{
			/*get mk(j) of seq2*/
			mutantSequence(seq2, &t_mk[tid], j);

			/*for all offsets options*/
			for (int i = 0; i < max_offset - 1; i++)
			{
				/*calculate score of Sequence Alignment seq1 and seq2 when seq2 after mk(j) and offset = i */
				d = compereSeq1AndSeq2(seq1, &t_mk[tid], fd, i); // compare seq1 and current MS(j) and offset = i
				if (d > temp_res[tid]->score) /*if new high score update  thread result*/
				{
					temp_res[tid]->score = d;
					temp_res[tid]->offset = i;
					temp_res[tid]->k = j;
				}
			}
			resetCharsArray(l2 + 1, &t_mk[tid]);  /*set all char at mk to '\0'*/
		}
#pragma omp for
		/*set res to the best result of all threads.*/
		for (int j = 0; j < nthreads; j++)
		{
			if (temp_res[j]->score > res->score)
			{
				res->score = temp_res[j]->score;
				res->offset = temp_res[j]->offset;
				res->k = temp_res[j]->k;
			}

		}

#pragma omp for
		/*Release memory */
		for (int j = 0; j < nthreads; j++)
			free(t_mk[j]);

#pragma omp single
		{
			/*Release memory */

			free(temp_res);

		}

	}
}

// /*This function checks all options (offset, MS (k)) and finds the best result*/
// struct ms_results* findBest(char **seq1, char **seq2, struct file_data fd,
// 		struct ms_results *res, char **mk, int start, int end)
// {
// 	float f;
// 	int i, j;
// 	int l1 = strlen(*seq1);
// 	int l2 = strlen(*seq2);
// 	int max_offset = l1 - l2;

// 	for (j = start; j < end; j++)
// 	{  //all MS(k) options

// 		mutantSequence(seq2, mk, j);

// 		for (i = 0; i < max_offset - 1; i++)
// 		{   //all offsets options

// 			f = compereSeq1AndSeq2(seq1, mk, fd, i); // compare seq1 and current MS(j) and offset = i
// 			if (f > res->score)
// 			{
// 				res->score = f;
// 				res->offset = i;
// 				res->k = j;
// 			}
// 		}
// 		resetArray(l2 + 1, mk);  //set all char at mk to '\0'
// 	}

// 	return res;
// }

/*compare 2 sequence. return score as double*/
double compereSeq1AndSeq2(char **seq1, char **seq2, struct file_data fd,
		int offset)
{
	/*get length of seq1 and seq2*/
	int l1 = strlen(*seq1);
	int l2 = strlen(*seq2);



	int numberOfStars = 0, numberOfColons = 0, numberOfPoints = 0,
			numberOfSpaces = 0, i; /*init var for calculate the score*/

	if (l1 == l2)  /*in case seq1 and seq2 have the same length*/
	{
		for (i = 0; i < l1; i++)
		{
			switch ( compareTowChars((*seq2)[i], (*seq1)[i]))  /* get sing for each pair of chars */
			{
			case ' ':
				numberOfSpaces++;
				break;
			case '.':
				numberOfPoints++;
				break;
			case ':':
				numberOfColons++;
				break;
			case '*':
				numberOfStars++;
				break;
			default:
				break;
			}
		}
	}
	else if (l1 > l2) /*in case the length of seq1 > seq2*/
	{
		int max_offset = l1 - l2; /* calculate the maximum offset for this compare*/
		if (offset > max_offset)
		{
			fprintf(stderr, "offset =  %d > max offset = %d\n", offset,
					max_offset);  /* invalid offset*/
			exit(1);
		}
		else
		{

			for (i = 0; i < l2; i++)
			{
				switch ((*seq2)[i] == HYPHEN ?
						HYPHEN :
						compareTowChars((*seq2)[i], (*seq1)[offset + i])) /* get sing for each pair of chars with specific offset*/
				{
				case ' ':
					numberOfSpaces++;
					break;
				case '.':
					numberOfPoints++;
					break;
				case ':':
					numberOfColons++;
					break;
				case '*':
					numberOfStars++;
					break;
				default:
					break;
				}
			}
		}
	}
	else
	{
		return LENGTH_ERROR; /*in case the length of seq1 < seq2 : invalid length*/
	}
	return calculateSimilarity(numberOfStars, numberOfColons, numberOfPoints,
			numberOfSpaces, fd); /*calculate score of this compare*/
}

/*compare pair*/
char compareTowChars(char a, char b)
{
	if (a == b)
	{
		return '*';
	}
	else
		return conservativeGroups(a, b);
}

/* compare pair in conservative groups*/
char conservativeGroups(char a, char b)
{
	int i;
	for (i = 0; i < 9; i++)
	{
		if ((strchr(CONSERVATIVE_GROUPS[i], (int) a) != NULL)
				&& (strchr(CONSERVATIVE_GROUPS[i], (int) b) != NULL))
			return ':';
	}
	return semiConservativeGroups(a, b);
}
/* compare pair in conservative groups*/
char semiConservativeGroups(char a, char b)
{
	int i;
	for (i = 0; i < 11; i++)
	{
		if ((strchr(SEMI_CONSERVATIVE_GROUPS[i], (int) a) != NULL)
				&& (strchr(SEMI_CONSERVATIVE_GROUPS[i], (int) b) != NULL))
			return ':';
	}
	return ' ';
}

/* calculate the score*/
double calculateSimilarity(int numberOfStars, int numberOfColons,
		int numberOfPoints, int numberOfSpaces, struct file_data fd)
{
	return (fd.w1 * numberOfStars) - (fd.w2 * numberOfColons)
			- (fd.w3 * numberOfPoints) - (fd.w4 * numberOfSpaces);
}

/*create mutant sequence  from a string/ copy the arr to des and inert '-' at the index i*/
void mutantSequence(char **from, char **des, int i)
{
	/*get len of original string*/
	int len = strlen(*from);
	if (len < i)
	{
		fprintf(stderr, "index %d out of bound. string length = %d\n", i, len);
		exit(1);
	}
/*copy all char from the first until the HYPHEN location */
	memcpy(*des, *from, i);
	/*insert HYPHEN*/
	strcpy(&(*des)[i], &HYPHEN);
	/*copy the rest of chars*/
	strcpy(&(*des)[i + 1], &(*from)[i]);

}

