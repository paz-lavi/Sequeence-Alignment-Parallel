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

const char *CONSERVATIVE_GROUPS[] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK",
		"NHQK", "FYW", "HY", "MILF" };
const char *SEMI_CONSERVATIVE_GROUPS[] = { "SAG", "ATV", "CSA", "SGND", "STPA",
		"STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };
const char HYPHEN = '-';

//This function checks all options (offset, MS (k)) and finds the best result
struct ms_results* findBest(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk ,int start , int end) {
	float f;
	int i, j;
	int l1 = strlen(*seq1);
	int l2 = strlen(*seq2);
	int max_offset = l1 - l2;

	for (j = start; j < end; j++) {  //all MS(k) options

		mutantSequence(seq2, mk, j);

		for (i = 0; i < max_offset - 1; i++) {   //all offsets options

			f = compereSeq1AndSeq2(seq1, mk, fd, i); // compare seq1 and current MS(j) and offset = i
			if (f > res->score) {
				res->score = f;
				res->offset = i;
				res->k = j;
			}
		}
		resetArray(l2 + 1, mk);  //set all char at mk to '\0'
	}

	return res;
}

//compare 2 sequence. return score as float
float compereSeq1AndSeq2(char **seq1, char **seq2, struct file_data fd,
		int offset) {

	int l1 = strlen(*seq1);
	int l2 = strlen(*seq2);

	//printf("l2 #2 = %d\n",l2);
	//printf("hi from compereSeq1AndSeq2\n");

	int numberOfStars = 0, numberOfColons = 0, numberOfPoints = 0,
			numberOfSpaces = 0, i; //init var for calculate the score

	if (l1 == l2) {
		for (i = 0; i < l1; i++) {
			switch (compareTowChars((*seq2)[i], (*seq1)[i])) { //in case seq1 and seq2 have the same length
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

			}
		}
	} else if (l1 > l2) {                    //in case the length of seq1 > seq2
		int max_offset = l1 - l2; // calculate the maximum offset for this compare
		if (offset > max_offset) {
			fprintf(stderr, "offset =  %d > max offset = %d\n", offset,
					max_offset);  // invalid offset
			exit(1);
		} else {

			for (i = 0; i < l2; i++) {
				switch (compareTowChars((*seq2)[i], (*seq1)[offset + i])) { // get sing for each pair of chars with specific offset
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

				}
			}
		}
	} else {
		return LENGTH_ERROR; //in case the length of seq1 < seq2 - invalid length
	}
	return calculateSimilarity(numberOfStars, numberOfColons, numberOfPoints,
			numberOfSpaces, fd); //calculate score of this compare
}

//compare pair
char compareTowChars(char a, char b) {
	if (a == b) {
		//printf("*\n");
		return '*';
	} else
		return conservativeGroups(a, b);
}

// compare pair in conservative groups
char conservativeGroups(char a, char b) {
	int i;
	for (i = 0; i < 9; i++) {
		if ((strchr(CONSERVATIVE_GROUPS[i], (int) a) != NULL)
				&& (strchr(CONSERVATIVE_GROUPS[i], (int) b) != NULL))
			return ':';
	}
	return semiConservativeGroups(a, b);
}
// compare pair in conservative groups
char semiConservativeGroups(char a, char b) {
	int i;
	for (i = 0; i < 11; i++) {
		if ((strchr(SEMI_CONSERVATIVE_GROUPS[i], (int) a) != NULL)
				&& (strchr(SEMI_CONSERVATIVE_GROUPS[i], (int) b) != NULL))
			return ':';
	}
	return ' ';
}

// calculate the score
float calculateSimilarity(int numberOfStars, int numberOfColons,
		int numberOfPoints, int numberOfSpaces, struct file_data fd) {
	return (fd.w1 * numberOfStars) - (fd.w2 * numberOfColons)
			- (fd.w3 * numberOfPoints) - (fd.w4 * numberOfSpaces);
}

//create mutant sequence  from a string/ copy the arr to des and inert '-' at the index i
void mutantSequence(char **from, char **des, int i) {
	int len = strlen(*from);
	if (len < i) {
		fprintf(stderr, "index %d out of bound. string length = %d\n", i, len);
		exit(1);
	}

	memcpy(*des, *from, i);
	strcpy(&(*des)[i], &HYPHEN);
	strcpy(&(*des)[i + 1], &(*from)[i]);

}

