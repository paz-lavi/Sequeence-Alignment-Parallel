/*
 * comparition.h
 *
 *  Created on: 12 Jul 2020
 *      Author: paz
 */
#pragma once
#ifndef COMPARITION_H_
#define COMPARITION_H_
#include "main_func.h"
#define LENGTH_ERROR -99996546546.999
//const char* CONSERVATIVE_GROUPS[] = {"NDEQ" , "NEQK" , "STA","MILV" , "QHRK","NHQK","FYW","HY","MILF"};
//const char* SEMI_CONSERVATIVE_GROUPS[] = {"SAG" , "ATV" , "CSA","SGND" , "STPA","STNK","NEQHRK","NDEQHK","SNDEQK","HFY","FVLIM"};
void findBestOpenMP(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk, int start, int end);
struct ms_results* findBest(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk, int start, int end);
double compereSeq1AndSeq2(char **seq1, char **seq2, struct file_data fd,
		int offset);
char compareTowChars(char a, char b);
char conservativeGroups(char a, char b);
char semiConservativeGroups(char a, char b);
double calculateSimilarity(int numberOfStars, int numberOfColons,
		int numberOfPoints, int numberOfSpaces, struct file_data fd);
void mutantSequence(char **from, char **des, int i);
#endif /* COMPARITION_H_ */
