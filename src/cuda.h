/*
 * cuda.h
 *
 *  Created on: 15 Jul 2020
 *      Author: paz
 */

#ifndef CUDA_H_
#define CUDA_H_

#define THREADS 1000

void findBestCuda(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk, int start, int end);

#endif /* CUDA_H_ */
