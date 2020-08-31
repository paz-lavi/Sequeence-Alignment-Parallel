/*
 * cuda.c
 *
 *  Created on: 15 Jul 2020
 *      Author: paz
 */
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "helpers.h"
#include <stdio.h>

/*check if string contains specific char  on GPU. return 1 for true 0 for false . original strchr from string.h not works on GPU*/
__device__ int _strchr(const char *str, char c)
{
	int i = 0;
	do
	{
		if (str[i] == c)
			return 1;
	} while (str[++i] != '\0');
	return 0;
}

/*calculate length of string on GPU. original strlen from string.h not works on GPU*/
__device__ int _strlen(const char *str)
{
	int i = 0;
	while (str[i++] != '\0')
		;
	return i;
}

/* calculate the score*/
__device__ double calculateSimilarityCuda(int numberOfStars, int numberOfColons,
		int numberOfPoints, int numberOfSpaces, double w1, double w2, double w3,
		double w4)
{
	return (w1 * numberOfStars) - (w2 * numberOfColons) - (w3 * numberOfPoints)
			- (w4 * numberOfSpaces);
}

/* compare pair in conservative groups*/
__device__ char semiConservativeGroupsCuda(char a, char b)
{
	/* define SEMI CONSERVATIVE_GROUPS */
	const char *SEMI_CONSERVATIVE_GROUPS[] =
	{ "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK",
			"HFY", "FVLIM" };
	int i;
	for (i = 0; i < 11; i++)
	{
		if ((_strchr(SEMI_CONSERVATIVE_GROUPS[i], a))
				&& (_strchr(SEMI_CONSERVATIVE_GROUPS[i], b)))
			return ':';
	}
	return ' ';
}

/* compare pair in conservative groups*/
__device__ char conservativeGroupsCuda(char a, char b)
{
	/*define CONSERVATIVE_GROUPS */
	const char *CONSERVATIVE_GROUPS[] =
	{	"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	int i;
	for (i = 0; i < 9; i++)
	{
		if ((_strchr(CONSERVATIVE_GROUPS[i], a))
				&& (_strchr(CONSERVATIVE_GROUPS[i], b)))
			return ':';
	}
	return semiConservativeGroupsCuda(a, b);
}

/*compare pair*/
__device__ char compareTowCharsCuda(char a, char b)
{
	if (a == b)
	{

		return '*';
	}
	else
		return conservativeGroupsCuda(a, b);
}

/*compare 2 sequence. return score as double*/
__device__ double compereSeq1AndSeq2Cuda(char *seq1, char *seq2, int offset,
		double w1, double w2, double w3, double w4)
{

	int l1 = _strlen(seq1);
	int l2 = _strlen(seq2);

	/*init var for calculate the score*/
	int numberOfStars = 0, numberOfColons = 0, numberOfPoints = 0,
			numberOfSpaces = 0, i;

	if (l1 == l2)
	{
		for (i = 0; i < l1; i++) /*in case seq1 and seq2 have the same length*/
		{
			switch (compareTowCharsCuda((seq2)[i], (seq1)[i]))
			/* get sing for each pair of chars with specific offset*/
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
			return -99999999.999;
		}
		else
		{

			for (i = 0; i < l2; i++)
			{
				switch ( compareTowCharsCuda((seq2)[i], (seq1)[offset + i]))
				/* get sing for each pair of chars with specific offset*/
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
		return -99999999.999; /*in case the length of seq1 < seq2 - invalid length*/
	}
	return calculateSimilarityCuda(numberOfStars, numberOfColons,
			numberOfPoints, numberOfSpaces, w1, w2, w3, w4); /*calculate score of this compare*/
}

/*create mutant sequence  from a string/ copy the arr to des and inert '-' at the index i*/
__device__ void mutantSequenceCuda(char *from, char *des, int i)
{
	int j = 0, k = 0;
	do
	{
		if (k == i)
		{
			des[j++] = '-';
		}

		des[j++] = from[k]; /*copy char*/
	} while (from[k++] != '\0');

}

__global__ void addCalculateKernel(char *seq1, char *seq2,
		struct ms_results **temp_res, double w1, double w2, double w3,
		double w4, int max_offset, double *scores, int *offsetsArr, int *kArr,
		int start)
{
	int tid = threadIdx.x;
	double d;
	scores[tid] = -999999999;
	offsetsArr[tid] = -1;
	kArr[tid] = -1;
	char t[SEQ2_MAX_LENGTH + 1];

	/*get mk(tid + start) of seq2*/
	mutantSequenceCuda(seq2, t, tid + start);

	/*for all offsets options*/
	for (int i = 0; i < max_offset - 1; i++)
	{
		d = compereSeq1AndSeq2Cuda(seq1, t, i, w1, w2, w3, w4); /* compare seq1 and current MS(tid + start) and offset = i*/
		if (d > scores[tid])
		{
			scores[tid] = d;
			offsetsArr[tid] = i;
			kArr[tid] = start + tid;
		}
	}

}
/*set final_res to the best result of all mk(k).*/

__global__ void findBest(struct ms_results *final_res,
		struct ms_results **temp_res, int size, double *scores, int *offsetsArr,
		int *kArr)
{
	final_res->score = -1;
	final_res->offset = -1;
	final_res->k = -1;
	int i;
	for (i = 0; i < size; i++)
	{
		if (final_res->score < scores[i])
		{

			(final_res->score) = scores[i];
			(final_res->offset) = offsetsArr[i];
			(final_res->k) = kArr[i];
		}
	}

}

/*copy result from ms_results to pointers that will copy from GPU to host*/
__global__ void getRes(struct ms_results *final_res, double *score, int *offset,
		int *k)
{
	*score = final_res->score;
	*offset = final_res->offset;
	*k = final_res->k;
}

void cudaCalculeteAll(struct file_data fd, struct ms_results **res, int proc)
{
	char *s1, *s2;
	int i, l2, max_offset, start, end;
	struct ms_results **temp_res;
	struct ms_results *final_res;
	double *score;
	int *offset, *k;
	double *scores;
	int *offsetsArr, *kArr;

	int l1 = strlen(fd.seq1);

	int size = 1000;

	cudaError_t cudaStatus;
	/* Choose which GPU to run on, change this on a multi-GPU system.*/
	cudaStatus = cudaSetDevice(0);

	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &scores, size * sizeof(double));

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc scores failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &offsetsArr, size * sizeof(int));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc offsetsArr failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &kArr, size * sizeof(int));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc kArr failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void***) &temp_res,
			size * sizeof(struct ms_results*));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc temp_res failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &final_res, sizeof(struct ms_results*));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc final_res failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &score, sizeof(double));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc score failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &offset, sizeof(int));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc offset failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}

	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &k, sizeof(int));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc k failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &s1, SEQ1_MAX_LENGTH * sizeof(char));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc s1 failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
	/* space allocation*/
	cudaStatus = cudaMalloc((void**) &s2, SEQ2_MAX_LENGTH * sizeof(char));
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "malloc s2 failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}

	/* Copy seq1 vectors from host memory to GPU buffers.*/
	cudaStatus = cudaMemcpy(s1, fd.seq1, SEQ1_MAX_LENGTH * sizeof(char),
			cudaMemcpyHostToDevice);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "copy seq1 failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}

	/*for each seq2*/
	for (i = 0; i < fd.sizeOfSeq2Arr; i++)
	{
		/*get seq2 len*/
		l2 = strlen(fd.arrOfSeq2[i]);
		/*calculate maximum offset for this calculation*/
		max_offset = l1 - l2;

		/*calculate start and end range for this proc*/
		start = proc == MASTER ? 1 : strlen((fd.arrOfSeq2[i])) / 2;
		end = proc == MASTER ?
				strlen(fd.arrOfSeq2[i]) / 2 + 1 : strlen(fd.arrOfSeq2[i]);
		size = end - start + 1;

		/* Copy seq2 vectors from host memory to GPU buffers.*/
		cudaStatus = cudaMemcpy(s2, fd.arrOfSeq2[i],
		SEQ2_MAX_LENGTH * sizeof(char), cudaMemcpyHostToDevice);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "copy seq2 failed: %s\n",
					cudaGetErrorString(cudaStatus));
			goto Error;
		}

		/*calculate results on kernel. thread for each mk*/
		addCalculateKernel<<<1, size>>>(s1,s2,temp_res,fd.w1,fd.w2,fd.w3,fd.w4,max_offset ,scores,offsetsArr , kArr,start);

		/* Check for any errors launching the kernel*/
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "addKernel launch failed: %s\n",
					cudaGetErrorString(cudaStatus));
			goto Error;
		}

		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr,
					"cudaDeviceSynchronize returned error code %d after launching addCalculateKernel!\n%s\n",
					cudaStatus, cudaGetErrorString(cudaStatus));
			goto Error;
		}

		/*find the best score of seq2*/
		findBest<<<1, 1>>>(final_res,temp_res , size ,scores,offsetsArr , kArr);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "addKernel launch failed: %s\n",
					cudaGetErrorString(cudaStatus));
			goto Error;
		}

		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr,
					"cudaDeviceSynchronize returned error code %d after launching findBest!\n%s\n",
					cudaStatus, cudaGetErrorString(cudaStatus));
			goto Error;
		}

		/*copy results to pointers*/
		getRes<<<1,1>>>(final_res,score,offset,k);

		/* Copy best result from GPU buffer to host memory.*/
		cudaStatus = cudaMemcpy(&res[i]->score, score, sizeof(double),
				cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy  1 failed!\n %s\n",
					cudaGetErrorString(cudaStatus));
			goto Error;
		}
		cudaStatus = cudaMemcpy(&res[i]->k, k, sizeof(int),
				cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy 2  failed!\n %s\n",
					cudaGetErrorString(cudaStatus));
			goto Error;
		}
		cudaStatus = cudaMemcpy(&res[i]->offset, offset, sizeof(int),
				cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy 3 failed!\n %s\n",
					cudaGetErrorString(cudaStatus));
			goto Error;
		}

	}

	/*Release GPU memory */
	cudaFree(s1);
	cudaFree(s2);
	cudaFree(temp_res);
	cudaFree(final_res);
	cudaFree(kArr);
	cudaFree(offsetsArr);
	cudaFree(scores);
	cudaFree(k);
	cudaFree(offset);
	cudaFree(score);

	cudaDeviceReset();

	Error:
	/*Release GPU memory */
	cudaFree(s1);
	cudaFree(s2);
	cudaFree(temp_res);
	cudaFree(final_res);
	cudaFree(kArr);
	cudaFree(offsetsArr);
	cudaFree(scores);
	cudaFree(k);
	cudaFree(offset);
	cudaFree(score);
	cudaDeviceReset();

}

