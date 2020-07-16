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


__device__ int _strchr(const char* str , char c){
int i = 0;
do{
if(str[i] == c)
return 1;
}while(str[++i]!= '\0');
return 0;
}
__device__ int _strlen(const char* str ){
int i = 0;
while(str[i++]!= '\0');
return i;
}

// calculate the score
__device__ double calculateSimilarityCuda(int numberOfStars, int numberOfColons,
		int numberOfPoints, int numberOfSpaces, double w1,double w2,double w3,double w4	) {
	return (w1 * numberOfStars) - (w2 * numberOfColons)
			- (w3 * numberOfPoints) - (w4 * numberOfSpaces);
}


// compare pair in conservative groups
__device__ char semiConservativeGroupsCuda(char a, char b) {
const char *SEMI_CONSERVATIVE_GROUPS[] = { "SAG", "ATV", "CSA", "SGND", "STPA",
		"STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };
	int i;
	for (i = 0; i < 11; i++) {
		if ((_strchr(SEMI_CONSERVATIVE_GROUPS[i],  a))
				&& (_strchr(SEMI_CONSERVATIVE_GROUPS[i],  b)))
			return ':';
	}
	return ' ';
}


// compare pair in conservative groups
__device__ char conservativeGroupsCuda(char a, char b) {
const char *CONSERVATIVE_GROUPS[] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK",
		"NHQK", "FYW", "HY", "MILF" };
	int i;
	for (i = 0; i < 9; i++) {
		if ((_strchr(CONSERVATIVE_GROUPS[i],  a) )
				&& (_strchr(CONSERVATIVE_GROUPS[i],  b) ))
			return ':';
	}
	return semiConservativeGroupsCuda(a, b);
}


//compare pair
__device__ char compareTowCharsCuda(char a, char b) {
	if (a == b) {
		//printf("*\n");
		return '*';
	} else
		return conservativeGroupsCuda(a, b);
}


__device__ //compare 2 sequence. return score as double
double compereSeq1AndSeq2Cuda(char *seq1, char *seq2,
		int offset ,double w1,double w2,double w3,double w4	) {

	int l1 = _strlen(seq1);
	int l2 = _strlen(seq2);

	//printf("l2 #2 = %d\n",l2);
//	printf("%s\n",seq2);

	int numberOfStars = 0, numberOfColons = 0, numberOfPoints = 0,
			numberOfSpaces = 0, i; //init var for calculate the score

	if (l1 == l2) {
		for (i = 0; i < l1; i++) {
			switch ((seq2)[i] == '-' ?
					'-' : compareTowCharsCuda((seq2)[i], (seq1)[i])) { //in case seq1 and seq2 have the same length
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
	} else if (l1 > l2) {                //in case the length of seq1 > seq2
		int max_offset = l1 - l2; // calculate the maximum offset for this compare
		if (offset > max_offset) {
			 // invalid offset
		
		} else {

			for (i = 0; i < l2; i++) {
				switch ((seq2)[i] == '-' ?
						'-' :
						compareTowCharsCuda((seq2)[i], (seq1)[offset + i])) { // get sing for each pair of chars with specific offset
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
	} else {
		return -9999.999; //in case the length of seq1 < seq2 - invalid length
	}
	return calculateSimilarityCuda(numberOfStars, numberOfColons, numberOfPoints,
			numberOfSpaces, w1,w2,w3,w4); //calculate score of this compare
}


__device__ void mutantSequenceCuda(char *from, char *des, int i) {
int j = 0 , k=0 ;
//printf("from - %s\n",from);
do{
if(k == i){
//printf("**\n");
des[j++] = '-';
//printf("****\n");
}

des[j++] =from[k];
//printf("j = %d - from = %c , des = %c\n" ,j, from[k],des[j++]);
}while(from[k++]!= '\0');

}


__global__ void addCalculateKernel(char *seq1, char *seq2, char** t_mk,  struct ms_results **temp_res,
double w1,double w2,double w3,double w4 ,int max_offset	,double* scores, int* offsetsArr , int* kArr ,int start) {
	int tid = threadIdx.x;
	double d;
	scores[tid] = -999999999;
					offsetsArr[tid] = -1;
					kArr[tid] = -1;
char t[SEQ2_MAX_LENGTH +1];
t_mk[tid] = t;
mutantSequenceCuda(seq2, t, tid);
	for (int i = 0; i < max_offset - 1; i++) {   //all offsets options
				d = compereSeq1AndSeq2Cuda(seq1, t, i , w1,w2,w3,w4); // compare seq1 and current MS(j) and offset = i
				if (d > scores[tid]) {
					scores[tid] = d;
					offsetsArr[tid] = i;
					kArr[tid] = start + tid;
				}
			}
			
		
}




__global__ void findBest( struct ms_results* final_res,	 struct ms_results** temp_res ,int size ,double* scores, int* offsetsArr , int* kArr){
					final_res->score = -1;
					final_res->offset = -1;
					final_res->k = -1;
					int i;
					printf("size = %d",size);
					for( i = 0 ; i < size;i++){
					//printf("score = %lf , offset = %d , k = %d \n",scores[i],offsetsArr[i],kArr[i]);
					if(
					final_res->score < scores[i]){
					
					(final_res->score) = scores[i] ;
					(final_res->offset) = offsetsArr[i];
					(final_res->k) = kArr[i];
					printf("toppppp i = %d\n",i);
					}
					}
					printf("cuda best score = %lf , offset = %d ,k = %d\n",final_res->score,final_res->offset,final_res->k);
					
}






void findBestCuda(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk, int start, int end){
	 char** t_mk ;
	 char* s1=0 ,*s2=0;
		 struct ms_results **temp_res;
		 struct ms_results *final_res;
		 
		double* scores; 
		int* offsetsArr , *kArr;
int i;
int l1 = strlen(*seq1);
	int l2 = strlen(*seq2);
	int max_offset = l1 - l2;
	
	cudaError_t cudaStatus;


//seq2 total size is maximum 2000 -> this part size will be max 1000 less then the maximum threads number
int size = end-start +1;
		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(0);

		cudaStatus = cudaMalloc((void***) &t_mk, size * sizeof(char*));
		cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
		cudaStatus = cudaMalloc( (void**)&scores, size * sizeof(double));
		cudaStatus = cudaMalloc( (void**)&offsetsArr, size * sizeof(int));
		cudaStatus = cudaMalloc( (void**)&kArr, size * sizeof(int));
		cudaStatus = cudaMalloc((void***) &temp_res, size * sizeof(struct ms_results*));
		cudaStatus = cudaMalloc((void**) &final_res,  sizeof(struct ms_results*));

		cudaStatus = cudaMalloc( (void**)&s1, SEQ1_MAX_LENGTH * sizeof(char));
		cudaStatus = cudaMalloc( (void**)&s2, SEQ2_MAX_LENGTH * sizeof(char));
				
						//cudaStatus = cudaMalloc(t_mk[0], (l2+1) * sizeof(char));
				
				for( i = 0 ; i < size; i++){
			//	printf("%d\n",i);
			//temp_res[i];
			//cudaStatus = cudaMalloc((void**) &(temp_res[i]),  sizeof(struct ms_results));
		//cudaStatus = cudaMalloc((void**)&t_mk[i], SEQ2_MAX_LENGTH *sizeof(char));
		//printf("%d\n",i);
				}
				
			// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(s1, *seq1, SEQ1_MAX_LENGTH * sizeof(char),
			cudaMemcpyHostToDevice);
					
			cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "copy seq1 failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
			cudaStatus = cudaMemcpy(s2, *seq2, SEQ2_MAX_LENGTH * sizeof(char),
			cudaMemcpyHostToDevice);
			cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "copy seq2 failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
			
		printf("%lf,%lf,%lf,%lf,\n",fd.w1,fd.w2,fd.w3,fd.w4);
addCalculateKernel<<<1, size>>>(s1,s2,t_mk,temp_res,fd.w1,fd.w2,fd.w3,fd.w4,max_offset ,scores,offsetsArr , kArr,start);
// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}


	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr,
				"cudaDeviceSynchronize returned error code %d after launching addCalculateKernel!\n%s\n",
				cudaStatus,cudaGetErrorString(cudaStatus));
		goto Error;
	}
	
findBest<<<1, 1>>>(final_res,temp_res , size ,scores,offsetsArr , kArr);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr,
				"cudaDeviceSynchronize returned error code %d after launching findBest!\n%s\n",
				cudaStatus,cudaGetErrorString(cudaStatus));
		goto Error;
	}


int temp* = (int*)malloc(sizeof(int));

printf("try to copy\n");
	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(&res->score, &final_res->score,  sizeof(double),
			cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy  1 failed!\n %s\n",cudaGetErrorString(cudaStatus));
		goto Error;
	}
	
			cudaStatus = cudaMemcpy(temp, &final_res->k,  sizeof(int),
			cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 3  failed!\n %s\n",cudaGetErrorString(cudaStatus));
		goto Error;
	}
			cudaStatus = cudaMemcpy(temp, &final_res->offset,  sizeof(int),
			cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 4 failed!\n %s\n",cudaGetErrorString(cudaStatus));
		goto Error;
	}
	printf("copy ok\n");
						printf("after copy best score = %lf , offset = %d ,k = %d\n",res->score,res->offset,res->k);
	

			cudaFree(t_mk);
			cudaFree(s1);
			cudaFree(s2);
			cudaFree(temp_res);
			cudaFree(final_res);
			cudaFree(kArr);
			cudaFree(offsetsArr);
			cudaFree(scores);



	Error: 
cudaFree(t_mk);
			cudaFree(s1);
			cudaFree(s2);
			cudaFree(temp_res);
cudaFree(final_res);
	cudaFree(kArr);
			cudaFree(offsetsArr);
			cudaFree(scores);


	
			
}

