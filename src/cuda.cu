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


__device__ //compare 2 sequence. return score as float
double compereSeq1AndSeq2Cuda(char **seq1, char **seq2,
		int offset ,double w1,double w2,double w3,double w4	) {

	int l1 = _strlen(*seq1);
	int l2 = _strlen(*seq2);

	//printf("l2 #2 = %d\n",l2);
	//printf("hi from compereSeq1AndSeq2\n");

	int numberOfStars = 0, numberOfColons = 0, numberOfPoints = 0,
			numberOfSpaces = 0, i; //init var for calculate the score

	if (l1 == l2) {
		for (i = 0; i < l1; i++) {
			switch ((*seq2)[i] == '-' ?
					'-' : compareTowCharsCuda((*seq2)[i], (*seq1)[i])) { //in case seq1 and seq2 have the same length
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
				switch ((*seq2)[i] == '-' ?
						'-' :
						compareTowCharsCuda((*seq2)[i], (*seq1)[offset + i])) { // get sing for each pair of chars with specific offset
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
int j = 1000 , k=0 ;
printf("kkkk\n");
do{
if(k == i){
printf("**\n");
//des[j++] = '-';
printf("****\n");}

//des[j++] ='-'; //from[k];
printf("%c\n",from[k++]);
}while(j-- != 0);

}


__global__ void addCalculateKernel(char **seq1, char **seq2, char** t_mk,  struct ms_results **temp_res,
double w1,double w2,double w3,double w4 ,int max_offset	) {
	int tid = threadIdx.x;
	double d;
printf("ppppppp\n");
mutantSequenceCuda(*seq2, t_mk[tid], tid);
printf("%s\n",t_mk[tid]);
	for (int i = 0; i < max_offset - 1; i++) {   //all offsets options

				d = compereSeq1AndSeq2Cuda(seq1, &t_mk[tid], i , w1,w2,w3,w4); // compare seq1 and current MS(j) and offset = i
				if (d > temp_res[tid]->score) {
					temp_res[tid]->score = d;
					temp_res[tid]->offset = i;
					temp_res[tid]->k = tid;
				}
			}
			
		
}




__global__ void findBest( struct ms_results* final_res,	 struct ms_results** temp_res ,int size ){
					printf("wtf\n");
					final_res->score = -1;
					final_res->offset = -1;
					final_res->k = -1;
					int i;
					for( i = 0 ; i < size;i++){
					if(
					final_res->score >temp_res[i]->score ){
					
					final_res->score = (temp_res[i])->score ;
					final_res->offset = (temp_res[i])->offset;
					final_res->k = (temp_res[i])->k;
					}
					}
					
}






void findBestCuda(char **seq1, char **seq2, struct file_data fd,
		struct ms_results *res, char **mk, int start, int end){
		printf("cuda 1\n");
	 char** t_mk ;
	 char* s1 ,*s2;
		 struct ms_results **temp_res;
		 struct ms_results *final_res;
int i;
int l1 = strlen(*seq1);
	int l2 = strlen(*seq2);
	int max_offset = l1 - l2;
	
	cudaError_t cudaStatus;

		printf("cuda 2\n");

//seq2 total size is maximum 2000 -> this part size will be max 1000 less then the maximum threads number
int size = end-start +1;
		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(0);

		cudaStatus = cudaMalloc((void**) &t_mk, size * sizeof(char*));
		cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}
		
		cudaStatus = cudaMalloc((void**) &temp_res, size * sizeof(char*));
		cudaStatus = cudaMalloc((void**) &final_res,  sizeof(char*));

			cudaStatus = cudaMalloc((void**) &s1, SEQ1_MAX_LENGTH * sizeof(char**));
				cudaStatus = cudaMalloc((void**) &s2, SEQ2_MAX_LENGTH * sizeof(char**));
						printf("cuda 3\n");
				
						//cudaStatus = cudaMalloc((void**) &t_mk[0], (l2+1) * sizeof(char));
				
				for( i = 0 ; i < size; i++){
			//	printf("%d\n",i);
		cudaStatus = cudaMalloc(t_mk[i], (l2+1) * sizeof(char));
		//printf("%d\n",i);
				}
						printf("cuda 4\n");
				
			// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(s1, *seq1, SEQ1_MAX_LENGTH * sizeof(char*),
			cudaMemcpyHostToDevice);
					printf("cuda 5\n");
			
			cudaStatus = cudaMemcpy(s2, *seq2, SEQ2_MAX_LENGTH * sizeof(char*),
			cudaMemcpyHostToDevice);
					printf("cuda 6\n");
			
			
	//	printf("%s\n",s2);
addCalculateKernel<<<1, size>>>(&s1,&s2,t_mk,temp_res,fd.w1,fd.w2,fd.w3,fd.w4,max_offset);
// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n",
				cudaGetErrorString(cudaStatus));
		goto Error;
	}

		printf("cuda 7\n");

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr,
				"cudaDeviceSynchronize returned error code %d after launching addCalculateKernel!\n%s\n",
				cudaStatus,cudaGetErrorString(cudaStatus));
		goto Error;
	}
	
findBest<<<1, 1>>>(final_res,temp_res , size);

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

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(res, final_res,  sizeof(struct ms_results),
			cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
for(int i = 0 ; i < size; i++)
cudaFree(t_mk[i]);
			cudaFree(t_mk);
			cudaFree(s1);
			cudaFree(s2);
			cudaFree(temp_res);
			cudaFree(final_res);



	Error: 
	for(int i = 0 ; i < size; i++)
cudaFree(t_mk[i]);
cudaFree(t_mk);
			cudaFree(s1);
			cudaFree(s2);
			cudaFree(temp_res);
cudaFree(final_res);


	
			
}

