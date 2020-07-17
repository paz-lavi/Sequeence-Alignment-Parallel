build:
	mpicxx -fopenmp -c SequenceAlignmentParallel.c -o SequeenceAlignmentParallel.o  
	mpicxx -fopenmp -c comparition.c -o comparition.o 
	mpicxx -fopenmp -c main_func.c -o main_func.o 
	mpicxx -fopenmp -c helpers.c -o helpers.o 
	/usr/local/cuda-10.2/bin/nvcc -I./inc -c cuda.cu -o cuda.o
	mpicxx -fopenmp -o SequeenceAlignmentParallel  SequeenceAlignmentParallel.o comparition.o main_func.o helpers.o cuda.o  /usr/local/cuda-10.2/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./SequeenceAlignmentParallel

run:
	mpiexec -np 2 ./SequeenceAlignmentParallel input.txt

runOn2:
	mpiexec -np 2 -machinefile  /home/linuxu/hostsIP.txt -map-by node ./SequeenceAlignmentParallel input.txt
