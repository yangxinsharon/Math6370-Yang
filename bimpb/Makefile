###################################################################
#  Makefile for screened coulombic potential computation
#
#  Xin (Sharon) Yang
#  SMU Mathematics
#  Math 6370 Project: bimpb with Weihua Geng, Jiahui Chen
###################################################################

# compilers
# CC = clang # on Yang's MacOS
CC = gcc # on Linux, i.e., HPC-M2
# NVCC = nvcc
MPICC = mpicc

# flags
CFLAGS = -c -O2
# NVCCFLAGS = -c -O3 -arch=sm_61 -use_fast_math	# for -arch=sm_xx check http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list
OMPFLAGS = -fopenmp

# makefile targets
all: bimpb.exe bimpb_omp.exe #bimpb_mpi.exe #bimpb_cuda.exe 

bimpb.exe: main.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec.o pp_timer.o
	$(CC) -lm -p -o $@ $^
# bimpb_cuda.exe: main_cuda.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_cuda.o pp_timer.o
# 	$(NVCC) -o -lm -o $@ $^
bimpb_omp.exe: main.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_omp.o pp_timer.o
	$(CC) $(OMPFLAGS) -lm -p -o $@ $^
# bimpb_mpi.exe: main_mpi.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_mpi.o pp_timer.o
# 	$(MPICC) -lm -o $@ $^
pp_timer.o: pp_timer.c
	$(CC) $(CFLAGS) pp_timer.c
readin.o: readin.c
	$(CC) $(CFLAGS) readin.c
matvec.o: matvec.c
	$(CC) $(CFLAGS) matvec.c
main.o: main.c
	$(CC) $(CFLAGS) main.c
# matvec_cuda.o: cuda/matvec_cuda.cu
# 	$(NNCC) $(NVCCFLAGS) cuda/matvec_cuda.cu
# main_cuda.o: cuda/main_cuda.c
# 	$(NVCC) $(NVCCFLAGS) cuda/main_cuda.c
matvec_omp.o: matvec_omp.c
	$(CC) $(CFLAGS) $(OMPFLAGS) matvec_omp.c
# main_omp.o: main_omp.c
# 	$(CC) $(CFLAGS) $(OMPFLAGS) main_omp.c
# matvec_mpi.o: matvec_mpi.c
# 	$(MPICC) $(CFLAGS) matvec_mpi.c
# main_mpi.o: main_mpi.c
# 	$(MPICC) $(CFLAGS) main_mpi.c
d_sign.o: gmres/d_sign.c
	$(CC) $(CFLAGS) gmres/d_sign.c
daxpy.o: gmres/daxpy.c
	$(CC) $(CFLAGS) gmres/daxpy.c
dcopy.o: gmres/dcopy.c
	$(CC) $(CFLAGS) gmres/dcopy.c
ddot.o: gmres/ddot.c
	$(CC) $(CFLAGS) gmres/ddot.c
dgemv.o: gmres/dgemv.c
	$(CC) $(CFLAGS) gmres/dgemv.c
dnrm2.o: gmres/dnrm2.c
	$(CC) $(CFLAGS) gmres/dnrm2.c
drot.o: gmres/drot.c
	$(CC) $(CFLAGS) gmres/drot.c
drotg.o: gmres/drotg.c
	$(CC) $(CFLAGS) gmres/drotg.c
dscal.o: gmres/dscal.c
	$(CC) $(CFLAGS) gmres/dscal.c
dtrsv.o: gmres/dtrsv.c
	$(CC) $(CFLAGS) gmres/dtrsv.c
gmres.o: gmres/gmres.c
	$(CC) $(CFLAGS) gmres/gmres.c
clean: 
	\rm -rf *.o *.out
realclean: clean
	\rm -rf *.exe test_proteins/*.face test_proteins/*.vert *.face *.vert *~

####### End of Makefile #######
