###################################################################
#  Makefile for screened coulombic potential computation
#
#  Xin (Sharon) Yang
#  SMU Mathematics
#  Math 6370 Project bimpb with Weihua Geng, Jiahui Chen
###################################################################

# compilers
CC = gcc 
# CU = nvcc

# flags
CFLAGS = -c
# CUCFLAGSS = -c -O3 -arch=sm_61 -use_fast_math
# for -arch=sm_xx check http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list
all: bimpb_s.exe #bimpb_c.exe 

# bimpb_c.exe: main_cuda.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec_cuda.o pp_timer.o
# 	$(comp) -o md_cuda.exe *.o -lm

bimpb_s.exe: main.o readin.o d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o matvec.o pp_timer.o
	$(CC) -o md_serial.exe *.o -lm
matvec.o: matvec.c
	$(CC) $(CFLAGS) matvec.c
main.o: main.c
	$(CC) $(CFLAGS) main.c

# matvec_cuda.o: matvec_cuda.cu
# 	$(CC) $(CFLAGS) cuda/matvec_cuda.cu
# main_cuda.o: main_cuda.c
# 	$(CU) $(CFLAGS) cuda/main_cuda.c
pp_timer.o: pp_timer.c
	$(CC) $(CFLAGS) pp_timer.c
readin.o: readin.c
	$(CC) $(CFLAGS) readin.c
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
	\rm -rf *.exe *~
