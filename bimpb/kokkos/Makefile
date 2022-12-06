# makefile for screened coulombic potential computation with kokkos
KOKKOS_PATH = ${HOME}/Kokkos/kokkos
KOKKOS_DEVICES = "Cuda,OpenMP"
EXE_NAME = "bimpb_kokkos"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "Volta70"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "BDW"
endif

# CXX = g++
CC = gcc
flag = -c -O2
AR = ar rc
#-lstdc++
all: libgmres.a $(EXE)

# matvec_kokkos.o: matvec_kokkos.cpp
# 	$(CXX) $(flag) matvec_kokkos.cpp 
main_kokkos.o: main_kokkos.cpp
	$(CXX) $(flag) main_kokkos.cpp
# readin.o: readin.cpp
# 	$(CXX) $(flag) readin.cpp

matvec_kokkos.o: matvec_kokkos.c
	$(CC) $(flag) matvec_kokkos.c
# main_kokkos.o: main_kokkos.c
# 	$(CC) $(flag) main_kokkos.c
readin.o: readin.c
	$(CC) $(flag) readin.c

pp_timer.o: pp_timer.c
	$(CC) $(flag) pp_timer.c 

d_sign.o: gmres/d_sign.c
	$(CC) $(flag) gmres/d_sign.c
daxpy.o: gmres/daxpy.c
	$(CC) $(flag) gmres/daxpy.c
dcopy.o: gmres/dcopy.c
	$(CC) $(flag) gmres/dcopy.c
ddot.o: gmres/ddot.c
	$(CC) $(flag) gmres/ddot.c
dgemv.o: gmres/dgemv.c
	$(CC) $(flag) gmres/dgemv.c
dnrm2.o: gmres/dnrm2.c
	$(CC) $(flag) gmres/dnrm2.c
drot.o: gmres/drot.c
	$(CC) $(flag) gmres/drot.c
drotg.o: gmres/drotg.c
	$(CC) $(flag) gmres/drotg.c
dscal.o: gmres/dscal.c
	$(CC) $(flag) gmres/dscal.c
dtrsv.o: gmres/dtrsv.c
	$(CC) $(flag) gmres/dtrsv.c
gmres.o: gmres/gmres.c
	$(CC) $(flag) gmres/gmres.c

libgmres.a: d_sign.o daxpy.o dcopy.o ddot.o dgemv.o dnrm2.o drot.o drotg.o dscal.o dtrsv.o gmres.o
	$(AR) libgmres.a $@ $?
$(EXE): main_kokkos.o readin.o matvec_kokkos.o pp_timer.o
	$(CXX) -static -o $(EXE) *.o -lm -L. -lgmres

clean: 
	\rm -rf *.o
realclean: clean
	\rm -rf *.exe *.a ../test_proteins/*.face ../test_proteins/*.vert *.face *.vert *~
