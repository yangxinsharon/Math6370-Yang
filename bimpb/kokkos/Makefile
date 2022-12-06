# makefile for screened coulombic potential computation with cuda
CC = gcc
CXX = g++
flag = -c -O2
AR = ar rc
#-lstdc++
all: libgmres.a bimpb_kokkos.exe

matvec_kokkos.o: matvec_kokkos.cpp
	$(CXX) $(flag) matvec_kokkos.cpp 
main_kokkos.o: main_kokkos.cpp
	$(CXX) $(flag) main_kokkos.cpp
readin.o: readin.cpp
	$(CXX) $(flag) readin.cpp

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
bimpb_kokkos.exe: main_kokkos.o readin.o matvec_kokkos.o pp_timer.o
	$(CXX) -static -o bimpb_kokkos.exe *.o -lm -L. -lgmres

clean: 
	\rm -rf *.o
realclean: clean
	\rm -rf *.exe ../test_proteins/*.face ../test_proteins/*.vert *.face *.vert *~
