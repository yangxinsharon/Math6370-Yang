/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gl_variables.h"
#include "gl_constants.h"
#include "mpi.h"
#include <string.h>

/* Prototypes */
int *matvec(double *alpha, double *x, double *beta, double *y);
void comp_soleng_wrapper(double soleng);
void comp_source_wrapper();
void comp_pot(const double* xvct, double *atmchr, double *chrpos, double *ptl, 
	double *tr_xyz,double *tr_q, double *tr_area, int nface, int nchr);
void comp_source( double* bvct, double *atmchr, double *chrpos, 
	double *tr_xyz, double *tr_q, int nface, int nchr);

 
void matvecmul(const double *x, double *y, double *q, int nface, 
	double *tr_xyz, double *tr_q, double *tr_area, double alpha, double beta) {
	/* declarations for mpi */
	int is, ie, i, j;
	int ierr, numprocs, myid;

	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_size = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	printf(" MYID = %i\n",myid);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_rank = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

    double pre1=0.50*(1.0+eps); /* const eps=80.0 */
    double pre2=0.50*(1.0+1.0/eps);
	printf(" nface = %i\n",nface);
	/* determine this processor's interval */
  	is = ((int) (1.0*nface/numprocs))*myid;
  	ie = ((int) (1.0*nface/numprocs))*(myid+1);
  	if (myid == numprocs-1)  ie = nface;

	printf(" is = %i\n",is);
	printf(" ie = %i\n",ie);
	printf(" nface = %i\n",nface);


  	int chunk=ie-is;
  	printf(" chunk = %i\n",chunk);
  	double *sbuf_y1, *sbuf_y2;
  	sbuf_y1 = (double *) calloc(chunk, sizeof(double));
	sbuf_y2 = (double *) calloc(chunk, sizeof(double));

    for (i=is; i<ie; i++) { 

    	double tp[3] = {tr_xyz[3*i], tr_xyz[3*i+1], tr_xyz[3*i+2]};
		double tq[3] = {tr_q[3*i], tr_q[3*i+1], tr_q[3*i+2]};

		double peng[2] = {0.0, 0.0};
		for (j=0; j<nface; j++) {

        	if (j != i) {
				double sp[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
				double sq[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
				double r_s[3] = {sp[0]-tp[0], sp[1]-tp[1], sp[2]-tp[2]};
				double sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
				double rs = sqrt(sumrs);
				double irs = 1.0/sqrt(sumrs) ; //rsqrt(sumrs);
				double G0 = one_over_4pi;
				G0 = G0*irs;
				double kappa_rs = kappa*rs;
				double exp_kappa_rs = exp(-kappa_rs);
				double Gk = exp_kappa_rs*G0;
	
				double cos_theta = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
				double cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
	
				double tp1 = G0*irs;
				double tp2 = (1.0+kappa_rs) * exp_kappa_rs;
	
				double G10 = cos_theta0*tp1;
				double G20 = tp2*G10;
	
				double G1 = cos_theta*tp1;
				double G2 = tp2*G1;
	
				double dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
				double G3 = (dot_tqsq-3.0*cos_theta0*cos_theta) * irs*tp1;
				double G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
				double L1 = G1-eps*G2;								// K2
				double L2 = G0-Gk;									// K1
				double L3 = G4-G3;									// K4
				double L4 = G10-G20/eps; //fdivide(G20,eps);		// K3
	
				double peng_old[2] = {x[j], x[j+nface]};
				double area = tr_area[j];
				peng[0] = peng[0] + (L1*peng_old[0] + L2*peng_old[1]) * area;
				peng[1] = peng[1] + (L3*peng_old[0] + L4*peng_old[1]) * area;
        	}
		}
		// y[i] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
  		// y[nface+i] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
  		sbuf_y1[i-myid*chunk] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
  		sbuf_y2[i-myid*chunk] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
	}

	// double ftime = MPI_Wtime();
	// double looptime = ftime-stime;
	// printf("looptime = %f\n",looptime);
	// printf("chunk = %i\n",chunk);
	
	// double *send_buf1, *send_buf2;
	// send_buf1 = (double *) calloc(chunk, sizeof(double));
	// send_buf2 = (double *) calloc(chunk, sizeof(double));
  	// for (i=0; i<chunk; i++) {
	// 	send_buf1[i]=y[is+i];
	// 	send_buf2[i]=y[nface+is+i];
  	// }


	double *rece_buf1,*rece_buf2;
	rece_buf1 = (double *) calloc(nface, sizeof(double));
	rece_buf2 = (double *) calloc(nface, sizeof(double));

	int irecv[numprocs];
	for (i=0; i<numprocs-1; i++) {
		irecv[i] = ((int) (1.0*nface/numprocs));
	}
	irecv[numprocs-1]=nface-((int) (1.0*nface/numprocs))*(numprocs-1);

	int idisp[numprocs];
	int sumtmp=0;
	for (i = 0; i <numprocs; i++) {
		idisp[i] = 0 + sumtmp;
		sumtmp = sumtmp+irecv[i];
	}


	// ierr = MPI_Allgather(sbuf_y1, chunk, MPI_DOUBLE, rece_buf1, chunk, MPI_DOUBLE, MPI_COMM_WORLD);
	ierr = MPI_Allgatherv(sbuf_y1, chunk, MPI_DOUBLE, rece_buf1, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather1 = %i\n",ierr);
  	}

	// ierr = MPI_Allgather(sbuf_y2, chunk, MPI_DOUBLE, rece_buf2, chunk, MPI_DOUBLE, MPI_COMM_WORLD);
	ierr = MPI_Allgatherv(sbuf_y2, chunk, MPI_DOUBLE, rece_buf2, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather2 = %i\n",ierr);
  	}


	// ierr = MPI_Allgather(y+myid*chunk+nface, chunk, MPI_DOUBLE, rece_buf1+nface, chunk, MPI_DOUBLE, MPI_COMM_WORLD);
  	// if (ierr != MPI_SUCCESS) {
  	//    	printf("Error in MPI_Allgather2 = %i\n",ierr);
  	// } y+myid*chunk+nface

	// ierr = MPI_Allgather(send_buf2, chunk, MPI_DOUBLE, rece_buf2, chunk, MPI_DOUBLE, MPI_COMM_WORLD);
  	// if (ierr != MPI_SUCCESS) {
  	//    	printf("Error in MPI_Allgather2 = %i\n",ierr);
  	// }

	// ftime = MPI_Wtime();
	// double commtime = ftime-stime;
	// printf("commtime = %f\n",commtime);
	
	// printf("nface = %i\n",nface);


	// stime = MPI_Wtime();
  	for (i=0; i<nface; i++) {
  		y[i] = rece_buf1[i];
  		// y[nface+i] = rece_buf1[i+nface];
  		y[nface+i] = rece_buf2[i];
  	}

	// memcpy(y,rece_buf1,sizeof(y));
  	free(sbuf_y1);
  	free(sbuf_y2);
  	free(rece_buf1);
  	free(rece_buf2);

	// ftime = MPI_Wtime();
	// double cpytime = ftime -stime;
	// printf("cpytime = %f\n",cpytime);



  	// printf("scount %i\n",scount);
  	// printf("rece_buf[N*numprocs] = %f\n",rece_buf[N*numprocs]);
  	// printf("y[1000] before = %f\n",y[1000]);

  	// memcpy(y,rece_buf+N*numprocs*sizeof(double),sizeof(y));
  	// printf("y[1000] after = %f\n",y[1000]);
	// printf("y[scount] = %f\n",y[scount]);
}



/* This subroutine wraps the matrix-vector multiplication */
int *matvec(double *alpha, double *x, double *beta, double *y) {
    matvecmul(x, y, tr_q, nface, tr_xyz, tr_q, tr_area, *alpha, *beta);
    // 
    return NULL;
}

/* This subroutine wraps the solvation energy computation */
void comp_soleng_wrapper(double soleng) {
    int i;
	double *chrptl;
	double units_para = 2.0;
    units_para = units_para *units_coef;
    units_para = units_para*pi;

	if ((chrptl=(double *) malloc(nface*sizeof(double)))==NULL) {
		printf("error in allcating chrptl");
	}

	comp_pot(xvct, atmchr, chrpos, chrptl, tr_xyz, tr_q, tr_area, nface, nchr);
	soleng=0.0;
	for (i=0; i<nface; i++) soleng = soleng+chrptl[i];
	soleng = soleng*units_para;
	printf("solvation energy = %f kcal/mol\n",soleng);
}



/* This subroutine calculates the element-wise potential */
void comp_pot(const double* xvct, double *atmchr, double *chrpos, double *ptl, 
	double *tr_xyz, double *tr_q, double *tr_area, int nface, int nchr) {
	int i, j;
    double sumrs, irs, rs, G0, Gk, kappa_rs, exp_kappa_rs;
    double cos_theta, G1, G2, L1, L2, tp1, tp2;
	for (j=0; j<nface; j++) {
    	ptl[j] = 0.0;
		double r[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
		double v[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
    	for (i=0; i<nchr; i++) {
        	double s[3] = {chrpos[3*i], chrpos[3*i+1], chrpos[3*i+2]};
			double r_s[3] = {r[0]-s[0], r[1]-s[1], r[2]-s[2]};
			sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
			rs = sqrt(sumrs);
			irs = 1.0/sqrt(sumrs);

        	G0 = one_over_4pi;
        	G0 = G0*irs;
        	kappa_rs = kappa*rs;
        	exp_kappa_rs = exp(-kappa_rs);
        	Gk = exp_kappa_rs*G0;

        	cos_theta = (v[0]*r_s[0]+v[1]*r_s[1]+v[2]*r_s[2]) * irs;

        	tp1 = G0*irs;
        	tp2 = (1.0+kappa_rs)*exp_kappa_rs;

        	G1 = cos_theta*tp1;
        	G2 = tp2*G1;

        	L1 = G1-eps*G2;
        	L2 = G0-Gk;

      		ptl[j] = ptl[j] + atmchr[i] * (L1*xvct[j]+L2*xvct[nface+j]) * tr_area[j];
		}
    }
}

/* This subroutine wraps the solvation energy computation */
void comp_source_wrapper() {
    comp_source(bvct, atmchr, chrpos, tr_xyz, tr_q, nface, nchr);
}


/* This subroutine calculates the source term of the integral equation */
/* atmchr=atom charge   chrpos=charge position */
/* bvct be located at readin.c */
void comp_source( double* bvct, double *atmchr, double *chrpos, 
	double *tr_xyz,double *tr_q, int nface, int nchr) {
	int i, j;
	double sumrs, cos_theta, irs, G0, G1, tp1;
	for (i=0; i<nface; i++) {
        bvct[i] = 0.0;
        bvct[i+nface] = 0.0;
        for (j=0; j<nchr; j++) {
            double r_s[3] = {chrpos[3*j]-tr_xyz[3*i], chrpos[3*j+1]-tr_xyz[3*i+1], 
            	chrpos[3*j+2]-tr_xyz[3*i+2]};
			sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2]; 
            cos_theta = tr_q[3*i]*r_s[0] + tr_q[3*i+1]*r_s[1] + tr_q[3*i+2]*r_s[2];
			irs = 1.0/sqrt(sumrs) ;//rsqrt(sumrs);//returns reciprocal square root of scalars and vectors.
            cos_theta = cos_theta*irs;
            G0 = one_over_4pi;//constant
            G0 = G0*irs;
            tp1 = G0*irs;
            G1 = cos_theta*tp1;
            bvct[i] = bvct[i]+atmchr[j]*G0;
            bvct[nface+i] = bvct[nface+i]+atmchr[j]*G1;
        }
    }
}
