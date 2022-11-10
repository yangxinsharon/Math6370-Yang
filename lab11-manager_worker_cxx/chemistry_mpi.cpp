/* Daniel R. Reynolds; Sharon Yang
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "mpi.h"

// Prototypes
void chem_solver(double, double*, double*, double*, double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // declarations
  int n, tag, maxit, its;
  int ierr, numprocs, myid, numsent;
  double *Pbuf, *Sbuf;
  double lam, eps, res, stime, ftime, runtime;
  bool more_work;
  MPI_Status status;
  int iend, sender, numreceived, mycount, ansentry;


  // initialize MPI
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
     std::cerr << "Error in calling MPI_Init\n";
     return 1;
  }

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
     std::cerr << "Error in calling MPI_Comm_size\n";
     return 1;
  }

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
     std::cerr << "Error in calling MPI_Comm_rank\n";
     return 1;
  }

  // 1. set solver input parameters
  maxit = 1000000;
  lam = 1.e-2;
  eps = 1.e-10;
  Pbuf = new double [4];
  Sbuf = new double [3];


  // manager code
  if (myid == 0) { 

    // 2. input the number of intervals
    std::cout << "Enter the number of intervals (0 quits):\n";
    std::cin >> n;
    if (n < 1) {
      return 1;
    }

    // 3. allocate temperature and solution arrays
    double *T = new double[n];
    double *u = new double[n];
    double *v = new double[n];
    double *w = new double[n];

    // 4. set random temperature field, initial guesses at chemical densities
    for (int i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
    for (int i=0; i<n; i++)  u[i] = 0.35;
    for (int i=0; i<n; i++)  v[i] = 0.1;
    for (int i=0; i<n; i++)  w[i] = 0.5;
  
    // 5. start timer
    stime = MPI_Wtime();

    // a counter to know how much work it has yet to send out
    numsent = 0;

    // first send initial tasks to each of the worker node
    for (int i=0; i<n; i++) {
      // fill send buffer
      Pbuf[0] = T[i];
      Pbuf[1] = u[i];
      Pbuf[2] = v[i];
      Pbuf[3] = w[i];
      // send with tag as entry in temperature array
      ierr = MPI_Send(Pbuf, 4, MPI_DOUBLE, i+1, numsent, MPI_COMM_WORLD);
      if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Send = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }
      numsent++;
    } 



    // obtain the workers’ solutions
    for (int i=0; i<n; i++) {

      // receive answers from any process
      ierr = MPI_Recv(Sbuf, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Recv = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }

      // decode the sender and solution entry from status
      sender = status.MPI_SOURCE;
      tag = status.MPI_TAG;

      // store results
      u[tag] = Sbuf[0];
      v[tag] = Sbuf[1];
      w[tag] = Sbuf[2];
      if (numsent < n) {  // send another row
        Pbuf[0] = T[numsent];
        Pbuf[1] = u[numsent];
        Pbuf[2] = v[numsent];
        Pbuf[3] = w[numsent];
        ierr = MPI_Send(Pbuf, 4, MPI_DOUBLE, sender, numsent, MPI_COMM_WORLD);
        if (ierr != MPI_SUCCESS) {
          printf("Error in MPI_Send = %i\n",ierr);
          ierr = MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }
        numsent++;

      // tell senders that work is complete, by sending message of zero size
      } else {
        ierr = MPI_Send(Pbuf, 0, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD);
        if (ierr != 0) {
          printf("Error in MPI_Send = %i\n",ierr);
          ierr = MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }
      } // if numsent
           
    }

    // 7. stop timer
    ftime = MPI_Wtime();
    runtime = ftime - stime;
  
    // 8. output solution time
    std::cout << "     runtime = " << runtime << std::endl;
  
    // 9. free temperature and solution arrays
    delete[] T;
    delete[] u;
    delete[] v;
    delete[] w;




  // worker code
  } else {

    // set up workers' own local variables
    double T, u, v, w;

    // 6. call solver in a while loop
    more_work = true;
    while (more_work) {

      // receive from the manager
      ierr = MPI_Recv(Pbuf, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (ierr != MPI_SUCCESS) {
        std::cerr << "Error in calling MPI_Recv in worker\n";
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }

      tag = status.MPI_TAG;
      std::cerr << "tag is %i\n"<<tag<<std::endl;
      // // check if work is complete
      // ierr = MPI_Get_count(&status, MPI_DOUBLE, &numreceived);
      // if (numreceived == 0) {

      if (tag == 0) {
        more_work = false;
        return 1;
      } else {
        // values received from the manager
        T = Pbuf[0];
        u = Pbuf[1];
        v = Pbuf[2];
        w = Pbuf[3];
        chem_solver(T, &u, &v, &w, lam, eps, maxit, &its, &res); 
        if (res < eps) {
          std::cout << "    i = " << tag << "  its = " << its << std::endl;
        }
        else {
          std::cout << "    error: i=" << tag << ", its=%i" << its << ", res=" << res 
          << ", u=" << u << ", v=" << v << ", w=" << w << std::endl;
          ierr = MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }

        // pack the solution buffer and send it back to the manager
        Sbuf[0] = u;
        Sbuf[1] = v;
        Sbuf[2] = w;
        ierr = MPI_Send(Sbuf, 3, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        if (ierr != 0) {
          std::cerr << "Error in calling MPI_Send in worker\n";
          ierr = MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }
      }
    } // end while loop
  } // end worker


  // finalize MPI
  ierr = MPI_Finalize();

  return 0;
} // end main
