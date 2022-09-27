/* FILE: omp_environment.cpp
   DESCRIPTION:
      OpenMP Example - Get environment information - C++ Version
   AUTHOR: Blaise Barney  5/99
   UPDATED: Daniel R. Reynolds (updated to C++), 1/13/2013 
   UPDATED: Sharon Yang, 09/27/2022 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

  // local variables
  int nthreads, tid, procs, maxt, inpar, dynamic, nested;

  // Start parallel region
# pragma omp parallel private(nthreads, tid)
  {

    // Obtain thread number
    tid = omp_get_thread_num();

    // Only master thread does this
    if (tid == 0) {
      printf("Thread %i getting environment info...\n", tid);

      // Get environment information
      procs    = omp_get_num_procs();
      nthreads = omp_get_num_threads();
      maxt     = omp_get_max_threads();
      inpar    = omp_in_parallel();
      dynamic  = omp_get_dynamic();
      nested   = omp_get_nested();

      // Print environment information
      printf("The number of processors available [%i]\n", procs);
      printf("The number of threads being used [%i]\n", nthreads);
      printf("The maximum number of threads available [%i]\n", maxt);
      printf("If you are in a parallel region [%i]\n", inpar);
      printf("If dynamic threads are enabled [%i]\n", dynamic);
      printf("If nested parallelism is supported [%i]\n", nested);
    }

  }  // end parallel region  
  return 0;
}  // end main
