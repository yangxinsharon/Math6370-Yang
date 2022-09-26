/* Sharon Yang
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include "vec2d.hpp"

// prototype for Gram-Schmidt routine
int GramSchmidt2d(vec2d X[], int numvectors);


// Example routine to test the vec2d class
int main(int argc, char* argv[]) {

  // create some vecs of size 5x3, and set some entries
  int m = 5; // nrow
  int n = 3; // ncol
  vec2d a(m,n), b(m,n), c(m,n);
  for (int i=0; i<b.nRow(); i++) {
    for (int j=0; j<b.nCol(); j++)
      b[i*b.nRow()+j] = (j+1)*0.1 + n*i*0.1;
  }
  for (int i=0; i<c.nRow(); i++) {
    for (int j=0; j<c.nCol(); j++)
      c[i*c.nRow()+j] = (j+1) + n*i;
  }

  // output to screen
  std::cout << "writing array of zeros:" << std::endl;
  a.Write();
  std::cout << "writing array of 0.1, 0.2, ... , 0.15:" << std::endl;
  b.Write();
  std::cout << "writing array of 1,2,...,15:" << std::endl;
  c.Write();

  // verify that b has length 5x3=15
  if (b.Length() != 15)
    std::cerr << "error: incorrect vector length" << std::endl;

  // access a's data array, update entries, and write to file
  double *dat = a.GetData();
  for (int j=0; j<a.nCol(); j++) {
    dat[0] = 10.0;
    dat[1] = 15.0;
    dat[2] = 20.0;
    dat[3] = 25.0;
    dat[4] = 30.0;
  } 
  a.Write("a_data");
  std::cout << "the new file 'a_data' on disk should have same row entries 10, 15, 20, 25, 30" << std::endl << std::endl;

  // access each entry of a and write to screen
  std::cout << "entries of a, one at a time (via []): should give 10, 10, 10, 15, 15, 15, ..., 30, 30, 30" << std::endl;
  for (int i=0; i<a.nRow(); i++){
    for (int j=0; j<a.nCol(); j++)
      std::cout << "  " << a[i*a.nRow()+j];
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // update one entry of a
  std::cout << "updating the last entry of a to be 31 (via a[4][2]):" << std::endl;
  a[14] = 31.0;
  a.Write();
  a[14] = 30.0;  // reset to original

  // Test arithmetic operators
  std::cout << "Testing vector constant, should all be -1" << std::endl;
  b.Constant(-1.0);
  b.Write();

  std::cout << "Testing vector copy, should be 1, 2, ..., 15" << std::endl;
  a.Copy(c);
  a.Write();

  std::cout << "Testing scalar multiply, should be 5, 10, ..., 75" << std::endl;
  c.Scale(5.0);
  c.Write();

  // create a few 2d vecs of size 10x3
  vec2d X[] = {vec2d(10,3), vec2d(10,3), vec2d(10,3), vec2d(10,3), vec2d(10,3)};

  // fill in the vectors
  for (int i=0; i<10; i++) {
    for (int j=0; j<3; j++) {
      X[0][i*X[0].nRow()+j] = (1.0*j+3*i);
      X[1][i*X[1].nRow()+j] = -5.0 + 1.0*(j+3*i);
      X[2][i*X[2].nRow()+j] = 2.0 + 2.0*(j+3*i);
      X[3][i*X[3].nRow()+j] = 20.0 - 1.0*(j+3*i);
      X[4][i*X[4].nRow()+j] = -20.0 + 1.0*(j+3*i);
    }
  }

  // check the LinearSum routine
  X[0].LinearSum(-2.0, X[1], 1.0, X[2]);
  std::cout << "Testing LinearSum, should be all 12:" << std::endl;
  X[0].Write();

  // check the various scalar output routines
  std::cout << "Testing TwoNorm (should be 3.872983346207417):  " << TwoNorm(b) << std::endl;

  std::cout << "Testing RmsNorm (should be 45.460605656619521):  " << RmsNorm(c) << std::endl;

  std::cout << "Testing MaxNorm (should be 1):  " << MaxNorm(b) << std::endl;

  std::cout << "Testing Min (should be 1):  " << a.Min() << std::endl;

  std::cout << "Testing Max (should be 75):  " << c.Max() << std::endl;

  std::cout << "Testing Dot (should be 1240):  " << Dot(a, c) << std::endl;

  std::cout << "Testing Linspace, should be 0 1 2 ... 14" << std::endl;
  vec2d d = Linspace(0.0, 14.0, 5, 3);
  d.Write();

  std::cout << "Testing Random" << std::endl;
  vec2d f = Random(5,3);
  f.Write();


  /// performance/validity tests (Gram-Schmidt process)
  int vec_sizes[6][2]={{10000, 1000}, {1000, 10000}, {100, 100000}, {10, 1000000}, {100000, 100}, {1000000, 10}};
  for (int i=0; i<6; i++){
    int m = vec_sizes[i][0]; // nrow
    int n = vec_sizes[i][1]; // ncol
  
    std::cout << "Running GramSchmidt2d process" << std::endl;
    vec2d Y[] = {Random(m,n), Random(m,n), Random(m,n), Random(m,n), Random(m,n)};
    std::chrono::time_point<std::chrono::system_clock> stime = std::chrono::system_clock::now();
    if (GramSchmidt2d(Y,5))
      std::cerr << "GramSchmidt2d returned error" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> ftime = std::chrono::system_clock::now();
    std::chrono::duration<double> rtime = ftime-stime;
    std::cout << "  GramSchmidt2d time: " << rtime.count() << std::endl << std::endl;
  
    std::cout << "Resulting vectors should be orthonormal:" << std::endl;
    bool pass = true;
    double tolerance = 1e-12;
    for (int i=0; i<5; i++) {
      if (std::abs(Dot(Y[i], Y[i]) - 1.0) > tolerance) {
        pass = false;
        std::cout << "  <Y[" << i << "],Y[" << i << "]> = " << Dot(Y[i], Y[i]) << std::endl;
      }
      for (int j=i+1; j<5; j++) {
        if (std::abs(Dot(Y[i], Y[j])) > tolerance) {
          pass = false;
          std::cout << "  <Y[" << i << "],Y[" << j << "]> = " << Dot(Y[i], Y[j]) << std::endl;
        }
      }
    }
    if (pass)
      std::cout << "  passed orthonormality check" << std::endl << std::endl;
    else
      std::cout << "  failed orthonormality check" << std::endl << std::endl;
  }

  return 0;
} // end main
