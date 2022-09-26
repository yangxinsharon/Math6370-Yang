/* Sharon Yang
   SMU Mathematics
   Math 4370 / 6370 */

#ifndef vec2d_b_DEFINED__
#define vec2d_b_DEFINED__

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>


// This defines a simple arithmetic vector class
class vec2d_b {

 private:

  ///// Contents /////
  long int length; // nrow * ncol
  long int nrow;
  long int ncol;
  double *data;

 public:

  ///// General class routines /////
  // constructor (initializes values to 0.0)
  vec2d_b(long int m, long int n);

  // destructor
  ~vec2d_b();

  // write myself to stdout
  void Write() const;

  // write myself to a file
  void Write(const char* outfile) const;

  // returns my overall length
  long int Length() const {return length;};

  // returns my number of rows
  long int nRow() const {return nrow;};

  // returns my number of columns
  long int nCol() const {return ncol;};

  // access my data array
  double* GetData() const {return data;};

  // access my data -- access value by location in array,
  // returns handle to value so that it can be changed as well,
  // does not check for legal 'i'
  // double* operator[](long int i) const {return data[i];};
  double& operator()(long int i) {return data[i];};
  double operator[](long int i) const {return data[i];};


  ///// Arithmetic operations defined on a vec2d_b /////

  // in-place operations (x is the calling vec2d_b) -- 0=success, 1=fail
  void LinearSum(double a, const vec2d_b& y,      // x = a*y + b*z
                 double b, const vec2d_b& z);
  void Scale(double a);                         // x = x*a
  void Copy(const vec2d_b& y);                    // x = y
  void Constant(double a);                      // x = a

  // scalar quantites derived from vectors
  double Min();                                 // min x_i
  double Max();                                 // max x_i

};  // end vec2d_b


// independent constructor routines
vec2d_b Linspace(double a, double b, long int m, long int n);
vec2d_b Random(long int m, long int n);

// independent arithmetic routines
double Dot(const vec2d_b& x, const vec2d_b& y);     // sum_i (x_i * y_i)
double TwoNorm(const vec2d_b& x);                 // sqrt(sum_i x_i^2)
double RmsNorm(const vec2d_b& x);                 // sqrt(sum_i x_i^2 / n)
double MaxNorm(const vec2d_b& x);                 // max_i |x_i|

#endif
