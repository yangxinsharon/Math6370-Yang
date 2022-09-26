/* Sharon Yang
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include "vec2d_b.hpp"


// This file implements the operations defined in the vec2d_b class


///// General class routines /////

// constructor (initializes values to 0.0)
vec2d_b::vec2d_b(long int m, long int n) {
  // if m or n is illegal, create an empty vector
  if (m < 1 || n < 1) {
    this->length = 0;
    this->nrow = 0;
    this->ncol = 0;
    this->data = NULL;
  } else {
    this->length = m*n;
    this->nrow = m;
    this->ncol = n;
    this->data = new double [m*n];
  }
}


// destructor (frees space associated with a vec2d_b)
vec2d_b::~vec2d_b() {
  if (this->data!=NULL)
    delete[] this->data;
  this->length = 0;
  this->nrow = 0;
  this->ncol = 0;
}


// write myself to stdout
void vec2d_b::Write() const {
  // throw exception if data array isn't allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d_b::Write error, empty data array" );

  // print data to screen
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      std::cout << std::setprecision(17) << "  " << this->data[i*this->ncol+j];
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


// write myself to a file
void vec2d_b::Write(const char *outfile) const {
  // throw exception if data array isn't allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d_b::Write error, empty data array" );

  // throw exception if 'outfile' is empty
  if (strlen(outfile) < 1)
    throw std::invalid_argument( "vec2d_b::Write error, empty outfile" );

  // open output file
  std::ofstream fstr;
  fstr.open(outfile);
  if (!fstr.is_open())
    throw std::invalid_argument( "vec2d_b::Write error, unable to open file for writing" );

  // print data to file
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      fstr << std::setprecision(17) << "  " << this->data[i*this->ncol+j];
    fstr << std::endl;
  }
  fstr << std::endl;

  // close output file
  fstr.close();
}



///// Arithmetic operations defined on a given vec2d_b /////

// x = a*y + b*z
void vec2d_b::LinearSum(double a, const vec2d_b& y, double b, const vec2d_b& z) {
  // check that array sizes match
  if (y.nrow != this->nrow || z.nrow != this->nrow || y.ncol != this->ncol || z.ncol != this->ncol)
    throw std::invalid_argument( "vec2d_b::LinearSum error, vector sizes do not match" );

  // check that data is not NULL
  if (this->data == NULL || y.data == NULL || z.data == NULL)
    throw std::invalid_argument( "vec2d_b::LinearSum error: empty data array" );

  // perform operation
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      this->data[i*this->ncol+j] = a*y.data[i*this->ncol+j] + b*z.data[i*this->ncol+j];
  }
}


//   x = x*a  (scales my data by scalar a)
void vec2d_b::Scale(double a) {
  // check that data is not NULL
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d_b::Scale error: empty data array" );

  // perform operation
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      this->data[i*this->ncol+j] *= a;
  }
}


//   x = y  (copies y into x)
void vec2d_b::Copy(const vec2d_b& y) {
  // check that array sizes match
  if (y.nrow != this->nrow || y.ncol != this->ncol)
    throw std::invalid_argument( "vec2d_b::Copy error, vector sizes do not match" );

  // check that data is not NULL
  if (this->data == NULL || y.data == NULL)
    throw std::invalid_argument( "vec2d_b::Copy error: empty data array" );

  // perform operation
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      this->data[i*this->ncol+j] = y.data[i*this->ncol+j];
  }
}


//   x = a  (sets all entries of x to the scalar a)
void vec2d_b::Constant(double a) {
  // check that data is not NULL
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d_b::Copy error: empty data array" );

  // perform operation and return
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      this->data[i*this->ncol+j] = a;
  }
}


///// scalar quantities derived from vectors /////

// min x_i
double vec2d_b::Min() {
  // check that my data is allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d_b::Min error: empty data array" );

  // perform operation and return
  double mn = this->data[0];
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      mn = std::min(mn, this->data[i*this->ncol+j]);
  }
  return mn;
}


// max x_i
double vec2d_b::Max() {
  // check that my data is allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d_b::Max error: empty data array" );

  // perform operation and return
  double mx = this->data[0];
  for (long int i=0; i<this->nrow; i++) {
    for (long int j=0; j<this->ncol; j++)
      mx = std::max(mx, this->data[i*this->ncol+j]);
  }
  return mx;
}



///// independent constructor routines /////

// create a vector of linearly spaced data
vec2d_b Linspace(double a, double b, long int m, long int n) {
  vec2d_b *x = new vec2d_b(m,n);
  double *xd = x->GetData();
  double h = (b-a)/(m*n-1);
  for (long int i=0; i<m; i++) {
    for (long int j=0; j<n; j++)
      xd[i*n+j] = a + h*j + h*n*i;
  }
  return *x;
}


// create a vector of uniformly-distributed random data
vec2d_b Random(long int m, long int n) {
  vec2d_b *x = new vec2d_b(m,n);
  double *xd = x->GetData();
  for (long int i=0; i<m; i++) {
    for (long int j=0; j<n; j++)
      xd[i*n+j] = random() / (pow(2.0,31.0) - 1.0);
  }
  return *x;
}


///// independent arithmetic routines /////

// dot-product of x and y
double Dot(const vec2d_b& x, const vec2d_b& y) {
  // check that array sizes match
  if (y.Length() != x.Length() || y.nRow() != x.nRow())
    throw std::invalid_argument( "vec2d_b::Dot error, vector sizes do not match" );

  // perform operation and return
  double sum = 0.0;
  for (long int i=0; i<x.nRow(); i++) {
    for (long int j=0; j<x.nCol(); j++)
      sum += x[i*x.nCol()+j] * x[i*x.nCol()+j];
  }
  return sum;
}



// ||x||_2
double TwoNorm(const vec2d_b& x) {
  double sum = 0.0;
  for (long int i=0; i<x.nRow(); i++)
    for (long int j=0; j<x.nCol(); j++)
      sum += x[i*x.nCol()+j] * x[i*x.nCol()+j];
  return sqrt(sum);
}


// ||x||_RMS
double RmsNorm(const vec2d_b& x) {
  double sum = 0.0;
  for (long int i=0; i<x.nRow(); i++) {
    for (long int j=0; j<x.nCol(); j++)
      sum += x[i*x.nCol()+j]*x[i*x.nCol()+j];
  }
  return sqrt(sum/x.Length());
}


// ||x||_infty
double MaxNorm(const vec2d_b& x) {
  double mx = 0.0;
  for (long int i=0; i<x.nRow(); i++) {
    for (long int j=0; j<x.nCol(); j++)
      mx = std::max(mx, std::abs(x[i*x.nCol()+j]));
  }
  return mx;
}
