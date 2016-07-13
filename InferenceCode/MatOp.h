/*
  MatOp.cpp

  Allows for exponentiation of matrices by the squaring method
 
  Nuno Rocha Nene   

*/

#ifndef __MatrixPow__MatOp__
#define __MatrixPow__MatOp__

// #include <cblas.h> // NRN

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdio.h>

// GSL
// ===
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h> // for faster matrix multiplication
#include <gsl/gsl_cblas.h>


class MatOp{
public:
    
    // Calculate matrix powers with Open_BLAS routines
    // ==============================================
    void PowMatBLAS(int expo, int PopSize, gsl_matrix *I, gsl_matrix *m);
    
    // Copy matrix
    // ===========
    void Matcopy(int PopSize,double *A, double *B);
    
    // Read in transition matrix
    
    void ReadInData(std::string dataFile,int N,gsl_matrix* m);
    
};

 

#endif  /* defined(__MatrixPow__MatOp__) */


