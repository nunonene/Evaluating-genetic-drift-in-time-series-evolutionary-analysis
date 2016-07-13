//
//  TMatBinomial.h
//  MatrixPow
//
//  Created by Nuno Rocha Nene on 21/07/2014.
//  Copyright (c) 2014 ___nn276___. All rights reserved.
//

#ifndef __MatrixPow__TMatBinomial__
#define __MatrixPow__TMatBinomial__


#include "shared.h"


class TMatBinomial{
    
public:
    
    int PopSize; // Population size
    
    int GridSize; // Grid size for transition matrix
    
    int NumGen; // Number of generations
    
    double Sigma; // Selection coefficient
    
    int NumThreads; // Number of threads for OpenBLAS
    
    void GetParameters (int argc, const char **argv); // Import parameters
    
    double Binomial_pdf (const size_t K,vector<double> &p,vector<int>& n); // Binomial distribution
    
    
    double Binomial_lnpdf (const size_t K, vector<double> &p,vector<int>& n);
};



#endif /* defined(__MatrixPow__TMatBinomial__) */
