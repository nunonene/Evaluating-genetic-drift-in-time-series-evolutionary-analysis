//
//  MatOp.h
//  MatrixPow
//
//  Created by Nuno Rocha Nene on 22/07/2014.
//  Copyright (c) 2014 ___nn276___. All rights reserved.
//

#ifndef __MatrixPow__MatOp__
#define __MatrixPow__MatOp__



#include "shared.h"


class MatOp{
    
    
public:
    
    // Calculate matrix powers with Open_BLAS routines
    // ==============================================
    
    void PowMatOpenBLAS(int expo, int PopSize, double *I, double *m);
    
    // Copy matrix
    // ===========
    
    void Matcopy(int PopSize,double *A, double *B);
    
    
    // Save gsl_matrix to file
    // =======================
    
    void SaveToFile(FILE *stream,double *m,int PopSize,int GridSize);
    
};


#endif /* defined(__MatrixPow__MatOp__) */
