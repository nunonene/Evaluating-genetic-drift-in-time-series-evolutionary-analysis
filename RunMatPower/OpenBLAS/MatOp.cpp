//
//  MatOp.cpp
//  MatrixPow
//
//  Created by Nuno Rocha Nene on 22/07/2014.
//  Copyright (c) 2014 ___nn276___. All rights reserved.
//

#include "MatOp.h"
#include "shared.h"

void MatOp::PowMatOpenBLAS(int expo, int PopSize, double  *I, double *m){ // with BLAS matrix multiplication routines
    
    int n=expo;
    
    //cout << "Stage 1: Allocating auxiliary matrices\n\n";
    
    double *dest;// = (double*) malloc((PopSize+1)*(PopSize+1));
    
    double *src;// = (double*) malloc((PopSize+1)*(PopSize+1));
    
    double *tempMat;// = (double*) malloc((PopSize+1)*(PopSize+1));
    
    
    
    dest = (double*) mkl_malloc((PopSize+1)*(PopSize+1)*sizeof(double),64);
    
    src = (double*) mkl_malloc((PopSize+1)*(PopSize+1)*sizeof(double),64);
    
    tempMat = (double*) mkl_malloc((PopSize+1)*(PopSize+1)*sizeof(double),64);
    
    
    
    
    
    
    Matcopy (PopSize,dest, I);
    
    Matcopy(PopSize,src, m);
    
    
    //cout<< "Calculating Matrix Powers\n\n";
    
    
   
    
    do
    {
       cout<< n << "\n";
        
        if(n % 2 != 0){
            
             //clock_t startTime = clock();
            
         
            
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,PopSize+1,PopSize+1,PopSize+1,
                            1.0, dest,PopSize+1, src,PopSize+1,
                         0.0, tempMat,PopSize+1); // void cblas_dgemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST double                      alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc);
            
            
            //cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds.\n\n" << endl;
            
            Matcopy(PopSize,dest,tempMat);
            
        }
        
        //clock_t startTime = clock();

        
        cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,PopSize+1,PopSize+1,PopSize+1,
                     1.0, src,PopSize+1, src,PopSize+1,
                     0.0, tempMat,PopSize+1);
        
       // cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds.\n\n" << endl;
        
        Matcopy(PopSize,src,tempMat);
        
        n/=2;
        
    }while (n>0);
    
    
    Matcopy(PopSize,m,dest);
    
    mkl_free(tempMat);
    mkl_free(dest);
    mkl_free(src);
    
};

void MatOp::Matcopy(int PopSize,double* A,double *B)
{
    for (int i=0;i<((PopSize+1)*(PopSize+1));i++){
    
        A[i]=B[i];
        
    }
}

void MatOp::SaveToFile(FILE *stream, double *m,int PopSize,int GridSize)
{
    
   
    
    gsl_matrix *propagator = gsl_matrix_alloc ((GridSize+1),(GridSize+1));
    
    double dx=1.0/(double)GridSize;
    
    // Continue with original program
    
	double val = 0.0, norm=0;
	
	   
    vector<int> a;
    vector<int> b;
    
    double dist=0.0;
    
    double disti=0.0;
    
    double distj=0.0;
    
    double Pn=0.0;
    double Pd=0.0;
    
    double r=0.0;
    
	for (int i=0; i<=GridSize; i++){ // NRN for (int i=0+1; i<= grid-1; i++){
		
        norm = 0.0;
		
        for (int j=0; j<=GridSize; j++){
            //cout << i << " " << j << "\n\n";
            
            if(i==0 && j==0){
                
                r=m[(i*(PopSize+1))+j];
                
            }else if (i==GridSize && j==GridSize){
                
                r=m[(PopSize*(PopSize+1))+PopSize];
                
            }else{
                
                if (i==GridSize && j==0){
                    r=m[(PopSize*(PopSize+1))+j];
                } else if (i==0 && j==GridSize){
                    r=m[(i*(PopSize+1))+PopSize];
                    
                }else{
                    if(i==0){
                        a.push_back((int)floor(i*dx*PopSize));
                        
                    }else if (i==GridSize){
                        a.push_back(PopSize);
                    }else{
                        a.push_back((int)floor(i*dx*PopSize));
                        a.push_back((int)floor(i*dx*PopSize)+1);
                    }
                    
                    if(j==0){
                        b.push_back((int)floor(j*dx*PopSize));
                    }else if (j==GridSize){
                        b.push_back(PopSize);
                    }else{
                        b.push_back((int)floor(j*dx*PopSize));
                        b.push_back((int)floor(j*dx*PopSize)+1);
                    }
                    
                    dist=0.0;
                    Pn=0.0;
                    Pd=0.0;
                    
                    // Inverse distance interpolation
                    disti=0.0;
                    distj=0.0;
                    
                    int flag=0;
                    
                    for (int k=0;k<a.size();k++){
                        for (int l=0;l<b.size();l++){
                            if(!flag){
                                
                                // cout<< i << " " << j << " " << a[k] << " " << b[l] << "\n\n";
                                
                                val=m[(a[k]*(PopSize+1))+b[l]];
                                
                                //cout << (a[k]*(PopSize+1))+b[l]<< " " << val << "\n";
                                
                                disti=((double(a[k])/double(PopSize))-(i*dx))*((double(a[k])/double(PopSize))-(i*dx));
                                distj=((double(b[l])/double(PopSize))-(j*dx))*((double(b[l])/double(PopSize))-(j*dx));
                                
                                if(disti==0 && distj==0){
                                    r=val;
                                    
                                    //cout << val <<"\n";
                                    flag=1;
                                }else{
                                    dist=sqrt(disti+distj);
                                    
                                    Pn=Pn+(pow(static_cast<double>(1.0/(dist)),static_cast<int>(2))*val);
                                    
                                    Pd=Pd+(pow(static_cast<double>(1.0/(dist)),static_cast<int>(2)));
                                    
                                    //   cout<< val << " " << i << " " << j << " " << a[k] << " " << b[l] << " " << disti << " " << distj << " " << 1/dist << " " << Pn << " " << Pd << "\n\n";
                                }
                            }
                        }
                    }
                    if(!flag) r=Pn/Pd;
                    flag=0;
                }
            }
            
 // cout << "Here";
            gsl_matrix_set(propagator,i,j,r);
            
            //            if (norm <0.0 ||norm!=norm){
            //                cout<< i<<" "<<j<<" "<<"\n";
            //                exit(0);
            //            }
            
            //norm += (j==0 || j==grid) ? 0.5*val : r;
            
            //if(r!=r) cout << i << " " << j <<  " " << r << " "<< Nlocal << "\n";
            
           // cout << r << "\n";
            
            norm += r;
            
            a.clear();
            b.clear();
		
	}
        
        
        
       
		
       row = gsl_matrix_row(propagator,i);
        
        
		norm *= dx;
		
		gsl_vector_scale(&row.vector,1.0/norm);
	}
    
    
       
    gsl_matrix_fwrite (stream, propagator); // save gsl_matrix
 gsl_matrix_free(propagator);
    
};
   


