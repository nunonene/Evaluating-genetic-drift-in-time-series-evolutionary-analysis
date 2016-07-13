/*
  MatOp.cpp

  Allows for exponentiation of matrices by the squaring method
 
  Nuno Rocha Nene   

*/

#include "MatOp.h"


void MatOp::PowMatBLAS(int expo, int PopSize, gsl_matrix *I, gsl_matrix *m){ // with BLAS matrix multiplication routines
    
    int n=expo;
    
    
    gsl_matrix *dest = gsl_matrix_alloc (PopSize+1,PopSize+1);
    
    gsl_matrix *src = gsl_matrix_alloc (PopSize+1,PopSize+1);
    
    gsl_matrix *tempMat = gsl_matrix_alloc (PopSize+1,PopSize+1);
    
    gsl_matrix_memcpy (dest, I);
    
    gsl_matrix_memcpy (src, m);
    
    
    //std::cout<< "Calculating Matrix Powers\n";


     // Perform matrix exponentiation
    // ==============================
    do
    {
       // std::cout<< n << "\n";
        
        if(n % 2 != 0){
            
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                            1.0, dest, src,
                            0.0, tempMat); // gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
            
            gsl_matrix_memcpy(dest,tempMat);
            
        }
        
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, src, src,0.0, tempMat); // gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
        
        
        gsl_matrix_memcpy(src,tempMat);
        
        n/=2;
        
    }while (n>0);
    
    
    gsl_matrix_memcpy(m,dest);
    
    gsl_matrix_free(tempMat);
    gsl_matrix_free(dest);
    gsl_matrix_free(src);
    
};

void MatOp::ReadInData(std::string dataFile,int N, gsl_matrix* m){
    
  
    
    std::ostringstream convert;   // stream used for the conversion
    
    std::fstream file;
    
    std::string line;
    
    double inputValue;
    
    std::string inputValueAllele;
    
    std::vector<double> d;
    
    
    // Open file containing data generated
    // ===================================
    file.open(dataFile.c_str(), std::fstream::in);
    
    if(file.fail())
    {
        std::cout<<"Error opening data file:"<< dataFile <<std::endl;
        exit(1);
    }
    
    // Count the number of lines of data
    // =================================
    
    int col=0;
    int row=0;
    
    while(getline(file, line)){
        std::istringstream iss(line); // treat the string as an I/O stream

        while (iss){
            iss>>inputValue;
            gsl_matrix_set (m, row, col,inputValue);
           // std::cout<<  gsl_matrix_get(m, row, col) <<"\n";
        };
        iss.clear();
        line.clear();
        
        col++;
        
        if(col==N+1){
            col=0;
            row++;
        }
       
    };
    
    file.close();
}
