//
//  main.cpp
//  MatrixPow
//
//  Created by Nuno Rocha Nene on 22/07/2014.
//  Copyright (c) 2014 ___nn276___. All rights reserved.
//


#include "TMatBinomial.h"
#include "MatOp.h"


int main(int argc, const char **argv)
{
    TMatBinomial TM; // class for binomial distribution
    
    TM.GetParameters(argc,argv); // parse inputs
    
    
    string MatMultOption="OpenBLAS"; // for the sake of file name consistency with other files generated elsewhere
    
    
    double *I;
    double *m;
    
    I=(double*) mkl_malloc(((TM.PopSize)+1)*((TM.PopSize)+1)*sizeof(double),64);
    
    m=(double*) mkl_malloc(((TM.PopSize)+1)*((TM.PopSize)+1)*sizeof(double),64);
    
    // Set I
    // =====
    for (int i=0;i<(TM.PopSize+1);i++){
        for (int j=0;j<(TM.PopSize+1);j++){
        
            if(i==j){
            
                I[(((TM.PopSize)+1)*i)+j]=1;
            
            }else{
                
                I[(((TM.PopSize)+1)*i)+j]=0;
            
            }
        
        }
        
    }
    
     
    
    MatOp MP; // instantiate class with exponentiation method
    
    
    int K=2; // necessary for binomial function
    
    double temp;
    
    vector <double> ps;
    
    vector <int> nn;
    
    
    // Calculate original transition matrix
    // ====================================
    
    for (int i=0;i< ((TM.PopSize)+1);i++){
        
        double sum=0;
        
        for (int j=0;j<((TM.PopSize)+1);j++){
            
            
            
            ps.push_back((i/(double)(TM.PopSize))+((TM.Sigma)*((i/(double)(TM.PopSize)) * (1-(i/(double)(TM.PopSize)))))); // is equal to ps.push_back(i/(double) PopSize) if Sigma=0;

            nn.push_back(j);
            
            ps.push_back(1-ps[0]);
                    
            nn.push_back((TM.PopSize)-nn[0]);
            
            
            temp=TM.Binomial_pdf (K,ps,nn);
            
            
            if(isnan(temp)){
                temp=1;
            }
            
            
            m[(((TM.PopSize)+1)*i)+j]=temp;
            
            
            sum+=temp;
         
            ps.clear();
            nn.clear();
            
        }
        
        ps.clear();
        nn.clear();
   
        // Row-wise normalization
        // ======================
        
        for (int k=0;k<((TM.PopSize)+1);k++){
            
            m[(((TM.PopSize)+1)*i)+k]= m[(((TM.PopSize)+1)*i)+k]/sum;
            
        }
        
    }
    
    //cout << "TMatBinomial is ready\n\n";
    
        
    // Exponentiate transition matrix (squares method)
    // ==============================================
    
    cout << "Start matrix exponentiation\n\n";
    
  
    
    
    
    MP.PowMatOpenBLAS((TM.NumGen),(TM.PopSize),I,m); // use OpenBLAS
    
    
     mkl_free(I);
    
    cout << "Matrix exponentiation completed.\n\n";
    
    // Save gsl_matrix to file
    // ===================
    
    FILE  *out_file; // using pointer to FILE as demanded by gsl_matrix_fwrite (FILE *stream, const gsl_matrix *m) or gsl_matrix_fwrite(FILE *stream, const gsl_matrix *m, const char* type)
    
    ostringstream convert;   // stream used for the conversion
    
    convert << (TM.PopSize);      // insert the textual representation of the number in the characters in the stream
    
    string PopSizeString = convert.str(); // set to the contents of the stream
    
    convert.str(""); // required to reset the string to be empty;
    convert.clear(); // required to clear any error flags that may be set
    
     convert << (TM.GridSize);
    
    string GridSizeString = convert.str(); // set to the contents of the stream
    
    convert.str(""); // required to reset the string to be empty;
    convert.clear(); // required to clear any error flags that may be set
    
    convert<<TM.NumGen;
    
    string PowString= convert.str(); //
    
    convert.str("");
    convert.clear();
    
    convert<<TM.Sigma;
    
    string SigmaString= convert.str(); //
    
    convert.str("");
    convert.clear();

    
    std::string filename = std::string("./")+ std::string("MatPow_Pow")+PowString+std::string("_PopSize")+PopSizeString+std::string("_GridSize")+GridSizeString+std::string("_Sigma")+SigmaString+std::string("_")+MatMultOption+std::string(".dat"); // filename
    
    
    out_file=fopen(filename.c_str(),"wb");
    
    
    
    if (out_file!=NULL) // Check if opened and save
    
    {
         MP.SaveToFile(out_file, m,TM.PopSize,TM.GridSize);
        
    }else{
     
        cout << "Unable to open file\n";
        perror(filename.c_str());
        exit(0);
    }
   
    
    fclose(out_file);
    
    
    // Free matrices
    // =================
    
    mkl_free (m);
    

    
    return 0;
    
    
};
