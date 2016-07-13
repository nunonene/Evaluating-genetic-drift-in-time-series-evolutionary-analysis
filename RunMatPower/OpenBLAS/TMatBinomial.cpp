//
//  TMatBinomial.cpp
//  MatrixPow
//
//  Created by Nuno Rocha Nene on 21/07/2014.
//  Copyright (c) 2014 ___nn276___. All rights reserved.
//

#include "TMatBinomial.h"




double TMatBinomial::Binomial_pdf (const size_t K, vector<double> &p,vector<int>& n)
{
    return exp (Binomial_lnpdf (K, p, n));
};


double TMatBinomial::Binomial_lnpdf (const size_t K, vector<double> &p,vector<int>& n)
{
    size_t k;
    double N = 0;
    double log_pdf = 0.0;
    double norm = 0.0;
    
    for (k = 0; k < K; k++)
    {
        N += n[k];
    }
    
    for (k = 0; k < K; k++)
    {
        norm += p[k];
    }
    
    log_pdf = gsl_sf_lnfact (N);
    
    for (k = 0; k < K; k++)
    {
        log_pdf -= gsl_sf_lnfact (n[k]);
    }
    
    for (k = 0; k < K; k++)
    {
        log_pdf += log (p[k] / norm) * n[k];
    }
    
    return log_pdf;
};

void TMatBinomial::GetParameters (int argc, const char **argv) {
	
    string p_switch;
    
    int x=1;
    
   
    
//    if(**argv=='/'){
//        
//        cout << "Incorrect usage!!!\n\n";
//        cout << "Please provide --PopSize, --NumGen, --Sigma and --NumThreads\n\n ";
//        exit(1);
//        
//    }
    
    int NumParams=0;
	
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		cout << p_switch << "\n";
		if (p_switch.compare("--PopSize")==0) {
			x++;
			PopSize=atoi(argv[x]);
            
            NumParams+=1;
		}
        else if (p_switch.compare("--GridSize")==0) {
            x++;
            GridSize=atoi(argv[x]);
            NumParams+=1;
        }
		else if (p_switch.compare("--NumGen")==0) {
			x++;
			NumGen=atoi(argv[x]);
             NumParams+=1;
		}
        else if (p_switch.compare("--Sigma")==0) {
			x++;
			Sigma=atof(argv[x]);
             NumParams+=1;
		}
        else if (p_switch.compare("--NumThreads")==0) {
			x++;
			NumThreads=atoi(argv[x]);
            NumParams+=1;
		}
		else{
            
            cout << "Incorrect usage\n ";
            cout << "Please provide --PopSize,--GridSize, --NumGen, --Sigma and --NumThreads\n\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
    
    if( NumParams!=5){
        cout << "Incorrect usage!!!\n\n";
        cout << "Please provide the 5 required parameters :--PopSize,--GridSize, --NumGen, --Sigma and --NumThreads\n\n ";
        exit(1);
    
    }
    
    
}
