/*
  DMS.cpp

	// Adapted from original code developed for the following paper: Andrej Fischer, Ignacio Vazquez-Garcia, Christopher J.R. Illingworth and Ville Mustonen. High-definition reconstruction of subclonal composition in cancer. Cell Reports (2014), http://dx.doi.org/10.1016/j.celrep.2014.04.055

	// Additional features for time-resolved data: neutral exact Wright-Fisher and Gaussian propagation with absorbing boundaries drift model on a frequency grid with matrix exponentiation routine and 
 
  Input: a file with four columns (extra columns of "read depth" correspond to more replicates which are analysed independently from each other):
  
	locus time reads depth

  Each locus is treated as an independent sample from the same process (linkage is not taken into account).

  For population sizes above 1000 the routines reads precomputed matrix powers from folder ./MatPow located in the same directory
  

  Output:
  -- The posterior mean and standard deviation of the emission rate.
  -- On request (--dist), the posterior distribution at each point (may generate large files and be slow to compute!).

  Required Arguments:

  --data        The input observed data file

  --mode        Emission model:
                1 Binomial
		2 Beta-Binomial

  --GWF         0 Gaussian propagation with absorbing boundaries
		1 Wright-Fisher propagation by matrix exponentiation (2N<=1000) or pre-computed matrix powers (2N>1000)
		
  Optional Arguments:
  --pre                 The prefix to put before all output files
  --grid      [int]     The grid size for the distributions (partitions [0,1] into [grid] bins).
  --dist                Prints all the posterior distributions as well.
  --sigmai    [double]  Initial drift amplitude (default minimum value sigmai=0.01)
  --sigmaf    [double]  Final drift amplitude (default minimum value sigmaf=0.05)
  --sigmagrid [int]     number of initial points for drift parameter optimization
  --shape     [double]  To fix the shape parameter in mode 2

Example:

./DMS --data [file] --mode 1 --GWF 1 -pre [str] --sigmai 0.01 --sigmaf 0.05 --sigmagrid 5 --grid 400 
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#include <list>


// GSL headers...

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"

//own headers...

//#include "emission.h"
//#include "drift.h"

#include "minimization.h"
#include "common-functions.h"



using namespace std;

struct cmdl_opts{
  const char * data_fn;
 
  const char * pre; 
  int grid, dist, mode, seed, filter_pVal,filter_shortSeg, reflect;
  double sigma, shape,xmin, xmax;
    double sigma_i,sigma_f,sigmagrid, shape_i;
    int GWF; // WF: GWF = 1 ; G: GWF = 0;
};


// *** OWN FUNCTIONS ***

void get_opts( int argc, const char ** argv, cmdl_opts& opts);
void test_opts(cmdl_opts& opts);
void default_opts(cmdl_opts& opts);
void print_opts();

double Q( const gsl_vector * x, void * p);
double find_D_parameters(Drift * myD, cmdl_opts& opts);


//

void init_parameter(gsl_vector*& var, gsl_vector*& range, vector<int>& to_opt, Drift * myD, cmdl_opts& opts);

// *** MAIN START***

int main (int argc, const char * argv[]){

    cmdl_opts opts;
  
    get_opts( argc, argv, opts);
  
    srand(opts.seed);
    
    vector<int> Loci; // vector of loci in each file
    vector<int> sTimes; // vector of sampling instances for each loci
    int nRep; // number of replicates
    int keepzero=1;
    
    //printf("%s\n",opts.data_fn);
    
       
    get_dims( opts.data_fn, nRep, Loci, sTimes, keepzero);
    
    int nLoci = (int) Loci.size();
    int total_nsTimes=0;
  
    for (int s=0; s<nLoci; s++){
        total_nsTimes += sTimes[s];
    }
  
    //announce:
    printf("\n================================================================================================\n");
    printf("\nDMS: Fitting drift model to data\n");                          
   
    printf("\nAlgorithm: Backward-forward/predict-update HMM  \n");

    printf("\nDrift model: ");
	
    if(opts.GWF){
       printf("Wright-Fisher \n");
    }else{
       printf("Gaussian \n");
    }

    printf("\nEmission model: ");
    
    if(opts.mode == 1){
        printf("Binomial \n"); // only the binomial emission model was used in the paper; other options have to be provided at the command line
    }else if (opts.mode == 2){
        printf("Beta-binomial \n");
    }else if (opts.mode == 3){
        printf("Poisson \n");
    }else if (opts.mode == 4){
        printf("Negative-binomial \n");
  
    }
    printf("\n================================================================================================\n");

    printf("\nNumber of loci taken: %i\n", nLoci);
    printf("\nNumber of replicates: %i\n", nRep);
    //printf("\nOverall number of instances: %i\n ", total_nsTimes);   
   
    //the data and emission model object
    Emission myEmit;
    myEmit.mode    = opts.mode; // 1 in paper
    myEmit.get_log = 1;
    myEmit.connect = 1; // all data will be retained
    myEmit.set( nRep, Loci, sTimes, opts.grid);
  
    get_data( opts.data_fn, &myEmit);
  
   
   
    double ** mean = new double * [nLoci];
    double ** std  = new double * [nLoci];
    double ** goftime  = new double * [nLoci];
  
    for (int s=0; s<nLoci; s++){
        mean[s] = new double [sTimes[s]];
        std[s]  = new double [sTimes[s]];
        goftime[s]  = new double [sTimes[s]]; // 
    }
  
   
    
  
    int ** mask = NULL; // the mask is used to filter out loci 
  
    if (opts.filter_pVal || opts.filter_shortSeg > 0){
        mask = new int * [nLoci];
        for (int s=0; s<nLoci; s++){
            mask[s] = new int [sTimes[s]];
            for (int l=0; l<sTimes[s]; l++){
                mask[s][l] = 1;
            }
        }
    }


                           
   
    
  // *** Filtering of each replicate ***

  for (int t=0; t<nRep; t++){
    
      printf("\nFiltering replicate %i of %i:\n", t+1, nRep);
    
      Drift myD( &myEmit, t);
      
      myD.DiffProp_opt=opts.GWF; // WF: GWF=1; G: GWF=0;
    
      //find maximum-likelihood estimates of all parameters
    
      double llh = find_D_parameters( &myD, opts);
    
      //calculate means and standard deviations...
    
      int uidx  = opts.reflect ? int(0.5*double( myD.gridSize)) :  myD.gridSize;
    
      double mx = opts.reflect ? 0.5 : myD.myEmit->xmax;
    
      double mn = myD.myEmit->xmin;
    
      gsl_vector * post = gsl_vector_alloc(uidx+1);
    
      double crit = 10.0/double(total_nsTimes);
    
      char buff[1024];
    
      sprintf(buff,"%s.posterior-%i.txt", opts.pre, t+1);
    
      FILE * total_fp = fopen(buff,"w");
      
      //cout << "HERE" << "\n";;
    
      fprintf(total_fp, "#Locus Time mean std-dev");
    
      fprintf(total_fp, " posterior %.5e %.5e\n", myD.myEmit->xmin, myD.myEmit->xmax);
      
      
    double gof=0, xobs=0, gofNorm=0;

  
     // For each loci 
    for (int s=0; s < myD.nLoci; s++){ // get posterior distribution with the ML parameters
        myD.get_posterior(s);
        double mstd = 0.0;
        
      for (int l=0; l < myD.sTimes[s]; l++){
          if (opts.reflect){//distribution in lower half
              gsl_vector_view lower = gsl_matrix_subrow( myD.gamma[s], l, 0, uidx+1);
              gsl_vector_memcpy( post, &lower.vector);
              
              double norm = gsl_blas_dasum(post);
              norm -= 0.5*(post->data[0] + post->data[uidx]);
              norm *= myEmit.dx;
	  
              if (norm <= 0.0) abort();
              gsl_vector_scale(post,1.0/norm);
          }else{
              gsl_matrix_get_row( post, myD.gamma[s], l);
          }
          
          mean[s][l] = get_mean( post, mn, mx);
          std[s][l]  = sqrt(get_var(post,mn,mx,mean[s][l]));
          mstd += std[s][l];
	
         
          
          //goodness of fit...
	
          if ((mask==NULL || mask[s][l] == 1) && myEmit.depths[t][s][l] > 0){
              
              xobs = double(myEmit.reads[t][s][l]) / double(myEmit.depths[t][s][l]);
	  
              if (opts.reflect) xobs = min(xobs,1.0-xobs);
              double x=0,g=0,dg=0;
              double b =  1.0;
              
              for (int i=0; i<=uidx; i++){
                  x = myEmit.xgrid[i] * b;
                  dg = fabs(x - xobs) * post->data[i];
                  if (i==0||i==uidx) dg *= 0.5;
                  g += dg;
              }
	  	
	      goftime[s][l]=g * myEmit.dx; // 


              gof += g * myEmit.dx;
              gofNorm += 1.0;
          }
      
      }


      
        mstd /= double(myD.sTimes[s]);
      
        //filter out data points which are not compatible with the emission model
      
        if (opts.filter_pVal){
            for (int l=0; l < myD.sTimes[s]; l++){
                double pval = myEmit.get_pval( t, s, l, mean[s][l]);
                if (pval < crit) mask[s][l] = 0;
            }
        }
      
        if (opts.filter_shortSeg > 0){
            int last=0;
            for (int l=1; l < myD.sTimes[s]; l++){
                if (fabs(mean[s][l] - mean[s][l-1]) > 4.0*mstd){
                    if (l < last + opts.filter_shortSeg){
                        for (int i=last; i<l; i++) mask[s][i] = 0;
                    }
                    last = l;
                }
            }
        }
      
        // print posterior information to file
      
        for (int l=0; l < myD.sTimes[s]; l++){
            fprintf(total_fp, "%i %6i %.2e %.2e",Loci[s], myD.times[s][l], mean[s][l], std[s][l]);
            if(opts.dist==1){// full posterior distribution? LARGE!
                for (int i=0; i <= myD.gridSize; i++){
                    fprintf(total_fp, " %.2e", gsl_matrix_get( myD.gamma[s], l, i));
                }
            }
            fprintf(total_fp,"\n");
        }
        gsl_matrix_free(myD.gamma[s]);
        myD.gamma[s] = NULL;
    }
   
    
      fclose(total_fp);
      
      gsl_vector_free(post);
    
      // *** PROCLAIM RESULTS ***

        // Export gof results at each time point and loci

	char goft[1024];
 	sprintf(goft,"%s.goftime-%i.txt", opts.pre, t+1);
	
	 FILE * goftime_estimates = fopen(goft,"w");

	for (int s=0; s < myD.nLoci; s++){//get posterior distribution with the ML parameters
     
      		for (int l=0; l < myD.sTimes[s]; l++){

		fprintf(goftime_estimates, "%i %6i %.2e",Loci[s], myD.times[s][l],goftime[s][l]/gofNorm);
		fprintf(goftime_estimates,"\n");

		}
	}

	fclose(goftime_estimates);

   
 
      // export ML parameters
      
      char buffparams[1024];
      
      sprintf(buffparams,"%s.paramsML-%i.txt", opts.pre, t+1);
      
      FILE * params_estimates = fopen(buffparams,"w");
      
      if(1/(myD.sigma*myD.sigma)<1000){ // above 2N=1000 all pre computed matrices were generated every 100

	if(opts.GWF){

          fprintf(params_estimates,"Filtered sample %i of %i replicates: llh = %.5e, gof = %.5e, --N %.3e",t+1, nRep, llh,  gof/gofNorm, 0.5 /(myD.sigma*myD.sigma));
         
          printf("Best model llh = %.5e, gof = %.5e, --N %.3e", llh,  gof/gofNorm, 0.5 /(myD.sigma*myD.sigma));

	}else{

	  fprintf(params_estimates,"Filtered sample %i of %i replicates: llh = %.5e, gof = %.5e, --Sigma %.3e",t+1, nRep, llh,  gof/gofNorm, myD.sigma);
         
          printf("Best model llh = %.5e, gof = %.5e, --Sigma %.3e", llh,  gof/gofNorm, myD.sigma);

	}
      }else{

         if(opts.GWF){    
   
	  fprintf(params_estimates,"Filtered sample %i of %i replicates: llh = %.5e, gof = %.5e, --N %.3e",t+1, nRep, llh,  gof/gofNorm, floor(((0.5/(myD.sigma*myD.sigma))+50)/100)*100);

          printf("Best model llh = %.5e, gof = %.5e, --N %.3e", llh,  gof/gofNorm, floor(((0.5/(myD.sigma*myD.sigma))+50)/100)*100);
         
          }else{

 	  fprintf(params_estimates,"Filtered sample %i of %i replicates: llh = %.5e, gof = %.5e, --Sigma %.3e",t+1, nRep, llh,  gof/gofNorm, myD.sigma);
         
          printf("Best model llh = %.5e, gof = %.5e, --Sigma %.3e", llh,  gof/gofNorm, myD.sigma);
         
         }
      }



      fclose(params_estimates);
      
      cout<<endl;
    
         
      // *** RESET ***
      myEmit.delete_old_Emit();
      myEmit.range_set = 0;
      myEmit.reset_mask();
  
  }
  
    //print filtered
  
    if (opts.filter_pVal || opts.filter_shortSeg > 0){
        char buff[1024];
        sprintf(buff,"%s.filtered.txt", opts.pre);
        FILE * filtered_fp = fopen(buff,"w");
        
        for (int s=0; s < myEmit.nLoci; s++){
            for (int l=0; l < myEmit.sTimes[s]; l++){
                if( mask[s][l] == 1 ){
                    fprintf( filtered_fp, "%i %6i", Loci[s], myEmit.times[s][l]);
                    for (int t=0; t<nRep; t++)
                        fprintf( filtered_fp, " %3i %3i", myEmit.reads[t][s][l], myEmit.depths[t][s][l]);
                    fprintf( filtered_fp, "\n");
                }
            }
            delete [] mask[s];
        }
        delete [] mask;
        fclose(filtered_fp);
    }
    
    for (int s=0; s < myEmit.nLoci; s++){
        delete [] mean[s];
        delete [] std[s];
	delete [] goftime[s];
    }
  
    delete [] mean;
    delete [] std;
    delete [] goftime;

    //done
    return (0);

}

// *** MAIN END ***


void default_opts(cmdl_opts& opts){
    opts.data_fn         = NULL ;
    opts.pre             = NULL ;
    opts.grid            = 400  ;
    opts.sigmagrid       = 5    ; // Test starting seed values of drift at 5 points from sigma_i to sigma_f
    opts.dist            = 0    ;
    opts.sigma           = -1.0 ; // not to be changed  
    opts.shape           = -1.0 ;
    opts.sigma_i         = 0.01 ; // seed value for starting drift amplitude tested in the likelihood optimization function
    opts.sigma_f         = 0.05 ; // seed value for starting drift amplitude tested in the likelihood optimization function
    opts.sigmagrid        = 5    ; // number of starting points for drift parameter optimization
    opts.shape_i         = -1.0 ;
    opts.mode            = 1    ; // default option is the Binomial emission model
    opts.xmin            = -1.0 ;
    opts.xmax            = -1.0 ;
    opts.seed            = (int) time(NULL);
    opts.filter_pVal     = 0    ;
    opts.filter_shortSeg = 0    ;
    opts.reflect         = 0    ;
    opts.GWF             = 1    ; // drift model option: 0-Gaussian 1-Wright-Fisher

}

// get command line arguments...

void get_opts( int argc, const char ** argv, cmdl_opts& opts){
  default_opts(opts);
  int opt_idx = 1;
  string opt_switch;  
  while ( opt_idx < argc && (argv[opt_idx][0] == '-')){
    opt_switch = argv[opt_idx];
    if ( opt_switch.compare("--print-options") == 0){
      print_opts();
      exit(0);
    }
    opt_idx++;
    if (opt_idx==argc) break;
    if ( argv[opt_idx][0] == '-') continue;
    if ( opt_switch.compare("--data") == 0){ // the input data
      opts.data_fn = argv[opt_idx];
    }
    else if ( opt_switch.compare("--pre") == 0){  // the prefix for all output files
      opts.pre = argv[opt_idx];
    }
    else if ( opt_switch.compare("--grid") == 0){ // the size of the grid for the continuous distributions
      opts.grid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigmagrid") == 0){ // the size of the grid for the continuous distributions
      opts.sigmagrid = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--seed") == 0){ // random seed
      opts.seed = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--mode") == 0){//emission model: see above
      opts.mode = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigma") == 0){//diffusion constant
      opts.sigma = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--shape") == 0){//shape parameter for mode 2/4
      opts.shape = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigmai") == 0){ // drift parameter minimum grid value 
      opts.sigma_i = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--sigmaf") == 0){ // drift parameter maximum grid value 
      opts.sigma_f = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--shapei") == 0){ // shape parameter for mode 2/4
      opts.shape_i = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--xmin") == 0){ // shape parameter for mode 2/4
      opts.xmin = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--xmax") == 0){ // shape parameter for mode 2/4
      opts.xmax = atof(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--dist") == 0){ // whether to print posterior
      opts.dist = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--filter-pVal") == 0){ // whether to filter out some data points
      opts.filter_pVal = 1;
    }
    else if ( opt_switch.compare("--filter-shortSeg") == 0){ // whether to filter out some data points
      opts.filter_shortSeg = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--reflect") == 0){ // whether to filter out some data points
      opts.reflect = atoi(argv[opt_idx]);
    }
    else if ( opt_switch.compare("--GWF") == 0){ //  use Gaussian propagation (0) or WF propagation (1)
        opts.GWF = atoi(argv[opt_idx]);
    }else {
      cout << "Usage: filterHD --print-options"<<endl;
      exit(1);
    }
    opt_switch.clear();
    opt_idx++;
  }
  test_opts(opts);
}


void test_opts(cmdl_opts& opts){
  if (opts.mode==0){
    cout<<"ERROR: choose emission mode with --mode [1,2,3,4]\n";
    exit(1);
  }
  
  if(opts.reflect==1 && (opts.mode==3||opts.mode==4)){
    cout<<"ERROR: --reflect [0/1] can only be used in mode 1 and 2.\n";
    exit(1);
  }
 
}

void print_opts(){
  cout<<"filterHD --data [file] --mode [1,2,3,4] --pre [string:./out] --grid [int:100] --sigmagrid [int] --shape [double] --dist [0/1:0] --sigmai [double] --sigmaf [double]  --shapei [double] --xmin [double] --xmax [double]"<<endl;
  exit(0);
}



// *** LEARN PARAMETERS OF THE JUMP-DIFFUSION MODEL ***
double find_D_parameters(Drift * myD, cmdl_opts& opts){
  vector<int> to_opt;
  int nvar=0;
  double llh=0;
  if (opts.xmin >= 0.0 && opts.xmax >= 0.0){
    myD->myEmit->xmin = opts.xmin;
    myD->myEmit->xmax = opts.xmax;
    myD->myEmit->ymin = opts.xmin;
    myD->myEmit->ymax = opts.xmax;
    myD->myEmit->set_grid();
  }
  else{
    myD->myEmit->init_range(myD->Rep);//get the initial range of rates
  }
 
   if (opts.sigma < 0.0){ // diffusion constant
    to_opt.push_back(1);
    nvar++;
  }
  else{
    myD->sigma = opts.sigma;
   
      if(myD->DiffProp_opt==0){
          myD->get_DiffProp();
      }else{
          myD->get_DiffPropWF();
      }
      
      
      
  }
  if ( opts.shape < 0.0 && (opts.mode == 2 || opts.mode == 4)){ // shape parameter
    to_opt.push_back(3);
    nvar++;
  }
  else{
    myD->myEmit->shape = opts.shape;
  }

  
    if(nvar>0){
    
        gsl_vector * var   = gsl_vector_calloc(nvar);
        gsl_vector * range = gsl_vector_calloc(nvar);
    
        //set initial values

        init_parameter(var, range,  to_opt, myD, opts);
        fpar myfpar;

	myfpar.myD      = myD;
        myfpar.to_opt   = to_opt;
        myfpar.sigm_i   = opts.sigma_i; // 
        myfpar.sigm_f   = opts.sigma_f; // 
        myfpar.sigmgrid = opts.sigmagrid; //
    
        void * param = static_cast<void*>(&myfpar);
    
        // get the ML estimates and ML value

        int steps = 0;
        gsl_vector ** simplex = NULL;
        gsl_vector * lower    = NULL;
      
	for (int i=0; i<nvar; i++){
           if(to_opt[i] == 1){//diffusion constant
           //     printf("%-11s ", "sigma");
      
            }
            else if(to_opt[i] == 3){//shape parameter
             //   printf("%-11s ", "shape");
            }
        }
            
         //printf("%-11s \n", "llh");
           
        llh = - find_local_optimum( 0, simplex, lower, var, range,  param, &Q, 1.0e-3, steps, 1);
        
       

        
        //adapt the range if needed
        if ((opts.mode==3 || opts.mode==4) && (opts.xmin < 0.0 && opts.xmax < 0.0)){
            int redo = myD->adapt_range();
            if (redo==1){
                printf("Adapted range to %.3e < x < %.3e\n", myD->myEmit->xmin, myD->myEmit->xmax);

                //init_parameter( var, range,  to_opt, myD, opts);

                llh = - find_local_optimum( 0, simplex, lower, var, range, param, &Q, 1.0e-3, steps, 1);
            }
        }
    
        //set the ML values into the objects
        for (int i=0; i<nvar; i++){
            if(to_opt[i] == 1){//diffusion constant
               
                myD->sigma = var->data[i];
             
                if(myD->DiffProp_opt==0){
                    myD->get_DiffProp();
                }else{
                    myD->get_DiffPropWF();
                }
                
          
            }
            else if(to_opt[i] == 3){//shape parameter
                myD->myEmit->shape = var->data[i];
                myD->myEmit->set_EmitProb(myD->Rep);
            }
        }
        gsl_vector_free(var);
    
        gsl_vector_free(range);
  
    }
    else{
        llh = myD->get_total_llh();
        //adapt the range if needed
        if ( (opts.mode==3 || opts.mode==4) && (opts.xmin < 0.0 && opts.xmax < 0.0) ){
            myD->adapt_range();
            printf("Adapted range to %.3e < x < %.3e\n", myD->myEmit->xmin, myD->myEmit->xmax);
            llh = myD->get_total_llh();
        }
    }
    
    return(llh);
}

void init_parameter(gsl_vector*& var, gsl_vector*& range, vector<int>& to_opt, Drift * myD, cmdl_opts& opts){
  int nvar = (int) var->size;
  for (int i=0; i<nvar; i++){
    if(to_opt[i] == 1){//diffusion constant
      if (opts.sigma_i > 0.0){
	var->data[i]   = opts.sigma_i;
      }
      else{
	var->data[i]   = 0.1*myD->myEmit->dx / sqrt(myD->myEmit->median_dist);
      }
	range->data[i] = 0.0;
    }
    else if(to_opt[i] == 3){//shape parameter
      var->data[i]   = (opts.shape_i > 0.0) ? opts.shape_i : 100.0;
      range->data[i] = 0.0;
    }
  }
}

double Q( const gsl_vector * x, void * p){
  //Drift * myD = static_cast<Drift*> (p);
  fpar * myfpar = static_cast<fpar*> (p);
  int nvar = (int) (myfpar->to_opt).size();
  gsl_vector * var   = gsl_vector_alloc(nvar);
  gsl_vector * range = gsl_vector_alloc(nvar);
  for (int i=0; i<nvar; i++){
   if((myfpar->to_opt)[i] == 1){//diffusion constant in [0,\infty]
      range->data[i] = 0.0;
    }
   else if((myfpar->to_opt)[i] == 3){//shape parameter in [0,\infty]
      range->data[i] = 0.0;
    }
  }
  gsl_vector ** simplex = NULL;
  gsl_vector * lower    = NULL;
  int err = arg_unmap( x, 0, simplex, lower, var, range);
  if (err==1){
    gsl_vector_free(var);
    gsl_vector_free(range);
    return(1.0e20);
  }
    
   double curr_sigma=0.0;
    
  //set the ML values into the objects
  for (int i=0; i<nvar; i++){
   if((myfpar->to_opt)[i] == 1){//diffusion constant
      myfpar->myD->sigma = var->data[i];
    
        curr_sigma=myfpar->myD->sigma;

	

        if(((floor(((1.0/(curr_sigma*curr_sigma))+50)/100)*100)<=10000)){ // necessary in order to avoid pop sizes above values that we didnot compute the matrices for
        
            if(myfpar->myD->DiffProp_opt==0){
        
//                myfpar->myD->get_DiffProp();
        
            }else{
        
                myfpar->myD->get_DiffPropWF();
        
          }
        }
        
            
        if(myfpar->myD->DiffProp_opt==0){
            
         myfpar->myD->get_DiffProp();
            
       }
    }else if((myfpar->to_opt)[i] == 3){//shape parameter
      myfpar->myD->myEmit->shape = var->data[i];
      myfpar->myD->myEmit->set_EmitProb(myfpar->myD->Rep);
    }
  }
  // Do Fwd to get LLH
      
    curr_sigma=myfpar->myD->sigma;
    
    double llh=0.0;
     
    if(((floor(((1.0/(curr_sigma*curr_sigma))+50)/100)*100)>10000)){ // necessary in order to avoid population sizes above values that do not have computed the matrix powers
	if(myfpar->myD->DiffProp_opt==0){
        llh = (myfpar->myD)->get_total_llh(); 
        }else{
	llh=-1e25;
	}    
   }else{
    
        llh = (myfpar->myD)->get_total_llh();
    }
    
  gsl_vector_free(var);
  gsl_vector_free(range);


  return(-llh);
}

