//drift.h


// Adapted from http://dx.doi.org/10.1016/j.celrep.2014.04.055


#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multimin.h"


using namespace std;

#include "hypgeo.hpp"



class Drift{
public:
  Drift(Emission * emit, int rp);
  ~Drift();
  Emission * myEmit; 
  int nLoci;
  int Rep;
  int mode;
  double sigma, N;
  int Fwd_done, Bwd_done, wTotal, save_alpha;
  int gridSize;
    
  int gridSize_WF;
    
  int * sTimes;
  unsigned int ** dist;
  unsigned int ** times;
  unsigned int ** mask;  
  
  double ** bias;
  void get_EmitProb(int read, int depth, double * xgrid, gsl_vector * eprob);// emission probability
  gsl_vector * proposal;
  gsl_matrix ** alpha;
  gsl_matrix ** gamma;
  gsl_matrix ** total;
  void set_pstay();
  int pstay_set;
  double do_Fwd(int s);
  void do_Bwd(int s);
  //
  gsl_matrix ** DiffProp;
  int set_DiffProp(gsl_matrix * propagator, double variance);
  void get_DiffProp();
	
  int DiffProp_set;
    int DiffProp_opt;
  void reset_DiffProp();
	
	int set_DiffPropWF(gsl_matrix * propagator, double N, double t);
	void get_DiffPropWF();
	
	//int DiffProp_set;
	//void reset_DiffPropWF();
	
  
	
vector<int> is_identity;
  map<unsigned int,int> position;
  //
  //void set_DiffProp_Log(gsl_matrix * propagator, double variance);
  int predict(gsl_vector * prior, gsl_vector * post, gsl_matrix*& DiffProp, gsl_matrix**& DP_pt, int sampe, int site);
  double total_llh;
  double get_total_llh();
  void get_posterior(int s);
  int adapt_range();
};

