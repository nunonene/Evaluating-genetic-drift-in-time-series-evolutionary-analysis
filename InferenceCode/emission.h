//emission.h

// Adapted from original code developed for the following paper: Andrej Fischer, Ignacio Vazquez-Garcia, Christopher J.R. Illingworth and Ville Mustonen. High-definition reconstruction of subclonal composition in cancer. Cell Reports (2014), http://dx.doi.org/10.1016/j.celrep.2014.04.055

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <set>
//#include <unordered_map>
#include <vector>
#include <list>
#include <algorithm>

// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_psi.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_cdf.h"



using namespace std;

class Emission{
public:
  Emission();
  void set(int nRep, vector<int>& Loci, vector<int>& stimes, int grid);
  ~Emission();
  void clear();
  void delete_old_Emit();
  int is_set;
  void set_dist();
  int dist_set;
  int connect;
  double median_dist;
  map<unsigned int, int> dist_count;
  map<unsigned int, int> frequent_dist;
  int get_log, get_der, get_mv;
  unsigned int nmax, Nmax;
  
  map< unsigned int, map< unsigned int, gsl_vector*> > EmitProb;
  map< unsigned int, map< unsigned int, gsl_vector*> > EmitLog;
  double shape, log_shape;
  double minRate, maxRate;
  //
  int mode, reflect, log_space;
  void set_EmitProb(int rp);
  void binomial(int N, int n);
  void beta_binomial(int N, int n);
  void poisson(int N, int n);
  void negative_binomial(int N, int n);
  
 
  int EmitProb_set;
  //
  int nRep, nLoci;
  int gridSize;
  double dx,xmin,xmax;
  double dy,ymin,ymax;
  double * xgrid;
  double * ygrid;
   unsigned int *** reads;
  unsigned int *** depths;
  unsigned int ** times;
  unsigned int ** mask;
  unsigned int ** dist;
  unsigned int *** nObs;
  void get_nObs();
  int *sTimes;
  int * Locus;
  std::set<int> Loci;
  int * idx_of;
  int maxlocus;
  double ** bias;
  double ** log_bias;
 
 
 
  int total_stimes, total_events;
  unsigned int total_dist;
  void set_grid();
  void init_range(int rp);
  int range_set;
  void reset_mask();
  double get_pval(int rp, int s, int site, double mean);
 
  double *** mean_tcn; // mean total copy number
  double **** av_cn; // copy number availability
  
  int * nEvents;
  unsigned int ** Event_of_idx; // map from idx to cnv-event 
 
  
  unsigned int ** idx_of_event; // map from event to idx
  unsigned int ** event_of_idx; // map from idx to event
  int idx_to_Event_mapped;
 
};


bool value_comparer(std::map<int,double>::value_type &i1, std::map<int,double>::value_type &i2);
