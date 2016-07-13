//emission.cpp

// Adapted from original code developed for the following paper: Andrej Fischer, Ignacio Vazquez-Garcia, Christopher J.R. Illingworth and Ville Mustonen. High-definition reconstruction of subclonal composition in cancer. Cell Reports (2014), http://dx.doi.org/10.1016/j.celrep.2014.04.055

//own headers...
#include "emission.h"

#define PI 3.1415926

using namespace std;

// Constructor
Emission::Emission(){
  is_set = 0;
  dx=0;
  EmitProb_set = 0;
  dist_set = 0;
  range_set = 0;
  shape    = 1.0;
  reflect  = 0;
  mode     = 1; // 1: binomial, 2: beta-binomial, 3: Poisson, 4: negative binomial
  xmin     = 0;
  xmax     = 0;
  ymin     = 0;
  ymax     = 0;
  xgrid = NULL;
  ygrid = NULL;
  maxRate  = 0; 
  get_log  = 0; //get log of emission probability?
  get_der  = 0; //get derivatives?
  get_mv   = 0; //get mean and variance?
  bias     = NULL;
  log_bias = NULL;
  mask     = NULL;

  idx_of_event   = NULL;
  event_of_idx   = NULL;
  nEvents        = NULL;
  Event_of_idx   = NULL;
  log_space      = 0;
  
  idx_to_Event_mapped=0;
  connect=0;
  nObs = NULL;
}


// real constructor
void Emission::set(int nrp, vector<int>& Lci, vector<int>& nts, int grid){
  if ((int) Lci.size() != (int) nts.size()){
    cout<<"FATAL ERROR (use gdb to locate).\n";
    abort();
  }
  nLoci  = (int) nts.size();
  nRep    = nrp;
  sTimes  = new int [nLoci]; 
  Locus     = new int [nLoci];
  Loci.clear();
  total_stimes   = 0;
  maxlocus = 0;
  for (int s=0; s<nLoci; s++){
    Locus[s]     = Lci[s];
    Loci.insert(Locus[s]);
    maxlocus      = max( maxlocus, Locus[s]);
    sTimes[s]   = nts[s];
    total_stimes += sTimes[s];
  }
  idx_of = new int [maxlocus+1];
  for (int i=0; i<=maxlocus; i++)  idx_of[i] = -1;
  for (int s=0; s<nLoci; s++) idx_of[Locus[s]] = s;
  // set xgrid...
  gridSize = grid;
  xgrid  = new double [gridSize+1];
  ygrid  = new double [gridSize+1];
  reads  = new unsigned int ** [nRep];
  depths = new unsigned int ** [nRep];
  for (int t=0; t<nRep; t++){
    reads[t]  = new unsigned int * [nLoci];
    depths[t] = new unsigned int * [nLoci];
    for (int s=0; s<nLoci; s++){
      reads[t][s]    = new unsigned int [sTimes[s]];
      depths[t][s]   = new unsigned int [sTimes[s]];
    }
  }
  times = new unsigned int * [nLoci];
  mask = new unsigned int * [nLoci];
 
  for (int s=0; s<nLoci; s++){
    times[s] = new unsigned int [sTimes[s]];
    mask[s] = new unsigned int [sTimes[s]];
   
    
  }
  Emission::reset_mask();
 
  is_set=1;
}
/*
void Emission::get_nObs(){
  if (nObs!=NULL) abort();
  nObs = new unsigned int ** [nRep];
  for (int t=0; t<nRep; t++){
    nObs[t] = new unsigned int * [nLoci];
    for (int s=0; s<nLoci; s++){
      if (nEvents[s] == 0){
	nObs[t][s] = NULL; 
	continue;
      }
      nObs[t][s] = new unsigned int [nEvents[s]];
    }
  }
  for (int s=0; s<nLoci; s++){
    for (int evt=0; evt<nEvents[s]; evt++){
      int first = idx_of_event[s][evt];     
      int last 
	= (evt < nEvents[s]-1) 
	? idx_of_event[s][evt+1] - 1 
	: sTimes[s] - 1; 
      unsigned int n,N;
      for (int t=0; t<nRep; t++){
	nObs[t][s][evt] = 0;
	for (int idx=first; idx<=last; idx++){
	  n = reads[t][s][idx];
	  N = depths[t][s][idx];
	  if (N==0){
	    if (n>0) abort();
	    continue;
	  }
	  nObs[t][s][evt]++;
	}
      }
    }
  }
}
*/

//Destructor
Emission::~Emission(){
  
  if (mean_tcn != NULL){
    for (int t=0; t<nRep; t++){
      for (int s=0; s<nLoci; s++){
	if (mean_tcn[t][s] != NULL) delete [] mean_tcn[t][s];
      }
      delete [] mean_tcn[t];
    }
    delete [] mean_tcn;
  }
  if (av_cn != NULL){
    for (int t=0; t<nRep; t++){
      for (int s=0; s<nLoci; s++){
	if (av_cn[t][s] != NULL){
	  for (int evt=0; evt<nEvents[s]; evt++){
	    delete [] av_cn[t][s][evt];
	  }
	  delete [] av_cn[t][s];
	}
      }
      delete [] av_cn[t];
    }
    delete [] av_cn;
  }
  if (bias != NULL){
    for (int s=0; s<nLoci; s++) delete [] bias[s];
    delete [] bias;
  }
  if (log_bias != NULL){
    for (int s=0; s<nLoci; s++) delete [] log_bias[s];
    delete [] log_bias;
  }
  if (is_set==1){
    Emission::clear();
  }
}

void Emission::clear(){
  if (is_set==0) abort();
  for (int t=0; t<nRep; t++){
    for (int s=0; s<nLoci; s++){
      delete [] reads[t][s];
      delete [] depths[t][s];
    }
    delete reads[t];
    delete depths[t];
  }



  for (int s=0; s<nLoci; s++){
    delete [] times[s];
    delete [] mask[s];
    delete [] dist[s];
  }



  delete [] reads;
  delete [] depths;
  delete [] times;
  delete [] mask;
  delete [] dist;
  delete [] xgrid;
  delete [] ygrid;
  delete [] sTimes;
  delete [] nEvents;




 // for (int s=0; s<nLoci; s++){
//    delete [] event_of_idx[s];
//    if (idx_of_event[s] != NULL)  delete [] idx_of_event[s];
//  }



 // delete [] event_of_idx;
 // delete [] idx_of_event;


  if (EmitProb_set==1) Emission::delete_old_Emit();


  is_set   = 0;
  dist_set = 0;
}





void Emission::delete_old_Emit(){
  //unordered_map< unsigned int, unordered_map< unsigned int, gsl_vector*> >::iterator it1;
  map< unsigned int, map< unsigned int, gsl_vector*> >::iterator it1;
  for (it1 = EmitProb.begin(); it1 != EmitProb.end(); ++it1){
    //unordered_map< unsigned int, gsl_vector* >::iterator it2;
    map< unsigned int, gsl_vector* >::iterator it2;
    for (it2 = (EmitProb[it1->first]).begin(); it2 != (EmitProb[it1->first]).end(); ++it2){
      gsl_vector_free( (EmitProb[it1->first])[it2->first]);
    }
    (EmitProb[it1->first]).clear();
  }
  EmitProb.clear();
  //
  if (get_log==1){
    for (it1 = EmitLog.begin(); it1 != EmitLog.end(); ++it1){
      //unordered_map< unsigned int, gsl_vector* >::iterator it2;
      map< unsigned int, gsl_vector* >::iterator it2;
      for (it2 = (EmitLog[it1->first]).begin(); it2 != (EmitLog[it1->first]).end(); ++it2){
	gsl_vector_free( (EmitLog[it1->first])[it2->first]);
      }
      (EmitLog[it1->first]).clear();
    }
    EmitLog.clear();
  }
  EmitProb_set = 0;
}


void Emission::reset_mask(){
  for (int s=0; s<nLoci; s++){
     for (int l=0; l<sTimes[s]; l++){
       mask[s][l] = 1;
     }
  }
}



void Emission::set_dist(){
  if (times == NULL){
    cout<<"ERROR-1 in Emission::set_dist()\n";
    exit(1);
  }
  dist = new unsigned int * [nLoci];
  total_dist=0.0;
  vector<double> distances;
  for (int s=0; s<nLoci; s++){
    dist[s] = new unsigned int [sTimes[s]];
    dist[s][0] = 0;
    for (int l=1; l <sTimes[s]; l++){
      if (mask[s][l] == 1){
	int k=l;
	unsigned distance = 0;
	while (k>0){
	  distance += times[s][k] - times[s][k-1];
	  if (mask[s][k-1]==1) break;
	  k--;
	}
	if (k==0 && mask[s][k] == 0){
	  dist[s][l] = 0;
	}
	else{
	  dist[s][l] = distance;
	  distances.push_back(double(distance));
	  dist_count[distance] += 1;
	  total_dist += distance;
	}
      }
      else{
	dist[s][l] = 0;
      }
    }
  }
  gsl_sort(distances.data(),1,distances.size());
  median_dist = gsl_stats_quantile_from_sorted_data ( distances.data(), 1, distances.size(), 0.5);
  distances.clear();
  // get the distribution of distances...
  map<unsigned int,int>::iterator it;
  for (it = dist_count.begin(); it != dist_count.end(); ++it){
    if (it->second > 1){
      frequent_dist.insert(pair<unsigned int,int>(it->first,it->second));
    }
  }
  dist_set = 1;
}


void Emission::set_grid(){
  dx = (xmax-xmin) / double(gridSize);
  dy = (ymax-ymin) / double(gridSize);
  for (int i=0; i<= gridSize; i++){
    xgrid[i] = xmin + double(i) * dx;
    ygrid[i] = ymin + double(i) * dy;
  }
}


//get initial values for xmin, xmax etc.
void Emission::init_range(int rp){
  if (mode==1 || mode==2){
    xmin = 0.0;
    ymin = 0.0;
    xmax = 1.0;
    ymax = 1.0;
  }
  else if (mode==3 || mode==4){
    vector<double> allx;
    vector<double> ally;
    double x,y;
    maxRate=0;
    minRate=1.0e6;
    for (int t=0; t<nRep; t++){
      if (rp >= 0 && rp != t) continue;
      for (int s=0; s<nLoci; s++){
	for (int i=0; i < sTimes[s]; i++){
	  if (mask[s][i] == 0) continue;
	  if ( depths[t][s][i] > 0 ){
	    y = double(reads[t][s][i]) / double(depths[t][s][i]);
	    x = (bias != NULL) ? y / bias[s][i] : y;
	    ally.push_back( y );
	    if (bias != NULL)  allx.push_back( x );
	    maxRate = max( maxRate, x);
	    minRate = min( minRate, x);
	  }
	}
      }
    }
    //
    gsl_sort( ally.data(), 1, ally.size());
    double eps = 1.0e-3;
    ymin = 0.01 + gsl_stats_quantile_from_sorted_data ( ally.data(), 1, ally.size(), eps);
    ymax = gsl_stats_quantile_from_sorted_data ( ally.data(), 1, ally.size(), 1.0-eps);
    if (bias != NULL){
      eps = 1.0e-3;
      gsl_sort( allx.data(), 1, allx.size());
      xmin = 0.01 + gsl_stats_quantile_from_sorted_data ( allx.data(), 1, allx.size(), eps);
      xmax = gsl_stats_quantile_from_sorted_data ( allx.data(), 1, allx.size(), 1.0-eps);
    }
    else{
      xmin = ymin;
      xmax = ymax;
    }
  }
  Emission::set_grid();
  range_set=1;
}




//emission probability as a function of total freq x
void Emission::set_EmitProb(int rp){ //
  if (mode == 0){
    printf("ERROR-1 in Emission::set_EmitProb(): mode not set.\n");
    exit(1);
  }
  if (mode==3||mode==4) reflect=0;
  if (range_set == 0) Emission::init_range(rp);
  //delete old stuff
  if (EmitProb_set == 1) Emission::delete_old_Emit();
  //
  for (int t=0; t<nRep; t++){
    if(rp >= 0 && t != rp) continue;
    for (int s=0; s<nLoci; s++){
      unsigned int n1,n2,N;
      for (int i=0; i < sTimes[s]; i++){
	//if (mask[s][i] == 0) continue;
	N  = depths[t][s][i];
	n1 = reads[t][s][i];
	n2 = (reflect == 1) ? N-n1 : n1;	
	//test
	if ( N==0 ){//no observation!
	  if (n1==0) continue;
	  abort();
	}
	else if (EmitProb.count(N) == 1 && EmitProb[N].count(n1) == 1){
	  continue;
	}
	else{//allocate
	  EmitProb[N][n1] = gsl_vector_alloc(gridSize+1);
	  if (n1!=n2) EmitProb[N][n2] = gsl_vector_alloc(gridSize+1);
	  if (get_log == 1){//logarithm of likelihood?
	    EmitLog[N][n1] = gsl_vector_alloc(gridSize+1);
	    if (n1!=n2) EmitLog[N][n2] = gsl_vector_alloc(gridSize+1);
	  }
	}//now get the actual values...
	if (mode == 1){
	  Emission::binomial(N,n1);
	  if (n1!=n2) Emission::binomial(N,n2);
	}
	else if (mode == 2){
	  Emission::beta_binomial(N,n1);
	  if (n1!=n2) Emission::beta_binomial(N,n2);
	}
	else if (mode==3){
	  Emission::poisson(N,n1);
	}
	else if (mode==4){
	  Emission::negative_binomial(N,n1);
	}
	else{
	  printf("ERROR-5 in Emission::set_EmitProb(): mode (%i) not set.\n", mode);
	  exit(1);
	}
      }
    }
  }
  EmitProb_set = 1;
}

//*** BINOMIAL ***
void Emission::binomial(int N, int n){
  double x,p;
  double pre = gsl_sf_lngamma(double(N+1));
  pre -= gsl_sf_lngamma(double(n+1));
  pre -= gsl_sf_lngamma(double(N-n+1));
  for (int j=0; j<=gridSize; j++){
    x = double(j)*dx;
    if (x==0.0){
      p = (n==0) ? 1.0 : 0.0;
    }
    else if (x==1.0){
      p = (n==N) ? 1.0 : 0.0;
    }
    else if (N==0){
      if (n==0){
	p = 1.0;
      }
      else{
	printf("ERROR-2 in Emission::set_EmitProb()\n");
	exit(1);
      }
    }
    else{
      p = gsl_ran_binomial_pdf(n, x, N);
    }
    if (p<0.0 || p!= p){
      printf("ERROR-3 in Emission::set_EmitProb(): %e\n", p);
    }
    gsl_vector_set( EmitProb[N][n], j, p);
    //logarithm of emission probability...
    if (get_log==1){
      if (x>0.0 && x< 1.0){
	gsl_vector_set(EmitLog[N][n], j, pre + double(n)*log(x) + double(N-n)*log(1.0-x));
      }
      else if ( (x==0.0 && n==0) || (x==1.0 && n==N) ){
	gsl_vector_set(EmitLog[N][n], j, 0.0);
      }
      else{
	gsl_vector_set(EmitLog[N][n], j, -1.0e6);
      }
    }	  
  }
}


//*** BETA-BINOMIAL ***
void Emission::beta_binomial(int N, int n){
  double x,p;
  double p0 = gsl_sf_lngamma(double(N+1)) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(double(N-n+1)); 
  for (int j=0; j<=gridSize; j++){
    x = double(j)*dx;
    if (x==0.0){
      gsl_vector_set(EmitProb[N][n], j, n==0 ? 1.0 : 0.0);
      if (get_log == 1) gsl_vector_set(EmitLog[N][n], j, n==0 ? 0.0 : -1.0e6);
    }
    else if (x==1.0){
      gsl_vector_set(EmitProb[N][n], j, n==N ? 1.0 : 0.0);
      if (get_log == 1) gsl_vector_set(EmitLog[N][n], j, n==N ? 0.0 : -1.0e6);
    }
    else{
      p  = p0 + gsl_sf_lnbeta(double(n) + shape*x, double(N-n) + shape*(1.0-x));
      p -= gsl_sf_lnbeta( shape*x,  shape*(1.0-x));
      if (p > 0.0) abort();
      if (get_log==1){
	gsl_vector_set(EmitLog[N][n], j, p);
      }
      gsl_vector_set( EmitProb[N][n], j, exp(p));
    } 
  }	
}

//POISSON
void Emission::poisson(int N, int n){
  double y,p;
  double pre = gsl_sf_lngamma( double(n) + 1.0);
  if (pre != pre) abort();
  for (int j=0; j<=gridSize; j++){
    y = ygrid[j];
    if (y>0.0){
      p  = -y*double(N) + n*(log(y)+log(double(N))) - pre;
    }
    else{
      p = (n==0) ? 0.0 : -1.0e6;
    }
    if (get_log==1){
      gsl_vector_set(EmitLog[N][n], j, p);
    }
    gsl_vector_set(EmitProb[N][n], j, exp(p));
  }
}



//negative binomial emission model
void Emission::negative_binomial(int N, int n){
  //TBD: N==0, x==0
  double y,p;
  double r = double(N)*shape;
  double pre = gsl_sf_lngamma( double(n) + r) - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(r);
  pre += r * log(r);
  //printf("n=%i N=%i shape=%e\n", n, N, shape);
  for (int j=0; j<=gridSize; j++){
    y = ygrid[j];
    if (y>0.0){
      p  = pre - (r + double(n))*log(r+y*double(N)) + double(n)*log(y*double(N));
    }
    else{
      p = (n==0) ? 0.0 : -1.0e6;
    }
    if (get_log==1){
      gsl_vector_set(EmitLog[N][n], j, p);
    }
    gsl_vector_set(EmitProb[N][n], j, exp(p));
    //printf("%.4f %.3e\n", y, p);
  }
  //exit(0);
}

double Emission::get_pval(int rp, int l, int st, double mean){ // rp:replicate ; st: sampling time
  unsigned int n = reads[rp][l][st];
  unsigned int N = depths[rp][l][st];
  if (reflect) n = min(n,N-n);
  double pval=0;
  if (mode==1){//binomial
    if (mean==0.0){
      pval = (n==0) ? 1.0 : 0.0;
    }
    else if (mean < 1.0){
      double p = gsl_cdf_binomial_P(n,mean,N);
      double q = gsl_cdf_binomial_Q(n,mean,N);
      q += gsl_ran_binomial_pdf(n,mean,N);
      pval = min(p,q);
    }
    else{
      pval = (n==N) ? 1.0 : 0.0;
    }
  }
  else if (mode==2){//beta-binomial
    double pre = gsl_sf_lngamma(double(N+1)) - gsl_sf_lnbeta( shape*mean,  shape*(1.0-mean));
    double pn = pre - gsl_sf_lngamma(double(n+1)) - gsl_sf_lngamma(double(N-n+1)); 
    pn += gsl_sf_lnbeta(double(n) + shape*mean, double(N-n) + shape*(1.0-mean));
    double q=0.0,p=0.0;
    if (n >= int(0.5*double(N))){
      for (unsigned int k=n+1;k<=N;k++){
	double y = pre - gsl_sf_lngamma(double(k+1)) - gsl_sf_lngamma(double(N-k+1));
	y += gsl_sf_lnbeta(double(k) + shape*mean, double(N-k) + shape*(1.0-mean));
	q += exp(y);
      }
      p = 1.0 - q - pn;
    }
    else{
      for (unsigned int k=n-1;k>=0;k--){
	double y  = pre - gsl_sf_lngamma(double(k+1)) - gsl_sf_lngamma(double(N-k+1));
	y += gsl_sf_lnbeta(double(k) + shape*mean, double(N-k) + shape*(1.0-mean));
	p += exp(y);
      }
      q = 1.0 - p - pn;
    }
    pval = min(p,q);
  }
  else if (mode==3){//poisson
    double p = gsl_cdf_poisson_P(n,mean*double(N));
    double q = gsl_cdf_poisson_Q(n,mean*double(N));
    q += gsl_ran_poisson_pdf(n,mean*double(N));
    pval = min(p,q);
  }
  else if (mode==4){//negative binomial
    double p = gsl_cdf_negative_binomial_P( n, shape / (mean+shape), double(N)*shape);
    double q = gsl_cdf_negative_binomial_Q( n, shape / (mean+shape), double(N)*shape);
    q += gsl_ran_negative_binomial_pdf( n, shape / (mean+shape), double(N)*shape);
    pval = min(p,q);
  }
  return(pval);
}






bool value_comparer(std::map<int,double>::value_type &i1, std::map<int,double>::value_type &i2){
  return i1.second < i2.second;
}

