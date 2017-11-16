//common-functions.cpp

// Adapted from http://dx.doi.org/10.1016/j.celrep.2014.04.055

#include "common-functions.h"
#include "emission.h"

#define PI 3.1415926
#define LOG2 0.693147

using namespace std;


//get general dimensions of a data set for...
void get_dims( const char * data_fn, 
	       int& nRep,
	       vector<int>& Loci,
	       vector<int>& sTimes,
	       int keep
	       ){
  ifstream data_ifs;
  string line;
  stringstream line_ss;
  data_ifs.open( data_fn, ios::in);
  if (data_ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", data_fn);
    exit(1);
  }
  sTimes.clear();
  Loci.clear();
  int ct=0,l,r,d;
  int locs=0,old=-1,nT=0;
  while( data_ifs.good()){
    line.clear();
    getline( data_ifs, line);
    if (line.empty()) break;
    if (line[0] == '#') continue;
    line_ss.clear();
    line_ss.str(line);

    //check first entry for nRep
    if (old == -1 && ct == 0){
      line_ss >> locs >> l; 

      while(line_ss >> r >> d){
	nT++;
      }
      line_ss.clear();
      line_ss.str(line);      
    }

    line_ss >> locs >> l;
    if (locs != old ){//new locus encounter     
      if (ct>0){
	sTimes.push_back(ct);
	Loci.push_back(old);
      }
      ct=0;
    }
    old=locs;
    r = 0;
    for( int t=0;t<nT; t++){
      line_ss >> r >> d;
      if (r>0) break;
    }
    if (keep || r>0) ct++;
  }
  if (ct>0){
    sTimes.push_back(ct);
    Loci.push_back(old);
  }
  nRep = nT;
  data_ifs.close();
}


// read in data: expects columns to be "locus time (depth reads)^x"
void get_data( const char * data_fn, Emission * myEmit){
  ifstream data_ifs;
  string line;
  stringstream line_ss;
  data_ifs.open( data_fn, ios::in);
  if (data_ifs.fail()){
    printf("ERROR: file %s cannot be opened.\n", data_fn);
    exit(1);
  }
  int ct=0,l;
  int locs=0,old=-1, sample=0;
  int d,r, keep=0;
 
    //now collect all data...
  
    while( data_ifs.good()){
    
      line.clear();
    
      getline( data_ifs, line);
    
      if (line.empty()) break;
    
      if (line[0] == '#') continue;
    
      line_ss.clear();
      line_ss.str(line);
      line_ss >> locs >> l;//locus and sampling time
    
      if (locs != old){
          if (myEmit->Loci.count(locs) == 0){
              cout<<"ERROR in get_data(): locus "<<locs<<" was not expected"<<endl<<line<<endl;
	
              for (int s=0; s<myEmit->nLoci; s++){
                  printf("sample %i = locus %i, idx = %i\n",s+1, myEmit->Locus[s], myEmit->idx_of[myEmit->Locus[s]]);
              }
              exit(1);
          }
      
          sample = myEmit->idx_of[locs];
          ct  = 0;
          old = locs;
      }
    
      if (ct >= myEmit->sTimes[sample]) continue;
    
      keep = 0;
    
      for (int t=0; t<myEmit->nRep; t++){
      
          myEmit->times[sample][ct] = l;//set locus
      
          line_ss >> r >> d;//get read and depth
      
          if (d == 0 && r > 0){
	
              printf("ERROR: depth = 0 in locus %i time %i\n", locs, l);
	
              cout<<line<<endl;
	
              exit(1);
      
          }
     
          //set read and depth
      
          myEmit->reads[t][  myEmit->idx_of[locs] ][ct] = r;
      
          myEmit->depths[t][ myEmit->idx_of[locs] ][ct] = d;
      
          if (r>0) keep=1;
    
      }
    
      if (keep || myEmit->connect) ct++;
  
  }
  
    data_ifs.close();
  
    // set the distances between loci
  
    for (int t=0; t<myEmit->nRep; t++){
    
        myEmit->set_dist();
  
    }

}

//***Mean and Variance function***
double get_mean(gsl_vector * dist, double xmin, double xmax){
  double mean=0.0,P1,P2;
  int n = (int) dist->size;
  double dx = (xmax - xmin) / double(n-1);
  for (int i=0; i < n-1; i++){
    P1 = gsl_vector_get(dist,i);
    P2 = gsl_vector_get(dist,i+1);
    mean += 3.0*(P1+P2)*(xmin+double(i)*dx) + (P1+2.0*P2)*dx;
  }
  mean = mean * dx / 6.0;
  return(mean);
}

double get_var(gsl_vector * dist, double xmin, double xmax, double mean){
  double var=0.0, P1, P2,dev;
  int n = (int) dist->size;
  double dx = (xmax - xmin) / double(n-1);
  for (int i=0; i<n-1; i++){
    P1 = gsl_vector_get(dist,i);
    P2 = gsl_vector_get(dist,i+1);
    dev = xmin + double(i)*dx - mean;
    var += (P1+3.0*P2)*dx*dx + 4.0*(P1+2.0*P2)*dev*dx + 6.0*(P1+P2)*pow(dev,2);
  }
  var = var*dx/12.0;
  return(var);
}
