//minimization.cpp

// Adapted from http://dx.doi.org/10.1016/j.celrep.2014.04.055

#include "minimization.h"



#define PI 3.1415926
#define LOG2 0.693147

using namespace std;


double find_local_optimum(
			  int nSimplex,
			  gsl_vector**& simplex,
			  gsl_vector * lower,
			  gsl_vector * other,
			  gsl_vector * range,
			  void * params,
			  double (*obj_fn)( const gsl_vector * x, void * p),
			  double prec,
			  int& steps,
			  int verbose
			  ){

    // Here starts the minimizing step...
    
    int iter,max_iter = 200,mint_iter=50; // max_iter max no. iterations for the simplex method; min_iter number of iterations necessary for evaluating additional convergence criteria
    int status=0,ct=0;
    double f1=0, f2=0;
    const gsl_multimin_fminimizer_type * T = NULL;
    gsl_multimin_fminimizer * s = NULL;
    gsl_multimin_function my_func;
    gsl_vector *  x=NULL;

    double size;
    
    // Sigma grid

    fpar * myfpar = static_cast<fpar*> (params);

    double sgm_init = gsl_vector_get(other,0);

    double sgm_final =(double)(myfpar->sigm_f); //

    int   sigmagrid=(double)(myfpar->sigmgrid); //  	

    double step=(double)(sgm_final-sgm_init)/(double)sigmagrid;

    //double sigma_test []={0.0051,0.012,0.0155,0.022,0.0253,0.059}; //

    double sigma_test [sigmagrid];
   

   for (int r=0;r<sigmagrid;r++){

      sigma_test[r]=sgm_init+r*step;

   }

    double fbest=1.0e20; 
    
    double ftest=0.0;
    
    int nvar;
         
    // map the arguments to a single vector
    
    arg_map( nSimplex, simplex, lower, other, range, &x);
    nvar = (int) x->size;
    
    // Allocate best estimate vector
    gsl_vector *  xbest = gsl_vector_alloc(nvar); 
    
    
    // Give the elements of the function-to-minimize object...
    my_func.n       = nvar;
    my_func.f       = obj_fn;
    my_func.params  = params;
    gsl_vector * dx = gsl_vector_alloc(nvar);
    
    
    // Define type of minimization procedure...
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc( T, x->size);
    
    
    // initial displacement vector...
    for (int i=0; i<nvar; i++){
        gsl_vector_set( dx, i, 0.1*fabs( gsl_vector_get( x, i))+0.01);
       
    }

    gsl_vector_memcpy(xbest,x); // copy first value into xbest
    
    
    for (int r=0;r<(sizeof(sigma_test)/sizeof(*sigma_test));r++){
        
        iter=0;
        steps=0;
	ct=0;
        
	

        if(r>0){ // reassign for
        
            // Set ML params once again
            x = NULL;
        
            simplex=NULL;
            lower=NULL;
        
            gsl_vector_set(other,0,sigma_test[r]);
      
            // map the arguments to a single vector
            arg_map( nSimplex, simplex, lower, other, range, &x);
        
            // initial displacement vector...
        
            for (int i=0; i<nvar; i++){
            
                gsl_vector_set( dx, i, 0.1*fabs( gsl_vector_get( x, i))+0.001);
            
               
            }
        }

	
    
    
        gsl_multimin_fminimizer_set( s, &my_func, x, dx);
    
        // Now iterate to find the minimum...
    
        do{
            iter++;
            steps++;

            status = gsl_multimin_fminimizer_iterate(s);
		size = gsl_multimin_fminimizer_size (s);
            if (status) break;
            status = gsl_multimin_test_size( gsl_multimin_fminimizer_size(s), prec);
            // stop-criterion
            f2 = gsl_multimin_fminimizer_minimum(s);
            if (iter > mint_iter && f1 != f2){
                if (f1-f2 > 0.0 && f1-f2 < 1.0){
                    ct++;
                    if (ct==mint_iter) status = 1;
                }
                else{
                    ct=0;
                }
                f1=f2;
            } else if(iter > mint_iter && f1 == f2){

 		ct++;

                    if (ct==mint_iter) status = 1;
			f1=f2;

	    }

            if (verbose==1 && iter % mint_iter == 0){
                printf("%4i ", iter);
                arg_unmap( s->x, nSimplex, simplex, lower, other, range);
                if (simplex != NULL){
                    for (int i=0; i< nSimplex; i++){
                        for (int j=0; j< (int) (simplex[i])->size; j++){
                            printf("%.5e ", gsl_vector_get( simplex[i],j) );
                        }
                        printf("/ ");
                    }
                    printf("/ ");
                }
                if (other != NULL){
                    for (int i=0;i<(int)other->size;i++) printf("%.5e ", other->data[i]);
                }
                printf("%.10e ", s->fval);
		printf("%.10e", size);
                cout<<endl;
            }
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        
        
        ftest = gsl_multimin_fminimizer_minimum(s);
        
        if(ftest<fbest){ // test if new value is better
            
            gsl_vector_memcpy(xbest,gsl_multimin_fminimizer_x(s));
            
            fbest=ftest;
        }

	//cout <<  sigma_test[r] << " " << -ftest <<"\n";
            
    }
        
   
    
    //back-transformation...
    //int err = arg_unmap( gsl_multimin_fminimizer_x(s), nSimplex, simplex, lower, other, range);
    int err = arg_unmap( xbest, nSimplex, simplex, lower, other, range);
    
    //cout << "Other:" << other->data[0] << "\n";
    
    
    if (err==1) abort();
   
    // cleanup...
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(dx);
    gsl_vector_free(xbest); //

   
    
    return(fbest);
}

void spherical_random_step_uniform( double ri, double& rf, double lower,
				    const gsl_vector * anglei, gsl_vector*& anglef,
				    double eps){
  int n = (anglei!=NULL) ? (int) anglei->size + 1 : 1;
  double p=0,na=0,r=0;
  for (int i=0; i<n-1; i++){
    p = double(rand()) / double(RAND_MAX);
    p = (p - 0.5) * eps * 0.5*PI;
    //na = anglei->data[i];
    na = anglei->data[i] + p;
    na = max( na, -na);
    na = min( na, PI-na);
    anglef->data[i] = na;
  }
  p = double(rand()) / double(RAND_MAX);
  p = (p - 0.5)*eps*(1.0-lower);
  r = ri+p;
  r=max(r,2*lower-r);
  rf=min(r,2-r);
}

void simplex_random_step_uniform(const gsl_vector*simplexi, gsl_vector*& simplexf,
				 double lower, double eps){
  int n = (int) simplexi->size;
  gsl_vector * anglei = NULL,* anglef = NULL;
  double ri=0,rf=0,L=0,l=0,p=0;
  if (n>1){
    anglei = gsl_vector_alloc(n-1);
    anglef = gsl_vector_alloc(n-1);
  }
  simplex_to_spherical( simplexi, ri, anglei);  
  L = spherical_to_simplex( ri, anglei, simplexf, 1);
  int acc=0;
  double low=sqrt(lower);
  while (acc==0){
    spherical_random_step_uniform( ri, rf, low, anglei, anglef, eps);
    l = spherical_to_simplex( rf, anglef, simplexf, 1);
    p = double(rand())/double(RAND_MAX);
    if ( l>L || p < exp(l-L)) acc=1;
  }
  if (anglei!=NULL) gsl_vector_free(anglei);
  if (anglef!=NULL) gsl_vector_free(anglef);
}




// map the points on/in simplex to spherical coordinates and then to REAL
void arg_map( 
	     int nSimplex, 
	     gsl_vector**& simplex, 
	     gsl_vector * lower, 
	     const gsl_vector * other, 
	     const gsl_vector * range, 
	     gsl_vector ** x
	      ){
  if (nSimplex == 0 && other == NULL) abort();
  int m=0;
  if (other != NULL) m = (int) other->size;
  int nvar=m;
  if (nSimplex>0){
    if (lower==NULL || (int) lower->size != nSimplex) abort();
    for (int t=0; t<nSimplex; t++){
      int n = (simplex[t])->size;
      nvar += (lower->data[t] < 1.0) ? n : n-1;//inside or on the simplex?
    }
  }
  if (nvar == 0) abort();
  if (*x != NULL) gsl_vector_free(*x);
  *x = gsl_vector_alloc(nvar);//allocate here
  gsl_vector * angle = NULL;
  double radial=0;
  // set transformed variables...
  int ct=0;
  for (int t=0; t<nSimplex; t++){    
    int n = (simplex[t])->size;
    if (n==1 && lower->data[t] < 1.0){//only one freq in [lower_t,1]
      radial = gsl_vector_get(simplex[t],0);
      gsl_vector_set( *x, ct, sqrt(radial));
      ct++;
    }
    else if (n>1){
      if (angle != NULL) gsl_vector_free(angle);
      angle = gsl_vector_alloc(n-1);
      simplex_to_spherical( simplex[t], radial, angle);
      if (lower->data[t] < 1.0){//inside
	gsl_vector_set( *x, ct, radial);
	ct++;
      }
      for (int j=0; j<n-1; j++){
	gsl_vector_set( *x, ct, gsl_vector_get(angle,j));
	ct++;
      }
    }
  }
  if (angle != NULL) gsl_vector_free(angle);
  if (ct != nvar - m) abort();
  // remaining variables...
  for (int k=0; k<m; k++){
    double R = gsl_vector_get(range,k);
    if (R > 0.0){
      gsl_vector_set( *x, ct, logify( gsl_vector_get(other,k), R));
    }
    else if (R==0.0){
      gsl_vector_set( *x, ct, log( gsl_vector_get(other,k) ));
    }
    else{
      gsl_vector_set( *x, ct, gsl_vector_get(other,k) );
    }
    ct++;
  }
}


int arg_unmap( 
	      const gsl_vector * x, 
	      int nSimplex, 
	      gsl_vector**& simplex, 
	      gsl_vector * lower,
	      gsl_vector * other,
	      const gsl_vector * range
	       ){
  if ( nSimplex > 0 ){
    if (lower == NULL || (int) lower->size != nSimplex) abort();
  } 
  // test dimensionality...
  int test=0;
  for (int t=0; t<nSimplex; t++){
    int n = (int) (simplex[t])->size;
    test += (lower->data[t] < 1.0) ? n : n-1;//inside or on simplex?
  }
  if (other != NULL) test += (int) other->size;
  int nvar = (int) x->size;
  if (test != nvar) abort();
  //...passed test
  double radial,val;
  gsl_vector * angle = NULL;
  int ct=0;
  for (int t=0; t<nSimplex; t++){ 
    if (lower->data[t] < 1.0){//inside simplex
      radial = x->data[ct];
      if (radial < sqrt(lower->data[t]) || radial > 1.0) return(1);//test radial bounds
      ct++;
    }
    else{//or on the surface
      radial = 1.0;
    }
    int n = (int) (simplex[t])->size;
    if (n>1){
      if (angle != NULL) gsl_vector_free(angle);
      angle = gsl_vector_alloc(n-1);
      for (int j=0; j<n-1; j++){
	gsl_vector_set( angle, j, x->data[ct]);
	if (angle->data[j] < 0.0 || angle->data[j] > 0.5*PI){//angle range test
	  gsl_vector_free(angle);
	  return(1);
	}
	ct++;
      }
      spherical_to_simplex( radial, angle, simplex[t], 0);
    }
    else{// n2==1
      gsl_vector_set( simplex[t], 0, pow(radial,2));
    }
  }
  if (angle != NULL) gsl_vector_free(angle);
  if (other == NULL) return(0);
  for (int k=0; k<(int) other->size; k++){
    double R = gsl_vector_get(range,k);
    val = x->data[ct];
    if (R>0.0){
      if ( val < -20.0 || val > 20.0) return(1);
      other->data[k] = delogify( val, R);
    }
    else if ( R == 0.0){
      if (val<-20.0 || val > 20.0) return(1);
      other->data[k] = exp(val);
    }
    else{
      other->data[k] = val;
    }
    ct++;
  }
  return(0);
}


double logify( double x, double R){
  if (x==0.0){
    return(-20.0);
  }
  else if (x<0.5*R){
    return( log(x) - log(0.5*R));
  }
  else if (x<R){
    return( log(0.5*R) - log(R-x));
  }
  else{
    return(20.0);
  }
}

double delogify(double y, double R){
  if (y<0.0){
    return( 0.5*R*exp(y) );
  }
  else{
    return( R-0.5*R*exp(-y) );
  }
}



double spherical_to_simplex( double radial, const gsl_vector * angle, gsl_vector *& simplex, int getLJD){
  int n = simplex->size;
  if (angle != NULL && n-1 != (int) angle->size) abort();
  double LJD=0.0;
  if (getLJD) LJD = double(n)*log(2.0) + (2.0*double(n)-1.0)*log(radial);
  if (n>1){
    double * sinV = new double[n-1];
    double * cosV = new double[n-1];
    for ( int i =0; i<n-1; i++){
      sinV[i] = sin(gsl_vector_get(angle,i));
      cosV[i] = cos(gsl_vector_get(angle,i));
      if (getLJD){
	LJD += log(cosV[i]);
	LJD += (2.0*double(n-i)-1.0)*log(sinV[i]);
      }
    }
    double mem=1;
    for ( int i=0; i<n-1; i++){
      gsl_vector_set( simplex, i, mem * cosV[i]);
      mem = mem*sinV[i];
    }
    gsl_vector_set( simplex, n-1, mem);
    for ( int i=0; i<n; i++) simplex->data[i] = pow(simplex->data[i],2);
    double norm = gsl_blas_dasum(simplex);
    if (norm<=0.0) abort();
    gsl_vector_scale( simplex, pow(radial,2) / norm);
    delete [] sinV;
    delete [] cosV;
  }
  else{
    simplex->data[0] = pow(radial,2);
  }
  return(LJD);
}




void simplex_to_spherical( const gsl_vector * simplex, double& radial, gsl_vector*& angle){
  int n = simplex->size;
  if (angle != NULL && n-1 != (int) angle->size) abort();
  if (n>1){
    double partSum = gsl_vector_get( simplex, n-1) + gsl_vector_get( simplex, n-2);
    double val=0;
    double x = gsl_vector_get( simplex, n-1);
    if (x > 0.0){
      val =  (sqrt(partSum) + sqrt(gsl_vector_get( simplex, n-2))) / sqrt(x);
      gsl_vector_set( angle, n-2,  PI - 2.0*atan(val));
    }
    else if (x == 0.0){
      gsl_vector_set( angle, n-2, 0.0);
    }
    else abort();
    for ( int i=3; i<n+1; i++){
      val = sqrt( gsl_vector_get(simplex,n-i) / partSum);
      val = PI/2.0 - atan(val);
      if ( val<0 || val>0.5*PI) abort();
      gsl_vector_set( angle, n-i, val);
      partSum += gsl_vector_get( simplex, n-i);
    }
    radial = sqrt(partSum);
  }
  else{
    radial = sqrt(simplex->data[0]);
  }
}

