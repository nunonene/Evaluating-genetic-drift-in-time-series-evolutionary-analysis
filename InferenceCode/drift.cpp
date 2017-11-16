//drift.cpp

// Adapted from http://dx.doi.org/10.1016/j.celrep.2014.04.055

//own headers...
#include "emission.h"
#include "log-space.h"
#include "drift.h"
#include "common-functions.h"
#include "MatOp.h" 

// GSL headers
#include <gsl/gsl_matrix.h>


#define PI 3.1415926

// ========== Constructor =============
Drift::Drift( Emission * emit, int rp){
  myEmit    = emit;
  Rep       = rp;
  nLoci     = myEmit->nLoci;
  sTimes    = myEmit->sTimes;
  dist      = myEmit->dist;
  times     = myEmit->times;
  mask      = myEmit->mask;
  mode      = myEmit->mode;
  gridSize  = myEmit->gridSize; 
  sigma     = -1.0;
  wTotal    = 0;
 
  // matrices...
  alpha  = new gsl_matrix * [nLoci]; // For combined forward-backward/predict-update 
  gamma  = new gsl_matrix * [nLoci]; // For combined forward-backward/predict-update 
  proposal = gsl_vector_alloc(gridSize+1);
  gsl_vector_set_all( proposal, 1.0);// uniform proposal distribution on [0,1]
  
  for (int s=0; s<nLoci; s++){
   
    alpha[s] = NULL;
    gamma[s] = NULL;
  }
  
  pstay_set  = 0;
  save_alpha = 0;
  DiffProp = NULL;
  DiffProp_set = 0;
    
  DiffProp_opt=1; // default option is for neutral Wright-Fisher drift model propagation  
}

Drift::~Drift(){
  for (int s=0; s<nLoci; s++){
   
    if (alpha[s] != NULL) gsl_matrix_free(alpha[s]);
    if (gamma[s] != NULL) gsl_matrix_free(gamma[s]);
  }
  delete [] alpha;
 
  
  delete [] gamma;
  Drift::reset_DiffProp();
}

double Drift::get_total_llh(){
  save_alpha  = 0;
  if (myEmit->EmitProb_set == 0){
    myEmit->set_EmitProb(Rep);
  }
  gsl_vector_set_all(proposal, 1.0/(myEmit->xmax - myEmit->xmin));
  	    
    if (DiffProp_set == 0) {
        if(DiffProp_opt==0){
            Drift::get_DiffProp();
        }else{
            Drift::get_DiffPropWF();
        }
    }
    
  int nl;
  total_llh   = 0.0;

  // All loci in parallel if structure is available:

#ifdef _OPENMP
  int nt = min( nLoci, omp_get_max_threads());
#pragma omp parallel for schedule( dynamic, 1) default(shared) num_threads(nt)
#endif
  for (nl=0; nl<nLoci; nl++){
    double llh = Drift::do_Fwd(nl);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      total_llh += llh;
    }
  }
  return(total_llh);
}

void Drift::get_posterior(int nl){
  if (alpha[nl] != NULL) gsl_matrix_free(alpha[nl]);
  if (gamma[nl] != NULL) gsl_matrix_free(gamma[nl]);
  alpha[nl] = gsl_matrix_alloc(sTimes[nl],gridSize+1);
  gamma[nl] = gsl_matrix_alloc(sTimes[nl],gridSize+1);
  save_alpha  = 1;
  if (myEmit->EmitProb_set == 0){
    myEmit->set_EmitProb(Rep);
  }
  gsl_vector_set_all(proposal, 1.0/(myEmit->xmax - myEmit->xmin));
 

    if (DiffProp_set == 0) {
        if(DiffProp_opt==0){
            Drift::get_DiffProp();
        }else{
            Drift::get_DiffPropWF();
        }
    }
  Drift::do_Fwd(nl);
  Drift::do_Bwd(nl);
  gsl_matrix_free(alpha[nl]);
  alpha[nl] = NULL;
}

void Drift::reset_DiffProp(){
  if (DiffProp != NULL){
    for (int i=0; i< (int) myEmit->frequent_dist.size(); i++){
      gsl_matrix_free(DiffProp[i]);
    }
    delete [] DiffProp;
    DiffProp = NULL;
  }
}

void Drift::get_DiffProp(){
  if (myEmit->dist_set == 0){
    cout<<"ERROR-1 in Drift::set_DiffProp()\n";
    exit(1);
  }
  // allocate...
  map<unsigned int,int>::iterator it;
  if (DiffProp == NULL){
    DiffProp = new gsl_matrix * [myEmit->frequent_dist.size()];
    int i=0;
    for (it = myEmit->frequent_dist.begin(); it != myEmit->frequent_dist.end(); ++it){
      DiffProp[i] = gsl_matrix_alloc(gridSize+1,gridSize+1);
      position.insert(pair<unsigned int,int>(it->first,i));
      i++;
    }
  }
  // now set the matrices...
  is_identity.clear();
  int i=0;
  for (it=myEmit->frequent_dist.begin(); it != myEmit->frequent_dist.end(); it++){
    int is_id = Drift::set_DiffProp( DiffProp[i],  sigma*sqrt(double(it->first)));
    is_identity.push_back(is_id);
    i++;
  }
  DiffProp_set = 1;
}

void Drift::get_DiffPropWF(){
	if (myEmit->dist_set == 0){
		cout<<"ERROR-1 in Drift::set_DiffPropWF()\n";
		exit(1);
	}
	// allocate...
	map<unsigned int,int>::iterator it;
	if (DiffProp == NULL){
		DiffProp = new gsl_matrix * [myEmit->frequent_dist.size()];
		int i=0;
		for (it = myEmit->frequent_dist.begin(); it != myEmit->frequent_dist.end(); ++it){
			DiffProp[i] = gsl_matrix_alloc(gridSize+1,gridSize+1);
			position.insert(pair<unsigned int,int>(it->first,i));
			i++;
		}
	}
	// now set the matrices...
	is_identity.clear();
	int i=0;
       
    
	for (it=myEmit->frequent_dist.begin(); it != myEmit->frequent_dist.end(); it++){
        
  
		int is_id = Drift::set_DiffPropWF( DiffProp[i], 1/(sigma*sigma), double(it->first));
		is_identity.push_back(is_id);
		i++;		
	}
	
	DiffProp_set = 1;
}


// ======== Gaussian drift model with absorbing boundaries at 0 and 1 ===============

int Drift::set_DiffProp( gsl_matrix * propagator, double sd){

  
    double dx = myEmit->dx;
  
    if (3.0*sd <= dx) return(1);// with the current resolution, propagator is Dirac-Delta!
    
    int range = 3 * ceil(sd / dx);// this means, only fluctuations up to 3 sigma are possible
    
     
    if (2*range > gridSize+1){
        sd = dx * double(gridSize) / 6.0;
        range = 3 * ceil(sd / dx);
    }
  
    gsl_vector * gauss = gsl_vector_alloc(2*range+1);
 
    for (int i=0; i<2*range+1; i++){
      gsl_vector_set(gauss,i, gsl_ran_gaussian_pdf( double(i-range)*dx, sd));
    }
 
    gsl_matrix_set_zero(propagator);
 
    double val = 0, norm=0, val_boundary = 0;
  
    gsl_vector_view row;
  
    for (int i=0; i<= gridSize; i++){
        norm = 0.0;

	val_boundary=0; 

	if(i==0 || i==gridSize){ 
		gsl_matrix_set(propagator,i,i,1);
		norm=1;
		
	}else{

      	  for (int j=i-range; j<=i+range; j++){

      		  
			if ( j<0 || j>gridSize){ 

			val_boundary =val_boundary+ gsl_vector_get( gauss, j-i+range);

      			continue;

			}else{

        	    	val = gsl_vector_get( gauss, j-i+range);


			}
      
        	  
            
		if(j==0){
		
        	    gsl_matrix_set(propagator,i,j,val+val_boundary); 
        	
		}else if(j>0 || j<=gridSize){
        	    gsl_matrix_set(propagator,i,j,val);
		}
           
            
        	norm += (j==0) ? val+val_boundary : val; 

		if(j==0) val_boundary=0;

        	
	  }

		gsl_matrix_set(propagator,i,gridSize,gsl_matrix_get(propagator,i, gridSize)+val_boundary);
      
	
 		norm +=val_boundary;
        }
	row = gsl_matrix_row(propagator,i);
        norm *= dx;
    
        if (norm <=0.0 ||norm!=norm){
            cout<<"ERROR in Drift::set_DiffProp()\n";
            printf("sd=%e dx=%e\n",sd,dx);
            exit(1);
        }
        gsl_vector_scale(&row.vector,1.0/norm);
    }
  
    gsl_vector_free(gauss);
    return(0);
}


// ==================== Wright-Fisher drift model ==========================

void setWFMatpropagator(gsl_matrix *FullMat,double s, double t, int Nlocal){
    
    gsl_matrix *m = gsl_matrix_alloc (Nlocal+1,Nlocal+1); // every matrix is a gsl_matrix
    gsl_matrix *I = gsl_matrix_alloc (Nlocal+1,Nlocal+1);
    gsl_matrix_set_identity (I);
        
    
    MatOp MP; // instantiate class with exponentiation method
    
    double temp;
    double ps;
    int nn;
    
    
    // Calculate original transition matrix
    // ====================================
    for (int i=0;i< (Nlocal+1);i++){
        
        double sum=0;
        
        for (int j=0;j<(Nlocal+1);j++){
            
            ps=(i/(double)(Nlocal))+((s)*((i/(double)(Nlocal)) * (1-(i/(double)(Nlocal)))));
            
            nn=j;
            
            temp=gsl_ran_binomial_pdf (nn,ps,Nlocal);
            
            
            if(isnan(temp)){
                temp=1;
            }
            
            gsl_matrix_set (m, i, j,temp);
            
            
            sum+=temp;
            
        }
        
        // Row-wise normalization
        // ======================
        
        for (int k=0;k<(Nlocal+1);k++){
            
            gsl_matrix_set(m,i,k,gsl_matrix_get(m,i,k)/sum);
            
        }
        
        
        
    }    
    
    // Exponentiate transition matrix (squares method)
    // ==============================================
    
    
    MP.PowMatBLAS(t,Nlocal,I,m); // use gsl_BLAS
    
    gsl_matrix_memcpy(FullMat,m);
    
    // Free matrices
    gsl_matrix_free (m);
    gsl_matrix_free (I);
}


int Drift::set_DiffPropWF( gsl_matrix * propagator, double N, double dt){
    
    double dx = myEmit->dx;
    
    
    gsl_vector * theta=  gsl_vector_alloc(1);
    gsl_vector_set(theta,0,N);
    gsl_matrix_set_zero(propagator);
    
    
    if(N>1000){ // Use pre computed matrices (with OpenBlas or MKL) if 2N> 1000
           
        gsl_matrix* FullMat=gsl_matrix_alloc(gridSize+1,gridSize+1);
      
        int Nlocal=(int) floor(N);
    
        double s=0.0;
    
        // Read gsl_matrix from file
    
        ostringstream convert;   // stream used for the conversion
    
        convert << fixed <<(int)Nlocal;     
        std::string PopSizeString = convert.str(); 
        convert.str(""); 
        convert.clear(); 
    
        convert<<(int)dt;
        std::string PowString= convert.str(); 
        convert.str("");
        convert.clear();
    
        convert<<(int)s;
        std::string SigmaString= convert.str(); 
        convert.str("");
        convert.clear();
    
        std::string MatMultOption="OpenBLAS"; 
    
	// read matrices from MatPow folder in current directory

         std::string filename = std::string("./MatPow/MatPow_Pow")+PowString+std::string("_PopSize")+PopSizeString+std::string("_GridSize400_Sigma")+SigmaString+std::string("_")+MatMultOption+std::string(".dat");
        
        FILE  *in_file=fopen(filename.c_str(),"rb");
    
    
    
        if(in_file==NULL)
        {
            Nlocal=(int)floor((N+50.00)/100.00)*100; // above 10000 only every 100 matrix powers are available
            convert <<Nlocal;   
            std::string PopSizeString = convert.str(); 
            convert.str("");
            convert.clear(); 
       
       // read matrices from MatPow folder in current directory

	    std::string filename = std::string("./MatPow/MatPow_Pow")+PowString+std::string("_PopSize")+PopSizeString+std::string("_GridSize400_Sigma")+SigmaString+std::string("_")+MatMultOption+std::string(".dat");
    

	 in_file=fopen(filename.c_str(),"rb");
        
        }
           
        gsl_matrix_fread (in_file, FullMat);
    
        fclose(in_file);
    
        // Continue with original program
        double norm=0.0;
        double r=0.0;
    
        gsl_vector_view row;
    
    
        for (int i=0; i<=gridSize; i++){ 
        
            norm = 0.0;
        
            for (int j=0; j<=gridSize; j++){
            
                r=gsl_matrix_get(FullMat,i,j);
            
                norm += r;
            
            }
        
            row = gsl_matrix_row(FullMat,i);
            norm *= dx;
        
            if (norm <=0.0 ||norm!=norm){
                cout<<"ERROR in Drift::set_DiffPropWFmat()\n";
                printf("Z=%e N=%d xinit=%e dx=%e dt=%e\n",norm,Nlocal,i*dx,dx,dt);
                exit(1);
            }
            gsl_vector_scale(&row.vector,1.0/norm);
        }
    
        gsl_vector_free(theta);
    
        gsl_matrix_memcpy(propagator,FullMat);
    
        gsl_matrix_free(FullMat);
    
    }else{ // Exponentiate matrices on the fly using gsl_blas if 2N<=1000
        
        double dx = myEmit->dx;
        
        gsl_vector * theta=  gsl_vector_alloc(1);
        gsl_vector_set(theta,0,N);
        gsl_matrix_set_zero(propagator);
        
        int Nlocal=(int)floor(N+1);
        
        gsl_matrix* FullMat=gsl_matrix_alloc(Nlocal+1,Nlocal+1);
        
       
        
        double s=0.0;
        
               
        setWFMatpropagator(FullMat,s,(int)dt,Nlocal); 
        
        
        // Continue with original program
        
        double val = 0, norm=0;
        
        gsl_vector_view row;
        
        vector<int> a;
        vector<int> b;
        
        double dist=0.0;
        
        double disti=0.0;
        
        double distj=0.0;
        
        double Pn=0.0;
        double Pd=0.0;
        
        double r=0.0;
        
        for (int i=0; i<=gridSize; i++){ 
            
            norm = 0.0;
            
            for (int j=0; j<=gridSize; j++){
               
                
                if(i==0 && j==0){
                    
                    r=gsl_matrix_get(FullMat,i,j);
                    
                }else if (i==gridSize && j==gridSize){
                    
                    r=gsl_matrix_get(FullMat,Nlocal,Nlocal);
                    
                }else{
                    
                    if (i==gridSize && j==0){
                        r=gsl_matrix_get(FullMat,Nlocal,j);
                    } else if (i==0 && j==gridSize){
                        r=gsl_matrix_get(FullMat,i,Nlocal);
                        
                    }else{
                        if(i==0){
                            a.push_back((int)floor(i*dx*Nlocal));
                            
                        }else if (i==gridSize){
                            a.push_back(Nlocal);
                        }else{
                            a.push_back((int)floor(i*dx*Nlocal));
                            a.push_back((int)floor(i*dx*Nlocal)+1);
                        }
                        
                        if(j==0){
                            b.push_back((int)floor(j*dx*Nlocal));
                        }else if (j==gridSize){
                            b.push_back(Nlocal);
                        }else{
                            b.push_back((int)floor(j*dx*Nlocal));
                            b.push_back((int)floor(j*dx*Nlocal)+1);
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
                                    
                                  
                                    
                                    val=gsl_matrix_get(FullMat,a[k],b[l]);
                                    
                                    disti=((double(a[k])/double(Nlocal))-(i*dx))*((double(a[k])/double(Nlocal))-(i*dx));
                                    distj=((double(b[l])/double(Nlocal))-(j*dx))*((double(b[l])/double(Nlocal))-(j*dx));
                                    
                                    if(disti==0 && distj==0){
                                        r=val;
                                        flag=1;
                                    }else{
                                        dist=sqrt(disti+distj);
                                        
                                        Pn=Pn+(pow(static_cast<double>(1.0/(dist)),static_cast<int>(2))*val);
                                        
                                        Pd=Pd+(pow(static_cast<double>(1.0/(dist)),static_cast<int>(2)));
                                        
                                 
                                    }
                                }
                            }
                        }
                        if(!flag) r=Pn/Pd;
                        flag=0;
                    }
                }
                
                gsl_matrix_set(propagator,i,j,r);
           
                norm += r;
                
                a.clear();
                b.clear();
            }
            
            row = gsl_matrix_row(propagator,i);
            norm *= dx;
            if (norm <=0.0 ||norm!=norm){
                cout<<"ERROR in Drift::set_DiffPropWFmat()\n";
                printf("Z=%e N=%d xinit=%e dx=%e dt=%e\n",norm,Nlocal,i*dx,dx,dt);
                exit(1);
            }
            gsl_vector_scale(&row.vector,1.0/norm);
        }
        gsl_vector_free(theta);
        gsl_matrix_free(FullMat);
        return(0);

        
        
    }
    
    return(0);
}




int Drift::predict(gsl_vector * prior, gsl_vector * post, 
			   gsl_matrix*& DP, gsl_matrix**& DP_pt, int s, int l){
  int is_id=0;

  //check whether DiffProp is pre-computed
  if (  (myEmit->frequent_dist).count(dist[s][l]) == 0 ){ // not found
	  
      if(DiffProp_opt==0){

          is_id = Drift::set_DiffProp( DP, sigma*sqrt(double(dist[s][l])));
      }else{
          is_id = Drift::set_DiffPropWF( DP, 1/(sigma*sigma), double(dist[s][l]));
      }
       
    DP_pt = &DP;
  }
  else{//exists
    if ( is_identity[ position[dist[s][l]] ] == 1 ){
      is_id = 1;
      DP_pt = &DP;
    }
    else{
      is_id = 0;
      DP_pt = &( DiffProp[ position[dist[s][l]] ] );
    }
  }




  //apply  propagator
  if (is_id == 0){
    gsl_vector * mem = gsl_vector_alloc(gridSize+1);
    post->data[0]        *= 0.5;
    post->data[gridSize] *= 0.5;
    gsl_blas_dgemv(CblasTrans, myEmit->dx, *DP_pt, post, 0.0, mem);


    gsl_vector_memcpy( post, mem);
    gsl_vector_free(mem);  
  }
  else{//increase variance by hand (exact for gaussian distribution)

 	if(DiffProp_opt==0){
    		double mn  = get_mean( post, myEmit->xmin, myEmit->xmax);
   		double var = get_var(  post, myEmit->xmin, myEmit->xmax, mn);
    		double e = var / ( var + 100*pow(sigma,2)*double(dist[s][l]));

    		for (int i=0; i<=gridSize; i++){
     			 double p = post->data[i];
     		 if (p>0.0) post->data[i] = pow(p,e);
   		 }
   		 double norm = gsl_blas_dasum(post) - 0.5*(post->data[0] + post->data[gridSize]);
    		norm *= myEmit->dx;
    		gsl_vector_scale(post,1.0/norm);
   	 }else{

		cout << " Analytical increase of variance not yet implemented for WF\n";
		exit(0);

	}
  }

  gsl_vector_scale( prior, 0.0);
  gsl_vector_scale( post,  1.0);
  gsl_vector_add( prior, post);

  return(is_id);
}

// ===================== Forward step ========================

double Drift::do_Fwd(int s){
  // prepare fwd...
  gsl_vector * eprob = gsl_vector_alloc(gridSize+1);
  gsl_vector * prior = gsl_vector_alloc(gridSize+1);
  gsl_vector * post  = gsl_vector_alloc(gridSize+1);
  //gsl_vector * mem   = gsl_vector_alloc(gridSize+1);
  gsl_matrix * DP    = gsl_matrix_alloc(gridSize+1,gridSize+1);
  gsl_matrix ** DP_pt = NULL;
  double norm;
  double llh=0.0;
  int get_log=0;

//cout<<"Rep:" << Rep<<"\n";
  // Forward Pass
  for (int l=0; l<sTimes[s]; l++){
    if (mask[s][l] == 0) continue;
    int N = myEmit->depths[Rep][s][l];
    int n = myEmit->reads[Rep][s][l];




    gsl_vector_memcpy( prior, proposal);

    if (dist[s][l] > 0) Drift::predict( prior, post, DP, DP_pt, s, l);


    //emission probability
    if (N>0){//if observation
    
          gsl_vector_memcpy( eprob, myEmit->EmitProb[N][n] );
          if (myEmit->reflect && n != N-n){
              gsl_vector_add( eprob, myEmit->EmitProb[N][N-n]);
             
          }
    
        gsl_vector_mul( prior, eprob); // at this time it is the posterior!
        norm  = gsl_blas_dasum(prior);


        norm -= 0.5*(prior->data[0] + prior->data[gridSize]);
        norm *= myEmit->dx;
        if (norm <=0 || norm != norm){
            cout<<"ERROR\n";
            abort();
        }
    
        gsl_vector_scale( prior, 1.0 / norm);
        llh += log(norm);// get part of the total log-likelihood
        
    }
    
      gsl_vector_memcpy(post,prior);
    
      if (save_alpha == 1){
      
          gsl_matrix_set_row(alpha[s], l, post);// save forward variable
    
      }

  }

  
    // clean up...
  gsl_vector_free(eprob);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_matrix_free(DP);


  return(llh);
}


// ==================== Backward Pass ==========================
void Drift::do_Bwd(int s){
  // prepare bwd...
  gsl_vector * eprob = gsl_vector_alloc(gridSize+1);
  gsl_vector * prior = gsl_vector_alloc(gridSize+1);
  gsl_vector * post  = gsl_vector_alloc(gridSize+1);
  gsl_vector * beta  = gsl_vector_alloc(gridSize+1);
  gsl_vector * last_beta  = gsl_vector_alloc(gridSize+1);
  gsl_vector * last_eprob = gsl_vector_alloc(gridSize+1);
  gsl_vector * mem   = gsl_vector_alloc(gridSize+1);
  gsl_vector * mem2  = gsl_vector_alloc(gridSize+1);
  gsl_matrix * DP    = gsl_matrix_calloc(gridSize+1,gridSize+1);
  gsl_matrix ** DP_pt = NULL;
  double x,y,norm;
  int get_log=0;
  int is_id=1;
  int last=-1;
  for (int l = sTimes[s]-1; l>=0; l--){
    if (mask[s][l]==0) continue;
    gsl_vector_memcpy(prior,proposal);
    if (last>0) is_id = Drift::predict( prior, post, DP, DP_pt, s, last);
    gsl_vector_memcpy( beta, prior);
    // get gamma, i.e. the total posterior probability vector
    gsl_vector_view alph = gsl_matrix_row(alpha[s],l);
    gsl_vector_memcpy( post, beta);
    gsl_vector_mul( post, &alph.vector);
    norm  = gsl_blas_dasum(post);
    norm -= 0.5*(post->data[0] + post->data[gridSize]);
    norm *= myEmit->dx;
    gsl_vector_scale( post, 1.0 / norm);
    // posterior on-site sojourn probability.
    gsl_matrix_set_row( gamma[s], l, post);
    // emission probability
    int N = myEmit->depths[Rep][s][l];
    int n = myEmit->reads[Rep][s][l];
    if (N>0){
    
          gsl_vector_memcpy( eprob, myEmit->EmitProb[N][n] );
          if (myEmit->reflect && n != N-n){
              gsl_vector_add( eprob, myEmit->EmitProb[N][N-n]);
          }
     
    
      // get posterior update for the next step...
      gsl_vector_mul( prior, eprob);// now it is the posterior!
      norm = gsl_blas_dasum(prior);
      norm -= 0.5*(prior->data[0] + prior->data[gridSize]);
      norm *= myEmit->dx;
      gsl_vector_scale(prior, 1.0 / norm);
    }
    gsl_vector_memcpy( post, prior);
    

    gsl_vector_memcpy(last_beta,beta);
    gsl_vector_memcpy(last_eprob,eprob);
    last = l;
  }
 
  //clean up...
  gsl_vector_free(eprob);
  gsl_vector_free(prior);
  gsl_vector_free(post);
  gsl_vector_free(beta);
  gsl_vector_free(mem);
  gsl_vector_free(mem2);
  gsl_vector_free(last_beta);
  gsl_vector_free(last_eprob);
  gsl_matrix_free(DP);
}



int Drift::adapt_range(){
  gsl_vector * cum = gsl_vector_calloc(gridSize+1);
  
  double mn = myEmit->xmax;
  double mx = myEmit->xmin;
  int low=0;
  int redo=0;
  double eps = 1.0e-4;
  for (int s=0; s<nLoci; s++){
    Drift::get_posterior(s);
    for (int l=0; l<sTimes[s]; l++){
      if (mask[s][l] == 0) continue;
      low=0;
      cum->data[0] = 0.0;
      for (int i=1; i<=gridSize;i++){
	cum->data[i]  = cum->data[i-1];
	cum->data[i] += 0.5*myEmit->dx*(gsl_matrix_get(gamma[s],l,i-1) + gsl_matrix_get(gamma[s],l,i));
	if (low==0 && cum->data[i] >= eps){
	  mn = min( mn, myEmit->xgrid[i-1]);
	  low=1;
	}
	if (cum->data[i] >= 1.0-eps){
	  mx = max( mx, myEmit->xgrid[i]);
	  break;	  
	}		
      }
    }
    gsl_matrix_free(gamma[s]);
    gamma[s] = NULL;
  }
  gsl_vector_free(cum);
  if (mx-mn < 10.0){
    double xc = 0.5*(mn+mx);
    mn = max(0.01,xc - 5.0);
    mx = xc + 5.0;
  }
  if (mn > myEmit->xmin || mx < myEmit->xmax){
    redo = 1;
    myEmit->xmin = mn;
    myEmit->xmax = mx;
    if (myEmit->bias == NULL){
      myEmit->ymin = mn;
      myEmit->ymax = mx;
    }
    myEmit->set_grid();
    myEmit->set_EmitProb(Rep);
    Fwd_done=0;
    Bwd_done=0;
  }
  return(redo);
}


