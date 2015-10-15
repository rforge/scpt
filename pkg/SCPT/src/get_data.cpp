#include <Rcpp.h>
using namespace Rcpp;
RcppExport SEXP get_data(SEXP lb, SEXP values, SEXP fpc, SEXP n, SEXP nfpc){
	
    NumericMatrix FPC(fpc);
    NumericVector sqval(values);
    NumericMatrix Lb(lb);//low triangular matrix;
    
    int N = as<int>(n); //number of spatial locations
    int NFPC = as<int>(nfpc); //number of components
    int L = 31;//number of days
    NumericMatrix Data(L,N);//matrix with data
    RNGScope scope;
    
    NumericMatrix scores(N, NFPC);
    NumericVector sc_tem(N);
    
    unsigned int i,j,k,l;
    
    for(i = 0; i < NFPC; i++){
        for(k = 0; k < N; k++){
            sc_tem(k) = (::Rf_rnorm(0, 1));
            scores(k,i) = 0;
        }
        //correlating
        for(k = 0; k < N; k++){
            for(l = 0; l <= k; l++){
                scores(k,i) += Lb(k,l)*sc_tem(l);
            }
        }
    }
    
    for(k = 0; k < N; k++){
        for(j = 0; j < L; j++){
            Data(j,k) = 0;
        }
    }
    
    
    for(k = 0; k < N; k++){
        for(i = 0; i < NFPC; i++){
            for(j = 0; j < L; j++){
                Data(j,k) += sqval(i)*scores(k,i)*FPC(j,i);
            }
        }
    }
    return(wrap(Data));
}

//bb<-cxxfunction(signature(m="integer", k="integer", lb="numeric"),body=src,plugin="Rcpp")
