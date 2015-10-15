#include "RcppArmadillo.h"


using namespace Rcpp;
RcppExport SEXP get_tvar_CL(SEXP vw, SEXP ch, SEXP DataArray, SEXP Dims){
    NumericVector w(vw);
    NumericVector chv(ch);
    
    NumericVector vecArray(DataArray);
    IntegerVector arrayDims(Dims);
    int T =  arrayDims[0];
    int K =  arrayDims[1]; //number of spatial locations
    int N =  arrayDims[2]; //number of years
    arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
    
    NumericMatrix CT_hat(T, T);
    NumericMatrix CT_sim(T, T);
    
    unsigned int i,j,k,n;
    double t1,t2;
    
    for(i = 0; i < T; i++){
        for(j = i; j < T; j++){
            CT_hat(i,j) = 0; CT_hat(j,i) = 0;
            CT_sim(i,j) = 0; CT_sim(j,i) = 0;
            for(k = 0; k < K; k++){
                t1 = 0; t2 = 0;
                for(n = 0; n < N; n++){
                    t1 += cubeArray(i,k,n)*cubeArray(j,k,n)/((double) N-1);
                }
                CT_hat(i,j) += w(k)*t1/(chv(k)*chv(k));
                CT_sim(i,j) += t1/(chv(k)*chv(k))/((double) K);
            }
            CT_hat(j,i)=CT_hat(i,j);
            CT_sim(j,i)=CT_sim(i,j);
        }
    }
    
    return wrap(List::create(CT_hat, CT_sim));
}
//bb<-cxxfunction(signature(m="integer", k="integer", lb="numeric"),body=src,plugin="Rcpp")
