#include "RcppArmadillo.h"


using namespace Rcpp;
RcppExport SEXP get_l2stat(SEXP DataArray, SEXP Dims, SEXP w){
    NumericVector vecArray(DataArray); // vector with data
    IntegerVector arrayDims(Dims);
    NumericVector weights(w);
    
    double tt, temp, sum, s0, s1, s2;
    int T =  arrayDims[0];
    int K =  arrayDims[1]; //number of spatial locations
    int N =  arrayDims[2]; //number of years
    double step = 1/((double) T - 1);
    arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
    
    NumericMatrix mt_hat(T,K);//total mean matrix
    NumericMatrix mt_hss(T,K);//partial mean matrix
    NumericVector r2_stat(N);//test statistics
    //NumericVector tst(N);
    //tst(0) = 0;
    
    
    unsigned int i,j,k,l; //counters
    
    for(k = 0; k < K; k++){
        for(i = 0; i < T; i++){
            mt_hat(i,k) = 0;
            mt_hss(i,k) = 0;
            for(j = 0; j < N; j++){
                mt_hat(i,k) += cubeArray(i,k,j)/((double) N);
            }
        }
    }
    
    for(j = 0; j < N; j++){
        temp = 0;
        for(k = 0; k < K; k++){
            for(i = 0; i < T; i++){
                mt_hss(i,k) += cubeArray(i,k,j);
            }
            tt = 0;
            r2_stat(j) = 0;
            //integration part
            sum = 0; s0 = 0; s1 = 0; s2 = 0;
            for (l = 1; l < T-1; l += 2){
                s0 += (mt_hss(l,k)-(j+1)*mt_hat(l,k))*(mt_hss(l,k)-(j+1)*mt_hat(l,k));
                s1 += (mt_hss(l-1,k)-(j+1)*mt_hat(l-1,k))*(mt_hss(l-1,k)-(j+1)*mt_hat(l-1,k));
                s2 += (mt_hss(l+1,k)-(j+1)*mt_hat(l+1,k))*(mt_hss(l+1,k)-(j+1)*mt_hat(l+1,k));
            }
            sum = step*(s1+4*s0+s2)/3;
            /* If n is even, add the last slice separately */
            if (T % 2 == 0){
                sum += step*(5*(mt_hss(T-1,k)-(j+1)*mt_hat(T-1,k))*(mt_hss(T-1,k)-(j+1)*mt_hat(T-1,k))+
                             8*(mt_hss(T-2,k)-(j+1)*mt_hat(T-2,k))*(mt_hss(T-2,k)-(j+1)*mt_hat(T-2,k))-
                             (mt_hss(T-3,k)-(j+1)*mt_hat(T-3,k))*(mt_hss(T-3,k)-(j+1)*mt_hat(T-3,k)))/12;
                //printf("Here!");
            }
            temp += weights(k)*sum;
        }
        r2_stat(j) = temp/((double) N);
    }
    return wrap(r2_stat);
}
//bb<-cxxfunction(signature(m="integer", k="integer", lb="numeric"),body=src,plugin="Rcpp")
