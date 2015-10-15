#include "RcppArmadillo.h"


using namespace Rcpp;
RcppExport SEXP get_svar_CL(SEXP locArray, SEXP DataArray, SEXP Dims){
    NumericMatrix loc(locArray);
    //int nrows = loc.nrow();
    //int ncolumns = loc.ncol();
    
    //double temp = 0;
    
    NumericVector vecArray(DataArray);
    IntegerVector arrayDims(Dims);
    int T =  arrayDims[0];
    int K =  arrayDims[1]; //number of spatial locations
    int N =  arrayDims[2]; //number of years
    arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
    NumericVector c_hat(K);//covariance at a single location;
    NumericMatrix dist(K, K);
    NumericMatrix rho(K, K);
    double temp;
    unsigned int i,j,l,n;//indicators
    
    //printf("%f %d %d %d", cubeArray(3-1,1-1,10-1), T, K, N);
    
    //for numerical integration:
    double h = 1/((double) T-1);
    double s0,s1,s2, sum;
    
    for(i = 0; i < K; i++){
        temp = 0;
        for (j = i; j < K; j++) {
            temp = 2*sqrt((sin((loc(i,0)-loc(j,0))/2))*(sin((loc(i,0)-loc(j,0))/2))+
                          cos(loc(i,0))*cos(loc(j,0))*(sin((loc(i,1)-loc(j,1))/2))*(sin((loc(i,1)-loc(j,1))/2)));
            dist(i,j)=temp; dist(j,i)=temp;
        }
        
        temp = 0;
        for(j = 0; j < N; j++){
            sum = 0; s0 = 0; s1 = 0; s2 = 0;
            for (l = 1; l < T-1; l += 2){
                s0 += cubeArray(l,i,j)*cubeArray(l,i,j);
                s1 += cubeArray(l-1,i,j)*cubeArray(l-1,i,j);
                s2 += cubeArray(l+1,i,j)*cubeArray(l+1,i,j);
            }
            sum = h*(s1+4*s0+s2)/3;
            /* If n is even, add the last slice separately */
            if (T % 2 == 0){
                sum += h*(5*cubeArray(T-1,i,j)*cubeArray(T-1,i,j)+
                          8*cubeArray(T-2,i,j)*cubeArray(T-2,i,j)-
                          cubeArray(T-3,i,j)*cubeArray(T-3,i,j))/12;
                //printf("Here!");
            }
            
            temp += sum/((double) N - 1);
        }
        c_hat(i) = sqrt(temp);
    }
    
    for(i = 0; i < K; i++){
        for(j = i; j < K; j++){
            temp = 0;
            for(n = 0; n < N; n++){
                sum = 0; s0 = 0; s1 = 0; s2 = 0;
                for (l = 1; l < T-1; l += 2){
                    s0 += cubeArray(l,i,n)*cubeArray(l,j,n);
                    s1 += cubeArray(l-1,i,n)*cubeArray(l-1,j,n);
                    s2 += cubeArray(l+1,i,n)*cubeArray(l+1,j,n);
                }
                sum = h*(s1+4*s0+s2)/3;
                /* If n is even, add the last slice separately */
                if (T % 2 == 0){
                    sum += h*(5*cubeArray(T-1,i,n)*cubeArray(T-1,j,n)+
                              8*cubeArray(T-2,i,n)*cubeArray(T-2,j,n)-
                              cubeArray(T-3,i,n)*cubeArray(T-3,j,n))/12;
                    //printf("Here!");
                }
                
                temp += sum/((double) N-1)/(c_hat(i)*c_hat(j));
            }
            rho(i,j) = temp;
            rho(j,i) = temp;
        }
    }
    
    return wrap(List::create(dist, c_hat, rho));
}
//bb<-cxxfunction(signature(m="integer", k="integer", lb="numeric"),body=src,plugin="Rcpp")
