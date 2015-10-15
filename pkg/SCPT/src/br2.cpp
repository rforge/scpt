#include <Rcpp.h>
using namespace Rcpp;
RcppExport SEXP br2_corr(SEXP m, SEXP k, SEXP lb){
	
	int M = as<int>(m);
	int K = as<int>(k);
	NumericMatrix LL(lb);//must be low triangular

	NumericVector W(M);//Weiner process
	NumericVector B(M);//Brownian Bridge
	NumericMatrix MAT(K,M);
	NumericMatrix BCR(K,M);
	NumericVector BR2(K);
	double sqdt = 1/sqrt((double) M - 1);
	double step = 1/((double) M - 1);
	double l = 0, sum = 0, s0 = 0, s1 = 0, s2 = 0;

	unsigned int i,j,t;
	RNGScope scope;

	for(j = 0; j < K; j++){
		W(0) = 0;
		for(i = 1; i < M; i++)
			W(i) = W(i-1)+ sqdt*(::Rf_rnorm(0, 1));

		l = W(M-1);
		for(i = 0; i < M; i++)
			MAT(j,i) = W(i) - i*step*l;
	}

	for(j = 0; j < M; j++){
		for(i = 0; i < K; i++){
			l = 0;
			for(t = 0; t <= i; t++){
				l += LL(i,t)*MAT(t,j);
			}
			BCR(i,j) = l;
		}
	}
	//Numerical Integration

	for(i = 0; i < K; i++){
		sum = 0; s0 = 0; s1 = 0; s2 = 0;
        for (j = 1; j < M-1; j += 2){
           	s0 += BCR(i,j)*BCR(i,j);
           	s1 += BCR(i,j-1)*BCR(i,j-1);
           	s2 += BCR(i,j+1)*BCR(i,j+1);
        }
        sum = step*(s1+4*s0+s2)/3;
        /* If n is even, add the last slice separately */
        if (M % 2 == 0){
            sum += step*(5*BCR(i,M-1)*BCR(i,M-1)+
						 8*BCR(i,M-2)*BCR(i,M-2)-
						   BCR(i,M-3)*BCR(i,M-3))/12;
					//printf("Here!");
        }
		BR2(i) = sum;
	}
	return wrap(BR2);
}

//bb<-cxxfunction(signature(m="integer", k="integer", lb="numeric"),body=src,plugin="Rcpp")
