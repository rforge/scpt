#' @title MC Sampling
#' @name gn_st
#' @details  This function is used to carry out Monte Carlo sampling to compute p-values
#' @param par Number of MC samples
#' @param Lb Cholesky decomposition of spatial covariance
#' @param val Temporal eigenvalues
#' @param w Spatial weights, test statistic pools across space using these weights.
#' @param L_per_curve Number of observations per curve


gn_st<-function(par=250,Lb,val,w,L_per_curve){
	#dyn.load("br2.so")
	K<-length(Lb[,1])
	N<-L_per_curve
	np<-length(val)#nfpc to include
	stat_l1.lvec<-c()
	stat_l2.lvec<-c()
	for(rrr in 1:par){
		BB<-c()
		for(i in 1:np)
			BB<-cbind(BB, t(as.matrix(.Call("br2_corr", N, K, Lb)))%*%w)
		stat1<-array(0,np)
		stat2<-array(0,np)
		stat1[1]<-BB[1,1]
		stat2[1]<-val[1]*BB[1,1]
		for(i in 2:np){
			stat1[i]<-stat1[i-1]+BB[1,i]
			stat2[i]<-stat2[i-1]+val[i]*BB[1,i]
		}
		stat_l1.lvec<-rbind(stat_l1.lvec, stat1)
		stat_l2.lvec<-rbind(stat_l2.lvec, stat2)
	}
	return(list(L1=stat_l1.lvec, L2=stat_l2.lvec))
}

