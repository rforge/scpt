#' @name SCPT
#' @title Caries out a change-point test for sequences of
#' functional data observed at different spatial locations.  Three tests are carried out
#' and three pvalues are returned.
#'
#' @param DATA3 Data for testing.  This should be an 3D array with the
#' first coordinate indicating the temporal location on a curve, the
#' second the spatial location, and the third the repitition number/time series iteration.
#' @param loc Two column matrix of spatial locations.
#' @param varprop Proportion of temporal variance explained in first two
#' tests (third test uses all explained variance)
#' @param ncores Number of cores to be used for monte carlo sampling.
#' @param reps_per_core Number of monte carlo samples per core.
#'
#' @details Carries out three different change point tests as outlined in Gromenko, Kokoszka,
#' and Reimherr (2015).  P-values are generated using Monte-Carlo.  The total MC sample size
#' is given by ncores*reps_per_core.  Be warned that the code will not try to use more cores
#' than are available, so users should choose the ncores value appropriately.  By default, the
#' code will use all available cores.  This will cause a computer to be slugish in other tasks
#' and so the user should be careful when using the default.
#'
#' Data3[i,j,k] indicates an the ith observation at the jth location and on the kth curve.  It is
#' important that this order is accurate.  The method assumes that the spatio-temporal
#' covariance is separable.  Methods such as Constantinou, Kokoszka, and Reimherr (2015) can be
#' used to check this assumption.  This allows the spatial and temporal covariances to be
#' estimated separately.  The spatial covariance is assumed to be stationary and isotropic.  A
#' tapered bsplines estimator is then used to estiamte it nonparametrically.  No assumptions are
#' made on the form of the temporal covariance function.  It is estimated by pooling appropriately
#' across space.
#'
#' The testing procedures are all based on integrated CUSUM test statistics.  The primary difference
#' between the three is in how they handle the temporal covariance.  The first statistic uses
#' FPCA to reduce the dimension of the data and normalizes by the corresponding eigenvalues.  The
#' second does the same, but does not normalize by the eigenvalue.  The third procedure does not use
#' FPCA and does not normalize by the eigenvalues.  It can be thought of as the same as procedure 2,
#' but with all dimensions included.  The first procedure should be used with caution as it can
#' be unstable for smaller sample sizes/larger numbers of FPCs.
#'
#' @return Returns vector of three values for three change-point tests.  The first normalizes each FPC, while the second two do not.  The only difference between the second and third is that the third uses all explained variance.
#'
#' @example 
#' mean(1:4)
#'
#' @import inline
#' @import Rcpp
#' @import RcppArmadillo
#' @import splines
#' @import parallel
#' @useDynLib scpt
#' @export

spatial_change_test<-function(DATA3,loc,varprop=.85,ncores='all',reps_per_core=100){
	#require(parallel); require(RcppArmadillo); require(splines); require(inline)
	#source("Simpson.R"); source("Fast_ST.R"); source("Inline.R"); source("Inline_L2.R"); source("Auxiliary_Func.R")

	L_per_curve<-dim(DATA3)[1]
	N<-dim(DATA3)[2] # number of spatial locations
	NY<-dim(DATA3)[3] # number of curves per location
	L<-L_per_curve*NY # total temporal observations
	DATA.DIFF3<-DATA3[,,2:NY]-DATA3[,,1:(NY-1)]

	#c.hat.pr <- get_svar(as.matrix(loc*pi/180), DATA.DIFF3, dim(DATA.DIFF3))
	c.hat.pr <- .Call('get_svar',PACKAGE='scpt',as.matrix(loc*pi/180), DATA.DIFF3, dim(DATA.DIFF3))
	dd<-c.hat.pr[[1]]
	H.hat<-c.hat.pr[[3]]
	dd.v<-as.vector(dd)
	H.hat.v<-as.vector(H.hat)
	#removing zeroes
	H.hat.vn<-H.hat.v[dd.v!=0]
	dd.vn<-dd.v[dd.v!=0]

	#noparametric BSpline estimation
	fm1 <- lm(H.hat.vn ~ bs(dd.vn, df = 9))
	#truncation
	#a<-0.25
	#b<-0.27
	a<-quantile(dd.vn,.75)
	b<-quantile(dd.vn,.95)
	slope<-predict(fm1, data.frame(dd.vn = a))/(b-a)

	H.final<-matrix(1,N,N)
	sp.cov<-predict(fm1, data.frame(dd.vn = dd.vn))
	for(i in 1:length(dd.vn)){
		H.final[dd.vn[i]==dd]<-sp.cov[i]
		if(dd.vn[i] > a & dd.vn[i] < b){
			H.final[dd.vn[i]==dd] <- slope*(b-a)-(dd.vn[i]-a)*slope;
		}
		if(dd.vn[i] > b){
			H.final[dd.vn[i]==dd] <- 0;
		}
	}

	H.vec<-as.vector(H.final)
	H.vec.n<-H.vec[dd!=0]
	dd.n<-dd[dd!=0]
	cov.dt<-cbind(dd.n, H.vec.n)
	cov.dt<-cov.dt[order(cov.dt[,1]),]

	#caclulating weights
	c.hm<-matrix(0,N,N)
	diag(c.hm)<-c.hat.pr[[2]]
	rho.hat.f<-c.hm%*%H.final%*%c.hm
	one<-matrix(1, ncol=1, nrow=N)
	w<-solve(rho.hat.f^2)%*%one
	w<-w/sum(w[,1])

	#make_cov_fig()
	#temporal covariance
	#tvar<-get_tvar(w[,1], c.hat.pr[[2]], DATA.DIFF3, dim(DATA.DIFF3))#may take some time
	tvar<-.Call('get_tvar',PACKAGE='scpt',w[,1], c.hat.pr[[2]], DATA.DIFF3, dim(DATA.DIFF3))
	FPC<-eigen(tvar[[1]], symmetric=T, only.values = FALSE, EISPACK = FALSE)
	val<-FPC$values[FPC$values>0]
	val<-val/sum(val)

	cvar<-val[1]
	for(i in 2:length(val))
		cvar<-c(cvar, cvar[i-1]+val[i])

	#NFPC<-length(val[cvar<.78])
	NFPC<-length(val)
	v.hat<-c()
	for(i in 1:NFPC)
		v.hat<-cbind(v.hat, FPC$vectors[,i]/sqrt(simpson_r(FPC$vectors[,i]^2)))

	rho_e<-eigen(rho.hat.f)
	rho.hat.f.clip<-(rho_e$vectors[,rho_e$values>0])%*%diag(rho_e$values[rho_e$values>0])%*%t(rho_e$vectors[,rho_e$values>0])
	suppressWarnings(Lb<-matrix(t(chol(rho.hat.f.clip,pivot=TRUE)),nrow=N))
	cores_all<-detectCores()
	if(ncores=='all'||ncores>=cores_all){
		XXX<-vector("list",cores_all)
		for(i in 1:cores_all)
			XXX[[i]]<-reps_per_core
		n.cores<-cores_all
	}else{XXX<-vector("list",ncores)
		for(i in 1:ncores)
			XXX[[i]]<-reps_per_core
		n.cores<-ncores
		}
	#XXX<-vector("list",n.cores)
	#for(i in 1:n.cores)
	#	XXX[[i]]<-reps_per_core
	#	#XXX[[i]]<-250

	#tmp<-gn_st(20,Lb,val,w,L_per_curve)
	#return(tmp)

	#------------------------------------
	#generating MC distribution
	#------------------------------------
	#!!! will use all available cores !!!
	#system.time(
	#result<-mcmapply(XXX, FUN=gn_st, mc.cores=n.cores,MoreArgs = list(Lb=Lb,val=val,w=w))
	result<-mclapply(XXX, FUN=gn_st, mc.cores=n.cores,Lb=Lb,val=val,w=w,L_per_curve = L_per_curve)
	#)

	#gathering data
	SL1<-c()
	SL2<-c()

	for(i in 1:n.cores){
		SL1<-rbind(SL1, result[[i]]$L1)
		SL2<-rbind(SL2, result[[i]]$L2)
	}

	lst<-length(SL2[,1])

	nnn<-NFPC
	#aaa<-get_stat(DATA3, dim(DATA3), as.matrix(v.hat[,1:nnn]), w[,1], val[1:nnn], nnn)
	aaa<-.Call('get_stat',PACKAGE='scpt',DATA3, dim(DATA3), as.matrix(v.hat[,1:nnn]), w[,1], val[1:nnn], nnn)
	r1<-aaa[[1]][-1,]
	r2<-aaa[[2]][-1,]

	L1<-rep(0,nnn)
	L2<-rep(0,nnn)
	pv1<-rep(0,nnn)
	pv2<-rep(0,nnn)
	for(i in 1:nnn){
		L1[i]<-simpson_r(r1[i,])
		L2[i]<-simpson_r(r2[i,])
		pv1[i]<-sum((rep(1, lst))[SL1[,i] >= L1[i]])/lst;
		pv2[i]<-sum((rep(1, lst))[SL2[,i] >= L2[i]])/lst;
	}
	spt_for_ret<-which(cumsum(val)>varprop*sum(val))[1]
	for_return<-data.frame(pv1[spt_for_ret],pv2[spt_for_ret],pv2[nnn])
	colnames(for_return)<-c("Pvalue_Test1", "Pvalue_Test2", "Pvalue_Test2_all")

	return(for_return)
}
