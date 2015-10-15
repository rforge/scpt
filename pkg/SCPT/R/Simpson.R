#' @title Simpsons Rule
#' @name Simpson Rule
#'
#' @details  Numeric integration using simpsons rule. Uses the simpson dynamic library

# simpson function which perfoms numerical integration
# of a time series using simpson's formula
# \int_{x_{i}}^{x_{i+2}}f(x)dx ~ \Delta x /3 (f(x_i)+4f(x_{i+1})+f(x_{i+1}) )

# important! for correct integration length
# must be odd
# parameters: series - initial time seies
#
#


#dyn.load("simpson.so")

simpson_r<-function(series){
	len<-length(series);
	sum<-0;
	sum<-.C("simpson", as.double(series), as.double(len-1), as.double(sum))[[3]];
	return(sum);
}

#  Example:
#  area<-simpson(c(data[,2]*data[,1]))
#

#
# Method of Trapezoids
#temp=(data[1,1]+ data[length(data[,1]),1])/2
#for (i in 2:(length(data[,1])-1)){
#temp=temp + data[i,1]
#}

# method of Rectangle
#temp=0
#for (i in 1:(length(data[,1])-1)){
#temp=temp + data[i,1]
#}


