dinvgauss <-
function(y, m, s, logcond=FALSE){
#	Inverse Gaussian distribution - Density 
#
#	Package: 	rmutil
#	Version: 	1.0
#	Title: 		Utilities for Nonlinear Regression and Repeated Measurements Models
#	Author: 	Jim Lindsey <jlindsey@luc.ac.be>
if(any(y<=0)){stop("y must contain positive values")}
if(any(m<=0)){stop("m must be positive")}
if(any(s<=0)){stop("s must be positive")}
tmp <- -(y-m)^2/(2*y*s*m^2)-(log(2*pi*s)+3*log(y))/2
if(!logcond){tmp <- exp(tmp)}
return(tmp)
}
