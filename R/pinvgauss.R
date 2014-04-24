pinvgauss <-
function(q, m, s){
#	Inverse Gaussian distribution - Probability 
#
#	Package: 	rmutil
#	Version: 	1.0
#	Title: 		Utilities for Nonlinear Regression and Repeated Measurements Models
#	Author: 	Jim Lindsey <jlindsey@luc.ac.be>
if(any(q<=0)){stop("q must contain positive values")}
if(any(m<=0)){stop("m must be positive")}
if(any(s<=0)){stop("s must be positive")}
t <- q/m
v <- sqrt(q*s)
out <- pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)
return(out)
}
