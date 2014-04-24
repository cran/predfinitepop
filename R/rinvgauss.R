rinvgauss <-
function(n=1, m, s){
#	Inverse Gaussian distribution - Simulation
#
#	Package: 	rmutil
#	Version: 	1.0
#	Title: 		Utilities for Nonlinear Regression and Repeated Measurements Models
#	Author: 	Jim Lindsey <jlindsey@luc.ac.be>
temp <- qinvgauss(runif(n),m=m,s=s)
return(temp)
}
