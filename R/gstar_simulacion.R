gstar_simulacion <-
function(rho,ystar,nSim){
#
#	This function generates Monte Carlo samples from the discrete component of \eqn{\hat{G}_j} for the Dirichlet process prior in the planned domain \eqn{j}.
#	
#	Input:
#		rho						-	Ux2 dimensional vector with weights associated with 'ystar'
#		ystar					-	Sample ties 'y^*_i' in 'S_j'
#		phi						-	Probability weight associated with 'G_{j0}' (the continuous part of 'Ghat_j')
#
#	Output:
#		gstar_simulacion_sal	- (nSim x 1) dimensional vector with simulated data
#

#	Normalising the weights "rho"
rho_norm <- rho / sum(rho) 

#-----------------------------------
#   Verification:
unicstar <- cbind(ystar,rho_norm)

#-----------------------------------
# Generating  the final sample discrete component (with replacement)
gstar_simulacion_sal <- sample(unicstar[,"ystar"], nSim, replace = TRUE, prob = unicstar[,"rho"])

#	output
salida <- gstar_simulacion_sal
return(salida)

#
#	--	END of gstar_simulacion.R --
}
