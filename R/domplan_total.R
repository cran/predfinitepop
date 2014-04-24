domplan_total <-
function(datos_j,rho,ystar,phi,g0_licitacion_sal,N_Stil_j,nSim){
#
#	Simulates Monte Carlo samples of size '\eqn{nSim}' of the final distribution for the total for the given planned domain \eqn{\mathcal{P}_j}.
#
#	Input:
#		datos_j						-	Data matrix with features and number of individuals in the sample 'S_j'
#		rho							-	(Ux2)-dimensional vector with weights associated with the sample ties 'ystar'
#		ystar						-	Sample ties 'y^*_i' in the sample 'S_j'
#		phi							-	Probability weight associated with 'G_{j0}' (the continuous part of Ghat_j)
#		g0_licitacion_sal			-	Object list with the continuous component in 'Ghat'
#		N_Stil_sal					-	Composition of individuals out of 'S_j' 
#		nSim						-	Number of Monte Carlo simulated replicates of the predictive distribution
#
#	Output:
#	Object list with two entries:
#		T_S_j		 				- 	Composition of 'T_j' across planned domains in the sample 'S_j'
#		domplan_total_sim			-	(1 x nSim) matrix with samples from the predictive distribution of the total 'T_Stil_j'
#

#	Repository:	Total non-sampling domain j -- "T_Stil_domnoplan_j"
T_Stil_domnplan_j <- array(NaN,c(1,nSim))
rownames(T_Stil_domnplan_j) <- c("T_Stil_domnplan_j")
colnames(T_Stil_domnplan_j) <- c(1:nSim)

# 	Computing the total "T_S_j" compositional in "S_j" (within the sample)
#	Vector with the totals of each group
T_S_j_gpo <- datos_j[,"n_i"] * datos_j[,"y_i"]
datos_j <- cbind(datos_j,T_S_j_gpo)

#	Composition of the total between domains on the sample unplanned "S_j"
T_S_j <- sum(datos_j[,"T_S_j_gpo"])

#	Simulations
sim <- 1
for(sim in 1:nSim){
	#	Computing the sample of size "N_Stil_j" of "ghat"
	y_i_aux <- ghat_simulacion(rho,ystar,phi,g0_licitacion_sal,N_Stil_j)
		
	#	Simulation subtotal "T_Stil_j"
	T_Stil_domnplan_j[sim] <- sum(y_i_aux)
	}

#	Output
salida <- list(T_S_j,T_Stil_domnplan_j)
return(salida)

#
#	--	END of  domplan_total.R	--
}
