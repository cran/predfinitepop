domnoplan_totalcomp <-
function(datos_j,rho,ystar,phi,g0_licitacion_sal,N_Stil_domnoplan_j_sim,nSim,colid_D){
#
#   This function simulates Monte Carlo samples from the predictive distribution of the vector \eqn{\mathbold{T}_j} across the \eqn{D} unplanned domains in a given planned domain \eqn{j}.
#   
#	Input:
#		datos_j						-	Data matrix with features and number of individuals in the sample 'S_j'
#		rho							-	(Ux2)-dimensional vector with weights associated with the sample ties 'ystar'
#		ystar						-	Sample ties 'y^*_i' in the sample 'S_j'
#		phi							-	Probability weight associated with "G_{j0}" (the continuous part of Ghat_j)
#		g0_licitacion_sal			-	Object list with the continuous component in 'Ghat'
#		N_Stil_domnoplan_j_sim		-	(1 x D x nSim) matrix with Monte Carlo samples of the predictive distribution of the composition of 'Stil_j' 
#		nSim						-	Number of Monte Carlo simulated replicates of the predictive distribution
#		colid_D						-	D-dimensional matrix array with the columns in ´datos´ that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of 'datos')
#
#	Output:
#	Object list with three entries:
#		T_S_domnoplan 				- 	Composition of 'T_j' for the planned domain 'j' in sample 'S_j'
#		T_Stil_domnoplan_j 			- 	(1xDxnSim)-dimensional array with simulated samples of the convolution for 'T_Stil_j'
#		domnoplan_totalcomp_sim 	- 	(1xDxnSim)-dimensional array with samples from the predictive distribution of the composition of 'Stil_j'
#

#	Number of unplanned domains
D <- length(colid_D)
dim_D_j <- dim(datos_j)

#	Repository: Unplanned domains "T_Stil_domnoplan_d_j"
domnoplan_totalcomp_sim <- array(NaN,c(D,nSim))

#	Repository: Non-sampled total the domain j - "T_Stil_domnoplan_j"
T_Stil_domnoplan_j <- array(NaN,c(1,nSim))
rownames(T_Stil_domnoplan_j) <- c("T_Stil_domnoplan_j")
colnames(T_Stil_domnoplan_j) <- c(1:nSim)

# Computing sample totals
d <- 1
for(d in 1:D){
	#	vector with number of members of the group for each unplanned domains
	T_S_domnoplan <- datos_j[,"n_i"] * datos_j[,"y_i"] * datos_j[,colid_D[d]]
	datos_j <- cbind(datos_j,T_S_domnoplan)
	}

#	Composition of the total across domains on the sample unplanned "S_j"
T_S_domnoplan <- colSums(datos_j[,c((dim_D_j[2]+1):(dim_D_j[2]+D))])

#	Simulations
sim <- 1
d <- 1
for(sim in 1:nSim){
	for(d in 1:D){
		#	Extract the auxiliary size N_Stil_domnoplan1_d
		N_d_aux <- N_Stil_domnoplan_j_sim[d,sim]
		
		if(N_d_aux > 0){
			#	simulating a sample of size "N_d_aux" of "ghat"
			y_d_aux <- ghat_simulacion(rho,ystar,phi,g0_licitacion_sal,N_d_aux)
			
			#	Simulation subtotal "T_Stil_domnoplan1_d"
			domnoplan_totalcomp_sim[d,sim] <- sum(y_d_aux)
		}else{
			#	Sum of the subtotal "T_Stil_domnoplan1_d" zero
			domnoplan_totalcomp_sim[d,sim] <- 0
		}
	}
	T_Stil_domnoplan_j[sim] <- sum(domnoplan_totalcomp_sim[,sim])
}

#	Output
salida <- list(T_S_domnoplan,T_Stil_domnoplan_j,domnoplan_totalcomp_sim)
return(salida)

#
#	--	END of domnoplan_totalcomp.R	--
}
