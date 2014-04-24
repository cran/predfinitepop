domnoplan_composicion <-
function(datos_j,N_j,N_S_j,N_Stil_j,colid_D,alpha_D,nSim){
#
#	This function simulates Monte Carlo samples from the predictive distribution of the vector \eqn{\mathbold{N}_j} across the \eqn{D} unplanned domains in a given planned domain \eqn{j}.
#
#	Input:
#		datos_j						-	Data matrix with features and number of individuals in the sample 'S_j'
#		N_j							-	Number of individuals in 'P_j' (jth planned domain)
#		N_S_j						-	Number of individuals in the sample 'S_j'
#		N_Stil_j					-	Number of individuals out the sample 'S_j'
#		colid_D						-	D-dimensional matrix array with the columns in ´datos´ that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of 'datos')
#		alpha_D						-	D-dimensional array with positive entries for the parameters of the multinomial-Dirichlet component for the composition across unplanned domains (NOTE: this one makes reference to a single planned domain)
#		nSim						-	Number of Monte Carlo simulated replicates of the predictive distribution
#
#	Output:
#		N_S_domnoplan				-	Composition of the number of individuals in sample 'S_j' ('N_S_j') across the 'D' unplanned domains
#		domnoplan_composicion_sim	-	(1 x D x nSim) matrix with Monte Carlo samples of the predictive distribution of the composition of 'Stil_j' 
#

#	Size of unplanned domains
D <- length(colid_D)
dim_D_j <- dim(datos_j)

#	Sample Repository
domnoplan_composicion_sim <- array(0,c(D,nSim))
colnames(domnoplan_composicion_sim) <- c(1:nSim)

# 	Computing sample counts
d <- 1
for(d in 1:D){
	#	Creating the table with members of the group number for each unplanned domain
	N_S_domnoplan <- datos_j[,"n_i"] * datos_j[,colid_D[d]]
	datos_j <- cbind(datos_j,N_S_domnoplan)
	}

#	Composition of unplanned domains in the sample "S_j"
N_S_domnoplan <- colSums(datos_j[,c((dim_D_j[2]+1):(dim_D_j[2]+D))])

if(N_Stil_j >0){
	#	Parameters of the Dirichlet distribution for the proportions of "domnoplan"
	alpha_D_new <- alpha_D + N_S_domnoplan

	#	Predictive simulation of multinomial-Dirichlet model
	sim <- 1
	for(sim in 1:nSim){
		#	Sample of the composition "p_domplan_j"
		p_domnoplan_j <- randDirichlet(alpha_D_new,1)
		
		#	Sample of the composition "N_Stil_domplan_j"
		N_Stil_domnoplan_j <- t(rmultinom(1, N_Stil_j, p_domnoplan_j))
		
		#	Repository
		domnoplan_composicion_sim[, sim] <- N_Stil_domnoplan_j
		}
}

#	Output
salida <- list(N_S_domnoplan,domnoplan_composicion_sim)
return(salida)

#
#	--	End of domnoplan_composicion.R	--
}
