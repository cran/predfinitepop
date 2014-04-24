ghat_simulacion <-
function(rho,ystar,phi,g0_licitacion_sal,N_Stil_j){
#
#	This function simulates Monte Carlo samples from 'Ghat_j'.
#	
#	Input:
#		rho						-	Ux2 dimensional vector with weights associated with 'ystar'
#		ystar					-	Sample ties 'y^*_i' in 'S_j'
#		phi						-	Probability weight associated with 'G_{j0}' (the continuous part of 'Ghat_j')
#		g0_licitacion_sal		-	Object list with three elements (produced with 'g0_licitacion'):
#									a)	String object for 'distribution'	
#											i) Gamma
#										   ii) Weibull
#										  iii) Log-normal
#										   iv) Inverse-Gaussian 
#									b) theta			-	Parameters associated with 'distribution'
#									c) mu				-	Expected value of 'distribution'
#		N_Stil_j				-	Number of samples to simulate (number of individuals in 'Stil_j') out of the sample
#
#	Output:
#		ghat_simulacion_sal	-	(N_Stil_j x 1) dimensional vector with simulated data
#

#	Repository
ghat_simulacion_sal <- matrix(NaN,N_Stil_j,1)
colnames(ghat_simulacion_sal) <- c("ghat_simulacion")
#	Simulation
sim <- 1
if(N_Stil_j > 0){
	for(sim in 1:N_Stil_j){
		#	Simulated component "g_0"
		#	1) Discrete "ystar"
		#	0) Continuous "g_0"
		ghat_comp_sim <- rbinom(1,1,sum(rho))
		
		#	Simulation "Y_jl"
		if( ghat_comp_sim == 1){
			#	Simulating the discrete component
			ghat_simulacion_sal[sim] <- gstar_simulacion(rho,ystar,1)
			}else{
			#	Simulating the discrete component
			ghat_simulacion_sal[sim] <- g0_simulacion(g0_licitacion_sal,1)
			     }
		                   }
}else{warning("The sample size is incorrect.")}
#	Output
salida <- ghat_simulacion_sal
return(salida)
#	--	END of ghat_simulacion.R --
}
