g0_simulacion <-
function(g0_licitacion_sal,nSim){
#
#	This function computes Monte Carlo samples from a given continuous \eqn{G_{j0}}.
#	
#	Input:
#		g0_licitacion_sal	-	Object list with three elements (produced with 'g0_licitacion'):
#								a)	String object for 'distribution'
#									  i) Gamma
#									 ii) Weibull
#									iii) Log-normal
#									iv)  Inverse-Gaussian 
#
#								b) theta			-	Vector of parameters associated with 'distribution'
#								c) mu				-	Expected value of 'distribution'		
#		nSim				-	Number of Monte Carlo simulations
#
#	Output:
#		g0_simulacion_sal	-	(nSim x 1) dimensional array with simulated data.  
#

#	Simulations case-particularl
if(g0_licitacion_sal[[1]]=='Lognormal'){
	g0_simulacion_sal <- as.matrix(rlnorm(nSim, meanlog =  g0_licitacion_sal[[2]][1], sdlog =  g0_licitacion_sal[[2]][2]))
	colnames(g0_simulacion_sal) <- c("ind_sim")
}else if(g0_licitacion_sal[[1]]=='Gamma'){
	g0_simulacion_sal <- as.matrix(rgamma(nSim, shape = g0_licitacion_sal[[2]][1], scale = g0_licitacion_sal[[2]][2]))
	colnames(g0_simulacion_sal) <- c("ind_sim")
}else if(g0_licitacion_sal[[1]]=='Weibull'){
	g0_simulacion_sal <- as.matrix(rweibull(nSim, shape = g0_licitacion_sal[[2]][1], scale = g0_licitacion_sal[[2]][2]))
	colnames(g0_simulacion_sal) <- c("ind_sim")
}else if (g0_licitacion_sal[[1]]=='Inverse-Gaussian'){
	g0_simulacion_sal <- as.matrix(rinvgauss(nSim, g0_licitacion_sal[[2]][1], g0_licitacion_sal[[2]][2]))
	colnames(g0_simulacion_sal) <- c("ind_sim")
}

#	Output
salida <- g0_simulacion_sal
return(salida)

#
#	--	END of g0_simulacion.R --
}
