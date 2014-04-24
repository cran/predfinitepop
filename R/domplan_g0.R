domplan_g0 <-
function(datos,domplan_N,alphaDP,nSim,g0_licitacion_sal){
#            
#	Generates Monte Carlo samples of the predictive distribution of totals of a finite population segmented in planned domains, using a predefined set of \eqn{(G_{j0})_{j=1}^{J}}.
#
#	Input:
#		datos		        -	(Mxp)-dimensional array with positive entries for {S}_j
#		domplan_N			-	Matrix array with counts of individuals in each planned domain
#		alphaDP				-   J-dimensional array with positive entries for the parameters of the Dirichlet process for F_j (with J being the number of planned domains)
#		nSim				-   Number of Monte Carlo simulated replicates of the predictive distribution
#   g0_licitacion_sal  		-   (Jx1) object list, each entry is another object list itself associated with each G_{j0} for the J planned domains. The first element for arch G_{j0} should be the name of the chosen distribution(see details below for alternatives), the second element should be a vector object with the parameters associated with distribution, and the third element should be its associated expectation.
#
#	Details:
#			- datos : Represents the data sample of the target population, unplanned domains labelled.
#  					    It should contain the following columns:
#							"domplan" - Categories fot planned domains.
#							"y_i"     - Actual individual measurements/outcomes (for the moment, they must be positive) for the group of observation.
#							"n_i"     - Number of individuals in the group (if the unit of observation in the sample are individuals, then "n_i" must be equal to 1)
#			
#			- domplan_N : Represents counts (or reference population) of the target population, divided by the planned domains.
#						Tagged data must be labelled by domains planned. It should contain the following columns:
#							"domplan" - Categories for planned domains.
#							"N_j"     - Number of individuals in each planned domain.
#
#			- g0_licitacion_sal: Chose one and only one of the distribution:
#							              i)	Gamma
#						 			     ii)	Weibull
#									    iii)	Lognormal
#								         iv)	Inverse-Gaussian 			
#
#						Parameterization for distribution:
#
#						 i) Gamma    with parameters theta = c(alpha>0 , beta>0) and density function
#							  
#								f(x)= x^{alpha-1} exp{-x/beta}
#								
#								where alpha is the shape parameter, and beta is the scale parameter.									
#
#						ii)	Weibull  with parameters theta = c(alpha>0 , beta>0) and density function
#							   
#								f(x)=(x/beta)^{alpha-1} exp{-(x/beta)^alpha}
#									
#								where alpha is the shape parameter, and beta is the scale parameter.
#
#					   iii)	Lognormal with parameters theta = c(alpha>0 , beta>0) and density function
#							     
#								f(x)= exp{-(log(x)-alpha)^2/(2*beta^2)}
#								
#								where alpha is the mean, and beta is the standard deviation of the logarithm.
#
#					    iv)	Inverse-Gaussian  with parameters theta = c(alpha>0 , beta>0) and density function
#							  
#								f(x)= .....
#									
#								where alpha is the shape parameter, and beta is the scale parameter.
#
#
#	Output:
#		total_domplan_sim	- "(J x 3 x nSim)" Matrix array dimension with the predictions of the planned domains
#									Column 1 - Indicator of the planned domains
#									Column 2 - T_j (totals of the planned domains)
#									Column 3 - N_j (composition of the planned domains)
#

#      Definition consulted domains			
cualdomplan <- as.matrix(domplan_N[,"domplan"])  
J <- length(cualdomplan) 
#		Repository totals in "unplanned domains" with three categories
total_domplan_sim <- array(NaN,c(J,3,nSim))

#----------------------------	
#	Validation
#----------------------------	

#	A.	Parameter alphaDP
if(any(alphaDP <= (0*alphaDP))){
	stop("Error in the specification of 'alphaDP'!!!")}
	
#	B.	Parameter nSim
if(nSim <= 0){
	stop("Error in the specification of 'nSim'!!!")}

#	-----------------------------------------------
#	Scanning planned domains
#	-----------------------------------------------
j <- 1
for(j in 1:J){

    #	Extraction of indexes in "data" (ie sample S_j)
	P_j <- which(datos[,"domplan"]==cualdomplan[j])
	#	Extraction of the data in "data"
	datos_j <- datos[P_j, ]

	#	Population size of \mathcal{P}_{j}
	P_j <- which(domplan_N[,"domplan"]==cualdomplan[j])
	#	Extraction of the data in "data"
	N_j <- domplan_N[P_j,"N_j"]

	#	Name the planned domain grooming and constant data impute "n_j"
	for(sim in 1:nSim){
		total_domplan_sim[j,1,sim] <- cualdomplan[j]
		total_domplan_sim[j,3,sim] <- N_j
		}
	
	#	``conteo''
	conteo_sal <- conteo(datos_j,N_j)

	#	Estimating the size of the population outside of the sample (in "Stil_j")
	N_S_j <- conteo_sal[[2]]

	#	Estimating the size of the population outside of the sample (in "Stil_j")
	N_Stil_j <- conteo_sal[[3]]

	if(N_Stil_j >0){
		#	``unicstar''
		unicstar_sal <- unicstar(datos_j)
		ystar <- unicstar_sal[[2]][,"ystar"]
	
		#	``pesos''
		pesos_sal <- pesosDP(unicstar_sal,alphaDP[j])
		rho <- pesos_sal[[1]]
		phi <- pesos_sal[[2]]
																								
	#-----------------------------------	
	#	Validation ``g0_licitacion_sal''	
	#-----------------------------------	
	#    i)	Gamma
	#	ii)	Weibull
	#  iii)	Log-normal
	#   iv)	Inverse-Gaussian 	 

	dis_sel <- list("Gamma","Weibull","Lognormal","Inverse-Gaussian")
	
	#	Validation of distribution
	
 	if(g0_licitacion_sal[[j]][[1]] == dis_sel[[1]] ||g0_licitacion_sal[[j]][[1]] == dis_sel[[2]] ||g0_licitacion_sal[[j]][[1]] == dis_sel[[3]] ||g0_licitacion_sal[[j]][[1]] == dis_sel[[4]]){
  	}else{
		stop("Error in the specification of 'distribution'!!!")}
			
	#	Validation of parameters
	if(g0_licitacion_sal[[j]][[2]][[1]]>0 & g0_licitacion_sal[[j]][[2]][[2]]>0){
	}else{
		stop("Error in the specification of 'parameters'!!!")}  

		#	Simulation the T_Stil_j
		domplan_total_sim <- domplan_total(datos_j,rho,ystar,phi,g0_licitacion_sal[[j]],N_Stil_j,nSim)
	
		#	Simulating data grouped "n_j" distributed "domnoplan1"
		T_S_domplan_j <- domplan_total_sim[[1]]
		T_Stil_domplan_j <- domplan_total_sim[[2]]
		total_domplan_sim[j,2,] <- matrix(t(T_S_domplan_j),1,nSim) + T_Stil_domplan_j
	}else if(N_Stil_j ==0){
		# 	Calculating the total "T_S_j" compositional in S_j (within the sample)
		T_S_j_gpo <- datos_j[,"n_i"] * datos_j[,"y_i"]

		#	Composition of the total between domains on the sample unplanned "S_j"
		T_S_domplan_j <- sum(T_S_j_gpo)

		#	Simulating data store "n_j" distributed "domnoplan1"
		total_domplan_sim[j,2,] <- matrix(t(T_S_domplan_j),1,nSim)
		}
		
	}

#----------------------------	
#	Output
#----------------------------	
return(total_domplan_sim)

#
#	--	END of	"domplan_g0.R"--
}
