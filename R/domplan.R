domplan <-
function(datos,datos_ant,domplan_N,alphaDP,inter,part,nSim){
#
#	Generates Monte Carlo samples of the predictive distribution of totals of a finite population segmented in planned domains.
#
#	Input:
#		datos		-	p-dimensional vector with positive entries		
#		datos_ant	-	Reference data for calibration of G_0
#		domplan_N	-	Free dimension matrix array represents counts
#		alphaDP		-   J-dimensional positive parameter for the Dirichlet process (with J being the number of planned domains)
#		inter		-   Tuning parameter for model comparison and selection
#		part		-	Number of partitions for predictive cross-validation
#		nSim		-   Number of simulated replicates
#
#	Details:
#			- datos : Represents the data sample of the target population, unplanned domains labelled.
# 					    It should contain the following columns:
#							"domplan" - Planned domains categories.
#							"y_i"     - Actual measurements of individual positive and the sample group.
#							"n_i"     - Number of members in the group (if the unit of observation in the sample are individuals, then "n_i" must be equal to 1)
#			
#			- datos_ant : Represents the data of reference used to calibrate G_{0}.
# 						The data must be labelled by domains planned. It should contain the following columns:
#							"domplan" - Planned categories of domains
#							"y_i"     - Positive real and individual measurements of each group in the sample (when the units of observation are groups, "y_i" should be a per capita measurement)
# 
#            - domplan_N : Represents counts (or population frame) of the target population, divided by the planned domains.
#						Tagged data must be labelled by domains planned. It should contain the following columns:
#							"domplan" - Planned domains categories
#							"N_j"     - Number of individuals in each population planned domain
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
	
#	C.	Parameter inter
if(inter <= 0){
	stop("Error in the specification of 'inter'!!!")}

#	-----------------------------------------------
#	Scanning planned domains
#	-----------------------------------------------
j <- 1
for(j in 1:J){
	#	Extraction of indexes in "datos_ant"
	P_j_ant <- which(datos_ant[,"domplan"]==cualdomplan[j])
	#	Extraction of the data in "datos_ant"
	datos_j_ant <- datos_ant[P_j_ant, ]

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
		
		#	``g0_licitacion''
		g0_licitacion_sal <- g0_licitacion(datos_j_ant,inter,part)
		#g0_licitacion_sal <- g0_licitacion(datos_j_ant,inter)

		#	Simulation the T_Stil_j
		domplan_total_sim <- domplan_total(datos_j,rho,ystar,phi,g0_licitacion_sal,N_Stil_j,nSim)
	
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
#	--	END of	"domplan.R"--
}
