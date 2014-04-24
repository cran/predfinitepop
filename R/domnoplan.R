domnoplan <-
function(datos,datos_ant,domplan_N,alphaDP,colid_D,alpha_D,inter,part,nSim){
#
#	Generates Monte Carlo samples of the predictive distribution of totals of a finite population 
#	segmented in planned and unplanned domains, along with simulations of the predictive distribution 
#	for the composition of the population between the unplanned domains.
#
#	Input:
#		datos		-	(Mxp)-dimensional array with positive entries for S_j
#		datos_ant	-	(Mxp)-dimensional reference array for calibration of G_0
#		domplan_N	-	Matrix array with counts of individuals in each planned domain
#		alphaDP		-   J-dimensional array with positive entries for the parameters of the Dirichlet process for F_j (with J being the number of planned domains)
#		colid_D		-   D-dimensional matrix array with the columns in ´datos´ that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of 'datos')
#		alpha_D		-   (JxD)-dimensional array with positive entries for the parameters of the multinomial-Dirichlet component for the composition across unplanned domains.
#						Note: Each one of the J rows is a vector of composition for P_j divided across the D unplanned domains
#		inter		-   Tuning parameter for model comparison and selection (related to calibration of G_0)
#		part		-	Number of partitions for predictive cross-validation (related to calibration of G_0)
#		nSim		-   Number of Monte Carlo simulated replicates of the predictive distribution
#
#	Details:
#			- datos : Represents the data sample of the target population, unplanned domains labelled.
#  					    It should contain the following columns:
#							"domplan" - Planned domains categories.
#							"y_i"     - Actual measurements of individual positive and the sample group.
#							"n_i"     - Number of members in the group (if the unit of observation in the sample are individuals, then "n_i" must be equal to 1)
#			
#			- datos_ant : Represents the data reference used to calibrate G_{j0}
# 						The data must be labelled by domains planned. It should contain the following columns:
#							"domplan" - Planned categories of domains
#							"y_i"     - Positive real and individual measurements of each group in the sample (when the units of observation are the groups, "y_i" should be per capita measurement)
#
#			- domplan_N : Represents counts (or reference population) of the target population, divided by the planned domains.
#						Tagged data must be labelled by domains planned. It should contain the following columns:
#							"domplan" - Planned domains categories
#							"N_j"     - Number of individuals in each population planned domain
#
#   Output:
#		total_domnoplan_sim		-	Matrix array of dimension "J x (3 + 2*D) x nSim" with predictions for relevant quantities of planned and unplanned domains
#									Column 1 - Indicator of the planned domains
#									Column 2 - T_j (totals of the planned domains)
#									Column 3 - N_j (composition of the planned domains)
#									Column 4 to (4+D-1) - T^d_j (totals of unplanned domains, such that T_j = sum  T^d_j (over d))
#									Column (4+D) to (3 + 2*D) - N^d_j (composition of unplanned domains, such that N_j = sum  N^d_j (over d))
#

#      Consulted domains
cualdomplan <- as.matrix(domplan_N[,"domplan"])
J <- length(cualdomplan)

#	Size unplanned domains
D <- length(colid_D)

#		Repository totals in "unplanned domains" with three categories
total_domnoplan_sim <- array(NaN,c(J,3+2*D,nSim))

#----------------------------	
#	validation
#----------------------------	

#	A.	Parameter alphaDP
if(any(alphaDP <= (0*alphaDP))){
	stop("Error in the specification of 'alphaDP'!!!")}
	
#	B.	Parameters alpha_D and colid_D
if(ncol(alpha_D) != ncol(colid_D)){
	stop("The dimenciones of 'alpha_D' and 'colid_D' are different!!!")}

#	C.	parameter alpha_D
if(any(alpha_D <= (0*alpha_D))){
	stop("Error in the specification of 'alpha_D'!!!")}

#	D.	Parameter nSim
if(nSim <= 0){
	stop("Error in the specification of 'nSim'!!!")}
	
#	E.	Parameter inter
if(inter <= 0){
	stop("Error in the specification of 'inter'!!!")}
	
#	F. Total domains unplanned
if(sum(datos[,colid_D]) != nrow(datos)) {
   stop("Error in total of unplanned domains")}

if(any(rowSums(datos[,colid_D]) != 1)) {
   stop("Error in total of unplanned domains")}
   
#	-----------------------------------------------
#	Scanning planned domains
#	-----------------------------------------------
j <- 1
for(j in 1:J){
	#	Extraction of indexes in "datos_ant"
	P_j_ant <- which(datos_ant[,"domplan"]==cualdomplan[j])
	#	Extraction of the data in "datos_ant"
	datos_j_ant <- datos_ant[P_j_ant, ]

    #	Extraction of indexes in "data" (i.e. sample S_j)
	P_j <- which(datos[,"domplan"]==cualdomplan[j])
	#	Extraction of the data in "data"
	datos_j <- datos[P_j, ]

	#	Population size in \mathcal{P}_{j}
	P_j <- which(domplan_N[,"domplan"]==cualdomplan[j])
	#	Extraction of the data in "data"
	N_j <- domplan_N[P_j,"N_j"]

	#	domain name the planned arrangement and constant data impute "n_j"
	sim <- 1
	for(sim in 1:nSim){
		total_domnoplan_sim[j,1,sim] <- cualdomplan[j]
		total_domnoplan_sim[j,3,sim] <- N_j
		}
	
	#	``conteo''
	conteo_sal <- conteo(datos_j,N_j)

	#	Estimating the size of the population out of the sample (in "Stil_j")
	N_S_j <- conteo_sal[[2]]

	#	Estimating the size of the population out of the sample (in "Stil_j")
	N_Stil_j <- conteo_sal[[3]]

	#	Simulating of the composition of "N_Stil_j" between domains unplanned
	domnoplan_composicion_sim <- domnoplan_composicion(datos_j,N_j,N_S_j,N_Stil_j,colid_D,alpha_D[j,],nSim)

	#	Simulating data"N_j" distributed in "domnoplan1"
	N_S_domnoplan_j <- domnoplan_composicion_sim[[1]]
	N_Stil_domnoplan_j_sim <- domnoplan_composicion_sim[[2]]
	total_domnoplan_sim[j,c((3+D+1):(3+2*D)),] <- matrix(t(N_S_domnoplan_j),D,nSim) + N_Stil_domnoplan_j_sim
	
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

	#	Simulating the composition of "T_Stil_j" between domains unplanned
	domnoplan_totalcomp_sim <- domnoplan_totalcomp(datos_j,rho,ystar,phi,g0_licitacion_sal,N_Stil_domnoplan_j_sim,nSim,colid_D)
		
	#	Simulating data grouped "n_j" distributed "domnoplan1"
	T_S_domnoplan_j <- domnoplan_totalcomp_sim[[1]]
	T_Stil_domnoplan_j <- domnoplan_totalcomp_sim[[2]]
	T_Stil_domnoplan_j_sim <- domnoplan_totalcomp_sim[[3]]
	total_domnoplan_sim[j,c((3+1):(3+D)),] <- matrix(t(T_S_domnoplan_j),D,nSim) + T_Stil_domnoplan_j_sim

	#	Totals for planned domains
	total_domnoplan_sim[j,c(2),] <- matrix(sum(T_S_domnoplan_j),1,nSim) + T_Stil_domnoplan_j
	
	}

#----------------------------	
#	Output
#----------------------------	
return(total_domnoplan_sim)

#
#	--	End of	"domnoplan.R"--
}
