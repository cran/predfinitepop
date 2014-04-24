domnoplan_g0 <-
function(datos,domplan_N,alphaDP,colid_D,alpha_D,nSim,g0_licitacion_sal){
#
#	Generates Monte Carlo samples of the predictive distribution of totals of a finite population 
#	segmented in planned and unplanned domains, using a predefined set of \eqn{(G_{j0})_{j=1}^{J}}, 
#	along with simulations of the predictive distribution for the composition of the population 
#	between the unplanned domains.
#					
#	Input:              
#		datos		-	(Mxp)-dimensional array with positive entries for S_j
#		domplan_N	-	Matrix array with counts of individuals in each planned domain
#		alphaDP		-   J-dimensional array with positive entries for the parameters of the Dirichlet process for F_j (with J being the number of planned domains)
#		colid_D		-   D-dimensional matrix array with the columns in ´datos´ that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of 'datos')
#		alpha_D		-   (JxD)-dimensional array with positive entries for the parameters of the multinomial-Dirichlet component for the composition across unplanned domains.
#						Note: Each one of the J rows is a vector of composition for P_j segmented across D unplanned domains
#       nSim		-   Number of Monte Carlo simulated replicates of the predictive distribution
# g0_licitacion_sal -   (Jx1) object list, each entry is another object list itself associated with each G_{j0} for the J planned domains. The first element for arch G_{j0} should be the name of the chosen ´distribution´(see details below for alternatives), the second element should be a vector object with the parameters associated with ´distribution´, and the third element should be its associated expectation 	
#
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
	  
	#	Simulating the composition of "T_Stil_j" between domains unplanned
	domnoplan_totalcomp_sim <- domnoplan_totalcomp(datos_j,rho,ystar,phi,g0_licitacion_sal[[j]],N_Stil_domnoplan_j_sim,nSim,colid_D)
		
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
#	--	End of	"domnoplan_g0.R"--
}
