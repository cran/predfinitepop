g0_licitacion <-
function(datos_ant,inter,part){
#   This function computes the prior elicitation of \eqn{G_{j0}} using reference data for the planned domain \eqn{j}. 
#	\eqn{G_{j0}} is used by the SSM as the baseline function. The distribution is elicited using a predictive cross-validation 
#	procedure for model comparison and selection among the following alternatives: Gamma, Weibull, Lognormal and Inverse-Gaussian. 
#	The function also computes the expectation of the chose distribution.
#             
#	input:
#		     datos_ant    -    Reference data for calibration of G_0
#		     inter		  -	   Tuning parameter for model comparison and selection
#							         
#
#	output:  List object:
#		     'distribution'	-   String object for 'distribution'
#							          i) Gamma
#									 ii) Weibull
#									iii) Log-normal
#								     iv) Inverse-Gaussian 
#                        
#		     theta		-	 Vector of parameters associated with 'distribution'
#		     mu			-	 Expected value of 'distribution'	
#

y <- datos_ant[,"y_i"] 
n <- length(y)
n1 <- n/20
# Vector with assessments of the 20 groups with missing observations
  g0LogN20 <- rep(NA,20)
  g0Gam20  <- rep(NA,20)
  g0Wei20  <- rep(NA,20)
  g0IGau20 <- rep(NA,20)
# Vector with assessments of each observation using the estimated parameters with the observations that were omitted
  g0LogNN <- rep(NA,n)
  g0GamN  <- rep(NA,n)
  g0WeiN  <- rep(NA,n)
  g0IGauN <- rep(NA,n)
# initial parameters
  sigmaLN0c <- rnorm(1,0,1)
  muLN0c    <- rnorm(1,0,1)
  nu1G0c    <- rnorm(1,0,1)
  nu2G0c    <- rnorm(1,0,1)
  nu1W0c    <- rnorm(1,0,1)
  nu2W0c    <- rnorm(1,0,1)
  nu1IG0c   <- rnorm(1,0,1)
  nu2IG0c   <- rnorm(1,0,1)
# Randomization group
  y <- sample(y)
# Vector with 20 parameters and parameter groups of the entire group
  muLN0  <- rep(NA,21)
  sigLN0 <- rep(NA,21) 
  shG0   <- rep(NA,21)
  scG0   <- rep(NA,21)
  shW0   <- rep(NA,21)
  scW0   <- rep(NA,21)
  shIG0  <- rep(NA,21)
  scIG0  <- rep(NA,21)
  
i <- 1
for (i in 1:20) 
{   
# Ranges to form 20 groups
    i5 <- round((i-1)*n1)+1                                        
    i6 <- round(i*n1)
    i7 <- i6+1
    if (i6 <= n)
    {              
    for (j in (i5:i6)) 
      {
      # Missing observations according to previous ranks
        y[j] <- NA                                                 
      }
    }else{for (j in i7:n){y[j] <- NA}}
    y <- na.omit(y) 
    
  # Prior elicitation for the Weibull distribution                                   
    nOmiOrd <- length(sort(y))                                     
    x00 <- c(1:nOmiOrd)                                                  
    x0  <- log(log(1/(1-(x00/(nOmiOrd+1)))))
    x1  <- (1/nOmiOrd)*sum(x0)
    y1  <- mean(log(sort(y)))
    
  # ESTIMATES IN EACH GROUP:
             
  # Mu de X ~ Normal, where X=log(Y)                                                       
    muLN0[i]  <- sum(log(y))/length(y)       
  # Sigma of X ~ Normal, where X=log(Y)                       
    sigLN0[i] <- sqrt(sum((log(y)-muLN0[i])^2)/length(y))          
  # Form parameter of Gamma 
    shG0[i]   <- 0.5/(log(mean(y))-mean(log(y)))            
  # Scale parameter of Gamma       
    scG0[i]   <- mean(y)/shG0[i]                                    
  # Form parameter of Weibull
    shW0[i]   <- ((nOmiOrd*(sum((log(sort(y)))*(x0))))-(sum(x0)*sum(log(sort(y)))))/((nOmiOrd*sum(log(sort(y))^2))-((sum(log(sort(y))))^2)) 
  # Scale parameter of Weibull
    scW0[i]   <- exp(y1-(x1/shW0[i]))                               
  # Parameter "mu" of the Inverse Gaussian (not a location parameter)
    shIG0[i]  <- mean(y)                                            
  # Parameter "sigma" of the Inverse Gaussian (not a scale parameter)
    scIG0[i]  <- (sum((1/y)-(1/shIG0[i])))/length(y)                

    y <- datos_ant[,"y_i"]                                 
    
    g0LogNPaso <- rep(NA,n)
    g0GamPaso  <- rep(NA,n)
    g0WeiPaso  <- rep(NA,n)
    g0IGauPaso <- rep(NA,n)

    for (j in (i5:i6)) 
    {
      # Conditional densities
        g0LogNPaso[j] <- (pnorm(log(y[j]+inter),mean=muLN0[i],sd=sigLN0[i])-pnorm(log(max(y[j]-inter,0)),mean=muLN0[i],sd=sigLN0[i]))/(2*inter)       
        g0GamPaso[j]  <- (pgamma(y[j]+inter,shape=shG0[i],scale=scG0[i])-pgamma(max(y[j]-inter,0),shape=shG0[i],scale=scG0[i]))/(2*inter)            
        g0WeiPaso[j]  <- (pweibull(y[j]+inter,shape=shW0[i],scale=scW0[i])-pweibull(max(y[j]-inter,0),shape=shW0[i],scale=scW0[i]))/(2*inter)         
        g0IGauPaso[j] <- (pinvgauss(y[j]+inter,shIG0[i],scIG0[i])-pinvgauss(max(y[j]-inter,0.000001),shIG0[i],scIG0[i]))/(2*inter)                          
      # Vectors with assessments of the "n" observations
        g0LogNN[j]  <- g0LogNPaso[j]                                 
        g0GamN[j]   <- g0GamPaso[j]
        g0WeiN[j]   <- g0WeiPaso[j]
        g0IGauN[j]  <- g0IGauPaso[j]
      # Vectors with assessments of the observations of the group "i" (20 groups)
        g0LogNsinNA  <- na.omit(g0LogNPaso)                          
        g0GamsinNA   <- na.omit(g0GamPaso)
        g0WeisinNA   <- na.omit(g0WeiPaso)
        g0IGausinNA  <- na.omit(g0IGauPaso)
    }
  # Vectors with the sum of the logarithms of the evaluations of the observations of the group "i"
    g0LogN20[i] <- sum(log(g0LogNsinNA))                            
    g0Gam20[i]  <- sum(log(g0GamsinNA))
    g0Wei20[i]  <- sum(log(g0WeisinNA))
    g0IGau20[i] <- sum(log(g0IGausinNA))
} 
# Average of the 20 sums of logarithms
criLN <- mean(g0LogN20)                                            
criG  <- mean(g0Gam20)
criW  <- mean(g0Wei20)  
criIG <- mean(g0IGau20)

# ESTIMATE FOR THE WHOLE SAMPLE :

# Previous calculations for Weibull parameter estimation                                             
x00  <- c(1:n)                                                
x0   <- log(log(1/(1-(x00/(n+1)))))
x1   <- (1/n)*sum(x0)  
y1   <- mean(log(y))
    
# Mu of Y ~ Lognormal (mu_{y})
muLN0[21]  <- sum(log(y))/length(y)   
# Sigma of Y ~ Lognormal (sigma_{y})                             
sigLN0[21] <- sqrt(sum((log(y)-muLN0[21])^2)/length(y))          
# Form parameter of Gamma  
shG0[21]   <- 0.5/(log(mean(y))-mean(log(y)))                      
# Scale parameter Gamma
scG0[21]   <- mean(y)/shG0[21]                                     
#Form parameter of Weibull
shW0[21]   <- ((n*(sum((log(sort(y)))*(x0))))-(sum(x0)*sum(log(sort(y)))))/((n*sum(log(sort(y))^2))-((sum(log(sort(y))))^2)) 
# Scale parameter of Weibull
scW0[21]   <- exp(y1-(x1/shW0[21]))                                
# Parameter "mu" of the Inverse Gaussian (not a location parameter)
shIG0[21]  <- mean(y)                                             
# Parameter "sigma" of the Inverse Gaussian (not a scale parameter) 
scIG0[21]  <- (sum((1/y)-(1/shIG0[21])))/length(y)                 

# COMPARISON CRITERIA:
     
# Including all criteria
#MaxCri <- max(criLN, criG, criW, criIG)  
# We removed the Gamma distribution criteria        
MaxCri <- max(criLN, criW, criIG)   
# Selected Distribution  
DenSel <- NA                          
# Parameters of the selected distribution
ParSel <- 0                           
# Expectation of selected distribution
MuSel <- 0                            
if((MaxCri==criLN)) {DenSel="Lognormal"
ParSel <- c(muLN0[21],sigLN0[21])
MuSel <- exp(muLN0[21]+(0.5*(sigLN0[21]^2)))}else{if((MaxCri==criW)) {DenSel="Weibull"
ParSel <- c(shW0[21],scW0[21])
MuSel <- shG0[21]*scG0[21]}else{if(MaxCri==criG) {DenSel="Gamma"
ParSel <- c(shG0[21],scG0[21])
MuSel <- shG0[21]*scG0[21]}else{if(MaxCri==criIG) {DenSel="Inverse-Gaussian"
ParSel <- c(shIG0[21],scIG0[21])
MuSel <- shIG0[21]}else{DenSel="Adjustment is not available"}}}}

g0_licitacion_sal <- list(DenSel,ParSel,MuSel)

#	Output:
salida <- g0_licitacion_sal
return(salida)

#	END of g0_licitacion.R
}
