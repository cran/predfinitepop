\name{domplan_g0}
\alias{domplan_g0}
\title{Finite population inference on unplanned domains using subjective \eqn{G_{j0}}.}
\description{Generates Monte Carlo samples of the predictive distribution of totals of a finite population segmented in planned domains, using a predefined set of baseline distributions \eqn{(G_{j0})_{j=1}^{J}}.}
\usage{domplan_g0(datos,domplan_N,alphaDP,nSim,g0_licitacion_sal)}
\arguments{
            \item{datos}{(\eqn{M} X \eqn{p})-dimensional array with positive entries for \eqn{\mathcal{S}_j}}
            \item{domplan_N}{Matrix array with counts of individuals in each planned domain}
            \item{alphaDP}{\eqn{J}-dimensional array with positive entries for the parameters of the Dirichlet process for \eqn{F_j} (with \eqn{J} being the number of planned domains)}
            \item{nSim}{Number of Monte Carlo simulated replicates of the predictive distribution}
			\item{g0_licitacion_sal}{(\eqn{J} X 1) object list, each entry is itself an object list associated with each \eqn{G_{j0}} for the \eqn{J} planned domains. The first element for arch \eqn{G_{j0}} should be the name of the chosen distribution (see details below for the admissible alternatives), the second element should be a vector object with the parameters associated with "distribution", and the third element should be its associated expectation}
		  }                          
		  
\details{
			- datos : Represents the data sample of the target population, unplanned domains labeled.
			
  					    It should contain the following columns:
						
						"domplan" - Categories fot planned domains.
						
						"y_i"  - Actual individual measurements/putcomes (for the moment, they must be positive) for the group of observation.
						
						"n_i"  - Number of individuals in the group (if the unit of observation in the sample are individuals, then "\eqn{n_{i}}" must be equal to 1)
							
            - domplan_N : Represents the counts (or number of individuals) in the target population, segmented in the planned domains.
			
						Data must be labeled by the planned domain. It should contain the following columns:
						
						"domplan" - Categories fot planned domains.
						
						"N_j" - Number of individuals in each planned domain.
							
			- g0_licitacion_sal: Chose one and only one of the distribution,
			
						i)	Gamma
						
						ii)	Weibull

					   iii) Lognormal

						iv)	Inverse-Gaussian 			
						
							Parameterizations for distribution:
							
						 i)	\eqn{\tt Gamma} distribution, with parametetrs \eqn{\boldsymbol{\theta} = c(\alpha>0 , \beta>0)} and density function     
							
							\eqn{f(x) \propto x^{\alpha-1} \exp\{-x/\beta\},}

							where \eqn{\alpha} is the shape parameter, and \eqn{\beta} is the scale parameter.

						ii)	\eqn{\tt Weibull} distribution, with parametetrs \eqn{\boldsymbol{\theta} = c(\alpha>0 , \beta>0)} and density function     
							      
							\eqn{f(x) \propto (x/\beta)^{\alpha-1} \exp\{-(x/\beta)^\alpha\},}
							
							where \eqn{\alpha} is the shape parameter, and \eqn{\beta} is the scale parameter.

					   iii)	\eqn{\tt Lognormal} distribution, with parametetrs \eqn{\boldsymbol{\theta} = c(\alpha>0 , \beta>0)} and density function     
							  
							\eqn{f(x) \propto \exp\{-(\mbox{log}x-\alpha)^2/2\beta^2\},}
																																		
							where \eqn{\alpha} is the mean, and \eqn{\beta} is the standard deviation of the logarithm.

					    iv)	\eqn{\tt Inverse-Gaussian } distribution, with parametetrs \eqn{\boldsymbol{\theta} = c(\alpha>0 , \beta>0)} and density function 
							    
							\eqn{f(x) \propto x^{-3/2} \exp\{-\alpha(x-\beta)^2/2x\beta^2\},}
																					
							where \eqn{\alpha} is the shape parameter, and \eqn{\beta} is the mean parameter.

		}
		  
\value{The function 'domplan' produces:
  
  \item{total_domnoplan_sim}{Matrix array of dimension "\eqn{J} * (3 +2 \eqn{D})* nSim" with predictions for the planned domains.
  
									Column 1 - Indicator of the planned domains
									
									Column 2 - \eqn{T_{j}} (totals of the planned domains)
									
									Column 3 - \eqn{N_{j}} (composition of the planned domains)
							}

     }

\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
		Sergio I. Olivares-Guzman, Adriana Roldan-Rodriguez, Juan Carlos Martinez-Ovando
		}

\keyword{Population Total}
\keyword{Convolutions}
\keyword{Planned domains}

\seealso{
		\code{\link[predfinitepop:domplan_total]{domplan_total}}
		}
\examples{
#--------------------------------------------------------
#	Uses of library: 	"predfinitepop"
#--------------------------------------------------------

#             Read files with functions in R
rm(list=ls())

#library("predfinitepop")

set.seed(12345)

#-------------------------------------------------------------
#		DATA
#-------------------------------------------------------------
data(datasample)
data(popcomp)

#		Sample (with 10% of the population)

datos <- datasample[[1]]
domplan_N <- popcomp[[1]]

#-------------------------------------------------------------
#		MODEL PARAMETERS
#-------------------------------------------------------------

#	a.	Parameter "alpha" of the Dirichlet process
#		(this is for J02 planned domains)
alphaDP <- matrix(0.03,2,1)   

#	b. vector with indicator variables for unplanned domains
colid_D	<- t(as.matrix(c(6:8)))	

#	c. Matrix with parameters for unplanned domains
#		(this is for J=2 planned domains)
alpha_D <- matrix(0.3,2,length(colid_D))

#-------------------------------------------------------------
#		CODE PARAMETERS AND INCORPORATION OF THE DISTRIBUTION
#-------------------------------------------------------------
#  Incorporate each element ("g0_licitacion_in_J") according to the number of domains planned
#    g0_licitacion_in_1 <- list("distribution",teta1,esp(teta1))
# 
# In this case we have two domains planned 
g0_licitacion_in_1 <- list("Lognormal",c(6.198,0.840),.5)  
g0_licitacion_in_2 <- list("Gamma",c(6.198,0.840),.5)  
				
g0_licitacion_sal <- list(g0_licitacion_in_1,g0_licitacion_in_2)
	
#	B. Number of Monte Carlo replicates of the predictive distributions
nSim <- 10				
 
#-------------------------------------------------------------
#		TESTING "domplan.R"
#-------------------------------------------------------------
date.ini <- date()
domplan_sim <- domplan_g0(datos,domplan_N,alphaDP,nSim,g0_licitacion_sal)
date.fin <- date()

#
#	--	END	--
}
