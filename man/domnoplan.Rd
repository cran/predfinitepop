\name{domnoplan}
\alias{domnoplan}
\title{Finite population inference on unplanned domains using reference data.}
\description{Generates Monte Carlo samples from the predictive distribution of totals of a finite population segmented in planned and unplanned domains, along with simulations of the predictive distribution for the composition of the population between the unplanned domains.}
\usage{
       domnoplan(datos,datos_ant,domplan_N,alphaDP,colid_D,alpha_D,inter,part,nSim)
      }
       
\arguments{
          \item{datos}{(\eqn{M} X \eqn{p})-dimensional array with positive entries for \eqn{\mathcal{S}_j}}
          \item{datos_ant}{(\eqn{M} X \eqn{p})-dimensional reference array for calibration of \eqn{G_{0}}}
          \item{domplan_N}{Matrix array with counts of individuals in each planned domain}
          \item{alphaDP}{\eqn{J}-dimensional array with positive entries for the parameters of the Dirichlet process for \eqn{F_j} (with \eqn{J} being the number of planned domains)}
          \item{colid_D}{\eqn{D}-dimensional matrix array with the columns in datos that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of "datos")}
          \item{alpha_D}{(\eqn{J} X \eqn{D})-dimensional array with positive entires for the parameters of the multinomial-Dirichlet component for the composition across unplanned domains. Note: Each one of the \eqn{J} rows is a vector of composition for \eqn{\mathcal{P}_j} divided across the \eqn{D} unplanned domains}
          \item{inter}{Tuning parameter for model comparison and selection (related to calibration of \eqn{G_{0}})}
		  \item{part}{Number of partitions for predictive cross-validation (related to calibration of \eqn{G_{0}})}
          \item{nSim}{Number of Monte Carlo simulated replicates of the predictive distribution}
          }

\details{
			- datos : Represent the data sample of the target population, unplanned domains labeled.
			
  					    It should contain the following columns:
						
						"domplan" - Categories for planned domains.
							
						"y_i" - Actual individual measurements (for the moment, they must be positive and continuous) for the group of observation.
							
						"n_i" - Number of individuals in the group (if the unit of observation in the sample are individuals, then "\eqn{n_{i}}" must be equal to 1)
			
			- datos_ant : Represents the data reference used to calibrate \eqn{G_{j0}}
			
  						The data must be labeled by planned domains. It should contain the following columns:
						
						"domplan" - Categories of planned domains
						
						"y_i" - Positive real and individual measurements of each group in the sample (when the units of observation are the groups, "\eqn{y_i}" should be per capita measurement)

            - domplan_N : Represents counts (or reference population) of the target population, divided by the planned domains.
			
						Tagged data must be labeled by domains planned. It should contain the following columns:
						
						"domplan" - Categories fot planned domains.
						
						"N_j" - Number of individuals in each planned domain.
		}		  
		  
\value{The function 'domnoplan' produces:
  
  \item{total_domnoplan_sim}{Matrix array of dimension "\eqn{J} * (3 +2 \eqn{D})* nSim" with predictions for the planned domains.
  
								Column 1 - Indicator of the planned domains
	
								Column 2 - \eqn{T_{j}} (totals of the planned domains)
								
								Column 3 - \eqn{N_{j}} (composition of the planned domains)

								Column 4 - \eqn{(4+D-1) - T_{d_j}} (totals of unplanned domains, such that  
											\eqn{T_{j}=} \eqn{\displaystyle\sum_{j=1}^J T_{d_j}}(across d))
																
								Column 5 - \eqn{(4+D)} to \eqn{(3 + 2*D) - N_{d_j}} (composition of unplanned domains, such  
											that \eqn{N_{j}= }\eqn{\displaystyle\sum_{j=1}^J N_{d_j}} (across \eqn{d}))
							}
       }

\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
}   
\author{
	Setgio I. Olivares-Guzman, Adriana Roldan-Rodriguez, Juan Carlos Martinez-Ovando
	   }

\seealso{
  \code{\link[predfinitepop:domnoplan_composicion]{domnoplan_composicion,}}
  \code{\link[predfinitepop:domnoplan_totalcomp]{domnoplan_totalcomp.}}
        }

\examples{
#--------------------------------------------------------
#	Uses of library: 	"predfinitepop"
#--------------------------------------------------------

#             Read files with functions in R
rm(list=ls())

library("predfinitepop")

set.seed(12345)

#-------------------------------------------------------------
#		DATA
#-------------------------------------------------------------
data(dataref)
data(datasample)
data(popcomp)

#		Sample (with 10% of the population)

datos <- datasample[[1]]
datos_ant <- dataref
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
#		CODE PARAMETERS
#-------------------------------------------------------------

#	A. Interval for calibration of G_0
inter <- 3 

#	B. Number of groups for predictive validation and model comparison
part <- 20

#	C. Number of Monte Carlo replicates of the predictive distributions
nSim <- 10			

#-------------------------------------------------------------
#		TESTING "domnoplan.R"
#-------------------------------------------------------------
date.ini <- date()
domnoplan_sim <- domnoplan(datos,datos_ant,domplan_N,alphaDP,colid_D,alpha_D,inter,part,nSim)
date.fin <- date()

#
#	--	END --
		  }

\keyword{Total Population}
\keyword{Convolutions}
\keyword{Domains planned}
\keyword{Compositional analysis}
\keyword{Domains unplanned}