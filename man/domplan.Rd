\name{domplan}
\alias{domplan}
\title{Finite population inference on unplanned domains using reference data.}
\description{Generates Monte Carlo samples of the predictive distribution of totals of a finite population segmented in planned domains.}
\usage{
		domplan(datos,datos_ant,domplan_N,alphaDP,inter,part,nSim)
	}
\arguments{
            \item{datos}{(\eqn{M} X \eqn{p})-dimensional array with positive entries for \eqn{\mathcal{S}_j}}
            \item{datos_ant}{(\eqn{M} X \eqn{p})-dimensional refenrence array for calibration of \eqn{G_0}}
            \item{domplan_N}{Matrix array with counts of individuals in each planned domain}
            \item{alphaDP}{\eqn{J}-dimensional array with positive entries for the parameters of the Dirichlet process for \eqn{F_j} (with \eqn{J} being the number of planned domains)}
            \item{inter}{Tunning parameter for model comparison and selection (related to calibration of \eqn{G_{0}})}
			\item{part}{Number of partitions for predictive cross-validation (related to calibration of \eqn{G_{0}})}
            \item{nSim}{Number of Monte Carlo simulated replicates of the predictive distribution}
          }
		  
\details{
			- datos : Represents the data sample of the target population, unplanned domains labeled.
			
  					    It should contain the following columns:
						
						"domplan" - Categories fot planned domains.
						
						"y_i" - Actual individual measurements(for the moment, they must be positive) for the group of observation.
						
						"n_i" - Number of individuals in the group (if the unit of observation in the sample are individuals, then "\eqn{n_{i}}" must be equal to 1)
			
			- datos_ant : Represents the data of reference used to calibrate \eqn{G_{0}}.
			
  						The data must be labeled by domains planned. It should contain the following columns:
			
						"domplan" - Planned categories of domains
						
						"y_i" - Positive real and individual measurements of each group in the sample (when the units of observation are groups, "\eqn{y_{i}}" should be a per capita measurement)
  
            - domplan_N : Represents counts (or reference population) of the target population, divided by the planned domains.
						
						Tagged data must be labeled by domains planned. It should contain the following columns:
						
						"domplan" - Categories fot planned domains.
						
						"N_j" - Number of individuals in each planned domain.
		}
		  
\value{The function 'domplan' produces:
  
  \item{total_domnoplan_sim}{Matrix array of dimension "J * (3 +2 D)* nSim" with predictions for the planned domains.
  
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

\keyword{Total Population}
\keyword{Convolutions}
\keyword{Domains planned}

\seealso{
		\code{\link[predfinitepop:domplan_total]{domplan_total}}
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
#		TESTING "domplan.R"
#-------------------------------------------------------------
date.ini <- date()
domplan_sim <- domplan(datos,datos_ant,domplan_N,alphaDP,inter,part,nSim)
date.fin <- date()

#
#	--	END	--
}

