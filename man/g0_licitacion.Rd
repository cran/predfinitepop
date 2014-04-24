\name{g0_licitacion}
\alias{g0_licitacion}
\title{Prior elicitation of \eqn{G_{j0}} using reference data.}
\description{This function computes the prior elicitation of \eqn{G_{j0}} using reference data for the planned domain \eqn{j}. \eqn{G_{j0}} is used by the SSM as the baseline function. The distribution is elicited using a predictive cross-validation procedure for model comparison and selection among the following alternatives: Gamma, Weibull, Lognormal and Inverse-Gaussian. The function also computes the expectation of the chose distribution.}
\usage{
		g0_licitacion(datos_ant,inter,part)
      }
\arguments{
  \item{datos_ant}{Reference data for calibration of \eqn{G_{j0}}}
  \item{inter}{Tunning parameter for model comparison and selection}
  \item{part}{Number of groups for predictive validation. The larger the number of groups the fewer observations are to be removed for calibration}

          }
\details{
			datos_ant: Data collected in a previous time or from another source for each planned domain. It must contain two columns, labeled by 'domplan' and 'y_i'
				'domplan' - Labels of planned domains
				'y_i' - individual outcome (for the sake of implementation, this must be positive)
		}	
\value{The function 'g0_licitacion' produces an object list with three elements:
		
		  \item{distribution }{List arrangement:
									String object for 'distribution'
							              i)  Gamma
							             ii)  Weibull
							            iii)  Log-normal
							             iv)  Inverse Gaussian
							  }
		  \item{theta }{Vector of parameters associated with 'distribution'}
		  \item{mu }{Expected value of 'distribution'}

}
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
	Adriana Roldan-Rodriguez
	}
\keyword{Model comparison and selection}
\keyword{Prior elicitation} 
\keyword{Calibration}