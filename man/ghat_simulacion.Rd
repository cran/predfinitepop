\name{ghat_simulacion}
\alias{ghat_simulacion}
\title{Simulates Monte Carlo samples from \eqn{\hat{G}_j}.}
\description{This function simulates Monte Carlo samples from \eqn{\hat{G}_j}.}
\usage{ghat_simulacion(rho,ystar,phi,g0_licitacion_sal,N_Stil_j)}
\arguments{
              \item{rho}{(\eqn{U} X 2)-dimensional vector with weights associated with 'ystar'}
              \item{ystar}{Sample ties \eqn{y^*_i} in \eqn{\mathcal{S}_j}} 
              \item{phi}{Probability weight asociated with \eqn{G_{j0}} (the continuos part of \eqn{\hat{G}_j})}
              \item{g0_licitacion_sal}{Object list with three elements (produced with 'g0_licitacion'):
			                               \tabular{ll}{
                          						 \tab a)String object for 'distribution'\cr	
															\tab	i) Gamma\cr
															\tab   ii) Weibull\cr
															\tab  iii) Log-normal\cr
															\tab  iv) Inverse-Gaussian\cr
                        						 \tab b) theta			-	Parameters associated with 'distribution'\cr
                        						 \tab c) mu				-	Expected value of 'distribution'\cr
												        }
										}				
              \item{N_Stil_j}{Number of samples to simulate (number of individuals in \eqn{\tilde{\mathcal{S}}_j}) out of the sample}
           }
\value{'ghat_simulacion' produces:		
		  \item{ghat_simulacion_sal}{(N_Stil_j X 1) dimensional vector with simulated data}
		  
	  }
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
		Juan Carlos Martinez-Ovando
		}
\seealso{
		\code{\link[predfinitepop:g0_simulacion]{g0_simulacion,}}
		\code{\link[predfinitepop:gstar_simulacion]{gstar_simulacion.}}
		}
\keyword{Prediction}
\keyword{Monte Carlo simulation}