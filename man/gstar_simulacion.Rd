\name{gstar_simulacion}
\alias{gstar_simulacion}
\title{Generates Monte Carlo samples from the discrete component of \eqn{\hat{G}_j}}
\description{This function generates Monte Carlo samples from the discrete component of \eqn{\hat{G}_j} for the Dirichlet process prior in the planned domain \eqn{j}.}
\usage{                                     
		gstar_simulacion(rho,ystar,nSim)
	  }

\arguments{
              \item{rho}{(\eqn{U} X 2)-dimensional vector with weights associated with 'ystar'}
              \item{ystar}{Sample ties \eqn{y^*_i} in \eqn{\mathcal{S}_j}} 
              \item{nSim}{Number of Monte Carlo simulations}
          }
\value{The function 'gstar_simulacion' produces:
		
		  \item{gstar_simulacion_sal}{(nSim X 1) dimensional vector with simulated data}
		  
	  }

\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
		Juan Carlos Martinez-Ovando
		}
\seealso{
		\code{\link[predfinitepop:ghat_simulacion]{ghat_simulacion}}
		}
\keyword{Prediction}
\keyword{Monte Carlo simulation}