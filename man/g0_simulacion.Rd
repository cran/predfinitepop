\name{g0_simulacion}
\alias{g0_simulacion}
\title{Computes Monte Carlo samples from a given \eqn{G_{j0}}.}
\description{This function computes Monte Carlo samples from a given continuous \eqn{G_{j0}}.}

\usage{
		g0_simulacion(g0_licitacion_sal,nSim)
	  }

\arguments{
            \item{g0_licitacion_sal}{Object list with three elements (produced with 'g0_licitacion'):
			\tabular{ll}{
								\tab a)	String object for 'distribution'\cr
									 \tab 		i)	Gamma\cr
									 \tab 		ii) Weibull\cr
									 \tab 		iii) Log-normal\cr
									 \tab 		iv)  Inverse-Gaussian\cr 
								\tab b) theta	-	Vector of parameters associated with 'distribution'\cr
								\tab c) mu		-	Expected value of 'distribution'\cr
						}
						             }		
            \item{nSim}{Number of Monte Carlo simulations}
          }
\value{This function 'g0_simulacion' produces:
		
		  \item{g0_simulacion_sal}{(nSim X 1) dimensional array with simulated data}
		  
	  }
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
		Adriana Roldan-Rodriguez, Juan Carlos Martinez-Ovando
		}
\seealso{
		\code{\link[predfinitepop:ghat_simulacion]{ghat_simulacion}}
		}
\keyword{Prediction}
\keyword{Monte Carlo simulation}