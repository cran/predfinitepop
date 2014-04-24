\name{domnoplan_composicion}
\alias{domnoplan_composicion}
\title{Generates Monte Carlo samples from the predictive distribution for the vector \eqn{\boldmath{N}_j} across the \eqn{D} unplanned domains within a given planned domain.}
\description{This function simulates Monte Carlo samples from the predictive distribution of the vector \eqn{\boldmath{N}_j} across the \eqn{D} unplanned domains within a given planned domain \eqn{j}.}

\usage{
      domnoplan_composicion(datos_j,N_j,N_S_j,N_Stil_j,colid_D,alpha_D,nSim)
      }

\arguments{
			  \item{datos_j}{Data matrix with features and number of individuals in the sample \eqn{\mathcal{S}_j}}
			  \item{N_j}{Number of individuals in \eqn{\mathcal{P}_j} (\eqn{j}th planned domain)}
			  \item{N_S_j}{Number of individuals in the sample \eqn{\mathcal{S}_j}}
			  \item{N_Stil_j}{Number of individuals out the sample \eqn{\mathcal{S}_j}}
			  \item{colid_D}{\eqn{D}-dimensional matrix array with the columns in 'datos' that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of 'datos')}
			  \item{alpha_D}{\eqn{D}-dimensional array with positive entires for the parameters of the multinomial-Dirichlet component for the composition across unplanned domains (NOTE: this one makes reference to a single planned domain)}
			  \item{nSim}{Number of Monte Carlo simulated replicates of the predictive distribution}
          }

\value{This function 'domnoplan_composicion' produces:

			\item{N_S_domnoplan }{Composition of the number of individuals in saample \eqn{\mathcal{S}_j} (\eqn{N^\mathcal{S}_j}) across the \eqn{D} unplanned domains}
			\item{domnoplan_composicion_sim }{(1 X \eqn{D} X nSim) matrix with Monte Carlo samples of the predictive distribution of the composition of \eqn{\tilde{\mathcal{S}}_j} }
	  }
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
	}
\author{Juan Carlos Martinez-Ovando}
\seealso{
		\code{\link[predfinitepop:domnoplan]{domnoplan}}
		}

