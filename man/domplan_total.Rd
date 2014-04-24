\name{domplan_total}
\alias{domplan_total}
\title{Monte Carlo simulation of total for a given planned domain.}
\description{Simulates Monte Carlo samples of size '\eqn{nSim}' of the final distribution for the total for the given planned domain \eqn{\mathcal{P}_j}.}
\usage{domplan_total(datos_j,rho,ystar,phi,g0_licitacion_sal,N_Stil_j,nSim)}
\arguments{
            \item{datos_j}{Data matrix with features and number of individuals in the sample '\eqn{\mathcal{S}_j}'}
            \item{rho}{(\eqn{U} X 2)-dimensional vector with weights associated with the sample ties '\eqn{ystar}'}
            \item{ystar}{Sample ties \eqn{y^*_i}, in the sample \eqn{\mathcal{S}_j}}
            \item{phi}{Probability weight asociated with '\eqn{G_{j0}}'(the continuos part of '\eqn{\hat{G}_j}')}
            \item{g0_licitacion_sal}{Object list with the continuous component in '\eqn{\hat{G}_j}'}
            \item{N_Stil_j}{Composition of individuals out of '\eqn{\mathcal{S}_j}'}
            \item{nSim}{Number of Monte carlosimulated replicates of the predictive distribution}
           }
\value{The function 'domplan_total' produces and object list with two entries:
  
  \item{T_S_j}{Composition of '\eqn{T_j}' of the given planned domains in the sample '\eqn{S_j}'}
  \item{domplan_total_sim}{(1 * nSim) matrix with samples from the predictive distribution of the total 'T_Stil_j'}
  
}

\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
	Juan Carlos Martinez-Ovando
	}
\seealso{
		\code{\link[predfinitepop:domnoplan]{domnoplan}}
		}

