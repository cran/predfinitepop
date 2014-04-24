\name{domnoplan_totalcomp}
\alias{domnoplan_totalcomp}
\title{Simulates Monte Carlo samples from the predictive distribution of the vector \eqn{\boldmath{T}_j} across the \eqn{D} unplanned domains within a given planned domain.}
\description{This function simulates Monte Carlo samples from the predictive distribution of the vector \eqn{\boldmath{T}_j} across the \eqn{D} unplanned domains in a given planned domain \eqn{j}.}
\usage{
domnoplan_totalcomp(datos_j,rho,ystar,phi,g0_licitacion_sal,N_Stil_domnoplan_j_sim,
					nSim,colid_D)
		}
\arguments{
              \item{datos_j}{Data matrix with features and number of individuals in the sample \eqn{\mathcal{S}_j}}
              \item{rho}{(\eqn{U} X 2)-dimensional vector with weights associated with the sample ties 'ystar'} 
              \item{ystar}{Sample ties \eqn{y^*_i} in the sample \eqn{\mathcal{S}_j}} 
              \item{phi}{Probability weight asociated with \eqn{G_{j0}} (the continuos part of \eqn{\hat{G}_j})}
              \item{g0_licitacion_sal}{Object list with the continuous component in \eqn{\hat{G}_j}}
              \item{N_Stil_domnoplan_j_sim}{(1 X \eqn{D} X nSim) matrix with Monte Carlo samples of the predictive distribution of the composition of \eqn{\tilde{\mathcal{S}_j}}}
              \item{nSim}{Number of Monte Carlo simulated replicates of the predictive distribution}
              \item{colid_D}{\eqn{D}-dimensional matrix array with the columns in 'datos' that correspond to the indicator variables of the planned domains (those indicator variables represent a partion of 'datos')}
          }
\value{The function 'domnoplan_totalcomp' produces an object list with three entries:
		\item{T_S_domnoplan}{Composition of \eqn{T_j} for the planned domain \eqn{j} in sample \eqn{\mathcal{S}_j}}
		\item{T_Stil_domnoplan_j}{(1 X \eqn{D} X nSim)-dimensional array with simulated samples of the convolution for \eqn{T^{\tilde{\mathcal{S}}_j}}}
		\item{domnoplan_totalcomp_sim}{(1 X \eqn{D} X nSim)-dimensional array with samples from the predictive distribution of the composition of \eqn{\tilde{\mathcal{S}}_j}}
	  }
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
		   }
\author{
		Sergio I. Olivares-Guzman, Adriana Roldan-Rodriguez, Juan Carlos Martinez-Ovando
	   }
\seealso{\code{\link[predfinitepop:domnoplan]{domnoplan.}}}

