\name{conteo}
\alias{conteo}
\title{Counts the number of individuals in \eqn{\mathcal{S}_j} and \eqn{\widetilde{\mathcal{S}}_j}.}
\description{This function counts the number of individuals in \eqn{\mathcal{S}_j} and \eqn{\widetilde{\mathcal{S}}_j}, with respect to the total number of individuals in a given planned domain \eqn{\mathcal{P}_j}.}
\usage{conteo(datos_j, N_j)}
\arguments{
		  \item{datos_j}{Data matrix with individual measurements for \eqn{\mathcal{S}_j}}
		  \item{N_j}{Number of individuals in \eqn{\mathcal{P}_j} (the \eqn{jth} planned domain) (\eqn{N_j} is assumed known)}
		  }          
\details{ 
            "datos_j"       -	This object should contain two columns labeled "\eqn{n_i}" and "domplan".  
			
								"domplan" - Column vector with categories for the planned domains.
								
								"n_i"- Number of individuals in the \eqn{ith} group of individuals, and 
								
             "N_j"  		-	This object should include two colums labeled "\eqn{N_j}" and "domplan".

					        	"domplan" - Column vector with categories for the planned domains.

						        "N_j" - Number of individuals in the \eqn{jth} planned domain, and 
		}

\value{The function 'conteo' produces an output list with three elements:
  \item{M}{Number of groups in "\eqn{\mathcal{S}_j}"} 
  \item{CardS}{Number of individuals in "\eqn{\mathcal{S}_j}" (if "\eqn{n_i}"= 1, for any "\eqn{j}", then \eqn{Cards=M=nrow(datos_j))}}  
  \item{CardNoS}{Number of individuals out of sample}
}
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.
}
\author{
	Sergio I. Olivares-Guzman, Adriana Roldan-Rodriguez, Juan Carlos Martinez-Ovando
		}

