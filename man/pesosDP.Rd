\name{pesosDP}
\alias{pesosDP}
\title{Computes the weights associated with the predictive distribution \eqn{\hat{G}_j} for a given planned domain \eqn{j} under the Dircihlet processes prior.}
\description{This function computes the weights associated with the predictive distribution \eqn{\hat{G}_j} for a given planned domain \eqn{j}, using the sample ties \eqn{(y^{*}_j)}, under the Dircihlet processes prior.}
											
\usage{
		pesosDP(unicos_sal,alphaDP)
	  }

\arguments{                 
		  \item{unicos_sal}{(\eqn{U} X 2) matrix with sample ties and asociated frequencies (produced with 'unicos')}
		  \item{alphaDP}{Scale parameter of the Dirichlet process}
		  }
\details{
		The matrix with sample ties 'unicos_sal' must include a column named 'm_k' for the frequencies of the sample ties.
		}
\value{The function 'pesosDP' produces an object list with two elements:
                               
	          \item{rho}{(\eqn{U} X 2) dimensional vector with weights associated with 'ystar' (the sample ties)}
              \item{phi}{Probability weight asociated with \eqn{G_{j0}} (the continuos part of \eqn{\hat{G}_j})}
               
      }
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
   }
		   
\author{
	Sergio I. Olivares-Guzman, Juan Carlos Martinez-Ovando
	}