\name{unicstar}
\alias{unicstar}
\title{Identifies the sample ties in \eqn{\mathcal{S}_j} and computes their associated sample frequencies.}
\description{This function identifies the sample ties \eqn{(y^{*}_k)} in \eqn{\mathcal{S}_j} and computes their associated sample frequencies.}
\usage{
		unicstar(datos_j)
	  }

\arguments{
		\item{datos_j}{Data matrix with individual measurements in "\eqn{\mathcal{S}_j}"}
		  }

\details{The column with individual measurements in "datos_j" should be named "y_i"}
         
\value{The function 'unicstar' produces an object list with two elements:
	\item{U}{Number of sample ties in "\eqn{\mathcal{S}_j}"}
	\item{unicstar}{Data matrix with sample ties '\eqn{y_i}' and asociated frequencies }
      }
\references{
    - "A Bayesian nonparametric framework to inference on totals of finite populations," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2013. Contributions of Young Researchers to Bayesian Statistics (eds. E. Lazarone & F. Ieva), Chapter 15 Springer.

	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
   }
\author{
	Sergio I. Olivares-Guzman, Adriana Roldan-Rodriguez, Juan Carlos Martinez-Ovando
       }

