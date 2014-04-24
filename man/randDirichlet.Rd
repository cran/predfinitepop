\name{randDirichlet}
\alias{randDirichlet}
\title{Simulates samples from the Dirichlet distribution with (px1) vector parameter \eqn{\alpha}.}
\description{This function simulates samples from the Dirichlet distribution with (px1) vector parameter \eqn{\alpha}.}
\usage{randDirichlet(alpha,n)}
\arguments{
		\item{alpha}{\eqn{p}-dimensional vector with positive entries}
		\item{n}{Number of simulated replicates}
          }        
\value{This function produces:
	\item{randDir }{(\eqn{n} X \eqn{p}) matrix with 'n' simulated replicates}
      }
\references{
	-	"Non-Uniform Random Variate Generation", Luc Devroy (1986), Berlin: Springer-Verlag.
			}
\author{
	Juan Carlos Martinez-Ovando
       }