\name{dataref}
\alias{dataref}
\docType{data}
\title{Reference data set of a simulated population segemented in planned and unplanned domains.}
\description{This data set corresponds to a simulated population of 500 individuals grouped in two planned domains. Each planned domain of the population has been segmented, as well, into the three unplanned domains.  Within the planned domain A, the segmentation is given as domain I (200), domain II (400) and domain III (200), whereas within planned domanin B, the segmentation  is given as domain I (780), domain II (400) and domain III (20). Individual measurements were simulated from log-normal and Weibull distributions with different parameterizations for each combination of planned and unplanned domains. See Martinez-Ovando et al. (2014) for further explanations.}
\usage{data(dataref)}
\format{
  A data frame with 500 observations on the following 3 variables.
  \describe{
    \item{\code{id}}{Identification code for simulated data.}
    \item{\code{y_i}}{Simulated individual ourcomes.}
    \item{\code{domplan}}{Identification for planned domains (domain A and B).}
  }
}
\source{Simulated data. See Martinez-Ovando et al. (2014) for further details.}
\references{
	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
}
\examples{data(dataref)}
\keyword{datasets}