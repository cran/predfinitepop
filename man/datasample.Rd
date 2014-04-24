\name{datasample}
\alias{datasample}
\docType{data}
\title{
	Data set of a simulated population segemented in planned and unplanned domains.
	}
\description{
	This is an object list with 20 simulated samples from the simulated population of 2,000 individuals 
	grouped in two planned domains. Planned domain A is formed of 800 individuals, and domain B is composed 
	of 1.2 thousand individuals. Each planned domain of the population has been segmented, as well, into 
	the three unplanned domains. 
	}
\usage{
	data(datasample)
	}
\details{
	Within the planned domain A, the segmentation is given as domain I (200), 
	domain II (400) and domain III (200), whereas within planned domanin B, the segmentation  is given as 
	domain I (780), domain II (400) and domain III (20). Individual measurements were simulated from log-normal 
	and Weibull distributions with different parameterizations for each combination of planned and unplanned domains. 
	The original sample is of size 5 percent. The upcoming samples were extrated increasing the sample size in 5 percent at a time. 
	See Martinez-Ovando et al. (2014) for further explanations.
	}
\format{
	List of lenght 20, with entries:
	\describe{
		\item{l1.}{'data.frame':	100 obs. of  8 variables with a sample of size 5 percent}
		\item{l2.}{'data.frame':	200 obs. of  8 variables with a sample of size 10 percent}
		\item{l3.}{'data.frame':	300 obs. of  8 variables with a sample of size 15 percent}
		\item{l4.}{'data.frame':	400 obs. of  8 variables with a sample of size 20 percent}
		\item{l5.}{'data.frame':	500 obs. of  8 variables with a sample of size 25 percent}
		\item{l6.}{'data.frame':	600 obs. of  8 variables with a sample of size 30 percent}
		\item{l7.}{'data.frame':	700 obs. of  8 variables with a sample of size 35 percent}
		\item{l8.}{'data.frame':	800 obs. of  8 variables with a sample of size 40 percent}
		\item{l9.}{'data.frame':	900 obs. of  8 variables with a sample of size 45 percent}
		\item{l10.}{'data.frame':	1000 obs. of  8 variables with a sample of size 50 percent}
		\item{l11.}{'data.frame':	1100 obs. of  8 variables with a sample of size 55 percent}
		\item{l12.}{'data.frame':	1200 obs. of  8 variables with a sample of size 60 percent}
		\item{l13.}{'data.frame':	1300 obs. of  8 variables with a sample of size 65 percent}
		\item{l14.}{'data.frame':	1400 obs. of  8 variables with a sample of size 70 percent}
		\item{l15.}{'data.frame':	1500 obs. of  8 variables with a sample of size 75 percent}
		\item{l16.}{'data.frame':	1600 obs. of  8 variables with a sample of size 80 percent}
		\item{l17.}{'data.frame':	1700 obs. of  8 variables with a sample of size 85 percent}
		\item{l18.}{'data.frame':	1800 obs. of  8 variables with a sample of size 90 percent}
		\item{l19.}{'data.frame':	1900 obs. of  8 variables with a sample of size 95 percent}
		\item{l20.}{'data.frame':	2000 obs. of  8 variables with a sample of size 100 percent}
		}
	List of variables:
	\describe{
		\item{id}{Individual identification code}
		\item{y_i}{Group measurenment}
		\item{domnoplan}{Identification code for unplanned domains}
		\item{n_i}{Group size}
		\item{domplan}{Identification code for planned domains}
		\item{1}{Indicator variable for unplanned domain I}
		\item{2}{Indicator variable for unplanned domain II}
		\item{3}{Indicator variable for unplanned domain III}
		}
	}
\source{
	Simulated data. See Martinez-Ovando et al. (2014) for further details.
	}
\references{
	- "Predictive inference on finite populations segmented in planned and unplanned domains," Martinez-Ovando, J. C., Olivares-Guzman, S. I., Roldan-Rodriguez, A., 2014. Submitted. Also available as Banco de Mexico Working Paper 2014-04. 	
	}
\examples{data(datasample)}
\keyword{datasets}