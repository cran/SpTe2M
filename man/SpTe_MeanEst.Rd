\name{SpTe_MeanEst}
\alias{SpTe_MeanEst}
\title{Estimate the spatio-temporal mean function}
\description{The function \code{SpTe_MeanEst} provides a major tool for estimating the spatio-temporal mean function nonparametrically (cf., Yang and Qiu 2018 and 2022).}

\usage{
SpTe_MeanEst(y,st,ht=NULL,hs=NULL,cor=FALSE,stE=NULL)
}

\arguments{
	\item{y}{A vector of spatio-temporal observations.}
	\item{st}{A three-column matrix specifying the spatial locations and times for all the spatio-temporal observations
	in \code{y}.}
	\item{ht}{The temporal kernel bandwidth \code{ht}; default is \code{NULL} and it will be chosen by the modified cross-validation \code{Mod_CV} if \code{ht=NULL}. }
	\item{hs}{The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will be chosen by the function \code{Mod_CV} if \code{hs=NULL}. }
	\item{cor}{A logical indicator where \code{cor=FALSE} implies that the covariance is not taken into account
	and the local linear kernel smoothing 
		procedure is used for estimating the mean function (cf., Yang and Qiu 2018) and \code{cor=TRUE} implies
		that the covariance is accommodated and the three-step local smoothing approach is used to estimate 
		the mean function (cf., Yang and Qiu 2022). Default is FALSE.}
	\item{stE}{A three-column matrix specifying the spatial locations and times where we want to calculate the estimate of the mean. Default is NULL, and \code{stE=st} if \code{stE=NULL}.}
}

\value{
\item{bandwidth}{The bandwidths (\code{ht}, \code{hs}) used in the estimation procedure.}
\item{stE}{Same as the one in the arguments.}
\item{muhat}{The estimated mean values at the spatial locations and times specified by \code{stE}.}
}

\author{ Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}

\examples{
library(SpTe2M)
data(SimDat)
y <- SimDat$y; st <- SimDat$st
IDs <- 1:500; y.sub <- y[IDs]; st.sub <- st[IDs,]
mu.est <- SpTe_MeanEst(y.sub,st.sub) 
}

\references{
Yang, K. and Qiu, P. (2018). Spatio-Temporal Incidence Rate Data Analysis by Nonparametric Regression. Statistics in Medicine, \strong{37}, 2094-2107.

Yang, K. and Qiu, P. (2022). A Three-Step Local Smoothing Approach for Estimating the Mean and Covariance Functions of Spatio-Temporal Data. \emph{Annals of the Institute of Statistical Mathematics}, \strong{74}, 49-68.
}

\keyword{SpTe_MeanEst}
