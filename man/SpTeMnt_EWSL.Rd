\name{SpTeMnt_EWSL}
\alias{SpTeMnt_EWSL}
\title{Spatio-temporal process monitoring using exponentially weighted spatial LASSO}
\description{Implementation of the online spatio-temporal process monitoring procedure 
described in Qiu and Yang (2022), in which spatial locations with the detected shifts are 
guaranteed to be small clustered spatal regions by the exponentially weighted spatial LASSO.}

\usage{
SpTeMnt_EWSL(y,st,type,ARL0=200,lambda=0.1,B=1000,bs=5,T=1,ht=NULL,hs=NULL,
gt=NULL,gs=NULL)
}

\arguments{
	\item{y}{A vector of \eqn{N} spatio-temporal observations.}
	\item{st}{An \eqn{N\times 3} matrix specifying the spatial locations and times for all the spatio-temporal observations in \code{y}.}
  	\item{type}{A vector of \eqn{N} characters specifying the types of the observations. Here, \code{type} could be \code{IC1}, \code{IC2} or \code{Mnt}, where \code{type='IC1'} denotes the in-control (IC) observations used to perform the block bootstrap procedure to determine the control limit of the CUSUM chart, \code{type='IC2'} denotes the IC observations used to estimate the spatio-temporal mean and covariance functions by \code{SpTe_MeanEst} and \code{SpTe_CovEst}, and \code{type='Mnt'} denotes the observations used for online process monitoring (cf., Yang and Qiu 2020). If there are only data points with either \code{type='IC1'} or \code{type='IC2'}, then these data points 
  		will be used to estimate the model and conduct the bootstrap procedure as well. This function will return an error if there are no observations with \code{type='IC1'} or \code{type='IC2'}.}
  		\item{ARL0}{The pre-specified IC average run length. Default is 200.}
  		\item{lambda}{The pre-specified weighting parameter \eqn{\lambda}. Default is 0.1.}
  	\item{B}{The bootstrap sizes used in the block bootstrap procedure for determining the control limit. Default value is 1,000.}
  		\item{bs}{The block size of the block bootstrap procedure. Default value is 5.}
	\item{T}{The period of the spatio-temporal mean and covariance. Default value is 1.}
		\item{ht}{The temporal kernel bandwidth \code{ht}; default is \code{NULL} and it will be chosen by the modified cross-validation \code{Mod_CV} if \code{ht=NULL}. }
	\item{hs}{The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will be chosen by the function \code{Mod_CV} if \code{hs=NULL}. }
		\item{gt}{The temporal kernel bandwidth \code{gt}; default is \code{NULL} and it will be chosen by minimizing the mean squared prediction error via \code{CV_MSPE} if \code{gt=NULL}. }
	\item{gs}{The spatial kernel bandwidth \code{gs}; default is \code{NULL},  and it will be chosen by the function \code{CV_MSPE} if \code{gs=NULL}.}
}

\value{
\item{ARL0}{Same as the one in the arguments.}
\item{lambda}{Same as the one in the arguments.}
\item{Cstat}{The charting statistics which can be used to make a plot 
for the control chart.}
\item{CL}{The control limit that is determined by the block bootstrap.}
\item{signal_time}{The signal time (i.e., the first time point when the charting statistic \code{Cstat} exceeds the control limit \code{CL}).}
}

\author{Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}

\examples{
library(SpTe2M)
data(ILI)
n <- 365; m <- 67
y <- ILI$Rate; st <- ILI[,3:5]
type <- rep(c('IC1','IC2','Mnt'),c(m*(n+1),(m*n),(m*n)))
IDs <- c(1:(5*m),((n+1)*m+1):(m*(n+6)),((2*n+1)*m+1):(m*(2*n+6)))
y.sub <- y[IDs]; st.sub <- st[IDs,]; type.sub <- type[IDs]
\donttest{ILI.EWSL <- SpTeMnt_EWSL(y.sub,st.sub,type.sub,B=5,ht=0.05,hs=6.5,gt=0.25,gs=1.5)}
}

\references{
Qiu, P. and Yang, K. (2022). Spatio-Temporal Process Monitoring Using Exponentially Weighted Spatial LASSO. \emph{Journal of Quality Technology}, in press.
}

\keyword{SpTeMnt_EWSL}
