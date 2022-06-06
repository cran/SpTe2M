\name{SpTe_DeCor}
\alias{SpTe_DeCor}
\title{Decorrelate the spatio-temporal data}
\description{The function \code{SpTe_DeCor} uses the estimated spatio-temporal mean and covariance to decorrelate the 
observed spatio-temporal data. After data decorrelation, each decorrelated observation should
have asymptotic mean of 0 and asymptotic variance of 1, and the decorrelated data should be 
asymptotically uncorrelated with each other.}

\usage{SpTe_DeCor(y,st,y0,st0,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL)}

\arguments{
	\item{y}{A vector of \eqn{N} spatio-temporal observations to decorrelate.}
	\item{st}{A three-column matrix specifying the spatial locations and observation times of the observations
	to decorrelate.}
   \item{y0}{A vector of \eqn{N_0} in-control (IC) spatio-temporal observations from which the IC spatio-temporal mean and 
   covariance functions can be estimated
   	via \code{SpTe_MeanEst} and \code{SpTe_CovEst}.}
	\item{st0}{A three-column matrix specifying the spatial locations and times for all the spatio-temporal observations in \code{y0}.}
	\item{T}{The period of the spatio-temporal mean and covariance. Default value is 1.}
		\item{ht}{The temporal kernel bandwidth \code{ht}; default is \code{NULL}, and it will be chosen by the modified cross-validation \code{Mod_CV} if \code{ht=NULL}. }
	\item{hs}{The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will be chosen by the function \code{Mod_CV} if \code{hs=NULL}. }
		\item{gt}{The temporal kernel bandwidth \code{gt}; default is \code{NULL}, and it will be chosen by minimizing the mean squared prediction error via \code{CV_MSPE} if \code{gt=NULL}. }
	\item{gs}{The spatial kernel bandwidth \code{gs}; default is \code{NULL},  and it will be chosen by the function \code{CV_MSPE} if \code{gs=NULL}. }
}

\value{
\item{st}{Same as the one in the arguments.}
\item{std.res}{The decorrelated data.}
}

\author{ Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}

\examples{
library(SpTe2M)
data(SimDat)
y <- SimDat$y; st <- SimDat$st
IDs <- 1:500; y.sub <- y[IDs]; st.sub <- st[IDs,]
decor <- SpTe_DeCor(y.sub,st.sub,y0=y.sub,st0=st.sub)
}

\references{
Yang, K. and Qiu, P. (2020). Online Sequential Monitoring of Spatio-Temporal Disease Incidence Rates. \emph{IISE Transactions}, \strong{52}, 1218-1233.
}

\keyword{SpTe_DeCor}
