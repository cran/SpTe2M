\name{SpTe_CovEst}
\alias{SpTe_CovEst}
\title{Estimate the spatio-temporal covariance function}
\description{ The function \code{SpTe_CovEst} is developed to estimate the spatio-temporal covariance \eqn{V(t,t';s,s')=\mbox{Cov}(y(t,s),y(t',s'))} by the weighted moment estimation procedure (cf., Yang and Qiu 2019). It should be noted that the estimated covariance from \code{SpTe_CovEst} may
not be positive semidefinite and thus it may not be a legitimate covariance function. In such cases, the projection-based modification needs to be used to make it positive semidefinite (cf., Yang and Qiu 2019).}
\usage{SpTe_CovEst(y,st,gt=NULL,gs=NULL,stE1=NULL,stE2=NULL)}
\arguments{\item{y}{A vector of length \eqn{N} containing data of the observed response.}
	\item{st}{An \eqn{N \times 3} matrix specifying the spatial locations and times for all the spatio-temporal observations in \code{y}.}
	\item{gt}{The temporal kernel bandwidth \code{gt}; default is \code{NULL},  and it will be chosen by minimizing the mean squared 
	prediction error via \code{CV_MSPE} if \code{gt=NULL}. }
	\item{gs}{The spatial kernel bandwidth \code{gs}; default is \code{NULL},  and it will be chosen by the function \code{CV_MSPE} if \code{gs=NULL}. }
	\item{stE1}{An \eqn{N_1 \times 3} matrix specifying the spatial locations \eqn{s} and times \eqn{t}. Default value is NULL, and \code{stE1=st} if \code{stE1=NULL}.}
	\item{stE2}{An \eqn{N_2 \times 3} matrix specifying the spatial locations \eqn{s'}
	and times \eqn{t'}. Default value is NULL, and \code{stE2=st} if \code{stE2=NULL} .}}
\value{\item{stE1}{Same as the one in the arguments.}
\item{stE2}{Same as the one in the arguments.}
\item{bandwidth}{The bandwidths \code{(gt, gs)} used in the weighted moment estimation procedure.}
\item{covhat}{An \eqn{N_1 \times N_2} covariance matrix estimate.}}
\author{ Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}
\examples{library(SpTe2M)
data(SimDat)
y <- SimDat$y; st <- SimDat$st
IDs <- 1:500; y.sub <- y[IDs]; st.sub <- st[IDs,]
cov.est <- SpTe_CovEst(y.sub,st.sub)
}
\references{Yang, K. and Qiu, P. (2019). Nonparametric Estimation of the Spatio-Temporal Covariance Structure. \emph{Statistics in Medicine}, \strong{38}, 4555-4565.}
\keyword{SpTe_CovEst}
