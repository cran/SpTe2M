\name{CV_MSPE}
\alias{CV_MSPE}
\title{Cross-validation mean squared prediction error}
\description{The spatio-temporal covariance function is estimated by the weighted moment estimation method in Yang and Qiu (2019). The function \code{CV_MSPE} is developed to select the bandwidths \code{(gt,gs)} used in the estimation of the spatio-temporal covariance function.}
\usage{CV_MSPE(y,st,gt=NULL,gs=NULL)}
\arguments{\item{y}{A vector of length \eqn{N} containing data of the observed response \eqn{y(t,s)}, where \eqn{N} is the total number of observations over space and time.}
	\item{st}{An \eqn{N\times 3} matrix specifying the spatial locations 
	(i.e., (\eqn{s_u},\eqn{s_v})) and times (i.e., \eqn{t}) for all the observations 
	in \code{y}. The three columns of \code{st} correspond to \eqn{s_u}, \eqn{s_v} and \eqn{t}, respectively.}
	\item{gt}{A sequence of temporal kernel bandwidth \code{gt} provided by users; default is \code{NULL},  and \code{CV_MSPE} will choose its own sequence if \code{gt=NULL}.}
	\item{gs}{A sequence of spatial kernel bandwidth \code{gs} provided by users; default is \code{NULL},  and \code{CV_MSPE} will choose its own sequence if \code{gs=NULL}.}}
\value{\item{bandwidth}{A matrix containing all the bandwidths (\code{gt}, \code{gs}) provided by users.}
\item{mspe}{The mean squared prediction errors for all the bandwidths provided by users.}
\item{bandwidth.opt}{The bandwidths \code{(gt, gs)} that minimizes the mean squared prediction error.}
\item{mspe.opt}{The minimal mean squared prediction error.}}
\author{Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}
\examples{library(SpTe2M)
data(SimDat)
y <- SimDat$y; st <- SimDat$st
gt <- seq(0.3,0.4,0.1); gs <- seq(0.3,0.4,0.1)
IDs <- 1:500; y.sub <- y[IDs]; st.sub <- st[IDs,]
mspe <- CV_MSPE(y.sub,st.sub,gt,gs)
}
\references{
Yang, K. and Qiu, P. (2019). Nonparametric Estimation of the Spatio-Temporal Covariance Structure. \emph{Statistics in Medicine}, \strong{38}, 4555-4565.}
\keyword{CV_MSPE}
