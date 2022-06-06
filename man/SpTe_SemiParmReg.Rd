\name{SpTe_SemiParmReg}
\alias{SpTe_SemiParmReg}
\title{Fit the semiparametric spatio-temporal model}
\description{The function \code{SpTe_SemiParmReg} fits the semiparametric spatio-temporal model to study the relationship between the response \eqn{y} and covariates \eqn{\bm{x}} by the method discussed in Qiu and Yang (2021), in which an iterative algorithm is used to compute the estimated regression coefficients.}

\usage{
SpTe_SemiParmReg(y,st,x,ht=NULL,hs=NULL,maxIter=1000,tol=10^(-4),stE=NULL)
}

\arguments{
	\item{y}{A vector of length \eqn{N} containing the data of the spatio-temporal response \eqn{y(t,s)}.}
	\item{st}{An \eqn{N \times 3} matrix specifying the spatial locations and times for all the spatio-temporal 
	observations in \code{y}.}
	\item{x}{An \eqn{N \times p} matrix containing the data of the \eqn{p} covariates.}
	\item{ht}{The temporal kernel bandwidth \code{ht}; default is \code{NULL} and it will be chosen by the modified cross-validation \code{Mod_CV} if \code{ht=NULL}. }
	\item{hs}{The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will be chosen by the function \code{Mod_CV} if \code{hs=NULL}. }
	\item{maxIter}{A positive integer specifying the maximum number of iterations allowed. Default value is 1,000.}
	\item{tol}{A positive numeric value specifying the tolerance level for the convergence criterion. Default value is 0.0001.}
	\item{stE}{A three-column matrix specifying the spatial locations and times where we want to calculate the estimate of the mean. Default is NULL, and \code{stE=st} if \code{stE=NULL}.}
}

\value{
\item{bandwidth}{The bandwidths (\code{ht}, \code{hs}) used in the estimation procedure.}
\item{stE}{Same as the one in the arguments.}
\item{muhat}{The estimated mean values at spatial locations and times specified by \code{stE}.}
\item{beta}{The vector of the estimated regression coefficient vector.}
}

\author{ Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}

\examples{
library(SpTe2M)
data(SimDat)
y <- SimDat$y; st <- SimDat$st; x <- SimDat$x
IDs <- 1:500; y.sub <- y[IDs]; st.sub <- st[IDs,]; x.sub <- x[IDs]
semi.est <- SpTe_SemiParmReg(y.sub,st.sub,x.sub,maxIter=2)   
}

\references{
Qiu, P. and Yang, K. (2021). Effective Disease Surveillance by Using Covariate Information. \emph{Statistics in Medicine}, \strong{40}, 5725-5745.
}

\keyword{SpTe_SemiParmReg}
