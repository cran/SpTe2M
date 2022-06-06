\name{Mod_CV}
\alias{Mod_CV}
\title{Modifed cross-validation for bandwidth selection}
\description{The spatio-temporal mean function can be estimated by the local linear 
kernel smoothing procedure (cf., Yang and Qiu 2018). The function \code{Mod_CV} provides a 
	reliable tool for selecting bandwidths \code{(ht, hs)} used in the local linear kernel
	smoothing procedure in cases when data are spatio-temporally correlated.}
	\usage{Mod_CV(y,st,ht=NULL,hs=NULL,eps=0.1)}
	\arguments{\item{y}{A vector of the spatio-temporal response \eqn{y(t,s)}.}
	\item{st}{A three-column matrix specifying the spatial locations and times for
	all the spatio-temporal observations in \code{y}.}
	\item{ht}{A sequence of temporal kernel bandwidth \code{ht} provided by users; default is \code{NULL}, and \code{Mod_CV} chooses its own sequence if \code{ht=NULL}. }
	\item{hs}{A sequence of temporal kernel bandwidth \code{hs} provided by users; default is \code{NULL}, and \code{Mod_CV} chooses its own sequence if \code{hs=NULL}. }
	\item{eps}{The value of this parametric is between 0 and 1. Default is 0.1.  The 
	following bimodal kernel function (cf., Yang and Qiu 2018) is used when calculting the modified 
	cross-validation score:
	\deqn{K_{\epsilon}(x) = \frac{4}{4-3\epsilon-\epsilon^3}
   \left\{ \begin{array}{ll}
           \frac{3}{4}(1-x^2)\mbox{I}(|x|\leq 1), & \mbox{ if } |x| \geq \epsilon, \\
           \frac{3(1-\epsilon^2)}{4\epsilon}|x|, & \mbox{ otherwise}.
           \end{array}
   \right.}
		The argument \code{eps} represents the parameter \eqn{\epsilon} in the above bimodal kernel, which controls the closeness of the bimodal kernel to the Epanechnikov kernel \eqn{K_e(x)=0.75(1-x^2)\mbox{I}(|x|\leq 1)}. The smaller the value, the closer the two kernels.}}
	\value{\item{bandwidth}{A matrix containing all the bandwidths (\code{ht}, \code{hs}) provided by users.}
\item{mcv}{The modified cross-validation scores for all the bandwidths 
provided by users.}
\item{bandwidth.opt}{The selected bandwidths \code{(ht, hs)} by the modified cross-validation.}
\item{mcv.opt}{The modified cross-validation score of the selected bandwidths.}}
\author{Kai Yang \email{kayang@mcw.edu} and Peihua Qiu}
\examples{library(SpTe2M)
data(SimDat)
y <- SimDat$y; st <- SimDat$st
ht <- seq(0.10,0.15,0.05); hs <- seq(0.20,0.30,0.10)
IDs <- 1:500; y.sub <- y[IDs]; st.sub <- st[IDs,]
mcv <- Mod_CV(y.sub,st.sub,ht,hs,eps=0.1)
}
\references{Yang, K. and Qiu, P. (2018). Spatio-Temporal Incidence Rate Data Analysis by Nonparametric Regression. \emph{Statistics in Medicine}, \strong{37}, 2094-2107.}
\keyword{Mod_CV}
