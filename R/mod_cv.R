#' @title 
#' Modifed cross-validation for bandwidth selection
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The spatio-temporal mean function can be estimated by the local linear kernel 
#' smoothing procedure (cf., Yang and Qiu 2018). The function \code{mod_cv} 
#' provides a reliable tool for selecting bandwidths \code{(ht, hs)} used in 
#' the local linear kernel smoothing procedure in cases when data are 
#' spatio-temporally correlated.
#' 
#' @param y 
#' A vector of the spatio-temporal response \eqn{y(t,s)}.
#' @param st 
#' A three-column matrix specifying the spatial locations and times for all 
#' the spatio-temporal observations in \code{y}.
#' @param ht
#' A sequence of temporal kernel bandwidth \code{ht} provided by users; default 
#' is \code{NULL}, and \code{mod_cv} chooses its own sequence if \code{ht=NULL}. 
#' @param hs
#' A sequence of temporal kernel bandwidth \code{hs} provided by users; default 
#' is \code{NULL}, and \code{mod_cv} chooses its own sequence if \code{hs=NULL}. 
#' @param eps
#' The value of this parametric is between 0 and 1. Default is 0.1. The 
#' following bimodal kernel function (cf., Yang and Qiu 2018) is used when 
#' calculting the modified cross-validation score:
#' \deqn{K_{\epsilon}(x) = \frac{4}{4-3\epsilon-\epsilon^3}
#' \left\{ \begin{array}{ll}
#' \frac{3}{4}(1-x^2)\mbox{I}(|x|\leq 1), & \mbox{ if } |x| \geq \epsilon, \\
#' \frac{3(1-\epsilon^2)}{4\epsilon}|x|, & \mbox{ otherwise}.
#' \end{array}
#' \right.}
#' The argument \code{eps} represents the parameter \eqn{\epsilon} in the above 
#' bimodal kernel, which controls the closeness of the bimodal kernel to the 
#' Epanechnikov kernel \eqn{K_e(x)=0.75(1-x^2)\mbox{I}(|x|\leq 1)}. The smaller 
#' the value, the closer the two kernels.
#' 
#' @return
#' \item{bandwidth}{A matrix containing all the bandwidths (\code{ht}, \code{hs}) 
#' provided by users.}
#' \item{mcv}{The modified cross-validation scores for all the bandwidths 
#' provided by users.}
#' \item{bandwidth.opt}{The selected bandwidths \code{(ht, hs)} by the modified 
#' cross-validation.}
#' \item{mcv.opt}{The modified cross-validation score of the selected 
#' bandwidths.}
#' 
#' @export
#'
#' @references
#' Yang, K. and Qiu, P. (2018). Spatio-Temporal Incidence Rate Data Analysis by 
#' Nonparametric Regression. \emph{Statistics in Medicine}, \strong{37}, 2094-2107.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords mod_cv
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(sim_dat)
#' y <- sim_dat$y; st <- sim_dat$st
#' ht <- seq(0.10,0.15,0.05); hs <- seq(0.20,0.30,0.10)
#' ids <- 1:500; y.sub <- y[ids]; st.sub <- st[ids,]
#' mcv <- mod_cv(y.sub,st.sub,ht,hs,eps=0.1)

mod_cv <- function(y,st,ht=NULL,hs=NULL,eps=0.1) 
{   # generate the default bandwidth sequence
    if(is.null(ht)==TRUE || is.null(hs)==TRUE) {
        sx.rg <- range(st[,1]); sy.rg <- range(st[,2])
        s.rg <- sqrt((sx.rg[2]-sx.rg[1])^2+(sy.rg[2]-sy.rg[1])^2)
        t.rg <- range(st[,3]); t.rg <- t.rg[2]-t.rg[1]
        ht <- seq(from=0.1*t.rg,to=0.6*t.rg,by=0.5*t.rg/9)
        hs <- seq(from=0.2*s.rg,to=0.8*s.rg,by=0.6*s.rg/9)
    }
    t <- sort(unique(st[,3])); n <- length(t)
    m <- rep(0,n)
    for(i in 1:n) {
        m[i] <- length(which(st[,3]==t[i]))
    }
    MAXm <- max(m); Y <- sx <- sy <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st[,3]==t[i])
        len <- length(ids)
        Y[i,1:len] <- y[ids]
        sx[i,1:len] <- st[ids,1]
        sy[i,1:len] <- st[ids,2]
    }
    Nbw <- length(ht); NUM <- dim(st)[1]
    # call Fortran code to perform the modified cross-validation
    modCV <- .Fortran("ModCV", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), NUM=as.integer(NUM), ht=as.vector(ht), 
                       hs=as.vector(hs), Nbw=as.integer(Nbw), eps=as.double(eps),
                       mcv=rep(as.double(0),Nbw))
    mcv <- modCV$mcv; results <- list()
    bw <- matrix(0,2,Nbw); bw[1,] <- ht; bw[2,] <- hs 
    row.names(bw) <- c('ht','hs'); results$bandwidth <- bw 
    results$mcv <- mcv; id <- which.min(mcv)
    results$bandwidth.opt <- bw[,id]; results$mcv.opt <- mcv[id]
    return(results)
}
