#' @title 
#' Cross-validation mean squared prediction error
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The spatio-temporal covariance function is estimated by the 
#' weighted moment estimation method in Yang and Qiu (2019). The function 
#' \code{cv_mspe} is developed to select the bandwidths \code{(gt,gs)} used 
#' in the estimation of the spatio-temporal covariance function.
#' 
#' @param y 
#' A vector of length \eqn{N} containing data of the observed response 
#' \eqn{y(t,s)}, where \eqn{N} is the total number of observations over space 
#' and time.
#' @param st 
#' An \eqn{N\times 3} matrix specifying the spatial locations 
#' (i.e., (\eqn{s_u},\eqn{s_v})) and times (i.e., \eqn{t}) for all the 
#' observations in \code{y}. The three columns of \code{st} correspond to 
#' \eqn{s_u}, \eqn{s_v} and \eqn{t}, respectively.
#' @param gt 
#' A sequence of temporal kernel bandwidth \code{gt} provided by users; 
#' default is \code{NULL},  and \code{cv_mspe} will choose its own sequence 
#' if \code{gt=NULL}.
#' @param gs 
#' A sequence of spatial kernel bandwidth \code{gs} provided by users; 
#' default is \code{NULL},  and \code{cv_mspe} will choose its own sequence 
#' if \code{gs=NULL}.
#' 
#' @return
#' \item{bandwidth}{A matrix containing all the bandwidths 
#' (\code{gt}, \code{gs}) provided by users.}
#' \item{mspe}{The mean squared prediction errors for all the bandwidths 
#' provided by users.}
#' \item{bandwidth.opt}{The bandwidths \code{(gt, gs)} that minimizes the 
#' mean squared prediction error.} 
#' \item{mspe.opt}{The minimal mean squared prediction error.}
#' 
#' @export
#'
#' @references
#' Yang, K. and Qiu, P. (2019). Nonparametric Estimation of the Spatio-Temporal 
#' Covariance Structure. \emph{Statistics in Medicine}, \strong{38}, 4555-4565.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords cv_mspe
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(sim_dat)
#' y <- sim_dat$y; st <- sim_dat$st
#' gt <- seq(0.3,0.4,0.1); gs <- seq(0.3,0.4,0.1)
#' ids <- 1:500; y.sub <- y[ids]; st.sub <- st[ids,]
#' mspe <- cv_mspe(y.sub,st.sub,gt,gs)

cv_mspe <- function(y,st,gt=NULL,gs=NULL) 
{   # generate the default bandwidth sequence
    if(is.null(gt)==TRUE || is.null(gs)==TRUE) {
        sx.rg <- range(st[,1]); sy.rg <- range(st[,2])
        s.rg <- sqrt((sx.rg[2]-sx.rg[1])^2+(sy.rg[2]-sy.rg[1])^2)
        t.rg <- range(st[,3]); t.rg <- t.rg[2]-t.rg[1]
        gt <- seq(from=t.rg*0.4,to=0.8*t.rg,by=0.4*t.rg/10)
        gs <- seq(from=s.rg*0.4,to=0.8*s.rg,by=0.4*s.rg/10)
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
    Nbw <- length(gt); NUM <- dim(st)[1]
    # run the modified cross-validation to select (ht, hs)
    mcv <- mod_cv(y,st,eps=0.1)
    ht.opt <- mcv$bandwidth.opt[1]; hs.opt <- mcv$bandwidth.opt[2]
    # call Fortran code to estimate mean by the local linear kernel smoothing
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                      sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                      ht=as.double(ht.opt), hs=as.double(hs.opt), stE=as.matrix(st), 
                      NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    muhat <- mu.est$muhat
    res <- y-muhat; RES <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st[,3]==t[i])
        len <- length(ids)
        RES[i,1:len] <- res[ids]
    }
    sx.rg <- range(st[,1]); sy.rg <- range(st[,2])
    s.rg <- sqrt((sx.rg[2]-sx.rg[1])^2+(sy.rg[2]-sy.rg[1])^2)
    t.rg <- range(st[,3]); t.rg <- t.rg[2]-t.rg[1]
    dt <- t.rg/n*3; ds <- 2*s.rg/sqrt(MAXm)
    # call Fortran code to compute the cross-validation MSPE
    cvMSPE <- .Fortran("CVMSPE", res=as.matrix(RES), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), NUM=as.integer(NUM), gt=as.vector(gt), 
                       gs=as.vector(gs), Nbw=as.integer(Nbw), dt=as.double(dt),
                       ds=as.double(ds), mspe=rep(as.double(0),Nbw))
    mspe <- cvMSPE$mspe; results <- list()
    bw <- matrix(0,2,Nbw); bw[1,] <- gt; bw[2,] <- gs 
    row.names(bw) <- c('gt','gs'); results$bandwidth <- bw 
    results$mspe <- mspe; id <- which.min(mspe)
    results$bandwidth.opt <- bw[,id]; results$mspe.opt <- mspe[id]
    return(results)
}
