#' @title 
#' Estimate the spatio-temporal covariance function
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The function \code{spte_covest} is developed to estimate the spatio-temporal 
#' covariance \eqn{V(t,t';s,s')=\mbox{Cov}(y(t,s),y(t',s'))} by the weighted 
#' moment estimation procedure (cf., Yang and Qiu 2019). It should be noted that 
#' the estimated covariance from \code{spte_covest} may not be positive 
#' semidefinite and thus it may not be a legitimate covariance function. In such 
#' cases, the projection-based modification needs to be used to make it positive 
#' semidefinite (cf., Yang and Qiu 2019).
#' 
#' @param y
#' A vector of length \eqn{N} containing data of the observed response.
#' @param st 
#' An \eqn{N \times 3} matrix specifying the spatial locations and times for all 
#' the spatio-temporal observations in \code{y}.
#' @param gt 
#' The temporal kernel bandwidth \code{gt}; default is \code{NULL},  and it will 
#' be chosen by minimizing the mean squared prediction error via \code{cv_mspe} 
#' if \code{gt=NULL}. 
#' @param gs The spatial kernel bandwidth \code{gs}; default is \code{NULL},  
#' and it will be chosen by the function \code{cv_mspe} if \code{gs=NULL}. 
#' @param stE1 An \eqn{N_1 \times 3} matrix specifying the spatial locations 
#' \eqn{s} and times \eqn{t}. Default value is NULL, and \code{stE1=st} if 
#' \code{stE1=NULL}.
#' @param stE2
#' An \eqn{N_2 \times 3} matrix specifying the spatial locations \eqn{s'}
#' and times \eqn{t'}. Default value is NULL, and \code{stE2=st} if 
#' \code{stE2=NULL}.
#' 
#' @return
#' \item{stE1}{Same as the one in the arguments.}
#' \item{stE2}{Same as the one in the arguments.}
#' \item{bandwidth}{The bandwidths \code{(gt, gs)} used in the weighted moment 
#' estimation procedure.}
#' \item{covhat}{An \eqn{N_1 \times N_2} covariance matrix estimate.}
#' 
#' @export
#'
#' @references
#' Yang, K. and Qiu, P. (2019). Nonparametric Estimation of the Spatio-Temporal 
#' Covariance Structure. \emph{Statistics in Medicine}, \strong{38}, 4555-4565.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords spte_covest
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(sim_dat)
#' y <- sim_dat$y; st <- sim_dat$st
#' ids <- 1:500; y.sub <- y[ids]; st.sub <- st[ids,]
#' cov.est <- spte_covest(y.sub,st.sub)

spte_covest <- function(y,st,gt=NULL,gs=NULL,stE1=NULL,stE2=NULL) 
{   
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
    # run the modified cross-validation to select (ht, hs)
    NUM <- dim(st)[1]; mcv <- mod_cv(y,st,eps=0.1)
    ht.opt <- mcv$bandwidth.opt[1]; hs.opt <- mcv$bandwidth.opt[2]
    # call Fortran code to estimate mean by the local linear kernel smoothing
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                       ht=as.double(ht.opt), hs=as.double(hs.opt), stE=as.matrix(st), 
                       NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    muhat <- mu.est$muhat
    # compute the residuals
    res <- y-muhat; RES <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st[,3]==t[i])
        len <- length(ids)
        RES[i,1:len] <- res[ids]
    }
    if(is.null(gt)==TRUE || is.null(gs)==TRUE) {
      # select (gt, gs) by the cross-validation MSPE
       mspe <- cv_mspe(y,st)
       gt <- mspe$bandwidth.opt[1]; gs <- mspe$bandwidth.opt[2]
    }
    if(is.null(stE1)==TRUE) {
       stE1 <- st
    }
    if(is.null(stE2)==TRUE) {
        stE2 <- st
    }
    NUM1 <- dim(stE1)[1]; NUM2 <- dim(stE2)[1]
    # call Fortran code to estimate covariance by the weighted moment estimation
    cov.est <- .Fortran("SpTeWME", res=as.matrix(RES), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), gt=as.double(gt), gs=as.double(gs), 
                       stE1=as.matrix(stE1), NUM1=as.integer(NUM1),
                       stE2=as.matrix(stE2), NUM2=as.integer(NUM2),
                       covhat=matrix(as.double(0),NUM1,NUM2))
    covhat <- cov.est$covhat; results <- list()
    results$stE1 <- stE1; results$stE2 <- stE2
    bandwidth <- matrix(c(gt,gs),2,1); row.names(bandwidth) <- c('gt','gs')
    results$bandwidth <- bandwidth[,1]; results$covhat <- covhat
    return(results)
}
