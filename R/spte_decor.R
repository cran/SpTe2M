#' @title 
#' Decorrelate the spatio-temporal data
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The function \code{spte_decor} uses the estimated spatio-temporal mean and 
#' covariance to decorrelate the observed spatio-temporal data. After data 
#' decorrelation, each decorrelated observation should have asymptotic mean of 
#' 0 and asymptotic variance of 1, and the decorrelated data should be 
#' asymptotically uncorrelated with each other.
#' 
#' @param y 
#' A vector of \eqn{N} spatio-temporal observations to decorrelate.
#' @param st
#' A three-column matrix specifying the spatial locations and observation 
#' times of the observations to decorrelate.
#' @param y0
#' A vector of \eqn{N_0} in-control (IC) spatio-temporal observations from 
#' which the IC spatio-temporal mean and covariance functions can be estimated
#' via \code{spte_meanest} and \code{spte_covest}.
#' @param st0
#' A three-column matrix specifying the spatial locations and times for all 
#' the spatio-temporal observations in \code{y0}.
#' @param T 
#' The period of the spatio-temporal mean and covariance. Default value is 1.
#' @param ht
#' The temporal kernel bandwidth \code{ht}; default is \code{NULL}, and it will 
#' be chosen by the modified cross-validation \code{mod_cv} if \code{ht=NULL}. 
#' @param hs 
#' The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will 
#' be chosen by the function \code{mod_cv} if \code{hs=NULL}. 
#' @param gt
#' The temporal kernel bandwidth \code{gt}; default is \code{NULL}, and it will 
#' be chosen by minimizing the mean squared prediction error via \code{cv_mspe} 
#' if \code{gt=NULL}. 
#' @param gs
#' The spatial kernel bandwidth \code{gs}; default is \code{NULL},  and it will 
#' be chosen by the function \code{cv_mspe} if \code{gs=NULL}. 
#' 
#' @return
#' \item{st}{Same as the one in the arguments.}
#' \item{std.res}{The decorrelated data.}
#' 
#' @export
#'
#' @references
#' Yang, K. and Qiu, P. (2020). Online Sequential Monitoring of Spatio-Temporal 
#' Disease Incidence Rates. \emph{IISE Transactions}, \strong{52}, 1218-1233.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords spte_decor
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(sim_dat)
#' y <- sim_dat$y; st <- sim_dat$st
#' ids <- 1:500; y.sub <- y[ids]; st.sub <- st[ids,]
#' decor <- spte_decor(y.sub,st.sub,y0=y.sub,st0=st.sub)

spte_decor <- function(y,st,y0,st0,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
{   
    y.ic <- y0; st.ic <- st0
    st.original <- st
    NUM <- dim(st.ic)[1]
    for(i in 1:NUM) {
        st.ic[i,3] <- st.ic[i,3] %% T
    }
    NUM <- dim(st)[1]
    for(i in 1:NUM) {
        st[i,3] <- st[i,3] %% T
    }
    
    t <- sort(unique(st.ic[,3])); n <- length(t)
    m <- rep(0,n)
    for(i in 1:n) {
        m[i] <- length(which(st.ic[,3]==t[i]))
    }
    MAXm <- max(m); Y <- sx <- sy <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st.ic[,3]==t[i])
        len <- length(ids)
        Y[i,1:len] <- y.ic[ids]
        sx[i,1:len] <- st.ic[ids,1]
        sy[i,1:len] <- st.ic[ids,2]
    }
    NUM <- dim(st.ic)[1]
    if(is.null(ht)==TRUE || is.null(hs)==TRUE) {
        # use the modified cross-validation to select (ht, hs)
        mcv <- mod_cv(y.ic,st.ic,eps=0.1)
        ht <- mcv$bandwidth.opt[1]; hs <- mcv$bandwidth.opt[2]
    } 
    # call Fortran code to estimate mean by the local linear kernel smoothing
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                       ht=as.double(ht), hs=as.double(hs), stE=as.matrix(st.ic), 
                       NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    muhat <- mu.est$muhat
    # compute the residuals
    res <- y.ic-muhat; RES <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st.ic[,3]==t[i])
        len <- length(ids)
        RES[i,1:len] <- res[ids]
    }
    if(is.null(gt)==TRUE || is.null(gs)==TRUE) {
      # use the cross-validation MSPE to select (gt, gs)
       mspe <- cv_mspe(y.ic,st.ic)
       gt <- mspe$bandwidth.opt[1]; gs <- mspe$bandwidth.opt[2]
    }
    
    NUM <- dim(st)[1]
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                       ht=as.double(ht), hs=as.double(hs), stE=as.matrix(st), 
                       NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    muhat <- mu.est$muhat
    # estimate the covariance by the weighted moment estimation
    cov.est <- .Fortran("SpTeWME", res=as.matrix(RES), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), gt=as.double(gt), gs=as.double(gs), 
                       stE1=as.matrix(st), NUM1=as.integer(NUM),
                       stE2=as.matrix(st), NUM2=as.integer(NUM),
                       covhat=matrix(as.double(0),NUM,NUM))
    covhat <- cov.est$covhat
    # start the data decorrelation procedure
    res <- y-muhat
    t <- sort(unique(st[,3])); n <- length(t)
    m <- rep(0,n)
    for(i in 1:n) {
        m[i] <- length(which(st[,3]==t[i]))
    }
    MAXm <- max(m); res.std <- rep(0,NUM)
    # decorrelate the data sequentially
    # Time=1
    Sig <- covhat[1:m[1],1:m[1]]; eig <- eigen(Sig) 
    M <- eig$vectors%*%diag(eig$values^(-0.5))%*%t(eig$vectors)
    res.std[1:m[1]] <- M%*%res[1:m[1]]
    # Time>1
    resid <- c(res[1:m[1]]); SIGMA <- Sig; num <- m[1]
    for(i in 2:n) {
        num.old <- sum(m[1:(i-1)]); num.new <- sum(m[1:i])
        Sigma <- covhat[(1+num.old-num):num.old,(1+num.old):num.new]
        SIG <- covhat[(1+num.old):num.new,(1+num.old):num.new]
        Sig <- SIG-t(Sigma)%*%solve(SIGMA)%*%Sigma
        M <- eig$vectors%*%diag(eig$values^(-0.5))%*%t(eig$vectors)
        res.std[(1+num.old):num.new] <- M%*%(res[(1+num.old):num.new]-t(Sigma)%*%solve(SIGMA)%*%resid)
        ## update
        resid <- c(resid,res[(1+num.old):num.new])
        SIGMA <- rbind(cbind(SIGMA,Sigma), cbind(t(Sigma),SIG))
        num <- num+m[i]
        if(i>5) {
            resid <- resid[-(1:m[i-5])]
            SIGMA <- SIGMA[(m[i-5]+1):dim(SIGMA)[1],(m[i-5]+1):dim(SIGMA)[1]]
            num <- num-m[i-5]
        }
    }
    results <- list(); results$st <- st.original; results$std.res <- res.std
    return(results)
}
