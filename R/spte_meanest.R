#' @title 
#' Estimate the spatio-temporal mean function
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The function \code{spte_meanest} provides a major tool for estimating the 
#' spatio-temporal mean function nonparametrically 
#' (cf., Yang and Qiu 2018 and 2022).
#' 
#' @param y
#' A vector of spatio-temporal observations.
#' @param st 
#' A three-column matrix specifying the spatial locations and times for all the 
#' spatio-temporal observations in \code{y}.
#' @param ht
#' The temporal kernel bandwidth \code{ht}; default is \code{NULL} and it will 
#' be chosen by the modified cross-validation \code{mod_cv} if \code{ht=NULL}.
#' @param hs
#' The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will 
#' be chosen by the function \code{mod_cv} if \code{hs=NULL}. 
#' @param cor
#' A logical indicator where \code{cor=FALSE} implies that the covariance is not 
#' taken into account and the local linear kernel smoothing procedure is used for 
#' estimating the mean function (cf., Yang and Qiu 2018) and \code{cor=TRUE} 
#' implies that the covariance is accommodated and the three-step local smoothing 
#' approach is used to estimate the mean function (cf., Yang and Qiu 2022). 
#' Default is FALSE.
#' @param stE 
#' A three-column matrix specifying the spatial locations and times where 
#' we want to calculate the estimate of the mean. Default is NULL, and 
#' \code{stE=st} if \code{stE=NULL}.
#' 
#' @return
#' \item{bandwidth}{The bandwidths (\code{ht}, \code{hs}) used in the estimation 
#' procedure.}
#' \item{stE}{Same as the one in the arguments.}
#' \item{muhat}{The estimated mean values at the spatial locations and times 
#' specified by \code{stE}.}
#' 
#' @export
#'
#' @references
#' Yang, K. and Qiu, P. (2018). Spatio-Temporal Incidence Rate Data Analysis by 
#' Nonparametric Regression. Statistics in Medicine, \strong{37}, 2094-2107.
#' @references 
#' Yang, K. and Qiu, P. (2022). A Three-Step Local Smoothing Approach for 
#' Estimating the Mean and Covariance Functions of Spatio-Temporal Data. 
#' \emph{Annals of the Institute of Statistical Mathematics}, \strong{74}, 49-68.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords spte_meanest
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(sim_dat)
#' y <- sim_dat$y; st <- sim_dat$st
#' ids <- 1:500; y.sub <- y[ids]; st.sub <- st[ids,]
#' cov.est <- spte_meanest(y.sub,st.sub)

spte_meanest <- function(y,st,ht=NULL,hs=NULL,cor=FALSE,stE=NULL) 
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
    if(is.null(ht)==TRUE || is.null(hs)==TRUE) {
        # select (ht, hs) by the modified cross-validation
        mcv <- mod_cv(y,st,eps=0.1)
        ht <- mcv$bandwidth.opt[1]; hs <- mcv$bandwidth.opt[2]
    }
    if(is.null(stE)==TRUE) {
        stE <- st
    }
    NUM <- dim(stE)[1]; NUM0 <- dim(st)[1]
    if(cor==FALSE) {
        # call Fortran code to estimate mean by the local linear kernel smoothing
        mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                           sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                           ht=as.double(ht), hs=as.double(hs), stE=as.matrix(stE), 
                           NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    } else {
        # select (gt, gs) by the modified cross-validation
        mspe <- cv_mspe(y,st)
        gt0 <- mspe$bandwidth.opt[1]; gs0 <- mspe$bandwidth.opt[2]
        # call Fortran code to estimate mean by the weighted local smoothing
        mu.est <- .Fortran("SpTeWLS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx), sy=as.matrix(sy), 
                           n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), NUM0=as.integer(NUM0),
                           ht=as.double(ht), hs=as.double(hs), ht0=as.double(ht), hs0=as.double(hs), 
                           gt0=as.double(gt0), gs=as.double(gs0), stE=as.matrix(stE), 
                           NUM=as.integer(NUM), muhat=rep(as.double(0),NUM))  
    }
    results <- list(); bandwidth <- matrix(c(ht,hs),2,1)
    row.names(bandwidth) <- c('ht','hs'); results$bandwidth <- bandwidth[,1]
    results$stE <- stE; results$muhat <- mu.est$muhat
    return(results)
}
