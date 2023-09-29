#' @title 
#' Fit the semiparametric spatio-temporal model
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The function \code{spte_semiparmreg} fits the semiparametric spatio-temporal 
#' model to study the relationship between the response \eqn{y} and covariates 
#' \eqn{\bm{x}} by the method discussed in Qiu and Yang (2021), in which an 
#' iterative algorithm is used to compute the estimated regression coefficients.
#' 
#' @param y
#' A vector of length \eqn{N} containing the data of the spatio-temporal 
#' response \eqn{y(t,s)}.
#' @param st
#' An \eqn{N \times 3} matrix specifying the spatial locations and times for all 
#' the spatio-temporal observations in \code{y}.
#' @param x
#' An \eqn{N \times p} matrix containing the data of the \eqn{p} covariates.
#' @param ht
#' The temporal kernel bandwidth \code{ht}; default is \code{NULL} and it will 
#' be chosen by the modified cross-validation \code{mod_cv} if \code{ht=NULL}. 
#' @param hs 
#' The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will 
#' be chosen by the function \code{mod_cv} if \code{hs=NULL}. 
#' @param maxIter
#' A positive integer specifying the maximum number of iterations allowed. 
#' Default value is 1,000.
#' @param tol
#' A positive numeric value specifying the tolerance level for the convergence 
#' criterion. Default value is 0.0001.
#' @param stE
#' A three-column matrix specifying the spatial locations and times where we 
#' want to calculate the estimate of the mean. Default is NULL, and 
#' \code{stE=st} if \code{stE=NULL}.
#' 
#' @return
#' \item{bandwidth}{The bandwidths (\code{ht}, \code{hs}) used in the estimation 
#' procedure.}
#' \item{stE}{Same as the one in the arguments.}
#' \item{muhat}{The estimated mean values at spatial locations and times 
#' specified by \code{stE}.}
#' \item{beta}{The vector of the estimated regression coefficient vector.}
#' 
#' @export
#'
#' @references
#' Qiu, P. and Yang, K. (2021). Effective Disease Surveillance by Using 
#' Covariate Information. \emph{Statistics in Medicine}, \strong{40}, 5725-5745.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords spte_semiparmreg
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(sim_dat)
#' y <- sim_dat$y; st <- sim_dat$st; x <- sim_dat$x
#' ids <- 1:500; y.sub <- y[ids]; st.sub <- st[ids,]; x.sub <- x[ids]
#' semi.est <- spte_semiparmreg(y.sub,st.sub,x.sub,maxIter=2)   

spte_semiparmreg <- function(y,st,x,ht=NULL,hs=NULL,maxIter=1000,tol=10^(-4),stE=NULL) 
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
    NUM <- dim(st)[1]
    if(is.null(ht)==TRUE || is.null(hs)==TRUE) {
        # select (ht, hs) by the modified cross-validation
        mcv <- mod_cv(y,st,eps=0.1)
        ht <- mcv$bandwidth.opt[1]; hs <- mcv$bandwidth.opt[2]
    }
    if(is.null(stE)==TRUE) {
        stE <- st
    }
    if(is.matrix(x)==TRUE) {
        p <- dim(x)[2] 
    } else {
        p <- 1
    }
    # implement the iterative estimation procedure
    k <- 0; beta.old <- rep(1,p); beta.new <- rep(0,p)
    while(k<maxIter & sum(abs(beta.old-beta.new))/sum(abs(beta.old))>tol) {
        beta.old <- beta.new
        if(p>1) {
            z <- y-x%*%beta.old
        } else {
            z <- y-x*beta.old
        }
        Z <- matrix(0,n,MAXm)
        for(i in 1:n) {
            ids <- which(st[,3]==t[i])
            len <- length(ids)
            Z[i,1:len] <- z[ids]
        }
        # call Fortran code to implement the local linear kernel smoothing
        mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Z), t=as.vector(t), sx=as.matrix(sx),
                           sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                           MAXm=as.integer(MAXm), ht=as.vector(ht), hs=as.vector(hs), 
                           stE=as.matrix(st), NUM=as.integer(NUM), eps=as.double(0),
                           muhat=rep(as.double(0),NUM))
        muhat <- mu.est$muhat; z.new <- y-muhat
        beta.new <- as.numeric(coef(lm(z.new~x-1)))
        k <- k+1
    }
    NUM1 <- dim(stE)[1]
    # use the local linear kernel smoothing to compute the final estimate
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Z), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), ht=as.vector(ht), hs=as.vector(hs), 
                       stE=as.matrix(stE), NUM=as.integer(NUM1), eps=as.double(0),
                       muhat=rep(as.double(0),NUM1))
    results <- list(); bandwidth <- matrix(c(ht,hs),2,1)
    row.names(bandwidth) <- c('ht','hs'); results$bandwidth <- bandwidth[,1]
    results$stE <- stE; results$muhat <- mu.est$muhat
    results$beta <- beta.new
    return(results)
}
