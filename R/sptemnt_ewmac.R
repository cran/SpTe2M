#' @title 
#' Spatio-temporal process monitoring using covariate information
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The function \code{sptemnt_ewmac} is developed to solve the spatio-temporal 
#' process montoring problems in cases when the information in covariates 
#' needs to be used. Please refer to Qiu and Yang (2021) for more details of 
#' the method.
#' 
#' @param y
#' A vector of \eqn{N} spatio-temporal observations.
#' @param x
#' An \eqn{N\times p} matrix containing the data of \eqn{p} covariates.
#' @param st
#' An \eqn{N\times 3} matrix specifying the spatial locations and times for 
#' all the spatio-temporal observations in \code{y}.
#' @param type
#' A vector of \eqn{N} characters specifying the types of the observations. 
#' Here, \code{type} could be \code{IC1}, \code{IC2} or \code{Mnt}, where 
#' \code{type='IC1'} denotes the in-control (IC) observations used to perform 
#' the block bootstrap procedure to determine the control limit of the CUSUM 
#' chart, \code{type='IC2'} denotes the IC observations used to estimate the 
#' spatio-temporal mean and covariance functions by \code{spte_meanest} and 
#' \code{spte_covest}, and \code{type='Mnt'} denotes the observations used 
#' for online process monitoring (cf., Yang and Qiu 2021). If there are only 
#' data points with either \code{type='IC1'} or \code{type='IC2'}, then these 
#' data points will be used to estimate the model and conduct the bootstrap 
#' procedure as well. This function will return an error if there are no 
#' observations with \code{type='IC1'} or \code{type='IC2'}.
#' @param ARL0 
#' The pre-specified IC average run length. Default is 200.
#' @param ARL0.z
#' The pre-specified IC average run length for the covariate chart.
#' Default is 200. Usually, set \code{ARL0.z=ARL0}.
#' @param lambda
#' The pre-specified weighting parameter in the EWMAC chart. Default is 0.1.
#' @param B
#' The bootstrap sizes used in the block bootstrap procedure for determining 
#' the control limit. Default value is 1,000.
#' @param bs 
#' The block size of the block bootstrap procedure. Default value is 5.
#' @param T
#' The period of the spatio-temporal mean and covariance. Default value is 1.
#' @param ht
#' The temporal kernel bandwidth \code{ht}; default is \code{NULL} and it will 
#' be chosen by the modified cross-validation via \code{mod_cv} if \code{ht=NULL}. 
#' @param hs
#' The spatial kernel bandwidth \code{hs}; default is \code{NULL},  and it will 
#' be chosen by the function \code{mod_cv} if \code{hs=NULL}. 
#' @param gt
#' The temporal kernel bandwidth \code{gt}; default is \code{NULL} and it will 
#' be chosen by minimizing the mean squared prediction error via \code{cv_mspe} 
#' if \code{gt=NULL}. 
#' @param gs
#' The spatial kernel bandwidth \code{gs}; default is \code{NULL},  and it will 
#' be chosen by the function \code{cv_mspe} if \code{gs=NULL}. 
#' 
#' @return
#' \item{ARL0}{Same as the one in the arguments.}
#' \item{lambda}{Same as the one in the arguments.}
#' \item{cstat}{The charting statistics which can be used to make 
#' a plot for the control chart.}
#' \item{cl}{The control limit that is determined by the block bootstrap.}
#' \item{signal_time}{The signal time (i.e., the first time point when the 
#' charting statistic \code{cstat} exceeds the control limit \code{cl}).}
#' 
#' @export
#'
#' @references
#' Qiu, P. and Yang, K. (2021). Effective Disease Surveillance by Using Covariate 
#' Information. \emph{Statistics in Medicine}, \strong{40}, 5725-5745. 
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords sptemnt_ewmac
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(ili_dat)
#' n <- 365; m <- 67
#' y <- ili_dat$Rate; x <- as.matrix(ili_dat[,7:8]); st <- ili_dat[,3:5]
#' type <- rep(c('IC1','IC2','Mnt'),c(m*(n+1),(m*n),(m*n)))
#' ids <- c(1:(5*m),((n+1)*m+1):(m*(n+6)),((2*n+1)*m+1):(m*(2*n+6)))
#' y.sub <- y[ids]; x.sub <- x[ids,]; st.sub <- st[ids,]; type.sub <- type[ids]
#' ili.ewmac <- sptemnt_ewmac(y.sub,x.sub,st.sub,type.sub,ht=0.05,hs=6.5,gt=0.25,gs=1.5)

sptemnt_ewmac <- function(y,x,st,type,ARL0=200,ARL0.z=200,lambda=0.1,B=1000,bs=5,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
{   # split the data into 3 parts: IC1, IC2, Mnt
    id1 <- which(type=='IC1'); id2 <- which(type=='IC2')
    id3 <- which(type=='Mnt')
    if(length(id1)!=0 || length(id2)!=0) {
        if(length(id1)==0) {
            id1 <- id2
        }
        if(length(id2)==0) {
            id2 <- id1
        }
    } else {
        stop("IC data must be provided!")
    }
    y.ic1 <- y[id1]; st.ic1 <- st[id1,]
    y.ic2 <- y[id2]; st.ic2 <- st[id2,]
    y <- y[id3]; st <- st[id3,]
    if(is.matrix(x)==TRUE) {
        x.ic1 <- x[id1,]
        x.ic2 <- x[id2,]; x <- x[id3,]
    } else {
        x.ic1 <- x[id1]
        x.ic2 <- x[id2];  x <- x[id3]
    }
    
    for(i in 1:dim(st.ic1)[1]) {
        st.ic1[i,3] <- st.ic1[i,3] %% T
    }
    for(i in 1:dim(st.ic2)[1]) {
        st.ic2[i,3] <- st.ic2[i,3] %% T
    }
    for(i in 1:dim(st)[1]) {
        st[i,3] <- st[i,3] %% T
    }
    
    ## semiparametric modeling
    beta <- spte_semiparmreg(y.ic2,st.ic2,x.ic2,ht,hs,maxIter=100)$beta
    z.ic1 <- x.ic1%*%beta; z.ic2 <- x.ic2%*%beta
    z <- x%*%beta
    
    ## block bootstrap for selecting the bandwidth
    res.y.ic <- spte_decor(y.ic1,st.ic1,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    res.z.ic <- spte_decor(z.ic1,st.ic1,z.ic2,st.ic2,T,ht,hs,gt,gs)$std.res

    t <- sort(unique(st.ic1[,3])); n <- length(t)
    m <- rep(0,n) 
    for(i in 1:n) {
        m[i] <- length(which(st.ic1[,3]==t[i]))
    }
    MAXm <- max(m); RES.y.ic <- RES.z.ic <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st.ic1[,3]==t[i])
        len <- length(ids)
        RES.y.ic[i,1:len] <- res.y.ic[ids]
        RES.z.ic[i,1:len] <- res.z.ic[ids]
    }
    LEN <- 2*bs*n; Ct.z <- matrix(0,B,LEN)
    data.y.b <- array(0,c(B,LEN,MAXm)); m.b <- matrix(0,B,LEN)
    set.seed(12345)
    for(b in 1:B) {
        id <- sample(1:(n-bs+1),2*n,replace=TRUE)
        data.z.b <- matrix(0,LEN,MAXm)
        for(i in 1:(2*n)) {
            for(j in 1:bs) {
                m.b[b,(i-1)*bs+j] <- m[id[i]+j-1]
                data.z.b[(i-1)*bs+j,] <- RES.z.ic[id[i]+j-1,]
                data.y.b[b,(i-1)*bs+j,] <- RES.y.ic[id[i]+j-1,]
            }
        }
        Ct.z[b,1] <- lambda*sum(data.z.b[1,1:m.b[b,1]])/sqrt(m.b[b,1])
        for(i in 2:LEN) {
            Ct.z[b,i] <- lambda*sum(data.z.b[i,1:m.b[b,i]])/sqrt(m.b[b,i])+(1-lambda)*Ct.z[b,i-1]
        }
    }
    hl <- 0; hu <- 100; ARL <- rep(0,B)
    while(abs(mean(ARL)-ARL0.z)>0.001) {
        h <- (hl+hu)/2
        for(b in 1:B) {
            id <- which(Ct.z[b,]>h); id <- c(id,LEN)
            ARL[b] <- min(id)
        }
        if(mean(ARL)>ARL0.z) { hu <- h }
        if(mean(ARL)<ARL0.z) { hl <- h }
        if(abs(hu-hl)<0.001) { break() }
    }
    kappa <- h
    Ct.y <- matrix(0,B,LEN)
    for(b in 1:B) {
        wt <- ifelse(Ct.z[b,1]>kappa,min(1,lambda+Ct.z[b,1]/kappa-1),lambda)
        Ct.y[b,1] <- wt*sum(data.y.b[b,1,1:m.b[b,1]])/sqrt(m.b[b,1])
        for(i in 2:LEN) {
            wt <- ifelse(Ct.z[b,i]>kappa,min(1,lambda+Ct.z[b,i]/kappa-1),lambda)
            Ct.y[b,i] <- wt*sum(data.y.b[b,i,1:m.b[b,i]])/sqrt(m.b[b,i])+(1-wt)*Ct.y[b,i-1]
        }
    }
    # the bi-section algorithm for determining the control limit
    hl <- 0; hu <- 100; ARL <- rep(0,B)
    while(abs(mean(ARL)-ARL0)>0.001) {
        h <- (hl+hu)/2
        for(b in 1:B) {
            id <- which(Ct.y[b,]>h); id <- c(id,LEN)
            ARL[b] <- min(id)
        }
        if(mean(ARL)>ARL0) { hu <- h }
        if(mean(ARL)<ARL0) { hl <- h }
        if(abs(hu-hl)<0.001) { break() }
    }
    CL <- h # the control limit
    
    ## process monitoring
    res.y <- spte_decor(y,st,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    res.z <- spte_decor(z,st,z.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    t <- sort(unique(st[,3])); n <- length(t)
    m <- rep(0,n) 
    for(i in 1:n) {
        m[i] <- length(which(st[,3]==t[i]))
    }
    MAXm <- max(m); RES.y <- RES.z <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st[,3]==t[i])
        len <- length(ids)
        RES.y[i,1:len] <- res.y[ids]
        RES.z[i,1:len] <- res.z[ids]
    }
    # compute the charting statistics
    Ct.y <- Ct.z <- wt <- rep(0,n)
    Ct.z[1] <- lambda*sum(RES.z[1,1:m[1]])/sqrt(m[1])
    wt[1] <- ifelse(Ct.z[1]>kappa,min(1,lambda+Ct.z[1]/kappa-1),lambda)
    Ct.y[1] <- wt[1]*sum(RES.y[1,1:m[1]])/sqrt(m[1])
    for(i in 2:n) {
        Ct.z[i] <- lambda*sum(RES.z[i,1:m[i]])/sqrt(m[i])+(1-lambda)*Ct.z[i-1]
        wt[i] <- ifelse(Ct.z[i]>kappa,min(1,lambda+Ct.z[i]/kappa-1),lambda)
        if(i< n/10) { wt[i] <- min(wt[i],1.5*lambda) }
        Ct.y[i] <- wt[i]*sum(RES.y[i,1:m[i]])/sqrt(m[i])+(1-wt[i])*Ct.y[i-1]
    }                   
    results <- list()
    results$ARLO <- ARL0; results$lambda <- lambda; results$wt <- wt
    results$cstat <- Ct.y; results$cl <- CL
    ids <- which(results$cstat>CL)
    if(length(ids)==0) {
      id <- n+1
    } else {
      id <- min(ids)
    }
    # calculate the signal time
    results$signal_time <- max(1,id-1)
    return(results)
}
