#' @title 
#' Online spatio-temporal process monitoring by a CUSUM chart
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' The function \code{sptemnt_cusum} implements the sequential online monitoring 
#' procedure described in Yang and Qiu (2020).
#' 
#' @param y
#' A vector of \eqn{N} spatio-temporal observations.
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
#' for online process monitoring (cf., Yang and Qiu 2020). If there are only 
#' data points with either \code{type='IC1'} or \code{type='IC2'}, then these 
#' data points will be used to estimate the model and conduct the bootstrap 
#' procedure as well. This function will return an error if there are no 
#' observations with \code{type='IC1'} or \code{type='IC2'}.
#' @param ARL0 
#' The pre-specified IC average run length. Default is 200.
#' @param gamma 
#' The pre-specified allowance constant in the CUSUM chart. Default is 0.1.
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
#' \item{gamma}{Same as the one in the arguments.}
#' \item{cstat}{The charting statistics which can be used to make 
#' a plot for the control chart.}
#' \item{cl}{The control limit that is determined by the block bootstrap.}
#' \item{signal_time}{The signal time (i.e., the first time point when the 
#' charting statistic \code{cstat} exceeds the control limit \code{cl}).}
#' 
#' @export
#'
#' @references
#' Yang, K. and Qiu, P. (2020). Online Sequential Monitoring of Spatio-Temporal 
#' Disease Incidence Rates. \emph{IISE Transactions}, \strong{52}, 1218-1233.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords sptemnt_cusum
#' 
#' @import glmnet MASS ggplot2 maps mapproj knitr rmarkdown
#' @importFrom stats coef lm median
#' 
#' @examples
#' library(SpTe2M)
#' data(ili_dat)
#' n <- 365; m <- 67
#' y <- ili_dat$Rate; st <- ili_dat[,3:5]
#' type <- rep(c('IC1','IC2','Mnt'),c(m*(n+1),(m*n),(m*n)))
#' ids <- c(1:(5*m),((n+1)*m+1):(m*(n+6)),((2*n+1)*m+1):(m*(2*n+6)))
#' y.sub <- y[ids]; st.sub <- st[ids,]; type.sub <- type[ids]
#' ili.cusum <- sptemnt_cusum(y.sub,st.sub,type.sub,ht=0.05,hs=6.5,gt=0.25,gs=1.5)

sptemnt_cusum <- function(y,st,type,ARL0=200,gamma=0.1,B=1000,bs=5,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
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
    
    for(i in 1:dim(st.ic1)[1]) {
        st.ic1[i,3] <- st.ic1[i,3] %% T
    }
    for(i in 1:dim(st.ic2)[1]) {
        st.ic2[i,3] <- st.ic2[i,3] %% T
    }
    for(i in 1:dim(st)[1]) {
        st[i,3] <- st[i,3] %% T
    }
    
    ## block bootstrap
    res.ic <- spte_decor(y.ic1,st.ic1,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    t <- sort(unique(st.ic1[,3])); n <- length(t)
    m <- rep(0,n) 
    for(i in 1:n) {
        m[i] <- length(which(st.ic1[,3]==t[i]))
    }
    MAXm <- max(m); RES.ic <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st.ic1[,3]==t[i])
        len <- length(ids)
        RES.ic[i,1:len] <- res.ic[ids]
    }
    set.seed(12345)
    LEN <- 2*bs*n; Ct <- matrix(0,B,LEN)
    for(b in 1:B) {
        id <- sample(1:(n-bs+1),2*n,replace=TRUE)
        data.b <- matrix(0,LEN,MAXm); m.b <- rep(0,LEN)
        for(i in 1:(2*n)) {
            for(j in 1:bs) {
                m.b[(i-1)*bs+j] <- m[id[i]+j-1]
                data.b[(i-1)*bs+j,] <- RES.ic[id[i]+j-1,]
            }
        }
        Ct[b,1] <- max(0,(sum(data.b[1,1:m.b[1]]^2)-m.b[1])/sqrt(2*m.b[1])-gamma)
        for(i in 2:LEN) {
            Ct[b,i] <- max(0,Ct[b,i-1]+(sum(data.b[i,1:m.b[i]]^2)-m.b[1])/sqrt(2*m.b[i])-gamma)
        }
    }
    # the bi-section search algorithm
    hl <- 0; hu <- 100; ARL <- rep(0,B)
    while(abs(mean(ARL)-ARL0)>0.001) {
        h <- (hl+hu)/2
        for(b in 1:B) {
            id <- which(Ct[b,]>h); id <- c(id,LEN)
            ARL[b] <- min(id)
        }
        if(mean(ARL)>ARL0) { hu <- h }
        if(mean(ARL)<ARL0) { hl <- h }
        if(abs(hu-hl)<0.001) { break() }
    }
    CL <- h # the control limit
    
    ## process monitoring
    res <- spte_decor(y,st,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    t <- sort(unique(st[,3])); n <- length(t)
    m <- rep(0,n) 
    for(i in 1:n) {
        m[i] <- length(which(st[,3]==t[i]))
    }
    MAXm <- max(m); RES <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st[,3]==t[i])
        len <- length(ids)
        RES[i,1:len] <- res[ids]
    }
    # compute the charting statistics
    Cstat <- rep(0,n)
    Cstat[1] <- max(0,(sum(RES[1,1:m[1]]^2)-m[1])/sqrt(2*m[1])-gamma)
    for(i in 2:n) {
        Cstat[i] <- max(0,Cstat[i-1]+(sum(RES[i,1:m[i]]^2)-m[i])/sqrt(2*m[i])-gamma)
    }
    results <- list()
    results$ARLO <- ARL0; results$gamma <- gamma
    results$cstat <- Cstat; results$cl <- CL
    ids <- which(results$cstat>CL)
    if(length(ids)==0) {
      id <- n+1
    } else {
      id <- min(ids)
    }
    # compute the signal time
    results$signal_time <- max(1,id-1)
    return(results)
}
