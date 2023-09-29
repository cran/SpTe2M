#' @title 
#' Spatio-temporal process monitoring using exponentially weighted spatial LASSO
#' 
#' @author 
#' Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#'
#' @description 
#' Implementation of the online spatio-temporal process monitoring procedure 
#' described in Qiu and Yang (2023), in which spatial locations with the detected 
#' shifts are guaranteed to be small clustered spatal regions by the exponentially 
#' weighted spatial LASSO.
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
#' for online process monitoring (cf., Yang and Qiu 2021). If there are only 
#' data points with either \code{type='IC1'} or \code{type='IC2'}, then these 
#' data points will be used to estimate the model and conduct the bootstrap 
#' procedure as well. This function will return an error if there are no 
#' observations with \code{type='IC1'} or \code{type='IC2'}.
#' @param ARL0 
#' The pre-specified IC average run length. Default is 200.
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
#' Qiu, P. and Yang, K. (2023). Spatio-Temporal Process Monitoring Using 
#' Exponentially Weighted Spatial LASSO. \emph{Journal of Quality Technology}, 
#' \strong{55}, 163-180.
#'
#' @useDynLib SpTe2M,.registration = TRUE
#' 
#' @keywords sptemnt_ewsl
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
#' \donttest{ili.ewsl <- sptemnt_ewsl(y.sub,st.sub,type.sub,ht=0.05,hs=6.5,gt=0.25,gs=1.5)}

sptemnt_ewsl <- function(y,st,type,ARL0=200,lambda=0.1,B=1000,bs=5,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
{ # split the data into 3 parts: IC1, IC2, Mnt
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
  
  ## decorrelation
  res.ic <- spte_decor(y.ic1,st.ic1,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
  res.mt <- spte_decor(y,st,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
  
  ## bandwidth selection
  h <- mod_cv(y.ic2,st.ic2)$bandwidth.opt[2]
  
  ## EWKS: lasso weights
  t <- sort(unique(st[,3])); n <- length(t); m <- rep(0,n) 
  for(i in 1:n) {
    m[i] <- length(which(st[,3]==t[i]))
  }
  MAXm <- max(m); RES <- sx <- sy <- matrix(0,n,MAXm)
  for(i in 1:n) {
    ids <- which(st[,3]==t[i])
    len <- length(ids); RES[i,1:len] <- res.mt[ids]
    sx[i,1:len] <- st[ids,1]; sy[i,1:len] <- st[ids,2]
  }
  mu.ewsl <- .Fortran("SpTeEWKS", y=as.matrix(RES), t=as.vector(t), sx=as.matrix(sx),
                      sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                      MAXm=as.integer(MAXm), lambda=as.double(lambda), 
                      h=as.double(h), LOO=as.integer(0), 
                      muhat=matrix(as.double(0),n,MAXm))$muhat
  varpi1 <- 1/pmax(abs(mu.ewsl),1/n); varpi2 <- matrix(0,n,MAXm)
  for(i in 1:n) {
    for(j in 1:m[i]) {
      dist <- rep(0,m[i])
      for(jj in 1:m[i]) {
        dist[jj] <- sqrt((sx[i,j]-sx[i,jj])^2+(sy[i,j]-sy[i,jj])^2)
      }
      tmp <- dist/h; w <- pmax(0,0.75*(1-tmp^2))
      varpi2[i,j] <-  mu.ewsl[i,j]-sum(w*mu.ewsl[i,1:m[i]])/sum(w)
    }
  }
  varpi2 <- pmin(1/pmax(abs(varpi2),1/n),30)/10
  ## process monitoring
  gamma1 <- gamma2 <- c(0,0.1,0.2); LEN <- length(gamma1)
  # compute the charting statistics
  Ct <- array(0,c(n,LEN)); mu.ewsl <- array(0,c(LEN,n,MAXm))
  for(I in 1:LEN) {
    
    gam1 <- gamma1[I]; gam2 <- gamma2[I]
    ## time = 1
    X <- matrix(0,m[1]^2,m[1]); Y <- rep(0,m[1]^2)
    dist <- matrix(0,m[1],m[1])
    for(i in 1:m[1]) {
      for(j in 1:m[1]) {
        dist[i,j] <- sqrt((sx[1,i]-sx[1,j])^2+(sy[1,i]-sy[1,j])^2)
      }
    }    
    for(i in 1:m[1]) {
      tmp <- dist[i,]/h
      w <- pmax(0,0.75*(1-tmp^2))
      X[((i-1)*m[1]+1):(i*m[1]),i] <- sqrt(w)
      Y[((i-1)*m[1]+1):(i*m[1])] <- RES[1,1:m[1]]*sqrt(w)
    }
    index <- which(Y==0); X <- X[-index,]; Y <- Y[-index]
    D <- matrix(0,2*m[1],m[1]); rowsum.D <- rep(0,m[1]) 
    for(i in 1:m[1]) {
      D[i,i] <- gam1*varpi1[1,i]
      rowsum.D[i] <- sum(D[i,]^2)
      vec <- rep(0,m[1]); vec[i] <- 1
      w <- pmax(0,0.75*(1-(dist[i,]/h)^2))
      D[i+m[1],] <- gam2*varpi2[1,i]*(vec-w/sum(w))
      rowsum.D[i+m[1]] <- sum(D[i+m[1],]^2)
    }
    D <- D[which(rowsum.D!=0),]
    
    if(dim(D)[1]==0) {
      mu.ewsl[I,1,] <- c(solve(t(X)%*%X)%*%t(X)%*%Y)
    } else {
      X1 <- X%*%solve(t(D)%*%D)%*%t(D)
      # run the path algorithm
      coef <- as.numeric(coef(glmnet::glmnet(x=X1,y=Y,intercept=FALSE,
                                             family='gaussian',lambda=1/sqrt(length(Y)))))
      coef <- coef[-1]
      mu.ewsl[I,1,] <- solve(t(D)%*%D)%*%t(D)%*%coef
    }
    INDEX <- rep(1,length(Y))
    
    ## time > 1
    for(it in 2:n) {
      dist <- matrix(0,m[it],m[it])
      for(ii in 1:m[it]) {
        for(jj in 1:m[it]) {
          dist[ii,jj] <- sqrt((sx[it,ii]-sx[it,jj])^2+(sy[it,ii]-sy[it,jj])^2)
        }
      }    
      X0 <- X; Y0 <- Y; len <- dim(X0)[1]
      X <- matrix(0,len+m[it]^2,m[it]); Y <- rep(0,len+m[it]^2)
      X[1:len,] <- X0*sqrt(1-lambda);  Y[1:len] <- Y0*sqrt(1-lambda)
      for(i in 1:m[it]) {
        tmp <- dist[i,]/h; w <- pmax(0,0.75*(1-tmp^2))
        X[(len+(i-1)*m[it]+1):(len+i*m[it]),i] <- sqrt(w)
        Y[(len+(i-1)*m[it]+1):(len+i*m[it])] <- RES[it,1:m[it]]*sqrt(w)
      }
      INDEX <- c(INDEX,rep(it,m[it]^2)); index <- which(Y==0) 
      X <- X[-index,]; Y <- Y[-index]; INDEX <- INDEX[-index]
      if(it>5) { 
        index <- which(INDEX==min(INDEX))
        X <- X[-index,]; Y <- Y[-index]; INDEX <- INDEX[-index]
      }
      D <- matrix(0,2*m[it],m[it]); rowsum.D <- rep(0,m[it]) 
      for(i in 1:m[it]) {
        D[i,i] <- gam1*varpi1[it,i]
        rowsum.D[i] <- sum(D[i,]^2)
        vec <- rep(0,m[it]); vec[i] <- 1
        w <- pmax(0,0.75*(1-(dist[i,]/h)^2))
        D[i+m[it],] <- gam2*varpi2[it,i]*(vec-w/sum(w))
        rowsum.D[i+m[it]] <- sum(D[i+m[it],]^2)
      }
      D <- D[which(rowsum.D!=0),]
      if(dim(D)[1]==0) {
        mu.ewsl[I,it,] <- c(solve(t(X)%*%X)%*%t(X)%*%Y)
      } else {
        X1 <- X%*%solve(t(D)%*%D)%*%t(D)
        # run the path algorithm
        coef <- as.numeric(coef(glmnet::glmnet(x=X1,y=Y,intercept=FALSE,
                                               family='gaussian',lambda=1/sqrt(length(Y)))))
        coef <- coef[-1]
        mu.ewsl[I,it,] <- solve(t(D)%*%D)%*%t(D)%*%coef
      }
    }
    
    for(it in 1:n) { 
      Ct[it,I] <- sum(mu.ewsl[I,it,1:m[it]]^2)
    }
  }
  Ct <- apply(Ct,1,max)
  ## bootstrap
  t <- sort(unique(st.ic1[,3])); n <- length(t); m <- rep(0,n) 
  for(i in 1:n) {
    m[i] <- length(which(st.ic1[,3]==t[i]))
  }
  MAXm <- max(m); varpi1 <- matrix(median(varpi1),n,MAXm)
  varpi2 <- matrix(median(varpi2),n,MAXm)
  RES <- sx <- sy <- matrix(0,n,MAXm)
  for(i in 1:n) {
    ids <- which(st.ic1[,3]==t[i])
    len <- length(ids); RES[i,1:len] <- res.ic[ids]
    sx[i,1:len] <- st.ic1[ids,1]; sy[i,1:len] <- st.ic1[ids,2]
  }
  Ct.ic <- array(0,c(n,LEN)); mu.ewsl.ic <- array(0,c(LEN,n,MAXm))
  for(I in 1:LEN) {
    
    gam1 <- gamma1[I]; gam2 <- gamma2[I]
    ## time = 1
    X <- matrix(0,m[1]^2,m[1]); Y <- rep(0,m[1]^2)
    dist <- matrix(0,m[1],m[1])
    for(i in 1:m[1]) {
      for(j in 1:m[1]) {
        dist[i,j] <- sqrt((sx[1,i]-sx[1,j])^2+(sy[1,i]-sy[1,j])^2)
      }
    }    
    for(i in 1:m[1]) {
      tmp <- dist[i,]/h
      w <- pmax(0,0.75*(1-tmp^2))
      X[((i-1)*m[1]+1):(i*m[1]),i] <- sqrt(w)
      Y[((i-1)*m[1]+1):(i*m[1])] <- RES[1,1:m[1]]*sqrt(w)
    }
    index <- which(Y==0); X <- X[-index,]; Y <- Y[-index]
    D <- matrix(0,2*m[1],m[1]); rowsum.D <- rep(0,m[1]) 
    for(i in 1:m[1]) {
      D[i,i] <- gam1*varpi1[1,i]
      rowsum.D[i] <- sum(D[i,]^2)
      vec <- rep(0,m[1]); vec[i] <- 1
      w <- pmax(0,0.75*(1-(dist[i,]/h)^2))
      D[i+m[1],] <- gam2*varpi2[1,i]*(vec-w/sum(w))
      rowsum.D[i+m[1]] <- sum(D[i+m[1],]^2)
    }
    D <- D[which(rowsum.D!=0),]
    if(dim(D)[1]==0) {
      mu.ewsl.ic[I,1,] <- c(solve(t(X)%*%X)%*%t(X)%*%Y)
    } else {
      X1 <- X%*%solve(t(D)%*%D)%*%t(D)
      # the path algorithm
      coef <- as.numeric(coef(glmnet::glmnet(x=X1,y=Y,intercept=FALSE,
                                             family='gaussian',lambda=1/sqrt(length(Y)))))
      coef <- coef[-1]
      mu.ewsl.ic[I,1,] <- solve(t(D)%*%D)%*%t(D)%*%coef
    }
    INDEX <- rep(1,length(Y))
    
    ## time > 1
    for(it in 2:n) {
      dist <- matrix(0,m[it],m[it])
      for(ii in 1:m[it]) {
        for(jj in 1:m[it]) {
          dist[ii,jj] <- sqrt((sx[it,ii]-sx[it,jj])^2+(sy[it,ii]-sy[it,jj])^2)
        }
      }    
      X0 <- X; Y0 <- Y; len <- dim(X0)[1]
      X <- matrix(0,len+m[it]^2,m[it]); Y <- rep(0,len+m[it]^2)
      X[1:len,] <- X0*sqrt(1-lambda);  Y[1:len] <- Y0*sqrt(1-lambda)
      for(i in 1:m[it]) {
        tmp <- dist[i,]/h; w <- pmax(0,0.75*(1-tmp^2))
        X[(len+(i-1)*m[it]+1):(len+i*m[it]),i] <- sqrt(w)
        Y[(len+(i-1)*m[it]+1):(len+i*m[it])] <- RES[it,1:m[it]]*sqrt(w)
      }
      INDEX <- c(INDEX,rep(it,m[it]^2)); index <- which(Y==0) 
      X <- X[-index,]; Y <- Y[-index]; INDEX <- INDEX[-index]
      if(it>5) { 
        index <- which(INDEX==min(INDEX))
        X <- X[-index,]; Y <- Y[-index]; INDEX <- INDEX[-index]
      }
      D <- matrix(0,2*m[it],m[it]); rowsum.D <- rep(0,m[it]) 
      for(i in 1:m[it]) {
        D[i,i] <- gam1*varpi1[it,i]
        rowsum.D[i] <- sum(D[i,]^2)
        vec <- rep(0,m[it]); vec[i] <- 1
        w <- pmax(0,0.75*(1-(dist[i,]/h)^2))
        D[i+m[it],] <- gam2*varpi2[it,i]*(vec-w/sum(w))
        rowsum.D[i+m[it]] <- sum(D[i+m[it],]^2)
      }
      D <- D[which(rowsum.D!=0),]
      if(dim(D)[1]==0) {
        mu.ewsl.ic[I,it,] <- c(solve(t(X)%*%X)%*%t(X)%*%Y)
      } else {
        X1 <- X%*%solve(t(D)%*%D)%*%t(D)
        # the path algorithm 
        coef <- as.numeric(coef(glmnet::glmnet(x=X1,y=Y,intercept=FALSE,
                                               family='gaussian',lambda=1/sqrt(length(Y)))))
        coef <- coef[-1]
        mu.ewsl.ic[I,it,] <- solve(t(D)%*%D)%*%t(D)%*%coef
      }
    }
    
    for(it in 1:n) { 
      Ct.ic[it,I] <- sum(mu.ewsl.ic[I,it,1:m[it]]^2)
    }
  }
  Ct.ic <- apply(Ct.ic,1,max); mu <- mean(Ct.ic); sd <- sd(Ct.ic)
  Ct[1:max(floor(n/10),1)] <- Ct[1:max(floor(n/10),1)]/sqrt(n/10)
  Ct.ic <- Ct.ic/sqrt(n/2.25)
  set.seed(12345)
  LEN <- 2*bs*n; Ct.b <- matrix(0,B,LEN)
  for(b in 1:B) {
    id <- sample(1:(n-bs+1),2*n,replace=TRUE)
    for(i in 1:(2*n)) {
      Ct.b[b,((i-1)*bs+1):(i*bs)] <- Ct.ic[id[i]:(id[i]+bs-1)]
    }
  }
  # the bi-section algorithm for selecting the control limit
  hl <- 0; hu <- 100; ARL <- rep(0,B)
  while(abs(mean(ARL)-ARL0)>0.001) {
    h <- (hl+hu)/2
    for(b in 1:B) {
      id <- which(Ct.b[b,]>h); id <- c(id,LEN)
      ARL[b] <- min(id)
    }
    if(mean(ARL)>ARL0) { hu <- h }
    if(mean(ARL)<ARL0) { hl <- h }
    if(abs(hu-hl)<0.001) { break() }
  }
  CL <- h # the control limit
  
  results <- list()
  results$ARLO <- ARL0; results$lambda <- lambda
  results$cstat <- Ct; results$cl <- CL
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
