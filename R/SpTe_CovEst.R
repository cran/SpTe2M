SpTe_CovEst <- function(y,st,gt=NULL,gs=NULL,stE1=NULL,stE2=NULL) 
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
    NUM <- dim(st)[1]; mcv <- Mod_CV(y,st,eps=0.1)
    ht.opt <- mcv$bandwidth.opt[1]; hs.opt <- mcv$bandwidth.opt[2]
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
    if(is.null(gt)==TRUE || is.null(gs)==TRUE) {
       mspe <- CV_MSPE(y,st)
       gt <- mspe$bandwidth.opt[1]; gs <- mspe$bandwidth.opt[2]
    }
    if(is.null(stE1)==TRUE) {
       stE1 <- st
    }
    if(is.null(stE2)==TRUE) {
        stE2 <- st
    }
    NUM1 <- dim(stE1)[1]; NUM2 <- dim(stE2)[1]
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
