SpTe_SemiParmReg <- function(y,st,x,ht=NULL,hs=NULL,maxIter=1000,tol=10^(-4),stE=NULL) 
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
        mcv <- Mod_CV(y,st,eps=0.1)
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
