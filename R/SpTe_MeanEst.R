SpTe_MeanEst <- function(y,st,ht=NULL,hs=NULL,cor=FALSE,stE=NULL) 
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
        mcv <- Mod_CV(y,st,eps=0.1)
        ht <- mcv$bandwidth.opt[1]; hs <- mcv$bandwidth.opt[2]
    }
    if(is.null(stE)==TRUE) {
        stE <- st
    }
    NUM <- dim(stE)[1]; NUM0 <- dim(st)[1]
    if(cor==FALSE) {
        mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                           sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                           ht=as.double(ht), hs=as.double(hs), stE=as.matrix(stE), 
                           NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    } else {
        mspe <- CV_MSPE(y,st)
        gt0 <- mspe$bandwidth.opt[1]; gs0 <- mspe$bandwidth.opt[2]
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
