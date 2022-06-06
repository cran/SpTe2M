CV_MSPE <- function(y,st,gt=NULL,gs=NULL) 
{   
    if(is.null(gt)==TRUE || is.null(gs)==TRUE) {
        sx.rg <- range(st[,1]); sy.rg <- range(st[,2])
        s.rg <- sqrt((sx.rg[2]-sx.rg[1])^2+(sy.rg[2]-sy.rg[1])^2)
        t.rg <- range(st[,3]); t.rg <- t.rg[2]-t.rg[1]
        gt <- seq(from=t.rg*0.4,to=0.8*t.rg,by=0.4*t.rg/10)
        gs <- seq(from=s.rg*0.4,to=0.8*s.rg,by=0.4*s.rg/10)
    }
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
    Nbw <- length(gt); NUM <- dim(st)[1]
    mcv <- Mod_CV(y,st,eps=0.1)
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
    sx.rg <- range(st[,1]); sy.rg <- range(st[,2])
    s.rg <- sqrt((sx.rg[2]-sx.rg[1])^2+(sy.rg[2]-sy.rg[1])^2)
    t.rg <- range(st[,3]); t.rg <- t.rg[2]-t.rg[1]
    dt <- t.rg/n*3; ds <- 2*s.rg/sqrt(MAXm)
    cvMSPE <- .Fortran("CVMSPE", res=as.matrix(RES), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), NUM=as.integer(NUM), gt=as.vector(gt), 
                       gs=as.vector(gs), Nbw=as.integer(Nbw), dt=as.double(dt),
                       ds=as.double(ds), mspe=rep(as.double(0),Nbw))
    mspe <- cvMSPE$mspe; results <- list()
    bw <- matrix(0,2,Nbw); bw[1,] <- gt; bw[2,] <- gs 
    row.names(bw) <- c('gt','gs'); results$bandwidth <- bw 
    results$mspe <- mspe; id <- which.min(mspe)
    results$bandwidth.opt <- bw[,id]; results$mspe.opt <- mspe[id]
    return(results)
}
