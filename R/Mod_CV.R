Mod_CV <- function(y,st,ht=NULL,hs=NULL,eps=0.1) 
{   
    if(is.null(ht)==TRUE || is.null(hs)==TRUE) {
        sx.rg <- range(st[,1]); sy.rg <- range(st[,2])
        s.rg <- sqrt((sx.rg[2]-sx.rg[1])^2+(sy.rg[2]-sy.rg[1])^2)
        t.rg <- range(st[,3]); t.rg <- t.rg[2]-t.rg[1]
        ht <- seq(from=0.1*t.rg,to=0.6*t.rg,by=0.5*t.rg/9)
        hs <- seq(from=0.2*s.rg,to=0.8*s.rg,by=0.6*s.rg/9)
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
    Nbw <- length(ht); NUM <- dim(st)[1]
    modCV <- .Fortran("ModCV", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), NUM=as.integer(NUM), ht=as.vector(ht), 
                       hs=as.vector(hs), Nbw=as.integer(Nbw), eps=as.double(eps),
                       mcv=rep(as.double(0),Nbw))
    mcv <- modCV$mcv; results <- list()
    bw <- matrix(0,2,Nbw); bw[1,] <- ht; bw[2,] <- hs 
    row.names(bw) <- c('ht','hs'); results$bandwidth <- bw 
    results$mcv <- mcv; id <- which.min(mcv)
    results$bandwidth.opt <- bw[,id]; results$mcv.opt <- mcv[id]
    return(results)
}
