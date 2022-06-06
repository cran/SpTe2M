SpTeMnt_EWMAC <- function(y,x,st,type,ARL0=200,ARL0.z=200,lambda=0.1,B=1000,bs=5,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
{  
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
    beta <- SpTe_SemiParmReg(y.ic2,st.ic2,x.ic2,ht,hs,maxIter=100)$beta
    z.ic1 <- x.ic1%*%beta; z.ic2 <- x.ic2%*%beta
    z <- x%*%beta
    
    ## block bootstrap
    res.y.ic <- SpTe_DeCor(y.ic1,st.ic1,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    res.z.ic <- SpTe_DeCor(z.ic1,st.ic1,z.ic2,st.ic2,T,ht,hs,gt,gs)$std.res

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
    CL <- h
    
    ## process monitoring
    res.y <- SpTe_DeCor(y,st,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
    res.z <- SpTe_DeCor(z,st,z.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
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
    results$Cstat <- Ct.y; results$CL <- CL
    ids <- which(results$Cstat>CL)
    if(length(ids)==0) {
      id <- n+1
    } else {
      id <- min(ids)
    }
    results$signal_time <- max(1,id-1)
    return(results)
}
