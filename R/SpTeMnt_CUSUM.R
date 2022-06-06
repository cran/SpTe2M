SpTeMnt_CUSUM <- function(y,st,type,ARL0=200,gamma=0.1,B=1000,bs=5,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
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
    res.ic <- SpTe_DeCor(y.ic1,st.ic1,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
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
    CL <- h
    
    ## process monitoring
    res <- SpTe_DeCor(y,st,y.ic2,st.ic2,T,ht,hs,gt,gs)$std.res
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
    Cstat <- rep(0,n)
    Cstat[1] <- max(0,(sum(RES[1,1:m[1]]^2)-m[1])/sqrt(2*m[1])-gamma)
    for(i in 2:n) {
        Cstat[i] <- max(0,Cstat[i-1]+(sum(RES[i,1:m[i]]^2)-m[i])/sqrt(2*m[i])-gamma)
    }
    results <- list()
    results$ARLO <- ARL0; results$gamma <- gamma
    results$Cstat <- Cstat; results$CL <- CL
    ids <- which(results$Cstat>CL)
    if(length(ids)==0) {
      id <- n+1
    } else {
      id <- min(ids)
    }
    results$signal_time <- max(1,id-1)
    return(results)
}
