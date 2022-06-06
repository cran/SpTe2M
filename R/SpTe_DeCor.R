SpTe_DeCor <- function(y,st,y0,st0,T=1,ht=NULL,hs=NULL,gt=NULL,gs=NULL) 
{   
    y.ic <- y0; st.ic <- st0
    st.original <- st
    NUM <- dim(st.ic)[1]
    for(i in 1:NUM) {
        st.ic[i,3] <- st.ic[i,3] %% T
    }
    NUM <- dim(st)[1]
    for(i in 1:NUM) {
        st[i,3] <- st[i,3] %% T
    }
    
    t <- sort(unique(st.ic[,3])); n <- length(t)
    m <- rep(0,n)
    for(i in 1:n) {
        m[i] <- length(which(st.ic[,3]==t[i]))
    }
    MAXm <- max(m); Y <- sx <- sy <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st.ic[,3]==t[i])
        len <- length(ids)
        Y[i,1:len] <- y.ic[ids]
        sx[i,1:len] <- st.ic[ids,1]
        sy[i,1:len] <- st.ic[ids,2]
    }
    NUM <- dim(st.ic)[1]
    if(is.null(ht)==TRUE || is.null(hs)==TRUE) {
        mcv <- Mod_CV(y.ic,st.ic,eps=0.1)
        ht <- mcv$bandwidth.opt[1]; hs <- mcv$bandwidth.opt[2]
    } 
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                       ht=as.double(ht), hs=as.double(hs), stE=as.matrix(st.ic), 
                       NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    muhat <- mu.est$muhat
    res <- y.ic-muhat; RES <- matrix(0,n,MAXm)
    for(i in 1:n) {
        ids <- which(st.ic[,3]==t[i])
        len <- length(ids)
        RES[i,1:len] <- res[ids]
    }
    if(is.null(gt)==TRUE || is.null(gs)==TRUE) {
       mspe <- CV_MSPE(y.ic,st.ic)
       gt <- mspe$bandwidth.opt[1]; gs <- mspe$bandwidth.opt[2]
    }
    
    NUM <- dim(st)[1]
    mu.est <- .Fortran("SpTeLLKS", y=as.matrix(Y), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m), MAXm=as.integer(MAXm), 
                       ht=as.double(ht), hs=as.double(hs), stE=as.matrix(st), 
                       NUM=as.integer(NUM), eps=as.double(0), muhat=rep(as.double(0),NUM))
    muhat <- mu.est$muhat
    cov.est <- .Fortran("SpTeWME", res=as.matrix(RES), t=as.vector(t), sx=as.matrix(sx),
                       sy=as.matrix(sy), n=as.integer(n), m=as.integer(m),
                       MAXm=as.integer(MAXm), gt=as.double(gt), gs=as.double(gs), 
                       stE1=as.matrix(st), NUM1=as.integer(NUM),
                       stE2=as.matrix(st), NUM2=as.integer(NUM),
                       covhat=matrix(as.double(0),NUM,NUM))
    covhat <- cov.est$covhat
    
    res <- y-muhat
    t <- sort(unique(st[,3])); n <- length(t)
    m <- rep(0,n)
    for(i in 1:n) {
        m[i] <- length(which(st[,3]==t[i]))
    }
    MAXm <- max(m); res.std <- rep(0,NUM)
    
    ## Time=1
    Sig <- covhat[1:m[1],1:m[1]]; eig <- eigen(Sig) 
    M <- eig$vectors%*%diag(eig$values^(-0.5))%*%t(eig$vectors)
    res.std[1:m[1]] <- M%*%res[1:m[1]]
    ## Time>1
    resid <- c(res[1:m[1]]); SIGMA <- Sig; num <- m[1]
    for(i in 2:n) {
        num.old <- sum(m[1:(i-1)]); num.new <- sum(m[1:i])
        Sigma <- covhat[(1+num.old-num):num.old,(1+num.old):num.new]
        SIG <- covhat[(1+num.old):num.new,(1+num.old):num.new]
        Sig <- SIG-t(Sigma)%*%solve(SIGMA)%*%Sigma
        M <- eig$vectors%*%diag(eig$values^(-0.5))%*%t(eig$vectors)
        res.std[(1+num.old):num.new] <- M%*%(res[(1+num.old):num.new]-t(Sigma)%*%solve(SIGMA)%*%resid)
        ## update
        resid <- c(resid,res[(1+num.old):num.new])
        SIGMA <- rbind(cbind(SIGMA,Sigma), cbind(t(Sigma),SIG))
        num <- num+m[i]
        if(i>5) {
            resid <- resid[-(1:m[i-5])]
            SIGMA <- SIGMA[(m[i-5]+1):dim(SIGMA)[1],(m[i-5]+1):dim(SIGMA)[1]]
            num <- num-m[i-5]
        }
    }
    results <- list(); results$st <- st.original; results$std.res <- res.std
    return(results)
}
