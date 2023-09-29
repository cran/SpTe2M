## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(SpTe2M)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("SpTe2M")
#  library(SpTe2M)

## -----------------------------------------------------------------------------
data(ili_dat)
head(ili_dat)

## ----fig.height = 6, fig.width = 7, fig.align = "center", eval=FALSE----------
#  library(maps)
#  library(ggplot2)
#  library(mapproj)
#  # turn shape files in the maps packages into a data frame
#  FL <- map_data('county','florida')
#  names(FL)[6] <- 'County'
#  # Only plot maps on Jun 15 and Dec 15
#  subdat <- subset(ili_dat,Date%in%c('6/15/2012','6/15/2013','6/15/2014',
#                            '12/15/2012','12/15/2013','12/15/2014'))
#  subdat$Date <- factor(subdat$Date,levels=c('6/15/2012','6/15/2013','6/15/2014',
#                                         '12/15/2012','12/15/2013','12/15/2014'))
#  # merge ILI data with map data
#  mydat <- merge(FL,subdat)
#  # make maps using ggplot
#  maps <- ggplot(data=mydat,aes(long,lat,group=group,fill=Rate))+geom_polygon()+
#    facet_wrap(~Date,ncol=3)+geom_path(colour='grey10',lwd=0.5)+
#    scale_fill_gradient2('',low='cyan',mid='white',high='navy',
#                         guide='colorbar',limits=c(0,0.0003),na.value='orange',
#                         breaks=c(0,0.0001,0.0002,0.0003),
#                         labels=c('0','1e-4','2e-4','3e-4'))+
#    guides(fill=guide_colorbar(barwidth=25,barheight=1,direction='horizontal'),
#           cex=1)+theme_bw(base_size=15)+xlab('longitude')+ylab('latitude')+
#    theme(legend.position='bottom',
#          axis.ticks=element_blank(),
#          line=element_blank(),
#          axis.text=element_blank(),
#          panel.border=element_rect(color='black',linewidth=1.2),
#          axis.line=element_line(colour='black'),
#          legend.margin=margin(-10,0,0,0),
#          legend.box.margin=margin(0,0,0,0))
#  # display the maps
#  maps

## ----eval=FALSE---------------------------------------------------------------
#  n <- 365; m <- 67; N1 <- (n+1)*m; N2 <- n*m
#  # extract the observed ILI data in year 2013
#  ili13 <- ili_dat[(N1+1):(N1+N2),]
#  y13 <- ili13$Rate; st13 <- ili13[,c('Lat','Long','Time')]

## ----eval=FALSE---------------------------------------------------------------
#  # estimate the mean
#  mu.est <- spte_meanest(y=y13,st=st13)

## ----fig.height = 6, fig.width = 7, fig.align = "center", eval=FALSE----------
#  mu <- mu.est$muhat; mu <- t(matrix(mu,m,n))
#  obs <- t(matrix(y13,m,n))
#  ids <- c(6,34,52,57) # IDs for Broward, Lake, Pinellas, Seminole
#  par(mfrow=c(2,2),mgp=c(2.4,1,0))
#  par(mar=c(3.5,3.5,1.5,0.5))
#  plot(1:365,mu[,ids[1]],type='l',lty=1,lwd=1.5,xaxt='n',
#       ylim=c(0,8e-5),main='Broward',xlab='Time',ylab='Incidence rate',
#       cex=1.2,cex.lab=1.3,cex.axis=1.2,cex.main=1.3)
#  points(1:365,obs[,ids[1]],cex=1)
#  axis(1,cex.axis=1.2,at=c(1+c(1,62,123,184,245,306, 367)),
#       label=c('Jan','Mar','May','July','Sep','Nov','Jan'))
#  par(mar=c(3.5,3,1.5,1))
#  plot(1:365,mu[,ids[2]],type='l',lty=1,lwd=1.5,xaxt='n',
#       ylim=c(0,8e-5),main='Lake',xlab='Time',ylab='',
#       cex=1.2,cex.lab=1.3,cex.axis=1.2,cex.main=1.3)
#  points(1:365,obs[,ids[2]],cex=1)
#  axis(1,cex.axis=1.2,at=c(1+c(1,62,123,184,245,306, 367)),
#       label=c('Jan','Mar','May','July','Sep','Nov','Jan'))
#  par(mar=c(3.5,3.5,1.5,0.5))
#  plot(1:365,mu[,ids[3]],type='l',lty=1,lwd=1.5,xaxt='n',
#       ylim=c(0,8e-5),main='Pinellas',xlab='Time',ylab='Incidence rate',
#       cex=1.2,cex.lab=1.3,cex.axis=1.2,cex.main=1.3)
#  points(1:365,obs[,ids[3]],cex=1)
#  axis(1,cex.axis=1.2,at=c(1+c(1,62,123,184,245,306, 367)),
#       label=c('Jan','Mar','May','July','Sep','Nov','Jan'))
#  par(mar=c(3.5,3,1.5,1))
#  plot(1:365,mu[,ids[4]],type='l',lty=1,lwd=1.5,xaxt='n',
#       ylim=c(0,8e-5),main='Seminole',xlab='Time',ylab=' ',
#       cex=1.2,cex.lab=1.3,cex.axis=1.2,cex.main=1.3)
#  points(1:365,obs[,ids[4]],cex=1)
#  axis(1.2,at=c(1+c(1,62,123,184,245,306, 367)),
#       label=c('Jan','Mar','May','July','Sep','Nov','Jan'))

## ----eval=FALSE---------------------------------------------------------------
#  # estimate the covariance
#  cov.est <- spte_covest(y=y13,st=st13)

## ----eval=FALSE---------------------------------------------------------------
#  y <- ili_dat$Rate; st <- ili_dat[,c('Lat','Long','Time')]
#  # specify the argument type
#  # data with type=IC1 are used to determine the control limit
#  # data with type=IC2 are used to estimate the regular pattern
#  # data with type =Mnt are used for sequential process monitoring
#  type <- rep(c('IC1','IC2','Mnt'),c(N1,N2,N2))

## ----eval=FALSE---------------------------------------------------------------
#  ili.cusum <- sptemnt_cusum(y,st,type,ht=0.05,hs=6.5,gt=0.25,gs=1.5)

## ----fig.height = 3, fig.width = 7, fig.align = "center", eval=FALSE----------
#  cstat <- ili.cusum$cstat; cl <- ili.cusum$cl
#  par(mfrow=c(1,2),mgp=c(2.1,1,0))
#  # plot the CUSUM chart
#  par(mar=c(3.5,3.5,1.5,0.5))
#  plot(1:n,cstat,type="l",xlab="Time",xaxt="n",ylab="Charting statistic",
#       main='CUSUM',cex=0.6,lwd=1.5,cex.lab=0.7,
#       cex.axis=0.7,cex.main=0.8)
#  abline(h=cl,lty=2,lwd=1.5)
#  axis(1,cex.axis=0.7,at=c(1,59,121,182,243,304,365),
#       labe= c('Jan','Mar','May','July','Sep','Nov','Jan'))
#  # plot the zoom-in part
#  par(mar=c(3.5,3,1.5,1))
#  plot(259:305,cstat[259:305],type="b",xlab="Time",xaxt="n",ylab=" ",
#       main='First signal time: 10/16/2014',cex=0.6,lwd=1.5,cex.lab=0.7,
#       cex.axis=0.7,cex.main=0.8)
#  abline(h=cl,lty=2,lwd=1.5)
#  axis(1,cex.axis=0.7,at=c(259,274,289,304),
#       label=c('Sep 15','Sep 30','Oct 15','Oct 31'))

## ----eval=FALSE---------------------------------------------------------------
#  x <- as.matrix(ili_dat[,c('Temp','RH')]) # the covariates
#  ili.ewmac <- sptemnt_ewmac(y,x,st,type,ht=0.05,hs=6.5,gt=0.25,gs=1.5)

## ----fig.height = 3, fig.width = 7, fig.align = "center", eval=FALSE----------
#  cstat <- ili.ewmac$cstat; cl <- ili.ewmac$cl
#  par(mfrow=c(1,2),mgp=c(2.1,1,0))
#  # plot the EWMAC chart
#  par(mar=c(3.5,3.5,1.5,0.5))
#  plot(1:n,cstat,type="l",xlab="Time",xaxt="n",ylab="Charting statistic",
#       main='EWMAC',cex=0.6,lwd=1.5,cex.lab=0.7,
#       cex.axis=0.7,cex.main=0.8)
#  abline(h=cl,lty=2,lwd=1.5)
#  axis(1,cex.axis=0.7,at=c(1,59,121,182,243,304,365),
#       labe= c('Jan','Mar','May','July','Sep','Nov','Jan'))
#  # plot its zoom-in part
#  par(mar=c(3.5,3,1.5,1))
#  plot(259:305,cstat[259:305],type="b",xlab="Time",xaxt="n",ylab=" ",
#       main='First signal time: 9/23/2014',cex=0.6,lwd=1.5,cex.lab=0.7,
#       cex.axis=0.7,cex.main=0.8)
#  abline(h=cl,lty=2,lwd=1.5)
#  axis(1,cex.axis=0.7,at=c(259,274,289,304),
#       label=c('Sep 15','Sep 30','Oct 15','Oct 31'))

## ----eval=FALSE---------------------------------------------------------------
#  ili.ewsl <- sptemnt_ewsl(y,st,type,ht=0.05,hs=6.5,gt=0.25,gs=1.5)

## ----fig.height = 3, fig.width = 7, fig.align = "center", eval=FALSE----------
#  cstat <- ili.ewsl$cstat; cl <- ili.ewsl$cl
#  par(mfrow=c(1,2),mgp=c(2.1,1,0))
#  # plot the EWSL chart
#  par(mar=c(3.5,3.5,1.5,0.5))
#  plot(1:n,cstat,type="l",xlab="Time",xaxt="n",ylab="Charting statistic",
#       main='EWSL',cex=0.6,lwd=1.5,cex.lab=0.7,
#       cex.axis=0.7,cex.main=0.8)
#  abline(h=cl,lty=2,lwd=1.5)
#  axis(1,cex.axis=0.8,at=c(1,59,121,182,243,304,365),
#       labe= c('Jan','Mar','May','July','Sep','Nov','Jan'))
#  # plot its zoom-in part
#  par(mar=c(3.5,3,1.5,1))
#  plot(259:305,cstat[259:305],type="b",xlab="Time",xaxt="n",ylab=" ",
#       main='First signal time: 10/6/2014',cex=0.6,lwd=1.5,cex.lab=0.7,
#       cex.axis=0.7,cex.main=0.8)
#  abline(h=cl,lty=2,lwd=1.5)
#  axis(1,cex.axis=0.7,at=c(259,274,289,304),
#       label=c('Sep 15','Sep 30','Oct 15','Oct 31'))

