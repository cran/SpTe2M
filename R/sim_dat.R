#' @title A simulated spatio-temporal dataset
#' @name sim_dat
#' @docType data
#' @author  Kai Yang \email{kayang@mcw.edu} and Peihua Qiu
#' @description 
#' This simulated dataset is saved as a list, and it contains the following 
#' three elements:
#' \describe{
#' \item{y}{A vector of length \eqn{N}; it contains the data of the observed 
#' response variable \eqn{y}.}
#' \item{x}{A vector of length \eqn{N}; it contains the data of the covariate 
#' \eqn{x}. }
#' \item{st}{An \eqn{N\times 3} matrix containing the spatial locations and 
#' times for all the observations in the dataset.}}
#' @format 
#' A list containing \eqn{N=10,000} observations.
#' @import MASS ggplot2 knitr rmarkdown
#' @keywords sim_dat
#' @usage data(sim_dat)
#' @examples 
#' library(MASS)
#' set.seed(100)
#' n <- 100; m <- 100; N <- n*m
#' t <- rep(seq(0.01,1,0.01),each=m)
#' su <- sv <- seq(0.1,1,0.1)
#' su <- rep(su,each=10); sv <- rep(sv,10)
#' su <- rep(su,n); sv <- rep(sv,n)
#' st <- matrix(0,N,3)
#' st[,1] <- su; st[,2] <- sv; st[,3] <- t
#' mu <- rep(0,N)
#' for(i in 1:N) {
#'   mu[i] <- 2+sin(pi*su[i])*sin(pi*sv[i])+sin(2*pi*t[i]) 
#' }
#' dist <- matrix(0,m,m) # distance matrix
#' for(i in 1:m) {
#'   for(j in 1:m) {
#'     dist[i,j] <- sqrt((su[i]-su[j])^2+(sv[i]-sv[j])^2)
#'   }
#' }
#' cov.s <- matrix(0,m,m) # spatial correlation
#' for(i in 1:m) {
#'   for(j in 1:m) {
#'     cov.s[i,j] <- 0.3^2*exp(-30*dist[i,j]) 
#'   }
#' }
#' noise <- matrix(0,n,m)
#' noise[1,] <- MASS::mvrnorm(1,mu=rep(0,m),Sigma=cov.s) 
#' for(i in 2:n) {
#'   noise[i,] <- 0.1*noise[i-1,]+sqrt(1-0.1^2)*
#'     MASS::mvrnorm(1,mu=rep(0,m),Sigma=cov.s)
#' }
#' noise <- c(t(noise)); x <- rnorm(N,0,0.3) 
#' beta <- 0.5; y <- mu+x*beta+noise
#' sim_dat <- list(); sim_dat$y <- y
#' sim_dat$x <- x; sim_dat$st <- st
NULL