\name{distsampOpen}
\alias{distsampOpen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Open population model for distance sampling data
}
\description{
  Fit the model of Dail and Madsen (2011) and Hostetler and Chandler
  (2015) with a distance sampling
  observation model. 
}
\usage{
distsampOpen(lambdaformula, gammaformula, omegaformula, sigmaformula, data, mixture = c("P", "NB", "ZIP"), K, 
dynamics = c("constant", "autoreg", "notrend", "trend"), fix = c("none", "gamma", "omega"), starts, method = "BFGS", 
se = TRUE, nintervals=10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambdaformula}{
        Right-hand sided formula for initial abundance
}
\item{gammaformula}{
  Right-hand sided formula for recruitment rate (when dynamics is "constant",
    "autoreg", or "notrend") or population growth rate (when dynamics is
    "trend", "ricker", or "gompertz")
}
\item{omegaformula}{
  Right-hand sided formula for apparent survival probability
    (when dynamics is "constant", "autoreg", or "notrend") or equilibrium
    abundance (when dynamics is "ricker" or "gompertz")
}
\item{sigmaformula}{
      Right-hand sided formula for detection probability
}
\item{data}{
      An object of class \code{\link{unmarkedFramePCO}}. See details
XXXXX will change soon XXXXX 
}
\item{mixture}{
  character specifying mixture: "P", "NB", or "ZIP" for
    the Poisson, negative binomial, and zero-inflated Poisson
    distributions.
}
\item{K}{
   Integer defining upper bound of discrete integration. This
    should be higher than the maximum observed count and high enough
    that it does not affect the parameter estimates. However, the higher
    the value the slower the compuatation.
}
\item{dynamics}{
  Character string describing the type of population
    dynamics. "constant" indicates that there is no relationship between
    omega and gamma. "autoreg" is an auto-regressive model in which
    recruitment is modeled as gamma*N[i,t-1]. "notrend" model gamma as
    lambda*(1-omega) such that there is no temporal trend. "trend" is
    a model for exponential growth, N[i,t] = N[i,t-1]*gamma, where gamma
    in this case is finite rate of increase (normally referred to as
    lambda). "ricker" and "gompertz" are models for density-dependent
    population growth.  "ricker" is the Ricker-logistic model, N[i,t] =
    N[i,t-1]*exp(gamma*(1-N[i,t-1]/omega)), where gamma is the maximum
    instantaneous population growth rate (normally referred to as r) and
    omega is the equilibrium abundance (normally referred to as K).  "gompertz"
    is a modified version of the Gompertz-logistic model, N[i,t] =
    N[i,t-1]*exp(gamma*(1-log(N[i,t-1]+1)/log(omega+1))), where the
    interpretations of gamma and omega are similar to in the Ricker model.
}
\item{fix}{
  If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0.
}
\item{starts}{
      vector of starting values
}
\item{method}{
      Optimization method used by \code{\link{optim}}.
}
\item{se}{
   logical specifying whether or not to compute standard errors.
 }
   \item{immigration}{
    logical specifying whether or not to include an immigration term (iota) in
    population dynamics.
  }
  \item{iotaformula}{
    Right-hand sided formula for average number of immigrants to a site
    per time step
  }
  \item{nintervals}{
    number of intervals used in computing the  integral of the distance function by the rectangular rule in order
    to compute multinomial cell probabilities.  This is NOT the number of distance bands but rather the number of intervals that each band is divided into.
    Default is 10. 
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Dail, D. and L. Madsen (2011) Models for Estimating Abundance from
Repeated Counts of an Open Metapopulation. \emph{Biometrics}. 67, pp
577-587.

Hostetler, J. A. and R. B. Chandler (2015) Improved State-space Models for
Inference about Spatial and Temporal Variation in Abundance from Count Data.
\emph{Ecology} 96:1713-1723.
}
\author{
Richard Chandler, Jeff Hostetler, Andy Royle 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pcount}, \link{unmarkedFramePCO}}
}
\examples{

sim7 <- function(lambda=1, gamma=0.5, omega=0.8, sigma=40, scale=NULL, M=100, T=5, J=4,type="line",
keyfun="halfnorm")
    {
    y <- array(NA, c(M, J, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
if(keyfun=="halfnorm"){
    if(type=="point")
       g <- function(x, sig) exp(-x^2/(2*sig^2))*x
    if(type=="line")
       g <- function(x, sig) exp(-x^2/(2*sig^2))
}
if(keyfun=="exp"){
    if(type=="point")
       g <- function(x, sig) exp(-x/sig)*x
    if(type=="line")
       g <- function(x, sig) exp(-x/sig)
}
if(keyfun=="hazard"){
    shape<- sigma
    if(type=="point")
       g <- function(x, shape, scale) (1 - exp(-(r/shape)^-scale)) * x
    if(type=="line")
       g <- function(x, sig) (1 - exp(-(x/shape)^-scale) )
}
if(keyfun=="uniform"){
     if(type=="point")
       g <-  function(x,sig) x
    if(type=="line")
       g <- function(x,sig) 1
}

    cp <- u <- a <- numeric(J)
if(keyfun!="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value
    }
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

}
if(keyfun=="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
   for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
 
    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]
 
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
 
    }
    }
    u <- a / sum(a)
    cp <- a
    cp[j+1] <- 0

}
    for(i in 1:M) {
    N[i,1] <- rpois(1, lambda)

    y[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]

    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
    }
    cp <- array(cp, c(J, M, T))
    cp <- matrix(aperm(cp, c(2,1,3)), M)
    cat("max(N) =", max(N), "\n")
    return(list(y=matrix(y, M),N=N))
}
library(unmarked)
set.seed(711)
lambda <- 4
gamma <- 2
omega <- 0.5
sigma <- 20
T <- 20
 
# half-normal examples 
#NEED TO ADD WARNING ABOUT TLENGTH!!!!!
y.sim7 <- sim7(lambda, gamma, omega, sigma=15, M=200, T=T,type="line", keyfun="halfnorm")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="halfnorm",
                 #starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 starts=c(0,0,0,2.21),
                 se=TRUE)

fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="halfnorm",
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma), 0),
                 mixture="NB",
                 se=TRUE)

y.sim7 <- sim7(lambda, gamma, omega, sigma=10, M=50, T=T,type="line", keyfun="halfnorm")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="halfnorm",
                 #starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 starts=c(0,0,0,2.204),
                 se=TRUE)
# exponential example 
y.sim7 <- sim7(lambda, gamma, omega, sigma=25, M=200, T=T,type="line", keyfun="exp")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))

fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="exp",
                 #starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 starts=c(0,0,0,2),
                   se=TRUE)

# uniform example 
y.sim7 <- sim7(lambda, gamma, omega, sigma=NULL, M=200, T=T,type="line", keyfun="uniform", scale=1)$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,
                   starts=rep(0,3),
                   keyfun="uniform", se=TRUE)

fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,
                   starts=rep(0,4),
                   mixture="NB",
                   keyfun="uniform", se=TRUE)

# hazard example
y.sim7 <- sim7(lambda, gamma, omega, sigma=20, scale = 1, M=200, T=T,type="line", keyfun="hazard")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="hazard",
                 #starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma), -1),
                 starts=c(0,0,0,0,0),
                 se=TRUE)




##
##
## do not run
###
##
sim7 <- function(lambda=1, gamma=0.5, omega=0.8, sigma=40, M=100, T=5, J=4,type="line")
    {
    y <- array(NA, c(M, J, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
    if(type=="point")
    g <- function(x, sig) exp(-x^2/(2*sig^2))*x
    if(type=="line")
    g <- function(x, sig) exp(-x^2/(2*sig^2)) 
    cp <- u <- a <- numeric(J)
  if(type=="point"){ 
   a[1] <- pi*db[2]^2
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    }
if(type=="line"){ 
   L <-  1
   a[1] <- L*db[2]
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value  
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value 
    }
    }

    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)
    for(i in 1:M) {
    N[i,1] <- rpois(1, lambda)
    y[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]
    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
    }
    cp <- array(cp, c(J, M, T))
    cp <- matrix(aperm(cp, c(2,1,3)), M)
    cat("max(N) =", max(N), "\n")
    return(list(y=matrix(y, M),N=N))
}
library(unmarked)
set.seed(711)
lambda <- 12
gamma <- 2
omega <- 0.5
sigma <- 20
T <- 20
simout<- matrix(NA,nrow=100,ncol=4)
for(sim in 1:100){
y.sim7 <- sim7(lambda, gamma, omega, sigma, M=200, T=T,type="line")$y
colSums(y.sim7)
umf7 <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200) )

fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7, K=150,
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE, nintervals = 12)
simout[sim,]<- fm@opt$par

}
colMeans(simout)

# results for lambda = 4 seem to work
# results of using nintervals=5
#> apply(simout,2,mean)
#[1] 1.37376613 0.65539517 0.05327946 2.99762181
 


# results using the analytic solution to the cell probabilities
# apply(simout,2,mean)
#[1] 1.37633799 0.65631764 0.05651401 2.99502916
# exp(1.376)
#[1] 3.959034
# exp(0.656)
#[1] 1.927069
#> plogis(0.056)
#[1] 0.5139963
#> exp(2.995)
#[1] 19.98536
 


library(unmarked)

set.seed(12)
y.sim7 <- sim7(lambda=6, gamma, omega, sigma=50, M=200, T=T,type="line")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T, 
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))

# want to test this out
 
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=120,keyfun="half",
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE, nintervals=5)

set.seed(123)
y.sim7 <- sim7(lambda=12, gamma, omega, sigma=50, M=200, T=T,type="point")
umf7b <- unmarkedFrameDSO(y = y.sim7$y, numPrimary=T, 
    dist.breaks = c(0, 25, 50, 75, 100), survey="point", unitsIn="m" )

# want to test this out
 
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=100,keyfun="half",
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE, nintervals=5)

victory<-  ranef(fm)
Nhat.mn<- bup(victory)
Nhat.md<- bup(victory,"mode")

 N<- y.sim7$N
n<- matrix(NA,nrow=200,ncol=T)
mn<- matrix(NA,nrow=200,ncol=T)
idx<- sort(rep(1:20,4))
for(s in 1:20){
 tmpy<- y.sim7$y[s,]
 n[s,]<- tapply(tmpy,idx,sum)
 pp<- victory@post[s,,]
 mn[s,]<- apply( (row(pp)-1)*pp,2,sum)
 }




## Now simulate the point transect case
set.seed(711)
lambda <- 4
gamma <- 2
omega <- 0.5
sigma <- 20
T <- 20
simout<- matrix(NA,nrow=100,ncol=4)
for(sim in 1:100){
y.sim7 <- sim7(lambda, gamma, omega, sigma, M=200, T=T,type="point")
colSums(y.sim7)
umf7 <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="point", unitsIn="m" )

fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7, K=40,
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE, nintervals = 5)
simout[sim,]<- fm@opt$par

}


# results:
#> apply(simout,2,mean)
#[1]  1.37105167  0.64604188 -0.09652492  2.99713509
#> exp(1.37)
#[1] 3.935351
#> exp(0.646)
#[1] 1.907894
#> plogis(-0.0965)
#[1] 0.4758937
#> exp(2.997)
#[1] 20.02537


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

