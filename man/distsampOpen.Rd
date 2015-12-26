\name{distsampOpen}
\alias{distsampOpen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Open population model for distance sampling data
}
\description{
  Fit the model of Dail and Madsen (2011) with a multinomial instead of
  a binomial observation model.
}
\usage{
distsampOpen(lambdaformula, gammaformula, omegaformula, sigmaformula, data, mixture = c("P", "NB", "ZIP"), K, dynamics = c("constant", "autoreg", "notrend", "trend"), fix = c("none", "gamma", "omega"), starts, method = "BFGS", se = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambdaformula}{
%%     ~~Describe \code{lambdaformula} here~~
}
  \item{gammaformula}{
%%     ~~Describe \code{gammaformula} here~~
}
  \item{omegaformula}{
%%     ~~Describe \code{omegaformula} here~~
}
  \item{sigmaformula}{
%%     ~~Describe \code{sigmaformula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{mixture}{
%%     ~~Describe \code{mixture} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{dynamics}{
%%     ~~Describe \code{dynamics} here~~
}
  \item{fix}{
%%     ~~Describe \code{fix} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


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
        y[i,1:J,t+1] <- rmultinom(1, N[i,t], cp)[1:J]
        }
    }
    cp <- array(cp, c(J, M, T))
    cp <- matrix(aperm(cp, c(2,1,3)), M)
    cat("max(N) =", max(N), "\n")
    return(matrix(y, M))
}


lambda <- 4
gamma <- 2
omega <- 0.5
sigma <- 20
T <- 20

# Do one case
y.sim7 <- sim7(lambda, gamma, omega, sigma, M=200, T=T,type="line")
colSums(y.sim7)
umf7 <- unmarkedFramePCO(y = y.sim7, numPrimary=T)
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7, K=40,
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE)

# Do a simulation study
simout<- matrix(NA,nrow=100,ncol=4)
for(sim in 1:100){
y.sim7 <- sim7(lambda, gamma, omega, sigma, M=200, T=T,type="line")
colSums(y.sim7)
umf7 <- unmarkedFramePCO(y = y.sim7, numPrimary=T)
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7, K=40,
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE)
simout[sim,]<- fm@opt$par

}

#
#  > apply(simout,2,mean)
#  [1] 1.38359217 0.67011133 0.03584008 2.99542523
#  > exp(1.3835)
#  [1] 3.988838
#  > exp(0.67)
#  [1] 1.954237
#  > plogis(0.036)
#  [1] 0.508999
#  > exp(2.995)
#  [1] 19.98536






}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

