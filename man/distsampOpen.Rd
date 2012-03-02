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

sim7 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5, J=4)
    {
    y <- matrix(NA, M, J*T)
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, p)
    return(y)
}

lambda <- 1
gamma <- 0.5
omega <- 0.8
p <- 0.7
T <- 5
y.sim7 <- sim7(lambda, gamma, omega, p, T=T)
umf7 <- unmarkedFramePCO(y = y.sim7, numPrimary=T)
m7 <- distsampOpen(~1, ~1, ~1, ~1, umf7, K=15,
                 starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
                 se=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

