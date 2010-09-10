\name{gmm}
\alias{gmm}
\title{Generalized multinomial mixture model}
\description{
A three level hierarchical model with a multinomial observation process
}
\usage{
gmm(lambdaformula, phiformula, pformula, data, mixture = c("P", "NB"), K, 
    starts, method = "BFGS", control = list(), se = TRUE)
}
\arguments{
    \item{lambdaformula}{Righthand side (RHS) formula describing abundance 
        covariates}
    \item{phiformula}{RHS formula describing availability covariates}
    \item{pformula}{RHS formula describing detection covariates}
    \item{data}{An object of class unmarkedFrameGMM}
    \item{mixture}{Either "P" or "NB" for Poisson and Negative Binomial mixing 
        distributions.}
    \item{K}{The upper bound of integration}
    \item{starts}{Starting values}
    \item{method}{Optimization method used by \code{\link{optim}}}
    \item{control}{List of control arguments passed to \code{\link{optim}}}
    \item{se}{Logical. Should standard errors be calculated?}  
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
Richard Chandler and Andy Royle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

n <- 50   # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals availabe for detection
    
for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)
    Nleft1 <- N[i,] - y[i,,1]
    y[i,,2] <- rbinom(T, Nleft1, p)
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }
    
apply(y, 3, sum)    

y.MxTJ <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])

umf <- unmarkedFrameGMM(y=y.MxTJ, numPrimary=T, piFun="removalPiFun")

summary(m <- gmm(~1, ~1, ~1, data=umf))
e <- coef(m)
# Should be 3, 0.5, 0.3
round(c(exp(e[1]), plogis(e[2]), plogis(e[3])), 2)

}
\keyword{ model }
