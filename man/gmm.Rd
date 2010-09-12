\name{gmm}
\alias{gmm}
\title{Generalized multinomial-mixture model}
\description{
A three level hierarchical model with a multinomial observation process. This 
model is useful when both availability and detection probability influence 
observed counts. A robust survey design is required with both primary 
and secondary sampling periods.
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
The latent transect-level super-population abundance distribution 
\eqn{f(M | \mathbf{\theta})}{f(M | theta)} can be set as either a Poisson or 
a negative binomial random variable, depending on the setting of the 
\code{mixture} argument. \code{mixture = "P"} or \code{mixture = "NB"} select 
the Poisson or negative binomial distribution respectively.  The mean of 
\eqn{M_i} is \eqn{\lambda_i}{lambda_i}.  If \eqn{M_i \sim NB}{M_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The number of individuals available for detection at time \eqn{\t}{t} 
is a modeled as binomial: 
\eqn{N_{it} \sim Binomial(M_i, \mathbf{\phi_ij})}{M_i ~ Binomial(M_i, phi_ij)}.

The detection process is modeled as multinomial: 
\eqn{y_{ij} \sim 
Multinomial(N_it, pi_{itj})}{y_itj ~ 
Multinomial(N_it, pi_it1, pi_it2, ..., pi_itJ)}, 
where \eqn{pi_itj} is the multinomial cell probability for plot i at time t on 
occasion j. 

Cell probabilities are computed via a user-defined function related to the 
sampling design. Alternatively, the default functions \code{\link{removalPiFun}}
or \code{\link{doublePiFun}} can be used for equal-interval removal sampling or
double observer sampling. Note that the function for computing cell probabilites
is specified when setting up the data using \code{\link{unmarkedFrameGMM}}.

Parameters \eqn{\lambda}{lambda}, \eqn{\phi}{phi} and \eqn{p}{p} can be 
modeled as linear functions of covariates using the log, logit and logit links
respectively.}

}
\value{
An object of class unmarkedFitGMM.
}
\references{
Royle, J. A. (2004) Generalized estimators of avian abundance from count 
survey data. \emph{Animal Biodiversity and Conservation} 27, pp. 375--386.
}
\author{
Richard Chandler \email{richard.chandlers@gmail.com} and Andy Royle
}
\note{
Three types of covariates can be supplied, site-level, site-by-year-level, and 
observation-level. These must be formatted correctly when organizing the data
with \code{\link{unmarkedFrameGMM}}
}

\seealso{
\code{\link{unmarkedFrameGMM}} for setting up the data and metadata. 
\code{\link{multinomPois}} for surveys where no secondary sampling periods were
used. Example functions to calculate multinomial cell probabilities are 
described \code{\link{piFuns}} 
}
\examples{

n <- 100  # number of sites
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


system.time(m <- gmm(~1, ~1, ~1, data=umf))
e <- coef(m)
# Should be 3, 0.5, 0.3
round(c(exp(e[1]), plogis(e[2]), plogis(e[3])), 2)




# An example with covariates


set.seed(37)

n <- 50   # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

sc <- rnorm(n)
ysc <- rnorm(n*T)
ysc <- matrix(ysc, n, T)
yr <- factor(rep(1:T, n))
oc <- rnorm(n*J*T)
oc <- array(oc, c(n, J, T))
pi <- array(NA, c(n, J, T))

lam <- exp(-1 + 1*sc)
phi <- plogis(2 + -2*ysc)
p <- plogis(1 + -1*oc)

y <- array(NA, c(n,J,T))
M <- rpois(n, lam)
N <- matrix(NA, n, T)

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi[i,])
    y[i,1,] <- rbinom(T, N[i,], p[i,1,])
    Nleft1 <- N[i,] - y[i,1,]
    y[i,2,] <- rbinom(T, Nleft1, p[i,2,])
    Nleft2 <- Nleft1 - y[i,2,]
    y[i,3,] <- rbinom(T, Nleft2, p[i,3,])
    }
    
apply(y, 2, sum)


y[1, 1, 1] <- NA


umf <- unmarkedFrameGMM(y=matrix(y, nrow=n), 
    siteCovs = data.frame(sc=sc), 
    obsCovs=list(oc=matrix(oc, nrow=n)),
    yearlySiteCovs=data.frame(ysc=as.numeric(ysc), yr=yr), 
    numPrimary=T, piFun="removalPiFun")
    
summary(umf)

system.time(m2 <- gmm(~sc, ~ysc, ~oc, umf)) # 6.77
summary(m2)


}
\keyword{ model }
