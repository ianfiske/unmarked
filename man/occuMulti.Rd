\name{occuMulti}

\alias{occuMulti}

\title{Fit the Rota et al. (2016) Multi-species Occupancy Model}

\usage{occuMulti(detformulas, stateformulas, data, starts, method="BFGS",
    se=TRUE, engine=c("R"), ...)}

\arguments{
    \item{detformulas}{List of formulas for the detection models, one per species.}
    \item{stateformulas}{List of formulas for the natural (*f*) parameters.}
    \item{data}{An \code{\link{unmarkedFrameOccuMulti}} object}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Either "C" to use fast C++ code or "R" to use native R
      code during the optimization.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits the multispecies occupancy model of Rota et al (2016).}

\details{

See \code{\link{unmarkedFrame}} and \code{\link{unmarkedFrameOccuMulti}} for a
description of how to supply data to the \code{data} argument.

\code{occuMulti} fits the multispecies occupancy model based on Rota et al. (2016).

}

\value{unmarkedFitOccuMulti object describing the model fit.}

\references{
Rota, C.T., et al. 2016. A multi-species occupancy model for two or more
  interacting species. Methods in Ecology and Evolution 7: 1164-1173.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuMulti}}}


\examples{

#Simulate data
N <- 100
nspecies <- 2
J <- 5

occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('par',1:10,sep='')

det_covs <- as.data.frame(matrix(rnorm(N*J*3),ncol=3))
names(det_covs) <- paste('par',1:3,sep='')

#True vals
beta <- c(0.5,0.5,-0.5,1.2,-1)
f1 <- beta[1] + beta[2]*occ_covs$par1
f2 <- beta[3] + beta[4]*occ_covs$par2
f3 <- beta[5]
f <- cbind(f1,f2,f3)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

psi <- exp(f %*% t(dm))
psi <- psi/rowSums(psi)

#True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
  ztruth[i,] <- as.matrix(z[sample(4,1,prob=psi[i,]),])
}

p_true <- c(0.6,0.7)

# fake y data

y1 <- y2 <- matrix(NA,N,J)
for (i in 1:N){
  for (j in 1:J){
    y1[i,j] <- rbinom(1,1,ztruth[i,1]*p_true[1])
    y2[i,j] <- rbinom(1,1,ztruth[i,2]*p_true[2])
}}

y <- list(y1,y2)

#Create data object
data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

# Formulas
occFormulas <- list('~par1','~par2','~1')
detFormulas <- list('~1','~1')

#Fit model
fit <- occuMulti(detFormulas,occFormulas,data,se=T)

fit
}

\keyword{models}
