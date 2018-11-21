\name{occuMulti}

\alias{occuMulti}

\title{Fit the Rota et al. (2016) Multi-species Occupancy Model}

\usage{occuMulti(detformulas, stateformulas, data, starts, method="BFGS",
    se=TRUE, engine=c("C","R"), silent=FALSE, ...)}

\arguments{
    \item{detformulas}{List of formulas for the detection models, one per species.}
    \item{stateformulas}{List of formulas for the natural parameters.}
    \item{data}{An \code{\link{unmarkedFrameOccuMulti}} object}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Either "C" to use fast C++ code or "R" to use native R
      code during the optimization.}
    \item{silent}{Boolean; if \code{TRUE}, suppress warnings.}
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
\dontrun{

#Simulate 3 species `data
N <- 1000
nspecies <- 3
J <- 5

occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('par',1:10,sep='')

det_covs <- as.data.frame(matrix(rnorm(N*J*3),ncol=3))
names(det_covs) <- paste('par',1:3,sep='')

#True vals
beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
f1 <- beta[1] + beta[2]*occ_covs$par1
f2 <- beta[3] + beta[4]*occ_covs$par2
f3 <- beta[5] + beta[6]*occ_covs$par3
f4 <- beta[7]
f5 <- beta[8]
f6 <- beta[9]
f7 <- beta[10]
f <- cbind(f1,f2,f3,f4,f5,f6,f7)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
colnames(z) <- paste('sp',1:nspecies,sep='')
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

psi <- exp(f %*% t(dm))
psi <- psi/rowSums(psi)

#True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
  ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
}

p_true <- c(0.6,0.7,0.5)

# fake y data
y <- list()

for (i in 1:nspecies){
  y[[i]] <- matrix(NA,N,J)
  for (j in 1:N){
    for (k in 1:J){
      y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
    }
  }
}
y[[1]][1,5] <- y[[2]][1,5] <- y[[3]][1,5] <- NA
names(y) <- c('coyote','tiger','bear')

##############

data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

#Summary
summary(data)
plot(data)

# Look at f parameter design matrix
data@fDesign

# Formulas

# Length should match number/order of columns in fDesign
occFormulas <- list('~par1','~par2','~par3','~1','~1','~1','~1')

#Length should match number/order of species in data@ylist
detFormulas <- list('~1','~1','~1')

fit <- occuMulti(detFormulas,occFormulas,data,se=T)

#Look at output
fit
plot(fit)

#Compare with known values
cbind(c(beta,log(p_true/(1-p_true))),fit@opt$par)

#predict method
lapply(predict(fit,'state'),head)
lapply(predict(fit,'det'),head)

#marginal occupancy
head(predict(fit,'state',species=2))
head(predict(fit,'state',species='bear'))
head(predict(fit,'det',species='coyote'))

#conditional occupancy
head(predict(fit,'state',species=2,cond=3)) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='bear')) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='-bear')) #bear absent
head(predict(fit,'state',species='tiger',cond=c('coyote','-bear')))

#residuals (by species)
lapply(residuals(fit),head)

#ranef (by species)
ranef(fit, species='coyote')

#parametric bootstrap
bt <- parboot(fit,nsim=30)

#update model
occFormulas <- list('~par1','~par2','~par2+par3','~1','~1','~1','~1')
fit2 <- update(fit,stateformulas=occFormulas)

#List of fitted models
fl <- fitList(fit,fit2)
coef(fl)

#Model selection
modSel(fl)

}

}

\keyword{models}
