\name{occuMS}

\alias{occuMS}

\title{Fit Single-Season Multi-State Occupancy Models}

\usage{occuMS(detformulas, stateformulas, data, parameterization="multinomial",
    starts, method="BFGS", se=TRUE, engine=c("C","R"), silent=FALSE, ...)}

\arguments{
    \item{detformulas}{Character vector of formulas for detection probabilities.
      See details for a description of how to order these formulas.}
    \item{stateformulas}{Character vector of formulas for occupancy probabilities. 
      See details for a description of how to order these formulas.}
    \item{data}{An \code{\link{unmarkedFrameOccuMS}} object}
    \item{parameterization}{Either \code{"multinomial"} for the multinomial 
      parameterization (MacKenzie et al. 2009) which allows an arbitrary
      number of occupancy states, or \code{"condbinom"} for the conditional 
      binomial parameterization (Nichols et al. 2007) which requires exactly
      3 occupancy states. See details.}
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

\description{This function fits single-season multi-state occupancy models with both the multinomial and conditional binomial parameterizations.}

\details{

Traditional occupancy models fit data with exactly two states: detection and
non-detection (MacKenzie et al. 2002).
The \code{occuMS} function fits models to occupancy data for which there are 
greater than 2 states (Nichols et al 2007, MacKenzie et al. 2009). For example, 
detections may be further divided into multiple biologically relevant categories, 
e.g. breeding vs. non-breeding, or some/many individuals present. As with
detection status, classification of these additional occupancy states is likely
to be imperfect.

Multiple parameterizations for multi-state occupancy models have been proposed.
The \code{occuMS} function fits two at present: the "conditional binomial" 
parameterization of Nichols et al. (2007), and the more general "multinomial"
parameterization of MacKenzie et al. (2009). At present, only single-season
models are possible with \code{occuMS}, although dynamic models have been
developed and we hope to add them later (MacKenzie et al. 2009).

The conditional binomial parameterization 
(\code{parameterization = 'condbinom'}) models occupancy and the presence or
absence of an additional biological state of interest given the species
is present (typically breeding status). Thus, there should be exactly 3 occupancy 
states in the data: 0 (non-detection); 1 (detection, no evidence of breeding);
or 2 (detection, evidence of breeding). 

Two state parameters are estimated:
\eqn{\psi}, the probability of occupancy, and \eqn{R}, the probability of 
successful reproduction given an occupied state (although this could be some
other binary biological condition). Covariates (in \code{siteCovs}) can be 
supplied for either or both of these parameters with the \code{stateformulas} 
argument, which takes a character vector of R-style formulas with length = 2, 
with formulas in the order (\eqn{\psi}, \eqn{R}). For example, to fit a model 
where \eqn{\psi} varies with a landcover covariate and \eqn{R} is constant,
\code{stateformulas = c('~landcover','~1')}. 

There are three detection parameters associated with the
conditional binomial parameterization: \eqn{p_1}, the probability of 
detecting the species given true state 1; \eqn{p_2}, the probability of detecting
the species given true state 2; and \eqn{\delta}, the probability of detecting 
state 2 (i.e., breeding), given that the species has been detected.
See MacKenzie et al. (2009), pages 825-826 for more details.
As with occupancy, covariates (in \code{obsCovs}) can be supplied for these 
detection probabilities with the \code{detformulas} argument, which takes a 
character vector of formulas with length = 3 in the order 
(\eqn{p_1}, \eqn{p_2}, \eqn{\delta}). So, to fit a model where \eqn{p_1} varies 
with temperature and the other two parameters are constant, 
\code{detformulas = c('~temp','~1','~1')}.

The multinomial parameterization (\code{parameterization = "multinomial"}) is
more general, allowing an arbitrary number of occupancy states \eqn{S}.
\eqn{S} - 1 occupancy probabilities \eqn{\psi} are estimated. Thus, if there
are \eqn{S} = 4 occupancy states (0, 1, 2, 3), \code{occuMS} estimates \eqn{\psi_1},
\eqn{\psi_2}, and \eqn{\psi_3} (the probability of state 0 can be obtained by 
subtracting the others from 1). Covariates can be supplied for each occupancy
probability with a character vector with length \eqn{S-1}, e.g. 
\code{stateformulas =} \code{c('~landcover','~1','~1')} where \eqn{\psi_1} varies with
landcover and \eqn{\psi_2} and \eqn{\psi_3} are constant.

The number of detection probabilities estimated quickly expands as \eqn{S}
increases, equal to \eqn{S \times (S-1) / 2}. In the simplest case
(when \eqn{S} = 3), there are 3 detection probabilities: \eqn{p_{11}}, 
the probability of detecting state 1 given true state 1; \eqn{p_{12}}, 
the probability of detecting state 1 given true state 2; and \eqn{p_{22}}, 
the probability of detecting state 2 given true state 2. 
Covariates can be supplied for any or all of these detection probabilities with 
the \code{detformulas} argument, which takes a character vector of formulas 
with length = 3 in the order (\eqn{p_{11}}, \eqn{p_{12}}, \eqn{p_{22}}). So, 
to fit a model where \eqn{p_{11}} varies with temperature and the other two detection 
probabilities are constant, \code{detformulas = c('~temp','~1','~1')}.
If there were \eqn{S} = 4 occupancy states, there are 6 estimated detection
probabilities and the order is (\eqn{p_{11}}, \eqn{p_{12}}, \eqn{p_{13}},
\eqn{p_{22}}, \eqn{p_{23}}, \eqn{p_{33}}), and so on. See MacKenzie et al. (2009)
for a more detailed explanation.
 
See \code{\link{unmarkedFrame}} and \code{\link{unmarkedFrameOccuMS}} for a
description of how to supply data to the \code{data} argument.
}

\value{unmarkedFitOccuMS object describing the model fit.}

\references{
MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

  MacKenzie, D. I., Nichols, J. D., Seamans, M. E., and R. J. Gutierrez, 2009. 
  Modeling species occurrence dynamics with multiple states and imperfect 
  detection. Ecology 90: 823-835.

Nichols, J. D., Hines, J. E., Mackenzie, D. I., Seamans, M. E., and 
  R. J. Gutierrez. 2007. Occupancy estimation and modeling with multiple states 
  and state uncertainty. Ecology 88: 1395-1400.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuMS}}}


\examples{

\dontrun{

#Simulate data

#Parameters
N <- 500; J <- 5; S <- 3
site_covs <- matrix(rnorm(N*2),ncol=2)
obs_covs <- matrix(rnorm(N*J*2),ncol=2)
a1 <- -0.5; b1 <- 1; a2 <- -0.6; b2 <- -0.7

##################################
## Multinomial parameterization ##
##################################

p11 <- -0.4; p12 <- -1.09; p22 <- -0.84
truth <- c(a1,b1,a2,b2,p11,0,p12,p22)

#State process
lp <- matrix(NA,ncol=S,nrow=N)
for (n in 1:N){
  lp[n,2] <- exp(a1+b1*site_covs[n,1])
  lp[n,3] <- exp(a2+b2*site_covs[n,2])
  lp[n,1] <- 1  
}
psi_mat <- lp/rowSums(lp)

z <- rep(NA,N)
for (n in 1:N){
  z[n] <- sample(0:2, 1, replace=T, prob=psi_mat[n,])
}

probs_raw <- matrix(c(1,0,0,1,exp(p11),0,1,exp(p12),exp(p22)),nrow=3,byrow=T)
probs_raw <- probs_raw/rowSums(probs_raw)
  
y <- matrix(0,nrow=N,ncol=J)
for (n in 1:N){

  probs <- switch(z[n]+1,
                  probs_raw[1,],
                  probs_raw[2,],
                  probs_raw[3,])
  if(z[n]>0){
    y[n,] <- sample(0:2, J, replace=T, probs)
  }
}

#Construct unmarkedFrame
umf <- unmarkedFrameOccuMS(y=y,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

#Formulas

#3 states, so detformulas is a character vector of formulas of 
#length 3 in following order:
#1) p[11]: prob of detecting state 1 given true state 1
#2) p[12]: prob of detecting state 1 given true state 2
#3) p[22]: prob of detecting state 2 given true state 2
detformulas <- c('~V1','~1','~1')
#If you had 4 states, it would be p[11],p[12],p[13],p[22],p[23],p[33] and so on

#3 states, so stateformulas is a character vector of length 2 in following order:
#1) psi[1]: probability of state 1
#2) psi[2]: probability of state 2
#You can get probability of state 0 (unoccupied) as 1 - psi[1] - psi[2]
stateformulas <- c('~V1','~V2')

#Fit model
fit <- occuMS(detformulas, stateformulas, umf, parameterization="multinomial")

#Look at results
fit
#Compare with truth
cbind(truth=truth,estimate=coef(fit))

#Generate predicted values
lapply(predict(fit,type='state'),head)
lapply(predict(fit,type='det'),head)

#Fit a null model
detformulas <- rep('~1',3)
stateformulas <- rep('~1',2)
fit_null <- occuMS(detformulas, stateformulas, umf, parameterization="multinomial")

#Compare fits
modSel(fitList(fit,fit_null))

###########################################
## Conditional binomial parameterization ##
###########################################

p11 <- 0.4; p12 <- 0.6; p22 <- 0.8
truth_cb <- c(a1,b1,a2,b2,qlogis(p11),0,qlogis(c(p12,p22)))

#Simulate data

#State process
psi_mat <- matrix(NA,ncol=S,nrow=N)
for (n in 1:N){
  psi_mat[n,2] <- plogis(a1+b1*site_covs[n,1])
  psi_mat[n,3] <- plogis(a2+b2*site_covs[n,2])
}
psi_bin <- matrix(NA,nrow=nrow(psi_mat),ncol=ncol(psi_mat))
psi_bin[,1] <- 1-psi_mat[,2]
psi_bin[,2] <- (1-psi_mat[,3])*psi_mat[,2]
psi_bin[,3] <- psi_mat[,2]*psi_mat[,3]
z <- rep(NA,N)
for (n in 1:N){
  z[n] <- sample(0:2, 1, replace=T, prob=psi_bin[n,])
}

#Detection process
y_cb <- matrix(0,nrow=N,ncol=J)
for (n in 1:N){
  #p11 = p1; p12 = p2; p22 = delta
  probs <- switch(z[n]+1,
                  c(1,0,0),
                  c(1-p11,p11,0),
                  c(1-p12,p12*(1-p22),p12*p22)) 
  if(z[n]>0){
    y_cb[n,] <- sample(0:2, J, replace=T, probs)
  }
}

#Build unmarked frame
umf2 <- unmarkedFrameOccuMS(y=y_cb,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

#Formulas

#detformulas is a character vector of formulas of length 3 in following order:
#1) p[1]: prob of detecting species given true state 1
#2) p[2]: prob of detecting species given true state 2
#3) delta: prob of detecting state 2 (eg breeding) given species was detected
detformulas <- c('~V1','~1','~1')

#stateformulas is a character vector of length 2 in following order:
#1) psi: probability of occupancy
#2) R: probability state 2 (eg breeding) given occupancyc
stateformulas <- c('~V1','~V2')

#Fit model
fit_cb <- occuMS(detformulas, stateformulas, umf2, parameterization='condbinom')

#Look at results
fit_cb
#Compare with truth
cbind(truth=truth_cb,estimate=coef(fit_cb))

#Generate predicted values
lapply(predict(fit_cb,type='state'),head)
lapply(predict(fit_cb,type='det'),head)
}

}

\keyword{models}
