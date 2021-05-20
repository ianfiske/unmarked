\name{occuTTD}
\alias{occuTTD}
\title{Fit Single-Season and Dynamic Time-to-detection Occupancy Models}

\usage{occuTTD(psiformula= ~1, gammaformula =  ~ 1, epsilonformula = ~ 1,
    detformula = ~ 1, data, ttdDist = c("exp", "weibull"), 
    linkPsi = c("logit", "cloglog"), starts, method="BFGS", se=TRUE, 
    engine = c("C", "R"), ...)}

\arguments{
  \item{psiformula}{Right-hand sided formula for the initial probability of 
    occupancy at each site.}
  \item{gammaformula}{Right-hand sided formula for colonization probability.}
  \item{epsilonformula}{Right-hand sided formula for extinction probability.}
  \item{detformula}{Right-hand sided formula for mean time-to-detection.}
  \item{data}{\code{unmarkedFrameOccuTTD} object that supplies the data 
    (see \code{\link{unmarkedFrameOccuTTD}}).}
  \item{ttdDist}{Distribution to use for time-to-detection; either
    \code{"exp"} for the exponential, or \code{"weibull"} for the Weibull,
    which adds an additional shape parameter \eqn{k}.}
  \item{linkPsi}{Link function for the occupancy model. Options are  
    \code{"logit"} for the standard occupancy model or \code{"cloglog"} 
    for the complimentary log-log link, which relates occupancy
    to site-level abundance.}
  \item{starts}{optionally, initial values for parameters in the optimization.}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{engine}{Either "C" or "R" to use fast C++ code or native R
    code during the optimization.}
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}
}

\description{Fit time-to-detection occupancy models of Garrard et al. 
  (2008, 2013), either single-season or dynamic. Time-to-detection can be 
  modeled with either an exponential or Weibull distribution.} 

\value{unmarkedFitOccuTTD object describing model fit.}

\details{

Estimates site occupancy and detection probability from time-to-detection 
(TTD) data, e.g. time to first detection of a particular bird species 
during a point count or time-to-detection of a plant species while searching 
a quadrat (Garrard et al. 2008). Time-to-detection can be modeled 
as an exponential (\code{ttdDist="exp"}) or Weibull (\code{ttdDist="weibull"}) 
random variable with rate parameter \eqn{\lambda} and, for the Weibull, 
an additional shape parameter \eqn{k}. Note that \code{occuTTD} puts covariates
on \eqn{\lambda} and not \eqn{1/\lambda}, i.e., the expected time between events.

In the case where there are no detections before the maximum sample time at
a site (\code{surveyLength}) is reached, we are not sure if the site is 
unoccupied or if we just didn't wait long enough for a detection. We therefore 
must censor the exponential or Weibull distribution at the maximum survey 
length, \eqn{Tmax}. Thus, assuming true site occupancy at site \eqn{i} is 
\eqn{z_i}, an exponential distribution for the TTD \eqn{y_i}, and that 
\eqn{d_i = 1} indicates \eqn{y_i} is censored (Kery and Royle 2016): 

\deqn{d_i = z_i * I(y_i > Tmax_i) + (1 - z_i)}

and

\deqn{y_i|z_i \sim Exponential(\lambda_i), d_i = 0}
\deqn{y_i|z_i = Missing, d_i = 1}

Because in \code{unmarked} values of \code{NA} are typically used to indicate 
missing values that were a result of the sampling structure (e.g., lost data), 
we indicate a censored \eqn{y_i} in \code{occuTTD} instead by setting 
\eqn{y_i = Tmax_i} in the \code{y} matrix provided to 
\code{\link{unmarkedFrameOccuTTD}}. You can provide either a single value of 
\eqn{Tmax} to the \code{surveyLength} argument of \code{unmarkedFrameOccuTTD}, 
or provide a matrix, potentially with a unique value of \eqn{Tmax} for each 
value of \code{y}. Note that in the latter case the value of \code{y} that will 
be interpreted by \code{occuTTD} as a censored observation (i.e., \eqn{Tmax}) 
will differ between observations!

Occupancy and detection can be estimated with only a single survey per site, 
unlike a traditional occupancy model that requires at least two replicated 
surveys at at least some sites. However, \code{occuTTD} also supports 
multiple surveys per site using the model described in Garrard et al. (2013). 
Furthermore, multi-season dynamic models are supported, using the same basic 
structure as for standard occupancy models (see \code{\link{colext}}). 
    
When \code{linkPsi = "cloglog"}, the complimentary log-log link 
function is used for \eqn{psi} instead of the logit link. The cloglog link
relates occupancy probability to the intensity parameter of an underlying
Poisson process (Kery and Royle 2016). Thus, if abundance at a site is 
can be modeled as \eqn{N_i ~ Poisson(\lambda_i)}, where 
\eqn{log(\lambda_i) = \alpha + \beta*x}, then presence/absence data at the 
site can be modeled as \eqn{Z_i ~ Binomial(\psi_i)} where 
\eqn{cloglog(\psi_i) = \alpha + \beta*x}. 

}

\references{

Garrard, G.E., Bekessy, S.A., McCarthy, M.A. and Wintle, B.A. 2008. When have 
  we looked hard enough? A novel method for setting minimum survey effort 
  protocols for flora surveys. Austral Ecology 33: 986-998.

Garrard, G.E., McCarthy, M.A., Williams, N.S., Bekessy, S.A. and Wintle, 
  B.A. 2013. A general model of detectability using species traits. Methods in 
  Ecology and Evolution 4: 45-52.

Kery, Marc, and J. Andrew Royle. 2016. \emph{Applied Hierarchical Modeling in
  Ecology}, Volume 1. Academic Press. 
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuTTD}}}

\keyword{models}

\examples{

\dontrun{

### Single season model
N <- 500; J <- 1

#Simulate occupancy
scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                    forest=runif(N,0,1),
                    wind=runif(N,0,1))

beta_psi <- c(-0.69, 0.71, -0.5)
psi <- plogis(cbind(1, scovs$elev, scovs$forest) \%*\% beta_psi)
z <- rbinom(N, 1, psi)

#Simulate detection
Tmax <- 10 #Same survey length for all observations
beta_lam <- c(-2, -0.2, 0.7)
rate <- exp(cbind(1, scovs$elev, scovs$wind) \%*\% beta_lam)
ttd <- rexp(N, rate)
ttd[z==0] <- Tmax #Censor at unoccupied sites
ttd[ttd>Tmax] <- Tmax #Censor when ttd was greater than survey length

#Build unmarkedFrame
umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

#Fit model
fit <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind, data=umf)

#Predict psi values
predict(fit, type='psi', newdata=data.frame(elev=0.5, forest=1))

#Predict lambda values
predict(fit, type='det', newdata=data.frame(elev=0.5, wind=0))

#Calculate p, probability species is detected at a site given it is present
#for a value of lambda. This is equivalent to eq 4 of Garrard et al. 2008
lam <- predict(fit, type='det', newdata=data.frame(elev=0.5, wind=0))$Predicted
pexp(Tmax, lam)

#Estimated p for all observations
head(getP(fit))

### Dynamic model

N <- 1000; J <- 2; T <- 2
scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                    forest=runif(N,0,1),
                    wind=runif(N,0,1))

beta_psi <- c(-0.69, 0.71, -0.5)
psi <- plogis(cbind(1, scovs$elev, scovs$forest) \%*\% beta_psi)
z <- matrix(NA, N, T)
z[,1] <- rbinom(N, 1, psi)

#Col/ext process
ysc <- data.frame(forest=rep(scovs$forest, each=T), 
                  elev=rep(scovs$elev, each=T))
c_b0 <- -0.4; c_b1 <- 0.3
gam <- plogis(c_b0 + c_b1 * scovs$forest)
e_b0 <- -0.7; e_b1 <- 0.4
ext <- plogis(e_b0 + e_b1 * scovs$elev)

for (i in 1:N){
  for (t in 1:(T-1)){
    if(z[i,t]==1){
      #ext
      z[i,t+1] <- rbinom(1, 1, (1-ext[i]))
    } else {
      #col
      z[i,t+1] <- rbinom(1,1, gam[i])
    }
  }
}

#Simulate detection
ocovs <- data.frame(obs=rep(c('A','B'),N*T))
Tmax <- 10
beta_lam <- c(-2, -0.2, 0.7)
rate <- exp(cbind(1, scovs$elev, scovs$wind) \%*\% beta_lam)
#Add second observer at each site
rateB <- exp(cbind(1, scovs$elev, scovs$wind) \%*\% beta_lam - 0.5)
#Across seasons
rate2 <- as.numeric(t(cbind(rate, rateB, rate, rateB)))
ttd <- rexp(N*T*2, rate2)
ttd <- matrix(ttd, nrow=N, byrow=T)
ttd[ttd>Tmax] <- Tmax
ttd[z[,1]==0,1:2] <- Tmax
ttd[z[,2]==0,3:4] <- Tmax
  
umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength = Tmax, 
                            siteCovs = scovs, obsCovs=ocovs,
                            yearlySiteCovs=ysc, numPrimary=2) 

dim(umf@y) #num sites, (num surveys x num primary periods)

fit <- occuTTD(psiformula=~elev+forest,detformula=~elev+wind+obs,
               gammaformula=~forest, epsilonformula=~elev, 
               data=umf,se=T,engine="C")

truth <- c(beta_psi, c_b0, c_b1, e_b0, e_b1, beta_lam, -0.5)

#Compare to truth
cbind(coef(fit), truth)

}
}
