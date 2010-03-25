\name{parboot}
\alias{parboot}
\alias{plot,parboot,missing-method}
\alias{show,parboot-method}
\title{Parametric bootstrap method for fitted models.}
\description{Evaluate goodness-of-fit of a fitted model.}
\arguments{
	\item{object}{a fitted model of class "unmarkedFit"}
	\item{statistic}{a function returning a vector of fit-statistics. 
        First argument must be the fitted model. 
        Default is sum of squared residuals.}
   	\item{nsim}{number of bootstrap replicates}
	\item{report}{print fit statistic every 'report' iterations during resampling}
    \item{...}{Additional arguments to be passed to statistic}
    }
\examples{
data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 
lengths <- linetran$Length

ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths*1000, survey = "line", unitsIn = "m")
    })

(fm <- distsamp(~ area ~habitat, ltUMF))

(pb <- parboot(fm, fitStats, nsim=25))
plot(pb, main="")

}
