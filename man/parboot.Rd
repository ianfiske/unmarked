\name{parboot}
\alias{parboot}
\alias{parboot-class}
\alias{show,parboot-method}
\title{Evaluate goodness-of-fit of a fitted model.}
\description{Evaluate goodness-of-fit of a fitted model.}
\arguments{
	\item{object}{a fitted model of class "unmarkedFit"}
	\item{nsim}{number of bootstrap replicates}
	\item{report}{print fit statistic every 'report' iterations during resampling}}
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

(pb <- parboot(fm))
plot(pb)

}
