\name{parboot}
\alias{parboot}
\title{Evaluate goodness-of-fit for a fitted distance-sampling model...}
\description{Evaluate goodness-of-fit for a fitted distance-sampling model}
\arguments{\item{object}{a fitted model of class "umDistsampFit"}
\item{R}{number of bootstrap replicates}
\item{report}{print fit statistic every 'report' iterations during resampling}
\item{label}{a label for the model}}
\examples{data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

lengths <- linetran$Length

ltUMF <- with(linetran, {
unmarkedFrame(y = cbind(dc1, dc2, dc3, dc4), 
siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
tlength = lengths*1000, survey = "line", unitsIn = "m")
})

(fm <- distsamp(~ area ~habitat, ltUMF))

(pb <- parboot(fm))
plot(pb)}
