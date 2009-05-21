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
(fm <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~habitat, linetran, 
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m"))

(pb <- parboot(fm))
plot(pb)}
