\name{linetran}
\Rdversion{1.1}
\alias{linetran}
\docType{data}
\title{
Simulated line transect data
}
\description{
Response matrix of animals detected in four distance classes plus transect 
	lengths and two covariates.
}
\usage{data(linetran)}
\format{
  A data frame with 12 observations on the following 7 variables.
  \describe{
    \item{\code{dc1}}{Counts in distance class 1 [0-5 m)}
    \item{\code{dc2}}{Counts in distance class 2 [5-10 m)}
    \item{\code{dc3}}{Counts in distance class 3 [10-15 m)}
    \item{\code{dc4}}{Counts in distance class 4 [15-20 m)}
    \item{\code{Length}}{Transect lengths in km}
    \item{\code{area}}{Numeric covariate}
    \item{\code{habitat}}{a factor with levels \code{A} and \code{B}}
  }
}
\examples{
data(linetran)
linetran

# Format for distsamp()
ltUMF <- with(linetran, {
        unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
        siteCovs = data.frame(Length, area, habitat), 
        dist.breaks = c(0, 5, 10, 15, 20),
        tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
        })

}
\keyword{datasets}
