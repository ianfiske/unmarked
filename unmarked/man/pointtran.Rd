\name{pointtran}
\Rdversion{1.1}
\alias{pointtran}
\docType{data}
\title{
Simulated point-transect data
}
\description{
Response matrix of animals detected in five distance classes plus two covariates.
}
\usage{data(pointtran)}
\format{
  A data frame with 30 observations on the following 7 variables.
  \describe{
    \item{\code{dc1}}{Counts in distance class 1 [0-5 m)}
    \item{\code{dc2}}{Counts in distance class 2 [5-10 m)}
    \item{\code{dc3}}{Counts in distance class 3 [10-15 m)}
    \item{\code{dc4}}{Counts in distance class 4 [15-20 m)}
    \item{\code{dc5}}{Counts in distance class 5 [20-25 m)}
    \item{\code{area}}{a numeric vector}
    \item{\code{habitat}}{a factor with levels \code{A} \code{B} \code{C}}
  }
}
\examples{
data(pointtran)
pointtran

# Format for distsamp()
ptUMF <- with(pointtran, {
        unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4, dc5), 
        siteCovs = data.frame(area, habitat), 
        dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
        })
}
\keyword{datasets}
