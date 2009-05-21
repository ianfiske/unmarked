\name{backTransform}
\alias{backTransform}
\title{Transform an unmarkedEstimate object to it's natural scale.}
\description{Transform an unmarkedEstimate object to it's natural scale.}
\details{The transformation is determined by the invlink and invlinkGrad slots
in \code{obj}.  These slots specify the name of a one-to-one function and
its gradient respectively.  The delta method is used to compute the transformed
estimate.}
\value{an unmarkedEstimate object representing the transformed estimate.
This object has invlink and invlinkGrad slots as the identity function.}
\arguments{\item{unmarkedEstimate}{object to be transformed}}
