\name{backTransform}
\alias{backTransform}
\title{Transform an unMarkedEstimate object to it's natural scale.}
\description{Transform an unMarkedEstimate object to it's natural scale.}
\details{The transformation is determined by the invlink and invlinkGrad slots
in \code{obj}.  These slots specify the name of a many-to-one function and
its gradient respectively.  The delta method is used to compute the transformed
estimate.}
\value{an unMarkedEstimate object representing the transformed estimate.
This object has invlink and invlinkGrad slots as the identity function.}
\arguments{\item{unMarkedEstimate}{object to be transformed}}
