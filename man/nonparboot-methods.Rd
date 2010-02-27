\name{nonparboot-methods}
\docType{methods}
\alias{nonparboot}
\alias{nonparboot-methods}
\alias{nonparboot,unmarkedFit-method}
\alias{nonparboot,unmarkedFitColExt-method}
\alias{nonparboot,unmarkedFitDS-method}
\alias{nonparboot,unmarkedFitMPois-method}
\alias{nonparboot,unmarkedFitOccu-method}
\alias{nonparboot,unmarkedFitOccuRN-method}
\alias{nonparboot,unmarkedFitPCount-method}
\title{ Nonparametric bootstrapping in unmarked }
\description{
Call \code{nonparboot} on an unmarkedFit to obtain non-parametric
bootstrap samples.  These can then be used by \code{vcov} in order to
get bootstrap estimates of standard errors.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFit")}}{ Obtain nonparametric
  bootstrap samples for a general unmarkedFit.  }

\item{\code{signature(object = "unmarkedFitColExt")}}{ Obtain nonparametric
  bootstrap samples for colext fits.  }

\item{\code{signature(object = "unmarkedFitDS")}}{ Obtain nonparametric
  bootstrap samples for a distsamp fits. }

\item{\code{signature(object = "unmarkedFitMPois")}}{ Obtain nonparametric
  bootstrap samples for a distsamp fits. }

\item{\code{signature(object = "unmarkedFitOccu")}}{  Obtain nonparametric
  bootstrap samples for a occu fits. }

\item{\code{signature(object = "unmarkedFitOccuRN")}}{  Obtain nonparametric
  bootstrap samples for a occuRN fits. }

\item{\code{signature(object = "unmarkedFitPCount")}}{  Obtain nonparametric
  bootstrap samples for a pcount fits. }
}}
\details{
  Calling \code{nonparboot} on an unmarkedFit returns the original
  unmarkedFit containing the bootstrap samples.  Subsequent calls to
  \code{vcov} with the argument \code{method="nonparboot"} will
  then use these bootstrap samples.
}
\examples{

}
\keyword{methods}

