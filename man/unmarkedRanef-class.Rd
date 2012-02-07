\name{unmarkedRanef-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedRanef-class}
\alias{bup}
\alias{bup,unmarkedRanef-method}
\alias{confint,unmarkedRanef-method}
\alias{plot,unmarkedRanef,missing-method}
\alias{show,unmarkedRanef-method}

\title{Class \code{"unmarkedRanef"}}

\description{
Stores the estimated posterior distributions of
the latent abundance or occurrence variables.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{ranef}}.
}

\section{Slots}{
  \describe{
    \item{\code{post}:}{A \code{\link{array}} with nSites rows and Nmax
  (K+1) columns and nPrimaryPeriod slices }
  }
}

\section{Methods}{
  \describe{
    \item{bup}{\code{signature(object = "unmarkedRanef")}: Extract the
    posterior means or modes}
    \item{confint}{\code{signature(object = "unmarkedRanef")}: Compute
      confidence intervals.}
    \item{plot}{\code{signature(x = "unmarkedRanef", y = "missing")}:
      Plot the posteriors using \code{\link{xyplot}}}
    \item{show}{\code{signature(object = "unmarkedRanef")}: Display the
    modes and confidence intervals}
}
}

\seealso{
\code{\link{ranef}}
}
\examples{
showClass("unmarkedRanef")
}
\keyword{classes}
