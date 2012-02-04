\name{unmarkedRanef1-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedRanef1-class}
\alias{postMean}
\alias{postMode}
\alias{postMean,unmarkedRanef1-method}
\alias{postMode,unmarkedRanef1-method}
\alias{confint,unmarkedRanef1-method}
\alias{plot,unmarkedRanef1,missing-method}
\alias{show,unmarkedRanef1-method}

\title{Class \code{"unmarkedRanef1"}}

\description{
Stores the empirical Bayes estimates of the posterior distributions of
the latent abundance or occurrence variables.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{ranef}}.
}

\section{Slots}{
  \describe{
    \item{\code{bup}:}{A \code{\link{array}} with nSites rows and Nmax
  (K+1) columns and nPrimaryPeriod slices }
  }
}

\section{Methods}{
  \describe{
    \item{postMode}{\code{signature(object = "unmarkedRanef1")}: Extract the
    posterior modes}
    \item{postMean}{\code{signature(object = "unmarkedRanef1")}: Extract the
    posterior means}
    \item{confint}{\code{signature(object = "unmarkedRanef1")}: Compute
      confidence intervals.}
    \item{plot}{\code{signature(x = "unmarkedRanef1", y = "missing")}:
      Plot the posteriors using \code{\link{xyplot}}}
    \item{show}{\code{signature(object = "unmarkedRanef1")}: Display the
    modes and confidence intervals}
}
}

\seealso{
\code{\link{ranef}}
}
\examples{
showClass("unmarkedRanef1")
}
\keyword{classes}
