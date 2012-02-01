\name{unmarkedRanef1-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedRanef1-class}
\alias{coef,unmarkedRanef1-method}
\alias{confint,unmarkedRanef1-method}
\alias{plot,unmarkedRanef1,missing-method}
\alias{show,unmarkedRanef1-method}

\title{Class \code{"unmarkedRanef1"}}

\description{
Stores the empirical Bayes estimates of the posterior distributions of
the latent abundance or occurrence variables. This class works for
closed population models fitted using functions such as
\code{\link{occu}}, \code{\link{pcount}}, \code{\link{multinomPois}},
and \code{\link{distsamp}}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{ranef}}.
}

\section{Slots}{
  \describe{
    \item{\code{bup}:}{A \code{"matrix"} with nSites rows and Nmax (K+1) columns }
    \item{\code{modes}:}{A vector of length nSites containing the
      posterior modes.}
  }
}

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "unmarkedRanef1")}: Extract the
    posterior modes}
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
