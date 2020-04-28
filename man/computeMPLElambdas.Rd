\name{computeMPLElambda}

\alias{computeMPLElambda}

\title{Compute the penalty weight for the MPLE penalized likelihood method}

\usage{computeMPLElambda(formula, data, knownOcc=numeric(0), starts,
method="BFGS",engine=c("C", "R"))}

\arguments{
    \item{formula}{Double right-hand side formula describing covariates of
        detection and occupancy in that order.}
    \item{data}{An \code{\link{unmarkedFrameOccu}} object}
    \item{knownOcc}{Vector of sites that are known to be occupied. These
    should be supplied as row numbers of the y matrix, eg, c(3,8) if
    sites 3 and 8 were known to be occupied a priori.}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{engine}{Either "C" or "R" to use fast C++ code or native R
      code during the optimization.}
  }

\description{This function computes the weight for the MPLE penalty of Moreno & Lele (2010). }

\details{

See \code{\link{occuPEN}} for details and examples.

}

\value{The computed lambda.}

\references{

Moreno, M. and S. R. Lele. 2010. Improved estimation of site occupancy
using penalized likelihood. Ecology 91: 341-346.

}

\author{Rebecca A. Hutchinson}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccu}},
    \code{\link{occu}}, \code{\link{occuPEN}}, \code{\link{occuPEN_CV}}, 
    \code{\link{nonparboot}}}


\keyword{models}
