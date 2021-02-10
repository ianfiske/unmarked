\name{sigma}
\alias{sigma}
\alias{sigma-methods}
\alias{sigma,unmarkedFit-method}
\alias{sigma,unmarkedEstimate-method}

\title{Extract estimates of random effect standard deviations}

\description{
  Extract estimates and summary statistics of random effect standard
  deviations from an \code{unmarkedFit} model or an \code{unmarkedEstimate}.
}

\usage{
\S4method{sigma}{unmarkedEstimate}(object, level=0.95, ...)
\S4method{sigma}{unmarkedFit}(object, type, level=0.95, ...)
}

\arguments{
  \item{object}{An object inheriting class \code{unmarkedEstimate} or 
    \code{unmarkedFit}}
  \item{level}{Significance level to use for confidence interval}
  \item{type}{If provided, return only random effect SDs from the chosen
    submodel type (as a character string)}
  \item{...}{Other arguments}
} 

\value{\code{data.frame} containing estimates, SEs, and confidence intervals
 for random effect standard deviations in the model. 
}

\author{Ken Kellner \email{contact@kenkellner.com}}
