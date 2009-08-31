\name{unmarkedMultFrame}
\alias{unmarkedMultFrame}
\title{Create an unmarkedMultFrame.}
\usage{unmarkedMultFrame(y, siteCovs, obsCovs, numPrimary, yearlySiteCovs,
    plotArea)}
\description{This function constructs an unmarkedMultFrame object.}
\details{unmarkedMultFrame objects are used in the fitting functions \code{\link{colext}} and \code{hmm}.

For a study with \emph{M} sites, \emph{T} years, and a maximum of \emph{J} observations per site-year, the data
are shaped as follows.  \code{y} is an \eqn{M \times TJ}{M by TJ} matrix, with each row corresponding to a site.
\code{siteCovs} is a data frame with \eqn{M} rows.
\code{yearlySiteCovs} is a data frame with \eqn{MT} rows which are in site-major, year-minor order.
\code{obsCovs} is a data frame with \eqn{MTJ} rows, which are ordered by site-year-observation, so that a 
column of \code{obsCovs} corresponds to \command{as.vector(t(y))}, element-by-element.  
The number of years must be specified in \code{numPrimary}.

If the data are in long format, the convenience function \code{\link{formatMult}} is useful for creating
the unmarkedMultFrame.}
\value{an unmarkedFrame object}
\arguments{\item{y}{A matrix of the observed data.}
\item{siteCovs}{Data frame of covariates that vary at the site level.}
\item{obsCovs}{Data frame of covariates that vary within site-year-observation level.}
\item{numPrimary}{Number of primary time periods (seasons in the multiseason model).}
\item{yearlySiteCovs}{Data frame containing covariates at the site-year level.}
\item{plotArea}{optional vector of plot areas.}}
