\name{unmarked}
\alias{unmarked-package}
\alias{unmarked}
\docType{package}
\title{
Models for Data from Unmarked Animals
}
\description{
  
This package fits several integrated likelihood models for
abundance and occurrence.  These models incorporate covariates on both the detection and state (occurrence or abundance)
process.

}
\details{
\tabular{ll}{
Package: \tab unmarked\cr
Type: \tab Package\cr
Version: \tab 0.1-4\cr
Date: \tab 2007-07-24\cr
License: \tab GPL version 2 or newer\cr
}

\strong{Overview of Model-fitting Functions:}  Unmarked provides
several functions for fitting integrated likelihood models for wildlife
abundance and occurrence to replicated survey data.  \command{\link{occu}}
fits occurrence models with no linkage between abundance and detection
(MacKenzie et al. 2006).  \command{\link{occuRN}} fits abundance models to
presence/absence data by exploiting the link between detection
probability and abundance (Royle and Nichols 2003).
\command{\link{pcount}} fits N-mixture models for repeated point count data
(Royle 2004, Kery and Royle 2005).  \command{\link{mnMix}} fits the
multinomial mixture model to categorical observations (Royle and Link
2005).  All of these functions allow the user to specify covariates that
affect the detection process and several also allow covariates for the
state process.

\strong{Data:}  All model-fitting functions analyze data that has
been collected from M sites with up to J samples taken at each site,
assuming closure within sites.  The only exception is \command{markovMN}
which uses data from multiple time periods, explicitly modeling the lack
of closure as a dynamic process.  See the documentation for
\command{\link{markovMN}} for its data requirements.  In general, survey
data consist of two fundamental components: \dfn{observation} and
\dfn{covariate} data.

\dfn{Observation} data are the observed responses, such as 0's and 1's
for occurrence data or non-negative integers for point count data.
Observation data is passed to unmarked functions as the \code{y} argument.
For a study with \eqn{M} sites and at most \eqn{J} observations per
site, \code{y} is an M x J matrix, where the \eqn{[i,j]^{th}} element is
the \eqn{j^{th}} observation for the \eqn{i^{th}} site.

\dfn{Covariate} data may be specified at either the observation level or
at the site level.  Observation-level covariates are passed to unmarked
functions as the \code{covdata.obs} argument.  \code{covdata.obs} may
either be a list containing M x J matrices for each covariate, or an M x
J x (number of covariates) three-dimensional array.  In either form,
\code{covdata.obs} must appropriately name its covariates.  In the list
format, each matrix must have a name that matches those used in
\code{detformula}.  In the array format, the dimnames for the 3rd
dimension must similarly match references in \code{detformula}.
Site-level covariates may be passed as the \code{covdata.site} argument
or with the observation-level covariates in \code{covdata.obs}.  In the
former form, \code{covdata.site} is a conventional data frame that
stores covariate data at the site level.  In the latter form, each
site-level covariate may be represented by an M x J matrix where each
column is identical.  This matrix can part of \code{covdata.obs} and
appropriately named.

\strong{Model Specification:}  Most of \emph{unmarked}'s
model-fitting functions allow specification of covariates for both the
state process and the detection process.  Covariates for the state
process (at the site level) and the detection process (at the site or
observation level) are specified with the right-hand sided formulas
passed to \code{siteformula} and \code{detformula} respectively.  Such a
formula looks like \eqn{~ x_1 + x_2 + \ldots + x_n}{~ x1 + x2 + ... +
x3} where \eqn{x_1}{x1} through \eqn{x_n}{xn} are additive covariates of
the process of interest.  The meaning of these covariates or, what they
model, is full described in the help files for the individual functions
and is not the same for all functions.

\strong{Utility Functions:}  \emph{Unmarked} contains several utility
functions for massaging data into the form required by its model-fitting
functions.  \command{\link{csvToData}} converts an appropriately
formated comma-separated values (.csv) file to a list containing the
components required by model-fitting functions.

}
\author{
Ian Fiske <ianfiske@gmail.com>
}
\references{

MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.
Amsterdam: Academic Press.
  
Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
pp. 777--790.

Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Kery, M. and Royle, J. A. (2005) Modeling Avaian Abundance from
Replicated Counts Using Binomial Mixture Models. \emph{Ecological
  Applications} 15(4), pp. 1450--1461.

Royle, J. A. and Link W. A. (2005) A general class of multinomial
mixture models for anuran calling survey data. \emph{Ecology}, \bold{86}(9),
pp. 2505--2512.

}
\keyword{package}
\seealso{

  \command{\link{occu}}, \command{\link{occuRN}},
  \command{\link{pcount}}, \command{\link{mnMix}},
  \command{\link{csvToData}}

}