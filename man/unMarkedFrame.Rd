\name{unMarkedFrame}
\title{Create an unMarkedFrame.}
\usage{unMarkedFrame(y, siteCovs, obsCovs, obsNum=ncol(y), primaryNum)}
\description{Constructor for unMarkedFrames.}
\details{unMarkedFrame is the S4 class that holds data structures to be passed to the model-fitting functions in unMarked.

An unMarkedFrame contains the observations (\code{y}), covariates measured at the observation level (\code{obsCovs}), and covariates measured at the site level (\code{siteCovs}).
For a data set with M sites and J observations at each site, y is an M x J matrix.
\code{obsCovs} and \code{siteCovs} are both data frames (see \link{data.frame}).  \code{siteCovs} has M rows so that each row contains the covariates for the corresponding sites.
\code{obsCovs} has M*obsNum rows so that each covariates is ordered by site first, then observation number.  Missing values are coded with \code{NA} in any of y, siteCovs, or obsCovs.

Additionally, unMarkedFrames contain metadata, obsNum and primaryNum.  obsNum is the number of observations measured at each site. primaryNum is the number of seasons in a robust design sampling scheme.
Typically, these can be automatically determined by the constructor.  If not specified, obsNum is taken to be the number of columns in y and primaryNum is taken to be 1.
However, for certain situations, these must be supplied.  For example, double observer sampling, y has 3 columns corresponding the observer 1, observer 2, and both, but there were only two independent observations.
In this situation, y has 3 columns, but obsNum must be specified as 2.  This flexibility is currenty only used in the function \link{multinomPois}.

For convenience, \code{obsCovs} can be a list of M x obsNum matrices, with each one corresponding to an observation level covariate.

All site-level covariates are automatically copied to obsCovs so that site level covariates are available at the observation level.}
\alias{unMarkedFrame}
\alias{obsCovs}
\alias{siteCovs}
\value{an unMarkedFrame object}
\arguments{\item{y}{A matrix of the observed measured data.}
\item{obsCovs}{Dataframe of covariates that vary within sites.}
\item{siteCovs}{Dataframe of covariates that vary at the site level.}
\item{obsNum}{Number of independent observations.}
\item{primaryNum}{Number of primary time periods (seasons in the multiseason model).}}
\examples{data(mallard)
mallardUMF <- unMarkedFrame(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)
obsCovs(mallardUMF)
obsCovs(mallardUMF, matrices = TRUE)}
