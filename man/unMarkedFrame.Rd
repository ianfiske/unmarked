\name{unMarkedFrame}
\alias{unMarkedFrame}
\title{Constuctor function to create an unmarkedFrame.}
\usage{unMarkedFrame(y, siteCovs, obsCovs, obsNum=ncol(y))}
\description{Constuctor function to create an unmarkedFrame.}
\details{This function takes observations (y) and covariates at the site level
(siteCovs) and observations level (obsCovs) and generates an
unMarkedFrame object.  This object can then be passed to any of the
unMarked functions.}
\value{an unmarkedFrame object}
\arguments{\item{y}{A matrix of the observed measured data.}
\item{obsCovs}{Dataframe of covariates that vary within sites.}
\item{siteCovs}{Dataframe of covariates that vary at the site level.}
\item{obsNum}{Number of independent observations.}}
