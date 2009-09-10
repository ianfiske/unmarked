\name{umHMM}
\alias{umHMM}
\title{Fit the general multistate multiseason occupancy model.}
\usage{umHMM(formula=~1, umf, detconstraint, phiconstraint, psiconstraint,
    phiMat, inits, B=0, max.bad=3, min.good=1, fit.stats=TRUE,
    method=BFGS, control=list())}
\description{Estimate parameters of the general multiseason multistate occpancy model.}
\details{Site level covariates are currently not implemented for markovMN and so stateformula is ignored.
See \link{unmarked} for a discussion of detformula.  See \link{unmarkedFrame} for a description of how to create
an unmarkedFrame for supplying data to the argument \code{umf}.

\code{K} and \code{phiMatrix} together determine the model form.  Multiple phi matrices
may be possible for each K.  Options for phi are:


TODO:  describe phi matrices in detail.

Currently, selection of the appropriate detection matrix for a given \code{K}
is given by \code{arDet}.  If \code{arDet} is TRUE, then the autoregressive flavored
detection matrix (reduced form) is chose.

TODO:  describe detection matrices in detail.


Each freely varying detection parameter can be modeled as a linear function of observation level covariates.
\code{detconstraint} is a matrix with number of rows equal to the detection matrix parameters and number of
columns equal to the number of covariates (plus 1 for the intercept).  Each column is a constraint vector beginning with 1 that indicates
which covariates are restricted to have the same effect on a given matrix parameter.}
\arguments{\item{formula}{right-hand side formula describing covariates of detection.}
\item{umf}{unmarkedFrame object that supplies the data (see \link{unmarkedFrame})..}
\item{detconstraint}{matrix to describe constraints on detection parameters}
\item{phiconstraint}{vector to describe phi constraints}
\item{psiconstraint}{vector to describe psi constraints}
#\item{K}{integer maximum level of y}
\item{phiMat}{character vector describing which phi specification to use (see Details)}
\item{inits}{optionally, initial values for parameters in the optimization}
\item{B}{number of bootstrap interations (0 indicates none)}
\item{min.good}{minimum number of successful model fits before returning}
\item{max.bad}{maximum number of unsuccessful model fits before aborting}
\item{fit.stats}{logical indicating whether or not to return the fit statistics}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}
}