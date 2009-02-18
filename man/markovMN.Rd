\name{markovMN}
\alias{markovMN}
\title{Fit the general multistate multiseason occupancy model.}
\usage{markovMN(stateformula=~1, detformula=~1, umf, detconstraint,
    phiconstraint, psiconstraint, K, phiMatrix=cumlogit, EM=FALSE,
    psi.init, phiParms.init, detParms.init, get.inits=TRUE,
    bootstrap.se=FALSE, B=50, trace=FALSE, arDet=FALSE, phiMatFun,
    detVecFun, n.starts=1, homo.det=FALSE, max.bad=3, min.good=1)}
\description{Estimate parameters of the general multiseason multistate occpancy model.}
\details{Site level covariates are currently not implemented for markovMN and so stateformula is ignored.
See \link{unmarked} for a discussion of detformula.  See \link{unMarkedFrame} for a description of how to create
an unMarkedFrame for supplying data to the argument \code{umf}.

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
\arguments{\item{stateformula}{right-hand side formula describing covariates of occurence.}
\item{detformula}{right-hand side formula describing covariates of detection.}
\item{umf}{unMarkedFrame object that supplies the data (see \link{unMarkedFrame})..}
\item{detconstraint}{matrix to describe constraints on detection parameters}
\item{phiconstraint}{vector to describe phi constraints}
\item{psiconstraint}{vector to describe psi constraints}
\item{K}{integer maximum level of y}
\item{phiMatrix}{character describing which phi to use}
\item{EM}{logical specifying to use EM (TRUE) or BFGS (FALSE)}
\item{psi.init}{initial values for psi}
\item{phiParms.init}{initial values for phi parameters}
\item{detParms.init}{initial values for detection parameters}
\item{get.inits}{logical, specifying whether or not to search for better initial values}
\item{booststrap.se}{logical. TRUE calculates standard errors using bootstrap}
\item{B}{number of bootstrap interations}
\item{trace}{logical, TRUE show BFGS steps}
\item{arDet}{TRUE chooses the autoregressive type detection structure}}
