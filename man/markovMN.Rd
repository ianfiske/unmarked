\name{markovMN}
\alias{markovMN}
\title{This fits the general multistate multiseason occupancy model of Fiske and Royle...}
\usage{markovMN(stateformula=~1, detformula=~1, umf, detconstraint,
    phiconstraint, psiconstraint, K, phiMatrix=cumlogit, EM=FALSE,
    psi.init, phiParms.init, detParms.init, get.inits=TRUE,
    bootstrap.se=FALSE, B=50, trace=FALSE, arDet=FALSE)}
\description{This fits the general multistate multiseason occupancy model of Fiske and Royle}
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
