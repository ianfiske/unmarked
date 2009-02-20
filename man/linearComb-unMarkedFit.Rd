\name{linearComb-unMarkedFit}
\title{Compute linear combinations of estimates in unMarkedFit objects.}
\description{Compute linear combinations of estimates in unMarkedFit objects.}
\details{This function computes the linear combination of parameter estimates in
\code{obj} given by the coefficient vector.  The user must
specify whether the detection or state covariates should be considered via the
\code{whichEstimate} argument.  This may be useful to estimate
quantities of interest from a fitted model or to test a hypothesis.}
\alias{linearComb,unMarkedFit-method}
\value{an unMarkedEstimate object}
\arguments{\item{obj}{an unMarkedFit object}
\item{coefficients}{vector of same length as obj}
\item{whichEstimate}{character, either "state" or "det"}}
