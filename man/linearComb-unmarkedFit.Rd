\name{linearComb-unmarkedFit}
\title{Compute linear combinations of estimates in unmarkedFit objects.}
\description{Compute linear combinations of estimates in unmarkedFit objects.}
\details{This function computes the linear combination of parameter estimates in
\code{obj} given by the coefficient vector.  The user must
specify whether the detection or state covariates should be considered via the
\code{whichEstimate} argument.  This may be useful to estimate
quantities of interest from a fitted model or to test a hypothesis.}
\alias{linearComb,unmarkedFit-method}
\value{an unmarkedEstimate object}
\arguments{\item{obj}{an unmarkedFit object}
\item{coefficients}{vector of same length as obj}
\item{type}{character, either "state" or "det"}}
