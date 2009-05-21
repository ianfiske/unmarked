\name{removalPiFun}
\alias{removalPiFun}
\title{Compute the cell probabilities for the observation classes...}
\usage{removalPiFun(p)}
\description{Compute the cell probabilities for the observation classes
in removal sampling.}
\details{Both p and the returned matrix are M x J for M sites and J sampling occasions.

This function and \link{doublePiFun} provide example functions for computing multinomial
cell probabilities for \link{multinomPois}.}
\value{matrix of cell probabilties for multinomial cells in removal sampling observation classes.}
\arguments{\item{p}{matrix of detection probabilities at each occasion for each site}}
