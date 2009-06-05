\name{obsCovs}
\alias{obsCovs}
\title{Extractor for observation level covariates...}
\usage{obsCovs(umf, matrices=FALSE)}
\description{Extractor for observation level covariates}
\value{either a data frame (default) or a list of matrices (if matrices = TRUE).}
\arguments{\item{umf}{an unmarkedFrame}
\item{matrices}{logical indicating whether to return the M * obsNum row data frame (default)
or a list of M x obsNum matrices (matrices = TRUE).}}
