\name{formatMult}
\alias{formatMult}
\title{Create unmarkedMultFrame from Long Format Data Frame}
\usage{formatMult(df.in)}
\description{This convenience function converts multi-year data in long format to unmarkedMultFrame Object.  See Details for more information.}
\details{\code{df.in} is a data frame with columns formatted as follows:

Column 1 = year number \cr
Column 2 = site name or number \cr
Column 3 = julian date or chronological sample number during year \cr
Column 4 = observations (y) \cr
Column 5 -- Final Column = covariates 

Note that if the data is already in wide format, it may be easier to create an unmarkedMultFrame object
directly with a call to \code{\link{unmarkedMultFrame}}.}
\value{unmarkedMultFrame object}
\arguments{\item{df.in}{a data.frame appropriately formatted (see Details).}}
