\name{csvToUMF}
\alias{csvToUMF}
\title{Convert .CSV File to an unmarkedFrame}
\usage{csvToUMF(filename, long=FALSE, species)}
\description{This function converts an appropriatedly formated comma-separated
values file (.csv) to a format usable by \emph{unmarked}'s fitting
functions (see \emph{Details}).}
\details{This function provides a quick way to take a .csv file with headers
named as described below and provides the data required and returns of
data in the format required by the model-fitting functions in
\code{\link{unmarked}}.  The .csv file can be in one of 2 formats: long or
wide.  See the first 2 lines of the \emph{examples} for what these
formats look like.

The .csv file is formatted as follows:
\itemize{
\item col 1 is site labels.
\item if data is in long format, col 2 is date of observation.
\item next J columns are the observations (y) - counts or 0/1's.
\item next is a series of columns for the site variables (one column
per variable).  The column header is the variable name.
\item next is a series of columns for the observation-level variables.
These are in sets of J columns for each variable, e.g., var1-1 var1-2
var1-3 var2-1 var2-2 var2-3, etc.  The column header of the first
variable in each group must indicate the variable name.
}}
\value{an unmarkedFrame object}
\keyword{utilities}
\author{Ian Fiske \email{ianfiske@gmail.com}}
\arguments{\item{filename}{string describing filename of file to read in}
\item{long}{\code{FALSE} if file is in long format or \code{TRUE} if
file is in long format (see \emph{Details})}
\item{species}{if data is in long format with multiple species, then
this can specify a particular species to extract if there is a
column named "species".}}
\examples{# examine a correctly formatted long .csv
read.csv(system.file("csv","frog2001pcru.csv", package="unmarked"))

# examine a correctly formatted wide .csv
read.csv(system.file("csv","widewt.csv", package="unmarked"))

# convert them!
dat1 <- csvToUMF(system.file("csv","frog2001pcru.csv", package="unmarked"), long = TRUE)
dat2 <- csvToUMF(system.file("csv","frog2001pfer.csv", package="unmarked"), long = TRUE)
dat3 <- csvToUMF(system.file("csv","widewt.csv", package="unmarked"), long = FALSE)}
