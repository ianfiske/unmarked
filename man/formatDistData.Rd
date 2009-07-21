\name{formatDistData}
\alias{formatDistData}
\title{Convert individual-level distance data to the transect-level format required...}
\usage{formatDistData(distData, distCol, transectNameCol, cutPoints)}
\description{Convert individual-level distance data to the transect-level format required 
by distsamp()}
\details{This function creates a site (M) by distance interval (J) response matrix
from a data.frame containing the detection distances for each individual and
the transect names.}
\note{It is very important that the factor containing transect names contains
levels for all the transects surveyed. This includes those where birds were
not detected. See the example for how to add levels to a factor.}
\seealso{unmarkedFrameDS distsamp}
\arguments{\item{distData}{data.frame where each row is a detected individual. Must have
at least 2 columns. One for distances and the other for transect names.}
\item{distCol}{character, the column name containing distances}
\item{transectNameCol}{character, column name containing transect names}
\item{cutPoints}{numeric vector of distance interval cutpoints. Length must
equal J+1.}}
\examples{# Create a data.frame containing distances of animals detected
# along 4 transects.
dat <- data.frame(transect=gl(4,5, labels=letters[1:4]), 
distance=rpois(20, 10))
dat

# Look at your transect names.
levels(dat$transect)

# Suppose that you also surveyed a transect named "e" where no animals were
# detected. You must add it to the levels of dat$transect
levels(dat$transect) <- c(levels(dat$transect), "e")
levels(dat$transect)

# Distance cut points defining distance intervals
cp <- c(6, 8, 10, 12, 14, 18)

# Create formated response data.frame
yDat <- formatDistData(dat, "distance", "transect", cp) 
yDat

# Now you could merge yDat with transect-level covariates and 
# then use unmarkedFrameDS to prepare data for distsamp}
