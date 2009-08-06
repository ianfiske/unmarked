\name{calcAreas}
\alias{calcAreas}
\title{Prepare area argument for distsamp().}
\usage{calcAreas(dist.breaks, tlength, survey, output, M, J, unitsIn,
    unitsOut)}
\description{Prepare area argument for distsamp(). This is primarily for internal use
but see details. Caution should be used because the returned matrix has 
different interpretations for different survey and output types.}
\value{An M x J numeric matrix. 

If output == "density" and survey == "line"
then the values are the areas of each distance class for each transect. If
output == "density" and survey == "point" then the values are the the radii
of each point transect. Currently, radii cannot vary.

If survey == "point" and output == "abund" a matrix of 1s is returned. If
survey == "line" and output == "abund" a matrix of transect lengths is
returned because transect lengths must be taken into account even if density
is not of interest.}
\note{This function might be useful if some distance classes for some 
transects were not surveyed. In such a case, missing values could be added
to the output of calcAreas() and the modified matrix could be supplied
to the plotArea argument of unmarkedFrameDS.}
\arguments{\item{dist.breaks}{numeric vector of distance class break poings}
\item{tlength}{numeric vector of transect lengths for line transects}
\item{M}{number of transects}
\item{J}{number of distance classes}
\item{survey}{either "line" or "point"}
\item{output}{either "abund" or "density"}
\item{unitsIn}{either "m" or "km" for units of both dist.breaks and tlength.}
\item{unitsOut}{either "ha" or "kmsq"}}
\examples{data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 
lengths <- linetran$Length * 1000

calcAreas(dbreaksLine, lengths, "line", "density", M=nrow(linetran), 
J = length(dbreaksLine) - 1, "m", "ha")

data(pointtran)
(dbreaksPt <- seq(0, 25, by=5))

calcAreas(dbreaksPt, survey="point", output="density", M=nrow(pointtran),
J = length(dbreaksPt) - 1, unitsIn = "m", unitsOut = "ha")}
