\name{unmarkedFrame}
\title{Create an unmarkedFrame, or one of its child classes.}
\alias{unmarkedFrame}
\alias{obsCovs}
\alias{siteCovs}
\usage{unmarkedFrame(y, siteCovs=NULL, obsCovs=NULL, mapInfo, plotArea, obsToY)}
\description{Constructor for unmarkedFrames.}
\arguments{
\item{y}{A matrix of the observed measured data.}
\item{siteCovs}{\code{\link{data.frame}} of covariates that vary at the site 
level.}
\item{obsCovs}{list of \code{\link{data.frame}}s of covariates that vary within
sites.}
\item{obsToY}{optional matrix specifying relationship between observation-level
covariates and response matrix}
\item{mapInfo}{geographic coordinate information}
\item{plotArea}{Vector of length M of plot areas. If specified, density rather 
than abundance will be modeled...for abundance models. Ignored by occupancy 
models.}
}
\details{unmarkedFrame is the S4 class that holds data structures to be passed 
to the model-fitting functions in unmarked.

An unmarkedFrame contains the observations (\code{y}), covariates measured at 
the observation level (\code{obsCovs}), and covariates measured at the site 
level (\code{siteCovs}).
For a data set with M sites and J observations at each site, y is an M x J matrix.
\code{obsCovs} and \code{siteCovs} are both data frames (see \link{data.frame}).  
\code{siteCovs} has M rows so that each row contains the covariates for the
corresponding sites.
\code{obsCovs} has M*obsNum rows so that each covariates is ordered by site first,
then observation number.  Missing values are coded with \code{NA} in any of y, 
siteCovs, or obsCovs.

Additionally, unmarkedFrames contain metadata: obsToY, plotArea, mapInfo. 
ObsToY is a matrix describing relationship between response matrix and 
observation-level covariates. Generally this does not need to be supplied by
the user; however, it may be needed when using \code{\link{multinomPois}}. For
example, double observer sampling, y has 3 columns corresponding the observer 1,
observer 2, and both, but there were only two independent observations.
In this situation, y has 3 columns, but obsToY must be specified.

Several child classes of \code{unmarkedFrame} require addional metadata. For
example, \code{unmarkedFrameDS} is used to organize distsance sampling data
for the \code{\link{distsamp}} function, and it has arguments dist.breaks,
tlength, survey, and unitsIn, which specify the distance interval cut points,
transect lengths, "line" or "point" transect, and units of measure, respectively.

All site-level covariates are automatically copied to obsCovs so that site level
 covariates are available at the observation level.
 }
\value{an unmarkedFrame object, see \code{\link{unmarkedFrame-class}}}
\examples{

# Set up data for pcount()
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
	obsCovs = mallard.obs)
summary(mallardUMF)


# Set up data for occu()
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)


# Set up data for distsamp()
data(linetran)
ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
	siteCovs = data.frame(Length, area, habitat), 
	dist.breaks = c(0, 5, 10, 15, 20),
	tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
	})
summary(ltUMF)


# Set up data for multinomPois()
data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
	siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])), 
	type = "removal")
summary(ovenFrame)	


# Set up data for colext()	
data(frogs)
frogUMF <- formatMult(masspcru)	
summary(frogUMF)

}
