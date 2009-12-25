\name{ovendata}
\alias{ovendata}
\alias{ovendata.list}
\docType{data}
\title{
Removal data for the Ovenbird
}
\description{
Removal sampling data collected for the Ovenbird (Seiurus aurocapillus).
}
\usage{data(ovendata)}
\format{
  The format is:
  chr "ovendata.list" which consists of
  	\describe{
	\item{data}{matrix of removal counts}
	\item{covariates}{data frame of site-level covariates}
	}

}
\source{
J.A. Royle (see reference below)
}
\references{
Royle, J. A. (2004). Generalized estimators of avian abundance from count survey data. Animal Biodiversity and Conservation, 27(1), 375-386.
}
\examples{
data(ovendata)
str(ovendata.list)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])), type = "removal")
}
\keyword{datasets}
