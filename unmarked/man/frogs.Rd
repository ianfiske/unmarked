\name{frogs}
\alias{frogs}
\alias{pcru.data}
\alias{pcru.bin}
\alias{pcru.y}
\alias{pfer.data}
\alias{pfer.bin}
\alias{pfer.y}
\alias{frog2001pfer}
\alias{frog2001pcru}
%\alias{delaware_bfow_2001-2006}
\docType{data}
\title{ 2001 Delaware North American Amphibian Monitoring Program Data }
\description{
frogs contains NAAMP data for Pseudacris feriarum (pfer) and Pseudacris
crucifer (pcru) in 2001.
}
\usage{data(frogs)}
\format{
	\describe{
	\item{pcru.y}{matrix of observed calling indices for pcru}
	\item{pcru.bin}{matrix of detections for pcru}
	\item{pcru.data}{array of covariates measured at the
	  observation-level for pcru}
	\item{pfer.y}{matrix of observed calling indices for pfer}
	\item{pfer.bin}{matrix of detections for pfer}
	\item{pfer.data}{array of covariates measured at the
	  observation-level for pfer}
	}
}
\details{
The rows of pcru.y, pcru.bin, pfer.y, and pfer.bin correspond to sites
and columns correspond to visits to each site.  The first 2 dimensions
of pfer.data and pcru.data are matrices of covariates that correspond to
the observation matrices (sites \eqn{\times} observation), with the 3rd
dimension corresponding to separate covariates.
}
\source{
https://www.pwrc.usgs.gov/naamp/
}
\references{
Mossman MJ, Weir LA. North American Amphibian Monitoring Program (NAAMP). Amphibian Declines: the conservation status of United States species. University of California Press, Berkeley, California, USA. 2005:307-313.
}
\examples{
data(frogs)
str(pcru.data)
}
\keyword{datasets}
