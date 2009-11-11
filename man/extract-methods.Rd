\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{[,unmarkedEstimateList,ANY,ANY,ANY-method}
\alias{[,unmarkedFit,ANY,ANY,ANY-method}
\alias{[,unmarkedFrame,missing,numeric,missing-method}
\alias{[,unmarkedFrame,numeric,missing,missing-method}
\alias{[,unmarkedFrame,numeric,numeric,missing-method}
\alias{[,unmarkedMultFrame,missing,numeric,missing-method}
\title{Methods for bracket extraction [ in Package `unmarked'}
\description{
Methods for bracket extraction \code{[} in Package `unmarked'
}
\usage{
\S4method{[}{unmarkedEstimateList,ANY,ANY,ANY}(x, i, j, drop)
\S4method{[}{unmarkedFit,ANY,ANY,ANY}(x, i, j, drop)
\S4method{[}{unmarkedFrame,numeric,numeric,missing}(x, i, j)
\S4method{[}{unmarkedMultFrame,missing,numeric,missing}(x, i, j)
}
\arguments{
	\item{x}{Object of appropriate S4 class}
	\item{i}{Row numbers}
	\item{j}{Observation numbers (eg occasions, distance classes, etc...)}
	\item{drop}{Not currently used}
	}
\section{Methods}{
\describe{

\item{x = "unmarkedEstimateList", i = "ANY", j = "ANY", drop = "ANY"}{Extract a
	unmarkedEstimate object from an unmarkedEstimateList by name (either 'det' 
	or 'state')}

\item{x = "unmarkedFit", i = "ANY", j = "ANY", drop = "ANY"}{Extract a
	unmarkedEstimate object from an unmarkedFit by name (either 'det' 
	or 'state')}

\item{x = "unmarkedFrame", i = "missing", j = "numeric", drop = "missing"}{
	Extract observations from an unmarkedFrame.}

\item{x = "unmarkedFrame", i = "numeric", j = "missing", drop = "missing"}{ 
	Extract rows from an unmarkedFrame }

\item{x = "unmarkedFrame", i = "numeric", j = "numeric", drop = "missing"}{ 
	Extract rows and observations from an unmarkedFrame }

\item{x = "unmarkedMultFrame", i = "missing", j = "numeric", drop = "missing"}{
	Extract primary sampling periods from an unmarkedMultFrame}
}}
\examples{

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
	obsCovs = mallard.obs)
summary(mallardUMF)

mallardUMF[1:5,]
mallardUMF[,1:2]
mallardUMF[1:5, 1:2]

}
\keyword{methods}

