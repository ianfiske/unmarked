\name{backTransform-methods}
\docType{methods}
\alias{backTransform}
\alias{backTransform-methods}
\alias{backTransform,unmarkedEstimate-method}
\alias{backTransform,unmarkedFit-method}
\alias{backTransform,unmarkedLinComb-method}
\alias{show,unmarkedBackTrans-method}

\title{Methods for Function backTransform in Package `unmarked'}
\description{Methods for function \code{backTransform} in Package `unmarked'. 
	This converts from link-scale to original-scale}
\usage{
\S4method{backTransform}{unmarkedFit}(obj, type)
\S4method{backTransform}{unmarkedEstimate}(obj)
}
\arguments{
	\item{obj}{Object of appropriate S4 class}
	\item{type}{one of names(obj), eg 'state' or 'det'}
	}
\section{Methods}{
\describe{

\item{obj = "unmarkedEstimate"}{Typically done internally}

\item{obj = "unmarkedFit"}{Back-transform a parameter from a fitted model. Only
	possible if no covariates are present. Must specify argument type 
	as one of the values returned by names(obj).}

\item{obj = "unmarkedLinComb"}{Back-transform a predicted value created by 
	\code{linearComb}. This is done internally by \code{\link{predict}} but
	can be done explicitly by user.}
}}
\keyword{methods}
