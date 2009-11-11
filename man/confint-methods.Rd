\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\alias{confint,unmarkedBackTrans-method}
\alias{confint,unmarkedEstimate-method}
\alias{confint,unmarkedFit-method}
\alias{confint,unmarkedLinComb-method}
\title{Methods for Function confint in Package `unmarked'}
\description{
Methods for function \code{confint} in Package `unmarked'
}
\usage{
\S4method{confint}{unmarkedBackTrans}(object, parm, level)
\S4method{confint}{unmarkedEstimate}(object, parm, level)
\S4method{confint}{unmarkedLinComb}(object, parm, level)
\S4method{confint}{unmarkedFit}(object, parm, level, type, method)
}
\arguments{
	\item{object}{Object of appropriate S4 class}
	\item{parm}{Name of parameter(s) of interest}
	\item{level}{Level of confidence}
	\item{type}{Either "state" or "det"}
	\item{method}{Either "normal" or "profile"}
	}
\value{A vector of lower and upper confidence intervals. These are
asymtotic unless method='profile' is used on unmarkedFit objects in which case
they are profile likelihood intervals.}
\seealso{\code{\link{unmarkedFit-class}}}
\keyword{methods}

