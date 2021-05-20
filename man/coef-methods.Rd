\name{coef-methods}
\docType{methods}
\alias{coef-methods}
\alias{coef,linCombOrBackTrans-method}
\alias{coef,unmarkedEstimate-method}
\alias{coef,unmarkedFit-method}
\title{Methods for Function coef in Package `unmarked'}
\description{
Extract coefficients
}
\usage{
\S4method{coef}{unmarkedFit}(object, type, altNames = TRUE, fixedOnly=TRUE)
\S4method{coef}{unmarkedEstimate}(object, altNames = TRUE, fixedOnly=TRUE, ...)
\S4method{coef}{linCombOrBackTrans}(object)
}
\arguments{
\item{object}{Object of appropriate S4 class}
\item{type}{Either 'state' or 'det'}
\item{altNames}{Return specific names for parameter estimates?}
\item{fixedOnly}{Return only fixed effect parameters?}
\item{...}{Further arguments. Not currently used}
}
\section{Methods}{
\describe{

\item{object = "linCombOrBackTrans"}{Object from linearComb}

\item{object = "unmarkedEstimate"}{unmarkedEstimate object}

\item{object = "unmarkedFit"}{Fitted model}
}}
\value{A named numeric vector of parameter estimates.}
\keyword{methods}

