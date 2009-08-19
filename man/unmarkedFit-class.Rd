\name{unmarkedFit-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedFit}

\title{Class "unmarkedFit" }
\description{ Contains fitted model information }
\section{Slots}{
	 \describe{
    \item{\code{fitType}:}{Object of class \code{"character"} }
    \item{\code{call}:}{Object of class \code{"call"} }
    \item{\code{formula}:}{Object of class \code{"formula"} }
    \item{\code{data}:}{Object of class \code{"unmarkedFrame"} }
    \item{\code{sitesRemoved}:}{Object of class \code{"numeric"} }
    \item{\code{estimates}:}{Object of class \code{"unmarkedEstimateList"} }
    \item{\code{AIC}:}{Object of class \code{"numeric"} }
    \item{\code{opt}:}{Object of class \code{"list"} containing results from \code{\link{optim}} }
    \item{\code{negLogLike}:}{Object of class \code{"numeric"} }
    \item{\code{nllFun}:}{Object of class \code{"function"} }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "unmarkedFit", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{backTransform}{\code{signature(obj = "unmarkedFit")}: ... }
    \item{coef}{\code{signature(object = "unmarkedFit")}: ... }
    \item{confint}{\code{signature(object = "unmarkedFit")}: ... }
    \item{fitted}{\code{signature(object = "unmarkedFit")}: ... }
    \item{getData}{\code{signature(object = "unmarkedFit")}: ... }
    \item{getP}{\code{signature(object = "unmarkedFit")}: ... }
    \item{hessian}{\code{signature(object = "unmarkedFit")}: ... }
    \item{linearComb}{\code{signature(obj = "unmarkedFit", coefficients = "matrixOrVector")}: ... }
    \item{mle}{\code{signature(object = "unmarkedFit")}: ... }
    \item{names}{\code{signature(x = "unmarkedFit")}: ... }
    \item{nllFun}{\code{signature(object = "unmarkedFit")}: ... }
    \item{parboot}{\code{signature(object = "unmarkedFit")}: ... }
    \item{plot}{\code{signature(x = "unmarkedFit", y = "missing")}: ... }
    \item{predict}{\code{signature(object = "unmarkedFit")}: ... }
    \item{profile}{\code{signature(fitted = "unmarkedFit")}: ... }
    \item{residuals}{\code{signature(object = "unmarkedFit")}: ... }
    \item{sampleSize}{\code{signature(object = "unmarkedFit")}: ... }
    \item{SE}{\code{signature(obj = "unmarkedFit")}: ... }
    \item{show}{\code{signature(object = "unmarkedFit")}: ... }
    \item{summary}{\code{signature(object = "unmarkedFit")}: ... }
    \item{update}{\code{signature(object = "unmarkedFit")}: ... }
    \item{vcov}{\code{signature(object = "unmarkedFit")}: ... }
	 }
}
\note{ This is a superclass with child classes for each fit type }
\examples{
showClass("unmarkedFit")
}
\keyword{classes}
