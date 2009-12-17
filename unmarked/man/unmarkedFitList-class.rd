\name{unmarkedFitList-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedFitList-class}
\alias{modSel,unmarkedFitList-method}
\alias{predict,unmarkedFitList-method}

\title{Class "unmarkedFitList"}
\description{Class to hold multiple nested fitted models from one of 
	\code{unmarked}'s fitting functions}
\section{Objects from the Class}{
Objects can be created by using the \code{\link{fitList}} function. 
}
\section{Slots}{
	 \describe{
    \item{\code{fits}:}{A \code{"list"} of models.}
  }
}
\section{Methods}{
  \describe{
    \item{modSel}{\code{signature(object = "unmarkedFitList")}: 
		Model selection}
    \item{predict}{\code{signature(object = "unmarkedFitList")}: 
		Model-averaged prediction}
	 }
}
\note{Model-averaging regression coefficients is intentionally not implemented.}
\seealso{
	\code{\link{fitList}},
	\code{\linkS4class{unmarkedFit}}
}
\examples{
showClass("unmarkedFitList")
}
\keyword{classes}
