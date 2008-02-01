\name{Mh}
\alias{Mh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fits Model Mh -- a closed population model with heterogeneity }
\description{
 This fits the logit-normal version of Model Mh in which the logit of
individual detection probability is assumed to have a normal distribution with 
parameters  mu and sigma. These parameters and the unobserved number of
individuals, n0, are estimated by integrated likelihood.
}
\usage{
Mh(nx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nx}{ the detection frequency vector, of length J where J is the number of sample occasions }
}
\details{
  At this time there are no further details
}
\value{
% # ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Just the nlm output so far}
  \item{comp2 }{Description of 'comp2'}
}
\references{ Coull and Agresit (1999) Dorazio and Royle (2003) }
\author{ Andy Royle }
\note{ ~~further notes~~ 
 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ closed population }
%\keyword{ model Mh }
\keyword{models}
% __ONLY ONE__ keyword per line
