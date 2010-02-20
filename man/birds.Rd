\name{birds}
\alias{birds}
\alias{catbird}
\alias{catbird.bin}
\alias{woodthrush}
\alias{woodthrush.bin}
\docType{data}
\title{ BBS Point Count and Occurrence Data from 2 Bird Species}
\description{
  Data frames for 2 species from the breeding bird survey (BBS).  Each
  data frame has a row for each site and columns for each sampling
  event.  There is a point count and occurrence--designated by
  \code{.bin}--  version for each species.
}
\usage{data(birds)}
\format{
  \describe{
    \item{\code{catbird}}{A data frame of point count observations for
      the catbird.}
    \item{\code{catbird.bin}}{A data frame of occurrence observations for
      the catbird.}
    \item{\code{woodthrush}}{A data frame of point count observations for
      the wood thrush.}
    \item{\code{woodthrush.bin}}{A data frame of point count observations for
      the wood thrush.}
  }
}
%\details{
 % ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Royle J. N-mixture models for estimating population size from spatially
replicated counts. Biometrics. 2004.  60(1):108--115.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(birds)
}
\keyword{datasets}
