\name{markovMN}
\alias{markovMN}
\title{ Fit Dynamic Multinomial Models }
\description{
%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

markovMN(stateformula = ~1, detformula = ~1, data = list(y = y,
covdata.site = covdata.site, covdata.obs = covdata.obs), y,
covdata.site = NULL, covdata.obs = NULL,detconstraint = NULL,
 phiconstraint = NULL, J, K)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stateformula}{ ~~Describe \code{stateformula} here~~ }
  \item{detformula}{ ~~Describe \code{detformula} here~~ }
  \item{data}{ provide the data}
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{covdata.site}{ ~~Describe \code{sitedata} here~~ }
  \item{covdata.obs}{ ~~Describe \code{obsdata} here~~ }
  \item{detconstraint}{ ~~Describe \code{detconstraint} here~~ }
  \item{phiconstraint}{ ~~Describe \code{phiconstraint} here~~ }
  \item{J}{ ~~Describe \code{J} here~~ }
  \item{K}{ Observations come from \{0,1,\ldots,K\}}
}
\details{
 % ~~ If necessary, more details than the description above ~~
}
\value{
  %~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Ian Fiske}
%\notes{ ~~further notes~~ 

 %~Make other sections like Warning with \section{Warning }{....} ~
%}
%~\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
