\name{gf}
\alias{gf}
\alias{gf.data}
\alias{gf.obs}
\docType{data}
\title{ Green frog count index data}
\description{
Multinomial calling index data.
}
\usage{data(gf)}
\format{
A list with 2 components
	\describe{
	\item{gf.data}{220 x 3 matrix of count indices}
	\item{gf.obs}{list of covariates}
	}
}
\references{
Royle, J. Andrew, and William A. Link. 2005. A General Class of Multinomial
 Mixture Models for Anuran Calling Survey Data. Ecology 86, no. 9: 2505--2512.}
\examples{
data(gf)
str(gf.data)
str(gf.obs)
}
\keyword{datasets}
