\name{lambda2psi}
\Rdversion{1.1}
\alias{lambda2psi}
\title{
Convert Poisson mean (lambda) to probability of occurrence (psi).
}
\description{
Abundance and occurrence are fundamentally related.
}
\usage{
lambda2psi(lambda)
}
\arguments{
  \item{lambda}{
	Numeric vector with values >= 0
}
}
\value{
	A vector of psi values of the same length as lambda.
}
\seealso{
\code{\link{pcount}}, \code{\link{multinomPois}}, \code{\link{distsamp}}
}
\examples{
lambda2psi(0:5)
}
