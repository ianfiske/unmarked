\name{imputeMissing}
\alias{imputeMissing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to impute missing entries in continuous obsCovs
}
\description{
This function uses an ad-hoc averaging approach to impute missing
entries in obsCovs.  The missing entry is replaced by an average of the
average for the site and the average for the visit number.

}
\usage{
imputeMissing(umf, whichCovs = seq(length=ncol(obsCovs(umf))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{umf}{
The data set who's obsCovs are being imputed.
}
  \item{whichCovs}{
    An integer vector giving the indices of the covariates to be imputed.
    This defaults to all covariates in \code{obsCovs}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A version of \code{umf} that has the requested obsCovs imputed.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Ian Fiske
}

\examples{

data(frogs)
pcru.obscovs <- data.frame(MinAfterSunset=as.vector(t(pcru.data[,,1])),
     Wind=as.vector(t(pcru.data[,,2])),
     Sky=as.vector(t(pcru.data[,,3])),
     Temperature=as.vector(t(pcru.data[,,4])))
pcruUMF <- unmarkedFrameOccu(y = pcru.bin, obsCovs = pcru.obscovs)
pcruUMF.i1 <- imputeMissing(pcruUMF)
pcruUMF.i2 <- imputeMissing(pcruUMF, whichCovs = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
