\name{gpcount}
\alias{gpcount}
\title{
Generalized N-mixture model for repeated count data
}
\description{
Fit the model of Chandler et al. (2011).
}
\usage{
gpcount(lambdaformula, phiformula, pformula, data, mixture = c("P", "NB"), K, starts, method = "BFGS", se = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambdaformula}{
%%     ~~Describe \code{lambdaformula} here~~
}
  \item{phiformula}{
%%     ~~Describe \code{phiformula} here~~
}
  \item{pformula}{
%%     ~~Describe \code{pformula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{mixture}{
%%     ~~Describe \code{mixture} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Richard Chandler
}
\note{
The interpretation of the parameter phi
}


\seealso{
\link{\code{gmultmix}}, \link{\code{gdistsamp}}
}
\examples{
set.seed(54)

nSites <- 50
nVisits <- 4
nReps <- 3

lambda <- 5
phi <- 0.6
p <- 0.4

M <- rpois(nSites, lambda) # super-population size

N <- matrix(NA, nSites, nVisits)
y <- array(NA, c(nSites, nReps, nVisits))
for(i in 1:nVisits) {
    N[,i] <- rbinom(nSites, M, phi) # population available during vist j
}
colMeans(N)

for(i in 1:nSites) {
    for(j in 1:nVisits) {
        y[i,,j] <- rbinom(nReps, N[i,j], p)
    }
}

ym <- matrix(y, nSites)
ym[1,] <- NA
ym[2, 1:nReps] <- NA
ym[3, (nReps+1):(nReps+nReps)] <- NA
umf <- unmarkedFrameGPC(y=ym, numPrimary=nVisits)

fmu <- gpcount(~1, ~1, ~1, umf, K=40, control=list(trace=TRUE, REPORT=1))

backTransform(fmu, type="lambda")
backTransform(fmu, type="phi")
backTransform(fmu, type="det")

}
