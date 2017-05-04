\name{jay}
\alias{jay}
\docType{data}
\title{
European Jay data from the Swiss Breeding Bird Survey 2002
}
\description{
The Swiss breeding bird survey ("Monitoring Haufige Brutvogel" MHB) has monitored the populations of 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol, where each square is surveyed up to three times during the breeding season (only twice above the tree line). Surveys are conducted along a transect that does not change over the years.

The list \code{jay} has the data for European Jay territories for 238 sites surveyed in 2002.
}

\usage{data("jay")}

\format{
\code{jay} is a list with 3 elements:

\describe{
  \item{\bold{caphist }}{a data frame with rows for 238 sites and columns for each of the observable detection histories. For the sites visited 3 times, these are \code{"100", "010", "001", "110", "101", "011", "111"}. Sites visited twice have \code{"10x", "01x", "11x"}. 

Each row gives the number of territories with the corresponding detection history, with NA for the detection histories not applicable: sites visited 3 times have NAs in the last 3 columns while those visited twice have NAs in the first 7 columns.
}
\describe{
  \item{\bold{sitescovs }}{a data frame with rows for 238 sites, and the following columns:

    \enumerate{
      \item elev : the mean elevation of the quadrat, m.
      \item length : the length of the route walked in the quadrat, km.
      \item forest : percentage forest cover.
    }}}
\describe{
  \item{\bold{covinfo }}{a data frame with rows for 238 sites, and the following columns:

    \enumerate{
      \item x, y : the coordinates of the site.
      \item date1, date2, date3 : the Julian date of the visit, with 1 April = 1. Sites visited twice have NA in the 3rd column.
      \item dur1, dur2, dur3 : the duration of the survey, mins. For 10 visits the duration is not available, so there are additional NAs in these columns.
    }}}
} }

\note{
In previous versions, \code{jay} had additional information not required for the analysis, and a data frame with essentially the same information as the \code{Switzerland} data set.
}
\source{
Swiss Ornithological Institute
}
\references{
Royle, J.A., Kery, M., Gauthier, R., Schmid, H. (2007) Hierarchical spatial models of abundance and occurrence from imperfect survey data. \emph{Ecological Monographs}, 77, 465-481.

Kery & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 7.9
}
\examples{
data(jay)
str(jay)

# Carry out a simple analysis, without covariates:
# Create a customised piFun (see ?piFun for details)
crPiFun <- function(p) {
   p1 <- p[,1] # Extract the columns of the p matrix, one for 
   p2 <- p[,2] #   each of J = 3 sample occasions
   p3 <- p[,3]
   cbind(      # define multinomial cell probabilities:
      "100" = p1 * (1-p2) * (1-p3),
      "010" = (1-p1) * p2 * (1-p3),
      "001" = (1-p1) * (1-p2) * p3,
      "110" = p1 * p2 * (1-p3),
      "101" = p1 * (1-p2) * p3,
      "011" = (1-p1) * p2 * p3,
      "111" = p1 * p2 * p3,
      "10x" = p1*(1-p2),
      "01x" = (1-p1)*p2,
      "11x" = p1*p2)
}
# Build the unmarkedFrame object
mhb.umf <- unmarkedFrameMPois(y=as.matrix(jay$caphist),
  obsToY=matrix(1, 3, 10), piFun="crPiFun")
# Fit a model
( fm1 <- multinomPois(~1 ~1, mhb.umf) )
}
\keyword{datasets}
