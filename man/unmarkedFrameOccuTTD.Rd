\name{unmarkedFrameOccuTTD}

\alias{unmarkedFrameOccuTTD}

\title{Create an unmarkedFrameOccuTTD object for the time-to-detection model
       fit by occuTTD}

\usage{
  unmarkedFrameOccuTTD(y, surveyLength, siteCovs=NULL, obsCovs=NULL, 
                           numPrimary=1, yearlySiteCovs=NULL)
}

\description{Organizes time-to-detection occupancy data along with covariates. 
  This S4 class is required by the data argument of \code{\link{occuTTD}}}

\arguments{
    \item{y}{An MxR matrix of time-to-detection data for a species, 
        where M is the number of sites and R is the maximum number of 
        observations per site (across all primary periods and observations, if 
        you have multi-season data). Values in \code{y} should be positive.}
    \item{surveyLength}{The maximum length of a survey, in the same units as 
        \code{y}. You can provide either a single value (if all surveys had
        the same max length), or a matrix matching the dimensions of \code{y}
        (if surveys had different max lengths).}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of 
        covariates that vary within sites, or a \code{\link{data.frame}} with 
        MxR rows in the ordered by site-observation (if single-season) or 
        site-primary period-observation (if multi-season).}
    \item{numPrimary}{Number of primary time periods (e.g. seasons) for the 
        dynamic or multi-season version of the model. There should be
        an equal number of secondary periods in each primary period.}
    \item{yearlySiteCovs}{A data frame with one column per covariate that varies 
        among sites and primary periods (e.g. years). It should have MxT rows
        where M is the number of sites and T the number of primary periods,
        ordered by site-primary period. These covariates only used for dynamic 
        (multi-season) models.}
  
}

\details{
    unmarkedFrameOccuTTD is the S4 class that holds data to be passed 
    to the \code{\link{occuTTD}} model-fitting function.
}

\note{
  If the time-to-detection values in \code{y} are very large (e.g., because
  they are expressed as numbers of seconds) you may have issues fitting models.
  An easy solution is to convert your units (e.g., from seconds to decimal
  minutes) to keep the values as close to 0 as possible.
}

\value{an object of class unmarkedFrameOccuTTD}

\author{Ken Kellner \email{contact@kenkellner.com}}

\examples{
  
  # For a single-season model
  N <- 100 #Number of sites
  psi <- 0.4 #Occupancy probability
  lam <- 7 #Parameter for exponential distribution of time to detection
  Tmax <- 10 #Maximum survey length

  z <- rbinom(N, 1, psi) #Simulate occupancy
  y <- rexp(N, 1/lam) #Simulate time to detection
  y[z==0] <- Tmax
  y[y>Tmax] <- Tmax
  
  sc <- as.data.frame(matrix(rnorm(N*2),ncol=2)) #Site covs
  oc <- as.data.frame(matrix(rnorm(N*2),ncol=2)) #obs covs

  umf <- unmarkedFrameOccuTTD(y=y, surveyLength=Tmax, siteCovs=sc, obsCovs=oc)
  
}
