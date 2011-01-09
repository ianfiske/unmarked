\name{unmarkedFrame}
\title{Create an object of class unmarkedFramePCO that contains data used by pcountOpen.}
\alias{unmarkedFramePCO}
\usage{unmarkedFramePCO(y, siteCovs=NULL, obsCovs=NULL, dates)}
\description{Organizes repeated count data along with the covariates and possibly the dates on which each survey was conducted. This S4 class is required by the data argument of \code{\link{pcountOpen}}}
\arguments{
	\item{y}{An MxT matrix of the repeated count data, where M is the number 
        of sites and T is the maximum number of survey occasions per site.}
	\item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have M rows and one column per covariate}
	\item{obsCovs}{Either a named list of \code{\link{data.frame}}s of covariates 
        that vary within sites, or a \code{\link{data.frame}} with MxT rows in 
        site-major order.}
  \item{dates}{integers indicating the dates of the surveys.}
}
\details{unmarkedFramePCO is the S4 class that holds data to be passed 
to the \code{\link{pcountOpen}} model-fitting function.

The unmarkedFramePCO class is similar to the unmarkedFramePCount class except
that it contains the dates for each survey, which needs to be supplied .

}

\value{an object of class unmarkedFramePCO}
\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}}, 
    \code{\link{pcountOpen}}}
\examples{

# Repeated count data
y <- matrix(c(
    0, 2, 3, 2, 0,
    2, 2, 3, 1, 1,
    1, 1, 0, 0, 3,
    0, 0, 0, 0, 0), nrow=4, ncol=5, byrow=TRUE)

# Site-specific covariates
sc <- data.frame(x1 = 1:4, x2 = c('A','A','B','B'))

# Observation-specific covariates
oc <- list(
    x3 = matrix(1:5, nrow=4, ncol=5, byrow=TRUE),
    x4 = matrix(letters[1:5], nrow=4, ncol=5, byrow=TRUE))

# Dates of surveys    
dates <- matrix(c(
    1, 2, 5, 7, 8,
    1, 2, 3, 4, 5,
    1, 2, 4, 5, 6,
    1, 3, 5, 6, 7), nrow=4, ncol=5, byrow=TRUE)
mode(dates) <- "integer"
     

# Create the unmarkedFrame    
umf <- unmarkedFramePCO(y=y, siteCovs=sc, obsCovs=oc, dates=dates)    

# Take a look
umf

summary(umf)


# Subset
umf[,1] # first survey occasion
umf[1,] # first row
umf[1:2, 2:3]


}
