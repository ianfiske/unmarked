\name{crossVal}
\alias{crossVal}
\alias{crossVal-methods}
\alias{crossVal,unmarkedFit-method}
\alias{crossVal,unmarkedFitList-method}
\alias{unmarkedCrossVal-class}
\alias{unmarkedCrossValList-class}
\alias{show,unmarkedCrossVal-method}
\alias{show,unmarkedCrossValList-method}

\title{Cross-validation methods for fitted unmarked models and fit lists}

\usage{crossVal(object, method=c("kfold","holdout","leaveoneout"),
  folds=10, holdoutPct=0.25)}

\description{Test predictive accuracy of fitted models using several cross-validation approaches}
\arguments{
  \item{object}{A fitted model inheriting class \code{unmarkedFit} or a list of fitted models with class \code{unmarkedFitList}}
  \item{method}{Cross validation method to use as string. Valid options are \code{"kfold"}, \code{"holdout"}, or \code{"leaveoneout"}}
  \item{folds}{Number of folds to use for k-fold cross validation}
  \item{holdoutPct}{Proportion of dataset (value between 0-1) to use as the "holdout" or "test" set, for the holdout method}
} 

\value{\code{unmarkedCrossVal} or \code{unmarkedCrossValList} object containing estimates of root mean square error and mean absolute error}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{fitList}},
  \code{\link{unmarkedFit}}
}

\examples{

#Get data
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))    
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

#Fit occupancy model
fm <- occu(~ obsvar1 ~ 1, pferUMF)

#k-fold cross validation with 10 folds
(kfold = crossVal(fm, method="kfold", folds=10))

#holdout method with 25% of data in test set and 75% in train set
(holdout = crossVal(fm,method='holdout', holdoutPct=0.25))

#Leave-one-out method
(leave = crossVal(fm, method='leaveoneout'))

#Fit a second model and combine into a fitList
fm2 <- occu(~1 ~1, pferUMF)
fl <- fitList(fm2,fm)

#Cross-validation for all fits in fitList using holdout method
(cvlist <- crossVal(fl, method='holdout'))

}
