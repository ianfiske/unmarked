\name{optimizePenalty-methods}
\docType{methods}
\alias{optimizePenalty}
\alias{optimizePenalty-methods}
\alias{optimizePenalty,unmarkedFitOccuMulti-method}

\title{ Identify Optimal Penalty Parameter Value }
\description{
  Identify the optimal value of the penalty term for \code{unmarked} models that
  support penalized likelihood. For each potential value of the penalty term, 
  K-fold cross validation is performed. Log-likelihoods for the test data in 
  each fold are calculated and summed. The penalty term that maximizes the sum 
  of the fold log-likelihoods is selected as the optimal value. Finally, the 
  model is re-fit with the full dataset using the selected penalty term.
  Right now only Bayes-inspired penalty of Hutchinson et al. (2015) is supported.
}

\usage{
\S4method{optimizePenalty}{unmarkedFitOccuMulti}(
  object, penalties = c(0, 2^seq(-4, 4)), k = 5, ...)
}

\arguments{
  \item{object}{A fitted model inheriting class \code{unmarkedFit}}
  \item{penalties}{Vector of possible penalty values, all of which must be >= 0}
  \item{k}{Number of folds to use for k-fold cross validation}
  \item{...}{Other arguments, currently ignored}
} 

\value{\code{unmarkedFit} object of same type as input, with the optimal
  penalty value applied.
}

\references{

Hutchinson, R. A., J. V. Valente, S. C. Emerson, M. G. Betts, and
  T. G. Dietterich. 2015. Penalized Likelihood Methods Improve Parameter
  Estimates in Occupancy Models. Methods in Ecology and Evolution. DOI:
  10.1111/2041-210X.12368

}


\author{Ken Kellner \email{contact@kenkellner.com}}
