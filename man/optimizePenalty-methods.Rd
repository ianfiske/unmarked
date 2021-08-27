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

  Currently the only fitting function that supports \code{optimizePenalty} is
  \code{occuMulti} for multispecies occupancy modeling; see Clipp et al. (2021).
}

\usage{
\S4method{optimizePenalty}{unmarkedFitOccuMulti}(
  object, penalties = c(0, 2^seq(-4, 4)), k = 5, boot = 30, ...)
}

\arguments{
  \item{object}{A fitted model inheriting class \code{unmarkedFit}}
  \item{penalties}{Vector of possible penalty values, all of which must be >= 0}
  \item{k}{Number of folds to use for k-fold cross validation}
  \item{boot}{Number of bootstrap samples to use to generate the 
    variance-covariance matrix for the final model.}
  \item{...}{Other arguments, currently ignored}
} 

\value{\code{unmarkedFit} object of same type as input, with the optimal
  penalty value applied.
}

\references{

Clipp, H. L., Evans, A., Kessinger, B. E., Kellner, K. F., and C. T. Rota.
  2021. A penalized likelihood for multi-species occupancy models improves 
  predictions of species interactions. Ecology.

Hutchinson, R. A., J. V. Valente, S. C. Emerson, M. G. Betts, and
  T. G. Dietterich. 2015. Penalized Likelihood Methods Improve Parameter
  Estimates in Occupancy Models. Methods in Ecology and Evolution. DOI:
  10.1111/2041-210X.12368

}


\author{Ken Kellner \email{contact@kenkellner.com}}
