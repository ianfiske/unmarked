\name{vif}
\alias{vif}
\title{Compute Variance Inflation Factors for an unmarkedFit Object.}
\description{
Compute the variance inflation factors (VIFs) for covariates in one level of the 
model (i.e., occupancy or detection). Calculation of VIFs follows the approach 
of function \code{vif} in package \code{car}, using the correlation matrix of 
fitted model parameters.
}
\usage{
vif(mod, type)
}
\arguments{
  \item{mod}{An unmarked fit object.}
  \item{type}{Level of the model for which to calculate VIFs (for example,
              \code{'state'})}
}
\value{
A named vector of variance inflation factor values for each covariate.}
