#' Estimate wildlife abundance or occupancy using hierarchical models.
#'
#' \tabular{ll}{
#' Package: \tab unmarked\cr
#' Type: \tab Package\cr
#' Version: \tab 0.5-6\cr
#' Date: \tab 2009-2-18\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
#'
#' Unmarked estimates wildlife parameters for many popular sampling methods including occupancy and point count data.
#'
#' \strong{Overview of Model-fitting Functions:}  Unmarked provides
#' several functions for fitting integrated likelihood models for wildlife
#' abundance and occurrence to replicated survey data.  \command{\link{occu}}
#' fits occurrence models with no linkage between abundance and detection
#' (MacKenzie et al. 2006).  \command{\link{occuRN}} fits abundance models to
#' presence/absence data by exploiting the link between detection
#' probability and abundance (Royle and Nichols 2003).
#' \command{\link{pcount}} fits N-mixture models for repeated point count data
#' (Royle 2004, Kery and Royle 2005).  \command{\link{mnMix}} fits the
#' multinomial mixture model to categorical observations (Royle and Link
#' 2005).  All of these functions allow the user to specify covariates that
#' affect the detection process and several also allow covariates for the
#' state process.
#'
#' \strong{Data:} All data is passed to unMarked's estimation functions
#'  as a formal S4 class called an unMarkedFrame.
#' See \code{\link{unMarkedFrame}} for a detailed description of unMarkedFrames
#' and how to create them using the convenient constructor \code{unMarkedFrame()}.
#'
#' \strong{Model Specification:}  Most of \emph{unmarked}'s
#' model-fitting functions allow specification of covariates for both the
#' state process and the detection process.  Covariates for the state
#' process (at the site level) and the detection process (at the site or
#' observation level) are specified with the right-hand sided formulas
#' passed to \code{siteformula} and \code{detformula} respectively.  Such a
#' formula looks like \eqn{~ x_1 + x_2 + \ldots + x_n}{~ x1 + x2 + ... +
#' x3} where \eqn{x_1}{x1} through \eqn{x_n}{xn} are additive covariates of
#' the process of interest.  The meaning of these covariates or, what they
#' model, is full described in the help files for the individual functions
#' and is not the same for all functions.
#'
#' \strong{Utility Functions:}  \emph{unMarked} contains several utility
#' functions for massaging data into the form required by its model-fitting
#' functions.  \code{\link{csvToUMF}} converts an appropriately
#' formated comma-separated values (.csv) file to a list containing the
#' components required by model-fitting functions.
#'
#' @references
#' MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.
#' Amsterdam: Academic Press.
#'
#' Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
#' Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
#' pp. 777--790.
#'
#' Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from
#' Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.
#'
#' Kery, M. and Royle, J. A. (2005) Modeling Avaian Abundance from
#' Replicated Counts Using Binomial Mixture Models. \emph{Ecological
#'   Applications} 15(4), pp. 1450--1461.
#'
#' Royle, J. A. and Link W. A. (2005) A general class of multinomial
#' mixture models for anuran calling survey data. \emph{Ecology}, \bold{86}(9),
#' pp. 2505--2512.
#' @name unmarked-package
#' @aliases unmarked
#' @docType package
#' @title Models for Data from Unmarked Animals
#' @author Ian Fiske \email{ijfiske@@ncsu.edu}
#' @keywords package
roxygen()