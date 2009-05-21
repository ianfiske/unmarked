#' @include unmarkedEstimate.R
#' @include unmarkedFrame.R
#' @include classes.R
roxygen()

# Class to store unMarked model fit information
#
# slot fitType Name of the model that was fit.
# @slot stateformula The abundance/occupancy formula.
# @slot detformula The formula governing the detection process.
# @slot data The unMarkedFrame containing the data that was fit.
# @slot stateMLE The MLE of the abundance/occupancy parameters
# @slot stateSE The standard errors of the state MLE's.
# @slot detMLE MLE of the detection parameters.
# @slot detSE Standard errors of the detection MLE's.
# @slot AIC The AIC of this model fit.
# @slot negLogLike The negative log likelihood of the fitted model.
# A Class to store fit results from unMarkedFrames.
#' @export
setClass("unmarkedFit",
    representation(fitType = "character",
        call = "call",
        data = "unMarkedFrame",
        estimates = "unmarkedEstimateList",
        AIC = "numeric",
        hessian = "matrix",
        negLogLike = "numeric"))



# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call,
    data, estimates, AIC, hessian) {
  umfit <- new("unmarkedFit", fitType = fitType,
      call = call, data = data,
      estimates = estimates, AIC = AIC,
      hessian = hessian)

  return(umfit)
}

# @export
setMethod("show", "unmarkedFit",
    function(object) {
      cat("\nCall:\n")
#      cat(object@fitType,"(stateformula = ~ ",
#          as.character(object@stateformula)[2],
#          ", detformula = ~ ",
#          as.character(object@detformula)[2],")\n\n", sep = "")
      print(object@call)

      cat("\n")

      show(object@estimates)

      #      show(object@stateEstimates)
#
#      cat("\n")
#
#      show(object@detEstimates)

#            stateEsts <- object@stateEstimates@estimates
#            stateSE <- SE()
#            stateZ <- stateEsts/stateSE
#            stateP <- 2*pnorm(abs(stateZ), lower.tail = FALSE)
#
#
#            detEsts <- object@detEstimates@estimates
#            detSE <- standardError(object@detEstimates)
#            detZ <- detEsts/detSE
#            detP <- 2*pnorm(abs(detZ), lower.tail = FALSE)

#            cat("\nState coefficients:\n")
#            stateDF <- data.frame(Estimate = stateEsts,
#                                  SE = stateSE,
#                                  z = stateZ,
#                                  "p-val" = stateP)
#            show(round(stateDF,3))
#
#            cat("\nDetection coefficients:\n")
#            detDF <- data.frame(Estimate = detEsts,
#                                  SE = detSE,
#                                  z = detZ,
#                                  p.val = detP)
#            show(round(detDF,3))

      cat("\nAIC:", object@AIC,"\n")
    })


#' Compute linear combinations of estimates in unmarkedFit objects.
#'
#' This function computes the linear combination of parameter estimates in
#' \code{obj} given by the coefficient vector.  The user must
#' specify whether the detection or state covariates should be considered via the
#' \code{whichEstimate} argument.  This may be useful to estimate
#' quantities of interest from a fitted model or to test a hypothesis.
#'
#' @name linearComb-unmarkedFit
#' @aliases linearComb,unmarkedFit-method
#' @param obj an unmarkedFit object
#' @param coefficients vector of same length as obj
#' @param whichEstimate character, either "state" or "det"
#' @return an unmarkedEstimate object
setMethod("linearComb",
    signature(obj = "unmarkedFit", coefficients = "numeric"),
    function(obj, coefficients, whichEstimate) {
      stopifnot(!missing(whichEstimate))
      stopifnot(whichEstimate %in% names(obj))
      estimate <- obj@estimates[whichEstimate]
#      if(whichEstimate == "det") {
#        lc <- linearComb(obj@detEstimates, coefficients)
#      } else {
#        lc <- linearComb(obj@stateEstimates, coefficients)
#      }
#      lc
      linearComb(estimate, coefficients)
    })

setMethod("[",
    "unmarkedFit",
    function(x, i, j, drop) {
      x@estimates[i]
    })

setMethod("names", "unmarkedFit",
    function(x) {
      names(x@estimates)
    })