#' @include classes.R
roxygen()

# Class to store actual parameter estimates
#' @export
setClass("unMarkedEstimate",
    representation(name = "character",
        estimates = "numeric",
        covMat = "matrix",
        invlink = "character",
        invlinkGrad = "character"),
    validity = function(object) {
      errors <- character(0)
      if(nrow(object@covMat) != length(object@estimates)) {
        errors <- c(errors, "Size of covMat does not match length of estimates.")
      }
      if(length(errors) > 0)
        errors
      else
        TRUE
    })

setClass("unMarkedEstimateLinearComb",
    representation(originalEstimate = "unMarkedEstimate",
        coefficients = "numeric"),
    contains = "unMarkedEstimate")

setClass("unMarkedEstimateBackTransformed",
    representation(transformation = "character"),
    contains = "unMarkedEstimateLinearComb")

setClass("unMarkedEstimateList",
    representation(estimates = "list"),
    validity = function(object) {
      errors <- character(0)
      for(est in object@estimates) {
        if(!is(est, "unMarkedEstimate")) {
          errors <- c("At least one element of unMarkedEstimateList is not an unMarkedEstimate.")
          break
        }
      }
      if(length(errors) == 0) {
        return(TRUE)
      } else {
        return(errors)
      }
    })

setMethod("show", "unMarkedEstimateList",
    function(object) {
      for(est in object@estimates) {
        show(est)
        cat("\n")
      }
    })

unMarkedEstimateList <- function(l) {
  new("unMarkedEstimateList", estimates = l)
}

#' @export
unMarkedEstimate <- function(name, estimates, covMat, invlink, invlinkGrad) {

  new("unMarkedEstimate",
      name = name,
      estimates = estimates,
      covMat = covMat,
      invlink = invlink,
      invlinkGrad = invlinkGrad)

}

setMethod("show",
    signature(object = "unMarkedEstimate"),
    function(object) {
      ests <- object@estimates
      SEs <- SE(object)
      Z <- ests/SEs
      p <- 2*pnorm(abs(Z), lower.tail = FALSE)

      if(is(object, "unMarkedEstimateLinearComb")) {
        printRowNames <- FALSE
      } else {
        printRowNames <- TRUE
      }

      cat(object@name,":\n", sep="")
      outDF <- data.frame(Estimate = ests,
          SE = SEs,
          z = Z,
          "P(>|z|)" = p,
          check.names = FALSE)
      print(outDF, row.names = printRowNames, digits = 3)
    })

setMethod("show",
    signature(object = "unMarkedEstimateLinearComb"),
    function(object) {
      coefTable <- data.frame(Estimate = object@originalEstimate@estimates,
          Coefficients = object@coefficients)

      callNextMethod(object)

      cat("\n")

      print(coefTable, digits = 3)

    })

setMethod("show",
    signature(object = "unMarkedEstimateBackTransformed"),
    function(object) {
      callNextMethod(object)
      cat("\nTransformation:", object@transformation,"\n")
    })


#' Compute linear combinations of estimates in unMarkedEstimate objects.
#'
#' This function computes the linear combination of parameter estimates in
#' \code{obj} given by the coefficient vector.  This may be useful to estimate
#' quantities of interest from a fitted model or to test a hypothesis.
#'
#' @name linearComb-unMarkedEstimate
#' @aliases linearComb,unMarkedEstimate-method
#' @param obj an unMarkedEstimate object
#' @param coefficients vector of same length as obj
#' @return an unMarkedEstimate object
setMethod("linearComb",
    signature(obj = "unMarkedEstimate", coefficients = "numeric"),
    function(obj, coefficients) {
      stopifnot(length(coefficients) == length(obj@estimates))
      e <- as.numeric(t(coefficients) %*% obj@estimates)
      v <- t(coefficients) %*% obj@covMat %*% coefficients
      umelc <- new("unMarkedEstimateLinearComb",
          name = paste("Linear combination of",obj@name,"estimate(s)"),
          estimates = e, covMat = v,
          invlink = obj@invlink, invlinkGrad = obj@invlinkGrad,
          originalEstimate = obj, coefficients = coefficients)
      umelc
    })



#' Transform an unMarkedEstimate object to it's natural scale.
#'
#' The transformation is determined by the invlink and invlinkGrad slots
#' in \code{obj}.  These slots specify the name of a one-to-one function and
#' its gradient respectively.  The delta method is used to compute the transformed
#' estimate.
#'
#' @param unMarkedEstimate object to be transformed
#' @return an unMarkedEstimate object representing the transformed estimate.
#' This object has invlink and invlinkGrad slots as the identity function.
setMethod("backTransform",
    signature(obj = "unMarkedEstimate"),
    function(obj) {
      stopifnot(length(obj@estimates) == 1)
      e <- eval(call(obj@invlink,obj@estimates))
      v <- (eval(call(obj@invlinkGrad,obj@estimates)))^2 * obj@covMat

      if(is(obj, "unMarkedEstimateLinearComb")) {
        coef <- obj@coefficients
        orig <- obj@originalEstimate
      } else {
        coef <- 1
        orig <- obj
      }

      umebt <- new("unMarkedEstimateBackTransformed",
          name = paste(obj@name,"transformed to native scale"),
          estimates = e, covMat = v,
          invlink = "identity", invlinkGrad = "identity",
          originalEstimate = orig, coefficients = coef,
          transformation = obj@invlink)
      umebt
    })

#' Compute standard error of an unMarkedEstimate object.
#'
#' This function computes the large-sample standard error from the inverse of the
#' hessian matrix.
#'
#' @name SE-unMarkedEstimate
#' @aliases SE,unMarkedEstimate-method
#' @param obj unMarkedEstimate whose standard error is returned
#' @return vector of the standard error(s) of estimates in obj
setMethod("SE",
    signature(obj = "unMarkedEstimate"),
    function(obj) {
      sqrt(diag(obj@covMat))
    })


setMethod("[",
    signature("unMarkedEstimateList"),
    function(x, i, j, drop) {
      x@estimates[[i]]
    })

setMethod("names", "unMarkedEstimateList",
    function(x) {
      names(x@estimates)
    })