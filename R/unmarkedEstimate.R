#' @include classes.R
roxygen()

# Class to store actual parameter estimates
#' @export
setClass("unmarkedEstimate",
    representation(name = "character",
				short.name = "character",
        estimates = "numeric",
        covMat = "matrix",
        invlink = "character",
        invlinkGrad = "character",
				backTransformed = "logical",
				coefficients = "matrix"),
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

#setClass("unmarkedEstimateLinearComb",
#    representation(originalEstimate = "unmarkedEstimate",
#        coefficients = "numeric"),
#    contains = "unmarkedEstimate")
#
## coefficients are now a matrix where each row is a 
#setClass("unmarkedEstimateMultLinearCombs",
#		representation(originalEstimate = "unmarkedEstimate",
#				coefficients = "matrix"),
#		contains = "unmarkedEstimate")

#setClass("unmarkedEstimateBackTransformed",
#    representation(transformation = "character"),
#    contains = "unmarkedEstimateLinearComb")

setClass("unmarkedEstimateList",
    representation(estimates = "list"),
    validity = function(object) {
      errors <- character(0)
      for(est in object@estimates) {
        if(!is(est, "unmarkedEstimate")) {
          errors <- c("At least one element of unmarkedEstimateList is not an unmarkedEstimate.")
          break
        }
      }
      if(length(errors) == 0) {
        return(TRUE)
      } else {
        return(errors)
      }
    })

setMethod("show", "unmarkedEstimateList",
    function(object) {
      for(est in object@estimates) {
        show(est)
        cat("\n")
      }
    })

#' @export
setGeneric("estimates",
    function(object) {
      standardGeneric("estimates")
    })

setMethod("estimates", "unmarkedEstimate",
    function(object) {
      object@estimates
    })

unmarkedEstimateList <- function(l) {
  new("unmarkedEstimateList", estimates = l)
}

#' @export
unmarkedEstimate <- function(name, short.name, estimates, covMat, invlink, invlinkGrad, backTransformed = FALSE,
		coefficients = matrix(1, 1, length(estimates))) {

  new("unmarkedEstimate",
      name = name,
			short.name = short.name,
      estimates = estimates,
      covMat = covMat,
      invlink = invlink,
      invlinkGrad = invlinkGrad,
			backTransformed = backTransformed,
			coefficients = coefficients)

}

setMethod("show",
    signature(object = "unmarkedEstimate"),
    function(object) {
      ests <- object@estimates
      SEs <- SE(object)
      Z <- ests/SEs
      p <- 2*pnorm(abs(Z), lower.tail = FALSE)

      if(!all(object@coefficients == 1)) {
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

#setMethod("show",
#    signature(object = "unmarkedEstimateLinearComb"),
#    function(object) {
#      coefTable <- data.frame(Estimate = object@originalEstimate@estimates,
#          Coefficients = object@coefficients)
#
#      callNextMethod(object)
#
#      cat("\n")
#
#      print(coefTable, digits = 3)
#
#    })
#
#setMethod("show",
#    signature(object = "unmarkedEstimateBackTransformed"),
#    function(object) {
#      callNextMethod(object)
#      cat("\nTransformation:", object@transformation,"\n")
#    })


#' Compute linear combinations of estimates in unmarkedEstimate objects.
#'
#' This function computes the linear combination of parameter estimates in
#' \code{obj} given by the coefficient vector.  This may be useful to estimate
#' quantities of interest from a fitted model or to test a hypothesis.
#'
#' @name linearComb-unmarkedEstimate
#' @aliases linearComb,unmarkedEstimate-method
#' @param obj an unmarkedEstimate object
#' @param coefficients vector of same length as obj
#' @return an unmarkedEstimate object
setMethod("linearComb",
    signature(obj = "unmarkedEstimate", coefficients = "numeric"),
    function(obj, coefficients) {
      stopifnot(length(coefficients) == length(obj@estimates))
      e <- as.vector(t(coefficients) %*% obj@estimates)
      v <- t(coefficients) %*% obj@covMat %*% coefficients
      umelc <- new("unmarkedEstimate",
          name = paste("Linear combination of",obj@name,"estimate(s)"),
          estimates = e, covMat = v,
          invlink = obj@invlink, invlinkGrad = obj@invlinkGrad,
          coefficients = t(coefficients))
      umelc
    })


# do the right thing if a matrix of coefficients supplied
setMethod("linearComb",
		signature(obj = "unmarkedEstimate", coefficients = "matrix"),
		function(obj, coefficients) {
			stopifnot(ncol(coefficients) == length(obj@estimates))
			e <- as.vector(coefficients %*% obj@estimates)
			v <- coefficients %*% obj@covMat %*% t(coefficients)
			umelc <- new("unmarkedEstimate",
					name = paste("Linear combinations of",obj@name,"estimate(s)"),
					short.name = obj@short.name,
					estimates = e, covMat = v,
					invlink = obj@invlink, invlinkGrad = obj@invlinkGrad,
					coefficients = coefficients, backTransformed = FALSE)
			umelc
		})


#' Transform an unmarkedEstimate object to it's natural scale.
#'
#' The transformation is determined by the invlink and invlinkGrad slots
#' in \code{obj}.  These slots specify the name of a one-to-one function and
#' its gradient respectively.  The delta method is used to compute the transformed
#' estimate.
#'
#' @param unmarkedEstimate object to be transformed
#' @return an unmarkedEstimate object representing the transformed estimate.
#' This object has invlink and invlinkGrad slots as the identity function.
#setMethod("backTransform",
#    signature(obj = "unmarkedEstimate"),
#    function(obj) {
#      stopifnot(length(obj@estimates) == 1)
#      e <- eval(call(obj@invlink,obj@estimates))
#      v <- (eval(call(obj@invlinkGrad,obj@estimates)))^2 * obj@covMat
#
#      if(is(obj, "unmarkedEstimateLinearComb")) {
#        coef <- obj@coefficients
#        orig <- obj@originalEstimate
#      } else {
#        coef <- 1
#        orig <- obj
#      }
#
#      umebt <- new("unmarkedEstimateBackTransformed",
#          name = paste(obj@name,"transformed to native scale"),
#          estimates = e, covMat = v,
#          invlink = "identity", invlinkGrad = "identity",
#          originalEstimate = orig, coefficients = coef,
#          transformation = obj@invlink)
#      umebt
#    })


# TODO: unify unmarkedEstimateMultLinearCombs and unmarkedEstimateLinearComb
setMethod("backTransform",
		signature(obj = "unmarkedEstimate"),
		function(obj) {
			
			stopifnot(!obj@backTransformed)
			## In general, MV delta method is Var=J*Sigma*J^T where J is Jacobian
			## In this case, J is diagonal with elements = gradient
			## This reduces to scaling the rows and then columns of Sigma by the gradient
			e <- eval(call(obj@invlink,obj@estimates))
			grad <- eval(call(obj@invlinkGrad,obj@estimates))
			v <- diag(grad) %*% obj@covMat %*% diag(grad) 
			
#			if(is(obj, "unmarkedEstimateLinearComb")) {
				coef <- obj@coefficients
#				orig <- obj@originalEstimate
#			} else {
#				coef <- 1
#				orig <- obj
#			}
			
			umebt <- new("unmarkedEstimate",
					name = paste(obj@name,"transformed to native scale"),
					short.name = obj@short.name,
					estimates = e, covMat = v,
					invlink = "identity", invlinkGrad = "identity",
					coefficients = coef,
					backTransformed = TRUE)
			
			umebt
		})

#' Compute standard error of an unmarkedEstimate object.
#'
#' This function computes the large-sample standard error from the inverse of the
#' hessian matrix.
#'
#' @name SE-unmarkedEstimate
#' @aliases SE,unmarkedEstimate-method
#' @param obj unmarkedEstimate whose standard error is returned
#' @return vector of the standard error(s) of estimates in obj
setMethod("SE",
    signature(obj = "unmarkedEstimate"),
    function(obj) {
      sqrt(diag(obj@covMat))
    })


setMethod("[",
    signature("unmarkedEstimateList"),
    function(x, i, j, drop) {
      x@estimates[[i]]
    })

setMethod("names", "unmarkedEstimateList",
    function(x) {
      names(x@estimates)
    })

setMethod("coef", "unmarkedEstimate",
		function(object, altNames = TRUE, ...) {
			coefs <- object@estimates
			names(coefs)[names(coefs) == "(Intercept)"] <- "Int"
			if(altNames) {
				names(coefs) <- paste(object@short.name, "(", names(coefs), ")", sep="")
			}
			if(object@backTransformed | !all(object@coefficients == 1)) names(coefs) <- NULL
			coefs
		})

setMethod("vcov", "unmarkedEstimate",
		function(object,...) {
			object@covMat
		})
#    })
#    
#    




# #' @exportMethod coef
#setMethod("coef", "unmarkedFit", function(object, type=NULL)
#		{
#			if(is.null(type)) {
#				eap <- object@estimates["state"]@estimates
#				edp <- object@estimates["det"]@estimates
#				e <- c(eap, edp)
#			}	else {
#				e <- object@estimates[type]@estimates
#			}
#			return(e)
#		})




# #' @exportMethod vcov
#setMethod("vcov", "unmarkedFit", function(object, type=NULL, drop=F)
#		{
#			if(is.null(type)) {
#				vc <- solve(object@hessian)
#			} 
#			else {
#				if(!is.null(type)) {
#					vc <- object@estimates[type]@covMat
#				}
#			}
#			return(vc)
#		})

