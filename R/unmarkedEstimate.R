#' @include classes.R
{}

setClassUnion("matrixOrVector", c("matrix","numeric"))

# Class to store actual parameter estimates
#' @export
setClass("unmarkedEstimate",
    representation(name = "character",
				short.name = "character",
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

setMethod("estimates", "unmarkedEstimateList",
		function(object) {
			object@estimates
		})

unmarkedEstimateList <- function(l) {
  new("unmarkedEstimateList", estimates = l)
}

#' @export
unmarkedEstimate <- function(name, short.name, estimates, covMat, invlink, invlinkGrad) {

  new("unmarkedEstimate",
      name = name,
			short.name = short.name,
      estimates = estimates,
      covMat = covMat,
      invlink = invlink,
      invlinkGrad = invlinkGrad)

}

setMethod("show",
    signature(object = "unmarkedEstimate"),
    function(object) {
      ests <- object@estimates
      SEs <- SE(object)
      Z <- ests/SEs
			p <- 2*pnorm(abs(Z), lower.tail = FALSE)
			
			printRowNames <- 
					!(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(as.numeric(names(ests)), 1))
			
			cat(object@name,":\n", sep="")
			outDF <- data.frame(Estimate = ests,
					SE = SEs,
					z = Z,
          "P(>|z|)" = p,
          check.names = FALSE)
      print(outDF, row.names = printRowNames, digits = 3)
			
    })




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
		signature(obj = "unmarkedEstimate", coefficients = "matrixOrVector"),
		function(obj, coefficients) {
			if(!is(coefficients, "matrix")) coefficients <- t(as.matrix(coefficients))
			stopifnot(ncol(coefficients) == length(obj@estimates))
			e <- as.vector(coefficients %*% obj@estimates)
			v <- coefficients %*% obj@covMat %*% t(coefficients)
			umelc <- new("unmarkedLinComb",
					parentEstimate = obj,
					estimate = e, covMat = v,
					coefficients = coefficients)
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



# backTransform is only valid for an unmarkedEstimate of length = 1.
# can backtranform a fit directly if it has length 1
# o.w. give error
setMethod("backTransform", "unmarkedEstimate",
		function(obj) {
			if(length(obj@estimates) == 1) {
				lc <- linearComb(obj, 1)
				return(backTransform(lc))
			} else {
				stop("Cannot directly back-transform an unmarkedEstimate with length > 1.")
			}
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
			coefs
		})

setMethod("vcov", "unmarkedEstimate",
		function(object,...) {
			v <- object@covMat
			rownames(v) <- colnames(v) <- names(coef(object))
			v
		})
 
setMethod("confint", "unmarkedEstimate", 
		function(object, parm, level = 0.95) {
			if(missing(parm)) parm <- 1:length(object@estimates)
			ests <- object@estimates[parm]
			ses <- SE(object)[parm]
			z <- qnorm((1-level)/2, lower.tail = FALSE)
			lower.lim <- ests - z*ses
			upper.lim <- ests + z*ses
			ci <- as.matrix(cbind(lower.lim, upper.lim))
			rownames(ci) <- names(ests)
			colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
			ci
		})





