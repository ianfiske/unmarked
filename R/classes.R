#' @import methods

#' @exportClass optionalDataFrame
setClassUnion("optionalDataFrame", c("data.frame","NULL"))

#' @exportClass optionalMatrix
setClassUnion("optionalMatrix", c("matrix","NULL"))

validUnMarkedFrame <- function(object) {
	errors <- character(0)
	M <- nrow(object@y)
	if(!is.null(object@siteCovs))
		if(nrow(object@siteCovs) != M)
			errors <- c(errors, "siteCovData does not have same size number of sites as y.")
	if(!is.null(object@obsCovs))
		if(nrow(object@obsCovs) != M*object@obsNum)
			errors <- c(errors, "obsCovData does not have M*obsNum rows.")
	if(length(errors) == 0)
		TRUE
	else
		errors
}

## not used in roxygen currently.
# @slot y A matrix of the observed measured data.
# @slot obsCovData Dataframe of covariates that vary within sites.
# @slot siteCovData Dataframe of covariates that vary at the site level.
# @slot obsNum Number of observations per site. For most models, this
# can be taken to be the number of columns in y.  But this is not always
# the case.  For example, double observer: y has 3 columns, but only 2
# independent observations were taken at each site.
# @slot primaryNum integer number of primary seasons for multiseason data only

#' Class to hold data for analyses in unmarked.
#'
#' @export
setClass("unMarkedFrame",
         representation(y = "matrix",
                        obsCovs = "optionalDataFrame",
                        siteCovs = "optionalDataFrame",
                        obsNum = "numeric",
                        primaryNum = "numeric"),
         validity = validUnMarkedFrame)

#' Constructor for unMarkedFrames.
#'
#' unMarkedFrame is the S4 class that holds data structures to be passed to the model-fitting functions in unMarked.
#'
#' An unMarkedFrame contains the observations (\code{y}), covariates measured at the observation level (\code{obsCovs}), and covariates measured at the site level (\code{siteCovs}).
#' For a data set with M sites and J observations at each site, y is an M x J matrix.
#' \code{obsCovs} and \code{siteCovs} are both data frames (see \link{data.frame}).  \code{siteCovs} has M rows so that each row contains the covariates for the corresponding sites.
#' \code{obsCovs} has M*obsNum rows so that each covariates is ordered by site first, then observation number.  Missing values are coded with \code{NA} in any of y, siteCovs, or obsCovs.
#'
#' Additionally, unMarkedFrames contain metadata, obsNum and primaryNum.  obsNum is the number of observations measured at each site. primaryNum is the number of seasons in a robust design sampling scheme.
#' Typically, these can be automatically determined by the constructor.  If not specified, obsNum is taken to be the number of columns in y and primaryNum is taken to be 1.
#' However, for certain situations, these must be supplied.  For example, double observer sampling, y has 3 columns corresponding the observer 1, observer 2, and both, but there were only two independent observations.
#' In this situation, y has 3 columns, but obsNum must be specified as 2.  This flexibility is currenty only used in the function \link{multinomPois}.
#'
#' For convenience, \code{obsCovs} can be a list of M x obsNum matrices, with each one corresponding to an observation level covariate.
#'
#' All site-level covariates are automatically copied to obsCovs so that site level covariates are available at the observation level.
#'
#' @title Create an unMarkedFrame.
#' @param y A matrix of the observed measured data.
#' @param obsCovs Dataframe of covariates that vary within sites.
#' @param siteCovs Dataframe of covariates that vary at the site level.
#' @param obsNum Number of independent observations.
#' @param primaryNum Number of primary time periods (seasons in the multiseason model).
#' @return an unMarkedFrame object
#' @examples
#' data(mallard)
#' mallardUMF <- unMarkedFrame(mallard.y, siteCovs = mallard.site,
#'                            obsCovs = mallard.obs)
#' obsCovs(mallardUMF)
#' obsCovs(mallardUMF, matrices = TRUE)
#' @export
unMarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL,
                          obsNum = ncol(y), primaryNum = NULL) {
  if(class(obsCovs) == "list") {
    obsVars <- names(obsCovs)
    for(i in seq(length(obsVars))) {
      if(class(obsCovs[[i]]) != "matrix")
        stop("At least one element of obsCovs is not a matrix.")
      if(ncol(obsCovs[[i]]) != ncol(y) | nrow(obsCovs[[i]]) != nrow(y))
        stop("At least one matrix in obsCovs has incorrect number of dimensions.")
    }
    if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
    obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
  }

  if(("data.frame" %in% class(y)) |
			("cast_matrix" %in% class(y))) y <- as.matrix(y)

  ## add obsCov for the observation number (sampling occasion)
  ## name it obs
  obs = data.frame(obs = rep(1:obsNum, nrow(y)))
  if(!is.null(obsCovs))
    obsCovs <- as.data.frame(cbind(obsCovs,obs))
  else
    obsCovs <- obs

  if(is.null(primaryNum)) primaryNum <- 1

  umf <- new("unMarkedFrame", y = y, obsCovs = obsCovs,
             siteCovs = siteCovs, obsNum = obsNum,
             primaryNum = primaryNum)

  ## copy siteCovs into obsCovs
  if(!is.null(siteCovs)) {
    umf@obsCovs <- as.data.frame(cbind(umf@obsCovs,
                                       sapply(umf@siteCovs, rep,
                                              each = umf@obsNum)))
  }
  return(umf)
}

#' @export
setMethod("show", "unMarkedFrame",
          function(object) {
            ## print y
            cat("Observation matrix:\n")
            print(object@y)

            ## site covariates
            if(!is.null(object@siteCovs)) {
              cat("\nSite level covariates:\n")
              print(object@siteCovs)
            }

            if(!is.null(object@obsCovs)) {
              cat("\nWithin-site covariates:\n")
              obsCovs <- object@obsCovs
              M <- nrow(object@y)
              J <- object@obsNum
              for(i in seq(length=ncol(obsCovs))) {
                cat("\n",colnames(obsCovs)[i],":\n", sep="")
                print(matrix(obsCovs[,i], M, J, byrow = TRUE))
              }
            }
          })

#' Extractor for site level covariates
#' @param umf an unMarkedFrame
#' @return a data frame containing the site level covariates.
#' @export
siteCovs <- function(umf) {
  return(umf@siteCovs)
}

#' Extractor for observation level covariates
#' @param umf an unMarkedFrame
#' @param matrices logical indicating whether to return the M * obsNum row data frame (default)
#'  or a list of M x obsNum matrices (matrices = TRUE).
#' @return either a data frame (default) or a list of matrices (if matrices = TRUE).
#' @export
obsCovs <- function(umf, matrices = FALSE) {
  M <- nrow(umf@y)
  J <- umf@obsNum
  if(matrices) {
    value <- list()
    for(i in seq(length=length(umf@obsCovs))){
      value[[i]] <- matrix(umf@obsCovs[,i], M, J, byrow = TRUE)
    }
    names(value) <- names(umf@obsCovs)
  } else {
    value <- umf@obsCovs
  }
  return(value)
}

# Class to store actual parameter estimates
#' @export
setClass("unMarkedEstimate",
    representation(name = "character",
        estimates = "numeric",
        covMat = "matrix",
        invlink = "character",
        invlinkGrad = "character"))

setClass("unMarkedEstimateLinearComb",
    representation(originalEstimate = "unMarkedEstimate",
        coefficients = "numeric"),
    contains = "unMarkedEstimate")

setClass("unMarkedEstimateBackTransformed",
    representation(transformation = "character"),
    contains = "unMarkedEstimateLinearComb")


#setClass("unMarkedLogitEstimate",
#    representation(),
#    contains = "unMarkedEstimate")

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

##' @export
#setGeneric("linearComb",
#    function(obj, contrast) {
#      standardGeneric("linearComb")
#    })
#
##' @export
#setMethod("linearComb",
#    signature(obj = "unMarkedEstimate", contrast = "missing"),
#    function(obj, contrast) {
#
#    })


### Compute standard error on the unconstrained (transformed) scale.
### this gets all of the std errors
##' @export
#setMethod("standardError",
#    signature(obj = "unMarkedEstimate", contrast = "missing"),
#    definition = function(obj) {
#      SE <- sqrt(diag(obj@covMat))
#      names(SE) <- names(obj@estimates)
#      SE
#    })
#
#
### Provide a contrast and use the delta method to
### estimate the SE error on the constrained (natural) scale.
##' @export
#setMethod("standardError",
#    signature(obj = "unMarkedEstimate", contrast = "numeric"),
#    definition = function(obj, contrast) {
#      as.numeric(obj@invlinkGrad(t(contrast) %*% obj@estimates) *
#          sqrt(t(contrast) %*% obj@covMat %*% contrast))
#    })

## Compute contrasts and parameters and place on natural scale
#' @export
setGeneric("linearComb",
    function(obj, coefficients) {
      standardGeneric("linearComb")
    })

#' Compute linear combinations of estimates in unMarkedEstimate objects.
#'
#' This function computes the linear combination of parameter estimates in
#' \code{obj} given by the coefficient vector.  This may be useful to estimate
#' quantities of interest from a fitted model or to test a hypothesis.
#'
#' @param obj an unMarkedEstimate object
#' @param coefficients vector of same length as obj
#' @return an unMarkedEstimate object
#' @export
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

#' Transform an object to it's natural scale.
#'
#' @param obj object to be transformed
#' @export
setGeneric("backTransform",
    function(obj) {
      standardGeneric("backTransform")
    })

#' Transform an unMarkedEstimate object to it's natural scale.
#'
#' The transformation is determined by the invlink and invlinkGrad slots
#' in \code{obj}.  These slots specify the name of a many-to-one function and
#' its gradient respectively.  The delta method is used to compute the transformed
#' estimate.
#'
#' @param unMarkedEstimate object to be transformed
#' @return an unMarkedEstimate object representing the transformed estimate.
#' This object has invlink and invlinkGrad slots as the identity function.
#' @export
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

#' Compute standard error of an object.
#'
#' @param obj object whose standard error is computed.
#' @export
setGeneric("SE",
    def = function(obj) {
      standardGeneric("SE")
    })

#' Compute standard error of an unMarkedEstimate object.
#'
#' This function computes the large-sample standard error from the inverse of the
#' hessian matrix.
#'
#' @param obj unMarkedEstimate whose standard error is returned
#' @return vector of the standard error(s) of estimates in obj
#' @export
setMethod("SE",
    signature(obj = "unMarkedEstimate"),
    function(obj) {
      sqrt(diag(obj@covMat))
    })



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
#' A Class to store fit results from unMarkedFrames.
#' @export
setClass("unMarkedFit",
    representation(fitType = "character",
        stateformula = "formula",
        detformula = "formula",
        data = "unMarkedFrame",
        stateEstimates = "unMarkedEstimate",
        detEstimates = "unMarkedEstimate",
        AIC = "numeric",
        hessian = "matrix"))

# constructor for unMarkedFit objects
unMarkedFit <- function(fitType,stateformula, detformula,
                        data, stateEstimates,
                        detEstimates, AIC, hessian) {
  umfit <- new("unMarkedFit", fitType = fitType,
               stateformula = stateformula,
               detformula = detformula, data = data,
               stateEstimates = stateEstimates,
               detEstimates = detEstimates, AIC = AIC,
               hessian = hessian)

  return(umfit)
}

#' @export
setMethod("show", "unMarkedFit",
          function(object) {
            cat("\nCall:\n")
            cat(object@fitType,"(stateformula = ~ ",
                as.character(object@stateformula)[2],
                ", detformula = ~ ",
                as.character(object@detformula)[2],")\n\n", sep = "")

            show(object@stateEstimates)

            cat("\n")

            show(object@detEstimates)

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


# TODO: make parent class unMarkedEstimate
# TODO: make binomial detection child class
# TODO: make binomial occ child class
# TODO: make poisson abundance child class
# TODO: make show method for each of these classes.
# TODO: make unMarkedFit show that calls the respective children.
# TODO: separate unmarkedFit class for each model type?... would contain different estimate types