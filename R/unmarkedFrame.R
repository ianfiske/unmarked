#' @include classes.R
{}

validunmarkedFrame <- function(object) {
  errors <- character(0)
  M <- nrow(object@y)
  if(!is.null(object@siteCovs))
    if(nrow(object@siteCovs) != M)
      errors <- c(errors, "siteCovData does not have same size number of sites as y.")
  if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
    if(nrow(object@obsCovs) != M*obsNum(object))
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
# @slot obsToY matrix that describes how the observations relate to y (see Details).
# @slot primaryNum integer number of seasons (1 for single season).

#' Class to hold data for analyses in unmarked.
#'
#' @export
setClass("unmarkedFrame",
    representation(y = "matrix",
        obsCovs = "optionalDataFrame",
        siteCovs = "optionalDataFrame",
        obsToY = "optionalMatrix",
        primaryNum = "numeric"),
    validity = validunmarkedFrame)

#' Constructor for unmarkedFrames.
#'
#' unmarkedFrame is the S4 class that holds data structures to be passed to the model-fitting functions in unMarked.
#'
#' An unmarkedFrame contains the observations (\code{y}), covariates measured at the observation level (\code{obsCovs}), and covariates measured at the site level (\code{siteCovs}).
#' For a data set with M sites and J observations at each site, y is an M x J matrix.
#' \code{obsCovs} and \code{siteCovs} are both data frames (see \link{data.frame}).  \code{siteCovs} has M rows so that each row contains the covariates for the corresponding sites.
#' \code{obsCovs} has M*obsNum rows so that each covariates is ordered by site first, then observation number.  Missing values are coded with \code{NA} in any of y, siteCovs, or obsCovs.
#'
#' Additionally, unmarkedFrames contain metadata, obsNum and primaryNum.  obsNum is the number of observations measured at each site. primaryNum is the number of seasons in a robust design sampling scheme.
#' Typically, these can be automatically determined by the constructor.  If not specified, obsNum is taken to be the number of columns in y and primaryNum is taken to be 1.
#' However, for certain situations, these must be supplied.  For example, double observer sampling, y has 3 columns corresponding the observer 1, observer 2, and both, but there were only two independent observations.
#' In this situation, y has 3 columns, but obsNum must be specified as 2.  This flexibility is currenty only used in the function \link{multinomPois}.
#'
#' For convenience, \code{obsCovs} can be a list of M x obsNum matrices, with each one corresponding to an observation level covariate.
#'
#' All site-level covariates are automatically copied to obsCovs so that site level covariates are available at the observation level.
#'
#' @title Create an unmarkedFrame.
#' @aliases unmarkedFrame obsCovs siteCovs
#' @param y A matrix of the observed measured data.
#' @param obsCovs Dataframe of covariates that vary within sites.
#' @param siteCovs Dataframe of covariates that vary at the site level.
#' @param obsNum Number of independent observations.
#' @param primaryNum Number of primary time periods (seasons in the multiseason model).
#' @return an unmarkedFrame object
#' @examples
#' data(mallard)
#' mallardUMF <- unmarkedFrame(mallard.y, siteCovs = mallard.site,
#'                            obsCovs = mallard.obs)
#' obsCovs(mallardUMF)
#' obsCovs(mallardUMF, matrices = TRUE)
#' @export
unmarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL,
    obsToY, primaryNum = NULL) {

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

  if(is.null(primaryNum)) primaryNum <- 1
	
	## if no obsToY is supplied, assume y <-> obsCov
	#if(is.null(obsToY)) obsToY <- diag(ncol(y))  assuming the obsToY can be dangerous... keep as NULL if not supplied.
	if(missing(obsToY)) obsToY <- NULL
	
  umf <- new("unmarkedFrame", y = y, obsCovs = obsCovs,
      siteCovs = siteCovs, obsToY = obsToY,
      primaryNum = primaryNum)

  return(umf)
}

#' @export
setMethod("show", "unmarkedFrame",
		function(object) {
			cat("unmarkedFrame object\n\n")
			cat(nrow(object@y),"sites\n")
			obsPres <- !is.null(object@obsCovs)
			sitePres <- !is.null(object@siteCovs)
			if(obsPres)  {
				cat("observation-level covariates present.\n") 
			} else {
				cat("No observation-level covariates present.\n")
			}
			if(sitePres)  {
				cat("site-level covariates present.\n")
			} else {
				cat("No site-level covariates present.\n")
			}
		})

#' Extractor for site level covariates
#' @param umf an unmarkedFrame
#' @return a data frame containing the site level covariates.
#' @export
siteCovs <- function(umf) {
  return(umf@siteCovs)
}


#
##' Extractor for observation level covariates
##' @param umf an unmarkedFrame
##' @param matrices logical indicating whether to return the M * obsNum row data frame (default)
##'  or a list of M x obsNum matrices (matrices = TRUE).
##' @return either a data frame (default) or a list of matrices (if matrices = TRUE).

setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))

#' @export
setMethod("obsCovs", "unmarkedFrame", 
		function(object, matrices = FALSE) {
			M <- numSites(object)
			R <- obsNum(object)
			if(matrices) {
				value <- list()
				for(i in seq(length=length(object@obsCovs))){
					value[[i]] <- matrix(object@obsCovs[,i], M, R, byrow = TRUE)
				}
				names(value) <- names(object@obsCovs)
			} else {
				value <- object@obsCovs
			}
			return(value)
		})

#' @export
setMethod("summary","unmarkedFrame",
    function(object,...) {
      cat("unmarkedFrame Object\n\n")
      cat(nrow(object@y), "sites\n")
      cat("Maximum observations per site:",obsNum(object),"\n\n")
      cat("Distribution of observations per site:")
      stem(rowSums(!is.na(object@y)), scale=0.5)
      cat("Tabulation of y observations:")
      print(table(object@y, exclude=NULL))
      if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:\n")
        print(summary(object@siteCovs))
      }
      if(!is.null(object@obsCovs)) {
        cat("\nObservation-level covariates:\n")
        print(summary(object@obsCovs))
      }
    })

setGeneric("obsNum", function(object) standardGeneric("obsNum"))

#' @export 
setMethod("obsNum", "unmarkedFrame", function(object) nrow(object@obsToY))

setGeneric("numSites", function(object) standardGeneric("numSites"))

#' @export 
setMethod("numSites", "unmarkedFrame", function(object) nrow(object@y))

setGeneric("numY", function(object) standardGeneric("numY"))

#' @export 
setMethod("numY", "unmarkedFrame", function(object) ncol(object@y))

setGeneric("obsToY", function(object) standardGeneric("obsToY"))

#' @export 
setMethod("obsToY", "unmarkedFrame", function(object) object@obsToY)

setGeneric("obsToY<-", function(object, value) standardGeneric("obsToY<-"))

#' @export 
setReplaceMethod("obsToY", "unmarkedFrame", function(object, value) {
			object@obsToY <- value
			object
		})

setGeneric("y", function(object) standardGeneric("y"))

#' @export
setMethod("y", "unmarkedFrame", function(object) object@y)

setAs("data.frame", "unmarkedFrame", function(from) {
			umf <- formatWide(from)
			umf
		})