#' @include classes.R
{}

############# VALIDATION FUNCTIONS #############################################

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

############ DATA CLASSES #########################################################

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
        obsToY = "optionalMatrix"),
    validity = validunmarkedFrame)



## a class for multi-season data
setClass("unmarkedMultFrame",
		representation(numPrimary = "numeric"),
		contains="unmarkedFrame")



## a class for distance sampling data
#' @exportClass unmarkedFrameDS 
setClass("unmarkedFrameDS", 
		representation(
				dist.breaks = "numeric",
				tlength = "numeric",
				survey = "character",
				unitsIn = "character", 
				distClassAreas = "matrix"),
		contains = "unmarkedFrame",
		validity = function(object) 
		{
			J <- numY(object)
			db <- object@dist.breaks
			if(J == length(db) - 1 && db[1] == 0) TRUE
			else "ncol(y) must equal length(dist.breaks)-1 and 
						dist.breaks[1] must equal 0"
		})

#' @export
setClass("unmarkedFrameOccu",
		contains = "unmarkedFrame")

#' @export
setClass("unmarkedFramePCount",
		contains = "unmarkedFrame")

#' @export
setClass("unmarkedFrameMPois",
		representation(samplingMethod = "character",
				piFun = "function"),
		contains = "unmarkedFrame")

################### CONSTRUCTORS ##########################################################

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
    obsToY) {

  if(class(obsCovs) == "list") {
    obsVars <- names(obsCovs)
    for(i in seq(length(obsVars))) {
      if(!(class(obsCovs[[i]]) %in% c("matrix", "data.frame")))
        stop("At least one element of obsCovs is not a matrix or data frame.")
      if(ncol(obsCovs[[i]]) != ncol(y) | nrow(obsCovs[[i]]) != nrow(y))
        stop("At least one matrix in obsCovs has incorrect number of dimensions.")
    }
    if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
    obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
  }

  if(("data.frame" %in% class(y)) |
      ("cast_matrix" %in% class(y))) y <- as.matrix(y)

	## if no obsToY is supplied, assume y <-> obsCov
	#if(is.null(obsToY)) obsToY <- diag(ncol(y))  assuming the obsToY can be dangerous... keep as NULL if not supplied.
	if(missing(obsToY)) obsToY <- NULL
	
  umf <- new("unmarkedFrame", y = y, obsCovs = obsCovs,
      siteCovs = siteCovs, obsToY = obsToY)

  return(umf)
}


# Constructor
#' @export unmarkedFrameDS
unmarkedFrameDS <- function(y, siteCovs = NULL, dist.breaks, tlength, survey,
		unitsIn, distClassAreas = NULL)
{
	if(is.null(distClassAreas))
		distClassAreas <- matrix(NA)
	if(missing(tlength) & survey == "point")
		tlength <- numeric(0)
	umfds <- new("unmarkedFrameDS", y = y, obsCovs = NULL,
			siteCovs = siteCovs, dist.breaks = dist.breaks, tlength = tlength,
			survey = survey, unitsIn = unitsIn, distClassAreas = distClassAreas,
			obsToY = matrix(1, 1, ncol(y)))
	return(umfds)
}


#' @export
unmarkedFrameOccu <- function(y, siteCovs = NULL, obsCovs = NULL) {
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J))
	umf <- as(umf, "unmarkedFrameOccu")
	umf
}


#' @export
unmarkedFramePCount <- function(y, siteCovs = NULL, obsCovs = NULL) {
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J))
	umf <- as(umf, "unmarkedFramePCount")
	umf
}


#' @export
unmarkedFrameMPois <- function(y, siteCovs = NULL, obsCovs = NULL, obsToY, piFun) {
	if(missing(obsToY)) stop("obsToY is required for multinomial-Poisson data.")
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY)
	umf <- as(umf, "unmarkedFrameMPois")
	umf@piFun <- piFun
	umf@samplingMethod <- as.character(quote(piFun))
	umf
}

################ SHOW METHODS ###########################################################

#' @export
setMethod("show", "unmarkedFrame",
		function(object) {
			obsCovs <- obsCovs(object)
			siteCovs <- siteCovs(object)
			y <- y(object)
			colnames(y) <- paste("y",1:ncol(y),sep=".")
			if(is.null(obsToY(object))) {
				obsNum <- ncol(y)
			} else {
				obsNum <- obsNum(object)
			}
			if(is.null(siteCovs)) siteCovs <- matrix(0,nrow(y),0)
			if(is.null(obsCovs)) {
				obsCovs <- matrix(0,nrow(y),0)
			} else {
				obsCovs <- data.frame(lapply(obsCovs, function(x) matrix(x, nrow(y), obsNum,byrow=T)))
			}
			df <- data.frame(y, siteCovs, obsCovs)
			cat("Data frame representation of unmarkedFrame object.\n")
			print(df)
		})

############################ EXTRACTORS #########################################################

#' Extractor for site level covariates
#' @param umf an unmarkedFrame
#' @return a data frame containing the site level covariates.
#' @export
setGeneric("siteCovs", function(object,...) standardGeneric("siteCovs"))

setMethod("siteCovs", "unmarkedFrame",
		function(object) {
			return(object@siteCovs)
		})

#
##' Extractor for observation level covariates
##' @param umf an unmarkedFrame
##' @param matrices logical indicating whether to return the M * obsNum row data frame (default)
##'  or a list of M x obsNum matrices (matrices = TRUE).
##' @return either a data frame (default) or a list of matrices (if matrices = TRUE).

#' @export
setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))

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
setGeneric("obsNum", function(object) standardGeneric("obsNum"))

setMethod("obsNum", "unmarkedFrame", function(object) nrow(object@obsToY))

#' @export 
setGeneric("numSites", function(object) standardGeneric("numSites"))

setMethod("numSites", "unmarkedFrame", function(object) nrow(object@y))

#' @export 
setGeneric("numY", function(object) standardGeneric("numY"))

setMethod("numY", "unmarkedFrame", function(object) ncol(object@y))

#' @export 
setGeneric("obsToY", function(object) standardGeneric("obsToY"))

setMethod("obsToY", "unmarkedFrame", function(object) object@obsToY)

#' @export 
setGeneric("obsToY<-", function(object, value) standardGeneric("obsToY<-"))

setReplaceMethod("obsToY", "unmarkedFrame", function(object, value) {
			object@obsToY <- value
			object
		})

#' @export
setGeneric("y", function(object) standardGeneric("y"))

setMethod("y", "unmarkedFrame", function(object) object@y)


################################### SUMMARY METHODS #############################################

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


setAs("data.frame", "unmarkedFrame", function(from) {
			umf <- formatWide(from)
			umf
		})

################################# PLOT METHODS ############################################
# TODO:  come up with nice show/summary/plot methods for each of these data types.


################################# EXTRACTORS ###############################################

# i is the vector of sites to extract
#' @exportMethod "["
setMethod("[", c("unmarkedFrame","numeric"),
		function(x, i, j = "missing", drop = "missing") {  
			y <- y(x)[i,]
			if (length(i) == 1) {
				y <- t(y)
			}
			siteCovs <- siteCovs(x)[i,]
			obsCovs <- obsCovs(x)
			obs.site.inds <- rep(1:numSites(x), each = obsNum(x))
			obsCovs <- obsCovs[obs.site.inds %in% i,]
			umf <- x
			umf@y <- y
			umf@siteCovs <- siteCovs
			umf@obsCovs <- obsCovs
			umf
		})