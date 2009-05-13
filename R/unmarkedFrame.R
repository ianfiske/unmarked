#' @include classes.R
roxygen()

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
#' @aliases unMarkedFrame obsCovs siteCovs
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
  if(!is.matrix(y)) {
    stop("y must be a matrix.")
  }
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

##' @export
# setGeneric("summary")
#    function(object,...) {
#      standardGeneric("summary")
#    })

#' @export
setMethod("summary","unMarkedFrame",
    function(object,...) {
      cat("unMarkedFrame Object\n\n")
      cat(nrow(object@y), "sites\n")
      cat("Maximum observations per site:",object@obsNum,"\n\n")
      cat("Distribution of observations per site:")
      stem(rowSums(!is.na(umf@y)), scale=0.5)
      cat("Tabulation of y observations:")
      print(table(object@y, exclude=NULL))
      if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:")
        print(summary(object@siteCovs))
      }
      if(!is.null(object@obsCovs)) {
        cat("\nObservation-level covariates:\n")
        print(summary(object@obsCovs))
      }
    })