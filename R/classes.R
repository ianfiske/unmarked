setClassUnion("optionalDataFrame", c("data.frame","NULL"))

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

# Class to hold data for analyses in unmarked.
#
# @slot y A matrix of the observed measured data.
# @slot obsCovData Dataframe of covariates that vary within sites.
# @slot siteCovData Dataframe of covariates that vary at the site level.
# @slot obsNum Number of observations per site. For most models, this
# can be taken to be the number of columns in y.  But this is not always
# the case.  For example, double observer: y has 3 columns, but only 2
# independent observations were taken at each site.
#' @export
setClass("unMarkedFrame",
         representation(y = "matrix",
                        obsCovs = "optionalDataFrame",
                        siteCovs = "optionalDataFrame",
                        obsNum = "numeric"),
         validity = validUnMarkedFrame)

#' Constuctor function to create an unmarkedFrame.
#'
#' This function takes observations (y) and covariates at the site level
#' (siteCovs) and observations level (obsCovs) and generates an
#' unMarkedFrame object.  This object can then be passed to any of the
#' unMarked functions.
#'
#' @param y A matrix of the observed measured data.
#' @param obsCovs Dataframe of covariates that vary within sites.
#' @param siteCovs Dataframe of covariates that vary at the site level.
#' @param obsNum Number of independent observations.
#' @return an unmarkedFrame object
#' @examples
#' data(mallard)
#' mallardUMF <- unMarkedFrame(mallard.y, siteCovs = mallard.site,
#'                            obsCovs = mallard.obs)
#' obsCovs(mallardUMF)
#' obsCovs(mallardUMF, matrices = TRUE)
#' @export
unMarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL,
                          obsNum = ncol(y)) {

  ## if obsCovs is a list of matrices, convert to a dataframe
  if(class(obsCovs) == "list") {
    obsVars <- names(obsCovs)
    for(i in seq(length(obsVars))) {
      if(class(obsCovs[[i]]) != "matrix")
        stop("At least one element of obsCovs is not a matrix.")
      if(ncol(obsCovs[[i]]) != ncol(y) | nrow(obsCovs[[i]]) != nrow(y))
        stop("At least one matrix in obsCovs has incorrect number of dimensions.")
    }
    obsCovs <- data.frame(lapply(obsCovs, as.vector))
  }

  if(class(y) == "data.frame") y <- as.matrix(y)

  ## add obsCov for the observation number (sampling occasion)
  ## name it obs
  obs = data.frame(obs = factor(rep(1:obsNum, nrow(y))))
  if(!is.null(obsCovs))
    obsCovs <- as.data.frame(cbind(obsCovs,obs))
  else
    obsCovs <- obs

  umf <- new("unMarkedFrame", y = y, obsCovs = obsCovs,
             siteCovs = siteCovs, obsNum = obsNum)
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

#' @export
siteCovs <- function(umf) {
  return(umf@siteCovs)
}
          
#' @export
obsCovs <- function(umf, matrices = FALSE) {
  M <- nrow(umf@y)
  J <- ncol(umf@obsNum)
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
            
# Class to store unMarked model fit information
#
# @slot fitType Name of the model that was fit.
# @slot stateformula The abundance/occupancy formula.
# @slot detformula The formula governing the detection process.
# @slot data The unMarkedFrame containing the data that was fit.
# @slot stateMLE The MLE of the abundance/occupancy parameters
# @slot stateSE The standard errors of the state MLE's.
# @slot detMLE MLE of the detection parameters.
# @slot detSE Standard errors of the detection MLE's.
# @slot AIC The AIC of this model fit.
#' @export
setClass("unMarkedFit",
         representation(fitType = "character",
                        stateformula = "formula",
                        detformula = "formula",
                        data = "unMarkedFrame",
                        stateMLE = "numeric",
                        stateSE = "numeric",
                        detMLE = "numeric",
                        detSE = "numeric",
                        AIC = "numeric"))

unMarkedFit <- function(fitType,stateformula, detformula,
                        data, stateMLE, stateSE,
                        detMLE, detSE, AIC) {
  umfit <- new("unMarkedFit", fitType = fitType,
               stateformula = stateformula,
               detformula = detformula, data = data,
               stateMLE = stateMLE, stateSE = stateSE,
               detMLE = detMLE, detSE = detSE, AIC = AIC)

  return(umfit)
}

#' @export
setMethod("show", "unMarkedFit",
          function(object) {
            cat("\nCall:\n")
            cat(object@fitType,"(stateformula = ~ ",
                as.character(object@stateformula)[2],
                ", detformula = ~ ",
                as.character(object@detformula)[2],")\n", sep = "")

            cat("\nState coefficients:\n")
            stateDF <- data.frame(Estimate = object@stateMLE,
                                  SE = object@stateSE)
            show(stateDF)
            
            cat("\nDetection coefficients:\n")
            detDF <- data.frame(Estimate = object@detMLE,
                                  SE = object@detSE)
            show(detDF)

            cat("\nAIC:", object@AIC,"\n")
          })


