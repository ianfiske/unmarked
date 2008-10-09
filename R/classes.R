setClassUnion("optionalDataFrame", c("data.frame","NULL"))

validUnMarkedFrame <- function(object) {
  errors <- character(0)
  M <- nrow(object@y)
  J <- object@obsNum
  if(!is.null(object@siteCovs))
    if(nrow(object@siteCovs) != M)
      errors <- c(errors, "siteCovData does not have same size number of sites as y.")
  if(!is.null(object@obsCovs))
    if(nrow(object@obsCovs) != M*J)
      errors <- c(errors, "obsCovData does not have M*J rows.")
  if(length(errors) == 0)
    TRUE
  else
    errors
}

#' Class to hold data for analyses in unmarked.
#'
#' @slot y A matrix of the observed measured data.
#' @slot obsCovData Dataframe of covariates that vary within sites.
#' @slot siteCovData Dataframe of covariates that vary at the site level.
#' @slot obsNum Number of observations per site
setClass("unMarkedFrame",
         representation(y = "matrix",
                        obsCovs = "optionalDataFrame",
                        siteCovs = "optionalDataFrame",
                        obsNum = "numeric"),
         validity = validUnMarkedFrame)

#' Constuctor function to create an unmarkedFrame.
#'
#' @param y A matrix of the observed measured data.
#' @param obsCovData Dataframe of covariates that vary within sites.
#' @param siteCovData Dataframe of covariates that vary at the site level.
#' @export
unMarkedFrame <- function(y, obsCovs = NULL, siteCovs = NULL,
                          obsNum = NULL) {
  if(is.null(obsNum)) obsNum <- ncol(y)
  return(new("unMarkedFrame", y = y, obsCovs = obsCovs,
             siteCovs = siteCovs, obsNum = obsNum))
}

#' Summary statistics for an unMarkedFrame
#'
#' @param object An unMarkedFrame to summarize.
#' @export
setMethod("summary", "unMarkedFrame",
          function(object) {
            ## get number of sites
            M <- nrow(object@y)
            print(paste("Number of sites:",M))
          })

#' Print an unMarkedFrame
#'
#' @param object
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
                cat("\n",colnames(obsCovs)[i],":\n")
                print(matrix(obsCovs[,i], M, J, byrow = TRUE))
              }
            }
          })

#' Extract the site covariates from an object
#' @param object An object to extract site covariates from.
setGeneric("siteCovs",
           function(object) {
             standardGeneric("siteCovs")
           })

#' Extract site covariates from an unMarkedFrame
#' @param object The unMarkedFrame whose covariates need to be extracted.
#' @export
setMethod("siteCovs", "unMarkedFrame",
          function(object) {
            return(object@siteCovs)
          })
          
setGeneric("obsCovs",
           function(object, matrices = NULL) {
             standardGeneric("obsCovs")
           })

setMethod("obsCovs",
          signature(object = "unMarkedFrame"),
          function(object, matrices = FALSE) {
            M <- nrow(object@y)
            J <- ncol(object@y)
            if(matrices) {
              value <- list()
              for(i in seq(length=length(object@obsCovs))){
                value[[i]] <- matrix(object@obsCovs[,i], M, J, byrow = TRUE)
              }
              names(value) <- names(object@obsCovs)
            } else {
              value <- object@obsCovs
            }
            value
          })
            


#' Class to store unMarked model fit information
#'
#' @slot fitType Name of the model that was fit.
#' @slot stateformula The abundance/occupancy formula.
#' @slot detformula The formula governing the detection process.
#' @slot data The unMarkedFrame containing the data that was fit.
#' @slot stateMLE The MLE of the abundance/occupancy parameters
#' @slot stateSE The standard errors of the state MLE's.
#' @slot detMLE MLE of the detection parameters.
#' @slot detSE Standard errors of the detection MLE's.
#' @slot AIC The AIC of this model fit.
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


#' Constructor function for unMarkedFit objects
#'
#' @param fitType
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


setMethod("show", "unMarkedFit",
          function(object) {
            cat("\nCall:\n")
            cat(object@fitType,"(stateformula = ~ ",
                as.character(object@stateformula)[2],
                ", detformula = ~ ",
                as.character(object@detformula)[2],")\n", sep = "")

            cat("\nState coefficients:\n")
            show(object@stateMLE)

            cat("\nDetection coefficients:\n")
            show(object@detMLE)
          })


setMethod("summary", "unMarkedFit",
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
          })
