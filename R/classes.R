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

#' Class to hold data for analyses in unmarked.
#'
#' @slot y A matrix of the observed measured data.
#' @slot obsCovData Dataframe of covariates that vary within sites.
#' @slot siteCovData Dataframe of covariates that vary at the site level.
#' @slot obsNum Number of observations per site. For most models, this
#' can be taken to be the number of columns in y.  But this is not always
#' the case.  For example, double observer: y has 3 columns, but only 2
#' independent observations were taken at each site.
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
#' @param obsNum Number of independent observations. 
#' @export
unMarkedFrame <- function(y, obsCovs = NULL, siteCovs = NULL,
                          obsNum = ncol(y)) {
  umf <- new("unMarkedFrame", y = y, obsCovs = obsCovs,
             siteCovs = siteCovs, obsNum = obsNum)
  ## copy siteCovs into obsCovs
  umf@obsCovs <- as.data.frame(cbind(umf@obsCovs,
                                     sapply(umf@siteCovs, rep,
                                            each = umf@obsNum)))
  return(umf)
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

siteCovs <- function(umf) {
  return(umf@siteCovs)
}
          
obsCovs <- function(umf, matrices = FALSE) {
  M <- nrow(umf@y)
  J <- ncol(umf@y)
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

            cat("\nAIC:", object@AIC,"\n")
          })


