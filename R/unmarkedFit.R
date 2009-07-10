#' @include unmarkedEstimate.R
#' @include unmarkedFrame.R
#' @include classes.R
{}

# Class to store unMarked model fit information
#
# slot fitType Name of the model that was fit.
# @slot stateformula The abundance/occupancy formula.
# @slot detformula The formula governing the detection process.
# @slot data The unmarkedFrame containing the data that was fit.
# @slot stateMLE The MLE of the abundance/occupancy parameters
# @slot stateSE The standard errors of the state MLE's.
# @slot detMLE MLE of the detection parameters.
# @slot detSE Standard errors of the detection MLE's.
# @slot AIC The AIC of this model fit.
# @slot negLogLike The negative log likelihood of the fitted model.
# A Class to store fit results from unmarkedFrames.
#' @export
setClass("unmarkedFit",
    representation(fitType = "character",
        call = "call",
				formula = "formula",
        data = "unmarkedFrame",
        estimates = "unmarkedEstimateList",
        AIC = "numeric",
				opt = "list",
        negLogLike = "numeric",
				nllFun = "function"))



# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call, formula,
    data, estimates, AIC, opt, negLogLike, nllFun) {
  umfit <- new("unmarkedFit", fitType = fitType,
      call = call, formula = formula, data = data,
      estimates = estimates, AIC = AIC,
      opt = opt, negLogLike = negLogLike, nllFun = nllFun)

  return(umfit)
}

# @export
setMethod("show", "unmarkedFit",
    function(object) {
      cat("\nCall:\n")
#      cat(object@fitType,"(stateformula = ~ ",
#          as.character(object@stateformula)[2],
#          ", detformula = ~ ",
#          as.character(object@detformula)[2],")\n\n", sep = "")
      print(object@call)

      cat("\n")

      show(object@estimates)

      #      show(object@stateEstimates)
#
#      cat("\n")
#
#      show(object@detEstimates)

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


#' Compute linear combinations of estimates in unmarkedFit objects.
#'
#' This function computes the linear combination of parameter estimates in
#' \code{obj} given by the coefficient vector.  The user must
#' specify whether the detection or state covariates should be considered via the
#' \code{whichEstimate} argument.  This may be useful to estimate
#' quantities of interest from a fitted model or to test a hypothesis.
#'
#' @name linearComb-unmarkedFit
#' @aliases linearComb,unmarkedFit-method
#' @param obj an unmarkedFit object
#' @param coefficients vector of same length as obj
#' @param whichEstimate character, either "state" or "det"
#' @return an unmarkedEstimate object
setMethod("linearComb",
    signature(obj = "unmarkedFit", coefficients = "matrixOrVector"),
    function(obj, coefficients, whichEstimate) {
      stopifnot(!missing(whichEstimate))
      stopifnot(whichEstimate %in% names(obj))
      estimate <- obj@estimates[whichEstimate]
      linearComb(estimate, coefficients)
    })

setMethod("backTransform", "unmarkedFit",
		function(obj, whichEstimate) {
			est <- obj[whichEstimate]
			if(length(est@estimates) == 1) {
				lc <- linearComb(est, 1)
				return(backTransform(lc))
			} else {
				stop('Cannot directly backTransform an unmarkedEstimate with length > 1.')
			}
		})

setMethod("[",
    "unmarkedFit",
    function(x, i, j, drop) {
      x@estimates[i]
    })

setMethod("names", "unmarkedFit",
    function(x) {
      names(x@estimates)
    })
    

#' @name predict-unmarkedFit
#' @aliases predict-umDistsampFit, unmarkedFit-method
#' @examples
#' 
#' data(mallard)
#' mallardUMF <- unmarkedFrame(mallard.y, siteCovs = mallard.site,
#' obsCovs = mallard.obs)
#' fm.mallard <- pcount(~ length + elev + forest, ~ ivel + date, mallardUMF)
#'
#' # Predict from fitted model using original data
#' predict(fm.mallard, type="state", backTran=TRUE)
#' predict(fm.mallard, type="det", backTran=FALSE)
#' 
#' # Create newdata with only 'length' varying, then predict and plot
#' nd <- data.frame(length=rnorm(10), elev=1, forest=0, ivel=-2, date=3)
#'
#' (Elam <- predict(fm.mallard, type="state", newdata=nd))
#' 
#' with(data.frame(Elam), {
#' 	plot(nd$length, Predicted, ylim=c(0, 25))
#' 	segments(nd$length, Predicted+SE, nd$length, Predicted-SE)
#' 	})
#' 
#' @exportMethod predict
setMethod("predict", "unmarkedFit", 
		function(object, type, newdata=NULL, backTran=TRUE, ...) {
			if(is.null(newdata))
				newdata <- object@data
			formula <- object@formula
#			stateformula <- object@stateformula
#			detformula <- object@detformula
			cls <- class(newdata)
			switch(cls, 
					unmarkedFrame = {
						designMats <- getDesign2(formula, newdata)
						switch(type, 
								state = X <- designMats$X,
								det = X <- designMats$V)
					},
					data.frame = {
						switch(type, 
								state = {
									detformula <- as.formula(formula[[2]])
									stateformula <- as.formula(paste("~",formula[3],sep=""))
									Terms <- delete.response(terms(stateformula))
									mf <- model.frame(Terms, newdata)
									X <- model.matrix(Terms, mf)
									},
								det = X <- model.matrix(detformula, newdata))
					})
			out <- matrix(NA, nrow(X), 2, dimnames=list(NULL, c("Predicted", "SE")))
			lc <- linearComb(object, X, type)
			if(backTran) lc <- backTransform(lc)
			out[,1] <- coef(lc)
			out[,2] <- SE(lc)
			return(out)
		}
)


#' @importFrom graphics plot
#' @importFrom stats coef vcov predict update profile


#' @export
setMethod("coef", "unmarkedFit",
		function(object, type, altNames = TRUE) {
			if(missing(type)) {
				co <- lapply(object@estimates@estimates, function(x) coef(x, altNames=altNames))
				names(co) <- NULL
				co <- unlist(co)
			} else {
				co <- coef(object[type], altNames=altNames)
			}
			co
		})

#' @export 
setMethod("vcov", "unmarkedFit",
		function(object, type, altNames = TRUE) {
			if(missing(type)) {
				v <- solve(hessian(object))
				rownames(v) <- colnames(v) <- names(coef(object, altNames=altNames))
			} else {
				v <- vcov(object[type])
				rownames(v) <- colnames(v) <- names(coef(object, type, altNames=altNames))
			}
			v
		})

#' @export
setMethod("confint", "unmarkedFit",
		function(object, parm, level = 0.95, type, method = c("normal", "profile")) {
			method <- match.arg(method)
			if(missing(type)) stop(paste("Must specify type as one of (", paste(names(object),collapse=", "),").",sep=""))
			if(missing(parm)) parm <- 1:length(object[type]@estimates)
			if(method == "normal") {
				callGeneric(object[type],parm = parm, level = level)
			} else {
				nllFun <- nllFun(object)
				ests <- mle(object)
				nP <- length(parm)
				ci <- matrix(NA, nP, 2)

				## create table to match parm numbers with est/parm numbers
				types <- names(object)
				numbertable <- data.frame(type = character(0), num = numeric(0))
				for(i in seq(length=length(types))) {
					length.est <- length(object[i]@estimates)
					numbertable <- rbind(numbertable, data.frame(type = rep(types[i], length.est), num = seq(length=length.est)))
				}
				parm.fullnums <- which(numbertable$type == type & numbertable$num %in% parm)
				
				for(i in seq(length=nP)) {
					cat("Profiling parameter",i,"of",nP,"...")
					se <- SE(object[type])
					whichPar <- parm.fullnums[i]
					ci[i,] <- profileCI(nllFun, whichPar=whichPar, MLE=ests, interval=ests[whichPar] + 10*se[i]*c(-1,1), level=level)
					cat(" done.\n")
				}
				rownames(ci) <- names(coef(object[type]))[parm]
				colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
				return(ci)
			}
		})

#' @export
setMethod("profile", "unmarkedFit",
		function(fitted, type, parm, seq) {
			stopifnot(length(parm) == 1)
			MLE <- mle(fitted)
			SE(fitted[type])
			nPar <- length(mle(fitted))
			nll <- nllFun(fitted)
			
			## create table to match parm numbers with est/parm numbers
			types <- names(fitted)
			numbertable <- data.frame(type = character(0), num = numeric(0))
			for(i in seq(length=length(types))) {
				length.est <- length(fitted[i]@estimates)
				numbertable <- rbind(numbertable, data.frame(type = rep(types[i], length.est), num = seq(length=length.est)))
			}
			parm.fullnums <- which(numbertable$type == type & numbertable$num == parm)
			
			f <- function(value) {
				fixedNLL <- genFixedNLL(nll, parm.fullnums, value)
				mleRestricted <- optim(rep(0,nPar), fixedNLL)$value
				mleRestricted
			}
			prof.out <- sapply(seq, f)
			prof.out <- cbind(seq, prof.out)
			new("profile", prof = prof.out)
		})

setMethod("hessian", "unmarkedFit",
		function(object) {
			object@opt$hessian
		})


setMethod

setGeneric("nllFun", function(object) standardGeneric("nllFun"))

#' @export 
setMethod("nllFun", "unmarkedFit", function(object) object@nllFun)

setGeneric("mle", function(object) standardGeneric("mle"))

#' @export 
setMethod("mle", "unmarkedFit", function(object) object@opt$par)

setClass("profile",
		representation(prof = "matrix"))
setMethod("plot", c("profile", "missing"),
		function(x) {
			plot(x@prof[,1], x@prof[,2], type = "l")
		})