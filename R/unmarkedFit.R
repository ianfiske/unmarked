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
				sitesRemoved = "numeric",  # vector of indices of sites that were removed from analysis
        estimates = "unmarkedEstimateList",
        AIC = "numeric",
				opt = "list",
        negLogLike = "numeric",
				nllFun = "function"))



# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call, formula,
    data, sitesRemoved, estimates, AIC, opt, negLogLike, nllFun) {
  umfit <- new("unmarkedFit", fitType = fitType,
      call = call, formula = formula, data = data, sitesRemoved = sitesRemoved,
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
      
      if(!identical(object@opt$convergence, 0L))
    		warning("Model did not converge. Try providing starting values or
    			increasing maxit control argment.")
    			
    })
    
    
    
    
#' @export
setMethod("summary", "unmarkedFit", 
	function(object) {
		show(object)
		cat("\nSample size:", sampleSize(object))
		if(length(object@sitesRemoved) > 0)
			cat("\nSites removed:", object@sitesRemoved)
		cat("\noptim convergence code:", object@opt$convergence)
		cat("\noptim iterations:", object@opt$counts[1], "\n", "\n")
		if(!identical(object@opt$convergence, 0L))
    		warning("Model did not converge. Try providing starting values or
    			increasing maxit control argment.")
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
	function(object, type, newdata=NULL, backTran=TRUE, na.rm = TRUE, ...) 
	{
		if(is.null(newdata))
			newdata <- object@data
		formula <- object@formula
		if(inherits(newdata, "unmarkedFrame"))
			class(newdata) <- "unmarkedFrame"
		cls <- class(newdata)
		switch(cls, 
			unmarkedFrame = {
				designMats <- getDesign2(formula, newdata, na.rm = na.rm)
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
#' @importFrom stats coef vcov predict update profile simulate


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



#' @exportMethod update
setMethod("update", "unmarkedFit", 
	function(object, formula., ..., evaluate = TRUE) {
		call <- object@call
		if (is.null(call)) 
            stop("need an object with call slot")
        extras <- match.call(expand.dots = FALSE)$...
        if (!missing(formula.)) 
            call$formula <- update.formula(formula(object), formula.)
        if (length(extras) > 0) {
            existing <- !is.na(match(names(extras), names(call)))
            for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
            if (any(!existing)) {
                call <- c(as.list(call), extras[!existing])
                call <- as.call(call)
            	}
        	}
        if (evaluate) 
            eval(call, parent.frame())
        else call
		}
	)


setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))

#' @export 
setMethod("sampleSize", "unmarkedFit",
		function(object) {
			data <- getData(object)
			M <- numSites(data)
			M <- M - length(object@sitesRemoved)
			M
		})
	
setGeneric("getData", function(object) standardGeneric("getData"))	

#' @export 
setMethod("getData", "unmarkedFit",
		function(object) {
			object@data
		})


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




 		
		

################################# CHILD CLASSES ###############################


#' @exportClass unmarkedFitDS
setClass("unmarkedFitDS",
		representation(
				keyfun = "character",
				unitsOut = "character"),
		contains = "unmarkedFit")



#' @exportClass unmarkedFitPCount
setClass("unmarkedFitPCount", 
		representation(
				K = "numeric",
				mixture = "character"),
		contains = "unmarkedFit")
				
				


############################# CHILD CLASS METHODS ##############################

# Extract detection probs
setGeneric("getP", function(object, ...) standardGeneric("getP"))


#' @export
setMethod("getP", "unmarkedFitDS", function(object, na.rm = TRUE) 
{
	formula <- object@formula
	detformula <- as.formula(formula[[2]])
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	V <- designMats$V
	M <- nrow(y)
	J <- ncol(y)
	ppars <- coef(object, type = "det")
	d <- umf@dist.breaks
	survey <- umf@survey
	key <- object@keyfun
	switch(key, 
		halfnorm = {
			sigma <- exp(V %*% ppars)
			p <- sapply(sigma, function(x) cp.hn(d = d, s = x, survey = survey))
			}, 
		exp = {
			rate <- exp(V %*% ppars)
			p <- sapply(rate, function(x) cp.exp(d = d, r = x, survey = survey))
			}, 
		hazard = {
			shape <- exp(V %*% ppars[-length(ppars)])
			scale <- exp(ppars[length(ppars)])
			p <- sapply(sigma, function(x) cp.haz(d = d, shape = shape, 
				scale = scale, survey = survey))
			})
    p <- matrix(p, M, J, byrow = TRUE)
	return(p)
})





#' @export
setMethod("getP", "unmarkedFitPCount", function(object, na.rm = TRUE) 
{
	formula <- object@formula
	detformula <- as.formula(formula[[2]])
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	V <- designMats$V
	M <- nrow(y)
	J <- ncol(y)
    ppars <- coef(object, type = "det")
	p <- plogis(V %*% ppars)
	p <- matrix(p, M, J, byrow = TRUE)
	return(p)
})




setGeneric("simulate", function(object,...) {
			standardGeneric("simulate")
		})


#' @export
setMethod("simulate", "unmarkedFitDS", 
	function(object, nsim = 1, seed = NULL, na.rm=TRUE)
{
	formula <- object@formula
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	X <- designMats$X
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	lamParms <- coef(object, type = "state")
	lam <- as.numeric(exp(X %*% lamParms))
	lamvec <- rep(lam, each = J) * a
	pvec <- c(t(getP(object)))
	yList <- list()
	for(i in 1:nsim) {
		yvec <- rpois(M * J, lamvec * pvec)
		yList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
		}
	return(yList)
})


#' @export
setMethod("simulate", "unmarkedFitPCount", 
	function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
	formula <- object@formula
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	X <- designMats$X
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	allParms <- coef(object, altNames = FALSE)
	lamParms <- coef(object, type = "state")
	lam <- as.numeric(exp(X %*% lamParms))
	lamvec <- rep(lam, each = J) * a
	pvec <- c(t(getP(object)))
	mix <- object@mixture
	yList <- list()
	for(i in 1:nsim) {
		switch(mix, 
			P = yvec <- rpois(M * J, lamvec * pvec),
			NB = yvec <- rnbinom(M * J, size = exp(allParms["alpha"]), 
				mu = lambec * pvec)
			)
		yList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
		}
	return(yList)
})









############################## PARAMETRIC BOOTSTRAP ###########################

#' @exportMethod parboot
setGeneric("parboot",
    def = function(object, ...) {
      standardGeneric("parboot")
    })




#' @exportClass parboot
setClass("parboot",
    representation(fitType = "character",
        call = "call",
        t0 = "numeric",
        t.star = "numeric",
        label = "character")
        )




#' Evaluate goodness-of-fit of a fitted model. 
#'
#' @param object a fitted model of class "unmarkedFit"
#' @param nsim number of bootstrap replicates
#' @param report print fit statistic every 'report' iterations during resampling
#'
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#'
#' lengths <- linetran$Length
#'
#' ltUMF <- with(linetran, {
#' 	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
#' 		siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
#' 		tlength = lengths*1000, survey = "line", unitsIn = "m")
#' 	})
#'
#' (fm <- distsamp(~ area ~habitat, ltUMF))
#'
#' (pb <- parboot(fm))
#' plot(pb)
#'
#' @exportMethod parboot
setMethod("parboot", "unmarkedFit", function(object, nsim=10, report=2, ...)  
{
	call <- match.call(call = sys.call(-1))
	formula <- object@formula
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = FALSE)
	X <- designMats$X; V <- designMats$V; y <- designMats$y
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	yvec0 <- c(t(y))
	ests <- as.numeric(coef(object, altNames = TRUE))
	aMat <- matrix(a, M, J, byrow = TRUE)
	state0 <- predict(object, type = "state", backTran = T, na.rm = F)[,1]
	statevec0 <- rep(state0, each = J) * a
	p0 <- getP(object, na.rm = FALSE)
	pvec0 <- c(t(p0))
	expected0 <- statevec0 * pvec0
	rmse0 <- sqrt(sum((sqrt(yvec0) - sqrt(expected0))^2, na.rm = TRUE))
	cat("t.star =", rmse0, "\n")      
	rmse <- numeric(nsim)
	fits <- list()
	simdata <- umf
	simList <- simulate(object, nsim = nsim, na.rm = FALSE)
	for(i in 1:nsim) {
		y.sim <- simList[[i]]
		yvec <- c(t(y.sim))
		simdata@y <- y.sim
		fits[[i]] <- update(object, data = simdata, starts = ests, ...)
		state <- predict(fits[[i]], type = "state", backTran = T, na.rm = F)[,1]
		statevec <- rep(state0, each = J) * a
		p <- getP(fits[[i]], na.rm = FALSE)
		pvec <- c(t(p))
		expected <- statevec * pvec
		rmse[i] <- sqrt(sum((sqrt(yvec) - sqrt(expected))^2, na.rm = TRUE))
		if(nsim > report && i %in% seq(report, nsim, by=report))
			cat(paste(round(rmse[(i-(report-1)):i], 1), collapse=", "), fill=T)
		}
	out <- new("parboot", call=call, t0 = rmse0, t.star = rmse)
	return(out)
	})






#' @exportMethod show
setMethod("show", "parboot", function(object) 
{
t.star <- object@t.star
t0 <- object@t0
bias <- mean(t0 - t.star)
bias.se <- sd(t0 - t.star)
R <- length(t.star)
p.val <- sum(abs(t.star-1) > abs(t0-1)) / (1+R)
stats <- c("original" = t0, "bias" = bias, "Std. error" = bias.se, 
   "p.value" = p.val)
cat("\nCall:", deparse(object@call), fill=T)
cat("\nBootstrap Statistics:\n")
print(stats, digits=3)
cat("\nt quantiles:\n")
print(quantile(t.star, probs=c(0,2.5,25,50,75,97.5,100)/100))        
})





#' @exportMethod plot
setMethod("plot", signature(x="parboot", y="missing"), function(x, y, ...)
{
op <- par(mfrow=c(1, 2))
t.star <- x@t.star
t0 <- x@t0
t.t0 <- c(t.star, t0)
bias <- mean(t0 - t.star)
bias.se <- sd(t0 - t.star)
R <- length(t.star)
p.val <- sum(abs(t.star-1) > abs(t0-1)) / (1+R)
hist(t.star, xlim=c(min(floor(t.t0)), max(ceiling(t.t0))), 
	main=paste("P =", round(p.val, 3), "; R =", format(R)))
rug(t.star)
abline(v=t0, lty=2)
qqnorm(t.star)
qqline(t.star)
title(outer=T, ...)
par(op)
})



		