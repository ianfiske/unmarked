setClass("unmarkedFitList",
    representation(fits = "list"),
    validity = function(object) {
    	fl <- object@fits
    	d1 <- getData(fl[[1]])
    	tests <- sapply(fl, function(x) all.equal(d1, getData(x)))
    	all(tests)
    	}
    )


# constructor of unmarkedFitList objects
fitList <- function(..., fits) {
	if(missing(fits)) {
		fits <- list(...)
		isList <- sapply(fits, function(x) is.list(x))
		if(sum(isList) > 1)
			stop("Specify nested models as named objects, or use fits = 'mylist'")
		if(isList[1L]) {
			warning("If supplying a list of fits, use fits = 'mylist'")
			fits <- fits[[1L]] 	# This is allowed for back-compatability.
			}
		}	
	umfl <- new("unmarkedFitList", fits=fits)
	return(umfl)
	}
	


setMethod("summary", "unmarkedFitList", function(object) {
	fits <- object@fits
	for(i in 1:length(fits))
		summary(fits[[i]])
	})



setMethod("predict", "unmarkedFitList", function(object, type, newdata=NULL, 
	backTransform = TRUE, appendData = FALSE) {
		fitList <- object@fits
		ese <- lapply(fitList, predict, type = type, newdata = newdata, 
			backTransform = backTransform)
		E <- sapply(ese, function(x) x[,"Predicted"])
		SE <- sapply(ese, function(x) x[,"SE"])
		ic <- sapply(fitList, slot, "AIC")
		deltaic <- ic - min(ic)
		wts <- exp(-deltaic / 2)
		wts <- wts / sum(wts)
		parav <- as.numeric(E %*% wts)
		seav <- as.numeric((SE + (E - parav)^2) %*% wts) # Double check this
		out <- data.frame(Predicted = parav, SE = seav)
		if(appendData) {
			if(missing(newdata))
				newdata <- getData(object@fits[[1]])
			out <- data.frame(out, newdata)
			}
		return(out)
	})






# Condition number
cn <- function(object) {
   	ev <- eigen(hessian(object))$value
   	max(ev) / min(ev)
   	}



# R-squared index from Nagelkerke (1991)				  
nagR2 <- function(fit, nullfit)
{
	n <- sampleSize(fit)
	devI <- 2 * fit@negLogLike
	devN <- 2 * nullfit@negLogLike
	r2 <- 1 - exp((devI - devN) / n)
	r2max <- 1 - exp(-1 * devN / n)
	return(r2 / r2max)
}



setGeneric("modSel",
		def = function(object, ...) {
			standardGeneric("modSel")
			}
		)

setClass("unmarkedModSel", 
	representation(
		Estimates = "matrix", 
		SE = "matrix",
		Full = "data.frame"
		)
	)
	


# Model selection results from an unmarkedFitList
setMethod("modSel", "unmarkedFitList", 
	function(object, nullmod=NULL) 
{
	fits <- object@fits
	estList <- lapply(fits, coef, altNames=T)
	seList <- lapply(fits, function(x) sqrt(diag(vcov(x, altNames=T))))
	eNames <- sort(unique(unlist(sapply(estList, names))))
	seNames <- paste("SE", eNames, sep="")
	eseNames <- character(l <- length(c(eNames, seNames)))
	eseNames[seq(1, l, by=2)] <- eNames
	eseNames[seq(2, l, by=2)] <- seNames
	cNames <- c("formula", eseNames)
	out <- data.frame(matrix(NA, ncol=length(cNames), nrow=length(fits)))
	rownames(out) <- names(fits)
	colnames(out) <- cNames
	eMat <- seMat <- matrix(NA, length(fits), length(eNames), 
		dimnames=list(names(fits), eNames))
	out$formula <- sapply(fits, function(x) deparse(x@formula))
	for(i in 1:length(eNames)) {
		eMat[,eNames[i]] <- out[,eNames[i]] <- sapply(estList, function(x) 
			x[eNames[i]])
		seMat[,eNames[i]] <- out[,seNames[i]] <- sapply(seList, function(x) 
			x[eNames[i]])
		}
	out$Converge <- sapply(fits, function(x) x@opt$convergence)
	out$CondNum <- sapply(fits, function(x) cn(x))
	out$negLogLike <- sapply(fits, function(x) x@negLogLike)
	out$K <- sapply(fits, function(x) length(coef(x)))
	out$n <- sapply(fits, function(x) sampleSize(x))
	if(!identical(length(table(out$n)), 1L))
		warning("Models are not nested. AIC comparisons not valid")
	out$AIC <- sapply(fits, function(x) x@AIC)
	out$deltaAIC <- out$AIC - min(out$AIC)
	out$AICwt <- exp(-out$deltaAIC / 2)
	out$AICwt <- out$AICwt / sum(out$AICwt)
	out$Rsq <- NA
	if(!is.null(nullmod))
		out$Rsq <- sapply(fits, nagR2, nullmod)
	out <- out[order(out$AIC),]
	out$AICwtCum <- cumsum(out$AICwt)
	msout <- new("unmarkedModSel", Estimates = eMat, SE = seMat, Full = out)
	return(msout)
})




setMethod("show", "unmarkedModSel", 
	function(object) {
		out <- object@Full[,c("n", "K", "AIC", "deltaAIC", "AICwt", "Rsq", 
			"AICwtCum")]
		print(out, digits=5)
		})
	



