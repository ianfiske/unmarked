#' @include unmarkedFit.R
roxygen()

#' @export
setClass("unmarkedFitList",
    representation(fitList = "list"),
    validity = function(object) {
    	fl <- object@fitList
    	d1 <- fl[[1]]@data
    	tests <- sapply(fl, function(x) all.equal(d1, x@data))
    	all(tests)
    	}
    )


#' constructor of unmarkedFitList objects
#'
#' @examples
#' # Fit some N-mixture models
#' data(mallard)
#' mallardUMF <- unmarkedFrame(mallard.y, siteCovs = mallard.site,
#' 	obsCovs = mallard.obs)
#' 
#' fm1 <- pcount(~ length, ~ ivel, mallardUMF)
#' fm2 <- pcount(~ length, ~ 1, mallardUMF)
#' fm3 <- pcount(~ 1, ~ 1, mallardUMF)
#' 
#' # Create an unmarkedFitList with a named list of models
#' fmList <- unmarkedFitList(fitList=list(Global=fm1, Length.=fm2, Null=fm3))
#' fmList
#' 
#' # Model-averaged prediction
#' predict(fmList, type="state")
#' 
#' # Model selection
#' modSel(fmList, nullmod=fm3)
#'
#' @export
unmarkedFitList <- function(fitList) {
	umfl <- new("unmarkedFitList", fitList=fitList)
	return(umfl)
}

##TODO Fix predict-unmarkedFit so that it won't fail when type="det" for 2 parameter detection functions (eg distsamp(keyfun="hazard"))
#' @exportMethod predict
setMethod("predict", "unmarkedFitList", function(object, type, newdata=NULL, 
				backTran=TRUE) {
			fitList <- object@fitList
			ese <- lapply(fitList, predict, type=type, newdata=newdata, backTran=backTran)
			E <- sapply(ese, function(x) x[,"Predicted"])
			SE <- sapply(ese, function(x) x[,"SE"])
			ic <- sapply(fitList, slot, "AIC")
			deltaic <- ic - min(ic)
			wts <- exp(-deltaic / 2)
			wts <- wts / sum(wts)
			parav <- as.numeric(E %*% wts)
			seav <- rowSums((SE + (E - parav)^2) %*% wts)
			out <- data.frame(Predicted = parav, SE = seav)
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
n <- nrow(fit@data@y)
devI <- 2 * fit@negLogLike
devN <- 2 * nullfit@negLogLike
r2 <- 1 - exp((devI - devN) / n)
r2max <- 1 - exp(-1 * devN / n)
return(r2 / r2max)
}



#' @exportMethod modSel
setGeneric("modSel",
		def = function(object, ...) {
			standardGeneric("modSel")
			}
		)

#' @export
setClass("unmarkedModSel", 
	representation(
		Estimates = "matrix", 
		SE = "matrix",
		Full = "data.frame"
		)
	)
	


#' Model selection results from an unmarkedFitList
#'
#' @name modSel-unmarkedFitList
#' @aliases modSel modSel-methods
#' @param object an object of class "umDistsampFitList"
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' 
#' ## Half-normal detection function. Density output. No covariates. 
#' ## lineDat$Length is transect lengths in km, so it has to be converted.
#' (fm1 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, survey="line", unitsIn="m"))
#'
#' (fm2 <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, survey="line", unitsIn="m"))
#'
#' (fm3 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~area, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, survey="line", unitsIn="m"))
#'
#' fl <- unmarkedFitList(fitList = list(Null=fm1, A.=fm2, .A=fm2))
#' fl
#'
#' (ms <- modSel(fl, nullmod=fm1))
#'
#' ms@@Full
#' @exportMethod modSel
setMethod("modSel", "unmarkedFitList", function(object, nullmod=NULL) 
{
fits <- object@fitList
estList <- lapply(fits, coef, altNames=T)
seList <- lapply(fits, function(x) sqrt(diag(vcov(x, altNames=T))))
eNames <- sort(unique(unlist(sapply(estList, names))))
seNames <- paste("SE", eNames, sep="")
eseNames <- character(l <- length(c(eNames, seNames)))
eseNames[seq(1, l, by=2)] <- eNames
eseNames[seq(2, l, by=2)] <- seNames
cNames <- c("stateformula", "detformula", eseNames)
out <- data.frame(matrix(NA, ncol=length(cNames), nrow=length(fits)))
rownames(out) <- names(fits)
colnames(out) <- cNames
eMat <- seMat <- matrix(NA, length(fits), length(eNames), 
   dimnames=list(names(fits), eNames))
out$stateformula <- sapply(fits, function(x) deparse(x@stateformula))
out$detformula <- sapply(fits, function(x) deparse(x@detformula))
for(i in 1:length(eNames)) {
	eMat[,eNames[i]] <- out[,eNames[i]] <- sapply(estList, function(x) x[eNames[i]])
	seMat[,eNames[i]] <- out[,seNames[i]] <- sapply(seList, function(x) x[eNames[i]])
	}
#out$Converge <- sapply(fits, function(x) x@optout$convergence)
#out$CondNum <- sapply(fits, cn)
out$negLogLike <- sapply(fits, function(x) x@negLogLike)
out$K <- sapply(fits, function(x) length(coef(x)))
out$n <- sapply(fits, function(x) nrow(x@data@y))
out$AIC <- sapply(fits, function(x) x@AIC)
out$deltaAIC <- out$AIC-min(out$AIC)
out$AICwt <- exp(-out$deltaAIC/2)
out$AICwt <- out$AICwt/sum(out$AICwt)
out$Rsq <- NA
if(!is.null(nullmod))
	out$Rsq <- sapply(fits, nagR2, nullmod)
out <- out[order(out$AIC),]
out$AICwtCum <- cumsum(out$AICwt)
msout <- new("unmarkedModSel", Estimates = eMat, SE = seMat, Full = out)
return(msout)
})





#' @exportMethod show
setMethod("show", "unmarkedModSel", 
	function(object) {
		out <- object@Full[,c("n", "K", "AIC", "deltaAIC", "AICwt", "Rsq", 
			"AICwtCum")]
		print(out, digits=5)
		})
	



