#' @include distsamp.R
#' @include umDistsampFit.R

roxygen()


#' Class to store nested umDistsampFit model fits information
#' 
#' @param fitList a (preferably named) list of nested model fits.
#'
#' @exportClass umDistsampFitList
setClass("umDistsampFitList",
    representation(fitList = "list"),
    validity = function(object) {
    	fl <- object@fitList
    	d1 <- fl[[1]]@data
    	tests <- sapply(fl, function(x) all.equal(d1, x@data))
    	all(tests)
    	}
    )

#' @export
dsFitList <- function(fitList) {
	umdsfl <- new("umDistsampFitList", fitList=fitList)
	return(umdsfl)
	}



cn <- function(object) {
   ev <- eigen(object@hessian)$value
   max(ev) / min(ev)
   }




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
setClass("umModSel", 
	representation(
		Estimates = "matrix", 
		SE = "matrix",
		Full = "data.frame"
		)
	)
	


#' Model selection results from an umDistsampFitList
#'
#' @name modSel-umDistsampFitList
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
#' fl <- dsFitList(fitList = list(Null=fm1, A.=fm2, .A=fm2))
#' fl
#'
#' (ms <- modSel(fl, nullmod=fm1))
#'
#' ms@Full
#' @exportMethod modSel
setMethod("modSel", "umDistsampFitList", function(object, nullmod=NULL) 
{
fits <- object@fitList
estList <- lapply(fits, coef, altNames=T)
seList <- lapply(fits, function(x) sqrt(diag(vcov(x, altNames=T))))
eNames <- sort(unique(unlist(sapply(estList, names))))
seNames <- paste("SE", eNames, sep="")
eseNames <- character(l <- length(c(eNames, seNames)))
eseNames[seq(1, l, by=2)] <- eNames
eseNames[seq(2, l, by=2)] <- seNames
cNames <- c("stateFormula", "detFormula", eseNames)
out <- data.frame(matrix(NA, ncol=length(cNames), nrow=length(fits)))
eMat <- seMat <- matrix(NA, length(fits), length(eNames), 
   dimnames=list(names(fits), eNames))
rownames(out) <- names(fits)
colnames(out) <- cNames
out$stateFormula <- sapply(fits, function(x) deparse(x@stateFormula))
out$detFormula <- sapply(fits, function(x) deparse(x@detFormula))
for(i in 1:length(eNames)) {
	eMat[,eNames[i]] <- out[,eNames[i]] <- sapply(estList, function(x) x[eNames[i]])
	seMat[,eNames[i]] <- out[,seNames[i]] <- sapply(seList, function(x) x[eNames[i]])
	}
out$Converge <- sapply(fits, function(x) x@optout$convergence)
out$CondNum <- sapply(fits, cn)
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
msout <- new("umModSel", Estimates = eMat, SE = seMat, Full = out)
return(msout)
})





#' @exportMethod show
setMethod("show", "umModSel", 
	function(object) {
		out <- object@Full[,c("n", "K", "AIC", "deltaAIC", "AICwt", "Rsq", 
			"AICwtCum")]
		print(out, digits=5)
		})
	



#' @name predict-umDistsampFitList
#' @aliases predict,umDistsampFitList-method
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
#' fl <- dsFitList(fitList = list(Null=fm1, A.=fm2, .A=fm2))
#' fl
#'
#' predict(fl)
#' @exportMethod predict
setMethod("predict", "umDistsampFitList", function(object, newdata=NULL, 
   notconstant=NULL, type=c("state", "det"), link=c("log", "logit"))
{
type <- match.arg(type)
link <- match.arg(link)
fitList <- object@fitList
EV <- sapply(fitList, function(x) {
	tmp <- predict(x, newdata=newdata, type=type, link=link)
	E <- tmp$Est.
	V <- tmp$SE
	return(list(E, V))
	})
E <- sapply(EV[1,], as.numeric)
V <- sapply(EV[2,], as.numeric)
ic <- sapply(fitList, function(x) x@AIC)
deltaic <- ic - min(ic)
wts <- exp(-deltaic / 2)
wts <- wts / sum(wts)
parav <- as.numeric(E %*% wts)
seav <- rowSums((V + (E - parav)^2) %*% wts)
if(!is.null(notconstant))
	notconstant <- newdata[,notconstant]
else
	notconstant <- NA
out <- data.frame(Predictor=notconstant, Est. = parav, SE = seav)
rownames(out) <- NULL
return(out)
})


