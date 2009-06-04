#' @include unmarkedFit.R
roxygen()

#' Class to store nested unMarked model fits information
#' 
#' @slot fitList a (preferably named) list of nested model fits.
#'
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
#' @export
unmarkedFitList <- function(fitList) {
	umfl <- new("unmarkedFitList", fitList=fitList)
	return(umfl)
	}


#' @exportMethod predict
setMethod("predict", "unmarkedFitList", function(object, type, newdata=NULL, 
	backTran=TRUE)
{
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









