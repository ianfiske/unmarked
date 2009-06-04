#' @include distsamp.R
roxygen()


#' List of fitted distance sampling models, each of class "umDistsampFit"
#'
#' Useful for model selection or model-averging predictions.
#'
#' @usage 
#' # Model-averged prediction
#' predict(object, newdata, notconstant, type, link)@param object a fitted distance sampling model of class 'umDistsampFit'
#' # Model-selection
#' modSel(object, nullmod)
#' 
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' 
#' #Fit a model
#' (fmhnA.H <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#' 	unitsIn="m"))
#'
#' # Create new data.frame for prediction
#' newhabitat <- factor("A")
#' newarea <- seq(min(linetran$area), max(linetran$area), length=20)
#' (ndA <- data.frame(area=newarea, habitat=newhabitat))
#'
#' # Predict density based upon fitted model and 'new' data
#' (Elam.A <- predict(fmhnA.H, newdata=ndA, notconstant="area", type="state", 
#' 	link="log"))
#' with(Elam.A, { # Plot relationships between density and area
#' 	plot(Predictor, Est., ylim=c(0.5, 1.2), xlab="Area", ylab="Density")
#' 	segments(Predictor, Est.-SE, Predictor, Est.+SE)
#' 	})
#' 
#' # Same as above but for detection 
#' ndH <- data.frame(area=mean(linetran$area), habitat=factor(c("a", "b")))
#' ndH
#' 
#' (Ep.H <- predict(fmhnA.H, newdata=ndH, notconstant="habitat", type="det", 
#' 	link="log")
#'
#' with(Ep.H, {  # Plot relationship (lack of) difference in detectability between habitat types
#' 	bp <- barplot(Est., ylim=c(0, 15), names=Predictor, 
#' 		ylab="Sigma (half-normal shape parameter)") 
#' 	arrows(bp, Est., bp, Est.+SE, code=2, angle=90, length=0.2)
#' 	box()
#' 	})
#' # OR
#' with(Ep.H, {
#' 	plot(function(x) unmarked:::gxhn(x, Est.[1]), 0, 25, 
#' 		xlab="Distance (m)", ylab="Detection probability")
#' 	plot(function(x) unmarked:::gxhn(x, Est.[2]), 0, 25, add=T, lty=2, lwd=2)
#' 	legend(1, 20, c("Habitat a", "Habitat b"), lty=1:2)
#' 	})
#' @exportMethod predict
#' @importFrom msm deltamethod
setMethod("predict", "umDistsampFit", 
   function(object, type=c("state", "det"), 
   link=c("log", "logit", "identity"), newdata=NULL, notconstant=NULL, 
	alpha=NULL, ...) 
{
require(msm)
type <- match.arg(type)
link <- match.arg(link)
if(is.null(newdata)) 
	newdata <- object@data
switch(type, 
state = {
	par <- coef(object, type="state")
	xpar <- paste("x", 1:length(par), sep="")
	varcov <- vcov(object, type="state") 
	covs <- paste("covs", 1:length(par), sep="")
	X <- model.matrix(object@stateformula[-2], newdata)
	},
det = {
	par <- coef(object, type="det")
	xpar <- paste("x", 1:length(par), sep="")
	varcov <- vcov(object, type="det")
	covs <- paste("covs", 1:length(par), sep="")
	X <- model.matrix(object@detformula, newdata)
	})
SEs <- rep(NA, nrow(X))
switch(link, 
	log = { 
		form <- paste("~exp(", paste(xpar, "*", covs, sep="", 
			collapse=" + "), ")", sep="")
		fitted <- exp(X %*% par)
	},
	logit = {
		form <- paste("~exp(", paste(xpar, "*", covs, sep="",
			collapse=" + "), ") / (1 + exp(", paste(xpar, "*", covs, 
			sep="", collapse=" + "), "))", sep="")
		fitted <- plogis(X %*% par)
	},
	identity = {
		form <- paste("~", paste(xpar, "*", covs, sep="", collapse=" + "),
			sep="")
		fitted <- X %*% par
	})
for(i in 1:nrow(X)) {
	for(j in 1:length(par)) {
		assign(covs[j], X[i, j], pos=.GlobalEnv)
	}
SEs[i] <- deltamethod(as.formula(form), par, varcov, ...)
}
rm(list=covs, pos=.GlobalEnv)
if(!is.null(notconstant))
	notconstant <- newdata[,notconstant]
else
	notconstant <- NA
out <- data.frame("Predictor" = notconstant, "Est." = fitted, "SE" = SEs)
if(!is.null(alpha)) {
	tval <- qt(alpha/2, Inf, lower.tail=F)
	out$lower <- out$Est. - out$SE * tval
	out$upper <- out$Est. + out$SE * tval
	attr(out, "alpha.level") <- alpha
	}
rownames(out) <- NULL
return(out)
})






#' Class for a list of distance-sampling fits each of class "umDistsampFit"
#'
#' @param fitlist a list of fitted distance sampling models all of class "umDistsampFit"
#'
#' @usage
#' # Model-averaged predictions
#' predict(object, newdata, notconstant, type, link)
#' # Model selection info
#' modSel(object, nullmod)
#'
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#'  
#' #Fit a few models
#' (fmhnA.AH <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
#'     dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'     unitsIn="m"))
#' (fmhnA <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~1, linetran, 
#'     dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'     unitsIn="m"))
#' (fmhn.A <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~area, linetran, 
#'     dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'     unitsIn="m"))
#' (fmhn <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m"))
#' (fmexp <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, key="exp",
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m", starts=c(0,0)))
#' (fmunif <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, key="uniform",
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m"))
#' 
#' # Put models in a list and make it a "umDistsampFitList"
#' fl <- list(hnNull=fmhn, expNull=fmexp, hnA.AH=fmhnA.AH, hnA=fmhnA, 
#'    hn.A=fmhn.A)
#' fmList <- umDistsampFitList(fitlist=fl)
#' 
#' # Make a new data.frame and average predictions from the list
#' # Note, you shouldn't do type="det" in this case because multiple keyfuns used
#' newarea <- seq(min(linetran$area), max(linetran$area), length=20)
#' newhabitat <- factor("a")
#' nd <- data.frame(area=newarea, habitat=newhabitat)
#' 
#' (ElamA <- predict(fmList, newdata=nd, notconstant="area", type="state"))
#'
#' with(ElamA, {
#'	plot(Predictor, Est., ylim=c(0.6, 1.2), xlab="Area", 
#'		ylab="Density (objects / ha)")
#'	segments(Predictor, Est.-SE, Predictor, Est.+SE)
#'	})
#'
#'
#' @exportClass umDistsampFitList
setClass("umDistsampFitList", 
	representation(fitlist="list"), 
    validity = function(object) {
    	fl <- object@fitList
    	d1 <- fl[[1]]@data
    	tests <- sapply(fl, function(x) all.equal(d1, x@data))
    	all(tests)
    	}
)


# Constructor of umDistsampFitList
umDistsampFitList <- function(fitList) {
	umdfl <- new("umDistsampFitList", fitList=fitList)
	return(umdfl)
	}


# TODO: add modSel function



#' @exportMethod predict
setMethod("predict", "umDistsampFitList", function(object, newdata=NULL, 
   notconstant=NULL, type=c("state", "det"), link=c("log", "logit"))
{
type <- match.arg(type)
link <- match.arg(link)
fitlist <- object@fitlist
EV <- sapply(fitlist, function(x) {
	tmp <- predict(x, newdata=newdata, type=type, link=link)
	E <- tmp$Est.
	V <- tmp$SE
	return(list(E, V))
	})
E <- sapply(EV[1,], as.numeric)
V <- sapply(EV[2,], as.numeric)
ic <- sapply(fitlist, function(x) x@AIC)
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






