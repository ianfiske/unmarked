#' @include distsamp.R
#' @include dsUtils.R
roxygen()

#setGeneric("predict",
#		def = function(object, ...) {
#			standardGeneric("predict")		})


#' @exportMethod predict
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





#' Model-averaged prediction from a list of fitted models.
#'
#' @param object a list of fitted distance sampling models all of class "umDistsampFit"
#' @param newdata an optional data.frame with new predictor variables
#' @param notconstant the column name of a variable in newdata that is not constant
#' @param type the type of prediction to make. Either "state" or "det".
#' @param link the type of link function to use
#'
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#'  
#' #Fit a few models
#' (fmhnA.AH <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m"))
#' (fmhnA <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~1, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m"))
#' (fmhn.A <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~area, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m"))
#' (fmhn <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m"))
#' (fmexp <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, key="exp",
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'  	unitsIn="m", starts=c(0,0)))
#' 
#' # Put models in a list and make it a "umDistsampFitList"
#' fl <- list(fmhn, fmexp, fmhnA.AH, fmhnA, fmhn.A)
#' fmList <- new("umDistsampFitList", fitlist=fl)
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
#'	plot(Predictor, Est., ylim=c(0.7, 1.1), xlab="Area", 
#'		ylab="Density (object / ha)")
#'	segments(Predictor, Est.-SE, Predictor, Est.+SE)
#'	})
#'
#' @exportClass umDistsampFitList
setClass("umDistsampFitList", representation(fitlist="list"))



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
out <- data.frame(Predictor=newdata[,notconstant], Est. = parav, SE = seav)
rownames(out) <- NULL
return(out)
})






