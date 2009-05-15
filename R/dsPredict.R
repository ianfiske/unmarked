# Model-based predictions and standard errors (and CI if alpha is specified).
# If newdata is not specified, original data is used.
# This requires the deltamethod() function from package msm.
setGeneric("predict",
		def = function(object, ...) {
			standardGeneric("predict")		})

setMethod("predict", "umDistsampFit", 
   function(object, type=c("state", "det"), 
   link=c("log", "logit", "identity"), newdata=NULL, notconstant=NULL, 
	alpha=NULL, par=NULL, ...) 
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
	mf <- model.frame(object@stateformula, newdata) #, na.action=na.pass)
	X <- model.matrix(object@stateformula, mf)#[,names(par)]
	},
det = {
	par <- coef(object, type="det")
	xpar <- paste("x", 1:length(par), sep="")
	varcov <- vcov(object, type="det")
	covs <- paste("covs", 1:length(par), sep="")
	mf <- model.frame(object@detformula, newdata, na.action=na.pass)
	X <- model.matrix(object@detformula, mf)#[,names(par)]
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








setGeneric("modavg",
    def = function(object, ...) {
      standardGeneric("modavg")
    })





setClass("umDistsampFitList", representation(fitlist="list"))



# This functions makes model-averaged predictions from a list of fitted models.

setMethod("modavg", "umDistsampFitList", function(object, newdata=NULL, 
   notconstant=NULL, type=c("state", "p"), link=c("log", "logit"))
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
out <- data.frame(Est. = parav, SE = seav)
rownames(out) <- NULL
return(out)
})






