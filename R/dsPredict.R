# Model-based predictions and standard errors (and CI if alpha is specified).
# If newdata is not specified, original data is used.
# This requires the deltamethod() function from package msm.

predict.mixmod <- function(fit, type=c("state", "p", "identity"), 
	link=c("log", "logit", "identity"), newdata=NULL, notconstant=NULL, 
	alpha=NULL, par=NULL, ...) 
{
type <- match.arg(type)
link <- match.arg(link)
mixture <- fit$mixture
require(msm)
switch(type, 
state = {
	par <- coef(fit, type="state")
	xpar <- paste("x", 1:length(par), sep="")
	varcov <- vcov(fit, type="state") 
	covs <- paste("covs", 1:length(par), sep="")
	if(is.null(newdata)) 
		X <- fit$Xlam
	else {
		mf <- model.frame(fit$stateformula, newdata, na.action=na.pass)
		X <- model.matrix(fit$stateformula, mf)[,names(par)]
		}
	},
p = {
	par <- coef(fit, type="p")
	xpar <- paste("x", 1:length(par), sep="")
	varcov <- vcov(fit, type="p")
	covs <- paste("covs", 1:length(par), sep="")
	if(is.null(newdata))
		X <- fit$Xp
	else { 
		mf <- model.frame(fit$detformula, newdata, na.action=na.pass)
		X <- model.matrix(fit$detformula, mf)[,names(par)]
		}
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
class(out) <- c("mixmodPred", "data.frame")
return(out)
}












# This functions makes model-averaged predictions from a list of fitted models.

predict.mixmodList <- function(fits, newdata=NULL, notconstant=NULL, 
	type=c("state", "p"), link=c("log", "logit"))
{
type <- match.arg(type)
link <- match.arg(link)
EV <- sapply(fits, function(x) {
	tmp <- predict(x, newdata=newdata, type=type, link=link)
	E <- tmp$Est.
	V <- tmp$SE
	return(list(E, V))
	})
E <- sapply(EV[1,], as.numeric)
V <- sapply(EV[2,], as.numeric)
ic <- sapply(fits, aic.c)
deltaic <- ic - min(ic)
wts <- exp(-deltaic / 2)
wts <- wts / sum(wts)
parav <- as.numeric(E %*% wts)
seav <- rowSums((V + (E - parav)^2) %*% wts)
if(is.null(notconstant))
	notconstant <- NA
else
	notconstant <- newdata[,notconstant]
out <- data.frame(Predictor = notconstant, Est. = parav, SE = seav)
rownames(out) <- NULL
return(out)
}






