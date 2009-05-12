

print.mixmod <- function(model, nulldev=NULL, pval=T, digits=3, ...) 
{
if(!exists("model")) stop("This model doesn't exist")
est <- model$par
lamInd <- which(attr(est, "type") == "state")
if(identical(model$mixture, "NB")) 
	lamInd <- c(lamInd, which(names(est)=="alpha"))
pInd <- which(attr(est, "type") == "p")
if(class(model$SE) != "try-error") {
	se <- model$SE
	if(is.null(se)) {
		se <- NA
		}
	outL <- data.frame(Est. = est[lamInd], SE = se[lamInd])
	outP <- data.frame(Est. = est[pInd], SE = se[pInd])
	}
else {
	outL <- data.frame("Est." = est[lamInd], "SE" = NA)
	outP <- data.frame("Est." = est[pInd], "SE" = NA)
	}
if(pval) {
	outL$'z value' <- zL <- outL$Est. / outL$SE
	outL$'Pr(>|z|)' <- 2 * pnorm(-abs(zL))
	outP$'z value' <- zP <- outP$Est. / outP$SE
	outP$'Pr(>|z|)' <- 2 * pnorm(-abs(zP))
	}
cat("\n", "Call: ", deparse(model$call), "\n", sep="", fill=T)
cat("State coefficients:\n")
print(outL, digits=digits, ...)
cat("\nDetection coefficients:\n")
print(outP, digits=digits, ...)
cat("\nSample size:", model$n)
if(!is.null(nulldev)) cat("\nNull Deviance:", format(nulldev))
cat("\nResidual Deviance:", format(2 * model$value), "\n")
cat("AICc:", aic.c(model), "\n")
cat("\noptim convergence code:", model$convergence)
cat("\noptim iterations:", model$counts[1], "\n", "\n")
}







coef.mixmod <- function(model, type=c("all", "state", "p", 
	"dispersion"), altNames=F)
{
type <- match.arg(type)
ests <- model$par
if(altNames)
	names(ests) <- attr(ests, "altNames")
switch(type, 
	all = co <- ests,
	state = co <- ests[attr(ests, "type") == "state"],
	p = co <- ests[attr(ests, "type") == "p"],
	dispersion = co <- ests["alpha"])
attr(co, "type") <- attr(co, "altNames") <- NULL
return(co)
}




vcov.mixmod <- function(model, type=c("all", "state", "p", 
	"dispersion"), altNames=F)
{
type <- match.arg(type)
ests <- model$par
hess <- model$hessian
if(altNames)
	dimnames(hess) <- list(attr(ests, "altNames"), attr(ests, "altNames"))
switch(type, 
	all = vc <- hess,
	state = vc <- hess[attr(ests, "type") == "state", 
		attr(ests, "type") == "state"],
	p = vc <- hess[attr(ests, "type") == "p", attr(ests, "type") == "p"],
	dispersion = vc <- hess["alpha", "alpha"])
vc <- solve(vc)
return(vc)
}

