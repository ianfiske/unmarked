

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










print.parBoot <- function(boot.out)
{
t.star <- boot.out$t.star
t0 <- boot.out$t0
bias <- mean(t0 - t.star)
bias.se <- sd(t0 - t.star)
R <- length(t.star)
t.val <- abs(bias / bias.se)
df <- as.integer(R - 1)
#p.val <- 1 - pt(t.val, df = df)
p.val <- sum(abs(t.star-1) > abs(t0-1)) / (1+R)
stats <- c("original" = t0, "bias" = bias, "Std. error" = bias.se, 
	"t" = t.val, "p.value" = p.val)
cat("\nCall:", deparse(boot.out$call), fill=T)
cat("\nBootstrap Statistics:\n")
print(stats, digits=3)
cat("\nt quantiles:\n")
print(quantile(t.star, probs=c(0,2.5,25,50,75,97.5,100)/100))
}






plot.parBoot <- function(boot.out, ...)
{
op <- par(mfrow=c(1, 2))
t.star <- boot.out$t.star
t0 <- boot.out$t0
t.t0 <- c(t.star, t0)
bias <- mean(t0 - t.star)
bias.se <- sd(t0 - t.star)
R <- length(t.star)
t.val <- abs(bias / bias.se)
df <- R - 1
p <- 1 - pt(t.val, df = df)
# p <- sum(abs(t.star-1) > abs(t0-1)) / (1+R)
hist(t.star, xlim=c(min(floor(t.t0)), max(ceiling(t.t0))), 
	main=paste("P =", format(p), "; R =", format(R)))
rug(t.star)
abline(v=t0, lty=2)
qqnorm(t.star)
qqline(t.star)
title(outer=T, ...)
par(op)
}










