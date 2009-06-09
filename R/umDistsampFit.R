#' @include unmarkedFit.R
roxygen()


#' @exportClass umDistsampFit
setClass("umDistsampFit",
		representation(
				optout = "list",
				keyfun = "character",
				dist.breaks = "numeric",
				tlength = "numeric",
				area = "numeric",
				survey = "character",
				unitsIn = "character",
				unitsOut = "character"),
		contains = "unmarkedFit")







#' Model-based prediction 
#'
#' Predict expectations of a model with standard errors (and CI if alpha is specified).
#' If newdata is not specified, original data is used.
#' Requires the deltamethod() function from package msm.
#'
#' @name predict-umDistsampFit
#' @aliases predict umDistsampFit-method
#' @param object a fitted distance sampling model of class 'umDistsampFit'
#' @param type the component of the model to predict. Either "state" or "det"
#' @param link the link function to be used
#' @param newdata an optional data.frame including variables necessary to predict from fitted model
#' @param notconstant an optional column name of a variable in newdata which is not constant. Only useful for plotting expected relationships.
#' @param alpha alpha level for creating confidence intervals
#' 
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' 
#' #Fit a model
#' (fmhnA.H <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
#' 		dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#' 		unitsIn="m"))
#'
#' # Create new data.frame for prediction
#' newhabitat <- factor("A")
#' newarea <- seq(min(linetran$area), max(linetran$area), length=20)
#' (ndA <- data.frame(area=newarea, habitat=newhabitat))
#'
#' # Predict density based upon fitted model and 'new' data
#' (Elam.A <- predict(fmhnA.H, newdata=ndA, notconstant="area", type="state", 
#' 		link="log"))
#'
#' with(Elam.A, { # Plot relationships between density and area
#' 		plot(Predictor, Est., ylim=c(0.5, 1.2), xlab="Area", ylab="Density")
#' 		segments(Predictor, Est.-SE, Predictor, Est.+SE)
#' 		})
#' 
#' # Same as above but for detection 
#' ndH <- data.frame(area=mean(linetran$area), habitat=factor(c("a", "b")))
#' ndH
#' 
#' (Ep.H <- predict(fmhnA.H, newdata=ndH, notconstant="habitat", type="det", 
#' 		link="log"))
#'
#' with(Ep.H, {  # Plot (lack of) difference in detectability between habitat types
#' 		bp <- barplot(Est., ylim=c(0, 15), names=Predictor, 
#' 		ylab="Sigma (half-normal shape parameter)") 
#' 	arrows(bp, Est., bp, Est.+SE, code=2, angle=90, length=0.2)
#' 	box()
#' 	})
#' # OR
#' with(Ep.H, {
#' 		plot(function(x) unmarked:::gxhn(x, Est.[1]), 0, 25, 
#' 		xlab="Distance (m)", ylab="Detection probability")
#' 		plot(function(x) unmarked:::gxhn(x, Est.[2]), 0, 25, add=TRUE, lty=2, lwd=2)
#' 		legend(1, 20, c("Habitat a", "Habitat b"), lty=1:2)
#' 		})
#'
#' @exportMethod predict
#' @importFrom msm deltamethod
setMethod("predict", "umDistsampFit", 
   function(object, type=c("state", "det"), 
   link=c("log", "logit", "identity"), newdata=NULL, notconstant=NULL, 
	alpha=NULL) 
{
require(msm)
type <- match.arg(type)
link <- match.arg(link)
if(is.null(newdata)) {
	newdata <- object@data
	newdata <- data.frame(newdata@y, newdata@siteCovs)
	}
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
SEs[i] <- deltamethod(as.formula(form), par, varcov)
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













#' @exportMethod parboot
setGeneric("parboot",
    def = function(object, ...) {
      standardGeneric("parboot")
    })




#' @exportClass umParboot
setClass("umParboot",
    representation(fitType = "character",
        call = "call",
        t0 = "numeric",
        t.star = "numeric",
        label = "character")
        )


#' @exportMethod update
setMethod("update", "umDistsampFit", 
	function(object, formula., ..., evaluate = TRUE) {
		call <- object@call
		if (is.null(call)) 
            stop("need an object with call slot")
        extras <- match.call(expand.dots = FALSE)$...
        if (!missing(formula.)) 
            call$formula <- update.formula(formula(object), formula.)
        if (length(extras) > 0) {
            existing <- !is.na(match(names(extras), names(call)))
            for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
            if (any(!existing)) {
                call <- c(as.list(call), extras[!existing])
                call <- as.call(call)
            	}
        	}
        if (evaluate) 
            eval(call, parent.frame())
        else call
		}
	)


#' Evaluate goodness-of-fit for a fitted distance-sampling model
#'
#' @param object a fitted model of class "umDistsampFit"
#' @param R number of bootstrap replicates
#' @param report print fit statistic every 'report' iterations during resampling
#' @param label a label for the model
#'
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' (fm <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~habitat, linetran, 
#'	   dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#'	   unitsIn="m"))
#'
#' (pb <- parboot(fm))
#' plot(pb)
#'
#' @exportMethod parboot
setMethod("parboot", "umDistsampFit", function(object, R=10, report=2, 
   label=character(0), ...)  
{
call <- match.call()
stateformula <- object@stateformula
detformula <- object@detformula
data <- object@data
data <- data.frame(data@y, data@siteCovs)
mf <- model.frame(stateformula, data)
Y <- model.response(mf)
ynames <- colnames(Y)
yvec0 <- c(t(Y))
Xlam <- model.matrix(stateformula, data)
Xp <-  model.matrix(detformula, data)
n <- nrow(Y)
J <- ncol(Y)
ests <- coef(object, altNames=T) 
tl <- object@tlength
a <- object@area
aMat <- matrix(a, n, J, byrow=T)
lam.pars0 <- coef(object, type="state")
ppars0 <- coef(object, type="det")
d <- object@dist.breaks
lam0 <- as.numeric(exp(Xlam %*% lam.pars0))
growlam0 <- rep(lam0, each=J)
key <- object@keyfun
survey <- object@survey
switch(key, 
"halfnorm" = {
	sigma <- exp(Xp %*% ppars0)
	p0 <- sapply(sigma, function(x) cp.hn(d=d, s=x, survey=survey))
	},
"exp" = {
	rate <- exp(Xp %*% ppars0)
	p0 <- sapply(rate, function(x) cp.exp(d=d, rate=x, survey=survey))
	},
"hazard" = {
	shape <- exp(Xp %*% ppars0[-length(ppars0)])
	scale <- exp(ppars0[length(ppars0)])
	p0 <- sapply(shape, function(x) cp.haz(d=d, shape=shape, scale=scale, 
		survey=survey))
	})
p0vec <- c(p0)
expected0 <- growlam0 * a * p0vec
sse0 <- sqrt(sum((sqrt(yvec0) - sqrt(expected0))^2, na.rm=T))
cat("t.star =", sse0, "\n")
SSEs <- numeric(R)
fits <- list()
for(i in 1:R) 
{
y.sim <- matrix(NA, n, J)
for(k in 1:n)
	y.sim[k,] <- rpois(J, lam0[k] * p0[,k] * aMat[k,])
simdata <- data[rownames(Y),]
simdata[,ynames] <- y.sim
fits[[i]] <- update(object, data=simdata, starts=as.numeric(coef(object)), ...)
lampars.i <- coef(fits[[i]], type="state")
lam <- as.vector(exp(Xlam %*% lampars.i))
growlam <- rep(lam, each=J)
ppars.i <- coef(fits[[i]], type="det")
switch(key, 
"halfnorm" = {
	sigma <- exp(Xp %*% ppars.i)
	p <- sapply(sigma, function(x) cp.hn(d=d, s=x, survey=survey))
	},
"exp" = {
	rate <- exp(Xp %*% ppars.i)
	p <- sapply(rate, function(x) cp.exp(d=d, r=x, survey=survey))
	},
"hazard" = {
	shape <- exp(Xp %*% ppars.i[-length(ppars.i)])
	scale <- exp(ppars.i[length(ppars.i)])
	p <- sapply(sigma, function(x) cp.haz(d=d, shape=shape, scale=scale, 
		survey=survey))
	})
pvec <- c(p)
expected <- growlam * a * pvec
yvec <- c(t(simdata[,ynames]))
sqrt.sim <- sqrt(yvec)
sqrt.model <- sqrt(expected)
SSEs[i] <- sqrt(sum((sqrt.sim - sqrt.model)^2, na.rm=T))
if(R > report && i %in% seq(report, R, by=report))
	cat(paste(round(SSEs[(i-(report-1)):i], 1), collapse=", "), fill=T)
}
out <- new("umParboot", call=call, t0 = sse0, t.star = SSEs, label=label)
return(out)
})













#' @exportMethod show
setMethod("show", "umParboot", function(object) 
{
t.star <- object@t.star
t0 <- object@t0
bias <- mean(t0 - t.star)
bias.se <- sd(t0 - t.star)
R <- length(t.star)
p.val <- sum(abs(t.star-1) > abs(t0-1)) / (1+R)
stats <- c("original" = t0, "bias" = bias, "Std. error" = bias.se, 
   "p.value" = p.val)
cat("\nCall:", deparse(object@call), fill=T)
cat("\nBootstrap Statistics:\n")
print(stats, digits=3)
cat("\nt quantiles:\n")
print(quantile(t.star, probs=c(0,2.5,25,50,75,97.5,100)/100))        
})





#' @exportMethod plot
setMethod("plot", signature(x="umParboot", y="missing"), function(x, y, ...)
{
op <- par(mfrow=c(1, 2))
t.star <- x@t.star
t0 <- x@t0
t.t0 <- c(t.star, t0)
bias <- mean(t0 - t.star)
bias.se <- sd(t0 - t.star)
R <- length(t.star)
p.val <- sum(abs(t.star-1) > abs(t0-1)) / (1+R)
hist(t.star, xlim=c(min(floor(t.t0)), max(ceiling(t.t0))), 
	main=paste("P =", round(p.val, 3), "; R =", format(R)))
rug(t.star)
abline(v=t0, lty=2)
qqnorm(t.star)
qqline(t.star)
title(outer=T, ...)
par(op)
})


  







