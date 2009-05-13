
parboot.distsamp <- function(data, model, R=10, report=2, sp.=NULL, ...)
{
stateformula <- model$stateformula
detformula <- model$detformula
Y <- model$y
ynames <- colnames(Y)
yvec0 <- c(t(Y))
Xlam <- model$Xlam
Xp <-  model$Xp
n <- model$n
J <- model$J
ests <- coef(model, altNames=T)
tl <- model$tlength
a <- model$area
saMat <- matrix(a, n, J, byrow=T)
lam.pars0 <- coef(model, type="state")
ppars0 <- coef(model, type="p")
d <- model$dist.breaks
lam0 <- as.numeric(exp(Xlam %*% lam.pars0))
growlam0 <- rep(lam0, each=J)
key <- model$key
survey <- model$survey
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
	y.sim[k,] <- rpois(J, lam0[k] * p0[,k] * saMat[k,])
simdata <- data[rownames(Y),]
simdata[,ynames] <- y.sim
fits[[i]] <- update(model, data=simdata, starts=coef(model), ...)
lampars.i <- coef(fits[[i]], type="state")
lam <- as.vector(exp(Xlam %*% lampars.i))
growlam <- rep(lam, each=J)
ppars.i <- coef(fits[[i]], type="p")
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
out <- list(t0 = sse0, t.star = SSEs, Call = call, fits = fits, species=sp.)
class(out) <- "parboot"
return(out)
}



