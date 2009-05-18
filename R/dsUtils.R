#' include classes.R
roxygen()


#' Vectorized integrate() accepting vector inputs for 'lower' and 'upper'
vIntegrate <- Vectorize(integrate, c("lower", "upper"))

# Detection functions for perpendicular (x) and radial (r) distances
gxhn <- function(x, sigma=1) exp(-x^2/(2 * sigma^2))
gxexp <- function(x, rate) exp(-x / rate)
gxhaz <- function(x, shape=1, scale=1)  1 - exp(-(x/shape)^-scale)
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
grexp <- function(r, rate) exp(-r / rate) * r
grhaz <- function(r, shape=1, scale=1)  (1 - exp(-(r/shape)^-scale)) * r




# Multinomial cell probabilities for line or point transects under half-normal model
cp.hn <- function(d, s, survey=c("line", "point")) 
{
survey <- match.arg(survey)
switch(survey, 
line = {
	strip.widths <- diff(d)
	f.0 <- 2 * dnorm(0, 0, sd=s)
	int <- 2 * (pnorm(d[-1], 0, sd=s) - pnorm(d[-length(d)], 0, sd=s))
	cp <- int / f.0 / strip.widths 
	},
point = {
	W <- max(d)
	int <- as.numeric(vIntegrate(grhn, d[-length(d)], d[-1], sigma=s)["value",])
	cp <- 2 / W^2 * int
	})
return(cp)
}



cp.exp <- function(d, rate, survey=c("line", "point")) 
{
survey <- match.arg(survey)
switch(survey, 
line = {
	strip.widths <- diff(d)
	f.0 <- dexp(0, rate=rate)
	int <- pexp(d[-1], rate=rate) - pexp(d[-length(d)], rate=rate)
	cp <- int / f.0 / strip.widths
	},
point = {
	W <- max(d)
	int <- as.numeric(vIntegrate(grexp, d[-length(d)], d[-1], 
		rate=rate)["value",])
	cp <- 2 / W^2 * int
	})
return(cp)
}



cp.haz <- function(d, shape, scale, survey=c("line", "point"))
{
survey <- match.arg(survey)
switch(survey, 
line = {
	strip.widths <- diff(d)
	int <- as.numeric(vIntegrate(gxhaz, d[-length(d)], d[-1], shape=shape, 
		scale=scale)["value",])
	cp <- int / strip.widths
	},
point = {
	W <- max(d)
	int <- as.numeric(vIntegrate(grhaz, d[-length(d)], d[-1], shape=shape, 
		scale=scale)["value",])
	cp <- 2 / W^2 *int
	})
return(cp)
}






ll.halfnorm <- function(param, Y, Xlam, Xp, J, d, a, nAP, nP, survey) 
{
sigma <- as.numeric(exp(Xp %*% param[(nAP+1):nP]))
lambda <- as.numeric(exp(Xlam %*% param[1:nAP]))
pvec <- c(sapply(sigma, function(x) cp.hn(d=d, s=x, survey=survey)))
growlam <- rep(lambda, each=J)
datavec <- c(t(Y))
ll <- dpois(datavec, growlam * pvec * a, log=T)
-sum(ll)
}




ll.exp <- function(param, Y, Xlam, Xp, K, J, a, d, nAP, nP, survey=survey)
{
rate <- as.numeric(exp(Xp %*% param[(nAP+1):nP]))
lambda <- as.numeric(exp(Xlam %*% param[1:nAP]))
pvec <- c(sapply(rate, function(x) cp.exp(d=d, rate=x, survey=survey)))
growlam <- rep(lambda, each=J)
datavec <- c(t(Y))
ll <- dpois(datavec, growlam * pvec * a, log=T)
-sum(ll)
}


ll.hazard <- function(param, Y, Xlam, Xp, J, a, d, nAP, nP, survey)
{
shape <- as.numeric(exp(Xp %*% param[(nAP+1):(nP-1)]))
scale <- as.numeric(exp(param[nP]))
lambda <- as.numeric(exp(Xlam %*% param[1:nAP]))
pvec <- c(sapply(shape, function(x) cp.haz(d=d, shape=x, scale=scale, 
	survey=survey)))
growlam <- rep(lambda, each=J)
datavec <- c(t(Y))
ll <- dpois(datavec, growlam * a * pvec, log=T)
-sum(ll)
}




ll.uniform <- function(param, Y, Xlam, Xp, J, a)
{
lambda <- exp(Xlam %*% param)
growlam <- rep(lambda, each=J)
datavec <- c(t(Y))
ll <- dpois(datavec, growlam * a, log=T)
-1*sum(ll)
}












# Condition number
cn <- function(fit) 
{
ev <- eigen(fit$hessian)$values
max(ev)/min(ev)
}



