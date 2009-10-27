setClass("unmarkedFit",
    representation(fitType = "character",
        call = "call",
				formula = "formula",
        data = "unmarkedFrame",
				sitesRemoved = "numeric",  # vector of indices of sites that were removed from analysis
        estimates = "unmarkedEstimateList",
        AIC = "numeric",
				opt = "list",
        negLogLike = "numeric",
				nllFun = "function"))

# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call, formula,
		data, sitesRemoved, estimates, AIC, opt, negLogLike, nllFun) {
	umfit <- new("unmarkedFit", fitType = fitType,
			call = call, formula = formula, data = data, sitesRemoved = sitesRemoved,
			estimates = estimates, AIC = AIC,
			opt = opt, negLogLike = negLogLike, nllFun = nllFun)
	
	return(umfit)
}

################################# CHILD CLASSES ###############################


setClass("unmarkedFitDS",
		representation(
				keyfun = "character",
				unitsOut = "character"),
		contains = "unmarkedFit")


setClass("unmarkedFitPCount", 
		representation(
				K = "numeric",
				mixture = "character"),
		contains = "unmarkedFit")


setClass("unmarkedFitOccu", 
		representation(knownOcc = "logical"),
		contains = "unmarkedFit")			


setClass("unmarkedFitMPois", 
		contains = "unmarkedFit")			


setClass("unmarkedFitOccuRN", 
		contains = "unmarkedFit")			


setClass("unmarkedFitColExt",
		representation(phi = "matrix",
				projected = "matrix"),
		contains = "unmarkedFit")

################################################################################



setMethod("show", "unmarkedFit",
    function(object) {
      cat("\nCall:\n")
      print(object@call)
      cat("\n")
      show(object@estimates)
      cat("\nAIC:", object@AIC,"\n")
      if(!identical(object@opt$convergence, 0L))
    		warning("Model did not converge. Try providing starting values or
    			increasing maxit control argment.")
    })
    
    
    
    
setMethod("summary", "unmarkedFit", 
	function(object) 
{
      	cat("\nCall:\n")
		print(object@call)
		cat("\n")
        summary(object@estimates)      	
		cat("Sample size:", sampleSize(object))
		if(length(object@sitesRemoved) > 0)
			cat("\nSites removed:", object@sitesRemoved)
		cat("\noptim convergence code:", object@opt$convergence)
		cat("\noptim iterations:", object@opt$counts[1], "\n", "\n")
		if(!identical(object@opt$convergence, 0L))
    		warning("Model did not converge. Try providing starting values or
    			increasing maxit control argment.")
})



# Compute linear combinations of estimates in unmarkedFit objects.

setMethod("linearComb",
    signature(obj = "unmarkedFit", coefficients = "matrixOrVector"),
    function(obj, coefficients, type) {
      stopifnot(!missing(type))
      stopifnot(type %in% names(obj))
      estimate <- obj@estimates[type]
      linearComb(estimate, coefficients)
    })

setMethod("backTransform", "unmarkedFit",
		function(obj, whichEstimate) {
			est <- obj[whichEstimate]
			if(length(est@estimates) == 1) {
				lc <- linearComb(est, 1)
				return(backTransform(lc))
			} else {
				stop('Cannot directly backTransform an unmarkedEstimate with length > 1.')
			}
		})

setMethod("[",
    "unmarkedFit",
    function(x, i, j, drop) {
      x@estimates[i]
    })

setMethod("names", "unmarkedFit",
    function(x) {
      names(x@estimates)
    })
    

# Prediction
setMethod("predict", "unmarkedFit", 
	function(object, type, newdata=NULL, backTran=TRUE, na.rm = TRUE, ...) 
	{
		if(is.null(newdata))
			newdata <- object@data
		formula <- object@formula
		detformula <- as.formula(formula[[2]])
		stateformula <- as.formula(paste("~",formula[3],sep=""))
		if(inherits(newdata, "unmarkedFrame"))
			class(newdata) <- "unmarkedFrame"
		cls <- class(newdata)
		switch(cls, 
			unmarkedFrame = {
				designMats <- getDesign2(formula, newdata, na.rm = na.rm)
					switch(type, 
						state = X <- designMats$X,
						det = X <- designMats$V)
				},
			data.frame = {
				switch(type, 
					state = {
						Terms <- delete.response(terms(stateformula))
						mf <- model.frame(Terms, newdata)
						X <- model.matrix(Terms, mf)
						},
					det = X <- model.matrix(detformula, newdata))
				})
		out <- matrix(NA, nrow(X), 2, dimnames=list(NULL, c("Predicted", "SE")))
		lc <- linearComb(object, X, type)
		if(backTran) lc <- backTransform(lc)
		out[,1] <- coef(lc)
		out[,2] <- SE(lc)
		return(out)
	}
)

setMethod("coef", "unmarkedFit",
		function(object, type, altNames = TRUE) {
			if(missing(type)) {
				co <- lapply(object@estimates@estimates, 
					function(x) coef(x, altNames=altNames))
				names(co) <- NULL
				co <- unlist(co)
			} else {
				co <- coef(object[type], altNames=altNames)
			}
			co
		})



setMethod("vcov", "unmarkedFit",
		function(object, type, altNames = TRUE, method = "hessian") {
			if(missing(type)) {
				switch(method,
						hessian = v <- solve(hessian(object)),
						nonparboot = {
							fmList <- nonparboot(object)
							coefs <- t(sapply(fmList, function(x) coef(x)))
							v <- cov(coefs)
						}
				)
				rownames(v) <- colnames(v) <- names(coef(object, altNames=altNames))
			} else {
				v <- vcov(object[type])
				rownames(v) <- colnames(v) <- names(coef(object, type, altNames=altNames))
			}
			v
		})

setMethod("SE", "unmarkedFit", 
		function(obj) {
			sqrt(diag(vcov(obj)))
		})



setMethod("confint", "unmarkedFit",
	function(object, parm, level = 0.95, type, method = c("normal", "profile")) {
		method <- match.arg(method)
		if(missing(type)) stop(paste("Must specify type as one of (", paste(names(object),collapse=", "),").",sep=""))
		if(missing(parm)) parm <- 1:length(object[type]@estimates)
		if(method == "normal") {
			callGeneric(object[type],parm = parm, level = level)
		} else {
			nllFun <- nllFun(object)
			ests <- mle(object)
			nP <- length(parm)
			ci <- matrix(NA, nP, 2)
    
			## create table to match parm numbers with est/parm numbers
			types <- names(object)
			numbertable <- data.frame(type = character(0), num = numeric(0))
			for(i in seq(length=length(types))) {
				length.est <- length(object[i]@estimates)
				numbertable <- rbind(numbertable, data.frame(type = 
					rep(types[i], length.est), num = seq(length=length.est)))
			}
			parm.fullnums <- which(numbertable$type == type & 
				numbertable$num %in% parm)
				
			for(i in seq(length=nP)) {
				cat("Profiling parameter",i,"of",nP,"...")
				se <- SE(object[type])
				whichPar <- parm.fullnums[i]
				ci[i,] <- profileCI(nllFun, whichPar=whichPar, MLE=ests, 
					interval=ests[whichPar] + 10*se[i]*c(-1,1), level=level)
				cat(" done.\n")
				}
			rownames(ci) <- names(coef(object[type]))[parm]
			colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
			return(ci)
			}
	})



setMethod("fitted", "unmarkedFit",
		function(object, na.rm = FALSE) {
			data <- object@data
			des <- getDesign2(object@formula, data, na.rm = na.rm)
			X <- des$X; V <- des$V; a <- des$plotArea
			state <- do.call(object['state']@invlink, 
				list(X %*% coef(object, 'state')))
			state <- as.numeric(state) * a  ## E(X) for most models
			p <- getP(object, na.rm = na.rm) # P(detection | presence)
			fitted <- state * p  # true for models with E[Y] = p * E[X]
			fitted
		})		



setMethod("fitted", "unmarkedFitOccu",
		function(object, na.rm = FALSE) {
			data <- object@data
			des <- getDesign2(object@formula, data, na.rm = na.rm)
			X <- des$X; V <- des$V; a <- des$plotArea
			state <- plogis(X %*% coef(object, 'state'))
			state <- as.numeric(state)  ## E(X) for most models
			state[object@knownOcc] <- 1
			p <- getP(object, na.rm = na.rm) # P(detection | presence)
			fitted <- state * p  # true for models with E[Y] = p * E[X]
			fitted
		})



setMethod("fitted", "unmarkedFitPCount",
		function(object, K, na.rm = FALSE) {
			data <- object@data
			des <- getDesign2(object@formula, data, na.rm = na.rm)
			X <- des$X; V <- des$V; a <- des$plotArea
			y <- des$y	# getY(data) ... to be consistent w/NA handling?
			M <- nrow(X)
			J <- ncol(y)
			state <- exp(X %*% coef(object, 'state')) * a
			p <- getP(object, na.rm = na.rm)
			mix <- object@mixture
			switch(mix,
					P = {
						fitted <- as.numeric(state) * p 
					},
					NB = {
						if(missing(K)) K <- max(y, na.rm = TRUE) + 20 
						k <- 0:K
						k.ijk <- rep(k, M*J)
						state.ijk <- state[rep(1:M, each = J*(K+1))]
						alpha <- exp(coef(object['alpha']))
						prob.ijk <- dnbinom(k.ijk, mu = state.ijk, size = alpha)
						all <- cbind(rep(as.vector(t(p)), each = K + 1), k.ijk, prob.ijk)
						prod.ijk <- rowProds(all)
						fitted <- colSums(matrix(prod.ijk, K + 1, M*J))
						fitted <- matrix(fitted, M, J, byrow = TRUE)
					})
			fitted
		})



setMethod("fitted", "unmarkedFitOccuRN", 
		function(object, K, na.rm = FALSE) {
			data <- object@data
			des <- getDesign2(object@formula, data, na.rm = na.rm)
			X <- des$X; V <- des$V; a <- des$plotArea
			y <- des$y	# getY(data) ... to be consistent w/NA handling?
			y <- truncateToBinary(y)
			M <- nrow(X)
			J <- ncol(y)
			lam <- exp(X %*% coef(object, 'state')) * a
			r <- plogis(V %*% coef(object, 'det'))
			if(missing(K)) K <- max(y, na.rm = TRUE) + 20 
			
			lam <- rep(lam, each = J)
			
			fitted <- 1 - exp(-lam*r) ## analytical integration.
			
			return(matrix(fitted, M, J, byrow = TRUE))
		})


setMethod("fitted", "unmarkedFitColExt", 
		function(object, na.rm = FALSE) {
			data <- object@data
			M <- numSites(data)
			nY <- data@numPrimary
			J <- obsNum(data)/nY
			psi <- plogis(coef(object, 'psi'))
			detParms <- coef(object, 'det')
			colParms <- coef(object, 'col')
			extParms <- coef(object, 'ext')
			designMats <- unmarked:::getDesign3(formula = object@formula, 
				object@data)
			V.itj <- designMats$V; X.it <- designMats$X
			
			detP <- plogis(V.itj %*% detParms)
			colP <- plogis(X.it  %*% colParms)
			extP <- plogis(X.it %*% extParms)
			
			detP <- array(detP, c(J, nY, M))
			colP <- matrix(colP, M, nY, byrow = TRUE)
			extP <- matrix(extP, M, nY, byrow = TRUE)
			
			## create transition matrices (phi^T)
			phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
			for(i in 1:M) {
				for(t in 1:(nY-1)) {
					phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], 
						extP[i,t], 1-extP[i,t])) 
				}
			}
			
			## first compute latent probs
			x <- array(NA, c(2, nY, M))
			x[,1,] <- rep(c(1-psi,psi), M)
			for(i in 1:M) {
				for(t in 2:nY) {
					x[,t,i] <- (phis[,,t-1,i] %*% x[,t-1,i])
				}
			}

			## then compute obs probs
			fitted <- array(NA, c(J, nY, M))
			for(i in 1:M) {
				for(t in 1:nY) {
					for(j in 1:J) {
						fitted[j,t,i] <- (x[,t,i] %*% matrix(c(1, 1 - detP[j,t,i], 0, detP[j,t,i]), 2, 2))[2]
					}
				}
			}	
				
			return(matrix(fitted, M, J*nY, byrow = TRUE))
			
		})



setMethod("profile", "unmarkedFit",
		function(fitted, type, parm, seq) {
			stopifnot(length(parm) == 1)
			MLE <- mle(fitted)
			SE(fitted[type])
			nPar <- length(mle(fitted))
			nll <- nllFun(fitted)
			
			## create table to match parm numbers with est/parm numbers
			types <- names(fitted)
			numbertable <- data.frame(type = character(0), num = numeric(0))
			for(i in seq(length=length(types))) {
				length.est <- length(fitted[i]@estimates)
				numbertable <- rbind(numbertable, 
				data.frame(type = rep(types[i], length.est), 
				num = seq(length=length.est)))
			}
			parm.fullnums <- which(numbertable$type == type & 
				numbertable$num == parm)
			
			f <- function(value) {
				fixedNLL <- genFixedNLL(nll, parm.fullnums, value)
				mleRestricted <- optim(rep(0,nPar), fixedNLL)$value
				mleRestricted
			}
			prof.out <- sapply(seq, f)
			prof.out <- cbind(seq, prof.out)
			new("profile", prof = prof.out)
		})

setMethod("hessian", "unmarkedFit",
		function(object) {
			object@opt$hessian
		})




setMethod("update", "unmarkedFit", 
	function(object, formula., ..., evaluate = TRUE) 
{
	call <- object@call
	origFormula <- formula(call)
	if (is.null(call)) 
		stop("need an object with call slot")
	extras <- match.call(expand.dots = FALSE)$...
	if (!missing(formula.)) {
		detformula <- as.formula(origFormula[[2]])
		stateformula <- as.formula(paste("~", origFormula[3], sep=""))
		newDetformula <- as.formula(formula.[[2]])
		upDetformula <- update.formula(detformula, newDetformula)
		newStateformula <- as.formula(paste("~", formula.[3], sep=""))
		upStateformula <- update.formula(stateformula, newStateformula)
		call$formula <- as.formula(paste(deparse(upDetformula), 
			deparse(upStateformula)))
		}
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
})


setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))
setMethod("sampleSize", "unmarkedFit",
		function(object) {
			data <- getData(object)
			M <- numSites(data)
			M <- M - length(object@sitesRemoved)
			M
		})
	
setGeneric("getData", function(object) standardGeneric("getData"))	
setMethod("getData", "unmarkedFit",
		function(object) {
			object@data
		})


setGeneric("nllFun", function(object) standardGeneric("nllFun"))
setMethod("nllFun", "unmarkedFit", function(object) object@nllFun)

setGeneric("mle", function(object) standardGeneric("mle"))
setMethod("mle", "unmarkedFit", function(object) object@opt$par)

setClass("profile", representation(prof = "matrix"))


setMethod("plot", c("profile", "missing"),
		function(x) {
			plot(x@prof[,1], x@prof[,2], type = "l")
		})


setMethod("residuals", "unmarkedFit", function(object, ...)
	{
		y <- getY(object@data)
		e <- fitted(object, na.rm = FALSE)	 
		r <- y - e
		return(r)
	})
	
setMethod("residuals", "unmarkedFitOccu", function(object, ...)
	{
		y <- getY(object@data)
		y <- truncateToBinary(y)
		e <- fitted(object, na.rm = FALSE)	 
		r <- y - e
		return(r)
	})
	
setMethod("residuals", "unmarkedFitOccuRN", function(object, ...)
	{
		y <- getY(object@data)
		y <- truncateToBinary(y)
		e <- fitted(object, na.rm = FALSE)	 
		r <- y - e
		return(r)
	})




setMethod("plot", c(x = "unmarkedFit", y = "missing"), 
	function(x, y, ...)
{
	r <- residuals(x)
	e <- fitted(x, na.rm = FALSE)
	plot(e, r, ylab = "Residuals", xlab = "Predicted values")
	abline(h = 0, lty = 3, col = "gray")
})
	
	
 		
		



############################# CHILD CLASS METHODS ##############################

# Extract detection probs
setGeneric("getP", function(object, ...) standardGeneric("getP"))


setMethod("getP", "unmarkedFit", function(object, na.rm = TRUE) 
	{
		formula <- object@formula
		detformula <- as.formula(formula[[2]])
		umf <- object@data
		designMats <- getDesign2(formula, umf, na.rm = na.rm)
		y <- designMats$y
		V <- designMats$V
		M <- nrow(y)
		J <- ncol(y)
		ppars <- coef(object, type = "det")
		p <- plogis(V %*% ppars)
		p <- matrix(p, M, J, byrow = TRUE)
		return(p)
	})




setMethod("getP", "unmarkedFitDS", function(object, na.rm = TRUE) 
{
	formula <- object@formula
	detformula <- as.formula(formula[[2]])
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	V <- designMats$V
	M <- nrow(y)
	J <- ncol(y)
	ppars <- coef(object, type = "det")
	d <- umf@dist.breaks
	survey <- umf@survey
	key <- object@keyfun
	switch(key, 
		halfnorm = {
			sigma <- exp(V %*% ppars)
			p <- sapply(sigma, function(x) cp.hn(d = d, s = x, survey = survey))
			}, 
		exp = {
			rate <- exp(V %*% ppars)
			p <- sapply(rate, function(x) cp.exp(d = d, r = x, survey = survey))
			}, 
		hazard = {
			shape <- exp(V %*% ppars)
			scale <- exp(coef(object, type="scale"))
			p <- sapply(shape, function(x) cp.haz(d = d, shape = x, 
				scale = scale, survey = survey))
			})
    p <- matrix(p, M, J, byrow = TRUE)
	return(p)
})






setMethod("getP", "unmarkedFitMPois", function(object, na.rm = TRUE) 
	{
		formula <- object@formula
		detformula <- as.formula(formula[[2]])
		piFun <- object@data@piFun
		umf <- object@data
		designMats <- getDesign2(formula, umf, na.rm = na.rm)
		y <- designMats$y
		V <- designMats$V
		M <- nrow(y)
		J <- ncol(y)
		ppars <- coef(object, type = "det")
		p <- plogis(V %*% ppars)
		p <- matrix(p, M, J, byrow = TRUE)
		pi <- do.call(piFun, list(p = p))
		return(pi)
	})





setMethod("simulate", "unmarkedFitDS", 
	function(object, nsim = 1, seed = NULL, na.rm=TRUE)
{
	formula <- object@formula
	umf <- object@data
	designMats <- getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	X <- designMats$X
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	lamParms <- coef(object, type = "state")
	lam <- as.numeric(exp(X %*% lamParms))
	lamvec <- rep(lam, each = J) * a
	pvec <- c(t(getP(object, na.rm = na.rm)))
	simList <- list()
	for(i in 1:nsim) {
		yvec <- rpois(M * J, lamvec * pvec)
		simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
		}
	return(simList)
})




setMethod("simulate", "unmarkedFitPCount", 
	function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
	formula <- object@formula
	umf <- object@data
	designMats <- unmarked:::getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	X <- designMats$X
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	allParms <- coef(object, altNames = FALSE)
	lamParms <- coef(object, type = "state")
	lam <- as.numeric(exp(X %*% lamParms)) * a
	lamvec <- rep(lam, each = J)
	pvec <- c(t(getP(object, na.rm = na.rm)))
	mix <- object@mixture
	simList <- list()
	for(i in 1:nsim) {
		switch(mix, 
			P = yvec <- rpois(M * J, lamvec * pvec),
			NB = {
				N <- rnbinom(M, size = exp(coef(object["alpha"])), mu = lam)
				yvec <- rbinom(M * J, size = rep(N, each = J), prob = pvec)
				}
			)
		simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
		}
	return(simList)
})



setMethod("simulate", "unmarkedFitMPois", 
	function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
	formula <- object@formula
	umf <- object@data
	designMats <- unmarked:::getDesign2(formula, umf, na.rm = na.rm)
	y <- designMats$y
	X <- designMats$X
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	lamParms <- coef(object, type = "state")
	lam <- as.numeric(exp(X %*% lamParms))
	lamvec <- rep(lam, each = J) * a
	pivec <- as.vector(t(getP(object, na.rm = na.rm)))
	simList <- list()
	for(i in 1:nsim) {
		yvec <- rpois(M * J, lamvec * pivec)
		simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
		}
	return(simList)
})




setMethod("simulate", "unmarkedFitOccu", 
		function(object, nsim = 1, seed = NULL, na.rm = TRUE)
		{
			formula <- object@formula
			umf <- object@data
			designMats <- getDesign2(formula, umf, na.rm = na.rm)
			y <- designMats$y
			X <- designMats$X
			M <- nrow(y)
			J <- ncol(y)
			allParms <- coef(object, altNames = FALSE)
			psiParms <- coef(object, type = "state")
			psi <- as.numeric(plogis(X %*% psiParms))
			p <- c(t(getP(object,na.rm = na.rm)))
			simList <- list()
			for(i in 1:nsim) {
				Z <- rbinom(M, 1, psi)
				Z[object@knownOcc] <- 1
				yvec <- rep(Z, each = J)*rbinom(M * J, 1, prob = p)
				simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
			}
			return(simList)
		})


setMethod("simulate", "unmarkedFitColExt",
		function(object, nsim = 1, seed = NULL, na.rm = TRUE) {
			data <- object@data
			psi <- plogis(coef(object, 'psi'))
			detParms <- coef(object, 'det')
			colParms <- coef(object, 'col')
			extParms <- coef(object, 'ext')
			designMats <- unmarked:::getDesign3(formula = object@formula, object@data)
			V.itj <- designMats$V; X.it <- designMats$X

			M <- nrow(X.it)
			nY <- data@numPrimary
			J <- obsNum(data)/nY

			detP <- plogis(V.itj %*% detParms)
			colP <- plogis(X.it  %*% colParms)
			extP <- plogis(X.it %*% extParms)
			
			detP <- array(detP, c(J, nY, M))
			detP <- aperm(detP, c(3, 1, 2))
			colP <- matrix(colP, M, nY, byrow = TRUE)
			extP <- matrix(extP, M, nY, byrow = TRUE)
			
			simList <- list()
			for(s in 1:nsim) {
				## generate first year's data
				x <- matrix(0, M, nY)
				x[,1] <- rbinom(M, 1, psi) 
				
				## create transition matrices (phi^T)
				phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
				for(i in 1:M) {
					for(t in 1:(nY-1)) {
						phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t], 1-extP[i,t])) 
					}
				}
				
				## generate latent years 2:T
				for(i in 1:M) {
					for(t in 2:nY) {
						x[i,t] <- rbinom(1, 1, phis[2,x[i,t-1]+1,t-1,i])
					}
				}
				
				## generate observations
				y <- array(NA, c(M, J, nY))
				for(t in 1:nY) {
					y[,,t] <- rbinom(M*J, 1, x[,t]*detP[,,t])
				}
				
				y.mat <- y[,,1]
				for(i in 2:dim(y)[3]) {
					y.mat <- cbind(y.mat,y[,,i])
				}
				simList[[s]] <- y.mat
			}
			
			return(simList)
			
		})
		
		


setMethod("simulate", "unmarkedFitOccuRN",
		function(object, nsim = 1, seed = NULL, na.rm = TRUE) {
			formula <- object@formula
			umf <- object@data
			designMats <- unmarked:::getDesign2(formula, umf, na.rm = na.rm)
			y <- designMats$y; X <- designMats$X; V <- designMats$V
			M <- nrow(y)
			J <- ncol(y)
			
			detParms <- coef(object, 'det')
			r.ij <- plogis(V %*% detParms)
			r <- matrix(r.ij, M, J, byrow = TRUE)
			
			lamParms <- coef(object, 'state')
			lambda <- exp(X %*% lamParms)

			simList <- list()
			for(s in 1:nsim) {
				N.i <- rpois(M, lambda)
				N.ij <- rep(N.i, each = J)
				y <- matrix(NA, M, J)
				for(i in 1:J) {
					y[,i] <- rbinom(M, N.i, r[,i])
				}
				simList[[s]] <- ifelse(y > 0, 1, 0)
			}
			return(simList)
		})



############################## PARAMETRIC BOOTSTRAP ###########################

setGeneric("parboot",
    def = function(object, ...) {
      standardGeneric("parboot")
    })


setClass("parboot",
    representation(fitType = "character",
        call = "call",
        t0 = "numeric",
        t.star = "numeric")
        )


# Evaluate goodness-of-fit of a fitted model. 
setMethod("parboot", "unmarkedFit", function(object, nsim=10, report=2, ...)  
{
	call <- match.call(call = sys.call(-1))
	formula <- object@formula
	umf <- object@data
	designMats <- unmarked:::getDesign2(formula, umf, na.rm = FALSE)
	y <- designMats$y
	if(class(object) %in% c("unmarkedFitOccu", "unmarkedFitOccuRN"))
		y <- truncateToBinary(y)
	yvec0 <- c(t(y))
	ests <- as.numeric(coef(object, altNames = TRUE))
	expected0 <- as.vector(t(fitted(object, na.rm = FALSE))) 
	rmse0 <- sqrt(sum((sqrt(yvec0) - sqrt(expected0))^2, na.rm = TRUE))
	cat("t0 =", rmse0, "\n")      
	rmse <- numeric(nsim)
	fits <- list()
	simdata <- umf
	simList <- simulate(object, nsim = nsim, na.rm = FALSE)
	for(i in 1:nsim) {
		y.sim <- simList[[i]]
		is.na(y.sim) <- is.na(y)
		yvec <- c(t(y.sim))
		simdata@y <- y.sim
		fits[[i]] <- update(object, data = simdata, starts = ests, se = F, ...)
		expected <- as.vector(t(fitted(fits[[i]], na.rm = FALSE)))
		rmse[i] <- sqrt(sum((sqrt(yvec) - sqrt(expected))^2, na.rm = TRUE))
		if(nsim > report && i %in% seq(report, nsim, by=report))
			cat(paste(round(rmse[(i-(report-1)):i], 1), collapse=", "), fill=T)
		}
	out <- new("parboot", call=call, t0 = rmse0, t.star = rmse)
	return(out)
})






setMethod("show", "parboot", function(object) 
		{
			t.star <- object@t.star
			t0 <- object@t0
			bias <- mean(t0 - t.star)
			bias.se <- sd(t0 - t.star)
			nsim <- length(t.star)
			p.val <- sum(abs(t.star - 1) > abs(t0 - 1)) / (1 + nsim)
			stats <- c("original" = t0, "bias" = bias, "Std. error" = bias.se, 
					"p.value" = p.val)
			cat("\nCall:", deparse(object@call), fill=T)
			cat("\nBootstrap Statistics:\n")
			print(stats, digits=3)
			cat("\nt quantiles:\n")
			print(quantile(t.star, probs=c(0,2.5,25,50,75,97.5,100)/100))        
		})


## nonparboot return entire list of fits... they will be processed by vcov, confint, etc.
setGeneric("nonparboot", function(object, B = 50, ...) {standardGeneric("nonparboot")})

setMethod("nonparboot", "unmarkedFit",
		function(object, B = 50, ...) {
			data <- object@data
			formula <- object@formula
			designMats <- getDesign2(formula, data)  # bootstrap only after removing sites
			removed.sites <- designMats$removed.sites
			data <- data[-removed.sites,]
			M <- numSites(data)
			boot.iter <- function(x) {
				sites <- sort(sample(1:M, M, replace = TRUE))
				fm <- update(object, data = data[sites,])
			}
			fmList <- lapply(1:B, boot.iter)
			fmList
		})


setMethod("plot", signature(x="parboot", y="missing"), function(x, y, ...)
		{
			op <- par(mfrow=c(1, 2))
			t.star <- x@t.star
			t0 <- x@t0
			t.t0 <- c(t.star, t0)
			bias <- mean(t0 - t.star)
			bias.se <- sd(t0 - t.star)
			nsim <- length(t.star)
			p.val <- sum(abs(t.star - 1) > abs(t0 - 1)) / (1 + nsim)
			hist(t.star, xlim=c(min(floor(t.t0)), max(ceiling(t.t0))), 
				main=paste("P =", round(p.val, 3), "; nsim =", format(nsim)), 
				xlab="t*")
			rug(t.star)
			abline(v=t0, lty=2)
			qqnorm(t.star, ylab="t*")
			qqline(t.star)
			title(outer=T, ...)
			par(op)
		})



		