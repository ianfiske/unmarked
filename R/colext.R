#' @include utils.R
#' @include unmarkedFit.R
#' @include unmarkedEstimate.R
#' @include markovMN.R
{}

#' Estimate parameters of the colonization-extinction model, including covariate-dependent rates and detection process.
#'
#' 
#' This function fits the colonization-extinction model of MacKenziet et al (2003).  The colonization and extinction rates
#' can be modeled with covariates that vary yearly at each site using a logit link.  These covariates
#' are supplied by special unmarkedMultFrame \code{yearlyCovs} slot.  These parameters are specified using the second formula in the \code{formula} argument.
#' 
#' The conditional detection rate can also be modeled as a function of covariates that vary at the secondary sampling
#' period (ie., repeat visits).  These covariates are specified by the first part of the \code{formula} argument and
#' the data is supplied via the usual \code{obsCovs} slot.
#' 
#' @title Fit the colonization-extinction model.
#' @param formula double right-hand side formula describing covariates of detection and occupancy in that order.
#' @param data unmarkedMultFrame object that supplies the data (see \link{unmarkedFrame}).
#' @param starts optionally, initial values for parameters in the optimization.
#' @param B number of bootstrap interations (the default 0 indicates no bootstrapping).
#' @param method Optimization method used by \code{\link{optim}}.
#' @param control Other arguments passed to \code{\link{optim}}.
#' @param se logical specifying whether or not to compute standard errors.
#' @return unmarkedFitColExt object describing model fit.
#' @examples
#' data(frogs)
#' umf <- formatMult(masspcru)
#' obsCovs(umf) <- scale(obsCovs(umf))
#' (fm <- colext(~ JulianDate + I(JulianDate^2) ~ 1, umf, control = list(trace=1, maxit=1e4)))
#' @keywords models
#' @references 
#' MacKenzie, D.I. et al. (2002) Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology, 83(8), 2248-2255.
#' MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008).
#' @useDynLib unmarked
#' @export
colext <-
		function(formula = ~ 1 ~ 1, data,	starts,	B = 0, method = "BFGS", control=list(), se = TRUE)
{
	
	K <- 1
	
	## truncate at K
	data@y[data@y > K] <- K
	
	y <- getY(data)
	J <- numY(data) / data@numPrimary
	
	M <- nrow(y)
	nY <- ncol(y)/J
	n.det <- sum(apply(y > 0, 1, any, na.rm = TRUE))
	
	fc <- match.call()
	fc[[1]] <- as.name("colext.fit")
	fc$formula <- as.name("formula")
	fc$bootstrap.se <- fc$covdata.site <- fc$covdata.obs <- fc$data <-  
			fc$B <- NULL
	fc$data <- as.name("data")
	fc$J <- as.name("J")
	fc$method <- as.name("method")
	fc$control <- as.name("control")
	fc$getHessian <- as.name("se")
	fc$se <- NULL
	if(missing(starts)) {
		fc$starts <- NULL
	} else {
		fc$starts <- eval(starts)
	}
	
	fm <- eval(fc)
	
	if(B > 0) {  # find bootstrap SEs?
		
		smooth.b <- projected.b <- array(NA, c(K + 1, nY, B))
		psi.b <- matrix(NA, B, K + 1)
		phi.b <- array(NA, c(K+1,K+1,B))
		ss.b <- matrix(NA, B, K + 1)
		mle.b <- matrix(NA, B, nrow(fm$mle))
		
		obsCovs.siteInd <- matrix(1:(M*J*nY),M,J*nY,byrow=T)
		
		b <- 1
		sd <- 0
		bad.boots <- 0
		while(b <= B) {
			sd <- sd + 1
			set.seed(sd)
			samp <- sort(sample(1:M, M, replace = TRUE))
			sites.unique <- unique(samp)
			sites.wts <- table(samp)
			
			data.b <- data
			data.b@y <- data@y[sites.unique,]
			obsCovsInd.b <- as.vector(t(obsCovs.siteInd[sites.unique,]))
			data.b@obsCovs <- data.b@obsCovs[obsCovsInd.b,]
			fc$data <- quote(data.b)
			fc$getHessian <- FALSE
			fc$starts <- fm$mle$value
			fc$wts <- sites.wts
			
			fm.b <- tryCatch(eval(fc),
					error = function(x) {
						bad.boots <<- bad.boots + 1
						cat("Caught failed bootstrap iteration.  Seed =",sd,"\n")
						cat(bad.boots, "bad boot iterations.\n")
						0  ## if fit breaks, then re-enter loop
					})
			if(identical(fm.b, 0)) next
			
			if(!(TRUE %in% is.nan(smooth.b[,,b]))) {
				smooth.b[,,b] <- fm.b$smooth
			} else {
				cat("smooth.b contained NaN.\n")
			}
			psi.b[b,] <- fm.b$psi
			phi.b[,,b] <- fm.b$phi
			mle.b[b,] <- fm.b$mle$value
			ss.b[b,] <- fm.b$ss
			projected.b[,,b] <- fm.b$projected
			
			cat(paste("Bootstrap iteration",b,"completed.\n"))
			b <- b + 1
		}
		
		mle.var <- apply(mle.b, 2, var)
		smooth.var <- apply(smooth.b, 1:2, var)
		projected.var <- apply(projected.b, 1:2, var)
		phi.var <- apply(phi.b, 1:2, var)
		
		## also get the time-series style covariances for K=1 only.
		if(identical(K,1)) {
			smooth.mat <- t(smooth.b[2,,])
			fm$smooth.covmat <- cov(smooth.mat, use="complete.obs")
			fm$smooth.cormat <- cor(smooth.mat, use="complete.obs")
		}
		
		fm$smooth.var <- smooth.var
		fm$projected.var <- projected.var
		fm$mle.var <- mle.var
		fm$phi.var <- phi.var
		
	}
		
	fm$n.det <- n.det
	opt <- fm$opt
	nP <- fm$nP;	M <- fm$M; nDP <- fm$nDP; nPhiP <- fm$nPhiP
	
	if(se) {
		tryCatch(covMat <- solve(opt$hessian),
				error=function(x) simpleError("Hessian is not invertible.  Try using fewer covariates."))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
	ests <- opt$par
	names(ests) <- fm$mle$names
	fmAIC <- 2 * opt$value + 2 * nP + 2*nP*(nP + 1)/(M - nP - 1)
	
	psiParam <- ests[1]
	colParams <- ests[2:(1+nPhiP)]
	extParams <- ests[(2 + nPhiP) : (1 + 2*nPhiP)]
	detParams <- ests[(2 + 2*nPhiP) : nP]
	
	psi <- unmarkedEstimate(name = "Initial", short.name = "psi",
			estimates = psiParam,
			covMat = as.matrix(covMat[1,1]), invlink = "logistic",
			invlinkGrad = "logistic.grad")
	
	col <- unmarkedEstimate(name = "Colonization", short.name = "col",
			estimates = colParams,
			covMat = as.matrix(covMat[2:(1+nPhiP),2:(1+nPhiP)]), invlink = "logistic",
			invlinkGrad = "logistic.grad")
	
	ext <- unmarkedEstimate(name = "Extinction", short.name = "ext",
			estimates = extParams,
			covMat = as.matrix(covMat[(2 + nPhiP) : (1 + 2*nPhiP), (2 + nPhiP) : (1 + 2*nPhiP)]), invlink = "logistic",
			invlinkGrad = "logistic.grad")
	
	det <- unmarkedEstimate(name = "Detection", short.name = "p",
			estimates = detParams,
			covMat = as.matrix(covMat[(2 + 2*nPhiP) : nP, (2 + 2*nPhiP) : nP]), invlink = "logistic",
			invlinkGrad = "logistic.grad")
	
	
	estimateList <- unmarkedEstimateList(list(psi = psi, col = col, ext = ext, det=det))
	
	umfit <- new("unmarkedFitColExt", fitType = "colext",
			call = match.call(), formula = formula, data = data, sitesRemoved = fm$designMats$removed.sites, 
			estimates = estimateList,
			AIC = fmAIC, opt = opt, negLogLike = opt$value, nllFun = fm$nll)
	
	return(umfit)
}


colext.fit <- function(formula, data,		J,
		starts=NULL, method, control, getHessian = TRUE, wts)
{
	K <- 1
	
	designMats <- getDesign3(formula = formula, data)
	V.itjk <- designMats$V; X.it <- designMats$X
	
	detParms <- colnames(V.itjk)
	colextParms <- colnames(X.it)
	
	
	y <- designMats$y
	M <- nrow(y)
	nY <- ncol(y)/J
	if(missing(wts)) wts <- rep(1, M)
	
	stateformula <- as.formula(paste("~",formula[3],sep=""))
	
	## remove final year from X.it
	X.it <- as.matrix(X.it[-seq(nY,M*nY,by=nY),])
	
	nDP <- ncol(V.itjk)
	nDMP <-  1
	nPhiP <- length(colextParms)
	
	## create linked list of parameters
	theta.df <- data.frame(parameter = character(), start = numeric(),
			end = numeric(), stringsAsFactors = FALSE)
	
	theta.df <- addParm(theta.df, "phiParms", nPhiP)
	
	theta.df <- addParm(theta.df, "detParms", nDP)

	nP <- nDP + 1 + 2*nPhiP  # total number of parameters
	
	y.itj <- as.numeric(t(y))
	
	## replace NA's with 99 before passing to C++
	## TODO: need better missing data passing mechanism (maybe NaN of Inf?)
	y.itj[is.na(y.itj)] <- 99
	V.itjk[is.na(V.itjk)] <- 9999
	
	
	# get ragged array indices
	y.it <- matrix(t(y), nY*M, J, byrow = TRUE)
	J.it <- rowSums(!is.na(y.it))
	
	alpha <- array(NA, c(K + 1, nY, M))
	
	V.arr <- array(t(V.itjk), c(nDP, nDMP, J, nY, M))
	V.arr <- aperm(V.arr, c(2,1,5,4,3))
	
	y.arr <- array(y.itj, c(J, nY, M))
	y.arr <- aperm(y.arr, c(3:1))
	storage.mode(J.it) <- storage.mode(y.arr) <- storage.mode(K) <- "integer"
	
	forward <- function(detParms, phis, psi, storeAlpha = FALSE) {
		
		negloglike <- 0
		psiSite <- matrix(c(1-psi,psi), K + 1, M)
		
		mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
		for(t in 1:nY) {
			storage.mode(t) <- "integer"
			detVecs <- .Call("getDetVecs", y.arr, mp, J.it[seq(from = t,to = length(J.it)-nY+t, by=nY)], t, K,
					PACKAGE = "unmarked")
			psiSite <- psiSite * detVecs
			if(storeAlpha) alpha[,t,] <<- psiSite[,]
			if(t < nY) {
				for(i in 1:M) {
					psiSite[,i] <- phis[,,t,i] %*% psiSite[,i]
				}
			} else {
				negloglike <- negloglike - sum(wts*log(colSums(psiSite)))
			}
		}
		negloglike
	}
	
	X <- X.it %x% c(-1,1)
	phis <- array(NA,c(2,2,nY-1,M))
	nll <- function(params) {
		psi <- plogis(params[1])
		colParams <- params[2:(1+nPhiP)]
		extParams <- params[(2 + nPhiP) : (1 + 2*nPhiP)]
		detParams <- params[(2 + 2*nPhiP) : nP]

		phis[,1,,] <- plogis(X %*% colParams) # these are in site-major, year-minor order
		phis[,2,,] <- plogis(X %*% -extParams)
		
		forward(detParams, phis, psi) + 0.01*sqrt(sum(params^2))
	}
	
	if(is.null(starts)) starts <- rnorm(nP)
	fm <- optim(starts, nll, method=method,hessian = getHessian,	control=control)
	
	mle <- fm$par
	parm.names <- c("", colextParms, colextParms, detParms)
	mle.df <- data.frame(names = parm.names, value = mle)
	rownames(mle.df) <- paste(c("psi", rep("col", nPhiP), rep("ext", nPhiP), rep("det", nDP)),
			c("",1:nPhiP,1:nPhiP, 1:nDP))
	
	list(mle = mle.df, opt=fm, nP = nP, M = M, nDP = nDP, nPhiP = nPhiP, nllFun = nll, designMats = designMats)
}



