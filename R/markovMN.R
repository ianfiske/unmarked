
# TODO:  improve documentation!!
# TODO:  improve mechanism for choosing the detection matrix!

# Estimate parameters of the general multiseason multistate occpancy model.
umHMM <-
  function(formula = ~ 1, umf,
           detconstraint = NULL,
           phiconstraint = NULL, psiconstraint = NULL,
           phiMat, inits,
           B = 0,
           max.bad = 3, min.good = 1, fit.stats = TRUE, method = "BFGS", control=list())
{

	phiMat <- get(phiMat)
	stopifnot(is(phiMat, "phiFunSpec"))
	K <- phiMat@K

	## truncate at K
  umf@y[umf@y > K] <- K

  ## coerce vector detcon for K = 1.
  if(identical(K,1)) {
    if(class(detconstraint) %in% c("numeric","integer")) {
      detconstraint <-  t(as.matrix(detconstraint))
    }
  }


  ##################################
  ## section determines appropriate default detconstraint...
  ## needs more investigation.
  nDCP <- length(attr(terms(formula), "term.labels")) + 1
#  if(arDet)
#    nDMP <- K
#  else
    nDMP <-  K*(K+1)/2

  if(is.null(detconstraint))
    detconstraint <- matrix(rep(1:nDMP,nDCP),nDMP, nDCP)
    #detconstraint <- matrix(1:(nDMP*nDCP), nDMP, nDCP)
#################################


  ## only do NA checking for variables being used!
  ## create reduced detection formula from detconstraint
  to.rm <- which(apply(detconstraint == 0, 2, all))
  if(length(to.rm) != 0 & length(to.rm) != (ncol(detconstraint) - 1)) {
    detformula.red <- eval(parse(text=paste("~ ",
                               paste(attr(terms(formula), "term.labels")[-(to.rm-1)],
                                     collapse="+"))))
  } else if(length(to.rm) == (ncol(detconstraint) - 1)) {
    detformula.red <- ~1
  } else {
    detformula.red <- formula
  }

	# fix state formula to "~1"
	formula <- as.formula(paste("~",as.character(detformula.red[2]),"~1"))
  #umf <- handleNA2(formula, umf)
  y <- getY(umf)
  J <- numY(umf) / umf@numPrimary

  M <- nrow(y)
  nY <- ncol(y)/J
  n.det <- sum(apply(y > 0, 1, any, na.rm = TRUE))

  fc <- match.call()
  fc[[1]] <- as.name("umHMM.fit")
	fc$formula <- as.name("formula")
  fc$bootstrap.se <- fc$covdata.site <- fc$covdata.obs <- fc$data <- fc$phiMat <- 
    fc$B <- fc$fit.stats <- NULL
  fc$umf <- as.name("umf")
  fc$J <- as.name("J")
  fc$detconstraint <- as.name("detconstraint")
	fc$method <- as.name("method")
	fc$control <- as.name("control")

	fc$K <- K
	fc$phiMatFun <- phiMat@fun
	fc$nPhiP.un <- phiMat@nPhiP.un

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
			
			umf.b <- umf
			umf.b@y <- umf@y[sites.unique,]
			obsCovsInd.b <- as.vector(t(obsCovs.siteInd[sites.unique,]))
			umf.b@obsCovs <- umf.b@obsCovs[obsCovsInd.b,]
			fc$umf <- quote(umf.b)
			fc$getHessian <- FALSE
			fc$inits <- fm$mle$value
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

  ## check hessian for NaN's, Inf's, 0's, etc.
  ## Their presence should trigger "convergence <- 1"
  if(any(is.na(fm$hessian)) || any(is.infinite(fm$hessian)) ||
      identical(sum(abs(fm$hessian)), 0)) {
    fm$convergence <- 1
  }

  if(fit.stats) {
    fm$fit.stats <- computeFitStats(detMats = fm$detMats, smooth = fm$smooth.sites, y = y, 
			J.it = fm$J.it)
  }

  return(fm)
}


umHMM.fit <- function(formula = ~ 1, umf,
		detconstraint = NULL,
		phiconstraint = NULL, psiconstraint = NULL, J, K, nPhiP.un,
		inits,
		phiMatFun, max.bad = 3, min.good = 3, method, control, getHessian = TRUE, wts)
{
	
	designMats <- getDesign2(formula = formula, umf)
	V.itj <- designMats$V
	y <- designMats$y
	M <- nrow(y)
	nY <- ncol(y)/J
	if(missing(wts)) wts <- rep(1, M)
	
	nDCP <- ncol(V.itj)
	detParms <- colnames(V.itj)
	
	nDMP <-  K*(K+1)/2
	
	if(nrow(detconstraint) != nDMP)
		stop(paste("detconstraint has wrong number of rows.\n
								It should have",nDMP))
	
	## for performance, remove variables eliminated by detconstraint
	to.rm <- which(apply(detconstraint == 0, 2, all))
	if(length(to.rm) != 0) {
		detconstraint <- detconstraint[,-to.rm]
		V.itj <- V.itj[,-to.rm]
		detParms <- detParms[-to.rm]
		nDCP <- nDCP - length(to.rm)
	}
	if(is.null(dim(detconstraint)))
		detconstraint <- matrix(detconstraint,nDMP,nDCP)
	
	## create design matrix for each matrix parameter
	# TODO: make this do the right thing if there are constraints.
	V.itjk <- V.itj %x%  diag(nDMP)
	#get a better line here    if(is.null(detconstraint)) detconstraint <- 1:nDMP.un
	
	nSP.un <- K                # number of parameters for psi vector of initial
	
	if(is.null(phiconstraint)) phiconstraint <- 1:nPhiP.un
	if(is.null(psiconstraint)) psiconstraint <- 1:nSP.un
	
	nPhiP <- max(phiconstraint)
	nSP <- max(psiconstraint)
	
	## create linked list of parameters
	theta.df <- data.frame(parameter = character(), start = numeric(),
			end = numeric(), stringsAsFactors = FALSE)
	
	theta.df <- addParm(theta.df, "phiParms", nPhiP)
	
	## convert from new easy-input style detcon to computational format
	prev.col.max <- 0
	detcon.corr <- matrix(0, nDMP, nDCP)
	for(j in 1:nDCP) {
		detcon.corr[,j] <- ifelse(detconstraint[,j] != 0,
				detconstraint[,j] + prev.col.max, 0)
		prev.col.max <- max(c(detcon.corr[,j],prev.col.max))
	}
	detcon.vec <- as.vector(detcon.corr)
	nDP <- max(detcon.vec)
	nDP.un <- length(detcon.vec)
	
	
	H.det <- matrix(0, nDP.un, nDP)
	for(i in 1:length(detcon.vec)){
		H.det[i,detcon.vec[i]] <- 1
	}
	
	theta.df <- addParm(theta.df, "detParms", nDP)
	nP <- nDP + nSP + nPhiP  # total number of parameters
	
	## construct constrain matrix for phi paramters
	H.phi <- matrix(0, nPhiP.un, nPhiP)
	for(i in 1:nPhiP.un){
		H.phi[i, phiconstraint[i]] <- 1
	}
	
	H.psi <- matrix(0, nSP.un, nSP)
	for(i in 1:nSP.un){
		H.psi[i, psiconstraint[i]] <- 1
	}
	
	y.itj <- as.numeric(t(y))
	
	## replace NA's with 99 before passing to C++
	## TODO: need better missing data passing mechanism (maybe NaN of Inf?)
	y.itj[is.na(y.itj)] <- 99
	V.itjk[is.na(V.itjk)] <- 9999
	
	
	# get ragged array indices
	y.it <- matrix(t(y), nY*M, J, byrow = TRUE)
	J.it <- rowSums(!is.na(y.it))
	
	## OLD findMLE
	alpha <- array(NA, c(K + 1, nY, M))
	beta <- array(NA, c(K + 1, nY, M))
	gamma <- array(NA, c(K + 1, nY, M))
	
	V.arr <- array(t(V.itjk), c(nDP, nDMP, J, nY, M))
	V.arr <- aperm(V.arr, c(2,1,5,4,3))
	
	y.arr <- array(y.itj, c(J, nY, M))
	y.arr <- aperm(y.arr, c(3:1))
	storage.mode(J.it) <- storage.mode(y.arr) <- storage.mode(K) <- "integer"
	
	forward <- function(detParms, phi, psi, storeAlpha = FALSE) {
		
		negloglike <- 0
		psiSite <- matrix(psi, K + 1, M)
		
		mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
		for(t in 1:nY) {
			storage.mode(t) <- "integer"
			detVecs <- .Call("getDetVecs", y.arr, mp, J.it[seq(from = t,to = length(J.it)-nY+t, by=nY)], t, K,
					PACKAGE = "unmarked")
			psiSite <- psiSite * detVecs
			if(storeAlpha) alpha[,t,] <<- psiSite[,]
			if(t < nY) {
				psiSite <- phi %*% psiSite
			} else {
				negloglike <- negloglike - sum(wts*log(colSums(psiSite)))
			}
		}
		negloglike
	}
	
	backward <- function(detParams, phi, psi) {
		for(i in 1:M) {
			backP <- rep(1, K + 1)
			for(t in nY:1) {
				
				beta[, t, i] <<- backP
				
				detVec <- rep(1, K + 1)
				for(j in 1:J) {
					if(y.arr[i,t,j] != 99) {
						mp <- V.arr[,,i,t,j] %*% detParams
						detVecObs <- .Call("getSingleDetVec", y.arr[i,t,j], mp, K, PACKAGE = "unmarked")
						detVec <- detVec * detVecObs
					}
				}
				
				backP <- t(phi) %*% (detVec * backP)
				
			}
		}
	}
	
	getDetMats <- function(detParams, phi, psi) {
		detMats <- array(NA, c(K+1,K+1,J,nY,M))
		for(i in 1:M) {
			for(t in nY:1) {
				for(j in 1:J) {
					if(y.arr[i,t,j] != 99) {
						for(k in 0:K) {
							mp <- V.arr[,,i,t,j] %*% detParams
							detMats[,k+1,j,t,i] <- .Call("getSingleDetVec", k, mp, K, PACKAGE = "unmarked")
						}
					}
				}
			}
		}
		return(detMats)
	}
	
	
	nll <- function(params) {
		
		psiParams <- c(0, H.psi %*% params[1:nSP])
		psi <- exp(psiParams)
		psi <- psi/sum(psi)
		
		phi = do.call(phiMatFun, list(p=H.phi %*% params[(nSP + 1):(nSP + nPhiP)]))
		detParams = H.det %*% params[(nSP + nPhiP+1):(nSP + nDP + nPhiP)]
		forward(detParams, phi, psi) + 0.01*sqrt(sum(params^2))
		
	}
	
	fmList <- list()
	GF <- 0
	BF <- 0
	run <- 1
	
	get.starts <- function() {
		y <- matrix(y.itj, M, J*nY, byrow = TRUE)
		
		## first get initial estimate of phi assuming N = N_max
		
		
		# commented out the wts stuff b/c we don't get inits for bstrap.
#		if(is.null(wts)) wts <- rep(1,nrow(y))
#		y.row <- split(y,1:nrow(y))
#
#		y2 <- lapply(seq_along(y.row), function(r){
#					matrix(rep(y.row[[r]], each = wts[r]), wts[r], ncol(y))
#				})
#		
#		y3 <- matrix(NA,0,ncol(y))
#		for(r in seq_along(y2)) {
#			y3 <- rbind(y3, y2[[r]])
#		}
#		y <- y3
#		M <- nrow(y)
		
		y[y == 99] <- NA
		y.it <- matrix(t(y), M*nY, J, byrow=TRUE)
		y.it <- apply(y.it, 1, max, na.rm=T)
		y.it[is.infinite(y.it)] <- NA
		y.it <- data.frame(year = rep(1:nY,M), y = as.factor(y.it))
		
		trans <- data.frame(year = numeric(0), y.from=numeric(0), y.to=numeric(0), Freq=numeric(0))
		for(i in 1:(nY-1)) {
			y.temp <- subset(y.it, year==i | year==i+1)  ## 1,2,1,2 (odd-even)
			y.from <- y.temp$y[row(y.temp) %% 2 == 1]
			y.to <- y.temp$y[row(y.temp) %% 2 == 0]
			y.tab <- table(y.from,y.to)
			trans <- rbind(trans,cbind(year=i,as.data.frame(y.tab)))
		}
		
		trans.agg <- with(trans, aggregate(Freq, list(y.from = y.from, y.to = y.to), sum))
		trans.mat <- matrix(trans.agg$x, K+1, K+1)
		trans.mat <- trans.mat / rowSums(trans.mat)
		
		trans.penalty <- function(phiParms) {
			phiParms <- H.phi %*% phiParms
			phi <- do.call(phiMatFun, list(p=H.phi %*% phiParms))
			sum((phi - trans.mat)^2)
		}
		phi.opt <- optim(rep(0,nPhiP), trans.penalty)
		phi.start <- phi.opt$par 
		
		y.it.mat <- matrix(as.numeric(as.character(y.it[,2])), M, nY, byrow = TRUE)
		dParms.nll <- function(detParms) {
			negloglike <- 0
			detParms <- H.det %*% detParms
			mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
			for(t in 1:nY) {
				storage.mode(t) <- "integer"
				detVecs <- matrix(.Call("getDetVecs", y.arr, mp, J.it[seq(from = t,to = length(J.it)-nY+t, by=nY)], t, K, 
								PACKAGE = "unmarked"),
						K+1, M)
				nll.t <- sum(log(detVecs[matrix(c(y.it.mat[,t] + 1, 1:M), M, 2)]), na.rm = TRUE)
				negloglike <- negloglike - nll.t
			}
			negloglike
		}
		
		det.opt <- optim(rep(0, nDP), dParms.nll, method = "BFGS", hessian = FALSE)
		det.start <- det.opt$par
		
		psi.mom <- table(y.it.mat[,1])		
		psi.mom2 <- numeric(K+1)
		names(psi.mom2) <- 0:K
		psi.mom2[names(psi.mom)] <- psi.mom 
		psi.mom2 <- psi.mom2/sum(psi.mom2)
		psi.penalty <- function(psiParms) {
			psiParms <- c(0, H.psi %*% psiParms)
			psi <- exp(psiParms)
			psi <- psi/sum(psi)
			sum((psi - psi.mom2)^2)
		}		
		psi.opt <- optim(rep(0,nSP), psi.penalty)
		psi.start <- psi.opt$par
		
		cat("Starting values found.\n")
		return(c(psi.start, phi.start, det.start))
	}
	
	
	if(missing(inits)) {
		start <- get.starts()
	} else {
		start <- inits
	}
	
#  while(GF < min.good && BF < max.bad) {
#    starts <- matrix(rnorm(nP*20), 20, nP)  # sample space initially for better start
#    starts.nll <- apply(starts, 1, nll)
#    start <- starts[which.min(starts.nll),]
#
#    fmList[[run]] <- optim(start, nll, method=method,hessian = TRUE,
#        control=control)
#
#		GF <- GF + 1
#		run <- run + 1  ## TODO: track bad fits.
#  }
#  nlls <- sapply(fmList, function(x) x$value)
#  fm <- fmList[[which.min(nlls)]]
	
#	num.tries <- 0
#	while(num.tries < 3) {
#		fm <- tryCatch(optim(start, nll, method=method,hessian = getHessian,
#					control=control),
#				error = function(x) {
#					num.tries <- num.tries + 1
#					start <- start + 0.1*rnorm(length(start))
#					simpleError("optim hit error.")
#				})
#		if(!identical(as.character(fm), "Error: optim hit error.\n")) break
#	}
#	if(identical(as.character(fm),"Error: optim hit error.\n")) stop("optim hit error.")
	fm <- optim(start, nll, method=method,hessian = getHessian,	control=control)
	opt<-fm
	
	mle <- fm$par
	psiParams <- c(0, H.psi %*% mle[1:nSP])
	psi <- exp(psiParams)
	psi <- psi/sum(psi)
	
	phiParams <- H.phi %*% mle[(nSP+1):(nSP + nPhiP)]
	phi = t(do.call(phiMatFun, list(p=phiParams)))
	detParams = as.vector(H.det %*% mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)])
	
	## smoothing
	forward(detParams, phi, psi, storeAlpha = TRUE)
	backward(mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)], phi, psi)
	beta[,,1]
	for(i in 1:M) {
		for(t in 1:nY) {
			gamma[,t,i] <- alpha[,t,i] * beta[,t,i] / sum(alpha[,t,i] * beta[,t,i])
		}
	}
	
	## get expected detection matrices
#  mp <- array(V.itjk %*% detParams, c(nDMP, J, nY, M))
#  expectedDetMats <- .Call("getDetMats", y.arr, mp, K)
#  expectedDetMats <- array(expectedDetMats, c(K+1,K+1,J,nY,M))
	detMats <- getDetMats(mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)], phi, psi)
	
	## call viterbi here using detMats
	
	## compute the projected trajectory
	projected <- matrix(NA, K + 1, nY)
	projected[,1] <- psi
	for(year in 2:nY) {
		projected[,year] <- t(phi) %*%
				projected[,year-1]
	}
	
	
	fm <- list(nll = fm$value - 0.01*sqrt(sum(mle^2)), phi = phi,
			psiParms = mle[1:nSP], detParms = mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)],
			phiParms = mle[(nSP+1):(nSP + nPhiP)],
			smooth = gamma, projected = projected,
			hessian = fm$hessian,
			convergence = fm$convergence, detMats = detMats, opt = opt)
	
	#### END OLD findMLE
	
	
	phiEsts <- H.phi %*% fm$phiParms
	detEsts <- H.det %*% fm$detParms
	detEsts <- matrix(detEsts, nDMP, nDCP)
	colnames(detEsts) <- detParms
	colnames(detconstraint) <- detParms
	phi <- fm$phi
	psiEsts <- H.psi %*% fm$psiParms
	psi <- exp(c(0,psiEsts))/sum(exp(c(0,psiEsts)))
	ss <- getSS(phi)
	hessian <- fm$hessian
	
	mle <- c(fm$psiParms, fm$phiParms, fm$detParms)
	parm.names <- c(rep("psi", nSP), rep("phi", nPhiP), rep("det", nDP))
	parm.names <- paste(parm.names, c(1:nSP, 1:nPhiP, 1:nDP), sep="")
	mle.df <- data.frame(names = parm.names, value = mle)
	rownames(mle.df) <- 1:nP
	
	## get epsi SE
	###     psi.ind <- grep("psi", rawfit.names)
	###     psi.theta <- ests[psi.ind]
	###     psi.hessian <- fm$hessian[psi.ind,psi.ind]
	###     epsi.se <- sqrt(deltaVar(psi.theta, psi.hessian, meanPsi))
	
	###     ## get ess SE
	###     phi.ind <- grep("phi", rawfit.names)
	###     phi.theta.con <- ests[phi.ind]
	###     phi.hessian <- fm$hessian[phi.ind, phi.ind]
	###     ess.se <- sqrt(deltaVar(phi.theta.con, phi.hessian, meanSS,
	###                             phiconstraint = phiconstraint,
	###                             phifun = phiMatrix))
	
	###     ## smoothing SE
	
#  smooth.mean <- apply(fm$smooth, c(2,3), meanstate)
#  smooth.overall.mean <- rowMeans(smooth.mean)
	
	smooth <- apply(fm$smooth, c(1,2), function(x) {
				sum(wts*x) / sum(wts) # weighted avg
			})
	
	smooth.sites <- fm$smooth
	
	list(mle = mle.df, psiEsts = psiEsts, phiEsts = phiEsts,
			detEsts = detEsts, phi = phi,
			AIC = 2*fm$nll + 2*nP,
			NegLogLike = fm$nll,
			psi = psi, epsi = meanstate(psi), #epsi.se = epsi.se,
			ss = ss, ess = meanstate(ss), #ess.se = ess.se,
			smooth = smooth, projected = fm$projected,
			#arDet = arDet,
			detform = formula,
			psiconstraint = psiconstraint,
			phiconstraint = phiconstraint,
			detconstraint = detconstraint,
			K = K, hessian = hessian,
			n = M, convergence = fm$convergence, detMats = fm$detMats,
			smooth.sites = smooth.sites, J.it = J.it, opt=opt)
}


#findMLE <-
#    function(y.itj, V.itjk, J.it, nDMP, nDP, nSP, nSP.un, nPhiP, nP, nDP.un,
#        nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
#        smooth.only = FALSE,
#        inits,
#        phiMatFun, get.inits = TRUE, 
#        max.bad, min.good, method, control, getHessian, wts) {
#
#
#}


phiLogit2 <- function(p) {
  p0 <- p[1]
  p1 <- p[2]
  matrix(plogis(c(-p0,p0,-p1,p1)),2,2)
}


## fun returns t(phi)
setClass("phiFunSpec",
		representation(fun = "function",
				nPhiP.un = "numeric", K = "numeric"
))

phiLogit4fun <- function(p) {
  p0 <- p[1]
  p1 <- p[2]
  p2 <- p[3]
  p3 <- p[4]
  p4 <- p[5]
  p5 <- p[6]
  p6 <- p[7]
  p7 <- p[8]
  p8 <- p[9]
  p9 <- p[10]
  p10 <- p[11]
  p11 <- p[12]

  matrix(c(1/(1 + exp(p0) + exp(p1) + exp(p2)),
          exp(p0)/(1 + exp(p0) + exp(p1) + exp(p2)),
          exp(p1)/(1 + exp(p0) + exp(p1) + exp(p2)),
          exp(p2)/(1 + exp(p0) + exp(p1) + exp(p2)),
          1/(1 + exp(p3) + exp(p4) + exp(p5)),
          exp(p3)/(1 + exp(p3) + exp(p4) + exp(p5)),
          exp(p4)/(1 + exp(p3) + exp(p4) + exp(p5)),
          exp(p5)/(1 + exp(p3) + exp(p4) + exp(p5)),
          1/(1 + exp(p6) + exp(p7) + exp(p8)),
          exp(p6)/(1 + exp(p6) + exp(p7) + exp(p8)),
          exp(p7)/(1 + exp(p6) + exp(p7) + exp(p8)),
          exp(p8)/(1 + exp(p6) + exp(p7) + exp(p8)),
          1/(1 + exp(p9) + exp(p10) + exp(p11)),
          exp(p9)/(1 + exp(p9) + exp(p10) + exp(p11)),
          exp(p10)/(1 + exp(p9) + exp(p10) + exp(p11)),
          exp(p11)/(1 + exp(p9) + exp(p10) + exp(p11))),4,4) # TODO:  make sure this is transposed correctly
}



phi4RedFun <- function (p) 
{
	p <- exp(-p)
	
	phi <- matrix(c(1, p[1] ^ (1:3),
					p[4], 1, p[2] ^ (1:2),
					p[5] ^ (2:1), 1, p[3],
					p[6] ^ (3:1), 1), 4, 4)
	
	phi / rep(colSums(phi),each=4)
}

phi4RandFun <- function (p) 
{
	p <- exp(p)
	den <- 1 + sum(p)
	matrix(c(1,p[1],p[2],p[3])/den, 4, 4)
}

## TODO: make phi generator.
phiGenerator <- function(K) {
	nPhiP.un <- K*(K-1)
	fun <- function(p) {
		
	}
	phiSpec <- new("phiFunSpec", fun = fun, K = K-1, nPhiP.un = nPhiP.un)
	return(phiSpec)
}

phi4 <- new("phiFunSpec",
		fun = phiLogit4fun,
		K = 3,
		nPhiP.un = 12)

phi4Red <- new("phiFunSpec",
		fun = phi4RedFun,
		K = 3,
		nPhiP.un = 6) 
		
phi4Rand <- new("phiFunSpec",
		fun = phi4RandFun,
		K = 3,
		nPhiP.un = 3)


phi2 <- new("phiFunSpec",
		fun = phiLogit2,
		K = 1,
		nPhiP.un = 2) 


computeFitStats.max <- function(detMats, smooth, y, J.it) {

  dims <- dim(detMats)
  K <- dims[1] - 1
  J <- dims[3]
  nY <- dims[4]
  M <- dims[5]

  ## computes expected cell probabilities for each i,t,j
  computeExpectedObs <- function(detMats, smooth) {
    expProbs <- array(NA,c(K+1,J,nY,M))
    q <- 1
    for(i in 1:M) {
      for(t in 1:nY) {
        psi.t <- smooth[,t,i]
        for(j in seq(length=J.it[q])) {
          expProbs[,j,t,i] <- as.vector(t(psi.t) %*% detMats[,,j,t,i])
        }
        q <- q + 1
      }
    }
    return(expProbs)
  }
  eo <- computeExpectedObs(detMats, smooth)

# given matrix of expected probilities for site i year t to compute the CDF of max
  maxCDF.siteyear <- function(mat) {
    mat[is.na(mat)] <- 0
    CDF.itj <- apply(mat, 2, cumsum)
    CDF.itj[CDF.itj == 0] <- 1
    CDF <- rowProds(CDF.itj)
    if(CDF[1] == 1) CDF <- rep(NA,length(CDF))  ## replace columns of all 1's with NA
    CDF
  }
  maxCDF <- apply(eo, 3:4, maxCDF.siteyear)

  ## convert a discrete CDF to a PDF
  CDFtoPDF <- function(CDF) {
    CDF2 <- c(0,CDF[-length(CDF)])
    CDF - CDF2
  }
  maxPDF <- apply(maxCDF, 2:3, CDFtoPDF)

  MaxDist <- apply(maxPDF, 1:2, sum, na.rm=TRUE)

  ## create y array (multinomial cell style)
  y.arr <- array(t(y), c(J, nY, M))#<- array(NA, c(K+1, J, nY, M))
  y.arr <- apply(y.arr, 1:3, function(x) {
        if(!is.na(x)) {
          v <- numeric(K+1)
          v[x+1] <- 1
          v
        } else {
          rep(NA, K+1)
        }
      })

  X2 <- sum((y.arr - eo)^2 / eo, na.rm = TRUE)


  return(list(MaxDist = MaxDist, eo = eo, X2 = X2))
}

computeFitStats <- function(detMats, smooth, y, J.it) {

	dims <- dim(detMats)
	K <- dims[1] - 1
	J <- dims[3]
	nY <- dims[4]
	M <- dims[5]
		
	## computes expected cell probabilities for each i,t,j
	computeExpectedObs <- function(detMats, smooth) {
		expProbs <- array(NA,c(K+1,J,nY,M))
		q <- 1
		for(i in 1:M) {
			for(t in 1:nY) {
				psi.t <- smooth[,t,i]
				for(j in seq(length=J.it[q])) {
					expProbs[,j,t,i] <- as.vector(t(psi.t) %*% detMats[,,j,t,i])
				}
				q <- q + 1
			}
		}
		return(expProbs)
	}
	eo <- computeExpectedObs(detMats, smooth)
	
	e.counts <- apply(eo, c(1,3), sum, na.rm = TRUE)
	
	y.arr <- array(t(y), c(J, nY, M))
	o.counts <- apply(y.arr, 2, function(x) {
				tab <- table(x, useNA = "no")
				fulltab <- numeric(K+1)
				names(fulltab) <- 0:K
				fulltab[names(tab)] <- tab
				fulltab
			})
	
	X2 <- sum((e.counts - o.counts)^2 / e.counts)
	
	list(e.counts = e.counts, o.counts = o.counts, chisq = X2)
}

viterbi <- function(y.site) {
	delta <- psi.v <- q.star <- matrix(NA, K+1, nY)
		
	## initialize
	for (i in 1:(K + 1)) {
		b.i1 <- rep(1, K + 1) 
		for (j in 1:J) {
			b.i1 <- b.i1 * detMats[,y.arr[i,1,j]+1,j,i,1]
		}
		delta[i,1] <- b.i1 * psi
	}
}