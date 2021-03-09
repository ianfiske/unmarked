# Fit the multinomial-Poisson abundance mixture model.

multinomPois <- function(formula, data, starts, method = "BFGS",
    se = TRUE, engine = c("C","R"), ...)
{
    check_no_support(split_formula(formula))

    if(!is(data, "unmarkedFrameMPois"))
		    stop("Data is not a data frame or unmarkedFrame.")
    engine <- match.arg(engine, c("C", "R"))
    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    if (is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
        }
    J <- ncol(y)
    R <- obsNum(data)
    M <- nrow(y)
    piFun <- data@piFun

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)
    nP <- nDP + nAP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(y)
    navec <- is.na(yvec)

    nll_R <- function(parms) {
        lambda <- exp(X %*% parms[1 : nAP] + X.offset)
        p <- plogis(V %*% parms[(nAP + 1) : nP] + V.offset)
        p.matrix <- matrix(p, M, R, byrow = TRUE)
        pi <- do.call(piFun, list(p = p.matrix))
        logLikeSite <- dpois(y, matrix(lambda, M, J) * pi, log = TRUE)
        logLikeSite[navec] <- 0
        -sum(logLikeSite)
        }

    nll_C <- function(params) {
        .Call("nll_multinomPois",
            params,piFun,
            X, X.offset, V, V.offset,
            yC, navecC, nP,nAP,
            PACKAGE = "unmarked")
    }

    if(engine=="R"){
      nll <- nll_R
    }else{
      yC <- as.numeric(t(y))
      navecC <- is.na(yC)
      nll <- nll_C
      if(!piFun%in%c('doublePiFun','removalPiFun','depDoublePiFun')){
        warning("Custom pi functions are not supported by C engine. Using R engine instead.")
        nll <- nll_R
      }
    }

    if(missing(starts))
        starts <- rep(0, nP)
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    fmAIC <- 2 * fm$value + 2 * nP
    names(ests) <- c(lamParms, detParms)

    stateName <- "Abundance"

    stateEstimates <- unmarkedEstimate(name = stateName,
                                       short.name = "lambda",
                                       estimates = ests[1:nAP],
                                       covMat = as.matrix(
                                           covMat[1:nAP,1:nAP]),
                                       invlink = "exp",
                                       invlinkGrad = "exp")

    detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
                                     estimates = ests[(nAP + 1) : nP],
                                     covMat = as.matrix(covMat[(nAP + 1) :
                                                 nP, (nAP + 1) : nP]),
                                     invlink = "logistic",
                                     invlinkGrad = "logistic.grad")

    estimateList <- unmarkedEstimateList(list(state=stateEstimates,
        det=detEstimates))

    umfit <- new("unmarkedFitMPois", fitType = "multinomPois",
        call = match.call(), formula = formula, data = data,
        estimates = estimateList, sitesRemoved = designMats$removed.sites,
        AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll)

    return(umfit)
}
