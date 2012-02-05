



# ----------------- Empirical Bayes Methods ------------------------------



setGeneric("ranef",
    function(object, ...) standardGeneric("ranef"))


setClass("unmarkedRanef1",
    representation(bup = "array"))



setMethod("ranef", "unmarkedFitPCount",
    function(object, ...)
{
    lam <- predict(object, type="state")[,1] # Too slow
    R <- length(lam)
    p <- getP(object)
    K <- object@K
    N <- 0:K
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    bup <- array(NA_real_, c(R, length(N), 1))
    colnames(bup) <- N
    mix <- object@mixture
    for(i in 1:R) {
        switch(mix,
               P  = f <- dpois(N, lam[i]),
               NB = {
                   alpha <- exp(coef(object, type="alpha"))
                   f <- dnbinom(N, mu=lam[i], size=alpha)
               },
               ZIP = {
                   psi <- plogis(coef(object, type="psi"))
                   f <- (1-psi)*dpois(N, lam[i])
                   f[1] <- psi + (1-psi)*exp(-lam[i])
               })
        g <- rep(1, K+1)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(p[i,j]))
                next
            g <- g * dbinom(y[i,j], N, p[i,j])
        }
        fudge <- f*g
        bup[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef1", bup=bup)
})






setMethod("ranef", "unmarkedFitOccu",
    function(object, ...)
{
    psi <- predict(object, type="state")[,1]
    R <- length(psi)
    p <- getP(object)
    z <- 0:1
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    bup <- array(0, c(R,2,1))
    colnames(bup) <- z
    for(i in 1:R) {
        f <- dbinom(z, 1, psi[i])
        g <- rep(1, 2)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(p[i,j]))
                next
            g <- g * dbinom(y[i,j], 1, z*p[i,j])
        }
        fudge <- f*g
        bup[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef1", bup=bup)
})








setMethod("ranef", "unmarkedFitOccuRN",
    function(object, K, ...)
{
    if(missing(K)) {
        warning("You did not specify K, the maximum value of N, so it was set to 50")
        K <- 50
    }
    lam <- predict(object, type="state")[,1] # Too slow
    R <- length(lam)
    r <- getP(object)
    N <- 0:K
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    bup <- array(NA_real_, c(R, length(N), 1))
    colnames(bup) <- N
    for(i in 1:R) {
        f <- dpois(N, lam[i])
        g <- rep(1, K+1)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(r[i,j]))
                next
            p.ijn <- 1 - (1-r[i,j])^N
            g <- g * dbinom(y[i,j], 1, p.ijn)
        }
        fudge <- f*g
        bup[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef1", bup=bup)
})





setMethod("ranef", "unmarkedFitMPois",
    function(object, K, ...)
{
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    if(missing(K)) {
        warning("You did not specify K, the maximum value of N, so it was set to max(y)+50")
        K <- max(y, na.rm=TRUE)+50
    }

    lam <- predict(object, type="state")[,1]
    R <- length(lam)
    cp <- getP(object)
    cp <- cbind(cp, 1-rowSums(cp))
    N <- 0:K
    bup <- array(0, c(R, K+1, 1))
    colnames(bup) <- N
    for(i in 1:R) {
        f <- dpois(N, lam[i])
        g <- rep(1, K+1)
        if(any(is.na(y[i,])) | any(is.na(cp[i,])))
            next
        for(k in 1:(K+1)) {
            yi <- y[i,]
            ydot <- N[k] - sum(yi)
            if(ydot<0) {
                g[k] <- 0
                next
            }
            yi <- c(yi, ydot)
            g[k] <- g[k] * dmultinom(yi, size=N[k], prob=cp[i,])
        }
        fudge <- f*g
        bup[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef1", bup=bup)
})







setMethod("ranef", "unmarkedFitDS",
    function(object, K, ...)
{
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    if(missing(K)) {
        warning("You did not specify K, the maximum value of N, so it was set to max(y)+50")
        K <- max(y, na.rm=TRUE)+50
    }
    lam <- predict(object, type="state")[,1]
    R <- length(lam)
    J <- ncol(y)
    survey <- object@data@survey
    tlength <- object@data@tlength
    db <- object@data@dist.breaks
    w <- diff(db)
    unitsIn <- object@data@unitsIn
    unitsOut <- object@unitsOut
    if(identical(object@output, "density")) {
        a <- matrix(NA, R, J)
        switch(survey, line = {
            for (i in 1:R) {
                a[i, ] <- tlength[i] * w
            }
        }, point = {
            for (i in 1:R) {
                a[i, 1] <- pi * db[2]^2
                for (j in 2:J)
                    a[i, j] <- pi * db[j + 1]^2 - sum(a[i, 1:(j - 1)])
            }
        })
        switch(survey, line = A <- rowSums(a) * 2, point = A <- rowSums(a))
        switch(unitsIn, m = A <- A/1e+06, km = A <- A)
        switch(unitsOut, ha = A <- A * 100, kmsq = A <- A)
        lam <- lam*A
    }
    cp <- getP(object)
    cp <- cbind(cp, 1-rowSums(cp))
    N <- 0:K
    bup <- array(0, c(R, K+1, 1))
    colnames(bup) <- N
    for(i in 1:R) {
        f <- dpois(N, lam[i])
        g <- rep(1, K+1)
        if(any(is.na(y[i,])) | any(is.na(cp[i,])))
            next
        for(k in 1:(K+1)) {
            yi <- y[i,]
            ydot <- N[k] - sum(yi)
            if(ydot<0) {
                g[k] <- 0
                next
            }
            yi <- c(yi, ydot)
            g[k] <- g[k] * dmultinom(yi, size=N[k], prob=cp[i,])
        }
        fudge <- f*g
        bup[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef1", bup=bup)
})







setMethod("ranef", "unmarkedFitGMM",
    function(object, ...)
{
    stop("method not written yet")
})



setMethod("ranef", "unmarkedFitGDS",
    function(object, ...)
{
    stop("method not written yet")
})





setMethod("ranef", "unmarkedFitColExt",
    function(object, ...)
{

#    stop("method not written")

    data <- object@data
    M <- numSites(data)
    nY <- data@numPrimary
    J <- obsNum(data)/nY
    psiParms <- coef(object, 'psi')
    detParms <- coef(object, 'det')
    colParms <- coef(object, 'col')
    extParms <- coef(object, 'ext')
    formulaList <- list(psiformula=object@psiformula,
        gammaformula=object@gamformula,
        epsilonformula=object@epsformula,
        pformula=object@detformula)
    designMats <- unmarked:::getDesign(object@data, object@formula)
    V.itj <- designMats$V
    X.it.gam <- designMats$X.gam
    X.it.eps <- designMats$X.eps
    W.i <- designMats$W

#    yumf <- getY(object@data)
#    yumfa <- array(yumf, c(M, J, nY))

    y <- designMats$y
    ya <- array(y, c(M, J, nY))

    psiP <- plogis(W.i %*% psiParms)
    detP <- plogis(V.itj %*% detParms)
    colP <- plogis(X.it.gam  %*% colParms)
    extP <- plogis(X.it.eps %*% extParms)

    detP <- array(detP, c(J, nY, M))
    colP <- matrix(colP, M, nY, byrow = TRUE)
    extP <- matrix(extP, M, nY, byrow = TRUE)

    ## create transition matrices (phi^T)
    phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
    for(i in 1:M) {
        for(t in 1:(nY-1)) {
            phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t],
                1-extP[i,t]))
            }
        }

    ## first compute latent probs
    x <- array(NA, c(2, nY, M))
    x[1,1,] <- 1-psiP
    x[2,1,] <- psiP
    for(i in 1:M) {
        for(t in 2:nY) {
            x[,t,i] <- (phis[,,t-1,i] %*% x[,t-1,i])
            }
        }

    z <- 0:1
    bup <- array(NA_real_, c(M, 2, nY))
    colnames(bup) <- z

    for(i in 1:M) {
        for(t in 1:nY) {
            g <- rep(1, 2)
            for(j in 1:J) {
                if(is.na(ya[i,j,t]) | is.na(detP[j,t,i]))
                    next
                g <- g * dbinom(ya[i,j,t], 1, z*detP[j,t,i])
            }
            tmp <- x[,t,i] * g
            bup[i,,t] <- tmp/sum(tmp)
        }
    }

    new("unmarkedRanef1", bup=bup)
})







setMethod("ranef", "unmarkedFitPCO",
    function(object, ...)
{
#    browser()
    dyn <- object@dynamics
    formlist <- object@formlist
    formula <- as.formula(paste(unlist(formlist), collapse=" "))
    D <- unmarked:::getDesign(object@data, formula)
    delta <- D$delta
    deltamax <- max(delta, na.rm=TRUE)

    lam <- predict(object, type="lambda")[,1] # Too slow, use D$Xlam instead
    om <- predict(object, type="omega")[,1]
    R <- length(lam)
    T <- object@data@numPrimary
    p <- getP(object)
    K <- object@K
    N <- 0:K
    y <- getY(getData(object))
    J <- ncol(y)/T
    if(dyn != "notrend") {
        gam <- predict(object, type="gamma")[,1]
        gam <- matrix(gam, R, T-1, byrow=TRUE)
    }
    om <- matrix(om, R, T-1, byrow=TRUE)
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    ya <- array(y, c(R, J, T))
    pa <- array(p, c(R, J, T))
    bup <- array(NA_real_, c(R, length(N), T))
    colnames(bup) <- N
    mix <- object@mixture
    if(dyn=="notrend")
        gam <- lam*(1-om)

    if(dyn %in% c("constant", "notrend")) {
        tp <- function(N0, N1, gam, om) {
            c <- 0:min(N0, N1)
            sum(dbinom(c, N0, om) * dpois(N1-c, gam))
        }
    } else if(dyn=="autoreg") {
        tp <- function(N0, N1, gam, om) {
            c <- 0:min(N0, N1)
            sum(dbinom(c, N0, om) * dpois(N1-c, gam*N0))
        }
    } else if(dyn=="trend") {
        tp <- function(N0, N1, gam, om) {
            dpois(N1, gam*N0)
        }
    }

    P <- matrix(NA_real_, K+1, K+1)

    for(i in 1:R) {
        switch(mix,
               P  = g2 <- dpois(N, lam[i]),
               NB = {
                   alpha <- exp(coef(object, type="alpha"))
                   g2 <- dnbinom(N, mu=lam[i], size=alpha)
               },
               ZIP = {
                   psi <- plogis(coef(object, type="psi"))
                   g2 <- (1-psi)*dpois(N, lam[i])
                   g2[1] <- psi + (1-psi)*exp(-lam[i])
               })
        g1 <- rep(1, K+1)
        for(j in 1:J) {
            if(is.na(ya[i,j,1]) | is.na(pa[i,j,1]))
                next
            g1 <- g1 * dbinom(ya[i,j,1], N, pa[i,j,1])
        }
        g1g2 <- g1*g2
        bup[i,,1] <- g1g2 / sum(g1g2)
        for(t in 2:T) {
            for(n0 in N) {
                for(n1 in N) {
                    P[n0+1, n1+1] <- tp(n0, n1, gam[i,t-1], om[i,t-1])
                }
            }
            delta.it <- delta[i,t-1]
            if(delta.it > 1) {
                P1 <- P
                for(d in 2:delta.it) {
                    P <- P %*% P1
                }
            }
            g1 <- rep(1, K+1)
            for(j in 1:J) {
                if(is.na(ya[i,j,t]) | is.na(pa[i,j,t]))
                    next
                g1 <- g1 * dbinom(ya[i,j,t], N, pa[i,j,t])
            }
            g <- colSums(P * bup[i,,t-1]) * g1
            bup[i,,t] <- g / sum(g)
        }
    }
    new("unmarkedRanef1", bup=bup)
})











setGeneric("postMode",
    function(object, ...) standardGeneric("postMode"))
setMethod("postMode", "unmarkedRanef1", function(object)
{
    bup <- object@bup
    N <- as.integer(colnames(bup))
    modes <- apply(bup, c(1,3), function(x) N[which.max(x)])
    modes <- drop(modes) # convert to vector if T=1
    return(modes)
})



setGeneric("postMean",
    function(object, ...) standardGeneric("postMean"))
setMethod("postMean", "unmarkedRanef1", function(object)
{
    bup <- object@bup
    dims <- dim(bup)
    N <- as.integer(colnames(bup))
    means <- apply(bup, c(1,3), function(x) sum(N*x))
    means <- drop(means) # convert to vector if T=1
    return(means)
})





setMethod("confint", "unmarkedRanef1", function(object, parm, level=0.95)
{
    if(!missing(parm))
        warning("parm argument is ignored. Did you mean to specify level?")
    bup <- object@bup
    N <- as.integer(colnames(bup))
    R <- nrow(bup)
    T <- dim(bup)[3]
    CI <- array(NA_real_, c(R,2,T))
    alpha <- 1-level
    c1 <- alpha/2
    c2 <- 1-c1
    colnames(CI) <- paste(c(c1,c2)*100, "%", sep="")
    for(i in 1:R) {
        for(t in 1:T) {
            pr <- bup[i,,t]
            ed <- cumsum(pr)
            lower <- N[which(ed >= c1)][1]
            upper <- N[which(ed >= c2)][1]
            CI[i,,t] <- c(lower, upper)
        }
    }
    CI <- drop(CI) # Convert to matrix if T==1
    return(CI)
})






setMethod("show", "unmarkedRanef1", function(object)
{
    bup <- object@bup
    dims <- dim(bup)
    T <- dims[3]
    if(T==1)
        print(cbind(#Mean=postMean(object),
                    Mode=postMode(object), confint(object)))
    else if(T>1) {
#        means <- postMean(object)
        modes <- postMode(object)
        CI <- confint(object)
        out <- array(NA_real_, c(dims[1], 3, T))
        dimnames(out) <- list(NULL,
                              c(#"Mean",
                                "Mode", "2.5%", "97.5%"),
                              paste("Year", 1:T, sep=""))
        for(t in 1:T) {
            out[,,t] <- cbind(#means[,t],
                              modes[,t], CI[,,t])
        }
        print(out)
    }
})



setAs("unmarkedRanef1", "array", function(from) {
    bup <- from@bup
    dims <- dim(bup)
    R <- dims[1]
    T <- dims[3]
    dimnames(bup) <- list(1:R, colnames(bup), 1:T)
    bup <- drop(bup)
    return(bup)
})


setAs("unmarkedRanef1", "data.frame", function(from) {
    bup <- from@bup
    dims <- dim(bup)
    R <- dims[1]
    lN <- dims[2]
    T <- dims[3]
    N <- as.integer(colnames(bup))
    N.ikt <- rep(rep(N, each=R), times=T)
    site <- rep(1:R, times=lN*T)
    year <- rep(1:T, each=R*lN)
    dat <- data.frame(site=site, year=year, N=N.ikt, p=as.vector(bup))
    dat <- dat[order(dat$site),]
    if(T==1)
        dat$year <- NULL
    return(dat)
})



setMethod("plot", c("unmarkedRanef1", "missing"), function(x, y, ...)
{
    bup <- x@bup
    T <- dim(bup)[3]
    N <- as.integer(colnames(bup))
    xlb <- ifelse(length(N)>2, "Abundance", "Occurrence")
    ylb <- "Probability"
    dat <- as(x, "data.frame")
    site.c <- as.character(dat$site)
    nc <- nchar(site.c)
    mc <- max(nc)
    dat$site.c <- paste("site", sapply(site.c, function(x)
         paste(paste(rep("0", mc-nchar(x)), collapse=""), x, sep="")),
         sep="")
    if(T==1)
        xyplot(p ~ N | site.c, dat, type="h", xlab=xlb, ylab=ylb, ...)
    else if(T>1) {
        year.c <- as.character(dat$year)
        nc <- nchar(year.c)
        mc <- max(nc)
        dat$year.c <- paste("year", sapply(year.c, function(x)
            paste(paste(rep("0", mc-nchar(x)), collapse=""), x, sep="")),
            sep="")
        xyplot(p ~ N | site.c+year.c, dat, type="h", xlab=xlb, ylab=ylb, ...)
    }
})












