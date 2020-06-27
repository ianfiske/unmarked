test.gdistsamp.covs <- function() {
    set.seed(343)
    R <- 30 # number of transects
    T <- 3  # number of replicates
    strip.width <- 50
    transect.length <- 100
    breaks <- seq(0, 50, by=10)

    lambda <- 5 # Abundance
    phi <- 0.6  # Availability
    sigma <- 30 # Half-normal shape parameter

    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    #Check that error thrown when length(tlength)!=nrow(y)
    checkException(unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, (R-1)), numPrimary=T))

    # Organize data
    umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    # Fit the model
    fm1 <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE)

    checkEqualsNumeric(coef(fm1), c( 1.71894803, -0.03744387, 3.54452329))

    re1 <- ranef(fm1)
    checkEqualsNumeric(bup(re1, "mode")[1:7], c(3,5,3,5,5,2,5))

}

test.gdistsamp.halfnorm <- function(){
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    #When output = density
    fm_R <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE, engine="R")
    fm_C <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE, engine="C")

    checkEqualsNumeric(coef(fm_R),c(0.6584008,-0.4440278,3.4817094),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #When output = abundance
    fm_R <- gdistsamp(~1, ~1, ~1, umf, output="abund", se=FALSE, engine="R")
    fm_C <- gdistsamp(~1, ~1, ~1, umf, output="abund", se=FALSE, engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.35067,-0.44262,3.48149),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #With covariates
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE, engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.24510,0.54419,-1.28146,-0.109737,
                                    3.46295,-0.13228),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #Predict


    #Negative binomial
    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      mixture="NB", engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      mixture="NB", engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.41241,0.52442,-1.49024,-0.10546,
                                    3.46284,-0.129831,3.16892),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #With missing values
    yna <- y
    yna[1,c(1,6)] <- NA
    umf <- unmarkedFrameGDS(y = yna, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
                      se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
                      se=FALSE, engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.35065,0.52558,-1.39758,-0.10675,
                                    3.46283,-0.136344),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #With an entire session missing
    yna <- y
    yna[1,1:J] <- NA
    umf <- unmarkedFrameGDS(y = yna, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
                      se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
                      se=FALSE, engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.30815,0.53527,-1.35387,-0.11038,
                                    3.46293,-0.13458),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point
    set.seed(123)
    data(issj)
    covs <- issj[,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      engine="R", se=F)
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      engine="C", se=F)

    checkEqualsNumeric(coef(fm_R),c(-2.42178,-0.17874,4.38373,0.62425),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)
}

test.gdistsamp.uniform <- function(){
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_R <- gdistsamp(~par1, ~par2, ~1, umf, output="density",
                      keyfun="uniform", se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~1, umf, output="density",
                      keyfun="uniform", se=FALSE, engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.17120,0.54748,-1.60963,-0.13009),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point: doesn't work with this dataset, find another one?
    #OR maybe uniform just doesn't work with point?
    set.seed(123)
    data(issj)
    covs <- issj[,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    checkException(gdistsamp(~1, ~1, ~1, jayumf, output='density',
                      keyfun='uniform', engine="R", se=F))
    checkException(gdistsamp(~elevation, ~1, ~1, jayumf, output='density',
                      keyfun='uniform', engine="C", se=F))
}


test.gdistsamp.exp <- function(){
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    #With covariates
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="exp",engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="exp",engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.28243,0.54312,-1.16608,-0.101122,
                                    3.86666,-0.2492846),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point
    set.seed(123)
    data(issj)
    covs <- issj[,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="exp",engine="R", se=F)
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="exp",engine="C", se=F)

    checkEqualsNumeric(coef(fm_R),c(-1.531876,-0.2037537,3.870335,
                                    0.89988),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)
}

test.gdistsamp.hazard <- function(){
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    #With covariates
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="hazard",engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="hazard",engine="C")

    checkEqualsNumeric(coef(fm_R),c(1.29425,0.54233,-1.41658,-0.09267,
                                    3.45436,-0.19978,0.8270215),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point
    set.seed(123)
    data(issj)
    covs <- issj[,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="hazard",engine="R", se=F)
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="hazard",engine="C", se=F)

    checkEqualsNumeric(coef(fm_R),c(0.70099,-0.23473,1.38888,
                                    1.40786,0.44896),tol=1e-4)
    checkEqualsNumeric(coef(fm_R),coef(fm_C),tol=1e-3)
}

test.gdistsamp.predict <- function(){
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    ysc <- data.frame(fac_cov = factor(rep(letters[1:T], R)))

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=ysc,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_C <- gdistsamp(~1, ~fac_cov -1, ~1, umf,
                      keyfun="halfnorm", se=TRUE, engine="C")

    #lambda
    pr <- predict(fm_C, "lambda")
    checkTrue(inherits(pr, "data.frame"))
    checkEqualsNumeric(dim(pr), c(30,4))
    checkEqualsNumeric(pr[1,1], 3.767935, tol=1e-5)
    nd <- data.frame(par1=0)
    pr2 <- predict(fm_C, type='lambda', newdata=nd)
    checkTrue(inherits(pr2, "data.frame"))
    checkEqualsNumeric(dim(pr2), c(1,4))
    checkEqualsNumeric(pr2[1,1], 3.767935, tol=1e-5)

    #phi
    pr <- predict(fm_C, "phi")
    checkTrue(inherits(pr, "data.frame"))
    checkEqualsNumeric(dim(pr), c(90,4))
    checkEqualsNumeric(pr[1,1], 0.4461197, tol=1e-5)

    nd <- data.frame(fac_cov=factor(letters[1:3]))
    pr2 <- predict(fm_C, type="phi", newdata=nd)
    checkTrue(inherits(pr2, "data.frame"))
    checkEqualsNumeric(dim(pr2), c(3,4))
    checkEqualsNumeric(pr2[1,1], 0.4461197, tol=1e-5)

    #sigma
    pr <- predict(fm_C, "det")
    checkTrue(inherits(pr, "data.frame"))
    checkEqualsNumeric(dim(pr), c(90,4))
    checkEqualsNumeric(pr[1,1], 32.51537, tol=1e-5)

    nd <- data.frame(par3=0)
    pr2 <- predict(fm_C, type='det', newdata=nd)
    checkTrue(inherits(pr2, "data.frame"))
    checkEqualsNumeric(dim(pr2), c(1,4))
    checkEqualsNumeric(pr2[1,1], 32.51537, tol=1e-5)
}
