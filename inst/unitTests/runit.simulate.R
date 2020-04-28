test.simulate.GDS <- function(){

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
    
    covs$par1[2] <- NA
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm <- gdistsamp(~par1, ~1, ~1, umf, se=FALSE, engine="C")
    
    #This used to error due to rmultinom not accepting size=NA
    s <- simulate(fm, nsim=2, na.rm=FALSE)
    checkEqualsNumeric(length(s), 2)
    checkEqualsNumeric(dim(s[[1]]), c(30,15))
    checkTrue(!any(is.na(s[[1]][1,])))
    checkTrue(all(is.na(s[[1]][2,])))

    pb <- parboot(fm, nsim=3)
    checkTrue(inherits(pb, "parboot"))

}
