





test.nonparboot.colext <- function() {

    nSites <- 20
    nReps <- 3
    nYears <- 5

    epsilon <- 0.7
    psi <- 0.5
    gamma <- 0.4
    p <- 0.6

    y <- array(NA, c(nSites, nReps, nYears))
    Z <- matrix(NA, nSites, nYears)

    phi <- 1-epsilon

    Z[,1] <- rbinom(nSites, 1, psi)
    for(t in 2:nYears) {
        muZ <- Z[,t-1] * phi + (1 - Z[,t-1]) * gamma
        Z[,t] <- rbinom(nSites, 1, muZ)
    }
    for(j in 1:nReps)
        for(t in 1:nYears)
            y[,j,t] <- rbinom(nSites, 1, Z[,t]*p)

    y <- matrix(y, nSites, nReps*nYears)

    sc <- rnorm(nSites)
    ysc <- matrix(rnorm(nSites*nYears), nSites)
    oc <- matrix(rnorm(nSites*nReps*nYears), nSites)

    umf <- unmarkedMultFrame(y=y, siteCovs=data.frame(sc=sc),
                             yearlySiteCovs=list(ysc=ysc),
                             obsCovs=list(oc=oc), numPrimary=nYears)

    m1 <- colext(~1, ~1, ~1, ~1, umf)
    m1 <- nonparboot(m1, B=2)

    m2 <- colext(~sc, ~ysc, ~1, ~oc, umf)
    m2 <- nonparboot(m2, B=2)


}
