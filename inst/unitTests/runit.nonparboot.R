





test.nonparboot.colext <- function() {

    set.seed(343)
    nSites <- 10
    nReps <- 3
    nYears <- 2

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

    umf1 <- unmarkedMultFrame(y=y, siteCovs=data.frame(sc=sc),
                              yearlySiteCovs=list(ysc=ysc),
                              obsCovs=list(oc=oc), numPrimary=nYears)

    m1 <- colext(~1, ~1, ~1, ~1, umf1)
    m1 <- nonparboot(m1, B=2)

    checkEqualsNumeric(vcov(m1, method="nonparboot"),
                       matrix(c(
  1.5313035, -0.4177529,  0.8101477,  0.4844295,
 -0.4177529,  0.1139666, -0.2210153, -0.1321566,
  0.8101477, -0.2210153,  0.4286148,  0.2562911,
  0.4844295, -0.1321566,  0.2562911,  0.1532498), 4, byrow=TRUE),
                       tolerance=1e-6)

    m1.1 <- colext(~sc, ~ysc, ~1, ~oc, umf1)
    m1.1 <- nonparboot(m1.1, B=2)

    sc[1] <- NA
#    ysc[2,1] <- NA
    oc[3,4] <- NA

    umf2 <- unmarkedMultFrame(y=y, siteCovs=data.frame(sc=sc),
                              yearlySiteCovs=list(ysc=ysc),
                              obsCovs=list(oc=oc), numPrimary=nYears)

    m2 <- colext(~sc, ~ysc, ~1, ~oc, umf2)
    m2 <- nonparboot(m2, B=2)



}
