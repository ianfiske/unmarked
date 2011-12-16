
test.colext <- function()
{

    nsites <- 6
    nyr <- 4
    nrep <- 2
    y <- matrix(c(
        1,0, 1,1, 0,0, 0,0,
        1,1, 0,0, 0,0, 0,0,
        0,0, 0,0, 0,0, 0,0,
        0,0, 1,1, 0,0, 0,0,
        1,1, 1,0, 0,1, 0,0,
        0,0, 0,0, 0,0, 1,1), nrow=nsites, ncol=nyr*nrep, byrow=TRUE)

    umf1 <- unmarkedMultFrame(y=y, numPrimary=4)
    fm1 <- colext(~1, ~1, ~1, ~1, umf1)
    checkEqualsNumeric(coef(fm1),
        c(0.1422577, -1.4950576,  0.2100365,  1.1998444),
        tol=1e-6)

    set.seed(232)
    oc <- matrix(rnorm(length(y)), nrow(y), ncol(y))
    umf2 <- unmarkedMultFrame(y=y, obsCovs=list(oc=oc), numPrimary=nyr)
    fm2 <- colext(~1, ~1, ~1, ~oc, umf2, starts=c(coef(fm1), 0))
    checkEqualsNumeric(coef(fm2),
        c(0.1548709, -1.4697222, 0.1745153, 1.1154250, -0.3557752),
        tol=1e-6)

    y1 <- y
    y1[1,3] <- NA

    umf3 <- unmarkedMultFrame(y=y1, obsCovs=list(oc=oc), numPrimary=nyr)
    fm3 <- colext(~1, ~1, ~1, ~1, umf3, starts=coef(fm1))
    checkEqualsNumeric(coef(fm3),
        c(0.2058462, -1.5612409, 0.4320085, 0.9616805),
        tol=1e-6)

    oc1 <- oc
    oc1[is.na(y1)] <- NA
    umf4 <- unmarkedMultFrame(y=y1, obsCovs=list(oc=oc1), numPrimary=nyr)
    fm4 <- colext(~1, ~1, ~1, ~oc, umf4, starts=coef(fm2))
    checkEqualsNumeric(coef(fm4),
        c(0.2783699, -1.5034240, 0.3470471, 0.7425092, -0.5110785),
        tol=1e-6)

    y2 <- y
    y2[4,] <- NA

    umf5 <- unmarkedMultFrame(y=y2, numPrimary=nyr)
    fm5 <- colext(~1, ~1, ~1, ~1, umf5, starts=coef(fm1))
    checkEqualsNumeric(coef(fm5),
        c(0.50002469, -1.99947927, -0.03660814, 1.09667556),
        tol=1e-6)
    checkEqualsNumeric(fm5@sitesRemoved, 4)

    ysc <- matrix(1:nyr, nsites, nyr, byrow=TRUE)
    ysc[1,1] <- NA
    umf6 <- unmarkedMultFrame(y=y, yearlySiteCovs=list(ysc=ysc),
                              numPrimary=nyr)
    checkException(fm3.1 <- colext(~1, ~1, ~ysc, ~1, umf6))
    checkException(fm3.2 <- colext(~1, ~ysc, ~1, ~1, umf6))

    ysc <- matrix(1:3, nsites, nyr, byrow=TRUE)
    ysc[1,1] <- NA
    y4 <- y
    y4[1,1:2] <- NA
    ysc4 <- ysc
    ysc4[1,1] <- 1 # NA
    umf7 <- unmarkedMultFrame(y=y4, yearlySiteCovs=list(ysc=ysc4),
                              numPrimary=nyr)
    fm7.1 <- colext(~1, ~1, ~ysc, ~1, umf7)
    fm7.2 <- colext(~1, ~ysc, ~1, ~1, umf7)


}
