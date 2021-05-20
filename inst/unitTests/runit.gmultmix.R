



test.gmultmix.fit <- function()
{
    y <- matrix(0:3, 5, 4)
    siteCovs <- data.frame(x = c(0,2,3,4,1))
    siteCovs[3,1] <- NA
    obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
    yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C")

    checkEquals(fm_R@sitesRemoved, 3)
    coef_truth <- c(2.50638554, 0.06226627, 0.21787839, 6.46029769, -1.51885928,
            -0.03409375, 0.43424295)
    checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
    checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)

    obsCovs[10,2] <- NA
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C")

    checkEquals(fm_R@sitesRemoved, 3)
    checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
    checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)

    yrSiteCovs[2, 1] <- NA
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C")

    coef_truth <- c(1.17280104, 0.37694710, 2.38249795, 2.87354955, -0.83875134,
            -0.08446507, 1.88056826)
    checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
    checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)

    #Negative binomial
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, mixture="NB", K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, mixture="NB", K=23, engine="C")
    checkEqualsNumeric(coef(fm_R), coef(fm_C), tol=1e-5)

    #Check error when random effect in formula
    checkException(gmultmix(~(1|dummy), ~1, ~1, umf))
}



test.gmultmix.na <- function()
{
    y <- matrix(0:3, 8, 9)
    oc <- matrix(c(
        0,1,1,    0,0,0,    1,2,1,
        0,0,1,    0,1,2,    2,0,3,
        2,0,0,    0,1,0,    0,0,NA,
        2,0,NA,   0,0,NA,   0,0,NA,
        NA,0,0,   NA,0,0,   NA,0,0,
        1,NA,3,   0,NA,0,   0,NA,0,
        NA,NA,NA, 1,2,0,    0,0,0,
        NA,NA,NA, NA,NA,NA, NA,NA,NA), byrow=TRUE, nrow=8, ncol=9)
    o2y <- diag(3)
    o2y[upper.tri(o2y)] <- 1
    m <- matrix(0, 3, 3)

    kronecker(diag(3), o2y)

    o2y <- rbind(
        cbind(o2y, m, m),
        cbind(m, o2y, m),
        cbind(m, m, o2y))

    oc.na <- is.na(oc)
    oc.na %*% o2y

    kronecker(diag(3), matrix(1, 3, 2))


    }



