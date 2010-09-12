



test.occu.fit.NA <- function() 
{
    y <- matrix(rep(0:1,10),5,2)
    siteCovs <- data.frame(x = c(0,2,3,4,1))
    siteCovs[3,1] <- NA
    obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
    umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
    fm <- occu(~ o1 + o2 ~ x, data = umf)
    checkEquals(fm@sitesRemoved, 3)
    checkEqualsNumeric(coef(fm), 
        c(8.70123, 4.58255, 0.66243, -0.22862, 0.58192), 
        tol = 1e-5)

    obsCovs[10,2] <- NA
    umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
    fm <- occu(~ o1 + o2 ~ x, data = umf)
    checkEquals(fm@sitesRemoved, 3)
    checkEqualsNumeric(coef(fm), 
        c(8.91289, 1.89291, -1.42471, 0.67011, -8.44608), 
        tol = 1e-5)

}
