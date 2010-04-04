test.distsamp.fit.simple.1 <- function() {  
    y <- matrix(rep(1, 10), 5, 2)
    umf <- unmarkedFrameDS(y = y, dist.breaks=c(0, 5, 10)/1000, survey="line", 
        tlength=rep(1, 5), unitsIn="km")
    fm <- distsamp(~ 1 ~ 1, data = umf)

    lam <- fm['state']
    det <- fm['det']

    lam <- coef(backTransform(lam))
    checkEqualsNumeric(lam, 1, tolerance=1e-2)

    det <- coef(backTransform(det))
    checkEqualsNumeric(det, 0.05778171)

    bt <- backTransform(fm, type = 'state')
    checkEqualsNumeric(coef(bt), 1, tolerance=1e-2)

    bt <- backTransform(fm, type = 'det')
    checkEqualsNumeric(coef(bt), 0.05778171)
    }


test.distsamp.fit.covs <- function() {
    y <- matrix(rep(4:1, 10), 5, 2, byrow=TRUE)
    siteCovs <- data.frame(x = c(0, 2, 3, 4, 1))
    umf <- unmarkedFrameDS(y = y, siteCovs = siteCovs, 
        dist.breaks=c(0, 5, 10)/1000, survey="line", tlength=rep(1, 5), 
        unitsIn="km")
    fm <- distsamp(~ x ~ x, data = umf)
  
    lam <- fm['state']
    det <- fm['det']

    checkEqualsNumeric(coef(lam), c(1.4340999, -0.1102387), tolerance = 1e-4)
    checkEqualsNumeric(coef(det), c(-4.64686395, -0.09337832), tolerance = 1e-4)

    lam.lc <- linearComb(fm, type = 'state', c(1, 2))
    det.lc <- linearComb(fm, type = 'det', c(1, 2))
    
    checkEqualsNumeric(coef(lam.lc), 1.213623, tol = 1e-4)
    checkEqualsNumeric(coef(det.lc), -4.833621, tol = 1e-4)

    checkEqualsNumeric(coef(backTransform(lam.lc)), 3.365655, tol = 1e-4)
    checkEqualsNumeric(coef(backTransform(det.lc)), 0.007957658, tol = 1e-4)
    }

