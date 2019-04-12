test.pifun.gmultmix <- function(){

    set.seed(123)
    y <- matrix(0:3, 5, 4)
    siteCovs <- data.frame(x = c(0,2,3,4,1))
    siteCovs[3,1] <- NA
    obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
    yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

    #bad type
    checkException(unmarkedFrameGMM(y = y, siteCovs = siteCovs, 
                                           obsCovs = obsCovs, 
                                           yearlySiteCovs = yrSiteCovs, 
                                           type="fake", numPrimary=2))
    
    #type = "removal"
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs, 
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C")
    coef_truth <- c(2.50638554, 0.06226627, 0.21787839, 6.46029769, -1.51885928, 
            -0.03409375, 0.43424295) 
    checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5) 
    checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)
 
    #type = "double"
    set.seed(123)
    y <- matrix(0:3, 5, 6)
    siteCovs <- data.frame(x = c(0,2,3,4,1))
    siteCovs[3,1] <- NA
    obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/30)
    yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))
    
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs, 
        yearlySiteCovs = yrSiteCovs, type="double", numPrimary=2)
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C")
    coef_truth <- c(2.42178997506, 0.0790163156156549, -0.483051837753635, 
                    0.0869555282843506,0.748689310997745, -0.0458486142792914, 
                    -0.374467809850705)
    checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5) 
    checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)

    #type = "depDouble"
    set.seed(123)
    y <- matrix(0:3, 5, 4)
    siteCovs <- data.frame(x = c(0,2,3,4,1))
    siteCovs[3,1] <- NA
    obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
    yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))
    
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs, 
        yearlySiteCovs = yrSiteCovs, type="depDouble", numPrimary=2)
    
    fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
    fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C")
    coef_truth <- c(2.50638554076642, 0.0622662713385639, 0.217878385786452, 
                    6.46029769667278,-1.51885927831977, -0.0340937462852327, 
                    0.434242947016337) 
    checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5) 
    checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)
}

test.pifun.multinomPois <- function(){

     nSites <- 50
     lambda <- 10
     p1 <- 0.5
     p2 <- 0.3
     cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
     set.seed(9023)
     N <- rpois(nSites, lambda)
     y <- matrix(NA, nSites, 3)
     for(i in 1:nSites) {
       y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
     }
     
     # incorrect type
     observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
     #TODO: This should return a more informative error than it does
     checkException(unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="fake"))
    
     #Type "double"
     umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="double")
     fm_R <- multinomPois(~observer-1 ~1, umf, engine="R") 
     fm_C <- multinomPois(~observer-1 ~1, umf, engine="C")
     coef_truth <- c(2.2586622, 0.1739752, -0.5685933)
     checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
     checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)
      
     #Type "depDouble" 
     nSites <- 50
     lambda <- 10
     p1 <- 0.5
     p2 <- 0.3
     cp <- c(p1, p2*(1-p1))
     set.seed(9023)
     N <- rpois(nSites, lambda)
     y <- matrix(NA, nSites, 2)
     for(i in 1:nSites) {
       y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:2]
     } 
     # Fit model
     observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
     umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="depDouble")
     fm_R <- multinomPois(~observer-1 ~1, umf, engine="R")
     fm_C <- multinomPois(~observer-1 ~1, umf, engine="C")
     coef_truth <- c(2.0416086, 0.7430343, 0.4564236)
     checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
     checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)

     #Type "removal"
     umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="removal")
     fm_R <- multinomPois(~observer-1 ~1, umf, engine="R")
     fm_C <- multinomPois(~observer-1 ~1, umf, engine="C")
     coef_truth <- c(2.0416086, 0.7430343, 0.4564236)
     checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
     checkEqualsNumeric(coef(fm_C), coef_truth, tol = 1e-5)
}
