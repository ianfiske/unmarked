

test.pcountOpen.null <- function() 
{
  y <- matrix(c(
      3, 2, 1, 4,
      3, 4, 2, 1,
      0, 1, 2, 3,
      5, 3, 3, 4,
      2, 4, 3, 3), 5, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:20)
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10, 
      starts=c(1, 0, 0, 7))
  checkEqualsNumeric(coef(fm1), c(0.9565311, 0.2741022, 0.1352888, 7.0041290), 
      tol = 1e-5)

  fm2 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="gamma", K=10)
  checkEqualsNumeric(coef(fm2), c(1.8219364, 8.7430266, -0.2873533), 
      tol = 1e-5)

  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="omega", K=10)
  checkEqualsNumeric(coef(fm3), c(1.6937578, -1.4762351, -0.1649877), 
      tol = 1e-5)
}






test.pcountOpen.na <- function() 
{
  y1 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA, 
      2, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf1 <- unmarkedFramePCO(y = y1, siteCovs = siteCovs, obsCovs = obsCovs)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf1, se=FALSE, K=10, 
      starts=c(1, 0, 0, 7))
  checkEqualsNumeric(coef(fm1), c(0.9093738, 0.3925347, -0.4829625, 6.9977870), 
      tol = 1e-5)

  y2 <- matrix(c(
      1, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1, 
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  oc <- y1 + -2:3
      
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- list(o1 = oc)
  umf2 <- unmarkedFramePCO(y = y2, siteCovs = siteCovs, obsCovs = obsCovs)

  fm2 <- pcountOpen(~1, ~1, ~1, ~o1, data = umf2, se=FALSE, K=10, 
      starts=c(1.5, 1, -6, -1, 1))
  checkEqualsNumeric(coef(fm2), c(1.3000908, -9.5548411, 2.3858981, -0.6812891, 
      1.1473823), tol = 1e-5)

  y3 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA, 
      NA, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf3 <- unmarkedFramePCO(y = y3, siteCovs = siteCovs, obsCovs = obsCovs)
  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf3, se=FALSE, K=10, 
      starts=c(1.5, 1, -6, -1))
  checkEqualsNumeric(coef(fm3), c(0.9910793, 0.3977758, -0.4453013, 3.6724174),
      tol = 1e-5)
  checkEquals(fm3@sitesRemoved, 6)
  

}







test.pcountOpen.delta <- function()
{
    M <- 5
    T <- 4
    y <- matrix(c(
        NA, 2, 1, 4,
        3, NA, 2, 1,
        0, 1, NA, 3,
        5, 3, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    dates <- matrix(c(1,3,5,7), M, T, byrow=TRUE)
    delta <- unmarked:::formatDelta(dates, y)
    ans <- matrix(c(
        1, 2, 2, 2,
        1, 2, 4, 2,
        1, 2, 2, 4,
        1, 2, 2, 2,
        1, 2, 2, 2), M, T, byrow=TRUE)

    checkEquals(delta, ans) 
    
    dates2 <- matrix(c(
      2, 4, 6, 8, 
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8, 
      2, 4, 6, 8), M, T, byrow=TRUE)
    delta2 <- unmarked:::formatDelta(dates2, y)
    ans2 <- matrix(c(
        2, 3, 2, 2,
        1, 3, 5, 2,
        2, 2, 2, 4,
        1, 3, 2, 2,
        2, 3, 2, 2), M, T, byrow=TRUE)

    checkEquals(delta2, ans2) 
    
    dates3 <- matrix(c(
      2, NA, 6, 8, 
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8, 
      2, 4, 6, 8), M, T, byrow=TRUE)
    checkException(unmarkedFramePCO(y=y, dates=dates3))

    dates4 <- dates2
    dates4[is.na(y)] <- NA
    mode(dates4) <- "integer"
    delta4 <- unmarked:::formatDelta(dates4, y)
    umf <- unmarkedFramePCO(y=y, dates=dates4)
    fm <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10, starts=c(1.2, 0, 1.4, 1.2))
    checkEqualsNumeric(coef(fm), 
        c(1.19052243, -0.3701578, 1.3852273, 1.3713127), tol = 1e-5)
    
}
  





test.pcountOpen.tranProbs <- function()
{

    tranProbsR <- function(K, omega, gamma, delta, dynamics) {
        lk <- length(K)
        bpsum <- matrix(NA, lk, lk)
        for(i in 1:lk) {
        for(j in 1:lk) {
            cmin0 <- 0:min(K[i], K[j])
            if(identical(dynamics, "autoreg"))
                gamma2 <- gamma * K[j]
            else
                gamma2 <- gamma
            bpsum[i, j] <- sum(dbinom(cmin0, K[j], omega) * 
                dpois(K[i]-cmin0, gamma2))
            }}
        if(delta>1) {
            for(i in 2:delta) {
                bpsum <- bpsum %*% bpsum
                }
            }
        return(bpsum)
        }

    tp1 <- unmarked:::tranProbs(Kr=0:4, omegaR=0.5, gammaR=1, deltaR=1, 
        dynamicsR="constant")
    tp1r <- tranProbsR(K=0:4, omega=0.5, gamma=1, delta=1, 
        dynamics="constant")
    checkEquals(tp1, tp1r)
    
    tp2 <- unmarked:::tranProbs(Kr=0:4, omegaR=0.5, gammaR=1, deltaR=3, 
        dynamicsR="constant")
    tp2r <- tranProbsR(K=0:4, omega=0.5, gamma=1, delta=3, 
        dynamics="constant")
    checkEquals(tp2, tp2r)
    
    tp3 <- unmarked:::tranProbs(Kr=0:4, omegaR=0.5, gammaR=1, deltaR=1, 
        dynamicsR="autoreg")
    tp3r <- tranProbsR(K=0:4, omega=0.5, gamma=1, delta=1, 
        dynamics="autoreg")
    checkEquals(tp3, tp3r)
    
    tp4 <- unmarked:::tranProbs(Kr=0:4, omegaR=0.5, gammaR=1, deltaR=5, 
        dynamicsR="autoreg")
    tp4r <- tranProbsR(K=0:4, omega=0.5, gamma=1, delta=5, 
        dynamics="autoreg")
    checkEquals(tp4, tp4r)
    
    
    
}    
    
        
        
        