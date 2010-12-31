

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
      tol = 1e-4)

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
      starts=c(1.6, 0.24, 1.16, -0.268))
  checkEqualsNumeric(coef(fm1), c(1.6722231, 0.2441238, 1.1623724, -0.2683158), 
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
      starts=c(1.4, -1.3, 1.8, -1.1, 0.7))
  checkEqualsNumeric(coef(fm2), c(1.4294109, -1.2762155, 1.8922928, -1.1444297, 
      0.7388486), tol = 1e-4)

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
      starts=c(1.5, 0, 1, 0))
  checkEqualsNumeric(coef(fm3), c(1.6931097, 0.2841317, 1.0706006, -0.2417096),
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
        c(1.69726936, -0.31631624, 1.88162678, -0.08689978), tol = 1e-5)
    
}
  




      