

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
