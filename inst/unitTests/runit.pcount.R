test.pcount.offest <- function() {

  y <- matrix(rep(8,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10)
  umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- pcount(~ o1 ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm), structure(c(-0.76302, 7.94584, 2.99753),
                                         .Names = c("lam(Int)", "p(Int)", "p(o1)")), tol = 1e-5)

}
