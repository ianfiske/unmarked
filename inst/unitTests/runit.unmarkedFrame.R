test.emptyframe <- function() {
  checkException(umf <- unmarkedFrame())
}

test.frame <- function() {
  M <- 10
  J <- 3
  y <- matrix(rbinom(J * M, 1, 0.5), M, J)
  siteCovs <- data.frame(a = rnorm(M), b = factor(gl(2,5)))
  umf <- unmarkedFrame(y = y, siteCovs = siteCovs)
}
