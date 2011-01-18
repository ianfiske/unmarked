test.colext.getP <- function() {
    y <- matrix(c(
        0,0, 1,1, 0,1,
        1,1, 0,0, 0,0,
        0,0, 0,0, 0,0), nrow=3, ncol=6, byrow=TRUE)
    umf <- unmarkedMultFrame(y=y, numPrimary=3)
    fm <- colext(~1, ~1, ~1, ~1, umf)
    checkException(p <- getP(fm))
}

