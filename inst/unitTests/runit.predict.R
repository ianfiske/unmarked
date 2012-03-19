
test.occu <- function() {
    set.seed(55)
    R <- 20
    J <- 4
    x1 <- rnorm(R)
    x2 <- factor(c(rep('A', R/2), rep('B', R/2)))
    x3 <- matrix(rnorm(R*J), R, J)
    z <- rbinom(R, 1, 0.5)
    y <- matrix(rbinom(R*J, 1, z*0.6), R, J)
    x1[1] <- NA
    x3[2,1] <- NA
    x3[3,] <- NA
    umf1 <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x1=x1, x2=x2),
                              obsCovs=list(x3=x3))
    fm1 <- occu(~x3 ~x1+x2, umf1)
    E1.1 <- predict(fm1, type="state")
    E1.2 <- predict(fm1, type="det")

    nd1.1 <- data.frame(x1=0, x2=factor('A', levels=c('A','B')))
    nd1.2 <- data.frame(x3=0)
    E1.3 <- predict(fm1, type="state", newdata=nd1.1, appendData=TRUE)
    E1.4 <- predict(fm1, type="det", newdata=nd1.2)

    r1 <- raster(matrix(rnorm(100), 10))
    checkException(predict(fm1, type="state", newdata=r1))
    s1 <- stack(r1)
    checkException(predict(fm1, type="state", newdata=s1))
    layerNames(s1) <- c("x3")
    E1.5 <- predict(fm1, type="det", newdata=s1)
    E1.5 <- predict(fm1, type="det", newdata=s1, appendData=TRUE)


}
