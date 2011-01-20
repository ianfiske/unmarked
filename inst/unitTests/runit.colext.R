test.colext.getP <- function() 
{
    y <- matrix(c(
        0,0, 1,1, 0,1,
        1,1, 0,0, 0,0,
        0,0, 0,0, 0,0), nrow=3, ncol=6, byrow=TRUE)
    umf <- unmarkedMultFrame(y=y, numPrimary=3)
    fm <- colext(~1, ~1, ~1, ~1, umf)
    checkException(p <- getP(fm))
    
}




test.colext.na <- function() 
{

    y <- matrix(c(
        0,0, 1,1, 0,1,
        1,1, 0,0, 0,0,
        0,0, 0,0, 0,0), nrow=3, ncol=6, byrow=TRUE)
    
    y2 <- y
    y2[1,3] <- NA
    y2[2,5:6] <- NA
    
    umf2 <- unmarkedMultFrame(y=y2, numPrimary=3)
    fm2 <- colext(~1, ~1, ~1, ~1, umf2)
    
    
    oc <- y2 + -1:1
    umf3 <- unmarkedMultFrame(y=y, obsCovs=list(oc=oc), numPrimary=3)
    fm3 <- colext(~1, ~1, ~1, ~oc, umf3)


    ysc <- matrix(1:3, 3, 3, byrow=TRUE)
    ysc[1,1] <- NA    
    umf4 <- unmarkedMultFrame(y=y, yearlySiteCovs=list(ysc=ysc), numPrimary=3)
    fm4.1 <- colext(~1, ~1, ~ysc, ~1, umf4)
    fm4.2 <- colext(~1, ~ysc, ~1, ~1, umf4)

    ysc <- matrix(1:3, 3, 3, byrow=TRUE)
    ysc[1,1] <- NA
    y5 <- y
    y5[1,1:2] <- NA
    ysc2 <- ysc
    ysc2[1,1] <- 1    
    umf5 <- unmarkedMultFrame(y=y5, yearlySiteCovs=list(ysc=ysc2), numPrimary=3)
    fm5.1 <- colext(~1, ~1, ~ysc, ~1, umf5)
    fm5.2 <- colext(~1, ~ysc, ~1, ~1, umf5)




}