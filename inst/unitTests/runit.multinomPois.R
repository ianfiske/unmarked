

test.multinomPois.na <- function() {

    y <- matrix(c(
        5, 3, 2,
        3, 3, 1,
        2, 0, 0,
        0, 0, 0), nrow=4, ncol=3, byrow=TRUE)
    
    sc <- data.frame(x1 = c(NA, 2, 3, 4))
    oc <- list(x2 = matrix(c(
        1, 1, 1,
        3, NA, 1,
        0, 0, 1,
        NA, NA, NA), nrow=4, ncol=3, byrow=TRUE)) 
    
    umf1 <- unmarkedFrameMPois(y = y, siteCovs = sc, obsCovs = oc, 
        type="removal") 


    checkEquals(obsToY(umf1), diag(3))

    m1 <- multinomPois(~1 ~1, umf1)
    checkEqualsNumeric(coef(m1), c(1.7490576,  -0.2330901), tol=1e-5)


    m2 <- multinomPois(~x2 ~1, umf1)
    checkEqualsNumeric(coef(m2), c(1.9159845, 0.2248897, -0.1808144), tol=1e-5)
    checkEquals(m2@sitesRemoved, 4)
    
    
    # The old mapping matrix    
    o2y <- diag(3)
    o2y[upper.tri(o2y)] <- 1
    umf1@obsToY <- o2y
    
    m3 <- multinomPois(~x2 ~1, umf1)
    checkEqualsNumeric(m3@sitesRemoved, 4)
    

    
    }


      