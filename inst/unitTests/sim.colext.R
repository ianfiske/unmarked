
library(unmarked)
library(RUnit)


# ----------------------------- simulate --------------------------------------


sim <- function(nSites=100, nReps=5, nYears=5, psi=0.5, gamma=0.2, epsilon=0.8, 
    p=0.4)
{

    y <- array(NA, c(nSites, nReps, nYears))
    Z <- matrix(NA, nSites, nYears)

    phi <- 1-epsilon

    Z[,1] <- rbinom(nSites, 1, psi)
    for(t in 2:nYears) {
        muZ <- Z[,t-1] * phi + (1 - Z[,t-1]) * gamma
        Z[,t] <- rbinom(n, 1, muZ)
        }
    for(j in 1:nReps)
        for(t in 1:nYears)
            y[,j,t] <- rbinom(n, 1, Z[,t]*p)
    
    y <- matrix(y, nSites, nReps*nYears)
    return(y)
}

# sim()


# ------------------------------- unmarked ------------------------------------




set.seed(3)
sim1 <- sim()
umf <- unmarkedMultFrame(y = sim1, numPrimary = nYears)

(m <- colext(~1, ~1, ~1, ~1, umf, control=list(trace=T, REPORT=1)))

backTransform(m, type="psi")
backTransform(m, type="col")
backTransform(m, type="ext")
backTransform(m, type="det")

checkEqualsNumeric(coef(m), c(-0.1047112, -1.3000613, 1.5203993, -0.3634747),
    tol=1e-5)












