
library(unmarked)


# Test null model

set.seed(26)

n <- 100  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals availabe for detection
    
for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)
    Nleft1 <- N[i,] - y[i,,1]
    y[i,,2] <- rbinom(T, Nleft1, p)
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }
    
y.MxTJ <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])


umf1 <- unmarkedFrameGMM(y=y.MxTJ, numPrimary=T, type="removal")


system.time(m1 <- gmm(~1, ~1, ~1, data=umf1)) #2.9


SSE(m1)

(pb1 <- parboot(m1, nsim=20, report=5))
plot(pb1)