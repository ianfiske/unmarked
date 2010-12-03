

sim <- function(lambda=5, phi=0.5, sigma=20, R=100, T=3, radius=50, 
    breaks=seq(0, 50, by=10))
{    
    J <- length(breaks)-1
    A <- (2*radius)^2 / 10000     # Area (ha) of square containing circle
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda * A) # Individuals within the square
        N <- rbinom(T, M, phi)    # Individuals available at time t
        for(t in 1:T) {
            # coordinates of each individual
            xy <- cbind(x=runif(N[t], -radius, radius), 
                y=runif(N[t], -radius, radius))
        
            # Distances from point
            d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
            d <- d[d <= radius]

            # Detection process
            if(length(d)) {
                p <- exp(-d^2 / (2 * sigma^2)) # half-normal
                d <- d[rbinom(length(d), 1, p) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
                }
            }
        }
    y <- matrix(y, nrow=R)
    return(y)
}

breaks <- seq(0, 50, by=10)
T <- 3
umf <- unmarkedFrameGDS(y = sim(T=T, breaks=breaks), survey="point", 
    unitsIn="m", dist.breaks=breaks, numPrimary=T)
summary(umf)
    
m1 <- gdistsamp(~1, ~1, ~1, umf)
backTransform(m1, type="lambda")
backTransform(m1, type="phi")
backTransform(m1, type="det")


nsim <- 50
simout <- matrix(NA, nsim, 3)
colnames(simout) <- c('lambda', 'phi', 'p')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 3
    y <- sim(T=T, breaks=breaks)
    umf <- unmarkedFrameGDS(y = y, survey="point", 
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, rel.tol=1e-3)
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2:3]))
    }
    
hist(simout[,1]); abline(v=5, col=4)    
    
    
    

