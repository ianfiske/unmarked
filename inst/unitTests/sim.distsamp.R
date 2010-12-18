



# Simulate with half-normal detection function

simDSpt <- function(lambda=5, sigma=20, npts=100, radius=50,
    breaks=seq(0, 50, by=10))
{      
    A <- (2*radius)^2 / 10000 # Area (ha) of square containing circle
    y <- matrix(0, npts, length(breaks)-1)
    for(i in 1:npts) {
        M <- rpois(1, lambda * A) # Individuals within the square
        # coordinates of each individual
        xy <- cbind(x=runif(M, -radius, radius), y=runif(M, -radius, radius))
        
        # Distances from each point
        d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
        d <- d[d <= radius]

        # Detection process
        if(length(d)) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    return(y)
}

colSums(simDSpt())

set.seed(3)
umf1 <- unmarkedFrameDS(y = simDSpt(), survey="point", 
    dist.breaks=seq(0, 50, by=10), unitsIn="m")
(m1 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20))))
m2 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)), output="abund")


cxcheckEqualsNumeric(coef(m1), c(1.813819, 2.893771), tol=1e-5)
checkEquals(exp(coef(m1, type="state")), 
    exp(coef(m2, type="state")) / (pi * 50^2 / 10000), tol=0.01)


nsims <- 100
simout1 <- matrix(NA, nsims, 2)
lam <- 20
sig <- 30
for(i in 1:nsims) {
    cat("sim", i, "\n"); flush.console()
    umf <- unmarkedFrameDS(y = simDSpt(lambda=lam, sigma=sig), survey="point", 
        dist.breaks=seq(0, 50, by=10), unitsIn="m")
    m <- distsamp(~1 ~1, umf, starts=c(log(lam), log(sig)), rel.tol=1e-3)
    simout1[i,] <- exp(coef(m))
    }
hist(simout1[,1]); abline(v=lam, lwd=2, col=3)
hist(simout1[,2]); abline(v=sig, lwd=2, col=3)




