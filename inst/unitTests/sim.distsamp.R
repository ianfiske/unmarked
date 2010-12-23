



# Simulate with half-normal detection function

simPt <- function(lambda=5, sigma=20, npts=100, radius=50,
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

colSums(simPt())

set.seed(3)
umf1 <- unmarkedFrameDS(y = simPt(), survey="point", 
    dist.breaks=seq(0, 50, by=10), unitsIn="m")
(m1 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20))))
(m2 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)), output="abund"))


checkEqualsNumeric(coef(m1), c(1.813819, 2.893771), tol=1e-5)
checkEquals(exp(coef(m1, type="state")), 
    exp(coef(m2, type="state")) / (pi * 50^2 / 10000), tol=0.01)


set.seed(11)
nsims <- 50
simout1 <- matrix(NA, nsims, 2)
lam <- 20
sig <- 30
for(i in 1:nsims) {
    cat("sim", i, "\n"); flush.console()
    umf <- unmarkedFrameDS(y = simPt(lambda=lam, sigma=sig), survey="point", 
        dist.breaks=seq(0, 50, by=10), unitsIn="m")
    m <- distsamp(~1 ~1, umf, starts=c(log(lam), log(sig)), rel.tol=1e-3, 
        output="abund")
    simout1[i,] <- exp(coef(m))
    }
hist(simout1[,1]); abline(v=lam*pi*50^2/10000, lwd=2, col=3)
hist(simout1[,2]); abline(v=sig, lwd=2, col=3)






integrate(unmarked:::grhn, 0, 10, sigma=1000)$value * 2 * pi
integrate(unmarked:::grhn, 10, 20, sigma=1000)$value * 2 * pi












simLine <- function(lambda=5, sigma=20, npts=100,
    breaks=seq(0, 50, by=10))
{      
    W <- max(breaks)
    A <- 2*W*100/10000 # Area (ha) of rectangle containing 100m line
    y <- matrix(0, npts, length(breaks)-1)
    for(i in 1:npts) {
        N <- rpois(1, lambda * A) # Individuals within the square
        # distance from the line
        d <- runif(N, 0, W)

        # Detection process
        if(length(d) > 0) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    return(y)
}

simLine()




set.seed(7)
nsims <- 100
simout2 <- matrix(NA, nsims, 2)
lam <- 20
sig <- 30
for(i in 1:nsims) {
    cat("sim", i, "\n"); flush.console()
    y.sim <- simLine(lambda=lam, sigma=sig)
    umf <- unmarkedFrameDS(y = y.sim, survey="line", 
        dist.breaks=seq(0, 50, by=10), unitsIn="m", 
        tlength=rep(100, nrow(y.sim)))
    m <- distsamp(~1 ~1, umf, starts=c(log(lam), log(sig)), rel.tol=1e-3)
    simout2[i,] <- exp(coef(m))
    }
hist(simout2[,1]); abline(v=lam, lwd=2, col=3)
hist(simout2[,2]); abline(v=sig, lwd=2, col=3)

