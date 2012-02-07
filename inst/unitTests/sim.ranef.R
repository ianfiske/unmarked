

# Evaluate bias of posterior mode

# ----------------------------- pcount ----------------------------------


library(unmarked)

sim.nmix <- function(R=100, J=5, lambda=5, p=0.7) {
    N <- rpois(R, lambda)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, N, p)
    return(list(y=y, N=N))
}


nsim <- 200
out.nmix <- matrix(NA, nsim, 2)
set.seed(8345)
for(i in 1:nsim) {
    lambda <- 10
    sim.i <- sim.nmix(J=10, lambda=lambda, p=0.5)
    umf <- unmarkedFramePCount(y=sim.i$y)
    K <- 40
    fm <- pcount(~1 ~1, umf, K=K, se=FALSE)
    lam.hat <- exp(coef(fm, type="state"))
    re <- ranef(fm)
    N <- sum(sim.i$N)
    N.hat <- sum(modes <- postMode(re))
    bias <- N.hat - N
    ci <- colSums(confint(re))
    cover <- (N >= ci[1] & N <= ci[2])
    out.nmix[i,] <- c(mean(modes), lam.hat)
    cat("sim", i, "\n")
    cat("  bias =", mean(modes)-lambda, "\n")
}

hist(out.nmix[,1], breaks=20)

colMeans(out.nmix)

plot(re, layout=c(5,5))






# ------------------------------- occu ----------------------------------

library(unmarked)




sim.occu <- function(R=100, J=5, psi=0.5, p=0.4) {
    z <- rbinom(R, 1, psi)
    y <- matrix(NA, R, J)
    y[] <- rbinom(R*J, 1, z*p)
    return(list(y=y, z=z))
}


nsim <- 300
out.occu <- matrix(NA, nsim, 2)
set.seed(38845)
for(i in 1:nsim) {
    cat("sim", i, "\n")
    sim.i <- sim.occu(psi=0.8, p=0.4)
    umf <- unmarkedFrameOccu(y=sim.i$y)
    K <- 30
    fm <- occu(~1 ~1, umf, se=FALSE)
    re <- ranef(fm)
    pao <- sum(sim.i$z)
    pao.hat <- sum(modes <- postMode(re))
    bias <- pao.hat - pao
    ci <- colSums(confint(re))
    cover <- (pao >= ci[1] & pao <= ci[2])
    out.occu[i,] <- c(bias, mean(modes))
}

hist(out.occu[,1])

colMeans(out.occu)

















# ------------------------------ distsamp -------------------------------




lambda <- 10
sigma <- 30
npts <- 100
radius <- 50
breaks <- seq(0, 50, by=10)
A <- (2*radius)^2 / 10000 # Area (ha) of square containing circle
y <- matrix(0, npts, length(breaks)-1)
N <- integer(npts)
for(i in 1:npts) {
    M <- rpois(1, lambda * A) # Individuals within the square
    xy <- cbind(x=runif(M, -radius, radius), y=runif(M, -radius, radius))
    d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
    d <- d[d <= radius]
    N[i] <- length(d)
    if(length(d)) {
        p <- exp(-d^2 / (2 * sigma^2)) # half-normal
        d <- d[rbinom(length(d), 1, p) == 1]
        y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
    }
}

max(N)
mean(N)

set.seed(3)
umf1 <- unmarkedFrameDS(y = y, survey="point",
    dist.breaks=breaks, unitsIn="m")
(m1 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20))))
(m2 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)), output="abund"))

backTransform(m1, type="state")
backTransform(m1, type="det")


re1 <- ranef(m1, K=20)
plot(re1)

re2 <- ranef(m2, K=20)
plot(re2)

all(postMode(re1) == postMode(re2))
all(confint(re1) == confint(re2))


reM <- postMode(re1)
reCI <- confint(re1)

sum(N)
sum(reM)
colSums(reCI)











# ------------------------------ multinomPois ----------------------------








# Simulate independent double observer data
nSites <- 50
lambda <- 10
p1 <- 0.5
p2 <- 0.3
cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
set.seed(9023)
N <- rpois(nSites, lambda)
y <- matrix(NA, nSites, 3)
for(i in 1:nSites) {
  y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
}

# Fit model
observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
    type="double")
fm <- multinomPois(~observer-1 ~1, umf)

# Estimates of fixed effects
e <- postMode(fm)
exp(e[1])
plogis(e[2:3])

# Estimates of random effects
re <- ranef(fm, K=20)
ltheme <- canonical.theme(color = FALSE)
lattice.options(default.theme = ltheme)
plot(re, layout=c(10,5))


sum(postMode(re))
colSums(confint(re))
sum(N)






# ---------------------------- gmultmix ---------------------------------


library(unmarked)
n <- 100  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals available for detection

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)    # Observe some
    Nleft1 <- N[i,] - y[i,,1]         # Remove them
    y[i,,2] <- rbinom(T, Nleft1, p)   # ...
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }

y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")

(m1 <- gmultmix(~1, ~1, ~1, data=umf1, K=30))

re <- ranef(m1)
plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)







# ------------------------------ gdistsamp ------------------------------

set.seed(36837)
R <- 50 # number of transects
T <- 5  # number of replicates
strip.width <- 50
transect.length <- 60 # so that abund != density
breaks <- seq(0, 50, by=10)

lambda <- 10 # Abundance
phi <- 0.6   # Availability
sigma <- 30  # Half-normal shape parameter

J <- length(breaks)-1
y <- array(0, c(R, J, T))
for(i in 1:R) {
    M <- rpois(1, lambda) # Individuals within the 1-ha strip
    for(t in 1:T) {
        # Distances from point
        d <- runif(M, 0, strip.width)
        # Detection process
        if(length(d)) {
            cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
            d <- d[rbinom(length(d), 1, cp) == 1]
            y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
y <- matrix(y, nrow=R) # convert array to matrix

# Organize data
umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
    dist.breaks=breaks, tlength=rep(transect.length, R), numPrimary=T)
summary(umf)

# Fit the model
m1 <- gdistsamp(~1, ~1, ~1, umf, output="abund", K=50)
summary(m1)
m2 <- gdistsamp(~1, ~1, ~1, umf, output="density", K=50)
summary(m2)


re1 <- ranef(m1)
plot(re1, xlim=c(-1, 30))
re2 <- ranef(m2)
plot(re2, xlim=c(-1, 30))

all(postMode(re1) == postMode(re2))
all(confint(re1) == confint(re2))

cbind(postMode(re1), postMode(re2))


# ----------------------------- colext ----------------------------------











# ----------------------------- pcountOpen -------------------------------




library(unmarked)
set.seed(7)
M <- 100
J <- 3
T <- 10
lambda <- 5
gamma <- 0.4
omega <- 0.9
p <- 0.5
N <- matrix(NA, M, T)
y <- array(NA, c(M, J, T))
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
y[,,1] <- rbinom(M*J, N[,1], p)
for(t in 1:(T-1)) {
    S[,t] <- rbinom(M, N[,t], omega)
    G[,t] <- rpois(M, gamma)
    N[,t+1] <- S[,t] + G[,t]
    y[,,t+1] <- rbinom(M*J, N[,t+1], p)
}


colSums(N)
colMeans(N)


# Prepare data
umf <- unmarkedFramePCO(y = matrix(y, M), numPrimary=T)
summary(umf)


# Fit model and backtransform
(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20))

e <- postMode(m1)
(lam <- exp(e[1]))
(gam <- exp(e[2]))
(om <- plogis(e[3]))
(p <- plogis(e[4]))

re <- ranef(m1)

plot(re, layout=c(5,5), subset = site %in% 1:25 & year %in% 1,
     xlim=c(-1,20))

postMode(re)
confint(re)

N.hat <- colSums(postMode(re))
CI <- apply(confint(re), c(2,3), sum)
rbind(N=colSums(N), N.hat=N.hat)

plot(1:T, N.hat, ylim=c(0, 1000), cex=1.5)
points(1:T, colSums(N), pch=16, col="blue")
segments(1:T, CI[1,], 1:T, CI[2,])





(m1nt <- update(m1, dynamics="notrend"))

re <- ranef(m1nt)

plot(re, layout=c(5,5), subset = site %in% sites, xlim=c(-1,10))


