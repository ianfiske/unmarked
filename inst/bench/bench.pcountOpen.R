
library(unmarked)
library(rbenchmark)


set.seed(34537)
nSites <- 100
T <- 20
covariates <- data.frame(veght=rnorm(nSites),
    habitat=factor(c(rep('A', 50), rep('B', 50))))
lampars <- c(-1, 1, -1)
gampars <- c(-1,0, 0)
ompars <- c(2, -1, 0)
ppars <- c(1, -1, 0)
X <- model.matrix(~veght+habitat, covariates) # design matrix
lam <- exp(X %*% lampars)
gam <- exp(X %*% gampars)
om <- plogis(X %*% ompars)
p <- plogis(X %*% ppars)
y <- N <- matrix(NA, nSites, T)
N[,1] <- rnbinom(nSites, size=1, mu=lam)       # true abund
for(t in 2:T) {
    S <- rbinom(nSites, N[,t-1], om)
    G <- rpois(nSites, gam)
    N[,t] <- S+G
    }
y[] <- rbinom(nSites*T, N, p)
#y[1,5] <- NA
#covariates[2,] <- NA
umf <- unmarkedFramePCO(y = y, siteCovs = covariates, numPrimary=T)
head(umf)
summary(umf)

{
st1 <- system.time(fmR <- pcountOpen(~veght+habitat, ~1, ~veght, ~veght,
                                     umf,
                              K=20, control=list(trace=TRUE, REPORT=1),
                              se=FALSE, engine="R")) # 886
st2 <- system.time(fmC <- pcountOpen(~veght+habitat, ~1, ~veght, ~veght,
                                     umf,
                              K=20, control=list(trace=TRUE, REPORT=1),
                              se=FALSE, engine="C")) # 722
}

benchmark(pcountOpen(~veght+habitat, ~1, ~veght, ~veght, umf, K=20,
                     control=list(trace=TRUE, REPORT=1)),
          columns=c("test", "elapsed", "relative"),
          replications=50)



st1 <- system.time(fm2R <- pcountOpen(~1, ~1, ~veght, ~1, umf, K=20,
                               control=list(trace=TRUE, REPORT=1, maxit=2),
                               se=FALSE, engine="R"))

st2 <- system.time(fm2C <- pcountOpen(~1, ~1, ~veght, ~1, umf, K=20,
                               control=list(trace=TRUE, REPORT=1, maxit=2),
                               se=FALSE, engine="C"))


