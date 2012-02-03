





library(unmarked)

set.seed(4564)
R <- 20
J <- 5
N <- rpois(R, 3)
y <- matrix(NA, R, J)
y[] <- rbinom(R*J, N, 0.5)
y[1,] <- NA

K <- 15

umf <- unmarkedFramePCount(y=y)
fm <- pcount(~1 ~1, umf, K=K)


(re <- ranef(fm))
coef(re)
confint(re, level=0.9)
plot(re, xlim=c(-1,10))

sum(N)
sum(coef(re))
colSums(confint(re))




library(unmarked)

set.seed(320)
R <- 50
J <- 3
z <- rbinom(R, 1, 0.5)
y <- matrix(NA, R, J)
y[] <- rbinom(R*J, z, 0.3)

visit <- matrix(as.character(1:J), R, J, byrow=TRUE)
umf <- unmarkedFrameOccu(y=y, obsCovs=list(visit=visit))
(fm <- occu(~visit ~1, umf))



(re <- ranef(fm))
coef(re)
confint(re, level=0.9)
plot(re)

sum(z)
sum(coef(re))
colSums(confint(re))







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

e <- coef(m1)
(lam <- exp(e[1]))
(gam <- exp(e[2]))
(om <- plogis(e[3]))
(p <- plogis(e[4]))

re <- ranef(m1)

sites <- paste("site", 1:25, sep="")
years <- paste("year", 1:1, sep="")

plot(re, layout=c(5,5), subset = site %in% sites & year %in% years,
     xlim=c(-1,20))

coef(re)
confint(re)

N.hat <- colSums(coef(re))
CI <- apply(confint(re), c(2,3), sum)
rbind(N=colSums(N), N.hat=N.hat)

plot(1:T, N.hat, ylim=c(0, 1000), cex=1.5)
points(1:T, colSums(N), pch=16, col="blue")
segments(1:T, CI[1,], 1:T, CI[2,])





(m1nt <- update(m1, dynamics="notrend"))

re <- ranef(m1nt)

plot(re, layout=c(5,5), subset = site %in% sites, xlim=c(-1,10))


