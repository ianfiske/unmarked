





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







set.seed(3)
M <- 50
T <- 5
lambda <- 1
gamma <- 0.5
omega <- 0.8
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
    S[,t] <- rbinom(M, N[,t], omega)
    G[,t] <- rpois(M, gamma)
    N[,t+1] <- S[,t] + G[,t]
}
y[] <- rbinom(M*T, N, p)


# Prepare data
umf <- unmarkedFramePCO(y = y, numPrimary=T)
summary(umf)


# Fit model and backtransform
(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20))

re <- ranef(m1)
sites <- paste("site", 1:25, sep="")
plot(re, layout=c(5,5), subset = site %in% sites, xlim=c(-1,10))

coef(re)
confint(re)



(m1nt <- update(m1, dynamics="notrend"))

re <- ranef(m1nt)

plot(re, layout=c(5,5), subset = site %in% sites, xlim=c(-1,10))


