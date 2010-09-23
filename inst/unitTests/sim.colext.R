


## Multiseason occupancy

	
n <- 100
J <- 5
T <- 10
dbhMat <- matrix(rnorm(n*T), n, T, byrow=T)
yearMat <- matrix(1:T, n, T, byrow=T)
yearlyCovs <- data.frame(dbh=matrix(t(dbhMat), n*T, 1), 
	year=matrix(t(yearMat), ncol=1))
detCovs <- data.frame(wind=rnorm(n * J * T))
treatment <- gl(2, 50, labels = c("Burned", "Mechanical"))
year <- matrix(1:T, n, T, byrow=T)
wind <- array(rnorm(n*J*T), c(n, J, T))
y <- p <- array(NA, c(n, J, T))
p[] <- 0.5


Xpsi <- model.matrix(~treatment + yearlyCovs[1:n, "dbh"])
Xyearly <- model.matrix(~year + dbh, yearlyCovs)

psi <- 0.6 #plogis(Xpsi %*% c(-1, 1, 0.5))
gamma <- matrix(plogis(Xyearly %*% c(-1, -0.1, 0)), n, T, byrow=T)
phi <- matrix(plogis(Xyearly %*% c(3, -0.5, 0)), n, T, byrow=T)


colMeans(gamma)
colMeans(phi)

#gamma <- runif(T-1)
#phi <- runif(T-1)

Z <- matrix(NA, n, T)
Z[,1] <- rbinom(n, 1, psi)
for(t in 2:T) {
	muZ <- Z[,t-1] * phi[,t] + (1 - Z[,t-1]) * gamma[,t]
	Z[,t] <- rbinom(n, 1, muZ)
	}
str(Z)
colSums(Z)

for(j in 1:J)
for(t in 1:T)
	y[,j,t] <- rbinom(n, 1, Z[,t]*p[,j,t])


library(unmarked)

mowaUMF <- unmarkedMultFrame(y = matrix(y, n, J*T), 
	siteCovs = data.frame(treatment),
	obsCovs = detCovs,
	yearlySiteCovs = yearlyCovs, 
	numPrimary = T)


#null <- colext(~1~1, mowaUMF, control=list(trace=1))	
(m1 <- colext(~1 ~ year + dbh, mowaUMF, control=list(trace=T)))
