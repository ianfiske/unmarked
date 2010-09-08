



lik0 <- function(pars, data, K=20) 
{


ldm <- function(x, sx, probs, kvec)
    lgamma(kvec+1) - lgamma(kvec - sx+1) + sum(x*log(probs[1:3])) + 
        (kvec-sx)*log(probs[4])


lam <- exp(pars[1])
phi <- plogis(pars[2])
p1 <- plogis(pars[3])

nsite <- nrow(data)

ll <- rep(NA, nsite)
for(i in 1:nsite) 
{

f <- dpois(0:(K-1), lam)
f <- f/sum(f)

# Primary period counts
sx <- c(sum(data[i, 1:3]), sum(data[i, 4:6]), sum(data[i, 7:9]))
k <- max(sx, na.rm=T):(K-1)
F <- f[(max(sx,na.rm=T)+1):K]

p2min <- 1-(1-p1)^2
p3min <- 1-(1-p1)^3
p5min <- 1-(1-p1)^5

# Cell probabilities
cp <- c(phi*p5min, phi*(1-p5min)*p3min, phi*(1-p5min)*(1-p3min)*p2min)
cp <- c(cp, 1-sum(cp))

a1 <- ldm(data[i, 1:3], sx[1], cp, k)
a2 <- ldm(data[i, 4:6], sx[2], cp, k)
a3 <- ldm(data[i, 7:9], sx[3], cp, k)
A <- cbind(a1,a2,a3)

ll1 <- exp(rowSums(A))        
ll[i] <- sum(ll1 * F)

}
-sum(log(ll))  
}

