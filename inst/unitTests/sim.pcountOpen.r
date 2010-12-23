## Simulate no covariates, constant sampling period intervals	
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
# y[M-5, 1:4] <- y[M-4, 2:5] <- y[M-3, 2:4] <- y[M-2, c(2,4)] <- y[M-1, 1] <- y[M, T] <- NA

                            
# Prepare data                               
umf <- unmarkedFramePCO(y = y)

summary(umf)

# Fit model and backtransform
system.time(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, 
    K=10))  # 14s on 64bit. K might be too small

backTransform(m1, "lambda") # 0.71 initial abundance
backTransform(m1, "gamma")  # 0.45 recruitment rate
backTransform(m1, "omega")  # 0.751 survival rate
backTransform(m1, "det")    # 0.718 detection probability


# Fix omega at 1 (no losses)
m1.fixo <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10, fix="omega")

## Simulate covariate model with constant intervals
set.seed(33)
M <- 50
T <- 5
veght <- rnorm(M)
isolation <- matrix(rnorm(M*T), M, T)
time <- matrix(rnorm(M*T, 1), M, T)
y <- p <- N <- gamma <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
lambda <- exp(-1)   # + 0.5*veght)
gamma[] <- exp(-1 + -1*isolation)
omega <- 0.8
p[] <- plogis(-1 + 1*time)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega)
	G[,t] <- rpois(M, gamma[,t])
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)


# Prepare data
umfC <- unmarkedFramePCO(y = y, siteCovs = data.frame(veght), 
	obsCovs = list(isolation=isolation, time=time))

# Fit covariate model
(m2 <- pcountOpen(~1, ~isolation, ~1, ~time, umfC, K=15, se=FALSE,
	control=list(maxit=30, trace=TRUE, REPORT=1), starts=c(-1,-1,-1,1.5,-1,1)))


## Simulate uneven sampling period intervals
set.seed(333)
M <- 25
T <- 5
date <- matrix(NA, M, T)
date[,1] <- rpois(M, 2)
for(t in 2:T) {
    date[,t] <- date[,t-1] + rpois(M, 10)
    }
datediff <- t(apply(date, 1, diff))    
lambda <- 4
gamma <- 0.05   # Daily survival rate
omega <- 0.95   # Daily recruitment rate
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega^datediff[,t])   
	G[,t] <- rpois(M, gamma*datediff[,t])           
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)


# Prepare data
umfO <- unmarkedFramePCO(y = y, delta = date)
umfO

# Fit model
(m3 <- pcountOpen(~1, ~1, ~1, ~1, umfO, K=15, se=FALSE, starts=c(2,-3,3,0.85),
    control=list(maxit=15, trace=T, REPORT=1)))
backTransform(m3, "lambda")
backTransform(m3, "gamma")
backTransform(m3, "omega")
backTransform(m3, "det")





# Auto-regressive model


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
	G[,t] <- rpois(M, gamma * N[,t])
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)


# Prepare data
umf4 <- unmarkedFramePCO(y = y)
umf4

# Fit model
(m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=15, dynamics="autoreg", se=FALSE, 
    starts=c(2,-3,3,0.85), control=list(maxit=15, trace=T, REPORT=1)))


       
