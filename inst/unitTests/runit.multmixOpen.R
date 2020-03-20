simData <- function(lambda=1, gamma=0.5, omega=0.8, p=0.5, M=100, T=5,
                    p2=NULL, type="removal")
    {
    y <- array(NA, c(M, 3, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
  
  if(type=='removal'){
    for(i in 1:M) {
      N[i,1] <- rpois(1, lambda)

      y[i,1,1] <- rbinom(1, N[i,1], p)    # Observe some
      Nleft1 <- N[i,1] - y[i,1,1]         # Remove them
      y[i,2,1] <- rbinom(1, Nleft1, p)   # ...
      Nleft2 <- Nleft1 - y[i,2,1]
      y[i,3,1] <- rbinom(1, Nleft2, p)

    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,1,t+1] <- rbinom(1, N[i,t+1], p)    # Observe some
        Nleft1 <- N[i,t+1] - y[i,1,t+1]         # Remove them
        y[i,2,t+1] <- rbinom(1, Nleft1, p)   # ...
        Nleft2 <- Nleft1 - y[i,2,t+1]
        y[i,3,t+1] <- rbinom(1, Nleft2, p)
        }
    }
  } else if(type == "double"){  
    cp <- c(p*(1-p2), p2*(1-p), p*p2)
    for(i in 1:M) {
      N[i,1] <- rpois(1, lambda)
      y[i,,1] <- rmultinom(1, N[i,1], c(cp, 1-sum(cp)))[1:3]

    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,,t+1] <- rmultinom(1, N[i,t+1], c(cp, 1-sum(cp)))[1:3]
        }
    }

  }

    return(list(y=matrix(y, M),N=N))
}

test.multmixOpen.removal <- function(){

  set.seed(123)
  simy <- simData(lambda=4, gamma=0.5, omega=0.8, p=0.5, 
                M=100, T=5)

  sc <- data.frame(x1=rnorm(100))

  umf <- unmarkedFrameMMO(y=simy$y, numPrimary=5, siteCovs=sc,
                        type="removal")

  fit <- multmixOpen(~x1, ~1, ~1, ~x1, K=30, data=umf)

  checkEqualsNumeric(coef(fit), c(1.38860,0.043406,-0.68448,
                                  1.40703,0.03174,-0.00235), tol=1e-5)

  #Check predict
  pr <- predict(fit, type='lambda')
  checkEqualsNumeric(as.numeric(pr[1,]), 
                     c(3.79942,0.298279,3.25808,4.43193), tol=1e-4)

  #Check getP
  pv <- getP(fit)
  checkEqualsNumeric(dim(pv), dim(umf@y))
  checkEqualsNumeric(pv[1,1:3], pv[1,4:6])
  checkEqualsNumeric(pv[1,1:3], c(0.5086598,0.2499250,0.1227982), tol=1e-5)
  
  #Check residuals
  r <- residuals(fit)
  checkEqualsNumeric(r[1,1:3], c(0.067122,-0.9497006,0.533337), tol=1e-4)

  #Check simulate

  #Check ranef

  #Check parboot

  #Check NA handling

}


test.multmixOpen.double <- function(){

  set.seed(123)
  simy <- simData(lambda=4, gamma=0.5, omega=0.8, p=0.5, p2=0.5,
                M=300, T=5, type="double")

  umf <- unmarkedFrameMMO(y=simy$y, numPrimary=5, 
                          siteCovs=data.frame(x1=rnorm(300)),
                        type="double")

  fit <- multmixOpen(~x1, ~1, ~1, ~x1, K=30, data=umf)
  

  checkEqualsNumeric(coef(fit), c(1.405123,-0.037941,-0.52361,
                                  1.321799,0.070564,-0.0150329), tol=1e-4)

  pv <- getP(fit)
  checkEqualsNumeric(dim(pv), dim(umf@y))
  checkEqualsNumeric(pv[1,1:3], pv[1,4:6])
  
  checkEqualsNumeric(pv[1,1:3], c(0.2497033,0.2497033,0.26752), tol=1e-5)

}

test.multmixOpen.NB <- function(){

  set.seed(123)
  simy <- simData(lambda=4, gamma=0.5, omega=0.8, p=0.5, 
                M=100, T=5)

  sc <- data.frame(x1=rnorm(100))

  umf <- unmarkedFrameMMO(y=simy$y, numPrimary=5, siteCovs=sc,
                        type="removal")

  fit <- multmixOpen(~x1, ~1, ~1, ~x1, K=30, mixture="NB", data=umf)

  checkEqualsNumeric(coef(fit), c(1.38861,0.0433983,-0.68451,1.40705,
                                  0.031728,-0.002354,9.81414),tol=1e-5)

}

test.multmixOpen.dynamics(){
  #TODO
}
