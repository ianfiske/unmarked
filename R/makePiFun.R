
# Factory functions to generate piFuns for use with 'unmarked' multinomial models
# ===============================================================================
# piFuns can only take 1 argument, so the only way to change the times for 'remPiFun' is
#  to bake it into the function environment. For the others, the factory approach is not
#  essential, but it does avoid re-evaluating all the preliminary stuff which depends only
#  on the number of occasions every time that the piFun is called - which is many times
#  for each fit.

# Mike Meredith, 2019-12-30

# Utility function to get capture histories
# nOcc : the number of survey occasions
getCombs <- function(nOcc) {
  combs <- matrix(NA, 2^nOcc, nOcc)
  for(i in 1:nOcc)
    combs[, i] <- rep(c(0, 1), each=2^(nOcc - i))
  combs <- combs[-1, ]  # Remove the all-zero capture history
  return(combs)
}

# Removal model with varying time intervals
# =========================================
# times : a vector of times for each interval, length(times)
#   is the number of survey occasions; can be all 1's if times the same.

# The function produced expects a sites x occasions matrix with probability of
#   detection *per unit time*.

makeRemPiFun <- function(times) {
  force(times)

  instRemPiFun <- function(p){
    nSites <- nrow(p)
    nOccs <- ncol(p)
    if(ncol(p) != length(times))
      stop("You have ", ncol(p), " occasions, but piFun expects ", length(times), ".", call.=FALSE)

    for(i in 1:nOccs)
      p[,i] <- 1 - (1 - p[,i])^times[i]
    pi <- p
    for(i in 2:nOccs)
      pi[,i] <- pi[, i - 1]/p[, i - 1] * (1 - p[, i - 1]) * p[, i]
    return(pi)
  }
  return(instRemPiFun)
}

# Basic capture-recapture model: M0, Mt, Mx
# =========================================
# nOcc : the number of survey occasions

# The function produced expects a sites x occasions matrix with probabilities of detection.

makeCrPiFun <- function(nOcc) {
  force(nOcc)
  combs <- getCombs(nOcc)
  histories <- apply(combs, 1, paste0, collapse="")

  crPiFun <- function(p) {
    if(ncol(p) != nOcc)
      stop("You have ", ncol(p), " occasions, but piFun expects ", nOcc, ".", call.=FALSE)
    out <- matrix(NA,  nrow(p), nrow(combs))
    colnames(out) <- histories
    q <- t(p)  # change to occasions x sites as multiplication below is done column-wise.
    for(i in 1:ncol(out)) {
      temp <- (1 - q) * (1 - combs[i, ]) + q * combs[i, ]
      out[, i] <- apply(temp, 2, prod)
    }
    return(out)
  }
  return(crPiFun)
}

# Models with a behavioural response: Mb
# ======================================
# Since p's do not change across occasions, we just need to sum up the captures for
#  naive and wise animals.
# nOcc : the number of survey occasions

# The function produced expects a sites x occasions matrix with probability of detection for naive
#   animals in column #1 and for wise animals in column #2 (!!); other entries are ignored.

makeCrPiFunMb <- function(nOcc) {
  force(nOcc)
  combs <- getCombs(nOcc)
  histories <- apply(combs, 1, paste0, collapse="")

  firstCap <- apply(combs, 1, function(x) which(x > 0)[1])
  naiveNon <- firstCap - 1                  # number of non-captures for naive animal
  wiseCap <- rowSums(combs) - 1             # number of captures for wise animal
  wiseNon <- nOcc - naiveNon - wiseCap - 1  # number of non-captures for wise animal
    #... and of course naiveCap == 1 everywhere

  crPiFunMb <- function(p) {
    if(ncol(p) != nOcc)
      stop("You have ", ncol(p), " occasions, but piFun expects ", nOcc, ".", call.=FALSE)
    pNaive <- p[,1]
    pWise <- p[,2]   # AHM1 p.355 uses column #3

    out <- matrix(NA,  nrow(p), nrow(combs))
    colnames(out) <- histories
    for(i in 1:ncol(out))
      out[, i] <- (1 - pNaive)^naiveNon[i] *
                   pNaive *
                   (1 - pWise)^wiseNon[i] *
                   pWise^wiseCap[i]
    return(out)
  }
  return(crPiFunMb)
}

# Models with individual heterogeneity: Mh
# ========================================
# Since p does not change across occasions, the order of captures/noncaptures does not
#  matter, only the number of each. So 001, 010, and 100 have 1 capture and will have
#  same value.
# Note that plogis(-x) == 1 - plogis(x) but with better precision when plogis(x) is close to 1.
# nOcc : the number of survey occasions

# The function produced expects a sites x occasions matrix with the mean probability of detection
#   in column #1; p[1, 2] is a value in [0, 1] which controls the scale parameter for the normal
#   distribution; other entries are ignored.

makeCrPiFunMh <- function(nOcc) {
  force(nOcc)
  combs <- getCombs(nOcc)
  histories <- apply(combs, 1, paste0, collapse="")

  nCaps <- rowSums(combs)  # number of captures

  crPiFunMh <- function(p) {
    if(ncol(p) != nOcc)
      stop("You have ", ncol(p), " occasions, but piFun expects ", nOcc, ".", call.=FALSE)
    mu <- stats::qlogis(p[,1])       # logit(p)
    sig <- exp(stats::qlogis(p[1,2]))
    nOcc <- ncol(p)
    nSite <- nrow(p)
    out <- matrix(NA,  nrow(p), nrow(combs))
    colnames(out) <- histories
    for(i in 1:nSite)
      for(j in 1:nOcc)
        out[i, nCaps==j] <- integrate( function(x) {
            plogis(mu[i]+x)^j * plogis(-(mu[i]+x))^(nOcc-j) * dnorm(x,0,sig)
            }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
    return(out)
  }
  return(crPiFunMh)
}
