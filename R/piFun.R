## Take an M x J matrix of detection probabilities and return a matrix
## of M x J observation probs
# Compute the cell probabilities for the observation classes
# in removal sampling.
#
# Both p and the returned matrix are M x J for M sites and J sampling occasions.

removalPiFun <- function(p){
  M <- nrow(p)
  J <- ncol(p)
  pi <- matrix(NA, M, J)
  pi[,1] <- p[,1]
  for(i in seq(from = 2, length = J - 1)) {
	  pi[, i] <- pi[,i-1] / p[,i-1] * (1-p[,i-1]) * p[,i]
  }
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# returns an M x 3 matrix of row=(1 not 2, 2 not 1, 1 and 2).
# Compute the cell probabilities for the observation classes
# in double observer sampling.

doublePiFun <- function(p){
  M <- nrow(p)
  pi <- matrix(NA, M, 3)
  pi[,1] <- p[,1] * (1 - p[,2])
  pi[,2] <- p[,2] * (1 - p[,1])
  pi[,3] <- p[,1] * p[,2]
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# returns an M x 2 matrix of row=(1, 2 not 1).
# Compute the cell probabilities for the observation classes
# in double observer sampling.

depDoublePiFun <- function(p){
  M <- nrow(p) 
  pi <- matrix(NA, M, 2) 
  pi[,1] <- p[,1] 
  pi[,2] <- p[,2]*(1-p[,1]) 
  return(pi) 
}
