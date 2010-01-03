\name{piFuns}
\alias{removalPiFun}
\alias{doublePiFun}
\alias{piFuns}
\title{Compute multinomial cell probabilities}
\usage{
removalPiFun(p)
doublePiFun(p)
}

\description{Compute the cell probabilities used in the multinomial-Poisson 
	model \link{multinomPois}.}

\arguments{\item{p}{
matrix of detection probabilities at each site for each observation}
}

\details{These two functions are provided as examples of 
possible functions to calculate multinomial cell probabilities. Users may 
write their own functions for specific sampling designs and use them in 
\link{multinomPois}. See the example for a user-defined function.
}

\value{
For removalPiFun, a matrix of cell probabilties for each site and sampling period.

For doublePiFun, a matrix of cell probabilities for each site and observer
combination. Column one is probability observer 1 but not observer 2 
detects the object, column two is probability that observer 2 but not 
observer 1 detects the object, and column 3 is probability of both 
detecting.
}
\examples{

(pRem <- matrix(0.5, nrow=3, ncol=3))	# Capture probabilities
removalPiFun(pRem)			# Cell probs

(pDouble <- matrix(0.5, 3, 2))		# Observer detection probs
doublePiFun(pDouble)			# Cell probs

# A user-defined piFun calculating removal probs when time intervals differ.
# Here 10-minute counts were divided into 2, 3, and 5 minute intervals.
# This function could be supplied to unmarkedFrameMPois along with the obsToY
# argument shown below.

instRemPiFun <- function(p) {
	M <- nrow(p)
	J <- ncol(p)
	pi <- matrix(NA, M, J)
	p[,1] <- pi[,1] <- 1 - (1 - p[,1])^2
	p[,2] <- 1 - (1 - p[,2])^3
	p[,3] <- 1 - (1 - p[,3])^5
	for(i in 2:J) {
		pi[,i] <- pi[, i - 1]/p[, i - 1] * (1 - p[, i - 1]) * p[, i]
		}
	return(pi)
	} 
	
instRemPiFun(pRem)

# Associated obsToY matrix required by unmarkedFrameMPois
o2y <- diag(3)
o2y[upper.tri(o2y)] <- 1
o2y

}		