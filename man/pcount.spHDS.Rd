\name{pcount.spHDS}
\alias{pcount.spHDS}
\encoding{UTF-8}

\title{
Fit spatial hierarchical distance sampling model
}
\description{
Function fits spatial hierarchical distance sampling model; \pkg{unmarked} model fitting function introduced in Section 9.8.4. See \code{\link[unmarked]{pcount}} for details.
}
\usage{
pcount.spHDS(formula, data, K, mixture = c("P", "NB", "ZIP"), starts,
  method = "BFGS", se = TRUE, engine = c("C", "R"), ...)
}
\arguments{
  \item{formula}{
Double right-hand side formula describing covariates of detection and abundance, in that order
}
  \item{data}{
an \code{unmarkedFramePCount} object supplying data to the model.
}
  \item{K}{
Integer upper index of integration for N-mixture. This should be set high enough so that it does not affect the parameter estimates. Note that computation time will increase with K.
}
  \item{mixture}{
character specifying mixture: Poisson (P), Negative-Binomial (NB), or Zero Inflated Poisson (ZIP).
}
  \item{starts}{
vector of starting values
}
  \item{method}{
Optimization method used by \code{\link{optim}}.
}
  \item{se}{
logical specifying whether or not to compute standard errors.
}
  \item{engine}{
Either "C" or "R" to use fast C++ code or native R code during the optimization.
}
  \item{\dots}{
Additional arguments to optim, such as lower and upper bounds
}
}

\value{
unmarkedFit object describing the model fit.
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.8.4
}
\author{
Kéry & Royle
}

\examples{

# From Kery and Royle (2016) AHM volume 1


# Simulator function for spatial distance sampling data
sim.spatialDS <- 
function(N=1000, beta = 1, sigma=1, keep.all=FALSE, B=B, model="halfnorm"){
# Function simulates coordinates of individuals on a square
# Square is [0,6] x [0,6], with a count location on the point (3,3)
#   N: total population size in the square
#   beta: coefficient of SOEMTHING on spatial covariate x
#   sigma: scale of half-normal detection function
#   B: circle radius
#   keep.all: return the data for y=0 individuals or not
library(raster)      # Load required packages
library(plotrix)

# Create coordinates for 30 x 30 grid
delta <- (2*B-0)/30                # '2D bin width'
grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
gr <- expand.grid(grx,grx)         # Create grid coordinates

# Create spatially correlated covariate x and plot it
V <- exp(-e2dist(gr,gr)/1)
x <- crossprod(t(chol(V)), rnorm(900) )
par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(gr,x)), col=topo.colors(10))
draw.circle(3, 3, B)
points(3, 3, pch="+", cex=3)
image.scale(x, col=topo.colors(10))

# Simulate point locations as function of habitat covariate x
probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
pixel.id <- sample(1:900, N, replace=TRUE, prob=probs)
# could simulate randomly within the pixel but it won't matter so place centrally
u1 <- gr[pixel.id,1]
u2 <- gr[pixel.id,2]
points(u1, u2, pch=20, col='black', cex = 1)  # plot points
title("This is so cool !")         # express your appreciation of all this

d <- sqrt((u1 - B)^2 + (u2-B)^2)   # distance to center point of square
#plot(u1, u2, pch = 1, main = "Point transect")
N.real <- sum(d<= B)               # Population size inside of count circle

# Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
# p <- ifelse(d< B, 1, 0) * exp(-d*d/(2*(sigma^2)))
# We do away with the circle constraint here.   
if(model=="hazard")
   p <- 1-exp(-exp(-d*d/(2*sigma*sigma)))
if(model=="halfnorm")
   p <- exp(-d*d/(2*sigma*sigma))
# Now we decide whether each individual is detected or not
y <- rbinom(N, 1, p)                                           # detected or not
#points(u1[d<= B], u2[d<= B], pch = 16, col = "black", cex = 1) # not detected
points(u1[y==1], u2[y==1], pch = 16, col = "red", cex = 1)     # detected

# Put all of the data in a matrix
if(!keep.all){
   u1 <- u1[y==1]
   u2 <- u2[y==1]
   d <- d[y==1]   }
# Output
return(list(model=model, N=N, beta=beta, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real, Habitat=x, grid=gr))
}

# Generate one data set
library(scrbook) # for e2dist   # Put this in Marc’s file!
# Simulate a data set and harvest the output
set.seed(1234)
str(tmp <- sim.spatialDS(N=200, beta=1, sigma=1.5, keep.all=FALSE, B=3))

# Harvest data
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Do data augmentation, including for pixel ID
M <- 400
nz <- M-nind
# We use discrete "pixel ID" here instead of "s"
pixel <- rep(NA, M)
y <- c(rep(1,nind), rep(0,nz))

# Pick some starting values and figure out the pixel of each observation
s <-cbind(u1,u2)
s <- rbind(s, matrix(NA,nrow=nz,ncol=2))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
  pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Bundle and summarize the data for BUGS
str(data <- list (B=B, nind=nind, y=y, nz=nz, Habitat=Habitat, Habgrid=Habgrid, G=G, pixel=pixel))



# Simulate a data set, N = 600 for the population size
set.seed(1234)
tmp <-sim.spatialDS(N=600, sigma=1.5, keep.all=FALSE, B=3, model= "hazard")

# Harvest stuff
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Find which pixel each observation belongs to
s <- cbind(u1,u2)
D <- e2dist(s[1:nind,], Habgrid)
pixel <-rep(NA,nind)
for(i in 1:nind){
   pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Create a vector of counts in each pixel and pad it with zeros
pixel.count <- rep(0, G)
names(pixel.count) <- 1:G
pixel.count[names(table(pixel))] <- table(pixel)
# Create a covariate: distance between observer and pixel center
dist <- sqrt( (Habgrid[,1]-3)^2 + (Habgrid[,2]-3)^2  )
# Construct an unmarkedFrame
umf <- unmarkedFramePCount(y=matrix(pixel.count,ncol=1), 
   siteCovs=data.frame(dist=dist,Habitat=Habitat))
summary(umf)
 
(fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K=20))


}
