
# data will need to be an unmarkedMultFrame
gmn <- function(lambdaformula, phiformula, pformula, data, mixture=c('P', 'NB'),
    K, starts, method = "BFGS", control = list(), se = TRUE)
{
if(!is(data, "unmarkedMultFrame"))
    stop("Data is not an unmarkedMultFrame.")

mixture <- match.arg(mixture)
k <- 0:K

D <- getDesign(data, formula)
Xlam <- D$Xlam
Xphi <- D$Xphi 
Xdet <- D$Xdet
y <- D$y

Xlam.offset <- D$X.offset
Xphi.offset <- D$Xphi.offset
Xdet.offset <- D$Xdet.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))
  
J <- ncol(y)
R <- obsNum(data)
M <- nrow(y)
piFun <- data@piFun

lamParms <- colnames(Xlam)
phiParms <- colnames(Xphi)
detParms <- colnames(Xdet)
nLP <- ncol(Xlam)
nPP <- ncol(Xphi)
nDP <- ncol(Xdet)
nP <- nLP + nPP + nDP + ifelse(mixture=='NB', 1, 0)

# y should probably be an MxRxJ array 
yvec <- as.numeric(y)
navec <- is.na(yvec)

# n.mr should be an MxR matrix of counts
n.mr <- apply(y.mrj, 1:2, sum)

ldm <- function(x, sx, pr, kvec) {
    
    lgamma(kvec+1) - lgamma(kvec-sx+1) + sum(x*log(pr[1:3])) + (kvec-sx)*log(pr[4])
    
    }


nll <- function(pars) 
{
		lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset) 
    phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
    p <- plogis(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)
    p.array <- array(p, c(M, R, J))
    
    f <- sapply(lambda, function(x) dpois(k, lambda))

    for(r in 1:R) {

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
    
    ll <- rowSums(ll1 * f)
    -sum(log(ll))  
}

