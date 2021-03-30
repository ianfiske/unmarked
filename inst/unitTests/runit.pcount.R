test.pcount.offest <- function()
{

  y <- matrix(c(
      8,7,
      6,7,
      8,8,
      8,6,
      7,7), nrow=5, ncol=2, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10)
  umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- pcount(~ o1 ~ offset(x), data = umf, K=30)
  checkEqualsNumeric(coef(fm), structure(c(-0.78814924, 2.62569034, -0.02578801),
      .Names = c("lam(Int)", "p(Int)", "p(o1)")), tol = 1e-5)

}



test.pcount.covs <- function()
{
  y <- matrix(c(
      8,7,7,8,
      6,7,7,5,
      8,8,7,8,
      4,5,5,5,
      4,4,3,3), nrow=5, ncol=4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = seq(-1, 1, length=length(y)))
  umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- pcount(~ o1 ~ x, data = umf, K=30)
  checkEqualsNumeric(coef(fm),
                     c(1.91984184, -0.02987393,  2.49421875, -0.23350448),
                     tol = 1e-5)

}

test.pcount.randomeffects <- function()
{

set.seed(35)
nSites <- 300
nVisits <- 3
x <- rnorm(nSites)               # a covariate
beta0 <- 0
beta1 <- 0.4

ran <- rnorm(100, 0, 1)
group <- factor(as.character(rep(1:100, each=3)))
ran_ind <- as.numeric(group)

lambda <- exp(beta0 + beta1*x +
              ran[ran_ind])   # expected counts at each site
N <- rpois(nSites, lambda)       # latent abundance
y <- matrix(NA, nSites, nVisits)
p <- c(0.3, 0.6, 0.8)            # detection prob for each visit
for(j in 1:nVisits) {
  y[,j] <- rbinom(nSites, N, p[j])
}

# Organize data
visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)

umf <- unmarkedFramePCount(y=y, siteCovs=data.frame(x=x,group=group),
         obsCovs=list(visit=visitMat))

fm <- pcount(~1~x, umf, K=50)
checkTrue(inherits(fm, "unmarkedFitPCount"))

fmr <- pcount(~visit~x+(1|group), umf, K=50)

checkEqualsNumeric(coef(fmr), c(0.05535, 0.3200, -0.8795, 1.3638, 2.07098),
                   tol=1e-4)

checkTrue(inherits(sigma(fmr), 'data.frame'))
checkEquals(sigma(fmr)$sigma, 1.060223, tol=1e-5)

pr <- predict(fmr, "state")
checkEqualsNumeric(as.numeric(pr[1,]),
                   c(1.0385, 0.5827, 0.3457, 3.1192), tol=1e-4)

pr2 <- predict(fmr, "state", re.form=NA)
checkEqualsNumeric(as.numeric(pr2[1,]),
                   c(1.4862, 0.2019, 1.1387, 1.9396), tol=1e-4)

pr3 <- predict(fmr, "det")
checkTrue(inherits(pr3, "data.frame"))

nd <- data.frame(x=siteCovs(umf)$x[c(1,4)], group=factor(c(1,2)))
pr4 <- predict(fmr, "state", newdata=nd)
checkEqualsNumeric(pr4$Predicted, pr$Predicted[c(1,4)])

# New group level
nd <- data.frame(x=c(0,1), group=factor(101))
checkException(predict(fmr, "state", newdata=nd))

nd <- data.frame(x=c(0,1))
checkException(predict(fmr, "state", newdata=nd))

pr5 <- predict(fmr, "state", newdata=nd, re.form=NA)
checkTrue(inherits(pr5, "data.frame"))

ft <- fitted(fmr)
checkEqualsNumeric(dim(ft), c(300,3))

r <- ranef(fmr)
checkTrue(inherits(r, "unmarkedRanef"))
b <- bup(r)
checkTrue(cor(N, b) > 0.95)

rt <- randomTerms(fmr)
checkTrue(inherits(rt, "data.frame"))
checkEqualsNumeric(dim(rt), c(100,8))
checkTrue(cor(ran, rt$Estimate) > 0.8)


# Multiple random effects
umf2 <- umf
siteCovs(umf2)$id <- sample(letters[1:3], 300, replace=T)

fmr2 <- pcount(~1~x+(1|group)+(1|id), umf2, K=50)

checkTrue(nrow(sigma(fmr2))==2)
rt2 <- randomTerms(fmr2)
checkTrue(all(rt2$Groups==c(rep("group",100), rep("id",3))))

# Check other distributions
fmnb <- pcount(~1~1, umf, engine="TMB", mixture="NB", K=50)
checkTrue(inherits(fmnb@TMB, "list"))
checkTrue(all(names(fmnb@estimates@estimates)==c("state","det","alpha")))

fmzip <- pcount(~1~1, umf, engine="TMB", mixture="ZIP", K=50)
checkTrue(inherits(fmnb@TMB, "list"))
checkTrue(all(names(fmnb@estimates@estimates)==c("state","det","alpha")))

}
