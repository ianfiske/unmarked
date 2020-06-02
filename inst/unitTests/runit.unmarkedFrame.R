test.emptyframe <- function() {
  checkException(umf <- unmarkedFrame())
}

test.frame <- function() {
  M <- 10
  J <- 3
  y <- matrix(rbinom(J * M, 1, 0.5), M, J)
  siteCovs <- data.frame(a = rnorm(M), b = factor(gl(2,5)))
  umf <- unmarkedFrame(y = y, siteCovs = siteCovs)
}


test.umfDS.args <- function() {
    y <- matrix(1, 1, 2)
    s <- "point"
    d <- 0:2
    uin <- "m"
    sc <- data.frame(1)
    oc <- matrix(1, 2)
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, obsCovs=oc,
        survey=s, dist.breaks=d, unitsIn=uin))
    umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, dist.breaks=d,
        unitsIn=uin)
    checkException(obsCovs(umf) <- oc)
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        dist.breaks=d))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, dist.breaks=d,
        unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        dist.breaks=0:3, unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s,
        dist.breaks=1:3, unitsIn=uin))

    }



test.obsToY <- function() {
    y <- matrix(c(
        1, 0, 0,
        2, 1, 0,
        1, 0, 1,
        2, 1, 2,
        1, 0, 3,
        1, 1, 1), nrow=6, ncol=3, byrow=TRUE)
    oc <- matrix(c(
        1, 0,
        2, 1,
        1, 1,
        NA, 0,
        1, NA,
        NA, NA), nrow=6, ncol=2, byrow=TRUE)
    umf <- unmarkedFrameMPois(y = y, obsCovs = list(x=oc), type="double")
    o2y <- obsToY(umf)

    checkEquals(o2y, matrix(1, 2, 3))
    oc.na <- is.na(oc)
    oc.na %*% o2y

    }


test.umf.yearlySiteCovs <- function() {

  n <- 50   # number of sites
  T <- 4    # number of primary periods
  J <- 3    # number of secondary periods
     
  site <- 1:50
  years <- data.frame(matrix(rep(2010:2013, each=n), n, T))
  years <- data.frame(lapply(years, as.factor))
  dummy <- matrix(rep(c('a','b','c','d'),n),nrow=n,byrow=T)
  occasions <- data.frame(matrix(rep(1:(J*T), each=n), n, J*T))   
  y <- matrix(0:1, n, J*T)
     
  umf <- unmarkedMultFrame(y=y,
        siteCovs = data.frame(site=site),
        obsCovs=list(occasion=occasions),
        yearlySiteCovs=list(year=years,dummy=dummy),
        numPrimary=T)

  as_df <- as(umf,'data.frame')

  checkEqualsNumeric(dim(as_df),c(50,33))
  checkTrue(all(names(as_df)[13:22] == c('site','year.1','year.2','year.3',
                                     'year.4','dummy.1','dummy.2','dummy.3',
                                     'dummy.4','occasion.1')))
  checkTrue(all(as_df$year.1==2010))
  checkTrue(all(as_df$dummy.1=='a'))


  umf2 <- unmarkedMultFrame(y=y,
        siteCovs = data.frame(site=site),
        obsCovs=list(occasion=occasions),
        numPrimary=T)

  as_df2 <- as(umf2,'data.frame')

  checkEqualsNumeric(dim(as_df2),c(50,25))

}

test.umf.char.to.factor <- function(){

  n <- 50   # number of sites
  T <- 4    # number of primary periods
  J <- 3    # number of secondary periods

  y <- matrix(0:1, n, J*T)

  #Site covs
  sc <- data.frame(x=rnorm(n), y=sample(letters, 50, replace=TRUE))
  checkEquals(sapply(sc, class), c(x="numeric", y="character"))
  
  options(warn=2)
  checkException(umf <- unmarkedFrame(y, siteCovs=sc))
  options(warn=0)
  umf <- unmarkedFrame(y, siteCovs=sc)
  checkEquals(sapply(siteCovs(umf), class), c(x="numeric", y="factor"))

  #Already factor
  sc2 <- data.frame(x=rnorm(n), y=factor(sample(letters, 50, replace=TRUE)))
  umf <- unmarkedFrame(y, siteCovs=sc2)
  checkEquals(sapply(siteCovs(umf), class), c(x="numeric", y="factor"))

  #Obs covs
  oc <- data.frame(x=rnorm(n*J*T), y=sample(letters, n*J*T, replace=TRUE))
  checkEquals(sapply(oc, class), c(x="numeric", y="character"))
  
  options(warn=2)
  checkException(umf <- unmarkedFrame(y, obsCovs=oc))
  options(warn=0)
  umf <- unmarkedFrame(y, obsCovs=oc)
  checkEquals(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  checkTrue(is.null(siteCovs(umf)))

  #as list
  oc <- list(x=matrix(oc$x, nrow=n), y=matrix(oc$y, nrow=n))
  options(warn=2)
  checkException(umf <- unmarkedFrameOccu(y, obsCovs=oc))
  options(warn=0)
  umf <- unmarkedFrameOccu(y, obsCovs=oc)
  checkEquals(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  checkTrue(is.null(siteCovs(umf)))
  
  #Check conversion
  df <- as(umf, "data.frame")
  checkEqualsNumeric(dim(df), c(50,36))

  #Yearly site covs
  ysc <- list(x=matrix(rnorm(n*T), nrow=n), 
             y=matrix(sample(letters, n*T, replace=TRUE), nrow=n))
  options(warn=2)
  checkException(umf <- unmarkedMultFrame(y, yearlySiteCovs=ysc, numPrimary=T))
  options(warn=0)
  umf <- unmarkedMultFrame(y, yearlySiteCovs=ysc, numPrimary=T)
  checkEquals(sapply(yearlySiteCovs(umf), class), c(x="numeric", y="factor"))
  checkTrue(is.null(siteCovs(umf)))

  #All
  options(warn=2)
  checkException(umf <- unmarkedMultFrame(y, yearlySiteCovs=ysc, obsCovs=oc, 
                                          siteCovs=sc, numPrimary=T))
  options(warn=0)
  umf <- unmarkedMultFrame(y, yearlySiteCovs=ysc, obsCovs=oc, 
                          siteCovs=sc, numPrimary=T)
  checkEquals(sapply(yearlySiteCovs(umf), class), c(x="numeric", y="factor"))
  checkEquals(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  checkEquals(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  
  df <- as(umf, "data.frame")
  checkEqualsNumeric(dim(df), c(50,46)) 
}
