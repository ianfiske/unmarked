
context("predict-related functions")

test_that("clean_up_covs works with dynamic model data",{

  # Dynamic data
  y <- matrix(c(
      3, 2, 1, 4,
      3, 4, 2, 1,
      0, 1, 2, 3
      ), 3, 4, byrow=TRUE)
  siteCovs <- data.frame(sc1 = 1:3)
  obsCovs <- data.frame(oc1 = 1:12)
  ysc <- data.frame(ysc1 = 1:6, ysc2=factor(rep(c("a","b"), 3)))
  #ysc <- data.frame(ysc1 = 1:6, ysc2=factor(rep(c("a","b"), each=3)))
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
                          yearlySiteCovs=ysc, numPrimary=2)

  dr <- unmarked:::clean_up_covs(umf, drop_final=TRUE)
  expect_equal(dr$site_covs, data.frame(sc1=1:3)) 
  expect_equal(dr$yearly_site_covs, data.frame(ysc1=c(1,NA,3,NA,5,NA),
                                                 ysc2=factor(c("a",NA,"a",NA,"a",NA)),
                                                 sc1=c(1,NA,2,NA,3,NA)))
  expect_equivalent(dr$obs_covs, data.frame(oc1=1:12, ysc1=rep(1:6, each=2),
                                       ysc2=factor(rep(c("a","b"), each=2)),
                                       sc1=rep(1:3, each=4)))

  no_drop <- unmarked:::clean_up_covs(umf)
  expect_equivalent(no_drop$yearly_site_covs, data.frame(ysc1=1:6,
                                                    ysc2=factor(rep(c("a","b"),3)),
                                                    sc1=rep(1:3, each=2)))

  umf <- unmarkedFramePCO(y=y, numPrimary=2)
  
  cc <- unmarked:::clean_up_covs(umf, drop_final=TRUE)
  expect_equivalent(cc$obs_covs,
                    data.frame(.dummy3=rep(1,12), .dummy2=rep(1,12), .dummy1=rep(1,12)))
})

test_that("clean_up_covs works with single-season models",{
  y <- matrix(c(0,1,0,0,1,1), nrow=3)
  umf <- unmarkedFrameOccu(y=y, siteCovs=data.frame(sc1=1:3),
                           obsCovs=data.frame(oc1=1:6))
  cc <- unmarked:::clean_up_covs(umf)
  expect_equal(names(cc), c("site_covs","obs_covs"))
  expect_equivalent(cc$site_covs, data.frame(sc1=1:3))
  expect_equivalent(cc$obs_covs, data.frame(oc1=1:6, sc1=rep(1:3, each=2)))
  cc2 <- unmarked:::clean_up_covs(umf, drop_final=TRUE)
  expect_equal(cc, cc2)
})

test_that("clean_up_covs works with models with no obs covs",{
  # single season
  ltUMF <- with(linetran, {
    unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
      siteCovs = data.frame(Length, area, habitat),
      dist.breaks = c(0, 5, 10, 15, 20),
      tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
    }) 
  ltUMF

  cc <- unmarked:::clean_up_covs(ltUMF)  
  expect_equal(names(cc), c("site_covs", "obs_covs"))
  expect_equal(dim(cc$obs_covs), c(12,4))
})

test_that("clean_up_covs works with models where length(y) != length(p)",{
  # double observer, etc
  nSites <- 3
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }
 
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  expect_warning(umf <- unmarkedFrameMPois(y=y,
         siteCovs <- data.frame(sc1=1:3),
         obsCovs=list(observer=observer),
         type="double"))

  cc <- unmarked:::clean_up_covs(umf)
  expect_equivalent(cc$site_covs, data.frame(sc=1:3))
  expect_equivalent(cc$obs_covs, data.frame(observer=factor(c(rep(c("A","B"), 3))),
                                            sc1=rep(1:3, each=2)))

})

