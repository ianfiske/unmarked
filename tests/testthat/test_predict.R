

library(unmarked)
library(testthat)
import_functions("unmarked")


  y <- matrix(c(
      3, 2, 1, 4, 3, 2, 1, 4,
      3, 4, 2, 1, 3, 4, 2, 1,
      0, 1, 2, 3, 0, 1, 2, 3
      ), 3, 8, byrow=TRUE)
  siteCovs <- data.frame(sc1 = 1:3)
  obsCovs <- data.frame(oc1 = 1:24)
  ysc <- data.frame(ysc1 = 1:6, ysc2=factor(rep(c("a","b"), 3)))
  #ysc <- data.frame(ysc1 = 1:6, ysc2=factor(rep(c("a","b"), each=3)))
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
                          yearlySiteCovs=ysc, numPrimary=2)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10,
                    control=list(maxit=1))

  clean_up_covs(fm1, drop_final=TRUE)
  clean_up_covs(fm1)

  umf <- unmarkedFramePCO(y = y, siteCovs = NULL, obsCovs = NULL,
                          yearlySiteCovs=NULL, numPrimary=2)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10,
                    control=list(maxit=1))

  clean_up_covs(fm1, drop_final=TRUE)
  clean_up_covs(fm1)

  umf <- unmarkedFramePCO(y = y, siteCovs = NULL, obsCovs = NULL,
                          yearlySiteCovs=ysc, numPrimary=2)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10,
                    control=list(maxit=1))

  clean_up_covs(fm1, drop_final=TRUE)
  clean_up_covs(fm1)

  # Models with no yearlySiteCovs

  # Models with no obsCovs

  # Models with R!=J (multinomPois)
