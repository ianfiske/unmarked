



test.crossVal.occu <- function() {

    set.seed(123)
    data(frogs)
    pferUMF <- unmarkedFrameOccu(pfer.bin)
    siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))    
    obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
     

    fm <- occu(~ obsvar1 ~ 1, pferUMF)

    kfold <- crossVal(fm, method='Kfold', folds=10)
    
    checkEqualsNumeric(nrow(kfold@stats),10)

    checkEqualsNumeric(as.numeric(kfold@stats[1,]), 
                       c(0.3213177,0.2159953), tolerance=1e-4)

    holdout <- crossVal(fm, method='holdout', holdoutPct=0.25)

    checkEqualsNumeric(as.numeric(holdout@stats[1,]), 
                       c(0.3695291,0.2414929), tolerance=1e-4)

    leave <- crossVal(fm, method='leaveOneOut')
    
    checkEqualsNumeric(nrow(leave@stats),130)
    checkEqualsNumeric(as.numeric(leave@stats[1,]),
                       c(0.6220418,0.4985790), tolerance=1e-4)

    #Check parallel
    set.seed(123)
    kfold <- crossVal(fm, method='Kfold', folds=10)
    set.seed(123)
    kfold_par <- crossVal(fm, method='Kfold', folds=10, parallel=TRUE)
    checkEqualsNumeric(kfold@stats, kfold_par@stats)

    #Check custom stat function
    checkException(crossVal(fm, statistic=function(x) "fake"))

    new_stat <- function(object){
      c(mean_res = mean(residuals(object),na.rm=T))
    }

    kfold_custom <- crossVal(fm, statistic=new_stat)
    checkEqualsNumeric(kfold_custom@stats[,1], rep(0,10), tol=0.05)

}

test.crossValList <- function(){

    set.seed(123)
    data(frogs)
    pferUMF <- unmarkedFrameOccu(pfer.bin)
    siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))    
    obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
     

    fm <- occu(~ obsvar1 ~ 1, pferUMF)
    fm2 <- occu(~1 ~1, pferUMF)

    fl <- fitList(fm2=fm2,fm=fm)

    cvlist <- crossVal(fl, method='Kfold')
    
    checkEqualsNumeric(length(cvlist@stats_list),2)
}


test.crossVal.multinomPois <- function(){

    set.seed(123)
    data(ovendata)
    ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")
    fm1 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)

    mout <- crossVal(fm1, method='Kfold')
    checkEqualsNumeric(as.numeric(mout@stats[1,]),
                       c(0.5521100,0.3335076), tolerance=1e-4)

}
