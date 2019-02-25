



test.crossVal.occu <- function() {

    set.seed(123)
    data(frogs)
    pferUMF <- unmarkedFrameOccu(pfer.bin)
    siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))    
    obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
     

    fm <- occu(~ obsvar1 ~ 1, pferUMF)

    kfold <- crossVal(fm, method='kfold', folds=10)
    
    checkEqualsNumeric(nrow(kfold@stats),10)

    checkEqualsNumeric(as.numeric(kfold@stats[1,]), 
                       c(0.289997,0.1984301), tolerance=1e-4)

    holdout <- crossVal(fm, method='holdout', holdoutPct=0.25)

    checkEqualsNumeric(as.numeric(holdout@stats[1,]), 
                       c(0.3313081,0.2210848), tolerance=1e-4)

    leave <- crossVal(fm, method='leaveoneout')
    
    checkEqualsNumeric(nrow(leave@stats),130)
    checkEqualsNumeric(as.numeric(leave@stats[1,]),
                       c(0.6220418,0.4985790), tolerance=1e-4)

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

    cvlist <- crossVal(fl, method='kfold')
    
    checkEqualsNumeric(length(cvlist@stats_list),2)
}


test.crossVal.multinomPois <- function(){

    set.seed(123)
    data(ovendata)
    ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")
    fm1 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)

    mout <- crossVal(fm1, method='kfold')
    checkEqualsNumeric(as.numeric(mout@stats[1,]),
                       c(0.5007878,0.3587654), tolerance=1e-4)

}
