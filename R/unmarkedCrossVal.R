
setGeneric("crossVal", function(object, method=c("kfold","holdout","leaveoneout"), 
    folds=10, holdoutPct=0.25) standardGeneric("crossVal"))

setClass("unmarkedCrossVal",
    representation(stats = "data.frame",
                   method = "character",
                   folds = "numeric",
                   holdoutPct = "numeric"),
    validity=function(object){
      errors <- character(0)
      hp <- object@holdoutPct
      if(hp<0|hp>1){
        errors <- c(errors,"holdoutPct must be between 0 and 1") 
      }
    }
)

#Constructor of crossVal objects
setMethod("crossVal", "unmarkedFit", 
          function(object, method=c("kfold","holdout","leaveoneout"), 
                   folds=10, holdoutPct=0.25){
  
  method <- match.arg(method, c('kfold','holdout','leaveoneout'))

  if(method=="kfold" & !is.integer(folds) & folds < 0){
    stop("folds must be a positive integer")
  }
  if(method=="holdout" & (holdoutPct>1 | holdoutPct<0)){
    stop("holdoutPct must be a proportion between 0 and 1")
  }

  partitions <- switch(method,
    kfold = partitionKfold(object, folds=folds),
    holdout = partitionHoldout(object, holdoutPct=holdoutPct),
    leaveoneout = partitionLeaveOneOut(object)
  )
  
  n_reps <- length(partitions)
  rmse <- mae <- numeric(n_reps)
  for (i in 1:n_reps){
    
    newfit <- update(object, data=partitions[[i]]$trainData)

    newfit@data <- partitions[[i]]$testData
    if(!is.null(attributes(newfit)$knownOcc)){
      newfit@knownOcc <- rep(FALSE,numSites(newfit@data))
    }

    res <- residuals(newfit)
    if(is.list(res)){
      res <- unlist(res)
    }

    mae[i] <- mean(abs(res),na.rm=T)
    rmse[i] <- sqrt(mean(res^2,na.rm=T))

  }

  stats <- data.frame(rmse=rmse,mae=mae)

  out <- new("unmarkedCrossVal", stats=stats, method=method,
             folds=folds, holdoutPct=holdoutPct)

  out
})

#Kfold partition function
partitionKfold <- function(object, folds){

  site_inds <- 1:numSites(object@data)
  shuf_site_inds <- sample(site_inds,numSites(object@data))
  fold_inds <- cut(site_inds, breaks=folds, labels=FALSE)

  fold_list <- vector(length=folds,"list")
  for (i in 1:folds){

    trainInds <- shuf_site_inds[fold_inds!=i]
    testInds <- shuf_site_inds[fold_inds==i]
  
    fold_list[[i]]$trainData <- object@data[trainInds,]
    fold_list[[i]]$testData <-  object@data[testInds,]
  }
  fold_list
}

#Holdout partition function
partitionHoldout <- function(object, holdoutPct){
  
  site_inds <- 1:numSites(object@data)
  shuf_site_inds <- sample(site_inds,numSites(object@data))

  splitInd <- round(numSites(object@data)*(1-holdoutPct))
  trainInds <- shuf_site_inds[1:splitInd]
  testInds <- shuf_site_inds[(splitInd+1):length(shuf_site_inds)]

  fold_list <- vector(length=1,"list")
  fold_list[[1]]$trainData <- object@data[trainInds,]
  fold_list[[1]]$testData <- object@data[testInds,]

  fold_list
}

#leave-one-out
partitionLeaveOneOut <- function(object){
  
  fold_list <- vector(length=numSites(object@data),"list")
  for (i in seq_along(fold_list)){
    fold_list[[i]]$trainData <- object@data[-i,]
    fold_list[[i]]$testData <- object@data[i,]
  }
  fold_list

}

setMethod("show", "unmarkedCrossVal", function(object)
{
  st <- object@stats

  if(object@method=='kfold'){
    cat(paste('Method: k-fold (',object@folds,' folds)\n\n',sep=''))
  } else if(object@method=='holdout'){
    cat(paste('Method: holdout (',round(object@holdoutPct*100),'% in test set)\n\n',sep=''))
  } else if(object@method=='leaveoneout'){
    cat('Method: leave-one-out\n\n')
  }

  cat('Root Mean Square Error:\n')
  print(data.frame(Estimate=mean(st$rmse),SD=sd(st$rmse)),row.names=FALSE,
        digits=3)
  cat('\n')
  cat('Mean Absolute Error:\n')
  print(data.frame(Estimate=mean(st$mae),SD=sd(st$mae)), row.names=FALSE,
        digits=3)

})

setClass("unmarkedCrossValList",
    representation(stats_list="list",
                   method = "character",
                   folds="numeric",
                   holdoutPct="numeric")
)

#CrossVal list constructor
setMethod("crossVal", "unmarkedFitList",
          function(object, method=c("kfold","holdout","leaveoneout"), 
                   folds=10, holdoutPct=0.25){

    stats <- lapply(object@fits, crossVal, method, folds, holdoutPct)

    out <- new("unmarkedCrossValList", stats_list=stats, method=method, folds=folds,
               holdoutPct=holdoutPct)

})


setMethod("show", "unmarkedCrossValList", function(object){

  sl <- object@stats_list
  nfits <- length(sl)

  if(object@method=='kfold'){
    cat(paste('Method: k-fold (',object@folds,' folds)\n\n',sep=''))
  } else if(object@method=='holdout'){
    cat(paste('Method: holdout (',round(object@holdoutPct*100),'% in test set)\n\n',sep=''))
  } else if(object@method=='leaveoneout'){
    cat('Method: leave-one-out\n\n')
  }

  rmse_means <- rmse_sds <- numeric(nfits)
  mae_means <- mae_sds <- numeric(nfits)
  mod_names <- names(sl)
  for (i in 1:nfits){
    rmse_means[i] <- mean(sl[[i]]@stats$rmse)
    rmse_sds[i] <- sd(sl[[i]]@stats$rmse)
    mae_means[i] <- mean(sl[[i]]@stats$mae)
    mae_sds[i] <- sd(sl[[i]]@stats$mae)
  }
  out_rmse <- data.frame(Estimate=rmse_means,SD=rmse_sds,
                         row.names=mod_names)
  out_mae <- data.frame(Estimate=mae_means,SD=mae_sds,
                        row.names=mod_names)


  cat('Root Mean Square Error:\n')
  print(out_rmse[order(out_rmse$Estimate),],digits=3)
  cat('\n')
  cat('Mean Absolute Error:\n')
  print(out_mae[order(out_mae$Estimate),],digits=3)

})
