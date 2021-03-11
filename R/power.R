setClass("unmarkedPower",
  representation(call="call", data="unmarkedFrame", M="numeric",
                 J="numeric", T="numeric", coefs="list", estimates="list",
                 alpha="numeric")
)

powerAnalysis <- function(object, coefs, design=NULL, alpha=0.05, nsim=100,
                          parallel=FALSE){

  stopifnot(inherits(object, "unmarkedFit"))

  submodels <- names(object@estimates@estimates)
  coefs <- check_coefs(coefs, object)
  fit_temp <- replace_estimates(object, coefs)

  T <- 1
  bdata <- NULL
  if(is.null(design)){
    sims <- simulate(fit_temp, nsim)
    M <- numSites(object@data)
    if(methods::.hasSlot(object@data, "numPrimary")){
      T <- object@data@numPrimary
    }
    J <- obsNum(object@data) / T
  } else {
    bdata <- bootstrap_data(fit_temp@data, nsim, design)
    sims <- lapply(bdata, function(x){
      fit_temp@data <- x
      #temporary workaround
      if(methods::.hasSlot(fit_temp, "knownOcc")){
        fit_temp@knownOcc <- rep(FALSE, design$M)
      }
      simulate(fit_temp, 1)[[1]]
    })
    M <- design$M
    J <- design$J
  }

  cl <- NULL
  if(parallel){
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(unmarked))
  }


  fits <- pbapply::pblapply(1:nsim, function(i, sims, fit, bdata=NULL){
    if(!is.null(design)) fit@data <- bdata[[i]]
    fit@data@y <- sims[[i]]
    update(fit, data=fit@data, se=TRUE)
  }, sims=sims, fit=object, bdata=bdata, cl=cl)

  sum_dfs <- lapply(fits, get_summary_df)

  new("unmarkedPower", call=object@call, data=object@data, M=M,
      J=J, T=T, coefs=coefs, estimates=sum_dfs, alpha=alpha)
}

bootstrap_data <- function(data, nsims, design){
  M <- design$M
  J <- design$J
  sites <- 1:numSites(data)
  if(!is.null(J) & methods::.hasSlot(data, "numPrimary")){
    stop("Can't automatically bootstrap observations with > 1 primary period", call.=FALSE)
  }
  if(J > obsNum(data)){
    stop("Can't currently bootstrap more than the actual number of observations", call.=FALSE)
  }
  obs <- 1:obsNum(data)

  if(M > numSites(data)){
    M_samps <- lapply(1:nsims, function(i) sample(sites, M, replace=TRUE))
  } else if(M < numSites(data)){
    M_samps <- lapply(1:nsims, function(i) sample(sites, M, replace=FALSE))
  } else {
    M_samps <- replicate(nsims, sites, simplify=FALSE)
  }

  if(J > obsNum(data)){
    J_samps <- lapply(1:nsims, function(i) sample(obs, J, replace=TRUE))
  } else if(J < obsNum(data)){
    J_samps <- lapply(1:nsims, function(i) sample(obs, J, replace=FALSE))
  } else {
    J_samps <- replicate(nsims, obs, simplify=FALSE)
  }

  lapply(1:nsims, function(i) data[M_samps[[i]], J_samps[[i]]])
}


check_coefs <- function(coefs, fit){
  required_coefs <- names(fit@estimates@estimates)
  required_lens <- sapply(fit@estimates@estimates, function(x) length(x@estimates))
  for (i in 1:length(required_coefs)){
    if(!required_coefs[i] %in% names(coefs)){
      stop(paste0("Missing entry '",required_coefs[i], "' in coefs list"), call.=FALSE)
    }
    if(length(coefs[[required_coefs[i]]]) != required_lens[i]){
      stop(paste0("Entry '",required_coefs[i], "' in coefs list must be length ",
                  required_lens[i]), call.=FALSE)
    }
  }
  coefs[required_coefs]
}

setMethod("summary", "unmarkedPower", function(object, ...){
  sum_dfs <- object@estimates
  npar <- nrow(sum_dfs[[1]])

  pow <- sapply(1:npar, function(ind){
    pcrit <- sapply(sum_dfs, function(x) x$`P(>|z|)`[ind]) < object@alpha
    direct <- sapply(sum_dfs, function(x) x$Estimate[ind]) * unlist(object@coefs)[ind]  > 0
    mean(pcrit & direct, na.rm=T)
  })

  out <- cbind(sum_dfs[[1]][,1:2], effect=unlist(object@coefs), power=pow)
  out <- out[out$param != "(Intercept)",,drop=FALSE]
  rownames(out) <- NULL
  names(out) <- c("Submodel", "Parameter", "Effect", "Power")
  out
})

setMethod("show", "unmarkedPower", function(object){
  cat(paste0("\nModel: ", deparse(object@call), "\n\n"))
  cat("Power Statistics:\n")
  print(summary(object), row.names=FALSE)
})

replace_estimates <- function(object, new_ests){
  for (i in 1:length(new_ests)){
    est <- object@estimates@estimates[[names(new_ests)[i]]]@estimates
    stopifnot(length(est) == length(new_ests[[i]]))
    object@estimates@estimates[[names(new_ests)[i]]]@estimates <- new_ests[[i]]
  }
  object
}

get_summary_df <- function(fit){
  n_est <- length(fit@estimates@estimates)
  #est_names <- unname(sapply(fit@estimates@estimates, function(x) x@name))
  est_names <- names(fit@estimates@estimates)
  all_est <- lapply(1:n_est, function(i){
    capture.output(out <- summary(fit@estimates@estimates[[i]]))
    out <- cbind(submodel=est_names[i], param=rownames(out), out)
    rownames(out) <- NULL
    out
  })
  do.call(rbind, all_est)
}

setClass("unmarkedPowerList", representation(powerAnalyses="list"))

setGeneric("unmarkedPowerList", function(object, ...){
             standardGeneric("unmarkedPowerList")})

setMethod("unmarkedPowerList", "list", function(object, ...){
  new("unmarkedPowerList", powerAnalyses=object)
})

setMethod("unmarkedPowerList", "unmarkedFit",
  function(object, coefs, design, alpha=0.05, nsim=100, parallel=FALSE, ...){

  ndesigns <- nrow(design)
  out <- lapply(1:ndesigns, function(i){
    cat(paste0("M = ",design$M[i],", J = ",scenarios$J[i],"\n"))
    powerAnalysis(object, coefs, as.list(design[i,]), alpha=alpha, nsim=nsim,
                  parallel=FALSE)
  })
  unmarkedPowerList(out)
})

setMethod("summary", "unmarkedPowerList", function(object, ...){
  out <- lapply(object@powerAnalyses, function(x){
    stats <- summary(x)
    cbind(M=x@M, T=x@T, J=x@J, stats)
  })
  out <- do.call(rbind, out)
  out$M <- factor(out$M)
  out$T <- factor(out$T)
  out$J <- factor(out$J)
  out
})

setMethod("show", "unmarkedPowerList", function(object){
  print(summary(object))
})

setMethod("plot", "unmarkedPowerList", function(x, beta=NULL, param=NULL, ...){
  dat <- summary(x)
  if(is.null(param)) param <- dat$Parameter[1]
  dat <- dat[dat$Parameter==param,,drop=FALSE]
  ylim <- range(dat$Power, na.rm=T)
  if(!is.null(beta)) ylim[2] <- max(1-beta, ylim[2])
  xlim <- range(as.numeric(as.character(dat$M)), na.rm=T)
  cols <- palette.colors(length(levels(dat$J)), palette="Dark 2")
  old_par <- par()[c("mfrow","mar")]
  nT <- length(levels(dat$T))
  mar <- old_par$mar
  if(nT == 1) mar <- c(5.1, 4.1, 2.1, 2.1)
  par(mfrow=c(length(levels(dat$T)),1), mar=mar)
  for (i in levels(dat$T)){
    plot_title <- ""
    if(nT > 1) plot_title <- paste0("T = ", i)
    tsub <- dat[dat$T==i,,drop=FALSE]
    Jlev <- levels(tsub$J)
    jsub <- tsub[tsub$J==Jlev[1],,drop=FALSE]
    plot(as.numeric(as.character(jsub$M)), jsub$Power, type="o",
        col=cols[1], ylim=ylim, xlim=xlim, xlab="Sites",
        ylab="Power", pch=19, main=plot_title)
    if(!is.null(beta)) abline(h=1-beta, lty=2)
    for (j in 2:length(Jlev)){
      jsub <- tsub[tsub$J==Jlev[j],,drop=FALSE]
      lines(as.numeric(as.character(jsub$M)), jsub$Power, type="o",
            col=cols[j], pch=19)
    }
    legend('bottomright', lwd=1, pch=19, col=cols, legend=Jlev, title="Observations")
  }
  par(mfrow=old_par)
})
