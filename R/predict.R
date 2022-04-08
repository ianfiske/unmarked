# General predict method-------------------------------------------------------

# Common predict function for all fit types
setGeneric("predict2", function(object,...) standardGeneric("predict2"))
setMethod("predict2", "unmarkedFit",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){

  # If no newdata, get actual data
  if(missing(newdata) || is.null(newdata)) newdata <- object@data

  # Check inputs
  check_predict_arguments(object, type, newdata)

  # Get model matrix (X) and offset
  # newdata is an unmarkedFrame, use getDesign
  is_raster <- FALSE
  if(inherits(newdata, "unmarkedFrame")){
    pred_inps <- predict_inputs_from_umf(object, type, newdata, na.rm, re.form)
  } else {
    # 1. Get original data and formula
    orig_data <- get_orig_data(object, type)
    orig_formula <- get_formula(object, type)

    # 2. If raster, get newdata from raster as data.frame
    if(inherits(newdata, c("RasterLayer","RasterStack"))){
      if(!require(raster)) stop("raster package required", call.=FALSE)
      is_raster <- TRUE
      orig_raster <- newdata
      newdata <- newdata_from_raster(newdata, all.vars(orig_formula))
    }

    # 3. Make model matrix and offset with newdata, informed by original data
    pred_inps <- make_mod_matrix(orig_formula, orig_data, newdata, re.form)
  }

  # Calculate predicted values in chunks (for speed) based on X and offset
  out <- predict_by_chunk(object, type, level, pred_inps$X, pred_inps$offset,
                          chunk_size = 70, backTransform, re.form)

  # Convert output to raster if newdata was raster
  if(is_raster){
    out <- raster_from_predict(out, orig_raster, appendData)
  } else if(appendData){
    # Append data if needed
    out <- data.frame(out, as(newdata, "data.frame"))
  }

  out

})

#Utility function to make model matrix and offset from newdata
make_mod_matrix <- function(formula, data, newdata, re.form=NULL){
  form_nobars <- lme4::nobars(formula)
  mf <- model.frame(form_nobars, data, na.action=stats::na.pass)
  X.terms <- stats::terms(mf)
  fac_cols <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_cols, levels)
  xlevs <- xlevs[names(xlevs) %in% names(mf)]
  nmf <- model.frame(X.terms, newdata, na.action=stats::na.pass, xlev=xlevs)
  #X <- model.matrix(X.terms, newdata, xlev=xlevs)
  X <- model.matrix(form_nobars, nmf)
  offset <- model.offset(nmf)
  if(is.null(re.form) & !is.null(lme4::findbars(formula))){
    Z <- get_Z(formula, data, newdata)
    X <- cbind(X, Z)
  }
  list(X=X, offset=offset)
}


# Fit-type specific methods----------------------------------------------------

# Fit type-specific methods to generate different components of prediction
# 1. check_predict_arguments(): Check arguments
# 2. predict_inputs_from_umf(): Generating inputs from an unmarked
#    frame (e.g. when no newdata),
# 3. get_formula: Get formula for submodel type
# 4. get_orig_data(): Get original dataset for use in building model frame
# 5. predict_by_chunk(): Take inputs and generate predictions
# Basic methods are shown below; fit-specific methods in their own sections

setGeneric("check_predict_arguments", function(object, ...){
  standardGeneric("check_predict_arguments")
})

setMethod("check_predict_arguments", "unmarkedFit",
  function(object, type, newdata, ...){
  # Check type
  check_type(object, type)

  # Check newdata
  if(!inherits(newdata, c("unmarkedFrame", "data.frame", "RasterLayer", "RasterStack"))){
    stop("newdata must be unmarkedFrame, data.frame, RasterLayer, or RasterStack", call.=FALSE)
  }
  invisible(TRUE)
})

# Check if predict type is valid
check_type <- function(mod, type){
  opts <- names(mod)
  if(type %in% opts) return(invisible(TRUE))
  stop("Valid types are ", paste(opts, collapse=", "), call.=FALSE)
}


# Get X and offset when newdata is umf
setGeneric("predict_inputs_from_umf", function(object, ...){
  standardGeneric("predict_inputs_from_umf")
})

setMethod("predict_inputs_from_umf", "unmarkedFit",
  function(object, type, newdata, na.rm, re.form){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  if(type == "state") list_els <- c("X","Z_state","X.offset")
  if(type == "det") list_els <- c("V","Z_det","V.offset")

  X <- designMats[[list_els[1]]]
  if(is.null(re.form)) X <- cbind(X, designMats[[list_els[2]]])
  offset <- designMats[[list_els[3]]]

  list(X=X, offset=offset)
})

# Get correct individual formula based on type
setGeneric("get_formula", function(object, type, ...){
  standardGeneric("get_formula")
})

setMethod("get_formula", "unmarkedFit", function(object, type, ...){
  if(type == "state") return(as.formula(paste("~", object@formula[3], sep="")))
  if(type == "det") return(as.formula(object@formula[[2]]))
  stop("Invalid type")
})

# When newdata is data.frame/raster, get original dataset
# For use in building correct model frame
setGeneric("get_orig_data", function(object, type, ...){
  standardGeneric("get_orig_data")
})

setMethod("get_orig_data", "unmarkedFit", function(object, type, ...){
  clean_covs <- clean_up_covs(object, drop_factor_levels=TRUE)
  datatype <- switch(type, state='site_covs', det='obs_covs')
  clean_covs[[datatype]]
})

# Deal with NULL covariate data frames and add site covs to ysc,
# ysc to obs covs, etc.
clean_up_covs <- function(object, drop_factor_levels=FALSE){
  M <- numSites(object@data)
  R <- ncol(object@data@y)
  T <- 1
  J <- R
  is_mult <- methods::.hasSlot(object@data, "numPrimary")
  if(is_mult){
    T <- object@data@numPrimary
    J <- R/T
  }

  sc <- siteCovs(object@data)
  if(is.null(sc)) sc <- data.frame(.dummy=rep(1,M))
  out <- list(site_covs=sc)

  if(is_mult){
    ysc <- yearlySiteCovs(object@data)
    if(is.null(ysc)) ysc <- data.frame(.dummy2=rep(1,M*T))
    if(!is.null(sc)) ysc <- cbind(ysc, sc[rep(1:M, each=T),,drop=FALSE])
  }

  if(methods::.hasSlot(object@data, "obsCovs")){
    oc <- obsCovs(object@data)
    if(is.null(oc)) oc <- data.frame(.dummy3=rep(1,M*T*J))
    if(is_mult){
      oc <- cbind(oc, ysc[rep(1:(M*T), each=J),,drop=FALSE])
    }
    out$obs_covs=oc
  }

  if(is_mult & (T > 1)){
    if(T > 1 & drop_factor_levels){
      # Drop factor levels only found in last year of data
      ysc <- droplevels_final_year(ysc, M, T)
    }
    out$yearly_site_covs <- ysc
  }

  out
}

#Remove data in final year of yearlySiteCovs
#then drop factor levels found only in that year
droplevels_final_year <- function(dat, nsites, nprimary){
  dat[seq(nprimary, nsites*nprimary, by=nprimary), ] <- NA
  dat <- lapply(dat, function(x) x[,drop = TRUE])
  as.data.frame(dat)
}


# Take inputs and generate prediction, done in chunks for speed
setGeneric("predict_by_chunk", function(object, ...){
  standardGeneric("predict_by_chunk")
})

setMethod("predict_by_chunk", "unmarkedFit",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){

  if(is.vector(xmat)) xmat <- matrix(xmat, nrow=1)
  nr <- nrow(xmat)
  ind <- rep(1:ceiling(nr/chunk_size), each=chunk_size, length.out=nr)

  # should find a way to keep xmat sparse, but it doesn't
  # work with linearComb
  #x_chunk <- lapply(split(as.data.frame(xmat), ind), as.matrix)
  x_chunk <- lapply(unique(ind),
                    function(i) as.matrix(xmat[ind==i,,drop=FALSE]))

  if(is.null(offsets)) offsets <- rep(0, nr)
  off_chunk <- split(offsets, ind)
  out <- mapply(function(x_i, off_i){
    has_na <- apply(x_i, 1, function(x_i) any(is.na(x_i)))
    x_i[has_na,] <- 0
    off_i[has_na] <- 0
    lc <- linearComb(object, x_i, type, offset=off_i, re.form=re.form)
    if(backTransform) lc <- backTransform(lc)
    se <- SE(lc)
    ci <- confint(lc, level=level)
    out <- data.frame(Predicted=coef(lc), SE=se,
                      lower=ci[,1], upper=ci[,2])
    out[has_na,] <- NA
    out
    }, x_chunk, off_chunk, SIMPLIFY=FALSE)
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
})

# Raster handling functions----------------------------------------------------


# Convert a raster into a data frame to use as newdata
newdata_from_raster <- function(rst, vars){
  nd <- raster::as.data.frame(rst)
  # Handle factors
  is_fac <- raster::is.factor(rst)
  rem_string <- paste(paste0("^",names(rst),"_"), collapse="|")
  names(nd)[is_fac] <- gsub(rem_string, "", names(nd)[is_fac])
  # Check if variables are missing
  no_match <- vars[! vars %in% names(nd)]
  if(length(no_match) > 0){
    stop(paste0("Variable(s) ",paste(no_match, collapse=", "), " not found in raster stack"),
         call.=FALSE)
  }
  return(nd)
}

# Convert predict output into a raster
raster_from_predict <- function(pr, orig_rst, appendData){
  new_rast <- data.frame(raster::coordinates(orig_rst), pr)
  new_rast <- raster::stack(raster::rasterFromXYZ(new_rast))
  crs(new_rast) <- raster::crs(orig_rst)
  if(appendData) new_rast <- raster::stack(new_rast, orig_rst)
  new_rast
}


# pcount methods---------------------------------------------------------------

setMethod("check_predict_arguments", "unmarkedFitPCount",
  function(object, type, newdata, ...){
  if(type %in% c("psi", "alpha")){
    stop(paste0(type, " is scalar. Use backTransform instead."), call.=FALSE)
  }
  methods::callNextMethod(object, type, newdata)
})

# Special predict stuff for ZIP distribution in pcount
setMethod("predict_by_chunk", "unmarkedFitPCount",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "state" & object@mixture == "ZIP"){

    out <- data.frame(matrix(NA, nrow(xmat), 4))
    names(out) <- c("Predicted", "SE", "lower", "upper")

    psi.hat <- plogis(coef(object, type="psi"))
    lamEst <- object["state"]
    psiEst <- object["psi"]
    fixedOnly <- !is.null(re.form)
    lam.mle <- coef(lamEst, fixedOnly=fixedOnly)
    lam_vcov <- vcov(lamEst, fixedOnly=fixedOnly)

    for(i in 1:nrow(xmat)) {
      if(nrow(xmat) > 5000) {
          if(i %% 1000 == 0)
              cat("  doing row", i, "of", nrow(xmat), "\n")
      }
      if(any(is.na(xmat[i,]))) next
      ## for the ZIP model the predicted values on the log scale have us
      ## add log(1-psi.hat) to the normal linear prediction
      out$Predicted[i] <-   xmat[i,] %*% lam.mle + offsets[i] + log(1 - psi.hat)
      ## to compute the approximate SE, I compute the variance of the usual
      ## linear part -- that is easy, and to that I add the variance of
      ## log(1-psi.hat) obtained by the delta approximation
      logit.psi<-coef(object,type="psi")
      #  To do that I took derivative of log(1-psi.hat) using application
      #  of chain rule.... hopefully correctly.
      delta.approx.2ndpart<-   ( ((1/(1-psi.hat))*(exp(logit.psi)/((1+exp(logit.psi))^2)))^2 ) * (SE(psiEst)^2)
      ## now the SE is the sqrt of the whole thing
      out$SE[i]<- sqrt( t(xmat[i,])%*% lam_vcov %*%xmat[i,] + delta.approx.2ndpart   )

      #From Mike Meredith
      alf <- (1 - level) / 2
      crit<-qnorm(c(alf, 1 - alf))
      ci <- out$Predicted[i] + crit * out$SE[i]
      out$lower[i]<- ci[1]
      out$upper[i]<- ci[2]
      if(backTransform){
        out$Predicted[i] <- exp(out$Predicted[i])
        ### If back-transform, delta approx says var = (exp(linear.predictor)^2)*Var(linear.predictor)
        ### also I exponentiate the confidence interval.....
        out$SE[i]<- out$Predicted[i]*out$SE[i]
        ci<-exp(ci)
      }
      out$lower[i] <- ci[1]
      out$upper[i] <- ci[2]
    }
    return(out)
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})


# colext methods---------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitColExt",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  list_el <- switch(type, psi="W", col="X.gam", ext="X.eps", det="V")
  # colext doesn't support offsets
  list(X=designMats[[list_el]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitColExt", function(object, type, ...){
  switch(type, psi=object@psiformula, col=object@gamformula,
                     ext=object@epsformula, det=object@detformula)
})

setMethod("get_orig_data", "unmarkedFitColExt", function(object, type, ...){
  clean_covs <- clean_up_covs(object)
  datatype <- switch(type, psi='site_covs', col='yearly_site_covs',
                     ext='yearly_site_covs', det='obs_covs')
  clean_covs[[datatype]]
})


# occuFP methods---------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitOccuFP",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@detformula, object@FPformula,
                          object@Bformula, object@stateformula, na.rm=na.rm)
  X_idx <- switch(type, state="X", det="V", fp="U", b="W")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=designMats[[off_idx]])
})

setMethod("get_formula", "unmarkedFitOccuFP", function(object, type, ...){
  switch(type, state=object@stateformula, det=object@detformula,
         b=object@Bformula, fp=object@FPformula)
})

setMethod("get_orig_data", "unmarkedFitOccuFP", function(object, type, ...){
  # Get obs data if fp, b, or det
  new_type <- ifelse(type %in% c("fp", "b"), "det", type)
  methods::callNextMethod(object, new_type, ...)
})


# Dail-Madsen model methods----------------------------------------------------

# Having problems with class unions, so defined regular functions common
# to all D-M models and then wrapped them in methods
# MMO inherits from PCO, so no need to write unique MMO methods here

check_predict_arg_dm <- function(object, type, newdata, ...){
  if(type %in% c("psi", "alpha", "scale")){
    stop(paste0(type, " is scalar. Use backTransform instead."), call.=FALSE)
  }
  dynamics <- object@dynamics
  immigration <- tryCatch(object@immigration, error=function(e) FALSE)
  if(identical(dynamics, "notrend") & identical(type, "gamma"))
    stop("gamma is a derived parameter for this model: (1-omega)*lambda")
  if(identical(dynamics, "trend") && identical(type, "omega"))
    stop("omega is not a parameter in the dynamics='trend' model")
  if(!immigration && identical(type, "iota"))
    stop("iota is not a parameter in the immigration=FALSE model")
}

setMethod("check_predict_arguments", "unmarkedFitPCO",
  function(object, type, newdata, ...){
  check_predict_arg_dm(object, type, newdata, ...)
  methods::callNextMethod(object, type, newdata)
})

setMethod("check_predict_arguments", "unmarkedFitDSO",
  function(object, type, newdata, ...){
  check_predict_arg_dm(object, type, newdata, ...)
  methods::callNextMethod(object, type, newdata)
})

pred_inp_umf_dm <- function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm=na.rm)
  X_idx <- switch(type, lambda="Xlam", gamma="Xgam", omega="Xom",
                  iota="Xiota", det="Xdet")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=designMats[[off_idx]])
}

setMethod("predict_inputs_from_umf", "unmarkedFitPCO",
  function(object, type, newdata, na.rm, re.form=NA){
  pred_inp_umf_dm(object, type, newdata, na.rm, re.form=NA)
})

setMethod("predict_inputs_from_umf", "unmarkedFitDSO",
  function(object, type, newdata, na.rm, re.form=NA){
  pred_inp_umf_dm(object, type, newdata, na.rm, re.form=NA)
})

get_formula_dm <- function(object, type, ...){
  fl <- object@formlist
  switch(type, lambda=fl$lambdaformula, gamma=fl$gammaformula,
         omega=fl$omegaformula, iota=fl$iotaformula, det=fl$detformula)
}

setMethod("get_formula", "unmarkedFitPCO", function(object, type, ...){
  get_formula_dm(object, type, ...)
})

setMethod("get_formula", "unmarkedFitDSO", function(object, type, ...){
  get_formula_dm(object, type, ...)
})

# This method differs between PCO and DSO
setMethod("get_orig_data", "unmarkedFitPCO", function(object, type, ...){
  clean_covs <- clean_up_covs(object)
  datatype <- switch(type, lambda='site_covs', gamma='yearly_site_covs',
                     omega='yearly_site_covs', iota='yearly_site_covs',
                     det='obs_covs')
  clean_covs[[datatype]]
})

setMethod("get_orig_data", "unmarkedFitDSO", function(object, type, ...){
  clean_covs <- clean_up_covs(object)
  datatype <- switch(type, lambda='site_covs', gamma='yearly_site_covs',
                     omega='yearly_site_covs', iota='yearly_site_covs',
                     det='yearly_site_covs')
  clean_covs[[datatype]]
})

# Special handling for ZIP distribution
pred_zip_dm <- function(object, type, level, xmat, offsets, chunk_size,
                        backTransform, re.form, ...){
  warning("Method to compute SE for ZIP model has not been written", call.=FALSE)
  out <- data.frame(matrix(NA, nrow(xmat), 4))
  names(out) <- c("Predicted", "SE", "lower", "upper")
  lam.mle <- coef(object, type="lambda")
  psi.hat <- plogis(coef(object, type="psi"))
  if(is.null(offsets)) offsets <- rep(0, nrow(xmat))
  out$Predicted <- as.numeric(xmat %*% lam.mle + offsets + log(1 - psi.hat))
  if(backTransform) out$Predicted <- exp(out$Predicted)
  out
}

setMethod("predict_by_chunk", "unmarkedFitPCO",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "lambda" & object@mixture == "ZIP"){
    return(pred_zip_dm(object, type, level, xmat, offsets, chunk_size, backTransform,
                  re.form, ...))
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})

setMethod("predict_by_chunk", "unmarkedFitDSO",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "lambda" & object@mixture == "ZIP"){
    return(pred_zip_dm(object, type, level, xmat, offsets, chunk_size, backTransform,
                  re.form, ...))
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})


# Temporary emigration models--------------------------------------------------

# All inherit from GMM so only one set of methods is required
# (except GDR which has its own predict method right now)

setMethod("predict_inputs_from_umf", "unmarkedFitGMM",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm=na.rm)
  X_idx <- switch(type, lambda="Xlam", phi="Xphi", det="Xdet")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=designMats[[off_idx]])
})

setMethod("get_formula", "unmarkedFitGMM", function(object, type, ...){
  fl <- object@formlist
  switch(type, lambda=fl$lambdaformula, phi=fl$phiformula, det=fl$pformula)
})

setMethod("get_orig_data", "unmarkedFitGMM", function(object, type, ...){
  clean_covs <- clean_up_covs(object, drop_factor_levels=FALSE)
  datatype <- switch(type, lambda='site_covs', phi='yearly_site_covs',
                     det='obs_covs')
  clean_covs[[datatype]]
})


# occuTTD----------------------------------------------------------------------

# Identical to colext

setMethod("predict_inputs_from_umf", "unmarkedFitOccuTTD",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  list_el <- switch(type, psi="W", col="X.gam", ext="X.eps", det="V")
  list(X=designMats[[list_el]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitOccuTTD", function(object, type, ...){
  switch(type, psi=object@psiformula, col=object@gamformula,
                     ext=object@epsformula, det=object@detformula)
})

setMethod("get_orig_data", "unmarkedFitOccuTTD", function(object, type, ...){
  clean_covs <- clean_up_covs(object)
  datatype <- switch(type, psi='site_covs', col='yearly_site_covs',
                     ext='yearly_site_covs', det='obs_covs')
  clean_covs[[datatype]]
})


# nmixTTD----------------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitNmixTTD",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  list_el <- switch(type, state="W", det="V")
  list(X=designMats[[list_el]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitNmixTTD", function(object, type, ...){
  switch(type, state=object@stateformula, det=object@detformula)
})

setMethod("get_orig_data", "unmarkedFitNmixTTD", function(object, type, ...){
  clean_covs <- clean_up_covs(object, drop_factor_levels=FALSE)
  datatype <- switch(type, state='site_covs', det='obs_covs')
  clean_covs[[datatype]]
})

