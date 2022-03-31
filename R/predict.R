
# Check if predict type is valid
check_type <- function(mod, type){
  opts <- names(mod)
  if(type %in% opts) return(invisible(TRUE))
  stop("Valid types are ", paste(opts, collapse=", "), call.=FALSE)
}

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

# General predict method - used by occu, occuRN, ...?
setMethod("predict", "unmarkedFit",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){

  # Check type
  check_type(object, type)

  # If no newdata, get actual data
  if(missing(newdata) || is.null(newdata)) newdata <- getData(object)

  # Check newdata
  if(!inherits(newdata, c("unmarkedFrame", "data.frame", "RasterStack"))){
    stop("newdata must be unmarkedFrame, data.frame, or RasterStack", call.=FALSE)
  }

  # Get model matrix (X) and offset
  # newdata is an unmarkedFrame, use getDesign
  is_raster <- FALSE
  if(inherits(newdata, "unmarkedFrame")){
    designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
    switch(type,
      state = {
        X <- designMats$X
        if(is.null(re.form)) X <- cbind(X, designMats$Z_state)
        offset <- designMats$X.offset
      },
      det = {
        X <- designMats$V
        if(is.null(re.form)) X <- cbind(X, designMats$Z_det)
        offset <- designMats$V.offset
      })

  } else {
    # When newdata is data.frame or RasterStack we need to do more work
    # 1. Get original data
    M <- numSites(getData(object))
    J <- obsNum(getData(object))
    sc <- siteCovs(getData(object))
    if(type == "state"){
      origform <- as.formula(paste("~", object@formula[3], sep=""))
      origdata <- sc
      if(is.null(origdata)) origdata <- data.frame(.dummy = rep(1, M))
    } else if(type == "det"){
      origform <- as.formula(object@formula[[2]])
      origdata <- obsCovs(getData(object))
      if(is.null(origdata)) origdata <- data.frame(.dummy = rep(1, M*J))
      # Add site covariates if necessary
      if(!is.null(sc)){
        origdata <- cbind(origdata, sc[rep(1:M, each=J),,drop=FALSE])
      }
    }

    # 2. If raster, get newdata from raster as data.frame
    if(inherits(newdata, "RasterStack")){
      if(!require(raster)) stop("raster package required", call.=FALSE)
      is_raster <- TRUE
      orig_raster <- newdata
      newdata <- newdata_from_raster(newdata, all.vars(origform))
    }

    # 3. Make model matrix and offset with newdata, informed by original data
    mm <- make_mod_matrix(origform, origdata, newdata, re.form)
    X <- mm$X
    offset <- mm$offset
  }

  # Calculate predicted values in chunks (for speed) based on X and offset
  out <- predict_by_chunk(object, type, level, X, offset, chunk_size = 70,
                          backTransform, re.form)

  # Convert output to raster if newdata was raster
  if(is_raster){
    out <- raster_from_predict(out, orig_raster, appendData)
  } else if(appendData){
    # Append data if needed
    out <- data.frame(out, as(newdata, "data.frame"))
  }

  out

})
