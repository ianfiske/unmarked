get_vars <- function(inp){
  if(is.list(inp)){
    out <- unique(unlist(lapply(inp, all.vars)))
  } else {
    out <- all.vars(inp)
  }
  names(out) <- out
  out
}

var_data <- function(var, guide, n){
  out <- rep(NA, n)
  gv <- guide[[var]]
  if(is.null(gv)){
    out <- rnorm(n, 0, 1)
  } else if(inherits(gv, "factor")){
    levs <- levels(gv)
    out <- factor(sample(levs, n, replace=TRUE), levels=levs)
  } else{
    gv$n <- n
    out <- do.call(gv$dist, gv[!names(gv)=="dist"])
  }
  out
}

generate_data <- function(formulas, guide, n){
  vars <- get_vars(formulas)
  if(length(vars)==0) return(NULL)
  as.data.frame(lapply(vars, var_data, guide=guide, n=n))
}

capitalize <- function(inp){
  paste0(toupper(substring(inp,1,1)),
         substring(inp,2,nchar(inp)))
}

parse_func_name <- function(inp){
  if(!is.character(inp)){
    stop("Argument must be a character string", call.=FALSE)
  }
  capitalize(inp)
}

blank_umFit <- function(fit_function){
  type <- parse_func_name(fit_function)
  type <- ifelse(type=="Pcount", "PCount", type)
  type <- ifelse(type=="MultinomPois", "MPois", type)
  type <- ifelse(type=="Distsamp", "DS", type)
  type <- ifelse(type=="Colext", "ColExt", type)
  type <- paste0("unmarkedFit", type)
  new(type)
}


setGeneric("get_umf_components", function(object, ...) standardGeneric("get_umf_components"))

setMethod("get_umf_components", "unmarkedFit",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$state, guide, design$M)
  oc <- generate_data(formulas$det, guide, design$J*design$M)
  yblank <- matrix(0, design$M, design$J)
  list(y=yblank, siteCovs=sc, obsCovs=oc)
})


setGeneric("simulate_fit", function(object, ...) standardGeneric("simulate_fit"))

setMethod("simulate_fit", "unmarkedFitOccu",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFrameOccu(y=parts$y, siteCovs=parts$siteCovs,
                           obsCovs=parts$obsCovs)
  occu(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
       data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("simulate_fit", "unmarkedFitPCount",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFramePCount(y=parts$y, siteCovs=parts$siteCovs,
                             obsCovs=parts$obsCovs)
  pcount(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
       data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("simulate_fit", "unmarkedFitOccuRN",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFrameOccu(y=parts$y, siteCovs=parts$siteCovs,
                           obsCovs=parts$obsCovs)
  occuRN(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
       data=umf, se=FALSE, control=list(maxit=1))
})


setMethod("simulate", "character",
  function(object, nsim=1, seed=NULL, formulas, coefs=NULL, design, guide=NULL, ...){
  model <- blank_umFit(object)
  fit <- suppressWarnings(simulate_fit(model, formulas, guide, design, ...))
  coefs <- check_coefs(coefs, fit)
  fit <- replace_estimates(fit, coefs)
  ysims <- simulate(fit, nsim)
  umf <- fit@data
  umfs <- lapply(ysims, function(x){
    umf@y <- x
    umf
  })
  if(length(umfs)==1) umfs <- umfs[[1]]
  umfs
})

setMethod("get_umf_components", "unmarkedFitMPois",
          function(object, formulas, guide, design, ...){
  args <- list(...)
  sc <- generate_data(formulas$state, guide, design$M)
  oc <- generate_data(formulas$det, guide, design$J*design$M)
  J <- ifelse(args$type=="double", 3, design$J)
  yblank <- matrix(0, design$M, J)
  list(y=yblank, siteCovs=sc, obsCovs=oc)
})

setMethod("simulate_fit", "unmarkedFitMPois",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFrameMPois(y=parts$y, siteCovs=parts$siteCovs,
                           obsCovs=parts$obsCovs, type=list(...)$type)
  multinomPois(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
               data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("get_umf_components", "unmarkedFitDS",
          function(object, formulas, guide, design, ...){
  #args <- list(...)
  sc <- generate_data(formulas$state, guide, design$M)
  sc2 <- generate_data(formulas$det, guide, design$M)
  sc <- cbind(sc,sc2)
  yblank <- matrix(0, design$M, design$J)
  list(y=yblank, siteCovs=sc)
})

setMethod("simulate_fit", "unmarkedFitDS",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(is.null(args$tlength)) args$tlength <- 0
  umf <- unmarkedFrameDS(y=parts$y, siteCovs=parts$siteCovs,
                         tlength=args$tlength, survey=args$survey, unitsIn=args$unitsIn,
                         dist.breaks=args$dist.breaks)
  keyfun <- ifelse(is.null(args$keyfun), "halfnorm", args$keyfun)
  output <- ifelse(is.null(args$output), "density", args$output)
  unitsOut <- ifelse(is.null(args$unitsOut), "ha", args$unitsOut)

  distsamp(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
               data=umf, se=FALSE, control=list(maxit=1), keyfun=keyfun,
          output=output, unitsOut=unitsOut)
})


setMethod("get_umf_components", "unmarkedFitColExt",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$psi, guide, design$M)
  ysc <- generate_data(list(formulas$col, formulas$ext), guide, design$M*design$T)
  oc <- generate_data(formulas$det, guide, design$J*design$M*design$T)
  yblank <- matrix(0, design$M, design$T*design$J)
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc, obsCovs=oc)
})


setMethod("simulate_fit", "unmarkedFitColExt",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedMultFrame(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           obsCovs=parts$obsCovs, numPrimary=design$T)
  colext(psiformula=formulas$psi, gammaformula=formulas$col,
         epsilonformula=formulas$ext,pformula=formulas$det,
         data=umf, se=FALSE, control=list(maxit=1))
})
