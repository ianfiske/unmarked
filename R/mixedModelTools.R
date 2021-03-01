model_frame <- function(formula, data, newdata=NULL){

  formula <- lme4::nobars(formula)
  mf <- model.frame(formula, data, na.action=stats::na.pass)

  if(is.null(newdata)){
    return(mf)
  }

  check_newdata(newdata, formula)
  model.frame(stats::terms(mf), newdata, na.action=stats::na.pass,
                        xlev=get_xlev(data, mf))
}

model_matrix <- function(formula, data, newdata=NULL){
  mf <- model_frame(formula, data, newdata)
  model.matrix(lme4::nobars(formula), mf)
}

model_offset <- function(formula, data, newdata=NULL){
  mf <- model_frame(formula, data, newdata)
  out <- model.offset(mf)
  if(is.null(out)) out <- rep(0, nrow(mf))
  out
}

get_xlev <- function(data, model_frame){
  fac_col <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_col, levels)
  xlevs[names(xlevs) %in% names(model_frame)]
}

get_reTrms <- function(formula, data, newdata=NULL){
  fb <- lme4::findbars(formula)
  mf <- model.frame(lme4::subbars(formula), data, na.action=stats::na.pass)
  if(is.null(newdata)) return(lme4::mkReTrms(fb, mf))
  new_mf <- model.frame(stats::terms(mf), newdata, na.action=stats::na.pass,
                        xlev=get_xlev(data, mf))
  lme4::mkReTrms(fb, new_mf, drop.unused.levels=FALSE)
}

get_Z <- function(formula, data, newdata=NULL){
  if(is.null(lme4::findbars(formula))){
    if(is.null(newdata)){
      return(Matrix::Matrix(matrix(0, nrow=nrow(data), ncol=0),sparse=TRUE))
    } else{
      return(Matrix::Matrix(matrix(0, nrow=nrow(newdata), ncol=0),sparse=TRUE))
    }
  }
  check_formula(formula, data)
  Zt <- get_reTrms(formula, data, newdata)$Zt
  Z <- t(as.matrix(Zt))
  Matrix::Matrix(Z, sparse=TRUE)
}

get_group_vars <- function(formula){
  rand <- lme4::findbars(formula)
  ifelse(is.null(rand), 0, length(rand))
}

get_nrandom <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(as.array(0))

  out <- sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    length(unique(data[[col_nm]]))
  })
  as.array(out)
}

has_random <- function(formula){
  length(lme4::findbars(formula)) > 0
}

sigma_names <- function(formula, data){
  if(!has_random(formula)) return(NA_character_)
  nms <- get_reTrms(formula, data)$cnms
  nms <- paste0(unlist(nms), "|", names(nms))
  nms <- gsub("(Intercept)", "1", nms, fixed=TRUE)
  #paste0("sigma [", nms, "]")
  nms
}

check_formula <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(invisible())

  char <- paste(deparse(formula))
  if(grepl(":|/", char)){
    stop("Nested random effects (using / and :) are not supported",
         call.=FALSE)
  }
  theta <- get_reTrms(formula, data)$theta
  if(0 %in% theta){
    stop("Correlated slopes and intercepts are not supported. Use || instead of |.",
         call.=FALSE)
  }
}

check_newdata <- function(newdata, formula){
  inp_vars <- names(newdata)
  term_vars <- all.vars(formula)
  not_found <- ! term_vars %in% inp_vars
  if(any(not_found)){
    stop(paste0("Required variables not found in newdata: ",
               paste(term_vars[not_found], collapse=", ")), call.=FALSE)
  }
}

split_formula <- function(formula){
  if(length(formula) != 3) stop("Double right-hand side formula required")
  p1 <- as.formula(formula[[2]])
  p2 <- as.formula(paste0(formula[[1]], deparse(formula[[3]])))
  list(p1, p2)
}

nobars_double <- function(form){
  spl <- split_formula(form)
  spl <- lapply(spl, lme4::nobars)
  spl <- paste(unlist(lapply(spl, as.character)),collapse="")
  as.formula(spl)
}

is_tmb_fit <- function(mod){
  if(!methods::.hasSlot(mod, "TMB")) return(FALSE)
  !is.null(mod@TMB)
}

get_b_vector <- function(tmb_report, type){
  #sdr <- TMB::sdreport(TMB)
  bname <- paste0("b_",type)
  bpar <- tmb_report$par.random
  bpar <- bpar[grepl(bname, names(bpar))]
  if(length(bpar)==0) return(NULL)
  bpar
}

get_joint_cov <- function(tmb_report, type=NULL, remove_sigma=TRUE){
  full <- tmb_report$jointPrecision
  if(is.null(full)){
    full <- tmb_report$cov.fixed
    colnames(full) <- rownames(full) <- names(tmb_report$par)
  } else {
    full <- solve(full)
  }

  if(is.null(type)) return(full)
  keep <- grepl(paste0("_",type), colnames(full))
  out <- full[keep,keep,drop=FALSE]
  if(!remove_sigma) return(out)
  keep2 <- !grepl(paste0("lsigma_",type), colnames(out))
  out <- out[keep2,keep2,drop=FALSE]
}

tmbfit_has_random <- function(mod, type){
  paste0(type,"RE") %in% names(mod@estimates@estimates)
}

use_tmb_bootstrap <- function(mod, type, re.form){
  is.null(re.form) && is_tmb_fit(mod) && tmbfit_has_random(mod, type)
}

# Gather information about grouping variables for a given submodel
get_randvar_info <- function(tmb_report, type, formula, data){
  ngv <- get_group_vars(formula)
  if(ngv == 0) return(list()) #Return blank list if there are no grouping variables

  sigma_type <- paste0("lsigma_",type)
  sigma_ind <- grepl(sigma_type, get_fixed_names(tmb_report))
  sigma_est <- tmb_report$par.fixed[sigma_ind]
  sigma_cov <- as.matrix(tmb_report$cov.fixed[sigma_ind,sigma_ind])
  re <- get_reTrms(formula, data)

  list(names=sigma_names(formula, data), estimates=sigma_est, covMat=sigma_cov,
       invlink="exp", invlinkGrad="exp", n_obs=nrow(data),
       n_levels=lapply(re$flist, function(x) length(levels(x))), cnms=re$cnms,
       levels=rownames(re$Zt))
}

get_fixed_names <- function(tmb_report){
  out <- names(tmb_report$par.fixed)
  if(is.null(out)) out <- 1:length(tmb_report$par.fixed)
  out
}

print_randvar_info <- function(object){
  group_info <- paste0(names(object$n_levels), ", ",
                       unlist(object$n_levels), collapse="; ")

  val <- do.call(object$invlink, list(object$estimates))

  disp <- data.frame(Groups=names(object$cnms), Name=unlist(object$cnms),
                     Variance=round(val^2,3), Std.Dev.=round(val,3))
  cat("Random effects:\n")
  print(disp, row.names=FALSE)
  #below needs to be corrected for missing values at some point
  #cat(paste0("Number of obs: ",object$n_obs,", groups: ",group_info,"\n"))
}

check_no_support <- function(formula_list){
  has_bars <- any(sapply(formula_list, function(x) !is.null(lme4::findbars(x))))
  if(has_bars){
    stop("This function does not support random effects", call.=FALSE)
  }
}

fit_TMB <- function(model, data, params, random,
                    starts, method, ...){

  fixed_sub <- names(params)[!names(params) %in% random]
  nfixed <- length(unlist(params[fixed_sub]))

  tmb_mod <- TMB::MakeADFun(data = c(model = model, data),
                            parameters = params,
                            random = random,
                            silent=TRUE,
                            DLL = "unmarked_TMBExports")

  if(is.null(starts)) starts <- rep(0, nfixed)
  if(length(starts) != nfixed){
    stop(paste("The number of starting values should be", nfixed))
  }

  opt <- optim(starts, fn=tmb_mod$fn, gr=tmb_mod$gr, method=method, ...)

  sdr <- TMB::sdreport(tmb_mod, getJointPrecision=TRUE)
  sdr$par <- tmb_mod$par

  AIC = 2 * opt$value + 2 * nfixed

  list(opt=opt, TMB=tmb_mod, sdr=sdr, AIC=AIC)
}

get_coef_info <- function(tmb_report, type, names, idx){
  no_sigma <- !grepl("lsigma", get_fixed_names(tmb_report))
  fixed <- tmb_report$par.fixed[no_sigma] #take out sigmas
  fixed <- fixed[idx]
  names(fixed) <- names
  rand <- get_b_vector(tmb_report, type)
  ests <- c(fixed, rand)
  covMat <- get_joint_cov(tmb_report, type)
  list(ests=ests, cov=covMat)
}

setMethod("sigma", "unmarkedEstimate", function(object, level=0.95, ...){
  rinf <- object@randomVarInfo
  if(length(rinf)==0){
    stop("No random effects in this submodel", call.=FALSE)
  }
  z <- qnorm((1-level)/2, lower.tail = FALSE)
  vals <- rinf$estimates
  ses <- sqrt(diag(rinf$covMat))
  lower <- vals - z*ses
  upper <- vals + z*ses
  Groups <- names(rinf$cnms)
  Name <- unlist(rinf$cnms)
  data.frame(Model=object@short.name, Groups=Groups, Name=Name, sigma=exp(vals),
             lower=exp(lower), upper=exp(upper))
})

setMethod("sigma", "unmarkedFit", function(object, type, level=0.95, ...){
  if(!missing(type)){
    return(sigma(object[type], level=level))
  }
  ests <- object@estimates@estimates
  has_rand <- sapply(ests, function(x) length(x@randomVarInfo)>0)
  if(!any(has_rand)){
    stop("No random effects in this model", call.=FALSE)
  }
  ests <- ests[has_rand]

  out_list <- lapply(ests, sigma, level=level)
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
})


setGeneric("randomTerms", function(object, ...) standardGeneric("randomTerms"))


setMethod("randomTerms", "unmarkedEstimate", function(object, level=0.95, ...){

  rv <- object@randomVarInfo
  if(length(rv)==0){
    stop("No random effects in this submodel", call.=FALSE)
  }

  Groups <- lapply(1:length(rv$cnms), function(x){
                  gn <- names(rv$cnms)[x]
                  rep(gn, rv$n_levels[[gn]])
          })
  Groups <- do.call(c, Groups)

  Name <-  lapply(1:length(rv$cnms), function(x){
                  gn <- names(rv$cnms)[x]
                  var <- rv$cnms[[x]]
                  rep(var, rv$n_levels[[gn]])
          })
  Name <- do.call(c, Name)

  rv_idx <- !1:length(object@estimates) %in% object@fixed
  b_var <- object@estimates[rv_idx]
  b_se <- sqrt(diag(object@covMat[rv_idx,rv_idx,drop=FALSE]))

  z <- qnorm((1-level)/2, lower.tail = FALSE)
  lower <- b_var - z*b_se
  upper <- b_var + z*b_se

  data.frame(Model=object@short.name, Groups=Groups, Name=Name, Level=rv$levels,
                    Estimate=b_var, SE=b_se, lower=lower, upper=upper)
})


setMethod("randomTerms", "unmarkedFit", function(object, type, level=0.95, ...){

  if(!missing(type)){
    return(randomTerms(object[type], level))
  }

  has_random <- sapply(object@estimates@estimates,
                       function(x) length(x@randomVarInfo) > 0)
  if(!any(has_random)){
    stop("No random effects in this model", call.=FALSE)
  }
  keep <- object@estimates@estimates[has_random]
  out <- lapply(keep, randomTerms, level=level)
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
})
