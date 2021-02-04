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
    #return(Matrix::Matrix(matrix(0,nrow=nrow(data),ncol=0),sparse=TRUE))
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

get_b_vector <- function(TMB, type){
  sdr <- TMB::sdreport(TMB)
  bname <- paste0("b_",type)
  bpar <- sdr$par.random
  bpar <- bpar[grepl(bname, names(bpar))]
  if(length(bpar)==0) return(NULL)
  bpar
}

get_joint_cov <- function(TMB, type=NULL, remove_sigma=TRUE){
  full <- solve(TMB::sdreport(TMB, getJointPrecision=TRUE)$jointPrecision)
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

get_randvar_info <- function(names, estimates, covMat, formula, data){
  re <- get_reTrms(formula, data)
  list(names=names, estimates=estimates, covMat=covMat, fixed=1:length(estimates),
       invlink="exp", invlinkGrad="exp", n_obs=nrow(data),
       n_levels=lapply(re$flist, function(x) length(levels(x))), cnms=re$cnms)
}

print_randvar_info <- function(object){
  group_info <- paste0(names(object$n_levels), ", ",
                       unlist(object$n_levels), collapse="; ")

  val <- do.call(object$invlink, list(object$estimates))

  disp <- data.frame(Groups=names(object$cnms), Name=unlist(object$cnms),
                     Variance=round(val,3), Std.Dev.=round(sqrt(val),3))
  cat("Random effects:\n")
  print(disp, row.names=FALSE)
  cat(paste0("Number of obs: ",object$n_obs,", groups: ",group_info,"\n"))
}
