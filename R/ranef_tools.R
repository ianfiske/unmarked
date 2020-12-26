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
    return(Matrix::Matrix(matrix(0,0,0),sparse=TRUE))
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
  spl <- unmarked:::split_formula(form)
  spl <- lapply(spl, lme4::nobars)
  spl <- paste(unlist(lapply(spl, as.character)),collapse="")
  as.formula(spl)
}

is_tmb_fit <- function(mod){
  if(!methods::.hasSlot(mod, "TMB")) return(FALSE)
  !is.null(mod@TMB)
}

get_b_vector <- function(mod, type){
  stopifnot(is_tmb_fit(mod))
  sdr <- TMB::sdreport(mod@TMB)
  bname <- paste0("b_",type)
  bpar <- sdr$par.random
  bpar <- bpar[grepl(bname, names(bpar))]
  if(length(bpar)==0) return(NULL)
  bpar
}

get_joint_cov <- function(mod, type=NULL){
  stopifnot(is_tmb_fit(mod))
  full <- solve(TMB::sdreport(mod@TMB, getJointPrecision=TRUE)$jointPrecision)
  if(is.null(type)) return(full)
  keep <- grepl(paste0("_",type), colnames(full))
  full[keep,keep]
}

tmbfit_has_random <- function(mod, type){
  paste0(type,"RE") %in% names(mod@estimates@estimates)
}

use_tmb_bootstrap <- function(mod, type, re.form){
  is.null(re.form) && is_tmb_fit(mod) && tmbfit_has_random(mod, type)
}

tmb_predict_bootstrap <- function(mod, type, newdata=NULL, backTransform=TRUE,
                                  level=0.95, nsim=100){
  cat("Bootstrapping confidence intervals with", nsim, "samples\n")
  qlow <- (1-level)/2
  qup <- level+qlow
  invlink <- mod@estimates@estimates[[type]]@invlink
  if(invlink=="logistic") invlink <- "plogis"

  beta <- coef(mod, type)
  b <- get_b_vector(mod, type)

  spl_forms <- split_formula(mod@formula)
  if(type=="state"){
    form <- spl_forms[[2]]
    dat <- siteCovs(getData(mod))
  } else if(type == "det"){
    form <- spl_forms[[1]]
    dat <- obsCovs(getData(mod))
  }

  X <- model_matrix(form, dat, newdata)
  Xoffset <- model_offset(form, dat, newdata)
  Z <- get_Z(form, dat, newdata)
  stopifnot(nrow(X)==nrow(Z))
  est <- as.vector(X %*% beta + Z %*% b + Xoffset)

  covMat <- get_joint_cov(mod, type)
  keep <- grepl("b_|beta_", colnames(covMat))
  covMat <- covMat[keep,keep]
  mu <- c(beta, b)
  samp <- MASS::mvrnorm(nsim, mu, covMat)
  which_beta <- grepl("beta_", colnames(covMat))
  samp_beta <- samp[,which_beta,drop=FALSE]
  samp_b <- samp[,!which_beta,drop=FALSE]

  lp_sim <- matrix(NA, nrow(X), nsim)
  for (i in 1:nsim){
    lp_sim[,i] <- as.vector(X %*% samp_beta[i,] + Z %*% samp_b[i,])
  }

  ci <- t(apply(lp_sim, 1, quantile, c(qlow,qup)))
  colnames(ci) <- c("lower","upper")

  if(backTransform){
    est <- do.call(invlink, list(est))
    ci <- do.call(invlink, list(ci))
    lp_sim <- do.call(invlink, list(lp_sim))
  }
  se <- apply(lp_sim, 1, sd)

  data.frame(Predicted=est, SE=se, ci)
}
