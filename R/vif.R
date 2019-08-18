#Variance inflation factors for unmarked models
#Calculated using var-cov matrix as with car::vif()

vif <- function(mod, type){

  v <- vcov(mod, type)
  type_options <- names(mod@estimates)
  if(!type%in%type_options){
    stop(paste("Possible types are",paste(type_options, collapse=', ')),
        call.=FALSE) 
  }
  est <- mod@estimates[type]@estimates
  labs <- names(est)
  ints <- grep("(Intercept)", labs)
  if(length(ints)>0){
    v <- v[-ints,-ints]
    labs <- labs[-ints]
  } else{
    stop("No intercept", call.=FALSE)
  }
  n.terms <- length(labs)
  if (n.terms < 2) stop("model contains fewer than 2 terms", call.=FALSE)

  R <- stats::cov2cor(v)
  detR <- det(R)

  result <- numeric(n.terms)
  names(result) <- labs
  for (i in 1:n.terms) {
    result[i] <- det(as.matrix(R[i, i])) * 
                    det(as.matrix(R[-i, -i])) / detR
  }
  result

}
