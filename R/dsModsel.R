# This is kind of ugly. It would probably be better to manage estimates and variances in seperate matrices and then make a print method that puts them togehter.

modSel.mixmodList <- function(fits, nullmod=NULL, response=NULL) 
{
classes <- table(sapply(fits, function(x) class(x)[1]))
if(length(names(classes)) > 1) 
	warning("More than 1 class of model in model in fits")
estList <- lapply(fits, coef, altNames=T)
vcovList <- lapply(fits, vcov, altNames=T)
vars <- sort(unique(unlist(sapply(estList, names))))
Vvars <- paste("VAR", vars, sep="")
varvars <- c(); l <- length(c(vars, Vvars))
varvars[seq(1, l, by=2)] <- vars
varvars[seq(2, l, by=2)] <- Vvars
cNames <- c("Response", "lamForm", "pForm", varvars, "Code", "CondNum", 
	"Deviance", "K", "n", "AICc", "AICcWt", "deltaAICc", "devExplnd", "NagR2")
sptable <- table(sapply(fits, function(x) substr(colnames(x$y)[1], 1, 4)))
if(length(names(sptable)) > 1) 
	warning("More than 1 response variable in list of models")
out <- data.frame(matrix(NA, ncol=length(cNames), nrow=length(fits)))
rownames(out) <- names(fits)
colnames(out) <- cNames
out$Response <- ifelse(is.null(response), toupper(names(sptable)), response)
out$lamForm <- sapply(fits, function(x) paste(x$stateformula, collapse=""))
out$pForm <- sapply(fits, function(x) paste(x$detformula, collapse=""))
for(i in 1:length(vars)) {
	out[,vars[i]] <- sapply(estList, function(x) x[vars[i]])
	out[,Vvars[i]] <- sapply(vcovList, function(x) diag(x)[vars[i]])
	}
out$Code <- sapply(fits, "[[", "convergence")
out$CondNum <- sapply(fits, cn)
out$Deviance <- 2 * sapply(fits, "[[", "value")
out$K <- sapply(fits, function(x) length(coef(x)))
out$n <- sapply(fits, "[[", "n")
out$AICc <- sapply(fits, aic.c)
if(is.null(nullmod) & !is.null(fits$Null)) {
	out$devExplnd <- 1 - out$Deviance / out["Null", "Deviance"]
	out$NagR2 <- sapply(fits, nag.r2, fits$Null) 
	}
else {
	out$devExplnd <- 1 - out$Deviance / (2 * nullmod$value)
	out$NagR2 <- sapply(fits, nag.r2, nullmod) 	
	}
out$deltaAICc <- out$AICc-min(out$AICc)
out$AICcWt <- exp(-out$deltaAICc/2)
out$AICcWt <- out$AICcWt/sum(out$AICcWt)
out <- out[order(out$AICc),]
out$AICcWtCum <- cumsum(out$AICcWt)
return(out)
}



