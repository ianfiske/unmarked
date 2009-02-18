#' @include classes.R
#' @import reshape
#' @import roxygen
roxygen()

genFixedNLL <- function(nll, whichFixed, fixedValues) {
  function(params) {
    params[whichFixed] <- fixedValues
    do.call(nll, list(params=params))
  }
}

# nll the original negative log likelihood function
# MLE the full vector of MLE values
profileCI <- function(nll, whichPar, MLE, interval){
  MLEnll <- nll(MLE)
  nPar <- length(MLE)
  f <- function(value) {
    fixedNLL <- genFixedNLL(nll, whichPar, value)
    mleRestricted <- optim(rep(0,nPar), fixedNLL)
    MLEnll - mleRestricted$value + 1.92
  }
## add some kind of try/catch block here for when interval is on boundary.
## first bnd should be (est-5*se, est+5*se)... catch this and expand?
  lower <- uniroot(f, c(interval[1],MLE[whichPar]))
  upper <- uniroot(f, c(MLE[whichPar], interval[2]))
  return(c(lower$root,upper$root))
}

## link functions and their gradients
logistic <- function(x) {
  1/(1 + exp(-x))
}

logistic.grad <- function(x) {
  exp(-x)/(exp(-x)+1)^2
}

log.grad <- function(x) { # duh! (but for clarity)
  1/x
}

identLink <- function(x) x
identLinkGrad <- function(x) 1

## use logarithms to vectorize row-wise products
## this speeds things up a LOT (vs. apply(x,1,prod))
rowProds <-
function(x, na.rm = FALSE)
{
  exp(rowSums(log(x), na.rm = na.rm))
}

# helper function to coerce an array of matrices to a list
arrToList <- function(x){
  nl <- list()
  for(i in 1:dim(x)[3]) {
    nl[[i]] <- x[,,i]
  }
  names(nl) <- dimnames(x)[[3]]
  nl
}

# compute estimated asymptotic variances of parameter estimates
# using the observed information matrix
sd.est <- function(fm) {
    sqrt(diag(solve(fm$hessian)))
}

# delta method for variance of proportion given variance of its logistic-
# transformed counterpart
sd.prop <- function(est,sd.est) {
    exp(-est)/(1 + exp(-est))^2 * sd.est
}

### track linked list of parameters using a data frame
### add row to linked list
addParm <- function(list.df, parm.name, parm.length) {
    if(parm.length > 0) {
        if(nrow(list.df) == 0) {
            last.ind <- 0
        } else {
            last.ind <- list.df$end[nrow(list.df)]
        }
        parm.df <- data.frame(parameter = parm.name, start = last.ind + 1,
                              end = last.ind + parm.length,
                              stringsAsFactors = FALSE)
        list.df <- rbind(list.df, parm.df)
    }
    return(list.df)
}

parmNames <- function(list.df) {
    npar <- list.df$end[nrow(list.df)]
    names <- character(npar)
    for(i in 1:npar) {
        which.par <- which(i >= list.df$start & i <= list.df$end)
        names[i] <- list.df$parameter[which.par]
    }
    return(names)
}


#' This function converts an appropriatedly formated comma-separated
#' values file (.csv) to a format usable by \emph{unmarked}'s fitting
#' functions (see \emph{Details}).
#'
#'   This function provides a quick way to take a .csv file with headers
#' named as described below and provides the data required and returns of
#' data in the format required by the model-fitting functions in
#' \code{\link{unmarked}}.  The .csv file can be in one of 2 formats: long or
#' wide.  See the first 2 lines of the \emph{examples} for what these
#' formats look like.
#'
#' The .csv file is formatted as follows:
#' \itemize{
#'   \item col 1 is site labels.
#'   \item if data is in long format, col 2 is date of observation.
#'   \item next J columns are the observations (y) - counts or 0/1's.
#'   \item next is a series of columns for the site variables (one column
#'     per variable).  The column header is the variable name.
#'   \item next is a series of columns for the observation-level variables.
#'   These are in sets of J columns for each variable, e.g., var1-1 var1-2
#'   var1-3 var2-1 var2-2 var2-3, etc.  The column header of the first
#'   variable in each group must indicate the variable name.
#' }
#' @title Convert .CSV File to an unMarkedFrame
#' @param filename string describing filename of file to read in
#' @param long \code{FALSE} if file is in long format or \code{TRUE} if
#'    file is in long format (see \emph{Details})
#' @param species if data is in long format with multiple species, then
#'    this can specify a particular species to extract if there is a
#'    column named "species".
#' @return an unMarkedFrame object
#' @keywords utilities
#' @examples
#' # examine a correctly formatted long .csv
#' read.csv(system.file("csv","frog2001pcru.csv", package="unmarked"))
#'
#' # examine a correctly formatted wide .csv
#' read.csv(system.file("csv","widewt.csv", package="unmarked"))
#'
#' # convert them!
#' dat1 <- csvToUMF(system.file("csv","frog2001pcru.csv", package="unmarked"), long = TRUE)
#' dat2 <- csvToUMF(system.file("csv","frog2001pfer.csv", package="unmarked"), long = TRUE)
#' dat3 <- csvToUMF(system.file("csv","widewt.csv", package="unmarked"), long = FALSE)
#' @author Ian Fiske \email{ianfiske@@gmail.com}
#' @export
csvToUMF <-
function(filename, long=FALSE, species = NULL)
{
  dfin <- read.csv(filename)

  if(long == TRUE) return(formatLong(dfin, species))
  else return(formatWide(dfin))
}

# utility function to create a variable that follows the dates as 1,2,3,...
# site id is first column
# julian date is second column
dateToObs <-
function(dfin)
{
  sitecol <- dfin[[1]]
  datecol <- dfin[[2]]


  # order by site, then obs date
  dfin <- dfin[order(sitecol,datecol),]
  sitecol <- dfin[[1]]
  datecol <- dfin[[2]]

  dTab <- table(datecol,sitecol)
  sites <- unique(sitecol)
  nSite <- length(sites)
  nStop <- colSums(dTab)
  nStop <- nStop[nStop > 0]  # get rid of the stops for sites with no stops

  obsNum <- numeric(length(sitecol))
  # for each site i, replace stops with 1:nStop[i]
  for(i in 1:nSite){
    stops <- which(sitecol == sites[i])
    obsNum[stops] <- 1:nStop[i]
  }

  dfout <- cbind(dfin,obsNum)
  dfout
}

# take long data set and return data list
# column names must be
# site names, first
# date, one column
# response, one column
# obs vars, one per column
#' @export
formatLong <-
function(dfin, species = NULL)
{

  ## copy dates to last column so that they are also a covdata var
  nc <- ncol(dfin)
  dfin[[nc+1]] <- dfin[[2]]
  names(dfin)[nc+1] <- "Date"

  if(!is.null(species)) {
    dfin$y <- ifelse(dfin$species == species, dfin$y, 0)
    dfin$y[is.na(dfin$y)] <- 0
    dfin$species = NULL
  }

  # sum up counts within time/site
  expr <- substitute(recast(dfin[,1:3], sv + dv ~ ..., id.var = 1:2,
                            fun.aggregate = sum),
                     list(sv = as.name(names(dfin)[1]),
                          dv = as.name(names(dfin)[2])))
  dfin2 <- eval(expr)
  dfin1 <- dfin[!duplicated(dfin[,1:2]),]

  dfin <- merge(dfin1,dfin2, by = 1:2)
  dfin[,3] <- dfin[,length(dfin)]
  dfin <- dfin[,-length(dfin)]
  names(dfin)[3] <- "y"

  dfin <- dateToObs(dfin)
  dfnm <- colnames(dfin)
  nV <- length(dfnm) - 1  # last variable is obsNum
  expr <- substitute(recast(dfin, newvar ~ obsNum + variable,
    id.var = c(dfnm[1],"obsNum"), measure.var = dfnm[3]),
                                list(newvar=as.name(dfnm[1])))
  y <- as.matrix(eval(expr)[,-1])
  attr(y,"class") <- "matrix"

  expr <- substitute(recast(dfin, newvar ~ obsNum ~ variable,
    id.var = c(dfnm[1],"obsNum"), measure.var = dfnm[4:nV]),
                                list(newvar=as.name(dfnm[1])))
  obsvars <- eval(expr)
  which.date <- which(dimnames(obsvars)$variable == "Date")
  dimnames(obsvars)$variable[which.date] <- "JulianDate"

  obsvars.matlist <- arrToList(obsvars)
  obsvars.veclist <- lapply(obsvars.matlist, function(x) as.vector(t(x)))
  obsvars.df <- data.frame(obsvars.veclist)

  unMarkedFrame(y, obsCovs = obsvars.df)
}

# column names must be
# site (optional, but if present, labeled "site")
# response: y.1, y.2, ..., y.J
# site vars: namefoo, namebar, ...
# obs vars: namefoo.1, namefoo.2, ..., namefoo.J, namebar.1, ..., namebar.J,...
#' @export
formatWide <-
function(dfin)
{
  # throw placeholder into sitedata
  sitedata <- data.frame(ones = rep(1,nrow(dfin)))

  obsdata <- list()

  if(identical(tolower(colnames(dfin))[1],"site")) dfin <- dfin[,-1]

  dfnm <- colnames(dfin)
  y <- grep("^y.",dfnm)
  J <- length(y)
  y <- as.matrix(dfin[,y])

  ncols <- length(dfnm)
  i <- J + 1
  while(i <= ncols) {     # loop through columns
    if(length(grep('.[[:digit:]]+$',dfnm[i]))) {  # check if this is obsdata
      nv <- sub('.[[:digit:]]+$','',dfnm[i])
      expr <- substitute(obsdata$newvar <- dfin[,i:(i+J-1)],
                list(newvar=as.name(nv)))
      eval(expr)
      i <- i + J
    }
    else {
      sitedata <- cbind(sitedata,dfin[,i])
      colnames(sitedata)[length(sitedata)] <- dfnm[i]
      i <- i + 1
    }
  }

  obsdata.veclist <- lapply(obsdata, function(x) as.vector(t(x)))
  obsdata.df <- data.frame(obsdata.veclist)
  sitedata$ones <- NULL

  unMarkedFrame(y, siteCovs = sitedata, obsCovs = obsdata.df)
}


# take a multiyear file and return correctly formated data
# formatted as above, but with year in first column
# col1 = year, c2 = site, c3 = juliandate or sample number
# c4 = y, c5 - cX = covariates
# add sample periods of NA to years with fewer samples
# to make balanced data... this eases future computations
#' @export
formatMult <-
function(df.in)
{
  years <- sort(unique(df.in[[1]]))
  nY <- length(years)
  df.obs <- list()
  nsamp <- numeric()
  maxsamp <- max(table(df.in[[1]], df.in[[2]])) # the maximum samples/yr
  for(t in 1:nY){
    df.t <- df.in[df.in[[1]] == years[t],] # subset for current year
    df.t <- df.t[,-1] # remove year column
    df.t <- dateToObs(df.t)
    nsamp <- max(df.t$obsNum)
    if(nsamp < maxsamp) {
      newrows <- df.t[1:(maxsamp - nsamp), ] # just a placeholder
      newrows[,"obsNum"] <- ((nsamp + 1) : maxsamp)
      newrows[,3 : (ncol(df.t) - 1)] <- NA
      df.t <- rbind(df.t, newrows)
    }
    df.obs <- rbind(df.obs,cbind(year = years[t],df.t))
  }
  dfnm <- colnames(df.obs)
  nV <- length(dfnm) - 1  # last variable is obsNum

  # create y matrix using reshape
  expr <- substitute(recast(df.obs, var1 ~ year + obsNum + variable,
    id.var = c(dfnm[2],"year","obsNum"), measure.var = dfnm[4]),
                                list(var1 = as.name(dfnm[2])))
  y <- as.matrix(eval(expr)[,-1])

  # create obsdata with reshape
  # include date (3rd col) and other measured vars
  expr <- substitute(recast(df.obs, newvar ~ year + obsNum ~ variable,
    id.var = c(dfnm[2],"year","obsNum"), measure.var = dfnm[c(3,5:nV)]),
                                list(newvar=as.name(dfnm[2])))
  obsvars <- eval(expr)

  rownames(y) <- dimnames(obsvars)[[1]]
  colnames(y) <- dimnames(obsvars)[[2]]
	umf <- unMarkedFrame(y, obsCovs = arrToList(obsvars), primaryNum = nY)
  return(umf)
}

# function to take data of form
# site  | species | count
# to
# site | spp1 | spp2 | ...
sppLongToWide <-
function(df.in)
{
  df.m <- melt(df.in, id = c("site", "spp"))
  df.out <- cast(df.m, site ~ spp, add.missing=T, fill = 0)
  df.out <- df.out[order(df.out$site),]
  df.out
}

# get estimated psi from rn fit
getPsi <-
function(lam)
{
  1-exp(-lam)
}

# get estimatd p from rn fit (only for a null type model so far)
getP.bar <-
function(lam, r)
{
  K = 30
  psi <- getPsi(lam)
  pN.k <- dpois(0:K,lam)
  pY.k <- 1 - (1 - r)^(0:30)
  sum(pY.k * pN.k)
}

handleNA <- function(stateformula, detformula, umf) {
  y <- umf@y
  # TODO: use J <- ncol(y) here and throughout instead of wrong use of obsNum?
  obsNum <- umf@obsNum
  M <- nrow(y)
  siteCovs <- umf@siteCovs
  obsCovs <- umf@obsCovs
  umf.clean <- umf

  ## set up obsCov indices
  sites <- rep(1:M, each = obsNum)
  obs <- rep(1:obsNum, M)

  ## assume that siteCovs have already been added to obsCovs
  X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
  V.mf <- model.frame(detformula, obsCovs, na.action = NULL)

  if(ncol(V.mf) > 0) {
    ## which sites have NA's in obsCovs included in detformula?
    V.NA <- apply(is.na(V.mf), 1, any)
    V.NA.obs <- cbind(sites[V.NA], obs[V.NA])
    V.NA.sites <- unique(sites[V.NA])
    umf.clean@y[V.NA.obs] <- NA
    if(any(!is.na(y[V.NA.obs]))) {
      warning(sprintf("NA(s) found in 'obsCovs' that were not in 'y' matrix.
Corresponding observation(s) 'y' were replaced with NA.
Observations removed from site(s) %s", paste(V.NA.sites,collapse=", ")))
    }
  }

  if(ncol(X.mf) > 0) {
    ## which sites have NA in site var included in stateformula?
    X.NA.sites <- unique(which(apply(is.na(X.mf), 1, any)))
    umf.clean@y[X.NA.sites,] <- NA
    if(length(X.NA.sites) > 0) {
      warning(sprintf("NA(s) found in 'siteCovs' that were not in 'y' matrix.
Corresponding site(s) in 'y' were replaced with NA: %s",
                      paste(X.NA.sites,collapse=", ")))
    }
  }

  ## which sites have all NA's in y?
  na.sites <- which(apply(is.na(umf.clean@y), 1, all))
  if(length(na.sites) > 0) {
    umf.clean@y <- umf.clean@y[-na.sites,]
    umf.clean@siteCovs <- subset(umf.clean@siteCovs,
                                 !seq(length=M) %in% na.sites)
    umf.clean@obsCovs <- umf.clean@obsCovs[!(sites %in% na.sites),]
  }

  return(umf.clean)
}



# function to move data around:
# converts array obsdata to a list
# copies site covariate info from obsdata to sitedata
# puts all site covariates back into obsdata
# needed because all site vars need to be available for both state and det models
arrangeData <-
function(data)
{
  require(abind)
  y <- data$y
  sitedata <- data$covdata.site
  obsdata <- data$covdata.obs

  J <- ncol(y)
  M <- nrow(y)
  nSV <- length(sitedata)

  # if not null, then just add "ones"
  if(!is.null(obsdata)) {
      if(class(obsdata) == "list") obsdata$ones <- matrix(1,M,J)
      if(class(obsdata) == "array") {
          obsdata <- abind(obsdata, ones = matrix(1,M,J))
      }
  }
  if(!is.null(sitedata)) sitedata <- cbind(ones = rep(1,M),sitedata)

  # if data components are null, create as just ones
  if(is.null(obsdata)) obsdata <- list(ones = matrix(1,M,J))
  if(is.null(sitedata)) sitedata=data.frame(ones = rep(1,M))

  # if obsdata is an array, coerce it to a list
  if(identical(class(obsdata),"array")) obsdata <- arrToList(obsdata)
  nOV <- length(obsdata)

  # move all site data (single vectors and matrices of J repeated vectors)
  # in obsdata into sitedata
  toDel <- numeric(0)
  nuniq <- function(x) length(as.numeric(na.omit(unique(x)))) # lil' helper fun
  for(i in 1:nOV){
    # test for equality across rows (matrix site data)
    eqRow <- as.numeric(apply(as.matrix(obsdata[[i]]), 1, nuniq) == 1)
    isRep <- as.logical(prod(eqRow)) # make sure all rows are
    # move into site data if (vector) or (repeated vector as matrix)
    if((dim(as.matrix(obsdata[[i]]))[2] == 1) || isRep){
      obsdata[[i]] <- matrix(obsdata[[i]],nrow = M, ncol = J)
      sitedata <- cbind(sitedata,obsdata[[i]][,1])
      colnames(sitedata)[length(sitedata)] <- names(obsdata[i])
      toDel <- c(toDel,i)
    }
    # ensure that obsdata is made of matrices rather than dataframes
    obsdata[[i]] <- as.matrix(obsdata[[i]])
  }
  if(length(toDel) > 0) {   #obsdata[[toDel]] <- NULL # remove sitedata from obsdata
    for(t in toDel) {
      obsdata[[t]] <- NULL
    }
  }
  if(length(obsdata) == 0) obsdata <- list(ones = matrix(1,M,J))
  nSV <- length(sitedata) # update nSV
  nOV <- length(obsdata)

  # make all site terms into obs terms by copying them to
  # obsdata (from vector to a matrix of repeated vectors)
  # needed if site variables are used to model detection.
  for(i in 1:nSV){
    obsdata[[nOV + i]] <- matrix(sitedata[[i]],nrow=M,ncol=J)
    names(obsdata)[nOV + i] <- colnames(sitedata)[i]
  }
  obsvars <- names(obsdata)
  nOV <- length(obsdata) # update length

  list(y = y, covdata.site = sitedata, covdata.obs = obsdata)
}

getDesign <- function(stateformula, detformula, umf) {

  M <- nrow(umf@y)

  ## Compute detection design matrix
  ## add site Covariates at observation-level
  if(!is.null(umf@obsCovs)) {
    V.mf <- model.frame(detformula, umf@obsCovs, na.action = NULL)
    V <- model.matrix(detformula, V.mf)
  } else {
    V <- matrix(1, M*umf@obsNum, 1)
    colnames(V) <- "(Intercept)"
  }

  ## Compute state design matrix
  if(!is.null(umf@siteCovs)) {
    X.mf <- model.frame(stateformula, umf@siteCovs, na.action = NULL)
    X <- model.matrix(stateformula, X.mf)
  } else {
    X <- matrix(1, M, 1)
    colnames(X) <- "(Intercept)"
  }
  return(list(X = X, V = V))
}


meanstate <- function(x) {
    K <- length(x) - 1
    sum(x*(0:K))
}

truncateToBinary <- function(y) {
  if(max(y, na.rm = TRUE) > 1) {
    y <- ifelse(y > 0, 1, 0)
    warning("Some observations were > 1.  These were truncated to 1.")
  }
  return(y)
}

getSS <- function(phi) {
	ev.length <- nrow(phi)
	ev <- tryCatch(eigen(t(phi))$vectors[,1],
			error = function(x) rep(NA, ev.length))
	ev/sum(ev)
}

#' @export
imputeMissing <- function(umf, whichCovs) {
# impute observation covariates
	if(!is.null(umf@obsCovs)) {
		obsCovs <- umf@obsCovs
		J <- umf@obsNum
		M <- nrow(obsCovs)/J
		obs <- obsCovs[,whichCovs]
		whichrows <- apply(obs, 1, function(x) any(!is.na(x)))
		if(sum(whichrows) == 0) return(obsCovs)
		whichels <- matrix(whichrows, M, J, byrow = TRUE)
		for(i in seq(length=length(whichCovs))) {
			obs.i <- obs[,i]
			obs.i.mat <- matrix(obs.i, M, J, byrow = TRUE) # get ith obsvar
			obs.i.missing <- is.na(obs.i.mat) & whichels
			obs.i.imputed <- obs.i.mat
			for(j in 1:M) {
				for(k in 1:J) {
					if(obs.i.missing[j,k])
						if(all(is.na(obs.i.mat[j,]))) {
							obs.i.imputed[j,k] <- mean(obs.i.mat[,k], na.rm = T)
						} else {
							obs.i.imputed[j,k] <- mean(c(mean(obs.i.mat[j,],na.rm = T),
											mean(obs.i.mat[,k], na.rm = T)))
						}
				}
			}
			obsCovs[,i] <- as.numeric(t(obs.i.imputed))
		}
		umf@obsCovs <- obsCovs
	}
	# TODO: impute site covariates
	return(umf)
}