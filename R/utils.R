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

# utility function to read a csv file and create data a required for other funs
# csv is formatted as follows:
# col 1 is site labels, must be named "site"

# if data is in long format, col 2 is date of observation
# then col for observations (y) and then all covariates
# long format may contain multiple species, where the column containing species names must
# be named "species"

# wide format:
# next J columns are the observations (y)
# next is a series of columns for the site variables (one column per variable)
#   the column header is the variable name
# next is a series of columns for the observation-level variables
#   these are in sets of J columns for each variable
#   e.g., var1-1 var1-2 var1-3 var2-1 var2-2 var2-3, etc.
#   the column header of the first variable in each group must indicate the
#   variable name.
csvToData <-
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
  require(reshape)

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
formatLong <-
function(dfin, species = NULL)
{
  require(reshape)

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
  expr <- substitute(recast(dfin, newvar ~ obsNum ~ variable,
    id.var = c(dfnm[1],"obsNum"), measure.var = dfnm[4:nV]),
                                list(newvar=as.name(dfnm[1])))
  obsvars <- eval(expr)

  which.date <- which(dimnames(obsvars)$variable == "Date")
  dimnames(obsvars)$variable[which.date] <- "JulianDate"

  list(ymat=y, covdata.obs = obsvars)
}

# column names must be
# site (optional, but if present, labeled "site")
# response: y.1, y.2, ..., y.J
# site vars: namefoo, namebar, ...
# obs vars: namefoo.1, namefoo.2, ..., namefoo.J, namebar.1, ..., namebar.J,...
formatWide <-
function(dfin)
{
  # throw placeholder into sitedata
  sitedata <- data.frame(ones = rep(1,nrow(dfin)))

  obsdata <- list()

  if(identical(colnames(dfin)[1],"site")) dfin <- dfin[,-1]

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

  list(ymat = y, covdata.site = sitedata, covdata.obs = obsdata)
}


# take a multiyear file and return correctly formated data
# formatted as above, but with year in first column
# col1 = year, c2 = site, c3 = juliandate or sample number
# c4 = y, c5 - cX = covariates
# add sample periods of NA to years with fewer samples
# to make balanced data... this eases future computations
formatMult <-
function(df.in, spp, state)
{
  require(reshape)
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
  list(ymat=y, covdata.obs = obsvars, J = ncol(y)/nY,
       species = spp, state = state)
}

# function to take data of form
# site  | species | count
# to
# site | spp1 | spp2 | ...
sppLongToWide <-
function(df.in)
{
  require(reshape)
  df.m <- melt(df.in, id = c("site", "spp"))
  df.out <- cast(df.m, site ~ spp, add.missing=T, fill = 0)
  df.out <- df.out[order(df.out$site),]
  df.out
}



# helper function to take siteformula, sitedata, detformula, obsdata
# and return design matrices, parameter names, numbers of parameters
getDesign <-
function(stateformula, detformula, y, sitedata, obsdata)
{
  if(is.null(dim(sitedata))) sitedata <- as.data.frame(sitedata)
  nSV <- length(sitedata)
  nOV <- length(obsdata)

  M = nrow(y)
  J = ncol(y)

  # get design matrix for occupancy variables
  mfOcc <- model.frame(formula = stateformula, data = sitedata)
  mtOcc <- attr(mfOcc, "terms")
  XOcc <- model.matrix(mtOcc, mfOcc)
  occParms <- c("psiconst",colnames(XOcc)[-1])
  nOP <- length(occParms)
  nOV <- length(obsdata)

  obsdata.ji.list <- list()
  for(i in 1:nOV) {
     obsdata.ji.list[[i]] <- as.numeric(obsdata[[i]])
  }
  names(obsdata.ji.list) <- names(obsdata)
  obsdata.ji <- do.call(cbind,obsdata.ji.list)
  detdf <- as.data.frame(obsdata.ji)
  mf <- model.frame(formula = detformula, data = detdf, na.action = na.pass)
  mt <- attr(mf,"terms")
  detParms <- c("pconst",attr(mt,"term.labels"))
  nDP <- length(detParms)
  XDet <- model.matrix(mt, mf)

  list(nOP = nOP, nDP = nDP, XDet = XDet, XOcc = XOcc,
    occParms = occParms, detParms = detParms)
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

# function to get rid of observations that have NA for any covariate
# and sites with NA's for all observations.
# Only throw out observations for which variables in varlist are missing.
# still need to implement for sitedata... unless sitedata is now vestigule
handleNA.old <-
function(data, stateformula, detformula)
{
  library(abind)
  y <- data$y
  sitedata <- data$covdata.site
  obsdata <- data$covdata.obs

  # ensure that data are names for giving informative warnings
  if(is.null(rownames(y))) rownames(y) <- 1:nrow(y)
##   if(!is.null(sitedata)) {
##     if(is.null(rownames(sitedata))) rownames(y) <- 1:nrow(sitedata)
##   }
##   if(!is.null(obsdata)) {
##     if(is.null(rownames(obsdata))) rownames(y) <- 1:nrow(obsdata)
##   }

  # get variables to handle, remove others
  # first, get names from formulas
  state.vars <- attr(terms(stateformula),"term.labels")
  obs.vars <- attr(terms(detformula),"term.labels")

  # ensure that "ones" is not called "one" (maybe deprecated)
  if(!is.null(sitedata$one)) {
    sitedata$ones <- sitedata$one
    sitedata$one <- NULL
  }
  if(!is.null(obsdata$one)) {
    names(obsdata)[names(obsdata) == "one"] <- "ones"
  }

  # coerce obsdata to an array for easier manipulation
  obsdata.ar <- abind(obsdata, along = 3)

  # if obsvar begins with I(), then remove obsvar before the following check
  # REALLY, I'd like to just parse out the variable name from the expression
  # and check it, but i'll do that later.
  if(length(grep("^I(",obs.vars,extended=FALSE)) > 0) {
    obs.vars2 <- obs.vars[-grep("^I(",obs.vars,extended=FALSE)]
  } else {
    obs.vars2 <- obs.vars
  }
  
  # check for formula specification that involves covariates not in obsdata
  # note that obs.vars2 is now here to ignore "I()" expressions
  if(any(!(obs.vars2 %in% dimnames(obsdata.ar)[[3]]))) {
    badvars <- !(obs.vars %in% dimnames(obsdata.ar)[[3]])
    badvars <- obs.vars[badvars]
    badvars <- paste(badvars, collapse = ", ")
    stop(sprintf("Detection covariate(s) %s were specified and were not present in the data.", badvars))
  }

  # remove variables that are not used so that NAs in unused variables
  # do not matter
  sitedata <- sitedata[,c("ones",state.vars)]
  obsdata.ar <- obsdata.ar[,,c("ones",obs.vars2)]  # NOTE obs.vars2 here too

  obsdata.NA <- is.na(obsdata.ar)
  obsdata.NA <- apply(obsdata.NA, c(1,2), any)

  # determine which sites have problem to give informative warnings.
  # NA in cov is problem only if it does not correspond to NA in y.
  whichsites <- names(which(apply(!is.na(y) & obsdata.NA, 1, any)))
  whichsites <- paste(whichsites, collapse = ", ")
  if(any(!is.na(y) & obsdata.NA)) {
    warning(sprintf("NA(s) found in 'covdata.obs' that were not in 'y' matrix.
Corresponding observation(s) 'y' were replaced with NA.
Observations removed from site(s) %s",whichsites))
  }
  is.na(y) <- (is.na(y) | obsdata.NA) # replace 'y' with NA if cov is NA

  # look for all NA's for a site
  whichsites <- names(which(apply(is.na(y),1,all)))
  whichsites <- paste(whichsites, collapse = ", ")
  if(any(apply(is.na(y), 1, all))) {
    warning(sprintf("Site(s) found with NA's for all observations in 'y'.
Site(s) %s cannot be analyzed and have been removed.",whichsites))
  }

  sitedata.na <- is.na(sitedata)
  if(!is.null(dim(sitedata.na))) {
    sitedata.na <- apply(sitedata.na, 1, any)
  }
  if(any(sitedata.na)) warning("NA(s) found in 'covdata.site'. site(s) were removed from 'y' and covariate data")

  # remove sites that have either all NA or an NA in any site covariate
  # from y matrix and all covariate data
  to.rm <- which(sitedata.na | apply(is.na(y), 1, all))
  sitedata <- as.matrix(sitedata)
  if(length(to.rm) > 0) {
    y <- y[ - to.rm, ]  # remove from y
    sitedata <- sitedata[ - to.rm, ]  # remove from sitedata
    for(i in 1:length(obsdata)) {      # remove from obsdata
      obsdata[[i]] <- obsdata[[i]][ - to.rm, ]
    }
  }
  sitedata <- as.data.frame(sitedata)
  list(y = y, covdata.site = sitedata, covdata.obs = obsdata)
}

handleNA <- function(stateformula, detformula, umf) {
  y <- umf@y
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
  
  ## which sites have NA in site var included in stateformula?
  X.NA.sites <- unique(which(apply(is.na(X.mf), 1, any)))
  umf.clean@y[X.NA.sites,] <- NA
  if(length(X.NA.sites) > 0) {
    warning(sprintf("NA(s) found in 'siteCovs' that were not in 'y' matrix.
Corresponding site(s) in 'y' were replaced with NA: %s",
                    paste(X.NA.sites,collapse=", ")))
  }

  ## which sites have all NA's in y?
  na.sites <- which(apply(is.na(umf.clean@y), 1, all))
  umf.clean@y <- umf.clean@y[-na.sites,]
  umf.clean@siteCovs <- umf.clean@siteCovs[-na.sites,]
  umf.clean@obsCovs <- umf.clean@obsCovs[!(sites %in% na.sites),]
  
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

meanstate <- function(x) {
    K <- length(x) - 1
    sum(x*(0:K))
}
