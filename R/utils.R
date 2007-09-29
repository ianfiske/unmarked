# use logarithms to vectorize row-wise products
# this speeds things up a LOT (vs. apply(x,1,prod))
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
formatLong <-
function(dfin, species)
{
  library(reshape)

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
  sitedata <- data.frame(one = rep(1,nrow(dfin)))
  
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
function(df.in)
{
  require(reshape)
  years <- sort(unique(df.in[[1]]))
  df.obs <- list()
  nsamp <- numeric()
  maxsamp <- max(table(df.in[[1]], df.in[[2]])) # the maximum samples/yr
  for(t in 1:length(years)){
    df.t <- df.in[df.in[[1]] == years[t],] # subset for current year
    df.t <- df.t[,-1] # remove year column
    df.t <- dateToObs(df.t)
    nsamp <- max(df.t$obsNum)
    if(nsamp < maxsamp) {
      newrows <- df.t[maxsamp - nsamp, ]
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
  list(ymat=y, covdata.obs = obsvars)
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
handleNA <-
function(data, stateformula, detformula)
{
  library(abind)
  y <- data$y
  sitedata <- data$covdata.site
  obsdata <- data$covdata.obs

  # get variables to handle, remove others
  state.vars <- attr(terms(stateformula),"term.labels")
  obs.vars <- attr(terms(detformula),"term.labels")
  
  obsdata.ar <- abind(obsdata, along = 3)

  # remove variables that are not used
  sitedata <- sitedata[,c("ones",state.vars)]   # CHECK THIS LINE
  obsdata.ar <- obsdata.ar[,,c("ones",obs.vars)]

  obsdata.NA <- is.na(obsdata.ar)
  obsdata.NA <- apply(obsdata.NA, c(1,2), any)

  # determine which sites have problem to give informative warnings.
  whichsites <- names(which(apply(!is.na(y) & obsdata.NA, 1, any)))
  whichsites <- paste(whichsites, collapse = ", ")

  if(any(!is.na(y) & obsdata.NA)) {
    warning("NA(s) found in 'covdata.obs' that were not in 'y' matrix; corresponding observations 'y' were replaced with NA")
    warning(sprintf("Sites were %s", whichsites))
  }
  is.na(y) <- (is.na(y) | obsdata.NA) # replace 'y' with NA if cov is NA
  if(any(apply(is.na(y), 1, all))) warning("site(s) found with NA's for all observations in 'y'; these sites cannot be analyzed and have been removed")

  sitedata.na <- is.na(sitedata)
  if(!is.null(dim(sitedata.na))) {
    sitedata.na <- apply(sitedata.na, 1, any)
  }
  if(any(sitedata.na)) warning("NA(s) found in 'covdata.site'. site(s) were removed from 'y' and covariate data")

  # remove sites that have either all NA or an NA in any site covariate
  # from y matrix and all covariate data
  to.rm <- which(sitedata.na | apply(is.na(y), 1, all))
  if(length(to.rm) > 0) {
    y <- y[ - to.rm, ]  # remove from y
    sitedata <- sitedata[ - to.rm, ]  # remove from sitedata
    for(i in 1:length(obsdata)){      # remove from obsdata
      obsdata[[i]] <- obsdata[[i]][ - to.rm, ]
    }
  }
  list(y = y, covdata.site = sitedata, covdata.obs = obsdata)
}

# function to move data around:
# converts array obsdata to a list
# copies site covariate info from obsdata to sitedata
# puts all site covariates back into obsdata
# needed because all site vars need to be available for both state and det models
arrangeData <-
function(data)
{

  y <- data$y
  sitedata <- data$covdata.site
  obsdata <- data$covdata.obs
  
  J <- ncol(y)
  M <- nrow(y)
  nSV <- length(sitedata)
  
  if(is.null(obsdata)) obsdata <- list(ones = matrix(1,M,J))

  # if obsdata is an array, coerce it to a list
  if(identical(class(obsdata),"array")) obsdata <- arrToList(obsdata)    
  nOV <- length(obsdata)

  # move all site data (single vectors and matrices of J repeated vectors)
  # in obsdata into sitedata
  if(is.null(sitedata)) sitedata=data.frame(ones = rep(1,M))
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
  }
  if(length(toDel) > 0) obsdata[[toDel]] <- NULL # remove sitedata from obsdata
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
