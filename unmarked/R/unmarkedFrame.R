
############# VALIDATION FUNCTIONS #############################################

validunmarkedFrame <- function(object) {
	errors <- character(0)
	M <- nrow(object@y)
	if(!is.null(object@siteCovs))
		if(nrow(object@siteCovs) != M)
			errors <- c(errors, 
				"siteCovData does not have same size number of sites as y.")
	if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
		if(nrow(object@obsCovs) != M*obsNum(object))
			errors <- c(errors, "obsCovData does not have M*obsNum rows.")
	if(!all(is.na(object@plotArea)) & any(object@plotArea < 0))
		errors <- c(errors, "plotArea cannot contain negative values.")	
	if(length(errors) == 0)
		TRUE
	else
		errors
}

############ DATA CLASSES ######################################################

# Class to hold data for analyses in unmarked.
setClass("unmarkedFrame",
    representation(y = "matrix",
        obsCovs = "optionalDataFrame",
        siteCovs = "optionalDataFrame",
		mapInfo = "optionalMapInfo",
		plotArea = "numeric",
        obsToY = "optionalMatrix"),
    validity = validunmarkedFrame)

## a class for multi-season data

setClass("unmarkedMultFrame",
		representation(numPrimary = "numeric",
				yearlySiteCovs = "optionalDataFrame"),  # a data frame in site-major, year-minor order describing site-level covariates
		contains="unmarkedFrame")

## a class for distance sampling data
setClass("unmarkedFrameDS", 
		representation(
				dist.breaks = "numeric",
				tlength = "numeric",
				survey = "character",
				unitsIn = "character"),
		contains = "unmarkedFrame",
		validity = function(object) 
		{
			J <- numY(object)
			db <- object@dist.breaks
			if(J == length(db) - 1 && db[1] == 0) TRUE
			else "ncol(y) must equal length(dist.breaks)-1 and 
						dist.breaks[1] must equal 0"
		})


setClass("unmarkedFrameOccu",
		contains = "unmarkedFrame")


setClass("unmarkedFramePCount",
		contains = "unmarkedFrame")


setClass("unmarkedFrameMPois",
		representation(
			samplingMethod = "character",
			piFun = "character"),
		contains = "unmarkedFrame")

################### CONSTRUCTORS ###############################################

# Constructor for unmarkedFrames.
unmarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo,
    plotArea, obsToY) {

  if(class(obsCovs) == "list") {
    obsVars <- names(obsCovs)
    for(i in seq(length(obsVars))) {
      if(!(class(obsCovs[[i]]) %in% c("matrix", "data.frame")))
        stop("At least one element of obsCovs is not a matrix or data frame.")
      if(ncol(obsCovs[[i]]) != ncol(y) | nrow(obsCovs[[i]]) != nrow(y))
        stop("At least one matrix in obsCovs has incorrect number of dimensions.")
    }
    if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
    obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
  }

  if(("data.frame" %in% class(y)) |
      ("cast_matrix" %in% class(y))) y <- as.matrix(y)

	if(missing(obsToY)) obsToY <- NULL
	if(missing(mapInfo)) mapInfo <- NULL
	if(missing(plotArea) || is.null(plotArea)) plotArea <- rep(1, nrow(y))
	
  umf <- new("unmarkedFrame", y = y, obsCovs = obsCovs,
      siteCovs = siteCovs, mapInfo = mapInfo, plotArea = plotArea, 
	  obsToY = obsToY)

  return(umf)
}


unmarkedFrameDS <- function(y, siteCovs = NULL, dist.breaks, tlength, survey,
		unitsIn, mapInfo = NULL, plotArea = NULL)
{
	if(is.null(plotArea))
		plotArea <- as.numeric(rep(NA, nrow(y)))
	else {
		if(is.matrix(plotArea))
			plotArea <- c(t(plotArea))
		else
			stop("plotArea must be NULL or an M x J matrix of plot areas in same units as unitsOut argument to be used in distsamp")
		}
	if(missing(tlength) & survey == "point")
		tlength <- numeric(0)
	umfds <- new("unmarkedFrameDS", y = y, obsCovs = NULL,
			siteCovs = siteCovs, dist.breaks = dist.breaks, tlength = tlength,
			survey = survey, unitsIn = unitsIn, plotArea = plotArea,
			obsToY = matrix(1, 1, ncol(y)))
	return(umfds)
}



unmarkedFrameOccu <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo) 
{
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J), 
		mapInfo = mapInfo)
	umf <- as(umf, "unmarkedFrameOccu")
	umf
}

# This function constructs an unmarkedMultFrame object.
unmarkedMultFrame <- function(y, siteCovs = NULL, obsCovs = NULL, numPrimary,
	yearlySiteCovs = NULL, plotArea = NULL) 
{
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J), 
		plotArea = plotArea)
	umf <- as(umf, "unmarkedMultFrame")
	umf@numPrimary <- numPrimary
	umf@yearlySiteCovs <- yearlySiteCovs
	umf
}



unmarkedFramePCount <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo,
 	plotArea = NULL) 
{
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J), 
		mapInfo = mapInfo, plotArea = plotArea)
	umf <- as(umf, "unmarkedFramePCount")
	umf
}



unmarkedFrameMPois <- function(y, siteCovs = NULL, obsCovs = NULL, type, obsToY, 
	mapInfo, piFun, plotArea = NULL) 
{
	if(!missing(type)) {
		switch(type,
			removal = {
				obsToY <- matrix(1, ncol(y), ncol(y))
				obsToY[col(obsToY) < row(obsToY)] <- 0
				piFun <- "removalPiFun"
			},
			double = {
				obsToY <- matrix(c(1, 0, 0, 1, 1, 1), 2, 3)
				piFun <- "doublePiFun"
			})
		}
	else {
		if(missing(obsToY)) 
			stop("obsToY is required for multinomial-Poisson data with no specified type.")
		type <- "userDefined"
		}
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY, 
		mapInfo = mapInfo, plotArea = plotArea)
	umf <- as(umf, "unmarkedFrameMPois")
	umf@piFun <- piFun
	umf@samplingMethod <- type
	umf
}

################ SHOW METHODS ##################################################


setMethod("show", "unmarkedFrame",
		function(object) {
			df <- as(object, "data.frame")
			cat("Data frame representation of unmarkedFrame object.\n")
			print(df)
		})

############################ EXTRACTORS ########################################

# Extractor for site level covariates
setGeneric("siteCovs", function(object,...) standardGeneric("siteCovs"))
setMethod("siteCovs", "unmarkedFrame",
		function(object) {
			return(object@siteCovs)
		})

setGeneric("yearlySiteCovs", 
	function(object,...) standardGeneric("yearlySiteCovs"))
setMethod("yearlySiteCovs", "unmarkedMultFrame",
		function(object) {
			return(object@yearlySiteCovs)
		})

setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))
setMethod("obsCovs", "unmarkedFrame", 
		function(object, matrices = FALSE) {
			M <- numSites(object)
			R <- obsNum(object)
			if(matrices) {
				value <- list()
				for(i in seq(length=length(object@obsCovs))){
					value[[i]] <- matrix(object@obsCovs[,i], M, R, byrow = TRUE)
				}
				names(value) <- names(object@obsCovs)
			} else {
				value <- object@obsCovs
			}
			return(value)
		})


setGeneric("obsNum", function(object) standardGeneric("obsNum"))
setMethod("obsNum", "unmarkedFrame", function(object) nrow(object@obsToY))


setGeneric("numSites", function(object) standardGeneric("numSites"))
setMethod("numSites", "unmarkedFrame", function(object) nrow(object@y))


setGeneric("numY", function(object) standardGeneric("numY"))
setMethod("numY", "unmarkedFrame", function(object) ncol(object@y))


setGeneric("obsToY", function(object) standardGeneric("obsToY"))
setMethod("obsToY", "unmarkedFrame", function(object) object@obsToY)


setGeneric("obsCovs<-", function(object, value) standardGeneric("obsCovs<-"))
setReplaceMethod("obsCovs", "unmarkedFrame", function(object, value) {
			object@obsCovs <- as.data.frame(value)
			object
		})


setGeneric("siteCovs<-", function(object, value) standardGeneric("siteCovs<-"))
setReplaceMethod("siteCovs", "unmarkedFrame", function(object, value) {
			object@siteCovs <- as.data.frame(value)
			object
		})


setGeneric("obsCovs<-", function(object, value) standardGeneric("obsCovs<-"))
setReplaceMethod("obsCovs", "unmarkedFrame", function(object, value) {
			object@obsCovs <- as.data.frame(value)
			object
		})


setGeneric("yearlySiteCovs<-", 
	function(object, value) standardGeneric("yearlySiteCovs<-"))
setReplaceMethod("yearlySiteCovs", "unmarkedMultFrame", function(object, value) {
			object@yearlySiteCovs <- as.data.frame(value)
			object
		})


setGeneric("obsToY<-", function(object, value) standardGeneric("obsToY<-"))
setReplaceMethod("obsToY", "unmarkedFrame", function(object, value) {
			object@obsToY <- value
			object
		})



setGeneric("getY", function(object) standardGeneric("getY"))
setMethod("getY", "unmarkedFrame", function(object) object@y)


setGeneric("coordinates", function(object) standardGeneric("coordinates"))
setMethod("coordinates", "unmarkedFrame",
		function(object) {
			object@mapInfo@coordinates
		})


setGeneric("projection", function(object) standardGeneric("projection"))
setMethod("projection", "unmarkedFrame",
		function(object) {
			object@mapInfo@projection
		})

################################### SUMMARY METHODS ############################


setMethod("summary", "unmarkedFrame",
	function(object,...) {
		cat("unmarkedFrame Object\n\n")
		cat(nrow(object@y), "sites\n")
		cat("Maximum number of observations per site:",obsNum(object),"\n")
			mean.obs <- mean(rowSums(!is.na(getY(object))))
		cat("Mean number of observations per site:",round(mean.obs,2),"\n")
		cat("Sites with at least one detection:", 
			sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))), 
				"\n\n")
		cat("Tabulation of y observations:")
		print(table(object@y, exclude=NULL))
		if(!is.null(object@siteCovs)) {
			cat("\nSite-level covariates:\n")
			print(summary(object@siteCovs))
			}
		if(!is.null(object@obsCovs)) {
			cat("\nObservation-level covariates:\n")
			print(summary(object@obsCovs))
			}
	})

	

setMethod("summary", "unmarkedFrameDS", 
	function(object, ...) 
{
	cat("unmarkedFrameDS Object\n\n")
	cat(object@survey, "-transect survey design", "\n", sep="")
	cat(paste("Distance class cutpoints (", object@unitsIn, "): ", sep=""), 
		object@dist.breaks, "\n")
	cat("plot area information supplied? :", !all(is.na(object@plotArea)), "\n\n")
	cat(nrow(object@y), "sites\n")
	cat("Maximum number of distance classes per site:", ncol(getY(object)), "\n")
		mean.dc <- mean(rowSums(!is.na(getY(object))))
	cat("Mean number of distance classes per site:", round(mean.dc, 2), "\n")
	cat("Sites with at least one detection:", 
		sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))), "\n\n")
	cat("Tabulation of y observations:")
	print(table(object@y, exclude=NULL))
	if(!is.null(object@siteCovs)) {
		cat("\nSite-level covariates:\n")
		print(summary(object@siteCovs))
		}
	if(!is.null(object@obsCovs)) {
		warning("Observation-level covariates cannot be used by distsamp()")
		}
})


################################# PLOT METHODS #################################
# TODO:  come up with nice show/summary/plot methods for each of these data types.

#setMethod("plot", c(x="unmarkedFrame", y="missing"),
#		function(x) {
#			M <- numSites(x)
#			J <- obsNum(x)
#			df <- data.frame(site = rep(1:M, each = J), 
#				obs = as.factor(rep(1:J, M)), y = as.vector(t(getY(x))))
#			df <- na.omit(df)
#			g <- ggplot(aes(x = site, y = obs, fill=y), data=df) + geom_tile() +
#				coord_flip()
#			g
#		})
#
#setMethod("plot", c(x="unmarkedFrameOccu", y="missing"),
#		function(x) {
#			M <- numSites(x)
#			J <- obsNum(x)
#			df <- data.frame(site = rep(1:M, each = J), 
#				obs = as.factor(rep(1:J, M)), 
#				y = as.ordered(as.vector(t(getY(x)))))
#			df <- na.omit(df)
#			g <- ggplot(aes(x = site, y = obs, fill=y), data=df) + geom_tile() +
#				coord_flip() + scale_fill_brewer(type="seq")
#			g
#		})
#
## this needs some work.
#setMethod("plot", c(x="unmarkedFramePCount", y="missing"),
#		function(x) {
#			## create lines for sites that have pos obs?
#			M <- numSites(x)
#			J <- obsNum(x)
#			df <- data.frame(site = rep(1:M, each = J), 
#				obs = as.factor(rep(1:J, M)), y = as.vector(t(getY(x))))
#			df <- na.omit(df)
#			g <- ggplot(aes(x = site, y = obs, colour=y,size=y), data=df) +
#				geom_point() + coord_flip() + theme_bw() + scale_fill_gradient()
#			g
#		})
#
#setMethod("plot", c(x="unmarkedMultFrame", y="missing"),
#		function(x) {
#			M <- numSites(x)
#			nY <- x@numPrimary
#			J <- obsNum(x)
#			obsNum <- J/nY
#			df <- data.frame(site = rep(1:M, each = J), 
#				obs = as.factor(rep(1:J, M)), 
#				y = as.ordered(as.vector(t(getY(x)))))
#			df <- na.omit(df)
#			yrs <- data.frame(x=rep(0,nY-1), xend=rep(M,nY-1), 
#				yrs=seq(obsNum, J-obsNum,by=obsNum)+0.5)
#			g <- ggplot(data=df) + geom_tile(aes(x = site, y = obs, fill=y)) +
#				coord_flip() + scale_fill_brewer(type="seq") +
#				geom_segment(aes(x=x,xend=xend,y=yrs,yend=yrs),data=yrs) +
#				theme_bw() + scale_y_discrete(breaks=seq(1,J-obsNum+1,by=obsNum))
#			g
#		})



setMethod("plot", c(x="unmarkedFrame", y="missing"),
	# plt <- 
	function (x, y, col=terrain.colors, zeroNAcolors=c(gray(0.5), gray(0)),
		xlab="Observation", ylab="Site", addgrid = FALSE, ...) 
{
	y <- getY(x)
	M <- nrow(y)
	J <- ncol(y)
	tab <- table(y, useNA = "ifany")
	vals <- as.numeric(names(tab))
	ymax <- max(y, na.rm=T)
	ncolors <- ymax + ifelse(any(vals==0), 1, 0)
	lt <- length(tab)
	op <- par(no.readonly = TRUE)
	laymat <- matrix(1, 5, 8)
	laymat[, 8] <- c(0, 2, 2, 2, 0)
	layout(laymat)
	if(is.function(col)) 
		color <- do.call(col, list(n=ncolors))
	else 
		color <- col
	color[which(vals==0)] <- zeroNAcolors[1]
	z <- t(y)
	z[!is.na(z)] <- 0
	z[is.na(z)] <- 1
	image(1:J, 1:M, z, col=zeroNAcolors, xaxt="n", xlab=xlab, ylab=ylab,
		...)
	image(1:J, 1:M, t(y), col=color, add=T)
	box()
	axis(1, at = 1:J, ...)
	if(addgrid)
		grid(J, M, lty = 1, lwd=0.1)
	par(mai = c(0.1, 0.2, 0.4, 0.4))
	zv <- vals
	zv[!is.na(zv)] <- 0
	zv[is.na(zv)] <- 1
	image(1:2, 1:lt, matrix(zv, 1), col=zeroNAcolors, xaxt = "n", 
		yaxt = "n", main = "Legend")
	image(1:2, 1:lt, matrix(vals, 1), col=color, add=T)
	box()
	if (any(is.null(names(tab)))) 
		laborder <- c(lt, 1:(lt - 1))
	else laborder <- 1:lt
	axis(2, at = 1:lt, labels = paste(names(tab))[laborder], las=1, ...)
	par(op)
})




setMethod("hist", "unmarkedFrameDS",
	function(x, ...)
{
	y <- getY(x)
	dbreaks <- x@dist.breaks
	nb <- length(dbreaks)
 	mids <- (dbreaks[-1] - dbreaks[-nb]) / 2 + dbreaks[-nb]
    distances <- unlist(mapply(rep, mids, each=colSums(y)))
	hist(distances, breaks=dbreaks, ...)
})



#setMethod("plot", c(x="unmarkedFrameOccu", y="missing"),
#		function(x) {
#			if(is.null(x@mapInfo)) stop("mapInfo is required to plot an unmarkedFrameOccu object.")
#			y <- getY(x)
#			## get sites w/ at least one pos
#			y <- as.factor(rowSums(y, na.rm = TRUE) > 0)
#			levels(y) <- c("non-detection", "detection")
#			siteCovs <- siteCovs(x)
#			coords <- coordinates(x)
#			if(is.null(x@mapInfo@projection)) {
#				proj <- list(x = coords[,1], y = coords[,2])	
#			} else {
#				proj <- mapproject(x = coords[,1], y = coords[,2], projection = x@mapInfo@projection,
#						parameters = x@mapInfo@parameters, orientation = x@mapInfo@orientation)
#			}
#			p <- qplot(x = proj$x, y = proj$y, colour = y, xlab = "longitude", ylab = "latitude")
#			if(!is.null(x@mapInfo@projection)) {
#				p + coord_map(project = x@mapInfo@projection)
#			} else {
#				p
#			}
#		})
#

#setMethod("plot", c(x="unmarkedFramePCount", y="missing"),
#		function(x) {
#			if(is.null(x@mapInfo)) stop("mapInfo is required to plot an unmarkedFramePCount object.")
#			y <- getY(x)
#			## plot maximums
#			y <- apply(y, 1, max, na.rm = TRUE)
#			siteCovs <- siteCovs(x)
#			coords <- coordinates(x)
#			if(is.null(x@mapInfo@projection)) {
#				proj <- list(x = coords[,1], y = coords[,2])	
#			} else {
#				proj <- mapproject(x = coords[,1], y = coords[,2], projection = x@mapInfo@projection,
#					parameters = x@mapInfo@parameters, orientation = x@mapInfo@orientation)
#			}
#			p <- qplot(x = proj$x, y = proj$y, colour = y, xlab = "longitude", ylab = "latitude")
#			if(!is.null(x@mapInfo@projection)) {
#				p + coord_map(project = x@mapInfo@projection)
#			} else {
#				p
#			}
#		})
################################# SELECTORS ####################################

# i is the vector of sites to extract

setMethod("[", c("unmarkedFrame", "numeric", "missing", "missing"),
	function(x, i) {  
		if(length(i) == 0) return(x)
		if(any(i < 0) && any(i > 0)) 
			stop("i must be all positive or all negative indices.")
		if(all(i < 0)) { # if i is negative, then convert to positive
			M <- numSites(x)
			i <- (1:M)[i]
			}
		y <- getY(x)[i,]
		if (length(i) == 1) {
			y <- t(y)
			}
		siteCovNames <- colnames(siteCovs(x))
		siteCovs <- as.data.frame(siteCovs(x)[i,])
		colnames(siteCovs) <- siteCovNames
		obsCovs <- obsCovs(x)
		obsCovsNames <- colnames(obsCovs(x))
		R <- obsNum(x)
		obs.site.inds <- rep(1:numSites(x), each = R)
		obs.site.sel <- rep(i, each = R)
		obsCovs <- as.data.frame(obsCovs[match(obs.site.sel,obs.site.inds),])
		colnames(obsCovs) <- obsCovsNames
		umf <- x
		umf@y <- y
		umf@siteCovs <- siteCovs
		umf@obsCovs <- obsCovs
		umf@plotArea <- umf@plotArea[i]
		umf
	})

## remove obs only
setMethod("[", c("unmarkedFrame", "missing", "numeric", "missing"),
		function(x, i, j) {  
			y <- getY(x)
			obsCovs <- obsCovs(x)
			obsToY <- obsToY(x)
			obs.remove <- rep(TRUE, obsNum(x))
			obs.remove[j] <- FALSE
			y.remove <- t(obs.remove) %*% obsToY > 0
			y <- y[,!y.remove, drop=FALSE]
			obsCovs <- obsCovs[!rep(obs.remove, numSites(x)),, drop=FALSE]
			x@obsCovs <- obsCovs
			x@y <- y
			x@obsToY <- obsToY[!obs.remove,!y.remove, drop=FALSE]
			x
			
		})

# i is as before and j is the obsNum to remove and corresponding y's
setMethod("[", c("unmarkedFrame","numeric", "numeric", "missing"),
		function(x, i, j) {  
			## first remove sites
			umf <- x[i,]
			umf <- umf[,j]
			umf
		})

## for multframes, must remove years at a time
setMethod("[", c("unmarkedMultFrame", "missing", "numeric", "missing"),
		function(x, i, j) {  
			J <- obsNum(x)/x@numPrimary
			obs <- rep(1:x@numPrimary, each = J)
			numPrimary <- length(j)
			j <- which(!is.na(match(obs,j)))
			u <- callNextMethod(x, i, j)
			u@numPrimary <- numPrimary
			u
		})


setMethod("head", "unmarkedFrame",
		function(x, n) {
			if(missing(n)) n <- 10
			umf <- x[1:n,]
			umf
		})
		
############################### COERCION #######################################

setAs("data.frame", "unmarkedFrame", function(from) {
			umf <- formatWide(from)
			umf
		})

setAs("unmarkedFrame", "data.frame", function(from) {
			obsCovs <- obsCovs(from)
			siteCovs <- siteCovs(from)
			y <- getY(from)
			colnames(y) <- paste("y",1:ncol(y),sep=".")
			if(is.null(obsToY(from))) {
				obsNum <- ncol(y)
			} else {
				obsNum <- obsNum(from)
			}
			if(is.null(siteCovs)) siteCovs <- matrix(0,nrow(y),0)
			if(is.null(obsCovs)) {
				obsCovs <- matrix(0,nrow(y),0)
			} else {
				obsCovs <- data.frame(lapply(obsCovs, 
					function(x) matrix(x, nrow(y), obsNum,byrow=T)))
			}
			df <- data.frame(y, siteCovs, obsCovs)
			df
		})
		
		
		


######################### Power analysis methods ##############################



setGeneric("powerAnalysis", function(formula, data, coefs, ...) 
	standardGeneric("powerAnalysis"))

setMethod("powerAnalysis", c("formula", "unmarkedFramePCount", "numeric"),
	function(formula, data, coefs, nsim=1, sig.level=0.05, mixture="P", 
		report=FALSE, ...) 
	{
	    mixture <- match.arg(mixture)
    	if (!is(data, "unmarkedFramePCount")) 
        	stop("Data is not an unmarkedFramePCount object.")
    	designMats <- unmarked:::getDesign2(formula, data, na.rm=FALSE)
    	X <- designMats$X
    	V <- designMats$V
    	y <- designMats$y
#		y[] <- NA	# Safety
		plotArea <- designMats$plotArea
    	J <- ncol(y)
    	M <- nrow(y)
   	    nDP <- ncol(V)
	    nAP <- ncol(X)
	    nP <- nAP + nDP + ifelse(identical(mixture, "NB"), 1, 0)
        lamParms <- coefs[(nDP+1):(nDP+nAP)]
    	detParms <- coefs[1:nDP]
		lambda <- exp(X %*% lamParms) * plotArea
		p <- plogis(V %*% detParms)
		umf <- data
		fits <- list()
		for(i in 1:nsim) {
			if(report)
				cat("sim", i, fill=T)
			switch(mixture, 
				P = {
					N <- rpois(M, lambda = lambda)
					N <- rep(N, each=J)
					},
				NB = {
					N <- rnbinom(M, mu = lambda, coefs[nP])
					N <- rep(N, each=J)
					}
				)
			yvec <- rbinom(M*J, N, p)
			umf@y <- matrix(yvec, M, J, byrow=TRUE)
			fits[[i]] <- pcount(formula, umf, mixture=mixture, ...)
			}
		cis <- lapply(fits, function(x) 
			rbind(confint(x, type = "state", level = 1-sig.level), 
				confint(x, type = "det", level = 1-sig.level)))
	 	tests <- sapply(cis, function(x) x[,1] <= 0 & 0 <= x[,2])
	 	test <- 1 - rowSums(tests) / ncol(tests)
	 	names(test) <- names(coef(fits[[1]], altNames=TRUE))
	 	return(list(call = match.call(call = sys.call(-1)), sites = M, 
		 	occasions = J, sig.level = sig.level, coefs = coefs, nsim = nsim, 
			power=test))
	})
			
		