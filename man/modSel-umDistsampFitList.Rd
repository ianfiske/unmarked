\name{modSel-umDistsampFitList}
\title{Model selection results from an umDistsampFitList...}
\description{Model selection results from an umDistsampFitList}
\alias{modSel}
\alias{modSel-methods}
\arguments{\item{object}{an object of class "umDistsampFitList"}}
\examples{data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

## Half-normal detection function. Density output. No covariates. 
## lineDat$Length is transect lengths in km, so it has to be converted.
(fm1 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, survey="line", unitsIn="m"))

(fm2 <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, survey="line", unitsIn="m"))

(fm3 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~area, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, survey="line", unitsIn="m"))

fl <- dsFitList(fitList = list(Null=fm1, A.=fm2, .A=fm2))
fl

(ms <- modSel(fl, nullmod=fm1))

ms}
