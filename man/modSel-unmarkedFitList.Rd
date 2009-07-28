\name{modSel-unmarkedFitList}
\title{Model selection results from an unmarkedFitList...}
\description{Model selection results from an unmarkedFitList}
\alias{modSel}
\alias{modSel-methods}
\arguments{\item{object}{an object of class "umDistsampFitList"}}
\examples{data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 
lengths <- linetran$Length * 1000

ltUMF <- with(linetran, {
unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
tlength = lengths, survey = "line", unitsIn = "m")
})

## Half-normal detection function. Density output. No covariates. 
## lineDat$Length is transect lengths in km, so it has to be converted.
(fm1 <- distsamp(~ 1, ~1, ltUMF))

(fm2 <- distsamp(~ area, ~1, ltUMF))

(fm3 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~area, ltUMF)

fl <- unmarkedFitList(fitList = list(Null=fm1, A.=fm2, .A=fm2))
fl

(ms <- modSel(fl, nullmod=fm1))

ms@Full}
