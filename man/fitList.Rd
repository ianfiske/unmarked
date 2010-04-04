\name{fitList}
\alias{fitList}
\title{constructor of unmarkedFitList objects}
\usage{fitList(..., fits)}
\description{Organize nested models for model selection or model-averaged prediction.}
\arguments{
\item{...}{Fitted nested models. Preferrably named.}
\item{fits}{An alternative way of providing the models. A (preferrably named) list of fitted models.}
}
\author{Richard Chandler \email{rchandler@nrc.umass.edu}}
\examples{# Fit some N-mixture models
data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 
lengths <- linetran$Length * 1000

ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths, survey = "line", unitsIn = "m")
	})

fm1 <- distsamp(~ 1 ~1, ltUMF)
fm2 <- distsamp(~ area ~1, ltUMF)
fm3 <- distsamp( ~ 1 ~area, ltUMF)

## Two methods of creating an unmarkedFitList using fitList()

# Method 1
fmList <- fitList(Null=fm1, .area=fm2, area.=fm3)

# Method 2. Note that the arugment name "fits" must be included in call.
models <- list(Null=fm1, .area=fm2, area.=fm3)
fmList <- fitList(fits = models)


# Model-averaged prediction
predict(fmList, type="state")

# Model selection
modSel(fmList, nullmod="Null")
}
