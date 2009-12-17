\name{fitList}
\alias{fitList}
\title{constructor of unmarkedFitList objects}
\usage{fitList(..., fits)}
\description{Organize nested models for model selection or model-averaged prediction.}
\arguments{
\item{...}{Fitted nested models. Preferrably named.}
\item{fits}{An alternative way of providing the models. A (preferrably named) list of fitted models.}
}
\examples{# Fit some N-mixture models
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)

fm1 <- pcount(~ ivel ~ length, mallardUMF)
fm2 <- pcount(~ ivel ~ 1, mallardUMF)
fm3 <- pcount(~ 1 ~ 1, mallardUMF)

## Two methods of creating an unmarkedFitList using fitList()

# Method 1
fmList <- fitList(Global=fm1, ivel.=fm2, Null=fm3)

# Method 2. Note that the arugment name "fits" must be included in call.
models <- list(Global=fm1, ivel.=fm2, Null=fm3)
fmList <- fitList(fits = models)


# Model-averaged prediction
predict(fmList, type="state")

# Model selection
modSel(fmList, nullmod=fm3)
}
