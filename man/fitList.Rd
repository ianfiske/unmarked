\name{fitList}
\alias{fitList}
\title{constructor of unmarkedFitList objects...}
\usage{fitList(fits)}
\description{constructor of unmarkedFitList objects}
\arguments{
\item{fits}{A (preferrably named) list of fitted models}
}
\examples{# Fit some N-mixture models
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)

fm1 <- pcount(~ ivel ~ length, mallardUMF)
fm2 <- pcount(~ ivel ~ 1, mallardUMF)
fm3 <- pcount(~ 1 ~ 1, mallardUMF)

# Create an unmarkedFitList with a named list of models
fmList <- fitList(fits=list(Global=fm1, ivel.=fm2, Null=fm3))
fmList

# Model-averaged prediction
predict(fmList, type="state")

# Model selection
modSel(fmList, nullmod=fm3)}
