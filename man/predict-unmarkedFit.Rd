\name{predict-unmarkedFit}
\title{predict-unmarkedFit}
\alias{predict-umDistsampFit,}
\alias{unmarkedFit-method}
\examples{data(mallard)
mallardUMF <- unmarkedFrame(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)
fm.mallard <- pcount(~ length + elev + forest, ~ ivel + date, mallardUMF)

# Predict from fitted model using original data
predict(fm.mallard, type="state", backTran=TRUE)
predict(fm.mallard, type="det", backTran=FALSE)

# Create newdata with only 'length' varying, then predict and plot
nd <- data.frame(length=rnorm(10), elev=1, forest=0, ivel=-2, date=3)

(Elam <- predict(fm.mallard, type="state", newdata=nd))

with(data.frame(Elam), {
plot(nd$length, Predicted, ylim=c(0, 25))
segments(nd$length, Predicted+SE, nd$length, Predicted-SE)
})}
