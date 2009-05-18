\name{umDistsampFit}
\alias{umDistsampFit}
\title{Model-based prediction...}
\description{Model-based prediction}
\details{Predict expectations of a model with standard errors (and CI if alpha is specified).

If newdata is not specified, original data is used.
Requires the deltamethod() function from package msm.}
\arguments{\item{object}{a fitted distance sampling model of class 'umDistsampFit'}
\item{type}{the component of the model to predict. Either "state" or "det"}
\item{link}{the link function to be used}
\item{newdata}{an optional data.frame including variables necessary to predict from fitted model}
\item{notconstant}{an optional column name of a variable in newdata which is not constant. Only useful for plotting expected relationships.}
\item{alpha}{alpha level for creating confidence intervals}}
\examples{data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

#Fit a model
(fmhnA.H <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m"))

# Create new data.frame for prediction
newhabitat <- factor("A")
newarea <- seq(min(linetran$area), max(linetran$area), length=20)
(ndA <- data.frame(area=newarea, habitat=newhabitat))

# Predict density based upon fitted model and 'new' data
(Elam.A <- predict(fmhnA.H, newdata=ndA, notconstant="area", type="state", 
link="log"))
with(Elam.A, { # Plot relationships between density and area
plot(Predictor, Est., ylim=c(0.5, 1.2), xlab="Area", ylab="Density")
segments(Predictor, Est.-SE, Predictor, Est.+SE)
})

# Same as above but for detection 
ndH <- data.frame(area=mean(linetran$area), habitat=factor(c("a", "b")))
ndH

(Ep.H <- predict(fmhnA.H, newdata=ndH, notconstant="habitat", type="det", 
link="log"))
with(Ep.H, {  # Plot relationship (lack of) difference in detectability between habitat types
bp <- barplot(Est., ylim=c(0, 15), names=Predictor, 
ylab="Sigma (half-normal shape parameter)") 
arrows(bp, Est., bp, Est.+SE, code=2, angle=90, length=0.2)
box()
})
# OR
with(Ep.H, {
plot(function(x) unmarked:::gxhn(x, Est.[1]), 0, 25, 
xlab="Distance (m)", ylab="Detection probability")
plot(function(x) unmarked:::gxhn(x, Est.[2]), 0, 25, add=T, lty=2, lwd=2)
legend(1, 20, c("Habitat a", "Habitat b"), lty=1:2)
})}
