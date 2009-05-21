\name{umDistsampFitList}
\alias{umDistsampFitList}
\title{Model-averaged prediction from a list of fitted models.}
\description{Model-averaged prediction from a list of fitted models.}
\arguments{\item{object}{a list of fitted distance sampling models all of class "umDistsampFit"}
\item{newdata}{an optional data.frame with new predictor variables}
\item{notconstant}{the column name of a variable in newdata that is not constant}
\item{type}{the type of prediction to make. Either "state" or "det".}
\item{link}{the type of link function to use}}
\examples{data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

#Fit a few models
(fmhnA.AH <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m"))
(fmhnA <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~1, linetran, 
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m"))
(fmhn.A <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~area, linetran, 
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m"))
(fmhn <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, 
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m"))
(fmexp <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, key="exp",
dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
unitsIn="m", starts=c(0,0)))

# Put models in a list and make it a "umDistsampFitList"
fl <- list(fmhn, fmexp, fmhnA.AH, fmhnA, fmhn.A)
fmList <- new("umDistsampFitList", fitlist=fl)

# Make a new data.frame and average predictions from the list
# Note, you shouldn't do type="det" in this case because multiple keyfuns used
newarea <- seq(min(linetran$area), max(linetran$area), length=20)
newhabitat <- factor("a")
nd <- data.frame(area=newarea, habitat=newhabitat)

(ElamA <- predict(fmList, newdata=nd, notconstant="area", type="state"))

with(ElamA, {
plot(Predictor, Est., ylim=c(0.7, 1.1), xlab="Area", 
ylab="Density (object / ha)")
segments(Predictor, Est.-SE, Predictor, Est.+SE)
})}
