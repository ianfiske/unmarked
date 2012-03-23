\name{crossbill}
\alias{crossbill.data}
\alias{crossbill}
\docType{data}
\title{
European crossbill data
}
\description{
Detection/non-detection data for the European crossbill in Switzerland,
and spatially-referenced covariates. The object crossbill.data is a list
with 2 components. The first is a data.frame with
detection/non-detection data and covariates at each of 267
square-kilometer sample locations. The second component is a list with
two matrices containing data on elevation and percent cover of forest at
each 1km-sq pixel in Switzerland. These can be easily converted to
rasters as shown in the example.
}
\usage{data(crossbill)}
\references{
Schmid, H., N. Zbinden, and V. Keller, 2004. Uberwachung der
bestandsentwicklung hauger brutvogel in der schweiz. Swiss
Ornithological Institute, Sempach Switzerland.
}
\examples{

\dontrun{
data(crossbill)

# Format for colext analysis
crossbill <- crossbill.data$crossbill
umf <- unmarkedMultFrame(y=as.matrix(crossbill[,5:31]),
    siteCovs=crossbill[,2:3], numPrimary=9)
(fm <- colext(~1, ~1, ~1, ~1, umf))

# Create rasters
library(raster)
elevation <- raster(crossbill.data$switzerland[[1]])
layerNames(elevation) <- "ele"
forest <- raster(crossbill.data$switzerland[[2]])
layerNames(forest) <- "forest"
rasters <- stack(elevation, forest)
plot(rasters)
}
}
\keyword{datasets}
