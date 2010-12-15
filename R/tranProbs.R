


tranProbs <- function(Kr, omegaR, gammaR) {
    .Call("tranProbs", 
        as.integer(Kr),
        as.double(omegaR),
        as.double(gammaR),
        PACKAGE = "unmarked")
    }