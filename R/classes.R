#' @import methods
{}

#' @exportClass optionalDataFrame
setClassUnion("optionalDataFrame", c("data.frame","NULL"))

#' @exportClass optionalMatrix
setClassUnion("optionalMatrix", c("matrix","NULL"))

#' @exportClass optionalNumeric
setClassUnion("optionalNumeric", c("numeric","NULL"))

#' @exportClass optionalCharacter
setClassUnion("optionalCharacter", c("character","NULL"))

#' Compute standard error of an object.
#'
#' @aliases SE-methods SE
#' @param obj object whose standard error is computed.
#' @exportMethod SE
setGeneric("SE",
    def = function(obj) {
      standardGeneric("SE")
    })

#' @importFrom stats confint fitted coef vcov predict update profile simulate residuals
#' @importFrom graphics plot barplot

#' @exportMethod confint
setGeneric("plot")
setGeneric("predict")
setGeneric("vcov")
setGeneric("coef")
setGeneric("summary")
setGeneric("update")
setGeneric("confint")
setGeneric("profile")
setGeneric("head")
setGeneric("fitted")
setGeneric("simulate")
setGeneric("residuals")
setGeneric("barplot")


#' Compute linear combinations of parameters.
#'
#' @aliases linearComb-methods linearComb
#' @seealso \code{\link{linearComb,unmarkedFit-method}}, \code{link{linearComb,unmarkedEstimate-method}}
#' @exportMethod linearComb
setGeneric("linearComb",
    function(obj, coefficients, ...) {
      standardGeneric("linearComb")
    })

#' Transform an object to it's natural scale.
#'
#' @aliases backTransform-methods
#' @param obj object to be transformed
#' @exportMethod backTransform
setGeneric("backTransform",
    function(obj, ...) {
      standardGeneric("backTransform")
    })

#' @exportMethod hessian
setGeneric("hessian",	function(object) standardGeneric("hessian"))

# TODO: make parent class unmarkedEstimate
# TODO: make binomial detection child class
# TODO: make binomial occ child class
# TODO: make poisson abundance child class
# TODO: make show method for each of these classes.
# TODO: make unmarkedFit show that calls the respective children.
# TODO: separate unmarkedFit class for each model type?... would contain different estimate types