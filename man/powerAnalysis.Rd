\name{powerAnalysis}
\alias{powerAnalysis}

\title{Conduct a power analysis on an unmarked model}

\description{
This function uses a simulation-based approach to estimate power for parameters 
in unmarked models. At a minimum, users must provide a fitted \code{unmarked} model object
(preferably fit with simulated data) which ensures the model has been properly
specified, a list of effect sizes for each parameter in the model (\code{coefs}), 
and the desired Type I error (\code{alpha}). It is also possible to get power
for a range of other sample sizes besides the sample size in the fitted
model object using the \code{design} argument to subsample within the
provided dataset. See the \code{unmarkedPower} vignette for more details and
examples.
}

\usage{
  powerAnalysis(object, coefs=NULL, design=NULL, alpha=0.05, nulls=list(),
                datalist=NULL, 
                nsim=ifelse(is.null(datalist), 100, length(datalist)), 
                parallel=FALSE)
}

\arguments{
  \item{object}{A fitted model inheriting class \code{unmarkedFit}. This
    could potentially be fit using real data, but ideally you would simulate
    an appropriate dataset using \code{simulate}}
  \item{coefs}{A list containing the desired effect sizes for which you want
    to estimate power. This list must follow a specific format. There is one
    named entry in the list per submodel (e.g., occupancy, detection). To
    get the required submodel names call \code{names(object)} on your fitted model.
    Then, each list entry is a named vector with the names corresponding to the
    parameter names for that submodel, and the values corresponding to the 
    desired effect sizes. It may be easier to leave \code{coefs=NULL}, which
    will generate an error message with a template that you can fill in.
  }
  \item{design}{An optional list of design/sample size parameters containing 
    at a minimum two named elements: \code{M}, the number of sites, and \code{J} 
    the number of observations per site. If this list is provided, \code{unmarked} 
    will subsample the provided dataset to the specified number of sites and
    observations, allowing you to test power for different designs. If
    your model has multiple primary periods you must also include \code{T},
    the number of periods, in the list.
  }
  \item{alpha}{Desired Type I error rate}
  \item{nulls}{If provided, a list matching the structure of \code{coefs} which
    defines the null hypothesis value for each parameter. By default the null
    is 0 for all parameters.
  }
  \item{datalist}{An optional list of previously-simulated datasets, in the form
    of \code{unmarkedFrames} matching the model type of \code{object}, which
    will be used for the power analysis simulations.
  }
  \item{nsim}{Number of simulations to conduct}
  \item{parallel}{If \code{TRUE}, run folds in parallel. This may speed up 
    the power analysis in some situations
  }
} 

\value{\code{unmarkedPower} object containing the results of the power analysis}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{unmarkedPowerList}}
}

\examples{

\dontrun{

# Simulate an occupancy dataset
# Covariates to include in simulation
forms <- list(state=~elev, det=~1)

# Covariate effects and intercept values
coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

# Study design
design <- list(M=300, J=8) # 300 sites, 8 occasions per site

# Simulate an unmarkedFrameOccu
occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)

# Fit occupancy model to simulated data
# This will contain all the model structure info powerAnalysis needs
# The estimates from the model aren't used
template_model <- occu(~1~elev, occu_umf)

# If we run powerAnalysis without specifying coefs we'll get a template list
powerAnalysis(template_model)

# Set desired effect sizes to pass to coefs
effect_sizes <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

# Run power analysis and look at summary
(pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05))

# Try a smaller sample size in the study design
(pa2 <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05,
                      design=list(M=100, J=2)))

}
}
