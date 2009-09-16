\name{detFuns}
\title{Detection functions used by distsamp}
\usage{
	gxhn(x, sigma)
	gxexp(x, rate)
	gxhaz(x, shape, scale)
	grhn(r, sigma)
	grexp(r, rate)
	grhaz(r, shape, scale)
	}
\alias{detFuns}
\alias{gxhn}
\alias{gxexp}
\alias{gxhaz}
\alias{grhn}
\alias{grexp}
\alias{grhaz}
\description{ These functions represent the currently available detection
	functions used for modeling line and point transect data with 
	\code{\link{distsamp}}
	}
\arguments{
	\item{x}{Perpendicular distance}
	\item{r}{Radial distance}
	\item{sigma}{Shape parameter of half-normal detection function}
	\item{rate}{Shape parameter of negative-exponential detection function}
	\item{shape}{Shape parameter of hazard-rate detection function}
	\item{scale}{Scale parameter of hazard-rate detection function}
	}
\seealso{\code{\link{distsamp}} for example of using these for plotting 
detection function}	