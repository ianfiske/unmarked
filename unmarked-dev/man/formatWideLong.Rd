\name{formatWideLong}
\Rdversion{1.1}
\alias{formatWide}
\alias{formatLong}

\title{
Convert between wide and long data formats.
}
\description{
Convert a data.frame between wide and long formats.
}
\usage{
formatWide(dfin, sep = ".", obsToY, type, ...)
formatLong(dfin, species = NULL, type)
}
\details{
  In order for these functions to work, the columns of \code{dfin} need to be in the
  correct order.  \code{formatLong} requires that the columns are in the
  following scheme:
  \enumerate{
    \item site name or number.
    \item date or observation number.
    \item response variable (detections, counts, etc).
    \item The remaining columns are observation-level covariates.
  }

  \code{formatWide} requires particular names for the columns.
  The column order for \code{formatWide} is
  \enumerate{
    \item (optional) site name, named \dQuote{site}.
    \item response, named \dQuote{y.1}, \dQuote{y.2}, \ldots, \dQuote{y.J}.
    \item columns of site-level covariates, each with a relevant name per column.
    \item groups of columns of observation-level covariates, each group
  having the name form \dQuote{someObsCov.1}, \dQuote{someObsCov.2},
  \ldots, 
  \dQuote{someObsCov.J}.
  }
}
\arguments{
  \item{dfin}{A data.frame to be reformatted.
}
  \item{sep}{A seperator of column names in wide format.
}
  \item{obsToY}{Optional matrix specifying relationship between covariate
  	column structure and response matrix structure.
}
  \item{type}{Type of unmarkedFrame to create?
}
  \item{species}{Character name of species response column
}
  \item{\dots}{Further arguments}
}
\value{A data.frame}

\seealso{
\code{\link{csvToUMF}}
}


