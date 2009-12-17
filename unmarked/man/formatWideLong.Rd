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


