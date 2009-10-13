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
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{species}{
%%     ~~Describe \code{species} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{csvToUMF}}
}


