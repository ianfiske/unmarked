\name{sight2perpdist}
\Rdversion{1.1}
\alias{sight2perpdist}
\title{
Convert sight distance and sight angle to perpendicular distance.
}
\description{
When distance data are collected on line transects using sight distances and
sight angles, they need to be converted to perpendicular distances before 
analysis. 
}
\usage{
sight2perpdist(sightdist, sightangle)
}
\arguments{
  \item{sightdist}{
	Distance from observer
}
  \item{sightangle}{
	Angle from center line. In degrees between 0 and 180.
}
}
\value{
Perpendicular distance
}
\seealso{
\code{\link{distsamp}}
}
\examples{
round(sight2perpdist(10, c(0, 45, 90, 135, 180)))
}

