\name{wtab}
\docType{data}
\alias{wtab}
\title{Real dynamic graph}
\description{
 A dataset reporting the cycle hires between the Santander stations of London, on September, 9, 2015 (one day).
 One interaction is a cycle hire from one station to another one. The stations are nodes in the dynamic graphs.
}
\usage{data(wtab)}

\format{ 
A data.frame with 3220 rows and 3 columns.

	\describe{
		\item{Start.Date}{: the time a bike was taken from a station, expressed in seconds from the midnight.}
		\item{StartStation.Id}{: the Id of the station the bike was taken from.}
		\item{EndStation.Id}{: the Id of the station the bike was given back.}
	}
}
\keyword{datasets; london; bike-sharing}