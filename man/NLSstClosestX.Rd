% $Id: NLSstClosestX.Rd,v 1.2 1999/06/04 13:31:40 bates Exp $
\name{NLSstClosestX}
\title{Find the x with the closest y}
\alias{NLSstClosestX}
\alias{NLSstClosestX.sortedXyData}
\usage{
NLSstClosestX(xy, yval)
}
\arguments{
 \item{xy}{a \code{sortedXyData} object}
 \item{yval}{the numeric value on the \code{y} scale to get close to}
}
\description{
  This function is used to determine the \code{x} value in a
  \code{sortedXyData} object that corresponds to the \code{y} value that
  is closest to \code{yval}.  The function is mostly used within the
  \code{initial} functions for a self-starting nonlinear regression
  models, which are in the \code{selfStart} class.
}
\value{
  A single numeric value which is one of the values of \code{x} in the
  \code{xy} object.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{sortedXyData}},\code{\link{NLSstLfAsymptote}},
  \code{\link{NLSstRtAsymptote}},\code{\link{selfStart}}
\examples{
library( lme )
data( DNase )
DNase.2 <- DNase[ DNase$Run == "2", ]
DN.srt <- sortedXyData( expression(log(conc)), expression(density), DNase.2 )
NLSstClosestX( DN.srt, 1.0 )
}
\keyword{data, manip}%-- one or more ...
