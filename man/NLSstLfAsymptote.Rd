% $Id: NLSstLfAsymptote.Rd,v 1.2 1999/06/04 13:31:40 bates Exp $
\name{NLSstLfAsymptote}
\title{Guess the horizontal asymptote on the left side}
\alias{NLSstLfAsymptote}
\alias{NLSstLfAsymptote.sortedXyData}
\usage{
NLSstLfAsymptote(xy)
}
\arguments{
 \item{xy}{a \code{sortedXyData} object}
}
\description{
  This function provides an initial guess at the horizontal asymptote on
  the left side (smaller values of \code{x}) of the graph of \code{y}
  versus \code{x} from the \code{xy} object.  It is mostly used within the
  \code{initial} functions for a self-starting nonlinear regression
  models, which are in the \code{selfStart} class.
}
\value{
  A single numeric value which is a guess at the \code{y} value that
  would be the asymptote for small \code{x}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{sortedXyData}},\code{\link{NLSstClosestX}},
  \code{\link{NLSstRtAsymptote}},\code{\link{selfStart}}
\examples{
library( lme )
data( DNase )
DNase.2 <- DNase[ DNase$Run == "2", ]
DN.srt <- sortedXyData( expression(log(conc)), expression(density), DNase.2 )
NLSstLfAsymptote( DN.srt )
}
\keyword{data, manip}
