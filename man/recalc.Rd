% $Id: recalc.Rd,v 1.3 1998/04/02 18:53:18 pinheiro Exp $
\name{recalc}
\title{Recalculate Condensed Linear Model Object}
\usage{
recalc(object, conLin)
}
\alias{recalc}
\arguments{
 \item{object}{any object which induces a recalculation of the condensed
 linear model object \code{conLin}.}
 \item{conLin}{a condensed linear model object, consisting of a list
   with components \code{"Xy"}, corresponding to a regression matrix
   (\code{X}) combined with a response vector (\code{y}), and
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying model.} 
 \item{\dots}{some methods for this generic function may require
   additional arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{corStruct}, \code{modelStruct},
  \code{reStruct}, and \code{varFunc}.
}
\value{
  the recalculated condensed linear model object.
}

\author{Jose Pinheiro and Douglas Bates}

\note{This function is only used inside model fitting functions
  which require recalculation of a condensed linear model object, like
  \code{lme} and \code{gls}.} 

\examples{
## see the method function documentation
}
\keyword{models}