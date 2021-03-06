% $Id: Variogram.corRatio.Rd,v 1.1 1999/05/13 15:44:54 pinheiro Exp $
\name{Variogram.corRatio}
\title{Calculate Semi-variogram for a corRatio Object}
\usage{
Variogram(object, distance, sig2, length.out)
}
\alias{Variogram.corRatio}
\arguments{
  \item{object}{an object inheriting from class \code{corRatio},
    representing an Rational Quadratic spatial correlation structure.}
  \item{distance}{an optional numeric vector with the distances at
    which the semi-variogram is to be calculated. Defaults to
    \code{NULL}, in which case a sequence of length \code{length.out}
    between the minimum and maximum values of
    \code{getCovariate(object)} is used.}
  \item{sig2}{an optional numeric value representing the process
    variance. Defaults to \code{1}.}
  \item{length.out}{an optional integer specifying the length of the
    sequence of distances to be used for calculating the semi-variogram,
    when \code{distance = NULL}. Defaults to \code{50}.}
}
\description{
  This method function calculates the semi-variogram values
  corresponding to the Rational Quadratic correlation model, using the
  estimated coefficients corresponding to \code{object}, at the
  distances defined by \code{distance}.
}
\value{
  a data frame with columns \code{variog} and \code{dist} representing,
  respectively, the semi-variogram values and the corresponding
  distances. The returned value inherits from class \code{Variogram}. 
}
\references{
  Cressie, N.A.C. (1993), "Statistics for Spatial Data", J. Wiley & Sons.
}  
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{corRatio}}, \code{\link{plot.Variogram}}}
\examples{
library(lme)
data(BodyWeight)
cs1 <- corRatio(7, form = ~ Time | Rat)
cs1 <- initialize(cs1, BodyWeight)
Variogram(cs1)[1:10,]
}
\keyword{models}
