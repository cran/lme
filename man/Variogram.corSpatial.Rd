% $Id: Variogram.corSpatial.Rd,v 1.1 1999/05/13 15:44:54 pinheiro Exp $
\name{Variogram.corSpatial}
\title{Calculate Semi-variogram for a corSpatial Object}
\usage{
Variogram(object, distance, sig2, length.out, FUN)
}
\alias{Variogram.corSpatial}
\arguments{
  \item{object}{an object inheriting from class \code{corSpatial},
    representing spatial correlation structure.}
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
  \item{FUN}{a function of two arguments, the distance and the range
    corresponding to \code{object}, specifying the semi-variogram
    model.} 
}
\description{
  This method function calculates the semi-variogram values
  corresponding to the model defined in \code{FUN}, using the estimated
  coefficients corresponding to \code{object}, at the distances defined
  by \code{distance}.
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

\seealso{\code{\link{corSpatial}}, \code{\link{Variogram.default}},
  \code{\link{Variogram.corExp}}, \code{\link{Variogram.corGaus}},
  \code{\link{Variogram.corLin}}, \code{\link{Variogram.corRatio}},
  \code{\link{Variogram.corSpher}}, \code{\link{plot.Variogram}}}
\examples{
library(lme)
data(BodyWeight)
cs1 <- corExp(3, form = ~ Time | Rat)
cs1 <- initialize(cs1, BodyWeight)
Variogram.corSpatial(cs1, FUN = function(x, y) (1 - exp(-x/y)))[1:10,]
}
\keyword{models}
