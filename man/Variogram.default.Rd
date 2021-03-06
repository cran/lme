% $Id: Variogram.default.Rd,v 1.2 1999/06/28 17:22:25 bates Exp $
\name{Variogram.default}
\title{Calculate Semi-variogram}
\usage{
Variogram(object, distance)
}
\alias{Variogram.default}
\arguments{
  \item{object}{a numeric vector with the values to be used for
    calculating the semi-variogram, usually a residual vector from a
    fitted model.}
  \item{distance}{a numeric vector with the pairwise distances
    corresponding to the elements of \code{object}. The order of the
    elements in \code{distance} must correspond to the pairs
    \code{(1,2), (1,3), ..., (n-1,n)}, with \code{n} representing the
    length of \code{object}, and must have length \code{n(n-1)/2}.}
}
\description{
  This method function calculates the semi-variogram for an arbitrary
  vector \code{object}, according to the distances in \code{distance}.
  For each pair of elements \eqn{x,y} in \code{object}, the
  corresponding semi-variogram is \eqn{(x-y)^2/2}. The semi-variogram is
  useful for identifying and modeling spatial correlation structures in
  observations with constant expectation and constant variance.   
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

\seealso{\code{\link{Variogram.gls}}, \code{\link{Variogram.lme}},
  \code{\link{plot.Variogram}}} 
\examples{
library(lme)
data(Ovary)
fm1 <- lm(follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time), Ovary,
          subset = Mare == 1)
\dontrun{
Variogram(resid(fm1), dist(1:29))[1:10,]
}
}
\keyword{models}
