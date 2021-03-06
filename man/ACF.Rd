% $Id: ACF.Rd,v 1.1 1999/02/01 19:06:52 pinheiro Exp $
\name{ACF}
\title{Autocorrelation Function}
\usage{
ACF(object, maxLag, ...)
}
\alias{ACF}
\arguments{
 \item{object}{any object from which an autocorrelation function can be
   obtained. Generally an object resulting from a model fit, from which
   residuals can be extracted.}
 \item{maxLag}{maximum lag for which the autocorrelation should be
   calculated.} 
 \item{\dots}{some methods for this generic require additional
   arguments.} 
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{gls} and  \code{lme}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\references{
  Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
  Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{ACF.gls}}, \code{\link{ACF.lme}}}
\examples{
## see the method function documentation
}
\keyword{models}
