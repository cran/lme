% $Id: getCovariate.Rd,v 1.3 1998/04/02 18:53:12 pinheiro Exp $
\name{getCovariate}
\title{Extract Covariate from an Object}
\usage{
getCovariate(object, form, data)
}
\alias{getCovariate}
\arguments{
 \item{object}{any object with a \code{covariate} component}
 \item{form}{an optional one-sided formula specifying the covariate(s)
   to be extracted. Defaults to \code{formula(object)}.}
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{form}. }
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{corStruct}, \code{corSpatial},
  \code{data.frame}, and \code{varFunc}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getCovariateFormula}}}
\examples{
## see the method function documentation
}
\keyword{models}
