% $Id: covariate<-.varFunc.Rd,v 1.3 1999/07/28 04:00:51 bates Exp $
\name{covariate<-.varFunc}
\title{Assign varFunc Covariate}
\usage{
covariate(object) <- value
}
\alias{covariate<-.varFunc}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{value}{a value to be assigned to the covariate associated with
   \code{object}.}
}
\description{
  The covariate(s) used in the calculation of the weights of the
  variance function represented by \code{object} is (are) replaced by
  \code{value}. If \code{object} has been initialized, \code{value} must
  have the same dimensions as \code{getCovariate(object)}.
}
\value{
  a \code{varFunc} object similar to \code{object}, but with its
  \code{covariate} attribute replaced by \code{value}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getCovariate.varFunc}}}

\examples{
library(lme)
data(Orthodont)
vf1 <- varPower(1.1, form = ~age)
covariate(vf1) <- Orthodont[["age"]]
}
\keyword{models}
