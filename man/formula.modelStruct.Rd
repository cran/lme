% $Id: formula.modelStruct.Rd,v 1.4 1998/04/02 18:53:12 pinheiro Exp $
\name{formula.modelStruct}
%\title{Extract Formula from \code{modelStruct} Object}
\title{Extract modelStruct Object Formula}
\usage{
formula(object)
}
\alias{formula.modelStruct}
\arguments{
 \item{object}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
}
\description{
  This method function extracts a formula from each of the 
  components of \code{object}, returning a list of formulas.
}
\value{
  a list with the formulas of each component of \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}}

\examples{
library(lme)
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   corStruct = corAR1(0.3))
formula(lms1)
}
\keyword{models}
