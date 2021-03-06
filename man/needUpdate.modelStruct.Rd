% $Id: needUpdate.modelStruct.Rd,v 1.3 1999/02/05 03:55:03 pinheiro Exp $
\name{needUpdate.modelStruct}
%\title{Check if \code{modelStruct} Object Needs Update}
\title{Check if a modelStruct Object Needs Updating}
\usage{
needUpdate(object)
}
\alias{needUpdate.modelStruct}
\arguments{
 \item{object}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
}
\description{
  This method function checks if any of the elements of \code{object}
  needs to be updated. Updating of objects usually takes place in
  iterative algorithms in which auxiliary quantities associated with the
  object, and not being optimized over, may change.
}
\value{
  a logical value indicating whether any element of \code{object} needs
  to be updated.  
}
\author{Jose Pinheiro and Douglas Bates}

\examples{
library(lme)
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   varStruct = varPower(form = ~age))
needUpdate(lms1)
}
\keyword{models}
