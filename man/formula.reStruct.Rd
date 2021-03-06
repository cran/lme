% $Id: formula.reStruct.Rd,v 1.3 1998/04/02 18:53:12 pinheiro Exp $
\name{formula.reStruct}
\title{Extract reStruct Object Formula}
\usage{
formula(object)
}
\alias{formula.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
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
rs1 <- reStruct(list(A = pdDiag(diag(2), ~age), B = ~1))
formula(rs1)
}
\keyword{models}
