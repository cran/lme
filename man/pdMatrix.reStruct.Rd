% $Id: pdMatrix.reStruct.Rd,v 1.3 1998/11/18 00:21:34 pinheiro Exp $
\name{pdMatrix.reStruct}
\title{Extract Matrix or Square-Root Factor from Components of an
  reStruct Object} 
\usage{
pdMatrix(object, fact)
}
\alias{pdMatrix.reStruct}

\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{fact}{an optional logical value. If \code{TRUE}, square-root
   factors of the positive-definite matrices represented by the elements
   of \code{object} are returned; else, if \code{FALSE}, the
   positive-definite matrices are returned. Defaults to \code{FALSE}.} 
}
\description{
  This method function extracts the positive-definite  matrices
  corresponding to the \code{pdMat} elements of \code{object}, or
  square-root factors of the positive-definite matrices. 
}
\value{
  a list with components given by the positive-definite matrices
  corresponding to the elements of \code{object}, or square-root factors
  of the positive-definite matrices.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.reStruct}}, \code{\link{reStruct}},
  \code{\link{pdMatrix.pdMat}}}

\examples{
library(lme)
data(Orthodont)
rs1 <- reStruct(pdSymm(diag(3), ~age+Sex, data = Orthodont))
pdMatrix(rs1)
}
\keyword{models}
