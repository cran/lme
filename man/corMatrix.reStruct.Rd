% $Id: corMatrix.reStruct.Rd,v 1.3 1998/11/18 00:21:34 pinheiro Exp $
\name{corMatrix.reStruct}
\title{Extract Correlation Matrix from Components of an reStruct Object}
\usage{
corMatrix(object)
}
\alias{corMatrix.reStruct}

\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
}
\description{
  This method function extracts the correlation matrices
  corresponding to the \code{pdMat} elements of \code{object}. 
}
\value{
  a list with components given by the correlation matrices
  corresponding to the elements of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.reStruct}}, \code{\link{reStruct}},
  \code{\link{pdMat}}}

\examples{
library(lme)
data(Orthodont)
rs1 <- reStruct(pdSymm(diag(3), ~age+Sex, data = Orthodont))
corMatrix(rs1)
}
\keyword{models}
