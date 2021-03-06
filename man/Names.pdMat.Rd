% $Id: Names.pdMat.Rd,v 1.6 1998/04/02 18:53:09 pinheiro Exp $
\name{Names.pdMat}
\title{Names of a pdMat Object}
\usage{
Names(object)
Names(object) <- value
}
\alias{Names.pdMat}
\alias{Names<-.pdMat}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat},
   representing a positive-definite matrix.}
 \item{value}{a character vector with the replacement values for the
   column and row names of the matrix represented by \code{object}. It
   must have length equal to the dimension of the matrix
   represented by \code{object} and, if names have been previously
   assigned to \code{object}, it must correspond to a permutation of the
   original names.}
}
\description{
  This method function returns the fist element of the \code{Dimnames}
  attribute of \code{object}, which contains the column names of the
  matrix represented by \code{object}. 
}
\value{
  if \code{object} has a \code{Dimnames} attribute then the first
  element of this attribute is returned; otherwise \code{NULL}. 
}
\author{Jose Pinheiro and Douglas Bates}
\section{SIDE EFFECTS}{
  On the left side of an assignment, sets the \code{Dimnames} attribute
  of \code{object} to \code{list(value, value)}.
}
\seealso{\code{\link{Names}}, \code{\link{Names.pdBlocked}}}

\examples{
library(lme)
data(Orthodont)
pd1 <- pdSymm(~age, data = Orthodont)
Names(pd1)
}
\keyword{models}
