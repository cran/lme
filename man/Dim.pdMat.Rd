% $Id: Dim.pdMat.Rd,v 1.6 1998/04/02 18:53:09 pinheiro Exp $
\name{Dim.pdMat}
\title{Dimensions of a pdMat Object}
\usage{
Dim(object)
}
\alias{Dim.pdMat}
\alias{Dim.pdCompSymm}
\alias{Dim.pdDiag}
\alias{Dim.pdIdent}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat},
   representing a positive-definite matrix.}
}
\description{
  This method function returns the dimensions of the matrix represented
  by \code{object}. 
}
\value{
  an integer vector with the number of rows and columns of the
  matrix represented by \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{Dim}}}
}
\examples{
library(lme)
Dim(pdSymm(diag(3)))
}
\keyword{models}
