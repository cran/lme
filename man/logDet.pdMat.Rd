% $Id: logDet.pdMat.Rd,v 1.4 1999/07/28 04:00:51 bates Exp $
\name{logDet.pdMat}
\title{Extract Log-Determinant from a pdMat Object}
\usage{
logDet(object)
}
\alias{logDet.pdMat}
\alias{logDet.pdBlocked}
\alias{logDet.pdCompSymm}
\alias{logDet.pdDiag}
\alias{logDet.pdIdent}
\alias{logDet.pdSymm}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.} 
}
\description{
  This method function extracts the logarithm of the determinant of a
  square-root factor of the positive-definite matrix represented by
  \code{object}.  
}
\value{
  the log-determinant of a square-root factor of the positive-definite
  matrix represented by \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{pdMat}}}

\examples{
library(lme)
pd1 <- pdSymm(diag(1:3))
logDet(pd1)
}
\keyword{models}
