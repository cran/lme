% $Id: corMatrix.pdMat.Rd,v 1.3 1998/04/02 18:53:11 pinheiro Exp $
\name{corMatrix.pdMat}
\title{Extract Correlation Matrix from a pdMat Object}
\usage{
corMatrix(object)
}
\alias{corMatrix.pdBlocked}
\alias{corMatrix.pdCompSymm}
\alias{corMatrix.pdDiag}
\alias{corMatrix.pdIdent}
\alias{corMatrix.pdMat}
\alias{corMatrix.pdSymm}

\arguments{
 \item{object}{an object inheriting from class code{pdMat}, representing
   a positive definite matrix.}
}

\description{
  The correlation matrix corresponding to the positive-definite matrix
  represented by \code{object} is obtained. 
}
\value{
  the correlation matrix corresponding to the positive-definite matrix
  represented by \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.pdMat}}, \code{\link{pdMatrix}}}

\examples{
library(lme)
pd1 <- pdSymm(diag(1:4))
corMatrix(pd1)
}
\keyword{models}
