% $Id: as.matrix.pdMat.Rd,v 1.6 1998/04/02 18:53:09 pinheiro Exp $
\name{as.matrix.pdMat}
\title{Matrix of a pdMat Object}
\usage{
as.matrix(x)
}
\alias{as.matrix.pdMat}
\arguments{
 \item{x}{an object inheriting from class \code{pdMat}, representing a
   positive-definite matrix.}
}
\description{
  This method function extracts the positive-definite matrix represented
  by \code{x}.
}
\value{
  a matrix corresponding to the positive-definite matrix represented by
  \code{x}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{pdMat}}, \code{\link{corMatrix}}}

\examples{
library(lme)
as.matrix(pdSymm(diag(4)))
}

\keyword{models}
