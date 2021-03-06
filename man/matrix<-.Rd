% $Id: matrix<-.Rd,v 1.4 1998/10/19 18:22:10 pinheiro Exp $
\name{matrix<-}
\title{Assign Matrix Values}
\usage{
matrix(object) <- value
}
\alias{matrix<-}
\arguments{
 \item{object}{any object to which \code{as.matrix} can be applied.}
 \item{value}{a matrix, or list of matrices, with the same dimensions as
   \code{as.matrix(object)} with the new values to be assigned to the
   matrix associated with \code{object}.} 
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{pdMat}, \code{pdBlocked}, and 
  \code{reStruct}.
}
\value{
  will depend on the method function; see the appropriate documentation.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix}}}

\examples{
## see the method function documentation
}
\keyword{models}
