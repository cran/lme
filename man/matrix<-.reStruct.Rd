% $Id: matrix<-.reStruct.Rd,v 1.3 1998/04/02 18:53:15 pinheiro Exp $
\name{matrix<-.reStruct}
\title{Assign reStruct Matrices}
\usage{
matrix(object) <-  value
}
\alias{matrix<-.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{value}{a matrix, or list of matrices, with the new values to be
   assigned to the matrices associated with the \code{pdMat} components
   of \code{object}.}  
}
\description{
  The individual matrices in \code{value} are assigned to each
  \code{pdMat} component of \code{object}, in the they are listed. The
  new matrices must have the same dimensions as the matrices they are
  meant to replace.
}
\value{
  an \code{reStruct} object similar to \code{object}, but with the
  coefficients of the individual \code{pdMat} components modified to
  produce the matrices listed in \code{value}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{reStruct}}, \code{\link{pdMat}}}

\examples{
library(lme)
data(Pixel)
rs1 <- reStruct(list(Dog = ~day, Side = ~1), data = Pixel)
matrix(rs1) <- list(diag(2), 3)
}
\keyword{models}