% $Id: Names.reStruct.Rd,v 1.7 1998/04/02 18:53:09 pinheiro Exp $
\name{Names.reStruct}
\title{Names of an reStruct Object}
\usage{
Names(object)
Names(object) <- value
}
\alias{Names.reStruct}
\alias{Names<-.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{value}{a list of character vectors with the replacement values
   for the names of the individual \code{pdMat} objects that form
   \code{object}. It must have the same length as \code{object}.}
}
\description{
  This method function extracts the column names of each of the
  positive-definite matrices represented the \code{pdMat}
  elements of \code{object}.
}
\value{
  a list containing the column names of each of the positive-definite
  matrices represented by the \code{pdMat} elements of \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\section{SIDE EFFECTS}{
  On the left side of an assignment, sets the \code{Names} of the
  \code{pdMat} elements of \code{object} to the corresponding element of
  \code{value}.  
}
\seealso{\code{\link{reStruct}}, \code{\link{pdMat}},
  \code{\link{Names.pdMat}}
} 
\examples{
library(lme)
data(Pixel)
rs1 <- reStruct(list(Dog = ~day, Side = ~1), data = Pixel)
Names(rs1)
}
\keyword{models}
