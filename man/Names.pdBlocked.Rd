% $Id: Names.pdBlocked.Rd,v 1.5 1998/04/02 18:53:09 pinheiro Exp $
\name{Names.pdBlocked}
\title{Names of a pdBlocked Object}
\usage{
Names(object, asList)
}
\alias{Names.pdBlocked}
\arguments{
 \item{object}{an object inheriting from class \code{pdBlocked}
   representing a positive-definite matrix with block diagonal
   structure} 
 \item{asList}{a logical value. If \code{TRUE} a \code{list} with the
   names for each block diagonal element is returned. If \code{FALSE}
   a character vector with all column names is returned. Defaults to
   \code{FALSE}.}
}
\description{
  This method function extracts the first element of the \code{Dimnames}
  attribute, which contains the column names, for each block diagonal
  element in the matrix represented by \code{object}.
}
\value{
  if \code{asList} is \code{FALSE}, a character vector with column names
  of the matrix represented by \code{object}; otherwise, if
  \code{asList} is \code{TRUE}, a list with components given by the
  column names of the individual block diagonal elements in the matrix
  represented by \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{Names}}, \code{\link{Names.pdMat}}}

\examples{
library(lme)
data(Orthodont)
pd1 <- pdBlocked(list(~Sex - 1, ~age - 1), data = Orthodont)
Names(pd1)
}
\keyword{models}
