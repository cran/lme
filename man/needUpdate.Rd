% $Id: needUpdate.Rd,v 1.4 1998/04/02 18:53:15 pinheiro Exp $
\name{needUpdate}
\title{Check if Update is Needed}
\usage{
needUpdate(object)
}
\alias{needUpdate}
\alias{needUpdate.default}
\alias{needUpdate.varComb}
\alias{needUpdate.varIdent}
\arguments{
 \item{object}{any object}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. By default, it tries to extract a
  \code{needUpdate} attribute of \code{object}. If this is \code{NULL}
  or \code{FALSE} it returns \code{FALSE}; else it returns \code{TRUE}.
  Updating of objects usually takes place in iterative algorithms in
  which auxiliary quantities associated with the object, and not being
  optimized over, may change.
}
\value{
  a logical value indicating whether \code{object} needs to be updated. 
}
\author{Jose Pinheiro and Douglas Bates}

\examples{
library(lme)
data(Orthodont)
vf1 <- varExpon()
vf1 <- initialize(vf1, data = Orthodont)
needUpdate(vf1)
}
\keyword{models}
