% $Id: isInitialized.Rd,v 1.4 1998/04/02 18:53:14 pinheiro Exp $
\name{isInitialized}
\title{Check if Object is Initialized}
\usage{
isInitialized(object)
}
\alias{isInitialized}
\alias{isInitialized.pdMat}
\alias{isInitialized.pdBlocked}
\alias{isInitialized.reStruct}
\arguments{
 \item{object}{any object requiring initialization.}
}
\description{
  Checks if \code{object} has been initialized (generally through a call
  to \code{initialize}), by searching for components and attributes
  which are modified during initialization.
}
\value{
  a logical value indicating whether \code{object} has been
  initialized. 
}
\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{initialize}}}

\examples{
library(lme)
pd1 <- pdDiag(~age)
isInitialized(pd1)
}
\keyword{models}
