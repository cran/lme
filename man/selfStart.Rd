% $Id: selfStart.Rd,v 1.5 1999/06/04 13:31:40 bates Exp $
\name{selfStart}
\title{Construct Self-starting Nonlinear Models}
\usage{
selfStart(model, initial, parameters, template)
}
\alias{selfStart}
\description{
  This function is generic; methods functions can be written to handle
  specific classes of objects. Available methods include
  \code{selfStart.default} and \code{selfStart.formula}.
  See the documentation on the appropriate method function.
}
\value{
  a function object of the \code{selfStart} class.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{selfStart.default}}, \code{\link{selfStart.formula}}
}
\examples{
## see documentation for the methods
}    
\keyword{models}
