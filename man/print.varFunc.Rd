% $Id: print.varFunc.Rd,v 1.2 1998/04/02 18:53:17 pinheiro Exp $
\name{print.varFunc}
\title{Print a varFunc Object}
\usage{
print(x, ...)
}
\alias{print.varFunc}
\alias{print.varComb}
\arguments{
 \item{x}{an object inheriting from class \code{varFunc}, representing a
 variance function structure.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  The class and the coefficients associated with \code{x} are printed.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{summary.varFunc}},
  \code{\link{print.summary.varFunc}}} 

\examples{
library(lme)
data(Orthodont)
vf1 <- varPower(0.3, form = ~age)
vf1 <- initialize(vf1, Orthodont)
print(vf1)
}
\keyword{models}
