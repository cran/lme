% $Id: print.corStruct.Rd,v 1.4 1998/04/02 18:53:17 pinheiro Exp $
\name{print.corStruct}
\title{Print a corStruct Object}
\usage{
print(x, ...)
}
\alias{print.corStruct}
\alias{print.corHF}
\alias{print.corSymm}

\arguments{
 \item{x}{an object inheriting from class \code{corStruct}, representing
   a correlation structure.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  If \code{x} has been initialized, its coefficients are printed in
  constrained form. 
}
\value{
  the printed coefficients of \code{x} in constrained form.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{print.default}}, \code{\link{coef.corStruct}}}

\examples{
library(lme)
cs1 <- corAR1(0.3)
print(cs1)
}
\keyword{models}
