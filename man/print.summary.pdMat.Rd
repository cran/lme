% $Id: print.summary.pdMat.Rd,v 1.4 1999/06/28 17:22:25 bates Exp $
\name{print.summary.pdMat}
\title{Print a summary.pdMat Object}
\usage{
print(x, sigma, rdig, Level, resid, ...)
}
\alias{print.summary.pdMat}
\arguments{
 \item{x}{an object inheriting from class \code{summary.pdMat},
   generally resulting from applying \code{summary} to an object
   inheriting from class \code{pdMat}.}
 \item{sigma}{an optional numeric value used as a multiplier for
   the square-root factor of the positive-definite matrix represented by
   \code{object} (usually the estimated within-group standard deviation
   from a mixed-effects model). Defaults to 1.} 
 \item{rdig}{an optional integer value with the number of significant
   digits to be used in printing correlations. Defaults to 3.}
 \item{Level}{an optional character string with a description of the
   grouping level associated with \code{object} (generally corresponding
   to levels of grouping in a mixed-effects model). Defaults to NULL.}
 \item{resid}{an optional logical value. If \code{TRUE} an extra row
   with the \code{"residual"} standard deviation given in \code{sigma}
   will be included in the output. Defaults to \code{FALSE}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  The standard deviations and correlations associated with the
  positive-definite matrix represented by \code{object} (considered as a
  variance-covariance matrix) are printed, together with the formula and
  the grouping level associated \code{object}, if any are present.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{summary.pdMat}},\code{\link{pdMat}}}

\examples{
library(lme)
data(Orthodont)
\dontrun{
pd1 <- pdCompSymm(3 * diag(3) + 1, form = ~age + age^2,
         data = Orthodont)
print(summary(pd1), sigma = 1.2, resid = TRUE)
}
}
\keyword{models}
