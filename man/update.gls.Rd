% $Id: update.gls.Rd,v 1.5 1999/07/28 04:00:51 bates Exp $
\name{update.gls}
\title{Update a gls Object}
\usage{
update(object, model, data, correlation, weights, subset, method,
       na.action, control)
}

\alias{update.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{other arguments}{defined as in \code{\link{gls}}. See
     that function's documentation for descriptions of and default values
     for these arguments.}
}
\description{
  The non-missing arguments in the call to the \code{update.gls} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{gls} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{gls} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}}

\examples{
library(lme)
data(Ovary)
\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
fm2 <- update(fm1, weights = varPower())
}
}
\keyword{models}
