% $Id: update.nlme.Rd,v 1.3 1999/07/28 04:00:51 bates Exp $
\name{update.nlme}
\title{Update an nlme Object}
\usage{
update(object, model, data, fixed, random, groups, start,
       correlation, weights, subset, method, na.action, naPattern,
       control, verbose)
}
\alias{update.nlme}
\arguments{
 \item{object}{an object inheriting from class \code{nlme}, representing
   a fitted nonlinear mixed-effects model.}
 \item{other arguments}{defined as in \code{\link{nlme}}. See that
     function's documentation for descriptions of and default values
     for these arguments.}
}
\description{
  The non-missing arguments in the call to the \code{update.nlme} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{nlme} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{nlme} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlme}}}

\examples{
library(lme)
data(Soybean)
\dontrun{
fm1 <- nlme(weight ~ SSlogis(Time, Asym, xmid, scal), data = Soybean,
            fixed = Asym + xmid + scal ~ 1, start = c(18, 52, 7.5))
fm2 <- update(fm1, weights = varPower())
}
}
\keyword{models}
