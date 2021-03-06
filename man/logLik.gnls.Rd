% $Id: logLik.gnls.Rd,v 1.2 1999/06/28 17:22:25 bates Exp $
\name{logLik.gnls}
\title{Log-Likelihood of a gnls Object}
\usage{
logLik(object)
}
\alias{logLik.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
}
\description{
  Returns the log-likelihood value of the nonlinear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\value{
  the log-likelihood of the linear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}}

\examples{
library(lme)
data(Soybean)
\dontrun{
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
logLik(fm1)
}
}
\keyword{models}
