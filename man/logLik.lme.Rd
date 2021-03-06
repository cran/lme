% $Id: logLik.lme.Rd,v 1.3 1999/07/28 04:00:51 bates Exp $
\name{logLik.lme}
\title{Log-Likelihood of an lme Object}
\usage{
logLik(object, REML)
}
\alias{logLik.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}. }
}
\description{
  If \code{REML=FALSE}, returns the log-likelihood value of the linear
  mixed-effects model represented by \code{object} evaluated at the
  estimated  coefficients; else, the restricted log-likelihood evaluated
  at the estimated coefficients is returned.
}
\value{
  the (restricted) log-likelihood of the linear mixed-effects model
  represented by \code{object} evaluated at the estimated coefficients. 
}
\references{
  Harville, D.A. (1974) "Bayesian Inference for Variance Components
  Using Only Error Contrasts", Biometrika, 61, 383-385.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lme(distance ~ Sex * age, Orthodont, random = ~ age, method = "ML")
logLik(fm1)
logLik(fm1, REML = TRUE)
}
\keyword{models}
