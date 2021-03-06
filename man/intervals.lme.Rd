% $Id: intervals.lme.Rd,v 1.5 1999/07/29 16:30:31 pinheiro Exp $
\name{intervals.lme}
\title{Confidence Intervals on lme Parameters}
\usage{
intervals(object, level, which)
}
\alias{intervals.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{level}{an optional numeric value with the confidence level for
   the intervals. Defaults to 0.95.} 
 \item{which}{an optional character string specifying the  subset
   of parameters for which to construct the confidence
   intervals. Possible values are \code{"all"} for all parameters,
   \code{"var-cov"} for the variance-covariance parameters only, and
   \code{"fixed"} for the fixed effects only. Defaults to \code{"all"}.}
}
\description{
  Approximate confidence intervals for the parameters in the linear
  mixed-effects model represented by \code{object} are obtained, using
  a normal approximation to the distribution of the (restricted)
  maximum likelihood estimators (the estimators are assumed to have a
  normal distribution centered at the true parameter values and with
  covariance matrix equal to the negative inverse Hessian matrix of the
  (restricted) log-likelihood evaluated at the estimated parameters).
  Confidence intervals are obtained in an unconstrained scale first,
  using the normal approximation, and, if necessary, transformed to the
  constrained scale. The \code{pdNatural} parametrization is used for
  general positive-definite matrices.
}
\value{
  a list with components given by data frames with rows corresponding to
  parameters and columns \code{lower}, \code{est.}, and \code{upper}
  representing respectively lower confidence limits, the estimated
  values, and upper confidence limits for the parameters. Possible
  components are:
  \item{fixed}{fixed effects, only present when \code{which} is not
    equal to \code{"var-cov"}.}
  \item{reStruct}{random effects variance-covariance parameters, only
    present when \code{which} is not equal to \code{"fixed"}.}
  \item{corStruct}{within-group correlation parameters, only
    present when \code{which} is not equal to \code{"fixed"} and a
    correlation structure is used in \code{object}.}
  \item{varFunc}{within-group variance function parameters, only
    present when \code{which} is not equal to \code{"fixed"} and a
    variance function structure is used in \code{object}.}
  \item{sigma}{within-group standard deviation.}
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}, \code{\link{print.intervals.lme}},
    \code{\link{pdNatural}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
intervals(fm1)
}
\keyword{models}

