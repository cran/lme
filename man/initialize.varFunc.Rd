% $Id: initialize.varFunc.Rd,v 1.4 1999/07/29 17:08:49 bates Exp $
\name{initialize.varFunc}
\title{Initialize varFunc Object}
\usage{
initialize(object, data, ...)
}
\alias{initialize.varFunc}
\alias{initialize.varComb}
\alias{initialize.varConstPower}
\alias{initialize.varExp}
\alias{initialize.varFixed}
\alias{initialize.varIdent}
\alias{initialize.varPower}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{data}{a data frame in which to evaluate the variables named in
   \code{formula(object)}. }
 \item{\dots}{this argument is included to make this method compatible
   with the generic.}
}
\description{
  This method initializes \code{object} by evaluating its associated
  covariate(s) and grouping factor, if any is present, in \code{data};
  determining if the covariate(s) need to be updated when the
  values of the coefficients associated with \code{object} change;
  initializing the log-likelihood and the weights associated with
  \code{object}; and assigning initial values for the coefficients in
  \code{object}, if none were present. The covariate(s) will only be
  initialized if no update is needed when \code{coef(object)} changes. 
}
\value{
  an initialized object with the same class as \code{object}
  representing a variance function structure.
}
\author{Jose Pinheiro and Douglas Bates}
\examples{
library(lme)
data(Orthodont)
vf1 <- varPower( form = ~ age | Sex )
vf1 <- initialize( vf1, Orthodont )
}
\keyword{models}
