% $Id: SSasympOff.Rd,v 1.3 1999/06/04 16:14:59 pinheiro Exp $
\name{SSasympOff}
\title{Alternative Asymptotic Regression Model with an Offset}
\usage{
SSasympOff(input, Asym, lrc, c0)
}
\alias{SSasympOff}
\arguments{
 \item{input}{a numeric vector of values at which to evaluate the model.}
 \item{Asym}{a numeric parameter representing the horizontal asymptote on
   the right side (very large values of \code{input}).}
 \item{lrc}{a numeric parameter representing the natural logarithm of
   the rate constant.}
 \item{c0}{a numeric parameter representing the \code{input} when
   response is zero.}
}
\description{
  This \code{selfStart} model evaluates the alternative asymptotic
  regression function and its gradient.  It has an \code{initial}
  attribute that will evaluate initial estimates of the parameters
  \code{Asym}, \code{lrc}, and \code{c0} for a given set of data.
}
\value{
  a numeric vector of the same length as \code{input}.  It is the value of
  the expression \code{Asym*(1 - exp(-exp(lrc)*(input - c0)))}.  If all of
  the arguments \code{Asym}, \code{lrc}, and \code{c0} are
  names of objects, as opposed to expressions or explicit numerical
  values, the gradient matrix with respect to these names is attached as
  an attribute named \code{gradient}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{nls}}, \code{\link{selfStart}}
}
\examples{
library( lme )
data( CO2 )
CO2.Qn1 <- CO2[CO2$Plant == "Qn1", ]
SSasympOff( CO2.Qn1$conc, 32, -4, 43 )  # response only
Asym <- 32
lrc <- -4
c0 <- 43
SSasympOff( CO2.Qn1$conc, Asym, lrc, c0 ) # response and gradient
}
\keyword{models}
