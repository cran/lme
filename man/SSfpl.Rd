% $Id: SSfpl.Rd,v 1.2 1999/03/31 15:41:22 pinheiro Exp $
\name{SSfpl}
\title{Four-parameter Logistic Model}
\usage{
SSfpl(input, A, B, xmid, scal)
}
\alias{SSfpl}
\arguments{
 \item{input}{a numeric vector of values at which to evaluate the model.}
 \item{A}{a numeric parameter representing the horizontal asymptote on
   the left side (very small values of \code{input}).}
 \item{B}{a numeric parameter representing the horizontal asymptote on
   the right side (very large values of \code{input}).}
 \item{xmid}{a numeric parameter representing the \code{input} value at the
   inflection point of the curve.  The value of \code{SSfpl} will be
   midway between \code{A} and \code{B} at \code{xmid}.}
 \item{scal}{a numeric scale parameter on the \code{input} axis.}
}
\description{
  This \code{selfStart} model evaluates the four-parameter logistic
  function and its gradient.  It has an \code{initial} attribute that
  will evaluate initial estimates of the parameters \code{A}, \code{B},
  \code{xmid}, and \code{scal} for a given set of data.
}
\value{
  a numeric vector of the same length as \code{input}.  It is the value of
  the expression \code{A+(B-A)/(1+exp((xmid-input)/scal))}.  If all of
  the arguments \code{A}, \code{B}, \code{xmid}, and \code{scal} are
  names of objects, as opposed to expressions or explicit numerical
  values, the gradient matrix with respect to these names is attached as
  an attribute named \code{gradient}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{nls}}, \code{\link{selfStart}}
}
\examples{
library( lme )
data( ChickWeight )
Chick.1 <- ChickWeight[ChickWeight$Chick == 1, ]
SSfpl( Chick.1$Time, 13, 368, 14, 6 )  # response only
A <- 13
B <- 368
xmid <- 14
scal <- 6
SSfpl( Chick.1$Time, A, B, xmid, scal ) # response and gradient
}
\keyword{models}
