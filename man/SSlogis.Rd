% $Id: SSlogis.Rd,v 1.2 1999/06/04 13:31:40 bates Exp $
\name{SSlogis}
\title{Logistic model}
\usage{
SSfpl(input, Asym, xmid, scal)
}
\alias{SSlogis}
\arguments{
 \item{input}{a numeric vector of values at which to evaluate the model.}
 \item{Asym}{a numeric parameter representing the asymptote.}
 \item{xmid}{a numeric parameter representing the \code{x} value at the
   inflection point of the curve.  The value of \code{SSlogis} will be
   \code{Asym/2} at \code{xmid}.}
 \item{scal}{a numeric scale parameter on the \code{input} axis.}
}
\description{
  This \code{selfStart} model evaluates the logistic
  function and its gradient.  It has an \code{initial} attribute that
  will evaluate initial estimates of the parameters \code{Asym},
  \code{xmid}, and \code{scal} for a given set of data.
}
\value{
  a numeric vector of the same length as \code{input}.  It is the value of
  the expression \code{Asym/(1+exp((xmid-input)/scal))}.  If all of
  the arguments \code{Asym}, \code{xmid}, and \code{scal} are
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
SSlogis( Chick.1$Time, 368, 14, 6 )  # response only
Asym <- 368
xmid <- 14
scal <- 6
SSlogis( Chick.1$Time, Asym, xmid, scal ) # response and gradient
}
\keyword{models}
