% $Id: formula.nlme.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{formula.nlme}
\title{Extract nlme Object Formula}
\usage{
formula(object)
}
\alias{formula.nlme}
\arguments{
 \item{object}{an object inheriting from class \code{nlme}, representing
   a fitted nonlinear mixed-effects model.}
}
\description{
  This method function extracts the nonlinear model formula
  associated with \code{object}. 
}
\value{
  a two-sided nonlinear formula specifying the model used to
  obtain \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlme}}}

\examples{
library(lme)
data(Soybean)
\dontrun{
fm1 <- nlme(weight ~ SSlogis(Time, Asym, xmid, scal), data = Soybean,
            fixed = Asym + xmid + scal ~ 1, start = c(18, 52, 7.5))
\dontrun{
formula(fm1)
}
}
}
\keyword{models}
