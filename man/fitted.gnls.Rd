% $Id: fitted.gnls.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{fitted.gnls}
\title{Extract gnls Fitted Values}
\usage{
fitted(object)
}
\alias{fitted.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlienar least squares fitted model.}
}
\description{
  The fitted values for the nonlinear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the nonlinear model represented by
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}, \code{\link{residuals.gnls}} }

\examples{
library(nlme)
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
fitted(fm1)
}
\keyword{models}
