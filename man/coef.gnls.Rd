% $Id: coef.gnls.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{coef.gnls}
\title{Extract gnls Coefficients}
\usage{
coef(object)
}
\alias{coef.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
}
\description{
  The estimated coefficients for the nonlinear model represented by
  \code{object} are extracted. 
}
\value{
  a vector with the estimated coefficients for the nonlinear model
  represented by \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}} 

\examples{
library(lme)
data(Soybean)
\dontrun{
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
coef(fm1)
}
}
\keyword{models}
