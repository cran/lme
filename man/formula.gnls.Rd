% $Id: formula.gnls.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{formula.gnls}
\title{Extract gnls Object Formula}
\usage{
formula(object)
}
\alias{formula.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
}
\description{
  This method function extracts the nonlinear model formula
  associated with \code{object}. 
}
\value{
  a two-sided formula specifying the nonlinear model used to
  obtain \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}}

\examples{
library(lme)
data(Soybean)
\dontrun{
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
formula(fm1)
}
}
\keyword{models}
