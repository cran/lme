% $Id: ranef.lmList.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{ranef.lmList}
\title{Extract lmList Random Effects}
\usage{
ranef(object)
}
\alias{ranef.lmList}
\alias{ranef.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The difference between the individual \code{lm} components
  coefficients and their average is calculated.
}
\value{
  a vector with the differences between the individual \code{lm}
  coefficients in \code{object} and their average. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{fixef.lmList}}}
 
\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
ranef(fm1)
}
\keyword{models}
