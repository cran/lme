% $Id: fixef.lmList.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{fixef.lmList}
\title{Extract lmList Fixed Effects}
\usage{
fixef(object)
}
\alias{fixef.lmList}
\alias{fixef.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The average of the coefficients corresponding to the \code{lm}
  components of \code{object} is calculated.
}
\value{
  a vector with the average of the individual \code{lm} coefficients in
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{ranef.lmList}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
fixef(fm1)
}
\keyword{models}
