% $Id: formula.lmList.Rd,v 1.3 1998/10/14 19:49:57 pinheiro Exp $
\name{formula.lmList}
\title{Extract lmList Object Formula}
\usage{
formula(object)
}
\alias{formula.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  This method function extracts the common linear model formula
  associated with each \code{lm} component of \code{object}. 
}
\value{
  a two-sided linear formula specifying the linear model used to
  obtain the \code{lm} components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
formula(fm1)
}
\keyword{models}
