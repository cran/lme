% $Id: getResponse.lmList.Rd,v 1.3 1998/10/14 19:49:57 pinheiro Exp $
\name{getResponse.lmList}
\title{Extract lmList Object Response}
\usage{
getResponse(object)
}
\alias{getResponse.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The response vectors from each of the \code{lm} components of
  \code{object} are extracted and combined into a single vector.  
}
\value{
  a vector with the response vectors corresponding to the \code{lm}
  components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
getResponse(fm1)
}
\keyword{models}
