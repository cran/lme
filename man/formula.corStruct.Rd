% $Id: formula.corStruct.Rd,v 1.5 1998/04/02 18:53:12 pinheiro Exp $
\name{formula.corStruct}
\title{Extract corStruct Object Formula}
\usage{
formula(object)
}
\alias{formula.corStruct}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct}
   representing a correlation structure.}
}
\description{
  This method function extracts the formula associated with a
  \code{corStruct} object, in which the covariate and the grouping
  factor, if any is present, are defined. 
}
\value{
  an object of class \code{formula} specifying the covariate and the
  grouping factor, if any is present, associated with \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}}

\examples{
library(lme)
cs1 <- corCAR1(form = ~ Time | Mare)
formula(cs1)
}
\keyword{models}
