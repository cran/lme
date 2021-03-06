% $Id: getGroups.corStruct.Rd,v 1.5 1998/04/02 18:53:12 pinheiro Exp $
\name{getGroups.corStruct}
\title{Extract corStruct Groups}
\usage{
getGroups(object, form, data, level)
}
\alias{getGroups.corStruct}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct}
   representing a correlation structure.}
 \item{form}{this argument is included to make the method function
   compatible with the generic. It will be assigned the value of
   \code{formula(object)} and should not be modified.}
 \item{data}{an optional data frame in which to evaluate the variables
   defined in \code{form}, in case \code{object} is not initialized and
   the grouping factor needs to be evaluated.}
 \item{level}{this argument is included to make the method function
   compatible with the generic and is not used.}
}
\description{
  This method function extracts the grouping factor associated with
  \code{object}, if any is present.
}
\value{
  if a grouping factor is present in the correlation structure
  represented by \code{object}, the function returns the corresponding
  factor vector; else the function returns \code{NULL}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getGroups}}}

\examples{
library(lme)
data(Orthodont)
cs1 <- corAR1(form = ~ 1 | Subject)
getGroups(cs1, data = Orthodont)
}
\keyword{models}
