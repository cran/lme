% $Id: lmList.groupedData.Rd,v 1.4 1999/07/28 04:00:51 bates Exp $
\name{lmList.groupedData}
\title{lmList Fit from a groupedData Object}
\usage{
lmList(object, data, level, na.action, pool)
}
\alias{lmList.groupedData}
\arguments{
    \item{object}{a data frame inheriting from class \code{groupedData}.}
    \item{data}{this argument is included for consistency with the generic
	function. It is ignored in this method function.}
    \item{other arguments}{identical to the arguments in the generic
	function call. See the documentation for \code{\link{lmList}}.}
}
\description{
  The response variable and primary covariate in \code{formula(object)}
  are used to construct the linear model formula. This formula
  and the \code{groupedData} object are passed as the \code{object} and
  \code{data} arguments to \code{lmList.formula}, together with any other
  additional arguments in the function call. See the documentation on
  \code{lmList.formula} for a description of that function.
}

\value{
  a list of \code{lm} objects with as many components as the number of
  groups defined by the grouping factor. Generic functions such as
  \code{coef}, \code{fixed.effects}, \code{lme}, \code{pairs},
  \code{plot}, \code{predict}, \code{random.effects}, \code{summary},
  and \code{update} have methods that can be applied to an \code{lmList}
  object. 
}

\seealso{
  \code{\link{groupedData}}, \code{\link{lm}}, \code{\link{lme.lmList}},
  \code{\link{lmList.formula}}
}
\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(Orthodont)
}
\keyword{models}
