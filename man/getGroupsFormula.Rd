% $Id: getGroupsFormula.Rd,v 1.4 1998/04/02 18:53:13 pinheiro Exp $
\name{getGroupsFormula}
\title{Extract Grouping Formula}
\usage{
getGroupsFormula(object, asList)
}
\alias{getGroupsFormula}
\alias{getGroupsFormula.default}
\arguments{
 \item{object}{any object from which a formula can be extracted.}
 \item{asList}{an optional logical value. If \code{TRUE} the returned
   value with be a list of formulas; else, if \code{FALSE} the returned
   value will be a one-sided formula. Defaults to \code{FALSE}.}
}
\description{
  The conditioning expression associated with \code{formula(object)}
  (i.e. an expression after a \code{|} operator) is returned either as
  a named list of one-sided formulas, or a single one-sided formula,
  depending  on the value of \code{asList}. The components of the
  returned list are ordered from outermost to innermost level and are
  named after the grouping factor expression.
}
\value{
  a one-sided formula, or a list of one-sided formulas, with the
  grouping structure associated with \code{formula(object)}. If no
  conditioning expression is present in \code{formula(object)} a
  \code{NULL} value is returned.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getGroups}}}

\examples{
library(lme)
getGroupsFormula(y ~ x | g1/g2)
}
\keyword{models}
