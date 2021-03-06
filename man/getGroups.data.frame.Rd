% $Id: getGroups.data.frame.Rd,v 1.3 1998/04/02 20:07:46 bates Exp $
\name{getGroups.data.frame}
\title{Extract Groups from a Data Frame}
\usage{
getGroups(object, form, level)
}
\alias{getGroups.data.frame}
\arguments{
 \item{object}{an object inheriting from class \code{data.frame}.}
 \item{form}{an optional formula with a conditioning expression on its
   right hand side (i.e. an expression involving the \code{|}
   operator). Defaults to \code{formula(object)}.}  
 \item{level}{a positive integer vector with the level(s) of grouping to
   be used when multiple nested levels of grouping are present. Defaults
   to all levels of nesting.} 
}
\description{
  Each variable named in the expression after the \code{|} operator on
  the right hand side of \code{form} is evaluated in \code{object}. If
  more than one variable is indicated in \code{level} they are combined
  into a data frame; else the selected variable is returned as a vector.
  When multiple grouping levels are defined in \code{form} and
  \code{level > 1}, the levels of the returned factor are obtained by
  pasting together the levels of the grouping factors of level greater
  or equal to \code{level}, to ensure their uniqueness. 
}
\value{
  either a data frame with columns given by the grouping factors
  indicated in \code{level}, from outer to inner, or, when a single
  level is requested, a factor representing the selected grouping
  factor. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{getGroupsFormula}}}
\examples{
library(lme)
data(Pixel)
getGroups(Pixel)
getGroups(Pixel, level = 2)
}
\keyword{models}
