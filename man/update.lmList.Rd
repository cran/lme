% $Id: update.lmList.Rd,v 1.3 1998/10/14 19:49:57 pinheiro Exp $
\name{update.lmList}
\title{Update an lmList Object}
\usage{
update(object, formula, data, level, subset, na.action, control, pool)
}

\alias{update.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} fitted objects.
 }
 \item{formula}{
   a two-sided linear formula with the common model for the individuals
   \code{lm} fits.
 }
 \item{other arguments}{defined as in \code{\link{lmList}}. See the
   documentation on that function for descriptions of and default values
   for these arguments.
 }
}
\description{
  The non-missing arguments in the call to the \code{update.lmList} method
  replace the corresponding arguments in the original call used to
  produce \code{object} and \code{lmList} is used with the modified call to
  produce an updated fitted object.
}
\value{
  an updated \code{lmList} object.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(Orthodont)
fm2 <- update(fm1, distance ~ I(age - 11))
}
\keyword{models}
