% $Id: fixef.lme.Rd,v 1.1 1998/10/14 19:49:57 pinheiro Exp $
\name{fixef.lme}
\title{Extract lme Fixed Effects}
\usage{
fixef(object)
}
\alias{fixef.lme}
\alias{fixed.effects.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
}
\description{
  The fixed effects estimates corresponding to the linear mixed-effects
  model represented by \code{object} are returned.
}
\value{
  a vector with the fixed effects estimates corresponding to
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{coef.lme}}, \code{\link{ranef.lme}} }

\examples{
library(lme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
fixef(fm1)
}
\keyword{models}
