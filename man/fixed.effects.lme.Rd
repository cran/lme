% $Id: fixed.effects.lme.Rd,v 1.5 1998/04/06 21:35:15 pinheiro Exp $
\name{fixed.effects.lme}
\title{Extract lme Fixed Effects}
\usage{
fixed.effects(object)
}
\alias{fixed.effects.lme}
\alias{fixef.lme}
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

\seealso{\code{\link{coef.lme}}, \code{\link{random.effects.lme}} }

\examples{
library(lme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
fixed.effects(fm1)
}
\keyword{models}
