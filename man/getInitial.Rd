% $Id: getInitial.Rd,v 1.1 1999/05/28 14:39:12 bates Exp $
\name{getInitial}
\title{get initial parameter estimates}
\alias{getInitial}
\alias{getInitial.default}
\alias{getInitial.formula}
\alias{getInitial.selfStart}
\usage{
getInitial(object, data, ...)
}
\arguments{
 \item{object}{a formula or a \code{selfStart} model that defines a
   nonlinear regression model}
 \item{data}{a data frame in which the expressions in the formula or
   arguments to the \code{selfStart} model can be evaluated}
 \item{\dots}{optional additional arguments}
}
\description{
  This function evaluates initial parameter estimates for a nonlinear
  regression model.  If \code{data} is a parameterized data frame or
  \code{pframe} object, its \code{parameters} attribute is returned.
  Otherwise the object is examined to see if it contains a call to a
  \code{selfStart} object whose \code{initial} attribute can be
  evaluated.
}
\value{
  A named numeric vector of starting estimates for the parameters.  The
  construction of many selfStart models is such that these "starting"
  estimates are in fact the converged parameter estimates.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{nls}},\code{\link{nlsList}},\code{\link{selfStart}},
  \code{\link{selfStart.default}},\code{\link{selfStart.formula}}
}
\examples{
library(lme)
data(Puromycin)
PurTrt <- Puromycin[ Puromycin$state == "treated", ]
\dontrun{
## check if the response should be vel or rate
getInitial( vel ~ SSmicmen( conc, Vm, K ), PurTrt )
}
}
\keyword{models,nonlinear,data}
