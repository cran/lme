% $Id: varComb.Rd,v 1.2 1998/04/02 18:53:20 pinheiro Exp $
\name{varComb}
\title{Combination of Variance Functions}
\usage{
varComb(...)
}
\alias{varComb}
\arguments{
 \item{\dots}{objects inheriting from class \code{varFunc} representing
   variance function structures.}
}
\description{
  This function is a constructor for the \code{varComb} class,
  representing a combination of variance functions. The corresponding
  variance function is equal to the product of the variance functions of
  the \code{varFunc} objects listed in \code{\dots}.
}
\value{
  a \code{varComb} object representing a combination of variance
  functions, also inheriting from class \code{varFunc}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{varWeights.varComb}}, \code{\link{coef.varComb}}}

\examples{
library(lme)
vf1 <- varComb(varIdent(form = ~1|Sex), varPower())
}
\keyword{models}
