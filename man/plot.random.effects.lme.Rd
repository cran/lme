% $Id: plot.random.effects.lme.Rd,v 1.3 1998/04/06 21:35:18 pinheiro Exp $
\name{plot.random.effects.lme}
\title{Plot a random.effects.lme Object}
\usage{
plot(object, outer, omitFixed, level,  ...)
}
\alias{plot.random.effects.lme}
\arguments{
 \item{object}{an object inheriting from class
   \code{random.effects.lme}, representing the estimated coefficients or
   estimated random effects for the \code{lme} object from which it was
   produced.} 
 \item{outer}{an optional one-sided formula specifying the grouping
   structure to be used for displaying the random effects
   (coefficients). Single factors (\code{~g}) or crossed factors
   (\code{~g1*g2}) are allowed. Default is \code{NULL}, in which case
   the row names of the random effects (coefficients) are used.}
 \item{omitFixed}{an optional logical value indicating whether
   columns with values that are constant across groups should be
   omitted. Default is \code{TRUE}.}
 \item{level}{an optional integer value giving the level of grouping
   to be used for \code{object}. Only used when \code{object} is a list
   with different components for each grouping level. Defaults to the
   highest or innermost level of grouping.}
 \item{\dots}{optional arguments passed to the Trellis \code{dotplot}
   function.} 
}
\description{
  A Trellis dot-plot of the random effects is generated, with a different
  panel for each random effect (coefficient). Rows in the dot-plot are
  determined by the \code{outer} argument (if not missing) or by the row
  names of the random effects (coefficients). If a single factor is
  specified in \code{outer}, its levels determine the dot-plot rows
  (with possibly multiple dots per row); otherwise, if \code{outer}
  specifies a crossing of factors, the dot-plot rows are determined by
  all combinations of the levels of the individual factors in the
  formula. The Trellis function \code{dotplot} is used in this method
  function. 
}
\value{
  a Trellis dot-plot of the estimated random-effects (coefficients).
}

\author{Jose Pinheiro and Douglas Bates}

\note{
  This function requires the \code{trellis} library.
}

\seealso{\code{\link{random.effects.lme}}, \code{\link{lme}},
  \code{\link{dotplot}}} 

\examples{
library(lme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
plot(ranef(fm1))
plot(ranef(fm1, aug = TRUE), outer = ~Sex)
}
\keyword{models}

