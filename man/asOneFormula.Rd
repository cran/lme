% $Id: asOneFormula.Rd,v 1.3 1998/04/02 18:53:09 pinheiro Exp $
\name{asOneFormula}
\title{Combine Formulas of a Set of Objects}
\usage{
asOneFormula(..., omit)
}
\alias{asOneFormula}
\arguments{
 \item{\dots}{objects, or lists of objects, from which a formula can be
   extracted.}
 \item{omit}{an optional character vector with the names of variables to
   be omitted from the returned formula. Defaults to c(".", "pi").}
}
\description{
  The names of all variables used in the formulas extracted from the
  objects defined in \code{...} are converted into a single linear
  formula, with the variables names separated by \code{+}. 
}
\value{
  a one-sided linear formula with all variables named in the formulas
  extracted from the objects in \code{...}, except the ones listed in
  \code{omit}. 
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}, \code{\link{all.vars}}}

\examples{
library(lme)
asOneFormula(y ~ x + z | g, list(~ w, ~ t * sin(2 * pi)))
}
\keyword{models}
