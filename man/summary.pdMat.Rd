% $Id: summary.pdMat.Rd,v 1.4 1999/07/28 04:00:51 bates Exp $
\name{summary.pdMat}
\title{Summarize a pdMat Object}
\usage{
summary(object, structName, noCorrelation)
}
\alias{summary.pdMat}
\alias{summary.pdBlocked}
\alias{summary.pdCompSymm}
\alias{summary.pdDiag}
\alias{summary.pdIdent}
\alias{summary.pdNatural}


\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.} 
 \item{structName}{an optional character string with a description of
   the \code{pdMat} class. Default depends on the method function:
   \code{"Blocked"} for \code{pdBlocked},
   \code{"Compound Symmetry"} for \code{pdCompSymm}, \code{"Diagonal"}
   for \code{pdDiag}, \code{"Multiple of an Identity"} for
   \code{pdIdent}, \code{"General Positive-Definite, Natural
     Parametrization"} for \code{pdNatural}, \code{"General
     Positive-Definite"} for \code{pdSymm}, and
     \code{data.class(object)} for \code{pdMat}.}
 \item{noCorrelation}{an optional logical value indicating whether
   correlations are to be printed in \code{print.summary.pdMat}. Default
   depends on the method function: \code{FALSE} for \code{pdDiag} and
   \code{pdIdent}, and \code{TRUE} for all other classes.}
}
\description{
  Attributes \code{structName} and \code{noCorrelation}, with the values
  of the corresponding arguments to the method function, are appended to
  \code{object} and its class is changed to \code{summary.pdMat}. 
}
\value{
  an object similar to \code{object}, with additional attributes
  \code{structName} and \code{noCorrelation}, inheriting from class
  \code{summary.pdMat}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{print.summary.pdMat}}, \code{\link{pdMat}}}

\examples{
library(lme)
summary(pdSymm(diag(4)))
}
\keyword{models}
