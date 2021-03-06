% $Id: pdClasses.Rd,v 1.1 1998/04/02 18:53:15 pinheiro Exp $
\name{pdClasses}
\title{Positive-Definite Matrix Classes}
\alias{pdClasses}

\description{
  Standard classes of positive-definite matrices (\code{pdMat})
  structures  available in the \code{lme} library. 
}
\value{
  Available standard classes:
  \item{pdSymm}{general positive-definite matrix, with no additional
    structure}
  \item{pdDiag}{diagonal}
  \item{pdIdent}{multiple of an identity}
  \item{pdCompSymm}{compound symmetry structure (constant diagonal and
    constant off-diagonal elements)}
  \item{pdBlocked}{block-diagonal matrix, with diagonal blocks of any
    "atomic" \code{pdMat} class}
  \item{pdNatural}{general positive-definite matrix in natural
    parametrization (i.e. parametrized in terms of standard deviations
    and correlations). The underlying coefficients are not unrestricted,
    so this class should NOT be used for optimization.}
}

\author{Jose Pinheiro and Douglas Bates}

\note{
  Users may define their own \code{pdMat} classes by specifying a
  \code{constructor} function and, at a minimum, methods for the
  functions \code{pdConstruct}, \code{pdMatrix} and \code{coef}. For
  examples of these functions, see the methods for classes \code{pdSymm}
  and \code{pdDiag}. 
}

\seealso{\code{\link{pdCompSymm}}, \code{\link{pdDiag}},
  \code{\link{pdIdent}}, \code{\link{pdNatural}}, \code{\link{pdSymm}}}

\keyword{models}
