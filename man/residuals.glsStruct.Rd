% $Id: residuals.glsStruct.Rd,v 1.2 1998/04/06 21:35:19 pinheiro Exp $
\name{residuals.glsStruct}
\title{Calculate glsStruct Residuals}
\usage{
fitted(object, glsFit)
}
\alias{fitted.glsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{glsStruct},
   representing a list of linear model components, such as
   \code{corStruct} and \code{varFunc} objects.}
 \item{glsFit}{an optional list with components \code{logLik}   
   (log-likelihood), \code{beta} (coefficients), \code{sigma} (standard
   deviation for error term), \code{varBeta} (coefficients' covariance
   matrix), \code{fitted} (residuals), and \code{residuals}
   (residuals). Defaults to \code{attr(object, "glsFit")}.}
}
\description{
  The residuals for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the residuals for the linear model represented by
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}
\note{
  This method function is generally only used inside \code{gls} and 
  \code{fitted.gls}.
}

\seealso{\code{\link{gls}}, \code{\link{residuals.gls}},
  \code{\link{fitted.glsStruct}} } 

\keyword{models}
