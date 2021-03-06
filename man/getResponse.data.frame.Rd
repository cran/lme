% $Id: getResponse.data.frame.Rd,v 1.4 1998/04/02 18:53:13 pinheiro Exp $
\name{getResponse.data.frame}
\title{Extract Response from a Data Frame}
\usage{
getResponse(object, form)
}
\alias{getResponse.data.frame}
\arguments{
 \item{object}{an object inheriting from class \code{data.frame}.}
 \item{form}{an optional formula specifying the response to be
   evaluated in \code{object}. Defaults to \code{formula(object)}.}
}
\description{
  The left hand side of \code{form} is evaluated in \code{object}.  
}
\value{
  the value of the left hand side of \code{form} evaluated in
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getResponseFormula}}}

\examples{
library(lme)
data(Orthodont)
getResponse(Orthodont)
}
\keyword{models}
