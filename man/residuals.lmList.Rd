% $Id: residuals.lmList.Rd,v 1.3 1998/10/14 19:49:57 pinheiro Exp $
\name{residuals.lmList}
\title{Extract lmList Residuals}
\usage{
residuals(object, type, subset, asList)
}
\alias{residuals.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{subset}{an optional character or integer vector naming the
   \code{lm} components of \code{object} from which the residuals 
   are to be extracted. Default is \code{NULL}, in which case all
   components are used.
 }
 \item{type}{an optional character string specifying the type of
   residuals to be extracted. Options include \code{"response"} for the
   "raw" residuals (observed - fitted), \code{"pearson"} for the
   standardized residuals (raw residuals divided by the estimated
   residual standard error) using different standard errors for each
   \code{lm} fit, and \code{"pooled.pearson"} for the standardized
   residuals using a pooled estimate of the residual standard
   error. Partial matching of arguments is used, so only the first 
   character needs to be provided. Defaults to \code{"response"}.
 }
 \item{asList}{an optional logical value. If \code{TRUE}, the returned
   object is a list with the residuals split by groups; else the
   returned value is a vector. Defaults to \code{FALSE}.
 }
}
\description{
  The residuals are extracted from each \code{lm} component of
  \code{object} and arranged into a list with as many components as
  \code{object}, or combined into a single vector.
}
\value{
  a list with components given by the residuals of each \code{lm}
  component of \code{object}, or a vector with the residuals for all
  \code{lm} components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{fitted.lmList}}}

\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
residuals(fm1)
}
\keyword{models}
