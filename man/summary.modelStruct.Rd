% $Id: summary.modelStruct.Rd,v 1.2 1998/03/10 22:12:37 bates Exp $
\name{summary.modelStruct}
%\title{Summarize \code{modelStruct} Object}
\title{Summarize a modelStruct Object}
\usage{
summary(object)
}
\alias{summary.modelStruct}
\arguments{
 \item{object}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
}
\description{
  This method function applies \code{summary} to each element of
  \code{object}. 
}
\value{
  a list with elements given by the summarized components of
  \code{object}. The returned value is of class
  \code{summary.modelStruct}, also inheriting from the same classes as
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{print.summary.modelStruct}}}

\examples{
library(lme)
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   corStruct = corAR1(0.3))
summary(lms1)
}
\keyword{models}
