% $Id: getData.lme.Rd,v 1.1 1999/03/31 15:41:22 pinheiro Exp $
\name{getData.lme}
\title{Extract lme Object Data}
\usage{
getData(object)
}
\alias{getData.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a linear mixed-effects fitted model.}
}
\description{
  If present in the calling sequence used to produce \code{object}, the
  data frame used to fit the model is obtained. 
}
\value{
  if a \code{data} argument is present in the calling sequence that
  produced \code{object}, the corresponding data frame (with
  \code{na.action} and \code{subset} applied to it, if also present in
  the call that produced \code{object}) is returned;
  else, \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}}

\examples{
library(lme)
data(Ovary)
fm1 <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
           random = ~ sin(2*pi*Time))
getData(fm1)
}
\keyword{models}
