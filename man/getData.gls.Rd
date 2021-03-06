% $Id: getData.gls.Rd,v 1.2 1999/06/28 17:22:25 bates Exp $
\name{getData.gls}
\title{Extract gls Object Data}
\usage{
getData(object)
}
\alias{getData.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
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

\seealso{\code{\link{gls}}}

\examples{
library(lme)
data(Ovary)
\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
getData(fm1)
}
}
\keyword{models}
