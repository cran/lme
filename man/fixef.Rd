% $Id: fixef.Rd,v 1.2 1999/07/28 04:00:51 bates Exp $
\name{fixef}
\title{Extract Fixed Effects}
\usage{
fixef(object, ...)
fixed.effects(object, ...)
}
\alias{fixef}
\alias{fixed.effects}
\arguments{
 \item{object}{any fitted model object from which fixed effects
   estimates can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{lmList} and \code{lme}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\seealso{\code{\link{fixef.lmList}}, \code{\link{fixef.lme}}} 
\examples{
## see the method function documentation
}
\keyword{models}
