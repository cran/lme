% $Id: selfStart.default.Rd,v 1.3 1998/10/14 19:49:57 pinheiro Exp $
\name{selfStart.default}
\title{Construct Self-starting Nonlinear Models}
\usage{
selfStart(model, initial, parameters, template)
}
\alias{selfStart.default}
\arguments{
 \item{model}{a function object defining a nonlinear model.}
 \item{initial}{a function object, with three arguments: \code{mCall},
   \code{data}, and \code{LHS}, representing, respectively, the
   expression on the right hand side of \code{model}, a data frame in
   which to interpret the variables in \code{mCall} and \code{LHS}, and
   a name, or expression, representing the variable to be used as the
   "response" in the initial values calculations. It should return
   initial values for the parameters on the right hand side of
   \code{model}. 

   NOTE: when a self-starting model is model function is used in
   \code{nls}, the  response vector (corresponding to the left hand side
   of the nonlinear formula in \code{nls}) is made available in
   \code{frame} 1 under the name \code{.nls.initial.response} and can be
   referenced from within \code{initial}, as in the example section
   below.} 
 \item{parameters, template}{these arguments are included to keep
   consistency with the call to the generic function, but are not used
   in the \code{default} method. See the documentation on
   \code{selfStart.formula}.}
}
\description{
  A method for the generic function `selfStart' for formula objects.
}
\value{
  a function object of class \code{selfStart}, corresponding to a
  self-starting nonlinear model function. An \code{initial} attribute
  (defined by the \code{initial} argument) is added to the function to
  calculate starting estimates for the parameters in the model
  automatically. 
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{selfStart.formula}}}
}
\examples{
library(lme)
# `first.order.log.model' is a function object defining a first order
# compartment model 
# `first.order.log.initial' is a function object which calculates initial
# values for the parameters in `first.order.log.model'

# self-starting first order compartment model
\dontrun{
SSfol <- selfStart(first.order.log.model, first.order.log.initial)
}
}
\keyword{models}


