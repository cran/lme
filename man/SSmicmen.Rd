% $Id: SSmicmen.Rd,v 1.4 1999/07/29 16:30:31 pinheiro Exp $
\name{SSmicmen}
\title{Michaelis-Menten model}
\usage{
SSmicmen(input, Vm, K)
}
\alias{SSmicmen}
\arguments{
 \item{input}{a numeric vector of values at which to evaluate the model.}
 \item{Vm}{a numeric parameter representing the maximum value of the response.}
 \item{K}{a numeric parameter representing the \code{input} value at
   which half the maximum response is attained.  In the field of enzyme
   kinetics this is called the Michaelis parameter.}
}
\description{
  This \code{selfStart} model evaluates the Michaelis-Menten model and
  its gradient.  It has an \code{initial} attribute that
  will evaluate initial estimates of the parameters \code{Vm} and \code{K}
}
\value{
  a numeric vector of the same length as \code{input}.  It is the value of
  the expression \code{Vm*input/(K+input)}.  If both
  the arguments \code{Vm} and \code{K} are
  names of objects, the gradient matrix with respect to these names is
  attached as an attribute named \code{gradient}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{nls}}, \code{\link{selfStart}}
}
\examples{
library( lme )
data( Puromycin )
PurTrt <- Puromycin[ Puromycin$state == "treated", ]
SSmicmen( PurTrt$conc, 200, 0.05 )  # response only
Vm <- 200
K <- 0.05
SSmicmen( PurTrt$conc, Vm, K ) # response and gradient
}
\keyword{models}
