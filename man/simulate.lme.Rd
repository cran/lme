% $Id: simulate.lme.Rd,v 1.1 1998/11/05 23:53:06 bates Exp $
\name{simulate.lme}
\title{simulate lme models}
\usage{
simulate.lme(m1, m2, Random.seed, nsim, sigma, niterEM)
}
\alias{simulate.lme}
\alias{plot.simulate.lme}
\arguments{
  \item{m1}{
    a list containing an lme model specification.  It should contain
    components \code{fixed}, \code{data}, and \code{random}
    with values suitable for a call to \code{\link{lme}}.  This model
    will be the null model.
  }
  \item{m2}{
    a list, like \code{m1} containing a second lme model specification.
    This model will be the alternate model.  Only those parts of the
    specification that change between model \code{m1} and \code{m2} need
    to be specified.
  }
  \item{Random.seed}{
    an option vector to seed the random number generator so as to
    reproduce a simulation.  This vector should be the same form as the
    \code{.Random.seed} object.
  }
  \item{nsim}{
    an optional positive integer specifying the number of simulations to
    perform.  Defaults to 1000.
  }
  \item{sigma}{
    an optional positive numeric value specifying the value of the
    residual standard deviation to use.  Defaults to the estimate of
    sigma from the fitted model \code{m1}.
  }
  \item{niterEM}{
    integer vector of length 2 giving the number of iterations of the EM
    algorithm to apply when fitting the \code{m1} and \code{m2} to each
    simulated set of data.
  }
}
\description{
  The model \code{m1} is fit to the data.  Using
  the fitted values of the parameters, \code{nsim} new data vectors from
  this model are simulated.  Both \code{m1} and \code{m2} are fit by
  maximum likelihood (ML) and by restricted maximum likelihood (REML) to
  each of the simulated data vectors.
}
\value{
  an object of class \code{Var.lme} with components \code{null} and
  \code{alt}.  Each of these has components \code{ML} and \code{REML}
  which are matrices.  An attribute called \code{Random.seed} contains
  the seed that was used for the random number generator.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{lme}},\code{\link{plot.simulate.lme}},
}
\examples{
library(lme)
data(Orthodont)
orthSim <-
 simulate.lme(m1 = list(fixed = distance ~ age, data = Orthodont,
                        random = ~ 1 | Subject),
              m2 = list(random = ~ age | Subject))
}
\keyword{models}
