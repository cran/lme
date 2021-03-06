% $Id: summary.nlsList.Rd,v 1.5 1999/06/28 17:22:26 bates Exp $
\name{summary.nlsList}
\title{Summarize an nlsList Object}
\usage{
summary(object, pool)
}
\alias{summary.nlsList}
\arguments{
 \item{object}{an object inheriting from class \code{nlsList},
   representing a list of \code{nls} fitted objects.
 }
 \item{pool}{an optional logical value indicating whether a pooled
   estimate of the residual standard error should be used. Default is
   \code{attr(object, "pool")}.
 }
}
\description{
  The \code{summary.nls} method is applied to each \code{nls} component of
  \code{object} to produce summary information on the individual fits,
  which is organized into a list of summary statistics. The returned
  object is suitable for printing with the \code{print.summary.nlsList}
  method. 
}
\value{
  a list with summary statistics obtained by applying \code{summary.nls}
  to the elements of \code{object}, inheriting from class
  \code{summary.nlsList}. The components of \code{value} are:
  \item{call}{a list containing an image of the \code{nlsList} call that
    produced \code{object}.
  }
  \item{parameters}{a three dimensional array with summary information
    on the \code{nls} coefficients. The first dimension corresponds to
    the names of the \code{object} components, the second dimension is
    given by   \code{"Value"}, \code{"Std. Error"}, \code{"t value"},
    and \code{"Pr(>|t|)"}, corresponding, respectively, to the
    coefficient estimates and their associated standard errors,
    t-values, and p-values. The third dimension is given by the
    coefficients names.
  }
  \item{correlation}{a three dimensional array with the 
    correlations between the individual \code{nls} coefficient
    estimates. The first dimension corresponds to the names of the
    \code{object} components. The third dimension is given by the
    coefficients names. For each coefficient, the rows of the associated
    array give the correlations between that coefficient and the
    remaining coefficients, by \code{nls} component.
  }
  \item{cov.unscaled}{a three dimensional array with the unscaled
    variances/covariances for the individual \code{lm} coefficient
    estimates (giving the estimated variance/covariance for the
    coefficients, when multiplied by the estimated residual errors). The
    first dimension corresponds to the names of the \code{object}
    components. The third dimension is given by the
    coefficients names. For each coefficient, the rows of the associated
    array give the unscaled covariances between that coefficient and the
    remaining coefficients, by \code{nls} component.
  }
  \item{df}{an array with the number of degrees of freedom for the model
    and for residuals, for each \code{nls} component.
  }
  \item{df.residual}{the total number of degrees of freedom for
    residuals, corresponding to the sum of residuals df of all \code{nls}
    components.
  }
  \item{pool}{the value of the \code{pool} argument to the function.
  }
  \item{RSE}{the pooled estimate of the residual standard error.}
  \item{sigma}{a vector with the residual standard error estimates for
    the individual \code{lm} fits.
  }
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlsList}}, \code{\link{summary.nls}}}

\examples{
library(lme)
data(Soybean)
\dontrun{
fm1 <- nlsList(weight ~ SSlogis(Time, Asym, xmid, scal) | Plot, Soybean)
summary(fm1)
}
}
\keyword{models}
