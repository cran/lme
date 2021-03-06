% $Id: plot.Variogram.Rd,v 1.4 1999/07/29 18:51:45 bates Exp $
\name{plot.Variogram}
\title{Plot a Variogram Object}
\usage{
plot(object, smooth, showModel, sigma, span, xlab, ylab, type, ylim, ...)
}
\alias{plot.Variogram}
\arguments{
 \item{object}{an object inheriting from class \code{Variogram},
   consisting of a data frame with two columns named \code{variog} and
   \code{dist}, representing the semi-variogram values and the corresponding
   distances.
 }
 \item{smooth}{an optional logical value controlling whether a
   \code{loess} smoother should be added to the plot. Defaults to
   \code{TRUE}, when \code{showModel} is \code{FALSE}.
 }
 \item{showModel}{an optional logical value controlling whether the
   semi-variogram corresponding to an \code{"modelVariog"} attribute of
   \code{object}, if any is present, should be added to the
   plot. Defaults to \code{TRUE}, when the \code{"modelVariog"}
   attribute is present.
 }
 \item{sigma}{an optional numeric value used as the height of a
   horizontal line displayed in the plot. Can be used to represent the
   process standard deviation. Default is \code{NULL}, implying that no
   horizontal line is drawn. 
 }
 \item{span}{an optional numeric value with the smoothing parameter for
   the \code{loess} fit. Default is 0.6.
 }
 \item{xlab,ylab}{optional character strings with the x- and y-axis
   labels. Default respectively to \code{"Distance"} and
   \code{"SemiVariogram"}. 
 }
 \item{type}{an optional character indicating the type of plot. Defaults
     to \code{"p"}.
 }
 \item{ylim}{an optional numeric vector with the limits for the
   y-axis. Defaults to \code{c(0, max(object$variog))}.
 }
 \item{\dots}{optional arguments passed to the Trellis \code{xyplot} function.}
}
\description{
  an \code{xyplot} of the semi-variogram versus the distances is
  produced. If \code{smooth = TRUE}, a \code{loess} smoother is added to
  the plot. If \code{showModel = TRUE} and \code{object} includes an
  \code{"modelVariog"} attribute, the corresponding semi-variogram
  is added to the plot.
}
\value{
  an \code{xyplot} Trellis plot.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{Variogram}}, \code{\link{xyplot}}, \code{\link{loess}}}
\examples{
library(lme)
data(Ovary)
\dontrun{
fm1 <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary)
plot(Variogram(fm1, form = ~ Time | Mare, maxDist = 0.7))
}
}
\keyword{models}
