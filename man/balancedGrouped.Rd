% $Id: balancedGrouped.Rd,v 1.3 1999/06/04 13:31:40 bates Exp $
\name{balancedGrouped}
\alias{balancedGrouped}
\title{Create a groupedData object from a matrix}
\usage{
balancedGrouped(form, data, labels=NULL, units=NULL)
}
\arguments{
 \item{form}{A formula of the form \code{y ~ x | g} giving the name of
   the response, the primary covariate, and the grouping factor.}
 \item{data}{A matrix or data frame containing the values of the
   response grouped according to the levels of the grouping factor
   (rows) and the distinct levels of the primary covariate (columns).
   The \code{dimnames} of the matrix are used to construct the levels of
   the grouping factor and the primary covariate.}
 \item{labels}{an optional list of character strings giving labels for
   the response and the primary covariate.  The label for the primary
   covariate is named \code{x} and that for the response is named
   \code{y}.  Either label can be omitted.}
 \item{units}{an optional list of character strings giving the units for
   the response and the primary covariate.  The units string for the
   primary covariate is named \code{x} and that for the response is
   named \code{y}.  Either units string can be omitted.}
}
\description{
  Create a \code{groupedData} object from a data matrix.  This function
  can only be used with balanced grouped data that will be representable
  as a matrix.  The opposite conversion (\code{groupedData} to
  \code{matrix}) is performed by \code{asTable}.
}
\value{
  A balanced \code{groupedData} object.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{groupedData}}, \code{\link{isBalanced}}, \code{\link{asTable}}}
\examples{
library(lme)
data(Orthodont)
OrthoMat <- asTable( Orthodont )
Orth2 <- balancedGrouped(distance ~ age | Subject, data = OrthoMat,
    labels = list(x = "Age",
                  y = "Distance from pituitary to pterygomaxillary fissure"),
    units = list(x = "(yr)", y = "(mm)"))
Orth2[ 1:10, ]        ## check the first few entries
}
\keyword{data}
