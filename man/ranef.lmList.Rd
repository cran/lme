% $Id: ranef.lmList.Rd,v 1.2 1999/07/28 04:00:51 bates Exp $
\name{ranef.lmList}
\title{Extract lmList Random Effects}
\usage{
ranef(object, augFrame, data, which, FUN, standard, omitGroupingFactor)
}
\alias{ranef.lmList}
\alias{random.effects.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.}
 \item{augFrame}{an optional logical value. If \code{TRUE}, the returned
   data frame is augmented with variables defined in the data frame used
   to produce \code{object}; else, if \code{FALSE}, only the random effects
   are returned. Defaults to \code{FALSE}.
 }
 \item{data}{an optional data frame with the variables to be used for
   augmenting the returned data frame when \code{augFrame =
     TRUE}. Defaults to the data frame used to fit \code{object}.}
 \item{which}{an optional positive integer or character vector
   specifying which columns of the data frame used to produce \code{object}
   should be used in the augmentation of the returned data
   frame. Defaults to all variables in the data.
 }
 \item{FUN}{an optional summary function or a list of summary functions
   to be applied to group-varying variables, when collapsing the data
   by groups.  Group-invariant variables are always summarized by the
   unique value that they assume within that group. If \code{FUN} is a
   single function it will be applied to each non-invariant variable by
   group to produce the summary for that variable.  If \code{FUN} is a
   list of functions, the names in the list should designate classes of
   variables in the frame such as \code{ordered}, \code{factor}, or
   \code{numeric}.  The indicated function will be applied to any
   group-varying variables of that class.  The default functions to be
   used are \code{mean} for numeric factors, and \code{Mode} for both
   \code{factor} and \code{ordered}.  The \code{Mode} function, defined
   internally in \code{gsummary}, returns the modal or most popular
   value of the variable.  It is different from the \code{mode} function
   that returns the S-language mode of the variable.
 }
 \item{standard}{an optional logical value indicating whether the
   estimated random effects should be "standardized" (i.e. divided by
   the corresponding estimated standard error). Defaults to
   \code{FALSE}.} 
 \item{omitGroupingFactor}{an optional logical value.  When \code{TRUE}
   the grouping factor itself will be omitted from the group-wise
   summary of \code{data} but the levels of the grouping factor will
   continue to be used as the row names for the returned data frame.
   Defaults to \code{FALSE}.
 }
}
\description{
  A data frame containing the differences between the coefficients of
  the individual \code{lm} fits and the average coefficients.
}
\value{
  a data frame with the differences between the individual \code{lm}
  coefficients in \code{object} and their average.  Optionally, the
  returned data frame may be augmented with covariates summarized
  over groups or the differences may be standardized.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{lmList}}, \code{\link{fixef.lmList}}}
\examples{
library(lme)
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
ranef(fm1)
ranef(fm1, standard = TRUE)
ranef(fm1, augFrame = TRUE)
}
\keyword{models}
