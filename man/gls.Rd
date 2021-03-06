% $Id: gls.Rd,v 1.8 1999/07/28 04:00:51 bates Exp $
\name{gls}
\title{Fit Linear Model Using Generalized Least Squares}
\usage{
gls(model, data, correlation, weights, subset, method, na.action,
    control, verbose)
}
\alias{gls}
\arguments{
  \item{model}{a two-sided linear formula object describing the
    model, with the response on the left of a \code{~} operator and the
    terms, separated by \code{+} operators, on the right.}
 \item{data}{an optional data frame containing the variables named in
   \code{model}, \code{correlation}, \code{weights}, and
   \code{subset}. By default the variables are taken from the
   environment from which \code{gls} is called.}
 \item{correlation}{an optional \code{corStruct} object describing the
   within-group correlation structure. See the documentation of
   \code{corClasses} for a description of the available \code{corStruct}
   classes. If a grouping variable is to be used, it must be specified in
   the \code{form} argument to the \code{corStruct}
   constructor. Defaults to \code{NULL}, corresponding to uncorrelated 
   errors.}  
 \item{weights}{an optional \code{varFunc} object or one-sided formula
   describing the within-group heteroscedasticity structure. If given as
   a formula, it is used as the argument to \code{varFixed},
   corresponding to fixed variance weights. See the documentation on
   \code{varClasses} for a description of the available \code{varFunc}
   classes. Defaults to \code{NULL}, corresponding to homoscesdatic
   errors.} 
 \item{subset}{an optional expression indicating which subset of the rows of
   \code{data} should  be  used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{method}{a character string.  If \code{"REML"} the model is fit by
   maximizing the restricted log-likelihood.  If \code{"ML"} the
   log-likelihood is maximized.  Defaults to \code{"REML"}.}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{gls} to print an error message and terminate if there are any
   incomplete observations.}
 \item{control}{a list of control values for the estimation algorithm to
   replace the default values returned by the function \code{glsControl}.
   Defaults to an empty list.}
 \item{verbose}{an optional logical value. If \code{TRUE} information on
   the evolution of the iterative algorithm is printed. Default is
   \code{FALSE}.}
}
\description{
  This function fits a linear model using generalized least
  squares. The errors are allowed to be correlated and/or have unequal
  variances.  
}
\value{
  an object of class \code{gls} representing the linear model
  fit. Generic functions such as \code{print}, \code{plot}, and 
  \code{summary} have methods to show the results of the fit. See
  \code{glsObject} for the components of the fit. The functions
  \code{resid}, \code{coef}, and \code{fitted} can be used to extract
  some of its components. 
}
\references{
 The different correlation structures available for the
 \code{correlation} argument are described in Box, G.E.P., Jenkins,
 G.M., and Reinsel G.C. (1994), Littel, R.C., Milliken, G.A., Stroup,
 W.W., and Wolfinger, R.D. (1997), and Venables, W.N. and Ripley,
 B.D. (1997). The use of variance functions for linear 
 and nonlinear models is presented in detail in Carroll, R.J. and Ruppert,
 D. (1988) and Davidian, M. and Giltinan, D.M. (1995).  

 Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
 Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 

 Carroll, R.J. and Ruppert, D. (1988) "Transformation and Weighting in
 Regression", Chapman and Hall.

 Davidian, M. and Giltinan, D.M. (1995) "Nonlinear Mixed Effects Models
 for Repeated Measurement Data", Chapman and Hall.

 Littel, R.C., Milliken, G.A., Stroup, W.W., and Wolfinger, R.D. (1997)
 "SAS Systems for Mixed Models", SAS Institute.

 Venables, W.N. and Ripley, B.D. (1997) "Modern Applied Statistics with
 S-plus", 2nd Edition, Springer-Verlag.
}

\author{Jose Pinheiro and Douglas Bates}
\seealso{
  \code{\link{glsControl}}, \code{\link{glsObject}},
  \code{\link{varFunc}}, \code{\link{corClasses}}, \code{\link{varClasses}}
}
\examples{
library(lme)
data(Ovary)
\dontrun{
# AR(1) errors within each Mare
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
# variance increases as a power of the absolute fitted values
fm2 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           weights = varPower())
}
}
\keyword{models}
