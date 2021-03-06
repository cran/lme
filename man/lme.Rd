% $Id: lme.Rd,v 1.10 1999/07/28 04:00:51 bates Exp $
\name{lme}
\title{Linear Mixed-Effects Models}
\usage{
lme(fixed, data, random, correlation, weights, subset, method,
    na.action, control)
}
\alias{lme}
\alias{lme.formula}
\arguments{
  \item{fixed}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right, an \code{lmList} object, or a \code{groupedData}
    object. The method functions \code{lme.lmList} and
    \code{lme.groupedData} are documented separately.}
 \item{data}{an optional data frame containing the variables named in
   \code{fixed}, \code{random}, \code{correlation}, \code{weights}, and
   \code{subset}.  By default the variables are taken from the
   environment from which \code{lme} is called.}
 \item{random}{optionally, any of the following: (i) a one-sided formula
   of the form \code{~x1+...+xn | g1/.../gm}, with \code{x1+...+xn}
   specifying the model for the random effects and \code{g1/.../gm} the
   grouping structure (\code{m} may be equal to 1, in which case no
   \code{/} is required). The random effects formula will be repeated
   for all levels of grouping, in the case of multiple levels of
   grouping; (ii) a list of one-sided formulas of the form
   \code{~x1+...+xn | g}, with possibly different random effects models
   for each grouping level. The order of nesting will be assumed the
   same as the order of the elements in the list; (iii) a one-sided
   formula of the form \code{~x1+...+xn}, or a \code{pdMat} object with
   a formula (i.e. a non-\code{NULL} value for \code{formula(object)}),
   or a list of such formulas or \code{pdMat} objects. In this case, the
   grouping structure formula will be derived from the data used to
   fit the linear mixed-effects model, which should inherit from class
   \code{groupedData}; (iv) a named list of formulas or \code{pdMat}
   objects as in (iii), with the grouping factors as names. The order of
   nesting will be assumed the same as the order of the order of the
   elements in the list; (v) an \code{reStruct} object. See the
   documentation on \code{pdClasses} for a description of the available
   \code{pdMat} classes. Defaults to a formula consisting of the right
   hand side of \code{fixed}.}  
 \item{correlation}{an optional \code{corStruct} object describing the
   within-group correlation structure. See the documentation of
   \code{corClasses} for a description of the available \code{corStruct}
   classes. Defaults to \code{NULL},
   corresponding to no within-group correlations.}
 \item{weights}{an optional \code{varFunc} object or one-sided formula
   describing the within-group heteroscedasticity structure. If given as
   a formula, it is used as the argument to \code{varFixed},
   corresponding to fixed variance weights. See the documentation on
   \code{varClasses} for a description of the available \code{varFunc}
   classes. Defaults to \code{NULL}, corresponding to homocesdatic
   within-group errors.} 
 \item{subset}{an optional expression indicating the subset of the rows of
   \code{data} that should be used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{method}{a character string.  If \code{"REML"} the model is fit by
   maximizing the restricted log-likelihood.  If \code{"ML"} the
   log-likelihood is maximized.  Defaults to \code{"REML"}.}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{lme} to print an error message and terminate if there are any
   incomplete observations.}
 \item{control}{a list of control values for the estimation algorithm to
   replace the default values returned by the function \code{lmeControl}.
   Defaults to an empty list.}
}
\description{
  This generic function fits a linear mixed-effects model in the
  formulation described in Laird and Ware (1982) but allowing for nested
  random effects. The within-group errors are allowed to be correlated
  and/or have unequal variances. 
}
\value{
  an object of class \code{lme} representing the linear mixed-effects
  model fit. Generic functions such as \code{print}, \code{plot} and
  \code{summary} have methods to show the results of the fit. See
  \code{lmeObject} for the components of the fit. The functions
  \code{resid}, \code{coef}, \code{fitted}, \code{fixed.effects}, and
  \code{random.effects}  can be used to extract some of its components.
}
\references{
 The computational methods are described in Bates, D.M. and Pinheiro
 (1998) and follow on the general framework of Lindstrom, M.J. and Bates,
 D.M. (1988). The model formulation is described in Laird, N.M. and Ware,
 J.H. (1982).  The variance-covariance parametrizations are described in
 <Pinheiro, J.C. and Bates., D.M.  (1996).   The different correlation
 structures available for the \code{correlation} argument are described
 in Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994), Littel, R.C.,
 Milliken, G.A., Stroup, W.W., and Wolfinger, R.D. (1997), and Venables,
 W.N. and Ripley, B.D. (1997). The use of variance functions for linear
 and nonlinear mixed effects models is presented in detail in Davidian,
 M. and Giltinan, D.M. (1995). 

 Bates, D.M. and Pinheiro, J.C. (1998) "Computational methods for
 multilevel models" available in PostScript or PDF formats at
 http://franz.stat.wisc.edu/pub/NLME/
 
 Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
 Analysis: Forecasting and Control", 3rd Edition, Holden-Day. 

 Davidian, M. and Giltinan, D.M. (1995) "Nonlinear Mixed Effects Models
 for Repeated Measurement Data", Chapman and Hall.

 Laird, N.M. and Ware, J.H. (1982) "Random-Effects Models for
 Longitudinal Data", Biometrics, 38, 963-974.  

 Lindstrom, M.J. and Bates, D.M. (1988) "Newton-Raphson and EM
 Algorithms for Linear Mixed-Effects Models for Repeated-Measures
 Data", Journal of the American Statistical Association, 83,
 1014-1022. 

 Littel, R.C., Milliken, G.A., Stroup, W.W., and Wolfinger, R.D. (1997)
 "SAS Systems for Mixed Models", SAS Institute.

 Pinheiro, J.C. and Bates., D.M.  (1996) "Unconstrained
 Parametrizations for Variance-Covariance Matrices", Statistics and
 Computing, 6, 289-296.

 Venables, W.N. and Ripley, B.D. (1997) "Modern Applied Statistics with
 S-plus", 2nd Edition, Springer-Verlag.
}

\author{Jose Pinheiro and Douglas Bates}
\seealso{
  \code{\link{lmeControl}}, \code{\link{lme.lmList}},
  \code{\link{lme.groupedData}}, \code{\link{lmeObject}},
  \code{\link{lmList}}, \code{\link{reStruct}}, \code{\link{reStruct}},
  \code{\link{varFunc}}, \code{\link{pdClasses}},
  \code{\link{corClasses}}, \code{\link{varClasses}}
}
\examples{
library(lme)
data(Orthodont)
fm1 <- lme(distance ~ age, data = Orthodont) # random is ~ age
fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
}
\keyword{models}
