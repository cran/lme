### $Id: newMethods.q,v 1.41 1999/07/29 21:14:34 pinheiro Exp $
###
###      Methods for generics from newGenerics.q for some standard classes
###
### Copyright 1997, 1999 Jose C. Pinheiro <jcp$research.bell-labs.com>,
###                      Douglas M. Bates <bates$stat.wisc.edu>
###
### This file is part of the nlme library for S and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


##*## Methods for some of the generics in newGenerics.q for standard classes

AIC.logLik <-
  ## AIC for logLik objects
  function(object)
{
  -2 * (c(object) - attr(object, "df"))
}

AIC.lm <- AIC.nls <- 
  ## AIC for various fitted objects
  function(object, ...) 
{
  if((rt <- nargs()) > 1) {
    object <- list(object, ...)
    val <- lapply(object, logLik)
    val <- 
      as.data.frame(t(sapply(val, function(el) c(attr(el, "df"), AIC(el)))))
    names(val) <- c("df", "AIC")
    row.names(val) <- as.character(match.call()[-1])
    val
  } else {
    AIC(logLik(object))
  }
}

BIC.logLik <-
  ## BIC for logLiks objects
  function(object)
{
  -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
}

BIC.lm <- BIC.nls <- 
  ## BIC for various fitted objects
  function(object, ...) 
{
  if((rt <- nargs()) > 1) {
    object <- list(object, ...)
    val <- lapply(object, logLik)
    val <- 
      as.data.frame(t(sapply(val, function(el) c(attr(el, "df"), BIC(el)))))
    names(val) <- c("df", "BIC")
    row.names(val) <- as.character(match.call()[-1])
    val
  } else {
    BIC(logLik(object))
  }
}

Dim.default <- function(object) dim(object)

fitted.nls <- 
  function(object, ...)
{
  val <- as.vector(object$m$fitted())
  lab <- "Fitted values"
  if (!is.null(aux <- attr(object, "units")$y)) {
    lab <- paste(lab, aux)
  }
  attr(val, "label") <- lab
  val
}


## may need formula.nls <- function(object) unlist( x$m$formula() ) for R?
formula.nls <- function(object) as.formula(object$call$formula)

getCovariate.data.frame <-
  function(object, form = formula(object), data)
{
  ## Return the primary covariate
  if (!(inherits(form, "formula"))) {
    stop("\"Form\" must be a formula")
  }
  aux <- getCovariateFormula(form)
  if (length(all.vars(aux)) > 0) {
    eval(aux[[2]], object)
  } else {
    rep(1, dim(object)[1])
  }
}

getGroups.data.frame <-
  ## Return the groups associated with object according to form for level
  function(object, form = formula(object), level, data, sep = "/")
{
  if (!missing(data)) {
    stop( "data argument to data.frame method for getGroups doesn't make sense" )
  }
  if (inherits(form, "formula")) {
    grpForm <- getGroupsFormula(form, asList = TRUE, sep = sep)
    if (is.null(grpForm)) {
      ## will use right hand side of form as the group formula
      grpForm <- splitFormula(asOneSidedFormula(form[[length(form)]]),
                              sep = sep)
      names(grpForm) <-
        unlist( lapply( grpForm, function(el) deparse( el[[ length(el) ]] ) ) )
    }
    if (any(unlist(lapply(grpForm,
                          function(el) length(el[[length(el)]]))) != 1)) {
      stop("Invalid formula for groups")
    }
    form <- grpForm
  } else if (data.class(form) == "list") {
    if (!all(unlist(lapply(form, function(el) inherits(el, "formula"))))) {
      stop("Form must have all components as formulas") 
    }
  } else {
    stop("Form can only be a formula, or a list of formulas")
  }
  vlist <- lapply(form,
                  function(x, dat, N) {
                    val <- eval(x[[length(x)]], dat)
                    if (length(val) == 1) {             # repeat groups
                      return(as.factor(rep(val, N)))
                    } else {
                      return(pruneLevels(as.factor(val)))
                    }
                  }, dat = object, N = nrow(object))
  if (length(vlist) == 1) return(vlist[[1]]) # ignore level - only one choice
  ## make the list into a data frame with appropriate names
  value <- do.call("data.frame", vlist)
  if (missing(level)) return(value)
  if (is.character(level)) {
    nlevel <- match(level, names(flist))
    if (any(aux <- is.na(nlevel))) {
      stop(paste("Level of", level[aux],"does not match formula \"",
		 deparse(as.vector(form)), "\""))
    }
  } else {
    nlevel <- as.numeric(level)
    if (any(aux <- is.na(match(nlevel, 1:ncol(value))))) { 
      stop(paste("level of ", level[aux]," does not match formula \"", 
	       deparse(as.vector(form)), "\""))
    }
  }
  if (length(nlevel) > 1)  return(value[, nlevel]) # multicolumn selection
  if (nlevel == 1)         return(value[, 1])     # no need to do more work
  value <- value[, 1:nlevel]
  val <- as.factor(do.call("paste", c(lapply(as.list(value),
					     as.character), sep = sep)))
  if (inherits(value[, 1], "ordered")) {
    value <- value[do.call("order", value),]
    aux <- unique(do.call("paste", c(lapply(as.list(value), 
					    as.character), sep = sep)))
    return(ordered(val, aux))
  } else {
    return(ordered(val, unique(as.character(val))))
  }
}

getGroupsFormula.default <-
  ## Return the formula(s) for the groups associated with object.
  ## The result is a one-sided formula unless asList is TRUE in which case
  ## it is a list of formulas, one for each level.
  function(object, asList = FALSE, sep = "/")
{
  form <- formula(object)
  if (!inherits(form, "formula")){
    stop("\"Form\" argument must be a formula")
  }
  form <- form[[length(form)]]
  if (!((length(form) == 3) && (form[[1]] == as.name("|")))) {
    ## no conditioning expression
    return(NULL)
  } 
  ## val <- list( asOneSidedFormula( form[[ 3 ]] ) )
  val <- splitFormula(asOneSidedFormula(form[[3]]), sep = sep)
  names(val) <- unlist(lapply(val, function(el) deparse(el[[2]])))
#  if (!missing(level)) {
#    if (length(level) == 1) {
#      return(val[[level]])
#    } else {
#      val <- val[level]
#    }
#  } 
  if (asList) val
  else eval(parse(text = paste("~",  paste(names(val), collapse = sep))))
}

getResponse.data.frame <-
  function(object, form = formula(object))
{
  ## Return the response, the evaluation of the left hand side of a formula
  ## on object
  if (!(inherits(form, "formula") && (length(form) == 3))) {
    stop("\"Form\" must be a two sided formula")
  }
  eval(form[[2]], object)
}


logLik.lm <-
  ## log-likelihood for lm objects
  function(object, REML = FALSE)
{
  res <- resid(object)
  p <- object$rank
  N <- length(res) 
  if(is.null(w <- object$weights)) {	
    w <- rep(1, N)
  } else {
    excl <- w == 0			# eliminating zero weights
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  
  N <- N - p * REML
  val <- (sum(log(w)) -N * (log(2 * pi) + 1 - log(N) + log(sum(w*res^2))))/2 -
    REML * sum(log(abs(diag(object$qr$qr)[1:p])))
  attr(val, "nall") <- N + REML * p
  attr(val, "df") <- p + 1
  attr(val, "nobs") <- N
  class(val) <- "logLik"
  val
}

logLik.nls <- 
  function(object, REML = FALSE)
{
  if (REML) {
    stop("Cannot calculate REML log-likelihood for nls objects")
  }
  res <- resid(object)
  N <- length(res)
  if(is.null(w <- object$weights)) {	
    w <- rep(1, N)
  }
  val <-  -N * (log(2 * pi) + 1 - log(N) - sum(log(w)) + log(sum(w*res^2)))/2
  attr(val, "df") <- length(object[["parameters"]]) + 1
  attr(val, "nobs") <- attr(val, "nall") <- N
  class(val) <- "logLik"
  val
}

Names.formula <-
  function(object, data = list(), exclude = c("pi", "."))
{
  if (!is.list(data)) { return(NULL) }  # no data to evaluate variable names
  allV <- all.vars(object)
  allV <- allV[is.na(match(allV, exclude))]

  if (length(allV) == 0) {
    if (attr(terms(object), "intercept")) { return("(Intercept)") }
    return(NULL)
  }

  if (any(is.na(match(allV, names(data))))) { return(NULL) }
  dimnames(model.matrix(object, model.frame(object, data)))[[2]]
}

Names.listForm <-
  function(object, data = list(), exclude = c("pi", "."))
{
  pnames <- as.character(unlist(lapply(object, "[[", 2)))
  nams <- lapply(object, function(el, data, exclude) {
    Names(getCovariateFormula(el), data, exclude)
    }, data = data, exclude = exclude)
  if (is.null(nams[[1]])) return(NULL)
  val <- c()
  for(i in seq(along = object)) {
    if ((length(nams[[i]]) == 1) && (nams[[i]] == "(Intercept)")) {
      val <- c(val, pnames[i])
    } else {
      val <- c(val, paste(pnames[i], nams[[i]], sep = "."))
    }
  }
  val
}

needUpdate.default <-
  function(object)
{
  val <- attr(object, "needUpdate")
  if (is.null(val) || !val) FALSE
  else TRUE
}

pairs.compareFits <-
  function(object, subset, key = TRUE, ...)
{

  if(!missing(subset)) {
    object <- object[subset,,]
  }
  dims <- dim(object)
  if(dims[3] == 1) {
    stop("At least two coefficients are needed.")
  }
  dn <- dimnames(object)
  coefs <- array(c(object), c(dims[1]*dims[2], dims[3]),
		 list(rep(dn[[1]], dims[2]), dn[[3]]))
  if(dims[3] > 2) {			# splom
    tt <- list(coefs = coefs,
	       grp = ordered(rep(dn[[2]], rep(dims[1], dims[2])), 
		   levels  = dn[[2]]))
    args <- list(formula = ~ coefs,
		  data = tt,
		  groups = tt$grp,
		  panel = function(x, y, subscripts, groups, ...)
		  {
		    panel.superpose(x, y, subscripts, groups)
		    aux <- groups[subscripts]
		    aux <- aux == unique(aux)[1]
		    segments(x[aux], y[aux], x[!aux], y[!aux], 
			     lty = 2, lwd = 0.5)
		  })
  } else {
    tt <- list(x = coefs[,1], y = coefs[,2],
	       grp = ordered(rep(dn[[2]], rep(dims[1], dims[2])),
		   levels = dn[[2]]))
    args <- list(formula = y ~ x,
		  data = tt,
		  groups = tt$grp,
		  panel = function(x, y, subscripts, groups, ...)
		  {
		    panel.grid()
		    panel.superpose(x, y, subscripts, groups)
		    aux <- groups[subscripts]
		    aux <- aux == unique(aux)[1]
		    segments(x[aux], y[aux], x[!aux], y[!aux], 
			     lty = 2, lwd = 0.5)
		  }, xlab = dn[[3]][1], ylab = dn[[3]][2])
  }
  dots <- list(...)
  args[names(dots)] <- dots
  if(is.logical(key)) {
    if(key && length(unique(tt$grp)) > 1) {
      args[["key"]] <- 
	list(points = Rows(trellis.par.get("superpose.symbol"), 1:2),
	     text = list(levels = levels(tt$grp)), columns = 2)
    }
  } else {
    args[["key"]] <- key
  }
  if(dims[3] > 2) do.call("splom", args) else do.call("xyplot", args)
}

##*## Plot method for ACF objects
plot.ACF <-
  function(object, alpha = 0, xlab = "Lag", ylab = "Autocorrelation",
           grid = FALSE, ...)
{
  ylim <- range(object$ACF)
  if (alpha) {
    assign("stdv",  qnorm(1-alpha/2)/sqrt(attr(object,"n.used")),
           frame = 1)
    stMax <- max(stdv)
    ylim <- c(min(c(-stMax, ylim[1])), max(c(ylim[2], stMax)))
  }
  assign("alpha", as.logical(alpha))
  assign("grid", grid)
  xyplot(ACF ~ lag, object, ylim = ylim,
         panel = function(x, y) {
           if (grid) panel.grid()
           panel.xyplot(x, y, type = "h")
           panel.abline(0, 0)
           if (alpha) {
             lines(x, stdv, lty = 2)
             lines(x, -stdv, lty = 2)
           }
         }, xlab = xlab, ylab = ylab, ...)
}

plot.augPred <-
  function(x, key = TRUE, grid = FALSE, ...)
{
  labels <- list(xlab = paste(attr(x, "labels")$x, attr(x, "units")$x),
		 ylab = paste(attr(x, "labels")$y, attr(x, "units")$y))
  labels <- labels[unlist(lapply(labels, function(el) length(el) > 0))]
  args <- c(list(formula = attr(x, "formula"),
		 groups = as.name(".type"),
		 data = x,
		 strip = function(...) strip.default(..., style = 1),
		 panel = if (length(levels(x[[".type"]])) == 2) {
                   ## single prediction level
                   function(x, y, subscripts, groups, ...) {
                     if (grid) panel.grid()
                     orig <- groups[subscripts] == "original"
                     panel.xyplot(x[orig], y[orig], ...)
                     panel.xyplot(x[!orig], y[!orig], ..., type = "l")
                   }
                 } else {             # multiple prediction levels
                   function(x, y, subscripts, groups, ...) {
                     if (grid) panel.grid()
                     orig <- groups[subscripts] == "original"
                     panel.xyplot(x[orig], y[orig], ...)
                     panel.superpose(x[!orig], y[!orig], subscripts[!orig],
                                     groups, ..., type = "l")
                   }
                 }), labels)
  ## perhaps include key argument allowing logical values
  dots <- list(...)
  args[names(dots)] <- dots
  if (is.logical(key) && key) {
    levs <- levels(x[[".type"]])
    if ((lLev <- length(levs)) > 2) {	# more than one levels
      lLev <- lLev - 1
      levs <- levs[1:lLev]
      aux <- !is.na(match(substring(levs, 1, 8), "predict."))
      if (sum(aux) > 0) {
	levs[aux] <- substring(levs[aux], 9)
      }
      args[["key"]] <- 
	list(lines = c(Rows(trellis.par.get("superpose.line"), 1:lLev),
		       list(size = rep(3, lLev))),
	     text = list(levels = levs), columns = lLev)
    } 
  } else {
    args[["key"]] <- key
  }
  assign("grid", grid)
  do.call("xyplot", args)
}

plot.compareFits <-
  function(object, subset, key = TRUE, mark = NULL, ...)
{

  if(!missing(subset)) {
    object <- object[subset,,]
  }
  dims <- dim(object)
  dn <- dimnames(object)
  assign("mark", rep(mark, rep(dims[1] * dims[2], dims[3])))
  tt <- data.frame(group = ordered(rep(dn[[1]], dims[2] * dims[3]),
		       levels = dn[[1]]),
		   coefs = as.vector(object),
		   what = ordered(rep(dn[[3]],
		       rep(dims[1] * dims[2], dims[3])), levels = dn[[3]]),
		   grp = ordered(rep(rep(dn[[2]], rep(dims[1], dims[2])), 
		       dims[3]), levels = dn[[2]]))
  args <- list(formula = group ~ coefs | what,
	       data = tt,
	       scales = list(x=list(relation="free")),
	       strip = function(...) strip.default(..., style = 1),
	       xlab = "",
	       groups = tt$grp,
	       panel = function(x, y, subscripts, groups, ...)
	       {
		 dot.line <- trellis.par.get("dot.line")
		 panel.abline(h = y, lwd = dot.line$lwd, 
			      lty = dot.line$lty, col = dot.line$col)
		 if(!is.null(mark)) {
		   panel.abline(v = mark[subscripts][1], lty = 2)
		 }
		 panel.superpose(x, y, subscripts, groups)
	       })
  dots <- list(...)
  args[names(dots)] <- dots
  if(is.logical(key)) {
    if(key && length(unique(tt$grp)) > 1) {
      args[["key"]] <- 
	list(points = Rows(trellis.par.get("superpose.symbol"), 1:2),
	     text = list(levels = levels(tt$grp)), columns = 2)
    }
  } else {
    args[["key"]] <- key
  }
  do.call("dotplot", args)
}

plot.Variogram <-
  function(object, smooth, showModel, sigma = NULL, span = 0.6,
           xlab = "Distance", ylab = "Semivariogram", type = "p", ylim,
           grid = FALSE, ...)
{
  trlLin <- trellis.par.get("superpose.line")
  coll <- attr(object, "collapse")
  modVrg <- attr(object, "modelVariog")
  lineT <- 1
  if (!is.na(match(type, c("l","o","b")))) {
    lineT <- lineT + 1
  }
  if (missing(showModel)) {
    showModel <- !is.null(modVrg)
  }
  if (showModel) {
    if (is.null(modVrg)) {
      stop("No model variogram available, with showModel = TRUE")
    }
    assign("ltyM", trlLin$lty[lineT])
    assign("colM", trlLin$col[lineT])
    assign("modVrg", modVrg)
    lineT <- lineT + 1
  }
  if (missing(smooth)) {
    smooth <- !showModel 
  }
  if (smooth) {
    assign("ltyS", trlLin$lty[lineT])
    assign("colS", trlLin$col[lineT])
  }
  assign("smooth", smooth)
  assign("showModel", showModel)
  assign("span", span)
  assign("type", type)
  assign("sigma", sigma)
  assign("grid", grid)
  if (missing(ylim)) {
    ylim <- c(0, max(object$variog))
  }
  xyplot(variog ~ dist, object, ylim = ylim, 
         panel = function(x, y, ...) {
           if (grid) panel.grid()
           panel.xyplot(x, y, type = type, ...)
           if (showModel) {
             panel.xyplot(modVrg$dist, modVrg$variog, type = "l",
                          col = colM, lty = ltyM, ...)
           }
           if (smooth) {
             panel.loess(x, y, span = span, col = colS, lty = ltyS, ...)
           }
           if (!is.null(sigma)) {
             panel.abline(c(sigma, 0), lty = 2)
           }
         }, xlab = xlab, ylab = ylab, ...)
}

print.compareFits <-
  function(x, ...)
{			# Will need to be changed for S4!
  print(unclass(x), ...)
}

print.correlation <-
  ## Print only the lower triangle of a correlation matrix
  function(x, title = " Correlation:", rdig = 3, ...)
{
  p <- dim(x)[2]
  if (p > 1) {
    cat(title, "\n")
    ll <- lower.tri(x)
    x[ll] <- format(round(x[ll], digits = rdig))
    x[!ll] <- ""
    if (!is.null(colnames(x))) {
      colnames(x) <- abbreviate(colnames(x), minlength = rdig + 3)
    }
   print(x[-1,  - p, drop = FALSE], ..., quote = FALSE)
  }
  invisible(x)
}

print.logLik <- 
  function(x, ...) print(c(x), ...)

pruneLevels.factor <-
  function(object)
{
  levs <- levels(object)
  factor(as.character(object),
         levels = levs[!is.na(match(levs, as.character(object)))])
}

pruneLevels.ordered <-
  function(object)
{
  levs <- levels(object)
  ordered(as.character(object),
          levels = levs[!is.na(match(levs, as.character(object)))])
}  

residuals.nls <- 
  function(object, type = c("response", "pearson"), ...)
{
  type <- match.arg(type)
  val <- as.vector(object$m$resid())
  if (type == "pearson") {
    std <- sqrt(sum(val^2)/(length(val) - length(coef(object))))
    val <- val/std
    attr(val, "label") <- "Standardized residuals"
  } else {
    lab <- "Residuals"
    if (!is.null(aux <- attr(object, "units")$y)) {
      lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
  }
  val
}

Variogram.default <-
  function(object, distance)
{
  ld <- length(distance)
  lo <- length(object)
  if (ld != round(lo*(lo-1)/2)) {
    stop("Distance and object have incompatible lengths")
  }
  val <- outer(object, object, function(x,y) ((x - y)^2)/2)
  val <- val[lower.tri(val)]
  val <- data.frame(variog = val, dist = distance)
  class(val) <- c("Variogram", "data.frame")
  val
}

## Local Variables:
## mode:S
## End:
