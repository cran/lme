### $Id: zzModels.q,v 1.14 1999/06/04 16:06:38 pinheiro Exp $
###
###       Individual selfStarting nonlinear regression models
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

##*## SSasymp - asymptotic regression model

SSasymp <- selfStart(~ Asym + (R0 - Asym) * exp(-exp(lrc) * input),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if (nrow(xy) < 3) {
    stop("Too few distinct input values to fit a asymptotic regression model")
  }
  if(nrow(xy) > 3) {
    xy$ydiff <- abs(xy$y - NLSstRtAsymptote(xy))
    lrc <- log( - coef(lm(log(ydiff) ~ x, data = xy))[2])
    names(lrc) <- NULL	
    ## This gives an estimate of the log (rate constant).  Use that
    ## with a partially linear nls algorithm
    pars <- coef(nls(y ~ cbind(1 - exp( - exp(lrc) * x),
				  exp(- exp(lrc) * x)),
		     data = xy,
		     start = list(lrc = lrc),
		     algorithm = "plinear"))
  }
  else {
    ydiff <- diff(xy$y)
    if(prod(ydiff) <= 0) {
      stop("Can't fit an asymptotic regression model to these data")
    }
    avg.resp <- xy$y
    frac <- (avg.resp[3] - avg.resp[1])/(avg.resp[2] - avg.resp[1])
    xunique <- unique(xy$x)
    xdiff <- diff(xunique)
    if(xdiff[1] == xdiff[2]) {	# equal spacing - can use a shortcut
      expmRd <- frac - 1
      rc <-  - log(expmRd)/xdiff[1]
      lrc <- log(rc)
      expmRx1 <- exp( - rc * xunique[1])
      bma <- ydiff[1]/(expmRx1 * (expmRd - 1))
      Asym <- avg.resp[1] - bma * expmRx1
      pars <- c(lrc = lrc, Asym = Asym, R0 = bma + Asym)
    }
    else {
      stop("Too few observations to fit an asymptotic regression model")
    }
  }
  names(pars) <- NULL
  val <- list(pars[2], pars[3], pars[1])
  names(val) <- mCall[c("Asym", "R0", "lrc")]
  val
},
c("Asym", "R0", "lrc"))

##*## SSasympOff - alternate formulation of asymptotic regression model
##*## with an offset

SSasympOff <- selfStart(~ Asym *( 1 - exp(-exp(lrc) * (input - c0) ) ),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if (nrow(xy) < 4) {
    stop("Too few distinct input values to fit the asympOff model")
  }
  xy$ydiff <- abs(xy$y - NLSstRtAsymptote(xy))
  lrc <- log( - coef(lm(log(ydiff) ~ x, data = xy))[2]) # log( rate constant)
  pars <- as.vector(coef(nls(y ~ cbind(1, exp(- exp(lrc) * x)),
			     data = xy, algorithm = "plinear",
			     start = list(lrc = lrc))))
  val <- list(pars[2], pars[1], exp(-pars[1]) * log(-pars[3]/pars[2]))
  names(val) <- mCall[c("Asym", "lrc", "c0")]
  val
},
c("Asym", "lrc", "c0"))

##*## SSasympOrig - exponential curve through the origin to an asymptote 

SSasympOrig <- selfStart(~ Asym * (1 - exp(-exp(lrc) * input)),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if (nrow(xy) < 3) {
    stop("Too few distinct input values to fit the asympOrig model")
  }
  ## get a preliminary estimate for A
  A0 <- NLSstRtAsymptote(xy)
  ## get a least squares estimate for log of the rate constant
  lrc <- log(abs(mean(log(1 - xy$y/A0)/xy$x)))
  ## use the partially linear form to converge quickly
  pars <- as.vector(coef(nls(y ~ 1 - exp(-exp(lrc)*x),
			     data = xy,
			     start = list(lrc = lrc),
			     algorithm = "plinear")))
  value <- c(pars[2], pars[1])
  names(value) <- mCall[c("Asym", "lrc")]
  value
},
c("Asym", "lrc"))

##*## SSbiexp - linear combination of two exponentials

SSbiexp <-
  selfStart(~ A1 * exp(-exp(lrc1)*input) + A2 * exp(-exp(lrc2) * input),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if (nrow(xy) < 5) {
    stop("Too few distinct input values to fit a biexponential")
  }
  ndistinct <- nrow(xy)
  nlast <- max(3, ndistinct/2)		# take at least half the data
  dlast <- xy[(ndistinct + 1 - nlast):ndistinct, ]
  pars2 <- coef(lm(log(y) ~ x, data = dlast))
  lrc2 <- log(abs(pars2[2]))		# log of the slope
  xy[["res"]] <- xy[["y"]] - exp(pars2[1]) * exp(-exp(lrc2)*xy[["x"]])
  dfirst <- xy[1:(ndistinct - nlast), ]
  pars1 <- coef(lm(log(res) ~ x, data = dfirst))
  lrc1 <- log(abs(pars1[2]))
  pars <- coef(nls(y ~ cbind(exp(-exp(lrc1)*x), exp(-exp(lrc2)*x)),
		   data = xy,
		   start = list(lrc1 = lrc1, lrc2 = lrc2),
		   algorithm = "plinear"))
  value <- c(pars[3], pars[1], pars[4], pars[2])
  names(value) <- mCall[c("A1", "lrc1", "A2", "lrc2")]
  value
},
c("A1", "lrc1", "A2", "lrc2"))

##*## SSfol - first order compartment model with the log of the rates
##*##       and the clearence  
SSfol <- 
  selfStart(~Dose * exp(lKe + lKa - lCl) * (exp(-exp(lKe) * input) -
		    exp(-exp(lKa) * input))/(exp(lKa) - exp(lKe)),
function(mCall, data, LHS)
{
  resp <- eval(LHS, data)
  input <- eval(mCall[["input"]], data)
  Dose <- eval(mCall[["Dose"]], data)
  n <- length(resp)
  if(length(input) != n) {
    stop(paste("must have length of response =", 
	       "length of second argument to SSfol"))
  }
  if(n < 4) {
    stop("must have at least 4 observations to fit an SSfol model")
  }
  rmaxind <- order(resp)[n]

  lresp <- log(resp)
  if(rmaxind == n) {
    lKe <- -2.5
  } else {
    lKe <- log((lresp[rmaxind] - lresp[n])/(input[n] - input[rmaxind]))
  }
  cond.lin <- nls(resp ~ (exp(-input * exp(lKe))-exp(-input * exp(lKa))) * Dose,
                  data = list(resp = resp, input = input, Dose = Dose, lKe = lKe),
                  start = list(lKa = lKe + 1), 
		  algorithm = "plinear")
  pars <- clearNames(coef(cond.lin))
  cond.lin <- nls(resp ~ (Dose * (exp(-input*exp(lKe))-
		     exp(-input*exp(lKa))))/(exp(lKa) - exp(lKe)), 
		  data = data.frame(list(resp = resp, input = input, 
		      Dose = Dose)), 
		  start = list(lKa = pars[1],lKe = lKe), 
		  algorithm = "plinear")
  pars <- clearNames(coef(cond.lin))
  lKa <- pars[1]
  lKe <- pars[2]
  Ka <- exp(lKa)
  Ke <- exp(lKe)
  value <- list(lKe, lKa, log((Ke * Ka)/(pars[3])))
  names(value) <- as.character(mCall)[4:6]
  value
}, c("lKe", "lKa", "lCl"))

##*## SSfpl - four parameter logistic model

SSfpl <- selfStart(~ A + (B - A)/(1 + exp((xmid - input)/scal)),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if (nrow(xy) < 5) {
    stop("Too few distinct input values to fit a four-parameter logistic")
  }
  ## get a preliminary estimate for xmid
  xydata <- c(as.list(xy), xmid = NLSstClosestX(xy, mean(range(xy[["y"]]))))
  pars <- as.vector(coef(nls(y ~ cbind(1, 1/(1 + exp((xmid - x)/exp(lscal)))), 
			     data = xydata,
			     start = list(lscal = 0),  # scal = 1
			     algorithm = "plinear")))
  param(xy, "lscal") <- pars[[1]]
  pars <- as.vector(coef(nls(y ~ cbind(1, 1/(1 + exp((xmid - x)/exp(lscal)))), 
			     data = xy,
                             start = list(xmid = xydata$xmid, lscal = pars[1]),
			     algorithm = "plinear")))
  value <- c(pars[3], pars[3] + pars[4], pars[1], exp(pars[2]))
  names(value) <- mCall[c("A", "B", "xmid", "scal")]
  value
},
c("A", "B", "xmid", "scal"))

##*## SSlogis - logistic model for nonlinear regression

SSlogis <- selfStart(~ Asym/(1 + exp((xmid - input)/scal)),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if(nrow(xy) < 4) {
    stop("Too few distinct input values to fit a logistic")
  }
  z <- xy[["y"]]
  if (min(z) <= 0) { z <- z - 1.05 * min(z) } # avoid zeroes
  z <- z/(1.05 * max(z))		# scale to within unit height
  xy[["z"]] <- log(z/(1 - z))		# logit transformation
  aux <- coef(lm(x ~ z, xy))
  pars <- as.vector(coef(nls(y ~ 1/(1 + exp((xmid - x)/scal)),
                             data = xy,
                             start = list(xmid = aux[1], scal = aux[2]),
                             algorithm = "plinear")))
  value <- c(pars[3], pars[1], pars[2])
  names(value) <- mCall[c("Asym", "xmid", "scal")]
  value
},
c("Asym", "xmid", "scal"))

##*## SSmicmen - Michaelis-Menten model for enzyme kinetics.

SSmicmen <- selfStart(~ Vm * input/(K + input),
function(mCall, data, LHS)
{
  xy <- sortedXyData(mCall[["input"]], LHS, data)
  if (nrow(xy) < 3) {
    stop("Too few distinct input values to fit a Michaelis-Menten")
  }
  ## take the inverse transformation
  pars <- as.vector(coef(lm(1/y ~ I(1/x), data = xy)))
  ## use the partially linear form to converge quickly
  pars <- as.vector(coef(nls(y ~ x/(K + x),
			     data = xy,
			     start = list(K = abs(pars[2]/pars[1])),
			     algorithm = "plinear")))
  value <- c(pars[2], pars[1])
  names(value) <- mCall[c("Vm", "K")]
  value
},
c("Vm", "K"))

### Local variables:
### mode: S
### End:
