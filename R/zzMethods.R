### $Id: zzMethods.q,v 1.9 1999/06/04 02:52:17 pinheiro Exp $
###
###   Miscellaneous methods that must be defined last in the library
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

AIC.lme <- AIC.lmList <- AIC.gls <- AIC.lm
BIC.lme <- BIC.lmList <- BIC.gls <- BIC.lm

comparePred.lme <- comparePred.lmList <- comparePred.gls

getData.nls <- getData.lme <- getData.gls

getData.nlme <- getData.gnls

qqnorm.nls <- qqnorm.lm <- qqnorm.gls 

plot.nls <- plot.lme
