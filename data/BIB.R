### $Id: BIB.q,v 1.1 1998/04/10 23:02:00 bates Exp $
### A balanced incomplete block design for analysis of covariance.  Data set 5.4
### from "SAS System for Mixed Models".  Perhaps artificial data.  The Treatment
### is referred to as a "diet model" but no other information is given.
"BIB" <-
  structure(list(
  Block = structure(ordered(c(1, 1, 1, 2, 2, 2, 
    3, 3, 3, 6, 6, 6, 5, 5, 5, 7, 7, 7, 8, 8, 8, 4, 4, 4), levels=1:8),
    class = c("ordered", "factor"),
    .Label = c("1", "2", "3", "8", "5", "4", "6", "7")),
Treatment = structure(factor(c(1, 
  2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4, 1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 
  3, 4), levels=1:4), class = "factor",
  .Label = c("1", "2", "3", "4")),
y = c(31, 29, 31, 29, 34, 33, 31, 28, 34, 39, 35, 32, 
  33, 35, 38, 35, 31, 42, 42, 43, 42, 27, 37, 29),
x = c(20, 18, 11, 37, 37, 39, 29, 12, 31, 37, 29, 28, 12, 19, 16, 31,
  13, 39, 38, 30, 25, 13, 39, 21),
Grp = structure(factor(c(1, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 1,
  1, 2, 2, 1, 1, 2, 2, 1, 2), levels=1:2), class = "factor",
  .Label = c("13", "24"))),
class = c("nfnGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = 1:24,
formula = y ~ x | Block,
FUN = function (x) max(x, na.rm = TRUE),
order.groups = TRUE)
