# $Id: PBIB.q,v 1.1 1998/04/10 23:02:01 bates Exp $
# A partially balanced incomplete blocked design.
# Data set 1.5.1 from "SAS System for Mixed Models".
"PBIB" <-
  structure(list(
response = c(2.4, 2.5, 2.6, 2, 2.7, 2.8, 2.4, 
    2.7, 2.6, 2.8, 2.4, 2.4, 3.4, 3.1, 2.1, 2.3, 4.1, 3.3, 3.3, 2.9, 
    3.4, 3.2, 2.8, 3, 3.2, 2.5, 2.4, 2.6, 2.3, 2.3, 2.4, 2.7, 2.8, 
    2.8, 2.6, 2.5, 2.5, 2.7, 2.8, 2.6, 2.6, 2.6, 2.3, 2.4, 2.7, 2.7, 
    2.5, 2.6, 3, 3.6, 3.2, 3.2, 3, 2.8, 2.4, 2.5, 2.4, 2.5, 3.2, 
    3.1),
Treatment = structure(factor(c(15, 9, 1, 13, 5, 7, 8, 1, 
  10, 1, 14, 2, 15, 11, 2, 3, 6, 15, 4, 7, 12, 4, 3, 1, 12, 14, 
  15, 8, 6, 3, 14, 5, 5, 4, 2, 13, 10, 12, 13, 6, 9, 7, 10, 3, 
  8, 6, 2, 9, 5, 9, 11, 12, 7, 13, 14, 11, 10, 4, 8, 11),
  levels=1:15), class = "factor",
  .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
    "12", "13", "14", "15")),
Block = structure(ordered(c(1, 1, 1, 1, 5, 5, 5, 
  5, 6, 6, 6, 6, 12, 12, 12, 12, 15, 15, 15, 15, 13, 13, 13, 13, 
  10, 10, 10, 10, 3, 3, 3, 3, 7, 7, 7, 7, 8, 8, 8, 8, 2, 2, 2, 
  2, 4, 4, 4, 4, 14, 14, 14, 14, 9, 9, 9, 9, 11, 11, 11, 11),
  levels=1:15),
  class = c("ordered", "factor"),
  .Label = c("1", "11", "8", "12", "2", "3", "9", "10", "14", "7",
    "15", "4", "6", "13", "5"))),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = response ~ Treatment | Block,
row.names = 1:60,
FUN = function (x) max(x, na.rm = TRUE), order.groups = TRUE)
