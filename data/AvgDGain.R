### $Id: AvgDGain.q,v 1.1 1998/05/19 19:45:40 bates Exp $
### Average daily weight gain of steers fed for 160 days.  The Treatment
### is the level of medicated feed additive in the diet.  It is given as both
### an ordered factor (Treatment) and a continuous covariate (Trt).  The initial
### weight is also recorded.  The steers were housed in 8 different barns,
### the Block factor.  This is data set 5.3 from "SAS System for Mixed Models".
"AvgDailyGain" <-
  structure(list(
  Id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 
    29, 30, 31, 32),
Block = structure(c(1, 1, 1, 1, 2, 2, 2, 2, 
  3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 
  8, 8, 8), .Label = c("1", "2", "3", "4", "5", "6", "7", "8"),
  class = c("ordered", "factor")),
Treatment = structure(c(1, 2, 3, 4, 1, 2, 3, 4, 1, 
  2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 
  3, 4), .Label = c("0", "10", "20", "30"),
  class = c("ordered", "factor")),
adg = c(1.03, 1.54, 1.82, 1.86, 1.31, 2.16, 2.13, 
  2.23, 1.59, 2.53, 2.33, 1.8, 2.09, 2.2, 2.21, 2.82, 1.66, 2.3, 
  2.65, 2.18, 1.42, 1.93, 1.58, 1.49, 1.41, 1.65, 1.08, 1.34, 0.18, 
  0.64, 0.76, 0.7),
InitWt = c(338, 477, 444, 370, 403, 451, 450, 
  393, 394, 499, 482, 317, 499, 411, 391, 396, 371, 418, 486, 333, 
  395, 325, 316, 311, 414, 313, 309, 323, 315, 376, 308, 439),
Trt = c(0, 
  10, 20, 30, 0, 10, 20, 30, 0, 10, 20, 30, 0, 10, 20, 30, 0, 10, 
  20, 30, 0, 10, 20, 30, 0, 10, 20, 30, 0, 10, 20, 30)),
row.names = c("1", 
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
"14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", 
"25", "26", "27", "28", "29", "30", "31", "32"),
class = c("nfnGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = adg ~ Trt | Block,
labels = list(x = "Level of medicated feed additive in diet",
y = "Average Daily Gain of steers fed for 160 days"),
FUN = function (x) max(x, na.rm = TRUE), order.groups = TRUE)
