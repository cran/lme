### $Id: DailyGain.q,v 1.1 1998/04/10 23:02:01 bates Exp $
### Average daily weight gain in steers fed for 160 days.  Treatment indicates
### the level of the medicated feed additive added to the base ration.  Data
### set 5.3 from "SAS System for Mixed Models"
"DailyGain" <-
  structure(list(
  Block = structure(ordered(c(3, 3, 3, 3, 5, 5, 
    5, 5, 6, 6, 6, 6, 8, 8, 8, 8, 7, 7, 7, 7, 4, 4, 4, 4, 2, 2, 2, 
    2, 1, 1, 1, 1), levels=1:8), class = c("ordered", "factor"),
    .Label = c("8", "7", "1", "6", "2", "3", "5", "4")),
Treatment = structure(ordered(c(1, 
  2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 
  3, 4, 1, 2, 3, 4, 1, 2, 3, 4), levels=1:4),
  class = c("ordered", "factor"),
  .Label = c("0", "10", "20", "30")),
AvgDailyGain = c(1.03, 
  1.54, 1.82, 1.86, 1.31, 2.16, 2.13, 2.23, 1.59, 2.53, 2.33, 1.8, 
  2.09, 2.2, 2.21, 2.82, 1.66, 2.3, 2.65, 2.18, 1.42, 1.93, 1.58, 
  1.49, 1.41, 1.65, 1.08, 1.34, 0.18, 0.64, 0.76, 0.7),
InitialWt = c(338, 
  477, 444, 370, 403, 451, 450, 393, 394, 499, 482, 317, 499, 411, 
  391, 396, 371, 418, 486, 333, 395, 325, 316, 311, 414, 313, 309, 
  323, 315, 376, 308, 439)),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = 1:32,
formula = AvgDailyGain ~ Treatment | Block,
labels = list(y = "Average Daily Weight Gain"),
FUN = function (x) max(x, na.rm = TRUE),
order.groups = TRUE)
