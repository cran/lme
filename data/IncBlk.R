### $Id: IncBlk.q,v 1.1 1998/04/10 23:02:01 bates Exp $
### An unbalanced incomplete block design.  Data set 5.5 from "SAS System
### for Mixed Models".  Probably artificial data.  The treatment is
### described as a "diet model".
"IncBlk" <-
  structure(list(
  Block = structure(ordered(c(12, 12, 9, 9, 10, 
    10, 4, 4, 8, 8, 11, 11, 6, 6, 1, 1, 3, 3, 5, 5, 7, 7, 2, 2),
    levels=1:12),
    class = c("ordered", "factor"),
    .Label = c("8", "12", "9", "4", "10", "7", "11", "5", 
      "2", "3", "6", "1")),
Treatment = structure(factor(c(1, 2, 1, 
  2, 1, 2, 1, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4),
  levels=1:4),
  class = "factor",
  .Label = c("1", "2", "3", "4")),
y = c(0.62, 0.91, 0.41, 0.48, 0.41, 0.49, 0.26, 0.28, 0.29, 
  0.37, 0.73, 0.72, 0.33, 0.31, 0.18, 0.18, 0.19, 0.25, 0.28, 0.32, 
  0.33, 0.27, 0.24, 0.23),
x = c(0.078, 0.01, 0.032, 0.05, 0, 0.015, 
  0.01, 0.016, 0.053, 0.069, 0.007, 0.062, 0.036, 0.068, 0.068, 
  0.057, 0.077, 0.09, 0.023, 0.039, 0.017, 0.062, 0.058, 0.082)),
class = c("nfnGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = 1:24,
formula = y ~ x | Block, inner = ~Treatment,
FUN = function (x) max(x, na.rm = TRUE),
order.groups = TRUE)
