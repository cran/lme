### $Id: Mississippi.q,v 1.1 1998/04/10 23:02:01 bates Exp $
### Nitrogen concentrations at several sampled sites of influents to the
### Mississippi River.  Data set 4.2 in "SAS System for Mixed Models"
"Mississippi" <-
  structure(list(
  influent = structure(ordered(c(4, 4, 4, 4, 4, 
    4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 
    5, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6), levels=1:6), 
    class = c("ordered", "factor"),
    .Label = c("3", "5", "2", "1", "4", "6")), 
y = c(21, 27, 29, 17, 19, 12, 
  29, 20, 20, 21, 11, 18, 9, 13, 23, 2, 20, 19, 20, 11, 14, 14, 
  24, 30, 21, 31, 27, 7, 15, 18, 4, 28, 41, 42, 35, 34, 30),
Type = structure(factor(c(2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 
  2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3), levels=1:3),
  class = "factor", .Label = c("1", "2", "3"))),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = 1:37,
outer = ~Type,
formula = y ~ 1 | influent,
labels = list(y = "Nitrogen concentration in Mississippi River"),
units = list(y = "(ppm)"),
FUN = function (x) max(x, na.rm = TRUE),
order.groups = TRUE)
