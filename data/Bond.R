# $Id: Bond.q,v 1.1 1998/04/10 23:02:00 bates Exp $
# Pressure to break a bond of two pieces of metal that used one of the metals
# as the bonding agent.  Example 1.2.4 in "SAS System for Mixed Models".
"Bond" <-
  structure(list(
  pressure = c(67, 71.9, 72.2, 67.5, 68.8, 66.4, 
    76, 82.6, 74.5, 72.7, 78.1, 67.3, 73.1, 74.2, 73.2, 65.8, 70.8, 
    68.7, 75.6, 84.9, 69),
Metal = structure(factor(c(3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3,
  2, 1, 3, 2, 1, 3, 2, 1), levels=1:3), .Label = c("c", "i", "n"),
  class = "factor"),
Ingot = structure(ordered(c(3, 3, 3, 1, 1, 1, 6, 6, 6, 5, 5, 5, 4,
  4, 4, 2, 2, 2, 7, 7, 7), levels=1:7), class = c("ordered", "factor"),
  .Label = c("2", "6", "1", "5", "4", "3", "7"))),
row.names = 1:21,
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
formula = pressure ~ 1 | Ingot,
inner = ~Metal,
labels = list(y = "Pressure required to break bond"),
FUN = function (x) mean(x, na.rm = TRUE),
order.groups = TRUE)
