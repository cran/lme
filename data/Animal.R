### $Id: Animal.q,v 1.1 1998/04/10 23:02:00 bates Exp $
### Animal breeding experiment with two dams each mated to five randomly
### selected sires.  Data set 6.4 from "SAS System for Mixed Models".
"Animal" <-
  structure(list(
  Sire = structure(factor(c(1, 1, 1, 1, 2, 2, 2, 
    2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5), levels=1:5),
    class = "factor", .Label = c("1", "2", "3", "4", "5")),
Dam = structure(factor(c(1, 1, 2, 2, 1, 
  1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2), levels=1:2),
  class = "factor", .Label = c("1", "2")),
AvgDailyGain = c(2.24, 1.85, 2.05, 2.41, 1.99, 1.93, 2.72, 
  2.32, 2.33, 2.68, 2.69, 2.71, 2.42, 2.01, 1.86, 1.79, 2.82, 2.64, 
  2.58, 2.56)),
class = c("nmGroupedData", "groupedData", "data.frame"),
row.names = 1:20,
formula = AvgDailyGain ~ 1 | Sire/Dam,
formulaList = list(Sire = ~Sire, Dam = ~Dam),
labels = list(y = "Average Daily Weight Gain"),
order.groups = list(Sire = TRUE, Dam = TRUE),
FUN = function (x) mean(x, na.rm = TRUE))
