### $Id: Cultivation.q,v 1.1 1998/04/10 23:02:01 bates Exp $
### Data Set 2.2(a) from SAS System for Mixed Models.  A split-plot design
### with Cult as the whole-plot treatment and Inoc as the sub-plot treatment
"Cultivation" <-
  structure(list(
  Block = structure(factor(c(1, 1, 1, 1, 1, 1, 2, 
    2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4), levels=1:4),
    class = "factor", .Label = c("1", "2", "3", "4")),
Cult = structure(factor(c(1, 1, 1, 2, 2, 2, 
  1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2), levels=1:2),
  .Label = c("a", "b"), class = "factor"),
Inoc = structure(factor(c(1, 2, 3, 1, 
  2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), levels=1:3),
  .Label = c("con", "dea", "liv"), class = "factor"),
drywt = c(27.4, 29.7, 34.5, 
  29.4, 32.5, 34.4, 28.9, 28.7, 33.4, 28.7, 32.4, 36.4, 28.6, 29.7, 
  32.9, 27.2, 29.1, 32.6, 26.7, 28.9, 31.8, 26.8, 28.6, 30.7)),
class = c("nmGroupedData", "groupedData", "data.frame"),
row.names = 1:24,
formula = drywt ~ 1 | Block/Cult,
formulaList = list(Block = ~Block, Cult = ~Cult),
labels = list(y = "Yield"),
inner = list(Cult = ~Inoc),
order.groups = list(Block = TRUE, Cult = TRUE),
FUN = function (x) max(x, na.rm = TRUE))
