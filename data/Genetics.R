### $Id: Genetics.q,v 1.1 1998/04/10 23:02:01 bates Exp $
### Randomized complete blocked design of wheat families assigned to blocks
### within locations.  Data set 4.5 from "SAS System for Mixed Models"
"Genetics" <-
  structure(list(
  Location = structure(factor(c(1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), levels=1:4), class = "factor",
    .Label = c("1", "2", "3", "4")),
Block = structure(factor(c(1, 2, 3, 1, 2, 3, 
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), levels=1:3), class = "factor",
  .Label = c("1", "2", "3")),
Family = structure(factor(c(1, 1, 1, 2, 2, 2, 3, 
  3, 3, 4, 4, 4, 5, 5, 5, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 
  5, 5, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1, 1, 1, 2, 
  2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5), levels=1:5), class = "factor",
  .Label = c("1", "2", "3", "4", "5")),
Yield = c(268, 279, 261, 242, 261, 258, 
  242, 245, 234, 225, 231, 219, 236, 260, 248, 238, 220, 243, 215, 
  192, 226, 198, 151, 191, 195, 182, 202, 201, 161, 196, 221, 216, 
  224, 208, 197, 201, 186, 173, 161, 207, 183, 186, 200, 207, 190, 
  194, 194, 197, 203, 191, 204, 177, 170, 180, 180, 195, 193, 199, 
  183, 208)),
class = c("nmGroupedData", "groupedData", "data.frame"),
row.names = 1:60,
formula = Yield ~ 1 | Location/Block,
formulaList = list(Location = ~Location, Block = ~Block),
labels = list(y = "Wheat yield"),
order.groups = list(Location = TRUE, Block = TRUE),
FUN = function (x) mean(x, na.rm = TRUE))

