### $Id: Demand.q,v 1.2 1998/07/02 21:48:44 bates Exp $
### Doubly repeated measures or "panel data" on per capita demand deposits.
### Example 3.6 from "SAS System for Mixed Models"
"Demand" <-
  structure(list(
State = structure(ordered(c(3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 
  7, 7, 7, 7, 7, 7, 7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1), levels=1:7),
  class = c("ordered", "factor"),
  .Label = c("WA", "FL", "CA", "TX", "IL", "DC", "NY")),
Year = structure(ordered(c(1, 
  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
  11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 
  9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5, 
  6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), levels=1:11),
  class = c("ordered", "factor"),
  .Label = c("1949", "1950", "1951", "1952", "1953", 
    "1954", "1955", "1956", "1957", "1958", "1959")),
d = c(533, 603, 669, 651, 609, 634, 665, 676, 642, 678, 714, 854, 1013, 
  1185, 1076, 1004, 1044, 1067, 1062, 1120, 1196, 1168, 408, 433, 
  469, 470, 464, 465, 545, 567, 531, 533, 522, 843, 860, 887, 914, 
  909, 928, 939, 944, 899, 919, 874, 1370, 1405, 1409, 1421, 1395, 
  1415, 1431, 1416, 1443, 1453, 1417, 573, 634, 679, 668, 666, 
  708, 722, 708, 675, 716, 703, 418, 501, 525, 519, 500, 537, 545, 
  525, 494, 521, 515),
y = c(1347, 1464, 1608, 1636, 1669, 1716, 
  1779, 1878, 1963, 2034, 2164, 1603, 1773, 2017, 1921, 1856, 1868, 
  1931, 1951, 2085, 2144, 2167, 1024, 1007, 1068, 1068, 1138, 1137, 
  1306, 1339, 1383, 1409, 1457, 1465, 1468, 1555, 1648, 1711, 1775, 
  1815, 1915, 1980, 2001, 2035, 1492, 1515, 1566, 1659, 1744, 1802, 
  1808, 1916, 2074, 2120, 2197, 995, 1052, 1154, 1176, 1228, 1285, 
  1335, 1358, 1416, 1457, 1520, 1146, 1324, 1433, 1481, 1531, 1602, 
  1649, 1656, 1711, 1754, 1809),
rd = c(0.343, 0.364, 0.367, 0.369, 
  0.41, 0.499, 0.496, 0.533, 0.63, 0.667, 0.664, 0.261, 0.267, 
  0.266, 0.267, 0.287, 0.308, 0.318, 0.322, 0.346, 0.36, 0.418, 
  0.354, 0.342, 0.335, 0.328, 0.354, 0.374, 0.378, 0.399, 0.447, 
  0.498, 0.523, 0.143, 0.146, 0.147, 0.144, 0.15, 0.164, 0.172, 
  0.183, 0.203, 0.214, 0.231, 0.112, 0.119, 0.119, 0.12, 0.134, 
  0.145, 0.146, 0.168, 0.189, 0.192, 0.203, 0.149, 0.147, 0.148, 
  0.147, 0.16, 0.182, 0.191, 0.208, 0.25, 0.278, 0.303, 0.358, 
  0.361, 0.365, 0.381, 0.414, 0.481, 0.529, 0.587, 0.681, 0.716, 
  0.73),
rt = c(1.114, 1.162, 1.493, 1.567, 1.594, 1.609, 1.637, 
  1.757, 2.641, 2.641, 2.648, 0.676, 0.662, 0.677, 0.729, 0.883, 
  1.5, 1.504, 1.598, 2.231, 2.1, 2.342, 0.909, 0.957, 1.002, 1.052, 
  1.118, 1.268, 1.339, 1.486, 2.42, 2.453, 2.489, 0.852, 0.847, 
  0.936, 1.059, 1.091, 1.13, 1.141, 1.354, 1.628, 1.737, 2.054, 
  0.687, 0.724, 0.795, 1.05, 1.241, 1.346, 1.406, 1.754, 2.231, 
  2.36, 2.521, 0.839, 0.836, 0.812, 1.07, 1.17, 1.328, 1.368, 1.544, 
  2.121, 2.241, 2.435, 0.937, 0.973, 1.039, 1.305, 1.342, 1.348, 
  1.77, 1.779, 2.313, 2.302, 2.495),
rs = c(2.905, 2.935, 3.093, 
  3.073, 3.357, 3.295, 3.451, 3.539, 3.93, 3.982, 4.047, 2.803, 
  2.877, 3.006, 2.975, 3.035, 3.083, 3.177, 3.25, 3.368, 3.457, 
  3.727, 2.314, 2.327, 2.428, 2.577, 2.625, 2.871, 2.882, 3.032, 
  3.338, 3.353, 3.575, 2.504, 2.448, 2.449, 2.568, 2.703, 2.748, 
  2.778, 2.932, 3.155, 3.402, 3.497, 2.099, 2.082, 2.218, 2.435, 
  2.477, 2.54, 2.655, 2.774, 2.957, 3.073, 3.223, 2.755, 2.74, 
  2.819, 2.88, 3.082, 3.093, 3.071, 3.068, 3.487, 3.413, 3.671, 
  2.068, 2.229, 2.367, 2.553, 2.848, 2.865, 2.907, 3.011, 3.252, 
  3.306, 3.507),
grp = structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1), .Label = "1", class = "factor")),
class = c("nffGroupedData", "nfGroupedData", "groupedData", "data.frame"),
row.names = as.character(1:77),
formula = d ~ Year | State,
labels = list(y = "per capita demand deposits"),
FUN = function (x) max(x, na.rm = TRUE),
order.groups = TRUE)
