### $Id: Wafer.q,v 1.3 1999/05/13 17:00:50 pinheiro Exp $
### Voltage versus intensity of current at 8 sites of 10 wafers
### from an experiment conducted at Lucent Technologies
"Wafer" <-
  structure(list
  (Wafer = structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
     2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
     4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
     5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 
     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 
     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
     10, 10, 10, 10),
     .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
     class = "factor"),
   Site = structure(.Data = c(1, 1, 1, 1, 1, 2, 
     2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 
     6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 1, 1, 1, 1, 
     1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 
     5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 1, 1, 
     1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 
     5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 
     1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 
     4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 
     8, 8, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 
     4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 
     8, 8, 8, 8, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 
     4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 
     7, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 
     3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 
     7, 7, 7, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 
     3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 
     7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
     2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 
     6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 2, 2, 
     2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 
     6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8),
     .Label = c("1", "2", "3", "4", "5", "6", "7", "8"),
     class = "factor"),
   voltage = c(0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2,
     1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8,
     1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4,
     0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2,
     2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6,
     2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2,
     1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8,
     1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4,
     0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2,
     2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6,
     2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2,
     1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8,
     1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4,
     0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2,
     2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6,
     2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2,
     1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8,
     1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4,
     0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2,
     2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6,
     2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2,
     1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8,
     1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4,
     0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2,
     2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6,
     2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2,
     1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8,
     1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4,
     0.8, 1.2, 1.6, 2, 2.4, 0.8, 1.2, 1.6, 2, 2.4),
   current = c(0.90088, 3.8682, 7.6406, 11.736, 15.934, 1.032, 4.1022,
     7.9316, 12.064, 16.294, 0.924, 3.908, 7.701, 11.818, 16.042,
     0.8816, 3.821, 7.5786, 11.668, 15.866, 0.93208, 3.909, 7.689,
     11.792, 16.002, 0.78162, 3.6044, 7.2764, 11.2924, 15.434, 0.8864,
     3.82, 7.5648, 11.64, 15.822, 0.91878, 3.8858, 7.6492, 11.734,
     15.924, 1.0216, 4.0374, 7.824, 11.896, 16.09, 1.249, 4.361,
     8.1752, 12.25, 16.434, 1.2824, 4.3766, 8.1534, 12.184, 16.328,
     0.95006, 3.903, 7.647, 11.688, 15.854, 1.2458, 4.3428, 8.1466,
     12.216, 16.396, 0.98324, 3.9484, 7.6942, 11.732, 15.892, 1.0834,
     4.093, 7.8468, 11.882, 16.038, 1.015, 4.0026, 7.7518, 11.788,
     15.946, 0.93266, 3.914, 7.6906, 11.786, 15.986, 0.98908, 3.9838,
     7.751, 11.83, 16.012, 1.0554, 4.1276, 7.9532, 12.082, 16.304,
     0.89568, 3.8454, 7.6058, 11.694, 15.886, 1.0648, 4.126, 7.9394,
     12.056, 16.272, 0.84492, 3.7396, 7.4564, 11.504, 15.664, 1.2228,
     4.3466, 8.172, 12.284, 16.482, 1.1904, 4.2582, 8.0294, 12.09,
     16.246, 0.96102, 3.9856, 7.7898, 11.904, 16.118, 1.0472, 4.1254,
     7.9546, 12.086, 16.31, 0.85842, 3.7592, 7.4678, 11.502, 15.644,
     0.84822, 3.7536, 7.488, 11.558, 15.742, 0.86412, 3.7524, 7.4538,
     11.486, 15.63, 0.8522, 3.7818, 7.5404, 11.632, 15.832, 0.9447,
     3.93, 7.7054, 11.8, 16, 0.96664, 3.9714, 7.757, 11.86, 16.062,
     0.9658, 3.9762, 7.7676, 11.874, 16.08, 1.1124, 4.1616, 7.938,
     12.014, 16.18, 1.1926, 4.3252, 8.17, 12.304, 16.522, 1.0442,
     4.099, 7.9144, 12.036, 16.254, 1.2096, 4.3498, 8.2024, 12.344,
     16.576, 0.90342, 3.85, 7.5978, 11.67, 15.848, 1.1834, 4.281,
     8.0934, 12.198, 16.396, 1.0908, 4.1304, 7.905, 11.982, 16.152,
     0.9582, 3.962, 7.7474, 11.846, 16.046, 0.96394, 3.9494, 7.7116,
     11.788, 15.964, 0.92248, 3.8868, 7.6508, 11.736, 15.922, 0.8499,
     3.755, 7.4888, 11.556, 15.734, 0.90988, 3.857, 7.61, 11.688,
     15.872, 0.82234, 3.6882, 7.391, 11.434, 15.592, 0.96714, 3.9624,
     7.7344, 11.822, 16.008, 0.94776, 3.9356, 7.7078, 11.796, 15.986,
     1.4158, 4.7124, 8.6864, 12.932, 17.256, 1.1388, 4.3252, 8.2552,
     12.482, 16.802, 1.68, 5.1116, 9.1526, 13.438, 17.786, 1.224,
     4.4692, 8.441, 12.704, 17.054, 1.2728, 4.5112, 8.46, 12.694,
     17.018, 1.0856, 4.227, 8.14, 12.362, 16.676, 1.0844, 4.194, 8.062,
     12.24, 16.516, 1.0606, 4.1536, 8.0124, 12.182, 16.456, 1.4158,
     4.776, 8.8104, 13.112, 17.488, 1.393, 4.7206, 8.7288, 13.01,
     17.368, 1.1104, 4.2704, 8.1884, 12.408, 16.726, 1.0346, 4.1582,
     8.0628, 12.28, 16.592, 1.047, 4.1514, 8.03, 12.222, 16.514,
     1.0194, 4.1246, 8.0188, 12.23, 16.542, 1.0346, 4.111, 7.9638,
     12.134, 16.404, 1.4284, 4.7316, 8.6998, 12.942, 17.258, 1.2358,
     4.4462, 8.3736, 12.59, 16.892, 1.2276, 4.473, 8.4444, 12.708,
     17.054, 1.1742, 4.3534, 8.271, 12.486, 16.794, 1.258, 4.5442,
     8.546, 12.834, 17.202, 1.1868, 4.3854, 8.3164, 12.542, 16.86,
     1.0848, 4.2498, 8.1914, 12.444, 16.792, 1.0652, 4.1588, 8.0194,
     12.192, 16.468, 1.4482, 4.8004, 8.8136, 13.092, 17.446, 1.0522,
     4.1844, 8.0974, 12.322, 16.644, 1.391, 4.7182, 8.726, 13.006,
     17.362, 1.3222, 4.6052, 8.5902, 12.856, 17.204, 1.1368, 4.34,
     8.3026, 12.566, 16.924, 1.175, 4.3784, 8.3282, 12.576, 16.916,
     1.0542, 4.1654, 8.0508, 12.25, 16.546, 1.3238, 4.5542, 8.4746,
     12.678, 16.968, 1.1512, 4.3224, 8.2322, 12.44, 16.744)),
.Names = c("Wafer",  "Site", "voltage", "current"),
row.names = as.character(1:400),
class = c("nmGroupedData", "groupedData", "data.frame"),
formula = current ~ voltage | Wafer/Site,
formulaList = structure(list(Wafer = ~Wafer, Site = ~Site),
  .Names = c("Wafer", "Site")),
labels = structure(list(x = "Voltage", y = "Intensity of Current"),
  .Names = c("x", "y")),
units = structure(list(x = "(V)", y = "(mA)"), .Names = c("x", "y")),
order.groups = structure(list(Wafer = TRUE, Site = TRUE),
  .Names = c("Wafer", "Site")),
FUN = function (x) max(x, na.rm = TRUE))
