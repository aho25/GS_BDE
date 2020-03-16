### Define function for data normalization
normalize <- function(x) {if (!all(is.na(x))) {(x + x/max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) + 1 - min(x, na.rm = TRUE))}}
