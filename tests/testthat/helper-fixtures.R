make_template <- function() {
  x <- terra::rast(nrows = 8, ncols = 10, xmin = 0, xmax = 10,
                   ymin = 0, ymax = 8, crs = "EPSG:3857")
  terra::values(x) <- 1
  x[4, 5] <- NA
  x[6, 8] <- NA
  x
}

make_domain <- function() {
  traverse_domain(make_template(), buffer = 2, buffer_units = "cells")
}
