#' @export
print.traverse_result <- function(x, ...) {
  cat("<traverse_result>\n")
  cat("  sources:", nrow(x$sources), "; targets:", nrow(x$targets),
      "; pairs:", nrow(x$pairs), "\n")
  cat("  aggregation:", x$model$aggregation, "\n")
  cat("  target raster:", terra::nrow(x$target_passage), "x",
      terra::ncol(x$target_passage), "cells\n")
  cat("  computational raster:", terra::nrow(x$computational_passage), "x",
      terra::ncol(x$computational_passage), "cells\n")
  invisible(x)
}

#' @export
plot.traverse_result <- function(x, ...) {
  terra::plot(x$target_passage, main = "Traverse target-domain passage", ...)
  invisible(x)
}
