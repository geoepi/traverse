#' Define target and computational domains
#'
#' A target domain is the area of scientific interpretation. The computational
#' domain expands that geometry by a buffer so terminals can be placed outside
#' the target and terminal artifacts can be masked from final outputs.
#'
#' @param template A single-layer [terra::SpatRaster]. Finite cells define the
#'   target domain; `NA` and non-finite cells are excluded.
#' @param buffer Non-negative buffer distance. A positive value is required for
#'   a computational domain larger than the target.
#' @param buffer_units Either `"map"` for map units or `"cells"` for cell
#'   counts. Map-unit distances require a projected CRS.
#' @return An object of class `traverse_domain`.
#' @export
#' @examples
#' template <- terra::rast(nrows = 5, ncols = 6, xmin = 0, xmax = 6,
#'   ymin = 0, ymax = 5, crs = "EPSG:3857")
#' terra::values(template) <- 1
#' domain <- traverse_domain(template, buffer = 2, buffer_units = "cells")
#' domain
traverse_domain <- function(template, buffer, buffer_units = c("map", "cells")) {
  traverse_assert_spatraster(template, "template")
  buffer_units <- match.arg(buffer_units)
  if (length(buffer) != 1L || !is.finite(buffer) || buffer <= 0) {
    traverse_stop("buffer must be a single positive finite value.")
  }
  if (buffer_units == "map") traverse_assert_projected(template, "template")

  target_values <- terra::values(template, mat = FALSE)
  if (!any(is.finite(target_values))) {
    traverse_stop("template must contain at least one finite target cell.")
  }
  target_mask <- template
  terra::values(target_mask) <- ifelse(is.finite(target_values), 1, NA_real_)

  resolution <- terra::res(template)
  padding <- if (buffer_units == "cells") {
    if (buffer != as.integer(buffer)) {
      traverse_stop("buffer must be an integer when buffer_units = 'cells'.")
    }
    rep.int(as.integer(buffer), 2L)
  } else {
    ceiling(buffer / resolution)
  }
  ext_values <- traverse_ext_values(template)
  computational_template <- terra::extend(template, padding, fill = NA)
  terra::values(computational_template) <- 1
  computational_extent <- terra::ext(computational_template)

  structure(
    list(
      target_template = template,
      target_mask = target_mask,
      target_extent = terra::ext(template),
      computational_extent = terra::ext(computational_template),
      computational_template = computational_template,
      buffer = buffer,
      buffer_units = buffer_units,
      buffer_map = c(x = padding[1L] * resolution[1L], y = padding[2L] * resolution[2L]),
      padding_cells = c(x = padding[1L], y = padding[2L]),
      crs = terra::crs(template, proj = TRUE),
      resolution = resolution
    ),
    class = "traverse_domain"
  )
}

#' @export
print.traverse_domain <- function(x, ...) {
  cat("<traverse_domain>\n")
  cat("  target:", terra::nrow(x$target_template), "x", terra::ncol(x$target_template),
      "cells;", sum(!is.na(terra::values(x$target_mask))), "included cells\n")
  cat("  computational:", terra::nrow(x$computational_template), "x",
      terra::ncol(x$computational_template), "cells\n")
  cat("  buffer:", x$buffer, x$buffer_units, "\n")
  cat("  CRS:", if (nzchar(x$crs)) x$crs else "none", "\n")
  invisible(x)
}

#' @export
plot.traverse_domain <- function(x, ...) {
  terra::plot(x$computational_template, main = "Target and computational domains", ...)
  terra::plot(x$target_mask, add = TRUE, legend = FALSE, col = "#D95F02", alpha = 0.35)
  terra::lines(terra::as.polygons(x$target_mask, dissolve = TRUE), col = "#D95F02", lwd = 2)
  invisible(x)
}
