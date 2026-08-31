#' Prepare a conductance or resistance surface
#'
#' Conductance expresses ease of movement (higher values are easier); resistance
#' expresses movement cost (higher values are harder). No normalization or
#' project-specific multiplication is performed automatically.
#' Finite zero values are retained as raw/prepared surface data. Before passage,
#' finite conductance values must be made strictly positive with an explicit
#' `minimum`, `transform`, or `NA` barrier decision; the backend never assigns a
#' conductance floor silently.
#'
#' @param x A single-layer [terra::SpatRaster].
#' @param domain Optional `traverse_domain` whose computational geometry should
#'   be returned.
#' @param type Either `"conductance"` or `"resistance"`.
#' @param transform Optional function applied to `x` before alignment.
#' @param minimum Optional lower bound applied to finite values.
#' @param outside_value Fill value used when extending a surface to the
#'   computational domain. It is required when extension is needed.
#' @return An object of class `traverse_surface` containing the aligned raster,
#'   its semantic type, and provenance metadata.
#' @export
traverse_surface <- function(x, domain = NULL,
                              type = c("conductance", "resistance"),
                              transform = NULL, minimum = NULL,
                              outside_value = NULL) {
  traverse_assert_spatraster(x, "x")
  type <- match.arg(type)
  if (!is.null(domain) && !inherits(domain, "traverse_domain")) {
    traverse_stop("domain must be a traverse_domain object or NULL.")
  }
  if (!is.null(transform)) {
    if (!is.function(transform)) traverse_stop("transform must be a function or NULL.")
    x <- transform(x)
    traverse_assert_spatraster(x, "the transformed surface")
  }
  if (!is.null(minimum)) {
    if (length(minimum) != 1L || !is.finite(minimum)) traverse_stop("minimum must be a finite scalar or NULL.")
    values_x <- terra::values(x, mat = FALSE)
    values_x[is.finite(values_x)] <- pmax(values_x[is.finite(values_x)], minimum)
    terra::values(x) <- values_x
  }

  if (!is.null(domain)) {
    comp <- domain$computational_template
    if (!traverse_grid_aligned(x, comp)) {
      traverse_stop("x and the computational domain must have matching CRS, resolution, and grid alignment.")
    }
    ex <- traverse_ext_values(x)
    ec <- traverse_ext_values(comp)
    inside_comp <- ex[1L] >= ec[1L] - 1e-8 && ex[2L] <= ec[2L] + 1e-8 &&
      ex[3L] >= ec[3L] - 1e-8 && ex[4L] <= ec[4L] + 1e-8
    covers_comp <- ex[1L] <= ec[1L] + 1e-8 && ex[2L] >= ec[2L] - 1e-8 &&
      ex[3L] <= ec[3L] + 1e-8 && ex[4L] >= ec[4L] - 1e-8
    if (traverse_grid_equal(x, comp)) {
      aligned <- x
    } else if (inside_comp) {
      if (is.null(outside_value)) {
        traverse_stop("outside_value must be supplied when extending x to the computational domain.")
      }
      aligned <- terra::extend(x, comp, fill = outside_value)
    } else if (covers_comp) {
      aligned <- terra::crop(x, comp, snap = "near")
    } else {
      traverse_stop("x must be contained by or contain the computational domain on an aligned grid.")
    }
    if (!traverse_grid_equal(aligned, comp)) {
      traverse_stop("Could not align x to the computational domain without changing its grid.")
    }
  } else {
    aligned <- x
  }
  structure(
    list(raster = aligned, type = type, domain = domain,
         minimum = minimum, outside_value = outside_value,
         transform = transform),
    class = "traverse_surface"
  )
}

#' @export
print.traverse_surface <- function(x, ...) {
  cat("<traverse_surface>", x$type, "\n")
  cat("  geometry:", terra::nrow(x$raster), "x", terra::ncol(x$raster), "cells\n")
  cat("  finite cells:", sum(is.finite(terra::values(x$raster))), "\n")
  invisible(x)
}
