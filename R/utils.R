traverse_stop <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

traverse_assert_spatraster <- function(x, name = "x") {
  if (!inherits(x, "SpatRaster")) {
    traverse_stop(sprintf("%s must be a terra::SpatRaster.", name))
  }
  if (terra::nlyr(x) != 1L) {
    traverse_stop(sprintf("%s must contain exactly one layer.", name))
  }
  invisible(x)
}

traverse_assert_projected <- function(x, name = "object") {
  crs_x <- terra::crs(x, proj = TRUE)
  if (!nzchar(crs_x) || isTRUE(terra::is.lonlat(x))) {
    traverse_stop(sprintf(
      "%s must use a projected CRS when map-unit distances are requested.",
      name
    ))
  }
  invisible(x)
}

traverse_ext_values <- function(x) {
  c(terra::xmin(x), terra::xmax(x), terra::ymin(x), terra::ymax(x))
}

traverse_crs_equal <- function(x, y) {
  isTRUE(terra::same.crs(x, y))
}

traverse_grid_equal <- function(x, y, tolerance = 1e-8) {
  if (!traverse_crs_equal(x, y)) return(FALSE)
  if (length(terra::res(x)) != 2L || any(abs(terra::res(x) - terra::res(y)) > tolerance)) {
    return(FALSE)
  }
  all(abs(traverse_ext_values(x) - traverse_ext_values(y)) <= tolerance)
}

traverse_grid_aligned <- function(x, y, tolerance = 1e-8) {
  if (!traverse_crs_equal(x, y)) return(FALSE)
  if (any(abs(terra::res(x) - terra::res(y)) > tolerance)) return(FALSE)
  ex <- traverse_ext_values(x)
  ey <- traverse_ext_values(y)
  offsets <- (ex[c(1L, 3L)] - ey[c(1L, 3L)]) / terra::res(y)
  all(abs(offsets - round(offsets)) <= tolerance)
}

traverse_point_input <- function(x, name = "points") {
  if (inherits(x, "sf") || inherits(x, "sfc")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      traverse_stop("sf input was supplied, but the sf package is not installed.")
    }
    x <- terra::vect(x)
  }
  if (!inherits(x, "SpatVector")) {
    traverse_stop(sprintf("%s must be a terra::SpatVector or sf object.", name))
  }
  if (!all(terra::geomtype(x) == "points")) {
    traverse_stop(sprintf("%s must contain point geometries.", name))
  }
  if (nrow(x) < 1L) traverse_stop(sprintf("%s must not be empty.", name))
  x
}

traverse_standardize_points <- function(x, prefix) {
  x <- traverse_point_input(x, prefix)
  ids <- if ("id" %in% names(x)) as.character(x$id) else {
    sprintf("%s%03d", prefix, seq_len(nrow(x)))
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    traverse_stop(sprintf("%s identifiers must be non-missing, non-empty, and unique.", prefix))
  }
  x$id <- ids
  x
}

traverse_point_crs_check <- function(points, reference, name) {
  if (!traverse_crs_equal(points, reference)) {
    traverse_stop(sprintf("%s and the computational raster must use the same CRS.", name))
  }
  invisible(points)
}

traverse_point_to_sp <- function(points) {
  pts <- traverse_point_input(points)
  sp::SpatialPoints(
    coords = terra::crds(pts),
    proj4string = sp::CRS(terra::crs(pts, proj = TRUE))
  )
}

traverse_surface_raster <- function(surface) {
  if (inherits(surface, "traverse_surface")) return(surface$raster)
  traverse_assert_spatraster(surface, "surface")
  surface
}

traverse_surface_type <- function(surface) {
  if (inherits(surface, "traverse_surface")) return(surface$type)
  "conductance"
}

traverse_rng_restore <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

traverse_allocate_counts <- function(n, side, allocation, edge_lengths = NULL) {
  k <- length(side)
  if (length(n) > 1L) {
    if (is.null(names(n)) || !all(side %in% names(n))) {
      traverse_stop("A vector n must be named for every requested side.")
    }
    counts <- as.integer(n[side])
  } else if (length(n) == 1L) {
    total <- as.integer(n)
    if (allocation == "equal") {
      counts <- rep.int(total %/% k, k)
      if (total %% k > 0L) counts[seq_len(total %% k)] <- counts[seq_len(total %% k)] + 1L
    } else {
      if (is.null(edge_lengths) || sum(edge_lengths) <= 0) {
        traverse_stop("Boundary lengths are unavailable for proportional allocation.")
      }
      raw <- total * edge_lengths / sum(edge_lengths)
      counts <- floor(raw)
      remainder <- total - sum(counts)
      if (remainder > 0L) {
        order_idx <- order(raw - counts, decreasing = TRUE)
        counts[order_idx[seq_len(remainder)]] <- counts[order_idx[seq_len(remainder)]] + 1L
      }
    }
  } else {
    traverse_stop("n must contain at least one value.")
  }
  if (any(!is.finite(counts)) || any(counts < 1L)) {
    traverse_stop("n must request at least one node per side.")
  }
  counts
}

traverse_backend_raster <- function(surface) {
  r <- surface$raster
  v <- terra::values(r, mat = FALSE)
  if (surface$type == "resistance") {
    if (any(is.finite(v) & v <= 0)) {
      traverse_stop("Resistance values must be strictly positive where finite.")
    }
    v[is.finite(v)] <- 1 / v[is.finite(v)]
    terra::values(r) <- v
  } else if (any(is.finite(v) & v <= 0)) {
    traverse_stop(paste(
      "The conductance surface contains finite non-positive values.",
      "Use minimum=, transform=, or NA barriers to define their",
      "movement interpretation before passage calculation."
    ))
  }
  if (!any(is.finite(v) & v > 0)) {
    traverse_stop("The movement surface must contain at least one positive finite value.")
  }
  raster::raster(r)
}

traverse_validate_model <- function(directions, theta, correction, workers,
                                    scale_correction = TRUE) {
  if (length(directions) != 1L || !is.finite(directions) ||
      directions != as.integer(directions) || !directions %in% c(4L, 8L, 16L)) {
    traverse_stop("directions must be one of 4, 8, or 16.")
  }
  if (length(theta) != 1L || !is.finite(theta) || theta <= 0) {
    traverse_stop("theta must be a single positive finite number.")
  }
  if (!is.null(correction) && (length(correction) != 1L || !is.character(correction) ||
                               !correction %in% c("c", "r"))) {
    traverse_stop("correction must be one of 'c', 'r', or NULL.")
  }
  if (length(scale_correction) != 1L || !is.logical(scale_correction) ||
      is.na(scale_correction)) {
    traverse_stop("scale_correction must be TRUE or FALSE.")
  }
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L ||
      workers != as.integer(workers)) {
    traverse_stop("workers must be a positive integer.")
  }
  invisible(TRUE)
}
