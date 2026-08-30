#' Calculate randomized-shortest-path passage
#'
#' Passage is reported as passage intensity or relative passage, not as a
#' probability. `theta` controls the randomized-shortest-path behavior and is
#' a scientific model parameter. It is intentionally configurable; no value is
#' claimed to be universally optimal. Resistance is converted internally to
#' conductance as `1 / resistance` for the `gdistance` backend.
#'
#' @param surface A `traverse_surface` or single-layer conductance
#'   [terra::SpatRaster].
#' @param sources Source points as a [terra::SpatVector] or `sf` object.
#' @param targets Target points as a [terra::SpatVector] or `sf` object.
#' @param pairs Optional output of [traverse_pairs()]. Defaults to all pairs.
#' @param directions Neighborhood size passed to `gdistance` (4, 8, or 16).
#' @param theta Positive randomized-shortest-path parameter.
#' @param correction Geographic correction type passed to `gdistance`, or
#'   `NULL` to omit correction.
#' @param aggregate Logical; if `TRUE`, return an aggregated raster. If
#'   `FALSE`, return pairwise rasters and metadata.
#' @param aggregation Either `"mean"` or `"sum"` when aggregating.
#' @param workers Explicit worker count reserved for future parallel execution.
#'   The safe default is one worker.
#' @return A [terra::SpatRaster] when `aggregate = TRUE`; otherwise a list with
#'   pairwise rasters and metadata.
#' @export
traverse_passage <- function(surface, sources, targets, pairs = NULL,
                             directions = 16, theta = 1, correction = "c",
                             aggregate = TRUE, aggregation = c("mean", "sum"),
                             workers = 1) {
  aggregation <- match.arg(aggregation)
  if (length(aggregate) != 1L || !is.logical(aggregate) || is.na(aggregate)) {
    traverse_stop("aggregate must be TRUE or FALSE.")
  }
  traverse_validate_model(directions, theta, correction, workers)
  raster_surface <- traverse_surface_raster(surface)
  surface_type <- traverse_surface_type(surface)
  surface_obj <- if (inherits(surface, "traverse_surface")) surface else
    traverse_surface(raster_surface, type = surface_type)
  sources <- traverse_standardize_points(sources, "source")
  targets <- traverse_standardize_points(targets, "target")
  traverse_point_crs_check(sources, raster_surface, "sources")
  traverse_point_crs_check(targets, raster_surface, "targets")
  if (is.null(pairs)) pairs <- traverse_pairs(sources, targets)
  required <- c("pair_id", "source_id", "target_id")
  if (!all(required %in% names(pairs))) traverse_stop("pairs must contain pair_id, source_id, and target_id columns.")
  if (nrow(pairs) < 1L) traverse_stop("pairs must contain at least one pair.")
  source_ids <- as.character(sources$id)
  target_ids <- as.character(targets$id)
  if (any(!as.character(pairs$source_id) %in% source_ids) ||
      any(!as.character(pairs$target_id) %in% target_ids)) {
    traverse_stop("pairs contains source_id or target_id values absent from the supplied points.")
  }
  cell <- terra::cellFromXY(raster_surface, terra::crds(rbind(sources, targets)))
  all_values <- terra::values(raster_surface, mat = FALSE)[cell]
  if (anyNA(cell) || any(!is.finite(all_values))) {
    traverse_stop("All source and target points must fall on finite movement-surface cells.")
  }

  transition_layer <- traverse_transition(surface_obj, directions, correction)
  source_sp <- traverse_point_to_sp(sources)
  target_sp <- traverse_point_to_sp(targets)
  pairwise <- vector("list", nrow(pairs))
  names(pairwise) <- as.character(pairs$pair_id)
  for (i in seq_len(nrow(pairs))) {
    source_index <- match(as.character(pairs$source_id[i]), source_ids)
    target_index <- match(as.character(pairs$target_id[i]), target_ids)
    raw <- gdistance::passage(
      transition_layer,
      origin = source_sp[source_index, , drop = FALSE],
      goal = target_sp[target_index, , drop = FALSE],
      theta = theta
    )
    pair_raster <- raster::raster(raw)
    raster::values(pair_raster) <- raster::getValues(raw)
    pairwise[[i]] <- terra::rast(pair_raster)
    terra::crs(pairwise[[i]]) <- terra::crs(raster_surface, proj = TRUE)
  }

  if (!aggregate) {
    return(list(pairwise = pairwise, pairs = pairs, sources = sources, targets = targets,
                aggregation = aggregation, directions = directions, theta = theta,
                correction = correction, workers = workers))
  }
  combined <- pairwise[[1L]]
  if (length(pairwise) > 1L) {
    for (i in 2:length(pairwise)) combined <- c(combined, pairwise[[i]])
  }
  fun <- if (aggregation == "mean") {
    function(values) if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
  } else {
    function(values) if (all(is.na(values))) NA_real_ else sum(values, na.rm = TRUE)
  }
  output <- terra::app(combined, fun = fun)
  attr(output, "traverse_passage") <- list(
    pairwise = pairwise, pairs = pairs, sources = sources, targets = targets,
    aggregation = aggregation, directions = directions, theta = theta,
    correction = correction, workers = workers
  )
  output
}
