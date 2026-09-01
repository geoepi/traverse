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
#' @param scale_correction Logical; whether to apply `gdistance` geographic
#'   correction scaling. `TRUE` reproduces the historical prototype behavior
#'   and interacts with the interpretation of `theta`. Scaling changes the
#'   numerical scale of corrected transition values and can be disabled with
#'   `FALSE`.
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
                             scale_correction = TRUE,
                             aggregate = TRUE, aggregation = c("mean", "sum"),
                             workers = 1) {
  aggregation <- match.arg(aggregation)
  if (length(aggregate) != 1L || !is.logical(aggregate) || is.na(aggregate)) {
    traverse_stop("aggregate must be TRUE or FALSE.")
  }
  traverse_validate_model(directions, theta, correction, workers, scale_correction)
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

  transition_layer <- traverse_transition(
    surface_obj, directions, correction, scale_correction
  )
  source_sp <- traverse_point_to_sp(sources)
  target_sp <- traverse_point_to_sp(targets)
  pairwise <- if (!aggregate) vector("list", nrow(pairs)) else NULL
  if (!is.null(pairwise)) names(pairwise) <- as.character(pairs$pair_id)
  running_sum <- NULL
  valid_count <- NULL
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
    raw_values <- raster::getValues(raw)
    raster::values(pair_raster) <- raw_values
    pair_result <- terra::rast(pair_raster)
    terra::crs(pair_result) <- terra::crs(raster_surface, proj = TRUE)
    if (!aggregate) {
      pairwise[[i]] <- pair_result
    } else if (i == 1L) {
      running_sum <- pair_result
      running_values <- terra::values(running_sum, mat = FALSE)
      valid_values <- !is.na(running_values)
      running_values[!valid_values] <- 0
      terra::values(running_sum) <- running_values
      valid_count <- pair_result
      terra::values(valid_count) <- as.numeric(valid_values)
    } else {
      pair_values <- terra::values(pair_result, mat = FALSE)
      sum_values <- terra::values(running_sum, mat = FALSE)
      count_values <- terra::values(valid_count, mat = FALSE)
      valid_values <- !is.na(pair_values)
      sum_values[valid_values] <- sum_values[valid_values] + pair_values[valid_values]
      count_values[valid_values] <- count_values[valid_values] + 1
      terra::values(running_sum) <- sum_values
      terra::values(valid_count) <- count_values
    }
  }

  if (!aggregate) {
    return(list(pairwise = pairwise, pairs = pairs, sources = sources, targets = targets,
                aggregation = aggregation, directions = directions, theta = theta,
                correction = correction, scale_correction = scale_correction,
                workers = workers, number_of_pairs = nrow(pairs)))
  }
  output_values <- terra::values(running_sum, mat = FALSE)
  count_values <- terra::values(valid_count, mat = FALSE)
  if (aggregation == "mean") {
    output_values[count_values > 0] <- output_values[count_values > 0] /
      count_values[count_values > 0]
  }
  output_values[count_values == 0] <- NA_real_
  finite_output <- output_values[is.finite(output_values)]
  if (length(finite_output) > 0L && all(finite_output == 0)) {
    warning(
      paste(
        "Passage returned all zeros despite a valid source-target setup;",
        "check theta, conductance scaling, and geographic-correction settings",
        "for numerical degeneration."
      ),
      call. = FALSE
    )
  }
  output <- running_sum
  terra::values(output) <- output_values
  attr(output, "traverse_passage") <- list(
    pairs = pairs, sources = sources, targets = targets,
    aggregation = aggregation, directions = directions, theta = theta,
    correction = correction, scale_correction = scale_correction,
    workers = workers, number_of_pairs = nrow(pairs)
  )
  output
}
