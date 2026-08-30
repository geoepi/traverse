#' Run the complete traverse workflow
#'
#' The workflow computes on the computational domain and interprets the final
#' result on the target domain. If either terminal set is omitted, it is
#' generated using the requested cardinal side and node settings.
#'
#' @param surface A `traverse_surface` or single-layer [terra::SpatRaster].
#' @param domain A `traverse_domain` object.
#' @param sources Optional source points.
#' @param targets Optional target points.
#' @param pairs Optional source-target pair table.
#' @param source_side Side used for generated sources.
#' @param target_side Side used for generated targets.
#' @param node_distance Map-unit distance for generated terminals.
#' @param n Node count for each generated terminal set.
#' @param placement Node placement method.
#' @param allocation Node allocation method.
#' @param seed Optional random seed.
#' @param directions Neighborhood size for the transition model.
#' @param theta Randomized-shortest-path parameter.
#' @param correction Geographic correction passed to `gdistance`, or `NULL`.
#' @param aggregation Pair aggregation method, `"mean"` or `"sum"`.
#' @param workers Explicit worker count; defaults to one.
#' @return A `traverse_result` object containing both target and computational
#'   passage rasters.
#' @export
traverse <- function(surface, domain, sources = NULL, targets = NULL, pairs = NULL,
                     source_side = "south", target_side = "north",
                     node_distance = NULL, n = NULL,
                     placement = "even", allocation = "equal", seed = NULL,
                     directions = 16, theta = 1, correction = "c",
                     aggregation = "mean", workers = 1) {
  if (!inherits(domain, "traverse_domain")) traverse_stop("domain must be a traverse_domain object.")
  if (is.null(sources) || is.null(targets)) {
    if (is.null(node_distance) || is.null(n)) {
      traverse_stop("node_distance and n are required when sources or targets are omitted.")
    }
    if (is.null(sources)) sources <- traverse_nodes(domain, side = source_side,
                                                    distance = node_distance, n = n,
                                                    placement = placement, allocation = allocation, seed = seed)
    if (is.null(targets)) targets <- traverse_nodes(domain, side = target_side,
                                                    distance = node_distance, n = n,
                                                    placement = placement, allocation = allocation, seed = seed)
  }
  surface_obj <- if (inherits(surface, "traverse_surface")) {
    if (!is.null(surface$domain) && !traverse_grid_equal(surface$raster, domain$computational_template)) {
      traverse_stop("surface does not match the supplied computational domain.")
    }
    if (traverse_grid_equal(surface$raster, domain$computational_template)) surface else {
      traverse_surface(
        surface$raster, domain = domain, type = surface$type,
        minimum = NULL, outside_value = surface$outside_value
      )
    }
  } else {
    traverse_surface(surface, domain = domain)
  }
  computational <- traverse_passage(
    surface_obj, sources, targets, pairs = pairs, directions = directions,
    theta = theta, correction = correction, aggregation = aggregation,
    workers = workers
  )
  metadata <- attr(computational, "traverse_passage")
  if (is.null(metadata)) traverse_stop("Passage metadata were not retained.")
  target <- traverse_mask(computational, domain)
  structure(
    list(
      target_passage = target,
      computational_passage = computational,
      sources = metadata$sources,
      targets = metadata$targets,
      pairs = metadata$pairs,
      domain = domain,
      surface = surface_obj,
      model = list(
        directions = directions, theta = theta, correction = correction,
        aggregation = aggregation, workers = workers
      )
    ),
    class = "traverse_result"
  )
}
