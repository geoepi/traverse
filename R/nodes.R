#' Generate external terminal nodes
#'
#' Nodes are generated from the finite cells of the target template, then
#' displaced outward into the computational domain. They are not derived from
#' arbitrary finite cells in a movement surface.
#'
#' @param domain A `traverse_domain` object.
#' @param side One or more of `"south"`, `"north"`, `"east"`, and `"west"`.
#' @param distance Positive distance in map units outside the target edge.
#' @param n A positive scalar total or a named count for each side.
#' @param placement Either `"even"` or `"random"`.
#' @param allocation For scalar `n` across multiple sides, either `"equal"`
#'   or `"boundary_length"`.
#' @param seed Optional random seed used only for random placement.
#' @return A point [terra::SpatVector] with `id`, `side`, and `distance`
#'   attributes.
#' @export
traverse_nodes <- function(domain, side = c("south", "north"), distance, n,
                            placement = c("even", "random"),
                            allocation = c("equal", "boundary_length"), seed = NULL) {
  if (!inherits(domain, "traverse_domain")) traverse_stop("domain must be a traverse_domain object.")
  side <- match.arg(side, c("south", "north", "east", "west"), several.ok = TRUE)
  placement <- match.arg(placement)
  allocation <- match.arg(allocation)
  traverse_assert_projected(domain$target_template, "domain")
  if (length(distance) != 1L || !is.finite(distance) || distance <= 0) {
    traverse_stop("distance must be a single positive finite map-unit distance.")
  }
  if (length(n) < 1L || any(!is.finite(n)) || any(n != as.integer(n)) || any(n < 1L)) {
    traverse_stop("n must contain positive integer node counts.")
  }
  if (!is.null(seed) && (length(seed) != 1L || !is.finite(seed))) {
    traverse_stop("seed must be a single finite value or NULL.")
  }

  target_values <- matrix(
    terra::values(domain$target_mask, mat = FALSE),
    nrow = terra::nrow(domain$target_mask),
    ncol = terra::ncol(domain$target_mask),
    byrow = TRUE
  )
  edge_cells <- lapply(side, function(s) {
    if (s %in% c("north", "south")) {
      rows <- if (s == "north") apply(target_values, 2L, function(z) {
        q <- which(!is.na(z)); if (length(q)) min(q) else NA_integer_
      }) else apply(target_values, 2L, function(z) {
        q <- which(!is.na(z)); if (length(q)) max(q) else NA_integer_
      })
      data.frame(row = rows, col = seq_along(rows))[!is.na(rows), , drop = FALSE]
    } else {
      cols <- if (s == "west") apply(target_values, 1L, function(z) {
        q <- which(!is.na(z)); if (length(q)) min(q) else NA_integer_
      }) else apply(target_values, 1L, function(z) {
        q <- which(!is.na(z)); if (length(q)) max(q) else NA_integer_
      })
      data.frame(row = seq_along(cols), col = cols)[!is.na(cols), , drop = FALSE]
    }
  })
  if (any(vapply(edge_cells, nrow, integer(1L)) == 0L)) {
    traverse_stop("At least one requested side has no finite target cells.")
  }
  edge_lengths <- vapply(edge_cells, nrow, numeric(1L)) *
    vapply(side, function(s) if (s %in% c("north", "south")) terra::res(domain$target_template)[1L] else terra::res(domain$target_template)[2L], numeric(1L))
  counts <- traverse_allocate_counts(n, side, allocation, edge_lengths)
  restore_rng <- traverse_rng_restore(seed)
  on.exit(restore_rng(), add = TRUE)

  comp_ext <- traverse_ext_values(domain$computational_template)
  target_ext <- traverse_ext_values(domain$target_template)
  resolution <- terra::res(domain$target_template)
  available <- c(
    south = target_ext[3L] - comp_ext[3L],
    north = comp_ext[4L] - target_ext[4L],
    west = target_ext[1L] - comp_ext[1L],
    east = comp_ext[2L] - target_ext[2L]
  )
  half_cell <- vapply(side, function(s) if (s %in% c("north", "south")) resolution[2L] else resolution[1L], numeric(1L)) / 2
  if (any(distance + half_cell >= available - 1e-8)) {
    traverse_stop(paste0(
      "Terminal distance must be smaller than the available computational buffer; ",
      " increase domain$buffer or reduce distance."
    ))
  }

  result <- lapply(seq_along(side), function(i) {
    s <- side[i]
    cells <- edge_cells[[i]]
    axis_values <- if (s %in% c("north", "south")) {
      terra::xyFromCell(domain$target_template, terra::cellFromRowCol(domain$target_template, cells$row, cells$col))[, 1L]
    } else {
      terra::xyFromCell(domain$target_template, terra::cellFromRowCol(domain$target_template, cells$row, cells$col))[, 2L]
    }
    axis_values <- unname(axis_values)
    edge_xy <- terra::xyFromCell(
      domain$target_template,
      terra::cellFromRowCol(domain$target_template, cells$row, cells$col)
    )
    edge_axis <- if (s %in% c("north", "south")) edge_xy[, 2L] else edge_xy[, 1L]
    edge_axis <- unname(edge_axis)
    edge_order <- order(axis_values)
    boundary_at <- function(positions) {
      stats::approx(
        x = axis_values[edge_order], y = edge_axis[edge_order],
        xout = positions, rule = 2
      )$y
    }
    positions <- if (counts[i] == 1L) mean(range(axis_values)) else if (placement == "even") {
      seq(min(axis_values), max(axis_values), length.out = counts[i])
    } else if (diff(range(axis_values)) == 0) {
      rep(axis_values[1L], counts[i])
    } else {
      sort(stats::runif(counts[i], min(axis_values), max(axis_values)))
    }
    if (s == "north") {
      x <- positions
      y <- boundary_at(positions) + distance + half_cell[i]
    } else if (s == "south") {
      x <- positions
      y <- boundary_at(positions) - distance - half_cell[i]
    } else if (s == "east") {
      x <- boundary_at(positions) + distance + half_cell[i]
      y <- positions
    } else {
      x <- boundary_at(positions) - distance - half_cell[i]
      y <- positions
    }
    data.frame(x = unname(x), y = unname(y), side = s, distance = distance,
               stringsAsFactors = FALSE)
  })
  xy <- do.call(rbind, result)
  points <- terra::vect(xy, geom = c("x", "y"), crs = terra::crs(domain$computational_template, proj = TRUE))
  points$id <- sprintf("node%03d", seq_len(nrow(points)))

  coords <- terra::crds(points)
  inside_comp <- coords[, 1L] > comp_ext[1L] & coords[, 1L] < comp_ext[2L] &
    coords[, 2L] > comp_ext[3L] & coords[, 2L] < comp_ext[4L]
  inside_target <- terra::extract(domain$target_mask, points)[, 2L]
  if (any(!inside_comp) || any(!is.na(inside_target))) {
    traverse_stop("Generated terminal positions are not outside the target and inside the computational domain.")
  }
  points
}
