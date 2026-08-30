#' Construct source-target pairs
#'
#' @param sources Source point [terra::SpatVector] or `sf` object.
#' @param targets Target point [terra::SpatVector] or `sf` object.
#' @param method Pairing method. The initial implementation supports `"all"`.
#' @return A data frame with stable pair, source, and target identifiers.
#' @export
traverse_pairs <- function(sources, targets, method = "all") {
  if (length(method) != 1L || !identical(method, "all")) {
    traverse_stop("The only supported pairing method is method = 'all'.")
  }
  sources <- traverse_standardize_points(sources, "source")
  targets <- traverse_standardize_points(targets, "target")
  source_id <- as.character(sources$id)
  target_id <- as.character(targets$id)
  grid <- expand.grid(source_index = seq_along(source_id), target_index = seq_along(target_id))
  data.frame(
    pair_id = sprintf("pair%04d", seq_len(nrow(grid))),
    source_id = source_id[grid$source_index],
    target_id = target_id[grid$target_index],
    source_index = grid$source_index,
    target_index = grid$target_index,
    stringsAsFactors = FALSE
  )
}
