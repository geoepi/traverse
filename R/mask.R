#' Restore passage to the target domain
#'
#' This masks the buffered computational result using the original target
#' template. It preserves internal `NA` holes and does not rely on rectangular
#' cropping alone.
#'
#' @param x A passage [terra::SpatRaster] covering the computational domain.
#' @param domain A `traverse_domain` object.
#' @return A [terra::SpatRaster] with exactly the target template geometry.
#' @export
traverse_mask <- function(x, domain) {
  traverse_assert_spatraster(x, "x")
  if (!inherits(domain, "traverse_domain")) traverse_stop("domain must be a traverse_domain object.")
  comp <- domain$computational_template
  if (!traverse_grid_equal(x, comp)) traverse_stop("x must match the computational-domain geometry.")
  cropped <- terra::crop(x, domain$target_template, snap = "near")
  masked <- terra::mask(cropped, domain$target_mask)
  if (!traverse_grid_equal(masked, domain$target_template)) {
    traverse_stop("Masking did not preserve the target raster geometry.")
  }
  masked
}
