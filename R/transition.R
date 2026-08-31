traverse_transition <- function(surface, directions = 16, correction = "c",
                                scale_correction = TRUE) {
  traverse_validate_model(
    directions, theta = 1, correction = correction, workers = 1,
    scale_correction = scale_correction
  )
  legacy_raster <- traverse_backend_raster(surface)
  raster::crs(legacy_raster) <- terra::crs(surface$raster, proj = TRUE)
  transition_layer <- gdistance::transition(
    legacy_raster,
    transitionFunction = mean,
    directions = directions
  )
  if (!is.null(correction)) {
    transition_layer <- gdistance::geoCorrection(
      transition_layer, type = correction, multpl = FALSE,
      scl = scale_correction
    )
  }
  transition_layer
}
