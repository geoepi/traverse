traverse_transition <- function(surface, directions = 16, correction = "c") {
  traverse_validate_model(directions, theta = 1, correction = correction, workers = 1)
  legacy_raster <- traverse_backend_raster(surface)
  raster::crs(legacy_raster) <- terra::crs(surface$raster, proj = TRUE)
  transition_layer <- gdistance::transition(
    legacy_raster,
    transitionFunction = mean,
    directions = directions
  )
  if (!is.null(correction)) {
    transition_layer <- gdistance::geoCorrection(transition_layer, type = correction)
  }
  transition_layer
}
