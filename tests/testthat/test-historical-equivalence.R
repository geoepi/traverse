test_that("historical and traverse workflows agree on a small fixture", {
  template <- terra::rast(
    nrows = 10, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 10,
    crs = "EPSG:3857"
  )
  terra::values(template) <- seq_len(terra::ncell(template)) / 100
  template[5, 7] <- NA_real_
  domain <- traverse_domain(template, buffer = 2, buffer_units = "cells")
  surface <- traverse_surface(
    template,
    domain = domain,
    type = "conductance",
    minimum = 0.01,
    outside_value = 0.01
  )
  sources <- traverse_nodes(domain, "south", distance = 1, n = 2)
  targets <- traverse_nodes(domain, "north", distance = 1, n = 2)
  pairs <- traverse_pairs(sources, targets)

  legacy_raster <- traverse:::traverse_backend_raster(surface)
  raster::crs(legacy_raster) <- terra::crs(surface$raster, proj = TRUE)
  legacy_transition <- gdistance::geoCorrection(
    gdistance::transition(
      legacy_raster, transitionFunction = mean, directions = 8
    ),
    type = "c", multpl = FALSE, scl = TRUE
  )
  source_sp <- traverse:::traverse_point_to_sp(sources)
  target_sp <- traverse:::traverse_point_to_sp(targets)
  legacy_pairwise <- lapply(seq_len(nrow(pairs)), function(i) {
    raw <- gdistance::passage(
      legacy_transition,
      origin = source_sp[pairs$source_index[i], , drop = FALSE],
      goal = target_sp[pairs$target_index[i], , drop = FALSE],
      theta = 1e-100
    )
    pair <- raster::raster(raw)
    raster::values(pair) <- raster::getValues(raw)
    pair
  })
  legacy_mean <- raster::calc(
    raster::stack(legacy_pairwise),
    fun = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  )
  legacy_target <- raster::mask(
    raster::crop(legacy_mean, raster::raster(domain$target_template)),
    raster::raster(domain$target_mask)
  )
  legacy_target <- terra::rast(legacy_target)
  terra::crs(legacy_target) <- terra::crs(template, proj = TRUE)

  modern <- traverse(
    surface, domain, sources = sources, targets = targets, pairs = pairs,
    directions = 8, theta = 1e-100, correction = "c",
    scale_correction = TRUE, aggregation = "mean"
  )
  expect_true(terra::compareGeom(
    legacy_target, modern$target_passage, stopOnError = FALSE
  ))
  expect_equal(
    terra::values(legacy_target, mat = FALSE),
    terra::values(modern$target_passage, mat = FALSE),
    tolerance = 1e-12
  )
})
