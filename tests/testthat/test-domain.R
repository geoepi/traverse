test_that("domain preserves target geometry and holes", {
  template <- make_template()
  domain <- make_domain()
  expect_true(terra::compareGeom(domain$target_template, template, stopOnError = FALSE))
  expect_true(terra::compareGeom(domain$target_mask, template, stopOnError = FALSE))
  expect_equal(terra::res(domain$computational_template), terra::res(template))
  expect_equal(terra::crs(domain$computational_template), terra::crs(template))
  expect_equal(terra::nrow(domain$computational_template), terra::nrow(template) + 4)
  expect_equal(terra::ncol(domain$computational_template), terra::ncol(template) + 4)
  expect_true(is.na(domain$target_mask[4, 5]))
  expect_true(is.na(domain$target_mask[6, 8]))
})

test_that("map-unit buffers require a projected CRS", {
  template <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                          ymin = 0, ymax = 2, crs = "EPSG:4326")
  terra::values(template) <- 1
  expect_error(traverse_domain(template, 1, "map"), "projected CRS")
})
