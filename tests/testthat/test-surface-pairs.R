test_that("surface semantics and alignment are explicit", {
  d <- make_domain()
  x <- make_template()
  s <- traverse_surface(x, d, type = "conductance", outside_value = 2,
                        minimum = 0.1)
  expect_s3_class(s, "traverse_surface")
  expect_equal(s$type, "conductance")
  expect_true(terra::compareGeom(s$raster, d$computational_template, stopOnError = FALSE))
  expect_equal(terra::values(s$raster)[1], 2)
  r <- traverse_surface(x, d, type = "resistance", outside_value = 1)
  expect_equal(r$type, "resistance")
  expect_error(traverse_surface(x, d, outside_value = NULL), "outside_value")
})

test_that("target-domain surface NA values are reported and preserved", {
  template <- terra::rast(
    nrows = 5, ncols = 6, xmin = 0, xmax = 6, ymin = 0, ymax = 5,
    crs = "EPSG:3857"
  )
  terra::values(template) <- 1
  surface <- template
  surface[3, 4] <- NA_real_
  domain <- traverse_domain(template, buffer = 1, buffer_units = "cells")

  expect_warning(
    prepared <- traverse_surface(surface, domain = domain, outside_value = 1),
    "Surface contains 1 NA cells inside the target domain"
  )
  target_cell <- terra::cellFromRowCol(surface, 3, 4)
  cell <- terra::cellFromXY(
    prepared$raster,
    terra::xyFromCell(surface, target_cell)
  )
  expect_true(is.na(terra::values(prepared$raster, mat = FALSE)[cell]))
})

test_that("all source-target pairs are produced", {
  d <- make_domain()
  sources <- traverse_nodes(d, "south", 1, 3)
  targets <- traverse_nodes(d, "north", 1, 4)
  pairs <- traverse_pairs(sources, targets)
  expect_equal(nrow(pairs), 12)
  expect_equal(nrow(unique(pairs[c("source_id", "target_id")])), 12)
  expect_equal(length(unique(pairs$pair_id)), 12)
})
