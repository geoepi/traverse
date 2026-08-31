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

test_that("all source-target pairs are produced", {
  d <- make_domain()
  sources <- traverse_nodes(d, "south", 1, 3)
  targets <- traverse_nodes(d, "north", 1, 4)
  pairs <- traverse_pairs(sources, targets)
  expect_equal(nrow(pairs), 12)
  expect_equal(nrow(unique(pairs[c("source_id", "target_id")])), 12)
  expect_equal(length(unique(pairs$pair_id)), 12)
})
