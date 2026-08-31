test_that("passage and masking work on a small synthetic regression case", {
  d <- traverse_domain(make_template(), buffer = 2, buffer_units = "cells")
  s <- traverse_surface(make_template(), d, outside_value = 1)
  sources <- traverse_nodes(d, "south", distance = 1, n = 2)
  targets <- traverse_nodes(d, "north", distance = 1, n = 2)
  pairs <- traverse_pairs(sources, targets)
  mean_passage <- traverse_passage(s, sources, targets, pairs = pairs,
                                   directions = 16, theta = 1,
                                   aggregation = "mean")
  sum_passage <- traverse_passage(s, sources, targets, pairs = pairs,
                                  directions = 16, theta = 1,
                                  aggregation = "sum")
  expect_s4_class(mean_passage, "SpatRaster")
  expect_true(terra::compareGeom(mean_passage, d$computational_template, stopOnError = FALSE))
  expect_true(all(terra::values(sum_passage) >= terra::values(mean_passage), na.rm = TRUE))
  expect_equal(terra::values(mean_passage), terra::values(
    traverse_passage(s, sources, targets, pairs = pairs, directions = 16, theta = 1)
  ))
  target <- traverse_mask(mean_passage, d)
  expect_true(terra::compareGeom(target, d$target_template, stopOnError = FALSE))
  expect_true(all(is.na(target[4, 5])))
  expect_true(all(is.na(target[6, 8])))
})

test_that("high-level traverse retains both interpretations", {
  d <- make_domain()
  s <- traverse_surface(make_template(), d, outside_value = 1)
  result <- traverse(s, d, node_distance = 1, n = 2, directions = 8)
  expect_s3_class(result, "traverse_result")
  expect_s4_class(result$target_passage, "SpatRaster")
  expect_s4_class(result$computational_passage, "SpatRaster")
  expect_true(terra::compareGeom(result$target_passage, d$target_template, stopOnError = FALSE))
  expect_true(terra::compareGeom(result$computational_passage, d$computational_template,
                                 stopOnError = FALSE))
  expect_equal(nrow(result$pairs), 4)
})
