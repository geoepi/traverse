test_that("nodes have exact counts, sides, and external placement", {
  domain <- make_domain()
  nodes <- traverse_nodes(domain, side = c("south", "north"), distance = 1,
                          n = 10, placement = "even")
  expect_s4_class(nodes, "SpatVector")
  expect_equal(nrow(nodes), 10)
  expect_setequal(unique(as.character(nodes$side)), c("south", "north"))
  target_values <- terra::extract(domain$target_mask, nodes)[, 2]
  expect_true(all(is.na(target_values)))
  comp <- terra::ext(domain$computational_template)
  xy <- terra::crds(nodes)
  expect_true(all(xy[, 1] > comp$xmin & xy[, 1] < comp$xmax &
                   xy[, 2] > comp$ymin & xy[, 2] < comp$ymax))
  expect_equal(length(unique(nodes$id)), 10)
})

test_that("random nodes are reproducible and all cardinal sides work", {
  domain <- make_domain()
  a <- traverse_nodes(domain, side = "north", distance = 1, n = 5,
                      placement = "random", seed = 42)
  b <- traverse_nodes(domain, side = "north", distance = 1, n = 5,
                      placement = "random", seed = 42)
  expect_equal(terra::crds(a), terra::crds(b))
  for (side in c("north", "south", "east", "west")) {
    nodes <- traverse_nodes(domain, side = side, distance = 1, n = 3)
    expect_equal(nrow(nodes), 3)
    expect_true(all(as.character(nodes$side) == side))
  }
})

test_that("terminal distance is checked against the computational buffer", {
  expect_error(
    traverse_nodes(make_domain(), side = "north", distance = 2, n = 2),
    "smaller than the available computational buffer"
  )
})

test_that("irregular target masks support cardinal edge placement", {
  x <- terra::rast(nrows = 7, ncols = 7, xmin = 0, xmax = 7, ymin = 0,
                   ymax = 7, crs = "EPSG:3857")
  terra::values(x) <- 1
  x[1, 1:2] <- NA
  x[7, 6:7] <- NA
  d <- traverse_domain(x, 2, "cells")
  for (side in c("north", "south", "east", "west")) {
    expect_equal(nrow(traverse_nodes(d, side, distance = 1, n = 4)), 4)
  }
})
