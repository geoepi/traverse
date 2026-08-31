test_that("scaled transition matches the direct gdistance calculation", {
  d <- make_domain()
  s <- traverse_surface(make_template(), d, outside_value = 1)
  actual <- traverse:::traverse_transition(
    s, directions = 8, correction = "c", scale_correction = TRUE
  )
  legacy <- traverse:::traverse_backend_raster(s)
  raster::crs(legacy) <- terra::crs(s$raster, proj = TRUE)
  expected <- gdistance::geoCorrection(
    gdistance::transition(legacy, transitionFunction = mean, directions = 8),
    type = "c", multpl = FALSE, scl = TRUE
  )
  expect_equal(actual@transitionMatrix, expected@transitionMatrix)
  expect_error(
    traverse:::traverse_transition(s, correction = "not-a-correction"),
    "one of 'c', 'r', or NULL"
  )
})

test_that("zero conductance is retained until backend construction", {
  x <- make_template()
  x[1, 1] <- 0
  s <- traverse_surface(x, type = "conductance")
  expect_equal(terra::values(s$raster)[1], 0)
  expect_error(
    traverse:::traverse_transition(s, correction = NULL),
    "finite non-positive values"
  )
  resolved <- traverse_surface(x, type = "conductance", minimum = 0.01)
  expect_s4_class(traverse:::traverse_transition(resolved, correction = NULL), "TransitionLayer")
})

test_that("aggregate passage streams values without pairwise metadata", {
  d <- make_domain()
  s <- traverse_surface(make_template(), d, outside_value = 1)
  sources <- traverse_nodes(d, "south", 1, 2)
  targets <- traverse_nodes(d, "north", 1, 2)
  explicit <- traverse_passage(s, sources, targets, aggregate = FALSE,
                               directions = 8, scale_correction = TRUE)
  streamed_mean <- traverse_passage(s, sources, targets, directions = 8,
                                    aggregation = "mean", scale_correction = TRUE)
  streamed_sum <- traverse_passage(s, sources, targets, directions = 8,
                                   aggregation = "sum", scale_correction = TRUE)
  stack <- explicit$pairwise[[1L]]
  if (length(explicit$pairwise) > 1L) {
    for (i in 2:length(explicit$pairwise)) stack <- c(stack, explicit$pairwise[[i]])
  }
  expected_mean <- terra::app(stack, fun = function(v) {
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  })
  expected_sum <- terra::app(stack, fun = function(v) {
    if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)
  })
  expect_equal(as.numeric(terra::values(streamed_mean)),
               as.numeric(terra::values(expected_mean)))
  expect_equal(as.numeric(terra::values(streamed_sum)),
               as.numeric(terra::values(expected_sum)))
  metadata <- attr(streamed_mean, "traverse_passage")
  expect_null(metadata$pairwise)
  expect_equal(metadata$number_of_pairs, 4)
})

test_that("top-level random terminal generation uses one stream and restores RNG", {
  d <- make_domain()
  s <- traverse_surface(make_template(), d, outside_value = 1)
  set.seed(123)
  caller_state <- .Random.seed
  first <- traverse(s, d, node_distance = 1, n = 1, placement = "random",
                    directions = 4, seed = 77)
  expect_equal(.Random.seed, caller_state)
  second <- traverse(s, d, node_distance = 1, n = 1, placement = "random",
                     directions = 4, seed = 77)
  expect_equal(terra::crds(first$sources), terra::crds(second$sources))
  expect_equal(terra::crds(first$targets), terra::crds(second$targets))

  restore <- traverse:::traverse_rng_restore(77)
  expected_sources <- traverse_nodes(d, "south", distance = 1, n = 1,
                                     placement = "random", seed = NULL)
  expected_targets <- traverse_nodes(d, "north", distance = 1, n = 1,
                                     placement = "random", seed = NULL)
  restore()
  expect_equal(terra::crds(first$sources), terra::crds(expected_sources))
  expect_equal(terra::crds(first$targets), terra::crds(expected_targets))
})
