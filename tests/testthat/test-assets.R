asset_file <- function(name) testthat::test_path("..", "..", "assets", name)

test_that("template asset preserves finite zero-valued target cells", {
  skip_if_not(file.exists(asset_file("template.tif")), "repository assets are not in the package check payload")
  template <- terra::rast(asset_file("template.tif"))
  domain <- traverse_domain(template, buffer = 150, buffer_units = "map")
  template_values <- terra::values(template, mat = FALSE)
  expect_equal(c(terra::nrow(template), terra::ncol(template)), c(202, 293))
  expect_true(any(is.finite(template_values) & template_values == 0))
  expect_equal(sum(!is.na(terra::values(domain$target_mask))), sum(is.finite(template_values)))
  expect_true(all(is.na(terra::values(domain$target_mask)) == is.na(template_values)))
  expect_true(terra::compareGeom(domain$target_template, template, stopOnError = FALSE))
  expect_true(terra::same.crs(domain$computational_template, template))
  expect_true(all(abs(terra::res(domain$computational_template) -
                         c(24.99502, 24.94366)) < 1e-4))
  expect_true(terra::nrow(domain$computational_template) > terra::nrow(template))
  expect_true(terra::ncol(domain$computational_template) > terra::ncol(template))
})

test_that("asset nodes use projected map-unit distances", {
  skip_if_not(file.exists(asset_file("template.tif")), "repository assets are not in the package check payload")
  template <- terra::rast(asset_file("template.tif"))
  domain <- traverse_domain(template, buffer = 150, buffer_units = "map")
  south <- traverse_nodes(domain, "south", distance = 100, n = 10)
  north <- traverse_nodes(domain, "north", distance = 100, n = 10)
  expect_equal(nrow(south), 10)
  expect_equal(nrow(north), 10)
  expect_true(all(as.character(south$side) == "south"))
  expect_true(all(as.character(north$side) == "north"))
  expect_true(all(is.na(terra::extract(domain$target_mask, south)[, 2])))
  expect_true(all(is.na(terra::extract(domain$target_mask, north)[, 2])))
  expect_gt(diff(range(terra::crds(south)[, 1])), 0)
  expect_gt(diff(range(terra::crds(north)[, 1])), 0)
  expect_equal(terra::crds(traverse_nodes(domain, "north", 100, 10, seed = 17)),
               terra::crds(traverse_nodes(domain, "north", 100, 10, seed = 17)))
})

test_that("asset surface preparation is explicit and backend-ready", {
  skip_if_not(file.exists(asset_file("template.tif")) && file.exists(asset_file("surface.tif")),
              "repository assets are not in the package check payload")
  template <- terra::rast(asset_file("template.tif"))
  raw_surface <- terra::rast(asset_file("surface.tif"))
  expect_true(terra::compareGeom(template, raw_surface, stopOnError = FALSE))
  template_values <- terra::values(template, mat = FALSE)
  raw_surface_values <- terra::values(raw_surface, mat = FALSE)
  expect_equal(
    sum(is.finite(template_values) & !is.finite(raw_surface_values)),
    0L
  )
  raw_prepared <- traverse_surface(raw_surface, type = "conductance")
  raw_idx <- which(is.finite(terra::values(raw_surface)))[1]
  expect_equal(terra::values(raw_prepared$raster)[raw_idx],
               terra::values(raw_surface)[raw_idx])
  domain <- traverse_domain(template, buffer = 150, buffer_units = "map")
  surface <- traverse_surface(
    raw_surface, domain = domain, type = "conductance",
    transform = function(x) x / 100, minimum = 0.01, outside_value = 0.01
  )
  expect_true(terra::compareGeom(surface$raster, domain$computational_template,
                                 stopOnError = FALSE))
  xy <- terra::xyFromCell(surface$raster, seq_len(terra::ncell(surface$raster)))
  ext <- terra::ext(template)
  outside <- xy[, 1] < ext$xmin | xy[, 1] > ext$xmax |
    xy[, 2] < ext$ymin | xy[, 2] > ext$ymax
  surface_values <- terra::values(surface$raster, mat = FALSE)
  expect_true(all(surface_values[outside] == 0.01))
  expect_true(all(surface_values[is.finite(surface_values)] > 0))
  idx <- which(is.finite(terra::values(raw_surface)) & terra::values(raw_surface) > 1)[1]
  comp_idx <- terra::cellFromXY(surface$raster, terra::xyFromCell(raw_surface, idx))
  expect_equal(surface_values[comp_idx], terra::values(raw_surface)[idx] / 100)
})
