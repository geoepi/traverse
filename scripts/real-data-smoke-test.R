#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else "scripts/real-data-smoke-test.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

library(terra)
library(traverse)

template <- rast(file.path(root, "assets", "template.tif"))
raw_surface <- rast(file.path(root, "assets", "surface.tif"))
domain <- traverse_domain(template, buffer = 150, buffer_units = "map")
surface <- traverse_surface(
  raw_surface, domain = domain, type = "conductance",
  transform = function(x) x / 100, minimum = 0.01, outside_value = 0.01
)
# The raw surface contains NA cells outside the target footprint inside its
# original rectangle. For this demonstration, explicitly assign the same
# positive terminal-zone value there while leaving target NA barriers intact.
target_on_computational <- terra::extend(domain$target_mask,
                                         domain$computational_template, fill = NA)
surface_values <- terra::values(surface$raster, mat = FALSE)
target_values <- terra::values(target_on_computational, mat = FALSE)
surface_values[is.na(surface_values) & is.na(target_values)] <- 0.01
terra::values(surface$raster) <- surface_values
sources <- traverse_nodes(domain, side = "south", distance = 100, n = 1)
targets <- traverse_nodes(domain, side = "north", distance = 100, n = 1)
result <- traverse(
  surface, domain, sources = sources, targets = targets,
  directions = 16, theta = 1, correction = "c", scale_correction = TRUE,
  aggregation = "mean"
)

stopifnot(
  terra::compareGeom(result$computational_passage, domain$computational_template,
                     stopOnError = FALSE),
  terra::compareGeom(result$target_passage, domain$target_template,
                     stopOnError = FALSE),
  all(is.na(terra::extract(domain$target_mask, sources)[, 2])),
  all(is.na(terra::extract(domain$target_mask, targets)[, 2])),
  any(is.finite(terra::values(result$target_passage)))
)

cat("real-data-smoke-test: OK\n")
cat("  target cells:", sum(!is.na(terra::values(domain$target_mask))), "\n")
cat("  pairs:", nrow(result$pairs), "\n")
cat("  computational dimensions:", terra::nrow(result$computational_passage), "x",
    terra::ncol(result$computational_passage), "\n")
