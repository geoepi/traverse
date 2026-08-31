#!/usr/bin/env Rscript

# Run from the traverse RStudio project with either:
#   source("scripts/local-demo.R")
# or from a terminal:
#   Rscript scripts/local-demo.R
# Add --real to run the larger asset-based smoke test after the fast demo.
# Add --save-plot to write local-demo-domain.png in the project root.

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Run this script inside the traverse project.")
    path <- parent
  }
}

root <- find_project_root()
args <- commandArgs(trailingOnly = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, quiet = TRUE)
} else if (requireNamespace("traverse", quietly = TRUE)) {
  library(traverse)
} else {
  stop("Install the package development dependencies or install traverse before running this demo.")
}

library(terra)

cat("Running fast synthetic traverse demo...\n")
template <- terra::rast(
  nrows = 12, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 12,
  crs = "EPSG:3857"
)
terra::values(template) <- 1
template[6, 8] <- NA_real_

domain <- traverse_domain(template, buffer = 2, buffer_units = "cells")
surface <- traverse_surface(
  template, domain = domain, type = "conductance", outside_value = 1
)
sources <- traverse_nodes(
  domain, side = "south", distance = 1, n = 3, placement = "even"
)
targets <- traverse_nodes(
  domain, side = "north", distance = 1, n = 3, placement = "even"
)
result <- traverse(
  surface, domain, sources = sources, targets = targets,
  directions = 8, theta = 1, correction = "c",
  scale_correction = TRUE, aggregation = "mean"
)

stopifnot(
  terra::compareGeom(result$computational_passage,
                     domain$computational_template, stopOnError = FALSE),
  terra::compareGeom(result$target_passage,
                     domain$target_template, stopOnError = FALSE),
  nrow(result$pairs) == 9L,
  all(is.na(terra::extract(domain$target_mask, sources)[, 2L])),
  all(is.na(terra::extract(domain$target_mask, targets)[, 2L])),
  is.na(result$target_passage[6, 8])
)

cat("local-demo: OK\n")
cat("  source nodes:", nrow(sources), "\n")
cat("  target nodes:", nrow(targets), "\n")
cat("  source-target pairs:", nrow(result$pairs), "\n")
cat("  target dimensions:", terra::nrow(result$target_passage), "x",
    terra::ncol(result$target_passage), "\n")

save_plot <- "--save-plot" %in% args
if (interactive() || save_plot) {
  if (save_plot) {
    grDevices::png(file.path(root, "local-demo-domain.png"),
                   width = 1000, height = 700, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  plot(domain)
  terra::plot(sources, add = TRUE, col = "#1B9E77", pch = 16,
              main = "Traverse local demo: external terminals")
  terra::plot(targets, add = TRUE, col = "#7570B3", pch = 16)
  legend("topright", legend = c("sources", "targets"),
         col = c("#1B9E77", "#7570B3"), pch = 16, bty = "n")
}

if ("--real" %in% args) {
  cat("Running the full asset-based smoke test...\n")
  sys.source(file.path(root, "scripts", "real-data-smoke-test.R"),
             envir = new.env(parent = globalenv()))
}

invisible(result)
