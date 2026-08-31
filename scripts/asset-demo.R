#!/usr/bin/env Rscript

# Presentation-quality asset demonstration for traverse.
#
# Run from the RStudio project with:
#   source("scripts/asset-demo.R")
# Or from a terminal with:
#   Rscript scripts/asset-demo.R
# Add --save-figures to save PNGs under demo-output/.

# ---- Configuration ----------------------------------------------------------

n_nodes <- 1
save_figures <- FALSE
if ("--save-figures" %in% commandArgs(trailingOnly = TRUE)) save_figures <- TRUE

# One source and one target produce one randomized-shortest-path calculation.
# Increasing n_nodes creates n_nodes^2 source-target passage calculations.
# For example, n_nodes = 3 produces 9 pairwise calculations and can be
# substantially slower on the full example raster.

# ---- Load packages and example data ----------------------------------------

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
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(root, quiet = TRUE)
}

library(terra)
library(ggplot2)
library(traverse)
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("The asset demo uses patchwork for the comparison figure; install it first.")
}
library(patchwork)

template <- terra::rast(file.path(root, "assets", "template.tif"))
raw_surface <- terra::rast(file.path(root, "assets", "surface.tif"))

# template.tif: finite value 0 means an included target study-area cell;
# NA means outside the target study area. The 0 is a mask value and has no
# movement meaning.
# surface.tif: values are approximately 0-1000; larger values indicate greater
# suitability/preference and therefore conceptually greater conductance.

template_values <- terra::values(template, mat = FALSE)
surface_values <- terra::values(raw_surface, mat = FALSE)
finite_surface <- surface_values[is.finite(surface_values)]
cat("Example raster summary\n")
cat("  template dimensions:", terra::nrow(template), "x", terra::ncol(template), "\n")
cat("  template resolution:", paste(signif(terra::res(template), 8), collapse = " x "), "\n")
cat("  template extent:", paste(signif(c(terra::xmin(template), terra::xmax(template),
                                         terra::ymin(template), terra::ymax(template)), 8),
                                  collapse = ", "), "\n")
cat("  template CRS:", terra::crs(template, proj = TRUE), "\n")
cat("  target cells:", sum(is.finite(template_values)), "\n")
cat("  raw surface range:", paste(signif(range(finite_surface), 8), collapse = " to "), "\n")

# Small helper for the repeated raster-to-ggplot conversion. It leaves values
# unchanged and retains NA rows so the plotting layer can decide how to display
# excluded cells.
raster_data <- function(x) {
  data <- terra::as.data.frame(x, xy = TRUE, na.rm = FALSE)
  data.frame(x = data[["x"]], y = data[["y"]], value = data[[3L]])
}

template_data <- raster_data(template)
surface_data <- raster_data(raw_surface)

# ---- Inspect the target study domain ---------------------------------------

target_cells <- template_data[is.finite(template_data$value), , drop = FALSE]
p_template <- ggplot() +
  geom_raster(data = target_cells, aes(x = x, y = y), fill = "#3B82A0") +
  coord_equal(expand = FALSE) +
  labs(
    title = "Target study domain",
    subtitle = "Finite cells define the area retained for interpretation",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid = element_blank())
print(p_template)

# ---- Inspect the movement/suitability surface ------------------------------

p_surface <- ggplot(surface_data, aes(x = x, y = y, fill = value)) +
  geom_raster(na.rm = FALSE) +
  scale_fill_viridis_c(na.value = "transparent", name = "Suitability") +
  coord_equal(expand = FALSE) +
  labs(
    title = "Input movement surface",
    subtitle = "Higher values indicate greater suitability or preference",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
print(p_surface)

# ---- Define the computational domain --------------------------------------

domain <- traverse_domain(
  template,
  buffer = 150,
  buffer_units = "map"
)

# These particular projected example rasters use kilometers as their map units,
# so this is a 150-km computational buffer. The traverse package itself does
# not assume kilometers.
print(domain)

# ---- Prepare conductance ----------------------------------------------------

# These are demonstration configuration choices, not universal scientific
# defaults. The explicit transformation puts the example surface on a smaller
# conductance scale, minimum resolves finite zeros, and outside_value fills
# cells newly added by the computational buffer.
surface <- traverse_surface(
  raw_surface,
  domain = domain,
  type = "conductance",
  transform = function(x) x / 100,
  minimum = 0.01,
  outside_value = 0.01
)

# The raw surface also contains NA cells outside the target footprint but inside
# the original raster rectangle. For this demonstration, explicitly assign the
# same low positive conductance to those non-target cells, while preserving NA
# barriers within the target study area. No hidden package behavior does this.
target_on_computational <- terra::extend(
  domain$target_mask, domain$computational_template, fill = NA
)
prepared_values <- terra::values(surface$raster, mat = FALSE)
target_mask_values <- terra::values(target_on_computational, mat = FALSE)
prepared_values[is.na(prepared_values) & is.na(target_mask_values)] <- 0.01
terra::values(surface$raster) <- prepared_values

cat("Prepared conductance range:", paste(signif(range(prepared_values, na.rm = TRUE), 8),
                                          collapse = " to "), "\n")

# ---- Place external source and target nodes --------------------------------

sources <- traverse_nodes(
  domain,
  side = "south",
  distance = 100,
  n = n_nodes,
  placement = "even"
)
targets <- traverse_nodes(
  domain,
  side = "north",
  distance = 100,
  n = n_nodes,
  placement = "even"
)

source_data <- data.frame(terra::crds(sources))
names(source_data) <- c("x", "y")
source_data$terminal <- "Source"
target_data <- data.frame(terra::crds(targets))
names(target_data) <- c("x", "y")
target_data$terminal <- "Target"
terminal_data <- rbind(source_data, target_data)

# ---- Visualize the passage experiment --------------------------------------

computational_extent <- data.frame(
  xmin = terra::xmin(domain$computational_template),
  xmax = terra::xmax(domain$computational_template),
  ymin = terra::ymin(domain$computational_template),
  ymax = terra::ymax(domain$computational_template)
)
target_extent <- data.frame(
  xmin = terra::xmin(template), xmax = terra::xmax(template),
  ymin = terra::ymin(template), ymax = terra::ymax(template)
)
p_geometry <- ggplot() +
  geom_rect(
    data = computational_extent,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, color = "#4D4D4D", linewidth = 0.8
  ) +
  geom_raster(data = target_cells, aes(x = x, y = y), fill = "#D9EAF0") +
  geom_rect(
    data = target_extent,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, color = "#D95F02", linewidth = 0.8
  ) +
  geom_point(
    data = terminal_data,
    aes(x = x, y = y, shape = terminal, color = terminal),
    size = 3
  ) +
  scale_shape_manual(values = c(Source = 17, Target = 15)) +
  scale_color_manual(values = c(Source = "#1B9E77", Target = "#7570B3")) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Traverse computational geometry",
    subtitle = "Terminals are outside the target but inside the computational buffer",
    x = NULL, y = NULL, shape = NULL, color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
print(p_geometry)

# ---- Run traverse -----------------------------------------------------------

# theta is a model parameter, not a computational tuning parameter. This value
# demonstrates the workflow and is not scientifically calibrated for a specific
# organism or process. correction = "c" and scale_correction = TRUE are also
# explicit model configuration choices.
message("Running traverse passage analysis...")
result <- traverse(
  surface,
  domain,
  sources = sources,
  targets = targets,
  directions = 16,
  theta = 1,
  correction = "c",
  scale_correction = TRUE,
  aggregation = "mean"
)
message("Passage analysis complete.")

# ---- Inspect the result -----------------------------------------------------

target_passage_values <- terra::values(result$target_passage, mat = FALSE)
finite_target_passage <- target_passage_values[is.finite(target_passage_values)]
cat("Result diagnostics\n")
cat("  source nodes:", nrow(result$sources), "\n")
cat("  target nodes:", nrow(result$targets), "\n")
cat("  source-target pairs:", nrow(result$pairs), "\n")
cat("  computational dimensions:", terra::nrow(result$computational_passage), "x",
    terra::ncol(result$computational_passage), "\n")
cat("  target dimensions:", terra::nrow(result$target_passage), "x",
    terra::ncol(result$target_passage), "\n")
cat("  finite target passage cells:", length(finite_target_passage), "\n")
cat("  target passage summary:", paste(signif(c(
  min = min(finite_target_passage),
  median = median(finite_target_passage),
  mean = mean(finite_target_passage),
  max = max(finite_target_passage)
), 8), collapse = ", "), "\n")
if (length(finite_target_passage) > 0L && all(finite_target_passage == 0)) {
  message(
    "Note: this full-resolution asset run returns exact zero passage at " ,
    "theta = 1. The source and target are connected, but the configured " ,
    "gdistance calculation underflows numerically on this raster; the raw " ,
    "result is shown unchanged below."
  )
}
print(result)

# ---- Plot computational-domain passage ------------------------------------

computational_passage_data <- raster_data(result$computational_passage)
p_computational <- ggplot(computational_passage_data,
                          aes(x = x, y = y, fill = value)) +
  geom_raster(na.rm = FALSE) +
  geom_rect(
    data = target_extent,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, color = "white", linewidth = 0.7
  ) +
  geom_point(
    data = terminal_data,
    aes(x = x, y = y, shape = terminal, color = terminal),
    inherit.aes = FALSE, size = 2.5
  ) +
  scale_fill_viridis_c(na.value = "transparent", name = "Passage intensity") +
  scale_shape_manual(values = c(Source = 17, Target = 15)) +
  scale_color_manual(values = c(Source = "#1B9E77", Target = "#7570B3")) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Computational passage surface",
    subtitle = "Terminal effects remain visible outside the target study domain",
    x = NULL, y = NULL, shape = NULL, color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
print(p_computational)

# ---- Plot final target-domain passage --------------------------------------

target_passage_data <- raster_data(result$target_passage)
p_final <- ggplot(target_passage_data, aes(x = x, y = y, fill = value)) +
  geom_raster(na.rm = FALSE) +
  scale_fill_viridis_c(na.value = "transparent", name = "Passage intensity") +
  coord_equal(expand = FALSE) +
  labs(
    title = "Target-domain passage intensity",
    subtitle = "External computational buffer and terminal artifacts removed",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
print(p_final)

# ---- Compare the complete traverse story ----------------------------------

# These variables have different meanings, so each panel keeps its own fill
# scale rather than forcing a common legend across suitability and passage.
p_comparison <- (p_surface | p_computational | p_final) +
  patchwork::plot_annotation(
    title = "Asset-based traverse demonstration",
    subtitle = "Target area + suitability + external terminals -> passage -> clean target result"
  )
print(p_comparison)

# ---- Optionally save demonstration figures ---------------------------------

if (save_figures) {
  output_dir <- file.path(root, "demo-output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, "01-target-domain.png"), p_template,
                  width = 7, height = 5, dpi = 160)
  ggplot2::ggsave(file.path(output_dir, "02-input-surface.png"), p_surface,
                  width = 7, height = 5, dpi = 160)
  ggplot2::ggsave(file.path(output_dir, "03-computational-geometry.png"), p_geometry,
                  width = 7, height = 5, dpi = 160)
  ggplot2::ggsave(file.path(output_dir, "04-computational-passage.png"), p_computational,
                  width = 7, height = 5, dpi = 160)
  ggplot2::ggsave(file.path(output_dir, "05-target-passage.png"), p_final,
                  width = 7, height = 5, dpi = 160)
  ggplot2::ggsave(file.path(output_dir, "06-traverse-demo.png"), p_comparison,
                  width = 16, height = 5, dpi = 160)
  message("Saved demonstration figures to ", output_dir)
}

invisible(result)
