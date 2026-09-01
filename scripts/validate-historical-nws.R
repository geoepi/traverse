#!/usr/bin/env Rscript

# Quantitatively validate traverse against the historical NWS passage workflow.
#
# This is a development/validation script, not package API. The reference
# implementation below deliberately uses raster, sp, and gdistance so that the
# old computational path remains visible and isolated from traverse itself.

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        dir.exists(file.path(path, "R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the traverse project root.", call. = FALSE)
    }
    path <- parent
  }
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
} else {
  file.path(getwd(), "scripts", "validate-historical-nws.R")
}
script_root <- dirname(script_path)
project_root <- if (dir.exists(file.path(script_root, "R")) &&
                    file.exists(file.path(script_root, "DESCRIPTION"))) {
  find_project_root(script_root)
} else {
  find_project_root(getwd())
}

load_traverse <- function(root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
  } else if (!requireNamespace("traverse", quietly = TRUE)) {
    stop(
      "Install traverse or the package development dependencies before running validation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

historical_reference <- function(surface, domain, sources, targets, pairs,
                                 directions, theta, correction,
                                 scale_correction) {
  # This section intentionally mirrors the historical raster/sp/gdistance
  # implementation. It does not call traverse() or traverse_passage().
  legacy_raster <- traverse:::traverse_backend_raster(surface)
  raster::crs(legacy_raster) <- terra::crs(surface$raster, proj = TRUE)

  transition_layer <- gdistance::transition(
    legacy_raster,
    transitionFunction = mean,
    directions = directions
  )
  if (!is.null(correction)) {
    transition_layer <- gdistance::geoCorrection(
      transition_layer,
      type = correction,
      multpl = FALSE,
      scl = scale_correction
    )
  }

  source_sp <- traverse:::traverse_point_to_sp(sources)
  target_sp <- traverse:::traverse_point_to_sp(targets)
  source_ids <- as.character(sources$id)
  target_ids <- as.character(targets$id)
  pairwise <- lapply(seq_len(nrow(pairs)), function(i) {
    source_index <- match(as.character(pairs$source_id[i]), source_ids)
    target_index <- match(as.character(pairs$target_id[i]), target_ids)
    raw <- gdistance::passage(
      transition_layer,
      origin = source_sp[source_index, , drop = FALSE],
      goal = target_sp[target_index, , drop = FALSE],
      theta = theta
    )
    pair_raster <- raster::raster(raw)
    raster::values(pair_raster) <- raster::getValues(raw)
    pair_raster
  })

  # Historical pairwise aggregation was performed after all pair rasters were
  # calculated. Keep this explicit rather than using traverse's streaming
  # aggregator.
  pairwise_stack <- raster::stack(pairwise)
  passage_overlay <- raster::calc(
    pairwise_stack,
    fun = function(x) {
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    }
  )

  # The old workflow cropped and masked the aggregate back to the target. The
  # mask is 1 on finite target cells and NA elsewhere, preserving holes.
  target_template <- raster::raster(domain$target_template)
  target_mask <- raster::raster(domain$target_mask)
  target_result <- raster::mask(
    raster::crop(passage_overlay, target_template, snap = "near"),
    target_mask
  )

  computational <- terra::rast(passage_overlay)
  target <- terra::rast(target_result)
  terra::crs(computational) <- terra::crs(surface$raster, proj = TRUE)
  terra::crs(target) <- terra::crs(domain$target_template, proj = TRUE)
  list(
    computational_passage = computational,
    target_passage = target,
    pairwise = pairwise
  )
}

repair_historical_surface <- function(template, raw_surface, domain,
                                      transform = function(x) x / 100,
                                      minimum = 0.01,
                                      outside_value = 0.01) {
  template_values <- terra::values(template, mat = FALSE)
  raw_values <- terra::values(raw_surface, mat = FALSE)
  if (sum(is.finite(template_values) & !is.finite(raw_values)) != 0L) {
    stop(
      "The historical-equivalent fixture has missing surface values inside finite target cells.",
      call. = FALSE
    )
  }

  surface <- traverse::traverse_surface(
    raw_surface,
    domain = domain,
    type = "conductance",
    transform = transform,
    minimum = minimum,
    outside_value = outside_value
  )

  target_on_computational <- terra::extend(
    domain$target_mask,
    domain$computational_template,
    fill = NA
  )
  surface_values <- terra::values(surface$raster, mat = FALSE)
  target_values <- terra::values(target_on_computational, mat = FALSE)
  fill_idx <- is.na(surface_values) & is.na(target_values)
  surface_values[fill_idx] <- outside_value
  terra::values(surface$raster) <- surface_values
  surface
}

raster_data <- function(x) {
  values <- terra::as.data.frame(x, xy = TRUE, na.rm = FALSE)
  data.frame(
    x = values[[1L]],
    y = values[[2L]],
    value = values[[3L]]
  )
}

summary_values <- function(x) {
  values <- terra::values(x, mat = FALSE)
  finite <- values[is.finite(values)]
  quantiles <- stats::quantile(
    finite,
    probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    names = FALSE
  )
  c(
    minimum = min(finite),
    q05 = quantiles[1L],
    q25 = quantiles[2L],
    median = quantiles[3L],
    q75 = quantiles[4L],
    q95 = quantiles[5L],
    mean = mean(finite),
    maximum = max(finite),
    finite_cells = length(finite),
    positive_cells = sum(finite > 0)
  )
}

compare_rasters <- function(legacy, candidate, absolute_tolerance,
                            correlation_tolerance,
                            normalized_rmse_tolerance) {
  legacy_values <- terra::values(legacy, mat = FALSE)
  candidate_values <- terra::values(candidate, mat = FALSE)
  both_finite <- is.finite(legacy_values) & is.finite(candidate_values)
  if (!all(both_finite == is.finite(legacy_values)) ||
      !all(both_finite == is.finite(candidate_values))) {
    stop("Raster comparison requires identical finite-cell masks.", call. = FALSE)
  }
  legacy_finite <- legacy_values[both_finite]
  candidate_finite <- candidate_values[both_finite]
  difference <- candidate_finite - legacy_finite
  reference_range <- diff(range(legacy_finite))
  stable_scale <- max(abs(legacy_finite))
  stable_reference <- abs(legacy_finite) > max(
    .Machine$double.eps,
    stable_scale * 1e-12
  )

  metrics <- c(
    pearson = stats::cor(legacy_finite, candidate_finite, method = "pearson"),
    spearman = stats::cor(legacy_finite, candidate_finite, method = "spearman"),
    rmse = sqrt(mean(difference^2)),
    mae = mean(abs(difference)),
    maximum_absolute_difference = max(abs(difference)),
    mean_absolute_difference = mean(abs(difference)),
    normalized_rmse = sqrt(mean(difference^2)) / reference_range,
    relative_rmse = sqrt(mean((difference[stable_reference] /
      legacy_finite[stable_reference])^2)),
    relative_maximum_difference = max(abs(difference[stable_reference] /
      legacy_finite[stable_reference]))
  )
  within <- abs(difference) <= absolute_tolerance
  metrics <- c(
    metrics,
    cells_within_tolerance = sum(within),
    percent_cells_within_tolerance = 100 * mean(within)
  )

  geometry_identical <- isTRUE(terra::compareGeom(
    legacy, candidate, stopOnError = FALSE
  ))
  list(
    metrics = metrics,
    geometry_identical = geometry_identical,
    finite_mask_identical = identical(
      is.finite(legacy_values),
      is.finite(candidate_values)
    ),
    acceptance = geometry_identical &&
      identical(
        is.finite(legacy_values),
        is.finite(candidate_values)
      ) &&
      isTRUE(metrics[["pearson"]] >= 1 - correlation_tolerance) &&
      isTRUE(metrics[["normalized_rmse"]] <= normalized_rmse_tolerance) &&
      isTRUE(metrics[["percent_cells_within_tolerance"]] == 100)
  )
}

make_diagnostic_plots <- function(validation) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Diagnostic figures require the ggplot2 and patchwork packages.",
      call. = FALSE
    )
  }
  ggplot2 <- asNamespace("ggplot2")
  target_data <- validation$target_data
  legacy_data <- validation$legacy_data
  traverse_data <- validation$traverse_data
  difference_data <- validation$difference_data
  terminals <- validation$terminals
  target_extent <- validation$target_extent
  computational_extent <- validation$computational_extent
  passage_limits <- range(
    c(legacy_data$value, traverse_data$value),
    na.rm = TRUE
  )
  difference_limit <- max(abs(difference_data$value), na.rm = TRUE)
  if (!is.finite(difference_limit) || difference_limit == 0) {
    difference_limit <- 1
  }
  map_theme <- ggplot2$theme_minimal(base_size = 10) +
    ggplot2$theme(
      panel.grid = ggplot2$element_blank(),
      axis.title = ggplot2$element_blank(),
      plot.title = ggplot2$element_text(face = "bold", size = 11),
      legend.position = "right"
    )
  passage_scale <- ggplot2$scale_fill_gradient(
    low = "#F7F4E9",
    high = "#063B4C",
    limits = passage_limits,
    na.value = "transparent",
    name = "Passage"
  )

  terminal_plot_data <- rbind(
    data.frame(x = terra::crds(validation$sources)[, 1L],
               y = terra::crds(validation$sources)[, 2L],
               terminal = "Source"),
    data.frame(x = terra::crds(validation$targets)[, 1L],
               y = terra::crds(validation$targets)[, 2L],
               terminal = "Target")
  )
  p_terminals <- ggplot2$ggplot() +
    ggplot2$geom_rect(
      data = computational_extent,
      ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, color = "#5B6668", linewidth = 0.55
    ) +
    ggplot2$geom_raster(
      data = target_data,
      ggplot2$aes(x = x, y = y),
      fill = "#DCE9E7"
    ) +
    ggplot2$geom_rect(
      data = target_extent,
      ggplot2$aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, color = "#A84B2A", linewidth = 0.75
    ) +
    ggplot2$geom_point(
      data = terminal_plot_data,
      ggplot2$aes(x = x, y = y, shape = terminal, color = terminal),
      size = 2.5
    ) +
    ggplot2$scale_shape_manual(values = c(Source = 17, Target = 15)) +
    ggplot2$scale_color_manual(values = c(Source = "#063B4C", Target = "#A84B2A")) +
    ggplot2$coord_equal(expand = FALSE) +
    ggplot2$labs(title = "Shared terminal geometry", shape = NULL, color = NULL) +
    map_theme

  p_legacy <- ggplot2$ggplot(legacy_data, ggplot2$aes(x = x, y = y, fill = value)) +
    ggplot2$geom_raster(na.rm = FALSE) + passage_scale +
    ggplot2$coord_equal(expand = FALSE) +
    ggplot2$labs(title = "A. Historical gdistance", fill = "Passage") + map_theme
  p_traverse <- ggplot2$ggplot(traverse_data, ggplot2$aes(x = x, y = y, fill = value)) +
    ggplot2$geom_raster(na.rm = FALSE) + passage_scale +
    ggplot2$coord_equal(expand = FALSE) +
    ggplot2$labs(title = "B. traverse", fill = "Passage") + map_theme
  p_difference <- ggplot2$ggplot(difference_data, ggplot2$aes(x = x, y = y, fill = value)) +
    ggplot2$geom_raster(na.rm = FALSE) +
    ggplot2$scale_fill_gradient2(
      low = "#557A98", mid = "#FBFAF7", high = "#A45A45",
      midpoint = 0, limits = c(-difference_limit, difference_limit),
      na.value = "transparent", name = "Difference"
    ) +
    ggplot2$coord_equal(expand = FALSE) +
    ggplot2$labs(title = "Absolute agreement: signed difference") + map_theme

  scatter_data <- data.frame(
    legacy = validation$legacy_finite,
    traverse = validation$traverse_finite
  )
  scatter_limits <- range(c(scatter_data$legacy, scatter_data$traverse), na.rm = TRUE)
  p_scatter <- ggplot2$ggplot(scatter_data, ggplot2$aes(x = legacy, y = traverse)) +
    ggplot2$geom_point(alpha = 0.22, size = 0.7, color = "#063B4C") +
    ggplot2$geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#A45A45") +
    ggplot2$coord_equal(xlim = scatter_limits, ylim = scatter_limits, expand = FALSE) +
    ggplot2$annotate(
      "text", x = scatter_limits[1L], y = scatter_limits[2L],
      hjust = 0, vjust = 1, size = 3.3,
      label = sprintf(
        "Pearson r = %.12g\nRMSE = %.6g",
        validation$comparison$metrics[["pearson"]],
        validation$comparison$metrics[["rmse"]]
      )
    ) +
    ggplot2$labs(
      title = "Cell-by-cell passage comparison",
      x = "Historical passage", y = "traverse passage"
    ) + map_theme

  list(
    terminals = p_terminals,
    passage = patchwork::wrap_plots(p_legacy, p_traverse, ncol = 2),
    difference = p_difference,
    scatter = p_scatter
  )
}

run_historical_validation <- function(
    root = project_root,
    n_nodes = 2L,
    absolute_tolerance = 1e-10,
    correlation_tolerance = 1e-10,
    normalized_rmse_tolerance = 1e-10,
    output_dir = file.path(root, "validation-output", "historical-equivalence"),
    write_figures = TRUE
) {
  required <- c("terra", "raster", "sp", "gdistance")
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing validation packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  template <- terra::rast(file.path(root, "assets", "template.tif"))
  raw_surface <- terra::rast(file.path(root, "assets", "surface.tif"))
  domain <- traverse::traverse_domain(template, buffer = 150, buffer_units = "map")
  surface <- repair_historical_surface(template, raw_surface, domain)

  sources <- traverse::traverse_nodes(
    domain, side = "south", distance = 100, n = n_nodes, placement = "even"
  )
  targets <- traverse::traverse_nodes(
    domain, side = "north", distance = 100, n = n_nodes, placement = "even"
  )
  pairs <- traverse::traverse_pairs(sources, targets)
  settings <- list(
    directions = 16L,
    theta = 1e-100,
    correction = "c",
    scale_correction = TRUE,
    aggregation = "mean",
    transform = "x / 100",
    minimum = 0.01,
    outside_value = 0.01
  )

  legacy_time <- system.time({
    legacy <- historical_reference(
      surface, domain, sources, targets, pairs,
      directions = settings$directions,
      theta = settings$theta,
      correction = settings$correction,
      scale_correction = settings$scale_correction
    )
  })
  traverse_time <- system.time({
    modern <- traverse::traverse(
      surface, domain, sources = sources, targets = targets, pairs = pairs,
      directions = settings$directions,
      theta = settings$theta,
      correction = settings$correction,
      scale_correction = settings$scale_correction,
      aggregation = settings$aggregation
    )
  })

  legacy_values <- terra::values(legacy$target_passage, mat = FALSE)
  traverse_values <- terra::values(modern$target_passage, mat = FALSE)
  both_finite <- is.finite(legacy_values) & is.finite(traverse_values)
  difference <- modern$target_passage - legacy$target_passage
  comparison <- compare_rasters(
    legacy$target_passage,
    modern$target_passage,
    absolute_tolerance = absolute_tolerance,
    correlation_tolerance = correlation_tolerance,
    normalized_rmse_tolerance = normalized_rmse_tolerance
  )

  target_data <- raster_data(template)
  target_data$value[!is.finite(target_data$value)] <- NA_real_
  legacy_data <- raster_data(legacy$target_passage)
  traverse_data <- raster_data(modern$target_passage)
  difference_data <- raster_data(difference)
  terminals <- data.frame(
    terminal = c(rep("Source", nrow(sources)), rep("Target", nrow(targets))),
    x = c(terra::crds(sources)[, 1L], terra::crds(targets)[, 1L]),
    y = c(terra::crds(sources)[, 2L], terra::crds(targets)[, 2L])
  )
  target_extent <- data.frame(
    xmin = terra::xmin(template), xmax = terra::xmax(template),
    ymin = terra::ymin(template), ymax = terra::ymax(template)
  )
  computational_extent <- data.frame(
    xmin = terra::xmin(domain$computational_template),
    xmax = terra::xmax(domain$computational_template),
    ymin = terra::ymin(domain$computational_template),
    ymax = terra::ymax(domain$computational_template)
  )

  settings_table <- data.frame(
    Parameter = c(
      "Target raster dimensions", "Resolution", "Target cells",
      "Source nodes", "Target nodes", "Source-target pairs", "Directions",
      "Theta", "Correction", "Scale correction", "Aggregation",
      "Conductance transform", "Minimum conductance", "Outside value"
    ),
    Value = c(
      sprintf("%d x %d", terra::nrow(template), terra::ncol(template)),
      paste(signif(terra::res(template), 8), collapse = " x "),
      format(sum(is.finite(terra::values(template, mat = FALSE))), big.mark = ","),
      nrow(sources), nrow(targets), nrow(pairs), settings$directions,
      format(settings$theta, scientific = TRUE), settings$correction,
      settings$scale_correction, settings$aggregation, settings$transform,
      settings$minimum, settings$outside_value
    ),
    stringsAsFactors = FALSE
  )
  summary_table <- data.frame(
    Metric = names(summary_values(legacy$target_passage)),
    Historical = unname(summary_values(legacy$target_passage)),
    Traverse = unname(summary_values(modern$target_passage)),
    row.names = NULL
  )
  metrics_table <- data.frame(
    Metric = names(comparison$metrics),
    Value = unname(comparison$metrics),
    row.names = NULL
  )
  validation <- list(
    root = root,
    template = template,
    raw_surface = raw_surface,
    surface = surface,
    domain = domain,
    sources = sources,
    targets = targets,
    pairs = pairs,
    settings = settings,
    settings_table = settings_table,
    legacy = legacy,
    modern = modern,
    legacy_time = legacy_time,
    traverse_time = traverse_time,
    legacy_finite = legacy_values[both_finite],
    traverse_finite = traverse_values[both_finite],
    target_data = target_data,
    legacy_data = legacy_data,
    traverse_data = traverse_data,
    difference_data = difference_data,
    terminals = terminals,
    target_extent = target_extent,
    computational_extent = computational_extent,
    comparison = comparison,
    summary_table = summary_table,
    metrics_table = metrics_table,
    tolerances = list(
      absolute = absolute_tolerance,
      correlation = correlation_tolerance,
      normalized_rmse = normalized_rmse_tolerance
    )
  )
  validation$plots <- make_diagnostic_plots(validation)

  if (isTRUE(write_figures)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(output_dir, "terminals.png"), validation$plots$terminals,
                    width = 8, height = 5.5, dpi = 160)
    ggplot2::ggsave(file.path(output_dir, "passage-comparison.png"), validation$plots$passage,
                    width = 11, height = 5.5, dpi = 160)
    ggplot2::ggsave(file.path(output_dir, "difference.png"), validation$plots$difference,
                    width = 8, height = 5.5, dpi = 160)
    ggplot2::ggsave(file.path(output_dir, "scatter.png"), validation$plots$scatter,
                    width = 7, height = 6, dpi = 160)
    validation$figure_dir <- normalizePath(output_dir, mustWork = TRUE)
  }
  if (!isTRUE(comparison$acceptance)) {
    stop(
      "Historical equivalence acceptance criteria failed. See comparison$metrics.",
      call. = FALSE
    )
  }
  validation
}

print_historical_validation <- function(validation) {
  cat("Historical NWS equivalence validation\n")
  cat("  source nodes:", nrow(validation$sources), "\n")
  cat("  target nodes:", nrow(validation$targets), "\n")
  cat("  source-target pairs:", nrow(validation$pairs), "\n")
  cat("  legacy elapsed seconds:", unname(validation$legacy_time[["elapsed"]]), "\n")
  cat("  traverse elapsed seconds:", unname(validation$traverse_time[["elapsed"]]), "\n")
  print(validation$metrics_table, row.names = FALSE)
  cat("  geometry identical:", validation$comparison$geometry_identical, "\n")
  cat("  finite masks identical:", validation$comparison$finite_mask_identical, "\n")
  cat("  acceptance criteria passed:", validation$comparison$acceptance, "\n")
  if (!is.null(validation$figure_dir)) {
    cat("  diagnostic figures:", validation$figure_dir, "\n")
  }
  invisible(validation)
}

# Set `options(traverse.skip_historical_validation = TRUE)` before sourcing
# this file from another analysis document. A direct non-interactive Rscript
# invocation therefore runs the validation, while the functions remain
# sourceable for the Quarto demonstration.
if (!isTRUE(getOption("traverse.skip_historical_validation", FALSE)) &&
    !interactive()) {
  load_traverse(project_root)
  validation <- run_historical_validation(root = project_root)
  print_historical_validation(validation)
}
