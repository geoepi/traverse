Traversing a gridded landscape
================

- <a href="#landscape-inputs" id="toc-landscape-inputs">Landscape
  inputs</a>
- <a href="#define-the-computational-domain"
  id="toc-define-the-computational-domain">Define the computational
  domain</a>
- <a href="#prepare-conductance" id="toc-prepare-conductance">Prepare
  conductance</a>
- <a href="#external-terminals-and-geometry"
  id="toc-external-terminals-and-geometry">External terminals and
  geometry</a>
- <a href="#passage-across-the-computational-domain"
  id="toc-passage-across-the-computational-domain">Passage across the
  computational domain</a>
- <a href="#target-domain-passage"
  id="toc-target-domain-passage">Target-domain passage</a>
- <a href="#the-traverse-workflow" id="toc-the-traverse-workflow">The
  traverse workflow</a>

`traverse` estimates relative passage across a gridded landscape between
source and destination sets. This compact demonstration uses the
package’s example rasters to show the full workflow: define a scientific
target domain, expand it for computation, place terminals outside the
target, and mask the result back to the area being interpreted.

## Landscape inputs

The demonstration separates two ideas. The **target domain** is the set
of cells for which passage will ultimately be interpreted. The
**movement surface** describes relative suitability or preference across
that landscape; higher values indicate easier or more preferred
movement. In `template.tif`, finite `0` values mark included cells and
`NA` marks cells outside the target. The zero is a mask value only; it
is not a movement cost. In `surface.tif`, larger continuous values
indicate greater suitability.

<details open>
<summary>Show code</summary>

``` r
p_template | p_surface
```

</details>

![The target mask defines the scientific area of interest; the
suitability raster supplies relative movement
preference.](asset-demo_files/figure-commonmark/input-maps-1.png)

`template.tif` and `surface.tif` share a projected grid whose map units
are kilometres in this example. The package does not assume that map
units are kilometres, however; the analyst supplies that interpretation
through the coordinate reference system and chosen distances.

<details open>
<summary>Show code</summary>

``` r
knitr::kable(input_summary)
```

</details>

| Property          | Value                                |
|:------------------|:-------------------------------------|
| Raster dimensions | 202 x 293 cells                      |
| Resolution        | 24.995 x 24.9437                     |
| Target cells      | 16,759                               |
| Surface range     | 0 to 927.268                         |
| Map units         | kilometres in this projected example |

## Define the computational domain

Passage values are naturally high near origins and destinations. Moving
terminals outside the scientific target area allows those terminal
effects to occur in a computational buffer rather than obscure the
landscape being interpreted. Here the buffer is 150 map units (150 km
for this example), while the final result will still be reported only on
the target mask.

<details open>
<summary>Show code</summary>

``` r
domain <- traverse_domain(
  template,
  buffer = 150,
  buffer_units = "map"
)
```

</details>

## Prepare conductance

The input surface is a suitability/preference surface, so this example
uses it as conductance: larger values mean easier movement. The
transformation and minimum below are demonstration settings, not
universal defaults. `traverse` intentionally leaves surface
transformation and minimum conductance to the analyst.

<details open>
<summary>Show code</summary>

``` r
surface <- traverse_surface(
  surface_raw,
  domain = domain,
  type = "conductance",
  transform = function(x) x / 100,
  minimum = 0.01,
  outside_value = 0.01
)

# The repaired asset has finite values across the intended target cells. Any
# remaining NA cells in the original rectangular raster are outside the target
# mask; assign them the same explicitly chosen low conductance so terminals
# can reach the landscape. The extended outer buffer is already 0.01.
target_on_computational <- terra::extend(
  domain$target_mask,
  domain$computational_template,
  fill = NA
)
surface_values <- terra::values(surface$raster, mat = FALSE)
target_values <- terra::values(target_on_computational, mat = FALSE)
surface_xy <- terra::xyFromCell(surface$raster, seq_len(terra::ncell(surface$raster)))
template_ext <- terra::ext(template)
inside_original <- surface_xy[, 1] >= template_ext$xmin &
  surface_xy[, 1] <= template_ext$xmax &
  surface_xy[, 2] >= template_ext$ymin &
  surface_xy[, 2] <= template_ext$ymax
fill_idx <- is.na(surface_values) & is.na(target_values) & inside_original
surface_values[fill_idx] <- 0.01
terra::values(surface$raster) <- surface_values
```

</details>

This distinguishes repaired target-domain data from intentionally low
conductance in the non-target computational area. The latter is a
boundary condition for this demonstration, not a claim that movement
outside the target is biologically uniform.

## External terminals and geometry

One source and one target keep this full-resolution example bounded:
there is one source-target passage calculation. Setting `n_nodes = 3`
would produce `3 x 3 = 9` pairwise calculations.

<details open>
<summary>Show code</summary>

``` r
n_nodes <- 1

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

computational_extent <- data.frame(
  xmin = terra::xmin(domain$computational_template),
  xmax = terra::xmax(domain$computational_template),
  ymin = terra::ymin(domain$computational_template),
  ymax = terra::ymax(domain$computational_template)
)
target_boundary <- as.data.frame(
  terra::geom(terra::as.polygons(domain$target_mask, dissolve = TRUE))
)
target_boundary$group <- interaction(
  target_boundary$geom, target_boundary$part, target_boundary$hole,
  drop = TRUE
)

p_geometry <- ggplot() +
  geom_rect(
    data = computational_extent,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#F5F6F3", color = "#5B6668", linewidth = 0.55
  ) +
  geom_raster(data = target_cells, aes(x = x, y = y), fill = "#DCE9E7") +
  geom_path(
    data = target_boundary,
    aes(x = x, y = y, group = group),
    color = "#A84B2A", linewidth = 0.8
  ) +
  geom_point(
    data = terminal_data,
    aes(x = x, y = y, shape = terminal, color = terminal),
    size = 2.7
  ) +
  scale_shape_manual(values = c(Source = 17, Target = 15)) +
  scale_color_manual(values = c(Source = "#063B4C", Target = "#A84B2A")) +
  coord_equal(expand = FALSE) +
  labs(title = "Computational geometry", shape = NULL, color = NULL) +
  map_theme
```

</details>
<details open>
<summary>Show code</summary>

``` r
p_geometry
```

</details>

![Source and target terminals are placed outside the target study area
but remain inside the buffered computational
domain.](asset-demo_files/figure-commonmark/geometry-figure-1.png)

## Passage across the computational domain

Randomized-shortest-path passage integrates many possible routes between
the terminal sets rather than returning one deterministic least-cost
line. `theta` controls that balance between shortest-path-like and
exploratory behaviour. We use the historical near-random-walk
configuration here; it is a model parameter, not a universal
calibration. Larger values can concentrate passage more strongly and may
be numerically difficult on large landscapes.

<details open>
<summary>Show code</summary>

``` r
result <- traverse(
  surface,
  domain,
  sources = sources,
  targets = targets,
  directions = 16,
  theta = 1e-100,
  correction = "c",
  scale_correction = TRUE,
  aggregation = "mean"
)

stopifnot(
  sum(is.finite(terra::values(result$target_passage))) ==
    sum(is.finite(terra::values(template)))
)

computational_passage_data <- raster_data(result$computational_passage)
target_passage_data <- raster_data(result$target_passage)

passage_scale <- scale_fill_gradient(
  low = "#F7F4E9",
  high = "#063B4C",
  trans = "sqrt",
  na.value = "transparent",
  name = "Passage intensity"
)

p_computational <- ggplot(computational_passage_data, aes(x = x, y = y, fill = value)) +
  geom_raster(na.rm = FALSE) +
  passage_scale +
  geom_path(data = target_boundary, aes(x = x, y = y, group = group), inherit.aes = FALSE,
            color = "#A84B2A", linewidth = 0.75) +
  geom_point(data = terminal_data,
             aes(x = x, y = y, shape = terminal, color = terminal),
             inherit.aes = FALSE, size = 2.2) +
  scale_shape_manual(values = c(Source = 17, Target = 15)) +
  scale_color_manual(values = c(Source = "#063B4C", Target = "#A84B2A")) +
  coord_equal(expand = FALSE) +
  labs(title = "C. Computational passage", shape = NULL, color = NULL) +
  map_theme

p_target <- ggplot(target_passage_data, aes(x = x, y = y, fill = value)) +
  geom_raster(na.rm = FALSE) +
  passage_scale +
  coord_equal(expand = FALSE) +
  labs(title = "D. Target-domain passage") +
  map_theme
```

</details>
<details open>
<summary>Show code</summary>

``` r
p_computational
```

</details>

![Passage is calculated across the expanded landscape.
Terminal-associated concentrations remain outside the area intended for
interpretation. The colour scale uses a square-root transformation for
visualisation
only.](asset-demo_files/figure-commonmark/computational-passage-figure-1.png)

## Target-domain passage

The final surface is the computational result masked back to the
scientific target domain. Pale cells have lower relative passage and
deep teal cells have higher relative passage; the analytical values are
unchanged by this plotting palette.

<details open>
<summary>Show code</summary>

``` r
p_target
```

</details>

![The computational buffer and external terminals are removed, leaving
passage intensity within the target study domain. The colour scale uses
a square-root transformation for visualisation
only.](asset-demo_files/figure-commonmark/target-passage-figure-1.png)

<details open>
<summary>Show code</summary>

``` r
target_values <- terra::values(result$target_passage, mat = FALSE)
finite_target <- target_values[is.finite(target_values)]
result_summary <- data.frame(
  Metric = c("Source nodes", "Target nodes", "Source-target pairs", "Target cells",
             "Minimum passage", "Median passage", "Mean passage", "Maximum passage"),
  Value = c(
    format(c(nrow(result$sources), nrow(result$targets), nrow(result$pairs),
             length(finite_target)), big.mark = ",", trim = TRUE),
    format(signif(c(min(finite_target), median(finite_target),
                    mean(finite_target), max(finite_target)), 6),
           scientific = TRUE, trim = TRUE)
  )
)
knitr::kable(result_summary)
```

</details>

| Metric              | Value       |
|:--------------------|:------------|
| Source nodes        | 1           |
| Target nodes        | 1           |
| Source-target pairs | 1           |
| Target cells        | 16,759      |
| Minimum passage     | 1.78886e-06 |
| Median passage      | 3.98580e-03 |
| Mean passage        | 1.00988e-02 |
| Maximum passage     | 2.53042e-01 |

## The traverse workflow

<details open>
<summary>Show code</summary>

``` r
(p_surface | p_geometry) / (p_computational | p_target)
```

</details>

![The complete workflow: suitability, buffered terminal geometry,
computational passage, and the masked target-domain
result.](asset-demo_files/figure-commonmark/overview-1.png)

The input suitability raster determines relative movement conductance.
Passage then integrates many possible routes between external source and
target terminals. Keeping those terminals outside the target reduces
boundary artifacts, and masking the final raster restores the scientific
interpretation area. Passage intensity is not a probability of movement;
it is a relative summary of the configured landscape and
randomized-shortest-path model.
