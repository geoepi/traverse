Validating traverse against the historical workflow
================

- <a href="#why-validate-equivalence"
  id="toc-why-validate-equivalence">Why validate equivalence?</a>
- <a href="#shared-landscape-and-parameters"
  id="toc-shared-landscape-and-parameters">Shared landscape and
  parameters</a>
- <a href="#shared-terminals" id="toc-shared-terminals">Shared
  terminals</a>
- <a href="#historical-implementation"
  id="toc-historical-implementation">Historical implementation</a>
- <a href="#traverse-implementation"
  id="toc-traverse-implementation"><code>traverse</code>
  implementation</a>
- <a href="#spatial-comparison" id="toc-spatial-comparison">Spatial
  comparison</a>
- <a href="#numerical-equivalence"
  id="toc-numerical-equivalence">Numerical equivalence</a>
- <a href="#cell-by-cell-comparison"
  id="toc-cell-by-cell-comparison">Cell-by-cell comparison</a>
- <a href="#interpretation" id="toc-interpretation">Interpretation</a>

`traverse` was extracted from a working landscape-connectivity workflow.
This note validates the documented historical configuration against that
original implementation before the package is generalized to broader
applications.

## Why validate equivalence?

This is software/model equivalence validation, not a claim that the
historical workflow was scientifically optimal. The comparison holds the
landscape, conductance values, terminal coordinates, pairings,
transition settings, passage settings, and target mask constant. The
historical implementation is run directly with `raster`, `sp`, and
`gdistance`; the new implementation keeps those classes behind the
`traverse` API.

## Shared landscape and parameters

The repaired repository assets provide a controlled proxy for the
historical application. Finite cells in `template.tif` define the target
study area, and the surface has no missing values in those finite target
cells. Missing cells outside the target are assigned the same low
conductance in both workflows.

<details open>
<summary>Show code</summary>

``` r
knitr::kable(
  validation$settings_table,
  caption = "Inputs and settings held constant between implementations."
)
```

</details>

| Parameter                | Value                 |
|:-------------------------|:----------------------|
| Target raster dimensions | 202 x 293             |
| Resolution               | 24.995019 x 24.943656 |
| Target cells             | 16,759                |
| Source nodes             | 2                     |
| Target nodes             | 2                     |
| Source-target pairs      | 4                     |
| Directions               | 16                    |
| Theta                    | 1e-100                |
| Correction               | c                     |
| Scale correction         | TRUE                  |
| Aggregation              | mean                  |
| Conductance transform    | x / 100               |
| Minimum conductance      | 0.01                  |
| Outside value            | 0.01                  |

Inputs and settings held constant between implementations.

The target raster is extended by a 150 map-unit computational buffer.
The surface is transformed as `x / 100`, given a finite minimum
conductance of `0.01`, and assigned `0.01` in the non-target
computational area. Two source and two target terminals create four
source-target passage pairs. The validation function accepts a larger
`n_nodes` value when a full 3 x 3 run is practical on the host system.

## Shared terminals

The terminal coordinates are generated once with `traverse_nodes()` and
passed unchanged to both implementations. This is important:
independently regenerating “equivalent” terminals could introduce a
geometry difference that would be mistaken for a model difference.

<details open>
<summary>Show code</summary>

``` r
validation$plots$terminals
```

</details>

![The same source and target coordinates are supplied to the historical
and traverse implementations. The red rectangle is the target extent;
the grey rectangle is the computational
extent.](historical-equivalence_files/figure-commonmark/terminals-figure-1.png)

## Historical implementation

The reference path converts the prepared conductance surface to a legacy
`RasterLayer`, creates a 16-neighbour transition layer, applies
`geoCorrection(type = "c", multpl = FALSE, scl = TRUE)`, calculates
passage for each of the four pairs with `theta = 1e-100`, averages those
pairwise rasters, and masks the aggregate back to the target.

<details open>
<summary>Show code</summary>

``` r
legacy_raster <- raster::raster(surface$raster)
raster::crs(legacy_raster) <- terra::crs(surface$raster, proj = TRUE)
legacy_transition <- gdistance::transition(
  legacy_raster, transitionFunction = mean, directions = 16
)
legacy_transition <- gdistance::geoCorrection(
  legacy_transition, type = "c", multpl = FALSE, scl = TRUE
)

source_sp <- sp::SpatialPoints(
  terra::crds(sources),
  proj4string = sp::CRS(terra::crs(sources, proj = TRUE))
)
target_sp <- sp::SpatialPoints(
  terra::crds(targets),
  proj4string = sp::CRS(terra::crs(targets, proj = TRUE))
)
legacy_pairwise <- lapply(seq_len(nrow(pairs)), function(i) {
  gdistance::passage(
    legacy_transition,
    origin = source_sp[pairs$source_index[i], , drop = FALSE],
    goal = target_sp[pairs$target_index[i], , drop = FALSE],
    theta = 1e-100
  )
})
legacy_mean <- raster::calc(
  raster::stack(legacy_pairwise),
  fun = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
)
legacy_target <- raster::mask(
  raster::crop(legacy_mean, raster::raster(domain$target_template)),
  raster::raster(domain$target_mask)
)
```

</details>

## `traverse` implementation

The equivalent `traverse` call receives the same prepared surface,
domain, terminals, pair table, and model settings. The package converts
to legacy backend classes only inside its transition and passage
adapters.

<details open>
<summary>Show code</summary>

``` r
traverse_result <- traverse(
  surface,
  domain,
  sources = sources,
  targets = targets,
  pairs = pairs,
  directions = 16,
  theta = 1e-100,
  correction = "c",
  scale_correction = TRUE,
  aggregation = "mean"
)
```

</details>

The development validation script runs both code paths and fails if the
core acceptance criteria are not met. It also records elapsed time for
each path; runtime is informational here and is not used to claim an
optimization.

## Spatial comparison

The historical and `traverse` maps use identical scale limits. The
difference map is signed (`traverse - historical`) and uses a restrained
diverging scale centered at zero.

<details open>
<summary>Show code</summary>

``` r
validation$plots$passage
```

</details>

![Historical and traverse target-domain passage surfaces shown on
identical colour
scales.](historical-equivalence_files/figure-commonmark/passage-comparison-1.png)

<details open>
<summary>Show code</summary>

``` r
validation$plots$difference
```

</details>

![Signed cell-by-cell difference between the traverse and historical
target-domain passage
surfaces.](historical-equivalence_files/figure-commonmark/difference-map-1.png)

## Numerical equivalence

Raster geometry includes dimensions, resolution, extent, CRS, and the
finite target mask. Values are compared only where both outputs are
finite, with relative metrics excluding reference values effectively
equal to zero.

<details open>
<summary>Show code</summary>

``` r
knitr::kable(
  validation$summary_table,
  digits = 8,
  caption = "Target-domain raster summaries."
)
```

</details>

| Metric         |   Historical |     Traverse |
|:---------------|-------------:|-------------:|
| minimum        | 3.354000e-05 | 3.354000e-05 |
| q05            | 3.206800e-04 | 3.206800e-04 |
| q25            | 1.245080e-03 | 1.245080e-03 |
| median         | 3.904660e-03 | 3.904660e-03 |
| q75            | 1.076351e-02 | 1.076351e-02 |
| q95            | 2.822189e-02 | 2.822189e-02 |
| mean           | 8.131840e-03 | 8.131840e-03 |
| maximum        | 4.729858e-01 | 4.729858e-01 |
| finite_cells   | 1.675900e+04 | 1.675900e+04 |
| positive_cells | 1.675900e+04 | 1.675900e+04 |

Target-domain raster summaries.

<details open>
<summary>Show code</summary>

``` r
metric_names <- c(
  "pearson", "spearman", "rmse", "mae",
  "maximum_absolute_difference", "mean_absolute_difference",
  "normalized_rmse", "relative_rmse", "relative_maximum_difference",
  "percent_cells_within_tolerance"
)
metric_table <- validation$metrics_table[
  match(metric_names, validation$metrics_table$Metric), , drop = FALSE
]
knitr::kable(
  metric_table,
  digits = 12,
  caption = "Cell-by-cell equivalence metrics."
)
```

</details>

|     | Metric                         | Value |
|:----|:-------------------------------|------:|
| 1   | pearson                        |     1 |
| 2   | spearman                       |     1 |
| 3   | rmse                           |     0 |
| 4   | mae                            |     0 |
| 5   | maximum_absolute_difference    |     0 |
| 6   | mean_absolute_difference       |     0 |
| 7   | normalized_rmse                |     0 |
| 8   | relative_rmse                  |     0 |
| 9   | relative_maximum_difference    |     0 |
| 11  | percent_cells_within_tolerance |   100 |

Cell-by-cell equivalence metrics.

<details open>
<summary>Show code</summary>

``` r
geometry_table <- data.frame(
  Check = c("Geometry identical", "Finite target mask identical", "Acceptance criteria"),
  Result = c(
    validation$comparison$geometry_identical,
    validation$comparison$finite_mask_identical,
    validation$comparison$acceptance
  )
)
knitr::kable(geometry_table, caption = "Structural equivalence checks.")
```

</details>

| Check                        | Result |
|:-----------------------------|:-------|
| Geometry identical           | TRUE   |
| Finite target mask identical | TRUE   |
| Acceptance criteria          | TRUE   |

Structural equivalence checks.

## Cell-by-cell comparison

The scatter plot provides a direct check of the one-to-one relationship
between the two output rasters. The dashed line is the identity line.

<details open>
<summary>Show code</summary>

``` r
validation$plots$scatter
```

</details>

![Historical passage values versus traverse passage values for every
finite target
cell.](historical-equivalence_files/figure-commonmark/scatter-1.png)

## Interpretation

<details open>
<summary>Show code</summary>

``` r
interpretation <- if (isTRUE(validation$comparison$acceptance)) {
  paste(
    "For this documented historical configuration, traverse reproduces the",
    "reference gdistance passage workflow within the reported numerical",
    "tolerance. Geometry and finite target masks are identical, and the",
    "cell-by-cell residuals are within tolerance. This conclusion is limited",
    "to the tested conductance preparation, terminal geometry, pair weighting,",
    "and gdistance settings."
  )
} else {
  paste(
    "The configured equivalence criteria were not met. Inspect the numerical",
    "metrics and signed difference map before changing package behavior."
  )
}
interpretation
```

</details>

    [1] "For this documented historical configuration, traverse reproduces the reference gdistance passage workflow within the reported numerical tolerance. Geometry and finite target masks are identical, and the cell-by-cell residuals are within tolerance. This conclusion is limited to the tested conductance preparation, terminal geometry, pair weighting, and gdistance settings."

The historical NWS application motivates this validation, but no
NWS-specific transformation is embedded in the package API. For the
documented historical configuration, the package architecture has now
been checked against the legacy implementation while retaining explicit
control over conductance, terminal placement, transition settings,
geographic correction, and pair aggregation.
