# traverse

`traverse` estimates spatial passage across gridded landscapes between user-defined source and destination sets. It separates the scientific target domain from a buffered computational domain, allowing terminal nodes to be placed outside the area of interpretation and removed cleanly from final outputs.

The package uses `terra::SpatRaster` and `terra::SpatVector` in its user-facing API. It currently uses `gdistance` as an internal randomized-shortest-path backend; legacy `raster`, `sp`, and `gdistance` classes are not required as user inputs.

## Core workflow

```r
library(traverse)

template <- terra::rast(
  nrows = 20, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 20,
  crs = "EPSG:3857"
)
terra::values(template) <- 1
template[10, 15] <- NA # an excluded cell remains a hole in the target mask

domain <- traverse_domain(template, buffer = 5, buffer_units = "cells")
surface <- traverse_surface(
  template, domain = domain, type = "conductance", outside_value = 1
)
sources <- traverse_nodes(domain, side = "south", distance = 1, n = 8)
targets <- traverse_nodes(domain, side = "north", distance = 1, n = 8)

result <- traverse(
  surface, domain, sources = sources, targets = targets,
  directions = 16, theta = 1, correction = "c", aggregation = "mean"
)

plot(result)                         # clean target-domain interpretation
result$computational_passage         # retained for diagnostics
```

## Important concepts

* **Target domain:** the area for which scientific interpretation and final output are desired. Finite cells in the template belong to this domain; `NA` cells are excluded.
* **Computational domain:** a larger raster containing the target plus an external buffer used for movement calculations and terminal placement.
* **Terminal zone:** locations outside the target where source and destination nodes are placed. Terminals can also be supplied directly as point `SpatVector` or `sf` objects.
* **Conductance and resistance:** higher conductance means easier movement; higher resistance means greater movement cost. Resistance is converted internally as `1 / resistance` only at the `gdistance` boundary.
* **Passage and connectivity:** the package reports passage intensity (or relative passage), not a probability unless a downstream analysis supplies an appropriate probabilistic interpretation.
* **Terminal artifacts:** path origins and destinations can have high passage values. External terminals and final masking keep those artifacts outside the target area; this is a methodological part of the workflow, not merely raster cleanup.
* **`theta`:** a randomized-shortest-path model parameter. It changes the balance between shortest-path-like and exploratory behavior, so it is explicit and should be justified for each application.
* **Pair aggregation:** `aggregation = "mean"` averages pairwise passage surfaces; `aggregation = "sum"` adds them. Equal weighting is the initial default and is not asserted to be universally appropriate.

`traverse` does not choose universal conductance transformations, environmental weights, resistance floors, or scientifically optimal `theta` and geographic-correction values. Those decisions remain configurable for the study design.

See [`docs/architecture.md`](docs/architecture.md) for the package data flow, object structure, backend boundary, and open methodological questions.

## Development

The package uses roxygen2 documentation and testthat tests. With R and the package dependencies installed:

```r
devtools::document()
devtools::test()
devtools::check()
```
