# traverse architecture

## Public API

The package keeps the public surface small:

* `traverse_domain()` defines a target template and a larger computational template.
* `traverse_nodes()` creates cardinal external terminals as a convenience layer.
* `traverse_surface()` validates movement-surface semantics and aligns a surface to a computational domain.
* `traverse_pairs()` creates stable source-target pair tables.
* `traverse_passage()` runs pairwise randomized-shortest-path passage and explicit mean or sum aggregation.
* `traverse_mask()` restores a computational raster to the exact target geometry and mask.
* `traverse()` orchestrates the workflow and returns a `traverse_result`.

The primary spatial classes are `terra::SpatRaster` and `terra::SpatVector`. `sf` point objects are accepted when the optional `sf` dependency is available. `raster`, `sp`, and `gdistance` objects are implementation details of the current backend.

## Object structure

### `traverse_domain`

This S3 object stores the original `target_template`, a binary `target_mask`, `target_extent`, `computational_extent`, a finite `computational_template`, buffer metadata, CRS, resolution, and cell padding. The target template is never replaced by the computational raster. This lets masking preserve internal holes and exact target geometry.

### `traverse_surface`

This S3 object stores an aligned `raster`, semantic `type` (`conductance` or `resistance`), and the optional transformation, minimum, and outside fill metadata. No 0–1 normalization or project-specific multiplication is applied by default.

### `traverse_result`

This S3 object stores `target_passage`, `computational_passage`, `sources`, `targets`, `pairs`, `domain`, `surface`, and explicit `model` parameters. The computational raster is retained so edge behavior and accessibility can be diagnosed rather than hidden.

## Data flow

```text
target template
       |
       v
traverse_domain() ---> target mask + computational template
       |                              |
       |                              +--> traverse_nodes()
       v
traverse_surface() ------------------+
                                      v
                            traverse_pairs()
                                      |
                                      v
                           internal gdistance adapter
                                      |
                                      v
                           pairwise passage surfaces
                                      |
                                      v
                            mean or sum aggregation
                                      |
                                      v
                            traverse_mask() ---> target result
```

The computational domain is intentionally larger than the target domain. Terminal cells can have high values because paths begin or end there. External placement and target masking reduce the risk that this terminal behavior obscures the landscape structure being interpreted.

## Backend isolation

`traverse_transition()` converts the prepared `SpatRaster` to a legacy `raster` object and creates a `gdistance::TransitionLayer` using a configurable neighborhood and geographic correction. Passage converts point coordinates to temporary `sp::SpatialPoints` only at the same boundary. The public functions do not require callers to construct either legacy class.

Resistance is converted as `conductance = 1 / resistance` for the backend. This is explicit in the implementation and documentation, but users should still consider whether it is scientifically appropriate for their application. Zero or non-finite resistance values are not silently repaired.

## Scientific assumptions and configurable decisions

* Finite template cells define the target; `NA` cells are excluded.
* Cardinal terminal generation uses the outermost finite target cell in each row or column. It is a practical first-release approximation, not a general polygon-normal algorithm.
* Node distance is measured in CRS map units and requires a projected CRS.
* `directions`, `theta`, and `correction` are model parameters. The default `theta = 1` is a neutral package starting point, not a claim of scientific optimality.
* Pairwise surfaces are equally weighted by default through a mean. Weighted pairs are intentionally not inferred from side labels.
* `workers = 1` is the safe default. No core function detects available cores or adds scheduler-specific behavior.

## Open methodological questions

1. What defaults, if any, are appropriate for randomized-shortest-path `theta` across different ecological and epidemiological applications?
2. Which geographic correction is appropriate for different `theta` regimes and spatial resolutions?
3. When is reciprocal resistance-to-conductance conversion defensible, and when should a user-supplied transformation or transition function be required?
4. Should source-target pairs be weighted by sampling effort, population size, boundary length, or another scientific quantity instead of equally averaged?
5. How should arbitrary boundary segments, line terminals, and polygon or patch terminals be represented beyond cardinal edge generation?
6. Would an alternative graph backend improve scalability or numerical behavior if `gdistance` becomes limiting?

These questions are left configurable or deferred rather than resolved by hard-coded values from a particular downstream project.
