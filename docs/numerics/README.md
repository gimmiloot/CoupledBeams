# Numerical Policies

This directory contains project-wide numerical workflow policies. These
policies define how calculations are organized, checked, resumed, and
reported; they do not define the governing equations of individual physical
models.

- [Frequency-map computation policy](frequency_map_computation_policy.md) --
  the canonical `frequency-map-v1` contract for ordinary maps, strict audits,
  and rendering from saved data.
- [`scripts/sweep_grid_policy.py`](../../scripts/sweep_grid_policy.py) -- the
  existing helper for primary parameter grids and explicitly requested local
  refinements.

Model-specific equations, root-quality thresholds, and physical assumptions
remain in their own theory, implementation, and validation documents.
