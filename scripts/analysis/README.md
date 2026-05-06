# Analysis scripts

Most historical analysis and audit scripts are currently kept at `scripts/*.py` because several of them import each other through `scripts.*`.
New standalone diagnostics can live here when they do not require broad import churn.

Current entry points:

- `compare_analytic_fem_tracked_descendant_shape.py`: analytic determinant-nullspace vs FEM local component comparison for one tracked bending descendant.

See `scripts/README.md` for the full logical inventory and recommended commands.
