# Internal helpers

This directory contains reusable helper modules that are not meant to be run directly.

Project-wide branch identity and diagnostic-tracking rules are summarized in
`../../docs/project_rules.md`.

- `analytic_branch_tracking.py` is the source-of-truth helper for analytic branch identity. It tracks branches in memory from `beta = 0`, `mu = 0` for each `epsilon`, separates stable `branch_id` from `current_sorted_index`, and treats low-MAC assignments as non-canonical unless a diagnostic caller explicitly allows them.
- `analytic_coupled_rods_shapes.py` provides determinant-nullspace reconstruction, endpoint diagnostics, normalization, and analytic arm-energy utilities used by analytic shape and tracking diagnostics.
- `thickness_mismatch_mac_tracking.py` provides diagnostic-only analytic shape
  reconstruction and adjacent-step MAC tracking for the mass-preserving
  thickness-mismatch eta model. It keeps nearest-frequency assignment only as a
  warning comparator, separates raw candidate assignments from accepted
  canonical sorted positions, and records diagnostic flags such as low MAC, low
  margin, unresolved assignments, sorted-position jumps, suspicious
  assignments, and refined-check requests.
- `thickness_mismatch_diagnostic_helpers.py` collects plotting/report helpers
  for thickness-mismatch diagnostics: fixed-eta descendant tracking wrappers,
  diameter-to-length validity summaries, solid/dashed applicability plotting,
  and the clamped-supported isolated-rod reference convention used in
  diagnostic `Lambda(mu)` plots.
- `tracked_bending_descendant_shapes.py` provides the shared tracked-state extraction, one-case normalization, one-case drawing, and output-path helpers used by both the single-shape and multi-panel tracked bending descendant commands.

Some historical helpers remain at root-level paths, especially `scripts/sweep_grid_policy.py`, because moving them would require broader import updates with no numerical benefit.
