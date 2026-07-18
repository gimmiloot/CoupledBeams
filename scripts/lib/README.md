# Internal helpers

This directory contains reusable helper modules that are not meant to be run directly.

Project-wide branch identity and diagnostic-tracking rules are summarized in
`../../docs/project_rules.md`.

- `analytic_branch_tracking.py` is the source-of-truth helper for analytic branch identity. It tracks branches in memory from `beta = 0`, `mu = 0` for each `epsilon`, separates stable `branch_id` from `current_sorted_index`, and treats low-MAC assignments as non-canonical unless a diagnostic caller explicitly allows them.
- `analytic_coupled_rods_shapes.py` provides determinant-nullspace reconstruction, endpoint diagnostics, normalization, and analytic arm-energy utilities used by analytic shape and tracking diagnostics.
- `in_plane_shape_geometry.py` is the shared display-only geometry helper for
  in-plane analytic mode shapes. It keeps determinant components separate from
  Cartesian plotting coordinates, provides the reflected Timoshenko bases
  `t1=(1,0)`, `n1=(0,-1)`, `t2=(cos(beta),sin(beta))`,
  `n2=(sin(beta),-cos(beta))`, and exposes the equivalent EB mapping for EB's
  opposite transverse-field sign convention. It must not own coupling
  equations, determinant transforms, root selection, or mode reconstruction.
- `diagnostic_common.py` provides small non-scientific utilities for diagnostic
  scripts: filename-safe number tokens, compact number text, inclusive grids,
  output-directory creation, finite-value coercion, CSV row writing, and simple
  float-list parsing. It must not own formulas, determinant entries, or root
  selection policy.
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
  and isolated-rod reference utilities/conventions used in diagnostic
  `Lambda(mu)` plots. Reference curves are interpretation aids and must state
  their boundary-condition family, such as clamped-supported / clamped-pinned
  (CS/CP) or clamped-clamped / fixed-fixed (CC/FF).
- `tracked_bending_descendant_shapes.py` provides the shared tracked-state extraction, one-case normalization, one-case drawing, and output-path helpers used by both the single-shape and multi-panel tracked bending descendant commands.

FEM comparison logic is intentionally split: reusable FEM model code stays in
`../../src/my_project/fem/python_fem.py`, while diagnostic comparison and
normalization notes remain local to the corresponding scripts/reports unless a
future task asks for a shared helper.

Some historical helpers remain at root-level paths, especially `scripts/sweep_grid_policy.py`, because moving them would require broader import updates with no numerical benefit.

Lightweight tests can be run with `python -m unittest discover -s tests`.
`pytest` is optional when it is available in the active interpreter.
