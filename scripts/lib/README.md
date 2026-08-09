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
- `family_inventory_local_repair.py` provides the diagnostic sorted-family
  missing-root detector, source-derived local-window inference, staged local
  matrix repair, multiplicity-aware merge, and isolated atomic cache used by
  `audit_family_inventory_local_repair.py`. It does not define descendant
  identity and does not call tracking, MAC, shapes, or strict verification. Its
  cache identity explicitly accepts only the isotropic circular coupled-rod
  EB/Timoshenko scope.
- `article_epsilon_family_inventory_integration.py` is the parent-process
  orchestration adapter for the article epsilon grid. It groups immutable
  pointwise sorted spectra by `(epsilon_0, mu, eta, theory)`, reuses the family
  detector and local matrix repair before an expensive-strict defer decision,
  and writes a separate shadow/provenance overlay. It contains no scientific
  matrices and does not import the rectangular-anisotropic research workflow.
- `article_epsilon_family_reconciliation.py` is the explicit zero-solve
  promotion layer for that verified shadow. It accepts only the isotropic
  circular EB/Timoshenko scope, validates source fingerprints and shadow gates,
  promotes only provenance-complete matrix-confirmed rows, preserves deferred
  rows with `N_true=NaN`, and writes the deterministic article-facing table and
  future resume plan without importing or calling a solver, matrix evaluator,
  detector, or local repair. The source point cache remains immutable.
- `article_epsilon_compact_certificates.py` is the zero-solve streaming
  migration layer for the same isotropic circular grid. It reads one full gzip
  trace at a time and writes a versioned per-case certificate containing only
  sorted roots through the scientific guard, `delta_f`, `N_true`, compact
  quality flags, and provenance. `article_epsilon_compact_poststage.py` groups
  only these compact records one beta-family at a time, permits narrow local
  matrix repair only for unresolved cases, and emits the scalar article-facing
  table plus a non-destructive raw-cache retention proposal. Neither module
  contains scientific matrices or imports the rectangular-anisotropic scope.
- `article_epsilon_targeted_resolution.py` is the target-only orchestration
  layer for deferred `epsilon_0=0.050` compact cases. It selects IDs from the
  unresolved table, reads at most one raw payload at a time, verifies only the
  required sorted prefix with four shifted local determinant/SVD phases and
  two refinement levels, and writes an immutable overlay plus a versioned
  finalization. It does not change the production matrices or tolerances; all
  target caches are isolated and force/full strict remains unused when T1 and
  the stored independent configuration agree.
- `yartsev_ch2_fast_beta_sweep.py` is the diagnostic-only generic coordinator
  for the Chapter-2 supervisor `Lambda(beta)` calculations. It provides
  sorted-frequency prediction windows, connected close-root clusters,
  mandatory global anchors and fallback, exact bounded transfer-matrix LRU
  caching, separate performance counters, and atomic family checkpoints. It
  contains no physical equations or boundary matrices, does not define modal
  descendants, and retains the existing global solver as oracle/fallback.
- `yartsev_ch2_monoclinic_rod.py` also exposes a narrow diagnostic propagation
  bridge from an arbitrary physical initial state and the same
  bending/shear/generalized-torsion energy components already used by the
  cantilever diagnostic. The supervisor workflow uses them only at `beta=0`
  to characterize sorted positions of the two-arm Chapter-2 spectrum; no
  determinant, material rotation, clamp, or joint equation is redefined.

FEM comparison logic is intentionally split: reusable FEM model code stays in
`../../src/my_project/fem/python_fem.py`, while diagnostic comparison and
normalization notes remain local to the corresponding scripts/reports unless a
future task asks for a shared helper.

Some historical helpers remain at root-level paths, especially `scripts/sweep_grid_policy.py`, because moving them would require broader import updates with no numerical benefit.

Lightweight tests can be run with `python -m unittest discover -s tests`.
`pytest` is optional when it is available in the active interpreter.
