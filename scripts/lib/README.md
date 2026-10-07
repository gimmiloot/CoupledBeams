# Internal helpers

This directory contains reusable helper modules that are not meant to be run directly.

- `weakly_nonlinear_planar_dynamics.py` restricts the audited quartic action
  to independent u,w,theta,c and discretizes it in an essential-BC Shen basis.
  Exact-degree quadrature, variable theta mass, both inertial terms and its
  analytic Jacobian support the [bounded time pilot](../../docs/theory/weakly_nonlinear_planar_time_pilot.md).
  Reversible M0 whitening is not mass replacement/modal truncation. Linear
  eigenpairs/time references, reconstruction, weak residual and energy checks
  are local numerical tools; no new physical closure or baseline API change.

- `weakly_nonlinear_spatial_rod.py` isolates the adopted seven-field reduced
  nonlinear model. Exact Rodrigues/Jr mass-form residuals and energies remain
  distinct from quartic-action/cubic polynomial expressions. Lazy exact
  Fraction derivation supplies independent A/B residuals and conjugate boundary
  covectors; no optional CAS install, mass inverse, spectrum solver or dynamics
  integration. [Model contract and audit](../../docs/theory/weakly_nonlinear_spatial_rod.md).

The [production Lambda geometry check](../../docs/theory/mindlin_herrmann_timoshenko_lambda_beta_large_checks.md)
composes `mindlin_herrmann_timoshenko_joint.arm_basis`, `arm_state`,
`arm_dynamic_stiffness` and the same `joint_matrix` for two separate section
objects. Physics helpers/APIs remain unchanged; section pairing and numeric
count-certified bracket prediction stay in the diagnostic CLI. No new
per-geometry solver/module is introduced.

- `coupled_longitudinal_comparators.py` is a bounded diagnostic comparator
  layer: existing H=0 elementary/planar Love axial basis plus identical
  verified Tim basis, six-DOF arms, six invariant joint rows and independent
  energy/Schur count with explicit pole-query conditioning. Love uses
  J=nu^2 rho Iy and N=(EA-J omega^2)u_x. No MH arm/joint API is changed.
  Geometric overlap and contraction norms support only the fixed-beta
  [screening contract](../../docs/theory/coupled_longitudinal_theory_hierarchy_screening.md).
  The subsequent [thickness screening](../../docs/theory/coupled_longitudinal_theory_thickness_screening.md)
  reuses this helper and the unchanged MH/Tim modules at each section;
  no separate thickness physics module or altered solver API was added.

Project-wide branch identity and diagnostic-tracking rules are summarized in
`../../docs/project_rules.md`.

- `mindlin_herrmann_timoshenko_joint.py` is the separate published reduced
  common-DOF rigid-joint helper: physical frame maps, endpoint signs, invariant
  eight-row residual and dual virtual work. Common `frame_boundary_matrix`
  supports project general beta and explicit swap/reflection frames; frozen
  beta0 API is a guard/wrapper of that same assembly. Local physics remains
  unchanged, global geometry mixes translations at nonzero beta. Exact
  arm Dirichlet-to-Neumann and nodal Schur matrices support the fixed pilot's
  energy count; bounded roots, full local modes/diagnostics and previous
  beta0/segmented-QR utilities stay here. No angle map or 3D joint theorem.
  Single/source modules do not
  import this helper. See [canonical joint note](../../docs/theory/mindlin_herrmann_timoshenko_rigid_joint.md).

- `mindlin_herrmann_longitudinal.py` is the isolated planar one-rectangle
  source-energy helper for M-H/Timoshenko diagnostics. Independent c,
  corrected source mass/stiffness, explicit factor inputs, natural boundary
  quantities, stable two-branch dispersion and analytic group velocities.
  It reuses the rectangular Timoshenko section coefficients and leaves that
  helper/API unchanged. Source factories have no implicit production defaults;
  no frame assembly or joint BC.
  Source fixtures: `data/input/mindlin_herrmann_timoshenko_sources.json`;
  [canonical audit](../../docs/theory/mindlin_herrmann_timoshenko_single_rod.md)
  retains variant-dependent source prescriptions and their historical
  unresolved coefficient audit. Source Jang requires explicit numeric kappa;
  Rucka factors are fitted.
  Subsequent selected project preset `project_jang_reduced_rectangular`
  uses the accepted rectangular K=5/6 in both shear-gradient terms, with
  unit contraction/rotary inertia. Source variants remain explicit and
  unchanged. Finite CC functions use bounded analytic columns, independent
  state-expm/QR, min-max count bounds and modal energy/mass checks; only the
  declared below-cutoff range is supported. No coupled-beam/joint solver.

- `bishop_longitudinal.py` is the isolated diagnostic longitudinal kernel:
  bounded analytic exponential/trigonometric bases, C/F/UP boundary and
  coaxial-interface assembly, separate H=0 equations, finite root search,
  mass/energy/ODE checks, Popov (13), and addressed high-precision Marais
  state-exponential verification. Source fixtures live in
  `data/input/bishop_literature_sources.json`; the
  [canonical report](../../docs/theory/bishop_literature_reproduction.md)
  defines signs, scope and qualifications. It does not import or alter the
  baseline bending/joint solver and has no runnable per-case scripts.

  The separate [single-rod kinematic audit](../../docs/theory/timoshenko_bishop_single_rod.md)
  reads this module and `isotropic_rectangular_timoshenko_coupled_beams.py`
  only as reference limits. It adds no combined helper/API: candidate fields
  have zero centered cross terms but incompatible raw bending self terms.

- `analytic_branch_tracking.py` is the source-of-truth helper for analytic branch identity. It tracks branches in memory from `beta = 0`, `mu = 0` for each `epsilon`, separates stable `branch_id` from `current_sorted_index`, and treats low-MAC assignments as non-canonical unless a diagnostic caller explicitly allows them.
- `analytic_coupled_rods_shapes.py` provides determinant-nullspace reconstruction, endpoint diagnostics, normalization, and analytic arm-energy utilities used by analytic shape and tracking diagnostics.
- `in_plane_shape_geometry.py` is the shared display-only geometry helper for
  in-plane analytic mode shapes. It keeps determinant components separate from
  Cartesian plotting coordinates, provides the reflected Timoshenko bases
  `t1=(1,0)`, `n1=(0,-1)`, `t2=(cos(beta),sin(beta))`,
  `n2=(sin(beta),-cos(beta))`, and exposes the equivalent EB mapping for EB's
  opposite transverse-field sign convention. It must not own coupling
  equations, determinant transforms, root selection, or mode reconstruction.
- `reddy_inplane_geometry.py` is the separate three-dimensional physical-frame
  helper for the RLB coordinate gate. Both local arm coordinates point from
  outer clamp to the future connection point; the helper defines the fixed
  global view basis, right-handed Reddy triads, and physical vector mappings
  used in virtual-work checks. It must remain separate from display-only
  geometry, production FEM transforms, connection equations, determinants,
  and root calculations.
- `reddy_symmetric_coupled_beams.py` is the narrow RLB-1 rigid-joint
  helper. It derives the canonical joint matrix from the physical maps,
  provides an independent closed-form comparator, reuses the verified
  single-beam state and transfer matrices, and exposes independent direct
  fixed--fixed and stepped references for the diagnostic `beta=0` pilot.
  It contains no material reduction, root finder, angle sweep, Ritz model,
  FEM, torsion, damping, or legacy coupled-rod imports.
- `reddy_symmetric_coupled_beams_ritz.py` is the independent RLB-1C
  two-arm constrained Rayleigh--Ritz helper. It assembles only the frozen
  physical energies and three endpoint kinematic constraints from reduced
  beam properties and `reddy_inplane_geometry.py`. It contains no transfer
  matrix, matrix exponential, determinant, root finder, force-equilibrium
  joint row, FEM, torsion, damping, or legacy coupled-rod import. The current
  `N=16` beta=0 guard does not pass the full first-13 bridge, so this module
  has not been used for a nonzero-angle spectral claim.
- `isotropic_rectangular_timoshenko_coupled_beams.py` is the independent
  closed-form rectangular Timoshenko comparator used by the finite
  four-equal-ply isotropic-limit audit. It supports mixed, exact-cutoff, and
  two-trigonometric spatial regimes and retains a circular-section backcompat
  path. It imports no RLB transfer, joint, laminate-reduction, Ritz, FEM, or
  Euler--Bernoulli module and does not contain a global root finder.
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
