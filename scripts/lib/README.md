# Internal helpers

- `two_springs_horizontal_combined_thresholds.py` closes HS-L1–HS-LA6:
  archived right-angle level and saturation coefficients produce two direct
  leading roots and independent two-term forecasts before full references
  are parsed. Its staged archival contract differs from the old full-threshold
  runner; frozen physics and source hashes remain unchanged. It reuses HS-L1
  Ritz and HS-LA4 prefix checks, with no full-model runner or new solver.


- `two_springs_horizontal_right_angle.py` is the bounded HS-LA5 right-angle
  limit workflow. It reuses HS-L1 integrated fields with the independently
  derived trace space and endpoint mass, verifies the congruence to the
  saturated inner stiffness, and checks exact leading-ODE roots and counts.
  Its different mass/boundary contract is kept separate to preserve source
  hashes of accepted HS-L1–HS-LA4 workflows. No full-model solver is called.


- `two_springs_horizontal_joint_stiffness.py` extends the unchanged HS-L1
  fields for the bounded HS-LA4 contract: exact zero-kernel deflation, exact
  J/D constraints, reaction coefficients and threshold sensitivities. Its
  independent EB endpoint-inertia count checks the lowest three positions.
  This is a reusable diagnostic entry point, not a full-model or universal
  solver. A separate module preserves old source hashes and avoids applying
  the positive-definite HS-LA2 contract to the zero mechanism. See section 12
  of the canonical spectral-structure report and `../README.md` for commands.


- `two_springs_horizontal_saturation.py` postprocesses HS-LA2 matrices and
  authoritative modes for the HS-LA3 second saturation coefficient. It has
  no eigensolver or fit; its small linear-response contract preserves the
  existing source-hashed HS-LA2 runner. See section 11 of the same report.
- `two_springs_horizontal_spectral_structure.py` is the bounded HS-LA2
  all-N secular/energy/reaction audit of the frozen HS-L1 matrices. It preserves
  high-mode arithmetic diagnostics, uses exact saturation and has no full-model
  spectral calls. See `docs/theory/isotropic_two_springs_horizontal_leading_spectral_structure.md`.
- `isotropic_two_springs_horizontal_leading.py` is the specialized HS-L1
  variational leading model with free external translation/rotation, analytic
  curvature integration, rank-one difference of endpoint displacements and
  exact saturation. An independent eight-condition EB ODE checks selected roots.
  `two_springs_horizontal_leading_validation.py` manages its bounded p20/24
  checks, four full-model threshold searches, archive reuse and CSV-only plots.
  See `docs/theory/isotropic_two_springs_horizontal_leading.md`.

- `two_springs_horizontal_pilot.py` implements HS-P1: the fixed seven-case
  input contract, verified archive reuse, three-position prefixes, triggered
  neighbour checks and analysis of saved first forms. It reuses the physical
  module and the frozen diagnostic `Run`, with no new spectral infrastructure.
  Its `plot_only` helper `two_springs_horizontal_pilot_figures.py` reads checked
  CSV and creates four PDF/PNG figures without solving eigenvalue problems.
  See `docs/theory/isotropic_two_springs_horizontal_pilot.md`.

- `isotropic_two_springs_horizontal_spring.py` is the separate physical
  one-horizontal-spring model with no outer guide or rotation spring.
  Reuses canonical RL_FULL/Timoshenko equations and node; adds the seven/six/five
  dimensional admissible trace spaces, exact zero kernel, low-frequency
  fundamental basis and independent semidefinite Ritz energy verification.
  Physical positions include zeros. Canonical verification is documented in
  `docs/theory/isotropic_two_springs_horizontal_spring.md`; the separate HS-P1
  workflow above uses this unchanged physical implementation.

- `two_springs_external_axial_mechanism.py` implements the bounded
  `--transition-mechanism` preset. It inherits the physical/count solver,
  source validation, sweep and range refinement; adds fixed-angle brackets,
  reflection checks, full-model curvature uncertainty and extrema evidence.
  Additional stiffnesses require the saved scientific trigger. No new
  physical assembly or universal spectral infrastructure is introduced.

- `two_springs_external_axial_transition.py` is the bounded small-stiffness
  preset of the axial CLI. It inherits source/physical/count validation,
  map traversal and angle-range refinement, adds the nine-angle preflight,
  signed growth, fixed-angle sign brackets and free-limit checks, and renders
  accepted data. No physical, leading or universal spectral solver is added.

- `two_springs_external_axial_high_ranges.py` implements the bounded
  `--high-frequency-ranges` preset of the existing axial CLI. It validates
  archive identity, reuses the axial solver/field gates and `high_ranges`,
  extends only missing tails and triggered brackets, retains nonterminal
  qualifying intervals, and renders four figures from accepted data.
  A temporary legacy column view reuses refinement logic without treating
  axial chi_u as transverse k_w. No physical or universal spectral solver.

- `two_springs_high_frequency_angle_ranges.py` implements the bounded
  `--high-frequency-ranges` workflow: verifies source archives, reuses existing
  full-model/count/physical gates, adds missing or reasoned refinement points,
  and determines all discrete 95%/98% ranges including terminal status.
  It also renders four figures from saved CSV. It is not a spectral solver;
  the existing transverse-support physics is unchanged.
  The `--high-frequency-ranges-90-95` preset reuses that completed archive,
  reproduces old maxima/thresholds, permits only 10%-boundary refinement,
  and handles zero-angle starts without treating them as missing values.
- `isotropic_two_springs_external_rotational_replacement.py` creates the
  existing analytic-curvature leading space and assembles K0 without Ke,
  retaining a free outer rotation. Full physics remains in the transverse library.
- `two_springs_external_rotational_replacement.py` reuses full/leading point
  validation, spectrum counting, I/O and saved-data guards for the bounded
  `--rotational-replacement` CLI mode; A/B/C/D comparisons and four figures.

This directory contains reusable helper modules that are not meant to be run directly.

- `two_springs_external_transverse_leading_resume.py` orchestrates the completed
  leading-only transverse continuation through the existing study CLI. It
  reuses `Leading` with analytic curvature, old archive validators and the
  frequency normalization helper; adds per-position trace/strong evidence and
  finite-Ritz-prefix inertia checks. Full-model spectral calls are explicitly
  forbidden. Rendering and finalization use saved results only.

- `isotropic_two_springs_external_transverse_spring.py` adds only the local
  transverse outer row, trace and spring energy to the full model, plus one
  constant w2 basis function in direct/leading Ritz spaces with consistent mass.
  `two_springs_external_transverse_spring.py` preserves the accepted full map
  and explicit leading strong-residual STOP; it forbids automatic resumption
  and renders only qualified full-model figures from saved data.

- `isotropic_two_springs_external_axial_spring.py` specializes only the outer
  local axial row, trace degree of freedom and spring energy. Private hooks reuse
  unchanged arm/count/root code; direct-energy Ritz adds one constant u2 function.
  `two_springs_external_axial_spring.py` provides the bounded experiment, archive
  comparisons, verified existing leading solutions and guarded saved-data plots.

- `two_springs_thickness_leading_comparison.py` compares eight archived full
  spectra with the unchanged individual-section leading Ritz model on matching
  xi nodes. It guards against full solves, reuses archive validation/writers,
  and records signed errors, thickness effects and addressed convergence checks.
  It contains no physical assembly or new spectral solver.

- `isotropic_two_springs_thickness_asymmetry.py` specializes the canonical
  RL_FULL/Timoshenko equations for individual rectangular sections at fixed
  volume. It reuses arm fundamentals and prefix isolation, while assembling
  individual energies, pole counts and weighted leading Ritz forms. Reference
  spring/frequency scales and original module globals remain unchanged.
- `two_springs_thickness_asymmetry.py` orchestrates the frozen 18-case study:
  algebraic and physical gates, archive reuse, conditional tails, local leading
  derivatives and saved-data figures. It reuses established stencils and checks;
  it is not a general spectral runner. See the dedicated CLI in `scripts/README.md`.

- `two_springs_length_asymmetry.py` is the bounded orchestration and plotting
  helper for the existing two-spring CLI's `length-asymmetry` preset. It reuses
  the canonical spectral/field checks and leading Ritz implementation, freezes
  the signed-length contract, and separates conditional tail/derivative checks
  from zero-solve rendering. It contains no new physical equations or solver.

- `nonsymmetric_laminated_kelvin_voigt.py` is the bounded NSL-3A complex
  adapter: unchanged coupled NSL H, centered propagation, analytic Frechet
  derivatives, full node and complete Hermitian energies including B/j1.
  Row/column scales are frozen at an elastic seed; SVD coefficients are
  transformed back as in NSL-2R. It reuses the old KV corrector and physical
  node rows, without importing the EB/RLB internal operator for NSL arms.
- `nonsymmetric_laminated_joint.py` is the isolated NSL-2 elastic angular
  adapter: NSL-1 physics, existing RLB physical node rows, covariant reference
  shifts, coupled exact-segment inertia and independent direct-energy Ritz.
  No material integration copy, damping or universal solver.
- `nonsymmetric_laminated_beam.py` is the isolated NSL-1 three-field FSDT
  model: labeled direct-source and coupled-ABD Schur reductions, m/j1/j2,
  energy-derived state operator and one CC propagation/count adapter.
  It imports the old layer integrator and essential Ritz basis; the new
  Ritz energy assembly includes B and j1 without calling a transfer API.
  Principal-axis 0/90 plies only; no angular joint or damped-spectrum code.
- `uniaxial_kelvin_voigt_split.py` extracts the saved seven-field quadratic
  action and implements the separate linear 14-condition one-axis KV joint.
  It checks physical vector/work duality, the exact block split and same-BC
  rigid OUT identity, then reuses finite Tim and scalar torsion bases for
  one G20 OUT prefix. Its three-trace Schur count has its own energy proof;
  old MH trace spaces and the Yartsev slope clamp are not substituted.
- `longitudinal_kelvin_voigt_comparison.py` provides the narrow elementary/RL
  complex continuation used by the fixed G20 comparison. It preserves the
  Tim/Jang signs, uses N=(EA+J*p^2)*U', differentiates a fixed bounded basis
  analytically and reuses existing assembly, physical masses and KV Newton.
  MH complex results and all elastic coefficients are read from verified
  completed data; this helper is not a universal solver or production change.

- `contraction_joint_sensitivity.py` is the isolated G20 FREE/CONTINUOUS/
  BLOCKED c diagnostic: variational trace maps, physical spring conditions,
  energy-index prefix isolation and local complex MH basis. It reuses the
  unchanged MH/Tim and reduced axial arms, DtN helpers and KV Newton corrector.
  It does not modify production dispatch or define a general sweep framework.

- Seven-field spatial verification reuses --spatial-verification on the
  existing continuation CLI. `weakly_nonlinear_spatial_dynamics.py` compiles
  the saved quartic action; `nlsp_spatial_verification_checks.py` retains exact
  planar/energy/parity evidence and separate strict strong/weak-flux checks.
  `nlsp_spatial_1d_program.py` owns bounded static/release/Radau cases;
  `nlsp_spatial_comparison.py` reads saved histories and correction controls.
  `nlsp_spatial_fem_protocol.py` and `nlsp_spatial_native_program.py` reuse the
  verified CCX generator, C3D10 quadrature, streaming native readers and
  section orientations without changing historical physics/helpers.
  [Scoped contract](../../docs/numerics/nlsp_spatial_nonlinear_3d_fem_verification.md).

- The saved-profile audit uses `--profile-audit` on the existing continuation
  CLI. `nlsp_spatial_profile_audit.py` owns source/hash selection, separate
  generated diagnostics and saved-data figures. `nlsp_profile_fem_diagnostics.py`
  reuses the historical section recovery and audited C3D10 quadrature for raw
  fits, unsmoothed interpolation, independent gradient/strain averages and
  one explicitly separate quadratic-fit probe. `nlsp_profile_1d_diagnostics.py`
  restores saved p48/p64 coordinates and derives contraction, symmetry,
  resultant and spectral-tail diagnostics; it never integrates motion.
  [Scoped postprocessing contract](../../docs/numerics/nlsp_spatial_profile_audit.md).

- FEM-3C reuses the same public continuation CLI through `--validation`.
  `nlsp_fem3c_validation.py` owns the six-case authorization/ledger, saved
  medium/fine source selection and conditional full-period gates. Its I/O
  adapter preserves each mesh-level preload and explicit dynamic output cadence.
  `nlsp_fem3c_diagnostics.py` transfers saved profiles onto the fixed 201-time
  grid, checks linear/PCHIP sensitivity, records fixed/updated metric scales,
  maps all seven fields and plots saved data. `nlsp_fem3c_1d.py` restores own
  static coordinates and complete linear factors, using the unchanged Radau
  runner for at most one p64 and one p48 nonlinear full-period solve. No new
  constitutive, element, mass, RHS/Jacobian, mesh or eigensolver implementation.
  [Scoped contract](../../docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification).

- FEM-3B `nlsp_fem3b_continuation.py` is a non-runnable scoped helper for the
  existing resume command. It records the separate parent/authorization,
  predeclared horizon decision, accepted Radau cubic output and actual native
  attempt ledger, reusing the unchanged four-field RHS/Jacobian and C3D10
  generation/recovery. Dense evaluation uses the accepted Radau polynomials,
  without a second ODE solve or an extra step multiplier. Explicit interpolation
  of differing 3D schedules stays separately marked with linear/PCHIP diagnostics.
  `nlsp_fem3b_diagnostics.py` reads saved profiles/native velocities for comparison,
  existing mass-weighted kinetic diagnostics and three figures; it never calls
  native solvers, ODE integration, equilibrium or eigenanalysis.
  [Contract](../../docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3b-nonlinear-correction-evolution-and-longer-horizon).

- FEM-3AR reuses the FEM-3A helpers/physics and adds no new library solver.
  A targeted repair in `nlsp_fem3a_transient_output.py` derives timestamp
  rounding intervals from native lexical precision, preserving raw tokens,
  unique actual increment association and rejection of inconsistent/ambiguous
  data. Correct successful files are reparsed with both snapshots retained;
  no interpolation, new dynamic constraint or changed scientific tolerance.
  Native energy/pseudo-time output semantics remain explicitly qualified.
  [Continuation evidence](../../docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3ar-controlled-continuation-after-native-elke-output-path-failure).

- FEM-3A narrow helpers `nlsp_fem3a_1d_reference.py` and
  `nlsp_fem3a_transient_output.py` restore frozen p64 static coordinates into the
  unchanged four-field mass/RHS/Jacobian/Radau implementation and stream actual
  CalculiX step/increment/time/field/energy output. They contain no new physical
  law or standalone FEM solver; static-origin metadata is not relabelled a
  dynamic frame. [Pilot contract](../../docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md).

- FEM-2R reuses the unchanged FEM-2 deck writer, static readers, reference-section
  polar recovery, consistent-bodyload reaction correction and report functions.
  Its scoped analysis CLI only orchestrates a separately authorized continuation
  over saved 1D profiles/meshes; it adds no library physics solver or dynamic
  helper change. The historical failed-attempt guard remains intact.
  [Continuation contract/result](../../docs/numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure).

- FEM-2 reuses unchanged `PlanarGalerkin.linear_stiffness`, `potential` gradient/
  Hessian and reconstruction for static stationarity of V_le4. A scoped analysis
  layer loads saved C3D10 meshes, builds static GRAV decks and contains final-
  increment-aware static parsing/recovery; no new library physics solver or
  dynamic helper change. Synthetic recovery does not substitute for the failed
  real medium equilibrium. [Static contract/result](../../docs/numerics/nlsp_nonlinear_static_3d_fem_validation.md).

- FEM-1R loads the immutable saved1D and three-grid reference, then reuses
  FEM-1 geometry/audit/modal/section/MAC functions for one .020 continuation.
  No new library/physics solver or modified historical helper is introduced.
  [Contract/result](../../docs/numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement).

- Full-family FEM-1 reuses the unchanged MH/Tim analytic finite-rod operators
  and Yartsev generalized isotropic torsional reduction. The new scoped solid
  orchestration/mesh/shape diagnostics live in one analysis CLI, not a new
  physics module. No GI_p substitution or book slope-clamp transfer.
  [Rectangular linear3D report](../../docs/numerics/nlsp_linear_rectangular_3d_fem_validation.md).

- The physical-sanity CLI reuses these unchanged prepared/planar helpers and
  runner. It reads saved leading stat/harm profiles directly and performs only
  one separately authorized half-amplitude ODE; reactions use canonical V4
  derivatives and explicit finite-p/rotation-truncation accounting. No second
  solver or new physical closure is added.
  [Bounded diagnostic contract](../../docs/theory/weakly_nonlinear_planar_physical_sanity_checks.md).

- One-T1 orchestration restores the already prepared q0 directly; it does
  not call any new numerical/MP projection in this helper. The existing Radau
  runner additionally accepts an optional writable float64 history buffer for
  memory-mapped output; its default numerical initial/stepping path is retained.
  [One-period contract](../../docs/theory/planar_prepared_initial_state.md#prepared-one-period-feasibility).

- `planar_prepared_initial_state.py` additionally loads the immutable saved
  common state without new BVP/eigenpair/Theta3, computes exact Gram/analytic
  moments at local high precision, includes quintic Theta3 directly and supplies
  a shared initial-only endpoint-constrained L2 representation. It retains the
  full later trial/test space and never changes the admission flag to bypass
  strict checks. The existing autonomous runner accepts explicit full q0
  coordinates as an optional input; its default initial projection is preserved.
  [Precision/feasibility scope](../../docs/theory/planar_prepared_initial_state.md#prepared-precision-feasibility).

- `planar_prepared_initial_state.py` restores constant/second-harmonic parts
  from saved full spectra and provides physical Legendre profiles, exact
  endpoint-trace checks and one quintic Hermite correction/common initial
  evaluator. An admission flag separates an algebraic candidate from an
  accepted initial state. It performs no eigensolve or time integration and
  does not modify the old dynamics helper or initial-condition workflow.
  [Preparation contract](../../docs/theory/planar_prepared_initial_state.md).

- `planar_second_order_axial_response.py` is a focused diagnostic companion
  for the leading forced u2,c2 response of the audited quartic action. It
  extracts the source terms, reuses the full linear M-H Shen space and evaluates
  stable exact-time kernels with analytic velocities/accelerations, including
  small-time and close-frequency limits. All 2(p−1) coordinates are retained;
  mass whitening and a full generalized eigendecomposition introduce no
  filtering or modal reduction. The prescribed continuous bending background
  is common to every p. This helper does not advance or change the existing
  nonlinear RHS, initial data or physical model.
  [Diagnostic contract](../../docs/theory/planar_second_order_axial_response.md).

- `weakly_nonlinear_planar_dynamics.py` restricts the audited quartic action
  to independent u,w,theta,c and discretizes it in an essential-BC Shen basis.
  Exact-degree quadrature, variable theta mass, both inertial terms and its
  analytic Jacobian support the [bounded time pilot](../../docs/theory/weakly_nonlinear_planar_time_pilot.md).
  Reversible M0 whitening is not mass replacement/modal truncation. Linear
  eigenpairs/time references, reconstruction, weak residual and energy checks
  are local numerical tools; no new physical closure or baseline API change.
  The targeted recovery splits energy-only, gradient and Hessian requests so
  ordinary RHS calls do not evaluate unused second derivatives. State-aware
  caches retain the same variable-mass solve and analytic Jacobian; numerical
  equivalence is checked against the preserved historical implementation.
  This optimization does not resolve the reported spatial-convergence gates.

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
- `two_springs_transverse_leading_diagnostic.py` evaluates saved transverse
  leading Ritz pairs with independent analytic derivatives and high precision,
  isolates assembly roundoff and checks a six-order single-case correction.
  Rational polynomial integrals evaluate the existing space; no new physical
  solver or frequency-map runner. The opt-in `analytic_curvature` flag in
  `isotropic_two_springs_external_transverse_spring.Leading` changes only the
  numerical evaluation of the curvature Gram matrix, preserving the old default.
