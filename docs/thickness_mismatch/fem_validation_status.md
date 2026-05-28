# FEM Validation Status

## Purpose

This document records the current diagnostic-only FEM validation status for
the thickness-mismatch / Timoshenko correction work. The FEM workflows are used
as independent checks of the one-dimensional models:

- Euler--Bernoulli beam references;
- diagnostic Timoshenko beam references;
- the planned thickness-mismatch / Timoshenko formulation.

The notes below are not article-level claims. They define what is currently
supported, what remains diagnostic, and what must be improved before using FEM
results in a paper.

## Current FEM Workflows

### `solid_fem_single_rod_fixed_fixed.py`

This workflow checks one circular rod with fixed-fixed ends using 3D solid FEM
through Gmsh and CalculiX. It generates the solid mesh, solves modal
frequencies when the external tools are available, parses displacement modes,
classifies bending/axial/torsion-like shapes, and pairs the circular-rod
bending doublets.

Current diagnostic conclusion: for the classified bending doublets, the
Timoshenko reference is closer to the 3D solid FEM result than the
Euler--Bernoulli reference.

### `solid_fem_coupled_equal_rods.py`

This workflow checks equal coupled rods using fused-cylinder geometry for
`beta = 15, 45, 90 deg`. It is useful as a geometric diagnostic because it
shows how the finite fused solid joint differs as the angle changes.

Current diagnostic conclusion: this workflow is not a direct validation of the
ideal one-dimensional point-joint model. The fused 3D joint has finite volume,
and at small angles the overlap/merger zone can be large relative to a rod
segment.

### `solid_fem_coupled_equal_rods_point_joint.py`

This workflow checks equal coupled rods using two separate solid cylinders and
a shared CalculiX `*RIGID BODY` reference joint at the inner end faces. This is
closer to the ideal one-dimensional point joint than the fused-cylinder model.
It also includes MAC-like centerline shape matching between parsed 3D solid
modes and reconstructed one-dimensional EB/Timoshenko shapes.

Current diagnostic conclusion: after MAC-like matching, more modes are closer
to Timoshenko than before. The current full beta/epsilon frequency-closeness
count after MAC-based mode selection is 23/36 closer to Timoshenko and 13/36
closer to Euler--Bernoulli. A first beta=15 mesh-convergence run has now been
completed for `epsilon = 0.025, 0.05` and mesh factors `1.5, 1.0, 0.75`;
the strong/moderate MAC rows were stable under the final refinement screen,
but weak rows and duplicate raw best-match warnings remain. An auxiliary
planar-constrained diagnostic has also been added for `beta = 15 deg` and
`epsilon = 0.025, 0.05`; it suppresses out-of-plane/torsion-like solid modes
to check the planar subspace corresponding to the 1D model. This is promising
diagnostic evidence, but it is not yet article-level validation.

## What Is Already Supported

- The single-rod 3D solid FEM workflow gives strong support for the
  Timoshenko correction in the fixed-fixed circular-rod benchmark.
- The fused-cylinder coupled FEM workflow should be treated as a geometry audit,
  not as direct validation or refutation of the point-joint one-dimensional
  model.
- The point-joint 3D solid FEM workflow is the current best validation path for
  coupled rods because it avoids the finite fused overlap volume.
- MAC-like centerline shape matching is preferred over sorted-frequency
  comparison and over simply taking the first classified in-plane modes.
- The current point-joint results are promising, and the first beta=15
  mesh-convergence diagnostic is stable for strong/moderate MAC rows, but the
  workflow remains diagnostic until weak/duplicate matches and visual
  mode-shape checks are resolved.
- The planar-constrained point-joint run is an auxiliary diagnostic of the
  planar subspace only; it must be labeled as such and must not replace the
  full 3D point-joint benchmark.
- Equal-thickness `Lambda(mu)` diagnostics must report Euler--Bernoulli
  diameter applicability and Timoshenko cut-off applicability separately:
  the EB `2*r_i/l_i <= 0.1` criterion does not dash or invalidate
  Timoshenko curves, whose diagnostic margin is `Omega/Omega_c`.
- The variable-length Timoshenko determinant used in the equal-thickness
  `Lambda(mu)` diagnostic is now diagnostic-verified for `eta=0` sorted-root
  limit checks: `mu=0` consistency with the older equal-rods Timoshenko
  diagnostic, the `epsilon -> 0` Euler--Bernoulli limit, and `beta=0`
  straight-rod `mu`-invariance all pass in
  `results/variable_length_timoshenko_limits_audit.md`.
- The tau-aware `eta != 0` variable-length Timoshenko extension is now
  diagnostic-verified for sorted-root limit checks in
  `results/variable_length_thickness_timoshenko_limits_audit.md`: eta=0
  regression, eta-to-zero continuity, rod-swap symmetry, the `epsilon -> 0`
  Euler--Bernoulli thickness-mismatch limit, and beta=0 straight composite-rod
  behavior all pass for the tested first six sorted roots. This does not verify
  descendant branch identity, energy plots, or FEM behavior.
- The joint-coupling-model audit at `epsilon=0.01`, `mu=0`, and `eta=0`
  now compares full rigid end-face coupling with central rigid patches for
  `beta = 0, 15, 45, 90 deg`. The `center_patch_0p5` variant is the current
  best diagnostic candidate: it passes the beta=0 straight sanity check and
  reduces the beta=15 mean best relative error from about 7.5% for full rigid
  end faces to about 4.7%. This is not article-ready validation.
- The joint-constraint bracketing audit at `epsilon=0.01`, `mu=0`, and
  `eta=0` now tests physically different constraint classes rather than using
  patch radius as a fit parameter. Translation-only central/average equation
  probes are classified as too soft and fail the beta=0 sanity check; the
  rotation-clamped reference reproduces the rigid-end-face beta=15 errors and
  remains incompatible with an analytic point joint whose rotation is not
  absolutely clamped; the attempted two-reference-node CalculiX deck remains
  unsupported; and `center_patch_0p5` remains only a diagnostic candidate
  needing mesh and shape review.
- The rigid end-face thickness trend audit now treats the full rigid end-face
  point-joint as the main physically clear 3D engineering model rather than
  tuning a patch radius. For `beta=15 deg`, `mu=0`, `eta=0`, and
  `epsilon=0.005,0.01,0.025,0.05`, all Gmsh/CalculiX cases complete without
  recorded solver warnings. Strong/moderate MAC rows show the mean
  Timoshenko/EB relative-error ratio decreasing from about `1.005` at
  `epsilon=0.005` to about `0.621` at `epsilon=0.05`, supporting the
  comparative claim that the Timoshenko reference becomes more favorable as
  thickness increases. This remains diagnostic evidence, not final article
  proof.

## Rigid-Joint Lambda(mu) Thickness-Trend Plots

The current main diagnostic visualization of the EB/Timoshenko/3D FEM
comparison along `mu` is the full rigid end-face `Lambda(mu)` epsilon sweep:

- `results/rigid_joint_equal_thickness_beta15_eps0p01_lambda_mu_eb_timoshenko_fem.png`
- `results/rigid_joint_equal_thickness_beta15_eps0p025_lambda_mu_eb_timoshenko_fem.png`
- `results/rigid_joint_equal_thickness_beta15_eps0p05_lambda_mu_eb_timoshenko_fem.png`
- `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem.csv`
- `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem_report.md`

These plots cover `beta = 15 deg`, `eta = 0`, equal thicknesses,
`epsilon = 0.01, 0.025, 0.05`, `mu = 0..0.9`, and the first six descendant
branches. Each plot shows Euler--Bernoulli descendant curves, Timoshenko
descendant curves, and MAC-filtered 3D FEM points.

The 3D FEM benchmark is a non-tuned full rigid end-face engineering joint: two
separate solid cylinders, clamped outer ends, and all inner end-face nodes
coupled by CalculiX `*RIGID BODY` to a common `JOINT_REF`. It uses no planar
constraint, no patch tuning, and no fused volume. Weak, ambiguous, and
duplicate FEM matches are retained in the CSV/report but are excluded from the
main plotted FEM points; the visible FEM points are strong/moderate
non-ambiguous MAC matches.

These plots should be interpreted comparatively. Exact coincidence with the
one-dimensional point-joint model is not expected, because the 3D benchmark is
a physically defined full rigid end-face joint rather than a tuned replica of
the analytic point joint. The current diagnostic result supports the
qualitative trend that Timoshenko becomes more appropriate than
Euler--Bernoulli as thickness increases. It is still not final article
validation until mesh convergence, representative mode-shape inspection, and a
reviewed assignment policy are completed.

Applicability rules are kept separate: Euler--Bernoulli uses the thin-rod
diameter criterion, while Timoshenko uses the `Omega/Omega_c` cut-off audit.
Do not add vertical applicability lines to these plots. Timoshenko curves must
not be dashed because of the Euler--Bernoulli thin-rod criterion; only EB
curves may be dashed where the EB diameter criterion fails.

## What Must NOT Be Claimed Yet

- "3D FEM fully validates the coupled Timoshenko model."
- "Fused-cylinder disagreement disproves the Timoshenko coupling."
- "Sorted FEM mode number equals analytic mode number."
- "All modes are reliably identified."

## Remaining Limitations

1. Mesh convergence has been started for the point-joint coupled FEM workflow:
   beta=15, `epsilon = 0.025, 0.05`, and mesh factors `1.5, 1.0, 0.75` are
   recorded in `results/coupled_equal_rods_point_joint_mesh_convergence_report.md`.
   Extension beyond this beta/epsilon diagnostic is still pending.
2. MAC matching is currently based on centerline displacement only.
3. Duplicate best-solid-mode matches exist.
4. Weak MAC rows remain.
5. Preliminary Hungarian assignment rows are now written for the convergence
   diagnostic, but the article workflow still needs reviewed unique assignment
   policy before final tables.
6. Representative strong, moderate, and weak mode shapes need visual
   inspection.
7. The planar-constrained diagnostic required an MPC-compatible fallback
   because direct `ALL_SOLID_NODES, 3, 3, 0` constraints conflict with rigid
   body dependent nodes in CalculiX; this constraint implementation must be
   described clearly if used.
8. The frequency effect of the CalculiX `*RIGID BODY` constraint should be
   checked, including possible artificial stiffness or mass effects.
9. Article-grade conclusions must remain separated from diagnostic
   observations.
10. The central-patch joint variants still need mesh refinement, visual
    mode-shape review, and a clearer mechanical interpretation of moment
    transfer before they can be promoted beyond diagnostic use.

## Article-Ready Validation Roadmap

### Step 1. Point-Joint Mesh Convergence

- Start with `beta = 15 deg`.
- Use `epsilon = 0.025` and `epsilon = 0.05`.
- Status: first diagnostic run completed for mesh factors `1.5`, `1.0`, and
  `0.75`; strong/moderate final-mesh rows passed the current stability screen.
- Remaining: inspect weak/duplicate cases and decide whether additional
  refinement or a stricter unique-assignment policy is needed before article
  use.

### Step 2. Unique MAC Assignment

- Add Hungarian assignment for analytic modes 1--6.
- Report strong, moderate, and weak matches separately.
- Keep duplicate raw best-match warnings visible as diagnostics.

### Auxiliary. Planar-Constrained Point-Joint Check

- Status: beta=15 diagnostic completed for `epsilon = 0.025` and `0.05`.
- Purpose: check whether suppressing out-of-plane/torsion-like 3D families
  clarifies comparison with the planar EB/Timoshenko model.
- Limitation: this is not a full 3D validation replacement and article use
  requires clear labeling as a constrained planar-subspace diagnostic.

### Auxiliary. Equal-Thickness `Lambda(mu)` FEM Overlay

- Status: corrected beta=15, `epsilon = 0.01`, `eta = 0` diagnostic written
  with analytic descendant branches and sparse planar point-joint FEM points.
- Purpose: inspect mode identity and 3D-vs-1D deviations across `mu` while
  keeping EB diameter validity, Timoshenko cut-off validity, mesh quality, and
  MAC strength separate.
- Limitation: large FEM discrepancies remain diagnostic; cases with negative
  Jacobian mesh warnings or weak MAC are excluded from the main FEM points.
- Follow-up audit: the `mu=0`, `epsilon=0.01` sanity audit confirms the
  `Lambda_FEM=sqrt(Omega/epsilon)` conversion and the 1D EB frame-FEM anchor,
  while leaving the planar 3D solid discrepancy unexplained pending
  constraint-isolation or full point-joint `epsilon=0.01` checks.
- Constraint-isolation follow-up:
  `audit_3d_fem_point_joint_mu0_eps0p01.py` now compares the current
  planar-constrained point-joint, the full point-joint without planar
  constraints, and a straight two-half-rod point-joint sanity case at
  `epsilon=0.01`. The discrepancy is already present in the full beta=15
  point-joint case, while the straight joint sanity case has small strong-MAC
  errors. This points away from the planar constraint as the sole source and
  away from a generic rigid-reference-joint failure; the angled point-joint
  idealization or its 3D rigid end-face coupling remains diagnostic-only and
  unresolved.
- Coupling-model follow-up:
  `audit_3d_fem_joint_coupling_models.py` tests full rigid end faces and
  central rigid patches of radius `0.25*r` and `0.5*r` at
  `beta = 0, 15, 45, 90 deg`. The `0.25*r` patch reduces the beta=15 error
  but fails the beta=0 sanity check; the `0.5*r` patch passes beta=0 and
  reduces the beta=15 error relative to full rigid end faces. No tested model
  is classified as article-ready.
- Constraint-bracketing follow-up:
  `audit_3d_fem_joint_constraint_bracketing.py` keeps the above models and
  adds translation-only central/average equation probes, a reference-node
  rotation clamp, and a two-reference-node rigid-body attempt. It explicitly
  records that this is not patch-radius calibration. The beta=15 mean best
  relative errors among strong/moderate first-six matches are about 7.5% for
  `rigid_end_faces`, 4.7% for `center_patch_0p5`, 3.3% for
  `center_displacement_only`, and 3.1% for `equation_average_probe`; however,
  the translation-only probes fail beta=0 with about 7.3% mean error and
  missing moment transfer, so they are too-soft brackets rather than analytic
  candidates. The `rotation_clamped_reference` run matches the rigid-end-face
  errors rather than improving them, and is incompatible with the analytic
  free point-joint rotation.
- Rigid-joint thickness-trend follow-up:
  `audit_3d_fem_rigid_joint_thickness_trend.py` now uses the full rigid
  end-face point-joint without planar constraints, fused volume, patch
  coupling, or fitted parameters. The validation question is comparative:
  whether Timoshenko improves over Euler-Bernoulli as thickness increases for
  this engineering joint. The completed run gives eligible MAC counts
  `8,7,6,7` for `epsilon=0.005,0.01,0.025,0.05`, with Timoshenko/EB mean
  error ratios approximately `1.005,0.991,0.970,0.621`. Mesh convergence,
  visual mode-shape review, and article-grade assignment policy remain
  pending. The companion equal-thickness frequency plot visualizes the same
  rigid-joint benchmark for `epsilon=0.01,0.025,0.05` using the existing
  trend-run CSV.
- Rigid-joint equal-thickness `Lambda(mu)` follow-up:
  `plot_rigid_joint_equal_thickness_beta15_lambda_mu_eps_sweep.py` extends the
  comparative benchmark to three separate descendant-branch `Lambda(mu)` plots
  for `epsilon=0.01,0.025,0.05`, using a 16-point full rigid end-face FEM
  `mu` grid and no patch tuning. The main plotted FEM rows are restricted to
  strong/moderate non-ambiguous MAC matches; thicker cases show more excluded
  ambiguity rows but the eligible mean Timoshenko/EB error ratio decreases from
  about `1.016` to `0.752`.
- Analytic limit follow-up: `audit_variable_length_timoshenko_limits.py`
  verifies the variable-length Timoshenko determinant for `eta=0` sorted-root
  limit checks, and
  `audit_variable_length_thickness_timoshenko_limits.py` verifies the
  tau-aware eta extension for sorted-root limit checks. These audits should
  still be kept separate from FEM validation and do not promote any article
  figure.

### Step 3. Shape Figures

- Export representative 3D mode shapes.
- Compare them with the corresponding one-dimensional centerline shapes.
- Use only clear strong-MAC cases in article figures.

### Step 4. Article Validation Table

- Use only modes with strong, or at least moderate, MAC support.
- Report Euler--Bernoulli, Timoshenko, and 3D solid FEM frequencies together.
- Mark weak or ambiguous modes as excluded or diagnostic-only.

### Step 5. Extend To Thickness Mismatch

- Extend to `eta != 0` only after equal-rod point-joint validation is stable.
- Keep thickness-mismatch FEM runs separate from the equal-rod validation
  baseline.

## Current Outputs Index

| Workflow | Output |
| --- | --- |
| Single fixed-fixed rod report | `results/single_rod_fixed_fixed_3d_solid_fem_report.md` |
| Single fixed-fixed rod bending doublets | `results/single_rod_fixed_fixed_3d_solid_fem_bending_doublet_comparison.csv` |
| Fused coupled rods report | `results/coupled_equal_rods_3d_solid_fem_report.md` |
| Fused coupled rods in-plane comparison | `results/coupled_equal_rods_3d_solid_fem_in_plane_comparison_by_beta.csv` |
| Point-joint coupled rods report | `results/coupled_equal_rods_point_joint_3d_solid_fem_report.md` |
| Point-joint MAC matches | `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matches.csv` |
| Point-joint EB MAC matrix | `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_eb.csv` |
| Point-joint Timoshenko MAC matrix | `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_timoshenko.csv` |
| Point-joint mesh convergence report | `results/coupled_equal_rods_point_joint_mesh_convergence_report.md` |
| Point-joint planar-constrained report | `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_report.md` |
| Equal-thickness corrected `Lambda(mu)` FEM overlay | `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected_report.md` |
| Equal-thickness `mu=0` FEM sanity audit | `results/coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.md` |
| Point-joint `mu=0`, `epsilon=0.01` discrepancy isolation audit | `results/3d_fem_point_joint_mu0_eps0p01_audit.md` |
| Point-joint coupling-model audit | `results/3d_fem_joint_coupling_models_audit.md` |
| Rigid end-face thickness trend audit | `results/3d_fem_rigid_joint_thickness_trend_report.md` |
| Rigid end-face equal-thickness frequency trend plot | `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend_report.md` |
| Rigid end-face equal-thickness `Lambda(mu)` epsilon sweep | `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem_report.md` |
| Variable-length Timoshenko limit audit | `results/variable_length_timoshenko_limits_audit.md` |
