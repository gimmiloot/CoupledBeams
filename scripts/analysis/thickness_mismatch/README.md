# Thickness-Mismatch Script Map

This directory is a navigation layer for the diagnostic-only
thickness-mismatch study. The runnable scripts still live in the flat
`scripts/analysis/` layout for compatibility. Moving them in this pass would
require broad import and output-path churn without changing any diagnostic
result.

Project-wide rules for descendant branch identity, sorted-position metadata,
low-MAC assignments, thin-rod applicability, and diagnostic/article separation
live in `../../../docs/project_rules.md`.

## Preferred Entry Points

| Task | Preferred command | Compatibility or historical scripts | Main outputs | Notes |
| --- | --- | --- | --- | --- |
| Single fixed-fixed rod EB/Timoshenko check | `python scripts/analysis/compare_single_rod_eb_timoshenko.py` | none | `results/single_rod_fixed_fixed_eb_timoshenko_frequencies.csv`, `results/single_rod_fixed_fixed_eb_timoshenko_lambda_vs_epsilon.png`, `results/single_rod_fixed_fixed_eb_timoshenko_report.md` | Preferred diagnostic for the one-rod Timoshenko check before coupled-beam use; parameterized by epsilon list and marked with the `2r/L <= 0.1` thin-rod criterion. |
| Single fixed-fixed rod 3D solid FEM workflow | `python scripts/analysis/solid_fem_single_rod_fixed_fixed.py` | optional `GMSH_EXE`, `CCX_EXE` | `docs/thickness_mismatch/solid_fem_audit.md`, `results/single_rod_fixed_fixed_3d_solid_fem_report.md`, `results/single_rod_fixed_fixed_3d_solid_fem_comparison.csv`, `results/single_rod_fixed_fixed_3d_solid_fem_mode_metrics.csv`, `results/single_rod_fixed_fixed_3d_solid_fem_bending_doublet_comparison.csv`, `results/solid_fem_single_rod/` | Diagnostic-only external solid-FEM workflow seed. Writes Gmsh/CalculiX input templates, builds fixed-end node sets for ccx, parses `.frd` displacement modes when available, classifies bending/axial/torsion-like modes, and compares classified bending doublets against the one-rod EB/Timoshenko references; external solvers remain optional. |
| Coupled equal rods multi-beta 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods.py` | `solid_fem_coupled_equal_rods_beta15.py` remains for reproducing the original beta=15 diagnostic outputs | `results/coupled_equal_rods_3d_solid_fem_report.md`, `results/coupled_equal_rods_3d_solid_fem_geometry_audit.csv`, `results/coupled_equal_rods_3d_solid_fem_sorted_comparison_by_beta.csv`, `results/coupled_equal_rods_3d_solid_fem_mode_metrics_by_beta.csv`, `results/coupled_equal_rods_3d_solid_fem_in_plane_comparison_by_beta.csv`, `results/solid_fem_coupled_equal_rods/` | Diagnostic-only beta sweep for equal rods at `beta=15,45,90 deg`, `mu=0`, `eta=0`; audits fused-joint overlap using `r/sin(beta)`, keeps coordinate-derived outer fixed node sets, parses `.frd` mode shapes, and compares classified in-plane solid modes with beta-specific EB/Timoshenko sorted references. |
| Coupled equal rods point-joint 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods_point_joint.py` | fused workflow remains the comparison baseline | Baseline mode (`MESH_SIZE_FACTORS=[1.0]`): `results/coupled_equal_rods_point_joint_3d_solid_fem_report.md`, `results/coupled_equal_rods_point_joint_3d_solid_fem_geometry_audit.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_sorted_comparison_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mode_metrics_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_in_plane_comparison_by_beta.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_eb.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_timoshenko.csv`, `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matches.csv`, `results/solid_fem_coupled_equal_rods_point_joint/`; mesh-convergence mode: `results/coupled_equal_rods_point_joint_mesh_convergence_report.md`, `results/coupled_equal_rods_point_joint_mesh_convergence_frequencies.csv`, `results/coupled_equal_rods_point_joint_mesh_convergence_mac.csv`, `results/coupled_equal_rods_point_joint_mesh_convergence_summary.csv`, `results/solid_fem_coupled_equal_rods_point_joint_mesh_convergence/`; planar-constrained mode: `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_report.md`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_sorted_comparison.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mode_metrics.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mac_matrix_eb.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mac_matrix_timoshenko.csv`, `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_mac_matches.csv`, `results/solid_fem_coupled_equal_rods_point_joint_planar/` | Diagnostic-only alternative joint idealization: two separate cylinder meshes whose inner end-face node sets are coupled to one shared CalculiX rigid reference node after a local `*RIGID BODY` capability probe. Adds centerline MAC-like matching against EB/Timoshenko analytic shape reconstructions, which is preferred over simple classified in-plane ordering for mode-identity diagnostics. Optional in-script partial-run filters `RUN_BETA_DEG_VALUES` and `RUN_EPSILON_VALUES` restrict reruns; multiple `MESH_SIZE_FACTORS` route cases into separate `mesh_<factor>` convergence folders; planar-constrained mode is an auxiliary planar-subspace diagnostic and does not replace full 3D validation. |
| Equal-thickness coupled rods `Lambda(mu)` EB/Timoshenko/FEM diagnostic | `python scripts/analysis/plot_coupled_equal_rods_beta15_eps0p01_eb_timoshenko_fem_vs_mu.py` | earlier uncorrected outputs remain historical diagnostics | `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected.png`, `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected.csv`, `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected_report.md`, `results/solid_fem_coupled_equal_rods_beta15_eps0p01_planar_mu/` | Corrected diagnostic-only `beta=15 deg`, `epsilon=0.01`, `eta=0` plot over `mu=0..0.9`; `mu` changes rod lengths, so this is equal-thickness rather than equal-length except at `mu=0`. Analytic EB and Timoshenko curves are descendant branches; sparse FEM points use planar point-joint 3D solid FEM and MAC-like branch assignment. EB thin-rod diameter validity and Timoshenko cut-off validity are reported separately; do not dash Timoshenko curves using the EB diameter criterion. |
| Equal-thickness full rigid-joint `Lambda(mu)` epsilon sweep | `python scripts/analysis/plot_rigid_joint_equal_thickness_beta15_lambda_mu_eps_sweep.py` | does not replace the earlier planar `epsilon=0.01` overlay | `results/rigid_joint_equal_thickness_beta15_eps0p01_lambda_mu_eb_timoshenko_fem.png`, `results/rigid_joint_equal_thickness_beta15_eps0p025_lambda_mu_eb_timoshenko_fem.png`, `results/rigid_joint_equal_thickness_beta15_eps0p05_lambda_mu_eb_timoshenko_fem.png`, `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem.csv`, `results/rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem_report.md`, `results/solid_fem_rigid_joint_equal_thickness_beta15_lambda_mu/` | Current main diagnostic visualization for comparing EB, Timoshenko, and full rigid end-face 3D FEM along `mu` at `beta=15 deg`, `eta=0`, and `epsilon=0.01,0.025,0.05`. Analytic curves are descendant branches over `mu=0..0.9`; FEM uses a 16-point `mu` grid with no planar constraint, no fused volume, and no patch tuning. Main plots include only strong/moderate non-ambiguous MAC rows; weak, ambiguous, and duplicate matches remain in the CSV/report. |
| Variable-length Timoshenko limit audit | `python scripts/analysis/audit_variable_length_timoshenko_limits.py` | reuses `scripts/lib/variable_length_timoshenko.py`, copied from the corrected `Lambda(mu)` diagnostic formulas | `results/variable_length_timoshenko_limits_audit.csv`, `results/variable_length_timoshenko_limits_audit.md` | Diagnostic-only sorted-root audit for the variable-length Timoshenko implementation. It checks `mu=0` consistency against the older equal-rods Timoshenko diagnostic, the `epsilon -> 0` EB limit for several beta/mu values, and `beta=0` straight-rod mu-invariance. Current run passes all three checks for `eta=0`; that audit does not verify `eta != 0` tau scaling, descendant branch identity, FEM validation, or high-frequency cut-off behavior. |
| Equal-thickness `mu=0` FEM sanity audit | `python scripts/analysis/audit_coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity.py` | consumes the corrected plot CSV and generated planar FEM case | `results/coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.md`, `results/coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.csv` | Diagnostic audit of the `beta=15 deg`, `epsilon=0.01`, `eta=0`, `mu=0` discrepancy. Records the current variable-length Timoshenko formulas, checks `Lambda_FEM=sqrt(Omega/epsilon)`, audits point-joint planar constraints/node sets, compares sorted and MAC-order planar FEM against EB/Timoshenko, and uses the existing 1D EB frame FEM as a sorted-frequency sanity anchor. |
| Point-joint 3D FEM discrepancy isolation audit | `python scripts/analysis/audit_3d_fem_point_joint_mu0_eps0p01.py` | reuses point-joint Gmsh/CalculiX helpers and writes only isolated audit outputs | `results/3d_fem_point_joint_mu0_eps0p01_audit.md`, `results/3d_fem_point_joint_mu0_eps0p01_audit.csv`, `results/solid_fem_point_joint_mu0_eps0p01_audit/` | Diagnostic-only isolation pass for `beta=15 deg`, `mu=0`, `eta=0`, `epsilon=0.01`. It compares the current planar-constrained point-joint, the full point-joint without planar constraints, and a straight two-half-rod point-joint sanity case. Current result: the discrepancy is already present in full point-joint, while the straight joint sanity case has small strong-MAC errors; this points away from the planar constraint as the sole source and away from a generic rigid-reference-joint failure. |
| Point-joint coupling-model audit | `python scripts/analysis/audit_3d_fem_joint_coupling_models.py` | writes isolated Gmsh/CalculiX cases only under the audit output directory | `results/3d_fem_joint_coupling_models_audit.md`, `results/3d_fem_joint_coupling_models_audit.csv`, `results/solid_fem_joint_coupling_models_audit/` | Diagnostic-only comparison of alternative joint coupling idealizations at `epsilon=0.01`, `mu=0`, `eta=0`, and `beta=0,15,45,90 deg`. It tests full rigid end faces and central rigid patches with radii `0.25*r` and `0.5*r`; the `0.5*r` patch is the current best diagnostic candidate because it passes the beta=0 sanity check and reduces the beta=15 error, but no model is article-ready. |
| Point-joint constraint bracketing audit | `python scripts/analysis/audit_3d_fem_joint_constraint_bracketing.py` | writes isolated Gmsh/CalculiX cases only under the bracketing audit output directory | `results/3d_fem_joint_constraint_bracketing_audit.md`, `results/3d_fem_joint_constraint_bracketing_audit.csv`, `results/solid_fem_joint_constraint_bracketing_audit/` | Diagnostic-only bracketing study at `epsilon=0.01`, `mu=0`, `eta=0`, and `beta=0,15,45,90 deg`. It keeps the rigid-end-face and center-patch probes, adds translation-only central/average equation probes, a rotation-clamped reference-node bracket, and a two-reference-node constraint attempt. Current result: translation-only probes are too soft and fail beta=0 sanity, the rotation clamp reproduces rigid-end-face errors while remaining physically incompatible with the analytic free point-joint rotation, the two-reference-node deck is unsupported by the attempted CalculiX constraint combination, and `center_patch_0p5` remains diagnostic promising but not article-ready or patch-radius calibration. |
| Rigid end-face thickness trend audit | `python scripts/analysis/audit_3d_fem_rigid_joint_thickness_trend.py` | writes isolated full rigid end-face Gmsh/CalculiX cases only under the trend output directory | `results/3d_fem_rigid_joint_thickness_trend_report.md`, `results/3d_fem_rigid_joint_thickness_trend.csv`, `results/3d_fem_rigid_joint_thickness_trend_error_ratio.png`, `results/solid_fem_rigid_joint_thickness_trend/` | Diagnostic-only comparative trend study for the physically clear full rigid end-face point-joint at `beta=15 deg`, `mu=0`, `eta=0`, and `epsilon=0.005,0.01,0.025,0.05`. It uses no planar constraint, no fused volume, and no patch tuning; MAC-matched strong/moderate rows show the mean Timoshenko/EB error ratio decreasing from about `1.005` at `epsilon=0.005` to about `0.621` at `epsilon=0.05`, supporting the comparative thickness trend while leaving mesh convergence and visual mode review pending. |
| Rigid end-face equal-thickness frequency trend plot | `python scripts/analysis/plot_rigid_joint_equal_thickness_beta15_eps_trend.py` | consumes the existing rigid end-face trend CSV; does not rerun Gmsh/CalculiX | `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend.png`, `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend.csv`, `results/rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend_report.md` | Diagnostic-only visualization of the `beta=15 deg`, `mu=0`, `eta=0` rigid-joint thickness benchmark for `epsilon=0.01,0.025,0.05`. It plots EB, Timoshenko, and MAC-eligible full rigid end-face 3D FEM frequencies for the first six branches, records weak/ambiguous exclusions, and keeps EB diameter applicability separate from the Timoshenko cut-off without vertical applicability lines or Timoshenko style changes based on the EB criterion. |
| Coupled equal rods beta=15 3D solid FEM workflow | `python scripts/analysis/solid_fem_coupled_equal_rods_beta15.py` | optional `GMSH_EXE`, `CCX_EXE` | `results/coupled_equal_rods_beta15_3d_solid_fem_report.md`, `results/coupled_equal_rods_beta15_3d_solid_fem_sorted_comparison.csv`, `results/coupled_equal_rods_beta15_3d_solid_fem_mode_metrics.csv`, `results/coupled_equal_rods_beta15_3d_solid_fem_in_plane_comparison.csv`, `results/solid_fem_coupled_equal_rods_beta15/` | Diagnostic-only 3D solid check for beta=15, mu=0, eta=0 equal rods. Generates a fused two-cylinder mesh, coordinate-derived outer fixed node sets, CalculiX modal runs, `.frd` shape classification, and classified in-plane FEM vs EB/Timoshenko sorted-frequency comparisons. |
| Coupled equal rods EB/Timoshenko check | `python scripts/analysis/compare_coupled_equal_rods_eb_timoshenko.py` | none | `results/coupled_equal_rods_beta15_eb_timoshenko_frequencies.csv`, `results/coupled_equal_rods_beta15_timoshenko_kappa_sensitivity.csv`, `results/coupled_equal_rods_beta15_eb_timoshenko_lambda_vs_epsilon.png`, `results/coupled_equal_rods_beta15_eb_timoshenko_report.md` | Preferred diagnostic for the first coupled-rod Timoshenko sanity check at `beta=15 deg`, `mu=0`, `eta=0`; compares sorted frequencies only, marks `4*epsilon <= 0.1`, records cut-off margins, and runs a small kappa sensitivity check. |
| Variable-length thickness-mismatch Timoshenko tau audit | `python scripts/analysis/audit_variable_length_thickness_timoshenko_limits.py` | reuses `scripts/lib/variable_length_timoshenko.py` and the EB eta determinant in `src/my_project/analytic/formulas_thickness_mismatch.py` | `results/variable_length_thickness_timoshenko_limits_audit.csv`, `results/variable_length_thickness_timoshenko_limits_audit.md` | Diagnostic-only sorted-root audit for tau-aware eta support. It checks eta=0 regression against the legacy variable-length Timoshenko path, eta-to-zero continuity, `(mu, eta) -> (-mu, -eta)` swap symmetry, the `epsilon -> 0` EB thickness-mismatch limit, and beta=0 straight composite-rod behavior. Current run passes these gates for the tested first six sorted roots; descendant branch identity and energy plots remain separate. |
| Timoshenko energy partition audit | `python scripts/analysis/audit_timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5.py` | reuses `scripts/lib/variable_length_timoshenko.py` reconstruction and energy helpers | `results/timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5.csv`, `results/timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5_report.md`, diagnostic PNGs in `results/` | Diagnostic-only descendant audit at `beta=15 deg`, `epsilon=0.01`, `eta=0,0.5`, and `mu=0..0.9`. It partitions Timoshenko potential energy into bending/shear/axial terms per rod, records EB diameter and Timoshenko cutoff margins, and compares matching EB descendants. Use it only to explain why EB/Timoshenko frequencies can remain close when the EB diameter criterion fails locally; it does not change article figures, FEM validation, old determinants, old solvers, or baseline results. |
| Eta-zero and swap checks | `python scripts/analysis/check_thickness_mismatch_eta_zero_limit.py` | none | `results/thickness_mismatch_eta_zero_roots_check.csv`, `results/thickness_mismatch_swap_symmetry_check.csv` | Sanity checks for the diagnostic eta extension. |
| Descendant `Lambda(mu)` eta sweep | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py` | `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py`, `track_lambda_mu_thickness_mismatch_eta_sweep.py` | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants.png` and report | Current clean presentation-style plot; tracking warnings stay in reports/data, not as x-markers. |
| Eta=0.5 global spectrum overview | `python scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py` | none | `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.*` | Shows sorted roots, descendant branches, canonical sorted positions, and unresolved candidate assignments. |
| Refined branch-identity audit | `python scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py` | none | `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.*` | Local descendant audit for branches 5--7 at eta=0.5. |
| Eta=0.5 crossing audit | `python scripts/analysis/audit_eta0p5_eps0p0025_lambda_mu_crossings.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.*` | Pairwise diagnostic for real descendant-frequency crossings versus finite gaps, sorted-position metadata, and tracking artifacts on the dense `mu=0..0.9` grid. |
| Eta=0.5 analytic/FEM overlay | `python scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py` | `fem_check_thickness_mismatch_eta_p0p5_beta15.py` | `results/thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5.*` | Lightweight plotting comparison. The older FEM check keeps the detailed CSV/MAC diagnostics. |
| Eta=0.5 isolated-rod references, beta=15 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods.*` | Uses the clamped-supported / clamped-pinned reference convention and eta-normalized `Lambda_rod`. |
| Eta=0.5 isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods.*` | Same reference convention as beta=15, different beta only. |
| Eta=0.5 fixed-fixed isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed.*` | Uses clamped-clamped / fixed-fixed reference roots. |
| Branch-shape eta overlay | `python scripts/analysis/plot_thickness_mismatch_branch_shapes_vs_eta.py` | `plot_thickness_mismatch_branch5_shapes.py` | `results/thickness_mismatch_shapes_branch5_beta15_mu0p0_eta_m0p5_0_p0p5.png`, `results/thickness_mismatch_shapes_branch5_beta15_mu0p15_eta_m0p5_0_p0p5.png` | Parameterized descendant-shape overlay for several eta values at target mu values, using shape-MAC tracking from `mu=0`. |
| Branch-shape diagnostics | `python scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` | none | `results/thickness_mismatch_branch5_shapes_*.png` | Shape-only diagnostic for tracked branch 5. |

## Current FEM Validation Status

The current diagnostic FEM validation path has three levels:

- `solid_fem_single_rod_fixed_fixed.py` is the one-rod fixed-fixed 3D solid FEM
  benchmark and currently gives the clearest support for the Timoshenko
  correction through classified bending-doublet comparisons.
- `solid_fem_coupled_equal_rods.py` is the fused-cylinder coupled-rod geometry
  audit across `beta=15,45,90 deg`; it is useful diagnostically but is not a
  direct validation of the ideal point-joint one-dimensional model.
- `solid_fem_coupled_equal_rods_point_joint.py` is the current best coupled-rod
  validation path, using separate solid cylinders with a CalculiX rigid
  reference joint and MAC-like centerline shape matching.

The status, limitations, and article-ready roadmap are tracked in
`../../../docs/thickness_mismatch/fem_validation_status.md`.

## Timoshenko Diagnostic Notes

The single-rod and coupled equal-rod EB/Timoshenko diagnostics now record the
Timoshenko cut-off frequency

```text
Omega_c = sqrt(kappa*G*A/(rho*I))
```

For the coupled equal-rod project normalization `Omega = epsilon*Lambda^2`,
with `r = 2*epsilon` and `L_segment = 1`, the report also records

```text
Lambda_c = (kappa/(2*(1 + nu)))**0.25/epsilon
```

The current diagnostic Timoshenko extension modifies only the flexural part:
- bending moment: `M = E I psi'`
- shear force: `Q = kappa G A (w' - psi)`

The axial part remains classical:
- `N = E A u'`

This is acceptable for the current low flexural-frequency diagnostics, but
higher-frequency longitudinal refinements may require Mindlin-Herrmann /
related rod models.

For `Lambda(mu)` diagnostics, keep applicability criteria separated:
Euler--Bernoulli diameter validity uses `2*r_i/l_i <= 0.1`, while
Timoshenko validity is audited with `Omega/Omega_c`. Do not dash or exclude
Timoshenko curves solely because the Euler--Bernoulli diameter criterion fails.

The variable-length Timoshenko implementation is considered
diagnostic-verified only when the limit audit passes all three gates:
`mu=0` consistency with the older equal-rods Timoshenko diagnostic, the
`epsilon -> 0` Euler--Bernoulli limit, and `beta=0` straight-rod
mu-invariance. The current audit output passes these gates for `eta=0` sorted
roots.

The tau-aware `eta != 0` extension in `scripts/lib/variable_length_timoshenko.py`
is considered diagnostic-verified only for the sorted-root checks in
`results/variable_length_thickness_timoshenko_limits_audit.md`: eta=0
regression, eta-to-zero continuity, rod-swap symmetry, the `epsilon -> 0`
Euler--Bernoulli thickness-mismatch limit, and beta=0 straight composite-rod
behavior. Descendant branch identity, energy plots, and FEM validation remain
separate.

## Historical And One-Off Diagnostics

- `plot_lambda_eta_thickness_mismatch.py` is the initial sorted-root
  `Lambda(eta)` diagnostic.
- `track_lambda_eta_thickness_mismatch.py` is the initial eta-continuation
  diagnostic.
- `track_lambda_mu_thickness_mismatch_eta_sweep.py` remains useful as an
  editable single-beta tracking sandbox, but the current clean eta sweep plot
  is `plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py`.
- `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py` is kept
  as the earlier large-eta diagnostic with warning artifacts.
- `fem_check_thickness_mismatch_eta_p0p5_beta15.py` is the heavier
  analytic-vs-FEM diagnostic with detailed CSV/MAC output.
- `plot_thickness_mismatch_branch5_shapes.py` is the older branch-5 shape
  diagnostic. Prefer `plot_thickness_mismatch_branch_shapes_vs_eta.py` when
  comparing several eta values on one shape plot or changing branch/mu/eta
  parameters.

## Helper Sources

- `src/my_project/analytic/formulas_thickness_mismatch.py` implements the
  diagnostic determinant and eta model helpers.
- `scripts/lib/thickness_mismatch_mac_tracking.py` implements diagnostic
  descendant tracking and unresolved-assignment flags.
- `scripts/lib/thickness_mismatch_diagnostic_helpers.py` collects shared
  plotting, validity, and isolated-rod reference utilities.
- `scripts/lib/variable_length_timoshenko.py` contains the diagnostic
  variable-length Timoshenko matrix helpers used by the eta=0 and tau-aware
  eta audits; it does not replace the old Euler--Bernoulli determinant or old
  solvers.

The CS/CP reference convention uses roots of `tan(alpha)=tanh(alpha)`. The
CC/FF reference convention uses roots of `cosh(alpha) cos(alpha)=1`. Both use
the eta-model normalization
`Lambda_rod_i = alpha_n*sqrt(tau_i)/(1 +/- mu)`.

## Future Refactor TODO

When a broader import cleanup is explicitly requested, the flat scripts can be
converted into parameterized entry points under this directory:

- `plot_descendant_eta_sweep.py`
- `plot_with_isolated_rods.py`
- `compare_fem.py`
- `audit_branch_identity.py`
- `check_eta_zero_limit.py`

The old script names should then remain as compatibility wrappers that print a
short deprecation note and call the new entry point with the old parameters and
output paths.
