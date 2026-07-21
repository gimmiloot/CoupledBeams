# Thickness-Mismatch Script Map

This directory is a navigation layer for the diagnostic-only
thickness-mismatch study. The runnable scripts still live in the flat
`scripts/analysis/` layout for compatibility. Moving them in this pass would
require broad import and output-path churn without changing any diagnostic
result.

Project-wide rules for descendant branch identity, sorted-position metadata,
low-MAC assignments, thin-rod applicability, and diagnostic/article separation
live in `../../../docs/project_rules.md`.

The current positive-gap status for frequency-crossing checks is summarized in
`../../../docs/thickness_mismatch/frequency_crossing_verification_status.md`.

The out-of-plane EB + Saint-Venant torsion subsystem is documented, covered by
lightweight determinant sanity tests, and now has a diagnostic 1D EB+torsion
FEM validation audit listed below.

## Directory Scaffold

New thickness-mismatch diagnostic scripts should be added under this directory
by diagnostic role:

```text
maps/
  eta-mu-beta frequency maps, heatmaps, sensitivity maps

shapes/
  mode-shape plotting, sorted-mode shapes, descendant-branch shapes

audits/
  local checks, sorted-vs-descendant audits, gap checks

lambda_mu/
  Lambda(mu) diagnostic plots, single-rod spectra overlays

postprocess/
  CSV-only summaries and global trend analysis
```

Existing flat scripts in `scripts/analysis/` remain supported for
compatibility. Future migration should use thin wrappers at the old paths, not
breaking moves.

Diagnostic scripts must explicitly state whether they use sorted frequencies
or descendant branches. Sorted frequencies are used for spectral maps and gap
metrics; descendant branches are used for tracking modal character.

## Smoke Mode Convention

Heavy diagnostic entry points should support `--smoke` or document an existing
`--quick` mode as the smoke mode. Smoke mode means a tiny grid, very fast
runtime, output under `results/_smoke/` or `results/smoke/`, and no overwrite
of normal diagnostic outputs unless a future script explicitly implements and
requires `--force`.

The current eta-mu-beta map, out-of-plane mode-character beta map, sorted
`Lambda(mu)` reference plot, descendant mode-shape beta scan, fixed sorted-mode
beta scan, EB/Timoshenko `Lambda(beta)` comparison, EB/Timoshenko `Lambda(mu)`
comparison, and full-spectrum analytic-vs-3D-FEM scaffold all expose `--smoke`.

New beta=0 EB/Timoshenko/3D FEM validation wrappers under `audits/` are
diagnostic-only stable entry points. They reuse the existing analytic and
single-cylinder/stepped-cylinder FEM helpers, and are intentionally separated
when the input/output contract differs: e.g. a uniform-cylinder first-8
low-spectrum plus bending-pair comparison is not the same workflow as a
two-radius stepped-cylinder audit.

## Mode-Shape Refactor TODO

`plot_mode_shapes_eta_beta_scan.py` and
`plot_mode_shapes_eta_beta_scan_sorted_modes.py` should eventually become one
parameterized entry point with `selection = sorted | descendant`. Shared
functionality to extract later includes beta grid construction, root solving,
coefficient extraction, geometry reconstruction, sign normalization,
grid/contact-sheet plotting, and CSV summary writing.

Example future descendant command:

```bash
python scripts/analysis/thickness_mismatch/shapes/plot_beta_scan.py \
  --eta 0.1 \
  --mu 0.3 \
  --epsilon 0.0025 \
  --beta-start 0 \
  --beta-end 11 \
  --selection descendant \
  --branch 2
```

Example future sorted-frequency command:

```bash
python scripts/analysis/thickness_mismatch/shapes/plot_beta_scan.py \
  --eta 0.1 \
  --mu 0.3 \
  --epsilon 0.0025 \
  --beta-start 0 \
  --beta-end 15 \
  --selection sorted \
  --sorted-modes 2 3
```

## Current Research Direction

The current direction is a conservative `K = 10` EB-based certificate for the
safe prefix of the sorted spectrum. The CSV-first postprocessor is implemented;
the engineering target, candidate EB-only indicators, held-out-geometry
validation policy, and operation-count requirements are recorded in the
[EB safe-spectrum-prefix research plan](../../../docs/thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md).

Generate the required K=10 source inputs without changing the legacy defaults:

```bash
python scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py \
  --n-reported-modes 10 \
  --n-candidate-roots 16

python scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py \
  --n-reported-modes 10 \
  --n-candidate-roots 20
```

Then run the CSV-only certification workflow:

```bash
python scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py \
  --stage1-dir results/eb_validity_vs_timoshenko_stage1 \
  --fixed-epsilon-dir results/eb_validity_fixed_epsilon_geometry_scan \
  --output-dir results/eb_safe_prefix_certification \
  --k-max 10
```

The postprocessor rejects K=8-only inputs with the two regeneration commands
above. It does not solve roots, run FEM, or use Timoshenko-derived quantities as
certificate features.

The completed 21-geometry pre-solution pilot is a separate, manifest-driven
diagnostic. It generates only the selected source points and then reruns all
calibration, fold, plot, and report work from CSV:

```bash
python scripts/analysis/thickness_mismatch/audits/run_eb_epsilon_apriori_pilot.py \
  --manifest scripts/analysis/thickness_mismatch/audits/data/eb_epsilon_apriori_pilot_cases.csv \
  --output-dir results/eb_epsilon_apriori_pilot \
  --k-max 10 \
  --n-candidate-roots 20 \
  --reuse-cache

python scripts/analysis/thickness_mismatch/postprocess/analyze_eb_epsilon_apriori_pilot.py \
  --pilot-dir results/eb_epsilon_apriori_pilot \
  --output-dir results/eb_epsilon_apriori_pilot/analysis \
  --k-max 10 \
  --frequency-error-threshold 0.10 \
  --candidate-grid-sizes 8 16 32
```

The pilot compares `epsilon_0` and helper-derived `epsilon_max` with Rules
A--D and Rule A-gap. It is diagnostic evidence from 21 selected geometries,
not a geometry-only replacement for the EB-based certificate or a
continuous-domain guarantee.

The refined research-step-2 straight-system baseline is generated with:

```bash
python scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_baseline_thresholds.py \
  --epsilon-min 0.005 \
  --epsilon-max 0.060 \
  --primary-epsilon-min 0.010 \
  --primary-epsilon-max 0.050 \
  --coarse-step 0.0005 \
  --k-max 10 \
  --n-spectrum-roots 12 \
  --n-candidate-roots 20 \
  --epsilon-tolerance 1e-6 \
  --delta-tolerance 1e-6 \
  --output-dir results/eb_epsilon_baseline_thresholds \
  --reuse-cache
```

It computes first-loss prefix thresholds only for
`(beta, mu, eta)=(0,0,0)`. The corrected source
`factorized_straight_spectrum_v2` retains the exact axial family and solves the
exact bending block extracted from the unchanged general 6x6 Timoshenko
matrix. The raw general sign scan is kept only as a completeness audit because
it can miss two close roots inside one scan interval. The workflow independently
force-recomputes resolved brackets, writes 14 audit CSVs, one Markdown report,
and six compact diagnostic plots, and supports `--smoke` and `--plot-only`.
Old pre-correction outputs/cache are preserved in named legacy directories and
cannot be loaded as the corrected cache. The corrected cache key includes the
algorithm version, geometry, requested/candidate root counts, and audit-local
solver settings. The completed audit has 23592/23592 passing factorized rows;
the 155 raw-general omissions are retained as evidence and all are confirmed
by local full-matrix SVD. Prefix 1 remains safe through the
upper scan endpoint; prefixes 2--10 have verified first-loss brackets. This
workflow is not a lower-envelope proof and does not implement step 3.

## Preferred Entry Points

Stage-1 universal-parameter post-processing is available after the EB
applicability scan has produced its CSV/cache outputs:

```bash
python scripts/analysis/thickness_mismatch/audits/analyze_universal_eb_validity_parameters_stage1.py
```

It writes `universal_parameter_*` CSV/report files and PNG diagnostics under
`results/eb_validity_vs_timoshenko_stage1/`, comparing `epsilon_max`,
`chi_max`, `chi_eff`, and the EB-mode estimate
`Pi_EB = Pi_shear + Pi_rotary`. The script reconstructs EB modes from saved
`Lambda_EB` values and reads Stage-1 cache data when critical-threshold rows
are available; it is diagnostic post-processing and does not recompute roots
when the Stage-1 data are present.

| Task | Preferred command | Compatibility or historical scripts | Main outputs | Notes |
| --- | --- | --- | --- | --- |
| Out-of-plane EB+torsion sorted `Lambda(beta)` maps | `python scripts/analysis/thickness_mismatch/maps/plot_out_of_plane_lambda_beta_mu_eta.py` | none | `results/out_of_plane_lambda_beta_mu_eta/out_of_plane_lambda_beta_mu_eta_sorted_roots.csv`, `results/out_of_plane_lambda_beta_mu_eta/out_of_plane_lambda_beta_gap_summary.csv`, low/extended eta-panel PNGs, overview PNG, and report | Diagnostic-only sorted-frequency beta maps for the out-of-plane Euler--Bernoulli plus Saint-Venant torsion determinant across selected `mu` and eta values. It includes low-spectrum first roots, extended plots reaching the first torsion-root region, and adjacent sorted-gap summaries. It does not use descendant tracking, crossing verification, 3D FEM, or article styling. |
| In-plane vs out-of-plane analytic `Lambda(beta)` maps | `python scripts/analysis/thickness_mismatch/maps/plot_in_plane_vs_out_of_plane_lambda_beta.py` | none | `results/in_plane_vs_out_of_plane_lambda_beta/in_plane_vs_out_of_plane_lambda_beta_sorted_roots.csv`, `in_plane_vs_out_of_plane_lambda_beta_case_summary.csv`, `in_plane_vs_out_of_plane_lambda_beta_spike_audit.csv`, `full_analytic_union_lambda_beta.csv`, per-case PNGs, overview PNG, and report | Diagnostic-only sorted-frequency comparison of the in-plane EB determinant roots against the out-of-plane EB+torsion determinant roots. In-plane curves are solid, out-of-plane curves are dashed, and spike handling uses CSV-first raw roots plus NaN plot segmentation without descendant tracking, smoothing, FEM, or article styling. |
| In-plane EB vs Timoshenko sorted `Lambda(beta)` cases | `python scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_beta_cases.py` | reuses `formulas_thickness_mismatch.py` and `scripts/lib/variable_length_timoshenko.py` | `results/eb_vs_timoshenko_lambda_beta_cases/lambda_beta_eb_vs_timo_*.png`, matching per-case CSVs, `eb_vs_timo_lambda_beta_case_summary.csv`, `eb_vs_timo_lambda_beta_spike_audit.csv`, `eb_vs_timo_lambda_beta_timing_report.csv`, cache NPZ files under `cache/`, `eb_vs_timo_lambda_beta_overview.png`, and `eb_vs_timo_lambda_beta_cases_report.md` | Diagnostic-only sorted-frequency comparison for `(mu, eta)=(0,0),(0,0.1),(0.3,0.1)` and `epsilon=0.0025,0.05` on `beta=0..90 deg`. Timoshenko curves are solid, Euler-Bernoulli curves are dashed, same sorted index uses the same color, and the optimized workflow supports cache reuse, plot-only regeneration, selected-case filters, timing output, continuation-assisted Timoshenko solves, local beta refinement, dense retry, raw/plot CSV columns, and targeted row-normalized smallest-singular-value recovery for missed close roots without descendant tracking, FEM/3D FEM, article styling, determinant changes, or shear-coefficient changes. |
| In-plane EB vs Timoshenko sorted `Lambda(mu)` eps scan | `python scripts/analysis/thickness_mismatch/maps/plot_eb_vs_timoshenko_lambda_mu_cases.py --beta-deg 45 --eta 0 --epsilon-values 0.005 0.01 0.02 0.03 0.04 --mu-min 0 --mu-max 0.7 --mu-step 0.005 --n-roots 6 --output-dir results/eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan` | reuses the EB-vs-Timoshenko `Lambda(beta)` diagnostic solve wrappers, `formulas_thickness_mismatch.py`, and `scripts/lib/variable_length_timoshenko.py` | `results/eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan/lambda_mu_eb_vs_timo_beta45_eta0_eps*.png`, matching per-case CSVs, `eb_vs_timo_lambda_mu_beta45_eta0_eps_scan_summary.csv`, `eb_vs_timo_lambda_mu_beta45_eta0_eps_scan_spike_audit.csv`, `timing_report.csv`, cache NPZ files under `cache/`, `eb_vs_timo_lambda_mu_beta45_eta0_eps_scan_overview.png`, and `eb_vs_timo_lambda_mu_beta45_eta0_eps_scan_report.md` | Diagnostic-only sorted-frequency comparison for `beta=45 deg`, `eta=0`, and `epsilon=0.005,0.01,0.02,0.03,0.04` on `mu=0..0.7`. Timoshenko curves are solid, Euler-Bernoulli curves are dashed, same sorted index uses the same color, and the workflow supports `--beta-deg`, cache reuse, plot-only regeneration, timing output, continuation-assisted Timoshenko solves, spike audit, local mu refinement only when needed, dense retry/local SVD fallback only on suspicious points, and requested raw/plot CSV columns without single-rod references, descendant tracking, FEM/3D FEM, article styling, determinant changes, or shear-coefficient changes. |
| Stage-1 EB applicability relative to Timoshenko | `python scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py` | reuses the existing in-plane EB root solver, Timoshenko root solver, sorted shape/energy reconstruction helpers, shape MAC helpers, and `thickness_mismatch_factors` | `results/eb_validity_vs_timoshenko_stage1/eb_timo_mode_level_metrics.csv`, `eb_timo_validity_summary.csv`, `eb_timo_critical_thickness_by_branch.csv`, `epsilon_branch_tracking_audit.csv`, `critical_threshold_refinement_audit.csv`, `chi_scaling_fit_summary.csv`, diagnostic PNGs, root cache under `cache/`, timing CSVs, and `eb_validity_vs_timoshenko_stage1_report.md` | Diagnostic-only `beta=45 deg`, `eta=0` applicability scan over `mu=0..0.7` and `epsilon_0=0.0025..0.06`. It keeps sorted-spectrum metrics separate from physical-branch metrics, seeds `base_branch_1..8` at the thin limit, continues EB and Timoshenko branches independently by shape MAC over candidate roots 1..12, computes energy character and chi metrics, refines 10% `delta_f` threshold brackets, and supports `--smoke`, `--reuse-cache`, `--force-recompute`, `--plot-only`, `--skip-critical-refinement`, grid overrides, and output-dir overrides. It does not run FEM/3D FEM/Gmsh/CalculiX, change formulas/determinants/root solvers/k', or touch article workspaces. |
| Fixed-epsilon EB applicability geometry scan | `python scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py` | reuses the existing EB/Timoshenko root wrappers, tau-aware geometry factors, analytic shape reconstruction, energy helpers, and shape-MAC utilities | `results/eb_validity_fixed_epsilon_geometry_scan/fixed_epsilon_mode_metrics.csv`, `fixed_epsilon_matching_audit.csv`, `fixed_epsilon_point_summary.csv`, predictor fit/CV/threshold/runtime CSVs, diagnostic PNGs, root cache under `cache/`, and `eb_validity_fixed_epsilon_geometry_scan_report.md` | Diagnostic-only follow-up scan at fixed `epsilon_0=0.02` over `beta=0..90 deg`, `mu=0..0.7`, and `eta=-0.1,0,0.1`. It keeps sorted metrics separate from per-point homologous EB/Timoshenko matching, uses mass-weighted local `[u,w]` MAC with one-to-one assignment and close-cluster subspace rows, computes EB-only `Theta_max_EB` and `Pi_EB` predictors, supports cache reuse, `--smoke`, `--plot-only`, optional local Timoshenko benchmarking, and grid overrides. It does not use descendant tracking, FEM/3D FEM/Gmsh/CalculiX, article styling, determinant changes, root-solver changes, or shear-coefficient changes. |
| K=10 EB safe-prefix certification | `python scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py --stage1-dir results/eb_validity_vs_timoshenko_stage1 --fixed-epsilon-dir results/eb_validity_fixed_epsilon_geometry_scan --output-dir results/eb_safe_prefix_certification --k-max 10` | consumes only the sorted rows from the two applicability CSV files and reuses fixed-epsilon EB reconstruction/predictor helpers | `results/eb_safe_prefix_certification/eb_safe_prefix_mode_audit.csv`, geometry, calibration, fold-prediction, validation, exclusion, overlap, predictor-consistency, and operation-count CSVs, plus `eb_safe_prefix_certification_report.md` | Diagnostic-only geometry-level certification with Rules A--D plus Rule A-gap, train-only deterministic threshold calibration, prefix-extrema-aware candidate priorities, complete-geometry held-out splits, explicit false-safe/conservative-loss separation, and raw operation counts. Requires complete K=10 inputs and performs no root recomputation, FEM, local Timoshenko solve, geometry-only prediction, or article work. |
| Manifest-driven EB epsilon a-priori pilot source generation | `python scripts/analysis/thickness_mismatch/audits/run_eb_epsilon_apriori_pilot.py --manifest scripts/analysis/thickness_mismatch/audits/data/eb_epsilon_apriori_pilot_cases.csv --output-dir results/eb_epsilon_apriori_pilot --k-max 10 --n-candidate-roots 20 --reuse-cache` | reuses the fixed-epsilon EB/Timoshenko root cache, point solve, mode reconstruction, predictor, MAC, cluster, and boundary-retry helpers plus the Stage-1 local-thickness helper | resolved manifest, K=10 mode and geometry CSVs, matching/exclusion/root-operation audits, cache files, and `epsilon_pilot_generation_report.md` under `results/eb_epsilon_apriori_pilot/` | Diagnostic-only solve of exactly 21 selected geometries. Root 11 is retained only as the right EB gap guard for target mode 10. No production tensor grid, FEM, solver modification, or article workflow is involved. |
| Geometry-only epsilon pilot post-processing | `python scripts/analysis/thickness_mismatch/postprocess/analyze_eb_epsilon_apriori_pilot.py --pilot-dir results/eb_epsilon_apriori_pilot --output-dir results/eb_epsilon_apriori_pilot/analysis --k-max 10 --frequency-error-threshold 0.10 --candidate-grid-sizes 8 16 32` | consumes only the pilot CSV files and imports production safe-prefix rule/calibration semantics | exact E0/Emax thresholds, fold predictions and summaries, reference/matched/same-epsilon audits, operation and convergence CSVs, six compact PNGs, and `epsilon_apriori_pilot_report.md` | CSV-only comparison of E0-ref, Emax-ref, E0-cal, Emax-cal, Rules A--D, and Rule A-gap. Complete geometries stay together in every fold; X-domain fallbacks and 16-to-32 search changes are explicit. The 21-case result does not establish a universal geometry-only criterion. |
| Straight-system epsilon baseline thresholds | `python scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_baseline_thresholds.py --epsilon-min 0.005 --epsilon-max 0.060 --primary-epsilon-min 0.010 --primary-epsilon-max 0.050 --coarse-step 0.0005 --k-max 10 --n-spectrum-roots 12 --n-candidate-roots 20 --epsilon-tolerance 1e-6 --delta-tolerance 1e-6 --output-dir results/eb_epsilon_baseline_thresholds --reuse-cache` | reuses the existing matrix/basis helpers and raw general root cache for comparison; corrected source uses `scripts/lib/straight_rod_factorized_spectrum.py` | 14 `baseline_*.csv` audit tables, `eb_epsilon_baseline_thresholds_report.md`, six diagnostic PNGs, versioned corrected cache, and preserved legacy outputs/cache under `results/eb_epsilon_baseline_thresholds/` | Research-step-2 baseline only at `(beta,mu,eta)=(0,0,0)`. Uses the exact EB family union and exact Timoshenko axial/4x4-bending factorization with multiplicity, while retaining the general 6x6 sign scan as an independent completeness audit. Computes corrected first-loss brackets, re-entry/order/multiplicity audits, mu-invariance, R1--R3 close-root regressions, full-matrix SVD confirmations, conservative display floors, and force-recompute verification. It does not change formulas/shared solvers, run FEM, prove a geometry lower envelope, or implement step 3. |
| Longitudinal-character audit for EB/Timoshenko `Lambda(mu)` suspect modes | `python scripts/analysis/thickness_mismatch/audits/audit_longitudinal_suspect_modes_eb_timo.py` | reuses the existing sorted EB/Timoshenko `Lambda(mu)` CSVs when available, the EB eta determinant for coefficients, and `scripts/lib/variable_length_timoshenko.py` for Timoshenko fields and section/kappa data | `results/eb_vs_timoshenko_longitudinal_suspect_modes/suspect_mode_energy_fractions.csv`, `control_mode_energy_fractions.csv`, `suspect_mode_shape_summary.csv`, `timoshenko_joint_continuity_audit.csv`, per-mode shape PNGs, `suspect_shapes_eps*_grid.png`, `longitudinal_suspect_modes_report.tex`, and optional PDF | Diagnostic-only sorted-frequency mode-shape and energy audit for the beta=45 deg, eta=0 suspect cases `epsilon=0.03` sorted 5 and `epsilon=0.05` sorted 4. It computes axial, bending, and Timoshenko shear energy fractions plus displacement fractions, audits Timoshenko joint compatibility before plotting, blocks Timoshenko shape figures with normalized kinematic joint gap above `1e-6`, includes neighboring sorted modes at `mu=0.35` as controls, and writes a Russian TeX report. It does not use descendant tracking, FEM/3D FEM/Gmsh/CalculiX, article outputs, determinant changes, root-solver changes, or shear-coefficient changes. |
| Timoshenko thin-limit shape-bug audit | `python scripts/analysis/thickness_mismatch/audits/audit_timoshenko_shape_bug_thin_limit.py` | reuses the sorted EB/Timoshenko root provider, analytic EB reconstruction, `scripts/lib/variable_length_timoshenko.py`, and `scripts/lib/in_plane_shape_geometry.py` | `results/timoshenko_shape_bug_audit/beta_unit_audit.md`, `display_transform_regression.csv`, `thin_limit_eb_timo_shape_mac.csv`, `timoshenko_null_vector_audit.csv`, `timoshenko_shape_residual_audit.csv`, three thin-limit full-centerline PNG grids, and `timoshenko_shape_bug_audit_report.md` | Diagnostic-only sorted-mode audit at `epsilon=0.0025` and comparison with `epsilon=0.03` for beta=45 deg, eta=0, mu=0,0.2,0.4,0.6, and modes 4-6. It checks explicit beta units, same-index local and corrected-display EB/Timoshenko MAC, rod-2-only display MAC, shear convergence, both smallest singular values, null residuals, clamp/joint residuals, and plot regularity. Corrected Timoshenko display normals are reflected relative to determinant components; older Timoshenko global plots are invalid for physical interpretation. It does not use descendant tracking, FEM/3D FEM/Gmsh/CalculiX, article outputs, determinant changes, root-solver changes, or shear-coefficient changes. |
| Clean corrected EB/Timoshenko full shape grids | `python scripts/analysis/thickness_mismatch/shapes/plot_eb_timo_full_mode_shapes_eps0p03_beta45_eta0_modes4_6.py --epsilon 0.03 --mode-indices 4 5 6 --mu-values 0 0.2 0.4 0.6 --clean-style --output-dir results/eb_timo_clean_mode_shapes` | reuses the sorted shape-construction audit and shared in-plane display geometry | `results/eb_timo_clean_mode_shapes/` clean EB, Timoshenko, and optional combined PNG grids, `clean_mode_shape_summary.csv`, and `clean_mode_shape_report.md` | Diagnostic-only full-displacement centerlines from sorted frequencies. The entrypoint also accepts `--epsilon 0.04 --mode-indices 3 4 5`, `--epsilon 0.02 --mode-indices 5 6 7`, `--epsilon 0.05 --mode-indices 3 4 5`, and other sorted-mode/mu selections. Scale 0.08 is checked across all EB+Timoshenko shapes for the requested epsilon; if any shape fails, the whole epsilon set checks 0.05, then 0.02, and uses one common scale. Panel titles are clean (`mu`, sorted `k`, `Lambda`) and the corrected in-plane display mapping is used. |
| Timoshenko modes 4-6 corrected vector/component diagnostics | `python scripts/analysis/thickness_mismatch/audits/audit_timoshenko_modes_4_6_shape_diagnostics.py` | shared sorted reconstruction and in-plane display geometry | `results/timoshenko_mode_shape_diagnostics/` corrected EB/Timoshenko displacement-vector grids, component grids, MAC/energy/reconstruction CSVs, conservative full-centerline grid, and report | Diagnostic-only component and vector-field interpretation of the same epsilon=0.03 grid; corrected vector plots use theory-specific determinant-to-display mappings and old Timoshenko vectors are invalidated. |
| Beta=90 eta-scan `Lambda(mu)` with single-rod references | `python scripts/analysis/thickness_mismatch/maps/plot_lambda_mu_beta90_eta_scan_with_single_rod_refs.py` | single-eta predecessor: `plot_lambda_mu_beta90_eta0p1_with_single_rod_refs.py`; helpers: `formulas_thickness_mismatch.py`, `FreqMuNet.py`, `solvers.py` | Per-eta folders under `results/lambda_mu_beta90_eta_scan_single_rod_refs/eta_p0p1/`, `eta_0p0/`, and `eta_m0p1/` containing sorted-system, reference, nearest-match, all-candidate, branch-family summary CSVs, readable/full PNGs, and reports; plus `eta_comparison_summary.csv` and `eta_scan_report.md` | Diagnostic-only sorted-frequency `Lambda(mu)` maps for `beta=90 deg`, `epsilon=0.0025`, `eta=0.1,0,-0.1`, and `mu=0..0.8`, overlaid with mu-dependent isolated-rod CP/FP and CC/FF bending references using `Lambda_ref=alpha*sqrt(tau_i)/L_i`. The readable PNGs default to the fixed range `0 <= Lambda <= 13` (`--readable-ymin`, `--readable-ymax`; `--auto-readable-ymax` restores the adaptive system-plus-long-reference limit), so high short-rod references do not force the scale; short-rod curves are clipped or omitted in the readable PNG while remaining in CSV and full-range PNGs. It reports nearest-reference families by sorted index only and does not use descendant tracking, FEM, axial references, article styling, crossing claims, strict gap verification, or determinant changes. |
| Beta=90 mu-scan `Lambda(eta)` with single-rod references | `python scripts/analysis/thickness_mismatch/maps/plot_lambda_eta_beta90_mu_scan_with_single_rod_refs.py` | generalized from the beta=90 `Lambda(mu)` reference workflow; helpers: `formulas_thickness_mismatch.py`, `plot_lambda_mu_beta90_eta0p1_with_single_rod_refs.py` | Per-mu folders under `results/lambda_eta_beta90_mu_scan_single_rod_refs/mu_0p0/`, `mu_0p3/`, and `mu_0p6/` containing sorted-system, reference, nearest-match, all-candidate, branch-family summary CSVs, readable/full PNGs, and reports; plus `mu_scan_summary.csv` and `mu_scan_report.md` | Diagnostic-only sorted-frequency `Lambda(eta)` maps for `beta=90 deg`, `epsilon=0.0025`, `mu=0,0.3,0.6`, and `eta=-0.5..0.5`, overlaid with eta-dependent isolated-rod CP/FP and CC/FF bending references using `Lambda_ref=alpha*sqrt(tau_i)/L_i`. The readable PNGs use the fixed range `0 <= Lambda <= 13`; short-rod curves may be clipped or omitted from that readable view while remaining in CSV and full-range PNGs. It reports nearest-reference families by sorted index only and does not use descendant tracking, FEM, axial references, article styling, crossing claims, strict gap verification, or determinant changes. |
| Beta=90 in-plane sorted roots vs isolated-rod references | `python scripts/analysis/thickness_mismatch/audits/compare_beta90_in_plane_first6_to_single_rod_refs.py` | `formulas_thickness_mismatch.py`, `FreqMuNet.py`, `solvers.py` | `results/beta90_in_plane_single_rod_reference_comparison/beta90_in_plane_first6_reference_matches.csv`, `beta90_in_plane_first6_reference_all_candidates.csv`, `beta90_single_rod_reference_frequencies.csv`, `beta90_in_plane_first6_reference_matches_report.md`, and optional diagnostic PNG | Single-case analytic-only proximity audit for `beta=90 deg`, `epsilon=0.0025`, `mu=0.3`, `eta=0.1`. It compares the first sorted in-plane EB roots with isolated-rod CP/FP and CC/FF bending references using `Lambda_ref=alpha*sqrt(tau_i)/L_i`. It does not use descendant tracking, axial references, FEM, article figures, crossing claims, or determinant changes. |
| Single-rod reference-curve correctness audit | `python scripts/analysis/thickness_mismatch/audits/audit_single_rod_reference_curves.py` | `formulas_thickness_mismatch.py`, `FreqMuNet.py`, `solvers.py`; inspects the beta=90 reference plotting/audit scripts and generated reference CSVs | `results/single_rod_reference_audit/eta_sign_length_thickness_audit.csv`, `boundary_condition_alpha_root_residuals.csv`, `single_rod_reference_curve_value_audit.csv`, and `single_rod_reference_curve_audit_report.md` | Diagnostic-only audit confirming that dashed isolated-rod reference families use `Lambda_ref=alpha*sqrt(tau_i)/L_i`, length labels come from `L_i`, eta-sign thickness labels come from `tau_i`, generated reference CSV values reproduce the formula, and CP/CC alpha roots satisfy the intended characteristic equations. It does not generate article figures, run FEM/3D FEM, make crossing claims, or change analytic formulas, determinants, old solvers, or baseline results. |
| Beta=90 first-six in-plane sorted mode shapes | `python scripts/analysis/thickness_mismatch/shapes/plot_in_plane_mode_shapes_beta90_mu0p3_eta0p1_first6.py` | `formulas_thickness_mismatch.py`, `analytic_coupled_rods_shapes.py`, `thickness_mismatch_mac_tracking.py`, `plot_mode_shapes_eta_beta_scan.py` geometry helpers | `results/in_plane_mode_shapes_beta90_mu0p3_eta0p1_first6/mode_shape_sorted1_beta90_mu0p3_eta0p1.png` through `mode_shape_sorted6_beta90_mu0p3_eta0p1.png`, `mode_shapes_sorted1_to_sorted6_beta90_mu0p3_eta0p1_grid.png`, `mode_shapes_beta90_mu0p3_eta0p1_first6_summary.csv`, and `mode_shapes_beta90_mu0p3_eta0p1_first6_report.md` | Diagnostic-only fixed-parameter in-plane mode-shape plots for the first six sorted roots at `beta=90 deg`, `epsilon=0.0025`, `mu=0.3`, and `eta=0.1`. Mode selection is fixed sorted index only, not descendant tracking. The visualization sign convention flips the whole eigenvector when the dominant global vertical lobe on the short rod points downward; fallback rules are recorded in the report. It does not use FEM, article styling, determinant changes, old solver changes, or baseline-result changes. |
| Out-of-plane EB+torsion mode-character beta maps | `python scripts/analysis/thickness_mismatch/maps/plot_out_of_plane_mode_character_beta.py` | none | `results/out_of_plane_mode_character_beta/out_of_plane_mode_character_beta.csv`, per-eta sorted-root mode-character PNGs, per-eta torsion-fraction PNGs, and report | Diagnostic-only sorted-mode character map for the same out-of-plane continuum, using the independent 1D FEM stiffness-energy split to add bending/torsion fractions. It is intended to show where torsion-dominated or mixed modes enter the sorted spectrum; it does not use descendant tracking, 3D FEM, or article styling. |
| Out-of-plane EB+torsion analytic vs 1D FEM validation | `python scripts/analysis/thickness_mismatch/audits/compare_out_of_plane_analytic_vs_1d_fem.py` | none | `results/out_of_plane_fem_validation/out_of_plane_analytic_vs_1d_fem_comparison.csv`, `results/out_of_plane_fem_validation/out_of_plane_analytic_vs_1d_fem_convergence.csv`, `results/out_of_plane_fem_validation/out_of_plane_analytic_vs_1d_fem_report.md` | Diagnostic-only sorted-frequency validation for the out-of-plane determinant against an independent one-dimensional Euler--Bernoulli plus Saint-Venant torsion FEM. It checks determinant signs, joint coupling, and Lambda scaling within the same continuum model; it does not implement or run 3D FEM and does not use descendant tracking. |
| 3D FEM environment smoke check | `python scripts/analysis/thickness_mismatch/audits/check_3d_fem_environment.py --gmsh-exe path/to/gmsh.exe --ccx-exe path/to/ccx.exe` | reuses `scripts/analysis/solid_fem_single_rod_fixed_fixed.py` straight-cylinder helpers | `results/_smoke/3d_fem_environment_check/3d_fem_environment_check_report.md`, `fem_3d_raw_modes.csv`, `fem_3d_environment_check_comparison.csv`, Gmsh/CalculiX logs, and `solid_fem_case/` | Diagnostic-only local solver smoke check for the beta=0, mu=0, eta=0 straight fixed-fixed cylinder at `epsilon=0.02`, `L=2`, `r0=2*epsilon`, and at least 40 requested modes. It resolves Gmsh/CalculiX paths, runs one small mesh/solve when both tools are available, parses the CalculiX RAD/TIME table to `Lambda=sqrt(omega/epsilon)`, and does not run beta, mu, eta, mesh-convergence, or coupled-angle 3D FEM scans. |
| Stepped beta=0 EB/Timoshenko/3D FEM validation | `python scripts/analysis/thickness_mismatch/audits/validate_eb_timo_3d_beta0_stepped.py --run-3d-fem` | reuses `formulas_thickness_mismatch.py`, `scripts/lib/variable_length_timoshenko.py`, and the straight fixed-fixed cylinder Gmsh/CalculiX helpers | `results/eb_vs_timoshenko_3d_validation/stepped_beta0_mu0p5_eta0p1/analytic_*.csv`, `fem_3d_raw_*.csv`, mode-by-mode and pair-average comparison CSVs, per-epsilon PNGs, combined PNG, and report | Diagnostic-only straight coaxial stepped-cylinder comparison for `beta=0`, `mu=0.5`, `eta=0.1`, and `epsilon=0.0025,0.01,0.05`. The 3D geometry is two fused coaxial cylinder segments with clamped ends, not a coupled-angle rod case; analytic formulas, determinants, root solvers, and the Timoshenko shear coefficient are not modified. |
| Uniform beta=0 EB/Timoshenko/3D FEM validation | `python scripts/analysis/thickness_mismatch/audits/validate_eb_timo_3d_beta0_uniform_eps0p05.py --run-3d-fem` | reuses `formulas_thickness_mismatch.py`, `scripts/lib/variable_length_timoshenko.py`, `scripts/analysis/compare_single_rod_eb_timoshenko.py`, and the straight fixed-fixed single-cylinder Gmsh/CalculiX helpers | `results/eb_vs_timoshenko_3d_validation/uniform_beta0_eps0p05/analytic_*.csv`, `fem_3d_raw_*.csv`, first-8 sorted and bending pair-average comparison CSVs, four PNGs, and report | Diagnostic-only straight uniform cylinder comparison for `beta=0`, `mu=0`, `eta=0`, and default `epsilon=0.05`, with optional `--also-run-eps0p1` stress testing. It compares the first eight sorted in-plane EB/Timoshenko roots to the first eight sorted 3D modes and separately compares bending roots to close 3D bending doublet pair averages. It does not run stepped-cylinder, beta-scan, or coupled-angle 3D FEM cases and does not modify formulas or solvers. |
| Experimental full-spectrum analytic-vs-3D-FEM scaffold | `python scripts/analysis/thickness_mismatch/audits/compare_full_spectrum_analytic_vs_3d_fem.py --smoke --skip-3d-fem`; straight-uniform smoke runner: `python scripts/analysis/thickness_mismatch/audits/run_full_spectrum_3d_fem_smoke_case.py --run-3d-fem`; extraction audit: `python scripts/analysis/thickness_mismatch/audits/audit_full_spectrum_3d_fem_smoke_extraction.py`; straight mesh convergence: `python scripts/analysis/thickness_mismatch/audits/run_straight_uniform_3d_mesh_convergence.py` | existing 3D solid FEM workflows remain separate baselines | Scaffold: `results/full_spectrum_analytic_vs_3d_fem/full_spectrum_analytic_union.csv`, placeholder raw/comparison FEM CSVs, nearest-match CSV, and report; smoke runner/audit: `results/full_spectrum_analytic_vs_3d_fem/smoke_straight_uniform/analytic_union.csv`, `fem_3d_raw_modes.csv`, `analytic_vs_3d_fem_matched.csv`, `nearest_match_analytic_vs_3d_fem.csv`, `fem_3d_raw_frequency_audit.csv`, `3d_fem_extraction_audit.md`, and `report.md`; convergence: `results/full_spectrum_analytic_vs_3d_fem/straight_uniform_mesh_convergence/mesh_convergence_summary.csv`, `mesh_convergence_modes.csv`, and `report.md` | Experimental audit scaffold only. The analytic full spectrum is the sorted union of in-plane EB thickness-mismatch roots and out-of-plane EB+torsion roots. The straight-uniform runner and mesh-convergence audit reuse the fixed-fixed straight-cylinder 3D helpers in dedicated result folders for the beta=0, mu=0, eta=0 total-length-2 smoke case. The extraction audit checks raw CalculiX frequency tables, parser inclusion, and nearest-frequency agreement; it does not modify production Gmsh/CalculiX workflows. |
| Coupled-angle full-spectrum 3D FEM smoke | `python scripts/analysis/thickness_mismatch/audits/run_coupled_mu0p3_eta0_beta15_3d_smoke.py` | existing coupled equal-rod solid FEM workflows remain separate baselines | `results/full_spectrum_analytic_vs_3d_fem/coupled_mu0p3_eta0_beta15/analytic_union.csv`, `fem_3d_raw_modes.csv`, `sorted_index_match.csv`, `nearest_frequency_match.csv`, and `report.md` | Diagnostic-only first coupled-angle smoke case for `beta=15 deg`, `mu=0.3`, `eta=0`, `epsilon=0.0025`, using variable-length fused cylinders and existing Gmsh/CalculiX extraction helpers. It is not article validation, does not scan beta/mu/eta, and does not modify analytic determinants or production FEM workflows. |
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
| Descendant `Lambda(mu)` eta sweep | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py --betas 7.5` | `plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py`, `track_lambda_mu_thickness_mismatch_eta_sweep.py` | `results/thickness_mismatch_lambda_mu_beta{beta}_eps0p0025_eta_m0p5_0_p0p5_descendants.*` | Current clean presentation-style plot; no args keeps the historical beta=15 output path, `--betas` writes one PNG/report/tracking-CSV/warning-CSV set per angle, and tracking warnings stay in reports/data, not as x-markers. |
| Eta=0.5 global spectrum overview | `python scripts/analysis/plot_thickness_mismatch_eta_p0p5_global_spectrum.py` | none | `results/thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.*` | Shows sorted roots, descendant branches, canonical sorted positions, and unresolved candidate assignments. |
| Lambda(beta) sorted/descendant identity audit | `python scripts/analysis/audit_lambda_beta_sorted_descendant_thickness_mismatch.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0p5_eps0p0025_mu0_lambda_beta_sorted_descendant_audit.*`, `results/eta0p5_eps0p0025_mu0_sorted5_identity_summary.csv`, plus eta=0 comparison outputs | Analytic-only EB beta sweep at `epsilon=0.0025`, `mu=0`, and `eta=0.5` with an eta=0 comparison. Descendants are seeded at `beta=0`, `mu=0`; sorted position remains metadata. The current audit finds sorted 5 at `beta=15 deg`, `eta=0.5` is descendant 5, while the eta=0 comparison has sorted 5 as descendant 6. |
| Eta scan first-six rearrangement audit | `python scripts/analysis/audit_eta_scan_first6_rearrangement.py` | `audit_lambda_beta_sorted_descendant_thickness_mismatch.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta_scan_eps0p0025_mu0_first6_rearrangement.csv`, `results/eta_scan_eps0p0025_mu0_first6_rearrangement_summary.csv`, `results/eta_scan_eps0p0025_mu0_first6_rearrangement_report.md`, `results/eta_scan_eps0p0025_mu0_first6_*.png` | Analytic-only EB eta audit over `eta=-0.5..0.5`, `beta=0..90 deg`, `epsilon=0.0025`, `mu=0`. Rearrangement means any descendant among 1..6 changes sorted position over beta. The current grid finds no rearrangement for `|eta| >= 0.018`, clean rearrangement near zero for `0.002 <= |eta| <= 0.016`, and the exact `eta=0` equal-thickness case as tracking-unreliable with pairwise sorted-position exchanges. |
| Eta=0 mu scan first-six rearrangement audit | `python scripts/analysis/audit_mu_scan_eta0_first6_rearrangement.py` | `audit_lambda_beta_sorted_descendant_thickness_mismatch.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/mu_scan_eta0_eps0p0025_first6_rearrangement.csv`, `results/mu_scan_eta0_eps0p0025_first6_rearrangement_summary.csv`, `results/mu_scan_eta0_eps0p0025_first6_rearrangement_report.md`, `results/mu_scan_eta0_eps0p0025_*.png` | Analytic-only EB mu audit over `mu=0..0.9`, `eta=0`, `beta=0..90 deg`, and `epsilon=0.0025`. Rearrangement means any descendant among 1..6 changes sorted position over beta after beta=0 seed degeneracy handling. The current grid finds true rearrangement only for `0.001 <= mu <= 0.002`, marks `mu=0` tracking-unreliable, and finds no rearrangement for `0.003 <= mu <= 0.9`. |
| Eta=0 fixed-mu `Lambda(beta)` rearrangement check | `python scripts/analysis/plot_lambda_beta_eta0_eps0p0025_mu_rearrangement_check.py` | `audit_mu_scan_eta0_first6_rearrangement.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/lambda_beta_eta0_eps0p0025_mu0p001_first6.png`, `results/lambda_beta_eta0_eps0p0025_mu0p002_first6.png`, `results/lambda_beta_eta0_eps0p0025_mu0p003_first6.png`, `results/lambda_beta_eta0_eps0p0025_mu_rearrangement_check.csv`, `results/lambda_beta_eta0_eps0p0025_mu_rearrangement_check_report.md` | Analytic-only EB refined beta-grid plot for `mu=0.001,0.002,0.003`, `eta=0`, and `epsilon=0.0025`. It uses the mu-scan shape-MAC descendant tracker, beta step `0.1 deg`, warning markers, and spike checks. Current result: no reproducible first-six rearrangement, no real crossings, and no spike artifacts for these three fixed-mu checks; the earlier `0.25 deg` narrow mu-scan signal is grid/tracking-sensitive diagnostic evidence. |
| Eta=0 positive-mu adjacent-gap verification | `python scripts/analysis/audit_eta0_mu_positive_gap_rearrangement_verification.py` | `audit_mu_scan_eta0_first6_rearrangement.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0_eps0p0025_mu_positive_gap_verification_summary.csv`, `results/eta0_eps0p0025_mu_positive_gap_verification_details.csv`, `results/eta0_eps0p0025_mu_positive_gap_verification_report.md`, `results/eta0_eps0p0025_mu_positive_gap_scaling.png`, `results/eta0_eps0p0025_mu_positive_gap_min_beta.png`, `results/eta0_eps0p0025_mu_positive_gap_lambda_beta_selected.png` | Strict analytic-only EB audit for positive `mu=1e-6..0.9` at `eta=0`, separating adjacent sorted-root gaps from descendant-label tracking. It refines beta-local gap minima for sorted pairs 1-2 through 6-7 and compares default/strict root settings. Current result: all 19 tested positive `mu` rows have resolved positive sorted gaps; the smallest is about `1.314e-4` at `mu=1e-6`, `beta=13.308 deg`, pair 5-6. Small-mu descendant swaps are fragile tracking diagnostics, not accepted true crossings. |
| Eta-parameter positive-gap verification | `python scripts/analysis/audit_eta_parameter_positive_gap_verification.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta_parameter_positive_gap_verification_summary.csv`, `results/eta_parameter_positive_gap_verification_details.csv`, `results/eta_parameter_positive_gap_verification_report.md`, `results/eta_parameter_positive_gap_min_gap_vs_eta.png`, `results/eta_parameter_positive_gap_min_gap_vs_beta_eta_scan.png`, `results/eta_parameter_positive_gap_min_gap_vs_mu_eta_scan.png`, `results/eta_parameter_positive_gap_classification_map.png`, plus selected `Lambda` plots | Strict analytic-only EB audit for eta-, beta-, and mu-parameter positive-gap checks at `epsilon=0.0025`. It solves the first 14 sorted roots, audits adjacent sorted pairs 1-2 through 6-7, tracks first-eight descendants only as label diagnostics, applies local strict root repairs where default sign scans miss close roots, and refines the special `mu` windows near `0.15864`, `0.379`, and `0.716`. Current result: all tested rows are resolved positive gaps with no true sorted-root crossings or unresolved cases; the smallest tested eta-nonzero gap is about `4.0747e-4` for `beta=5 deg`, `eta=0.5`, sorted pair 1-2. |
| Eta-mu-beta sorted-frequency maps | `python scripts/analysis/plot_diagnostic_eta_mu_beta_frequency_maps.py` | `formulas_thickness_mismatch.py` | `results/diagnostic_lambda_eta_beta0_mu0_eps0p0025.*`, `results/diagnostic_lambda_beta_eps0p0025_mu_eta_slices.csv`, `results/diagnostic_lambda_beta_eps0p0025_mu_eta_grid_overview.png`, `results/diagnostic_eta_mu_beta_heatmap_metrics_eps0p0025.*`, and six `results/diagnostic_heatmap_*_eta_mu_eps0p0025.png` heatmaps | Diagnostic-only sorted-frequency mapping for the eta/mu/beta landscape. It plots `Lambda(eta)` at `beta=0, mu=0`, sorted `Lambda(beta)` slices for `mu=0.1,0.2,0.3` and `eta=+-0.1,+-0.2,+-0.3`, and eta-mu heatmaps of minimum adjacent sorted gaps plus beta-sensitivity metrics. Current fallback heatmap grid is eta step `0.025`, mu step `0.05`, beta step `1 deg`, first 10 sorted roots; current run finds min `g_min=0.00838632600513` at `eta=0.5`, `mu=0.35`, `beta=23 deg`, pair 6-7, and the max `S_mean`/`S_max` location at `eta=-0.3`, `mu=0.6`, with no solver warnings or missing roots. These are diagnostic close-approach candidates, not crossing/no-crossing claims. |
| Eta=0.1, mu=0.3 branch-2 beta transition audit | `python scripts/analysis/audit_mode_shape_branch2_beta_transition_eta0p1_mu0p3.py` | `plot_mode_shapes_eta_beta_scan.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/mode_shape_branch2_transition_audit_eta_0p1_mu_0p3_eps_0p0025.csv`, `results/mode_shape_branch2_transition_audit_eta_0p1_mu_0p3_eps_0p0025.md`, `results/mode_shape_branch2_transition_*_eta_0p1_mu_0p3_eps_0p0025.png` | Local diagnostic-only audit for the apparent descendant-2 sorted-index change in the eta=0.1, mu=0.3 mode-shape beta scan. It checks the sorted 2-3 gap over beta=4..7 deg at step 0.05, tracks descendants 2 and 3 only as shape-label diagnostics from beta=0, and reports a positive sorted gap with local tracking ambiguity rather than a frequency crossing. |
| Eta=0.1, mu=0.3 fixed sorted 2/3 beta-scan mode shapes | `python scripts/analysis/plot_mode_shapes_eta_beta_scan_sorted_modes.py` | `plot_mode_shapes_eta_beta_scan.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/mode_shapes_sorted2_sorted3_eta_0p1_mu_0p3_eps_0p0025/` | Diagnostic-only analytic EB mode-shape plots for fixed sorted indices 2 and 3 over beta=0..15 deg. Selection is the sorted root position independently at each beta, not descendant tracking; the default `short-rod-up` sign convention orients plotted eigenvectors, while MAC and sorted 2-3 gap outputs are metadata for interpretation. |
| Eta=+-0.5 sorted `Lambda(mu)` with single-rod references | `python scripts/analysis/plot_lambda_mu_eta_m0p5_with_single_beam_refs.py --eta 0.5 --betas 0 5 10 15` | `formulas_thickness_mismatch.py`, `FreqMuNet.py`, `solvers.py` | `results/diagnostic_lambda_mu_eta_*eps0p0025_beta_*deg_single_beam_refs.*` | Diagnostic-only sorted-root plots for selected beta values at `epsilon=0.0025`. The coupled curves are sorted frequencies, not descendants. CP/FP and CC/FF single-rod reference families are variable `Lambda(mu)` curves computed separately as `alpha*sqrt(tau_i)/L_i`; they are not horizontal reference lines. For `eta=-0.5`, rod 1 is short/thick and rod 2 is long/thin; for `eta=0.5`, rod 1 is short/thin and rod 2 is long/thick, so reference styles follow the actual thin/thick role rather than rod number. Default PNGs zoom to the coupled spectrum and clip reference segments only for plotting; `--full-scale-refs` writes separate full-scale PNGs. |
| Eta-mu-beta global trend post-processing | `python scripts/analysis/summarize_diagnostic_eta_mu_beta_global_trends.py` | consumes the generated diagnostic map CSVs only | `results/diagnostic_eta_mu_beta_global_trends_by_mu_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_trends_by_eta_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_trends_eta_sign_abs_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_pair_transition_by_mu_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_beta_at_gmin_histogram_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_sensitivity_branch_dominance_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_surrogate_fits_eps0p0025.csv`, `results/diagnostic_eta_mu_beta_global_trends_eps0p0025_report.md`, and eight `results/diagnostic_global_*_eps0p0025.png` plots | CSV-only global post-processing of the sorted-frequency map. It computes by-mu/by-eta quantiles, eta-sign asymmetry, pair transitions, beta-at-minimum histograms, branch-sensitivity dominance, and low-order surrogate fits without recomputing roots or doing local refinement. Current result: mu is the stronger global driver of `g_min`, `S_mean`, and `S_mean_rel`; median `g_min` increases with mu, while median beta sensitivity decreases with mu. |
| Selected-eta `Lambda(beta)` descendant plots | `python scripts/analysis/plot_eta_scan_selected_lambda_beta_descendants.py` | `audit_lambda_beta_sorted_descendant_thickness_mismatch.py`, `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0.png`, `results/eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0.csv`, `results/eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0_report.md` | Analytic-only 2x3 plot of smooth tracked descendant `Lambda(beta)` curves for eta values `-0.016,-0.004,-0.002,0.002,0.004,0.016` at `epsilon=0.0025`, `mu=0`. Uses shape-MAC continuation from `beta=0`, accepts near-degenerate best-MAC steps as diagnostic metadata instead of falling back to sorted roots, and runs a spike-artifact check on the plotted curves. |
| Refined branch-identity audit | `python scripts/analysis/check_thickness_mismatch_branch_identity_eta_p0p5.py` | none | `results/thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.*` | Local descendant audit for branches 5--7 at eta=0.5. |
| Eta=0.5 crossing audit | `python scripts/analysis/audit_eta0p5_eps0p0025_lambda_mu_crossings.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py` | `results/eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.*` | Pairwise diagnostic for real descendant-frequency crossings versus finite gaps, sorted-position metadata, and tracking artifacts on the dense `mu=0..0.9` grid. |
| Article beta=5 desc1/desc2 crossing audit | `python scripts/analysis/audit_article_beta5_eta0p5_desc1_desc2_crossing.py` | `formulas_thickness_mismatch.py`, `thickness_mismatch_mac_tracking.py`, article local CSV metadata | `results/article_check_beta5_eta0p5_eps0p0025_mu_desc1_desc2_crossing.*` | Focused analytic-only check for the article `beta=5 deg`, `epsilon=0.0025`, `eta=0.5` `Lambda(mu)` figure. It resolves the apparent desc1/desc2 crossing with a local `mu` refinement and fine root scan, reporting that the refined descendants keep a finite positive gap and no sorted-position swap. |
| Descendant-5 corrected veering shape points | `python scripts/analysis/plot_thickness_mismatch_desc5_veering_mode_shapes.py` | `plot_thickness_mismatch_branch_shapes_vs_eta.py`, existing descendant-6 veering points table | `results/thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc5_near45_rank*.png`, `results/thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc5_near56_rank*.png`, `results/thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc5_veering_points.*` | Analytic-only corrected shape set for the eta=0.5 sorted/descendant identity audit: reuses the six mu/local-gap labels from the descendant-6 table but reconstructs descendant 5, because sorted 5 is descendant 5 for eta=0.5. Descendant-6 outputs are retained separately. |
| Descendant-6 veering shape points | `python scripts/analysis/plot_thickness_mismatch_desc6_veering_mode_shapes.py` | `plot_thickness_mismatch_branch_shapes_vs_eta.py`, existing crossing audit outputs | `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_near45_rank*.png`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_near56_rank*.png`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_veering_points.*` | Analytic-only TASK 2 wrapper: recomputes dense sorted adjacent gaps for eta=0.5, selects three distinct local minima each for sorted pairs 4-5 and 5-6, then plots descendant 6 shapes with ranked CSV/report metadata. |
| Descendant-6 veering localization audit | `python scripts/analysis/audit_desc6_veering_localization.py` | `plot_thickness_mismatch_branch_shapes_vs_eta.py`, `desc6_veering_points.csv` | `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_localization_audit.csv`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/desc6_localization_audit_report.md` | Analytic-only rod-localization audit for the six selected descendant-6 veering shapes; reports displacement, slope, curvature, and tau-weighted EB bending-energy fractions without selecting new mu points. |
| Eta=0.5 analytic/FEM overlay | `python scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py` | `fem_check_thickness_mismatch_eta_p0p5_beta15.py` | `results/thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5.*` | Lightweight plotting comparison. The older FEM check keeps the detailed CSV/MAC diagnostics. |
| Eta=0.5 isolated-rod references, beta=15 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods.*` | Uses the clamped-supported / clamped-pinned reference convention and eta-normalized `Lambda_rod`. |
| Eta=0.5 isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods.*` | Same reference convention as beta=15, different beta only. |
| Eta=0.5 fixed-fixed isolated-rod references, beta=45 | `python scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py` | none | `results/thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed.*` | Uses clamped-clamped / fixed-fixed reference roots. |
| Branch-shape eta overlay | `python scripts/analysis/plot_thickness_mismatch_branch_shapes_vs_eta.py --branch-index 5 --mus 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 --etas -0.5 0 0.5` | `plot_thickness_mismatch_branch5_shapes.py` | `results/thickness_mismatch_desc5_mode_shapes_beta15_eps0p0025/`, `results/thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5/` | Parameterized descendant-shape overlay for several eta values at target mu values, using shape-MAC tracking from `mu=0`; optional compact filenames and aggregate CSV/report outputs support descendant-5 mu sweeps and descendant-6 near-interaction shape diagnostics. |
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
