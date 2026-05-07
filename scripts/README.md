# Scripts guide

This directory keeps the historical script names working, but the preferred entry points are now in `scripts/run/`.
The original root-level scripts remain in place because many scripts import each other through `scripts.*`; moving them now would require broad import churn with no numerical benefit.

## Branch identity and current sorted index

Analytic branch identity is defined at the base point `beta = 0`, `mu = 0` for each `epsilon` independently. A `branch_id` such as `bending_desc_05` means "the branch seeded by base sorted index 5"; it does not mean "whatever root is currently fifth."

The current place of that branch in the sorted spectrum is `current_sorted_index`. It may change along `beta` or `mu`, so plots and CSV tables should label it separately from `branch_id`. Do not label a continued branch only by root number.

All analytic branch-selecting frequency and shape plots should use `scripts/lib/analytic_branch_tracking.py`. That helper performs in-memory shape-MAC tracking with a small frequency tie-breaker, adaptive step refinement for low-MAC transitions, and local smallest-singular-value candidate recovery when a sign-change root scan misses a branch near a close interaction. CSV/JSON tracking paths are optional debug artifacts only, disposable, and not a source of truth.

Branch tracking results used in figures must be regression-tested when they resolve a previous inconsistency. Low-MAC assignments are not canonical unless explicitly accepted in a diagnostic run, for example with `--allow-low-mac`.

## Main commands

### Beta sweep at `mu = 0`

- Task: compare analytic and FEM branches over `beta` at `mu = 0` for four radii.
- Command: `python scripts/run/run_beta_sweep_mu0_four_radii.py`
- Main parameters: fixed by the underlying presentation script; radii are `0.005, 0.01, 0.015, 0.02`.
- Results: `results/beta_sweep_mu0_four_radii_compare.png`, `results/beta_sweep_mu0_four_radii_compare.csv`, and per-radius CSV files when recomputed.
- Use when: you need the four-radius `beta` sweep figure for the low spectrum.
- Do not use when: you need a `mu` sweep, a full branchwise FEM audit, or custom branch descendant shapes.

### Mu sweep at `beta = 0`

- Task: compare analytic and FEM bending branches over `mu` at `beta = 0` for four radii.
- Command: `python scripts/run/run_mu_sweep_beta0_four_radii.py`
- Main parameters: fixed by the underlying presentation script; radii are `0.005, 0.01, 0.015, 0.02`.
- Results: `results/mu_sweep_beta0_four_radii_compare.png`, `results/mu_sweep_beta0_four_radii_compare.csv`, and per-radius CSV files.
- Use when: you need the baseline `beta = 0` `Lambda(mu)` comparison.
- Do not use when: the coupling angle is nonzero or you need an analytic-only comparison.

### Mu sweep at fixed nonzero beta

- Task: compare analytic and FEM bending branches over `mu` for fixed `beta = 7.5 deg` and `beta = 15 deg`.
- Command: `python scripts/run/run_mu_sweep_fixed_beta_four_radii.py`
- Main parameters: fixed by the underlying presentation script; beta values are `7.5` and `15.0` degrees.
- Results: `results/mu_sweep_beta7p5_four_radii_compare.png`, `results/mu_sweep_beta15_four_radii_compare.png`, and matching combined/per-radius CSV files.
- Use when: you need presentation figures for positive-beta `mu` sweeps across four radii.
- Do not use when: you need `beta = 0` only or a full 2D branchwise audit.

### Analytic-only mu sweep for selected betas

- Task: plot analytic-only `mu` sweeps at fixed radius for selected beta values.
- Command: `python scripts/run/run_mu_sweep_four_betas_analytic.py`
- Main parameters: `--betas`; default is `15 30 45 60`.
- Results: `results/mu_sweep_r015_selected_betas_analytic.png`, `results/mu_sweep_r015_selected_betas_analytic.csv`, and per-beta analytic CSV files.
- Use when: you need an analytic-only comparison of several beta values at the fixed presentation radius.
- Do not use when: FEM markers are required.

Branch-selecting Lambda(mu) plots should use the shared analytic branch tracker. The standalone `src/my_project/analytic/FreqMuNet.py` path now tracks coupled branches by the same `branch_id` / `current_sorted_index` rule as analytic shape plots.

### Fixed-beta analytic Lambda(mu) plot

- Task: plot the coupled-rods analytic `Lambda(mu)` curves at one fixed `beta` for an arbitrary number of canonical tracked branches.
- Command: `python scripts/run/run_lambda_mu_fixed_beta_analytic.py --beta 15 --epsilon 0.0025 --num-modes 7 --num-dashed-lines 7`
- Main parameters: `--beta`, `--epsilon`, `--num-modes`, `--num-dashed-lines`, `--mu-min`, `--mu-max`, `--mu-step`, `--y-max`, `--output`, `--csv-output`, `--show`, `--allow-low-mac`.
- Results: by default deterministic PNG/CSV paths such as `results/lambda_mu_beta15_eps0p0025_modes7_dash7.png` and `results/lambda_mu_beta15_eps0p0025_modes7_dash7.csv`.
- Use when: you need the same visual style as `FreqMuNet.py` but with a convenient script wrapper and a summary CSV for any chosen number of coupled modes.
- Do not use when: you need FEM comparison markers or one highlighted branch.

The solid lines are the first `--num-modes` canonical analytic branches seeded at `beta = 0`, `mu = 0`; dashed CS reference families use `--num-dashed-lines` roots per family. The CSV records the displayed line rank, canonical `branch_id`, `base_sorted_index`, `current_sorted_index`, and `Lambda` at each plotted `mu`. No full tracking debug CSV is saved by default. If the shared tracker detects a low-MAC assignment on a displayed branch, ordinary plotting fails fast; `--allow-low-mac` is only for exploratory diagnostics.

### Analytic coupled-rods mode shape

- Task: build one analytic mode shape directly from the determinant nullspace, without FEM comparison.
- Command: edit the `USER PARAMETERS` block in `scripts/run/run_analytic_coupled_rods_mode_shape_ru.py`, then run `python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py`.
- Main parameters: `--mode-number`, `--beta`, `--mu`, `--epsilon`, `--l-total`, `--plot-kind`, `--mode-scale`, `--normalize`, `--output`, `--dpi`, `--figsize`, `--show`, `--save-samples-csv`, `--save-diagnostics-csv`, `--print-diagnostics`, `--no-print-diagnostics`.
- Results: by default `results/analytic_mode_shape_{plot_kind}_mode..._beta..._mu..._eps..._scale..._ru.png`; optional samples and diagnostics CSV paths are supported.
- Use when: you need the analytic reconstructed local components or analytic deformed geometry for a selected determinant root.
- Do not use when: you need FEM descendant tracking or an analytic-vs-FEM component comparison.

`--mode-number` is the analytic root index in ascending `Lambda` order; it is not a FEM descendant label such as `bending_desc_05`. The right-arm output convention is joint-to-external, so the right component curves use `s/L = 0` at the joint and `s/L = 1` at the right external clamp. The diagnostics report the determinant null-vector residual, endpoint clamp residuals, and analytic arm-wise axial/bending energy shares.

Examples:

```bash
python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py
python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py --mode-number 5 --beta 15 --mu 0 --epsilon 0.0025 --plot-kind components
python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py --mode-number 5 --beta 30 --mu 0 --epsilon 0.0025 --plot-kind full --mode-scale 0.08
python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py --mode-number 5 --beta 30 --mu 0 --epsilon 0.0025 --plot-kind transverse
```

### Tracked bending descendant shapes

- Task: plot Russian-labeled mode shapes for a tracked bending descendant branch.
- Command: `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_01`
- Main parameters: `--branch-id`, `--beta`, `--target-mus`, `--radii`, `--output`, `--title-label`.
- Results: by default `results/tracked_bending_descendant_shapes_beta..._ru.png`; custom `--output` is supported.
- Use when: you need mode-shape panels for a tracked bending descendant such as `bending_desc_01`, `bending_desc_02`, or `bending_desc_04`.
- Do not use when: you need branchwise tables, energy summaries, analytic-vs-FEM frequency sweeps, or a single explicitly parameterized shape.

### Single tracked descendant shape

- Task: plot one Russian-labeled mode shape for one tracked bending descendant at explicit `branch`, `mu`, `epsilon`, and `beta`.
- Command: edit the `USER PARAMETERS` block in `scripts/run/run_tracked_bending_descendant_shape_ru.py`, then run `python scripts/run/run_tracked_bending_descendant_shape_ru.py`.
- Main parameters: `--branch-number`, `--branch-id`, `--mu`, `--epsilon`, `--beta`, `--l-total`, `--plot-kind`, `--mode-scale`, `--normalize`, `--output`, `--title-label`, `--dpi`, `--figsize`, `--show`, `--print-diagnostics`, `--save-diagnostics-csv`, `--diagnostics-level`.
- Results: by default `results/tracked_bending_descendant_shape_{plot_kind}_beta..._ru.png`; custom `--output` is supported.
- Use when: you need one final tracked descendant shape for a specific parameter point.
- Do not use when: you need comparison panels over several `mu` values or several radii.

Ordinary use:

1. Edit `USER PARAMETERS` near the top of `scripts/run/run_tracked_bending_descendant_shape_ru.py`.
2. Run:

```bash
python scripts/run/run_tracked_bending_descendant_shape_ru.py
```

CLI arguments remain available as optional overrides. For example:

`--branch-number` is the descendant number used to form `bending_desc_XX`; it is not the current sorted mode index printed as "место в спектре".

`--plot-kind full` shows the total modal displacement, analogous to a scaled total-deformation mode-shape view. `--plot-kind transverse` removes local axial displacement from the geometry and shows the local bending-deflection shape on each arm. `--plot-kind components` draws diagnostic local axial/transverse component curves for the two arms, normalized together so relative axial/transverse size remains visible. `--normalize auto` resolves to `max-full` for `full` and `components`, and to `max-transverse` for `transverse`.

`--mu` must be non-negative because the tracked continuation starts at `mu = 0`; values above `0.9` run with a warning because they are outside the usual analysis window.

`--diagnostics-level all` adds three checks for local components. The projection check reconstructs global components from local right-arm components. The derivative diagnostics estimate whether large local axial displacement also produces large axial strain-like slopes; these are nodal amplitude slopes, not exact element strain fields. The arm-energy diagnostics compares axial and bending participation by arm for the selected eigenvector amplitude using the current FEM local element stiffness matrices and the same local rotation convention as assembly.

Examples:

```bash
python scripts/run/run_tracked_bending_descendant_shape_ru.py --branch-number 5 --mu 0 --epsilon 0.0025 --beta 30 --plot-kind full --mode-scale 0.05
python scripts/run/run_tracked_bending_descendant_shape_ru.py --branch-number 5 --mu 0 --epsilon 0.0025 --beta 30 --plot-kind transverse
python scripts/run/run_tracked_bending_descendant_shape_ru.py --branch-number 5 --mu 0 --epsilon 0.0025 --beta 30 --plot-kind components --save-diagnostics-csv results/desc05_beta30_components_diag.csv
python scripts/run/run_tracked_bending_descendant_shape_ru.py --branch-number 5 --mu 0 --epsilon 0.0025 --beta 30 --plot-kind components --diagnostics-level all --save-diagnostics-csv results/desc05_beta30_all_diag.csv
```

### Analytic vs FEM mode-shape comparison

- Task: compare determinant-nullspace analytic local components with FEM local components for one tracked bending descendant branch.
- Command: `python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py --branch-number 5 --beta 30 --mu 0 --epsilon 0.0025`
- Main parameters: `--branch-number`, `--branch-id`, `--beta`, `--mu`, `--epsilon`, `--l-total`, `--num-roots`, `--root-window`, `--output-prefix`, `--right-coordinate`, `--use-best-orientation`, `--check-endpoint-consistency`, `--compare-energies`, `--check-fem-direct-energy`.
- Results: by default `results/analytic_fem_shape_compare_beta..._{branch_id}_mu..._eps..._components_overlay.png`, matching sampled component CSV, one-row diagnostics CSV, and `*_orientation_scan.csv`; with `--check-endpoint-consistency`, also writes `*_endpoint_consistency.csv` and `*_endpoint_consistency_summary.txt`; with `--compare-energies`, also writes `*_energy_comparison.csv` and appends the key energy fields to the diagnostics CSV; with `--check-fem-direct-energy`, also writes `*_fem_direct_energy_check.csv`.
- Use when: you need an independent analytic/FEM check of local axial and transverse mode-shape components for one tracked descendant.
- Do not use when: you need to change determinant entries, FEM baseline behavior, or branch tracking.

The default comparison convention is `s/L = 0` at each external clamp and `s/L = 1` at the joint. The FEM right-arm local arrays are reversed by default because the current FEM node order on that arm runs from the joint to the right clamp. CSV and PNG component values are max-normalized and analytic signs are aligned to the FEM vector. The orientation scan tries right-arm coordinate order and analytic component sign variants, reports the top five variants to stdout, and writes the full table to CSV. The main overlay keeps the current convention unless `--use-best-orientation` is supplied.

Endpoint consistency diagnostics:

```bash
python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py \
  --branch-number 5 \
  --beta 30 \
  --mu 0 \
  --epsilon 0.0025 \
  --check-endpoint-consistency
```

This checks analytic field reconstruction against the determinant matrix rows, checks basis-column consistency, reports external-clamp residuals, and compares FEM joint local kinematics under the available right-arm coordinate conventions. It is diagnostic postprocessing only and does not modify the model.

Energy comparison diagnostics:

```bash
python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py --branch-number 5 --beta 30 --mu 0 --epsilon 0.0025 --compare-energies
```

This computes arm-wise axial/bending diagnostic energies for the normalized reconstructed analytic shape and for the same FEM mode. FEM energies are evaluated from the current `fem.elem_K(Le)` local axial and bending stiffness blocks after rotating each element vector to local coordinates. The values are reported as diagnostic fractions and shares; the arbitrary modal scale is also recorded in the CSV.

Direct FEM element-stiffness energy check:

```bash
python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py --branch-number 5 --beta 30 --mu 0 --epsilon 0.0025 --check-fem-direct-energy
```

This independently recomputes FEM arm-wise axial/bending energies from the current `fem.elem_K(Le)` local stiffness blocks and compares them against the default `arm_energy_diagnostics(...)` path. It should now report agreement; it is diagnostic postprocessing only and does not modify the FEM baseline.

### Analytic shape in FEM residual diagnostic

- Task: embed a canonically tracked analytic determinant-nullspace shape directly into the baseline FEM DOF vector and compute `Kq - omega^2 Mq`.
- Command: `python scripts/analysis/check_analytic_shape_in_fem_residual.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025`
- Main parameters: edit the `USER PARAMETERS` block or use `--branch-id`, `--branch-number`, `--beta`, `--mu`, `--epsilon`, `--l-total`, `--num-analytic-roots`, `--output-prefix`, `--scan-theta-conventions`, `--localize-residual`, `--scan-bending-basis-pairings`, `--scan-axial-conventions`, `--scan-axial-scales`, `--scan-right-transform-conventions`, `--audit-frequency-scaling`, `--save-debug-csv`, `--no-save-debug-csv`.
- Results: by default `results/analytic_shape_in_fem_residual_beta..._{branch_id}_mu..._eps....csv` plus optional `*_nodes.csv`; with `--scan-theta-conventions`, also writes `*_theta_convention_scan.csv`; with `--localize-residual`, also writes `*_modal_projection.csv`, `*_residual_groups.csv`, and `*_element_residual_energy.csv`; with `--scan-bending-basis-pairings`, also writes `*_bending_basis_pairing_scan.csv` and `*_full_basis_column_audit.csv`; with `--scan-axial-conventions`, also writes `*_axial_convention_scan.csv`; with `--scan-axial-scales`, also writes `*_axial_scale_scan.csv`; with `--scan-right-transform-conventions`, also writes `*_right_transform_convention_scan.csv` and `*_right_transform_sanity.csv`; with `--audit-frequency-scaling`, also writes `*_frequency_scaling_audit.csv`.
- Use when: local component overlays are ambiguous and you need a direct FEM residual plus global DOF/MAC checks for the analytically reconstructed shape.
- Do not use when: you need to select analytic branch identity by FEM frequency; the script uses canonical analytic branch tracking for the analytic branch and FEM tracking only as a diagnostic comparison target.

The diagnostic reports residuals for both plausible right-arm rotation signs, records the selected FEM eigenvector self-residual as an assembly check, compares global/translational/rotational/arm/joint DOF vectors after optimal scalar alignment, and checks that analytic `q -> local components` reproduces the directly reconstructed analytic components. The optional theta-convention scan compares FEM theta DOFs against finite-difference local transverse derivatives and scans analytic embeddings using `dw/ds`, `dw/dxi`, and `dw/dz` sign/scale variants for each arm. The optional residual localization layer records the analytic Rayleigh quotient, FEM modal projection shares, residual DOF-group shares, and element-wise local residual/energy diagnostics. The optional bending-basis pairing scan audits the full `cos`, `sin`, `cosh`, `sinh` representation directly, without Krylov functions, and checks identity/swap/sign pairings plus right-arm `z` sign variants against determinant columns and FEM residuals. The optional axial-convention scan keeps bending fields fixed while changing only local axial components before FEM embedding, covering right-arm sign/reversal variants and left-arm sign flips, and reports Rayleigh, residual, axial-strain, and direct element axial-energy diagnostics. The optional axial-scale scan keeps bending fields fixed while scaling local axial components by physical candidate factors, L2-optimal local axial factors, and a coarse two-arm log grid, then reports residual, Rayleigh, MAC, residual-share, and direct element axial-energy metrics. The optional right-transform scan keeps analytic local components fixed while testing right-arm `T`, `T.T`, `R(+beta)`, `R(-beta)`, inverse-`T`, and FEM-energy-consistent local/global embedding variants; it also writes a unit-vector sanity table and explicitly states the `K_global = T.T K_local T` implication. The optional frequency-scaling audit checks whether `Lambda`, `Lambda^2`, `Lambda^4`, FEM `omega`, FEM `omega^2`, or the analytic Rayleigh factor is being used as the spectral multiplier in `Kq - factor Mq`, and verifies the selected FEM vector self-residual under the same choices. It is diagnostic postprocessing only and does not modify the determinant, formulas, solvers, FEM baseline, physical model, or eigenfrequencies.

### Joint dynamic stiffness analytic/FEM comparison

- Task: condense each clamped arm to the coupled joint and compare FEM Schur-complement dynamic stiffness against an independent full-basis analytic arm response.
- Command: `python scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025`
- Row audit command: `python scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025 --row-audit`
- Constrained row audit command: `python scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025 --constrained-row-audit`
- Force-row scaling audit command: `python scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025 --audit-force-row-scaling`
- Results: writes `results/joint_dynamic_stiffness_compare_...csv`, `*_summary.csv`, and `*_determinant_crosscheck.csv`; with `--row-audit`, also writes `*_determinant_row_audit.csv`, `*_determinant_row_scaling_audit.csv`, `*_determinant_force_row_contribution_audit.csv`, and `*_determinant_null_vector_row_audit.csv`; with `--constrained-row-audit`, also writes `results/constrained_determinant_force_row_audit_*.csv`, `*_summary.csv`, and `results/constrained_null_vector_force_check_*.csv`; with `--audit-force-row-scaling`, also writes `results/force_row_scaling_audit_*.csv`, `*_summary.csv`, and `*_diagnostic_matrix.csv`.
- Use when: direct analytic-in-FEM residuals are large and you need to determine whether the mismatch is in the arm dynamic stiffness operators or in the coupled determinant/joint force rows.

The FEM side assembles per-arm `K - omega^2 M` with the same `elem_K`, `elem_M`, and right-arm transform convention as `python_fem.py`, then eliminates internal DOFs and fixed external clamp DOFs by Schur complement. The analytic side uses full `cos`, `sin`, `cosh`, `sinh` bending fields and `sin`/`cos` axial fields to map joint `[ux, uy, theta]` to nodal residual `[Fx, Fy, M]`. The determinant cross-check records the scaled force rows assumed to correspond to `[M/Lambda^2, epsilon*Fy/Lambda^2, epsilon*Fx/Lambda^2]`. It is diagnostic postprocessing only and does not modify the determinant, formulas, solvers, FEM baseline, physical model, or eigenfrequencies.

The optional row audit compares the current determinant matrix against an independent full-basis endpoint-row matrix and several dynamic-stiffness-induced force-order/sign variants. It prints and records the assumed joint displacement order, generalized force order, determinant force-row mapping, per-row best scalar fits, force-row contribution terms, and determinant-null-vector residuals under each interpretation. This is an audit of conventions and row scaling only; it does not patch `formulas.py` or assert a corrected determinant.

The constrained row audit first restricts coefficient vectors to the SVD nullspace of the first three kinematic determinant rows, then compares only the force-row operators on that kinematically compatible subspace. It records the kinematic singular values, nullspace dimension, left/right joint DOF mismatch on the subspace, force-order/sign variants, row-scaling fits, and the projected determinant-null-vector force check. This avoids interpreting arbitrary coefficient basis columns as valid joint states.

The force-row scaling audit decomposes the determinant force rows into local moment, shear, and axial-force contributions and compares each contribution with the FEM-compatible nondimensional endpoint quantities `M`, `Q`, and `N` for `EI=1` and `EA=1/epsilon^2`. It records the inferred factor such as `1/Lambda^2` or `epsilon/Lambda^2`, including signs and arm-length candidates, and builds a diagnostic-only common physical-force row matrix at the same `Lambda`. This is a scaling/convention audit only and does not change the determinant or `formulas.py`.

### Single-beam exact-shape FEM residual sanity check

- Task: verify that the FEM residual normalization and exact-shape embedding pass simple single-member tests before interpreting coupled-rods analytic-in-FEM residuals.
- Command: `python scripts/analysis/check_single_beam_exact_shape_residual.py`
- Results: writes `results/single_beam_exact_shape_residual_bending.csv` and `results/single_beam_exact_shape_residual_axial.csv` for `N_ELEM = 20, 40, 80, 160`.
- Use when: coupled analytic-in-FEM residuals are unexpectedly large and you need to separate residual-diagnostic errors from coupled-operator or embedding mismatches.

The script uses the current `python_fem.py` `elem_K` and `elem_M` only as diagnostic matrix sources. It checks a clamped-free Euler-Bernoulli bending mode with the exact cantilever shape and a clamped-free axial bar mode sampled on the FEM nodes, then reports the same relative residual normalization, Rayleigh quotient error, and first-FEM-mode comparison in the CSV notes. It is diagnostic postprocessing only and does not modify the determinant, formulas, solvers, FEM baseline, physical model, or eigenfrequencies.

### Article shape validation at `beta = 15 deg`

- Task: run the analytic-vs-FEM component and energy diagnostics for the article small-angle cases.
- Command: `python scripts/analysis/validate_article_shape_cases_beta15.py`
- Main parameters: `--branch-ids`, `--mus`, `--beta`, `--epsilon`, `--right-coordinate`, `--use-best-orientation`, `--output`, `--output-prefix-root`.
- Results: by default `results/article_shape_validation_beta15_summary.csv`, `results/article_shape_validation_beta15_report.md`, plus one overlay PNG, samples CSV, diagnostics CSV, orientation scan CSV, energy comparison CSV, and direct-energy check CSV per branch/mu case.
- Use when: you need a compact review table for `bending_desc_01` and `bending_desc_02` at `beta = 15 deg`, `epsilon = 0.0025`, and the article `mu = 0, 0.1, 0.2` set.
- Do not use when: you need to change determinant entries, FEM baseline behavior, branch tracking, or root matching.

The summary table records the FEM tracked descendant, nearest analytic root index, `Lambda` mismatch, full local component errors, component-wise errors, transverse-only and axial-only shape metrics, analytic/FEM direct axial-energy fractions, right axial/bending shares, and review flags. The FEM energy fields are explicitly sourced from direct local element-stiffness diagnostics. The `article_use_flag` evaluates transverse plotting separately from axial-component mismatch and is a review aid only.

Examples:

```bash
python scripts/analysis/validate_article_shape_cases_beta15.py
python scripts/analysis/validate_article_shape_cases_beta15.py --branch-ids bending_desc_01 bending_desc_02 --mus 0 --beta 15 --epsilon 0.0025
```

### Desc05 analytic full-shape epsilon sweep

- Task: plot analytic full mode shapes for the fifth analytic bending descendant at `beta = 15 deg` over several `epsilon` values.
- Command: `python scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py`
- Main parameters: edit the script `USER PARAMETERS` block or use `--branch-number`, `--beta`, `--mus`, `--epsilons`, `--mode-scale`, `--output-dir`, `--beta-steps`, `--mu-steps`, `--max-refinement-depth`, `--min-beta-step`, `--min-mu-step`, `--shape-metric`, `--save-tracking-debug`, `--allow-low-mac`, `--check-known-contradiction`.
- Results: by default one PNG per configured `mu` value and `results/analytic_full_shapes_desc05_beta15_eps_sweep_summary.csv`.
- Use when: you need discussion figures for the analytic branch seeded as root/branch number 5 at `beta = 0`, `mu = 0`.
- Do not use when: you need FEM descendant tracking or FEM-based root selection.

For each `epsilon`, the script starts independently at `beta = 0`, `mu = 0`, selects the requested analytic base root number, then calls `scripts/lib/analytic_branch_tracking.py` to track by global assignment over sampled-shape MAC. No FEM Lambda is used for selecting analytic roots. The summary CSV records `branch_id`, `base_sorted_index`, `current_sorted_index`, Lambda, MAC diagnostics, analytic energy fractions, and review flags. Full tracking-path CSVs are not written by default; pass `--save-tracking-debug` to write disposable debug files under `results/debug/`. Low-MAC assignments fail fast unless `--allow-low-mac` is passed for an explicit diagnostic run; at larger `beta`, adaptive refinement runs first, and users can tune `--max-refinement-depth`, `--min-beta-step`, and `--min-mu-step`. If canonical tracking still fails, the helper writes a one-case `results/debug/analytic_tracking_failure_*.csv` when possible. The known-contradiction diagnostic is opt-in via `--check-known-contradiction` and is only relevant to the historical `beta=15`, `epsilon=0.0025`, `mu=0.8` case.

### Branchwise FEM audit

- Task: track FEM branch identities over the configured 2D `(beta, mu)` analysis grid.
- Command: `python scripts/run/run_branchwise_fem_audit.py`
- Main parameters: fixed by `scripts/analyze_branchwise_fem_spectrum.py` and `scripts/sweep_grid_policy.py`.
- Results: branchwise CSV tables and diagnostic plots under `results/branchwise_*`.
- Use when: you need a full FEM-only audit of branch identity, gaps, extrema, and axial fraction.
- Do not use when: you only need a presentation figure; this is heavier than the plotting wrappers.

## Analysis / audit scripts

| file | category | purpose | recommended command | keep/move/wrapper/archive |
| --- | --- | --- | --- | --- |
| `scripts/analyze_branchwise_fem_spectrum.py` | analysis/audit | Full FEM branchwise audit over beta/mu grid. | `python scripts/run/run_branchwise_fem_audit.py` | keep root; wrapper in `scripts/run/` |
| `scripts/analyze_flat_mu_bending_energy.py` | analysis/audit | Select flattening branches, compute bending-energy fractions, and plot follow-up shapes. | `python scripts/analyze_flat_mu_bending_energy.py` | keep root; imported by tracked-shape scripts |
| `scripts/analyze_flat_mu_branches.py` | analysis/audit | Build flat-mu candidate metrics, slot occupancy, and representative cases. | `python scripts/analyze_flat_mu_branches.py` | keep root; analysis dependency |
| `scripts/analyze_target_descendants_beta15_r5.py` | analysis/audit | Targeted beta=15, r=5 mm descendant analysis for selected bending branches. | `python scripts/analyze_target_descendants_beta15_r5.py` | keep root; specialized audit |
| `scripts/analysis/check_analytic_shape_in_fem_residual.py` | analysis/diagnostic | Embed a canonically tracked analytic shape into the FEM DOF vector and compute direct FEM residual/global DOF diagnostics. | `python scripts/analysis/check_analytic_shape_in_fem_residual.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025` | keep in `scripts/analysis/` |
| `scripts/analysis/check_single_beam_exact_shape_residual.py` | analysis/diagnostic | Sanity-check FEM residual/Rayleigh convergence on exact single-beam bending and axial-bar modes. | `python scripts/analysis/check_single_beam_exact_shape_residual.py` | keep in `scripts/analysis/` |
| `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py` | analysis/diagnostic | Compare per-arm FEM Schur-complement joint dynamic stiffness with full-basis analytic arm response, with optional determinant force-row audit. | `python scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py --branch-id bending_desc_05 --beta 15 --mu 0.2 --epsilon 0.0025` | keep in `scripts/analysis/` |
| `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py` | analysis/diagnostic | Compare analytic determinant-nullspace local components with one tracked FEM descendant shape. | `python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py --branch-number 5 --beta 30 --mu 0 --epsilon 0.0025` | keep in `scripts/analysis/` |
| `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py` | analysis/diagnostic | Plot analytic-only full shapes for the fifth base bending descendant at beta=15 over an epsilon sweep. | `python scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py` | keep in `scripts/analysis/` |
| `scripts/analysis/validate_article_shape_cases_beta15.py` | analysis/diagnostic | Batch analytic/FEM component and energy validation for article `beta=15` descendant cases. | `python scripts/analysis/validate_article_shape_cases_beta15.py` | keep in `scripts/analysis/` |
| `scripts/compare_beta0_analytic_vs_fem.py` | analysis/audit | Beta=0 analytic/FEM type-aware audit and shared helper source. | `python scripts/compare_beta0_analytic_vs_fem.py` | keep root; many imports depend on it |
| `scripts/compare_beta_positive_type_aware.py` | analysis/audit | Positive-beta type-aware analytic/FEM comparison. | `python scripts/compare_beta_positive_type_aware.py` | keep root; specialized audit |

## Presentation / figure scripts

| file | category | purpose | recommended command | keep/move/wrapper/archive |
| --- | --- | --- | --- | --- |
| `scripts/plot_beta_sweep_mu0_compare.py` | presentation/figure | Per-radius beta sweep support script. | `python scripts/plot_beta_sweep_mu0_compare.py` | keep root; imported by four-radii script |
| `scripts/plot_beta_sweep_mu0_four_radii_compare.py` | presentation/figure | Four-radius beta sweep at `mu = 0`. | `python scripts/run/run_beta_sweep_mu0_four_radii.py` | keep root; wrapper in `scripts/run/` |
| `scripts/plot_mu_sweep_beta0_four_radii_compare.py` | presentation/figure | Four-radius `mu` sweep at `beta = 0`. | `python scripts/run/run_mu_sweep_beta0_four_radii.py` | keep root; wrapper in `scripts/run/` |
| `scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py` | presentation/figure | Four-radius `mu` sweep for fixed positive beta values. | `python scripts/run/run_mu_sweep_fixed_beta_four_radii.py` | keep root; wrapper in `scripts/run/` |
| `scripts/plot_mu_sweep_radius_fixed_four_betas_analytic.py` | presentation/figure | Analytic-only `mu` sweeps for selected betas at fixed radius. | `python scripts/run/run_mu_sweep_four_betas_analytic.py --betas 15 30 45 60` | keep root; wrapper in `scripts/run/` |
| `scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` | main user-facing command | Build one analytic coupled-rods mode shape directly from the determinant nullspace. | `python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` | keep |
| `scripts/plot_tracked_bending_descendant_shapes_ru.py` | main user-facing command | Parameterized Russian mode-shape plot for tracked bending descendants. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_01` | keep root; wrapper in `scripts/run/` |
| `scripts/plot_freq_mu_vs_fem.py` | one-off/obsolete | Older beta=15 frequency comparison against `results/fem_spectrum.csv`. | `python scripts/plot_freq_mu_vs_fem.py` | archive candidate after manual review |

## Internal helpers

| file | category | purpose | recommended command | keep/move/wrapper/archive |
| --- | --- | --- | --- | --- |
| `scripts/lib/analytic_coupled_rods_shapes.py` | internal helper | Shared analytic null-vector reconstruction, endpoint diagnostics, normalization, and analytic arm-energy utilities. | none | keep in `scripts/lib/` |
| `scripts/lib/analytic_branch_tracking.py` | internal helper | Source-of-truth analytic branch identity tracking from `beta=0`, `mu=0` using shape-MAC assignment and `current_sorted_index` diagnostics. | none | keep in `scripts/lib/` |
| `scripts/lib/tracked_bending_descendant_shapes.py` | internal helper | Shared tracked-state extraction and one-case drawing for tracked bending descendant shape plots. | none | keep in `scripts/lib/` |
| `scripts/sweep_grid_policy.py` | internal helper | Shared presentation/analysis grid policy. | none | keep root; moving would require broad import updates |

## Legacy wrappers

These preserve old command paths and old output filenames.

| file | category | purpose | recommended command | keep/move/wrapper/archive |
| --- | --- | --- | --- | --- |
| `scripts/plot_flat_mu_bending_desc_01_mu0_0p1_0p2_ru.py` | legacy compatibility wrapper | Old `bending_desc_01`, `mu=0,0.1,0.2` shape command. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_01` | keep root wrapper |
| `scripts/plot_flat_mu_bending_desc_02_mu0_0p1_0p2_ru.py` | legacy compatibility wrapper | Old `bending_desc_02`, `mu=0,0.1,0.2` shape command. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_02` | keep root wrapper |
| `scripts/plot_flat_mu_bending_desc_04_mu0_0p1_0p2_ru.py` | legacy compatibility wrapper | Old `bending_desc_04`, `mu=0,0.1,0.2` shape command. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_04` | keep root wrapper |
| `scripts/plot_flat_mu_bending_desc_06_clean_ru.py` | legacy compatibility wrapper | Old `bending_desc_06`, `mu=0,0.5,0.9` clean shape command. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_06 --target-mus 0 0.5 0.9` | keep root wrapper |
| `scripts/plot_flat_mu_bending_desc_06_mu0_0p1_0p2_ru.py` | legacy compatibility wrapper | Old `bending_desc_06`, `mu=0,0.1,0.2` shape command. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_06 --target-mus 0 0.1 0.2` | keep root wrapper |

## One-off or obsolete scripts

| file | category | purpose | recommended command | keep/move/wrapper/archive |
| --- | --- | --- | --- | --- |
| `scripts/run_fem_case.py` | one-off/obsolete | Placeholder FEM entry point. | none | archive candidate after manual review |
| `scripts/run_program1.py` | one-off/obsolete | Placeholder entry point. | none | archive candidate after manual review |
| `scripts/run_program2.py` | one-off/obsolete | Placeholder entry point. | none | archive candidate after manual review |

## User-facing wrappers

| file | category | purpose | recommended command | keep/move/wrapper/archive |
| --- | --- | --- | --- | --- |
| `scripts/run/run_beta_sweep_mu0_four_radii.py` | main user-facing command | Friendly wrapper for four-radius beta sweep at `mu = 0`. | `python scripts/run/run_beta_sweep_mu0_four_radii.py` | keep |
| `scripts/run/run_mu_sweep_beta0_four_radii.py` | main user-facing command | Friendly wrapper for four-radius `mu` sweep at `beta = 0`. | `python scripts/run/run_mu_sweep_beta0_four_radii.py` | keep |
| `scripts/run/run_mu_sweep_fixed_beta_four_radii.py` | main user-facing command | Friendly wrapper for fixed-beta four-radius `mu` sweeps. | `python scripts/run/run_mu_sweep_fixed_beta_four_radii.py` | keep |
| `scripts/run/run_mu_sweep_four_betas_analytic.py` | main user-facing command | Friendly wrapper for analytic-only selected-beta `mu` sweeps. | `python scripts/run/run_mu_sweep_four_betas_analytic.py --betas 15 30 45 60` | keep |
| `scripts/run/run_lambda_mu_fixed_beta_analytic.py` | main user-facing command | Friendly wrapper for fixed-beta analytic `Lambda(mu)` plots with arbitrary mode count and summary CSV. | `python scripts/run/run_lambda_mu_fixed_beta_analytic.py --beta 15 --epsilon 0.0025 --num-modes 7 --num-dashed-lines 7` | keep |
| `scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` | main user-facing command | Friendly wrapper for one analytic determinant-nullspace mode shape. | `python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` | keep |
| `scripts/run/run_tracked_bending_descendant_shape_ru.py` | main user-facing command | Friendly wrapper for one tracked bending descendant shape at explicit `mu`, `epsilon`, `beta`, plot kind, scale, and diagnostics settings. | `python scripts/run/run_tracked_bending_descendant_shape_ru.py` | keep |
| `scripts/run/run_tracked_bending_descendant_shapes_ru.py` | main user-facing command | Friendly wrapper for tracked bending descendant shapes. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_01` | keep |
| `scripts/run/run_branchwise_fem_audit.py` | main user-facing command | Friendly wrapper for the branchwise FEM audit. | `python scripts/run/run_branchwise_fem_audit.py` | keep |
