# Scripts guide

This directory keeps the historical script names working, but the preferred entry points are now in `scripts/run/`.
The original root-level scripts remain in place because many scripts import each other through `scripts.*`; moving them now would require broad import churn with no numerical benefit.

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

`--diagnostics-level all` adds three checks for local components. The projection check reconstructs global components from local right-arm components. The derivative diagnostics estimate whether large local axial displacement also produces large axial strain-like slopes; these are nodal amplitude slopes, not exact FEM strain fields. The arm-energy diagnostics compares axial and bending participation by arm for the selected eigenvector amplitude.

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

This computes arm-wise axial/bending diagnostic energies for the normalized reconstructed analytic shape and for the same FEM mode through the existing arm-energy helper. The values are reported as diagnostic fractions and shares; the arbitrary modal scale is also recorded in the CSV.

Direct FEM element-stiffness energy check:

```bash
python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py --branch-number 5 --beta 30 --mu 0 --epsilon 0.0025 --check-fem-direct-energy
```

This recomputes FEM arm-wise axial/bending energies directly from the current `fem.elem_K(Le)` local stiffness blocks and compares them against `arm_energy_diagnostics(...)`. It is diagnostic postprocessing only and does not modify the FEM baseline.

### Article shape validation at `beta = 15 deg`

- Task: run the analytic-vs-FEM component and energy diagnostics for the article small-angle cases.
- Command: `python scripts/analysis/validate_article_shape_cases_beta15.py`
- Main parameters: `--branch-ids`, `--mus`, `--beta`, `--epsilon`, `--right-coordinate`, `--use-best-orientation`, `--output`, `--output-prefix-root`.
- Results: by default `results/article_shape_validation_beta15_summary.csv` plus one overlay PNG, samples CSV, diagnostics CSV, orientation scan CSV, and energy comparison CSV per branch/mu case.
- Use when: you need a compact review table for `bending_desc_01` and `bending_desc_02` at `beta = 15 deg`, `epsilon = 0.0025`, and the article `mu = 0, 0.1, 0.2` set.
- Do not use when: you need to change determinant entries, FEM baseline behavior, branch tracking, or root matching.

The summary table records the FEM tracked descendant, nearest analytic root index, `Lambda` mismatch, local component relative L2 errors, analytic/FEM axial-energy fractions, right axial shares, and a conservative `conclusion_flag`. The flag is a review aid only; it does not make a physical conclusion automatically.

Examples:

```bash
python scripts/analysis/validate_article_shape_cases_beta15.py
python scripts/analysis/validate_article_shape_cases_beta15.py --branch-ids bending_desc_01 bending_desc_02 --mus 0 --beta 15 --epsilon 0.0025
```

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
| `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py` | analysis/diagnostic | Compare analytic determinant-nullspace local components with one tracked FEM descendant shape. | `python scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py --branch-number 5 --beta 30 --mu 0 --epsilon 0.0025` | keep in `scripts/analysis/` |
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
| `scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` | main user-facing command | Friendly wrapper for one analytic determinant-nullspace mode shape. | `python scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` | keep |
| `scripts/run/run_tracked_bending_descendant_shape_ru.py` | main user-facing command | Friendly wrapper for one tracked bending descendant shape at explicit `mu`, `epsilon`, `beta`, plot kind, scale, and diagnostics settings. | `python scripts/run/run_tracked_bending_descendant_shape_ru.py` | keep |
| `scripts/run/run_tracked_bending_descendant_shapes_ru.py` | main user-facing command | Friendly wrapper for tracked bending descendant shapes. | `python scripts/run/run_tracked_bending_descendant_shapes_ru.py --branch-id bending_desc_01` | keep |
| `scripts/run/run_branchwise_fem_audit.py` | main user-facing command | Friendly wrapper for the branchwise FEM audit. | `python scripts/run/run_branchwise_fem_audit.py` | keep |
