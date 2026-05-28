# CHANGELOG

## 2026-05-28

- Documented the full rigid end-face `Lambda(mu)` EB/Timoshenko/3D FEM
  diagnostic plots as the current main comparative visualization along `mu`,
  recording that the non-tuned 3D engineering joint supports the trend that
  Timoshenko improves agreement as thickness increases while remaining
  diagnostic-only pending mesh convergence and mode-shape review. No
  calculations, Gmsh/CalculiX runs, or result regenerations were performed.
- Added a diagnostic-only full rigid end-face `Lambda(mu)` epsilon sweep for
  equal-thickness coupled rods at `beta=15 deg`, `eta=0`, and
  `epsilon=0.01,0.025,0.05`, writing three separate EB/Timoshenko/3D-FEM PNGs,
  one combined CSV, one Markdown report, and isolated FEM outputs under
  `results/solid_fem_rigid_joint_equal_thickness_beta15_lambda_mu/`. The
  workflow uses descendant analytic branches, a 16-point FEM `mu` grid, the
  full rigid end-face joint with no planar constraint or patch tuning, and
  MAC-based filtering that keeps weak, ambiguous, and duplicate matches out of
  the main plots while preserving them in the CSV/report.
- Added a diagnostic-only equal-thickness rigid-joint frequency trend plot for
  `beta=15 deg`, `mu=0`, `eta=0`, and `epsilon=0.01,0.025,0.05`, consuming
  the existing full rigid end-face 3D FEM thickness-trend CSV without rerunning
  Gmsh/CalculiX. The plot-data CSV, PNG, and report show EB/Timoshenko
  analytic frequencies with MAC-eligible 3D FEM points for the first six
  branches, record weak-MAC exclusions, and keep EB diameter applicability
  separate from the Timoshenko cut-off.
- Added a diagnostic-only rigid end-face 3D FEM thickness-trend audit for
  `beta=15 deg`, `mu=0`, `eta=0`, and
  `epsilon=0.005,0.01,0.025,0.05`, writing isolated full rigid end-face
  point-joint Gmsh/CalculiX outputs under
  `results/solid_fem_rigid_joint_thickness_trend/` plus CSV/Markdown/PNG
  summaries. The audit uses no planar constraint, fused volume, patch
  coupling, or fitted parameter; strong/moderate MAC rows show the mean
  Timoshenko/EB error ratio decreasing from about `1.005` at `epsilon=0.005`
  to about `0.621` at `epsilon=0.05`, while leaving mesh convergence and
  visual mode review pending.
- Added a diagnostic-only 3D FEM joint-constraint bracketing audit for
  `epsilon=0.01`, `mu=0`, `eta=0`, and `beta=0,15,45,90 deg`, writing
  isolated Gmsh/CalculiX outputs under
  `results/solid_fem_joint_constraint_bracketing_audit/` plus CSV/Markdown
  summaries. The study keeps the rigid-end-face and central-patch probes,
  adds translation-only central/average equation probes, a
  rotation-clamped-reference bracket, and a two-reference-node constraint
  attempt; it records that this is not patch-radius calibration, classifies
  translation-only probes as too soft, the rotation clamp as incompatible with
  the analytic free point-joint rotation, and `center_patch_0p5` as diagnostic
  promising but not article-ready.
- Added a diagnostic-only 3D FEM joint-coupling-model audit for
  `epsilon=0.01`, `mu=0`, `eta=0`, and `beta=0,15,45,90 deg`, writing
  isolated Gmsh/CalculiX outputs under
  `results/solid_fem_joint_coupling_models_audit/` plus CSV/Markdown
  summaries comparing full rigid end-face coupling with `0.25*r` and
  `0.5*r` central rigid patches. The `0.5*r` patch is recorded as the current
  best diagnostic candidate, while no tested model is promoted to article-ready
  validation.
- Added a diagnostic-only 3D FEM point-joint discrepancy isolation audit for
  `beta=15 deg`, `mu=0`, `eta=0`, and `epsilon=0.01`, writing isolated
  Gmsh/CalculiX outputs under
  `results/solid_fem_point_joint_mu0_eps0p01_audit/` plus CSV/Markdown
  summaries comparing the current planar-constrained point-joint, full
  point-joint without planar constraints, and a straight two-half-rod
  point-joint sanity case without changing article files, old determinants,
  `formulas.py`, old solvers, existing FEM physical models, baseline results,
  or presentation figures.
- Added a diagnostic-only tau-aware Timoshenko energy partition audit for
  coupled rods at `beta=15 deg`, `epsilon=0.01`, and `eta=0,0.5`, writing
  CSV/Markdown/PNG outputs that track descendants by beta-then-mu shape MAC,
  integrate per-rod bending/shear/axial energy fractions, record EB diameter
  and Timoshenko cutoff margins, and compare matching EB descendants without
  changing article files, old determinants, `formulas.py`, old solvers, FEM
  models, Gmsh/CalculiX workflows, baseline results, or presentation figures.

## 2026-05-27

- Extended the diagnostic-only variable-length Timoshenko helper with explicit
  tau-aware `eta != 0` section scaling and added a sorted-root thickness-
  mismatch audit that verifies eta=0 regression, eta-to-zero continuity,
  `(mu, eta) -> (-mu, -eta)` swap symmetry, the `epsilon -> 0`
  Euler--Bernoulli thickness-mismatch limit, and beta=0 straight composite-rod
  behavior, writing CSV/Markdown outputs without changing article files, old
  determinants, `formulas.py`, old solvers, FEM models, baseline results, or
  presentation figures.
- Added a diagnostic-only variable-length Timoshenko limit audit with a small
  reusable helper, writing CSV/Markdown outputs that verify `eta=0` sorted-root
  checks for `mu=0` consistency against the older equal-rods Timoshenko
  diagnostic, the `epsilon -> 0` Euler--Bernoulli limit, and `beta=0`
  straight-rod `mu`-invariance without changing article files, old
  determinants, `formulas.py`, old solvers, FEM models, baseline results, or
  presentation figures.

## 2026-05-26

- Added a diagnostic-only `mu=0` sanity audit for the equal-thickness
  `beta=15 deg`, `epsilon=0.01`, `eta=0` EB/Timoshenko/planar-FEM overlay,
  documenting the current variable-length Timoshenko formulas, checking
  `Lambda_FEM=sqrt(Omega/epsilon)`, auditing point-joint planar constraints,
  and comparing sorted/MAC-order planar 3D solid modes against analytic and
  1D EB frame-FEM references without changing the FEM model.
- Corrected the diagnostic-only equal-thickness `Lambda(mu)` plot at
  `beta=15 deg`, `epsilon=0.01`, and `eta=0`, writing new corrected outputs
  with a denser FEM `mu` grid, separate EB diameter and Timoshenko cut-off
  applicability audits, mesh-quality/MAC-based FEM point exclusions, and
  large-discrepancy diagnostics without changing article or baseline outputs.
- Added a diagnostic-only `Lambda(mu)` plot script for identical coupled
  circular rods at `beta=15 deg`, `epsilon=0.01`, and `eta=0`, combining EB
  and Timoshenko descendant branches with sparse planar point-joint 3D solid
  FEM points and reporting MAC matching, thin-rod applicability, and EB-vs-
  Timoshenko closeness without changing article or baseline outputs.
- Added a separate diagnostic-only planar-constrained mode for the point-joint
  coupled-rods 3D solid FEM workflow, writing planar-specific report, sorted
  comparison, mode metrics, MAC matrices, and MAC match outputs; the run uses
  an explicitly reported MPC-compatible fallback for CalculiX rigid-body
  dependent nodes and does not replace the full 3D validation benchmark.
- Implemented a separate diagnostic-only mesh-convergence mode for the
  point-joint coupled-rods 3D solid FEM workflow, running beta=15 at
  `epsilon=0.025,0.05` over mesh factors `1.5,1.0,0.75`, writing separate
  convergence report/frequency/MAC/summary outputs with preliminary Hungarian
  assignment rows while leaving baseline and fused outputs untouched.
- Refactored the diagnostic-only point-joint coupled-rods 3D solid FEM workflow
  to pass `beta_deg` explicitly through geometry-dependent helpers, added
  optional partial-run beta/epsilon filters, and documented the rerun controls
  without changing the physical FEM model or fused/article workflows.
- Added a thickness-mismatch FEM validation status document and short index
  links documenting the current 3D solid FEM validation state, article-ready
  roadmap, and candidate figure guidance without recalculating code or results.
- Extended the diagnostic-only point-joint 3D solid FEM workflow with
  centerline MAC-like shape matching between parsed CalculiX solid modes and
  EB/Timoshenko analytic mode reconstructions, writing EB/Timoshenko MAC
  matrices, match summaries, duplicate/low-MAC warnings, and a report section
  that prioritizes shape matching over simple classified-mode order.
- Added a diagnostic-only point-joint 3D solid FEM workflow for coupled equal
  rods using two separate cylinder meshes, a shared CalculiX rigid reference
  node for the inner end faces, a local `*RIGID BODY` capability probe, modal
  runs/classification by beta, and high-level comparison against the previous
  fused-cylinder diagnostic outputs without changing those fused outputs.
- Generalized the diagnostic-only coupled equal-rods 3D solid FEM workflow with
  a multi-beta sweep (`15`, `45`, `90` degrees), fused-joint overlap geometry
  audit, connected-mesh and fixed-node-set summaries, classified in-plane
  comparison by beta, and aggregate cut-off/thin-rod reporting while leaving
  the beta=15-specific workflow and existing outputs intact.
- Added a diagnostic-only beta=15 equal-coupled-rods 3D solid FEM workflow
  using portable Gmsh/CalculiX resolution, fused-cylinder mesh generation,
  coordinate-derived outer fixed node sets, CalculiX modal runs, `.frd`
  mode-shape classification, sorted and classified in-plane comparison CSVs,
  cut-off checks, and thin-rod applicability reporting.
- Added diagnostic-only post-processing for the single fixed-fixed circular-rod
  3D solid FEM workflow: CalculiX `.frd` displacement mode parsing, nodal
  slice-based mode-shape metrics, bending/axial/torsion-like classification,
  classified bending-doublet pairing, and separate mode-metric/doublet CSVs.
- Extended the diagnostic-only single fixed-fixed circular-rod 3D solid FEM
  workflow with portable CalculiX resolution through `CCX_EXE`, coordinate-
  derived fixed-end node sets, sanitized CalculiX modal inputs, actual ccx
  modal runs when available, and preliminary EB/Timoshenko/3D solid comparison
  output.
- Added mesh sanity summaries to the diagnostic-only single fixed-fixed
  circular-rod 3D solid FEM workflow, reporting node/solid-element counts,
  bounding boxes, inferred dimensions, end-face physical-group detection, and
  captured Gmsh warning/error lines before any CalculiX modal solve is used.
- Updated the diagnostic-only single fixed-fixed circular-rod 3D solid FEM
  workflow to resolve a portable Gmsh executable through `GMSH_EXE`, a
  user-editable script default, or PATH, reporting the resolved path/version
  and generating meshes without requiring a system PATH change.
- Added a diagnostic-only single fixed-fixed circular-rod 3D solid FEM workflow
  seed with local Gmsh/CalculiX/Code_Aster/Salome-Meca audit notes, optional
  Gmsh and CalculiX input generation, one-rod EB/Timoshenko reference
  comparison rows, and graceful no-solver reporting without making any
  external solver a project dependency.
- Extended the diagnostic-only single-rod and coupled equal-rod Timoshenko
  comparison scripts with cut-off frequency columns and report warnings,
  added a small coupled equal-rod kappa sensitivity CSV/report section, and
  documented that the current Timoshenko diagnostic changes only the flexural
  part while keeping the axial force law classical.

## 2026-05-25

- Added a diagnostic-only coupled equal-radius circular-rod comparison between
  the existing Euler--Bernoulli determinant and a trial explicit-function
  Timoshenko model at `beta=15 deg`, `mu=0`, and `eta=0`, with CSV/PNG/Markdown
  outputs, sorted-frequency warnings, the circular shear coefficient, and
  thin-rod validity marking without changing article files, figures,
  determinants, old solvers, formulas.py, baseline results, or the FEM model.
- Added a diagnostic-only single fixed-fixed circular-rod comparison between
  Euler--Bernoulli and exact state-space Timoshenko beam frequencies, including
  CSV/PNG/Markdown outputs, the Diaz-de-Anda circular shear coefficient, and
  thin-rod validity warnings without changing article files, figures,
  determinants, old solvers, formulas.py, results baselines, or the FEM model.
- Added local Timoshenko/shear-coefficient literature integration notes for
  eight imported PDFs, including circular-rod coefficient selection, critical
  frequency/second-spectrum cautions, frame-oriented Timoshenko sources, and
  bibliography/source-index updates without changing code, article files,
  figures, solvers, determinants, FEM model, or results.
- Added a parameterized diagnostic-only thickness-mismatch branch-shape overlay
  script for comparing descendant analytic mode shapes across eta values at
  selected mu points, reusing the eta determinant, analytic reconstruction
  helpers, and shape-MAC descendant tracking without changing article files,
  figures, old determinants, old solvers, results, or the FEM model.
- Added project-wide Script Proliferation Control guidance, linked it from the
  agent and script guides, and shortened duplicated thickness-mismatch script
  and global-rule documentation without changing code, scripts, article files,
  figures, determinants, solvers, results, or FEM model behavior.

## 2026-05-24

- Added a conservative documentation cleanup for the diagnostic
  thickness-mismatch series: a script audit/navigation map, a diagnostic-to-
  article workflow note, a thickness-mismatch / Timoshenko article-writing
  skill, and a skeleton `paper_thickness_mismatch_timoshenko/` article folder
  without moving old scripts or changing numerical models.
- Added three diagnostic-only thickness-mismatch plots for the beta-15,
  epsilon-0.0025 eta=0.5 study: a descendant `Lambda(mu)` eta sweep for
  `eta=-0.5, 0, 0.5`, an eta=0.5 analytic descendant vs FEM frequency
  overlay, and an eta=0.5 descendant plot with clamped-supported isolated-rod
  reference curves, all preserving the diameter-based applicability rule and
  branch-descendant convention.
- Added a diagnostic-only beta-45 eta=0.5 thickness-mismatch `Lambda(mu)` plot
  with descendant branches and the same clamped-supported isolated-rod
  reference convention, including the `sqrt(tau_i)` thickness-mismatch Lambda
  normalization and a cautious beta-45 vs beta-15 nearest-reference-grid
  observation.
- Added a separate diagnostic-only beta-45 eta=0.5 isolated-rod plot using
  clamped-clamped / fixed-fixed single-rod references, preserving the same
  descendant tracking, thin-rod applicability split, and `sqrt(tau_i)`
  thickness-mismatch Lambda normalization while writing distinct output files.
- Added `docs/project_rules.md` as the central project-wide rules document and
  linked it from the root README, theory assumptions, thickness-mismatch docs,
  veering docs, and scripts helper guides.
- Added a diagnostic-only large-eta thickness-mismatch `Lambda(mu)` plot for
  `beta=15 deg`, `epsilon=0.0025`, and `eta=-0.5, 0, 0.5`, splitting tracked
  branches into solid valid segments and dashed segments where the diameter
  criterion `2r_i/l_i <= 0.1` is violated.
- Added a diagnostic-only FEM check for the mass-preserving thickness-mismatch
  model at `eta=0.5`, `beta=15 deg`, and `epsilon=0.0025`, comparing analytic
  sorted/tracked `Lambda` values against a separate two-radius Euler-Bernoulli
  frame FEM assembly with a coarse/refined mesh check and branch-5/6 shape MAC
  diagnostics.
- Added common thickness-mismatch helpers for the diameter-to-length validity
  check and documented the diagnostic rule that curve segments must be dashed
  with a warning/report when `2r_i/l_i > 0.1` for either rod.
- Replaced nearest-frequency `Lambda(mu)` continuation in thickness-mismatch
  diagnostics with analytic shape-MAC tracking plus low-MAC warnings,
  added nearest-frequency disagreement warnings for the large-eta beta-15 plot,
  track ten branches internally while plotting the first seven for the
  large-eta beta-15 case, and updated the eta=0.5 FEM check to compare against
  MAC-tracked analytic branches without changing the determinant, solvers, FEM
  model, article files, or article figures.
- Documented the thickness-mismatch branch convention that branch numbers are
  descendant mode-shape identities while sorted positions are diagnostic
  metadata, added suspicious-assignment flags for large sorted-position jumps,
  and added a refined eta=0.5 branch-identity audit for descendants 5--7.
- Added an eta=0.5 global thickness-mismatch spectrum overview diagnostic that
  compares the first eight sorted roots with the first eight descendant
  branches, plots accepted canonical sorted positions, records unresolved
  candidate assignments, tracking warnings, and the diameter-based thin-rod
  validity status, and writes CSV/PNG/Markdown outputs.
- Updated thickness-mismatch MAC tracking so low-MAC, low-margin, or large-jump
  candidate assignments remain unresolved diagnostics instead of changing the
  accepted descendant identity or canonical sorted position.

## 2026-05-23

- Added a diagnostic-only mass-preserving thickness-mismatch analytic model with
  `eta=(r2-r1)/(r1+r2)`, a separate determinant implementation in
  `src/my_project/analytic/formulas_thickness_mismatch.py`, eta=0/root/swap
  checks, Lambda(eta) diagnostic plots, and
  `docs/thickness_mismatch/README.md`, leaving the baseline equal-radius
  determinant, FEM model, article files, and article figures unchanged.
- Added `scripts/analysis/track_lambda_eta_thickness_mismatch.py` to seed the
  first six thickness-mismatch branches at `eta=0`, continue them to positive
  and negative eta by unique nearest-root matching, and write tracked-vs-sorted
  CSV/PNG diagnostics plus `results/thickness_mismatch_lambda_eta_tracking_report.md`
  without changing the baseline determinant, solvers, FEM model, article files,
  or article figures.
- Added `scripts/analysis/track_lambda_mu_thickness_mismatch_eta_sweep.py` to
  track the first six mass-preserving thickness-mismatch `Lambda(mu)` branches
  from `mu=0` to `mu=0.9` for `eta=-0.1, 0, 0.1`, writing diagnostic CSV/PNG
  outputs and `results/thickness_mismatch_lambda_mu_eta_sweep_tracking_report.md`
  with eta=0 consistency, sorted-index switch, gap, jump, and eta-sensitivity
  summaries.
- Updated `scripts/analysis/track_lambda_mu_thickness_mismatch_eta_sweep.py`
  into a single-beta, user-parameterized plotting script with dynamic
  beta/epsilon output names and PNG-only output by default; CSV and Markdown
  report files are now optional debug outputs.
- Added `scripts/analysis/plot_thickness_mismatch_branch5_shapes.py` to plot
  diagnostic full deformed shapes for tracked thickness-mismatch branch 5 at
  `beta=15 deg`, `epsilon=0.0025`, comparing `eta=-0.1, 0, 0.1` for
  `mu=0` and `mu=0.1` without changing the baseline determinant, solvers, FEM
  model, article files, or article figures.
- Added `scripts/analysis/plot_article_fig3_with_fp_and_ff_refs.py` to build a
  results-only diagnostic version of article Figure 3 with the original
  clamped-pinned single-rod references and an additional clamped-clamped
  single-rod reference family, writing only
  `results/article_fig3_with_fp_and_ff_refs.png` and the matching CSV.

## 2026-05-18

- Restored automatic equation numbering and `\eqref`/`\ref` cross-references in
  `paper_dorofeev_style/main.tex` after the article file rename, preserving
  formula content, figure files, bibliography, and mathematical model.
- Converted manually tagged equations in
  `paper_dorofeev_style/dorofeev_article_journal_template.tex` to automatic
  LaTeX equation numbering with `\label`/`\eqref`, preserving formula content,
  order, notation, figures, bibliography, and article conclusions.
- Rebalanced article Figure 8 layout in
  `paper_dorofeev_style/generate_descendant5_vibration_shape_figures.py` by
  padding only panel view limits to a shared aspect before saving
  `Dorofeev.Fig.8.eps` and a preview PDF, keeping beta, mu, epsilon,
  mode-scale, branch tracking, reconstructed shapes, and beta=45 deg geometry
  unchanged.
- Added `paper_dorofeev_style/dorofeev_article_journal_template.tex` as a
  single-file journal-template article draft assembled from the current
  article text, final EPS figures, and manual literature/References lists, with
  missing required metadata marked as `to do` and no changes to formulas,
  plotted data, figures, determinant, solver logic, FEM model, or branch
  tracking.
- Finalized the current article figure set as journal-ready EPS files
  `Dorofeev.Fig.1.eps` through `Dorofeev.Fig.8.eps`, updated article LaTeX
  paths to those final names, removed preview/test/variant article figure
  outputs, and converted multipanel article figures into one EPS file per
  figure without changing plotted data, determinant, solver logic, FEM model,
  or branch tracking.
- Fixed Figure 2 marker styling by making the colored FEM markers edge-free in
  the final `Dorofeev.Fig.2.eps`, preserving marker size, point count,
  coordinates, analytic curves, plotted data, determinant, solver logic, FEM
  model, and branch tracking.
- Improved article Figure 2 EPS output quality by saving the primary EPS
  directly from matplotlib at 600 dpi and tightening Figure 2 line/tick/spine
  styling while keeping all plotted data, branch identities, FEM markers,
  determinant, solver logic, and branch tracking unchanged.
- Updated article Figure 2 generation to write `Dorofeev.Fig.2.eps`, use
  math-only axis labels `\beta` in degrees and `\Lambda`, and render the
  right-edge branch numbers in black italic without changing the analytic/FEM
  data, branch set, determinant, solver logic, FEM model, or branch tracking.

## 2026-05-13

- Updated the descendant-5 vibration-shape article generator to use the same
  axis-free panel styling as the local `target_descendants_beta15_r5` figures,
  and added an `--axis-off` styling flag to the shared desc05 full-shape plotter
  to hide axes, ticks, frame, labels, and grid without changing calculations or
  plotted shapes.
- Added `paper_dorofeev_style/generate_descendant5_vibration_shape_figures.py`
  to regenerate descendant-5 vibration-shape article figures into
  `paper_dorofeev_style/figures/` from a top-of-file `USER PARAMETERS` block,
  using the existing desc05 analytic full-shape workflow with article styling
  defaults of no in-figure title and no legend; the shared plotter now exposes
  `--no-title` and `--no-legend` styling flags, with no determinant,
  `formulas.py`, `solvers.py`, `python_fem.py`, FEM physical model, branch
  tracking, shape reconstruction, or mathematical changes.

## 2026-05-11

- Updated `scripts/run/run_analytic_coupled_rods_vibration_shapes_beta15_mu06_eps0025_001_ru.py`
  so beta, mu, epsilon/root-label cases, output path, title, mode scale, DPI,
  and tracking-grid parameters are edited through a top-of-file `USER
  PARAMETERS` block; no determinant, `formulas.py`, `solvers.py`,
  `python_fem.py`, FEM physical model, branch tracking, shape reconstruction,
  or mathematical changes were made.
- Added `scripts/run/run_analytic_coupled_rods_vibration_shapes_beta15_mu06_eps0025_001_ru.py`
  for the two-case descendant-5 vibration-shape overlay at beta=15 deg and
  `mu=0.6` (`epsilon=0.0025, root=6` and `epsilon=0.01, root=5`), and extended
  the existing desc05 analytic full-shape sweep plotter with optional
  case-root legend labels, single-PNG `--output`, title-prefix control, and
  geometry-legend hiding; the same analytic branch tracking, shape
  reconstruction, normalization, root search, and sampling grid are reused, with
  no determinant, `formulas.py`, `solvers.py`, `python_fem.py`, FEM physical
  model, or mathematical changes.
- Updated the article-local fixed-fixed reference figure generator so its beta
  values, sweep parameters, compact y-limit, and axis labels are edited through
  a top-of-file configuration block; the article x-axis label is now
  `Параметр μ`, with no determinant, `formulas.py`, `solvers.py`,
  `python_fem.py`, FEM physical model, branch tracking, or mathematical
  changes.
- Added `paper_dorofeev_style/generate_fixed_fixed_reference_spectral_figures.py`
  so the compact beta=15 deg and beta=45 deg fixed-fixed reference
  `Lambda(mu)` article figures can be regenerated directly into
  `paper_dorofeev_style/figures/`, reusing the existing analytic-only plotting
  workflow without determinant, `formulas.py`, `solvers.py`, `python_fem.py`,
  FEM physical model, branch tracking, or mathematical changes.
- Updated the compact fixed-fixed reference plotting workflow to generate
  separate article-style `Lambda(mu)` outputs for beta=15 deg and beta=45 deg,
  using Russian axis labels, no in-figure title, no legend, compact
  `ylim = [0, 11.2]`, solid coupled branches for `mu <= 0.6`, dashed coupled
  branches for `mu > 0.6`, and full-range dotted fixed-fixed bending
  references only; no determinant, `formulas.py`, `solvers.py`,
  `python_fem.py`, FEM physical model, branch tracking, or mathematical
  changes were made.
- Updated `scripts/run/run_lambda_mu_beta15_eps001_fixed_fixed_ref.py` to
  write separate compact article-style outputs with `ylim = [0, 11.2]`, keep
  fixed-fixed reference lines dotted over the full `mu` range, and draw the
  same tracked coupled branches solid for `mu <= 0.6` and dashed for
  `mu > 0.6`; the main CSV now records `region_style`, with no determinant,
  `formulas.py`, `solvers.py`, `python_fem.py`, FEM physical model, branch
  tracking, or mathematical changes.
- Added `scripts/run/run_lambda_mu_beta15_eps001_fixed_fixed_ref.py` to plot
  the first six beta=15 deg, epsilon=0.01 analytic `Lambda(mu)` branches with
  six horizontal fixed-fixed single-beam `L = 2 m` reference lines, writing
  PNG, main-curve CSV, and reference-root CSV outputs; the old CS
  variable-length reference families and existing fixed-beta defaults remain
  unchanged, with no determinant, `formulas.py`, `solvers.py`,
  `python_fem.py`, FEM physical model, or article figure changes.
- Added `scripts/analysis/plot_lambda_vs_beta_fixed_mu.py` to plot
  canonically tracked analytic `Lambda(beta)` branches at fixed `mu` values,
  writing sorted-root CSVs, focused branch-path CSVs, per-`mu` PNGs, and a
  gap/exchange summary report for beta-driven veering diagnostics; no
  determinant, `formulas.py`, `solvers.py`, `python_fem.py`, FEM physical
  model, article prose, or article figure styling changes were made.
- Updated the `mu -> 1` single-rod limit diagnostic to compute and plot 12
  coupled-system roots by default, matching the 12-mode FP/FF reference set and
  making the beta 10/45/90 plots visually comparable; the diagnostic root scan
  now uses a finer `scan_step=0.01` to recover the 12th displayed root in the
  high-mu cases without changing the determinant, `formulas.py`, `solvers.py`,
  `python_fem.py`, FEM physical model, article prose, or article figure styling.
- Expanded `scripts/analysis/check_mu_to_one_single_rod_limit.py` so the
  clamped-pinned and clamped-clamped single-rod reference roots are computed
  numerically through 12 modes, validated against known FP/FF prefixes, and
  reported with an explicit comparison against the earlier first-five-reference
  diagnostic; no determinant, `formulas.py`, `solvers.py`, `python_fem.py`,
  FEM physical-model, article prose, or article figure-styling changes were
  made.
- Added `scripts/analysis/check_mu_to_one_single_rod_limit.py` to compare the
  thin-rod `mu -> 1` coupled-system roots against clamped-pinned and
  clamped-clamped single-rod references using both finite `L2 = 1 + mu` and
  limiting `L = 2` lengths, with sorted-root CSVs, optional FEM cross-check,
  branch-focused CSV, diagnostic plots, and a summary report; no determinant,
  `formulas.py`, `solvers.py`, `python_fem.py`, FEM physical-model, article
  prose, or article figure-styling changes were made.
- Documented the public FEM right-arm transform convention in `AGENTS.md` and
  `scripts/README.md`, and kept generated/private article artifacts ignored;
  no determinant, `formulas.py`, `solvers.py`, FEM physical-model, article
  prose, or article figure-styling changes were made by this cleanup pass.

## 2026-05-06

- Corrected the production FEM right-arm transform convention in `src/my_project/fem/python_fem.py` to the local-to-global map `q_global = R(+beta) q_local`, assembling right-arm element matrices as `K_global = T K_local T.T` and `M_global = T M_local T.T`; updated related residual/joint/energy diagnostics and added `tests/test_fem_right_transform_convention.py`, with no analytic determinant, formulas, solvers, physical parameters, or analytic eigenfrequency changes.
- Added `scripts/analysis/check_fem_right_transform_variants.py` to reassemble diagnostic-only right-arm FEM transform variants and compare determinant-null shapes against each corrected-convention variant, writing variant audit and local-direction sanity CSVs with no `python_fem.py`, determinant, formula, solver, production-eigenfrequency, or physical-model changes.
- Added `scripts/analysis/check_beta0_single_rod_limit.py` to verify that the current determinant roots and SVD null-vector reconstruction reduce to the exact single clamped-clamped rod of total length 2 at `beta = 0`, writing root, shape, sample, and mu-invariance CSV diagnostics with no determinant, formula, solver, FEM baseline, physical-model, or eigenfrequency changes.
- Added `--derive-determinant-as-joint-system` to `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py`, deriving the current determinant force equations as a `q_joint=[ux,uy,theta]` operator on the determinant kinematic nullspace and comparing current/sign/order variants against analytic/FEM joint dynamic stiffness, with no determinant, formula, solver, FEM baseline, physical-model, or eigenfrequency changes.
- Added `--audit-determinant-nullspace-conditioning` to `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py`, comparing raw determinant SVD coefficients against coefficients reconstructed from the joint-displacement dynamic-stiffness null vector and auditing row, column, row+column, and axial physical column scalings, with no determinant, formula, solver, FEM baseline, physical-model, or eigenfrequency changes.
- Added `--compare-joint-nullspace` to `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py`, comparing `D_total_FEM` and `D_total_analytic` directly on `q_joint=[ux,uy,theta]`, writing matrix, summary, vector, pairwise, and reconstruction-residual CSV diagnostics while avoiding determinant internal-force row conventions, with no determinant, formula, solver, FEM baseline, physical-model, or eigenfrequency changes.
- Added `results/right_arm_force_projection_sign_audit.md`, deriving the right-arm mixed force projection signs from the FEM virtual-work convention `K_global = T.T K_local T`, identifying candidate future sign flips for `N2*sin(beta)` in the transverse row and `Q2*sin(beta)` in the axial row, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a force-row scaling audit behind `--audit-force-row-scaling` in `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py`, decomposing determinant force rows into FEM-compatible nondimensional `M`, `Q`, and `N` endpoint contributions and writing scaling, summary, and diagnostic-matrix CSVs, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a constrained determinant force-row audit behind `--constrained-row-audit` in `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py`, restricting coefficient comparisons to the SVD nullspace of the kinematic compatibility rows before comparing determinant and dynamic-stiffness force operators, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a determinant force-row audit behind `--row-audit` / `--audit-determinant-rows` in `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py`, writing row-level, row-scaling, force-contribution, and null-vector CSVs for current determinant rows, independent full-basis endpoint rows, and dynamic-stiffness-induced force-row variants, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `scripts/analysis/compare_joint_dynamic_stiffness_analytic_fem.py` to compare per-arm FEM Schur-complement joint dynamic stiffness with an independent full-basis analytic arm response and determinant force-row cross-checks, writing `joint_dynamic_stiffness_compare_*.csv`, `*_summary.csv`, and `*_determinant_crosscheck.csv`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a frequency-scaling audit behind `--audit-frequency-scaling` in `scripts/analysis/check_analytic_shape_in_fem_residual.py`, writing `*_frequency_scaling_audit.csv` to compare direct residuals under `Lambda`, `Lambda^2`, `Lambda^4`, FEM `omega`, FEM `omega^2`, and Rayleigh spectral factors, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `scripts/analysis/check_single_beam_exact_shape_residual.py` to sanity-check the FEM residual diagnostic on exact clamped-free single-beam bending and axial-bar modes using the current `python_fem.py` element matrices, writing `results/single_beam_exact_shape_residual_bending.csv` and `results/single_beam_exact_shape_residual_axial.csv`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a right-arm local/global transform convention scan behind `--scan-right-transform-conventions` in `scripts/analysis/check_analytic_shape_in_fem_residual.py`, writing `*_right_transform_convention_scan.csv` and `*_right_transform_sanity.csv` for `T`, `T.T`, `R(+beta)`, `R(-beta)`, inverse-`T`, and FEM-energy-consistent embedding variants, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added an axial amplitude scale scan behind `--scan-axial-scales` in `scripts/analysis/check_analytic_shape_in_fem_residual.py`, writing `*_axial_scale_scan.csv` with physical candidate scale factors, L2-optimal local axial factors, coarse two-arm log-grid variants, Rayleigh/residual/MAC metrics, residual-share diagnostics, and direct element axial-energy fractions, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added an axial-component convention scan behind `--scan-axial-conventions` in `scripts/analysis/check_analytic_shape_in_fem_residual.py`, writing `*_axial_convention_scan.csv` with right-arm axial sign/reversal variants, left-arm sign flips, Rayleigh/residual metrics, axial-strain comparisons, and direct element axial-energy diagnostics, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a full `cos`/`sin`/`cosh`/`sinh` bending-basis coefficient pairing scan behind `--scan-bending-basis-pairings` in `scripts/analysis/check_analytic_shape_in_fem_residual.py`, writing `*_bending_basis_pairing_scan.csv` and `*_full_basis_column_audit.csv` to audit identity/swap/sign pairings and right-arm `z` sign variants without Krylov functions, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added Rayleigh quotient, FEM modal projection, residual DOF-group, and element-wise residual/energy localization diagnostics behind `--localize-residual` in `scripts/analysis/check_analytic_shape_in_fem_residual.py`, writing `*_modal_projection.csv`, `*_residual_groups.csv`, and `*_element_residual_energy.csv`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added a rotational DOF convention scan to `scripts/analysis/check_analytic_shape_in_fem_residual.py`: the residual diagnostic now compares FEM theta DOFs against finite-difference transverse slopes, can scan analytic theta embeddings based on `dw/ds`, `dw/dxi`, and `dw/dz` sign/scale variants, and writes `*_theta_convention_scan.csv`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `scripts/analysis/check_analytic_shape_in_fem_residual.py` to embed a canonically tracked analytic determinant-nullspace shape directly into the baseline FEM DOF vector, compute direct FEM residuals, compare global/translational/rotational/arm/joint DOFs, and write summary plus node debug CSVs, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added adaptive refinement and failure diagnostics to `scripts/lib/analytic_branch_tracking.py`: low-MAC transitions now bisect beta/mu steps before failing, missed sign-change roots can be recovered locally by minimizing the determinant matrix smallest singular value for tracking candidates, failure reports include per-candidate MAC/frequency/cost CSVs under `results/debug/`, and the desc05 plotter exposes refinement controls, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Improved low-MAC diagnostics for analytic branch tracking and `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py`: beta/mu continuation steps are explicit user parameters, fail-fast messages now include branch id, beta, mu, MAC, current sorted index, Lambda, and a refinement recommendation, exploratory `--allow-low-mac` runs are labeled in stdout, and tracking debug CSVs include an explicit `Lambda` column under `results/debug/`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Fixed `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py` so the historical known-contradiction diagnostic is opt-in behind `--check-known-contradiction` and skips gracefully when the current tracking result does not contain `beta = 15`, `epsilon = 0.0025`, `mu = 0.8`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `tests/test_analytic_branch_tracking_regression.py` for the `bending_desc_05`, `beta = 15 deg`, `epsilon = 0.0025` tracking inconsistency, made low-MAC analytic assignments fail fast unless `allow_low_mac` / `--allow-low-mac` is explicitly requested, added branch-Lambda/current-index consistency assertions, and documented that low-MAC results are not canonical, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Added `scripts/run/run_lambda_mu_fixed_beta_analytic.py` as a convenient fixed-`beta` analytic `Lambda(mu)` plotting entrypoint for arbitrary `--num-modes`, preserving the existing `FreqMuNet.py` visual style, using shared canonical analytic branch tracking, and writing a PNG plus summary CSV such as `results/lambda_mu_beta15_eps0p0025_modes7_dash7.*`, with no determinant, FEM baseline, formula, solver, physical-model, or eigenfrequency changes.
- Introduced `scripts/lib/analytic_branch_tracking.py` as the shared source-of-truth helper for analytic branch identity from `beta = 0`, `mu = 0` for each `epsilon`, separated `branch_id` from `current_sorted_index`, made tracking CSVs optional debug artifacts only, synchronized the desc05 analytic shape sweep and `FreqMuNet.py` Lambda(mu) branch selection through the same helper, and documented the rule in `AGENTS.md` and `scripts/README.md`, with no determinant, FEM baseline, formula, solver, FEM tracking, or eigenfrequency changes.
- Updated `scripts/analysis/plot_desc05_full_shapes_beta15_eps_sweep.py` so the fifth-descendant epsilon-sweep figures use analytic-only branch tracking from `beta = 0`, `mu = 0` for each `epsilon` independently, select roots by sampled-shape MAC instead of FEM Lambda, label `current_sorted_index` instead of branch identity as a root number, and save full tracking-path CSV diagnostics only behind `--save-tracking-debug`, with no determinant, FEM baseline, formula, solver, FEM tracking, or eigenfrequency changes.
- Updated `scripts/analysis/validate_article_shape_cases_beta15.py` to use explicit corrected direct FEM element-stiffness energy fields, add transverse-only and axial-only shape metrics, split frequency/transverse/full-component/energy flags, and write `results/article_shape_validation_beta15_report.md`, with no determinant, FEM baseline, formula, solver, root, tracking, or eigenfrequency changes.
- Fixed FEM arm-wise energy diagnostics by making `arm_energy_diagnostics(...)` use direct local element-stiffness energy blocks for both arms, preserving the legacy nodal-gradient implementation as `legacy_arm_energy_diagnostics(...)` for reference and adding a regression smoke test for `bending_desc_05`, `beta = 30 deg`, `mu = 0`, `epsilon = 0.0025`, with no FEM baseline, determinant, formula, solver, root, tracking, or eigenfrequency changes.
- Added `scripts/analysis/validate_article_shape_cases_beta15.py` for batch analytic-vs-FEM component and energy validation of the article `beta = 15 deg`, `epsilon = 0.0025` small-angle cases (`bending_desc_01`, `bending_desc_02`, `mu = 0, 0.1, 0.2` by default), writing `results/article_shape_validation_beta15_summary.csv` plus per-case overlay/diagnostic files with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added `scripts/run/run_analytic_coupled_rods_mode_shape_ru.py` and `scripts/lib/analytic_coupled_rods_shapes.py` for analytic mode-shape reconstruction directly from the determinant nullspace, including Russian-labeled full/transverse/components plots, samples CSV, diagnostics CSV, endpoint residuals, and analytic arm-energy shares, with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added a direct FEM element-stiffness energy check to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_fem_direct_energy_check.csv` behind `--check-fem-direct-energy` and comparing `arm_energy_diagnostics(...)` against energies recomputed from `fem.elem_K(Le)` local axial/bending blocks, with no FEM baseline or `python_fem.py` changes.
- Added analytic-vs-FEM arm-wise axial/bending energy comparison diagnostics to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_energy_comparison.csv` behind `--compare-energies` and appending key energy fields to the diagnostics CSV, with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added endpoint/matrix consistency diagnostics to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_endpoint_consistency.csv` and an optional text summary for analytic null-vector rows, analytic basis columns, external clamps, and FEM joint kinematics, with no determinant, FEM, formula, solver, root, tracking, or eigenfrequency changes.
- Added an orientation/sign-convention scan to `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py`, writing `*_orientation_scan.csv`, printing the top matching variants, and keeping the main overlay on the current convention unless `--use-best-orientation` is explicitly supplied, with no determinant, FEM, formula, solver, root, or tracking changes.
- Added `scripts/analysis/compare_analytic_fem_tracked_descendant_shape.py` for analytic-vs-FEM local component comparison of tracked descendants; it reconstructs the analytic mode shape from the determinant nullspace and changes no determinant, FEM, formulas, root-finding, or branch tracking.

## 2026-05-05

- Added `--diagnostics-level all` to the single tracked bending descendant shape runner, with local projection reconstruction checks, strain-like local derivative diagnostics, and arm-wise axial/bending diagnostic energies for the selected mode shape, without changing determinant, FEM, formulas, root-finding, branch tracking, or eigenfrequencies.
- Added an editable `USER PARAMETERS` block to the single tracked bending descendant shape runner so no-argument runs use file-local working defaults while CLI arguments remain available as overrides, without changing determinant, FEM, formulas, root-finding, branch tracking, or eigenfrequencies.
- Expanded the single tracked bending descendant shape runner CLI with explicit branch defaults, `--l-total`, `--dpi`, `--figsize`, `--show`, `--normalize auto`, deterministic scale-aware output names, and optional one-row diagnostics CSV export, without changing determinant, FEM, formulas, root-finding, branch tracking, eigenfrequencies, or parameter meanings.
- Added `full`, `transverse`, and `components` plot modes plus `--mode-scale` and normalization controls to the single tracked bending descendant shape runner; the new diagnostics report local axial/transverse amplitudes without changing determinant, FEM, formulas, root-finding, branch tracking, eigenfrequencies, or baseline mode-shape normalization.
- Generalized Russian title-label inference for tracked bending descendants so any `bending_desc_NN` branch id renders as `потомок N-й изгибной ветви`, with a unit-style fallback check and no determinant, FEM, formula, root-finding, tracking, normalization, layout, or color changes.
- Cleaned up the single-case tracked bending descendant shape PNG layout by removing its top legend, shortening the title to the requested beta/mu/epsilon and Russian branch label metadata, and reducing the single-case figure height, without changing tracking, normalization, scale, colors, formulas, determinant, root-finding, or FEM logic.
- Split tracked bending descendant shape plotting into a reusable `scripts/lib/tracked_bending_descendant_shapes.py` core used by both the single-shape runner and the multi-panel builder; multi-panel figures now call shared in-process helpers instead of depending on single-case CLI behavior, with no determinant, FEM, formula, root-finding, or branch-tracking changes.
- Added `scripts/run/run_tracked_bending_descendant_shape_ru.py` as a single-shape user-facing runner for tracked bending descendant mode shapes; it exposes descendant branch number/id, `mu`, `epsilon`, and `beta`, without changing determinant, FEM, formula, root-finding, or branch-tracking logic.
- Introduced `scripts/plot_tracked_bending_descendant_shapes_ru.py` as the parameterized tracked bending-descendant mode-shape plotting entrypoint; converted the old `scripts/plot_flat_mu_bending_desc*.py` scripts into compatibility wrappers that preserve their branch ids and output PNG paths, without changing determinant, FEM, formula, or branch-tracking logic.
- Added a `scripts/README.md` script map, introduced clear user-facing wrappers in `scripts/run/`, and documented analysis/helper/legacy groupings while preserving old root-level command paths; mathematical model, determinant/root logic, and FEM baseline were unchanged.

## 2026-05-02

- Updated `src/my_project/analytic/FreqMuNet.py` as the standalone fixed-`beta`
  `Lambda(mu)` plotting entrypoint: added CLI controls for `--beta`,
  `--epsilon`, `--num-modes`, dashed CS reference-line counts, `mu` range, and
  output path, without changing determinant assembly, shared root-finding, or
  FEM/data files.

## 2026-04-22

- Integrated four new veering/localization literature sources into the project
  literature system: added notes, bibliography entries, and source-index blocks
  for Manconi--Mace 2017, Ehrhardt et al. 2018, Lacarbonara et al. 2005, and
  Fontanela et al. 2021; updated the veering assessment/terminology only to
  clarify mechanism-vs-geometry priority and the weak-coupling versus
  strong-coupling distinction, with no solver-code, script, test, data, figure,
  or manuscript edits.

## 2026-04-21

- Added `nair_1973_quasi_degeneracies` to the project literature system and
  recorded its use as a mechanism-level source for quasi-degeneracy / rapid
  modal-reorganization vocabulary; updated veering notes with a conservative
  `bending_desc_04` tracked-pair search from existing `beta = 15 deg`,
  `r = 5 mm` data only, without new solves, code changes, or manuscript edits.

- Синхронизирована проектная теория с текущим текстом статьи в `paper_dorofeev_style/`: в `docs/theory/main_note.md` и `docs/theory/assumptions.md` зафиксированы смысл параметра `\mu`, нормировка `\Lambda`, эталонная роль `\mu=0`, пунктирные CS reference-линии и необходимость интерпретировать `\Lambda(\mu)` вместе с формами колебаний без изменения формул, determinant, solver logic, FEM baseline или текста статьи.

## 2026-04-20

- Переписано введение статьи в `paper_dorofeev_style/sections/introduction.tex`: убраны неподтверждённые обобщения, встроены реальные ссылки по месту, а работы `perkins1986`, `pierre1988` и `liu2002` теперь используются как источники по общим спектральным механизмам `veering`/`mode localization`, а не как прямые аналоги текущей геометрии.
- Обновлена локальная библиография статьи в `paper_dorofeev_style/bib/references.bib`: добавлены записи `liu2002derivatives` и `bauer2025coupledrods`, причём ссылка на предыдущую нашу статью теперь привязана к DOI `10.24412/0136-4545-2025-3-73-81`.

## 2026-04-19

- Rebuilt article Figure 2 in `paper_dorofeev_style/generate_article_spectral_figures.py` around branch identity instead of the smoothed sorted-spectrum CSV: the figure now seeds working branch numbers from the `beta = 0` bending descendants, tracks FEM descendants across `0 <= beta <= 90 deg` by MAC+frequency continuity, matches analytic roots to those tracked descendants at each `beta`, and therefore preserves branch colors/numbers through real crossings without changing the determinant, shared solver logic, or baseline FEM.
- Updated `paper_dorofeev_style/sections/baseline_spectral_picture_mu0.tex` so the Figure 2 caption explicitly states that the right-edge numbers denote the working branch numbering and that color/number stay attached to a branch even when its current place in the spectrum changes.
- Reworked the article mode-shape presentation layer for Figures 4 and 5: `paper_dorofeev_style/generate_targeted_analysis_figures.py` now saves silent per-panel PNGs with only the undeformed geometry and the two deformed arm curves, while `paper_dorofeev_style/manuscript.tex` and `paper_dorofeev_style/sections/mode_shapes_and_frequency_shift.tex` assemble those panels through `subfigure` blocks and move the panel labels `а)`–`е)` / `а)`–`в)` into LaTeX without changing the tracked states, branch logic, or data.

- Fixed the article mode-shape sign presentation in `paper_dorofeev_style/generate_targeted_analysis_figures.py`: within each figure, each branch family now uses a reference panel plus scalar-product sign alignment, with the reference orientation anchored by the short-arm deflection direction, so panel-to-panel sign flips no longer create artificial visual "rearrangements" of the same tracked mode.

## 2026-04-18

- Added `scripts/analyze_target_descendants_beta15_r5.py` for targeted analysis-only tracking of `bending_desc_01`, `bending_desc_02`, and `bending_desc_04` at `beta = 15 deg`, `r = 0.005`, and `0 <= mu <= 0.9`: the script reuses the established FEM continuation and MAC-based descendant tracking, counts local half-wave order on the left and right arms from reconstructed beam-shape profiles, compares each tracked state against single-rod CS reference families, and saves compact CSV/PNG outputs in `results/` without changing the determinant, shared solver logic, or baseline FEM.
- Updated `paper_dorofeev_style/sections/analytic_model.tex` and the `Lambda(mu)` figure caption to document the compact coupled-system characteristic determinant and to identify the dashed reference curves as single-rod clamped-supported families with $l_1=l(1-\mu)$ and $l_2=l(1+\mu)$, without changing the verified mathematics or determinant code.

## 2026-04-15

- Added `scripts/plot_flat_mu_bending_desc_01_mu0_0p1_0p2_ru.py` to generate the clean Russian `2x3` mode-shape figure for the `beta = 15 deg` descendant branch `bending_desc_01` at `r = 0.005, 0.02` and `mu = 0, 0.1, 0.2`, reusing the established FEM MAC-tracking plotting workflow without changing the determinant, shared solver logic, or baseline FEM.
- Added `scripts/plot_flat_mu_bending_desc_02_mu0_0p1_0p2_ru.py` to generate the clean Russian `2x3` mode-shape figure for the `beta = 15 deg` descendant branch `bending_desc_02` at `r = 0.005, 0.02` and `mu = 0, 0.1, 0.2`, reusing the established FEM MAC-tracking plotting workflow without changing the determinant, shared solver logic, or baseline FEM.

## 2026-04-13

- Added `scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py` to generate shared `2x2` fixed-`beta` mu-sweep figures in `Lambda` for `beta = 7.5 deg` and `beta = 15 deg` across radii `r = 0.005, 0.01, 0.015, 0.02`, reusing the established positive-`beta` low-branch comparison workflow together with the same presentation style and muted CS reference lines as the existing `beta = 0` four-radii figure.
- Added `scripts/plot_mu_sweep_beta0_four_radii_compare.py` for a shared `2x2` analytic-vs-FEM `beta = 0` mu-sweep figure in `Lambda` at radii `r = 0.005, 0.01, 0.015, 0.02`, reusing the existing type-aware bending matching and muted CS single-rod dashed reference-line logic and saving one common PNG together with per-radius and combined CSV tables in `results/`.
- Добавлен отдельный сценарий `scripts/compare_beta0_analytic_vs_fem.py` для сравнения determinant-based аналитических частот с baseline FEM при `beta = 0` без изменения формул, determinant или solver logic.
- Новый сценарий использует существующий `results/fem_spectrum.csv`, если он согласуется с текущим baseline FEM при `mu = 0`; в противном случае пересчитывает `beta = 0` sweep из `src/my_project/fem/python_fem.py` в памяти.
- В сценарий добавлены воспроизводимые выгрузки в `results/`: таблица сравнения при `mu = 0`, таблица классификации bending/axial мод по FEM-формам при `mu = 0`, таблицы radius-sensitivity при `mu = 0` и comparison plot для первых `beta = 0` ветвей.
- Сценарий `scripts/compare_beta0_analytic_vs_fem.py` расширен type-aware сопоставлением мод: отдельный `mu`-график для bending-ветвей, отдельные графики первой axial-ветви по `mu` и по `r`, а также CSV-таблицы с `analytic_branch_id`, `fem_mode_id`, `mode_type`, частотами и относительной ошибкой.
- Добавлен сценарий `scripts/compare_beta_positive_type_aware.py` для multibeta type-aware comparison при `beta > 0`: ветви переносятся продолжением от проверенного случая `beta = 0`, FEM-потомки сопоставляются по MAC+frequency continuity, а в `results/` сохраняются multibeta CSV-таблицы и графики bending/axial ветвей и смешения по `fem_axial_fraction`.
- Добавлен сценарий `scripts/analyze_branchwise_fem_spectrum.py` для FEM-only branch-wise анализа полного отслеживаемого спектра на сетке `(beta, mu)`: сохраняются branch-id-таблица с `current_sorted_index`, `axial_fraction`, соседними gap-ами и численными `df/dmu`, `df/dbeta`, а также summary/extrema CSV и графики branch-colored spectrum, axial-fraction heatmaps и карты минимальных gap-ов.
- Добавлен сценарий `scripts/plot_beta_sweep_mu0_compare.py` для презентационного сравнения analytic vs FEM при `mu = 0` на beta-sweep `0..90°` и радиусах `r = 0.005` и `r = 0.01`: в `results/` сохраняются отдельные PNG-графики в переменной `Λ` и CSV-таблицы с `beta`, `branch_id`, `analytic_hz`, `fem_hz`, `relative_error` и matched `fem_mode_id`.
- Сценарий `scripts/plot_beta_sweep_mu0_compare.py` обновлён для сглаженного presentation sweep без не физических high-branch spikes: в проблемных окнах по `beta` добавлены локально уплотнённые точки, analytic acquisition расширен до 16 candidate roots, локальный `scan_step` по `Λ` уменьшен, а для тесных пар используется local continuation поверх существующего determinant/root-finding слоя без изменения формул, shared solver logic или baseline FEM.
- Problem-window smoothing из `scripts/plot_beta_sweep_mu0_compare.py` расширен на радиусы `r = 0.015` и `r = 0.02`, чтобы presentation beta-sweep оставался гладким и без branch-swap artifacts и для дополнительных high-branch close-pair зон.
- Добавлен сценарий `scripts/plot_beta_sweep_mu0_four_radii_compare.py` для общей `2x2` figure analytic vs FEM при `mu = 0` и радиусах `r = 0.005, 0.01, 0.015, 0.02`; в `results/` сохраняются общий PNG, combined CSV по всем 4 радиусам и per-radius CSV для новых случаев `r = 0.015` и `r = 0.02`.

## 2026-04-07

- Проведена проверка согласованности формул между `docs/theory/equations.tex` и аналитическими программами в `src/my_project/analytic/`.
- Зафиксировано замечание по `docs/literature/pdf/2003JSVb.pdf`: determinant-like запись матрицы `T` в формуле `(39)` требует ручной сверки знаков с условиями непрерывности и не должна переноситься в код без проверки.
- Выполнен структурный рефакторинг `src/my_project/analytic/FreqFromAngle.py` и `src/my_project/analytic/FreqFromMu.py` без изменения формул, determinant, порядка неизвестных и математических соглашений.
- Общий математический слой вынесен в `src/my_project/analytic/formulas.py` и `src/my_project/analytic/solvers.py`.
- Добавлен smoke test `tests/test_analytic_smoke.py` для проверки контрольных значений после рефакторинга.
- Добавлена baseline-программа `src/my_project/analytic/FreqMuNet.py` для `mu`-графика в переменной `Lambda` с дополнительными reference-кривыми одиночного стержня CS.
- `src/my_project/analytic/FreqFromMu.py` и `src/my_project/analytic/FreqMuNet.py` переведены на общий `mu`-sweep слой в `src/my_project/analytic/solvers.py` без изменения determinant, формул и численного смысла; для `FreqMuNet.py` сохранён baseline greedy branch tracking, для `FreqFromMu.py` — текущий shared tracking mode.
- `README.md` обновлён под новый analytic workflow и запуск `FreqMuNet.py`.
- Добавлен сценарий `scripts/plot_freq_mu_vs_fem.py` для наложения аналитических ветвей `FreqFromMu` на tracked FEM-частоты из `results/fem_spectrum.csv` и сохранения comparison PNG в `results/`.
- Обновлено оформление `scripts/plot_freq_mu_vs_fem.py`: FEM-маркеры теперь рисуются цветом соответствующей аналитической моды, при этом легенда comparison-графика остаётся компактной.
