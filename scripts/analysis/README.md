# Analysis scripts

Most historical analysis and audit scripts are currently kept at `scripts/*.py` because several of them import each other through `scripts.*`.
New standalone diagnostics can live here when they do not require broad import churn.

Current entry points:

- `postprocess_rl_full_rr_continuation.py`: diagnostic continuation from saved
  square RL/Timoshenko CSV/JSON/NPZ. Reuses the assignment helper with the full
  inertia product on a common material coordinate, midpoint/end-point metric
  checks and pair-subspace diagnostics. Never calls spectral or basis routines;
  unresolved descendants remain explicitly unconfirmed.

  Run `python scripts/analysis/postprocess_rl_full_rr_continuation.py --config <local.json> --output results/<new-directory>`.
  The local configuration selects existing maps and their manifests, grids,
  tracked counts and diagnostic thresholds. Use `--redraw-from <completed-output>
  --output results/<new-figures>` to reuse the derived tables without repeating
  assignment. Completed results are never overwritten. This separate saved-field
  input contract supports changing lengths without changing the symmetric-only
  publication renderer or the historical aspect-four audit. Focused tests:
  `tests/test_rl_full_rr_continuation.py`.

  The `p26_figure5` presentation profile accepts fixed-section beta paths at
  nonzero tau with K6/guard7, as well as legacy K8/guard9 data. Canonical section
  coefficients and all saved inertia components are checked; legacy mu paths
  remain restricted to tau=0. It reuses an accepted symmetric figure manifest,
  exports raw pair observations and principal angles, and draws unresolved
  frequencies as neutral markers without assigning branch colors. Branch IDs
  start at beta=0 separately for each asymmetric construction.
  Local preset: `private/dissertation/research/chapter1/chapter1_asymmetric_beta_publication_postprocess.json`.
  Use an explicit new `--output` directory; the delivery command and qualified
  results are in the adjacent private `.md` report. `--redraw-from` preserves
  derived assignments and never solves spectra. No solver fallback is allowed.
  The refined local preset is `chapter1_asymmetric_beta_publication_refined.json`
  in the same private directory. Its `refinement_window` compares retained-node
  steps and performs a separate reverse pass of the selected beta series, using
  only saved fields. Mode panels select confirmed descendants and use one dark
  solid color for deformed axes; unresolved fields receive only sorted labels.

  Presentation-only options for this profile: `render_figures` limits emitted
  PDF/PNG files to the listed figure names; `mode_label_style: "branch_only"`
  omits the sorted index for confirmed mode panels (unresolved modes keep their
  raw-position labels); `mode_label_fontsize` sets the mode-axis label size.
  Defaults preserve the diagnostic presentation. The redraw manifest records
  the actual configuration and renderer hash. The local manuscript preset is
  `chapter1_section132_presentation.json`, in the same private directory;
  use it with `--redraw-from` and a new output directory, without new tracking.

  A presentation may combine confirmed series from completed bundles through
  `continuation_sources` (series ID to bundle directory). This redraw mode
  checks source manifests and hashes, physical fields, frequencies and the
  selected six assignments at every node; it retains the original guard forms.
  `frequency_panel_series` orders the panels, `frequency_ylim` sets their common
  ordinate range, and `scene_angles` selects existing shape nodes. Empty angle
  lists suppress shape sheets. `mode_label_style: "ordinal_rows"` labels rows
  by mode number only when the confirmed descendant matches that sorted
  position; `mode_xlim` and `mode_ylim` control the shared viewing window.
  Unrequested component figures are not constructed. Existing presets retain
  their defaults. Example with three asymmetric cases and no new assignments:
  `results/rl_full_rr_maps/_section132_illustrations_20261110/config.json`;
  its adjacent `review.md` records the exact command and reproduction check.

- `audit_rl_full_rr_reflection_intersection.py`: saved-class audit and local
  intersection diagnostic for the accepted symmetric square RL/RR model.
  Uses the previously derived 3x3 restrictions, unchanged basis and full mass;
  verifies first roots with finite inertia queries using saved pole catalogues.
  It never runs a full prefix, builds catalogues, or invokes Ritz.

  Run `python scripts/analysis/audit_rl_full_rr_reflection_intersection.py`
  with `--stage a` for the read-only archive check, `--stage all --output <new>`
  for the two-stage diagnostic, or `--stage plot` for saved-data figures.
  Completed outputs are not overwritten. CSV, NPZ, checks and PNG/PDF remain
  in ignored `results/`; the audit report and provenance accompany the local run.

- `classify_rl_full_rr_reflection_modes.py`: read-only reflection classification
  of the symmetric square RL/RR archive. Reuses the existing reflection,
  full-mass integrators and physical trace/state operators; accepts saved
  fields only and blocks basis evaluation, spectral counts and root searches.
  A separate field/projection output contract avoids changing the historical
  aspect-four postprocessor or source-hashed spectrum workflow.

  Run `python scripts/analysis/classify_rl_full_rr_reflection_modes.py` with
  the accepted local defaults, or specify `--source`, `--inventory`, `--output`.
  Existing outputs are preserved; choose a fresh output directory to repeat.
  `--diagnostics-only` first writes distributions without assigning classes.
  Outputs are classification/class-prefix CSV, physical and projector checks,
  sampled sign-change intervals, PNG/PDF and a SHA256 manifest. A class rank
  is not a global tracked-branch identity. Results remain in ignored `results/`.

- `compare_analytic_fem_tracked_descendant_shape.py`: analytic determinant-nullspace vs FEM local component comparison for one tracked bending descendant.
- `plot_rl_full_rr_publication.py`: read-only publication presentation of saved
  RL/RR sorted BETA maps (six or eight positions from the configuration,
  local first-pair panels and overview; the guard stays hidden).
  Uses `map_rl_full_rr.saved_table` for CSV checksum verification; refuses
  unresolved bundles. Does not evaluate the spectrum. The shared geometry
  presenter also supplies line styles to diagnostic rendering.
  A separate entry point preserves the source-hashed solver/cache
  contract; no script is needed for each case or panel.

  Run `python scripts/analysis/plot_rl_full_rr_publication.py --bundle <bundle>
  --config <config>` (one line). The configuration must match the saved maps and
  include `grid_provenance.analysis_step`, `local_step` and `local_windows`
  keyed by map name. Windows only select saved nodes. `--dpi` and `--font-size`
  change presentation without solving. Requires Matplotlib and Times New Roman;
  output PNG/PDF, sampled gaps, gallery and hashed manifest are written to a
  separate fingerprinted directory under `results/rl_full_rr_maps/figures/`.

See `scripts/README.md` for the full logical inventory and recommended commands.

- `plot_rl_full_rr_reflection_publication.py`: publication figures from accepted
  reflection-class CSV and a saved crossing-field NPZ. Run
  `python scripts/analysis/plot_rl_full_rr_reflection_publication.py` with the
  accepted local defaults, or supply `--classification` and `--intersection`.
  Both directories must contain complete manifests with matching data hashes.
  `--dpi`, `--font-size`, `--displacement-fraction` affect presentation only.
  Outputs are two-panel frequency and mode PDF/PNG, plotted-frequency and axis-coordinate CSV,
  normalization checks and a manifest in a fresh fingerprinted directory under
  `results/rl_full_rr_maps/figures/`; existing outputs are never overwritten.

  This is a distinct input contract from `plot_rl_full_rr_publication.py`
  (sorted frequency maps): it reads class ranks and a separate verified crossing
  solution, preserving their provenance without changing the map workflow.
  It reuses the saved quadrature weights and physical-coordinate convention
  documented by the intersection diagnostic; standard NumPy/Matplotlib handle
  the arrays and rendering. The branch view imports the existing tracker's
  assignment helpers; none of its spectral or basis routines are called.
  Unit-inertia mode amplitudes are retained;
  one common multiplier and one sign per mode control their display. Styling
  follows the published frequency panels and black mode-axis illustrations.

  `--frequency-view sorted` selects a single panel of positions 1--8 from
  `--spectrum-bundle` / `--map-name` (defaults: the accepted refined symmetric
  square map). It checks all nine CSV positions against each hashed point JSON,
  the saved guard-count certificate and NPZ fields. The two already certified
  coincident roots are added only to positions 1 and 2; no other frequencies
  are generated. Class CSV is used only to verify the descriptions of the two
  mode panels. The default `--frequency-view classes` retains the earlier view.

  `--frequency-view branches` selects the eight lowest modes at beta=mu=0,
  then continues their saved shapes independently inside each reflection class
  with `analytic_branch_tracking.assignment_diagnostics` and
  `build_next_points_from_assignment`. The adapter embeds the full RL/Timoshenko
  inertia product in quadrature-weighted features (global translations, section
  rotation and axial derivative); it does not use the baseline EB mass metric.
  Assignment is checked against alternative candidates and frequency penalties.
  Ambiguous/low-overlap assignments stop the render. Class ranks and full sorted
  positions remain provenance fields, not branch IDs. The single-panel plot keeps
  one color per base branch, including at the saved first-pair crossing. The CSV
  records every branch's source shape and overlap; checks report coverage and
  sensitivity. All spectral values and local roots are read from existing archives.

## Geometric line styles

Since the author's decision of 2026-11-08, all RL/RR map, publication,
reflection and saved-continuation frequency renderers use
`scripts/lib/frequency_geometry_style.py`. For square sections h is the side:
solid means both L_i/h_i >= 10 (including equality), dashed means at least one
ratio is smaller. This is a preliminary geometric indicator, not a rigorous
boundary of theory validity or a numerical quality flag. Branch colors and
the existing sorted/tracked spectrum semantics are preserved.

The shared API supports beta, mu, tau and epsilon paths, increasing or decreasing,
and multiple panels. Use `FrequencyGeometry` once per path and
`plot_frequency_curve` for each main curve; do not assign dashes by rank or class.
For auxiliary CS/CC curves use `plot_partial_frequency` (neutral dotted lines
with explicit labels). Canonical `Case.lengths` and `arm.section.thickness`
provide the geometry. Geometric bisection never calls a spectral routine.
Any inserted boundary vertex is solely a linear display point on a saved
segment; original observations and their markers remain unchanged.

Each new presentation writes `geometry_nodes.csv` and `geometry_diagnostics.json`
and embeds the latter's node diagnostics in its manifest. They record both
individual ratios, their minimum, pass/fail, threshold 10, solid/dashed status,
geometric boundary positions and source hashes. Mixed-style figures receive
the author's explanatory caption; entirely solid figures do not. New figure
identities include geometry presenter hashes. Existing outputs are preserved.
Continuation redraw retains scientific-input checks and records both the old
manifest and the current presenter; it does not repeat assignments.

The common saved-map renderer accepts an EPSILON axis for presentation;
this does not add an epsilon computation mode to the spectral runner.
No archived dissertation figure is regenerated merely by updating this code.
Verification: `python -m pytest -q tests/test_frequency_geometry_style.py`;
the tests forbid spectral calls and use synthetic frequency values.

## Square prefix pilot

The existing `map_rl_full_rr.py` accepts exactly `K_plot/K_guard=6/7` or
`8/9`; omitted fields retain 8/9. Material epsilon may be .005 or .01.
All other material constraints and numerical tolerances remain unchanged.
No K10 mode or new solver was added. The canonical model recalculates all
section, inertia, stiffness, frequency-scale and regularity data for epsilon.

Verified higher prefixes supply lower ones only when the complete physical
and numerical contracts agree, including reference time and normalization.
New derived JSON/NPZ records retain the selected full modes, their mass checks
and the source record/NPZ hashes. The source archive is never truncated.
If a cluster straddles the guard, its completion is preserved and labelled;
extra ranks never become plotted branches. Source migrations require known,
hash-verified snapshots and the reviewed interface AST fingerprint; all
physical dependency hashes remain exact.

For the authorized pilot, `control_reference` and `control_reference_sha256`
select an existing PASS square control bundle. This prevents an implicit new
Ritz control run. Those archived independent controls cover epsilon=.005;
epsilon=.01 has the separately recorded coefficient and bounded pilot checks.
Omitting this option retains the older control-gate workflow, which can run
Ritz and must not be used where such work is not authorized.

The local configs under `results/rl_full_rr_maps/_pilot_20261109/` contain
six points in `pilot_eps0.005.json` and two in `pilot_eps0.01.json`.
Replay each with `python scripts/analysis/map_rl_full_rr.py --config <path>
--no-plots` (one line). Completed pilot bundles are reused with zero solves.
The runner's internal clamped-arm pole catalogues remain part of its energetic
count; these are distinct from standalone publication CS/CC spectrum studies.

The pilot is finished: 8/8 PASS, seven complete modes per point; one additional
bounded recomputation of an existing geometry passed against its archived K8.
Replay, CSV/JSON/NPZ hashes, saved-field quality and temporary plot-only checks
passed. The 1694-point main set and separate CS/CC studies were not started;
1686 planned new geometries remain after the eight pilot points.
Numerical PASS does not establish accuracy relative to three-dimensional elasticity.
Tests: `tests/test_rl_full_rr_pilot_interface.py` and the existing prefix/map suites.
