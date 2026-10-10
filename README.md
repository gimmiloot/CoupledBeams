# CoupledBeams

CoupledBeams is a research repository for frequency models and computations for coupled beams. The repository combines analytic frequency calculations, a baseline FEM implementation of the same problem, and the local theory, literature notes, and consistency checks used to support them.

## Seven-field spatial nonlinear verification

[The spatial report](docs/numerics/nlsp_spatial_nonlinear_3d_fem_verification.md)
uses the existing continuation CLI's `--spatial-verification` mode and a
[separate config](data/input/nlsp_spatial_nonlinear_3d_fem_verification.json).
Six p48/p64 1D cases and all four medium/fine 3D jobs reach 0.25T1.
Overall NUMERICAL_PARTIAL: full nonlinear w/v differences are 3.55%/3.03%,
but v correction/evolution differences remain 47.29%/38.80% on declared
common scales. Small rotational corrections, dynamic twist, c, velocities,
single-level temporal checks and energy remain qualified. Frozen physics,
load and historical results are preserved. Completed cache/report/plot/
postprocess replay makes zero scientific calls; no further study is started.

## Saved spatial-profile diagnostic

[The profile audit](docs/numerics/nlsp_spatial_profile_audit.md) examines the
saved u, theta, c and 3D effective-contraction profiles after FEM-3C. The existing
continuation CLI exposes `--profile-audit` with a separate
[postprocessing config](data/input/nlsp_spatial_profile_audit.json).
It reads saved nodal/coordinate histories, checks raw section recovery,
independent C3D10 strain averages and p48/p64 sensitivity, and performs no new
FEM, ODE, static, eigen or mesh calculation. The effective contraction remains
a diagnostic proxy; historical FEM-3C and spatial/energy qualifications remain
unchanged. No smoothing or phase/amplitude fitting is applied.

## FEM-3C: limited straight-rod verification

The [technical report](docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification)
and [dissertation summary](docs/numerics/nlsp_straight_rod_3d_fem_verification_summary.md)
combine the unchanged linear, static and dynamic comparisons. Explicit
`python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --validation --compute`
uses the [separate bounded config](data/input/nlsp_nonlinear_dynamic_validation.json)
and preserves the old default and `--long-horizon` workflows. Four quarter-period
medium/fine controls pass the predeclared robustness guide: the updated evolving
bending-correction difference is 7.76472% on its own scale. Full-period p64/p48
1D trajectories are complete; their all-eight spatial check remains PARTIAL
(only u,w displacement pass). Both medium 3D full-period jobs reach T1. Full nonlinear w/evolving-correction
max differences are 9.45%/10.07% on declared full-horizon scales; this comparison
is illustrative. The straight-rod verification is complete with qualifications. Native energy bookkeeping remains PARTIAL; no coefficients,
loads, initial equilibria or physics are fitted to FEM. Matching completed
cache/report/plot perform no new scientific calls. Execution stays within the
separate authorization; no automatic study follows this block.

## Project Layout

- [FEM-3B: evolution of the nonlinear dynamic correction](docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3b-nonlinear-correction-evolution-and-longer-horizon)
  -- `python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --long-horizon --compute`.
  The separately authorized stage diagnoses the old .05T1 data, compares saved-state
  p64/p48 dynamics to .5T1, fixes one longer horizon before 3D results and reuses
  the medium mesh for one linear/nonlinear pair. Initial static offsets and evolving
  corrections remain separate; no physics, load, mesh or temporal policy change.
  Both new native jobs complete .25T1 with502 frames each. The evolving-w
  model difference is7.30%; the signal exceeds observed output/recovery differences.
  Native energy and full1D spatial qualifications remain PARTIAL. One 3D mesh/time
  level does not certify nonlinear dynamic accuracy. Old workflows and their
  authorization-bound caches remain unchanged.

- [FEM-3AR: completed short preload/release continuation](docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3ar-controlled-continuation-after-native-elke-output-path-failure)
  -- `python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --run-pilot`.
  Two corrected medium 3D jobs and the saved-state p64 references reach .05T1.
  Preload/release gates pass; native energy bookkeeping remains PARTIAL.
  The small NL-minus-L response is largely the inherited static offset, not
  independent certification of nonlinear dynamic accuracy. Cached replay does
  zero new solves; historical FEM-3A failure remains unchanged.

- [FEM-3A: static preload and short free-motion pilot](docs/numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md)
  -- `python scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py --run-pilot`
  preserves the first medium native access-violation failure. Source/protocol/
  saved-state preflight pass; actual transfer and trajectories remain unverified.
  No NL/1D dynamic run or automatic retry; this blocked pilot does not establish
  dynamic convergence or authorize another/longer study.

- [FEM-2R: authorized nonlinear static continuation](docs/numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure)
  -- `python scripts/analysis/resume_nlsp_nonlinear_static_3d_fem.py --run-fem`
  reuses the preserved 1D equilibria and three saved meshes under the same frozen
  load. A separate authorization-bound ledger preserves the original failed
  FEM-2 attempt; cache/report/plot perform no new scientific calculations.
  This bounded static comparison does not authorize FEM-3.

- [FEM-2: nonlinear static comparison](docs/numerics/nlsp_nonlinear_static_3d_fem_validation.md)
  -- `python scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py --preflight`:
  frozen dead load and p48/p64 quartic equilibrium PASS. The first medium 3D
  linear attempt failed at input parsing; remaining jobs were not run. Overall
  PARTIAL, no 1D/3D nonlinear validation. Report/plot use preserved evidence;
  another solver attempt requires a separate decision after the hard-gate stop.


- [FEM-1R: one extra .020 mesh](docs/numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement)
  -- `python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py --check-source`
  and `--run-fem`: immutable parent/reference,one24-mode job,all8 preset mesh
  checks accepted. Matching cache/report/plot do zero new solver calls; no
  new1D roots,old-job repeats,fifth grid or nonlinear calculation.

- [Full-family rectangular linear3D FEM-1](docs/numerics/nlsp_linear_rectangular_3d_fem_validation.md)
  -- `python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py --preflight`
  and `--run-fem`: frozen thick diagnostic,3 audited C3D10 levels,24 modes each,
  both bending planes/MH/twist uniquely identified. Mesh convergence remains
  PARTIAL; no nonlinear validation or automatic refinement. Cached report/plot
  perform zero new solver/eigen calls.

- [Bounded nonlinear physical sanity checks](docs/theory/weakly_nonlinear_planar_physical_sanity_checks.md)
  -- `python scripts/analysis/check_weakly_nonlinear_planar_physics.py --compute`
  reuses the large one-T1 result and permits one half-amplitude p64 trajectory.
  Expected amplitude/leading-response patterns, formal stretching and reactions
  are checked with explicit finite-p/truncation caveats. Diagnostic complete with
  qualifications; report/plot/cache do zero new integration. No physical/3D
  validation claim or automatic next stage.

- [Prepared movement through one linear period](docs/theory/planar_prepared_initial_state.md#prepared-one-period-feasibility)
  -- `python scripts/analysis/prepare_planar_initial_state.py --compute --one-T1`
  reuses saved initial coordinates and completes3 exploratory runs toT1.
  Prefix/temporal/energy checks PASS; spatial comparison remains7/8 PARTIAL.
  Report/plot/cache reuse saved evidence without integration. The historical
  `--compute --feasibility` horizon remains0.1T1.

- [Prepared-state short feasibility](docs/theory/planar_prepared_initial_state.md#prepared-precision-feasibility)
  -- explicit `python scripts/analysis/prepare_planar_initial_state.py --compute --feasibility`
  completes three0.1T1 controls using the same frozen physical initial state.
  Initial representation passes; strict strong/weak remains PARTIAL. Runs are
  EXPLORATORY_NOT_CERTIFIED, temporal comparison PASS, spatial comparison PARTIAL.
  Report/plot-only and matching compute read cached evidence without new ODE.

- [Historical prepared planar initial-state audit](docs/theory/planar_prepared_initial_state.md)
  -- constant/second-harmonic profiles from saved spectra, a common O2 initial
  state and a checked cubic endpoint correction. The original zero-axial
  case remains separate. Preparation passes; the permitted nonlinear pairs
  fail initial-projection admission, so the bounded pilot is PARTIAL with
  zero new integrations. Run
  `python scripts/analysis/prepare_planar_initial_state.py --compute`;
  matching compute and report/plot-only reuse the saved cache.

- [Exact-time leading axial response diagnostic](docs/theory/planar_second_order_axial_response.md)
  -- leading u2,c2 forced by the same continuous first bending mode, with all
  Shen coefficients retained and no ODE integration. The bounded diagnostic
  is COMPLETE; spatial convergence remains PARTIAL and the old nonlinear
  pilot/recovery statuses are preserved. Run
  `python scripts/analysis/verify_planar_second_order_axial_response.py --compute`;
  `--report-only <bundle>` and `--plot-only <bundle>` reuse saved evidence.

- [Targeted planar solver diagnosis and recovery](docs/theory/weakly_nonlinear_planar_time_pilot.md)
  -- historical field-error localization, initial-compatibility audit and
  equivalent lazy energy/gradient/Hessian evaluation. Recovery remains
  PARTIAL: the guarded full p48 estimate exceeds the remaining fixed budget.
  Use `python scripts/analysis/diagnose_weakly_nonlinear_planar_rod.py --diagnose`
  for diagnosis without integration; `--report-only <bundle>` and
  `--plot-only <bundle>` read saved evidence. The previous pilot is retained.

- [Four-field planar nonlinear time pilot](docs/theory/weakly_nonlinear_planar_time_pilot.md)
  -- bounded free-motion calculation of the audited quartic action, two
  amplitudes, spatial/time controls and preserved per-case histories.
  Reproduce with `python scripts/analysis/simulate_weakly_nonlinear_planar_rod.py --compute`;
  use `--plot-only <bundle>` for redraw without integration. Read the note
  for fixed budgets, actual convergence statuses and model limitations.

- [Seven-field spatial nonlinear action audit](docs/theory/weakly_nonlinear_spatial_rod.md)
  -- isolated accepted reduced energy, exact/cubic expressions and bounded
  algebra/linear/amplitude checks. Reproduce with
  `python scripts/analysis/verify_weakly_nonlinear_spatial_rod.py --compute`;
  `--report-only <bundle>` reads saved evidence with zero roots/derivations.
  Nonlinear evolution and critical amplitudes have not been calculated.
- [Research index](docs/research_index.md) -- research directions, canonical documentation,
  and current scientific status.
- [Generated results index](docs/results_index.md) -- workflow map for ignored/generated results.
- [Longitudinal Rayleigh–Bishop literature checks](docs/theory/bishop_literature_reproduction.md)
  -- fixed circular Marais/Popov controls, source precision, reproducible CLI
  and separate numerical/print-match statuses; no angled-joint extension.
- [Single-rod Timoshenko--Bishop kinematic audit](docs/theory/timoshenko_bishop_single_rod.md)
  -- exact mixed-energy checks and an unresolved combined-kinematics gate;
  reproduce with `python scripts/analysis/audit_timoshenko_bishop_single_rod.py --compute`.
  This command performs algebra, not a combined spectrum calculation.
- [Mindlin–Herrmann + Timoshenko literature map](docs/literature/mindlin_herrmann_timoshenko_sources.md)
  -- current candidate for the combined in-plane model, with source-specific
  correction factors and reduced constitutive assumptions. The subsequent
  [single-rod source audit](docs/theory/mindlin_herrmann_timoshenko_single_rod.md)
  retains `MHTIM_VARIANT_DEPENDENT` for source variants. The separate
  [finite single-rod gate](docs/theory/mindlin_herrmann_timoshenko_single_rod.md#14-production-formulation-decision)
  selects the Jang reduced project closure and existing rectangular K=5/6;
  finite spectrum passes, hierarchy is qualified PARTIAL_PASS for contraction
  clamp interpretation. Source Jang kappa remains unstated.
  Run the finite gate with `python scripts/analysis/verify_mindlin_herrmann_timoshenko_single_rod.py --compute`.
  Reproduce with `python scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py --compute --case all --jang-kappa 5/6`;
  Here 5/6 is an explicit source-control input, not a recovered Jang value;
  the project preset belongs to the separate finite command.
  Bishop stays a standalone reference and retains its closed kinematics gate.
- [Scoped research memory](docs/memory/README.md) -- RLB-2I/RLB-2J and the EB/RLB rotational-spring theory stage, sources, and decisions.
- [M-H–Timoshenko reduced-joint beta0 gate](docs/theory/mindlin_herrmann_timoshenko_rigid_joint.md)
  -- qualified published common-DOF closure, three artificial splits reproduce
  one direct fixed--fixed rod including c/R and modal forms. Run
  `python scripts/analysis/verify_mindlin_herrmann_timoshenko_beta0_joint.py --compute`.
  No nonzero-angle spectrum or direct 3D joint derivation is claimed.
- [General-angle M-H–Timoshenko joint gate](docs/theory/mindlin_herrmann_timoshenko_rigid_joint.md#9-general-angle-geometry-действующий-project-contract)
  -- common project geometry/duality assembly, frozen beta0 and symmetry gates,
  then one fixed5/45/90 pilot with certified12+guard13. Run
  `python scripts/analysis/verify_mindlin_herrmann_timoshenko_general_beta_joint.py --compute`.
  Reduced contraction closure is retained; no 3D joint proof, angle map,
  applicability/hierarchy study or mode tracking.
- [Archive policy](docs/archive_policy.md) -- preservation and soft/hard archive rules.
- [Production M-H/Timoshenko Lambda(beta) checks](docs/theory/mindlin_herrmann_timoshenko_lambda_beta_large_checks.md)
  -- length asymmetry and section contrast, one fixed canonical Lambda
  reference,37 angles and certified sorted12+guard13. Run
  `python scripts/analysis/verify_mindlin_herrmann_timoshenko_lambda_beta_large_checks.py --compute`.
  Two PDF/PNG figures; no mode tracking or new physical model.
- [Coupled longitudinal-theory screening](docs/theory/coupled_longitudinal_theory_hierarchy_screening.md)
  -- one G20 control, eight fixed angles, elementary/planar Love/M-H with
  identical Timoshenko arms; independently sorted12+guard13 and fixed-beta
  geometric overlaps. Run
  `python scripts/analysis/screen_coupled_longitudinal_theory_hierarchy.py --compute`.
  No energy classification, mode continuation or applicability threshold.
- [Refactoring status](docs/refactoring/README.md) -- verified inventory and staged refactoring status.
- [Bounded thickness screening](docs/theory/coupled_longitudinal_theory_thickness_screening.md)
  -- five h/h0 values1–2, beta0/45/90 and three unchanged axial theories;
  coefficient scaling, certified12+guard13, direct beta0 profiles and
  fixed-case geometric overlaps. Run
  `python scripts/analysis/screen_coupled_longitudinal_theory_thickness.py --compute`.
  No fitted law, applicability threshold or across-thickness tracking.
- [Script status](scripts/STATUS.md) -- preferred, active, completed, historical, and
  compatibility workflows.
- `docs/project_rules.md` -- global project rules for branch identity,
  diagnostics, thin-rod applicability, model-extension checks, and FEM
  comparison conventions.
- [`docs/numerics/`](docs/numerics/README.md) -- project-wide numerical
  workflow policies that are independent of any one physical model.
- [Frequency-map computation policy](docs/numerics/frequency_map_computation_policy.md)
  -- the canonical `frequency-map-v1` contract for ordinary maps, strict
  audits, and rendering from saved data.
- `docs/writing/` -- diagnostic-to-article workflow notes.

- `docs/theory/` — verified local theory, equations, assumptions, and theory notes.
- `docs/literature/` — literature PDFs, source notes, and bibliography material.
- `src/my_project/analytic/` — analytic Python programs for the coupled-beam frequency problem.
- `src/my_project/fem/` — baseline FEM implementation.
- `tests/` — smoke tests and local verification helpers.
- `results/` — generated and ignored computational outputs and tables; see
  `docs/results_index.md` for the workflow map.

Article workspaces referenced by historical documentation are not part of the
tracked public checkout. The thickness-mismatch / Timoshenko article remains a
planned or local workflow whose authoritative workspace status requires manual
review; no absent `paper_*` directory is assumed to exist here.

## Research Directions

- baseline isotropic in-plane coupled beams and their analytic/FEM comparison;
- descendant branch tracking, veering, modal exchange, and localization;
- mass-preserving thickness mismatch;
- Euler--Bernoulli versus Timoshenko applicability and validation diagnostics;
- the completed and closed EB safe-prefix engineering study;
- out-of-plane Euler--Bernoulli bending plus Saint-Venant torsion;
- completed Chapter-2 single-rod source gates and a passing small elastic
  rigid-joint pilot for two monoclinic rods; no final coupled model exists.

See the [research index](docs/research_index.md) for canonical documents, implementations,
conclusions, and status definitions.

## Project Status

The isotropic analytic/FEM baseline is stable. Several diagnostic branches are
active, while the K=10 completeness, Step-3A, and exact Rules A/B/S stages have
completed records. The Rule-S engineering-selector path is closed after the
negative cost result `rule_S_cost_not_beneficial`; this does not refute Rule S
mathematically. The selected anisotropic Chapter-2 single-rod source line has
completed free-free and cantilever gates, and its first ideal rigid-joint
elastic pilot passed. A finite rectangular orthotropic endpoint 1D-FEM gate
has also been completed and remains `PARTIAL_PASS` after targeted
equal-element-length refinement; a final coupled model, production API,
unequal-thickness study, and 3D FEM validation have not been started.

## Analytic Layer

- `src/my_project/analytic/FreqFromAngle.py` — analytic scenario sweeping the coupling angle `beta`.
- `src/my_project/analytic/FreqFromMu.py` — analytic scenario sweeping the length-asymmetry parameter `mu` in frequency units, with tracked branches and optional close-pair diagnostics.
- `src/my_project/analytic/FreqMuNet.py` — baseline fixed-`beta` `mu`-sweep plot in dimensionless `Lambda`, with additional single-beam CS reference curves over the coupled-beam branches and CLI controls for `--beta`, `--epsilon`, `--num-modes`, `--num-dashed-lines`, and output path.
- `src/my_project/analytic/formulas.py` — shared matrix and determinant assembly extracted during refactoring.
- `src/my_project/analytic/solvers.py` — shared numerical solver logic extracted during refactoring.
- `scripts/README.md` — script guide with the main commands, analysis/audit scripts, internal helpers, legacy wrappers, outputs, and usage notes.
- [Script status](scripts/STATUS.md) — concise workflow-status and preferred-entry-point map.

The analytic refactoring did not change the formulas, determinant structure, unknown ordering, signs, or coefficients. It only extracted the common layer for reuse. `FreqFromMu.py` and `FreqMuNet.py` now share the same common mathematical layer and differ only in plotting/output behavior and in their preserved branch-tracking mode.

Run from the repository root:

```bash
python scripts/run/run_beta_sweep_mu0_four_radii.py
python scripts/run/run_mu_sweep_beta0_four_radii.py
python scripts/run/run_mu_sweep_fixed_beta_four_radii.py
python scripts/run/run_mu_sweep_four_betas_analytic.py --betas 15 30 45 60
python scripts/run/run_tracked_bending_descendant_shape_ru.py
python scripts/run/run_branchwise_fem_audit.py
```

For the single tracked descendant shape runner, ordinary runs use the editable `USER PARAMETERS` block at the top of `scripts/run/run_tracked_bending_descendant_shape_ru.py`; CLI arguments remain available as overrides.

See `scripts/README.md` for the full script inventory and legacy command map.

The diagnostic Chapter-2 single-rod reproduction is run with:

```bash
python scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py
```

It writes Git-ignored evidence under
`results/anisotropic_rods/yartsev_ch2_free_free/` and does not alter the
isotropic analytic/FEM baseline.

The completed cantilever source reproduction and its saved-data-only boundary
decision are documented in
[`docs/anisotropic_rods/yartsev_ch2_cantilever_reproduction.md`](docs/anisotropic_rods/yartsev_ch2_cantilever_reproduction.md).
The separate internal rigid angular-joint derivation and small elastic pilot
are complete with `PASS`; this is not a stable or final coupled-rod baseline.

## FEM Baseline

- Baseline file: `src/my_project/fem/python_fem.py`
- Dependencies: `numpy`, `scipy`
- Input files: none
- Output CSV: `results/fem_spectrum.csv`

Run from the repository root:

```bash
python src/my_project/fem/python_fem.py
```

## Theory And References

Base notation in the theory-facing materials is oriented to `docs/literature/pdf/Статья-Дорофеев-2025.pdf`.

When comparing against `docs/literature/pdf/2003JSVb.pdf`, account for the known sign issue in its determinant-like matrix record. The printed sign pattern from that source must not be copied blindly. For the current local implementation, the verified local theory and the corresponding local code are treated as the source of truth.

## Tests

The analytic smoke test is `tests/test_analytic_smoke.py`.

Run from the repository root:

```bash
python -m unittest discover -s tests -p "test_analytic_smoke.py"
```
