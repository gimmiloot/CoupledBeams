# CoupledBeams

CoupledBeams is a research repository for frequency models and computations for coupled beams. The repository combines analytic frequency calculations, a baseline FEM implementation of the same problem, and the local theory, literature notes, and consistency checks used to support them.

## Project Layout

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
