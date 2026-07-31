# CoupledBeams

CoupledBeams is a research repository for frequency models and computations for coupled beams. The repository combines analytic frequency calculations, a baseline FEM implementation of the same problem, and the local theory, literature notes, and consistency checks used to support them.

## Project Layout

- [Research index](docs/research_index.md) -- research directions, canonical documentation,
  and current scientific status.
- [Generated results index](docs/results_index.md) -- workflow map for ignored/generated results.
- [Archive policy](docs/archive_policy.md) -- preservation and soft/hard archive rules.
- [Refactoring status](docs/refactoring/README.md) -- verified inventory and staged refactoring status.
- [Script status](scripts/STATUS.md) -- preferred, active, completed, historical, and
  compatibility workflows.
- `docs/project_rules.md` -- global project rules for branch identity,
  diagnostics, thin-rod applicability, model-extension checks, and FEM
  comparison conventions.
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
- completed Chapter-2 source reproduction for one free-free and cantilever
  monoclinic rod from Yartsev; the coupled-rod model is the next research step
  and is not implemented.

See the [research index](docs/research_index.md) for canonical documents, implementations,
conclusions, and status definitions.

## Project Status

The isotropic analytic/FEM baseline is stable. Several diagnostic branches are
active, while the K=10 completeness, Step-3A, and exact Rules A/B/S stages have
completed records. The Rule-S engineering-selector path is closed after the
negative cost result `rule_S_cost_not_beneficial`; this does not refute Rule S
mathematically. The selected anisotropic Chapter-2 single-rod source line has
completed free-free and cantilever gates. Coupled anisotropic rods, a
production API, and FEM validation have not been started.

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
The next separate research gate is the derivation and verification of internal
rigid angular-joint conditions; no coupled-rod determinant is implemented.

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
