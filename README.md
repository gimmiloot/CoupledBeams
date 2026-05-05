# CoupledBeams

CoupledBeams is a research repository for frequency models and computations for coupled beams. The repository combines analytic frequency calculations, a baseline FEM implementation of the same problem, and the local theory, literature notes, and consistency checks used to support them.

## Project Layout

- `docs/theory/` — verified local theory, equations, assumptions, and theory notes.
- `docs/literature/` — literature PDFs, source notes, and bibliography material.
- `src/my_project/analytic/` — analytic Python programs for the coupled-beam frequency problem.
- `src/my_project/fem/` — baseline FEM implementation.
- `tests/` — smoke tests and local verification helpers.
- `results/` — generated computational outputs and tables.

## Analytic Layer

- `src/my_project/analytic/FreqFromAngle.py` — analytic scenario sweeping the coupling angle `beta`.
- `src/my_project/analytic/FreqFromMu.py` — analytic scenario sweeping the length-asymmetry parameter `mu` in frequency units, with tracked branches and optional close-pair diagnostics.
- `src/my_project/analytic/FreqMuNet.py` — baseline fixed-`beta` `mu`-sweep plot in dimensionless `Lambda`, with additional single-beam CS reference curves over the coupled-beam branches and CLI controls for `--beta`, `--epsilon`, `--num-modes`, `--num-dashed-lines`, and output path.
- `src/my_project/analytic/formulas.py` — shared matrix and determinant assembly extracted during refactoring.
- `src/my_project/analytic/solvers.py` — shared numerical solver logic extracted during refactoring.
- `scripts/README.md` — script guide with the main commands, analysis/audit scripts, internal helpers, legacy wrappers, outputs, and usage notes.

The analytic refactoring did not change the formulas, determinant structure, unknown ordering, signs, or coefficients. It only extracted the common layer for reuse. `FreqFromMu.py` and `FreqMuNet.py` now share the same common mathematical layer and differ only in plotting/output behavior and in their preserved branch-tracking mode.

Run from the repository root:

```bash
python scripts/run/run_beta_sweep_mu0_four_radii.py
python scripts/run/run_mu_sweep_beta0_four_radii.py
python scripts/run/run_mu_sweep_fixed_beta_four_radii.py
python scripts/run/run_mu_sweep_four_betas_analytic.py --betas 15 30 45 60
python scripts/run/run_tracked_bending_descendant_shape_ru.py --branch-number 4 --mu 0.2 --epsilon 0.0025 --beta 15
python scripts/run/run_branchwise_fem_audit.py
```

See `scripts/README.md` for the full script inventory and legacy command map.

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
