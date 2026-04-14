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
- `src/my_project/analytic/FreqMuNet.py` — baseline `mu`-sweep plot in dimensionless `Lambda`, with additional single-beam CS reference curves over the coupled-beam branches.
- `src/my_project/analytic/formulas.py` — shared matrix and determinant assembly extracted during refactoring.
- `src/my_project/analytic/solvers.py` — shared numerical solver logic extracted during refactoring.
- `scripts/plot_freq_mu_vs_fem.py` — comparison plot overlaying analytic `mu`-branches with FEM frequencies from `results/fem_spectrum.csv`.
- `scripts/compare_beta0_analytic_vs_fem.py` — focused `beta=0` audit with explicit `bending`/`axial` mode labeling, type-aware analytic/FEM matching, separate bending and axial plots, and saved comparison tables without changing the baseline formulas or FEM model.
- `scripts/compare_beta_positive_type_aware.py` — extends the type-aware audit to several `beta > 0` values using continuation from the validated `beta=0` branches, saves multibeta comparison tables, and plots bending descendants, the first axial descendant, and FEM axial-fraction mixing diagnostics versus `mu`.
- `scripts/analyze_branchwise_fem_spectrum.py` — FEM branch-continuation audit over a 2D `(beta, mu)` grid, seeded from the pure `beta=0` spectrum and saved as branch-wise tables, extrema summaries, gap maps, and branch-colored spectrum plots for axial-descendant and bending-descendant evolution.
- `scripts/plot_beta_sweep_mu0_compare.py` — presentation-oriented comparison of analytic and FEM beta-sweeps at `mu = 0` for selected radii, plotted in the dimensionless frequency parameter `Λ` with solid analytic branches, same-color FEM markers, and saved smoothed PNG/CSV outputs in `results/`; the script uses local beta-grid densification plus finer analytic root acquisition and local continuation in diagnosed high-branch spike windows without changing the shared formulas or FEM baseline.
- `scripts/plot_beta_sweep_mu0_four_radii_compare.py` — builds one shared `2x2` presentation figure for `mu = 0` and radii `r = 0.005, 0.01, 0.015, 0.02`, reusing the smoothed beta-sweep comparison style, saving one common PNG plus a combined CSV and per-radius CSV outputs for newly computed radii.
- `scripts/plot_mu_sweep_beta0_four_radii_compare.py` — builds one shared `2x2` presentation figure for the `beta = 0` mu-sweep at radii `r = 0.005, 0.01, 0.015, 0.02`, reusing the type-aware bending matching from `scripts/compare_beta0_analytic_vs_fem.py` and the muted CS single-rod dashed reference curves from `src/my_project/analytic/FreqMuNet.py`, and saving one common PNG plus per-radius and combined CSV outputs in `results/`.
- `scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py` — builds shared `2x2` presentation figures for fixed `beta` mu-sweeps at `beta = 7.5 deg` and `beta = 15 deg` across radii `r = 0.005, 0.01, 0.015, 0.02`, keeping the same four-radii style as the `beta = 0` figure while reusing the existing positive-`beta` low-branch comparison logic and the muted CS single-rod dashed reference curves.

The analytic refactoring did not change the formulas, determinant structure, unknown ordering, signs, or coefficients. It only extracted the common layer for reuse. `FreqFromMu.py` and `FreqMuNet.py` now share the same common mathematical layer and differ only in plotting/output behavior and in their preserved branch-tracking mode.

Run from the repository root:

```bash
python src/my_project/analytic/FreqFromAngle.py
python src/my_project/analytic/FreqFromMu.py
python src/my_project/analytic/FreqMuNet.py
python scripts/plot_freq_mu_vs_fem.py
python scripts/compare_beta0_analytic_vs_fem.py
python scripts/compare_beta_positive_type_aware.py
python scripts/analyze_branchwise_fem_spectrum.py
python scripts/plot_beta_sweep_mu0_compare.py
python scripts/plot_beta_sweep_mu0_four_radii_compare.py
python scripts/plot_mu_sweep_beta0_four_radii_compare.py
python scripts/plot_mu_sweep_beta_fixed_four_radii_compare.py
```

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
