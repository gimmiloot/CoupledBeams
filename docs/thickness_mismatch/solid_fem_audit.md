# 3D Solid FEM Audit For Single-Rod Diagnostics

## Scope

This is a diagnostic-only preparation note for an independent 3D solid FEM
check of the one-rod fixed-fixed Euler-Bernoulli / Timoshenko comparison. It
does not promote the model to article use and does not change the baseline
determinant, old solvers, `src/my_project/analytic/formulas.py`, or the
existing FEM physical model.

The first benchmark is intentionally a single straight circular rod, not the
coupled-rod geometry.

## Local Tool Audit

The local audit was run from `D:\PHD\CoupledBeams\CoupledBeams` on
2026-05-26.

Detected:

- Python with `numpy` and `scipy`.
- Existing in-repo 1D/beam FEM utilities:
  - `src/my_project/fem/python_fem.py`
  - `scripts/analysis/fem_check_thickness_mismatch_eta_p0p5_beta15.py`
  - `scripts/analysis/plot_thickness_mismatch_fem_comparison_eta_p0p5_beta15.py`

Not detected:

- `gmsh` command-line executable on PATH.
- `gmsh` Python API.
- `meshio` Python API.
- CalculiX `ccx` on PATH.
- Code_Aster `as_run` or `run_aster` on PATH.
- Salome-Meca commands `salome` or `salome_meca` on PATH.
- Matching Gmsh / CalculiX / Code_Aster / Salome-Meca directories in the
  checked common Windows install locations.

## Chosen Workflow

The preferred lightweight path is:

1. Generate a 3D cylinder mesh with Gmsh.
2. Export a Gmsh `.msh` file and an Abaqus/CalculiX-style `.inp` mesh.
3. Run a CalculiX modal eigenfrequency analysis if `ccx` is available.
4. Convert 3D angular frequencies to the single-rod project normalization:

```text
Lambda_FEM = sqrt(Omega)*(L/2)*(rho*A/(E*I))**0.25
```

Because no mesh generator or external solid solver was detected locally, the
current run writes reproducible input templates only. It does not create fake
3D frequency results.

## Benchmark Definition

- Geometry: straight circular cylinder.
- Length: `L = 1.0`.
- Radius: `r = epsilon*L`.
- Epsilons: `0.025`, `0.05`, `0.075`; the last value is a stress-test.
- Material: `E = 1`, `rho = 1`, `nu = 0.3`.
- Boundary conditions: both end faces fully fixed, `Ux = Uy = Uz = 0`.
- Requested solid modes in the template: `16`.

Circular 3D solid rods should show bending doublets because the two transverse
bending planes are equivalent. Small numerical splitting is expected. The
first workflow report compares sorted rows cautiously and does not claim
exact 1D-to-3D mode identity.

## Diagnostic Entry Point

```text
python scripts/analysis/solid_fem_single_rod_fixed_fixed.py
```

Current outputs:

- `results/single_rod_fixed_fixed_3d_solid_fem_report.md`
- `results/single_rod_fixed_fixed_3d_solid_fem_comparison.csv`
- `results/solid_fem_single_rod/eps_*/single_rod_fixed_fixed_eps*.geo`
- `results/solid_fem_single_rod/eps_*/single_rod_fixed_fixed_eps*_ccx_modal.inp`

If Gmsh is installed later, one manual mesh command for the first case is:

```text
gmsh results/solid_fem_single_rod/eps_0p025/single_rod_fixed_fixed_eps0p025.geo -3 -format inp -o results/solid_fem_single_rod/eps_0p025/single_rod_fixed_fixed_eps0p025_mesh.inp
```

Then install or expose CalculiX `ccx` on PATH and rerun the diagnostic script.

## Extension Notes

Before extending this workflow to coupled rods:

- verify single-cylinder mesh convergence for the first bending doublets;
- inspect 3D mode shapes rather than relying only on sorted frequency order;
- decide how doublet splitting should be summarized;
- keep Gmsh / CalculiX / Code_Aster as optional external tools, not mandatory
  project dependencies.
