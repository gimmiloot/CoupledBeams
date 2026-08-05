# Limited Chapter-2 3D FEM anchor — 3D-A0 design and readiness

## Status

```text
Unequal-thickness 1D validation phase:
COMPLETE

Overall unequal-thickness 1D validation:
PARTIAL_PASS

Limited 3D FEM anchor phase:
DESIGN_ONLY

3D-A0 design and readiness audit:
BLOCKED_MATERIAL_AND_SOLVER

3D mesh generation:
NOT_STARTED

3D eigenfrequency calculations:
NOT_STARTED
```

The 1D statuses are regressions, not inputs to the 3D readiness decision:
UT-0 is `PASS`, UT-1 is `PARTIAL_PASS`, UT-1a is
`INCONCLUSIVE_NUMERICAL_AUDIT`, UT-2 is `PASS`, and UT-3 is `PASS`.
3D-A0 changes none of those results.

## 1. Scope

3D-A0 is a design and local-readiness audit for a future, limited 3D solid
anchor. It freezes material, geometry, boundary, joint, parity, case, mesh,
and acceptance contracts before any scientific 3D solve. It is not a 3D
implementation and provides no frequency result.

The future anchor is intended to check low out-of-plane bending--torsion
frequencies of rectangular orthotropic arms with an ideal rigid point joint.
It would independently include 3D elasticity, local material orientation,
transverse shear without an assigned shear factor, rotary inertia,
Saint-Venant warping, rigid end-section compatibility, and thickness mismatch.
It will not validate a physical adhesive or mechanical joint zone, finite joint
compliance, local joint stresses, failure, damping, off-axis material, or a
parameter sweep.

The dedicated CLI is:

```powershell
python scripts/analysis/anisotropic_rods/audit_yartsev_ch2_limited_3d_fem_readiness.py
```

It accepts `--output-dir` and an optional provenance-bearing
`--material-json`. It does not extend the large unequal-thickness CLI and does
not add a production API.

## 2. Full 3D material requirement

A fully orthotropic solid requires nine independent engineering constants and
density:

\[
E_1,E_2,E_3,\quad G_{12},G_{13},G_{23},\quad
\nu_{12},\nu_{13},\nu_{23},\quad \rho.
\]

The current tracked HMS/DX-209 `BookMaterial` was inspected directly. Its
available real constants are:

| Constant | Value | Source |
| --- | ---: | --- |
| \(E_1\) | `191 GPa` | Chapter 1, Table 1.1, printed p.45 |
| \(E_2\) | `5 GPa` | Chapter 1, Table 1.1, printed p.45 |
| \(G_{12}\) | `3 GPa` | Chapter 1, Table 1.1, printed p.45 |
| \(G_{13}\) | `3 GPa` | Chapter 1, Table 1.1, printed p.45 |
| \(G_{23}\) | `2.5 GPa` | Chapter 1, Table 1.1, printed p.45 |
| \(\nu_{12}\) | `0.279` | Chapter 1, Table 1.1, printed p.45 |
| \(\rho\) | `1580 kg/m^3` | Chapter 1, Table 1.1, printed p.45 |

The dataclass has no fields for `E3`, `nu13`, or `nu23`; the HMS factory does
not set them elsewhere. The following required data therefore remain absent:

```text
E3
nu13
nu23
```

They are represented as missing values, never inferred values.

## 3. Local source audit

The local scan `docs/literature/pdf/Глава 1_compressed.pdf` was inspected only
on the relevant printed pages 40--46. It is a scanned fragment without a usable
text layer, so those pages were reviewed visually; no whole-book OCR was run.

Printed p.40, equation (1.50), states the reduced compliance convention

\[
S_{12}=-\frac{\nu_{12}}{E_1}
=-\frac{\nu_{21}}{E_2}.
\]

Printed pp.40--43 develop the reduced monoclinic-layer transformations used by
the current Chapter-2 line. Printed p.45, Table 1.1, lists the seven tracked
HMS/DX-209 quantities above, along with mechanical-loss coefficients and ply
thickness. It does not list `E3`, `nu13`, or `nu23`, and pp.40--46 do not state
transverse isotropy or equivalent relations that determine them.

Consequently the audit does not adopt

\[
E_3=E_2,\qquad \nu_{13}=\nu_{12},\qquad
\nu_{23}=\frac{E_2}{2G_{23}}-1.
\]

No typical internet value, other composite, memory-based value, or
positive-definiteness fit was used.

## 4. Compliance convention and tensor checks

For a future source-confirmed complete input, the CLI uses engineering strain
order

```text
[epsilon_11, epsilon_22, epsilon_33, gamma_23, gamma_13, gamma_12]
```

and stress order

```text
[sigma_11, sigma_22, sigma_33, tau_23, tau_13, tau_12].
```

The normal block is constructed from

\[
S_{11}=1/E_1,\quad S_{22}=1/E_2,\quad S_{33}=1/E_3,
\]

\[
S_{12}=-\nu_{12}/E_1,\quad
S_{13}=-\nu_{13}/E_1,\quad
S_{23}=-\nu_{23}/E_2,
\]

and the shear diagonal is

\[
S_{44}=1/G_{23},\quad S_{55}=1/G_{13},\quad S_{66}=1/G_{12}.
\]

The minor Poisson ratios are derived only for the reciprocity audit:

\[
\nu_{21}=\nu_{12}E_2/E_1,\quad
\nu_{31}=\nu_{13}E_3/E_1,\quad
\nu_{32}=\nu_{23}E_3/E_2.
\]

The CLI checks finiteness, positive `E_i`, `G_ij`, and density, compliance
symmetry, all three reciprocal relations, positive definiteness of compliance,
positive definiteness of `C=S^{-1}`, minimum eigenvalues, condition numbers,
and exact `theta=0` recovery of `E1`, `G12`, and `G13`. The numerical
symmetry/reciprocity/recovery tolerance is `1e-12`.

These checks were unit-tested with complete positive-definite and deliberately
non-positive-definite synthetic tensors. They were not applied to an alleged
full HMS/DX-209 tensor, because no such source-supported tensor was found.
Nothing was added to production `BookMaterial`.

An optional JSON can be material-ready only if all constants are present and
the following provenance fields are nonempty:

```text
material_name
units
source_file
source_printed_page
source_table_or_equation
poisson_convention
source_audit_status = SOURCE_CONFIRMED
```

An incomplete, ambiguous, or unconfirmed JSON remains blocked.

## 5. Solver-stack readiness

The audit used passive executable lookup for `abaqus`, `abq*`, `ansys`,
`mapdl`, `ccx`, `gmsh`, `ElmerSolver`, `salome`, `FreeCADCmd`, and `comsol`.
It used `importlib.util.find_spec`, without importing the package, for `gmsh`,
`meshio`, `pyvista`, `sfepy`, `dolfinx`, `fenics`, `petsc4py`, `slepc4py`, and
`ansys.mapdl.core`.

No listed executable or Python package was detected in the audited environment.
No version process, GUI, solver job, license checkout, download, or install was
attempted. The tracked `src/my_project/fem/python_fem.py` was found, but it is a
planar 1D axial/Euler--Bernoulli frame/beam baseline with `[u,v,theta]` nodal
DOFs. It is explicitly not a 3D solid solver and fails the 3D capability gate.

A suitable future solver must confirm all of:

1. headless scripted linear 3D solid elasticity;
2. generalized eigenfrequency extraction;
3. a full orthotropic tensor;
4. per-arm local material orientations;
5. quadratic solid elements;
6. a mass matrix;
7. exact six-DOF reference-point face coupling or equivalent MPC equations;
8. nodal eigenvector and eigenvalue export;
9. node/element/result metadata export and a reproducible parser;
10. sufficient local license/runtime readiness.

Because no candidate is present, all candidate-specific capabilities remain
`NOT_AUDITED_NOT_DETECTED`. No primary solver and no fallback solver are
selected. This is `BLOCKED_NO_CAPABLE_3D_SOLVER` as a readiness dimension.

## 6. Geometry and material axes

Each coupled case uses two independent complete rectangular prisms. For arm
`i`, local coordinates are

\[
x_i:\ \text{outer clamp to joint},\qquad
y_i\parallel\mathbf n_i,\qquad z_i\parallel\mathbf e_z,
\]

\[
y_i\in[-b_i/2,b_i/2],\qquad z_i\in[-a_i/2,a_i/2].
\]

Material axes are frozen as

\[
\mathbf e_1=\mathbf t_i,\qquad
\mathbf e_2=\mathbf n_i,\qquad
\mathbf e_3=\mathbf e_z.
\]

The project convention is

\[
\mathbf t_1=(1,0,0),\qquad
\mathbf n_1=\mathbf e_z\times\mathbf t_1,
\]

\[
\mathbf t_2=-\cos\beta\,\mathbf t_1
+\sin\beta\,\mathbf n_1,\qquad
\mathbf n_2=\mathbf e_z\times\mathbf t_2.
\]

The joint-face centers coincide at `r_J=(0,0,0)`; the outer-face center is
`-L_i t_i`. Design-level orientation matrices at 0, 30, and 90 degrees have a
maximum orthogonality residual `7.437084071669e-18` and determinant residual
`0`, below the frozen future `1e-12` orientation threshold.

## 7. Overlapping-part idealization and mass gate

Near the point joint the full prisms may geometrically overlap, but they remain
independent substructures. There is no Boolean union, node merge, contact,
shared element, overlap subtraction, miter, rigid cube, adhesive, fillet, or
finite joint zone. Each arm retains its complete volume and mass. The overlap
is therefore not a claim about physical double material; it is an abstraction
of two beam-like continua connected only through one ideal point joint.

Future coupled models must verify

\[
m_{3D}=\rho(A_1L_1+A_2L_2),
\]

and S0 must verify `m_3D=rho*A*L`, with relative residual `<=1e-10`. Reference
point mass and rotary inertia must be exactly zero.

## 8. Ideal joint contract

One massless six-DOF reference point is placed at the origin. Both inner faces
are coupled to that same point through exact rigid-face kinematics using all
three translations and rotations. The point is not fixed and the two arms do
not receive independent joint rotations. Penalty stiffness is forbidden.

Before any future eigenfrequency extraction, six algebraic rigid-motion tests
must prescribe three unit translations and three unit rotations and verify

\[
\mathbf u(\mathbf r)=\mathbf U_J+
\boldsymbol\Omega_J\times(\mathbf r-\mathbf r_J)
\]

at every coupled face node.

## 9. Outer clamps and 1D comparators

Every node on each outer face has all three solid translations fixed. This is
a complete solid-face clamp, not a single-node, centroid-only, elastic, remote
spring, or transverse-only support.

The source-facing 1D condition is

```text
book_slope_clamp: w=0, w'=0, Phi=0,
w'=psi+Q/(kappa*Gxz*A).
```

A fully fixed solid face corresponds in the 1D limit primarily to

```text
timoshenko_section_clamp: w=0, psi=0, Phi=0.
```

Therefore the primary future comparator is `state_corrected` Timoshenko with
`timoshenko_section_clamp` at the external ends and the unchanged internal
rigid joint. `state_corrected + book_slope_clamp` remains only a
boundary-sensitivity/source diagnostic. The rectangular EB comparator remains
diagnostic; because `w'=psi` in EB, its two clamp descriptions coincide.
Future section-clamped assembly must remain script-local. UT-0--UT-3 used the
source-facing clamp and are not recalculated or relabeled.

## 10. Target parity family

The solid spectrum also contains axial, in-plane bending, cross-sectional, and
other families. The first seven sorted 3D modes must not be compared directly
with the seven 1D out-of-plane/torsion roots.

For a mesh and model symmetric under `z -> -z`, define

\[
\mathcal R_z\mathbf u(x,y,z)=Q_z\mathbf u(x,y,-z),\qquad
Q_z=\operatorname{diag}(1,1,-1).
\]

The target family has `phi approximately -R_z phi`: `u_x,u_y` are odd in `z`
and `u_z` is even. For a mass-normalized mode,

\[
p=\frac{\phi^TM\mathcal R_z\phi}{\phi^TM\phi}.
\]

The frozen criterion is `p<=-0.99`, or equivalently the antisymmetric residual

\[
r_-=\frac{\lVert\phi+\mathcal R_z\phi\rVert_M}
{2\lVert\phi\rVert_M}\le10^{-3}.
\]

The mesh must have reflected node pairs and identical topology above and below
`z=0`; clamps, material axes, and joint coupling must preserve the reflection.
For an eigenvalue cluster with relative gap `<1e-3`, native vectors are not
classified separately. Instead the future implementation forms

\[
P_c=V_c^TM\mathcal R_zV_c,
\]

diagonalizes it, and classifies the resulting parity-pure basis. This is not
MAC, physical branch tracking, or a shape-result study.

Future extraction requests 24 modes, then 32 and 40 only if fewer than seven
target modes have been identified. Failure to find seven by 40 is
`TARGET_FAMILY_NOT_FOUND_WITHIN_LIMIT`.

## 11. Frozen anchor cases

| Case | Geometry | Dimensions and angle | Future role |
| --- | --- | --- | --- |
| `S0_uniform_5mm` | one straight prism, both end faces clamped | `a=5 mm`, `b=20 mm`, `L=0.4 m`, `theta=0` | material, orientation, element, parity, and mesh-convergence anchor |
| `J30_equal_5_5` | two independent prisms and ideal point joint | `a1=a2=5 mm`, `b1=b2=20 mm`, `L1=L2=0.2 m`, `beta=30 deg` | generic equal-thickness joint anchor |
| `J30_unequal_4_6` | two independent prisms and ideal point joint | `a1=4 mm`, `a2=6 mm`, common `b=20 mm`, `L1=L2=0.2 m`, `beta=30 deg` | primary unequal-thickness anchor |
| `J90_unequal_4_6` | two independent prisms and ideal point joint | `a1=4 mm`, `a2=6 mm`, common `b=20 mm`, `L1=L2=0.2 m`, `beta=90 deg` | demanding quarter-turn anchor |

There is no swapped `(6,4)` case, negative angle, intermediate angle, or
thickness/angle sweep in the initial implementation contract.

## 12. Frozen mesh policy

Quadratic structured hexahedra are preferred. Quadratic tetrahedra are allowed
only as a documented fallback with their own convergence audit. Linear
tetrahedra are not an acceptable primary anchor.

| Level | Along 200 mm | Along S0 400 mm | Across 20 mm | Through 4 mm | Through 5 mm | Through 6 mm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| M0 | 50 | 100 | 10 | 2 | 3 | 3 |
| M1 | 75 | 150 | 15 | 3 | 4 | 5 |
| M2 | 100 | 200 | 20 | 4 | 5 | 6 |

Counts may be revised only before scientific results and only for documented
topology or memory constraints, never after inspecting frequency agreement.

## 13. Staged future plan and frozen thresholds

The staged plan is:

- 3D-A1: `S0_uniform_5mm` at M0, M1, and M2;
- 3D-A2: `J30_equal_5_5` at M2 only;
- 3D-A3: `J30_unequal_4_6` at M2 and `J90_unequal_4_6` at M1/M2.

No stage was started in 3D-A0.

Future mesh-convergence thresholds are

\[
\max_{k\le3}\frac{|f_k^{M2}-f_k^{M1}|}{f_k^{M2}}\le5\times10^{-3},
\]

\[
\max_{k\le6}\frac{|f_k^{M2}-f_k^{M1}|}{f_k^{M2}}\le10^{-2}.
\]

Using the section-clamped Timoshenko spectrum as the primary reference, future
coupled-case thresholds are

\[
E_3^{3D/T}=\max_{k\le3}
\frac{|f_k^{3D}-f_k^T|}{f_k^{3D}}\le0.05,
\]

\[
E_6^{3D/T}=\max_{k\le6}
\frac{|f_k^{3D}-f_k^T|}{f_k^{3D}}\le0.10.
\]

Mode 7 must be present and classified as target parity but has no separate
frequency threshold. Book-slope Timoshenko and EB differences are saved only
as diagnostics; neither receives a hard 3D acceptance threshold.

## 14. Readiness classification

Readiness is deliberately separated:

| Dimension | Result |
| --- | --- |
| material readiness | `BLOCKED_INCOMPLETE_3D_CONSTITUTIVE_DATA` |
| solver readiness | `BLOCKED_NO_CAPABLE_3D_SOLVER` |
| geometry readiness | `DESIGN_READY` |
| boundary-condition readiness | `DESIGN_READY` |
| mode-extraction readiness | `DESIGN_READY` |
| implementation readiness | `NOT_STARTED` |

Because material and solver are both unavailable, the single 3D-A0 status is

```text
BLOCKED_MATERIAL_AND_SOLVER
```

This is a normal scientific readiness result, not a failure of the accepted 1D
models. Permissible next audits are to locate a full HMS/DX-209 datasheet or its
primary reference, select a different material with a complete source-backed
tensor, or separately authorize an isotropic solver-infrastructure smoke that
is explicitly not anisotropic validation. None is performed here.

## 15. Evidence and exclusions

The readiness CLI writes only:

```text
solver_inventory.csv
constitutive_audit.csv
anchor_case_manifest.csv
mesh_policy.csv
readiness_summary.json
readiness_report.md
```

under
`results/anisotropic_rods/yartsev_ch2_limited_3d_fem_anchor/design_readiness/`.
Generated evidence is ignored and may be absent from a fresh clone.

3D-A0 explicitly excludes mesh/geometry creation, solver input decks, binary
results, solver or GUI execution, eigenfrequencies, physical mode shapes,
plots, PDF, package installation, network download, constitutive assumptions,
production material/API changes, and execution of 3D-A1/A2/A3.
