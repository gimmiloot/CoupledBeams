# FEM Validation Status

## Purpose

This document records the current diagnostic-only FEM validation status for
the thickness-mismatch / Timoshenko correction work. The FEM workflows are used
as independent checks of the one-dimensional models:

- Euler--Bernoulli beam references;
- diagnostic Timoshenko beam references;
- the planned thickness-mismatch / Timoshenko formulation.

The notes below are not article-level claims. They define what is currently
supported, what remains diagnostic, and what must be improved before using FEM
results in a paper.

## Current FEM Workflows

### `solid_fem_single_rod_fixed_fixed.py`

This workflow checks one circular rod with fixed-fixed ends using 3D solid FEM
through Gmsh and CalculiX. It generates the solid mesh, solves modal
frequencies when the external tools are available, parses displacement modes,
classifies bending/axial/torsion-like shapes, and pairs the circular-rod
bending doublets.

Current diagnostic conclusion: for the classified bending doublets, the
Timoshenko reference is closer to the 3D solid FEM result than the
Euler--Bernoulli reference.

### `solid_fem_coupled_equal_rods.py`

This workflow checks equal coupled rods using fused-cylinder geometry for
`beta = 15, 45, 90 deg`. It is useful as a geometric diagnostic because it
shows how the finite fused solid joint differs as the angle changes.

Current diagnostic conclusion: this workflow is not a direct validation of the
ideal one-dimensional point-joint model. The fused 3D joint has finite volume,
and at small angles the overlap/merger zone can be large relative to a rod
segment.

### `solid_fem_coupled_equal_rods_point_joint.py`

This workflow checks equal coupled rods using two separate solid cylinders and
a shared CalculiX `*RIGID BODY` reference joint at the inner end faces. This is
closer to the ideal one-dimensional point joint than the fused-cylinder model.
It also includes MAC-like centerline shape matching between parsed 3D solid
modes and reconstructed one-dimensional EB/Timoshenko shapes.

Current diagnostic conclusion: after MAC-like matching, more modes are closer
to Timoshenko than before. The current full beta/epsilon frequency-closeness
count after MAC-based mode selection is 23/36 closer to Timoshenko and 13/36
closer to Euler--Bernoulli. A first beta=15 mesh-convergence run has now been
completed for `epsilon = 0.025, 0.05` and mesh factors `1.5, 1.0, 0.75`;
the strong/moderate MAC rows were stable under the final refinement screen,
but weak rows and duplicate raw best-match warnings remain. An auxiliary
planar-constrained diagnostic has also been added for `beta = 15 deg` and
`epsilon = 0.025, 0.05`; it suppresses out-of-plane/torsion-like solid modes
to check the planar subspace corresponding to the 1D model. This is promising
diagnostic evidence, but it is not yet article-level validation.

## What Is Already Supported

- The single-rod 3D solid FEM workflow gives strong support for the
  Timoshenko correction in the fixed-fixed circular-rod benchmark.
- The fused-cylinder coupled FEM workflow should be treated as a geometry audit,
  not as direct validation or refutation of the point-joint one-dimensional
  model.
- The point-joint 3D solid FEM workflow is the current best validation path for
  coupled rods because it avoids the finite fused overlap volume.
- MAC-like centerline shape matching is preferred over sorted-frequency
  comparison and over simply taking the first classified in-plane modes.
- The current point-joint results are promising, and the first beta=15
  mesh-convergence diagnostic is stable for strong/moderate MAC rows, but the
  workflow remains diagnostic until weak/duplicate matches and visual
  mode-shape checks are resolved.
- The planar-constrained point-joint run is an auxiliary diagnostic of the
  planar subspace only; it must be labeled as such and must not replace the
  full 3D point-joint benchmark.
- Equal-thickness `Lambda(mu)` diagnostics must report Euler--Bernoulli
  diameter applicability and Timoshenko cut-off applicability separately:
  the EB `2*r_i/l_i <= 0.1` criterion does not dash or invalidate
  Timoshenko curves, whose diagnostic margin is `Omega/Omega_c`.

## What Must NOT Be Claimed Yet

- "3D FEM fully validates the coupled Timoshenko model."
- "Fused-cylinder disagreement disproves the Timoshenko coupling."
- "Sorted FEM mode number equals analytic mode number."
- "All modes are reliably identified."

## Remaining Limitations

1. Mesh convergence has been started for the point-joint coupled FEM workflow:
   beta=15, `epsilon = 0.025, 0.05`, and mesh factors `1.5, 1.0, 0.75` are
   recorded in `results/coupled_equal_rods_point_joint_mesh_convergence_report.md`.
   Extension beyond this beta/epsilon diagnostic is still pending.
2. MAC matching is currently based on centerline displacement only.
3. Duplicate best-solid-mode matches exist.
4. Weak MAC rows remain.
5. Preliminary Hungarian assignment rows are now written for the convergence
   diagnostic, but the article workflow still needs reviewed unique assignment
   policy before final tables.
6. Representative strong, moderate, and weak mode shapes need visual
   inspection.
7. The planar-constrained diagnostic required an MPC-compatible fallback
   because direct `ALL_SOLID_NODES, 3, 3, 0` constraints conflict with rigid
   body dependent nodes in CalculiX; this constraint implementation must be
   described clearly if used.
8. The frequency effect of the CalculiX `*RIGID BODY` constraint should be
   checked, including possible artificial stiffness or mass effects.
9. Article-grade conclusions must remain separated from diagnostic
   observations.

## Article-Ready Validation Roadmap

### Step 1. Point-Joint Mesh Convergence

- Start with `beta = 15 deg`.
- Use `epsilon = 0.025` and `epsilon = 0.05`.
- Status: first diagnostic run completed for mesh factors `1.5`, `1.0`, and
  `0.75`; strong/moderate final-mesh rows passed the current stability screen.
- Remaining: inspect weak/duplicate cases and decide whether additional
  refinement or a stricter unique-assignment policy is needed before article
  use.

### Step 2. Unique MAC Assignment

- Add Hungarian assignment for analytic modes 1--6.
- Report strong, moderate, and weak matches separately.
- Keep duplicate raw best-match warnings visible as diagnostics.

### Auxiliary. Planar-Constrained Point-Joint Check

- Status: beta=15 diagnostic completed for `epsilon = 0.025` and `0.05`.
- Purpose: check whether suppressing out-of-plane/torsion-like 3D families
  clarifies comparison with the planar EB/Timoshenko model.
- Limitation: this is not a full 3D validation replacement and article use
  requires clear labeling as a constrained planar-subspace diagnostic.

### Auxiliary. Equal-Thickness `Lambda(mu)` FEM Overlay

- Status: corrected beta=15, `epsilon = 0.01`, `eta = 0` diagnostic written
  with analytic descendant branches and sparse planar point-joint FEM points.
- Purpose: inspect mode identity and 3D-vs-1D deviations across `mu` while
  keeping EB diameter validity, Timoshenko cut-off validity, mesh quality, and
  MAC strength separate.
- Limitation: large FEM discrepancies remain diagnostic; cases with negative
  Jacobian mesh warnings or weak MAC are excluded from the main FEM points.
- Follow-up audit: the `mu=0`, `epsilon=0.01` sanity audit confirms the
  `Lambda_FEM=sqrt(Omega/epsilon)` conversion and the 1D EB frame-FEM anchor,
  but marks the variable-length Timoshenko extension as incomplete and leaves
  the planar 3D solid discrepancy unexplained pending constraint-isolation or
  full point-joint `epsilon=0.01` checks.

### Step 3. Shape Figures

- Export representative 3D mode shapes.
- Compare them with the corresponding one-dimensional centerline shapes.
- Use only clear strong-MAC cases in article figures.

### Step 4. Article Validation Table

- Use only modes with strong, or at least moderate, MAC support.
- Report Euler--Bernoulli, Timoshenko, and 3D solid FEM frequencies together.
- Mark weak or ambiguous modes as excluded or diagnostic-only.

### Step 5. Extend To Thickness Mismatch

- Extend to `eta != 0` only after equal-rod point-joint validation is stable.
- Keep thickness-mismatch FEM runs separate from the equal-rod validation
  baseline.

## Current Outputs Index

| Workflow | Output |
| --- | --- |
| Single fixed-fixed rod report | `results/single_rod_fixed_fixed_3d_solid_fem_report.md` |
| Single fixed-fixed rod bending doublets | `results/single_rod_fixed_fixed_3d_solid_fem_bending_doublet_comparison.csv` |
| Fused coupled rods report | `results/coupled_equal_rods_3d_solid_fem_report.md` |
| Fused coupled rods in-plane comparison | `results/coupled_equal_rods_3d_solid_fem_in_plane_comparison_by_beta.csv` |
| Point-joint coupled rods report | `results/coupled_equal_rods_point_joint_3d_solid_fem_report.md` |
| Point-joint MAC matches | `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matches.csv` |
| Point-joint EB MAC matrix | `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_eb.csv` |
| Point-joint Timoshenko MAC matrix | `results/coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_timoshenko.csv` |
| Point-joint mesh convergence report | `results/coupled_equal_rods_point_joint_mesh_convergence_report.md` |
| Point-joint planar-constrained report | `results/coupled_equal_rods_point_joint_planar_3d_solid_fem_report.md` |
| Equal-thickness corrected `Lambda(mu)` FEM overlay | `results/coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected_report.md` |
| Equal-thickness `mu=0` FEM sanity audit | `results/coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.md` |
