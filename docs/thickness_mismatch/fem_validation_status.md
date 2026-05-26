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
to Timoshenko than before. The current frequency-closeness count after
MAC-based mode selection is 23/36 closer to Timoshenko and 13/36 closer to
Euler--Bernoulli. This is promising diagnostic evidence, but it is not yet
article-level validation.

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
- The current point-joint results are promising, but still diagnostic.

## What Must NOT Be Claimed Yet

- "3D FEM fully validates the coupled Timoshenko model."
- "Fused-cylinder disagreement disproves the Timoshenko coupling."
- "Sorted FEM mode number equals analytic mode number."
- "All modes are reliably identified."

## Remaining Limitations

1. No mesh convergence study has been completed for the point-joint coupled FEM
   workflow.
2. MAC matching is currently based on centerline displacement only.
3. Duplicate best-solid-mode matches exist.
4. Weak MAC rows remain.
5. A unique assignment step, for example Hungarian matching, is still needed.
6. Representative strong, moderate, and weak mode shapes need visual
   inspection.
7. The frequency effect of the CalculiX `*RIGID BODY` constraint should be
   checked, including possible artificial stiffness or mass effects.
8. Article-grade conclusions must remain separated from diagnostic
   observations.

## Article-Ready Validation Roadmap

### Step 1. Point-Joint Mesh Convergence

- Start with `beta = 15 deg`.
- Use `epsilon = 0.025` and `epsilon = 0.05`.
- Run at least 2 or 3 mesh levels.
- Compare both frequencies and MAC stability.

### Step 2. Unique MAC Assignment

- Add Hungarian assignment for analytic modes 1--6.
- Report strong, moderate, and weak matches separately.
- Keep duplicate raw best-match warnings visible as diagnostics.

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

