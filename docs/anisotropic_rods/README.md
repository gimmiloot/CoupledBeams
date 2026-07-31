# Anisotropic Rods — Chapter-2 source line consolidated

```text
Status: active-research
Source assimilation: completed for the selected Chapter-2 single-rod model
Primary source: registered
Free-free Chapter-2 reproduction: completed — PASS_WITHIN_GRAPH_RESOLUTION
Cantilever Chapter-2 reproduction: completed
Cantilever source clamp: BOOK_SLOPE_CLAMP_CONFIRMED
Canonical source-consistent formulation: state_corrected
Printed sign variant: diagnostic only
Rigid angular-joint gate: completed — PASS
Coupled-rod code: small elastic diagnostic pilot only — PASS
Final coupled-rod model or parameter study: not started
Rectangular orthotropic Euler--Bernoulli--Saint-Venant gate: PARTIAL_PASS after targeted refinement
Production anisotropic API: not started
Independent rectangular EB + torsion 1D FEM: completed for the finite gate
```

## Established single-rod baseline

The selected source line is one rectangular monoclinic Timoshenko rod with
generalized torsion. The internally consistent `state_corrected` formulation
restores `I_y` in equation (2.1). `eliminated_corrected` is its independent
equivalence check, while `eliminated_printed` preserves the signs printed
after equation (2.16) only for historical diagnosis.

The free-free Figure-2.2 benchmark is
`PASS_WITHIN_GRAPH_RESOLUTION`. The cantilever Figure-2.8 source check is
`BOOK_SLOPE_CLAMP_CONFIRMED`: the book's calculated curves use
`w(0)=0`, `w'(0)=0`, `Phi(0)=0`. This identifies the source boundary
condition; it does not make the alternative section-rotation clamp physically
invalid. The separate rigid-joint gate has now derived that internal joint
from physical vector compatibility and equilibrium.

## Rigid angular-joint pilot

The two-arm helper retains the book state
`[w, psi, Phi, Q, M, M_T]^T`, the project bases `e_z,t_i,n_i`, and the
accepted `state_corrected` arm equations. Its `6 x 12` joint matrix passed the
old-notation sign gate, rank and random-state checks, virtual work, requested
angle limits, orthotropic block separation, equal-arm exchange, and the
independent `beta=0` straight-rod spectrum check. The small real-elastic
HMS/DX-209 pilot at `beta=0,30,90 deg` is `PASS`; it is not a stable baseline
or a final parameter study.

## Rectangular orthotropic EB validation

The `theta=0` comparator reuses the existing rectangular `Cbar`, unchanged
`J_book`, book state, and project bases. Unequal-length invariance, exact
fixed--fixed bending/torsion families, the Timoshenko slender limit, FEM
structure, and the original fixed-elements-per-arm convergence checks passed.
The original gate is `PARTIAL_PASS` because 64 elements per arm gives unequal
physical element lengths and misses the strict first-three `1e-5` target
(`5.18e-5` worst case). The targeted proportional sequence through
`(N_1,N_2)=(64,192)` closes that accuracy threshold with the raw value
`6.18e-6`, without changing any model coefficient or threshold, but its
separate unchanged monotonic-refinement gate fails for bending mode 1 at the
dense-eigensolver conditioning floor. Its targeted status is
`FAIL_CONVERGENCE_ORDER`, so the overall status remains `PARTIAL_PASS`. This
is a project comparator, not a source reproduction, general anisotropic EB
model, final coupled baseline, unequal-thickness study, or 3D validation.

## Navigation

- [Yartsev source note and local scan map](source_note_yartsev_2024.md)
- [Literature source-index entry](../literature/source_index.md#yartsev_2024_coupled_composite_structures)
- [Free-free Figure-2.2 reproduction note](yartsev_ch2_single_rod_reproduction.md)
- [Cantilever Figure-2.8 reproduction note](yartsev_ch2_cantilever_reproduction.md)
- [Chapter-2 notation translation](yartsev_ch2_notation_translation.md)
- [Rigid angular-joint theory and gate](yartsev_ch2_rigid_angular_joint.md)
- [Rectangular orthotropic EB validation](yartsev_ch2_rectangular_eb_validation.md)
- [Free-free diagnostic CLI](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py)
- [Cantilever diagnostic CLI](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py)
- [Small coupled-joint pilot CLI](../../scripts/analysis/anisotropic_rods/pilot_yartsev_ch2_coupled_rods.py)
- [Rectangular EB validation CLI](../../scripts/analysis/anisotropic_rods/validate_yartsev_ch2_rectangular_eb.py)

Generated local evidence, which may be absent from a fresh clone, is kept in:

- [`results/anisotropic_rods/yartsev_ch2_free_free/`](../../results/anisotropic_rods/yartsev_ch2_free_free/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever/`](../../results/anisotropic_rods/yartsev_ch2_cantilever/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/)
- [`results/anisotropic_rods/yartsev_ch2_coupled_joint_pilot/`](../../results/anisotropic_rods/yartsev_ch2_coupled_joint_pilot/)
- [`results/anisotropic_rods/yartsev_ch2_rectangular_eb_validation/`](../../results/anisotropic_rods/yartsev_ch2_rectangular_eb_validation/)

## Current boundary

The rigid-point-joint pilot and finite `theta=0` rectangular EB/1D-FEM gate
are complete. The latter closed its original strict accuracy shortfall under
equal-element-length refinement but was not promoted because the separate
mode-1 monotonic-refinement gate failed. Any off-axis or broad spectrum study, final coupled-rod model,
damping calculation, unequal-thickness or 3D validation, or production API
requires a separately scoped next stage.

## Preservation rule

> No existing isotropic equation, determinant entry, sign convention, unknown
> ordering, root solver, or FEM model may be changed to support this research
> direction without a separate consistency audit.
