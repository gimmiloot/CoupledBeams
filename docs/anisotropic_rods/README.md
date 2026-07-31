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
Coupled-rod theory: not started
Coupled-rod code: not started
Euler--Bernoulli--Saint-Venant comparison: not started
Production anisotropic API: not started
FEM: not started
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
invalid and does not prescribe a future internal rigid joint.

## Navigation

- [Yartsev source note and local scan map](source_note_yartsev_2024.md)
- [Literature source-index entry](../literature/source_index.md#yartsev_2024_coupled_composite_structures)
- [Free-free Figure-2.2 reproduction note](yartsev_ch2_single_rod_reproduction.md)
- [Cantilever Figure-2.8 reproduction note](yartsev_ch2_cantilever_reproduction.md)
- [Free-free diagnostic CLI](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py)
- [Cantilever diagnostic CLI](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py)

Generated local evidence, which may be absent from a fresh clone, is kept in:

- [`results/anisotropic_rods/yartsev_ch2_free_free/`](../../results/anisotropic_rods/yartsev_ch2_free_free/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever/`](../../results/anisotropic_rods/yartsev_ch2_cantilever/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/)

## Next scientific gate

Derive and verify the rigid angular-joint conditions for two monoclinic rods
using the existing out-of-plane coordinate convention. In particular, the
future derivation must treat displacement compatibility, global physical
section-rotation compatibility, and vector force/moment equilibrium as an
internal-joint problem. No coupled-rod determinant or solver has been created.

## Preservation rule

> No existing isotropic equation, determinant entry, sign convention, unknown
> ordering, root solver, or FEM model may be changed to support this research
> direction without a separate consistency audit.
