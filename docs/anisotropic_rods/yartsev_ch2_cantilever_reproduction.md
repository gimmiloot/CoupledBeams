# Yartsev Chapter-2 cantilever reproduction

## 1. Scope

This source-reproduction gate concerns one straight rectangular monoclinic
rod, with Timoshenko bending coupled to generalized torsion. The source
comparison uses the book's complex moduli for HMS/DX-209. It does not model a
pair of coupled rods or an internal angular joint.

## 2. Source

The material record is taken from Chapter 1, Table 1.1 (printed page 45). The
cantilever formulation and results are in Chapter 2, Section 2.2.2 (printed
pages 64--68); Figure 2.8 is on printed page 65. Bibliographic data, the local
scan inventory, and source-specific corrections are recorded in the
[Yartsev source note](source_note_yartsev_2024.md). No book image or extended
passage is reproduced here.

## 3. Geometry and material

The implementation uses

| Quantity | Value |
| --- | ---: |
| `a` | `0.005 m` |
| `b` | `0.020 m` |
| `L` | `0.200 m` |
| `k` | `5/6` |
| `E1` | `191 GPa` |
| `E2` | `5 GPa` |
| `G12`, `G13` | `3 GPa`, `3 GPa` |
| `G23` | `2.5 GPa` |
| `nu12` | `0.279` |
| `rho` | `1580 kg/m^3` |

The loss factors used by the `book_complex` mode are the HMS/DX-209 values in
the canonical source record and implementation; they are not duplicated as a
second material table here.

## 4. Corrected Chapter-2 model

`state_corrected` is the canonical source-consistent formulation for further
work. It restores the dimensionally required `I_y` in the rotary-inertia term
of equation (2.1) and uses the internally consistent coefficient signs.
`eliminated_corrected` is an independent verification form of the same model.
In particular, the corrected eliminated form uses positive `d0` and `f0`.
`eliminated_printed` preserves the negative signs printed after equation
(2.16) only as historical diagnostic evidence; it is not the physical
baseline.

## 5. Two clamp variants

For the state ordering
`y = [w, psi_b, Phi_t, Q, M, M_T]^T`, both variants have the free-end
conditions `Q(L)=M(L)=M_T(L)=0`. Their clamped-end conditions are:

- `book_slope_clamp`: `w(0)=0`, `w'(0)=0`, `Phi(0)=0`, where
  `w' = psi_b + Q/(k Gxz A)`;
- `timoshenko_section_clamp`: `w(0)=0`, `psi_b(0)=0`, `Phi(0)=0`.

With `K_s = k Gxz A`, the implemented state-space basis matrices are

```text
B_book = [[0,       0, 0],
          [-1/K_s,  0, 0],
          [0,       0, 0],
          [1,       0, 0],
          [0,       1, 0],
          [0,       0, 1]]

B_section = [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 0],
             [1, 0, 0],
             [0, 1, 0],
             [0, 0, 1]]
```

Here `w'` is the centerline slope, whereas `psi_b` is the independent section
rotation. They need not coincide at finite shear deformation. The matrices
are implemented in
[`yartsev_ch2_monoclinic_rod.py`](../../scripts/lib/yartsev_ch2_monoclinic_rod.py).

## 6. Quick elastic boundary gate

The smoke gate used `theta = [0, 15, 45, 90] deg`; the completed quick gate
used `theta = [0, 15, 30, 45, 60, 75, 90] deg`. It retained the first five
positive frequencies and a sixth root as a completeness guard for both clamp
variants. The maximum normalized determinant residual was `4.344008e-14`,
the maximum relative singular residual was `1.411361e-13`, and halving the
scan step at `theta=15 deg` changed the first six roots by at most
`1.894311e-12` relatively. At `theta=0 deg`, the independent torsional roots
agreed with the analytical cantilever torsion formula to `1.899077e-13`
relatively.

The largest difference between the two spectra was `4.842132%` in frequency
and `9.449802%` in squared frequency, at `theta=0 deg`, sorted mode 5. Pure
torsional roots at that angle coincide; the main sensitivity is in the
bending modes and decreases rapidly with increasing `theta`. This preliminary
gate deliberately did not select a preferred clamp. The complete tables are
in the generated quick-gate CSV files; no possible sorted-root reordering was
flagged on this grid.

## 7. Figure 2.8 source check

The frequency panel alone gave `SOURCE_AMBIGUITY`: the difference between
clamps was graphically resolvable at only 1 of 35 digitized points. The
declared, pre-comparison frequency-reading uncertainty was `±0.08 kHz`.

Adding the calculated loss-factor curves gave
`BOOK_SLOPE_CLAMP_CONFIRMED`. The declared loss-reading uncertainty was
`±0.04` in `eta*1e2`. The two independent decisive points at `theta=0 deg`
were:

| Mode | Source `eta*1e2` | `book_slope_clamp` | `timoshenko_section_clamp` |
| ---: | ---: | ---: | ---: |
| 3 | 0.130 | 0.132935 | 0.202015 |
| 5 | 0.240 | 0.244444 | 0.329288 |

The book-slope results are within the predeclared graphical uncertainty in
both points; the section-clamp results are outside it. Existing saved complex
roots and continuity labels were reused, no parameters were fitted, and the
postprocessor recorded:

```text
determinant evaluations = 0
root solves = 0
complex continuations = 0
scientific matrix-exponential calls = 0
```

Figure 2.11 was not needed for the clamp decision.

## 8. Accepted source boundary condition

The source-faithful cantilever condition is

```text
w(0) = 0
w'(0) = 0
Phi(0) = 0
```

It is confirmed as the boundary condition used for the calculated curves of
Figure 2.8 within graph resolution. This statement does not make it the only
physically meaningful Timoshenko clamp: `timoshenko_section_clamp` remains a
valid boundary-sensitivity alternative, but it is not the source baseline.

## 9. Relation to the future coupled-rod model

An external source clamp and an internal rigid angular joint are different
objects. Figure 2.8 does not determine the future joint conditions. A separate
scientific gate must derive and verify displacement continuity, continuity of
the global physical section-rotation vector
`theta_i = Phi_i t_i + psi_b,i n_i`, force equilibrium, and moment-vector
equilibrium using the existing out-of-plane coordinate convention. These
conditions have not been derived here, and no coupled-rod determinant exists.

## 10. Reproduction commands

The actual CLI supports these entry modes:

```text
# ordinary full orientation and length workflow
python scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py --study both --clamp-variant both --material-mode both

# reduced smoke run
python scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py --smoke

# elastic boundary-sensitivity gate
python scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py --quick-boundary-gate

# saved-data-only Figure-2.8 boundary-source check
python scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py --postprocess-boundary-source-check
```

The postprocess-only path exits before the scientific root workflow is
imported or called.

## 11. Generated outputs

- [`results/anisotropic_rods/yartsev_ch2_cantilever/`](../../results/anisotropic_rods/yartsev_ch2_cantilever/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/)

These are generated local results and may be absent from a fresh public clone.
The postprocess source check supersedes the provisional clamp ambiguity in the
full-run generated report; it does not alter that preserved report.

## 12. Final status and limitations

```text
cantilever source reproduction: completed
source boundary condition: BOOK_SLOPE_CLAMP_CONFIRMED
physical section-clamp alternative: retained as diagnostic sensitivity case
coupled rods: not started
Euler--Bernoulli--Saint-Venant comparison: not started
FEM: not started
```
