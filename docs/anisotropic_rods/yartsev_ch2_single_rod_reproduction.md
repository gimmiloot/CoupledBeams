# Yartsev Chapter-2 Single-Rod Reproduction

## Scope and status

Status: `PASS_WITHIN_GRAPH_RESOLUTION` for the calculated curves in Figure 2.2.

This diagnostic gate reproduces one free-free monoclinic Timoshenko rod. It is
not a coupled-rod model, production anisotropic API, Euler--Bernoulli model, or
FEM validation. The seven reported tones are the first seven positive roots;
the eighth positive elastic root is retained as a completeness guard, and the
three rigid-body roots at zero are excluded.

## Sources and parameters

Only the registered local Chapter-1 and Chapter-2 source fragments and the
tracked source-consistency records were used. Chapter 1 equations (1.32)/(1.34) define
`M*=Re(M*)[1+i eta]`; equations (1.41)/(1.42) define the exact and small-loss
modal coefficients. Material rotation follows (1.52), the pure shear modulus
follows (1.56), and the one-dimensional model follows Chapter 2 equations
(2.1)--(2.18). Table 1.2 and the Figure-2.2 specimen description give:

| quantity | value |
| --- | ---: |
| `a` | `9.76e-3 m` |
| `b` | `2.524e-2 m` |
| `L` | `0.295 m` |
| `k` | `5/6` |
| `Re(E1)` | `28.4 GPa` |
| `Re(E2)` | `11.6 GPa` |
| `Re(G12)` | `4.4 GPa` |
| `Re(G13)` | `3.7 GPa` |
| `Re(G23)` | `2.6 GPa` |
| `nu12` | `0.22` |
| `rho` | `1650 kg/m^3` |
| `eta1, eta2` | `4.3e-3, 13.9e-3` |
| `eta12, eta13, eta23` | `2.0e-2, 2.74e-2, 3.1e-2` |
| `h_ply` | `3.5e-4 m` (manifest only; unused by the 1D equations) |

## Printed and corrected variants

- `state_corrected` restores the dimensionally required `I_y` in the rotary
  inertia term printed in equation (2.1).
- `eliminated_corrected` uses
  `d0=+rho*omega^2*I_p/Cbar` and
  `f0=+rho*omega^2*mu_xy_x/(2*k*Gxz)`.
- `eliminated_printed` retains the negative `d0` and `f0` printed after
  equation (2.16) and is diagnostic only.

The state and corrected eliminated formulations reproduce the same first eight
positive elastic roots to strict numerical tolerance and have small boundary
and field-mapping residuals. The printed signs suppress the independent
positive torsional family in the orthotropic endpoints and therefore produce a
different sorted spectrum; they are not repaired by parameter or sign fitting.

## Run and outputs

```text
python scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py
```

Useful options are `--material-mode elastic|book-complex|both`,
`--equation-variant corrected|printed|both`, `--theta-step-deg`,
`--num-positive-modes`, `--output-dir`, and `--smoke`.

The generated [full report](../../results/anisotropic_rods/yartsev_ch2_free_free/single_rod_reproduction_report.md),
[source manifest](../../results/anisotropic_rods/yartsev_ch2_free_free/source_manifest.md),
root CSV files, mode-continuity diagnostics, and reproduction figures remain
under the Git-ignored results directory. Root tables are sorted by `Re(f)` at
each angle and record the MAC-continuity label separately. The source-style
Figure-2.2 plot uses continuity labels through the two 6/7 sorted-position
exchanges.

## Figure comparison limits

The book contains no digital table for Figure 2.2. A separate CSV records a
manual, rounded reading of the calculated solid curves at
`theta=0,15,...,90 deg`, with estimated graphical uncertainties of `±0.06 kHz`
and `±0.06` in `eta*1e2`. All 98 calculated comparisons lie within those
declared reading uncertainties. This is not an exact percentage agreement or
a statistical confidence statement. Experimental markers were not digitized,
and digitized values were not used to tune the model.

## Related source gate

The subsequent cantilever source gate is documented separately in the
[Chapter-2 cantilever reproduction note](yartsev_ch2_cantilever_reproduction.md).
It confirms the boundary condition used for the calculated Figure-2.8 curves
and identifies the separate derivation of internal rigid-joint conditions as
the next scientific gate; this free-free note remains the canonical Figure-2.2
benchmark record.
