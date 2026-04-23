# Theorem Map

## Status labels

- `proved_here`
- `accepted_standard_background`
- `explicit_hypothesis`
- `needs_caution`
- `not_reached_yet`

## Consolidated step map

| id | name | statement | status | depends_on | what_it_gives | what_it_does_not_give |
| --- | --- | --- | --- | --- | --- | --- |
| PCL-SR | Simple-root packet construction | In the simple-root isolated-cluster setting, pointwise right/left null data can be assembled into local rank-2 packets `V(p),W(p)` once rank, regularity, and pairing hypotheses are supplied. | `proved_here` | regular local domain; isolated two-root cluster; simple roots; accepted standard background for local root/null-vector continuation | local retained right/left packet in the simple-root setting | complements; Schur block; normalization; symmetric form; `delta-kappa`; any non-simple-root theory |
| CR | Complement regularity | Given a local rank-2 packet, one can complete it to bases `[V Z]`, `[W U]`, but regular complement requires the additional hypothesis `det D(Lambda,p) != 0` on the claimed region. | `proved_here` | compatible right/left data; local basis completion; `D = U^T M Z` | regular complement and Schur-readiness | root capture; normalization; any claim that basis completion alone is enough |
| SRC | Schur root capture | If the basis extensions are invertible and `det D != 0` on the claimed region, then `det M` and `det S` have the same local zeros there. | `proved_here` | compatible packet; regular complement; standard Schur determinant identity | exact local reduced determinant for the retained cluster | normalized eigenvalue-style form; symmetric form; branch application |
| LN | Lambda-normalization | If `S` is analytic and `partial_Lambda S(Lambda0,p0)` is invertible, then locally `S = R((Lambda-Lambda0)I-K(Lambda,p))` with `det R != 0`. | `proved_here` | Schur root capture; analytic `S`; nonsingular base-point spectral derivative | exact local normalized block | a `p`-only reduced matrix; symmetric form; `delta-kappa` |
| FORM | First-order reduced model | Freezing `K(Lambda,p)` at `Lambda0` gives a `p`-only model `((Lambda-Lambda0)I-K0(p))` with matrix remainder `O(|Lambda-Lambda0|)`. | `proved_here` | exact normalized block; analytic `K` | controlled matrix-level first-order approximation | exact root capture for the frozen model; exact multiplicity transfer; parameterwise comparison |
| FRC | Fixed-`p` frozen-model root comparison | At fixed `p`, if the frozen determinant has one simple root in a local disc and the determinant perturbation is boundary-dominated, then the exact normalized determinant has one nearby root in that disc. | `proved_here` | exact normalized block; frozen model; determinant perturbation bound; accepted standard background via Rouche-type root comparison | first honest local spectral link between exact normalized and frozen models | global branch theorem in `p`; quantitative root shift; symmetric or veering interpretation |

## What is exact

Already exact in the current line:

- the full matrix object `M(Lambda,p)`;
- basis-completion and block-definition algebra once packets are supplied;
- the determinant relation between `M` and `M_tilde`;
- Schur elimination once `det D != 0` is assumed on the claimed region;
- exact local root capture by `det S` on the claimed region;
- exact local `Lambda`-normalization

  ```text
  S = R((Lambda-Lambda0)I-K(Lambda,p));
  ```

- the exact matrix remainder identity for first-order freezing;
- the exact determinant perturbation identity for fixed-`p` comparison.

## What is exact only under explicit hypotheses

These steps are exact as theorem statements only when their stated hypotheses
are carried along:

- packet construction in the simple-root setting;
- regular complement on the claimed local region;
- Schur root capture on the claimed local region;
- exact local `Lambda`-normalization near `(Lambda0,p0)`;
- fixed-`p` nearby-root existence in a chosen disc.

## What is only accepted standard background

Imported but not reproved in the package:

- simple-root corank-one consequences for finite-dimensional matrix problems;
- local continuation of simple roots;
- local right/left kernel representative selection;
- local regular basis completion after shrinking;
- Schur determinant formula for invertible `D`;
- analytic matrix difference-quotient decomposition;
- Rouche-type one-variable analytic root comparison.

## What is only approximate

Approximate, not exact:

- the frozen `p`-only reduced model

  ```text
  ((Lambda-Lambda0)I-K0(p));
  ```

- any spectral reading of that frozen model before extra fixed-`p` boundary
  hypotheses are imposed;
- any interpretation of that frozen model as symmetric, canonical, or already
  in `delta-kappa` form.

## Current next real obstruction

The next real obstruction is:

```text
turn the fixed-p spectral comparison into a controlled parameterwise
comparison in p, without skipping prematurely to symmetric normal form,
delta-kappa, or a final veering criterion.
```
