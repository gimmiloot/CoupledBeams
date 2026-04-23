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
| PMRC | Uniform parameterwise moving-disc root comparison | If a frozen-model root branch `Lambda_fr(p)` exists on a local interval, one uniform radius `r` works for all moving discs `D_r(p)`, and uniform boundary domination holds, then each disc contains exactly one exact normalized root. | `proved_here` | exact normalized block; frozen-model root branch; moving discs; determinant perturbation bound; accepted standard background via Rouche-type comparison | parameterwise nearby-root selection over a local interval | regular exact root branch without separate simple-root continuation; quantitative exact-vs-frozen distance; symmetric or veering interpretation |
| ERBC | Exact-root branch continuation | Starting from the moving-disc selected exact root at `p0`, if that exact root is simple, then a local regular exact root branch `Lambda_ex(p)` exists and is identified with the moving-disc selected root by uniqueness in each disc. | `proved_here` | parameterwise moving-disc comparison; joint analyticity; accepted standard background via implicit-function-style simple-root continuation | local exact root branch consistent with the moving-disc selection | quantitative branch separation from the frozen root; global continuation in `p`; symmetric or veering interpretation |
| QPBC | Quantitative parameterwise branch comparison | If exact branch continuation, moving-disc isolation, a uniform frozen simple-root lower bound, and determinant discrepancy control are already in place, then `|Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p)/m0`. | `proved_here` | exact-root branch continuation; uniform moving-disc isolation; determinant discrepancy `eta_loc(p)`; accepted standard background for simple-root factorization / derivative lower bound | a conditional local quantitative comparison estimate between the exact and frozen branches | derivative-level shift formulas; symmetric normal form; `delta-kappa`; final veering criterion |
| DDMR | Determinant-discrepancy bound from matrix remainder | If `E=F-G` on the moving discs, then in the local `2x2` Frobenius-norm setting the exact identity `det(G+E)-det(G)=tr(adj(G)E)+det(E)` yields `eta_loc(p) <= M_G(p) M_E(p) + (1/2) M_E(p)^2`; combined with `QPBC`, this gives a constructive conditional bound for `|Lambda_ex(p)-Lambda_fr(p)|`. | `proved_here` | first-order reduced model; moving discs; local matrix controls `M_G,M_E`; accepted standard background for the Frobenius determinant bound | a constructive route from matrix-level remainder control to determinant discrepancy and then to branch-distance control | derivative law, asymptotic shift formula, symmetry, `delta-kappa`, final veering criterion |
| MDC | Moving-disc closure | If the constructive branch bound `B_shift(p) = [M_G(p) M_E(p) + (1/2) M_E(p)^2]/m0` satisfies `B_shift(p) < r` or uniformly `B_* < r`, then the exact branch selected earlier is certified to lie inside the chosen moving discs. | `proved_here` | QPBC; DDMR; chosen radius `r` | self-consistency of the moving-disc family with the current quantitative branch bound | a new existence theorem, derivative law, symmetry, `delta-kappa`, final veering criterion |
| FOSF | First-order root-shift formula | If the exact branch already exists, the moving-disc closure condition holds, `partial_Lambda g_p(Lambda_fr(p))` has a positive lower bound, and local derivative controls for `f_p-g_p` and `g_p` hold on the moving disc, then `Lambda_ex(p)-Lambda_fr(p) = -(f_p(Lambda_fr(p))-g_p(Lambda_fr(p)))/partial_Lambda g_p(Lambda_fr(p)) + Rem_shift(p)` with an explicit local remainder bound. | `proved_here` | exact-root branch continuation; moving-disc closure; first-order derivative lower bound; derivative discrepancy control; second-derivative control of `g_p` | a conditional local first-order comparison law for the branch shift | a final asymptotic law, constructive derivative controls from the matrix model, symmetry, `delta-kappa`, final veering criterion |
| DCL | Derivative control from reduced-model objects | If the frozen model `G` and matrix remainder `E` are defined on the moving discs, then in the local `2x2` setting one has the exact identities `partial_Lambda g_p(Lambda)=tr G(Lambda,p)=2(Lambda-Lambda0)-tr K0(p)` and `partial_Lambda^2 g_p(Lambda)=2`, together with the constructive bound `sup_{D_r(p)} |partial_Lambda(f_p-g_p)| <= sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p)`. | `proved_here` | first-order reduced model; moving discs; local matrix controls `M_G,M_E,M_{E,1}`; accepted standard background for Frobenius trace inequalities | a constructive bridge from reduced-model matrices to the derivative quantities entering `FOSF` | a uniform positive derivative lower bound, a final asymptotic shift law, symmetry, `delta-kappa`, final veering criterion |
| NCDL | Frozen-model denominator control from noncoalescence | If a second local frozen root branch `Lambda_fr,2(p)` is available for the frozen quadratic and the frozen roots remain separated, then `partial_Lambda g_p(Lambda_fr(p)) = Lambda_fr(p) - Lambda_fr,2(p)` and any lower bound on frozen-root separation, or exclusion of the second frozen root from `D_r(p)`, gives the denominator lower bound needed in `FOSF`. | `proved_here` | frozen `2x2` model; chosen frozen root branch `Lambda_fr(p)`; local second frozen root branch `Lambda_fr,2(p)` | converts the remaining denominator hypothesis in `FOSF` into a frozen-root noncoalescence or second-root-exclusion condition | existence of the second frozen root branch for CoupledBeams, proof of a separation bound, asymptotic shift law, symmetry, `delta-kappa`, final veering criterion |
| FDPL | Frozen discriminant persistence and denominator control | If the frozen discriminant `Disc_fr(p0)` is nonzero at a base point, then after shrinking to a smaller interval `I0` and locally labeling the two frozen roots there, `Disc_fr` remains nonzero and there exists `c0 > 0` such that `|partial_Lambda g_p(Lambda_fr(p))| >= c0` for all `p in I0`. | `proved_here` | frozen `2x2` model; frozen discriminant; accepted standard background for local labeling of noncoalesced quadratic roots | a local coefficient-level route from base-point frozen noncoalescence to the denominator lower bound used in `FOSF` | proof that the needed discriminant lower bound holds for CoupledBeams, asymptotic shift law, symmetry, `delta-kappa`, final veering criterion |
| FDCC | Frozen discriminant persistence from coefficient control | If `T(p)=tr K0(p)` and `D0(p)=det K0(p)` satisfy explicit local bounds `|T(p)| <= M_T`, `|T(p)-T(p0)| <= eps_T`, `|D0(p)-D0(p0)| <= eps_D`, then `|Disc_fr(p)-Disc_fr(p0)| <= 2 M_T eps_T + 4 eps_D`; hence, if `|Disc_fr(p0)|` exceeds that margin, then `Disc_fr` stays nonzero on the interval, and after shrinking and local root labeling one gets the local denominator lower bound from `FDPL`. | `proved_here` | frozen discriminant notation; coefficient functions `T,D0`; FDPL downstream | a simple coefficient-level criterion for local frozen noncoalescence persistence | proof that the needed coefficient bounds hold on the actual CoupledBeams regime, asymptotic shift law, symmetry, `delta-kappa`, final veering criterion |

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
- the exact moving-disc zero count conclusion once the uniform domination
  hypotheses are imposed;
- local exact-root branch continuation once the selected exact root is simple.
- the quantitative bound

  ```text
  |Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p)/m0
  ```

  once the exact branch exists and the frozen branch satisfies a uniform
  simple-root lower bound on the moving discs.
- the constructive determinant-discrepancy estimate

  ```text
  eta_loc(p) <= M_G(p) M_E(p) + (1/2) M_E(p)^2
  ```

  in the local `2x2` Frobenius-norm setting.
- the moving-disc closure implication

  ```text
  B_shift(p) < r => Lambda_ex(p) in int D_r(p)
  ```

  and its uniform version with `B_*`.
- the first-order shift identity with remainder once the derivative lower
  bound and derivative regularity hypotheses are imposed.
- the exact frozen-derivative identities

  ```text
  partial_Lambda g_p(Lambda)=tr G(Lambda,p),
  partial_Lambda^2 g_p(Lambda)=2
  ```

  and the constructive derivative-discrepancy bound from `G`, `E`, and
  `partial_Lambda E`.
- whenever two local frozen root branches are available, the exact
  factorization

  ```text
  g_p(Lambda)
  =
  (Lambda-Lambda_fr(p))(Lambda-Lambda_fr,2(p))
  ```

  and the denominator identity

  ```text
  partial_Lambda g_p(Lambda_fr(p))
  =
  Lambda_fr(p)-Lambda_fr,2(p).
  ```
- the exact frozen quadratic formula and discriminant identity

  ```text
  g_p(Lambda)
  =
  (Lambda-Lambda0)^2 - tr K0(p) (Lambda-Lambda0) + det K0(p),
  Disc_fr(p) = (tr K0(p))^2 - 4 det K0(p).
  ```

## What is exact only under explicit hypotheses

These steps are exact as theorem statements only when their stated hypotheses
are carried along:

- packet construction in the simple-root setting;
- regular complement on the claimed local region;
- Schur root capture on the claimed local region;
- exact local `Lambda`-normalization near `(Lambda0,p0)`;
- fixed-`p` nearby-root existence in a chosen disc.
- uniform parameterwise moving-disc comparison on a chosen local interval;
- local exact-root branch continuation near a base point `p0`.
- local quantitative branch comparison on an interval where the frozen
  simple-root lower bound and determinant error indicator are controlled.
- local constructive comparison from matrix-level remainder to
  determinant-level discrepancy on the moving discs.
- local moving-disc closure under the constructive branch bound.
- local first-order root-shift comparison under explicit derivative lower
  bounds and derivative regularity hypotheses.
- local derivative control from reduced-model objects on the moving discs.
- local denominator control from frozen-root noncoalescence or second-root
  exclusion, when a second frozen root branch is available.
- local discriminant persistence and denominator control from a nonvanishing
  frozen discriminant at a base point.
- local discriminant persistence from explicit coefficient control on
  `T(p)=tr K0(p)` and `D0(p)=det K0(p)`.

## What is only accepted standard background

Imported but not reproved in the package:

- simple-root corank-one consequences for finite-dimensional matrix problems;
- local continuation of simple roots;
- local right/left kernel representative selection;
- local regular basis completion after shrinking;
- Schur determinant formula for invertible `D`;
- analytic matrix difference-quotient decomposition;
- Rouche-type one-variable analytic root comparison;
- simple-root branch continuation by an implicit-function-style theorem.
- local factorization of a simple analytic root and its equivalence with a
  derivative lower bound after shrinking.
- the bound `|det E| <= (1/2)||E||_F^2` for `2x2` matrices.
- one-variable Taylor expansion with integral remainder.
- Frobenius trace inequalities for `2x2` matrices.
- local labeling of both roots of a simple quadratic family away from
  coalescence.
- local analytic square-root selection for a nonvanishing discriminant after
  shrinking.
- continuity or mean-value shrinking arguments for preserving a nonzero scalar
  quantity on a sufficiently small interval.

## What is only approximate

Approximate, not exact:

- the frozen `p`-only reduced model

  ```text
  ((Lambda-Lambda0)I-K0(p));
  ```

- any spectral reading of that frozen model before extra fixed-`p` boundary
  hypotheses are imposed;
- any quantitative identification of the exact-root branch with the frozen
  root branch beyond the conditional determinant-error estimates recorded
  here;
- any first-order shift formula beyond the conditional remainder-controlled
  comparison recorded here;
- any interpretation of that frozen model as symmetric, canonical, or already
  in `delta-kappa` form.

## Detailed status for the new parameterwise steps

### PMRC. Uniform parameterwise moving-disc root comparison

- `proved_here`:
  once a frozen-model root branch, a uniform radius, and uniform boundary
  domination are supplied, the same Rouche-type comparison gives one exact
  normalized root in each moving disc.
- `accepted_standard_background`:
  one-variable Rouche-type equality of zero count inside each disc.
- `explicit_hypothesis`:
  existence of `Lambda_fr(p)`, one uniform radius `r`, one disc family
  inside the claimed normalization region, simple frozen root in each disc,
  and uniform boundary domination.
- `what_it_still_does_not_give`:
  regular exact branch in `p`, quantitative root-shift bounds, symmetric
  normal form, `delta-kappa`, or a veering criterion.

### ERBC. Exact-root branch continuation

- `proved_here`:
  once the moving-disc step has selected a unique exact root at each `p`, a
  simple-root continuation theorem yields a local exact branch near `p0`,
  and disc uniqueness identifies that branch with the selected moving-disc
  root.
- `accepted_standard_background`:
  implicit-function-style continuation of a simple analytic root.
- `explicit_hypothesis`:
  joint analyticity near `(Lambda_ex(p0),p0)` and nonvanishing
  `partial_Lambda f_{p0}(Lambda_ex(p0))`.
- `what_it_still_does_not_give`:
  global continuation, quantitative exact-vs-frozen distance control,
  symmetric normal form, `delta-kappa`, or a veering criterion.

### QPBC. Quantitative parameterwise branch comparison

- `proved_here`:
  once the exact branch `Lambda_ex(p)` is already available and the frozen
  determinant has the factorization

  ```text
  g_p(Lambda) = a_p(Lambda)(Lambda-Lambda_fr(p)),
  |a_p(Lambda)| >= m0 > 0
  ```

  on the moving discs, evaluating `f_p(Lambda_ex(p))=0` yields

  ```text
  |Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p)/m0.
  ```
- `accepted_standard_background`:
  simple-root factorization of an analytic determinant and the passage from a
  derivative lower bound to a nonvanishing factor on a sufficiently small
  neighborhood.
- `explicit_hypothesis`:
  existence of the exact branch `Lambda_ex(p)`, moving-disc isolation, a
  uniform frozen simple-root lower bound `m0`, and a controlled determinant
  error indicator `eta_loc(p)`.
- `what_it_still_does_not_give`:
  derivative-level comparison formulas, quantitative gap laws, symmetric
  normal form, `delta-kappa`, or a final veering criterion.

### DDMR. Determinant-discrepancy bound from matrix remainder

- `proved_here`:
  in the local `2x2` setting, the exact identity

  ```text
  det(G+E)-det(G)=tr(adj(G)E)+det(E)
  ```

  yields the scalar bound

  ```text
  |f_p-g_p| <= ||G||_F ||E||_F + (1/2)||E||_F^2,
  ```

  and hence

  ```text
  eta_loc(p) <= M_G(p) M_E(p) + (1/2) M_E(p)^2.
  ```

  Combined with `QPBC`, this gives a constructive conditional bound for the
  branch distance.
- `accepted_standard_background`:
  the Frobenius-norm estimate `|det E| <= (1/2)||E||_F^2`.
- `explicit_hypothesis`:
  the exact normalized block and frozen model are already defined on the
  moving discs, and the local matrix controls `M_G(p), M_E(p)` are finite.
- `what_it_still_does_not_give`:
  derivative laws, asymptotic shift formulas, symmetric normal form,
  `delta-kappa`, or a final veering criterion.

### MDC. Moving-disc closure

- `proved_here`:
  once `B_shift(p)` is available, the implication

  ```text
  B_shift(p) < r => |Lambda_ex(p)-Lambda_fr(p)| < r
  ```

  is immediate, so the current quantitative bound certifies that the exact
  branch stays inside the chosen moving disc.
- `accepted_standard_background`:
  none beyond elementary inequality use.
- `explicit_hypothesis`:
  the constructive branch bound is already available and the radius `r` is
  the one defining the moving-disc family.
- `what_it_still_does_not_give`:
  a new nearby-root theorem, a derivative law, symmetric normal form,
  `delta-kappa`, or a final veering criterion.

### FOSF. First-order root-shift formula

- `proved_here`:
  fixing `p`, expanding `f_p` at `Lambda_fr(p)`, and comparing the linear term
  with `partial_Lambda g_p(Lambda_fr(p))` yields the first-order shift
  identity with an explicit remainder.
- `accepted_standard_background`:
  one-variable Taylor expansion with integral remainder.
- `explicit_hypothesis`:
  existence of the exact branch, moving-disc closure, a positive lower bound
  for `|partial_Lambda g_p(Lambda_fr(p))|`, and local bounds for
  `partial_Lambda(f_p-g_p)` and `partial_Lambda^2 g_p` on the moving disc.
- `what_it_still_does_not_give`:
  a final asymptotic shift law, constructive derivative controls from the
  matrix model, symmetric normal form, `delta-kappa`, or a final veering
  criterion.

### DCL. Derivative control from reduced-model objects

- `proved_here`:
  in the local `2x2` setting, the frozen derivative quantities satisfy

  ```text
  partial_Lambda g_p(Lambda)=tr G(Lambda,p),
  partial_Lambda^2 g_p(Lambda)=2,
  ```

  and the derivative discrepancy obeys

  ```text
  sup_{D_r(p)} |partial_Lambda(f_p-g_p)|
  <=
  sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p).
  ```
- `accepted_standard_background`:
  Frobenius trace inequalities for `2x2` matrices.
- `explicit_hypothesis`:
  the frozen model and matrix remainder are already defined on the moving
  discs, and `M_G(p), M_E(p), M_{E,1}(p)` are finite there.
- `what_it_still_does_not_give`:
  a uniform positive lower bound for `|partial_Lambda g_p(Lambda_fr(p))|`,
  a final asymptotic shift law, symmetric normal form, `delta-kappa`, or a
  final veering criterion.

### NCDL. Frozen-model denominator control from noncoalescence

- `proved_here`:
  once the frozen quadratic is written in terms of two local frozen root
  branches,

  ```text
  g_p(Lambda)
  =
  (Lambda-Lambda_fr(p))(Lambda-Lambda_fr,2(p)),
  ```

  differentiating at `Lambda_fr(p)` gives

  ```text
  partial_Lambda g_p(Lambda_fr(p))
  =
  Lambda_fr(p)-Lambda_fr,2(p),
  ```

  so any lower bound on frozen-root separation gives the denominator lower
  bound used in `FOSF`.
- `accepted_standard_background`:
  local labeling of the second frozen root branch away from coalescence.
- `explicit_hypothesis`:
  the second local frozen root branch exists on the interval under
  discussion, and either a separation lower bound or second-root exclusion
  from the moving disc is assumed there.
- `what_it_still_does_not_give`:
  a proof that the second branch exists for the current CoupledBeams regime,
  a proof that frozen-root separation stays uniformly positive, a final
  asymptotic shift law, symmetry, `delta-kappa`, or a final veering
  criterion.

### FDPL. Frozen discriminant persistence and denominator control

- `proved_here`:
  the frozen determinant is the monic quadratic

  ```text
  g_p(Lambda)
  =
  (Lambda-Lambda0)^2 - T(p) (Lambda-Lambda0) + D0(p),
  ```

  so distinct frozen roots are equivalent to

  ```text
  Disc_fr(p) = T(p)^2 - 4 D0(p) != 0.
  ```

  After shrinking and locally labeling the frozen roots, one has

  ```text
  Disc_fr(p) = (Lambda_fr(p)-Lambda_fr,2(p))^2,
  |partial_Lambda g_p(Lambda_fr(p))| = |Disc_fr(p)|^(1/2).
  ```

  Therefore a nonzero frozen discriminant at `p0` yields, on a smaller
  interval `I0`, a constant `c0 > 0` such that

  ```text
  |partial_Lambda g_p(Lambda_fr(p))| >= c0
  ```

  for all `p in I0`.
- `accepted_standard_background`:
  local labeling of the two frozen roots, or equivalently local choice of a
  square root of a nonvanishing discriminant, after shrinking.
- `explicit_hypothesis`:
  `K0(p)` is defined on a local interval and `Disc_fr(p0) != 0` at the base
  point.
- `what_it_still_does_not_give`:
  a proof that the discriminant stays uniformly far from zero on the actual
  CoupledBeams regime of interest, a final asymptotic shift law, symmetry,
  `delta-kappa`, or a final veering criterion.

### FDCC. Frozen discriminant persistence from coefficient control

- `proved_here`:
  with

  ```text
  T(p)=tr K0(p),
  D0(p)=det K0(p),
  Disc_fr(p)=T(p)^2-4D0(p),
  ```

  the algebraic difference identity gives

  ```text
  |Disc_fr(p)-Disc_fr(p0)|
  <=
  2 M_T eps_T + 4 eps_D
  ```

  under the stated coefficient bounds. Therefore the explicit margin

  ```text
  |Disc_fr(p0)| > 2 M_T eps_T + 4 eps_D
  ```

  implies `Disc_fr(p) != 0` on the interval, and then `FDPL` yields the local
  denominator lower bound after shrinking and local labeling.
- `accepted_standard_background`:
  continuity-based shrinking, and the derivative-based variant via a mean-value
  estimate when bounds on `T'(p)` and `D0'(p)` are used instead.
- `explicit_hypothesis`:
  local coefficient bounds on `T(p)` and `D0(p)`, together with base-point
  nonvanishing of `Disc_fr(p0)` by a margin exceeding the explicit error
  bound.
- `needs_caution`:
  this does not yet prove project-specific frozen noncoalescence on the actual
  CoupledBeams regime; it only turns that question into a transparent
  coefficient-control check.
- `what_it_still_does_not_give`:
  a project-specific regime proof, a final asymptotic shift law, symmetry,
  `delta-kappa`, or a final veering criterion.

## Current next real obstruction

The next real obstruction is:

```text
turn the first-order shift formula into a sharper asymptotic or derivative-
level law by turning the coefficient-level frozen discriminant checks into an
actual project-specific regime verification, without skipping prematurely to
symmetric normal form, delta-kappa, or a final veering criterion.
```
