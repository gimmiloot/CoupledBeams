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
| RUFOSF | Ready-to-use local first-order shift corollary | If branch selection/continuation and moving-disc closure are already in force, `DDMR` supplies the constructive branch bound, `DCL` supplies derivative-discrepancy control through `M_G,M_E,M_{E,1}`, and `FDPL` or `FDCC` supplies a frozen-discriminant denominator lower bound through `T,D0,Disc_fr`, then the `FOSF` first-order shift formula applies with a constructive remainder bound and without separate abstract denominator or derivative-discrepancy hypotheses. | `proved_here` | PMRC; ERBC; QPBC; DDMR; MDC; DCL; FDPL or FDCC | a ready-to-use conditional local first-order shift estimate for `Lambda_ex(p)-Lambda_fr(p)` | project-specific coefficient bounds, a sharper asymptotic or derivative-level law, symmetric normal form, `delta-kappa`, final veering criterion, branch-case applications |
| PDC | Parameterwise derivative comparison | If the exact and frozen branch equations can be differentiated on the local interval, frozen discriminant persistence gives `|partial_Lambda g_p(Lambda_fr(p))| >= c0`, the exact denominator is protected by `C_Lambda(p)+2B_shift(p)<c0`, and the `p`-derivative determinant discrepancy `H_pDelta(p)` is controlled on the moving discs, then differentiating `f_p(Lambda_ex(p))=0` and `g_p(Lambda_fr(p))=0` gives exact branch derivative formulas and a conditional bound for `d/dp(Lambda_ex-Lambda_fr)`. | `proved_here` | ERBC; RUFOSF; DCL; FDPL or FDCC; explicit `p`-derivative discrepancy control | local exact/frozen branch derivative identities and a conditional derivative-difference estimate | project-specific `p`-derivative discrepancy bounds, an asymptotic branch-shift law, symmetric normal form, `delta-kappa`, final veering criterion, branch-case applications |
| PDCL | `p`-derivative discrepancy control from reduced-model objects | In the local `2x2` setting, differentiating the exact determinant discrepancy identity gives `partial_p(f_p-g_p)=tr(adj(partial_p G)E)+tr(adj(G)partial_pE)+tr(adj(E)partial_pE)`, hence `H_pDelta(p) <= M_{G,p}(p) M_E(p) + (M_G(p)+M_E(p)) M_{E,p}(p)` on the moving discs. | `proved_here` | exact normalized block; frozen model; matrix remainder `E=F-G`; local controls `M_G,M_E,M_{G,p},M_{E,p}` | a constructive reduced-model bound for the `H_pDelta` input in `PDC` | project-specific smallness of the `p`-derivative matrix controls, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion |
| RUPDC | Ready-to-use parameterwise derivative comparison | Combining `PDC` with `PDCL` replaces the abstract `H_pDelta` input by `M_{G,p} M_E + (M_G+M_E)M_{E,p}` in the derivative-difference bound. | `proved_here` | PDC; PDCL; RUFOSF; DCL; FDPL or FDCC | a local conditional derivative-comparison estimate stated directly in reduced-model matrix controls | project-specific bounds for `M_{G,p},M_{E,p},T',D0'`, an asymptotic branch-shift law, symmetry, `delta-kappa`, final veering criterion |
| FCDL | Frozen coefficient-derivative control | In the frozen local `2x2` model, `partial_p G=-K0'`, hence `M_{G,p}=||K0'||_F`, while `T'=tr K0'` and `D0'=tr(adj(K0)K0')`; Frobenius bounds give `|T'| <= sqrt(2)||K0'||_F` and `|D0'| <= ||K0||_F||K0'||_F`. | `proved_here` | frozen `2x2` model; differentiable `K0(p)`; Frobenius norm | coefficient-level control of the frozen derivative inputs in `RUPDC` | project-specific bounds for `||K0'||_F` or `M_{E,p}`, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion |
| RUPDC-FC | Frozen-coefficient ready-to-use derivative comparison | Combining `RUPDC` with `FCDL` substitutes the frozen coefficient bounds for `M_{G,p}`, `T'`, and `D0'`, leaving the derivative-difference estimate in terms of `K0`, `K0'`, `M_E`, `M_{E,p}`, the existing matrix controls, and `c0`. | `proved_here` | RUPDC; FCDL; RUFOSF; DCL; FDPL or FDCC | a more directly usable local conditional derivative-comparison estimate for the frozen part of the reduced model | project-specific control of `K0'`, `M_{E,p}`, existing matrix controls, frozen-discriminant margins, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion |
| ERPC | Exact-remainder `p`-derivative control | Since `E=-(K(Lambda,p)-K0(p))`, one has `partial_p E=-(partial_p K(Lambda,p)-partial_p K(Lambda0,p))`; if the mixed derivative `partial_{Lambda p}K` is controlled on the segments from `Lambda0` to the moving disc, then `M_{E,p}(p) <= rho_0(p) M_{K,Lambda p}(p)`. | `proved_here` | exact normalized block; frozen model; differentiable `K`; segment inclusion; mixed derivative control | replaces the abstract exact-remainder input `M_{E,p}` by a structural normalized-block mixed-derivative control | project-specific bounds for `M_{K,Lambda p}`, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion |
| RUPDC-ER | Exact-remainder ready-to-use derivative comparison | Combining `RUPDC-FC` with `ERPC` substitutes `M_{E,p} <= rho_0 M_{K,Lambda p}` into the derivative-difference bound. | `proved_here` | RUPDC-FC; ERPC; RUFOSF; DCL; FDPL or FDCC | a derivative-comparison estimate stated in terms of `K0`, `K0'`, the mixed derivative of `K`, existing matrix controls, and `c0` | project-specific control of `K0'`, `M_{K,Lambda p}`, existing matrix controls, frozen-discriminant margins, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion |
| QDR | CoupledBeams quantitative-regime derivative comparison | If the local regime template `Q_der(I_Q)` supplies uniform bounds for `K0'`, `M_{K,Lambda p}`, `M_G`, `M_E`, `M_{E,1}`, `K0`, the frozen branch offset, the moving-disc closure, and the frozen denominator margin, then the `RUPDC-ER` estimate becomes directly usable with one set of regime constants. | `proved_here` | RUPDC-ER; local regime `Q_der(I_Q)` | a single checkable CoupledBeams regime package under which the derivative-comparison estimate can be cited without restating every remaining bound | verification that `Q_der(I_Q)` holds for any actual branch, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion, branch-case applications |
| QVR | Branch-ready `Q_der` verification protocol | For one concrete candidate local regime, fill the verification record: candidate branch/packet/window data, local construction checks, frozen-model margins, matrix and derivative bounds, and a final pass/fail decision. If all checks pass, then `Q_der(I_Q)` is verified for that candidate and `QDR` applies; otherwise no `QDR` conclusion may be drawn. | `proved_here` | QDR; branch-ready verification record | a practical checklist for verifying `Q_der(I_Q)` without pretending it has already been verified | actual verification for any candidate, packet construction, asymptotic derivative law, symmetry, `delta-kappa`, final veering criterion, branch-case theorem |

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
- the ready-to-use local first-order shift formula once the earlier branch
  selection/closure, determinant-discrepancy, derivative-control, and frozen
  discriminant-persistence hypotheses are all in force.
- the exact implicit-differentiation identities for the exact and frozen
  branches once the branch equations are differentiable and the relevant
  denominators are nonzero.
- the exact `p`-derivative determinant-discrepancy identity in the local
  `2x2` setting once `G` and `E` are `p`-differentiable.
- the exact frozen coefficient-derivative identities

  ```text
  partial_p G(Lambda,p) = -K0'(p),
  T'(p) = tr K0'(p),
  D0'(p) = tr(adj(K0(p)) K0'(p))
  ```

  and the resulting Frobenius bounds for `M_{G,p}`, `T'`, and `D0'` in the
  local `2x2` frozen model.
- the exact-remainder identity

  ```text
  E(Lambda,p) = -(K(Lambda,p)-K0(p))
  ```

  and its `p`-derivative identity

  ```text
  partial_p E(Lambda,p)
  =
  -(partial_p K(Lambda,p)-partial_p K(Lambda0,p)).
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
- local ready-to-use first-order shift control under the earlier branch
  selection/closure, constructive matrix-control, and frozen discriminant
  persistence hypotheses.
- local parameterwise derivative comparison under explicit denominator
  closure and `p`-derivative discrepancy hypotheses.
- local `p`-derivative determinant-discrepancy control under explicit
  `p`-differentiability and finite matrix-control hypotheses.
- ready-to-use parameterwise derivative comparison after substituting the
  constructive `p`-derivative discrepancy bound into `PDC`.
- local frozen coefficient-derivative control under differentiability of
  `K0(p)`.
- frozen-coefficient ready-to-use derivative comparison after substituting
  the coefficient-level bounds into `RUPDC`.
- local exact-remainder `p`-derivative control under differentiability of
  `K`, segment inclusion from `Lambda0` to the moving discs, and finite mixed
  derivative control `M_{K,Lambda p}`.
- exact-remainder ready-to-use derivative comparison after substituting
  `M_{E,p} <= rho_0 M_{K,Lambda p}` into `RUPDC-FC`.
- CoupledBeams quantitative-regime derivative comparison after packaging the
  remaining `RUPDC-ER` inputs into the explicit local regime `Q_der(I_Q)`.
- branch-ready `Q_der` verification protocol after turning the regime
  assumptions into a concrete candidate record with pass/fail fields.

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
- Cauchy-Schwarz estimates for traces in the Frobenius inner product.
- local labeling of both roots of a simple quadratic family away from
  coalescence.
- local analytic square-root selection for a nonvanishing discriminant after
  shrinking.
- continuity or mean-value shrinking arguments for preserving a nonzero scalar
  quantity on a sufficiently small interval.
- chain-rule / implicit-differentiation formulas for regular scalar branch
  equations.
- the determinant derivative formula and the linearity of the adjugate map
  for `2x2` matrices.
- the fundamental theorem of calculus along straight line segments for
  matrix-valued functions.

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
- the ready-to-use first-order shift corollary is still conditional and local,
  and its remainder must still be checked against the intended scale;
- the parameterwise derivative comparison bound is still conditional on
  `p`-derivative discrepancy controls and is not an asymptotic law by itself;
- the ready-to-use derivative-comparison corollary is still conditional on
  local `p`-derivative matrix controls and coefficient-derivative controls;
- the frozen-coefficient ready-to-use derivative-comparison corollary is
  still conditional on usable bounds for `K0'`, `M_{E,p}`, the existing matrix
  controls, and the frozen-discriminant margin;
- the exact-remainder ready-to-use derivative-comparison corollary is still
  conditional on usable bounds for `K0'`, `M_{K,Lambda p}`, the existing
  matrix controls, and the frozen-discriminant margin;
- the CoupledBeams quantitative-regime derivative-comparison corollary is
  still conditional on verifying the local regime `Q_der(I_Q)` for the actual
  branch and interval under study;
- the branch-ready verification protocol is still only a protocol; it does
  not verify any candidate unless the record is actually filled and every
  required check passes;
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

### RUFOSF. Ready-to-use local first-order shift corollary

- `proved_here`:
  this step only composes earlier results. `DDMR` and `MDC` supply the
  constructive branch bound

  ```text
  |Lambda_ex(p)-Lambda_fr(p)| <= B_shift(p),
  ```

  `DCL` supplies

  ```text
  H_Delta(p)
  <=
  sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p),
  H_g2(p) = 2,
  ```

  and `FDPL` or `FDCC` supplies a local constant `c0 > 0` such that

  ```text
  |partial_Lambda g_p(Lambda_fr(p))| >= c0.
  ```

  Substituting these into `FOSF` gives the ready-to-use estimate with
  constructive remainder

  ```text
  |Rem_shift(p)|
  <=
  [
    (sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p)) B_shift(p)
    + B_shift(p)^2
  ] / c0.
  ```

  The final `B_shift(p)^2` term is exactly the `(1/2) H_g2(p) B_shift(p)^2`
  contribution from `FOSF` after using the `DCL` identity `H_g2(p)=2`.
- `accepted_standard_background`:
  only the same background already used by the component steps: Rouche-type
  selection, simple-root continuation, Taylor expansion, Frobenius trace
  inequalities, and local labeling or persistence of noncoalesced frozen
  roots.
- `explicit_hypothesis`:
  the earlier branch-selection and closure hypotheses, finiteness of the
  matrix controls `M_G,M_E,M_{E,1}`, and either a frozen discriminant
  persistence input or the coefficient-control margin on `T(p)` and `D0(p)`.
- `needs_caution`:
  the corollary removes two abstract inputs from `FOSF`, but it still depends
  on local constructive bounds and on a frozen discriminant margin that must
  be verified in the intended regime.
- `what_it_still_does_not_give`:
  a project-specific coefficient-bound proof, a final asymptotic or
  derivative-level branch-shift law, symmetry, `delta-kappa`, a final veering
  criterion, or branch-case applications.

### PDC. Parameterwise derivative comparison

- `proved_here`:
  differentiating the already selected branch equations gives

  ```text
  Lambda_ex'(p)
  =
  - partial_p f_p(Lambda_ex(p))
    / partial_Lambda f_p(Lambda_ex(p)),

  Lambda_fr'(p)
  =
  - partial_p g_p(Lambda_fr(p))
    / partial_Lambda g_p(Lambda_fr(p)).
  ```

  For the frozen quadratic this becomes

  ```text
  Lambda_fr'(p)
  =
  [T'(p)(Lambda_fr(p)-Lambda0)-D0'(p)]
  /
  [2(Lambda_fr(p)-Lambda0)-T(p)].
  ```

  The comparison identity for
  `Lambda_ex'(p)-Lambda_fr'(p)` is algebraic substitution of these two
  formulas. The displayed bound follows by using `RUFOSF` for
  `|Lambda_ex-Lambda_fr| <= B_shift`, `DCL` for the
  `partial_Lambda(f-g)` control and `partial_Lambda^2 g=2`, and the explicit
  `H_pDelta` hypothesis for `partial_p(f-g)`.
- `accepted_standard_background`:
  chain-rule / implicit differentiation for regular scalar branch equations.
- `explicit_hypothesis`:
  differentiability in `(Lambda,p)`, the frozen denominator lower bound `c0`,
  the exact-denominator closure condition

  ```text
  C_Lambda(p) + 2 B_shift(p) < c0,
  ```

  and a local `p`-derivative determinant-discrepancy control

  ```text
  H_pDelta(p)
  =
  sup_{Lambda in D_r(p)} |partial_p(f_p-g_p)(Lambda)| < infinity.
  ```
- `needs_caution`:
  this is the first derivative-comparison step, but its useful size depends
  on controlling `H_pDelta(p)`, `T'(p)`, and `D0'(p)` in the intended local
  regime.
- `what_it_still_does_not_give`:
  a project-specific proof of those `p`-derivative controls, an asymptotic
  branch-shift law, symmetry, `delta-kappa`, a final veering criterion, or
  branch-case applications.

### PDCL. `p`-derivative discrepancy control from reduced-model objects

- `proved_here`:
  in the local `2x2` setting, start from the exact identity

  ```text
  f_p-g_p = tr(adj(G)E) + det E.
  ```

  Since the `2x2` adjugate map is linear, differentiating at fixed `Lambda`
  gives

  ```text
  partial_p(f_p-g_p)
  =
  tr(adj(partial_p G)E)
  + tr(adj(G)partial_p E)
  + tr(adj(E)partial_p E).
  ```

  Frobenius trace bounds and `||adj(A)||_F=||A||_F` for `2x2` matrices then
  give

  ```text
  H_pDelta(p)
  <=
  M_{G,p}(p) M_E(p)
  + (M_G(p)+M_E(p)) M_{E,p}(p).
  ```
- `accepted_standard_background`:
  determinant derivative formula, linearity of the `2x2` adjugate map, and
  Frobenius trace inequalities.
- `explicit_hypothesis`:
  `G` and `E` are `p`-differentiable on the moving discs, and the local matrix
  controls `M_G,M_E,M_{G,p},M_{E,p}` are finite there.
- `needs_caution`:
  the estimate is constructive, but it is only useful if the new
  `p`-derivative matrix controls are bounded in the intended local regime.
- `what_it_still_does_not_give`:
  project-specific smallness of `M_{G,p}` or `M_{E,p}`, an asymptotic
  derivative law, symmetry, `delta-kappa`, a final veering criterion, or
  branch-case applications.

### RUPDC. Ready-to-use parameterwise derivative comparison

- `proved_here`:
  substituting `PDCL` into `PDC` replaces the abstract input

  ```text
  H_pDelta(p)
  ```

  by the reduced-model matrix-control expression

  ```text
  M_{G,p}(p) M_E(p)
  + (M_G(p)+M_E(p)) M_{E,p}(p).
  ```

  Thus the derivative-difference bound is now directly stated in terms of
  `M_G`, `M_E`, `M_{E,1}`, `M_{G,p}`, `M_{E,p}`, `T'`, `D0'`, `B_shift`, and
  the frozen-discriminant lower bound `c0`.
- `accepted_standard_background`:
  only the standard background already used by `PDC` and `PDCL`.
- `explicit_hypothesis`:
  all `PDC` hypotheses except the abstract `H_pDelta` bound, plus the finite
  `p`-derivative matrix controls from `PDCL`.
- `needs_caution`:
  this is genuinely more usable than `PDC`, but it still has to be paired with
  project-specific estimates for the new matrix controls and coefficient
  derivatives.
- `what_it_still_does_not_give`:
  an asymptotic or project-specific branch-shift derivative law, symmetry,
  `delta-kappa`, a final veering criterion, or branch-case applications.

### FCDL. Frozen coefficient-derivative control

- `proved_here`:
  in the frozen local `2x2` model

  ```text
  G(Lambda,p)=((Lambda-Lambda0)I-K0(p)),
  ```

  differentiating at fixed `Lambda` gives

  ```text
  partial_p G(Lambda,p) = -K0'(p),
  M_{G,p}(p)=||K0'(p)||_F.
  ```

  The coefficient derivatives are

  ```text
  T'(p)=tr K0'(p),
  D0'(p)=tr(adj(K0(p))K0'(p)).
  ```

  For `2x2` matrices, Frobenius Cauchy-Schwarz and
  `||adj(K0)||_F=||K0||_F` give

  ```text
  |T'(p)| <= sqrt(2)||K0'(p)||_F,
  |D0'(p)| <= ||K0(p)||_F ||K0'(p)||_F.
  ```
- `accepted_standard_background`:
  Frobenius Cauchy-Schwarz estimates, the `2x2` adjugate norm identity, and
  the standard determinant derivative formula.
- `explicit_hypothesis`:
  `K0(p)` is differentiable at the parameter value under discussion, and the
  coefficient norms are finite there.
- `needs_caution`:
  these bounds are intentionally coarse coefficient-level controls. They do
  not prove that `||K0'(p)||_F` is small in the CoupledBeams regime.
- `what_it_still_does_not_give`:
  control of the exact-remainder derivative `M_{E,p}`, an asymptotic
  derivative law, symmetry, `delta-kappa`, a final veering criterion, or
  branch-case applications.

### RUPDC-FC. Frozen-coefficient ready-to-use derivative comparison

- `proved_here`:
  substituting `FCDL` into `RUPDC` gives a derivative-difference bound whose
  frozen derivative inputs are controlled by `K0` and `K0'`. The substitutions
  are

  ```text
  M_{G,p}(p) = ||K0'(p)||_F,
  |T'(p)| <= sqrt(2)||K0'(p)||_F,
  ```

  and

  ```text
  |-T'(p)(alpha(p)-Lambda0)+D0'(p)|
  <=
  (sqrt(2)|alpha(p)-Lambda0|+||K0(p)||_F)||K0'(p)||_F.
  ```
- `accepted_standard_background`:
  only the standard background already used by `RUPDC` and `FCDL`.
- `explicit_hypothesis`:
  all `RUPDC` hypotheses, differentiability of `K0(p)`, and finite
  coefficient-level norms for `K0(p)` and `K0'(p)`.
- `needs_caution`:
  this is a useful substitution corollary, not a project-specific smallness
  theorem. The exact-remainder derivative control `M_{E,p}` remains explicit.
- `what_it_still_does_not_give`:
  project-specific bounds for `K0'` or `M_{E,p}`, an asymptotic
  derivative-level branch-shift law, symmetry, `delta-kappa`, a final veering
  criterion, or branch-case applications.

### ERPC. Exact-remainder `p`-derivative control

- `proved_here`:
  from the exact normalized and frozen blocks,

  ```text
  E(Lambda,p)
  =
  F(Lambda,p)-G(Lambda,p)
  =
  -(K(Lambda,p)-K0(p)).
  ```

  Therefore

  ```text
  partial_p E(Lambda,p)
  =
  -(partial_p K(Lambda,p)-partial_p K(Lambda0,p)).
  ```

  If the straight segments from `Lambda0` to the moving disc stay inside the
  normalization region, the integral identity

  ```text
  partial_p E(Lambda,p)
  =
  -(Lambda-Lambda0)
  integral_0^1
  partial_{Lambda p}K(Lambda0+t(Lambda-Lambda0),p) dt
  ```

  gives

  ```text
  M_{E,p}(p)
  <=
  rho_0(p) M_{K,Lambda p}(p),
  rho_0(p) = sup_{Lambda in D_r(p)} |Lambda-Lambda0|.
  ```
- `accepted_standard_background`:
  the fundamental theorem of calculus along a straight segment for
  matrix-valued functions.
- `explicit_hypothesis`:
  differentiability of `K` in `p`, existence and continuity of the mixed
  derivative `partial_{Lambda p}K` on the segment set, and inclusion of those
  segments in the local normalization region.
- `needs_caution`:
  the factor is `rho_0(p)`, not automatically the moving-disc radius `r`.
  Only in the special case where the moving disc is centered at `Lambda0` (or
  otherwise has `rho_0 <= r`) does the bound become an `r`-factor estimate.
- `what_it_still_does_not_give`:
  project-specific bounds for `M_{K,Lambda p}`, an asymptotic derivative law,
  symmetry, `delta-kappa`, a final veering criterion, or branch-case
  applications.

### RUPDC-ER. Exact-remainder ready-to-use derivative comparison

- `proved_here`:
  substituting `ERPC` into `RUPDC-FC` replaces

  ```text
  M_{E,p}(p)
  ```

  by

  ```text
  rho_0(p) M_{K,Lambda p}(p).
  ```

  Thus the derivative-difference bound no longer contains the abstract
  exact-remainder derivative input and is instead stated through `K0`,
  `K0'`, the mixed derivative of `K`, the existing matrix controls, and the
  frozen-discriminant lower bound.
- `accepted_standard_background`:
  only the standard background already used by `RUPDC-FC` and `ERPC`.
- `explicit_hypothesis`:
  all `RUPDC-FC` hypotheses plus the segment and mixed-derivative hypotheses
  from `ERPC`.
- `needs_caution`:
  this is a genuine structural substitution, but it still requires estimating
  `M_{K,Lambda p}` in the intended local CoupledBeams regime.
- `what_it_still_does_not_give`:
  project-specific bounds for `K0'` or `M_{K,Lambda p}`, an asymptotic
  derivative-level branch-shift law, symmetry, `delta-kappa`, a final veering
  criterion, or branch-case applications.

### QDR. CoupledBeams quantitative-regime derivative comparison

- `proved_here`:
  under the explicit local regime `Q_der(I_Q)`, the remaining inputs of
  `RUPDC-ER` are collected into one set of uniform constants:

  ```text
  G_*, E_*, E_{1,*}, K_{0,*}, K_{p,*}, K_{Lambda p,*}, A_*,
  m0, c0, r.
  ```

  With

  ```text
  B_Q = [G_* E_* + (1/2) E_*^2] / m0,
  C_Q = sqrt(2) E_* + (G_* + E_*) E_{1,*},
  D_Q = c0 - C_Q - 2 B_Q > 0,
  ```

  the derivative-comparison estimate from `RUPDC-ER` is available directly
  with these constants. This removes repeated restatement of the remaining
  bottleneck assumptions in later local regime arguments.
- `accepted_standard_background`:
  none beyond the background already used by `RUPDC-ER`; this is a
  substitution/packaging corollary.
- `explicit_hypothesis`:
  the local regime `Q_der(I_Q)`: canonical construction through `RUPDC-ER`,
  uniform bounds for `K0'`, `M_{K,Lambda p}`, `M_G`, `M_E`, `M_{E,1}`, `K0`,
  the frozen branch offset, the moving-disc closure, and a frozen denominator
  or discriminant margin.
- `needs_caution`:
  `QDR` does not verify `Q_der(I_Q)`. It only records that once those
  project-specific checks are supplied, the derivative-comparison machinery is
  ready to cite.
- `what_it_still_does_not_give`:
  a proof that any actual CoupledBeams branch satisfies the regime,
  an asymptotic derivative-level branch-shift law, symmetry, `delta-kappa`, a
  final veering criterion, or branch-case applications.

### QVR. Branch-ready `Q_der` verification protocol

- `proved_here`:
  a concrete candidate can be tested by filling a verification record with
  candidate selection data, local construction checks, frozen-model checks,
  matrix and derivative bounds, and a final pass/fail decision. If every
  field required by `Q_der(I_Q)` is present and every check passes, then
  `QDR` applies on that candidate interval.
- `accepted_standard_background`:
  none beyond `QDR`; this is a project-facing verification protocol.
- `explicit_hypothesis`:
  the proposed candidate branch or branch pair, parameter interval, base
  point, local spectral window, packet/construction data, frozen margin data,
  and all constants entering `Q_der(I_Q)`.
- `needs_caution`:
  an incomplete record is not partial verification. Missing packet data,
  missing complement regularity, missing derivative bounds, or a failed
  closure/margin check blocks `QDR`. The current C1 record now fixes a
  candidate-level base point, spectral-window target, and retained-pair
  target from repository evidence, but it is still not a verified
  `Q_der(I_Q)` instance.
- `what_it_still_does_not_give`:
  actual verification for a named candidate, packet construction, an
  asymptotic derivative-level law, symmetry, `delta-kappa`, final veering
  criterion, or branch-case theorem.

## Current next real obstruction

The next real obstruction is:

```text
turn the candidate-level C1 theorem-line object into a proved local retained
packet with certified spectral window, complement regularity, and passing
`Q_der` constants, or replace this route by a different reduction strategy.
```

At this boundary, the inputs already reduced are:

- `H_pDelta(p)`, replaced by reduced-model matrix controls in `PDCL`/`RUPDC`;
- `M_{G,p}(p)`, `T'(p)`, and `D0'(p)`, replaced by `K0(p)` and `K0'(p)` in
  `FCDL`/`RUPDC-FC`;
- `M_{E,p}(p)`, replaced by `rho_0(p) M_{K,Lambda p}(p)` in `ERPC`/`RUPDC-ER`.

The genuinely remaining independent bottlenecks are `K0'(p)`,
`M_{K,Lambda p}(p)`, the existing local matrix and moving-disc closure
constants, and the frozen-discriminant margin. The direction audit in
`proof_notes.md` shows that rewriting `K0'(p)` and `M_{K,Lambda p}(p)` through
`R^{-1}S` is exact algebra but not yet a reduction of abstraction. The
follow-up gauge audit shows that imposing `R_hat(Lambda0,p)=I` removes the
explicit `partial_p R_hat(Lambda0,p)` term from the base-point formula, but
only by hiding the same derivative input in `partial_p S_hat`, and it does not
remove the mixed `R_hat`-derivatives.

The new `QDR` corollary packages these remaining inputs into the local regime
`Q_der(I_Q)`, and `QVR` states how to check that regime on one concrete
candidate. The C1 record in `proof_notes.md` now fixes the candidate-level
construction target: `mu0 = 0.26624999999999999`, the sorted-position `2/3`
spectral target anchored by
`Lambda_pair = 5.4876154692852506` and
`Lambda_target = 6.2086633731319516`, and the retained-pair target
`bending_desc_01` / `bending_desc_04` on the documented `beta=15 deg`,
physical-radius-`5 mm`, `0.20 <= mu <= 0.35` window. This advances C1 from a
repo-grounded seed to a theorem-line candidate object, but current repo
evidence still leaves `Q_der(I_Q)` undecidable there: certified open spectral
endpoints, packet construction, complement regularity, Schur/normalization
construction, frozen margins, matrix constants, and derivative constants are
not yet available. Further progress needs those C1 construction and bound
checks, derivative control from a concrete normalization construction, or a
different reduction strategy, without skipping prematurely to symmetric
normal form, delta-kappa, a final veering criterion, or branch-case
application theorems.
