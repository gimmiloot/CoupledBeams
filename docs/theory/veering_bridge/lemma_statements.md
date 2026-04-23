# Lemma Statements

## Scope

This file collects clean statement-level versions of the current bridge lemmas.
Definitions and conventions are taken from [foundations.md](foundations.md).

## 1. Packet construction in the simple-root setting

Fix a regular local domain `Omega x I` and assume:

1. `Omega` is an isolated two-root cluster for `det M(Lambda,p)` on `I`;
2. both roots in `Omega` are simple for `p in I`;
3. the two roots admit locally regular right and left null-vector
   representatives;
4. the two right directions and the two left directions remain linearly
   independent;
5. the retained pairing is nondegenerate.

Then, after shrinking `I` if necessary, one obtains local rank-2 packets

```text
V(p), W(p) in C^{6x2}
```

associated with the isolated two-root cluster.

This is a simple-root conditional packet lemma only. It does not treat
defective roots, multiple roots, complements, or reduced-block normal forms.

## 2. Complement regularity lemma

Assume a local rank-2 packet `V(p),W(p)` has already been supplied on a
regular local domain. Choose complements `Z(p),U(p)` so that `[V Z]` and
`[W U]` are invertible local bases, and define

```text
D(Lambda,p) = U(p)^T M(Lambda,p) Z(p).
```

If

```text
det D(Lambda,p) != 0
```

on the claimed local region, then the complement is regular there and Schur
elimination is locally meaningful.

This statement separates algebraic basis completion from the extra
regular-complement hypothesis. Basis completion alone does not imply
regularity of `D`.

## 3. Schur root-capture lemma

Assume:

1. compatible right/left data `V,W` and complements `Z,U` are supplied on a
   regular local domain;
2. `[V Z]` and `[W U]` are invertible on the claimed local region;
3. with

   ```text
   M_tilde = [W U]^T M [V Z] = [A B; C D],
   ```

   one has `det D(Lambda,p) != 0` on that region.

Define the Schur block

```text
S(Lambda,p) = A(Lambda,p) - B(Lambda,p) D(Lambda,p)^(-1) C(Lambda,p).
```

Then `det M(Lambda,p)` and `det S(Lambda,p)` have the same local zeros on the
claimed region. For fixed `p`, the same one-variable zero multiplicities in
`Lambda` are preserved.

This is exact local root capture, but only on the region where the basis
factors and `det D` remain nonzero.

## 4. Lambda-normalization lemma

Fix a base point `(Lambda0,p0)` in the claimed local region and assume:

1. the Schur root-capture step has already supplied a local analytic reduced
   block `S(Lambda,p)`;
2. the spectral derivative is nonsingular at the base point:

   ```text
   det(partial_Lambda S(Lambda0,p0)) != 0.
   ```

Then, after shrinking the local region if necessary, there exist analytic
`2x2` matrices `R(Lambda,p)` and `K(Lambda,p)` such that

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p)),
det R(Lambda,p) != 0.
```

Consequently the exact normalized block

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p))
```

has the same local zeros as `S`, and therefore the same retained local zeros
as `M`.

## 5. First-order reduced model lemma

Assume the exact normalized block already exists on a local region:

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p)),
det R(Lambda,p) != 0,
```

and `K(Lambda,p)` is analytic in `Lambda` there.

Define the frozen `p`-only matrix

```text
K0(p) = K(Lambda0,p).
```

Then there exists an analytic matrix `K1(Lambda,p)` such that

```text
K(Lambda,p) = K0(p) + (Lambda-Lambda0) K1(Lambda,p),
```

and hence

```text
((Lambda-Lambda0) I - K(Lambda,p))
=
((Lambda-Lambda0) I - K0(p))
- (Lambda-Lambda0) K1(Lambda,p).
```

On a sufficiently small local region, the remainder is bounded by

```text
O(|Lambda-Lambda0|)
```

at matrix level.

This is a controlled first-order approximation statement, not an exact
spectral equivalence theorem.

## 6. Frozen-model root-comparison lemma

Fix one parameter value `p` and define

```text
f_p(Lambda) = det F(Lambda,p),
g_p(Lambda) = det G(Lambda,p),
```

where

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)),
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)).
```

Let `D` be a closed disc in the `Lambda`-plane inside the claimed local
normalization region. Assume:

1. `g_p` has exactly one zero in `D`, and that zero is simple;
2. `g_p` has no zeros on `partial D`;
3. on `partial D`,

   ```text
   sup_{partial D} |f_p - g_p|
   <
   inf_{partial D} |g_p|.
   ```

Then `f_p` has exactly one zero in `D`, counted with multiplicity. Hence the
exact normalized block has one nearby root in the same disc, unique and
simple in that disc.

This is a narrow fixed-`p` local comparison theorem. It is not a global
continuation statement in `p`.

## 7. Uniform parameterwise moving-disc root-comparison lemma

Let `I` be a local interval in the selected scalar sweep parameter `p`, and
assume:

1. there is a local frozen-model root branch `Lambda_fr(p)` on `I`;
2. one fixed radius `r > 0` is chosen for all `p in I`;
3. every disc

   ```text
   D_r(p) = { Lambda : |Lambda - Lambda_fr(p)| <= r }
   ```

   lies inside the claimed local normalization region;
4. for each `p in I`, the frozen determinant `g_p` has exactly one zero in
   `D_r(p)`, and that zero is simple;
5. for each `p in I`, `g_p` has no zeros on `partial D_r(p)`;
6. the uniform boundary-domination inequality

   ```text
   sup_{p in I} sup_{Lambda in partial D_r(p)} |f_p(Lambda) - g_p(Lambda)|
   <
   inf_{p in I} inf_{Lambda in partial D_r(p)} |g_p(Lambda)|
   ```

   holds.

Then, for each `p in I`, the exact normalized determinant `f_p` has exactly
one zero in `D_r(p)`, counted with multiplicity. Hence for each `p in I`, the
exact normalized block has one root in `D_r(p)`, unique and simple in that
disc.

This is a parameterwise family of local fixed-`p` comparison statements with
a uniform radius and uniform boundary control. It is not yet, by itself, a
regular exact root branch theorem.

## 8. Exact-root branch continuation lemma

Assume the hypotheses and conclusion of the uniform parameterwise moving-disc
root-comparison lemma. Fix one base parameter value `p0 in I`, and let
`Lambda_ex(p0)` be the unique exact normalized root selected in `D_r(p0)`.

Assume further that:

1. the exact normalized determinant is jointly analytic in `(Lambda,p)` on a
   neighborhood of `(Lambda_ex(p0),p0)`;
2. the selected exact root is simple:

   ```text
   partial_Lambda f_{p0}(Lambda_ex(p0)) != 0.
   ```

Then, after shrinking to a smaller local interval `I_ex` around `p0` if
necessary, there exists a local regular exact root branch

```text
Lambda_ex(p)
```

on `I_ex` such that

```text
f_p(Lambda_ex(p)) = 0
```

for `p in I_ex`.

Because each disc `D_r(p)` contains exactly one exact normalized root, this
local branch is identified with the moving-disc selected root for every
`p in I_ex`.

This is a local branch-continuation statement for the exact root. It does not
yet quantify the distance between `Lambda_ex(p)` and `Lambda_fr(p)`.

## 9. Quantitative parameterwise branch-comparison lemma

Assume the hypotheses and conclusion of the exact-root branch continuation
lemma on a local interval `I_ex`. Thus:

- `Lambda_fr(p)` is the local frozen-model root branch;
- `Lambda_ex(p)` is the local exact normalized root branch identified by the
  moving-disc selection;
- each moving disc

  ```text
  D_r(p) = { Lambda : |Lambda - Lambda_fr(p)| <= r }
  ```

  contains exactly one exact normalized root, namely `Lambda_ex(p)`.

Assume further that:

1. there exists a constant `m0 > 0` and analytic functions `a_p(Lambda)` on
   `D_r(p)` such that

   ```text
   g_p(Lambda) = a_p(Lambda) (Lambda - Lambda_fr(p)),
   |a_p(Lambda)| >= m0
   ```

   for all `Lambda in D_r(p)` and `p in I_ex`;
2. the local determinant error indicator

   ```text
   eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda) - g_p(Lambda)|
   ```

   is finite on `I_ex`.

Then, for every `p in I_ex`,

```text
|Lambda_ex(p) - Lambda_fr(p)| <= eta_loc(p) / m0.
```

Consequently, if

```text
eta_* = sup_{p in I_ex} eta_loc(p)
```

is finite, then

```text
sup_{p in I_ex} |Lambda_ex(p) - Lambda_fr(p)| <= eta_* / m0.
```

This is a local quantitative comparison theorem between the exact and frozen
branches. It is conditional on a uniform frozen simple-root lower bound and a
determinant-level error bound. It does not yet provide a sharper asymptotic
law, a derivative formula for the branch shift, or a veering criterion.

## 10. Determinant-discrepancy bound from matrix remainder lemma

Work on a local interval where the exact normalized block and frozen model are
already defined:

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)),
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
E(Lambda,p) = F(Lambda,p) - G(Lambda,p).
```

For one fixed `p`, define

```text
f_p(Lambda) = det F(Lambda,p),
g_p(Lambda) = det G(Lambda,p),
```

and use the Frobenius norm on `2x2` matrices.

Then for every `Lambda`,

```text
f_p(Lambda) - g_p(Lambda)
=
tr(adj(G(Lambda,p))) E(Lambda,p) + det(E(Lambda,p)),
```

and consequently

```text
|f_p(Lambda) - g_p(Lambda)|
<=
||G(Lambda,p)||_F ||E(Lambda,p)||_F + (1/2) ||E(Lambda,p)||_F^2.
```

If `D_r(p)` is a moving disc and

```text
M_G(p) = sup_{Lambda in D_r(p)} ||G(Lambda,p)||_F,
M_E(p) = sup_{Lambda in D_r(p)} ||E(Lambda,p)||_F,
```

then

```text
eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda) - g_p(Lambda)|
<=
M_G(p) M_E(p) + (1/2) M_E(p)^2.
```

Consequently, if the hypotheses of the quantitative parameterwise
branch-comparison lemma also hold, then

```text
|Lambda_ex(p) - Lambda_fr(p)|
<=
[M_G(p) M_E(p) + (1/2) M_E(p)^2] / m0.
```

This is a constructive local quantitative estimate derived from the matrix
remainder of the frozen model. It is still conditional on the earlier moving-
disc and frozen simple-root hypotheses, and it does not yet yield a
derivative-level or asymptotic shift law.

## 11. Moving-disc closure corollary

Assume the hypotheses of the determinant-discrepancy bound from matrix
remainder lemma and of the quantitative parameterwise branch-comparison lemma
on a local interval `I_ex`. Define

```text
B_shift(p) = [M_G(p) M_E(p) + (1/2) M_E(p)^2] / m0.
```

Then, for every `p in I_ex`,

```text
|Lambda_ex(p) - Lambda_fr(p)| <= B_shift(p).
```

Consequently:

1. if `B_shift(p) < r` for one fixed parameter value, then

   ```text
   Lambda_ex(p) in int D_r(p);
   ```
2. if

   ```text
   B_* = sup_{p in I_ex} B_shift(p) < r,
   ```

   then

   ```text
   Lambda_ex(p) in int D_r(p)
   ```

   for every `p in I_ex`.

This is a closure/self-consistency corollary for the chosen moving discs. It
does not add a new existence theorem; it shows that the current quantitative
bound is itself compatible with the disc family.

## 12. First-order root-shift formula lemma

Assume the hypotheses and conclusion of the exact-root branch continuation
lemma on a local interval `I_ex`, and assume the moving-disc closure
condition holds there so that

```text
Lambda_ex(p) in D_r(p)
```

for `p in I_ex`.

Fix one `p in I_ex`, and write

```text
delta_p = Lambda_ex(p) - Lambda_fr(p).
```

Assume further that:

1. the frozen determinant derivative is nonzero at the frozen root,

   ```text
   partial_Lambda g_p(Lambda_fr(p)) != 0;
   ```
2. there is a lower bound

   ```text
   |partial_Lambda g_p(Lambda_fr(p))| >= m1(p) > 0;
   ```
3. on the moving disc `D_r(p)`,

   ```text
   H_Delta(p) = sup_{Lambda in D_r(p)} |partial_Lambda(f_p-g_p)(Lambda)|
   ```

   is finite;
4. on the same disc,

   ```text
   H_g2(p) = sup_{Lambda in D_r(p)} |partial_Lambda^2 g_p(Lambda)|
   ```

   is finite.

Then there exists a remainder `Rem_shift(p)` such that

```text
Lambda_ex(p) - Lambda_fr(p)
=
- (f_p(Lambda_fr(p)) - g_p(Lambda_fr(p))) / partial_Lambda g_p(Lambda_fr(p))
+ Rem_shift(p),
```

and

```text
|Rem_shift(p)|
<=
[H_Delta(p) |delta_p| + (1/2) H_g2(p) |delta_p|^2] / m1(p).
```

Consequently, if the constructive bound `|delta_p| <= B_shift(p)` is already
available, then

```text
|Rem_shift(p)|
<=
[H_Delta(p) B_shift(p) + (1/2) H_g2(p) B_shift(p)^2] / m1(p).
```

This is a conditional local first-order comparison law between the exact and
frozen branches. It is not yet a final asymptotic law unless the remainder is
shown to be smaller than the leading term in a further argument.

## 13. Derivative-control lemma from reduced-model objects

Work on a local interval where the frozen model and the matrix remainder are
already defined:

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
E(Lambda,p) = F(Lambda,p) - G(Lambda,p).
```

Let `D_r(p)` be the current moving discs, and define

```text
M_G(p) = sup_{Lambda in D_r(p)} ||G(Lambda,p)||_F,
M_E(p) = sup_{Lambda in D_r(p)} ||E(Lambda,p)||_F,
M_{E,1}(p) = sup_{Lambda in D_r(p)} ||partial_Lambda E(Lambda,p)||_F.
```

Then, in the local `2x2` setting:

1. the frozen determinant derivative is given exactly by

   ```text
   partial_Lambda g_p(Lambda)
   =
   tr G(Lambda,p)
   =
   2 (Lambda-Lambda0) - tr K0(p);
   ```

   hence at the frozen root

   ```text
   |partial_Lambda g_p(Lambda_fr(p))|
   =
   |2 (Lambda_fr(p)-Lambda0) - tr K0(p)|.
   ```
2. the frozen second derivative is given exactly by

   ```text
   partial_Lambda^2 g_p(Lambda) = 2,
   ```

   so

   ```text
   sup_{Lambda in D_r(p)} |partial_Lambda^2 g_p(Lambda)| = 2;
   ```
3. the derivative discrepancy entering the first-order root-shift formula
   satisfies

   ```text
   sup_{Lambda in D_r(p)} |partial_Lambda(f_p-g_p)(Lambda)|
   <=
   sqrt(2) M_E(p) + (M_G(p) + M_E(p)) M_{E,1}(p).
   ```

Consequently, the derivative hypotheses used in the first-order root-shift
formula are reduced to:

- a lower bound on the explicit frozen-model scalar

  ```text
  |2 (Lambda_fr(p)-Lambda0) - tr K0(p)|;
  ```
- local matrix controls for `E` and `partial_Lambda E` on the moving discs.

This is a conditional constructive derivative-control lemma. It does not yet
prove that the derivative lower bound is uniform or separated away from zero;
that remains an explicit noncoalescence hypothesis on the frozen branch.

## 14. Frozen-model denominator-control lemma from noncoalescence

Work on a local interval where:

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p))
```

is the frozen `2x2` model, `Lambda_fr(p)` is the chosen local frozen root
branch, and a second local frozen root branch

```text
Lambda_fr,2(p)
```

is available.

Then, for each fixed `p` in that interval,

```text
g_p(Lambda)
=
det G(Lambda,p)
=
(Lambda - Lambda_fr(p)) (Lambda - Lambda_fr,2(p)).
```

Consequently,

```text
partial_Lambda g_p(Lambda_fr(p))
=
Lambda_fr(p) - Lambda_fr,2(p),
```

and therefore

```text
|partial_Lambda g_p(Lambda_fr(p))|
=
|Lambda_fr(p) - Lambda_fr,2(p)|
=
sep_fr(p).
```

Hence:

1. if the frozen roots satisfy the separation bound

   ```text
   |Lambda_fr(p) - Lambda_fr,2(p)| >= sigma_fr(p) > 0,
   ```

   then

   ```text
   |partial_Lambda g_p(Lambda_fr(p))| >= sigma_fr(p);
   ```
2. if a uniform frozen-root noncoalescence bound

   ```text
   |Lambda_fr(p) - Lambda_fr,2(p)| >= sigma_fr > 0
   ```

   holds on the interval, then

   ```text
   |partial_Lambda g_p(Lambda_fr(p))| >= sigma_fr
   ```

   there;
3. if the second frozen root stays outside the moving disc `D_r(p)`, then

   ```text
   |partial_Lambda g_p(Lambda_fr(p))| > r.
   ```

This is a conditional denominator-control lemma for the frozen model. It does
not prove that the second frozen root branch exists or that the frozen-root
separation bound holds for CoupledBeams; it shows that such a separation
hypothesis is exactly the denominator hypothesis needed in `FOSF`.

## 15. Frozen discriminant persistence and denominator-control lemma

Fix a base parameter value `p0` and work near it with the frozen `2x2` model

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)).
```

Define the frozen determinant and frozen discriminant by

```text
T(p) = tr K0(p),
D0(p) = det K0(p),

g_p(Lambda)
=
(Lambda-Lambda0)^2 - T(p) (Lambda-Lambda0) + D0(p),

Disc_fr(p) = T(p)^2 - 4 D0(p).
```

Then, for each fixed `p`:

1. the frozen roots are distinct if and only if

   ```text
   Disc_fr(p) != 0;
   ```
2. after shrinking to a sufficiently small interval on which two local frozen
   root branches `Lambda_fr(p), Lambda_fr,2(p)` can be labeled, one has

   ```text
   Disc_fr(p) = (Lambda_fr(p) - Lambda_fr,2(p))^2;
   ```
3. consequently,

   ```text
   |partial_Lambda g_p(Lambda_fr(p))|
   =
   |Lambda_fr(p) - Lambda_fr,2(p)|
   =
   |Disc_fr(p)|^(1/2).
   ```

Assume now that

```text
Disc_fr(p0) != 0.
```

Then, after shrinking to a smaller local interval `I0` around `p0` if
necessary:

1. the frozen discriminant remains nonzero on `I0`;
2. the two frozen roots remain distinct on `I0`;
3. there exist local frozen root branches `Lambda_fr(p), Lambda_fr,2(p)` on
   `I0`;
4. on `I0`, one has the local root-label identity

   ```text
   Disc_fr(p) = (Lambda_fr(p) - Lambda_fr,2(p))^2;
   ```
5. there exists a constant `c0 > 0` such that

   ```text
   |Disc_fr(p)| >= c0^2,
   |Lambda_fr(p) - Lambda_fr,2(p)| >= c0,
   |partial_Lambda g_p(Lambda_fr(p))| >= c0
   ```

   for all `p in I0`.

This is a local persistence lemma for frozen noncoalescence. It converts the
remaining denominator-control hypothesis into a coefficient-level frozen
discriminant condition near `p0`. It does not prove that the needed
discriminant lower bound holds in the actual CoupledBeams regime.

## 16. Frozen discriminant persistence from coefficient control lemma

Fix a base parameter value `p0` and a local interval `I` around it. For the
frozen `2x2` model

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
```

define

```text
T(p) = tr K0(p),
D0(p) = det K0(p),
Disc_fr(p) = T(p)^2 - 4 D0(p).
```

Assume that on `I` one has the coefficient-level bounds

```text
|T(p)| <= M_T,
|T(p) - T(p0)| <= eps_T,
|D0(p) - D0(p0)| <= eps_D.
```

Then, for every `p in I`,

```text
|Disc_fr(p) - Disc_fr(p0)|
<=
2 M_T eps_T + 4 eps_D.
```

Consequently, if

```text
|Disc_fr(p0)| > 2 M_T eps_T + 4 eps_D,
```

then

```text
Disc_fr(p) != 0
```

for all `p in I`.

Therefore, after possibly shrinking to a smaller local interval `I0` around
`p0` and locally labeling the two frozen roots there, the conclusions of the
frozen discriminant persistence and denominator-control lemma apply on `I0`.
In particular, there exists a constant `c0 > 0` such that

```text
|partial_Lambda g_p(Lambda_fr(p))| >= c0
```

for all `p in I0`.

This is a coefficient-level persistence lemma for frozen noncoalescence. It
does not require the second frozen root branch as a primary input. It does
not prove that the needed coefficient bounds hold on the actual CoupledBeams
regime.

## 17. Ready-to-use local first-order shift corollary

Work on a local interval `I0` on which the exact normalized block, frozen
model, matrix remainder, frozen root branch, and moving discs are already
defined:

```text
F(Lambda,p), G(Lambda,p), E(Lambda,p)=F(Lambda,p)-G(Lambda,p),
Lambda_fr(p), D_r(p).
```

Assume the previous branch-selection and closure steps have already supplied
a local exact branch `Lambda_ex(p)` on `I0` such that

```text
f_p(Lambda_ex(p)) = 0,
Lambda_ex(p) in D_r(p),
|Lambda_ex(p)-Lambda_fr(p)| <= B_shift(p).
```

Here `B_shift(p)` is the constructive branch-shift bound supplied by the
determinant-discrepancy and moving-disc closure steps.

Assume also that the determinant-discrepancy and derivative-control
quantities from the reduced-model objects are finite on the moving discs:

```text
M_G(p), M_E(p), M_{E,1}(p) < infinity.
```

Finally, assume the frozen discriminant persistence step has been verified on
`I0`, either directly through the nonvanishing of

```text
Disc_fr(p) = T(p)^2 - 4 D0(p),
T(p)=tr K0(p), D0(p)=det K0(p),
```

or through the coefficient-control criterion. Thus there is a constant
`c0 > 0` such that, for all `p in I0`,

```text
|partial_Lambda g_p(Lambda_fr(p))| >= c0.
```

Then, for every `p in I0`, the local first-order shift formula can be used in
the ready-to-use form

```text
Lambda_ex(p) - Lambda_fr(p)
=
- (f_p(Lambda_fr(p)) - g_p(Lambda_fr(p)))
  / partial_Lambda g_p(Lambda_fr(p))
+ Rem_shift(p),
```

with the constructive remainder bound

```text
|Rem_shift(p)|
<=
[
  (sqrt(2) M_E(p) + (M_G(p) + M_E(p)) M_{E,1}(p)) B_shift(p)
  + B_shift(p)^2
] / c0.
```

The coefficient `1` in front of `B_shift(p)^2` comes from the `FOSF`
remainder term

```text
(1/2) H_g2(p) |delta_p|^2
```

after using the frozen local `2x2` identity `H_g2(p)=2` from derivative
control.

In this corollary, the denominator lower bound is supplied by frozen
discriminant persistence through `Disc_fr(p)`, and the derivative-discrepancy
control is supplied by the reduced-model matrix controls
`M_G(p), M_E(p), M_{E,1}(p)`. No separate abstract denominator hypothesis or
separate abstract derivative-discrepancy hypothesis is being added here.

This is a local conditional corollary packaging the earlier steps. It does not
prove a project-specific discriminant margin, a sharper asymptotic law,
symmetric normal form, project-defined `delta-kappa`, a final veering
criterion, or any branch-case application theorem.

## 18. Parameterwise derivative comparison lemma

Work on a local interval `I0` where the ready-to-use local first-order shift
corollary applies. Write

```text
lambda(p) = Lambda_ex(p),
alpha(p) = Lambda_fr(p),
delta_p = lambda(p) - alpha(p).
```

Assume that `f(Lambda,p)=f_p(Lambda)` and `g(Lambda,p)=g_p(Lambda)` are
regular enough in `(Lambda,p)` to differentiate the two branch equations

```text
f_p(lambda(p)) = 0,
g_p(alpha(p)) = 0.
```

Assume also that the frozen denominator lower bound from frozen discriminant
persistence is available:

```text
|partial_Lambda g_p(alpha(p))| >= c0 > 0.
```

Let the derivative-control bound from reduced-model objects be

```text
C_Lambda(p)
=
sqrt(2) M_E(p) + (M_G(p) + M_E(p)) M_{E,1}(p).
```

Assume the exact-branch denominator is kept nonzero by the local closure
condition

```text
C_Lambda(p) + 2 B_shift(p) < c0.
```

Then

```text
|partial_Lambda f_p(lambda(p))|
>=
c0 - C_Lambda(p) - 2 B_shift(p)
> 0.
```

Consequently the exact branch satisfies

```text
lambda'(p)
=
- partial_p f_p(lambda(p)) / partial_Lambda f_p(lambda(p)).
```

The frozen branch satisfies

```text
alpha'(p)
=
- partial_p g_p(alpha(p)) / partial_Lambda g_p(alpha(p)).
```

Since

```text
g_p(Lambda)
=
(Lambda-Lambda0)^2 - T(p) (Lambda-Lambda0) + D0(p),
```

this frozen derivative formula is equivalently

```text
alpha'(p)
=
[T'(p) (alpha(p)-Lambda0) - D0'(p)]
/
[2 (alpha(p)-Lambda0) - T(p)].
```

The branch-shift derivative satisfies the exact comparison identity

```text
lambda'(p) - alpha'(p)
=
- [partial_p f_p(lambda(p)) - partial_p g_p(alpha(p))]
  / partial_Lambda f_p(lambda(p))
+
partial_p g_p(alpha(p))
  [partial_Lambda f_p(lambda(p)) - partial_Lambda g_p(alpha(p))]
  /
  [partial_Lambda f_p(lambda(p)) partial_Lambda g_p(alpha(p))].
```

If, in addition, the `p`-derivative determinant discrepancy is controlled on
the moving discs by

```text
H_pDelta(p)
=
sup_{Lambda in D_r(p)} |partial_p(f_p-g_p)(Lambda)|
< infinity,
```

then the derivative difference obeys the conditional bound

```text
|lambda'(p) - alpha'(p)|
<=
[H_pDelta(p) + |T'(p)| B_shift(p)]
/
[c0 - C_Lambda(p) - 2 B_shift(p)]
+
|-T'(p) (alpha(p)-Lambda0) + D0'(p)|
  [C_Lambda(p) + 2 B_shift(p)]
/
[
  (c0 - C_Lambda(p) - 2 B_shift(p)) c0
].
```

This is a local conditional derivative-comparison step. It proves the
implicit-differentiation identities and the displayed comparison bound under
explicit denominator and `p`-derivative discrepancy hypotheses. It does not
prove that `H_pDelta(p)` or the coefficient derivatives are small in the
actual CoupledBeams regime, and it does not yield an asymptotic branch-shift
law, symmetric normal form, project-defined `delta-kappa`, a final veering
criterion, or a branch-case application theorem.

## 19. `p`-derivative determinant-discrepancy control lemma

Work on a local interval where the exact normalized block, frozen model, and
matrix remainder are already defined:

```text
F(Lambda,p) = G(Lambda,p) + E(Lambda,p),
f_p(Lambda) = det F(Lambda,p),
g_p(Lambda) = det G(Lambda,p).
```

Assume that `G` and `E` are differentiable in `p` on the moving discs
`D_r(p)`. In the local `2x2` setting, differentiating the exact determinant
discrepancy identity gives, at fixed `Lambda`,

```text
partial_p(f_p-g_p)(Lambda)
=
tr(adj(partial_p G(Lambda,p)) E(Lambda,p))
+ tr(adj(G(Lambda,p)) partial_p E(Lambda,p))
+ tr(adj(E(Lambda,p)) partial_p E(Lambda,p)).
```

Consequently, using the Frobenius norm,

```text
|partial_p(f_p-g_p)(Lambda)|
<=
||partial_p G(Lambda,p)||_F ||E(Lambda,p)||_F
+
(||G(Lambda,p)||_F + ||E(Lambda,p)||_F)
||partial_p E(Lambda,p)||_F.
```

Define the local `p`-derivative matrix controls on `D_r(p)` by

```text
M_{G,p}(p) = sup_{Lambda in D_r(p)} ||partial_p G(Lambda,p)||_F,
M_{E,p}(p) = sup_{Lambda in D_r(p)} ||partial_p E(Lambda,p)||_F.
```

Then the `p`-derivative determinant discrepancy in `PDC` satisfies

```text
H_pDelta(p)
=
sup_{Lambda in D_r(p)} |partial_p(f_p-g_p)(Lambda)|
<=
M_{G,p}(p) M_E(p)
+ (M_G(p) + M_E(p)) M_{E,p}(p).
```

This is an exact local algebraic identity followed by a constructive
Frobenius-norm bound. It is conditional on the stated `p`-differentiability
and finiteness of the local matrix controls. It does not prove that those
controls are small in the actual CoupledBeams regime.

## 20. Ready-to-use parameterwise derivative-comparison corollary

Assume the hypotheses of the parameterwise derivative comparison lemma on a
local interval `I0`. Assume also that the `p`-derivative
determinant-discrepancy control lemma applies on the same moving discs.

With

```text
C_Lambda(p)
=
sqrt(2) M_E(p) + (M_G(p) + M_E(p)) M_{E,1}(p),
```

and the exact-denominator closure condition

```text
C_Lambda(p) + 2 B_shift(p) < c0,
```

the branch-shift derivative obeys the ready-to-use bound

```text
|lambda'(p) - alpha'(p)|
<=
[
  M_{G,p}(p) M_E(p)
  + (M_G(p) + M_E(p)) M_{E,p}(p)
  + |T'(p)| B_shift(p)
]
/
[c0 - C_Lambda(p) - 2 B_shift(p)]
+
|-T'(p) (alpha(p)-Lambda0) + D0'(p)|
  [C_Lambda(p) + 2 B_shift(p)]
/
[
  (c0 - C_Lambda(p) - 2 B_shift(p)) c0
].
```

Thus the abstract `H_pDelta(p)` input in `PDC` has been replaced by local
matrix controls for `G`, `E`, `partial_p G`, and `partial_p E` on the moving
discs.

This remains a local conditional corollary. It does not prove project-specific
smallness of `M_{G,p}(p)`, `M_{E,p}(p)`, `T'(p)`, or `D0'(p)`, and it does
not yield an asymptotic branch-shift law, symmetric normal form,
project-defined `delta-kappa`, a final veering criterion, or a branch-case
application theorem.

## 21. Frozen coefficient-derivative control lemma

Work in the frozen local `2x2` setting

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
T(p) = tr K0(p),
D0(p) = det K0(p),
```

and assume `K0(p)` is differentiable in `p`. Write

```text
K0(p) =
[ a(p)  b(p)
  c(p)  d(p) ].
```

Then the frozen matrix derivative is exact:

```text
partial_p G(Lambda,p) = -K0'(p).
```

Consequently, on every moving disc `D_r(p)`,

```text
M_{G,p}(p)
=
sup_{Lambda in D_r(p)} ||partial_p G(Lambda,p)||_F
=
||K0'(p)||_F.
```

The frozen coefficient derivatives are

```text
T'(p) = tr K0'(p) = a'(p) + d'(p),
```

and

```text
D0'(p)
=
tr(adj(K0(p)) K0'(p))
=
a'(p)d(p) + a(p)d'(p) - b'(p)c(p) - b(p)c'(p).
```

Using the Frobenius norm in the local `2x2` setting gives the coefficient
bounds

```text
|T'(p)| <= sqrt(2) ||K0'(p)||_F,
|D0'(p)| <= ||K0(p)||_F ||K0'(p)||_F.
```

For any frozen branch value `alpha(p)=Lambda_fr(p)`, the frozen numerator
appearing in the derivative-comparison estimate therefore satisfies

```text
|-T'(p)(alpha(p)-Lambda0) + D0'(p)|
<=
(sqrt(2) |alpha(p)-Lambda0| + ||K0(p)||_F) ||K0'(p)||_F.
```

This is exact coefficient algebra followed by Frobenius-norm bounds. It
removes `M_{G,p}(p)`, `T'(p)`, and `D0'(p)` as completely free inputs in the
frozen part of the derivative comparison, but it still requires
project-specific control of `K0'(p)` and of the exact-remainder derivative
control `M_{E,p}(p)`.

## 22. Frozen-coefficient ready-to-use derivative-comparison corollary

Assume the hypotheses of the ready-to-use parameterwise
derivative-comparison corollary. Assume also that the frozen
coefficient-derivative control lemma applies at the same parameter value.

Then, with

```text
C_Lambda(p)
=
sqrt(2) M_E(p) + (M_G(p) + M_E(p)) M_{E,1}(p),
```

and

```text
C_Lambda(p) + 2 B_shift(p) < c0,
```

the branch-shift derivative obeys the coefficient-level ready-to-use bound

```text
|lambda'(p) - alpha'(p)|
<=
[
  ||K0'(p)||_F M_E(p)
  + (M_G(p) + M_E(p)) M_{E,p}(p)
  + sqrt(2) ||K0'(p)||_F B_shift(p)
]
/
[c0 - C_Lambda(p) - 2 B_shift(p)]
+
[
  (sqrt(2) |alpha(p)-Lambda0| + ||K0(p)||_F) ||K0'(p)||_F
]
  [C_Lambda(p) + 2 B_shift(p)]
/
[
  (c0 - C_Lambda(p) - 2 B_shift(p)) c0
].
```

Thus the `RUPDC` inputs `M_{G,p}(p)`, `T'(p)`, and `D0'(p)` have been
replaced by coefficient-level frozen-model quantities involving only
`K0(p)` and `K0'(p)`, together with the already present exact-remainder
control `M_{E,p}(p)`.

This remains a local conditional corollary. The displayed estimate is not a
project-specific asymptotic law: it still needs useful bounds for
`||K0'(p)||_F`, `M_{E,p}(p)`, the existing matrix controls, and the frozen
discriminant margin in the intended CoupledBeams regime.

## 23. Exact-remainder `p`-derivative control lemma

Work on a local interval where the exact normalized block and frozen model are
defined:

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)),
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
K0(p)=K(Lambda0,p),
E(Lambda,p)=F(Lambda,p)-G(Lambda,p).
```

Then the exact-remainder identity is

```text
E(Lambda,p)
=
-(K(Lambda,p)-K0(p)).
```

Assume `K` is differentiable in `p` and `K0'(p)=partial_p K(Lambda0,p)`.
Then, at fixed `Lambda`,

```text
partial_p E(Lambda,p)
=
- (partial_p K(Lambda,p) - partial_p K(Lambda0,p)).
```

Assume further that, for every `Lambda in D_r(p)`, the straight segment

```text
Lambda0 + t (Lambda - Lambda0), 0 <= t <= 1,
```

lies in the local normalization region, and that the mixed derivative
`partial_{Lambda p} K` is continuous there. Define

```text
rho_0(p) = sup_{Lambda in D_r(p)} |Lambda - Lambda0|
```

and

```text
M_{K,Lambda p}(p)
=
sup_{Lambda in D_r(p), t in [0,1]}
||partial_{Lambda p} K(Lambda0 + t (Lambda - Lambda0), p)||_F.
```

Then

```text
partial_p E(Lambda,p)
=
- (Lambda-Lambda0)
  integral_0^1
  partial_{Lambda p} K(Lambda0 + t (Lambda-Lambda0), p) dt,
```

and consequently

```text
M_{E,p}(p)
=
sup_{Lambda in D_r(p)} ||partial_p E(Lambda,p)||_F
<=
rho_0(p) M_{K,Lambda p}(p).
```

Since

```text
rho_0(p) <= |Lambda_fr(p)-Lambda0| + r,
```

one also has the explicit moving-disc bound

```text
M_{E,p}(p)
<=
(|Lambda_fr(p)-Lambda0| + r) M_{K,Lambda p}(p).
```

This is exact normalized-block algebra followed by a conditional local
mixed-derivative bound. It reduces `M_{E,p}(p)` to a structural control of
the exact normalized matrix `K(Lambda,p)`, but it still requires verifying
the mixed derivative bound in the intended local region.

## 24. Exact-remainder ready-to-use derivative-comparison corollary

Assume the hypotheses of the frozen-coefficient ready-to-use
derivative-comparison corollary. Assume also that the exact-remainder
`p`-derivative control lemma applies at the same parameter value, so that

```text
M_{E,p}(p) <= rho_0(p) M_{K,Lambda p}(p).
```

Then, with

```text
C_Lambda(p)
=
sqrt(2) M_E(p) + (M_G(p) + M_E(p)) M_{E,1}(p),
```

and

```text
C_Lambda(p) + 2 B_shift(p) < c0,
```

the branch-shift derivative obeys the mixed-derivative ready-to-use bound

```text
|lambda'(p) - alpha'(p)|
<=
[
  ||K0'(p)||_F M_E(p)
  + (M_G(p) + M_E(p)) rho_0(p) M_{K,Lambda p}(p)
  + sqrt(2) ||K0'(p)||_F B_shift(p)
]
/
[c0 - C_Lambda(p) - 2 B_shift(p)]
+
[
  (sqrt(2) |alpha(p)-Lambda0| + ||K0(p)||_F) ||K0'(p)||_F
]
  [C_Lambda(p) + 2 B_shift(p)]
/
[
  (c0 - C_Lambda(p) - 2 B_shift(p)) c0
].
```

Equivalently, one may replace `rho_0(p)` by the explicit upper bound
`|Lambda_fr(p)-Lambda0|+r`.

Thus the `RUPDC-FC` input `M_{E,p}(p)` has been replaced by the structural
mixed-derivative control `M_{K,Lambda p}(p)` for the exact normalized block.

This remains a local conditional corollary. It does not prove that
`M_{K,Lambda p}(p)`, `||K0'(p)||_F`, the existing matrix controls, or the
frozen discriminant margin are small in the actual CoupledBeams regime.

## 25. CoupledBeams quantitative-regime derivative-comparison corollary

Assume a local CoupledBeams derivative-comparison quantitative regime
`Q_der(I_Q)` holds on a small parameter interval `I_Q` in the sense of
`foundations.md`. In particular, the canonical local construction through
`RUPDC-ER` is in force and the regime constants satisfy

```text
|a_p(Lambda)| >= m0 on D_r(p),
M_G(p) <= G_*,
M_E(p) <= E_*,
M_{E,1}(p) <= E_{1,*},
||K0(p)||_F <= K_{0,*},
||K0'(p)||_F <= K_{p,*},
M_{K,Lambda p}(p) <= K_{Lambda p,*},
|Lambda_fr(p)-Lambda0| <= A_*,
```

with frozen denominator margin

```text
|partial_Lambda g_p(Lambda_fr(p))| >= c0,
```

and with

```text
B_Q = [G_* E_* + (1/2) E_*^2] / m0,
C_Q = sqrt(2) E_* + (G_* + E_*) E_{1,*},
D_Q = c0 - C_Q - 2 B_Q > 0.
```

Then, writing

```text
lambda(p) = Lambda_ex(p),
alpha(p) = Lambda_fr(p),
```

the derivative comparison supplied by `RUPDC-ER` is directly usable on
`I_Q`: for every `p in I_Q`,

```text
|lambda'(p) - alpha'(p)|
<=
[
  K_{p,*} E_*
  + (G_* + E_*) (A_* + r) K_{Lambda p,*}
  + sqrt(2) K_{p,*} B_Q
]
/
D_Q
+
[
  (sqrt(2) A_* + K_{0,*}) K_{p,*}
]
  [C_Q + 2 B_Q]
/
[
  D_Q c0
].
```

This is only a conditional regime corollary. It packages the remaining
project-specific inputs into a single checkable regime `Q_der(I_Q)` and then
substitutes its uniform bounds into `RUPDC-ER`. It does not prove that the
regime holds for any actual CoupledBeams branch, does not give an asymptotic
law, and does not introduce symmetric normal form, project-defined
`delta-kappa`, a veering criterion, or branch-case application language.

## 26. Branch-ready `Q_der` verification protocol

For a proposed concrete CoupledBeams local regime, prepare a branch-ready
verification record with the fields specified in `foundations.md`:

1. candidate selection data;
2. local construction checks;
3. frozen-model checks;
4. matrix and derivative bounds;
5. final decision.

If the record supplies all data needed for `Q_der(I_Q)` and every listed
check passes, then the local regime `Q_der(I_Q)` is verified for that concrete
candidate and the CoupledBeams quantitative-regime derivative-comparison
corollary applies on `I_Q`.

If any required entry is absent, or if any construction, margin, derivative,
or closure check fails, then the record is incomplete and no `QDR` conclusion
may be drawn for that candidate.

This protocol is branch-ready but still conditional. It does not verify any
specific candidate by itself, does not construct a packet, does not prove the
frozen-discriminant margin, and does not establish an asymptotic law,
symmetric normal form, project-defined `delta-kappa`, a veering criterion, or
any branch-case theorem.
