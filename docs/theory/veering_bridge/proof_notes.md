# Proof Notes

## Scope

This file is the non-canonical working companion to the clean statements in
[lemma_statements.md](lemma_statements.md).
It records proof attempts, status distinctions, failure modes, and the main
reasons the bridge line can break.

### Canonical workflow reminder

New mathematics should be added by editing existing canonical files, not by
creating auxiliary notes. If a result cannot be summarized inside the current
five-file structure, it is probably not mature enough to enter the package yet.

## Status discipline

The package uses five labels:

- `proved_here`
- `accepted_standard_background`
- `explicit_hypothesis`
- `needs_caution`
- `not_reached_yet`

`proved_here` means the package contains the relevant local argument.
`accepted_standard_background` means the step is mathematically standard but
not reproved in full detail inside this package.

## 1. Packet construction in the simple-root setting

### What is actually proved here

- `det M(Lambda_j,p)=0` implies singularity of the square matrix
  `M(Lambda_j,p)`;
- once two right directions and two left directions are supplied at fixed `p`
  with rank two, they define rank-2 packets;
- once continuous representatives exist and the relevant determinants are
  nonzero at a base point, rank and pairing persist after shrinking.

### What is only accepted standard background

- simple determinant root implies one-dimensional right and left kernels;
- local continuation of simple roots;
- local representative selection for one-dimensional right and left kernel
  lines.

### Why this step can break

- roots cease to be simple;
- the two right or left directions lose independence;
- a continuous packet basis is not justified;
- the retained pairing degenerates.

### Why this is still useful

This is the first honest step from scalar root data to packet data. It is
deliberately narrow and does not yet touch complements.

## 2. Complement regularity

### The key separation

Two facts must not be conflated:

1. algebraic basis completion:

   ```text
   V,W -> [V Z], [W U];
   ```

2. regular complement:

   ```text
   det D(Lambda,p) != 0,
   D = U^T M Z.
   ```

The first is ordinary finite-dimensional linear algebra. The second is the
actual theorem-level regularity hypothesis.

### What is proved here

- basis completion exists at fixed `p`;
- the block matrix `M_tilde = [W U]^T M [V Z]` is well defined once the bases
  are supplied;
- basis completion does not imply invertibility of `D`;
- base-point invertibility of `D` gives only local invertibility after
  shrinking.

### Main failure modes

- `D` is already singular at the base point;
- `D` is invertible near the base point but loses invertibility on the claimed
  window;
- one silently treats basis completion as if it already implied regular
  complement.

## 3. Schur root capture

### Three layers

1. block determinant identity:

   ```text
   det M_tilde = det D det S;
   ```

2. basis-factor relation:

   ```text
   det M_tilde = det([W U]^T) det M det([V Z]);
   ```

3. local root capture:
   once all prefactors are nonzero on the claimed region, the local zeros of
   `det M` are exactly the zeros of `det S`.

### What is exact

The determinant factorization is exact once the packet, complements, and
regular complement hypothesis are supplied.

### Caution on multiplicity

Multiplicity is only discussed in the fixed-`p`, one-variable sense in
`Lambda`. The package does not use any multivariable multiplicity notion on
`Omega x I`.

### Main failure modes

- `det D` vanishes somewhere in the claimed region;
- basis determinants vanish;
- the local theorem is overstretched beyond the region where the nonzero
  factors are controlled.

## 4. Exact Lambda-normalization

### What is proved here

Given analytic `S(Lambda,p)` and the explicit hypothesis

```text
det(partial_Lambda S(Lambda0,p0)) != 0,
```

the exact local factorization

```text
S(Lambda,p) = R(Lambda,p) ((Lambda-Lambda0) I - K(Lambda,p))
```

is available after shrinking if needed.

### Why this matters

This is the first exact eigenvalue-style reduced form. It is still general:
`K` may depend on both `Lambda` and `p`.

### What it is not

It is not yet:

- a `p`-only model;
- a symmetric normal form;
- a `delta-kappa` representation.

### Main failure modes

- `partial_Lambda S(Lambda0,p0)` is singular;
- the invertible spectral coefficient is lost on a window larger than the
  region where the proof really works;
- a first-order model is mistaken for this exact normalization.

## 5. First-order freezing

### Exact versus approximate

The exact normalized block is

```text
F(Lambda,p) = ((Lambda-Lambda0) I - K(Lambda,p)).
```

The frozen model is

```text
G(Lambda,p) = ((Lambda-Lambda0) I - K0(p)),
K0(p) = K(Lambda0,p).
```

The package proves only the exact remainder formula

```text
F - G = -(Lambda-Lambda0) K1(Lambda,p),
```

with matrix-level size `O(|Lambda-Lambda0|)`.

### What this step does not prove

It does not prove exact root capture for the frozen model. It only provides a
controlled matrix approximation.

### Main failure modes

- `K` is not regular enough in `Lambda`;
- the local region is too large for the remainder to be informative;
- one silently transfers exact spectral statements from `F` to `G`.

## 6. Fixed-`p` root comparison

### The crucial distinction

Matrix-level approximation does not by itself imply spectral approximation.

To compare roots of `det F` and `det G`, the package fixes one parameter value
`p`, works with one-variable analytic determinants

```text
f_p = det F(.,p),
g_p = det G(.,p),
```

and adds a genuinely spectral local hypothesis on a chosen disc.

### What is proved here

- for `2x2` matrices with `F = G + E`,

  ```text
  det(G+E) - det(G) = tr(adj(G)E) + det(E);
  ```

- therefore the earlier matrix remainder gives a determinant perturbation
  bound on the boundary of a small disc;
- if the frozen determinant has one simple zero in that disc, no boundary
  zero, and the boundary domination inequality holds, then the exact
  normalized determinant has one nearby root in the same disc.

### What is only standard background

The step from boundary domination to equal zero count is imported as a
Rouche-type one-variable analytic root theorem.

### What this still does not give

- any global branch theorem in `p`;
- a quantitative root-shift bound;
- any veering interpretation.

## 7. Parameterwise moving-disc comparison and exact-branch continuation

### Pointwise-in-`p` Rouche versus branch regularity

Two layers must stay separate.

First layer:

```text
for each fixed p,
one Rouche-type comparison on one disc.
```

Second layer:

```text
as p varies,
those selected exact roots fit together into one regular branch.
```

The first layer is a parameterwise family of fixed-`p` arguments. The second
layer needs a separate simple-root continuation theorem.

### Shifted variable around the moving center

Let `Lambda_fr(p)` be a local frozen-model root branch of the frozen model and
set

```text
xi = Lambda - Lambda_fr(p).
```

Then the moving disc condition is simply

```text
|xi| <= r.
```

This makes clear that the real new object is not one fixed disc around one
fixed center, but one family of translated discs moving with the frozen root.

### What is genuinely new here

The real new ingredients are:

- one uniform radius `r > 0`;
- one local interval `I` in the scalar parameter `p`;
- one frozen-model root branch `Lambda_fr(p)`;
- one uniform boundary-domination estimate on all disc boundaries.

Pointwise domination at `p0` alone is not enough. The new theorem step needs a
single radius and a single strict domination margin that survive over a whole
local interval in `p`.

### Why "one root in each disc" is not yet a regular branch

Even if every moving disc contains exactly one exact root, that statement is
still set-theoretic in `p`: it gives a selected root value for each parameter,
but not yet a regular parametrization.

To obtain a local branch `Lambda_ex(p)`, one still needs:

- joint analyticity in `(Lambda,p)`;
- a simple exact root at the base point;
- an implicit-function-style continuation step.

### How uniqueness and continuation combine

The logical order is:

1. uniform moving-disc comparison gives one exact root in each disc;
2. at `p0`, if that exact root is simple, implicit-function-style continuation
   yields a local regular exact branch;
3. uniqueness in each disc prevents ambiguity, so the continued branch must be
   exactly the moving-disc selected root for nearby `p`.

This is the clean bridge from fixed-`p` nearby-root comparison to a local
parameterwise exact-root branch.

### Main failure modes

- one cannot choose a single uniform radius `r`;
- the frozen root is not isolated enough to support one moving disc family;
- domination holds at `p0` but not uniformly over the interval;
- the selected exact root loses simplicity;
- a moving disc captures competing exact or frozen roots.

### What this still does not give

This step still does not provide:

- a quantitative estimate for `|Lambda_ex(p)-Lambda_fr(p)|`;
- a global continuation theorem in `p`;
- a symmetric normal form;
- project-defined `delta-kappa`;
- any veering interpretation.

## 8. Quantitative parameterwise branch comparison

### Exact branch versus frozen branch

At this stage the bridge line has four distinct layers that must stay
separate:

1. exact normalized determinant:

   ```text
   f_p(Lambda) = det F(Lambda,p);
   ```
2. frozen-model determinant:

   ```text
   g_p(Lambda) = det G(Lambda,p);
   ```
3. existence and uniqueness in moving discs:
   for each `p`, one exact root lies in `D_r(p)`;
4. quantitative branch comparison:
   estimate the distance between `Lambda_ex(p)` and `Lambda_fr(p)`.

The new step concerns only layer 4. It must not be confused with the earlier
existence/uniqueness theorem.

### Where simple-root nondegeneracy is used

Two different simple-root inputs appear, and they play different roles.

- Exact simple-root nondegeneracy at the base point is used upstream, in the
  exact-root branch continuation step, to produce the local regular branch
  `Lambda_ex(p)`.
- Frozen simple-root nondegeneracy is used here, in quantitative form, to
  control how far `Lambda_ex(p)` can move away from `Lambda_fr(p)`.

The present estimate uses the second input directly.

### Preferred quantitative line

Work on the interval where the earlier moving-disc selection and exact-branch
continuation already hold. For each `p`, assume the frozen determinant admits
the factorization

```text
g_p(Lambda) = a_p(Lambda) (Lambda - Lambda_fr(p)),
|a_p(Lambda)| >= m0 > 0
```

on `D_r(p)`.

Define the determinant discrepancy by

```text
eta_loc(p) = sup_{Lambda in D_r(p)} |f_p(Lambda)-g_p(Lambda)|.
```

Since `Lambda_ex(p)` lies in `D_r(p)` and satisfies `f_p(Lambda_ex(p))=0`,
one has

```text
|g_p(Lambda_ex(p))|
=
|g_p(Lambda_ex(p)) - f_p(Lambda_ex(p))|
<=
eta_loc(p).
```

Using the frozen factorization at `Lambda_ex(p)` gives

```text
|a_p(Lambda_ex(p))| |Lambda_ex(p)-Lambda_fr(p)|
<=
eta_loc(p),
```

hence

```text
|Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p) / m0.
```

This is the cleanest honest quantitative statement currently available in the
package.

It is still only a pointwise-in-`p` estimate unless `eta_loc(p)` is itself
controlled uniformly. In the present notation, uniform smallness on an
interval is obtained only after controlling

```text
eta_* = sup_{p in I_ex} eta_loc(p).
```

### What is proved here and what is standard background

- `proved_here`:
  once the factorization and lower bound are assumed, the displacement
  estimate follows immediately by evaluating the determinant discrepancy at
  the exact root.
- `accepted_standard_background`:
  the passage from simple-root isolation or a derivative lower bound to the
  factorization

  ```text
  g_p(Lambda) = a_p(Lambda)(Lambda-Lambda_fr(p))
  ```

  with `a_p` nonvanishing on a small neighborhood.

If one prefers to formulate the hypothesis through
`|partial_Lambda g_p(Lambda_fr(p))| >= m0`, then converting that derivative
bound into the factorized lower bound above is standard local analytic-root
theory, not something reproved in full detail here.

### Why the derivative lower bound matters

The estimate divides by the frozen nondegeneracy constant. If the lower bound
for `|a_p|` deteriorates, then the comparison constant blows up and the
quantitative control becomes useless even if nearby-root selection still
holds.

This is exactly where the quantitative step is stricter than the earlier
existence theorem.

### Main failure modes

- the exact branch exists, but the frozen simple-root lower bound is too weak
  or collapses somewhere on the interval;
- the moving-disc radius is large enough that the factorized frozen
  coefficient is no longer uniformly bounded away from zero;
- the determinant discrepancy `eta_loc(p)` is finite but too large to give a
  useful comparison;
- one silently treats this quantitative estimate as if it were already a
  derivative formula or a gap law.

### What this still does not give

This step still does not provide:

- an asymptotic formula for `Lambda_ex(p)-Lambda_fr(p)`;
- a derivative identity for `d/dp (Lambda_ex-Lambda_fr)`;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- any veering criterion.

## 9. Determinant-discrepancy control from the matrix remainder

### Exact identity versus scalar bound

At the matrix level one has

```text
F(Lambda,p) = G(Lambda,p) + E(Lambda,p).
```

For `2x2` matrices, the determinant discrepancy is governed by the exact
identity

```text
det(G+E)-det(G)=tr(adj(G)E)+det(E).
```

This identity is exact. The scalar bound derived from it is a second step.

### Why the Frobenius norm is convenient here

In the current local `2x2` setting, the Frobenius norm gives a clean bound:

```text
|tr(adj(G)E)| <= ||adj(G)||_F ||E||_F = ||G||_F ||E||_F,
```

because for `2x2` matrices one has `||adj(G)||_F = ||G||_F`. Together with

```text
|det(E)| <= (1/2) ||E||_F^2,
```

this yields

```text
|f_p(Lambda)-g_p(Lambda)|
<=
||G(Lambda,p)||_F ||E(Lambda,p)||_F + (1/2)||E(Lambda,p)||_F^2.
```

### Moving-disc consequence

On a moving disc `D_r(p)`, define

```text
M_G(p)=sup_{Lambda in D_r(p)} ||G(Lambda,p)||_F,
M_E(p)=sup_{Lambda in D_r(p)} ||E(Lambda,p)||_F.
```

Then the determinant error indicator satisfies

```text
eta_loc(p)
<=
M_G(p) M_E(p) + (1/2) M_E(p)^2.
```

This is the missing constructive bridge from matrix-level freezing error to
the scalar discrepancy used in the earlier quantitative branch-comparison
lemma.

### Combined branch-distance estimate

Once the earlier `QPBC` hypotheses also hold, substitution into

```text
|Lambda_ex(p)-Lambda_fr(p)| <= eta_loc(p)/m0
```

gives the constructive bound

```text
|Lambda_ex(p)-Lambda_fr(p)|
<=
[M_G(p) M_E(p) + (1/2) M_E(p)^2] / m0.
```

This is still conditional. It does not yet produce a derivative law, because
it only controls the size of the branch shift through local matrix norms.

### Main failure modes

- `M_E(p)` is finite but not small enough for the bound to be useful;
- `M_G(p)` becomes large on the moving discs, weakening the estimate;
- the matrix bound is available only on discs that are too small to support
  the earlier moving-disc comparison cleanly;
- one confuses the exact determinant identity with the derived norm bound;
- one treats this constructive estimate as if it were already an asymptotic
  expansion.

### What this still does not give

This step still does not provide:

- a rate sharper than quadratic-linear control in the local matrix norms;
- a derivative formula for the branch shift;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- a veering criterion.

## 10. Moving-disc closure

### Why this is only a closure step

The constructive branch bound

```text
|Lambda_ex(p)-Lambda_fr(p)| <= B_shift(p)
```

does not by itself create a new exact root or a new branch. The earlier
moving-disc comparison and exact-branch continuation already do that work.

The only new point here is self-consistency: if

```text
B_shift(p) < r,
```

then the current quantitative estimate itself certifies that the selected
exact root lies inside the chosen moving disc.

### Uniform version

If

```text
B_* = sup_{p in I_ex} B_shift(p) < r,
```

then the same certification holds for every `p in I_ex`.

This is useful because the moving-disc framework is no longer justified only
by the earlier existence theorem; it is now also compatible with the current
quantitative estimate.

### Failure mode

If `B_shift(p) >= r`, then the existing constructive estimate is too weak to
certify disc containment. The earlier existence theorem may still guarantee
that the root lies in the disc, but the quantitative bound no longer proves
that fact by itself.

## 11. First-order root-shift formula

### What this step tries to add

The earlier quantitative step gives only a coarse distance bound. The next
local question is whether the branch shift admits a first-order comparison
formula around the frozen root.

Write, at fixed `p`,

```text
alpha = Lambda_fr(p),
lambda = Lambda_ex(p),
delta = lambda - alpha,
h_p(Lambda) = f_p(Lambda) - g_p(Lambda).
```

The target is

```text
delta
=
- h_p(alpha) / partial_Lambda g_p(alpha)
+ Rem_shift(p),
```

which is the same as

```text
Lambda_ex(p) - Lambda_fr(p)
=
- (f_p(Lambda_fr(p)) - g_p(Lambda_fr(p))) / partial_Lambda g_p(Lambda_fr(p))
+ Rem_shift(p).
```

### Proof line

Because `f_p(lambda)=0`, the one-variable Taylor formula with integral
remainder gives

```text
0 = f_p(alpha) + delta * A_p,
```

where

```text
A_p = integral_0^1 partial_Lambda f_p(alpha + t delta) dt.
```

Hence

```text
delta = - f_p(alpha) / partial_Lambda g_p(alpha)
        - delta * (A_p - partial_Lambda g_p(alpha))
          / partial_Lambda g_p(alpha).
```

Since `g_p(alpha)=0`, the leading term may also be written as

```text
- (f_p(alpha)-g_p(alpha)) / partial_Lambda g_p(alpha).
```

Define

```text
Rem_shift(p)
=
- delta * (A_p - partial_Lambda g_p(alpha))
  / partial_Lambda g_p(alpha).
```

This is an exact identity.

### Where the remainder bound comes from

To bound `A_p - partial_Lambda g_p(alpha)`, separate the difference into two
pieces:

1. derivative discrepancy between `f_p` and `g_p`;
2. variation of `partial_Lambda g_p` away from the frozen root.

If

```text
H_Delta(p) = sup_{Lambda in D_r(p)} |partial_Lambda(f_p-g_p)(Lambda)|
```

and

```text
H_g2(p) = sup_{Lambda in D_r(p)} |partial_Lambda^2 g_p(Lambda)|,
```

then

```text
|A_p - partial_Lambda g_p(alpha)|
<=
H_Delta(p) + (1/2) H_g2(p) |delta|.
```

Therefore, if

```text
|partial_Lambda g_p(alpha)| >= m1(p) > 0,
```

one obtains

```text
|Rem_shift(p)|
<=
[H_Delta(p) |delta| + (1/2) H_g2(p) |delta|^2] / m1(p).
```

Substituting the earlier constructive bound

```text
|delta| <= B_shift(p)
```

gives the stated disc-level remainder control.

### What is proved here and what is only background

- `proved_here`:
  the algebraic rearrangement from the Taylor identity to the first-order
  shift formula and the explicit remainder estimate.
- `accepted_standard_background`:
  the one-variable Taylor formula with integral remainder.
- `explicit_hypothesis`:
  exact branch existence, moving-disc closure, a lower bound for
  `|partial_Lambda g_p(Lambda_fr(p))|`, and disc-level bounds for
  `partial_Lambda(f_p-g_p)` and `partial_Lambda^2 g_p`.

### Why this is still not a final asymptotic law

The formula is exact, but the usefulness of the first-order term depends on
the size of the remainder. If the derivative discrepancy or second-derivative
bound is too large, then the remainder may be of the same order as the
leading term.

So this step gives a local first-order comparison law with a controlled
remainder, not yet a final asymptotic expansion.

### Main failure modes

- the lower bound for `|partial_Lambda g_p(Lambda_fr(p))|` deteriorates;
- the earlier closure condition fails, so the disc-level derivative bounds do
  not automatically apply at `Lambda_ex(p)`;
- `H_Delta(p)` is too large, so the discrepancy of derivatives dominates;
- `H_g2(p)` is too large, so the quadratic correction is not small;
- one mistakes this conditional first-order law for a fully constructive
  derivative formula.

## 12. Derivative control from reduced-model objects

### Why this is the real bridge to a model-linked law

The first-order root-shift formula introduced the quantities

```text
|partial_Lambda g_p(Lambda_fr(p))|,
sup_{D_r(p)} |partial_Lambda(f_p-g_p)|,
sup_{D_r(p)} |partial_Lambda^2 g_p|.
```

Up to that point they were still abstract scalar hypotheses. The present step
links them back to the concrete reduced-model objects `G`, `E`, and
`partial_Lambda E`.

### Exact frozen-derivative identities

Because

```text
G(Lambda,p) = ((Lambda-Lambda0)I - K0(p))
```

is affine in `Lambda` with slope `I`, the frozen determinant is a monic
quadratic polynomial in `Lambda`. In the local `2x2` setting one gets the
exact identities

```text
partial_Lambda g_p(Lambda) = tr G(Lambda,p)
                           = 2(Lambda-Lambda0) - tr K0(p),
partial_Lambda^2 g_p(Lambda) = 2.
```

Therefore the second-derivative quantity in `FOSF` is no longer abstract at
all, and the first-derivative denominator is reduced to an explicit scalar
expression in the frozen reduced model.

### Where the remaining noncoalescence hypothesis sits

This step does not prove a positive lower bound for

```text
|partial_Lambda g_p(Lambda_fr(p))|
=
|2(Lambda_fr(p)-Lambda0) - tr K0(p)|.
```

It only shows exactly what must stay away from zero. So the denominator
assumption in `FOSF` is now model-linked but still explicit.

### Derivative discrepancy bound

Differentiate the exact determinant-discrepancy identity

```text
f_p-g_p = tr(adj(G)E) + det(E).
```

Since `G` is affine in `Lambda` with slope `I`, one has for `2x2` matrices

```text
partial_Lambda adj(G) = I.
```

Hence

```text
partial_Lambda(f_p-g_p)
=
tr(E) + tr(adj(G) partial_Lambda E) + tr(adj(E) partial_Lambda E).
```

Now use Frobenius trace bounds:

```text
|tr(E)| <= sqrt(2) ||E||_F,
|tr(adj(G) partial_Lambda E)| <= ||G||_F ||partial_Lambda E||_F,
|tr(adj(E) partial_Lambda E)| <= ||E||_F ||partial_Lambda E||_F.
```

Taking suprema on the moving discs gives

```text
sup_{D_r(p)} |partial_Lambda(f_p-g_p)|
<=
sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p).
```

This is the first point where the derivative discrepancy used in `FOSF` is
controlled directly by reduced-model matrix quantities.

### What is exact and what is only bounded

- exact:
  `partial_Lambda g_p = tr G`,
  `partial_Lambda^2 g_p = 2`,
  and the differentiated discrepancy identity;
- bounded:
  the disc-level control of `partial_Lambda(f_p-g_p)` via the local Frobenius
  suprema.

### Main failure modes

- `M_{E,1}(p)` is not finite or is too large to make the derivative
  discrepancy bound useful;
- `|2(Lambda_fr(p)-Lambda0) - tr K0(p)|` becomes small, so the denominator in
  `FOSF` is still badly conditioned;
- the moving discs are too large for the local Frobenius suprema to remain
  informative;
- one mistakes exact derivative identities for a full noncoalescence theorem.

### What this still does not give

This step still does not provide:

- a constructive proof that the frozen derivative denominator is uniformly
  separated from zero;
- a final asymptotic or derivative-level shift law;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- a veering criterion.

## 13. Frozen-model denominator control from noncoalescence

### Why the derivative at a simple frozen root equals root separation

Once a second local frozen root branch `Lambda_fr,2(p)` is available, the
frozen determinant is a monic quadratic in `Lambda`, so at fixed `p` one has

```text
g_p(Lambda)
=
(Lambda-Lambda_fr(p))(Lambda-Lambda_fr,2(p)).
```

Differentiating gives

```text
partial_Lambda g_p(Lambda)
=
(Lambda-Lambda_fr,2(p)) + (Lambda-Lambda_fr(p)).
```

Evaluating at `Lambda = Lambda_fr(p)` yields

```text
partial_Lambda g_p(Lambda_fr(p))
=
Lambda_fr(p) - Lambda_fr,2(p).
```

So the denominator in `FOSF` is not an arbitrary derivative quantity anymore:
it is exactly the frozen-root separation from the selected root to the second
frozen root.

### Equivalent discriminant reading

Because `G(Lambda,p)=((Lambda-Lambda0)I-K0(p))`, the frozen determinant is the
characteristic polynomial of `K0(p)` in the shifted variable `Lambda-Lambda0`.
Hence the squared denominator may also be read as the frozen quadratic
discriminant:

```text
(partial_Lambda g_p(Lambda_fr(p)))^2
=
(tr K0(p))^2 - 4 det K0(p).
```

This is an exact identity, but the package does not yet use it as a concrete
regime test.

### Why this is the correct next bridge

The previous derivative-control step reduced the `FOSF` denominator to the
explicit scalar

```text
|partial_Lambda g_p(Lambda_fr(p))|
=
|2(Lambda_fr(p)-Lambda0) - tr K0(p)|.
```

The present step identifies that same quantity geometrically as

```text
|Lambda_fr(p)-Lambda_fr,2(p)|.
```

So the remaining denominator hypothesis is narrowed from a free derivative
lower bound to a frozen noncoalescence condition.

### Moving-disc interpretation

If the second frozen root lies outside `D_r(p)`, then

```text
|Lambda_fr(p)-Lambda_fr,2(p)| > r,
```

and therefore

```text
|partial_Lambda g_p(Lambda_fr(p))| > r.
```

This is the cleanest way to connect denominator control with the existing
moving-disc language.

### Where the argument can fail

- a second local frozen root branch is not available on the interval;
- the two frozen roots coalesce, so the denominator vanishes;
- the second frozen root enters the moving disc, destroying the useful lower
  bound;
- one confuses the exact denominator identity with a proof that the needed
  separation hypothesis actually holds for CoupledBeams.

### What is exact and what remains conditional

- exact:
  once both frozen root branches are available, the factorization of `g_p`
  and the identity

  ```text
  partial_Lambda g_p(Lambda_fr(p))
  =
  Lambda_fr(p)-Lambda_fr,2(p)
  ```

  are exact;
- conditional:
  any positive lower bound for the denominator still comes from an explicit
  frozen-root separation or second-root exclusion hypothesis.

## 14. Frozen discriminant persistence and denominator control

### Why this is the correct next bridge

The previous denominator-control step still phrased the useful hypothesis in
root language:

```text
|Lambda_fr(p)-Lambda_fr,2(p)| >= sigma_fr > 0.
```

The present step rewrites that same condition at coefficient level through the
frozen `2x2` model itself. Since

```text
G(Lambda,p) = ((Lambda-Lambda0)I-K0(p)),
```

one has the exact quadratic formula

```text
g_p(Lambda)
=
(Lambda-Lambda0)^2 - tr K0(p) (Lambda-Lambda0) + det K0(p),
```

so the discriminant is

```text
Disc_fr(p) = (tr K0(p))^2 - 4 det K0(p).
```

This is the right next bridge because it replaces a root-separation condition
by a coefficient-level condition on `K0(p)`.

### Distinct roots, discriminant, and denominator

For a monic quadratic, the following are equivalent at fixed `p`:

1. the two frozen roots are distinct;
2. `Disc_fr(p) != 0`;
3. `partial_Lambda g_p(Lambda_fr(p)) != 0`.

Once the roots are locally labeled on a sufficiently small interval after
shrinking,

```text
Disc_fr(p) = (Lambda_fr(p)-Lambda_fr,2(p))^2,
```

and therefore

```text
|partial_Lambda g_p(Lambda_fr(p))|
=
|Lambda_fr(p)-Lambda_fr,2(p)|
=
|Disc_fr(p)|^(1/2).
```

So the remaining denominator in `FOSF` is now tied directly to the frozen
discriminant.

### Local persistence near a base point

If

```text
Disc_fr(p0) != 0,
```

then continuity immediately gives a smaller local interval `I0` on which

```text
Disc_fr(p) != 0.
```

After shrinking, one may locally label the two frozen roots on `I0`. On that
smaller interval the roots remain distinct, and because `|Disc_fr|` is
continuous and nonzero on the compact closure of `I0`, there exists `c0 > 0`
such that

```text
|partial_Lambda g_p(Lambda_fr(p))| >= c0
```

for all `p in I0`.

So the honest local conclusion is:

```text
Disc_fr(p0) != 0
->
after shrinking to I0 and locally labeling the frozen roots,
there exists c0 > 0 with |partial_Lambda g_p(Lambda_fr(p))| >= c0 on I0.
```

### What is exact and what is only standard background

- exact:
  the explicit quadratic formula for `g_p`, the discriminant identity, and
  the relation

  ```text
  Disc_fr(p) = (Lambda_fr(p)-Lambda_fr,2(p))^2;
  ```
- `accepted_standard_background`:
  local labeling of the two roots, or equivalently local square-root choice
  for a nonvanishing discriminant, after shrinking.

### Where the argument can fail

- `Disc_fr(p0) = 0`, so the frozen roots coalesce already at the base point;
- the discriminant is nonzero at `p0` but becomes too small on the interval
  one wants to use, so the resulting denominator lower bound is weak;
- one mistakes local persistence near `p0` for a proof that the needed
  discriminant bound is already available for the actual CoupledBeams sweep;
- one overreads this coefficient-level control as a symmetry or veering
  statement.

### What this still does not give

This step still does not provide:

- a project-specific proof that `Disc_fr(p)` is bounded away from zero on the
  desired regime;
- a sharper asymptotic or derivative-level shift law;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- a veering criterion.

## 15. Frozen discriminant persistence from coefficient control

### Main coefficient-level bound

Write

```text
T(p) = tr K0(p),
D0(p) = det K0(p),
Disc_fr(p) = T(p)^2 - 4 D0(p).
```

Then

```text
Disc_fr(p) - Disc_fr(p0)
=
(T(p)-T(p0))(T(p)+T(p0)) - 4 (D0(p)-D0(p0)).
```

If on an interval `I` around `p0` one has

```text
|T(p)| <= M_T,
|T(p)-T(p0)| <= eps_T,
|D0(p)-D0(p0)| <= eps_D,
```

then, because `p0 in I` also implies `|T(p0)| <= M_T`,

```text
|T(p)+T(p0)| <= 2 M_T,
```

and therefore

```text
|Disc_fr(p)-Disc_fr(p0)|
<=
2 M_T eps_T + 4 eps_D.
```

This is the clean coefficient-level estimate that pushes frozen noncoalescence
back to trace/determinant control.

### Nonvanishing implication

If

```text
|Disc_fr(p0)| > 2 M_T eps_T + 4 eps_D,
```

then

```text
Disc_fr(p) != 0
```

throughout `I`.

After that, the previous frozen-discriminant persistence step applies: after
possibly shrinking to a smaller interval `I0` and locally labeling the frozen
roots, one obtains a constant `c0 > 0` such that

```text
|partial_Lambda g_p(Lambda_fr(p))| >= c0
```

for all `p in I0`.

So the new bridge is:

```text
coefficient control on T and D0
->
discriminant stays nonzero
->
local frozen-root labeling
->
uniform local denominator lower bound for FOSF.
```

### Convenient derivative-based variant

Sometimes it is easier to control derivatives of the frozen coefficients than
their direct oscillation. In that case,

```text
Disc_fr'(p) = 2 T(p) T'(p) - 4 D0'(p).
```

Hence, if

```text
|T(p)| <= M_T,
|T'(p)| <= M_{T,1},
|D0'(p)| <= M_{D,1},
```

then

```text
|Disc_fr'(p)| <= 2 M_T M_{T,1} + 4 M_{D,1}.
```

By a mean-value estimate, after shrinking around `p0`, nonvanishing of
`Disc_fr(p0)` persists on a small interval. This is a convenient corollary,
not the main theorem statement.

### What is exact and what is only standard background

- `proved_here`:
  the algebraic bound for `|Disc_fr(p)-Disc_fr(p0)|` from coefficient control,
  and the resulting explicit nonvanishing implication;
- `accepted_standard_background`:
  the continuity/mean-value shrinking step used to pass from pointwise
  nonvanishing of `Disc_fr(p0)` to local interval persistence;
- `explicit_hypothesis`:
  local coefficient bounds on `T` and `D0`, and a base-point discriminant
  margin larger than the explicit error term.

### Main failure modes

- the coefficient oscillation is too large relative to `|Disc_fr(p0)|`;
- `|T(p)|` is not controlled on the interval, so the simple bound on
  `|T(p)+T(p0)|` is unavailable;
- the discriminant margin is nonzero but too small to give a useful interval;
- one mistakes this local coefficient-level test for a proof that the actual
  CoupledBeams regime already satisfies it.

## 16. Ready-to-use local first-order shift corollary

### Why this is only a corollary

No new theorem mechanism is being introduced here. The point is to stop
leaving the usable first-order shift estimate scattered across five earlier
steps.

The input chain is:

1. moving-disc selection and exact-branch continuation supply
   `Lambda_fr(p)`, `Lambda_ex(p)`, and the disc family `D_r(p)`;
2. `DDMR`, together with the quantitative comparison and closure steps,
   supplies the constructive bound

   ```text
   |Lambda_ex(p)-Lambda_fr(p)| <= B_shift(p);
   ```

3. `DCL` supplies the derivative-discrepancy and frozen second-derivative
   controls

   ```text
   H_Delta(p)
   <=
   sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p),
   H_g2(p) = 2;
   ```

4. `FDPL`, or `FDCC` followed by `FDPL`, supplies the denominator lower bound
   from the frozen discriminant:

   ```text
   |partial_Lambda g_p(Lambda_fr(p))| >= c0.
   ```

Substituting these three displayed inputs into `FOSF` gives the corollary.
The abstract `m1(p)`, `H_Delta(p)`, and `H_g2(p)` slots in `FOSF` are replaced
by `c0`, `M_G(p)`, `M_E(p)`, `M_{E,1}(p)`, and the exact value `2`.

### What is ready-to-use now

For every `p` on the local interval where the prior inputs hold,

```text
Lambda_ex(p) - Lambda_fr(p)
=
- (f_p(Lambda_fr(p)) - g_p(Lambda_fr(p)))
  / partial_Lambda g_p(Lambda_fr(p))
+ Rem_shift(p),
```

with

```text
|Rem_shift(p)|
<=
[
  (sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p)) B_shift(p)
  + B_shift(p)^2
] / c0.
```

The coefficient `1` multiplying `B_shift(p)^2` is not a new estimate. It is
the old `FOSF` contribution

```text
(1/2) H_g2(p) B_shift(p)^2
```

after `DCL` has identified `H_g2(p)=2` for the frozen local `2x2` model.

So a later argument can cite this single corollary instead of separately
collecting the first-order formula, determinant-discrepancy bound,
derivative-discrepancy bound, and frozen-discriminant denominator control.

### Honest boundary of the result

The corollary is still local and conditional. It does not prove that the
CoupledBeams sweep satisfies the coefficient bounds on `T(p)` and `D0(p)`, it
does not prove that the matrix controls are small in a useful scale, and it
does not turn the first-order comparison into an asymptotic law.

It also does not introduce symmetric normal form, project-defined
`delta-kappa`, a veering criterion, or branch-case application language.

## 17. Parameterwise derivative comparison

### Exact branch derivative

Start from the already selected exact branch:

```text
f_p(Lambda_ex(p)) = 0.
```

If the branch equation is differentiable in `p` and

```text
partial_Lambda f_p(Lambda_ex(p)) != 0,
```

then ordinary chain-rule differentiation gives

```text
partial_p f_p(Lambda_ex(p))
+
partial_Lambda f_p(Lambda_ex(p)) Lambda_ex'(p)
=
0,
```

hence

```text
Lambda_ex'(p)
=
- partial_p f_p(Lambda_ex(p))
  / partial_Lambda f_p(Lambda_ex(p)).
```

The new issue is not the algebra. The issue is keeping the exact denominator
away from zero using quantities already controlled in the package.

### Exact denominator closure

Compare the exact denominator to the frozen one at the frozen root:

```text
partial_Lambda f_p(Lambda_ex(p))
-
partial_Lambda g_p(Lambda_fr(p)).
```

Split it as

```text
[partial_Lambda(f_p-g_p)(Lambda_ex(p))]
+
[partial_Lambda g_p(Lambda_ex(p))
 - partial_Lambda g_p(Lambda_fr(p))].
```

`DCL` bounds the first bracket by

```text
C_Lambda(p)
=
sqrt(2) M_E(p) + (M_G(p)+M_E(p)) M_{E,1}(p),
```

and the frozen local `2x2` identity `partial_Lambda^2 g_p=2` bounds the
second bracket by `2 |Lambda_ex(p)-Lambda_fr(p)|`. With the previous
constructive branch bound this gives

```text
|partial_Lambda f_p(Lambda_ex(p))
 - partial_Lambda g_p(Lambda_fr(p))|
<=
C_Lambda(p) + 2 B_shift(p).
```

Therefore the frozen denominator lower bound `c0` gives the exact denominator
lower bound if

```text
C_Lambda(p) + 2 B_shift(p) < c0.
```

This is the honest extra closure condition for differentiating the exact
branch in a controlled way.

### Frozen branch derivative

For the frozen branch,

```text
g_p(Lambda_fr(p)) = 0.
```

The same chain-rule calculation gives

```text
Lambda_fr'(p)
=
- partial_p g_p(Lambda_fr(p))
  / partial_Lambda g_p(Lambda_fr(p)).
```

Since

```text
g_p(Lambda)
=
(Lambda-Lambda0)^2 - T(p)(Lambda-Lambda0) + D0(p),
```

one has

```text
partial_p g_p(Lambda)
=
-T'(p)(Lambda-Lambda0) + D0'(p),
```

so

```text
Lambda_fr'(p)
=
[T'(p)(Lambda_fr(p)-Lambda0)-D0'(p)]
/
[2(Lambda_fr(p)-Lambda0)-T(p)].
```

This is exact for the frozen quadratic once the denominator is nonzero.

### Difference identity and bound

Let

```text
lambda = Lambda_ex(p),
alpha = Lambda_fr(p).
```

Subtracting the two branch-derivative formulas and using the exact
denominators gives the identity

```text
lambda' - alpha'
=
- [partial_p f_p(lambda) - partial_p g_p(alpha)]
  / partial_Lambda f_p(lambda)
+
partial_p g_p(alpha)
  [partial_Lambda f_p(lambda) - partial_Lambda g_p(alpha)]
  /
  [partial_Lambda f_p(lambda) partial_Lambda g_p(alpha)].
```

The second numerator is controlled by the denominator-closure estimate above.
For the first numerator, split

```text
partial_p f_p(lambda) - partial_p g_p(alpha)
=
partial_p(f_p-g_p)(lambda)
+
[partial_p g_p(lambda)-partial_p g_p(alpha)].
```

The new quantity is

```text
H_pDelta(p)
=
sup_{Lambda in D_r(p)} |partial_p(f_p-g_p)(Lambda)|.
```

The frozen part is explicit because

```text
partial_p g_p(lambda)-partial_p g_p(alpha)
=
-T'(p)(lambda-alpha).
```

Using `|lambda-alpha| <= B_shift(p)` gives the displayed conditional
derivative-difference bound in `lemma_statements.md`.

### Honest boundary

This step proves the local derivative identities and a conditional comparison
bound, but it does not yet prove that `H_pDelta(p)`, `T'(p)`, or `D0'(p)` are
small in a project-specific regime. That is now the real obstruction before
any asymptotic derivative-level shift law can be claimed.

It also does not introduce symmetric normal form, project-defined
`delta-kappa`, a veering criterion, or branch-case application language.

## 18. `p`-derivative discrepancy control from reduced-model objects

### Exact identity first

The starting point is still the exact local `2x2` determinant discrepancy
identity:

```text
f_p-g_p = tr(adj(G)E) + det E.
```

Differentiate this at fixed `Lambda` with respect to `p`. For `2x2` matrices
the adjugate map is linear, so

```text
partial_p adj(G) = adj(partial_p G).
```

Also,

```text
partial_p det E = tr(adj(E) partial_p E).
```

Therefore

```text
partial_p(f_p-g_p)
=
tr(adj(partial_p G)E)
+ tr(adj(G)partial_p E)
+ tr(adj(E)partial_p E).
```

This is the exact algebraic line. No smallness or approximation has entered.

### Frobenius bound

Use the Frobenius trace inequality and, for `2x2` matrices,

```text
||adj(A)||_F = ||A||_F.
```

Then

```text
|tr(adj(partial_p G)E)|
<=
||partial_p G||_F ||E||_F,
```

```text
|tr(adj(G)partial_p E)|
<=
||G||_F ||partial_p E||_F,
```

and

```text
|tr(adj(E)partial_p E)|
<=
||E||_F ||partial_p E||_F.
```

Taking suprema on the moving disc gives

```text
H_pDelta(p)
<=
M_{G,p}(p) M_E(p)
+ (M_G(p)+M_E(p)) M_{E,p}(p).
```

This is the direct `p`-analogue of the earlier constructive controls, but it
uses `partial_p G` and `partial_p E` rather than `partial_Lambda E`.

### Ready-to-use consequence for PDC

Substitute this bound into the `PDC` derivative-difference estimate. The only
change is that the abstract term

```text
H_pDelta(p)
```

is replaced by

```text
M_{G,p}(p) M_E(p)
+ (M_G(p)+M_E(p)) M_{E,p}(p).
```

This gives a genuinely more usable derivative-comparison corollary because a
later argument can work with matrix-level quantities from the reduced model
instead of a free scalar discrepancy hypothesis.

### Honest boundary

The new bound is constructive but not yet project-specific. The next
obstruction is to control `M_{G,p}`, `M_{E,p}`, `T'`, and `D0'` in the actual
local regime. Without those estimates, the derivative comparison remains a
local conditional theorem and not an asymptotic branch-shift law.

It also does not introduce symmetric normal form, project-defined
`delta-kappa`, a veering criterion, or branch-case application language.

## 19. Frozen coefficient-derivative control for RUPDC

### Exact frozen coefficient identities

In the frozen local `2x2` model,

```text
G(Lambda,p) = ((Lambda-Lambda0)I - K0(p)).
```

Since the `Lambda` part is independent of `p`,

```text
partial_p G(Lambda,p) = -K0'(p).
```

This is stronger than a generic moving-disc estimate: the right-hand side is
independent of `Lambda`. Therefore

```text
M_{G,p}(p)
=
sup_{Lambda in D_r(p)} ||partial_p G(Lambda,p)||_F
=
||K0'(p)||_F.
```

Writing

```text
K0(p) =
[ a(p)  b(p)
  c(p)  d(p) ],
```

the frozen trace and determinant coefficients satisfy

```text
T'(p) = a'(p)+d'(p) = tr K0'(p),
```

and

```text
D0'(p)
=
a'(p)d(p) + a(p)d'(p) - b'(p)c(p) - b(p)c'(p)
=
tr(adj(K0(p))K0'(p)).
```

These are exact coefficient identities in the frozen model, not asymptotic
approximations.

### Frobenius bounds

The trace bound is just Cauchy-Schwarz in the Frobenius inner product:

```text
|T'(p)|
=
|tr(I K0'(p))|
<=
||I||_F ||K0'(p)||_F
=
sqrt(2)||K0'(p)||_F.
```

For the determinant coefficient,

```text
|D0'(p)|
=
|tr(adj(K0(p))K0'(p))|
<=
||adj(K0(p))||_F ||K0'(p)||_F.
```

In the `2x2` setting, `||adj(K0)||_F=||K0||_F`, giving

```text
|D0'(p)| <= ||K0(p)||_F ||K0'(p)||_F.
```

Consequently the frozen numerator that appears in the `RUPDC` denominator
comparison obeys

```text
|-T'(p)(alpha(p)-Lambda0)+D0'(p)|
<=
(sqrt(2)|alpha(p)-Lambda0|+||K0(p)||_F)||K0'(p)||_F.
```

The constants here are not hidden theorem assumptions: `sqrt(2)` is
`||I||_F` in dimension two, and the determinant bound uses the exact `2x2`
identity for the adjugate norm.

### Ready-to-use substitution into RUPDC

Substituting the displayed estimates into `RUPDC` removes the frozen
coefficient derivative inputs

```text
M_{G,p}(p), T'(p), D0'(p)
```

as independent quantities. The resulting estimate still depends on

```text
||K0'(p)||_F, M_{E,p}(p),
M_G(p), M_E(p), M_{E,1}(p), B_shift(p), c0.
```

This is worth recording as a corollary because it is directly reusable:
future arguments can cite one bound rather than repeatedly substituting
`M_{G,p}=||K0'||_F` and the two coefficient inequalities by hand.

### Honest boundary

The new step controls only the frozen-model derivative inputs. It does not
control the exact-remainder derivative `M_{E,p}(p)`, does not prove
`||K0'(p)||_F` is small, and does not verify the frozen-discriminant margin
or the existing matrix controls in any project-specific regime.

It also does not introduce symmetric normal form, project-defined
`delta-kappa`, a veering criterion, or branch-case application language.

## 20. Exact-remainder `p`-derivative control

### Exact identities

The exact normalized and frozen blocks are

```text
F(Lambda,p) = ((Lambda-Lambda0)I - K(Lambda,p)),
G(Lambda,p) = ((Lambda-Lambda0)I - K0(p)),
K0(p)=K(Lambda0,p).
```

Therefore

```text
E(Lambda,p)
=
F(Lambda,p)-G(Lambda,p)
=
-(K(Lambda,p)-K0(p)).
```

Differentiating at fixed `Lambda` gives the exact identity

```text
partial_p E(Lambda,p)
=
-(partial_p K(Lambda,p)-partial_p K(Lambda0,p)).
```

So `M_{E,p}` is not an unrelated scalar input. It measures how much the
`p`-derivative of the exact normalized matrix changes between `Lambda0` and
the moving-disc point.

### Segment integral bound

Assume the straight segment from `Lambda0` to `Lambda` stays in the local
normalization region. Then

```text
partial_p K(Lambda,p)-partial_p K(Lambda0,p)
=
(Lambda-Lambda0)
integral_0^1
partial_{Lambda p}K(Lambda0+t(Lambda-Lambda0),p) dt.
```

Hence

```text
||partial_p E(Lambda,p)||_F
<=
|Lambda-Lambda0|
sup_{0<=t<=1}
||partial_{Lambda p}K(Lambda0+t(Lambda-Lambda0),p)||_F.
```

Taking the supremum over the moving disc gives

```text
M_{E,p}(p)
<=
rho_0(p) M_{K,Lambda p}(p),
```

where

```text
rho_0(p) = sup_{Lambda in D_r(p)} |Lambda-Lambda0|
```

and `M_{K,Lambda p}(p)` is the corresponding mixed-derivative supremum on the
segment set.

### Why the factor is not automatically `r`

The moving disc is centered at the frozen branch:

```text
D_r(p) = { Lambda : |Lambda-Lambda_fr(p)| <= r }.
```

The exact-remainder identity compares `K(Lambda,p)` with `K(Lambda0,p)`, so
the natural distance factor is `|Lambda-Lambda0|`, not
`|Lambda-Lambda_fr(p)|`. Thus

```text
rho_0(p) <= |Lambda_fr(p)-Lambda0| + r.
```

Only in a special normalization where the moving disc is centered at
`Lambda0`, or after another argument gives `rho_0(p) <= r`, does the bound
become an `r M_{K,Lambda p}` estimate.

### Ready-to-use consequence

Substituting

```text
M_{E,p}(p) <= rho_0(p) M_{K,Lambda p}(p)
```

into `RUPDC-FC` removes the last abstract exact-remainder derivative input
from the derivative comparison. The resulting corollary is genuinely reusable
because later work can estimate the mixed derivative of the exact normalized
matrix `K` rather than working with the free quantity `M_{E,p}`.

### Honest boundary

This is still local and conditional. It requires:

- the exact normalized matrix `K` to have the needed mixed derivative;
- the segment set from `Lambda0` to the moving discs to remain inside the
  local normalization region;
- a useful bound for `M_{K,Lambda p}`.

It does not prove project-specific smallness, a derivative-level asymptotic
law, symmetric normal form, project-defined `delta-kappa`, a veering
criterion, or branch-case application language.

## 21. Normalization-derivative audit after ERPC

### Exact algebra through `R` and `S`

The exact normalization gives

```text
S(Lambda,p)
=
R(Lambda,p) ((Lambda-Lambda0)I - K(Lambda,p)),
det R(Lambda,p) != 0.
```

Thus, wherever `R` is invertible,

```text
((Lambda-Lambda0)I - K(Lambda,p))
=
R(Lambda,p)^(-1) S(Lambda,p),
```

and therefore

```text
K(Lambda,p)
=
(Lambda-Lambda0)I - R(Lambda,p)^(-1) S(Lambda,p).
```

Set

```text
A(Lambda,p)=R(Lambda,p)^(-1).
```

Then

```text
partial_p A = -A (partial_p R) A.
```

Differentiating `K=(Lambda-Lambda0)I-AS` at fixed `Lambda` gives the exact
identity

```text
partial_p K
=
A (partial_p R) A S - A (partial_p S).
```

At the frozen point this becomes

```text
K0'(p)
=
A(Lambda0,p) (partial_p R)(Lambda0,p)
 A(Lambda0,p) S(Lambda0,p)
-
A(Lambda0,p) (partial_p S)(Lambda0,p).
```

This is exact algebra. It expresses `K0'` through normalization data, but it
requires controls for `R^{-1}`, `partial_p R`, `S`, and `partial_p S` at
`Lambda0`.

### Mixed derivative algebra

Differentiating the same `partial_p K` identity in `Lambda` gives

```text
partial_{Lambda p} K
=
partial_Lambda[
  A (partial_p R) A S - A (partial_p S)
].
```

Writing subscripts for derivatives, and using

```text
A_Lambda = -A R_Lambda A,
```

one obtains the expanded exact identity

```text
K_{Lambda p}
=
- A R_Lambda A R_p A S
+ A R_{Lambda p} A S
- A R_p A R_Lambda A S
+ A R_p A S_Lambda
+ A R_Lambda A S_p
- A S_{Lambda p}.
```

Consequently, a Frobenius product bound would involve terms such as

```text
||A||_F^3 ||R_Lambda||_F ||R_p||_F ||S||_F,
||A||_F^2 ||R_{Lambda p}||_F ||S||_F,
||A||_F^2 ||R_p||_F ||S_Lambda||_F,
||A||_F^2 ||R_Lambda||_F ||S_p||_F,
||A||_F ||S_{Lambda p}||_F.
```

Such a bound is formally valid under the corresponding finite suprema on the
local region, but those suprema are not currently among the constructive
objects controlled by the package.

### Direction-test verdict

This test gives exact identities, but it does not yet reduce abstraction in a
theorem-worthy way.

What improves:

- `K0'` and `M_{K,Lambda p}` are not mysterious independent objects; they are
  determined by the normalization data `R,S` and their derivatives.

What does not improve:

- the needed bounds are simply transferred to `R^{-1}`, derivatives of `R`,
  and derivatives of `S`;
- the package does not yet provide constructive control of those quantities
  from Schur data, packet data, complement regularity, or the verified
  CoupledBeams matrix;
- without a controlled normalization construction, the `R`-derivative inputs
  are at least as hard as the `K`-derivative inputs they replace.

Therefore this is classified as:

```text
only algebraic rewriting, not promoted.
```

It is useful as an audit result, but not as a new theorem step or
ready-to-use corollary. The next honest bottleneck is derivative control for
the normalization itself: one needs project-usable estimates for `R^{-1}`,
`partial_p R`, `partial_Lambda R`, `partial_{Lambda p}R`, and the
corresponding derivatives of `S`, or a different normalization construction
that eliminates those inputs.

## 22. Gauge-normalization audit after the `R^{-1}S` test

### Candidate gauge: normalize the left factor at `Lambda0`

Let

```text
R0(p) = R(Lambda0,p),
C(p) = R0(p)^(-1).
```

Define the gauge-renormalized Schur block and left factor by

```text
S_hat(Lambda,p) = C(p) S(Lambda,p),
R_hat(Lambda,p) = C(p) R(Lambda,p).
```

Then

```text
S_hat(Lambda,p)
=
R_hat(Lambda,p) ((Lambda-Lambda0)I - K(Lambda,p)),
```

and the exact normalized block `K(Lambda,p)` is unchanged. The gauge has the
base-point normalization

```text
R_hat(Lambda0,p) = I.
```

Consequently,

```text
partial_p R_hat(Lambda0,p) = 0
```

at fixed `Lambda=Lambda0`.

### What simplifies at the base point

In the hatted gauge,

```text
K(Lambda,p)
=
(Lambda-Lambda0)I - R_hat(Lambda,p)^(-1) S_hat(Lambda,p).
```

At `Lambda=Lambda0`, since `R_hat(Lambda0,p)=I`,

```text
S_hat(Lambda0,p) = -K0(p).
```

Differentiating gives the formally simple identity

```text
K0'(p) = - partial_p S_hat(Lambda0,p).
```

Thus the explicit `partial_p R_hat(Lambda0,p)` term has disappeared from the
base-point formula.

However, this is not a free reduction. Expanding the hatted block gives

```text
partial_p S_hat(Lambda0,p)
=
partial_p(R0(p)^(-1) S(Lambda0,p))
=
- R0(p)^(-1) R0'(p) R0(p)^(-1) S(Lambda0,p)
+ R0(p)^(-1) partial_p S(Lambda0,p).
```

Using `S(Lambda0,p)=-R0(p)K0(p)`, this is exactly

```text
partial_p S_hat(Lambda0,p)
=
R0(p)^(-1) R0'(p) K0(p)
+ R0(p)^(-1) partial_p S(Lambda0,p).
```

So the old `R0'(p)` input has not vanished; it is contained in the
`p`-derivative of the chosen gauge `C(p)`.

### Mixed derivative check

The gauge also simplifies some base-point mixed derivative terms, but not the
bottleneck. In hatted notation let

```text
A_hat = R_hat^(-1).
```

The general identity remains

```text
K_{Lambda p}
=
- A_hat R_hat_Lambda A_hat R_hat_p A_hat S_hat
+ A_hat R_hat_{Lambda p} A_hat S_hat
- A_hat R_hat_p A_hat R_hat_Lambda A_hat S_hat
+ A_hat R_hat_p A_hat S_hat_Lambda
+ A_hat R_hat_Lambda A_hat S_hat_p
- A_hat S_hat_{Lambda p}.
```

At `Lambda=Lambda0`, the identities `A_hat=I`, `R_hat_p=0`, and
`S_hat(Lambda0,p)=-K0(p)` reduce this to

```text
K_{Lambda p}(Lambda0,p)
=
- R_hat_{Lambda p}(Lambda0,p) K0(p)
+ R_hat_Lambda(Lambda0,p) partial_p S_hat(Lambda0,p)
- S_hat_{Lambda p}(Lambda0,p).
```

Equivalently, using `partial_p S_hat(Lambda0,p)=-K0'(p)`,

```text
K_{Lambda p}(Lambda0,p)
=
- R_hat_{Lambda p}(Lambda0,p) K0(p)
- R_hat_Lambda(Lambda0,p) K0'(p)
- S_hat_{Lambda p}(Lambda0,p).
```

Thus `R_hat_Lambda` and `R_hat_{Lambda p}` still enter. Moreover, away from
`Lambda0`, the condition `R_hat_p(Lambda0,p)=0` does not imply
`R_hat_p(Lambda,p)=0` on the moving discs.

### Direction-test verdict

The base-point gauge `R_hat(Lambda0,p)=I` is algebraically legitimate and
preserves exact root capture, because it multiplies the Schur block by the
invertible factor `C(p)` depending only on `p`.

It gives one useful bookkeeping simplification:

```text
K0'(p) = -partial_p S_hat(Lambda0,p).
```

But it does not reduce abstraction unless `partial_p S_hat` and
`S_hat_{Lambda p}` are themselves controlled directly. If those hatted
derivatives are expanded back in terms of the original Schur block, the
derivatives of the gauge factor bring back `R0'(p)`. For mixed derivatives,
the hard inputs `R_hat_Lambda` and `R_hat_{Lambda p}` remain.

Therefore this gauge test is classified as:

```text
gauge test failed to reduce abstraction; not promoted.
```

The honest bottleneck is not the absence of a base-point gauge convention.
It is the lack of project-usable derivative control for the normalization
construction itself, or a different reduction strategy that gives such
control without hiding it in the gauge derivative.

## 23. Abstract derivative-line ceiling decision test

### Remaining inputs by layer

After `RUPDC-ER`, the derivative-comparison line has already removed several
formerly abstract inputs.

Already reduced:

- `H_pDelta(p)` is no longer a free scalar hypothesis; `PDCL` and `RUPDC`
  replace it by local matrix controls involving `M_G`, `M_E`, `M_{G,p}`, and
  `M_{E,p}`.
- `M_{G,p}(p)`, `T'(p)`, and `D0'(p)` are no longer independent frozen
  derivative inputs; `FCDL` and `RUPDC-FC` express them through `K0(p)` and
  `K0'(p)`.
- `M_{E,p}(p)` is no longer an unrelated exact-remainder derivative input;
  `ERPC` and `RUPDC-ER` reduce it to the structural mixed-derivative control
  `rho_0(p) M_{K,Lambda p}(p)`.
- the denominator lower bound has a local frozen-discriminant route through
  `FDPL`/`FDCC`, but the actual discriminant margin remains a quantitative
  input to verify in the intended regime.
- the moving-disc and closure constants are organized by `DDMR`, `MDC`, and
  `RUFOSF`, but their useful size is still a local quantitative hypothesis.

Still genuinely independent at the current abstract level:

- `K0'(p)=partial_p K(Lambda0,p)`;
- `M_{K,Lambda p}(p)`, the mixed derivative control for the exact normalized
  matrix on the segment set from `Lambda0` to the moving discs;
- existing local matrix controls such as `M_G`, `M_E`, `M_{E,1}`,
  `B_shift`, `rho_0`, and the exact-denominator closure margin;
- the frozen-discriminant margin or an equivalent frozen-root separation
  lower bound.

If one tries to express `K0'` and `M_{K,Lambda p}` through the normalization
identity, the previous two audits show what happens: the inputs are rewritten
as controls on `R^{-1}`, derivatives of `R`, and derivatives of `S`, or on the
corresponding hatted gauge objects. That is exact algebra, but it is not a
reduction of abstraction unless the normalization construction itself supplies
new derivative bounds.

### Candidate direction 1: stronger abstract normalization gauge

The base-point gauge `R_hat(Lambda0,p)=I` is the natural abstract test case.
It preserves root capture and simplifies the displayed base-point identity to

```text
K0'(p) = -partial_p S_hat(Lambda0,p).
```

However, expanding `partial_p S_hat` brings back the derivative of the gauge
factor, and the mixed derivative formula still contains `R_hat_Lambda` and
`R_hat_{Lambda p}`. Thus a gauge imposed after an arbitrary normalization
does not remove hard derivative inputs; it moves them into the renormalized
block.

A genuinely useful normalization convention would have to come with a
controlled construction of the normalization itself. That is not supplied by
the current abstract package.

### Candidate direction 2: stop the abstract derivative line here

The current abstract line now gives a reusable conditional toolkit:

```text
branch derivative comparison
  -> determinant p-discrepancy from reduced matrices
  -> frozen derivative inputs from K0,K0'
  -> exact-remainder derivative input from M_{K,Lambda p}.
```

One more abstract step would be theorem-worthy only if it removed one of
`K0'`, `M_{K,Lambda p}`, the local matrix controls, or the denominator/closure
margins, or replaced them by objects already controlled earlier in the
package. The two normalization audits show that the obvious abstract
rewritings do not do that.

### Decision

The abstract theorem line has likely reached its natural ceiling here, before
project-specific regime assumptions or a different reduction strategy.

This is not a negative result about the CoupledBeams problem. It is a
classification of the package boundary: the current theorem line has produced
the conditional local derivative-comparison machinery, but it cannot by itself
prove that the remaining derivative and margin inputs have the useful size
needed for an asymptotic or project-specific branch-shift law.

The next honest work is therefore one of:

- impose and justify project-specific regime assumptions that control `K0'`,
  `M_{K,Lambda p}`, the local matrix controls, and the frozen-discriminant
  margin;
- derive those controls directly from a concrete normalization construction;
- use a different reduction strategy that avoids the same normalization
  derivative bottleneck.

This ceiling diagnosis is classified as:

```text
negative audit / ceiling diagnosis, not promoted.
```

It does not add a new lemma or ready-to-use corollary.

## 24. Project-specific quantitative regime setup after the ceiling

### A1 versus A2

The right move after the ceiling diagnosis is `A1`: promote a narrow
conditional regime corollary, not merely a checklist.

The reason is that `RUPDC-ER` is already the endpoint of the abstract
derivative line. The remaining inputs have a natural project-specific shape
and can be checked together on a small interval:

- derivative size of the frozen reduced matrix, through `||K0'(p)||_F`;
- mixed derivative size of the exact normalized matrix, through
  `M_{K,Lambda p}(p)`;
- existing local matrix controls `M_G`, `M_E`, `M_{E,1}`;
- moving-disc closure through `B_shift` and `r`;
- denominator persistence through a frozen denominator or discriminant margin.

Packaging these into `Q_der(I_Q)` is not a proof that the regime holds. It is
still a conditional theorem under explicit regime assumptions. But it is
theorem-useful because later project-specific work can cite one regime
package instead of restating every `RUPDC-ER` input.

### Why the regime corollary is only substitution

Under `Q_der(I_Q)`, the pointwise quantities in `RUPDC-ER` obey

```text
B_shift(p) <= B_Q,
C_Lambda(p) <= C_Q,
rho_0(p) <= |Lambda_fr(p)-Lambda0| + r <= A_* + r,
||K0'(p)||_F <= K_{p,*},
M_{K,Lambda p}(p) <= K_{Lambda p,*}.
```

The strict regime condition

```text
D_Q = c0 - C_Q - 2 B_Q > 0
```

protects the exact denominator in the derivative comparison. The frozen
numerator in the second term is bounded by

```text
(sqrt(2)|Lambda_fr(p)-Lambda0| + ||K0(p)||_F)||K0'(p)||_F
<=
(sqrt(2)A_* + K_{0,*})K_{p,*}.
```

Substituting these bounds into `RUPDC-ER` gives the statement in
`lemma_statements.md`.

### Honest boundary

The regime corollary reduces future repetition, but it does not reduce the
project-specific work itself. The next job is still to verify, estimate, or
numerically audit the regime quantities for a concrete CoupledBeams interval
and branch choice.

In particular, this step does not prove:

- that the bounds in `Q_der(I_Q)` hold for any actual branch;
- that the derivative-comparison estimate is asymptotically small;
- a symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- a final veering criterion or branch-case application theorem.

## 25. Branch-ready verification template after `QDR`

### Candidate-regime search result

The project materials do contain one sufficiently grounded candidate seed for
a future `Q_der` verification attempt, so the choice here is `C1`.

The candidate seed is:

```text
fixed project regime: beta = 15 deg, physical radius = 5 mm,
target tracked branch: bending_desc_04,
paired tracked branch seed: bending_desc_01,
candidate parameter window: 0.20 <= mu <= 0.35.
```

This is grounded in the current `docs/veering/` materials:

- `candidate_cases.md` identifies `bending_desc_04` as the strongest current
  candidate because it shows staged local modal-order reorganization;
- `paired_branch_candidates_beta15_r5.csv` lists
  `bending_desc_04` versus `bending_desc_01` on `mu=0.20..0.35` as the
  nearest available low tracked branch during the first staged local-order
  reorganization of `bending_desc_04`;
- `bending_desc_04_pair_search_beta15_r5.csv` records an interior gap minimum
  near `mu=0.26625`, with

  ```text
  Lambda_target = 6.2086633731319516,
  Lambda_pair   = 5.4876154692852506,
  gap_abs       = 0.72104790384670103,
  gap_rel       = 0.123295265709942.
  ```

- `strict_veering_decision_table_beta15_r5.csv` explicitly labels this as
  `exchange of modal character without strict veering`, with missing pairwise
  MAC / branch-pair exchange evidence.

Thus this is not a verified theorem-line packet, not a verified close
two-root cluster, and not evidence of strict veering. It is only the best
currently grounded local candidate seed for trying the `Q_der` verification
workflow.

The word `r` in the theorem-line moving disc must not be confused with the
physical beam radius `5 mm` in the candidate data.

### Branch-ready verification record

For this candidate, or for any later candidate replacing it, the verification
record should be filled as follows.

1. Candidate selection data.
   Record the candidate id, fixed project parameters, target branch, paired
   branch seed, interval `I_Q`, base parameter `p0`, base roots, and local
   spectral window. For the current C1 seed, `p0` would naturally be chosen
   near the documented interior gap minimum, but the actual theorem-line base
   point and spectral window still have to be fixed by the verifier.

2. Local construction checks.
   Construct or specify the retained two-root packet, verify complement
   regularity on the local spectral window, build the Schur block, verify the
   exact local `Lambda`-normalization, and ensure the moving discs and `ERPC`
   segment set remain inside the normalization region.

3. Frozen-model checks.
   Choose the frozen branch `Lambda_fr(p)` and local second frozen branch or
   discriminant route. Verify the frozen denominator margin `c0`, the frozen
   simple-root lower bound `m0`, the theorem-line moving-disc radius, and the
   closure condition.

4. Exact-vs-frozen discrepancy checks.
   Estimate `M_G`, `M_E`, and `M_{E,1}` on the moving discs. Then compute

   ```text
   B_Q = [G_* E_* + (1/2) E_*^2] / m0,
   C_Q = sqrt(2) E_* + (G_* + E_*) E_{1,*},
   D_Q = c0 - C_Q - 2 B_Q.
   ```

   The checks `B_Q < r` and `D_Q > 0` must both pass.

5. Derivative-level checks.
   Estimate `||K0'(p)||_F`, `M_{K,Lambda p}(p)`, `||K0(p)||_F`, and the
   offset `|Lambda_fr(p)-Lambda0|` uniformly on `I_Q`. These supply
   `K_{p,*}`, `K_{Lambda p,*}`, `K_{0,*}`, and `A_*`.

6. Final decision.
   Mark `Q_der(I_Q)` as verified only if every construction, margin,
   discrepancy, derivative, closure, and denominator check has passed. If any
   field is missing, the record should name the first missing or failed item
   and no `QDR` conclusion should be drawn.

### Filled C1 verification record from current repo evidence

Use the following labels:

```text
passed from repo evidence
partially grounded
not yet available
blocked by missing quantitative bound
blocked by missing construction
fails on current evidence
not checked
```

This record uses only the current repository evidence. It does not upgrade
tracked branch data into a theorem-ready packet.

1. Candidate selection data.

   - `candidate_id = C1`:
     `partially grounded`. `C1` is the package-local name for the candidate
     seed selected after `QDR`; it is not an external branch id in the data.
   - Fixed project parameters:
     `passed from repo evidence`. `docs/veering/README.md`,
     `docs/veering/candidate_cases.md`, and
     `docs/veering/data/dataset_index.md` all record the baseline analysis
     regime `beta = 15 deg`, physical radius `5 mm`, and `0 <= mu <= 0.9`.
   - Target tracked branch `bending_desc_04`:
     `passed from repo evidence` for project candidate selection.
     `candidate_cases.md` calls it the strongest current candidate because of
     staged local modal-order reorganization; `dataset_index.md` records full
     tracked-table coverage and selected states for this branch.
   - Paired tracked branch seed `bending_desc_01`:
     `partially grounded`. `paired_branch_candidates_beta15_r5.csv` lists
     `bending_desc_04` versus `bending_desc_01` on `mu=0.20..0.35`, but marks
     the status as `missing_paired_branch_evidence`.
   - Parameter interval `I_Q`:
     `partially grounded`. The candidate window `0.20 <= mu <= 0.35` is
     explicit in `paired_branch_candidates_beta15_r5.csv`,
     `bending_desc_04_pair_search_beta15_r5.csv`, and
     `strict_veering_decision_table_beta15_r5.csv`. It is grounded as a
     project scan window, not yet as a theorem-ready interval where all
     `Q_der` hypotheses are known.
   - Base parameter:
     `partially grounded`. `paired_gap_scan_beta15_r5.csv` gives the smallest
     recorded `bending_desc_04`/`bending_desc_01` gap at

     ```text
     mu = 0.26624999999999999.
     ```

     This is the C1 candidate-level anchor fixed by the construction test
     below, but it is not yet a theorem base point certified by a retained
     packet and local spectral window.
   - Base roots:
     `partially grounded`. At that grid point, the tracked data give

     ```text
     Lambda_target = 6.2086633731319516,
     Lambda_pair   = 5.4876154692852506,
     gap_abs       = 0.72104790384670103,
     gap_rel       = 0.12329526570994202.
     ```

     These are tracked-branch values, not yet roots of a constructed local
     Schur block or frozen model.
   - Local spectral window:
     `not yet available` at theorem-ready level. The construction test below
     records a candidate-level spectral-window target anchored by the tracked
     sorted-position `2/3` pair, but the repo does not yet specify certified
     open endpoints, guard margins, a moving-disc radius, or a verified
     retained two-root window for C1.

2. Local construction checks.

   - Isolated cluster / retained pair evidence:
     `partially grounded`. The paired gap scan gives adjacent sorted
     positions `3` and `2` at the minimum-gap grid point, but the gap is
     explicitly described as large and not supportive of a close avoided
     crossing. This is project candidate evidence only.
   - Packet construction:
     `blocked by missing construction`. No C1 right/left packet `V,W`, basis
     completion, or retained two-dimensional spectral packet is recorded in
     the allowed repo materials.
   - Complement / Schur-readiness:
     `blocked by missing construction`. No C1 complement choice, determinant
     lower bound for `D`, or Schur block `S` is recorded.
   - Exact local `Lambda`-normalization:
     `blocked by missing construction`. Without a C1 Schur block, the
     normalization objects `R` and `K` are not available.
   - Theorem-ready retained packet:
     `not yet available`. The current data support tracked candidate
     histories, not a theorem-ready retained packet.

3. Frozen-model checks.

   - Frozen branch `Lambda_fr(p)`:
     `not yet available`. The tracked branch `bending_desc_04` is not yet a
     frozen-model branch of a constructed `K0(p)`.
   - Second frozen root / discriminant / denominator margin:
     `not yet available`. No C1 frozen discriminant `Disc_fr`, second frozen
     root branch, or lower bound `c0` is recorded.
   - Frozen simple-root lower bound `m0`:
     `not yet available`.
   - Theorem-line moving-disc radius `r`:
     `not yet available`. The physical radius `5 mm` in the project data is
     not the theorem-line moving-disc radius.
   - Moving-disc closure:
     `blocked by missing quantitative bound`. `B_Q` cannot be computed until
     `m0`, `G_*`, and `E_*` are available.

4. Matrix and derivative bounds.

   - `M_G`, `M_E`, `M_{E,1}`:
     `blocked by missing construction`. These require the local exact
     normalized block, frozen model, matrix remainder, and moving discs.
   - `||K0(p)||_F` and `||K0'(p)||_F`:
     `blocked by missing construction`. `K0(p)` has not been constructed for
     C1.
   - `M_{K,Lambda p}`:
     `blocked by missing construction`. The normalized matrix
     `K(Lambda,p)` and its mixed derivative domain are not available for C1.
   - Offset bound `A_*`:
     `not yet available`. The tracked branch offset relative to a theorem-line
     `Lambda0` cannot be fixed before the base point and frozen branch are
     chosen.
   - `B_Q`, `C_Q`, and `D_Q`:
     `blocked by missing quantitative bound`. None of these constants can be
     computed from the current C1 project evidence.

5. Evidence that fails or weakens a stronger interpretation.

   - Close-pair / strict interaction support:
     `fails on current evidence`. `bending_desc_04_pair_assessment.md` and
     `strict_veering_decision_table_beta15_r5.csv` both say the
     `bending_desc_04`/`bending_desc_01` gap is large, with no pairwise MAC or
     branch-pair exchange evidence.
   - Pairwise MAC:
     `not yet available`. `pair_mac_beta15_r5.csv` records `NA` values and
     notes that pairwise MAC is not available from existing data.

6. Final decision.

   ```text
   Q_der(I_Q) cannot yet be decided for C1 from current repo evidence.
   ```

   The first blocking item is the missing theorem-ready retained packet and
   certified local spectral window for C1. Until that construction exists,
   the Schur block, normalized block, frozen model, discriminant margin,
   matrix-control constants, derivative constants, and closure checks cannot
   be evaluated.

The next concrete action is therefore not to claim `Q_der`, but to use the
candidate-level construction target below to build the retained two-root
packet, certify the local spectral window, and test complement regularity.
Only after that can the `Q_der` constants be estimated.

### Promotion decision

This template is promoted as `QVR` because it is not merely editorial. It is
the branch-ready gateway between the abstract regime package `Q_der(I_Q)` and
any future branch-case work: if the record passes, then `QDR` applies; if it
does not pass, the package says exactly where the verification failed.

It still does not verify the C1 candidate, construct its packet, prove a
discriminant margin, prove derivative bounds, or establish any branch-case
theorem.

### C1 theorem-line candidate construction test

Outcome: `T1`, but only at candidate-record level.

The repo contains enough grounded information to turn C1 from a tracked
candidate seed into a theorem-line candidate object. This is strictly weaker
than a packet proof, complement regularity, a certified local window, or
`Q_der` verification.

1. Candidate base parameter.

   The grounded candidate base point is

   ```text
   mu0 = 0.26624999999999999.
   ```

   This is not an arbitrary midpoint of `0.20 <= mu <= 0.35`. It is the grid
   point where `paired_gap_scan_beta15_r5.csv` records the smallest stored
   `bending_desc_04`/`bending_desc_01` gap in the C1 window. The supporting
   row gives

   ```text
   Lambda_pair   = 5.4876154692852506,
   Lambda_target = 6.2086633731319516,
   gap_abs       = 0.72104790384670103,
   gap_rel       = 0.12329526570994202,
   sorted_position_pair   = 2,
   sorted_position_target = 3.
   ```

   This fixes the candidate anchor for construction. It does not prove that
   the two roots form an isolated theorem packet.

2. Candidate spectral-window target.

   The grounded spectral target is the local two-root window around the
   sorted-position `2/3` pair at `mu0`, anchored by the tracked values

   ```text
   5.4876154692852506 <= Lambda <= 6.2086633731319516
   ```

   as the inner pair bracket. This is more specific than "near the
   interaction": it says which two tracked roots the future retained window
   must contain and at which base parameter.

   It is still not a certified theorem-line spectral window. The repo does
   not provide open lower/upper endpoints, guard margins to the adjacent
   sorted roots, a proof that exactly these two roots remain inside the
   window after shrinking, or the moving-disc radius used later in `Q_der`.

3. Retained two-root packet construction target.

   The retained packet construction target is the two tracked branches

   ```text
   bending_desc_01 at sorted position 2,
   bending_desc_04 at sorted position 3,
   ```

   over the documented C1 scan window `0.20 <= mu <= 0.35`, with the
   construction anchored at `mu0`. Around the anchor, the copied tracked table
   records both branches with local half-wave descriptor `(1,2)` and high
   along-track assignment MAC. This supports using the pair as a candidate
   retained two-root object.

   The support is still heuristic at theorem-line level.
   `pair_mac_beta15_r5.csv` records no pairwise MAC/cross-MAC data, and the
   assessment files say the C1 gap is large and not evidence of a close
   two-branch interaction. Thus the target is usable for the next
   construction attempt, but it is not a proved retained packet.

4. Candidate-object status.

   C1 is now recorded as a theorem-line candidate object:

   ```text
   parameter anchor: mu0 = 0.26624999999999999,
   spectral target: sorted-position 2/3 pair bracketed by
     Lambda_pair = 5.4876154692852506 and
     Lambda_target = 6.2086633731319516,
   retained-packet target:
     bending_desc_01 / bending_desc_04.
   ```

   The next honest step is to construct the right/left retained pair at this
   candidate object, choose certified open spectral endpoints, and prove the
   complement/Schur-readiness checks. Only after that can the C1 `QVR` record
   move from candidate target to actual `Q_der` verification.

## 26. Ideal-model algebra stays downstream

The earlier `ideal symmetric 2x2` algebra remains useful only as an ideal
model. It must not be transferred to the CoupledBeams bridge line until a
later symmetric or self-adjoint reduced normal form is actually justified.

Therefore the current package does not claim:

- project-defined `delta(p)` and `kappa(p)`;
- ideal gap formulas for the actual reduced block;
- a direct bridge from the present reduced block to a final veering
  criterion.

## 27. Constructive remarks

The most realistic constructive packet path remains:

1. isolate a two-root cluster;
2. build pointwise right/left root data in the simple-root setting;
3. promote those data to rank-2 packets;
4. test complement regularity;
5. only then use Schur reduction and normalization.

This is still a constructive program, not a theorem that any named packet in
CoupledBeams already satisfies the hypotheses.

## 28. Current next obstruction

The bridge line has reached exact local reduced-block statements, local exact
normalization, controlled freezing, local parameterwise nearby-root
selection, exact-branch continuation, a conditional local quantitative
comparison estimate between the exact and frozen branches, and a constructive
determinant-discrepancy bound from the matrix remainder, a closure check for
the moving discs, a conditional local first-order root-shift formula, and a
ready-to-use local first-order shift corollary combining the constructive
denominator and derivative controls, plus a conditional parameterwise
derivative-comparison step and constructive `p`-derivative discrepancy
control from reduced-model objects, and frozen coefficient-derivative control
for the `RUPDC` frozen inputs, plus exact-remainder `p`-derivative control
through a mixed derivative of `K`. The normalization-derivative audit above
shows that rewriting `K0'` and `M_{K,Lambda p}` through `R` and `S` is exact
but does not by itself decrease abstraction. The gauge-normalization audit
also shows that imposing `R_hat(Lambda0,p)=I` improves bookkeeping but hides
the same derivative difficulty in `partial_p S_hat` and leaves mixed
`R_hat`-derivatives in the problem. The ceiling test therefore leaves the
abstract derivative line at the following boundary:

```text
the abstract theorem line has likely reached its natural ceiling before
project-specific regime assumptions or a different reduction strategy.
```

The next real obstruction is to turn the conditional parameterwise derivative
comparison into a useful asymptotic or project-specific estimate by controlling
`K0'`, `M_{K,Lambda p}`, the existing local matrix and closure constants, and
the frozen-discriminant margin. The new regime package `Q_der(I_Q)` turns
this list into a single checkable local assumption set, and the `QDR`
corollary makes `RUPDC-ER` directly citeable under that set. The branch-ready
`QVR` protocol now says how to test one grounded candidate record, with the
current C1 seed being `bending_desc_04` versus `bending_desc_01` on
`beta=15 deg`, physical radius `5 mm`, and `0.20 <= mu <= 0.35`. None of this
verifies the regime. Trying to push the remaining inputs through `R^{-1}S` or
a base-point gauge requires controlling the normalization derivatives
themselves: `R^{-1}`, `partial_p R`, `partial_Lambda R`,
`partial_{Lambda p}R`, and the corresponding `S`-derivatives.

That step must be handled before any honest move to:

- symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- final veering criteria;
- branch-case application theorems.

## 29. Freeze decision

The current veering-bridge theorem line is temporarily frozen here, pending
broader regime evidence and additional examples.

### Why freeze now

The abstract derivative/reduction line has likely reached its natural ceiling.
The remaining bottlenecks are no longer plausibly resolved by one more
abstract lemma. They are candidate-dependent inputs: packet construction,
certified spectral windows, complement regularity, frozen margins, `K0'(p)`,
`M_{K,Lambda p}(p)`, and the local closure constants. Pushing the same strict
line further right now would mostly rename these unresolved inputs rather than
turn them into usable project evidence.

Cheap feasibility questions now look more informative than further abstract
tightening: which parameter regimes actually exhibit clean positive and
negative veering-like behavior, where `Q_der(I_Q)` might plausibly close, and
which examples fail early enough to rule out overinvestment in the wrong local
picture.

### What remains useful

This package is not abandoned. It remains useful as the canonical record of:

- the layer discipline for the bridge line;
- the separation between strong theorem claims, conditional mechanism claims,
  and candidate-level records;
- what is exact, approximate, conditional, or still unresolved;
- the framework to revisit later once better evidence exists.

### What would justify reopening

An honest return to this line would need one or more of:

- more parameter examples, with both positive and negative veering-like
  cases;
- clearer beta-regime evidence and at least one named candidate with enough
  data to decide `Q_der(I_Q)`;
- a cheap feasibility or predictor criterion that helps sort promising local
  regimes;
- a different reduction or normalization strategy that genuinely decreases the
  current derivative-control burden.
