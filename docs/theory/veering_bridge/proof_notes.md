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

## 16. Ideal-model algebra stays downstream

The earlier `ideal symmetric 2x2` algebra remains useful only as an ideal
model. It must not be transferred to the CoupledBeams bridge line until a
later symmetric or self-adjoint reduced normal form is actually justified.

Therefore the current package does not claim:

- project-defined `delta(p)` and `kappa(p)`;
- ideal gap formulas for the actual reduced block;
- a direct bridge from the present reduced block to a final veering
  criterion.

## 17. Constructive remarks

The most realistic constructive packet path remains:

1. isolate a two-root cluster;
2. build pointwise right/left root data in the simple-root setting;
3. promote those data to rank-2 packets;
4. test complement regularity;
5. only then use Schur reduction and normalization.

This is still a constructive program, not a theorem that any named packet in
CoupledBeams already satisfies the hypotheses.

## 18. Current next obstruction

The bridge line has reached exact local reduced-block statements, local exact
normalization, controlled freezing, local parameterwise nearby-root
selection, exact-branch continuation, a conditional local quantitative
comparison estimate between the exact and frozen branches, and a constructive
determinant-discrepancy bound from the matrix remainder, a closure check for
the moving discs, and a conditional local first-order root-shift formula. The
next real obstruction is:

```text
turn the first-order shift formula into a sharper asymptotic or derivative-
level law by upgrading the coefficient-level frozen discriminant checks into
an actual project-specific regime verification.
```

That step must be handled before any honest move to:

- symmetric or self-adjoint reduced normal form;
- project-defined `delta-kappa`;
- final veering criteria;
- branch-case application theorems.
