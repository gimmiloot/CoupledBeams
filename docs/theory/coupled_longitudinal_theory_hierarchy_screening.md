# Coupled longitudinal-theory hierarchy: bounded screening

Date: 2026-10-07. This diagnostic compares independently sorted spectral
positions of three linear frame models. Mindlin–Herrmann (MH) is the
**reference model for comparison**, not an exact solution or experimental truth.
No applicability threshold, safe prefix or production-model replacement is
inferred. The numerical tables below are evidence for one geometry only.

## Question, common inputs and model contracts

How much do the first 12 coupled frequencies change when elementary or
planar Rayleigh–Love axial motion is replaced by the selected Jang-type MH
block, while retaining the same Timoshenko bending model and frame geometry?

The accepted G20 contract is `tests/data/reddy_four_ply_isotropic_limit_cases.json`:
E=1, rho=1, nu=0.3, b=0.20 along y, h=0.05 along z, A=bh,
Iy=bh^3/12, total L=1, L1=L2=0.5, kappa=5/6. These are normalized diagnostic
inputs; frequency values are f*=f L/sqrt(E/rho). They are not steel frequencies.
The eight fixed beta values are 0,5,15,30,45,60,75,90 degrees.

All models use the **identical** verified Timoshenko basis, with
S=kappa GA, B=EIy and rotary mass rho Iy. There is no fitting or change of
section, shear coefficient, rotary inertia or physical centroidal clamp.

| Model | Axial energy, with the common factor 1/2 integral dx omitted | Independent DOFs |
| --- | --- | --- |
| E | T=rho A u_t^2; V=EA u_x^2 | u |
| RL | T=rho A u_t^2+J u_xt^2; V=EA u_x^2; J=nu^2 rho Iy, H=0 | u |
| MH | T=rho A u_t^2+rho Iy c_t^2; V=C(u_x^2+2nu u_x c+c^2)+H c_x^2; C=EA/(1-nu^2), H=kappa G Iy | u,c |

The common bending terms are T=rho A w_t^2+rho Iy theta_t^2 and
V=EIy theta_x^2+kappa GA(w_x-theta)^2. These self terms are unchanged.
The MH normal block is the selected reduced Jang closure; Ng/Fernandes's
Lamé variant and the Rucka/Jang source-specific variants are not substituted.
Bishop stays a standalone reference/diagnostic theory.

This extension has no new frozen block in `equations.tex` or
`src/my_project/analytic/`: their original determinants/signs/order are
unchanged. The reduced axial kernel is the independently checked H=0 branch
of `bishop_longitudinal.py`; the planar J uses Iy, **not** the polar Ip of
the circular Bishop literature cases. The existing single-rod hierarchy
uses the same planar specialization. Beta0 direct/split checks below establish
its consistency with the common frame assembly. No theory/code mismatch was
found. SymPy/Lean are unavailable in the working Python environment;
variation is shown directly and exact Fraction/array/limit tests are used.

## Rayleigh–Love resultant and boundary limitation

In delta(T−V), integrating the J term first in x, then in time, gives
the spatial endpoint contribution

\[
 [(-EAu_x-Ju_{xtt})\,\delta u]_0^L.
\]

Thus physical traction is N=EAu_x+Ju_xtt. With exp(i omega t),

\[
 (EA-J\omega^2)U''+\rho A\omega^2U=0,
 \qquad N=(EA-J\omega^2)U'.
\]

N has force units; J omega^2 has the units of EA. The natural endpoint
force is sigma N, sigma=−1 at x=0 and +1 at x=Li. There is no separate
u_x essential constraint. H=0 is handled as a second-order harmonic problem,
not as a fourth-order model with a small artificial H.

Both positive local coordinates run from external clamp to joint. For A/B
the external essential set is u=w=theta=0. C additionally imposes c=0.
These clamps share the physical centroid/bending constraints but **do not
have identical resolved cross-sectional contraction kinematics**. A/B also
omit the MH scalar c-continuity/R-balance pair at the joint. Consequently this
screening compares the declared reduced theories and their associated
boundary closures, not an isolated bulk-coefficient perturbation.

## Geometry and virtual work

The unchanged project convention is beta = signed angle of the
joint-to-right-clamp ray from +EX toward +EY. With local x toward the joint,

\[
 t_1=(1,0),\quad n_1=(0,-1),\qquad
 t_2=(-\cos\beta,-\sin\beta),\quad n_2=(-\sin\beta,\cos\beta).
\]

d_i=u_i t_i+w_i n_i. For reduced local q=(u,w,theta), global
qg=(dX,dY,theta), G=diag([t,n],1) and qg=G q.
The dual force is fg=G p, p=(N,Q,M), since
p^T delta q=(G p)^T delta qg. Both joint endpoint signs are +1.
The six comparator conditions are d1=d2, theta1=theta2,
sum F_node=0 and sum M_node=0. There is no fictitious c or R row.

MH uses the unchanged eight-condition common-frame closure, adding c1=c2
and R1+R2=0. **c1=c2 is variational reduced-frame closure, not a direct
3D-elasticity derivation for a finite welded joint region.** Local arms remain
block separated; nonzero-angle translations create structural coupling.
Reflection reverses local w,theta,Q,M but preserves c/R; no new convention
was introduced. Sources and their limited roles remain in the
[canonical joint note](mindlin_herrmann_timoshenko_rigid_joint.md).

## Solver, count certificate and gates

One small comparator layer assembles a 12x12 full boundary system from the
existing H=0 cosine/sine axial basis and the unchanged bounded Timoshenko
basis. The latter uses anchored exponentials, not unbounded hyperbolics or
long products of transfer matrices. MH retains its 16x16 production system.

The independent root count is the fixed-arm Dirichlet count J0 plus the
negative inertia of the energy-derived nodal Schur matrix. Tim arm-pole
catalogs are certified by saturation of the existing min-max count bound;
elementary/RL arm poles have exact fixed-fixed formulas. Positive fixed
congruences balance the Schur matrix without changing its inertia. Raw
symmetry, distance of its eigenvalues from zero and arm conditioning are
checked before using roundoff-only symmetrization.

For beta0 equal halves, a genuine global elementary/RL root can coincide
with an arm Dirichlet pole. A count query within 1e−7 max(1,omega) of a pole
is moved to its right by twice that exclusion distance. Requested and
effective nodes are saved. Counts and determinant brackets use the same
effective endpoints. **Counts on opposite sides are not forced equal**.
The full analytic boundary matrix remains regular at the pole and determines
the unshifted physical root. This is a query-conditioning measure, not a
root correction or exclusion. A dedicated test retains such roots.

Initial ceiling f*=3.75, 400 deterministic intervals, root xtol=1e−11,
rtol=1e−12. One count-guided bisection mechanism is bounded by depth20 and
2000 subdivisions. Only absence of guard13 permits one ceiling multiplier
1.5; no expansion was needed. No angle refinement is permitted. Certified
complete inventories can exceed the scientific prefix; only positions1–12
enter metrics, position13 is a guard. Near multiple/unresolved roots fail
locally rather than being duplicated or silently removed.

Pre-existing general/single-rod tolerances are retained: boundary/joint/ODE
1e−9, internal total-energy quotient 5e−8, mass Gram 5e−7, nonzero singular
condition <1e8. No physical-difference or overlap threshold is a PASS gate.
Internal energy verifies the eigenproblem and is not used for mode labels.

Before screening, A/B each pass direct fixed-fixed beta0 recovery on splits
0.50 and 0.35/0.65; same-frame duality, rank6, beta45 arm swap and reflection
including remapped kinematic profiles, root counts and residuals. All roots
of these bounded gate inventories are retained. MH reuses immutable
general-angle results at beta5/45/90 and beta0 roots; new MH inventories at
15/30/60/75 use its unchanged solver. Code/source/config/artifact hashes are
validated. Beta0 MH profiles are freshly evaluated by the unchanged basis,
without root search, and its first13 are independently count-certified.

The local `frequency-map-v1` instance uses `certified_audit` because the user
explicitly requested comparator gates, completeness and fixed-beta shape
overlaps. Its finite grid, guard, bounded repair and sorted semantics are in
the config. This does not redefine the ordinary `fast_plot` default.

## Comparison and shape diagnostics

For independently sorted k,
delta_E=(fE_k−fMH_k)/fMH_k and delta_RL=(fRL_k−fMH_k)/fMH_k.
Signed and absolute values are saved. Sorted index does not assert physical
mode identity across theories or angles. No roots are reassigned.

At each **fixed beta**, all three theory pairs have complete 12x12 matrices

\[
 O_d=\frac{|\sum_i\int d_i^a\cdot d_i^b dx|^2}
 {\sum_i\int|d_i^a|^2dx\;\sum_i\int|d_i^b|^2dx},\qquad
 O_\theta=\frac{|\sum_i\int\theta_i^a\theta_i^bdx|^2}
 {\sum_i\int|\theta_i^a|^2dx\;\sum_i\int|\theta_i^b|^2dx}.
\]

Geometric displacement and rotation are kept separate because their units
differ. Sign and amplitude are immaterial. Common Gauss order200 is used;
order300 convergence is tested at beta0/45/90 from saved profiles without
new root calculations. Theta norms below
128 eps_machine cond_nonzero ||d||/L receive SMALL_NORM and associated
overlaps are null/NOT_INFORMATIVE. Norms and numeric zero scales are saved.
All row/column argmax, diagonals, best/second/margins are saved; off-diagonal
argmax is a diagnostic, not a failure. There is no arbitrary MAC threshold,
assignment, continuation or across-beta shape comparison.

For MH only,

\[
 D_c=\frac{\|c+\nu u'\|}{\|c\|+\|\nu u'\|}.
\]

All three norms are retained. Zero denominator uses NOT_DEFINED_SMALL_FIELD,
not D_c=0. The derivative-field zero scale multiplies the theta roundoff
scale by max(1,L max sqrt(|k_spatial^2|)); this includes evanescent spatial
amplification and is not a physical classification threshold. Triangle
inequality gives 0<=D_c<=1 for a defined field. D_c measures departure from
quasistatic Poisson contraction; independent dynamics and end/joint layers
can contribute. It is not a modal type or an energy fraction.

Adjacent sorted gaps (f_{k+1}−f_k)/f_k, k1–11, are saved only as context.
They do not establish veering, avoided crossings or internal resonance.

## Numerical results and status

Canonical local bundle: `results/coupled_longitudinal_theory_hierarchy_screening/9a35ff23c3f43c66/`.

| Status | Result |
| --- | --- |
| ELEMENTARY_TIM_COUPLED_COMPARATOR | PASS |
| RAYLEIGH_LOVE_TIM_COUPLED_COMPARATOR | PASS |
| MHTIM_REFERENCE_REUSE | PASS |
| COUPLED_HIERARCHY_ROOT_INVENTORY | PASS |
| COUPLED_HIERARCHY_FIXED_BETA_SHAPE_OVERLAP | COMPLETE |
| COUPLED_HIERARCHY_SCREENING | COMPLETE |

Comparator direct/split maximum relative differences are 7.00e-12 (E)
and 7.00e-12 (RL); swap frequency differences <=2.23e-16, reflection 0,
remapped kinematic L2 <=2.83e-14. Geometry/virtual-work maximum2.23e-16;
rank6 at all algebraic test angles. Both gates pass before screening.

| beta (deg) | max abs delta_E (%) [k] | max abs delta_RL (%) [k] | median E / RL (%) | non-diagonal E / RL rows | max D_c [k] |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.141030 [5] | 0.175141 [11] | 0.000000 / 0.000000 | 0 / 0 | 0.058849 [5] |
| 5 | 0.140860 [5] | 0.175014 [11] | 0.000383 / 0.000403 | 0 / 0 | 0.085301 [7] |
| 15 | 0.139505 [5] | 0.173993 [11] | 0.003492 / 0.003698 | 0 / 0 | 0.085310 [7] |
| 30 | 0.134935 [5] | 0.170480 [11] | 0.013977 / 0.015114 | 0 / 0 | 0.085340 [7] |
| 45 | 0.127346 [5] | 0.164384 [11] | 0.025681 / 0.028535 | 0 / 0 | 0.085393 [7] |
| 60 | 0.116798 [5] | 0.155304 [11] | 0.023991 / 0.034166 | 0 / 0 | 0.085474 [7] |
| 75 | 0.103426 [5] | 0.142606 [11] | 0.030618 / 0.035726 | 0 / 0 | 0.085590 [7] |
| 90 | 0.087503 [5] | 0.125413 [11] | 0.038626 / 0.046567 | 0 / 0 | 0.085751 [7] |

Global maxima: E/MH0.1410296214% at beta0,k5; RL/MH0.1751407056%
at beta0,k11. Both signed deviations are negative. These maxima are over
the eight declared angles and first12 only; no continuous-angle extremum
was searched. RL does not reduce the maximum MH difference for this input.

All displacement row/column best correspondences of all three pairs are
diagonal. All informative theta best correspondences are also diagonal.
Thus there are **no non-diagonal cases to list** or to associate with close
gaps. Minimum diagonal O_d: E/MH0.9998862625, RL/MH0.9998031367.
Minimum adjacent gaps over the declared grid occur at beta75,k9->10:
E0.0128691744, RL0.0130729679, MH0.0126369418. These are diagnostic gaps,
not evidence of veering or an internally resonant pair.

Defined D_c ranges approximately0.04162--0.08575; the exact three norms and
ratios are in contraction.csv. At beta0, positions1,2,3,4,6,7,9,10,12 have
numerically negligible u/c fields and receive NOT_DEFINED_SMALL_FIELD.
No modal type is assigned. The largest defined D_c is0.0857507021 at
beta90,k7; this can reflect independent contraction and end/joint layers.
It does not imply a frequency error or an applicability failure.

### Representative independently sorted frequencies

All values below are normalized f*. Full precision, signed differences and
guard brackets are in the local JSON/CSV. Equal displayed decimals do not
imply bitwise equality.

| k | beta0 E / RL / MH | beta45 E / RL / MH | beta90 E / RL / MH |
| --- | --- | --- | --- |
| 1 | 0.050527603 / 0.050527603 / 0.050527603 | 0.136050322 / 0.136050318 / 0.136051160 | 0.134718488 / 0.134718465 / 0.134723377 |
| 2 | 0.136325767 / 0.136325767 / 0.136325767 | 0.155284476 / 0.155283886 / 0.155336543 | 0.185619381 / 0.185619176 / 0.185632402 |
| 3 | 0.260072546 / 0.260072546 / 0.260072546 | 0.307452527 / 0.307448151 / 0.307572760 | 0.382999535 / 0.382992492 / 0.383164056 |
| 4 | 0.416258846 / 0.416258846 / 0.416258846 | 0.409712693 / 0.409710412 / 0.409759369 | 0.396315445 / 0.396301324 / 0.396591527 |
| 5 | 0.500000000 / 0.499953743 / 0.500706144 | 0.500496085 / 0.500454253 / 0.501134260 | 0.501973417 / 0.501944523 / 0.502413043 |
| 6 | 0.599831008 / 0.599831008 / 0.599831008 | 0.586725005 / 0.586716876 / 0.586829718 | 0.559613924 / 0.559595393 / 0.559866141 |
| 7 | 0.806045922 / 0.806045922 / 0.806045922 | 0.808215283 / 0.808213020 / 0.808241880 | 0.817209868 / 0.817198857 / 0.817337288 |
| 8 | 1.000000000 / 0.999630095 / 1.001229213 | 0.934807248 / 0.934712804 / 0.935209840 | 0.905449217 / 0.905419729 / 0.905590601 |
| 9 | 1.030751272 / 1.030751272 / 1.030751272 | 1.152634627 / 1.152330521 / 1.153232146 | 1.241167515 / 1.241028385 / 1.241571900 |
| 10 | 1.270446989 / 1.270446989 / 1.270446989 | 1.265246752 / 1.265219077 / 1.265318983 | 1.270361657 / 1.270048723 / 1.270814792 |
| 11 | 1.500000000 / 1.498752436 / 1.501381967 | 1.491830619 / 1.490674059 / 1.493128522 | 1.463850642 / 1.463012640 / 1.464849755 |
| 12 | 1.522253495 / 1.522253495 / 1.522253495 | 1.519295732 / 1.519126554 / 1.519478634 | 1.511904575 / 1.511327018 / 1.512533709 |

| beta | guard13 E | guard13 RL | guard13 MH | certified inventory E / RL / MH |
| --- | --- | --- | --- | --- |
| 0 | 1.783833632 | 1.783833632 | 1.783833632 | 23 / 23 / 13 |
| 5 | 1.783886324 | 1.783885676 | 1.783887293 | 23 / 23 / 23 |
| 15 | 1.784310285 | 1.784304450 | 1.784319027 | 23 / 23 / 23 |
| 30 | 1.785773507 | 1.785750055 | 1.785808680 | 24 / 24 / 24 |
| 45 | 1.788326970 | 1.788273812 | 1.788406840 | 24 / 24 / 24 |
| 60 | 1.792156912 | 1.792061539 | 1.792300587 | 24 / 24 / 24 |
| 75 | 1.797554008 | 1.797403566 | 1.797781418 | 24 / 24 / 24 |
| 90 | 1.804946885 | 1.804728638 | 1.805278177 | 24 / 24 / 24 |

MH beta0 certifies the reused prefix13 up to its saved guard bracket,
rather than claiming a complete MH inventory to3.75 from the older shorter
bending catalog. Other cases certify all23/24 roots to3.75. All24 first12
prefixes and guards are complete; no failed numerical intervals, missing
guard, duplicated root, ceiling expansion or unresolved case remains.

### Numerical quality and provenance

| Scaled diagnostic, all accepted roots | Maximum |
| --- | --- |
| clamp_scaled_residual | 1.48697e-12 |
| equation_scaled_residual | 5.57257e-13 |
| singular_ratio | 1.09896e-12 |
| nonzero_singular_condition | 1060.29 |
| energy_relative_error | 6.53833e-12 |
| d_X | 2.19706e-12 |
| d_Y | 2.26184e-13 |
| theta | 2.17594e-12 |
| F_X | 3.00153e-11 |
| F_Y | 4.5933e-12 |
| M_node | 1.31164e-13 |
| mass_gram_max_error | 3.63059e-11 |
| c | 3.8406e-14 |
| R_node | 4.73534e-16 |

Comparator conditioning is somewhat higher than in the earlier MH pilot
(1060 vs about461), but passes the unchanged1e8 bound; it is not hidden.
Energy quotient is only an internal total-eigenproblem verification.
No energy fractions or classification fields are generated.

202 targeted/regression tests pass:32 new screening tests,166 unchanged
MH general/beta0/single/source, Bishop and prior kinematic checks, plus
4 rectangular Timoshenko regressions (7 unrelated tests deselected).
The new tests verify direct/split/swap/reflection, exact RL force, identical
Tim basis/geometry, pole coincidence retention, null handling, sign-invariant
overlaps, D_c triangle bound, provenance/cache identity and Gauss200/300
convergence. Full expensive research suite was not run.

Two initial implementation smoke issues were corrected without modifying
physics/tolerances: a dispersion dictionary key and a gate comparison that
included direct catalog poles above the gate ceiling. Earlier attempt
directories are retained; the latter has failure.json. Recomputed bundles
after adding mass-Gram/zero-scale diagnostics have distinct fingerprints.
No failed scientific interval has been removed from an accepted inventory.

### Reproduction from repository root

```powershell
python scripts/analysis/screen_coupled_longitudinal_theory_hierarchy.py --check-sources
python scripts/analysis/screen_coupled_longitudinal_theory_hierarchy.py --compute
python scripts/analysis/screen_coupled_longitudinal_theory_hierarchy.py --plot-only results/coupled_longitudinal_theory_hierarchy_screening/9a35ff23c3f43c66
python -m pytest tests/test_coupled_longitudinal_theory_hierarchy_screening.py -q
```

Working interpreter: `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`
(Python3.12.4, NumPy2.1.3, SciPy1.15.2, Matplotlib3.9.2). No packages installed.
Identity includes exact config/policy, source/PDF hashes, numerical versions,
code and frozen-reference hashes. Manifest retains Git HEAD/branch/status
and command; changed inputs/code never reuse the old fingerprint.
The immutable MH general bundle a36e72f715cb526b and beta0 bundle
3059d70b1b50ea2e are not overwritten. New profiles contain all state fields;
A/B have null c/R columns, not added contraction coordinates.

Three figures only: max_differences.png, position_differences.png and
contraction_diagnostic.png. Lines between declared points are visual guides,
not additional computed angles. Full fixed-beta overlap matrices, best/second
margins, theta norms and adjacent gaps remain tabular; no energy plots.
Plot-only regenerates from saved tables with zero roots and checks artifact
bytes. Matching compute rerun likewise checks/reuses data with zero roots.


## Scope and follow-up

No energy-based modal classification, across-beta tracking, spectral
reordering, applicability threshold, safe-prefix certification, beta
refinement/extremum search, mu/tau/thickness scan, fitting, FEM or nonlinear
derivation was performed. Historical MH gates/source formulations/Bishop
are retained unchanged. The four-field contraction closure remains qualified.
The fixed input/grid does not establish conclusions for other dimensions,
length ratios, materials, higher spectral prefixes or joint closures.

All correspondences in this screening are diagonal and model differences
are small. No position-exchange candidate was found. RL does not improve
the maximum MH difference here. Do not expand the study automatically;
retain this result as evidence for this G20 control only.
