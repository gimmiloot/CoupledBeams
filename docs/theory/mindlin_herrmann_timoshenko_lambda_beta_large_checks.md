# Production M-H/Timoshenko: large geometry implementation checks

2026-10-07. **MHTIM_LAMBDA_BETA_LARGE_CHECKS=PASS** for the declared
geometry/grid contract. This stage checks the existing production implementation,
not a new physical theory, scientific sensitivity study or applicability map.
Previous source, single-rod, joint, hierarchy and thickness evidence is retained.

## Canonical Lambda and the fixed reference

The governing project convention is the beam parameter in
`docs/theory/equations.tex` (frequency normalization, lines233–234),
`docs/project_rules.md` (Canonical Equal-Radius In-Plane Model Conventions)
and `docs/theory/main_note.md` (current article synchronization):

\[
 \Lambda^4 J/(S l^2)=\omega^2 l^2\rho/E.
\]

Here old geometric J/S becomes Iy/A for this rectangle. Thus

\[
 \Lambda^4=\frac{\rho A_{ref}\omega^2l_{ref}^4}{E I_{ref}},
 \qquad \Lambda=l_{ref}(\rho A_{ref}/(E I_{ref}))^{1/4}\sqrt{\omega}.
\]

There is a local notation difference in the earlier RLB-2B rectangular
note/scripts: their common `Lambda=omega*l^2*sqrt(rho*A/(E*Iy))` is the
**square of this canonical Lambda**. The subsequent RLB-2C
`rectangular_weakly_orthotropic_models_vs_beta_note.md` §2 explicitly restores
the mapping Omega=Lambda^2. This stage uses the verified baseline convention;
the local RLB-2B convention and all old formulas/data are left unchanged.
It neither silently relabels f_star nor silently corrects another branch.

One fixed reference is used for every curve:

\[
 l_{ref}=0.5,\quad b_{ref}=0.20,\quad h_{ref}=0.05,
 \quad A_{ref}=0.01,\quad I_{ref}=2.083333333333334\,10^{-6},
 \quad E=\rho=1.
\]

Consequently Lambda^4=300 omega^2, Lambda^2=sqrt(300) omega and

\[
 \Lambda=300^{1/4}\sqrt{2\pi f}
 =\sqrt{2\pi\sqrt{300}\,f_*},\qquad f_*=fL/\sqrt{E/\rho},\quad L=1.
\]

The exact reference factor300 is checked with Fraction arithmetic; float
conversion is checked against the fourth-power formula. The source of truth
is clear after distinguishing the scoped RLB-2B name from the baseline.
No per-arm/per-curve I or length enters this output conversion.

## Fixed production contract and geometry families

All arms use unchanged `project_jang_reduced_rectangular`, E=rho=1,
nu=.3,kappa5/6, C=EA/(1-nu^2),H=kappa GI,j=rho I. State order remains
(u,c,w,theta,N,R,Q,M). Both local x point external clamp→joint, project
t/n and beta convention are imported from the validated joint helper.
External u=c=w=theta=0 clamps and the eight common-DOF/dual joint conditions
are unchanged. Source variants, Bishop and old production arm APIs are untouched.

Length cases, total L=1 and identical b=.20,h=.05:

| mu | L1 | L2 |
| --- | --- | --- |
| 0 | .500 | .500 |
| .25 | .375 | .625 |
| .50 | .250 | .750 |

Thickness cases, L1=L2=.5,b1=b2=.20:

| thickness_contrast delta_h | h1 | h2 |
| --- | --- | --- |
| 0 | .050 | .050 |
| .20 | .040 | .060 |
| .40 | .030 | .070 |

The searched existing canonical notes/configs contained no conflicting
geometric delta_h definition. `thickness_contrast` is the machine name;
delta_h is localized to this note/figure. It is not another eta/tau parameter.
Both families preserve mass0.01: length cases preserve section and total
length; thickness cases preserve h1+h2=2h0 at equal lengths. Individual
A,I,EI,rho I,H vary physically; no mass/material rescaling is applied.

Main beta grid is exactly0:2.5:90 degrees,37 points. Spectra are independently
sorted at each point. Figures show Lambda1–12 only; root13 is a saved
completeness guard. Equal baseline appears in both figures but is solved once.

## Production composition and beta0 references

The old joint convenience API has one arm model shared by two arms. The
new diagnostic adapter composes two calls to the existing `joint.arm_basis`
with the **same** `joint.joint_matrix(joint.frames(beta))`. It copies no local
wave/constitutive equations or sign convention. The two-model composition
is bitwise equal to the old homogeneous full boundary matrix when models
are equal, including asymmetric lengths. The per-arm end dynamic stiffness
also comes from the unchanged verified helper.

At beta0, changing mu is an artificial partition of one uniform CC rod.
All three first12+guard13 inventories match the preserved direct-rod
frequencies. Full profiles/resultants additionally match a separately
assembled straight-coordinate representation of that same rod.

For unequal sections, the independent stepped reference uses two positive-X
segment states, first clamp at segment1 x=0, second clamp at segment2 x=L2,
and eight **literal state-continuity** equations at the interface. It never
calls a general-angle transformation. Its determinant is solved independently
inside the primary certified brackets. Its q/p profiles are compared using
the already verified reversal map for arm2, not signs chosen from spectra.
This is independent assembly of the same 1D theory, not experimental validation.

## Arm swap and signed rotation

Controls use parameter pairs ±.25/±.50 and ±.20/±.40 at beta0/22.5/45/67.5/90.
With canonical frames, swapping the parameter is physical reflection in the
angle bisector followed by arm relabeling:

\[
 S_\beta=\begin{bmatrix}-\cos\beta&-\sin\beta\\-\sin\beta&\cos\beta\end{bmatrix},
 \quad S_\beta t_1=t_2,\quad -S_\beta n_1=n_2.
\]

The same existing mirror state map applies:
(u,c,w,theta,N,R,Q,M)→(u,c,−w,−theta,N,R,−Q,−M).
Contraction c/R is scalar, theta/M is signed under reflection. Proper frame
rotations leave the two scalar nodal coordinates unchanged; force maps remain
dual to displacement maps. This remapping is a same-system symmetry test,
not across-beta mode tracking. Frequencies, kinematic fields and resultants
all pass; eigenvector global sign is immaterial.

## Root localization, completeness and cost

This is a fast_plot instance of frequency-map-v1 with the explicitly requested
normalization/reference/symmetry gates. There are5 unique main geometries,
185 main geometry×beta cases and20 sparse negative-parameter controls.
The figure/long-CSV contract has6×37=222 logical cases because baseline is
shared by both families. No additional beta points or geometry refinement occur.

Each geometry begins with a certified400-interval seed search. Next points
use previous **sorted frequencies only** to form midpoint brackets. Every
point recomputes the independent energy/Schur count, uses bounded count-guided
subdivision as needed, independently sorts the resulting roots and checks
the absolute count indices1–13. No forms, MAC, assignments or descendant
identities enter numerical acceleration. Failed local localization would
trigger one full-scan fallback for that point; none occurred.

The independent count is the sum of per-arm certified fixed-fixed poles plus
negative nodal Schur inertia. Pole catalogs use existing finite/min-max
functions. A fixed **reference** positive congruence balances the Schur matrix;
raw symmetry and conditioning are checked before roundoff symmetrization.
Queries too near a pole shift away and record requested/effective frequencies.
The physical determinant roots are not shifted. This handles possible
root/pole coincidences without forcing equal counts across the pole.

Initial configured f_star ceiling3.75 is retained, with one allowed1.5
expansion on missing guard. All prefixes are covered below the unchanged
analytic basis optical cutoffs; no expansion was used. No unsupported tail,
unresolved multiplicity, bad point interpolation or curve smoothing is accepted.

At every root: clamp, eight joint rows, force/moment, c/R, PDE, mass norm,
energy quotient, SVD and per-case mass Gram checks. Old thresholds remain:
boundary/joint/PDE1e-9, energy5e-8, Gram5e-7, nonzero condition1e8.
These are numerical quality gates, not applicability or spectral-difference
thresholds. The immutable production physics helpers and old bundles retain
their byte hashes. No new physics module or another general-beta theory exists.

## Curvature flags

Discrete second differences of each sorted Lambda sequence are inspected
with a median/MAD outlier diagnostic (multiplier12 plus machine-scale floor).
This is not a physical smoothness criterion or an extrema/veering search.
All68 flagged samples keep their certified brackets, counts and residuals.
Each is independently checked **at the same beta** by harmonic-state
short-step expm/positive-diagonal QR, retaining both displacements and forces
at each arm end. No ill-conditioned long transfer-matrix product is used.

All independent roots agree within the existing1e-10 frequency tolerance and
their projected boundary SVD gate. Flags are reviewed/retained, not erased,
smoothed or interpolated. This supports numerical correctness of the sampled
features without assigning a physical cause to them or claiming proof of
smooth derivatives/branch identity between samples.

## Results and implementation assessment

Canonical local bundle: `results/mindlin_herrmann_timoshenko_lambda_beta_large_checks/356c3be4953268f2/`.

| Status | Result |
| --- | --- |
| MHTIM_LAMBDA_NORMALIZATION | PASS |
| MHTIM_LAMBDA_MAP_BASELINE_REGRESSION | PASS |
| MHTIM_LENGTH_BETA0_COLLAPSE | PASS |
| MHTIM_LENGTH_ARM_SWAP | PASS |
| MHTIM_LENGTH_LAMBDA_BETA_MAP | PASS |
| MHTIM_THICKNESS_STEPPED_BETA0 | PASS |
| MHTIM_THICKNESS_ARM_SWAP | PASS |
| MHTIM_THICKNESS_LAMBDA_BETA_MAP | PASS |
| MHTIM_LAMBDA_BETA_ROOT_QUALITY | PASS |
| MHTIM_LAMBDA_BETA_LARGE_CHECKS | PASS |

| Family / parameter | 37-point prefixes | Max joint residual | Max nonzero condition | beta0 matching frequency difference | Swap frequency difference | Curvature flags |
| --- | --- | --- | --- | --- | --- | --- |
| length_asymmetry / 0 | PASS | 1.942e-11 | 10639 | 3.331e-15 | 0 | 25 |
| length_asymmetry / 0.25 | PASS | 1.135e-11 | 103.05 | 2.665e-15 | 6.147e-12 | 12 |
| length_asymmetry / 0.5 | PASS | 1.586e-11 | 92.512 | 1.554e-15 | 4.828e-12 | 8 |
| thickness_contrast / 0 | PASS | 1.942e-11 | 10639 | 3.331e-15 | 0 | 25 |
| thickness_contrast / 0.2 | PASS | 1.835e-11 | 292.66 | 2.22e-16 | 1.511e-12 | 12 |
| thickness_contrast / 0.4 | PASS | 1.421e-11 | 106.06 | 1.11e-15 | 6.054e-12 | 11 |

Length beta0 collapse against the preserved direct CC rod: maximum relative
frequency difference6.00498e-12 across all3 partitions and13 roots.
Positive-X uniform matching/profile comparison: max frequency3.33067e-15,
full kinematic/resultant L2 error8.20346e-14. Stepped unequal-section
reference: max frequency1.11023e-15, full profile/resultant L2 error4.04448e-13.

Sparse swaps: length frequency6.14731e-12, full-field L2 error4.48409e-11;
thickness frequency6.05450e-12, full-field L2 error3.35882e-11. All c/R,
signed theta/M and physical force rows pass; no signs were selected from
frequency agreement. Single eigenvector sign alignment is conventional.

Preserved baseline8-angle regression: maximum frequency3.98237e-12 and
Lambda1.99130e-12, consistent with the square-root conversion. Baseline
common points are freshly solved; previous results are not copied as new roots.

| Quality metric,185 main +20 control cases | Maximum |
| --- | --- |
| clamp_scaled_residual | 6.78768e-12 |
| equation_scaled_residual | 1.47958e-11 |
| singular_ratio | 2.26004e-12 |
| energy_relative_error | 9.55747e-12 |
| d_X | 1.08798e-12 |
| d_Y | 1.25199e-12 |
| c | 2.2961e-14 |
| theta | 4.55061e-12 |
| F_X | 1.94173e-11 |
| F_Y | 1.83482e-11 |
| R_node | 1.24548e-14 |
| M_node | 6.69163e-13 |
| mass_gram_max_error | 5.16933e-11 |

Nonzero condition range 2.84793--10639; no tolerance change.
All205 unique computed cases certify positions1?13. No duplicated/lost
root, NaN/Inf, unresolved point, failed quality gate, fallback or ceiling
expansion. Main data contain185 unique points and222 logical family points;
Lambda_spectra.csv has2886 rows including guards, not a tracked-branch table.

68 unique curvature flags: baseline25,length .25=12,length .50=8,
section .20=12,section .40=11. Baseline flags are shared by two figures.
Independent same-point QR maximum frequency difference7.20602e-12.
All flags remain in result.json with thresholds, certified brackets, quality
and independent diagnostics. They are reviewed numerical candidates; no
unresolved numerical spike was found and no physical extremum/crossing
claim is made. Curves are not smoothed.

| Performance counter | Final run |
| --- | --- |
| determinant_evaluations | 32196 |
| count_evaluations | 4879 |
| catalog_evaluations | 7862 |
| reference_determinant_evaluations | 583 |
| full_scans | 5 |
| local_continuations | 200 |
| fallback_full_scans | 0 |
| subdivisions | 74 |
| anomaly_independent_evaluations | 975 |
| ceiling_expansions | 0 |
| compute_runtime_seconds | 137.267874 |

The200 predictors comprise180 successive main-grid points and20
same-system sparse controls seeded by their positive counterpart. This
is numerical localization, not shape or branch continuation. Full400-node
seed scans occur only at beta0 for the5 unique geometries. Both prototype
and final audit runs are retained; the final adds mass Gram and independent
flag checks. Its137.27s includes all references/controls/QR reviews, not
just plotting. Plot-only and matching-cache checks perform zero root calls.

273 relevant tests PASS:36 new normalization/paired-assembly/geometry/
collapse/stepped/swap/quality/independent-review/cache/plot-only checks,
233 unchanged hierarchy/thickness/MH general/beta0/single/source/Bishop
checks, plus4 rectangular Tim regressions (7 unrelated tests deselected).
New tests verify paired-arm matrix equality, reflection signs, total mass,
canonical conversion and source/old-bundle identity; no full expensive
research suite was run. git diff --check and new-file syntax/whitespace
checks pass. Physics/source/previous-test files and all old bundles are
preserved by initial byte-hash comparison.

The production general-beta implementation is sufficiently checked for
these declared geometry/grid ranges under the adopted model. This is a
sampled implementation assessment, not a proof of physical3D accuracy or
a universal all-geometry/all-frequency validation.


## Figures, commands and provenance

Two figures only, each with3 panels, shared axes, fixed reference Lambda
scale, matched typography and the same color/dash identity for sorted k:

- `length_asymmetry_Lambda_beta.pdf` and `.png`;
- `thickness_contrast_Lambda_beta.pdf` and `.png`.

Caption for both: **Production M-H/Timoshenko; independently sorted spectral
positions k=1–12 at37 fixed angles, not tracked modal branches.** Root13 is
omitted from the plots. No energy labels, theory overlays or third science
figure are generated. PDF is vector with embedded fonts; PNG is300dpi.

```powershell
python scripts/analysis/verify_mindlin_herrmann_timoshenko_lambda_beta_large_checks.py --check-sources
python scripts/analysis/verify_mindlin_herrmann_timoshenko_lambda_beta_large_checks.py --compute
python scripts/analysis/verify_mindlin_herrmann_timoshenko_lambda_beta_large_checks.py --plot-only results/mindlin_herrmann_timoshenko_lambda_beta_large_checks/356c3be4953268f2
python -m pytest tests/test_mindlin_herrmann_timoshenko_lambda_beta_large_checks.py -q
```

Working interpreter: D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe
(Python3.12.4, NumPy2.1.3, SciPy1.15.2, Matplotlib3.9.2). No packages installed.
PDF/PNG plots were visually checked for labels/legend/axes; export is deterministic.


New CLI is justified by the paired-section validation/reference/accelerated-map
contract. It reuses production arm, dynamic-stiffness, geometry and source
helpers; old fixed-geometry CLIs and APIs are not broadened or invalidated.
Manifest records geometry/beta grids, normalization and its source hashes,
production/solver/code/version provenance, Git, command and artifact hashes.
Full q/p coefficient profiles, references, swap controls, counts/brackets,
quality/Gram/anomaly reviews, Lambda CSV and performance are retained.
Repeat compute checks/reuses the bundle with zero roots. Plot-only reads
saved data and also performs zero roots; deterministic PDF/PNG bytes are checked.

## Remaining modeling caveat

These checks support numerical stability, sampled continuity, normalization
and symmetry of the implementation across the declared large geometric
changes. They do not validate 1D physics against3D elasticity, establish
scientific sensitivities/veering/applicability or authorize another study.
**c1=c2 remains variational reduced-frame closure, not a direct 3D derivation
for the finite welded joint region.** Coefficients/kappa/clamps were not fitted
or changed. No across-beta tracking, energy classification, beta refinement,
damping, FEM, nonlinear or out-of-plane work was performed. Stop here.
