# FEM-1: full low-frequency spectrum of a thick rectangular rod

Calculations: 2026-10-08; documentation/verification finalized: 2026-10-09.
Main HEAD `7d1363cc47d9cea426515fa0783a0643b47a5361`.
**Bounded execution and all-four-family identification complete; mesh convergence
and quantitative all-family validation PARTIAL.** Three real solid meshes/jobs
were accepted; five of eight medium-to-fine changes exceed the predeclared0.1%
numerical convergence criterion. No fourth grid or nonlinear calculation followed.

## Scope and pre-FEM geometry decision

This is a linear check of the accepted seven-field model
q=(u,w,v,Phi,psi,theta,c), not validation of every nonlinear coefficient in V0.
The [spatial model](../theory/weakly_nonlinear_spatial_rod.md),
[Jang-type single rod](../theory/mindlin_herrmann_timoshenko_single_rod.md),
[rigid-joint note](../theory/mindlin_herrmann_timoshenko_rigid_joint.md),
[Yartsev reproduction](../anisotropic_rods/yartsev_ch2_single_rod_reproduction.md),
[old anchor design](../anisotropic_rods/yartsev_ch2_limited_3d_fem_anchor_design.md),
and [readiness audit](nlsp_3d_fem_environment_readiness.md) set the contracts.
Frozen physical files, old FEM scripts/results and pre-existing readiness edits
were preserved. Historical G20 h=.05 was not changed.

Chosen diagnostic: L=1,b=.20,h=.10,E=rho=1,nu=.3,kappa=5/6.
h/L=.10,b/L=.20,b/h=2. Width/mass were not rescaled. h=.10 passed the independent
1D selection before ANY new FEM result: first acoustic axial mode is position8
(<15). The only permitted backup h=.12 was NOT evaluated. Input config contains
the completed pre-FEM screen hash and was frozen before mesh generation.

## Frozen 1D comparator and torsion

The four independent blocks are MH(u,c), Tim(w,theta), Tim(v,psi), scalar twistPhi.
Both ends impose all seven field VALUES zero; no book_slope_clamp is used.
Existing analytic finite-rod root/basis solvers supply MH/Tim roots and profiles;
independent expm/QR boundary residuals, PDE and Rayleigh-energy checks pass.
Each block saturates an existing min-max count upper bound in the final window.
Torsion uses its exact scalar Dirichlet spectrum. This 1D computation is a
comparator/completeness check, not independent3D validation.

A0=.02,I_parallel=1.6666666667e-5,I_perp=6.6666666667e-5,
Ip=8.3333333333e-5. The h direction uses I_parallel; b direction uses I_perp.
m=.02,j_parallel=I_parallel,j_perp=I_perp,C=.0219780219780,
H=5.34188034188e-6,S=.00641025641026,B_parallel=I_parallel,B_perp=I_perp.

C_T=1.759089824002232e-5 comes from the existing Yartsev generalized section
reduction, book Geometry(a=projectb,b=projecth), isotropic input at theta=0.
Sbar16=0,C_T=Cbar;367 series terms give estimated relative tail<1e-12.
**G*Ip was not substituted**, no stiffness was fitted. Condensed sectional warping
in C_T does not introduce an independent dynamic warping coordinate.

Required maximum is first MH acoustic omega3.150123638896353. A fixed10% guard
sets omega_max=3.4651360027859885. Eight modes lie inside:3 in-plane,2 out-of-plane,
2 torsion,1 MH. Initial24 eigenpairs were specified before FEM; the single allowed
extension to36 was unnecessary. Contraction optical cutoffomega36.31365196 and
its associated branch are OUTSIDE_FEM_1_FREQUENCY_WINDOW. No hundreds-mode search. The hashed original pre-selection scratch JSON retains
draft keys ending_hz; these values mean cycles per normalized time. Active
preflight/CSV keys and the comparison use explicit normalized units, not Hz.

## Geometry, axes, solid BC and software

One monolithic Box, global X in[0,1],Y in[-.05,.05],Z in[-.10,.10].
B=[t,n,k]=diag(1,-1,-1); positive w is -Y, positive v is -Z.
Local rotation=(Phi,-psi,theta). The h direction is n, b direction k;
no interchange of inertias or new coordinate convention occurs.
Every node on BOTH complete X faces has Ux=Uy=Uz=0; lateral faces are free.
No internal interface, two-arm split, rigid body or angular joint is introduced.

Gmsh4.15.2: `D:\PHD\gmsh-4.15.2-Windows64\gmsh.exe`.
CCX2.22 single-threaded build: old cache
`results/_smoke/3d_fem_environment_check/calculix_2p22/CalculiX-2.22.0-win-x64/bin/ccx.exe`.
Current binaries match the readiness SHA256 evidence; sidecar DLLs remain in place.
No relocation, install or PATH change. Absolute paths recorded in input/jobs.
Python3.12.4,NumPy2.1.3,SciPy1.15.2,Matplotlib3.9.2. Processes use one CPU thread,
900s individual timeouts,4GiB memory ceiling and3600s total numerical budget.

Each input contains only linear ELASTIC/DENSITY/FREQUENCY/NODE FILE U.
No NLGEOM/STATIC/DYNAMIC. Units are consistent normalized G20 units:
omega=sqrt(eigenvalue)=2pi*f; no experimental-Hz or circular Lambda conversion.
The three independent DAT columns agree within8.97e-7 relative printed precision,
below the predeclared1e-5 format-consistency gate.

## Mesh quality and resources

C3D10 quadratic tetrahedra, target sizes h/2,h/3,h/4 fixed BEFORE FEM.
The reused helper generates two exports per level (MSH4.1 and INP); these are
three resolution levels, six Gmsh CLI invocations, not six refinement levels.
Geometry/mass, face node sets, face connectivity/conforming quadratic nodes,
unused nodes, straight midsides and sampled Jacobians are checked BEFORE CCX.
Straight affine midsides make J constant; corner/centroid/14-quadrature samples
confirm positive J. Quality is a normalized corner determinant, not a universal
FEM-error tolerance. No Gmsh or CCX warnings were recorded.

| Mesh | Target size | Nodes | C3D10 | Left/right fixed nodes | Min J | Min corner quality | Total stage s | Peak MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| coarse | 0.05000000 | 2092 | 1054 | 57/57 | 3.130048e-05 | 0.153313 | 2.819 | 29.30 |
| medium | 0.03333333 | 5649 | 3120 | 119/119 | 1.241477e-05 | 0.152998 | 7.190 | 76.93 |
| fine | 0.02500000 | 11553 | 6670 | 193/193 | 4.841763e-06 | 0.125578 | 15.857 | 165.40 |

All meshes have volume/mass.02 to arithmetic accuracy, one face-connected solid,
zero nonconforming quadratic faces, zero negative/zero Jacobian elements,
nondegenerate fully fixed face sets and no unconstrained rigid-body motion.
Actual corner-edge segment counts are h:2/3/4 and b:4/6/8. Midside nodes do not
count as additional elements. Nine interior rays count intersected tetrahedra
h:6-8/8-11/10-13 and b:12-15/19-24/25-32; these are cell intersections, NOT claims
of structured layers. Bbox and actual ray coverage are stored for each level.

Before fine generation, medium costs yielded a predeclared forecast31.18s and
280.8MiB; actual fine stage15.86s/165.4MiB. All three CCX modal jobs returned0
with24 positive eigenpairs and complete nodal vectors each; maximum clamp
relative residual0. Total staged numerical work26.305s/3600s (CCX calls1.21,
3.42,7.85s); no eigenpair extension or fourth mesh.

## Shape identification and FRD reader qualification

Comparison uses full mass-weighted displacement MAC at positive14-point C3D10
quadrature, plus rigid-section fits over41 disjoint axial slabs. The degree5
rule exactly integrates quadratic FE mass norms on affine tetrahedra; the
lifted1D fields are sampled from401-point analytic profiles. Their overlap is
numerical shape evidence, not an exact continuous identity.

The linear mass-kinematic lift is
U_local=(u-theta*eta-psi*zeta,w-Phi*zeta+c*eta,v+Phi*eta),
with eta=-Y,zeta=-Z. Fits recover D and Omega=(Phi,-psi,theta), allowing linear
axial variation inside each slab so normal axial variation is not labeled warping.
Family diagnostics are displacement mass-norm components, not stiffness-energy
classification. Residuals can retain local/cross-sectional or ambiguous character.

Shape-only one-to-one assignment retains independent-best duplicates/margins;
frequencies are NOT assignment inputs. Diagnostic criteria minMAC.70,margin.08,
family dominance.70,section residual fraction.35 were frozen before FEM, without
any physical model-discrepancy acceptance threshold. All eight matches are unique,
all family characters agree, minimum full-vector MAC.990059,margin.987935.
No ambiguous match or close-subspace identity was accepted silently.

A confirmed OLD FRD reader limitation merged fixed-width node IDs with positive
Ux. A scoped strict adapter fixes this only in the NEW workflow; old parser and
results are unchanged. It validates all IDs/components/blocks and E12.5E3 overflow.
Actual evidence:24 complete vectors of2092/5649/11553 nodes respectively, with
no missing vectors or imputed nodal values. Tests cover the historical format.

## Full bounded comparison

All angular frequencies below are rad per normalized time. Absolute difference
uses abs(omega1D-omegaFine)/omegaFine; signed values are negative for all eight.
Mesh change is separately abs(omegaFine-omegaMedium)/omegaFine.
Numerical mesh criterion was0.1% plus decreasing successive changes, selected
before seeing FEM. It is not a physical applicability threshold.

| Sorted1D | Family/local | 1D omega | Coarse omega | Medium omega | Fine omega | Model difference % | Fine MAC | Medium-fine % | Mesh status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | inplane_bending/1 | 0.605416730 | 0.6182542 | 0.6166689 | 0.6159927 | 1.71690 | 0.9998746 | 0.10977 | MESH_UNRESOLVED |
| 2 | outplane_bending/1 | 1.038923589 | 1.0552370 | 1.0543180 | 1.0536820 | 1.40065 | 0.9998622 | 0.06036 | PASS |
| 3 | torsion/1 | 1.443392698 | 1.5231980 | 1.5132500 | 1.5108350 | 4.46391 | 0.9970668 | 0.15985 | MESH_UNRESOLVED |
| 4 | inplane_bending/2 | 1.551535483 | 1.5863410 | 1.5810000 | 1.5788980 | 1.73301 | 0.9993243 | 0.13313 | MESH_UNRESOLVED |
| 5 | outplane_bending/2 | 2.378101731 | 2.4199780 | 2.4178940 | 2.4165000 | 1.58900 | 0.9993740 | 0.05769 | PASS |
| 6 | inplane_bending/3 | 2.804275820 | 2.8728140 | 2.8605930 | 2.8559300 | 1.80866 | 0.9981096 | 0.16327 | MESH_UNRESOLVED |
| 7 | torsion/2 | 2.886785396 | 3.0670850 | 3.0451820 | 3.0400340 | 5.04102 | 0.9901749 | 0.16934 | MESH_UNRESOLVED |
| 8 | axial_mh/1 | 3.150123639 | 3.1679230 | 3.1666160 | 3.1657330 | 0.49307 | 0.9972070 | 0.02789 | PASS |

All changes decrease and all matched3D frequencies decrease with refinement.
Only first axial and two out-of-plane bending modes meet the preset criterion;
the five remaining modes remain MESH_UNRESOLVED. Frequency differences therefore
are measured finite-grid discrepancies, not fully established continuum model errors.
Cross-mesh section-profile MAC is at least.9999877 coarse/medium and.9999977
medium/fine; identities are stable independently of raw sorted positions.

### Bending planes

In-plane differences are1.717-1.809%, medium/fine0.110-0.163% (unresolved criterion).
Out-of-plane differences1.401/1.589%, medium/fine0.0604/0.0577% (accepted numerical
comparison). Distinct inertias and section rotations were used. No source slope
clamp or circular bending-doublet averaging was imported.

### Axial MH and effective contraction

First MH acoustic3.15012364 matches fine3D3.165733 (raw FEM position8), difference
0.493073%, medium/fine0.0278924%. Its full-vector MAC.997207 and axial section
u-overlap.9999626 support the main displacement structure. Classical pi*sqrt(E/rho)/L
may explain a leading limit but was NOT substituted as the main comparator.

c_eff=<eta*Ures_n>/<eta²> is an effective thickness-strain diagnostic after rigid
section removal, not a Cartesian DOF or exact identity with MHc. Width strain is
recorded separately. Aligned unit-mass profiles give c shape MAC.998731,
maximum c difference.784150 versus reference characteristic9.16027 (about8.56%).
These are arbitrary LINEAR EIGENVECTOR scales, not finite-motion strains or a
safety violation. Slab sampling/fit and end restraint qualify this measurement;
no boundary zone is removed from the profile comparison. The 1D and solid-face
clamps constrain local deformation differently. No H/j/kappa/nu/C fitting.

### Torsion, warping and end restraint

First and second generalized-twist differences are4.46391/5.04102%, the largest
observed. Medium/fine changes0.15985/0.16934% remain numerically unresolved.
Strong full-vector MAC.997067/.990175 confirms a twist-like counterpart rather
than accidental frequency matching. Fine axial warping mass fractions are.001928
and.007224, total section residual fractions.001958 and.007456. Small squared
mass-norm residuals do not bound stiffness/frequency error.

The source C_T condenses static sectional warping; no independent dynamic warping
field exists. Solid-face clamps restrict axial warping, unlike only Phi=0 in1D.
These are plausible physical sources of persistent differences, NOT separately
proven causal attribution. Fine-grid changes are much smaller than4-5% observed
differences but still exceed the chosen accuracy criterion. CT was not corrected.

There are **zero additional3D modes** inside the selected window, zero unmatched
1D modes and zero ambiguous/duplicate assignments. Raw24-mode spectra/vectors
are retained; extra modes beyond the window are not claimed missing1D roots.
Contraction-dominated branch remains OUTSIDE_FEM_1_FREQUENCY_WINDOW; no extra search.

## Reproduction, artifacts, tests and stop

New [CLI](../../scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py),
[immutable config](../../data/input/nlsp_linear_rectangular_3d_fem.json),
[targeted tests](../../tests/test_nlsp_linear_rectangular_3d_fem.py).
One new workflow has a genuinely new full-seven-field rectangular shape-comparison
contract; it reuses existing MH/Tim/Yartsev physics and modal helpers, not a second
physics solver. Results: `results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/`.

```powershell
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py --preflight
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py --run-fem
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py --report-only results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py --plot-only results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad
```

Matching preflight/run-fem cache and report/plot perform zero new Gmsh/CCX,
1D eigenanalysis or BVP. Cache identity covers code/config, models, binary hashes
and dependency versions. Old outputs are never a stale-cleanup target.
The bundle retains frozen geometry/config,1D roots/count certificates/profiles,
GEO/MSH/INP/DAT/FRD, logs/warnings/commands, mesh/Jacobian/resolution data,
all nodal eigenvectors, section/MAC evidence, comparison/raw-spectrum CSV,
convergence and axial diagnostics, resource forecast/counters and hashes.

Three PDF+PNG figures: `combined_linear_spectrum`,
`frequency_difference_and_mesh_convergence`, `representative_all_family_shapes`.
Mode shapes are displayed at common unit displacement-mass normalization with
one overall sign; amplitudes are not physical nonlinear initial displacements.
No curves connect missing identities; no unresolved frequency point is concealed.

49 targeted tests PASS, including synthetic geometry/quadrature/Jacobian/FRD,
frozen coefficients/BC/CT/completeness, matching/cache/units and actual saved-data
checks. Tests start no3D jobs. Links, protected-file/index preservation and
`git diff --check` are checked separately; no historical test modified.

| Status | Outcome |
|---|---|
| NLSP_FEM1_GEOMETRY_SCREEN | PASS |
| NLSP_FEM1_LINEAR_1D_SPECTRUM | PASS |
| NLSP_FEM1_3D_MESH_QUALITY | PASS |
| NLSP_FEM1_3D_MODAL_EXECUTION | PASS |
| NLSP_FEM1_3D_MESH_CONVERGENCE | PARTIAL |
| NLSP_FEM1_MODE_IDENTIFICATION | PASS |
| NLSP_FEM1_ALL_FAMILY_COMPARISON | PARTIAL |

This bounded linear result supports identification of all four families and
quantifies differences; it does not validate all nonlinear terms in V0, prove
universal1D applicability or treat fine3D as exact truth. Mesh convergence PARTIAL
is not MODEL_INVALID. No fourth mesh or threshold fitting, h=.12, parameter map,
nonlinear static/dynamic, angular joint or out-of-plane stability calculation ran.

[NLSP-D09](../memory/decisions.md#nlsp-d09) and
[NLSP-K10](../memory/knowledge.md#nlsp-k10) record scope/result. LONG CLOSED,
EB/RLB-KV PAUSED, angular same-clamp UNAVAILABLE, prepared planar strict PARTIAL
and physical sanity DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS remain unchanged.
FEM-2/FEM-3 and further refinement require a separate explicit decision.


<a id="fem-1r-one-additional-mesh-refinement"></a>
## FEM-1R: one additional mesh refinement

2026-10-09; main HEAD `efb7966fbaf44ef956666477a89904c00731e033`.
This separately authorized continuation adds exactly one target-size
.020 C3D10 level to the historical coarse/medium/fine results above. The original
FEM-1 report, three-grid PARTIAL statuses, decisions and bundle remain history;
they are not rewritten as if the fourth grid had already existed.

The original geometry/material, axes, solid-face clamps, mesh generator/algorithm,
quadratic element order, mass, output, 24-mode request and matching policy remain
unchanged. L=1,b=.20,h=.10,E=rho=1,nu=.3,kappa=5/6. The same normalized time is used:
omega is rad per normalized time, with cyclic frequency=omega/(2*pi); no Hz or
Lambda conversion is introduced. V0, cubic equations and all 1D coefficients,
including C_T and H, are untouched.

### Frozen parent and continuation contract

Parent: `results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/`.
Its pinned manifest SHA256 is
`9f2d5139b84b2aa133b20d9a7cae806bac085c178fba506e087cae2941da1d90`.
The historical manifest and all listed artifacts are checked before a new job;
the saved eight 1D frequencies, mode profiles, section reduction, completeness
and original three mesh results are loaded directly. No new roots, eigensystems,
BVP, torsional-series evaluation or 1D profiles are computed. Current helper and
binary hashes must match the parent identities; a changed new CLI hash does not
invalidate or regenerate the historical bundle.

The [continuation config](../../data/input/nlsp_linear_rectangular_3d_fem_refinement.json)
and [CLI](../../scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py)
use a separate schema/output contract. This is a diagnostic parent-validated
orchestration layer, not a second FEM solver. It reuses `rectangular_geo`,
`audit_rectangular_mesh`, `inspect_modal`, strict FRD reading, section projection,
MAC/assignment and the historical Gmsh/CalculiX generator/runner. Altering the
frozen three-grid config/schema would change its cache semantics; this new
immutable-parent continuation contract is the Script Proliferation Control
reason for a separate reusable entry point.

### New mesh and modal job

| Level | Target size | Nodes | C3D10 | Left/right fixed-face nodes | Min Jacobian | Corner quality min/median/max | Volume/mass |
|---|---:|---:|---:|---:|---:|---|---:|
| refined | 0.020000 | 20752 | 12687 | 279/279 | 2.559649e-06 | 0.147842/0.430992/0.707107 | 0.020000 |

The bbox remains X=[0,1],Y=[-.05,.05],Z=[-.10,.10]. One face-connected solid,
consistent quadratic faces, straight midsides, no negative/zero Jacobian element,
and full end-face constraints pass before CalculiX. Volume relative discrepancy
is5.20e-16. Corner/centroid/positive14-point volume-quadrature samples confirm the
affine positive Jacobians. No unused or artificially restrained internal nodes.

Actual corner-edge segment counts, excluding midside nodes, are5 through h and10
through b. Nine interior rays intersect15-18 tetrahedra through h and28-40 through
b; these are cell intersections, not structured-layer counts. Ray coverage is
complete to arithmetic accuracy. Gmsh and CalculiX warning lists are empty.

Same Gmsh4.15.2/CalculiX2.22 binaries and one-thread policy were used, without
installation, relocation or PATH changes. Two deterministic Gmsh format exports
(MSH4.1 and INP) describe ONE additional resolution level. Exactly ONE CalculiX
modal job returned0 and produced24 positive elastic eigenpairs and24 complete
20752-node displacement blocks. Success also requires parsed DAT/FRD evidence;
exit status alone is not sufficient. Maximum fixed-face relative displacement
is0. The three printed DAT frequency columns agree within6.33e-7 relative, below
the inherited1e-5 format-consistency gate.

Mode24 has omega10.79142, so the original omega_max3.4651360027859885 window is
covered. No additional 3D mode occurs inside this window. No extension of the
24-mode request, fifth grid, backup h=.12 or nonlinear calculation was attempted.
The stage used32.8429s/1200s; CalculiX18.1133s and353.410MiB observed peak working
set. Individual timeout900s and memory limit4GiB were retained.

### Shape identity and full four-mesh frequencies

All eight modes are independently matched to the saved 1D lift using the same
positive mass-quadrature MAC and shape-only unique assignment. Frequency proximity
is not an assignment input. Family fractions, section residuals, axial warping,
independent-best duplicates/margins and all24 full vectors are retained. All
family dominance checks pass (minimum matched family fraction>.99999967), with
zero ambiguous/duplicate assignments. The old raw IDs1-8 re-emerge from the shape
checks; they were not copied as the assumed assignment.

Fine/refined identity is checked against all24 refined candidates with the
historical common241-point translation/rotation section mass metric; c_eff is
excluded from that metric. Minimum cross-mesh MAC=.999999142, minimum margin
=.966888. Minimum full refined1D/3D MAC=.9902147, assignment margin=.9881806.
In particular, both twists, in-plane mode3 and axial mode1 retain identity.

| Family | Mode | 1D omega | Coarse .050 | Medium .033333 | Fine .025 | Refined .020 |
|---|---:|---:|---:|---:|---:|---:|
| In-plane bending | 1 | 0.605416730 | 0.6182542 | 0.6166689 | 0.6159927 | 0.6156668 |
| Out-of-plane bending | 1 | 1.038923589 | 1.0552370 | 1.0543180 | 1.0536820 | 1.0533170 |
| Torsion | 1 | 1.443392698 | 1.5231980 | 1.5132500 | 1.5108350 | 1.5098780 |
| In-plane bending | 2 | 1.551535483 | 1.5863410 | 1.5810000 | 1.5788980 | 1.5779590 |
| Out-of-plane bending | 2 | 2.378101731 | 2.4199780 | 2.4178940 | 2.4165000 | 2.4157600 |
| In-plane bending | 3 | 2.804275820 | 2.8728140 | 2.8605930 | 2.8559300 | 2.8540220 |
| Torsion | 2 | 2.886785396 | 3.0670850 | 3.0451820 | 3.0400340 | 3.0379440 |
| Axial M-H acoustic | 1 | 3.150123639 | 3.1679230 | 3.1666160 | 3.1657330 | 3.1652110 |

All values above are ACTUALLY calculated FEM frequencies or immutable parent1D
values. No extrapolated continuum frequency appears in the table.

### Mesh criterion and updated model differences

For each transition use abs(omega_new-omega_old)/omega_new. The new main control
is abs(omega_.020-omega_.025)/omega_.020, with the unchanged preset criterion
<=1e-3 (0.1%). The inherited decreasing-successive-change qualification also
holds. All eight satisfy the criterion, with fine/refined changes.01649-.06880%.
Each is recorded as **MESH_ACCEPTED_AT_PRESET_TOLERANCE**, not
EXACT_3D_CONTINUUM_REFERENCE. Historical medium/fine statuses remain unchanged.

| Family | Mode | Coarse-medium % | Medium-fine % | Fine-refined % | Fine/refined section MAC |
|---|---:|---:|---:|---:|---:|
| In-plane bending | 1 | 0.25707 | 0.10977 | 0.05293 | 0.999999940 |
| Out-of-plane bending | 1 | 0.08717 | 0.06036 | 0.03465 | 0.999999964 |
| Torsion | 1 | 0.65739 | 0.15985 | 0.06338 | 0.999999816 |
| In-plane bending | 2 | 0.33782 | 0.13313 | 0.05951 | 0.999999768 |
| Out-of-plane bending | 2 | 0.08619 | 0.05769 | 0.03063 | 0.999999886 |
| In-plane bending | 3 | 0.42722 | 0.16327 | 0.06685 | 0.999999425 |
| Torsion | 2 | 0.71927 | 0.16934 | 0.06880 | 0.999999142 |
| Axial M-H acoustic | 1 | 0.04127 | 0.02789 | 0.01649 | 0.999999854 |

All eight frequencies decrease across all four grids, and their absolute changes
decrease at each refinement. No identity change or irregular-trend flag occurs.
This is a regular observed trend, not a proven asymptotic law or rigorous error
bound. Grid-size ratios differ; no constant-ratio Richardson formula, effective
order fit or extrapolation was needed or used.

Updated model difference is abs(omega_1D-omega_.020)/omega_.020; the signed quantity
(omega_1D-omega_.020)/omega_.020 is retained separately. Every signed value remains
negative. The last column below uses ACCEPTED as shorthand for
MESH_ACCEPTED_AT_PRESET_TOLERANCE.

| Family | Mode | Fine-refined % | Signed 1D/refined % | Absolute 1D/refined % | Change from previous model difference, percentage points | Refined shape MAC | Mesh status |
|---|---:|---:|---:|---:|---:|---:|---|
| In-plane bending | 1 | 0.05293 | -1.66487 | 1.66487 | -0.05203 | 0.9998768 | ACCEPTED |
| Out-of-plane bending | 1 | 0.03465 | -1.36648 | 1.36648 | -0.03417 | 0.9998637 | ACCEPTED |
| Torsion | 1 | 0.06338 | -4.40336 | 4.40336 | -0.06055 | 0.9970807 | ACCEPTED |
| In-plane bending | 2 | 0.05951 | -1.67454 | 1.67454 | -0.05848 | 0.9993307 | ACCEPTED |
| Out-of-plane bending | 2 | 0.03063 | -1.55886 | 1.55886 | -0.03015 | 0.9993788 | ACCEPTED |
| In-plane bending | 3 | 0.06685 | -1.74302 | 1.74302 | -0.06564 | 0.9981234 | ACCEPTED |
| Torsion | 2 | 0.06880 | -4.97569 | 4.97569 | -0.06533 | 0.9902147 | ACCEPTED |
| Axial M-H acoustic | 1 | 0.01649 | -0.47666 | 0.47666 | -0.01641 | 0.9972084 | ACCEPTED |

The additional grid reduces the measured model differences by.01641-.06564
percentage points, while leaving their family ordering and scale intact.
In-plane modes1-3 now differ by1.66487,1.67454,1.74302%; out-of-plane modes1-2
by1.36648,1.55886%. Both bending planes use the inherited distinct inertias,
section rotations and solid-face restraint; no shear coefficient is fitted.

Twist differences remain largest:4.40336/4.97569%, with mesh changes.06338/.06880%.
Refined axial-warp squared mass fractions.00192689/.00722057 and total section
residual fractions.00195705/.00745194 are close to the earlier diagnostics.
Condensed static sectional warping in C_T, absence of a separate dynamic warping
field and differing end restraint are plausible contributors, not proven unique
causes. Small mass fractions do not bound a stiffness/frequency effect; C_T is
not corrected from FEM.

First acoustic M-H difference is now.476662%, mesh change.0164918%, axial section
MAC=.9999644. The effective thickness-strain diagnostic c_eff has shape
MAC=.9988310 and maximum aligned difference.7799998 versus reference scale9.1602715
(about8.52%). These are arbitrary unit-mass LINEAR eigenvector amplitudes, not
finite-motion strains. The41-slab projection and full-end restraint qualification
remain; no boundary zone is omitted. c_eff is not an exact FEM M-H coordinate or
a direct Cartesian DOF. The contraction-dominated branch remains outside the
unchanged frequency window.

### Reliability and readiness for a bounded FEM-2 test

The observed model differences persist as the mesh changes become substantially
smaller. The new level strengthens their interpretation as systematic finite-grid
1D/3D discrepancies in this declared linear window. It does not uniquely separate
shear stiffness, section deformation, constrained warping, local face-clamp effects
and other reduced-model assumptions, or establish a continuum reference.

All eight meet the declared linear numerical criterion. The available four-grid
evidence is sufficient as the LINEAR baseline for a limited nonlinear static
verification, IF separately authorized; no fifth linear grid is required solely
to satisfy this completed preset criterion. Nonlinear static mesh convergence,
load/control and constitutive/BC comparison contracts would still require their
own evidence. FEM-1R neither executes FEM-2/FEM-3 nor proves physical accuracy of
nonlinear V0. No new scientific stage is started automatically.

### Reproduction, artifacts, figures and stop

Bundle: `results/nlsp_linear_rectangular_3d_fem_refinement/63d44daae533389c/`.
Old data are referenced by parent path/hashes rather than copied large arrays.
Continuation identity covers parent manifest, geometry/material, one mesh size,
24-mode/window/matching policy, source-helper and executable/runtime hashes,
new CLI/config hashes, Python/dependencies and resource/semantic policy.
New GEO/MSH/INP/DAT/FRD, logs, complete nodal vectors, mesh/section/assignment data,
`four_mesh_comparison.csv`, `full_refined_spectrum.csv`, axial diagnostic,
resource counters, execution-code copy, provenance and manifest are retained.

```powershell
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py --check-source
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py --run-fem
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py --report-only results/nlsp_linear_rectangular_3d_fem_refinement/63d44daae533389c
python scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py --plot-only results/nlsp_linear_rectangular_3d_fem_refinement/63d44daae533389c
```

Matching cached compute and report/plot start no Gmsh/CalculiX,1D roots or BVP.
A corrupt parent or unmanifested attempted job is a stop condition, not an
instruction to regenerate history or retry a new mesh automatically.
[Targeted tests](../../tests/test_nlsp_linear_rectangular_3d_fem_refinement.py)
use synthetic fixtures and saved evidence, with no extra real FEM jobs:59 tests
PASS (3.52s). Sixteen selected historical FEM1 regressions PASS (1.32s); the other
33 were not run to avoid new preflight roots. Cache replay of matching --run-fem,
report-only and plot-only confirms zero new Gmsh/CCX/1D/BVP calls and unchanged
figure hashes. All15 frozen source/model/config/test hashes, HEAD and staging
index are preserved. Link and whitespace verification are recorded separately.

Exactly two PDF+PNG figures: `four_mesh_convergence` shows all four actual grids;
`refined_model_difference_and_mesh_change` separates model differences from
successive numerical mesh changes. Historical mode-shape figures are reused;
no new identity-problem figure is needed.

| Status | Outcome |
|---|---|
| NLSP_FEM1R_SOURCE_PRESERVATION | PASS |
| NLSP_FEM1R_NEW_MESH_QUALITY | PASS |
| NLSP_FEM1R_MODAL_EXECUTION | PASS |
| NLSP_FEM1R_MODE_IDENTIFICATION | PASS |
| NLSP_FEM1R_MESH_CONVERGENCE | PASS |
| NLSP_FEM1R_ALL_FAMILY_COMPARISON | PASS |

[NLSP-D10](../memory/decisions.md#nlsp-d10) and
[NLSP-K11](../memory/knowledge.md#nlsp-k11) record this separately authorized
continuation and bounded result. LONG CLOSED, EB/RLB-KV PAUSED, angular same-clamp
UNAVAILABLE, prepared planar strict PARTIAL and physical sanity
DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS are preserved. Linear convergence PASS
is not universal1D applicability or nonlinear validation. No new1D calculation,
repeat of the original three modal jobs, fifth mesh, parameter fitting/map,
nonlinear static/dynamic or change to physics/BC occurred. Stop after FEM-1R.
