# Anisotropic Rods — Chapter-2 source line consolidated

```text
Status: active-research
Source assimilation: completed for the selected Chapter-2 single-rod model
Primary source: registered
Free-free Chapter-2 reproduction: completed — PASS_WITHIN_GRAPH_RESOLUTION
Cantilever Chapter-2 reproduction: completed
Cantilever source clamp: BOOK_SLOPE_CLAMP_CONFIRMED
Canonical source-consistent formulation: state_corrected
Printed sign variant: diagnostic only
Rigid angular-joint gate: completed — PASS
Coupled-rod code: small elastic diagnostic pilot only — PASS
Final coupled-rod model or parameter study: not started
Rectangular orthotropic Euler--Bernoulli--Saint-Venant gate: PARTIAL_PASS after targeted refinement
Unequal-thickness 1D phase: COMPLETE; overall: PARTIAL_PASS; UT-0: PASS; UT-1: PARTIAL_PASS; UT-1a: INCONCLUSIVE; UT-2: PASS; UT-3: PASS
Limited 3D FEM anchor: DESIGN_ONLY; 3D-A0: BLOCKED_MATERIAL_AND_SOLVER
Supervisor Chapter-2 rectangular-rod figures: completed — PASS
Fast Chapter-2 beta-sweep oracle validation: completed — PASS
Extended supervisor Figures 5–8: completed — PASS
Small-theta supervisor Figures 9–12: completed — PASS
AP-0 theta=2 small-grid spectral-applicability screening: completed — PASS
AP-1 theta=5 same-grid spectral-applicability screening: completed — PASS
AP-2 theta=3/4 intermediate-angle same-grid screening: completed — PASS
Production anisotropic API: not started
Independent rectangular EB + torsion 1D FEM: completed for the finite gate
```

## Established single-rod baseline

The selected source line is one rectangular monoclinic Timoshenko rod with
generalized torsion. The internally consistent `state_corrected` formulation
restores `I_y` in equation (2.1). `eliminated_corrected` is its independent
equivalence check, while `eliminated_printed` preserves the signs printed
after equation (2.16) only for historical diagnosis.

The free-free Figure-2.2 benchmark is
`PASS_WITHIN_GRAPH_RESOLUTION`. The cantilever Figure-2.8 source check is
`BOOK_SLOPE_CLAMP_CONFIRMED`: the book's calculated curves use
`w(0)=0`, `w'(0)=0`, `Phi(0)=0`. This identifies the source boundary
condition; it does not make the alternative section-rotation clamp physically
invalid. The separate rigid-joint gate has now derived that internal joint
from physical vector compatibility and equilibrium.

## Rigid angular-joint pilot

The two-arm helper retains the book state
`[w, psi, Phi, Q, M, M_T]^T`, the project bases `e_z,t_i,n_i`, and the
accepted `state_corrected` arm equations. Its `6 x 12` joint matrix passed the
old-notation sign gate, rank and random-state checks, virtual work, requested
angle limits, orthotropic block separation, equal-arm exchange, and the
independent `beta=0` straight-rod spectrum check. The small real-elastic
HMS/DX-209 pilot at `beta=0,30,90 deg` is `PASS`; it is not a stable baseline
or a final parameter study.

## Rectangular orthotropic EB validation

The `theta=0` comparator reuses the existing rectangular `Cbar`, unchanged
`J_book`, book state, and project bases. Unequal-length invariance, exact
fixed--fixed bending/torsion families, the Timoshenko slender limit, FEM
structure, and the original fixed-elements-per-arm convergence checks passed.
The original gate is `PARTIAL_PASS` because 64 elements per arm gives unequal
physical element lengths and misses the strict first-three `1e-5` target
(`5.18e-5` worst case). The targeted proportional sequence through
`(N_1,N_2)=(64,192)` closes that accuracy threshold with the raw value
`6.18e-6`, without changing any model coefficient or threshold, but its
separate unchanged monotonic-refinement gate fails for bending mode 1 at the
dense-eigensolver conditioning floor. Its targeted status is
`FAIL_CONVERGENCE_ORDER`, so the overall status remains `PARTIAL_PASS`. This
is a project comparator, not a source reproduction, general anisotropic EB
model, final coupled baseline, unequal-thickness study, or 3D validation.

## Unequal-thickness UT-0 through UT-3

The isolated UT-0 gate constructs independent `a=4,5,6 mm` rectangular
`RodPoint` objects at `theta=0` and verifies section properties, arm-specific
Timoshenko/EB state coefficients, and recovery of each physical transfer from
the existing scaling. For the only spectral case `(a_1,a_2)=(4,6) mm` at
`beta=0`, independent script-local stepped determinants agree with the coupled
first-three sorted spectra to `1.42e-13` relatively for Timoshenko and
`7.26e-14` for EB. UT-0 is `PASS`, while overall unequal-thickness validation
remains `IN_PROGRESS`; other angles, the swapped case, root 7, FEM, shapes,
tracking, and parameter studies were not run.

UT-1 extends the same script-local validation to seven roots for `(5,5)`,
`(4,6)`, and `(6,4) mm`. All 12 coupled/direct-stepped spectra, root quality,
the homogeneous equal-thickness regression, analytic stepped-order exchange,
mesh-64 FEM accuracy, and the fixed `(16,32,64)` aggregate convergence gate
pass. UT-1 is `PARTIAL_PASS` only because the independent mesh-64 FEM exchange
difference is `5.76e-8`, above its fixed `1e-8` threshold. No FEM coefficient,
continuum formula, or threshold was changed.

UT-1a isolates that FEM exchange difference at matrix level for only the two
`N_1=N_2=64`, `beta=0` systems `(4,6)` and `(6,4)`. The a-priori full and
signed-reduced arm permutations give exactly zero saved residual for
permutation, reduction-intertwining, full-matrix exchange, and forward/reverse
reduced congruence. The canonicalized spectra also agree exactly. The audit is
nevertheless `INCONCLUSIVE_NUMERICAL_AUDIT`: native/transported backward
residuals reach `7.28e-7 > 1e-8`, and the transported Rayleigh difference is
`1.48e-7 > 1e-10`. No assembly asymmetry was detected within the declared
matrix thresholds, but the failed numerical-quality gates prevent attributing
the native difference conclusively to DOF-order-sensitive eigensolver error.
UT-1 remains `PARTIAL_PASS` and the overall status remains `IN_PROGRESS`.

UT-2 is the separate `beta=+/-30 deg` angular-joint gate for `(5,5)`, `(4,6)`,
and `(6,4) mm`. All 12 seven-root Timoshenko/EB spectra pass the unchanged
root-quality checks. Reflection and arm-relabeling continuum spectra agree to
`0` and `3.00e-13`, respectively. Six independent mesh-64 EB FEM assemblies
pass the fixed first-three, first-six, and root-7 accuracy gates. The declared
full/reduced reflection and relabeling maps satisfy their endpoint, reduction,
and stiffness/mass identities within the fixed `1e-13` relative matrix
threshold; the largest relabeling relative matrix residual is `1.38e-19`.
UT-2 is `PASS`. Native FEM spectral relabeling reaches `2.01e-7`, but remains
diagnostic because the a-priori matrix congruence and FEM accuracy gates pass.
No beta-zero status is changed. At the UT-2 boundary the overall validation
remained `IN_PROGRESS`, pending the separately scoped quarter-turn gate.

UT-3 closes the one-dimensional phase at `beta=+/-90 deg`. Exact integer
quarter-turn joint matrices agree with unchanged `J_book(±pi/2)` to a maximum
entry residual of `6.12e-17`; 84 standard plus 42 independent exact-limit
continuum roots are accepted, and their maximum spectral difference is
`1.96e-14`. Six mesh-64 EB FEM systems pass structure and fixed accuracy
gates. Reflection matrix residuals are zero and the largest relative
relabeling residual is `6.12e-17`, below `1e-13`. UT-3 is `PASS`; the 1D
phase is `COMPLETE`, while its overall status is `PARTIAL_PASS` because UT-1
and the declared UT-1a numerical audit retain their historical limitations.

## Navigation

- [Yartsev source note and local scan map](source_note_yartsev_2024.md)
- [Literature source-index entry](../literature/source_index.md#yartsev_2024_coupled_composite_structures)
- [Free-free Figure-2.2 reproduction note](yartsev_ch2_single_rod_reproduction.md)
- [Cantilever Figure-2.8 reproduction note](yartsev_ch2_cantilever_reproduction.md)
- [Chapter-2 notation translation](yartsev_ch2_notation_translation.md)
- [Rigid angular-joint theory and gate](yartsev_ch2_rigid_angular_joint.md)
- [Rectangular orthotropic EB validation](yartsev_ch2_rectangular_eb_validation.md)
- [Unequal-thickness 1D validation (UT-0 through UT-3)](yartsev_ch2_unequal_thickness_validation.md)
- [Limited 3D FEM anchor design and 3D-A0 readiness](yartsev_ch2_limited_3d_fem_anchor_design.md)
- [Supervisor Chapter-2 rectangular-rod figures](yartsev_ch2_supervisor_figures.md)
- [AP-0/AP-1/AP-2 sampled-theta spectral-applicability screening](yartsev_ch2_spectral_applicability_screening.md)
- [Free-free diagnostic CLI](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py)
- [Cantilever diagnostic CLI](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py)
- [Small coupled-joint pilot CLI](../../scripts/analysis/anisotropic_rods/pilot_yartsev_ch2_coupled_rods.py)
- [Rectangular EB validation CLI](../../scripts/analysis/anisotropic_rods/validate_yartsev_ch2_rectangular_eb.py)
- [Unequal-thickness UT-0/UT-1/UT-1a/UT-2/UT-3 CLI](../../scripts/analysis/anisotropic_rods/validate_yartsev_ch2_unequal_thickness.py)
- [Limited 3D FEM 3D-A0 readiness CLI](../../scripts/analysis/anisotropic_rods/audit_yartsev_ch2_limited_3d_fem_readiness.py)
- [Supervisor-figure CLI](../../scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py)
- [AP-0/AP-1/AP-2 screening CLI](../../scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py)

Generated local evidence, which may be absent from a fresh clone, is kept in:

- [`results/anisotropic_rods/yartsev_ch2_free_free/`](../../results/anisotropic_rods/yartsev_ch2_free_free/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever/`](../../results/anisotropic_rods/yartsev_ch2_cantilever/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/)
- [`results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/`](../../results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/)
- [`results/anisotropic_rods/yartsev_ch2_coupled_joint_pilot/`](../../results/anisotropic_rods/yartsev_ch2_coupled_joint_pilot/)
- [`results/anisotropic_rods/yartsev_ch2_rectangular_eb_validation/`](../../results/anisotropic_rods/yartsev_ch2_rectangular_eb_validation/)
- [`results/anisotropic_rods/yartsev_ch2_unequal_thickness_validation/`](../../results/anisotropic_rods/yartsev_ch2_unequal_thickness_validation/)
- [`results/anisotropic_rods/yartsev_ch2_limited_3d_fem_anchor/design_readiness/`](../../results/anisotropic_rods/yartsev_ch2_limited_3d_fem_anchor/design_readiness/)
- [`results/anisotropic_rods/yartsev_ch2_supervisor_figures/`](../../results/anisotropic_rods/yartsev_ch2_supervisor_figures/)
- [`results/anisotropic_rods/yartsev_ch2_spectral_applicability_screening/`](../../results/anisotropic_rods/yartsev_ch2_spectral_applicability_screening/)

## Current boundary

The rigid-point-joint pilot, finite `theta=0` rectangular EB/1D-FEM gate,
and unequal-thickness UT-0 through UT-3 one-dimensional gates are complete.
The rectangular gate closed its original strict accuracy shortfall under
equal-element-length refinement but was not promoted because the separate
mode-1 monotonic-refinement gate failed. The unequal-thickness 1D phase is
`COMPLETE` with overall `PARTIAL_PASS`: UT-1 retains its historical native
FEM exchange failure and UT-1a remains inconclusive under its fixed numerical
quality gates. Other angles, off-axis material axes, and broad spectrum
angles beyond the explicitly reported `theta=5` and `15 deg` diagnostics, and
broad spectrum studies remain outside the completed scope. The separate 3D-A0
design/readiness audit is complete with `BLOCKED_MATERIAL_AND_SOLVER`: the
local HMS/DX-209 source lacks `E3`, `nu13`, and `nu23`, and no capable local
3D solver stack was detected. Geometry, ideal-joint, clamp, target-parity,
case, mesh, and future acceptance contracts are frozen, but mesh generation
and eigenfrequency calculations remain `NOT_STARTED`. No 3D implementation is
part of the completed 1D phase. The separate supervisor-figure workflow is a
presentation/report layer over the verified Chapter-2 rectangular-rod
builders; it does not alter these 1D or 3D statuses and is explicitly separate
from the circular-isotropic article workflow. Its diagnostic-only fast
beta-sweep reproduces all 6335 saved Figure-2--4 oracle roots within `1e-8`
and reduces the recorded sequential runtime from `860.775 s` to `139.714 s`;
Figures 5--8 then add the prescribed unequal-length, direct unequal-thickness,
weak-monoclinic approximation, and within-Chapter-2 material-rotation plots.
Figures 9--12 continue the Figure-7 diagnostic at local material-axis angles
`theta1=theta2=1,2,3,4 deg`, reusing the Figure-3 theta-zero rectangular EB
baseline exactly and retaining independently sorted spectral positions. The
separate AP-0 diagnostic then screens nine volume-preserving similar-section
geometries at `theta=2 deg` and `beta=0,5,...,90 deg`: all nine remain within
the finite-grid 10% Lambda criterion, while no exact theta boundary,
refinement, tracking, shapes, energy analysis, or FEM is attempted. AP-1
reuses exactly the same geometry/beta grid at `theta=5 deg`: all nine cases
then exceed the 10% criterion, with a descriptive same-pair gap analysis over
1026 observations. AP-2 then evaluates only the prescribed intermediate
sampled orientations `theta=3,4 deg`: all 171 configurations and all nine
families remain within 10% at both sampled angles, while theta=5 retains
66/171 pointwise within and 105/171 exceeding. All sampled point and family
sequences are nondecreasing, but no continuous monotonicity, exact crossing,
or critical angle is inferred.

## Preservation rule

> No existing isotropic equation, determinant entry, sign convention, unknown
> ordering, root solver, or FEM model may be changed to support this research
> direction without a separate consistency audit.
