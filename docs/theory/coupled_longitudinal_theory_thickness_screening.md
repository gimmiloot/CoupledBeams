# Coupled longitudinal theories: bounded thickness screening

2026-10-07. This is a diagnostic comparison of the adopted 1D theories,
with independently sorted positions at each fixed (beta,s_h). MH is a
comparison reference, not exact truth. The previous
[G20 hierarchy screening](coupled_longitudinal_theory_hierarchy_screening.md)
and its immutable bundle remain unchanged.

## Question and scope

Does increasing rectangular thickness make the longitudinal-theory
differences more pronounced for the same two-arm frame? Only five thicknesses
and beta0/45/90 are considered. There is no target percentage, applicability
threshold, fitted law, threshold-thickness search or extension beyond s_h=2.

Inputs: E=rho=1, nu=0.3, b=0.20 along y, h along z, total L=1,
L1=L2=a=0.5, kappa=5/6. Bending is in x-z, I=Iy=bh^3/12. s_h=h/h0,
h0=0.05. A=bh. Total mass rho bhL grows with h; section/material are not
rescaled to preserve mass. Solver mass normalization fixes eigenvector
amplitude within each case only, not physical mass across geometries.

Primary frequency remains f*=f L/sqrt(E/rho); no additional division by h.
The normalized G20 input does not represent dimensional steel frequencies.

## Analytic scaling audit

This is exact coefficient algebra within the already selected planar
formulation, not a new theory or a fitted spectral law. At fixed E,rho,nu,b,a,

\[
 A=bh,\quad I=bh^3/12,\quad I/A=h^2/12,
 \qquad q_h=I/(Aa^2)=h^2/(12a^2).
\]

q_h and derived slenderness lambda_h=a/h are reporting quantities, not new
solver parameters; lambda_h is distinct from the project's frequency Lambda.
Existing eta,tau,mu names are not reused for thickness.

| Quantity | h dependence | Relevant normalized ratio |
| --- | --- | --- |
| A,m=rho A,EA,C=EA/(1-nu^2),S=kappa GA | h | EA/m=E/rho; C/m=E/[rho(1-nu^2)]; S/m=kappa G/rho |
| I,J_RL=nu^2 rho I,j=rho I,H=kappa GI,B=EI,r=rho I | h^3 | I/A=h^2/12 |
| Planar Love gradient inertia | h^3 | J/m=nu^2 h^2/12 |
| MH lateral inertia | h^3 | j/m=h^2/12 |
| MH gradient/normal ratio | h^2 | H/C=kappa(1-nu)(I/A)/2 |
| MH gradient/inertia ratio | h^0 | H/j=kappa G/rho |
| MH normal/inertia ratio | h^-2 | C/j=E/[rho(1-nu^2)] A/I |
| Tim bending stiffness / translation mass | h^2 | B/m=(E/rho) I/A |
| Tim rotary/translation mass | h^2 | r/m=I/A |
| Tim bending/shear ratio | h^2 | B/S=E/(kappa G) I/A |
| Tim shear/rotary ratio | h^-2 | S/r=kappa G/rho A/I |

Both optical cutoffs therefore scale as h^-1:

\[
 f_{c,MH}=\frac{1}{2\pi}\sqrt{C/j},\qquad
 f_{c,Tim}=\frac{1}{2\pi}\sqrt{S/r}.
\]

For an explicit dimensionless check, set xi=x/a, tbar=t sqrt(E/rho)/a,
u=a ubar, w=a wbar. The source energies give

\[
 \bar u_{\bar t\bar t}=\bar u_{\xi\xi}\quad(E),
\]
\[
 \bar u_{\bar t\bar t}-\nu^2 q_h\bar u_{\xi\xi\bar t\bar t}
 -\bar u_{\xi\xi}=0\quad(RL),
\]
\[
 \bar u_{\bar t\bar t}=\frac{\bar u_{\xi\xi}+\nu c_\xi}{1-\nu^2},
 \quad q_h c_{\bar t\bar t}=\kappa\frac GE q_h c_{\xi\xi}
 -\frac{c+\nu\bar u_\xi}{1-\nu^2}\quad(MH),
\]
\[
 \bar w_{\bar t\bar t}=\kappa\frac GE(\bar w_{\xi\xi}-\theta_\xi),
 \quad q_h\theta_{\bar t\bar t}=q_h\theta_{\xi\xi}
 +\kappa\frac GE(\bar w_\xi-\theta)\quad(Tim).
\]

Thus q_h is a natural refined-inertia/gradient scale. This does **not** imply
that every finite coupled frequency difference is proportional to h^2.
In particular, static MH equilibrium has constant N and
H c''=EA c+nu N. Its homogeneous end layer has length

\[
 \ell_{MH}=\sqrt{H/(EA)}
 =h\sqrt{\kappa/[24(1+\nu)]}.
\]

The independently resolved c=0 clamp makes the small-q_h limit singular:
the bulk c=-nu u' condition generally cannot meet it without end layers.
Layer length scales linearly with h. This is an exact static reduction,
not a claimed asymptotic law for frame eigenvalues. Geometry, changing
sorted prefix and joint interaction can further modify observed differences.

Width b multiplies every local mass/stiffness coefficient by the same
factor. It cancels in the ratios above and, for equal homogeneous widths
and no added nodal mass/spring, in homogeneous force-equilibrium rows.
Therefore the relevant dimensionless spectra depend on h and span scales,
not on a common b factor within this planar 1D formulation. Keeping b fixed
isolates the selected thickness change; this is not a statement about width
independence of full 3D elasticity. Exact Fraction tests check cancellation
and q_h scaling; no new symbolic dependency was installed. Lean/SymPy are
unavailable, so algebra is derived directly and checked exactly/numerically.

## Geometry, models and boundary qualification

| s_h | h | h/a | lambda_h=a/h | b/h | q_h=I/(A a^2) |
| --- | --- | --- | --- | --- | --- |
| 1.00 | 0.0500 | 0.100 | 10.0000 | 4.0000 | 0.000833333 |
| 1.25 | 0.0625 | 0.125 | 8.0000 | 3.2000 | 0.001302083 |
| 1.50 | 0.0750 | 0.150 | 6.6667 | 2.6667 | 0.001875000 |
| 1.75 | 0.0875 | 0.175 | 5.7143 | 2.2857 | 0.002552083 |
| 2.00 | 0.1000 | 0.200 | 5.0000 | 2.0000 | 0.003333333 |

The same verified Timoshenko block is used by all three models at each h:
rho A, rho I, EI and kappa GA are identical. Production kappa remains5/6,
MH inertia factor1. Elementary has H=J=0; planar Love has H=0,
J=nu^2 rho Iy and harmonic N=(EA-J omega^2)U'. No polar Ip or additional
slope clamp is introduced. MH retains C=EA/(1-nu^2), H=kappa GI, j=rho I,
the Jang reduced normal block and the validated c/R joint closure.
Source Rucka/Jang/Ng/Fernandes variants and Bishop remain unchanged.

Common project beta and dual local/global maps are unchanged. A/B clamp
u=w=theta=0; MH additionally clamps c=0 and uses c1=c2,R1+R2=0 at the joint.
The centroid/bending clamps are shared, but cross-sectional contraction
constraints are not identical. c1=c2 is variational reduced-frame closure,
not direct 3D elasticity of a welded finite joint region.

Frozen equations.tex and src/my_project/analytic determinants are unchanged.
This diagnostic calls the previously verified source-energy and H=0 kernels.
No new mathematical baseline block was installed; coefficient and beta0
direct-profile checks found no theory/code mismatch.

## Numerical policy and available analytic domain

Configured initial frequency ceiling remains f*=3.75. The unchanged finite
basis is valid below optical cutoff. At h=0.0875 and0.10 the Tim cutoff falls
below3.75. Rather than introduce another above-cutoff solver, the workflow
certifies only the required prefix13 within the available domain. It does
not claim an inventory covering the unused tail up to3.75.

Predeclared catalogs use Tim simply-supported index n+0.9 and MH scalar
index n+0.5. Tim catalogs stay below0.99 optical cutoff; MH catalogs stay
below0.99 of the Young lower-form contraction cutoff, so that min-max count
is sharp here. A count query window ends at the smaller of configured
ceiling and0.98 of each catalog ceiling. These are representation/conditioning
guards, not physical applicability or spectral-difference thresholds.
Catalog scans can be doubled once on the same interval if count saturation
fails; every attempt is retained. No such retry was needed.

Full boundary determinants and independent energy/Schur counts, pole-query
conditioning, 400 deterministic intervals and bounded count-guided subdivision
come from the existing hierarchy/general helpers. A missing guard permits
one configured ceiling expansion by1.5. If the unchanged analytic domain
still cannot cover it, the case stays incomplete rather than changing physics.
No expansion or unresolved case occurred. Root brackets, counts, subdivisions
and failed-interval lists are saved for every case.

All accepted tolerances are unchanged: boundary/joint/PDE1e-9, total-energy
quotient5e-8, mass Gram5e-7, nonzero singular condition1e8. The latter is
saved even when a close pair worsens conditioning. Total energy remains
internal eigenproblem validation; no energy fractions or mode classes occur.

The local frequency-map-v1 instance declares certified_audit because this
task explicitly requests every root certificate, direct beta0 profiles and
fixed-case overlaps. This does not change the ordinary fast_plot default.
Same-config rerun checks hashes and computes zero roots; plot-only uses
saved tables. Baseline nine cases are reused read-only after verifying their
code/source/config/artifact identity.

## Metrics and interpretation

delta_E=(fE_k−fMH_k)/fMH_k and delta_RL=(fRL_k−fMH_k)/fMH_k are signed
model-to-model differences at independently sorted k. Absolute values,
medians, maxima/positions and secondary RL/E differences are retained.
Maxima are over first12 only; root13 is a completeness guard.

At the same (beta,s_h), geometric d=u t+w n and theta overlaps use the
unchanged separate 12x12 metrics, common Gauss200 quadrature and machine-scale
small-norm policy. Row/column best, second and margins are retained. No
across-h or across-beta shape comparison, assignment or root reordering is
performed. A non-diagonal best match would be a flag, not a failure.

For MH, D_c=||c+nu u'||/(||c||+||nu u'||), with three norms and
NOT_DEFINED_SMALL_FIELD handling, is unchanged. It describes departure
from quasistatic Poisson contraction, including dynamics and end/joint
layers, not a modal class. Adjacent gaps are context only.

For each h, all three beta0 first13 are compared with independently assembled
direct fixed-fixed rod boundary problems, including mass-normalized profiles
and transparent joint residuals. Direct operator-block indices are used only
for this separated beta0 benchmark, not for classifying coupled modes.
Fresh direct checks and Gauss200/300 overlap/D_c checks also cover the thickest
case. No new general solver or tracking machinery was created.

## Results

Canonical local bundle: `results/coupled_longitudinal_theory_thickness_screening/7be0fce968fd2b35/`.

| Status | Result |
| --- | --- |
| COUPLED_THICKNESS_SCALING_AUDIT | PASS |
| COUPLED_THICKNESS_ROOT_INVENTORY | PASS |
| COUPLED_THICKNESS_FIXED_CASE_OVERLAP | COMPLETE |
| COUPLED_THICKNESS_SCREENING | COMPLETE |

| beta | s_h | max abs delta_E (%) [k] | max abs delta_RL (%) [k] | max D_c [k] | min diagonal O_d E / RL | non-diagonal E / RL rows |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 1.00 | 0.141030 [5] | 0.175141 [11] | 0.058849 [5] | 0.999985200 / 0.999985200 | 0 / 0 |
| 0 | 1.25 | 0.174481 [5] | 0.227916 [11] | 0.065961 [5] | 0.999977344 / 0.999977344 | 0 / 0 |
| 0 | 1.50 | 0.207249 [4] | 0.284218 [10] | 0.072430 [4] | 0.999967940 / 0.999967940 | 0 / 0 |
| 0 | 1.75 | 0.239354 [4] | 0.343902 [10] | 0.078409 [4] | 0.999956990 / 0.999956990 | 0 / 0 |
| 0 | 2.00 | 0.270814 [4] | 0.479289 [12] | 0.086181 [12] | 0.999893434 / 0.999893434 | 0 / 0 |
| 45 | 1.00 | 0.127346 [5] | 0.164384 [11] | 0.085393 [7] | 0.999984392 / 0.999980643 | 0 / 0 |
| 45 | 1.25 | 0.089428 [4] | 0.163090 [11] | 0.096945 [7] | 0.999941920 / 0.999761713 | 0 / 0 |
| 45 | 1.50 | 0.177129 [4] | 0.234490 [10] | 0.101161 [7] | 0.999953799 / 0.999832665 | 0 / 0 |
| 45 | 1.75 | 0.219876 [4] | 0.314687 [10] | 0.102156 [8] | 0.999957718 / 0.999731846 | 0 / 0 |
| 45 | 2.00 | 0.253385 [4] | 0.376504 [10] | 0.109122 [5] | 0.999947380 / 0.999572563 | 0 / 0 |
| 90 | 1.00 | 0.087503 [5] | 0.125413 [11] | 0.085751 [7] | 0.999983038 / 0.999969837 | 0 / 0 |
| 90 | 1.25 | 0.100589 [4] | 0.127183 [9] | 0.096836 [7] | 0.999975901 / 0.999926474 | 0 / 0 |
| 90 | 1.50 | 0.146563 [3] | 0.205211 [9] | 0.101359 [7] | 0.999973561 / 0.999939119 | 0 / 0 |
| 90 | 1.75 | 0.179916 [3] | 0.254300 [9] | 0.103238 [7] | 0.999965750 / 0.999912571 | 0 / 0 |
| 90 | 2.00 | 0.205757 [3] | 0.284585 [9] | 0.112479 [5] | 0.999956617 / 0.999869136 | 0 / 0 |

Global E/MH maximum0.2708143008% at beta0,s_h2,k4 (signed negative);
global RL/MH maximum0.4792888457% at beta0,s_h2,k12 (signed negative).
These are discrete-grid maxima over first12; no extremum/threshold search.

### Descriptive growth ratios, without a fit

| s_h | q_h/q_h0=s_h^2 | beta0 E / RL | beta45 E / RL | beta90 E / RL |
| --- | --- | --- | --- | --- |
| 1.00 | 1.0000 | 1.00000 / 1.00000 | 1.00000 / 1.00000 | 1.00000 / 1.00000 |
| 1.25 | 1.5625 | 1.23719 / 1.30133 | 0.70224 / 0.99213 | 1.14955 / 1.01411 |
| 1.50 | 2.2500 | 1.46954 / 1.62280 | 1.39093 / 1.42648 | 1.67495 / 1.63628 |
| 1.75 | 3.0625 | 1.69719 / 1.96358 | 1.72660 / 1.91434 | 2.05612 / 2.02770 |
| 2.00 | 4.0000 | 1.92027 / 2.73659 | 1.98974 / 2.29039 | 2.35144 / 2.26918 |

The max statistic is not the same mode across thicknesses. Its sorted
position changes (e.g. beta0 E:k5->k4; RL:k11->k12). At s_h2, E growth
is1.9203/1.9897/2.3514 and RL growth2.7366/2.2904/2.2692 for beta0/45/90,
rather than4. beta45 has an initial decrease at1.25, followed by growth.
Thus the natural h^2 coefficient scale does not give a universal spectral
difference law over this finite prefix. Static layers and joint geometry
provide mechanisms for a different response; no fitted exponent is inferred.

### Shape correspondence, contraction and close gaps

All displacement row/column best matches for all three theory pairs are
diagonal. All informative theta best matches are diagonal as well. There
are no non-diagonal cases to list and no additional overlap figure was made.
Overall minimum diagonal O_d is0.9998934335 E/MH and0.9995725627 RL/MH.

Max defined D_c grows with h for all three angles. Endpoint ratios are
about1.46 (beta0),1.28 (beta45),1.31 (beta90), with changing sorted positions.
The global maximum0.1124793677 is at beta90,s_h2,k5. This is a modest
increase in departure from quasistatic contraction, not a modal type or
an error-to-truth criterion. The9 baseline null fields remain null through
s_h1.75;8 first12 fields are undefined at beta0,s_h2. Other cases retain
defined norms/ratios. Their exact values are in contraction.csv.

| Model | Smallest adjacent gap | beta | s_h | Positions |
| --- | --- | --- | --- | --- |
| elementary | 0.000259556362 | 90 | 1.25 | 9--10 |
| rayleigh_love_planar | 0.0002263357442 | 90 | 1.25 | 9--10 |
| mindlin_herrmann | 3.672664694e-05 | 90 | 1.25 | 9--10 |

The closest pair is beta90,s_h1.25,k9?10. MH gap=3.67266469e-5
(0.00367266%), E0.02595564%, RL0.02263357%. It causes the largest
nonzero SVD condition52194.6 (MH position9), but all fixed-case best shapes
remain diagonal. No crossing, avoided crossing or resonance is concluded.
Every per-case/model minimum gap and positions are also in summary.json.

By the maximum first12 difference, RL is farther from MH than E in all
15 cases, and this ordering persists with thickness. The stored180 rows
show no resolved individual improvement by RL beyond numerical noise;
this is limited to these inputs and the selected closure/BC, not a claim
of general superiority of elementary axial theory.

### Completeness, direct checks and numerical quality

| s_h | Tim cutoff f* | usable search ceiling f* | largest guard13 f* | root-count range over9 cases |
| --- | --- | --- | --- | --- |
| 1.00 | 6.242570465 | 3.750000000 | 1.805278177 | 13--24 |
| 1.25 | 4.994056372 | 3.750000000 | 1.990214663 | 22--23 |
| 1.50 | 4.161713644 | 3.750000000 | 2.118643153 | 22--22 |
| 1.75 | 3.567183123 | 3.294270785 | 2.202065301 | 19--20 |
| 2.00 | 3.121285233 | 2.801065334 | 2.251025337 | 16--16 |

All45 first12+guard13 inventories are complete in their saved ranges.
Baseline MH beta0 reuses its prefix13 certificate, not an older full3.75
tail claim. No guard expansion, catalog retry, failed numerical interval
or unresolved case remains. Count-guided subdivisions are saved (maximum9).
The shorter sub-cutoff ranges are not silent exclusions of requested roots:
their independent count certificates and guard13 prove the required prefixes.

All15 direct beta0 checks match13 direct fixed-fixed frequencies and forms.
Maximum relative frequency difference6.3465e-12, dimensionally scaled
kinematic L2 difference8.4010e-12. Joint force/moment and c/R transmission
pass separately. The direct separated-block inventory is a same-geometry
verification, not tracking between thicknesses.

| Scaled diagnostic, all accepted roots | Maximum |
| --- | --- |
| clamp_scaled_residual | 5.05352e-12 |
| equation_scaled_residual | 1.14439e-12 |
| singular_ratio | 2.57789e-12 |
| nonzero_singular_condition | 52194.6 |
| energy_relative_error | 7.19391e-12 |
| d_X | 3.17113e-12 |
| d_Y | 2.20542e-13 |
| theta | 6.13518e-12 |
| F_X | 6.38624e-12 |
| F_Y | 4.21548e-12 |
| M_node | 1.31164e-13 |
| mass_gram_max_error | 1.47777e-10 |
| c | 6.77208e-14 |
| R_node | 6.80039e-16 |

Accepted tolerances were not relaxed. The poorer conditioning of the
identified close pair is disclosed; its residuals, normalization and gaps
are retained. It is still below the existing1e8 nonzero-condition bound.

237 targeted/regression checks pass:35 new thickness tests,198 unchanged
hierarchy/MH general/beta0/single/source/Bishop checks, plus4 rectangular
Tim regressions (7 unrelated tests deselected). Exact moments/ratios and
width cancellation, shared Tim bases, no coefficient fit, direct recovery
at baseline/intermediate/thickest h and all saved checks, overlap sign/null
handling, D_c and Gauss200/300 convergence are verified. Full expensive
research suite was not run.

### Reproduction and artifacts

```powershell
python scripts/analysis/screen_coupled_longitudinal_theory_thickness.py --check-sources
python scripts/analysis/screen_coupled_longitudinal_theory_thickness.py --compute
python scripts/analysis/screen_coupled_longitudinal_theory_thickness.py --plot-only results/coupled_longitudinal_theory_thickness_screening/7be0fce968fd2b35
python -m pytest tests/test_coupled_longitudinal_theory_thickness_screening.py -q
```

Working interpreter: D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe
(Python3.12.4, NumPy2.1.3, SciPy1.15.2, Matplotlib3.9.2). No dependency installed.
The new CLI is justified by the scaling/direct-profile/trend output contract;
it calls existing root, mode, count and reporting helpers rather than
duplicating them. Only one CLI/config/test file and canonical note were added.

Manifest retains source/input/code/version/baseline hashes, Git and command.
The baseline hierarchy bundle9a35ff23c3f43c66 is read-only and unchanged.
Artifacts: scaling_audit.json, all45 inventories with coefficients/brackets/
counts/diagnostics, direct_beta0_checks.json,180-row frequencies.csv,
495 adjacent gaps, full fixed-case overlaps, D_c/norm data, five full-profile
CSVs and summary.json. No output field assigns energy-based modal types.
Three figures: thickness_max_differences.png, thickness_position_differences.png
(beta45) and thickness_contraction_diagnostic.png. Connected lines guide
the eye between five computed points; no intermediate h was computed.

### Answers to the bounded scientific questions

1. E/MH maxima generally grow toward s_h2; beta45 initially decreases.
2. RL/MH maxima also grow overall; beta45 initially nearly unchanged/slightly lower.
3. Growth ratios1.92?2.35 (E) and2.27?2.74 (RL) at2 differ from s_h^2=4.
4. Endpoint E/MH maximum decreases from beta0 to90; RL/MH also largest at0.
   Relative growth and max positions depend on angle; no angle extrema searched.
5. No non-diagonal fixed-case best correspondence appeared.
6. Max D_c increases modestly, to0.11248, with unchanged diagnostic semantics.
7. RL remains farther from MH by max difference for all15 cases.
8. Larger h strengthens refinement effects, but all max differences remain
   sub-percent and correspondence stays diagonal. This is insufficient reason
   to extend thickness automatically; the existing close pair is the more
   localized candidate if further analysis is scientifically needed.


## Limits and recommendation

Increasing h reduces slenderness: the final members have a/h=5 and b/h=2.
These results compare the adopted 1D models; they do not validate them
against 2D/3D elasticity. No source threshold was assigned to these ratios.
There was no energy classification, across-case tracking, automatic spectral
reordering, fitted coefficient/power law, applicability threshold, threshold
thickness search, beta refinement, unequal-arm/section/material scan, FEM
or nonlinear derivation. s_h>2 was not calculated.

The bounded result gives stronger but still sub-percent model differences
and preserves diagonal fixed-case correspondence. It does not compel an
automatic extension to thicker rods. If a later targeted investigation is
authorized, the existing beta90,s_h1.25 positions9–10 close pair is a more
specific numerical/structural question than an open-ended thickness search;
it is not presently a veering claim. Any extension beyond2 needs a separate
scientific decision and an explicit check of the available 1D/basis domain.
