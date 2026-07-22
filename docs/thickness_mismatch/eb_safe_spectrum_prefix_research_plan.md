# EB Safe-Spectrum-Prefix Research Plan

## Status and Scope

This document records the current diagnostic-only research direction for a
conservative Euler--Bernoulli (EB) certificate for the low sorted spectrum of
the coupled-rod system. It defines the engineering target, the evidence that
may be used, the next computational step, and the quality and cost metrics that
must be reported. It does not introduce a new analytic model or claim a
continuous-domain guarantee.

Project-wide conventions for sorted frequencies, descendant branches, model
applicability, and diagnostic-to-article promotion remain those in the
[project rules](../project_rules.md).

## Implementation Status

The first CSV-first certification postprocessor was implemented on
`2026-07-20` at
[`analyze_eb_safe_prefix_certification.py`](../../scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py).
It performs no root solve. It requires complete `K = 10` sorted rows from both
source applicability workflows, reconstructs canonical EB-only predictors
from saved `Lambda_EB` values, and implements:

- Rule A: one `Pi_EB` limit;
- Rule A-gap: Rule A plus the complete-cluster adjacent EB gap guard;
- Rule B: separate `Pi_shear_EB` and `Pi_rotary_EB` limits;
- Rule C: Rule B plus a complete-cluster adjacent EB spectral-gap guard;
- Rule D: Rule C plus the EB-only `mixed` modal-character guard.

Thresholds are recalibrated on each train fold, while complete geometries are
held out for cross-source, leave-one-mu, leave-one-epsilon, leave-one-beta, and
leave-one-eta checks. The workflow also writes incomplete/warning exclusions,
source-overlap and fixed-source predictor-consistency audits, and raw
postprocessor operation counts. Existing `first8` source-summary fields retain
their original eight-mode meaning; added `firstK` fields support explicit
`--n-reported-modes 10` runs without changing source defaults.

Remaining steps are computational and scientific rather than implementation
claims:

- regenerate both complete source datasets with `K = 10` and the documented
  candidate-root margins;
- review held-out false-safe cases, conservative loss, and cluster triggers on
  those complete datasets;
- decide whether Rules A--D and the Rule A-gap ablation retain enough safe frequencies before introducing
  any new physical composite criterion;
- extend evidence beyond the two current parameter slices only after the
  complete K=10 results have been reviewed;
- add source root-operation counters later if this can be done without changing
  solver behavior.

## Geometry-Only Epsilon Pilot

A deliberately small geometry-only pilot was implemented and run on
`2026-07-20` through the manifest-driven source runner
[`run_eb_epsilon_apriori_pilot.py`](../../scripts/analysis/thickness_mismatch/audits/run_eb_epsilon_apriori_pilot.py)
and the CSV-only postprocessor
[`analyze_eb_epsilon_apriori_pilot.py`](../../scripts/analysis/thickness_mismatch/postprocess/analyze_eb_epsilon_apriori_pilot.py).
The fixed manifest contains 21 selected geometries and compares two candidate
pre-solution indicators:

- `epsilon_0`, which is available directly from the base geometry;
- `epsilon_max`, computed through the verified local-thickness helper.

Both quantities remain possible low-cost screening variables, not replacements
for the EB-based certificate. The pilot kept the sorted `K = 10` target,
calibrated ten exact finite-observation thresholds for each geometry-only rule,
used an explicit predictor-domain fallback to `N_hat = 0`, and compared those
rules with Rules A--D plus a `Pi_EB`/spectral-gap ablation, Rule A-gap. The
production threshold candidate search now gives priority to the running prefix
extrema that actually control certification decisions.

The legacy 2026-07-20 run reported all 21 geometries without root or
candidate-boundary warnings. On
the baseline-reference transfer to the 14 non-baseline cases, the observed
results were:

| rule | false-safe geometries | retained safe EB frequencies | mean conservative loss | X-domain fallbacks |
|---|---:|---:|---:|---:|
| E0-ref | 0 | 0.7863 | 1.786 | 1 |
| Emax-ref | 0 | 0.4017 | 5.000 | 3 |
| Rule A | 0 | 0.8803 | 1.000 | 0 |
| Rule A-gap | 0 | 0.7179 | 2.357 | 0 |
| Rule B | 0 | 0.6752 | 2.714 | 0 |
| Rule C | 0 | 0.6752 | 2.714 | 0 |
| Rule D | 0 | 0.3504 | 5.429 | 0 |

The two matched-`epsilon_max` triplets had `N_true` ranges of two and six
modes, respectively. Thus `epsilon_max` did not collapse the selected pilot
geometries to a common safe prefix, and it performed worse than `epsilon_0` as
the transferred baseline indicator. Across the broader repeated held-out fold
evaluations, calibrated geometry-only rules produced observed false-safe cases.
In addition, many reduced searches for Rules A-gap through D changed between
the 16- and 32-candidate grids. Those rules therefore retain an explicit
calibration-convergence caveat in this pilot.

The current diagnostic recommendation is to retain a cascade rather than
replace the EB certificate: geometry-only screening may be tested as an early
conservative fallback, followed by EB frequencies and gaps, EB-only Pi guards,
and Timoshenko for the uncertified suffix. Any decision to promote `epsilon_0`
or `epsilon_max` requires a larger held-out geometry design. Zero observed
false-safe on these selected finite cases is not a continuous-domain guarantee.

## Refined Straight-System Epsilon Baseline

### Straight-baseline close-root correction

Research step 2 was implemented and run on `2026-07-21` through the stable
diagnostic entry point
[`audit_eb_epsilon_baseline_thresholds.py`](../../scripts/analysis/thickness_mismatch/audits/audit_eb_epsilon_baseline_thresholds.py).
It is restricted to the straight homogeneous system
`(beta, mu, eta) = (0, 0, 0)` and uses the dimensional-frequency discrepancy
based on squared `Lambda`. For every prefix `n=1,...,10`, it evaluates

\[
\Delta_n(\epsilon)=\max_{1\leq k\leq n}\delta_{f,k}(\epsilon)
\]

and defines first loss by the first increasing-`epsilon` transition from
`Delta_n <= 0.10` to `Delta_n > 0.10`. Equality remains on the safe side. The
reported `epsilon_certified_n` is the verified safe lower bracket endpoint,
not the midpoint estimate. If first loss is not reached, the certified value
is only the right-censored last verified scan point; it is not an estimate of
the unknown crossing.

The original step-2 run used the general 6x6 determinant sign scan as its
Timoshenko production spectrum. That scan can miss two close simple roots of
different families when both fall inside one `Lambda=0.01` interval and leave
the same endpoint sign. Therefore its prefix-5 through prefix-10 thresholds
are superseded. The original files and cache remain under
`legacy_pre_factorized_root_fix/` and
`cache_legacy_pre_factorized_root_fix/` for provenance only.

The corrected source `factorized_straight_spectrum_v2` uses an exact
special-limit factorization of the unchanged matrix. With unknown ordering
`(A1,B1,P1,A2,B2,P2)`, bending uses rows/columns
`[0,2,3,4]/[0,1,3,4]` and axial motion uses `[1,5]/[2,5]`. The EB source is
the exact axial plus fixed-fixed EB bending union; the Timoshenko source is the
same exact axial family plus roots of the extracted 4x4 bending block found by
sign changes and complementary SVD minima. Cross-family duplicates are kept;
only numerical duplicates inside one family are merged. The raw general 6x6
scan remains an independent completeness audit and the production mechanism
outside this special straight homogeneous limit.

The corrected full run used a `0.0005` coarse step on `0.005..0.060`, screened
the midpoint of every coarse interval, refined all detected
order/cluster/threshold events, kept the first 12 sorted roots for quality
control, and independently force-recomputed every critical bracket. The
resulting straight-system thresholds are:

| prefix n | status | epsilon_certified_n | epsilon_star_estimate | trigger |
|---:|---|---:|---:|---|
| 1 | not reached through 0.060 | 0.060000000 | -- | -- |
| 2 | resolved | 0.049140625 | 0.049141113 | sorted 2, EB bending |
| 3 | resolved | 0.037009766 | 0.037010254 | sorted 3, EB bending |
| 4 | resolved | 0.029705078 | 0.029705566 | sorted 4, EB bending |
| 5 | resolved | 0.029705078 | 0.029705566 | sorted 4, EB bending |
| 6 | resolved | 0.024823242 | 0.024823730 | sorted 6, EB bending |
| 7 | resolved | 0.021326172 | 0.021326660 | sorted 7, EB bending |
| 8 | resolved | 0.018695312 | 0.018695801 | sorted 8, EB bending |
| 9 | resolved | 0.016643555 | 0.016644043 | sorted 9, EB bending |
| 10 | resolved | 0.016643555 | 0.016644043 | sorted 9, EB bending |

Prefixes 4--5 and 9--10 form simultaneous transition groups within the
requested epsilon tolerance. All nine finite first-loss estimates lie in the primary
range `0.010..0.050`; prefix 1 is instead certified only through the upper
buffer endpoint. Independent force-recompute verification passed for all nine
critical brackets. Later behavior remains separate from first loss: the scan
recorded four unsafe-to-safe re-entries, four subsequent safe-to-unsafe returns,
53 family reorder events, and 721 evaluated points with at least one individual
late pass. None of these increases the conservative certificate.

All 1004 corrected quality rows are resolved. The factorized-spectrum audit
has 23592/23592 passing rows (11796 EB and 11796 Timoshenko), the nine R1--R3
and `epsilon +/- 1e-6` rows all retain both roots, and all 720 first-12
comparisons at `mu=0,0.3,0.7,0.9` pass with unchanged order, multiplicity, and
`N_true`. The independent raw general-6x6 comparison misses 155 factorized
roots over all primary, verification, and mu-audit scopes (91 EB and 64
Timoshenko); all 155 are independently confirmed by local full-matrix SVD
refinement. All 2754 required axial records pass their block/full-matrix
checks, including 27 Timoshenko axial records absent from the raw sign scan.
These omissions remain visible diagnostics but do not contaminate the
corrected threshold source.

The straight baseline defines

```text
N_certified_0(epsilon) = max { n : epsilon <= epsilon_certified_n }
```

only over the scanned buffer. It is not a global lower envelope over `beta`,
`mu`, or `eta`. Research step 3 remains pending and was not started by this
correction. A future targeted geometry search may probe the recorded
`epsilon_near_n = 0.999 epsilon_star_n` and
`epsilon_buffer_n = 0.99 epsilon_star_n` values while keeping first loss
distinct from re-entry.

## General-Spectrum Completeness Audit Before Step 3

Research step 2.5 was implemented and run on `2026-07-22` through
[`audit_eb_timo_general_spectrum_completeness.py`](../../scripts/analysis/thickness_mismatch/audits/audit_eb_timo_general_spectrum_completeness.py).
It is an offline audit of the unchanged general coupled Euler--Bernoulli and
Timoshenko `6x6` matrices. It does not replace the production root wrappers.
The candidate union combines unshifted and shifted determinant scans,
independent half-step scans, normalized smallest-singular-value and
singular-value-ratio valleys, adaptive local refinement, continuation seeds,
EB seed windows for Timoshenko, and the straight factorized oracle where its
`beta=0`, `eta=0` assumptions apply. A candidate is accepted only after a
row-normalized full-matrix SVD check. Lambda proximity alone does not merge
roots: self-MAC and compatible detection histories are also required, while
exact nullity and coalesced continuation tracks remain separate metadata.

The primary search retained 20 candidate roots and the independent
verification search retained 24, with the first 12 used for the `K=10`
completeness decision. Root 11 is the right spectral-gap guard and root 12 is
a candidate-boundary margin. The two searches use different steps and phases;
verification receives no primary candidate objects or primary seeds. The
versioned JSON cache records geometry, model, search bounds, steps, phases,
refinement, SVD settings, root tolerances, and
`general_complete_svd_v1`. Stale algorithm versions are rejected explicitly.

The finite full run produced the following evidence:

- the general `6x6` audit resolved both EB and Timoshenko for 17 of the 21
  pilot geometries; unresolved rows were B07 Timoshenko, G01 EB, G02 EB, and
  both M02 models;
- the straight oracle comparison passed 431 of 432 first-12 rows; the sole
  mismatch was G02 EB root 12, and 65 roots absent from raw sign-scan prefixes
  were recovered and checked in the full matrix;
- the requested R1--R3, `epsilon +/- 1e-6`, and
  `beta=0,0.1,1,5 deg` stress set plus pilot cluster controls contained 110
  rows. Twenty-six stress rows failed or remained unresolved, including 21
  unresolved model summaries. The minimum retained close-pair gap was
  `9.19871808003e-4`;
- across all 114 model/geometry audit rows, the candidate records contained 67
  close clusters and no asserted algebraic multiplicity. The run performed
  2,290,411 characteristic-matrix evaluations, 2,106,161 full `6x6` SVD
  calls, 8,040 adaptive subdivisions, and 254,040 Brent/bisection evaluations;
  this first full run had 114 cache misses and no cache hits;
- the corrected pilot's explicit auto-spectrum mode used the exact
  factorized source only in the straight special limit and the general audit
  elsewhere. It included 20/21 geometries and excluded M02 with
  `EB_spectrum_unresolved;Timo_spectrum_unresolved`; it never substituted
  legacy roots for an unresolved corrected result;
- all included first-ten roots and all 21 reported `N_true` values matched the
  legacy pilot. False-safe geometry counts were unchanged in all 53 rule
  comparison rows. Retention and conservative-loss summaries changed because
  M02 was excluded, with maximum absolute changes of `0.172414` and
  `1.333333`, respectively; thresholds were not retuned in response to these
  results.

This is a finite numerical audit, not a mathematical root-count proof. Its
readiness decision is `not_ready_for_step3`: the general pilot audit retains
unresolved cases, G02 retains a straight-oracle mismatch at root 12, and the
small-angle stress audit retains unresolved/failing rows. Step 3, any
four-dimensional lower-envelope search, FEM, and article work remain unrun.

## Branch-informed continuation gateway to step 3

Research step 2.5b was completed on `2026-07-22` with
[`audit_eb_timo_branch_continuation_gateway.py`](../../scripts/analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py).
It addresses the unresolved step-2.5 cases without changing either production
characteristic matrix. At `beta=0`, the existing EB and Timoshenko coefficient
orderings are used to extract separate exact axial and bending parent blocks;
the measured off-block maximum over the completed audit is zero. Cross-family
multiplicity is preserved, and the straight `eta=0` parent spectra retain their
factorized-oracle regression.

For `beta>0`, isolated roots are continued through adaptive local windows.
Close roots are continued as left/right null-subspace clusters: a reduced
matrix proposes candidates, but the unchanged full `6x6` matrix determines
stationarity, singularity quality, and the accepted root record. Seeds never
become terminal roots by themselves. Direct `seed_only` acceptance is
prohibited; a new record requires `seed_refined_to_new_root`, while convergence
to an already accepted stationary minimum is recorded separately. A cheap
global guard searches below root 11, and strict global search is used only as
a triggered fallback or in the separate force-verification scope.

The completed gateway resolved `K10_guard_resolved` for all 122 audited
model/geometry records: the R1--R3 base and `epsilon +/- 1e-6` cases at
`beta=0,0.1,1,5 deg`, B07/G01/G02/M02, and both models for all 21 pilot
geometries. The full-12 margin resolved for 103/122 records; this is reported
separately and is not substituted for the `K=10` decision. The branch-informed
pilot included 21/21 geometries and changed no first-ten root, `N_true`, or
first-failure result relative to the comparison pilots. All local-independent,
cluster, root-11 guard, straight-oracle, and requested force-global readiness
checks passed. The resulting decision is `ready_for_targeted_step3`.

Because that gate passed, the future-only compact manifest
[`eb_epsilon_lower_envelope_step3_cases.csv`](../../scripts/analysis/thickness_mismatch/audits/data/eb_epsilon_lower_envelope_step3_cases.csv)
was written. Its 28 unique geometries use only full-precision
`epsilon_near_n` and `epsilon_buffer_n` values from the corrected straight
baseline, deduplicate the prefix 4/5 and 9/10 thresholds, and select baseline,
small-angle, 45/90-degree, high-`mu`, signed-`eta`, and mixed probes. No root
calculation or lower-envelope search from that manifest was run in step 2.5b.

## Step 3A: targeted lower-envelope screening

Step 3A executes the fixed 28-row manifest produced by the branch gateway; it
does not expand it into a tensor grid or optimize over geometry. For sorted
mode `k`, prefix `n`, and the fixed threshold `tau=0.10`, the audited
quantities are

\[
\delta_{f,k}=\frac{|\Lambda_{\mathrm{EB},k}^{2}-
\Lambda_{\mathrm{Timo},k}^{2}|}{\Lambda_{\mathrm{Timo},k}^{2}},\qquad
\Delta_n=\max_{1\leq k\leq n}\delta_{f,k},
\]

\[
V_n=\Delta_n-\tau,\qquad M_n=\tau-\Delta_n=-V_n.
\]

Equality is safe. `N_true` is the continuous prefix ending before the first
failed mode; a later individual pass cannot repair it. The groups for
prefixes 4/5 and 9/10 generate distinct screening rows even though each pair
shares one corrected straight-baseline threshold.

Every manifest epsilon is checked at full CSV precision against the corrected
`factorized_straight_spectrum_v2` safe-threshold products before root
calculation. The scientific inclusion gate requires branch-informed EB and
Timoshenko roots 1--10 and resolved root 11; root 12 is retained only as a
diagnostic margin. The 14 straight controls are compared independently with
the corrected factorized oracle. Legacy roots never fill an unresolved case.

`provisional_counterexample` and near-boundary/cluster/quality triggers do not
determine the conclusion from one run. They cause a force-recomputed second
continuation with an independently fixed beta path, separate cache, unchanged
full-`6x6` residual checks, root-11 guard, and automatic strict fallback when
local continuation needs it. A `confirmed_counterexample` requires both
resolved runs, root/cluster agreement, a violation above the fixed numerical
tolerance, and failure of the target prefix or full corrected baseline
certificate. Sign instability is reported as
`numerically_indeterminate_near_threshold` rather than rounded into a pass or
failure.

The runner records one of `counterexample_found`,
`no_counterexample_in_28_case_screen`,
`inconclusive_due_to_unresolved_cases`, or
`inconclusive_due_to_numerical_boundary`. It also writes an unexecuted Step-3B
proposal pairing near/buffer epsilon values at the same selected geometry.
Adversarial near/buffer rows in Step 3A itself are usually different
geometries and are not interpreted as thickness trajectories. Even a complete
28-case result is finite numerical evidence, not a continuous-domain proof of
a global lower envelope.

The completed run resolved K10/root 11 for all 28 geometries and all 56
primary model spectra; 40/56 primary spectra also resolved the optional root
12. All 18 baseline-control/prefix rows passed the corrected factorized-oracle
pipeline. Automatic triggers force-recomputed 27 geometries with the separate
verification configuration, and every verification row passed its K10,
root-match, and cluster checks. Two candidates were confirmed by both runs:
`S3_12` at prefix 5 (`beta=90 deg`, `mu=0.7`, `eta=0`, verified
`V_5=1.73946990918e-2`) and `S3_14` at prefix 6 (`beta=45 deg`, `mu=0.5`,
`eta=-0.1`, verified `V_6=5.09348548033e-4`, triggered by sorted mode 5).
There were no unresolved or numerically indeterminate geometries. The Step-3A
decision is therefore `counterexample_found`. The 18-row paired Step-3B
proposal covers prefixes 2--10 and remains unexecuted; no local epsilon
refinement was performed.

## Engineering Objective

For each complete system parameter point, first solve the Euler--Bernoulli
model for the first `K = 10` sorted eigenfrequencies. Use those frequencies,
their EB mode shapes, and inexpensive EB-only indicators to estimate the
largest prefix length

```text
N_hat in {0, 1, ..., 10}
```

for which the first `N_hat` frequencies may be retained from the EB model.
Frequencies with sorted indices `N_hat + 1, ..., 10` are then to be computed
with the Timoshenko model. Thus `N_hat = 0` selects Timoshenko for all ten
frequencies, while `N_hat = 10` retains the complete requested EB spectrum.

The immediate objective is an **EB-based** certificate: an EB solve is already
available and supplies modal information that should be exploited. A predictor
based only on geometric parameters, without solving the EB eigenproblem, is a
possible later stage and is not the primary task now.

## Formal Target

Let `g` denote one complete geometry and material parameter point. At fixed
`g`, let `Lambda_EB,k` and `Lambda_Timo,k` be the independently sorted
dimensionless frequency parameters of the two one-dimensional theories. The
relative dimensional-frequency discrepancy is

\[
\delta_{f,k}(g) =
\frac{\left|\Lambda_{\mathrm{EB},k}^{2}(g)
- \Lambda_{\mathrm{Timo},k}^{2}(g)\right|}
{\Lambda_{\mathrm{Timo},k}^{2}(g)}.
\]

The square is required because, at one fixed physical parameter point, the
dimensional frequency is proportional to `Lambda^2` with the same dimensional
factor for EB and Timoshenko.

For the 10% acceptance threshold, define the true safe sorted prefix as

\[
N_{\mathrm{true}}(g) =
\max\left\{n \in \{0,\ldots,10\}\;:\;
\delta_{f,k}(g) \leq 0.10
\text{ for every } k=1,\ldots,n\right\}.
\]

The empty prefix `n = 0` is admissible, so `N_true = 0` when the first sorted
frequency fails. A later frequency that happens to satisfy the threshold does
not repair an earlier failure and does not extend the safe prefix.

`N_hat(g)` is the output of the conservative EB-only criterion. The dangerous
error, or **false-safe**, is

\[
N_{\mathrm{hat}}(g) > N_{\mathrm{true}}(g).
\]

The conservative error `N_hat(g) < N_true(g)` is acceptable, but its cost must
be measured as the number of lost safe EB modes,

\[
L_{\mathrm{conservative}}(g) =
\max\left(0, N_{\mathrm{true}}(g)-N_{\mathrm{hat}}(g)\right).
\]

The practical target is defined by the first ten **sorted frequencies**.
Homologous-mode MAC, one-to-one matching, close-cluster subspace checks, and
branch-character diagnostics remain essential controls for root quality,
matching ambiguity, and physical interpretation. They do not replace or
redefine the sorted-spectrum target.

## Validation Reference

The Timoshenko theory is the primary one-dimensional reference used to build
and test the certificate. The resulting quantity is the discrepancy between
Euler--Bernoulli and Timoshenko, or an EB applicability boundary relative to
the Timoshenko model; it is not an exact error relative to a real three-
dimensional body.

Existing 3D FEM diagnostics provide encouraging support for the Timoshenko
correction over the thickness range investigated so far. In particular, the
single-rod benchmark supports the correction for classified bending doublets,
and the non-tuned rigid-end-face coupled-rod diagnostics support the
comparative trend that Timoshenko becomes more favorable as thickness
increases. These results remain diagnostic: coupled-mode matching, mesh
convergence, joint idealization, and article-grade assignment limitations are
still open. No stronger validation claim is permitted here than the one stated
in the current [FEM validation status](fem_validation_status.md).

## Current Evidence and Parameter Slices

The current evidence comes from three complementary diagnostic workflows:

- [`audit_eb_validity_vs_timoshenko_stage1.py`](../../scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py)
  provides the thickness scan at `beta = 45 deg`, `eta = 0`, with varying
  `epsilon` and `mu`, while keeping sorted-spectrum and physical-branch
  comparisons separate.
- [`analyze_universal_eb_validity_parameters_stage1.py`](../../scripts/analysis/thickness_mismatch/audits/analyze_universal_eb_validity_parameters_stage1.py)
  post-processes the first scan and compares existing local-thickness,
  modal-thickness, and EB-mode indicators without recomputing roots.
- [`audit_eb_validity_fixed_epsilon_geometry_scan.py`](../../scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py)
  provides the complementary slice at `epsilon = 0.02`, with varying `beta`,
  `mu`, and `eta`, and includes individual-mode and close-cluster matching
  diagnostics.

Together these workflows cover two complementary slices of the parameter
space:

- `beta = 45 deg`, `eta = 0`, varying `epsilon` and `mu`;
- `epsilon = 0.02`, varying `beta`, `mu`, and `eta`.

They provide evidence about candidate indicators and transfer behavior, but
they do not prove that any criterion is universal. In addition, the existing
applicability datasets report eight modes; a `K = 10` source-data extension is
required before the proposed ten-frequency certificate can be calibrated.

## Candidate EB-Only Indicators

The first certification stage must test the existing EB-only quantities before
introducing a new arbitrary scalar predictor:

- `Pi_shear_EB`;
- `Pi_rotary_EB`;
- `Pi_EB = Pi_shear_EB + Pi_rotary_EB`;
- `chi_max_EB`;
- `chi_eff_EB`;
- `Theta_max_EB`;
- `epsilon_max`;
- EB axial and bending energy fractions.

An adjacent EB sorted spectral gap is to be added as an inexpensive protective
feature. Small gaps within a proposed prefix, or at its upper boundary, can
trigger a conservative fallback because close roots increase the risk of
root-order sensitivity, mode mixing, and unreliable individual-mode
interpretation. The gap must remain an explicit guard with a documented
normalization; it must not be hidden inside an unexplained composite scalar.

For fixed material constants, `Theta_max_EB` is a monotonic transformation of
`chi_max_EB`. Agreement between these two quantities therefore is not
independent physical confirmation and must not be counted as such.

## Immediate Next Computational Step

The recommended next stage is the following.

1. Extend both source applicability workflows to `K = 10` reported sorted
   frequencies. Retain a sufficient margin of candidate roots and preserve
   candidate-boundary, root-recovery, close-cluster, and matching diagnostics.

2. Use the separate CSV-only postprocessor at
   `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py`.
   The stable entry point is now implemented and reads only existing source
   CSV files; it does not recompute roots.

3. For every complete geometry, compute `N_true` and candidate `N_hat` values
   for nested rule families:

   - **Rule A:** `Pi_EB <= threshold`.
   - **Rule A-gap:** Rule A plus the minimum adjacent EB sorted spectral-gap
     guard.
   - **Rule B:** `Pi_shear_EB <= threshold_s` and
     `Pi_rotary_EB <= threshold_r`.
   - **Rule C:** Rule B plus a minimum adjacent EB sorted spectral-gap guard.
   - **Rule D:** Rule C plus EB-only modal-character guards or subgroup-specific
     thresholds based on EB axial and bending energy fractions.

   Each rule must return a prefix: isolated passing modes after the first
   rejected mode must not increase `N_hat`.

4. Calibrate each rule by maximizing the total number of accepted EB
   frequencies subject to `false_safe_geometry_count = 0`, that is, **zero
   observed false-safe** on the calibration set. Threshold selection, tie
   breaking, and any subgroup boundaries must be reproducible and reported.

5. Validate on held-out **complete geometries**, not on randomly split mode
   rows. The same geometry must never contribute modes to both calibration and
   validation. Required transfer checks hold out beta/eta groups, epsilon
   groups or ranges, and mu groups or ranges.

6. Introduce a new physically motivated composite criterion only if the
   existing indicators produce systematic false-safe cases or unacceptably low
   coverage, and only after the failure structure identifies the missing
   physical or spectral dependence.

## Required Quality Metrics

Every calibration and held-out validation summary must report at least:

- `false_safe_geometry_count`: number of geometries with
  `N_hat > N_true`;
- `false_safe_frequency_count`: total
  `sum(max(N_hat - N_true, 0))` over geometries;
- `mean N_true`;
- `mean N_hat`;
- `exact_N_match_rate`: fraction of geometries with `N_hat = N_true`;
- `mean conservative loss N_true - N_hat`, evaluated as the mean of
  `max(N_true - N_hat, 0)`;
- `maximum conservative loss`;
- `fraction of usable EB frequencies retained`, evaluated as
  `sum(min(N_hat, N_true)) / sum(N_true)` when the denominator is nonzero;
- `cluster-guard trigger count`;
- `out-of-calibration-domain count`.

Results must distinguish **zero observed false-safe on finite calibration or
validation data** from a mathematical guarantee over a continuous parameter
domain. The former is an empirical result conditional on the sampled
geometries; it must never be presented as the latter. Out-of-domain points must
be flagged rather than silently certified by extrapolation.

## Computational-Cost Policy

The primary efficiency evidence is operation counts or counts of algorithmic
primitives, not wall-clock time. Future source workflows and certification
audits should expose, aggregate, and preserve counters for:

- EB characteristic-matrix evaluations;
- EB determinant evaluations;
- Timoshenko characteristic-matrix evaluations;
- Timoshenko determinant evaluations;
- `6x6` SVD calls;
- root-bracketing intervals;
- Brent/function evaluations;
- EB mode-shape reconstructions;
- shape sample evaluations;
- quadrature-point evaluations;
- local Timoshenko root attempts;
- full Timoshenko fallbacks.

Wall-clock time may remain an auxiliary diagnostic for a stated machine,
software environment, cache state, and worker count. It is not the principal
evidence for computational savings because it mixes implementation details,
hardware, parallelism, and cache effects.

## Non-Goals

This research direction does not currently require:

- changing analytic formulas or determinants;
- changing root solvers;
- changing the Timoshenko shear coefficient;
- changing FEM geometry, constraints, or discretization;
- running a full four-dimensional tensor grid;
- editing article workspaces;
- introducing a geometry-only predictor before the EB-based certificate has
  been evaluated.

The plan also does not promote existing diagnostic outputs into article
results. Any later promotion remains subject to the repository's separate
diagnostic-to-article workflow.
