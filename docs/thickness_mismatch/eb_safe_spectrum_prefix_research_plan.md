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
- decide whether Rules A--D retain enough safe frequencies before introducing
  any new physical composite criterion;
- extend evidence beyond the two current parameter slices only after the
  complete K=10 results have been reviewed;
- add source root-operation counters later if this can be done without changing
  solver behavior.

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
