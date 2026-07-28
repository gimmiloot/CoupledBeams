# EB-only Safe-Prefix Study: Final Status and Stage Closure

## 1. Status

**Research status:** `closed after a negative engineering cost result`

Rule S survived the checked finite-sample safety tests but did not reduce the
cost of the reliable `K=10` computation relative to the direct Timoshenko
workflow.

This closes the present empirical-selector stage. It does **not** refute Rule S
mathematically, prove that every EB-only indicator must fail, or rule out every
possible hybrid algorithm.

## 2. Initial research objective

For a geometry and parameter point \(g\), the first \(K=10\) Euler–Bernoulli
(EB) frequencies were assumed to be available. The intended selector was to
use only those relatively inexpensive EB data to estimate the largest
continuous sorted-frequency prefix

\[
\widehat N(g)\leq 10
\]

that could be retained from EB. The suffix
\(\widehat N(g)+1,\ldots,10\) would then be computed with the Timoshenko model.

The reference safe prefix was

\[
N_{\mathrm{true}}(g)=
\max\left\{
n\leq 10:
\delta_{f,k}(g)\leq 0.10
\quad\text{for every }k=1,\ldots,n
\right\},
\]

where

\[
\delta_{f,k}(g)=
\frac{
\left|\Lambda_{\mathrm{EB},k}^{2}(g)-
\Lambda_{\mathrm{Timo},k}^{2}(g)\right|
}{
\Lambda_{\mathrm{Timo},k}^{2}(g)
}.
\]

The square is required because the dimensional cyclic frequency at one fixed
physical point is proportional to \(\Lambda^2\) with the same dimensional
factor for the two theories. The safety condition was

\[
\widehat N(g)\leq N_{\mathrm{true}}(g).
\]

A **false-safe** prediction, \(\widehat N>N_{\mathrm{true}}\), leaves at least
one frequency outside the 10% tolerance in the accepted EB prefix. It is more
serious than a conservative prediction, \(\widehat N<N_{\mathrm{true}}\),
which only sends additional safe modes to the more expensive reference
workflow. Retention was therefore optimized only after zero observed
false-safe had been imposed. The final objective was computational rather than
classification-only: a safe selector was useful as an engineering algorithm
only if EB selection plus the Timoshenko suffix cost less than a reliable
direct Timoshenko `K=10` calculation.

## 3. Parameters and investigated domain

The study used the mass-preserving two-arm geometry described in the
[thickness-mismatch model note](README.md):

- \(\epsilon_0=r_0/(2l)\) is the base thickness parameter formed with the
  equal-radius reference \(r_0\) and base length
  \(l=(l_1+l_2)/2\);
- \(\beta\) is the coupling angle;
- \(\mu=(l_2-l_1)/(l_1+l_2)\) controls the arm lengths
  \(l_1=l(1-\mu)\), \(l_2=l(1+\mu)\);
- \(\eta=(r_2-r_1)/(r_1+r_2)\) controls the mass-preserving radius mismatch.

The target was always the first ten **sorted frequencies**, not ten descendant
branches, with a 10% dimensional-frequency tolerance between Euler–Bernoulli
and Timoshenko results. The accepted set had to be a continuous prefix: a later
passing mode could not repair an earlier failure.

Close roots and clusters made ordinary sign scans and index-by-index
interpretation unsafe. Root multiplicity had to be preserved, and root 11 was
used as the right-hand completeness guard for the first ten roots. The
verified EB equations and notation remain in
[the theory equations](../theory/equations.tex) and
[the assumptions note](../theory/assumptions.md). The diagnostic Timoshenko
implementation is documented by
[`variable_length_timoshenko.py`](../../scripts/lib/variable_length_timoshenko.py),
while the EB-only predictors used here originate in
[`audit_eb_validity_fixed_epsilon_geometry_scan.py`](../../scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py).
This closure note does not revise any formula, determinant, unknown ordering,
solver, or model assumption.

## 4. Research chronology

### 4.1 Geometry-only screening

The first question was whether geometry could determine a conservative prefix
before any EB modal analysis. The 21-geometry pilot compared:

- \(\epsilon_0\), the base thickness;
- \(\epsilon_{\max}\), the maximum local thickness parameter;
- a straight symmetric baseline at \((\beta,\mu,\eta)=(0,0,0)\);
- the hypothesis that this baseline provided a conservative lower envelope
  over non-straight and asymmetric geometries.

\(\epsilon_{\max}\) was a weak predictor. Matched-
\(\epsilon_{\max}\) groups had materially different true prefixes, and the
quantity depended on the artificial internal-arm split used to describe even a
straight homogeneous geometry. It was therefore not accepted as a physically
reliable certificate.

\(\epsilon_0\) transferred better and remained useful as a screening
quantity, but the pilot did not make it a certificate. The straight-baseline
lower-envelope hypothesis still required an adversarial geometry check.

### 4.2 Repair of spectral completeness

The original general determinant sign scan could miss two close simple roots
inside one scan interval, particularly when axial and bending families nearly
coincided. That defect affected the sorted spectrum and could therefore
corrupt both \(N_{\mathrm{true}}\) and any learned threshold.

The reference-data path was repaired and audited through:

- exact axial/bending factorisation in the straight symmetric special case;
- preservation of cross-family multiplicity;
- branch-informed isolated-root and cluster continuation away from the
  special limit;
- root 11 as the `K=10` completeness guard;
- independent strict verification and full-matrix SVD quality checks;
- the explicit `K10_guard_resolved` status, kept separate from the optional
  full-12 margin.

The branch-informed gateway ultimately resolved the `K=10` guard for all
122 audited model/geometry records and for both models in all 21 pilot
geometries. This infrastructure was necessary to build a credible reference
dataset. It is not part of the final one-threshold Rule-S selector.

### 4.3 Step 3A and rejection of the epsilon lower envelope

Step 3A tested the fixed 28-case directed manifest against the corrected
straight baseline. It found two independently verified counterexamples.

For `S3_12`:

```text
epsilon_0 = 0.029408510742187498
beta = 90 degrees
mu = 0.7
eta = 0
straight-baseline prediction = 5
N_true = 4
delta_f,5 = 0.11739469909177262
```

For `S3_14`:

```text
epsilon_0 = 0.024798906738281248
beta = 45 degrees
mu = 0.5
eta = -0.1
baseline-minus-one prediction = 5
N_true = 4
delta_f,5 = 0.10050934854803291
```

Both cases had quality-approved roots 1–10 and a resolved root-11 guard in
the primary and independently triggered verification paths. Neither result
was unresolved or numerically indeterminate. The conclusion was that
\(\epsilon_0\) could remain a screening quantity but could not serve as the
proposed geometry-only certificate. Step 3B was not needed to answer that
main hypothesis and was not executed.

### 4.4 Historical Rules A–D

The initial EB-only families were:

- **Rule A:** one upper threshold for \(\Pi_{\mathrm{EB}}\);
- **Rule B:** separate upper thresholds for
  \(\Pi_{\mathrm{shear}}\) and \(\Pi_{\mathrm{rotary}}\);
- **Rule A-gap:** Rule A plus an adjacent EB sorted-gap condition;
- **Rule C:** Rule B plus the gap condition;
- **Rule D:** Rule C plus an EB modal-character guard.

The historical A-gap and C results did not isolate the effect of adding a
gap guard. Their base predictor thresholds and gap threshold were optimized
jointly, so the resulting rules were not strict nested ablations of frozen
Rule A or Rule B. Their retention changes could not be interpreted as a clean
gap effect. The A-gap, C, and D paths were therefore closed rather than
recalibrated. No new guard, machine-learning model, subgroup threshold, or
composite indicator was introduced.

## 5. Exact Rule A/B/S experiment and data lock

The deciding Phase-I experiment used four disjoint partitions:

| Partition | Geometries | Role |
|---|---:|---|
| development | 21 | exact threshold selection |
| Step-3A baseline controls | 14 | straight-baseline controls, never used for threshold selection |
| nonbaseline directed validation | 12 | directed validation, never used for threshold selection |
| locked adversarial holdout | 2 | `S3_12` and `S3_14`, opened only after freezing thresholds |

There were 49 unique, quality-approved geometries. Full-precision geometry
deduplication used independent \(10^{-12}\) tolerances on
\((\epsilon_0,\beta,\mu,\eta)\); no included geometry leaked between
partitions. Thresholds were selected only on the 21-case development set and
were not changed after directed validation or after opening the S3 holdout.
Phase I read saved spectra and computed no new EB or Timoshenko roots.

## 6. Rule A benchmark

For geometry \(g\) and prefix \(n\), define

\[
P_{g,n}=\max_{1\leq k\leq n}\Pi_{\mathrm{EB},g,k}.
\]

Rule A returns

\[
N_A(g)=\max\{n:P_{g,n}\leq T_A\}.
\]

The exact search evaluated 159 observed-development candidates, including the
reject-all endpoint, and selected

\[
T_A=0.20310844707256814.
\]

It retained 147 of 155 safe development frequencies,
\(147/155=0.948387\ldots\), and produced zero observed false-safe geometries
on all checked partitions. Its fixed status remains
`rule_A_no_false_safe_on_checked_sets_benchmark_only`. This is a useful
benchmark, not a production or continuous-domain certificate.

## 7. Rule B and its degeneration to Rule S

Define the prefix extrema

\[
S_{g,n}=\max_{1\leq k\leq n}\Pi_{\mathrm{shear},g,k},
\qquad
R_{g,n}=\max_{1\leq k\leq n}\Pi_{\mathrm{rotary},g,k}.
\]

Rule B returns

\[
N_B(g)=\max\left\{n:
S_{g,n}\leq T_s,\quad R_{g,n}\leq T_r
\right\}.
\]

The exact development search used 159 shear candidates, 158 rotary
candidates, and all \(159\times158=25{,}122\) Cartesian pairs. It found 32
equally optimal pairs and one Pareto-nondominated optimum. The corrected frozen
representative was

\[
T_s=0.16762413001084248,
\qquad
T_r=0.046719283392029604,
\]

selected by the declared minimum-\(T_s\), then minimum-\(T_r\) tie-break. Its
development objective was

\[
153/155=0.9870967741935484.
\]

All 32 equally optimal pairs were applied to five partition summaries. The
audit contains 160 partition-level sensitivity rows and 896 per-geometry
comparisons. Every prediction vector and every safety signature was identical
to the representative result, so the observed optimum did not depend on which
equally optimal pair was selected.

The independent shear-only search defined

\[
N_S(g)=\max\left\{n:S_{g,n}\leq T_s\right\}.
\]

It independently selected the same \(T_s\), retained the same 153 of 155
development frequencies, and produced \(N_S=N_B\) for all 49 geometries.
The rotary decision-binding count was zero.

> On all checked data, the optimal rectangular Rule B degenerates to the
> one-threshold shear-only Rule S; the rotary component is nonbinding.

This is an empirical statement for the checked finite datasets, not a theorem
over the continuous parameter domain.

## 8. Frozen validation results

| Partition | Geometries | False-safe | Rule-S retention | Mean conservative loss |
|---|---:|---:|---:|---:|
| development | 21 | 0 | 0.9870967741935484 | 0.09523809523809523 |
| baseline controls | 14 | 0 | 0.75 | 1.5714285714285714 |
| nonbaseline directed validation | 12 | 0 | 0.9324324324324325 | 0.4166666666666667 |
| locked holdout | 2 | 0 | 1.0 | 0.0 |

On the locked holdout:

```text
S3_12: N_true = N_S = 4
S3_14: N_true = N_S = 4
```

In both cases mode 5 was the first rejected mode, and the trigger was the
shear component. These are finite-sample validation results, not a proof of
global safety.

## 9. Threshold stability and margins

The leave-one-development-geometry-out audit recalibrated 21 train-20 folds.
No held-out fold produced a false-safe result. The selected thresholds ranged
from

```text
0.16666777788825218
to
0.16762413001084248
```

and the frozen full-development threshold was not changed.

The exact Rule-S threshold-separation audit was:

| Partition | Nearest accepted prefix; accepted margin | Nearest unsafe-above prefix; unsafe margin |
|---|---|---|
| development | `G05, n=8`; `0` | `G03, n=9`; `0.0019349017842357208` |
| directed validation | `S3_28, n=9`; `0.0040508828965558075` | `S3_22, n=9`; `0.0016242567340815917` |
| locked S3 holdout | `S3_12, n=4`; `0.016357859510945455` | `S3_14, n=5`; `0.00948000082193301` |

The zero development accepted margin means that `G05, n=8` lies exactly on
the frozen threshold; the finite development result is therefore stable under
the checked LOO recalibrations but does not contain a positive accepted-side
development buffer.

## 10. Finite-sample limitation of the monotone threshold class

The development prefix plane contained 155 safe prefix points and 55 unsafe
prefix points. Two safe points could not be accepted by the monotone
rectangular threshold class without also accepting an unsafe point:

```text
unsafe G03, n=9 dominates safe G03, n=8
unsafe G03, n=9 dominates safe G05, n=10
```

Here “dominates” means that the unsafe point is componentwise no larger in
both \((S,R)\), so any upper-threshold rectangle that includes the safe point
also includes that unsafe witness. Consequently, the maximum attainable
zero-false-safe development retention for this finite Rule-B class was

\[
153/155.
\]

This is a **finite-sample limitation of the monotone threshold class**. It is
not an impossibility result over the continuous parameter domain or over all
possible EB-only indicators.

## 11. Frozen five-case cost benchmark

Finite-sample safety did not by itself answer the original engineering
question. Phase II therefore compared two locked workflows.

### Direct workflow

The direct workflow computed the first ten Timoshenko frequencies with the
existing reliable `K=10` path.

### Hybrid workflow

The hybrid workflow used:

1. the existing EB roots;
2. sequential EB mode reconstruction;
3. frozen Rule S, stopping at the first rejected mode;
4. local Timoshenko solves only for the predicted suffix;
5. full `K=10` fallback when a local solve was incomplete, ambiguous, failed
   its reliability checks, or disagreed with post-solve verification.

Saved Timoshenko roots were withheld from online bracketing and were used only
after solving for verification.

The locked cases and their manifest roles were:

| Case | Manifest source | Partition | \(N_S\) | Suffix | Benchmark role |
|---|---|---|---:|---:|---|
| `B01` | `branch_informed_pilot` | development | 10 | 0 | all-EB/no-suffix control |
| `G04` | `branch_informed_pilot` | development | 9 | 1 | short-suffix development case |
| `S3_06` | `step3a` | directed validation | 4 | 6 | nonbaseline directed-validation suffix case |
| `S3_12` | `step3a` | locked adversarial holdout | 4 | 6 | first adversarial counterexample |
| `S3_14` | `step3a` | locked adversarial holdout | 4 | 6 | second adversarial counterexample |

The legacy manifest `purpose` value is
`future_direct_Timoshenko_K10_cost_break_even_benchmark` for every row; the
source, partition, and suffix columns above distinguish the case roles.

## 12. Canonical cost result

Only the final canonical benchmark outputs are used for the cost conclusion.

Canonical direct operation counts were:

```text
5 Timoshenko spectrum calls
56,611 Timoshenko determinant evaluations
```

Canonical hybrid operation counts were:

```text
35 EB reconstructions
35 EB 6x6 SVD calls
210 quadrature calls
19 Timoshenko calls
52,001 Timoshenko determinant evaluations
15 Timoshenko SVD calls
4 full K=10 fallbacks
```

All four geometries with a nonempty suffix required full `K=10` fallback.
`S3_06` failed an online local reliability check; `G04`, `S3_12`, and
`S3_14` triggered fallback after post-solve sorted-mode verification
mismatches. The local suffix path therefore failed to avoid the reliable full
spectrum calculation in every case where it was needed. The reduction of
4,610 Timoshenko determinant evaluations did not remove the EB selector,
local-solve, verification, and fallback overhead.

The median runtimes over all final case/repeat observations were

\[
t_{\mathrm{direct}}=1.3191177000007883\ \mathrm{s},
\qquad
t_{\mathrm{hybrid}}=1.4112544000017806\ \mathrm{s}.
\]

Runtime was supporting evidence; no arbitrary scalar weighting was imposed on
determinant evaluations, reconstructions, SVDs, quadratures, and time. The
final cost status is

```text
rule_S_cost_not_beneficial
```

## 13. Diagnostic runs and provenance

Before the final verification-mismatch fallback was introduced, a preliminary
diagnostic benchmark was run to inspect the local suffix behavior. Those
diagnostic outputs were overwritten by the final benchmark. They contributed
63 additional instrumented calls and do not enter the canonical operation or
runtime comparison. One still earlier command was stopped by a one-second
timeout before it could become a scientific benchmark result.

These diagnostic calls are execution provenance, not the cost of the
production algorithm. All scientific cost tables and the status above use
only the final benchmark repetitions recorded in the canonical Phase-II CSVs.

## 14. Scientific and engineering interpretation

### Positive scientific result

Rule S is simple and physically interpretable. It reached the finite-sample
153/155 ceiling of the checked monotone rectangular Rule-B class, produced no
observed false-safe result on development, baseline controls, directed
validation, or the locked S3 holdout, and was stable in the 21-fold LOO audit.
It also established that the rotary component supplied no additional decision
information in the checked Rule-B datasets. This is useful finite-sample
evidence and a useful diagnostic of EB validity.

### Negative engineering result

For the current production Timoshenko solver, quality contract, local suffix
strategy, \(K=10\), and the five locked benchmark cases, the Rule-S hybrid was
not cheaper than direct Timoshenko `K=10`.

> The shear-only Rule S is a useful finite-sample diagnostic of EB validity,
> but it is not a cost-effective engineering selector for the reliable
> first-ten-frequency workflow implemented in this repository.

The negative cost result does not erase the positive finite-sample safety and
separability findings; it answers the separate engineering question that
motivated the selector.

## 15. What is closed

The following directions are closed within the present research branch:

- retuning Rule A;
- retuning Rule B or Rule S;
- Rule A-gap;
- Rule C;
- Rule D;
- geometry caps;
- promotion of \(\epsilon_0\) to a certificate;
- Step 3B;
- additional S3 screens whose sole purpose is to rescue the selector;
- machine learning;
- a new fitted empirical composite predictor;
- post-hoc guards;
- a large cost grid under the unchanged algorithm.

The reason is cumulative: the best simple candidate was identified, its
finite-sample safety was checked, the limitation of its monotone threshold
class was exposed, and its engineering cost was measured. The original
cost-saving objective was not achieved.

## 16. What is not claimed

This stage does not claim that:

- Rule S is mathematically safe over the entire continuous parameter domain;
- all possible EB-only indicators are impossible;
- every possible local suffix solver is always uneconomic;
- no hybrid algorithm can be constructed;
- Timoshenko must always be faster for every (K), parameter domain, solver,
  or implementation.

The correct scope is:

> Negative engineering conclusion for the current `K=10` workflow, solver
> implementation, quality requirements, and checked benchmark cases.

## 17. Conditions for reopening

This branch should be reopened only after a material change in the problem,
for example:

- a new suffix solver can reliably find only the required high roots without
  a full `K=10` fallback;
- the production Timoshenko `K=10` workflow becomes substantially more
  expensive;
- \(K\) increases enough to change the relative suffix cost;
- a theoretically derived a posteriori majorant becomes available instead of
  another fitted threshold;
- the engineering error criterion or parameter domain changes materially;
- Rule S is needed as a diagnostic rather than as a cost-saving selector.

Merely adding ordinary validation geometries to the same empirical setup is
not, by itself, a sufficient reason to reopen the stage.

## 18. Reproducibility map

| Purpose | Canonical file |
|---|---|
| Research plan and historical audit trail | [`docs/thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md`](eb_safe_spectrum_prefix_research_plan.md) |
| Step-3A lower-envelope decision | [`results/eb_epsilon_lower_envelope_step3a/eb_epsilon_lower_envelope_step3a_report.md`](../../results/eb_epsilon_lower_envelope_step3a/eb_epsilon_lower_envelope_step3a_report.md) |
| Exact A/B/S postprocessor | [`scripts/analysis/thickness_mismatch/postprocess/analyze_eb_rule_ab_exact_pareto.py`](../../scripts/analysis/thickness_mismatch/postprocess/analyze_eb_rule_ab_exact_pareto.py) |
| Rule-S cost benchmark | [`scripts/analysis/thickness_mismatch/benchmarks/benchmark_rule_s_cost_break_even.py`](../../scripts/analysis/thickness_mismatch/benchmarks/benchmark_rule_s_cost_break_even.py) |
| Phase-I report | [`results/eb_rule_ab_exact_pareto/eb_rule_ab_exact_pareto_report.md`](../../results/eb_rule_ab_exact_pareto/eb_rule_ab_exact_pareto_report.md) |
| Phase-II report | [`results/eb_rule_s_cost_break_even/rule_S_cost_break_even_report.md`](../../results/eb_rule_s_cost_break_even/rule_S_cost_break_even_report.md) |
| Partition audit | [`results/eb_rule_ab_exact_pareto/data_partition_audit.csv`](../../results/eb_rule_ab_exact_pareto/data_partition_audit.csv) |
| Rule-S exact search | [`results/eb_rule_ab_exact_pareto/rule_S_exact_search.csv`](../../results/eb_rule_ab_exact_pareto/rule_S_exact_search.csv) |
| Rule-S predictions | [`results/eb_rule_ab_exact_pareto/rule_S_predictions.csv`](../../results/eb_rule_ab_exact_pareto/rule_S_predictions.csv) |
| Dominance witnesses | [`results/eb_rule_ab_exact_pareto/rule_B_dominance_witnesses.csv`](../../results/eb_rule_ab_exact_pareto/rule_B_dominance_witnesses.csv) |
| Validation summary | [`results/eb_rule_ab_exact_pareto/rule_S_validation_summary.csv`](../../results/eb_rule_ab_exact_pareto/rule_S_validation_summary.csv) |
| Cost summary | [`results/eb_rule_s_cost_break_even/cost_break_even_summary.csv`](../../results/eb_rule_s_cost_break_even/cost_break_even_summary.csv) |
| Phase-I operation counts | [`results/eb_rule_ab_exact_pareto/operation_counts.csv`](../../results/eb_rule_ab_exact_pareto/operation_counts.csv) |
| Phase-II operation counts | [`results/eb_rule_s_cost_break_even/aggregate_operation_counts.csv`](../../results/eb_rule_s_cost_break_even/aggregate_operation_counts.csv) |
| Targeted tests | [`tests/test_eb_rule_ab_exact_pareto.py`](../../tests/test_eb_rule_ab_exact_pareto.py) |
