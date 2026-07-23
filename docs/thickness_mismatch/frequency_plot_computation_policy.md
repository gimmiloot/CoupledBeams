# Frequency-Map Computation Policy

## Scope

This document defines the computation contract for future frequency-versus-
parameter maps in the thickness-mismatch study. It separates ordinary figure
production from research-grade spectral certification. It does not change any
characteristic matrix, determinant, unknown ordering, physical coefficient,
root-solver default, or tolerance, and it does not implement a new solver or
command-line interface.

The governing principle is:

> Spectrum calculation and figure rendering are separate stages.

The three calculation modes are `fast_plot`, `certified_audit`, and
`plot_only`. They have distinct purposes, quality requirements, operation
counts, and cache identities. A plotting entry point must state which mode
produced its data.

This policy complements the
[safe-spectrum-prefix research plan](eb_safe_spectrum_prefix_research_plan.md),
the [project rules](../project_rules.md), and the
[thickness-mismatch script map](../../scripts/analysis/thickness_mismatch/README.md).
The current branch machinery is documented by
[`branch_informed_spectrum_continuation.py`](../../scripts/lib/branch_informed_spectrum_continuation.py)
and the
[branch-continuation gateway](../../scripts/analysis/thickness_mismatch/audits/audit_eb_timo_branch_continuation_gateway.py).

## Mode Summary

| mode | purpose | spectrum work | expected use |
| --- | --- | --- | --- |
| `fast_plot` | Ordinary illustrative frequency maps | Sequential branch-informed continuation with roots 1--10 plus root 11 as a guard | Default for future map and figure scripts |
| `certified_audit` | Independent scientific verification | Expanded roots, global searches, independent configurations, and detailed cluster/multiplicity checks | Counterexamples, reference data, solver validation, and difficult clusters |
| `plot_only` | Figure rendering from saved data | None | Any change limited to visual presentation or output format |

The modes are not interchangeable labels for one cache. Their computational
contracts must remain distinguishable in both metadata and cache keys.

## `fast_plot`

`fast_plot` is the standard mode for ordinary illustrative frequency maps. It
must reuse the existing Euler--Bernoulli and Timoshenko matrices and the
branch-informed continuation machinery. A legacy determinant sign scan must
not be the sole evidence that the requested spectrum is complete.

For each fixed `(model, epsilon_0, mu, eta)` combination, the solver must:

1. solve the `beta=0` parent spectrum once;
2. sort the requested beta grid;
3. continue sequentially through `beta_1, beta_2, ..., beta_N`;
4. use adaptive intermediate steps only between adjacent requested nodes;
5. retain branch and cluster metadata for quality control while writing the
   current sorted positions as the plotting target.

The required structure is:

```text
solve beta=0 parent spectrum once

continue sequentially through:
    beta_1, beta_2, ..., beta_N
```

The following structure is prohibited:

```text
for each beta_target:
    restart continuation from beta=0
```

Strict fallback or recovery may start an explicitly recorded independent
search at a suspect target. Such recovery does not redefine the normal path.

### One-path performance contract

For each case/model combination:

```text
parent_spectrum_solve_count = 1
beta_path_restart_count = 0
```

The only permitted exception to `beta_path_restart_count = 0` is an explicitly
recorded strict fallback or recovery. Continuation cost should be approximately
`O(N_beta * K)`, not `O(N_beta^2 * K)`.

Every run must preserve at least these path operation counts:

```text
beta_target_count
continuation_step_count
adaptive_inserted_step_count
path_restart_count
parent_solve_count
```

If `path_restart_count > 0` outside recorded fallback or recovery, the report
must emit a performance warning.

### Root-count contract

The default fast-map target is:

```text
K_plot = 10
n_guard_roots = 11
```

Roots 1--10 are written as plot values. Root 11 is retained as the right
gap/completeness guard for sorted root 10 and is not plotted. Root 12 is not
computed by default, and `full12_resolved` is not a fast-mode requirement.

The acceptance condition is:

```text
fast_plot requires K10_guard_resolved,
not full12_resolved.
```

### Ordinary isolated branches

A well-separated branch should use:

1. a predictor from preceding beta points;
2. a small local `Lambda` window;
3. local root refinement;
4. one final full-matrix SVD residual check;
5. a root-order check against its neighbours.

A full global SVD scan of the complete `Lambda` axis must not be run for every
ordinary isolated branch.

### Close branches and clusters

Cluster/subspace continuation is required when any of the following occurs:

- a small adjacent gap;
- overlapping local windows;
- a low MAC margin;
- branch reorder;
- a change in cluster dimension.

The implementation must continue the cluster subspace, perform a local dense
SVD scan, preserve every resolved root, and invoke strict fallback only when
local cluster resolution fails. Individual vectors must not be continued as
separate canonical branches inside an unresolved cluster.

### Cheap global completeness guard

Fast mode must schedule a cheap global guard at:

- `beta=0`;
- the final beta endpoint;
- scientifically important anchor angles;
- a configurable periodic spacing, normally about 5 degrees;
- cluster or sorted-order-rearrangement events;
- points following a poor residual or poor MAC result;
- any point where the root count changes;
- the angle used as a scientific counterexample, when applicable.

The guard schedule may insert checks near an event. It must not perform the
complete shifted/half-step global audit at every `beta=0.5 deg` target.

### Strict fallback triggers

Strict fallback is allowed only for:

- a missing expected root;
- an unexplained low-sigma valley;
- an unresolved cluster;
- root-order disagreement;
- a bad full-matrix residual;
- a candidate-boundary warning that can affect roots 1--11;
- disagreement with an independent local refinement;
- explicit scientific-anchor verification.

Every fallback must be recorded in the audit CSV with its target, trigger,
result, and operation counts.

## `certified_audit`

`certified_audit` is intended for:

- confirmation of a potential counterexample;
- production of reference training or validation data;
- validation of a root-solver correction;
- investigation of an especially close spectral cluster;
- final independent scientific verification.

This mode may use roots 1--12 plus an additional candidate margin,
independent search configurations, shifted and half-step scans, a full global
SVD audit, strict force-recompute, and detailed multiplicity diagnostics. Its
cost can be substantially higher than an ordinary map.

`certified_audit` must not be the default for an ordinary map or plotting
script.

> A figure may use data previously generated by `certified_audit`, but
> changing its visual formatting must never rerun the certified audit.

The Step-3A screen and branch gateway are examples of scientific workflows
that need stronger evidence than ordinary illustration. Their detailed
operation counts must not be described as PDF-rendering cost.

## `plot_only`

`plot_only` reads a completed CSV and creates figure files. It must not:

- assemble characteristic matrices;
- search for roots;
- perform SVD calculations;
- modify a spectrum cache.

Changes to line style, line width, labels, font size, axis limits, DPI, or
PDF/PNG/SVG output format must use `plot_only` and should take seconds rather
than rerun spectral calculations.

The required invariant is:

```text
plot_only root calculation count must be zero.
```

CSV round-trip tests should verify that a saved spectrum can reproduce the
figure without importing or invoking the root-calculation path.

## Scientific Anchors

When a map includes an already confirmed counterexample, strict verification
is required at least at:

- `beta=0`;
- the counterexample beta;
- both beta endpoints;
- sorted-order-rearrangement points;
- minimum-gap points;
- points where continuation reports a warning.

For the current confirmed cases, the minimum named anchors are:

- `S3_12`: `beta=0` and `beta=90 deg`;
- `S3_14`: `beta=0` and `beta=45 deg`.

Other grid points may use fast continuation when all quality guards pass.
Anchor certification may reuse previously certified data; visual regeneration
does not require recertification.

## Unresolved-Point Policy

If fast mode cannot resolve roots 1--11, it must not substitute a legacy
spectrum or silently interpolate. It must:

1. write `NaN` for affected plotted values;
2. preserve the exact warning and failed quality status;
3. leave a visible line break through the unresolved interval;
4. recommend a `certified_audit` fallback.

Lines must not connect across an unresolved beta interval without explicit,
recorded scientific approval.

## Sorted Spectrum and Branch Lineage

For maps related to `N_true`, the plotting target is the first ten sorted
frequencies at each parameter point. Branch-informed continuation provides
computational stability and completeness evidence, but branch identity does
not replace the sorted-spectrum target. Physical branches may exchange sorted
positions. Branch, parent-family, cluster, MAC, and reorder metadata must be
retained in CSV audit fields.

This distinction follows the project-wide branch identity and sorted-position
rules in [`docs/project_rules.md`](../project_rules.md).

## Cache Policy

For each case/model combination, fast-mode cache data must represent the full
ordered beta path rather than a set of independent restarts from `beta=0`.

The cache identity must include:

- algorithm version;
- calculation mode;
- `epsilon_0`, `mu`, and `eta`;
- the complete beta grid;
- `K_plot`;
- guard-root count;
- continuation settings;
- global-guard spacing;
- strict-fallback settings.

`fast_plot`, `certified_audit`, and `plot_only` must not share an
indistinguishable cache contract. `plot_only` may read spectrum CSV data but
must not create or mutate a root cache.

## Performance Reporting

Every new frequency-map workflow must report:

```text
calculation_mode
n_cases
n_models
n_beta_targets
parent_solve_count
continuation_step_count
adaptive_insertions
path_restart_count
local_matrix_evaluations
local_SVD_calls
guard_matrix_evaluations
guard_SVD_calls
strict_fallback_count
strict_fallback_SVD_calls
cache_hits
cache_misses
plot_only_root_calculation_count
```

Wall-clock time is auxiliary because it depends on hardware, process count,
and cache state. Reports must distinguish spectrum-generation cost,
figure-rendering cost, and test/smoke cost. They must not state that PDF
plotting consumed the full root-calculation time.

## Recommended Future CLI Contract

Future map entry points should expose:

```text
--calculation-mode fast_plot
--calculation-mode certified_audit
--calculation-mode plot_only
```

The recommended default for map and figure scripts is `fast_plot`. Reproduction
of certified Step-3A figures should request `certified_audit`. Changes to the
formatting of saved figures should request `plot_only`.

This is a documented future contract only. The present documentation change
does not add these flags or modify an existing CLI.

## Current Counterexample Figures

The existing
[`S3_12_dimensional_frequency_vs_beta.pdf`](../../results/eb_timo_counterexample_dimensional_frequency_beta/S3_12_dimensional_frequency_vs_beta.pdf)
and
[`S3_14_dimensional_frequency_vs_beta.pdf`](../../results/eb_timo_counterexample_dimensional_frequency_beta/S3_14_dimensional_frequency_vs_beta.pdf)
were generated from certified branch-continuation data by
[`plot_counterexample_dimensional_frequency_beta.py`](../../scripts/analysis/thickness_mismatch/maps/plot_counterexample_dimensional_frequency_beta.py).

They are valid certified outputs, not failed or suspect figures. Their
generation time reflects spectral verification work, including local and
cluster continuation, guards, full-matrix checks, and strict recovery, rather
than the cost of rendering two PDF files. They must not be recomputed merely
to change presentation. The saved CSV is the source for any subsequent
`plot_only` rendering.
