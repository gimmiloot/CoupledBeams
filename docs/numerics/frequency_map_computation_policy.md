# Frequency-Map Computation Policy

**Policy identifier:** `frequency-map-v1`

## Scope and governing principle

This document defines the project-wide computation contract for
frequency-versus-parameter maps. It applies to isotropic coupled beams,
thickness-mismatch models, anisotropic rods, laminated beams, and future
eigenfrequency models in this repository.

The governing principle is:

> Ordinary frequency-map production is not the same task as research-grade
> spectral certification.

An ordinary map uses `calculation_mode: fast_plot` by default. A
`certified_audit` is not a prerequisite for producing a figure, a mandatory
validation of a new universal runner, or the standard calculation at every
grid point. It also does not automatically require a spectral tail beyond the
range needed by the figure.

The introduction of a new physical model is not, by itself, a reason to
develop a universal sweep runner, compute full inventories, certify every root
with independent global scans, or postpone the first map until an
infrastructure audit is complete. Any deviation from `fast_plot` must record a
`policy_override_reason` tied to an explicit scientific purpose or a concrete
numerical trigger.

This policy is a numerical orchestration contract. It does not change a
model's governing equations, determinant, state ordering, coordinate
conventions, root-quality thresholds, or physical parameters. Those remain
under the authority of the model-specific theory, implementation, and
validation documents.

## Calculation modes

| Mode | Purpose | Root calculations |
| --- | --- | --- |
| `fast_plot` | Ordinary frequency maps on a declared finite grid | Sequential, model-specific calculation of the plotted range and one completeness guard |
| `certified_audit` | Reference data, counterexamples, solver validation, or unresolved numerical events | Explicitly scoped independent or global checks required by the audit question |
| `plot_only` | Rendering from completed spectrum data | None |

These modes have distinct metadata and cache identities. A figure entry point
must state which mode produced its spectrum data.

## Required local policy instance

Every frequency map must declare a local instance of this policy. The instance
defines the map's semantics, grid, plotted range, guard, and recovery behavior.
It must include:

- `frequency_map_policy`;
- `calculation_mode`;
- `spectrum_semantics`;
- `sweep_parameter`;
- `parameter_grid`;
- `K_plot`;
- `K_guard`;
- `guard_root_role`;
- `neighbour_audit`;
- `local_repair_policy`;
- `strict_audit_default`.

A recommended template is:

```yaml
frequency_map_policy: frequency-map-v1
calculation_mode: fast_plot
spectrum_semantics: sorted_positions
sweep_parameter: mu
parameter_grid: 0.00:0.02:0.80
K_plot: 8
K_guard: 9
guard_root_role: completeness_only
neighbour_audit: enabled
local_repair_policy: triggered_only
strict_audit_default: false
```

The values in this template are examples, not fixed values for every model.
In the default contract, `K_guard = K_plot + 1`, but each workflow must state
both values explicitly. A local note may add model-specific fields without
redefining the project-wide modes.

If a workflow overrides this policy, its local instance must add
`policy_override_reason` and identify the scientific question or numerical
event that requires the override.

## Spectrum semantics

### `sorted_positions`

With `spectrum_semantics: sorted_positions`, each parameter point uses its own
independently sorted eigenfrequencies.

- Branch tracking is not performed by default.
- MAC is not computed by default.
- Mode shapes are not computed by default.
- A change in the physical character of neighbouring modes does not prevent
  construction of the sorted spectrum.
- The workflow checks completeness, ordering, root quality, and reasonable
  evolution relative to neighbouring parameter points.
- A sorted index must not be described as a modal descendant or branch
  identity.

### `tracked_branches`

With `spectrum_semantics: tracked_branches`, branch identity is a separate
scientific task. The workflow must follow the branch-identity rules in
[`docs/project_rules.md`](../project_rules.md). It may require mode shapes,
MAC, cluster or subspace treatment, and local grid refinement.

Tracked branches must not be enabled automatically for an ordinary
sorted-position plot. A local workflow must state its continuation start,
identity rule, acceptance diagnostics, and behavior near close or repeated
roots.

## `fast_plot` workflow

For an ordered parameter grid

```text
q_0 < q_1 < ... < q_M
```

`fast_plot` follows this contract:

1. Compute `q_0` with the existing reliable model-specific solver.
2. Traverse the requested grid sequentially.
3. Use one or two completed neighbouring points only to predict root
   locations, choose local search windows, and estimate an adequate upper
   search bound.
4. Never use interpolation or a predictor as a final eigenfrequency.
5. At every point, compute roots `1...K_plot` and root
   `K_guard = K_plot + 1`.
6. Do not compute an additional spectral tail by default.
7. Refine every accepted root with the existing model-specific refiner.
8. For plotted roots, verify sorting, all frozen residual and quality gates,
   absence of unresolved candidates below the plotted range, and consistency
   with neighbouring points.
9. After the primary sweep, perform a neighbour-jump audit.
10. If the audit marks a point as suspicious, recalculate only the affected
    point. The local repair may narrow or densify its root search and, when
    needed, add an intermediate parameter point.
11. Use strict or global recovery only after a recorded unresolved numerical
    event.

The predictor is a numerical locator. It does not define a final root or a
modal branch. A completed point must be written only after its local quality
contract passes.

## Neighbour audit and local repair

The neighbour audit compares each plotted sorted position with adjacent grid
points using a robust, predeclared jump criterion. It is a trigger for
verification, not a data transformation.

At a suspicious point, the workflow must not smooth the spectrum, replace a
value by interpolation, average neighbouring frequencies, or delete a point
manually. It must instead preserve the original diagnostic and run the declared
`local_repair_policy`. If local repair reproduces the feature, the value is
retained. If it corrects a locator failure, both the original and corrected
provenance must be recorded. An unresolved event remains visible as a failed
or qualified point according to the local output contract.

Close roots must remain distinct whenever the model-specific determinant or
singular-value diagnostics resolve distinct roots. Automatic merging based
only on a small frequency gap is prohibited.

## Plotted roots and the guard root

By default,

```text
K_guard = K_plot + 1
guard_root_role = completeness_only
```

The guard root:

- is not plotted;
- is not an additional scientific branch;
- confirms that exactly `K_plot` accepted roots lie below it;
- must preserve the correct ordering;
- must be separated from the right boundary of the search range;
- must have no unresolved candidate before it;
- must pass the model-specific residual and physical-quality gates required
  for a completeness guard.

The local instance must distinguish the quality criteria for plotted roots
from the completeness and ordering criteria for the guard root. If the guard
frequency is not a plotted result, it need not automatically inherit a stricter
cross-run frequency-agreement gate imposed on plotted roots. A guard-only
qualification must retain its diagnostics and the original strict warning; it
must not hide or relabel that warning as a strict pass.

## Strict recovery and `certified_audit`

Strict global recovery or `certified_audit` requires a recorded trigger, such
as:

- fewer than `K_guard` roots were found;
- an unexplained candidate lies below the guard root;
- sorted ordering is violated;
- a plotted root fails a residual or physical-quality gate;
- a close cluster is not resolved;
- a suspicious neighbour jump remains after local repair;
- two genuinely independent local refinements disagree for a plotted root;
- a reference or certification dataset is being constructed;
- a scientific counterexample is being checked;
- the user explicitly requests independent certification.

The following are not triggers by themselves:

- introduction of a new physical model;
- production of an ordinary illustrative map;
- absence of formal authorization for a universal runner;
- disagreement in the last digits of an auxiliary guard root when plotted
  completeness is already established;
- a presentation-only change to a completed figure.

A strict audit must declare its own finite scope, independent evidence, root
range, and acceptance thresholds. It does not silently redefine the ordinary
map contract.

## `plot_only`

`plot_only` reads completed spectrum CSV data and renders figures. It performs
none of the following:

- matrix assembly;
- determinant evaluation;
- SVD;
- root search or refinement;
- spectrum recomputation;
- cache creation or mutation.

Changes limited to labels, axis titles, fonts, line styles, DPI, or
PNG/PDF/SVG format must use `plot_only`. A test or manifest check should record
that the root-calculation count is zero.

## Performance, checkpoints, and reproducibility

An ordinary map uses one sequential parameter path. Completed points are not
recomputed; `missing-only` execution is allowed and preferred. Each point is a
transaction: calculate, validate the plotted range and guard, write a temporary
file, atomically publish it, and update the checkpoint.

Wide scans and roots beyond `K_guard` are not computed without a recorded
trigger. Caches must be bounded and their identity must include the model
contract, full parameter grid, calculation mode, policy version, `K_plot`,
`K_guard`, and relevant numerical settings.

Spectrum-generation time and figure-rendering time must be reported
separately. When runtime is unknown, measure a representative point before
estimating the complete job. An ETA must not be promised without a measured
basis, and a benchmark must not turn into development of new infrastructure.

The reproducible outputs are normally:

- spectrum CSV;
- figure;
- report;
- run manifest.

The run manifest must store the complete local policy instance, counts of
reused and newly calculated points, local repairs and strict recoveries, and
separate runtime for spectrum generation and rendering.

## Local model documentation

The frequency-map section of a model note records only the local instance and
any justified override. Physical contracts and model-specific quality gates
remain in their established local sections, but they do not reproduce the
universal workflow. The note must link to this document. Existing certified or
historical workflows may retain their original root ranges and provenance;
adoption of `frequency-map-v1` does not rewrite completed evidence.
