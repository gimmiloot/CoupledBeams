# Terminology policy

## Priority

Use terminology in this order:

1. the user's current instruction;
2. the latest approved manuscript;
3. current project theory, assumptions, source notes, and notation policy;
4. the primary source for the model;
5. this skill's corpus-derived phrases;
6. a neutral descriptive expression.

Older corpus wording never overrides an approved current term.

## Exact repetition is acceptable

Repeat a precise term whenever it denotes the same object. Do not substitute decorative alternatives such as:

- `собственные частоты` → `спектральные показатели`;
- `изменение спектра` → `спектральная реорганизация`;
- `сближение кривых` → `кластеризация ветвей`;
- `угол ориентации армирующих волокон` → `ориентационный параметр`.

## Before using a new term

Check whether the needed wording exists in the current manuscript, project documentation, corpus, or primary source. If not, prefer a plain description. If a technical term is genuinely needed, keep it outside the finished passage until the user approves it or introduce it explicitly with its basis.

## Keep models distinct

Do not mix terminology belonging to Euler--Bernoulli, Rayleigh, Timoshenko, anisotropic, three-dimensional, FEM, or reduced models. Do not describe one model as true, complete, or exact unless the current source of truth explicitly does so.

Distinguish:

- physical frequency, angular frequency, and a dimensionless frequency parameter;
- difference, relative difference, and error when no reference truth is defined;
- sorted spectral index, current sorted position, and continuation-defined `branch_id`;
- observed curve proximity, veering, crossing, localization, and tracked mode exchange.

## Modal identity, modal composition, and localization

- A sorted spectral position is not a tracked branch. Ordered frequencies alone establish only spectral order, convergence, or separation.
- Modal composition describes the components expressed in one form. Component fields, energy fractions, coupling coefficients, and partial-mode projections may characterize that composition, but do not determine modal ancestry.
- Localization describes where displacement or energy is concentrated. The dominance of total bending, shear, or torsion energy is not localization; require a spatial distribution between or along structural elements.
- In a multiparameter study, a branch label is defined only along an explicitly declared continuation path. Do not assume that individual modal identity is independent of the path in parameter space.
- Near close or repeated eigenvalues, individual shape matching may be ambiguous. The stable object may be the corresponding modal subspace rather than a particular eigenvector.
- Mutual mode transformation requires tracked forms or modal subspaces and evidence that their character changes. A change in modal composition without tracking does not establish an exchange of descendants.

For CoupledBeams, follow the continuation definition and matching requirements in the current project rules. This skill does not replace them with a universal continuation algorithm or a universal matching threshold.

## Source terms

Do not replace a term taken from a primary source with a more elegant synonym. When source terminology conflicts with the current project's approved terminology, report the mismatch rather than blending the terms silently.
