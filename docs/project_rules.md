# Project Rules

## Purpose

This document collects project-wide rules that apply across analytic scripts,
FEM diagnostics, theory notes, and exploratory studies. Detailed derivations and
case-specific parameters stay in their local documents; this file records the
global conventions that should not be rediscovered or redefined in each script.

## Branch Identity and Mode Descendants

- `sorted mode k` or `k-th sorted frequency` means the k-th eigenfrequency in
  increasing order at one fixed parameter point.
- `descendant branch k` means the continuation of the k-th mode shape selected
  at the tracking start point.
- In this project, `branch k` means a descendant branch by default, not the k-th
  sorted frequency at every parameter value.
- A `sorted_position` or `current_sorted_index` is diagnostic metadata: it says
  where the current frequency of a descendant branch sits in the sorted
  spectrum.
- A branch must not be renamed only because its frequency becomes closer to a
  different sorted root.
- Plots and reports must state whether they show sorted frequencies or
  descendant branches.

For the baseline analytic branch family, the default seed is
`beta = 0`, `mu = 0` for each `epsilon`. For fixed-eta
thickness-mismatch diagnostics, the local seed is usually `mu = 0`.

## Sorted Frequencies vs Descendant Branches

Nearest-frequency tracking is not reliable in veering windows, close
interactions, near-crossings, or quasi-degenerate regions. Use shape continuity,
MAC/correlation data, and local refined grids when branch identity matters.

MAC or shape correlation is evidence for tracking, but low-MAC assignments are
not canonical. If an assignment has low MAC, low margin, or a large sorted
position jump:

- mark the step as unresolved;
- keep the previous accepted canonical sorted position;
- record the raw candidate sorted position separately;
- require a refined local shape analysis before interpreting a branch exchange.

Jumps larger than one sorted position, such as `5 -> 7`, are suspicious until
confirmed by a local refined analysis of mode shapes.

## Applicability of the Thin-Rod Theory

For the current theory of thin circular rods, applicability is checked using
the diameter-to-length ratio for each rod:

```text
2*r_i/l_i <= 0.1,  i = 1, 2.
```

Here thickness means the circular-section diameter `2*r_i`. Do not replace
this with the radius criterion `r_i/l_i <= 0.1`.

If either rod violates the criterion, that part of the result is questionable.
Diagnostic plots must show such curve segments as dashed, and the script or
report must warn with:

- the relevant parameter values;
- the affected `mu` range;
- the rod number;
- the maximum value of `2*r_i/l_i` on the computed grid.

For the mass-preserving thickness-mismatch eta parameterization,

```text
l1 = l*(1 - mu)
l2 = l*(1 + mu)
r1 = r0*tau1
r2 = r0*tau2
epsilon = r0/(2*l)
```

so

```text
2*r1/l1 = 4*epsilon*tau1/(1 - mu)
2*r2/l2 = 4*epsilon*tau2/(1 + mu)
```

## Diagnostic Studies and Article Workflow

Diagnostic studies are exploratory and must not silently modify article
artifacts or baseline models.

- Diagnostic scripts write their outputs under `results/` unless a task
  explicitly says otherwise.
- Article figures under `paper_dorofeev_style/figures/` are not overwritten by
  diagnostics.
- `paper_dorofeev_style` is not edited in diagnostic tasks unless the user
  explicitly requests article work.
- New diagnostic models do not replace the baseline equal-radius determinant,
  old solvers, or the baseline FEM physical model.
- Diagnostic conclusions must be labeled as diagnostic until backed by the
  required consistency checks and evidence.
- Diagnostic results become article material only through the promotion
  workflow in `docs/writing/article_workflow.md`.

## Script Proliferation Control

Do not create a new script for every parameter set, plot variant, or output
path.

If only `beta`, `mu`, `epsilon`, `eta`, the number of branches, the reference
family, output path, FEM mesh, plotting style, or CSV/report generation changes,
parameterize the existing script or add a documented preset instead.

A new script is acceptable only when:

- it introduces a genuinely new workflow;
- the input/output contract changes;
- it will be a stable reusable entry point;
- merging it into an existing script would make the code unclear or unsafe.

One-off exploratory scripts should eventually be handled in one of these ways:

- documented as historical/one-off;
- replaced by a parameterized entry point;
- removed after a documented reproducible command exists.

Every proposed new script should answer:

1. why is this not a parameter or preset of an existing script?
2. what helper does it reuse?
3. is it article-facing, diagnostic-only, compatibility wrapper, or one-off?

## Consistency Checks for New Model Extensions

New analytic model extensions must preserve the verified baseline unless the
user explicitly asks for a mathematical revision.

- Do not silently change determinant entries, signs, coefficients, unknown
  ordering, or solver behavior.
- Check limiting cases such as `eta = 0`, `mu = 0`, `beta = 0`, symmetric
  limits, and relevant single-rod limits.
- For thickness mismatch, `eta = 0` must reproduce the baseline equal-radius
  model.
- For rod-exchange extensions, check the appropriate swap symmetry; in the
  mass-preserving eta model this is `(mu, eta) -> (-mu, -eta)`.
- Record assumptions and new notation in the appropriate theory or diagnostic
  documentation.

For theory-facing work, verified formulas in `docs/theory/equations.tex` have
priority over assumptions, literature notes, and code. If theory and code
disagree, record the mismatch first; do not silently "fix" either side.

## Eta-Model Rules

The thickness-mismatch eta model is diagnostic-only unless explicitly promoted
by a later task. Its global methodological rules are:

- fixed total length `l1 + l2 = 2l`;
- `mu = (l2 - l1)/(l1 + l2)`, so `l1 = l*(1 - mu)` and
  `l2 = l*(1 + mu)`;
- `eta = (r2 - r1)/(r1 + r2)`;
- positive `eta` means rod 2 is thicker than rod 1, and negative `eta` means
  rod 1 is thicker than rod 2;
- the mass-preserving factors `tau1`, `tau2` keep the total mass fixed relative
  to the equal-radius base radius `r0`;
- `eta = 0` must reduce to the baseline equal-radius model;
- swap symmetry should be checked when the two rods are exchanged.

Detailed formulas live in `docs/thickness_mismatch/README.md`.

## Veering and Modal-Exchange Claims

Strict veering claims require more than close sorted frequencies.

A strict veering claim must have:

- tracked descendant branch identities;
- a small local spectral gap between the tracked branches;
- evidence that this is not a true crossing;
- mode-shape evidence, such as pairwise MAC/correlation or equivalent shape
  diagnostics;
- evidence that modal-character change occurs in the same interaction window.

Sorted frequency closeness alone is not enough. A change in sorted position is
not by itself evidence of veering or modal exchange. Localization claims should
be supported by shape or energy metrics, not only by visual inspection.

When evidence is incomplete, use cautious labels such as modal-character
reorganization, quasi-degeneracy-like behavior, or localization tendency rather
than strict veering.

## Nondimensionalization and Parameter Conventions

Project notation follows the local theory and article-facing notation unless a
diagnostic document explicitly states a local extension.

- `l = (l1 + l2)/2` is the base length for the two-arm system.
- `mu = (l2 - l1)/(l1 + l2)`.
- For circular rods, `epsilon = r/(2*l)` in the equal-radius model.
- `Lambda` is the dimensionless frequency parameter used by the analytic
  determinant. Equivalently, for the baseline equal-radius scaling,
  `Lambda^4 J/(S*l^2) = omega^2*l^2*rho/E`.
- The external ends are clamped and the joint is treated as a rigid connection
  in the baseline model.
- Reference single-rod curves in Lambda(mu) diagnostics are interpretation
  aids, not replacements for branch-shape analysis.

## FEM Comparison Rules

FEM comparisons are validation and diagnostic tools unless a task explicitly
asks for FEM model changes.

- Do not change the FEM physical model, discretization, constitutive model, or
  boundary conditions as part of a diagnostic comparison.
- The production right-arm transform convention is local-to-global:
  `q_global = T @ q_local`, with right axial mapped to `(cos beta, sin beta)`
  and right transverse mapped to `(-sin beta, cos beta)`.
- Production assembly uses
  `K_global = T @ K_local @ T.T` and
  `M_global = T @ M_local @ T.T`.
- FEM frequency comparisons must state the normalization used to compare FEM
  angular frequencies with analytic `Lambda`.
- Shape comparisons must state coordinate orientation, sign alignment, and
  whether they compare sorted modes or tracked descendants.

## Where Detailed Rules Live

- `docs/thickness_mismatch/README.md` -- eta-model formulas, thin-rod
  applicability rule, diagnostic outputs, and thickness-mismatch scripts.
- `docs/veering/terminology.md` -- project terminology for veering,
  quasi-degeneracy, modal-character exchange, and localization.
- `docs/veering/strict_veering_assessment.md` -- evidence criteria and current
  limitations for strict veering claims.
- `docs/veering/literature_assessment.md` -- literature-specific cautions and
  how sources should and should not be used.
- `docs/theory/assumptions.md` -- active assumptions and notation notes.
- `docs/theory/AGENTS.md` -- detailed theory-edit rules, source-of-truth
  priority, and formula-edit discipline.
- `docs/theory/equations.tex` -- verified local theory formulas.
- `docs/writing/article_workflow.md` -- diagnostic-to-article promotion
  checklist, figure selection workflow, and TODO placeholder policy.
- `docs/literature/source_index.md` -- source-specific warnings, including the
  known sign issue in `2003JSVb.pdf`.
- `scripts/README.md` -- runnable script guide, branch-tracking conventions,
  FEM transform convention, diagnostic command inventory, and links to script
  proliferation rules.
- `scripts/analysis/thickness_mismatch/README.md` -- thickness-mismatch script
  audit, preferred diagnostic entry points, historical scripts, and future
  wrapper/refactor TODOs.
- `scripts/lib/README.md` -- internal helper responsibilities, including
  branch tracking and thickness-mismatch MAC diagnostics.
- `paper_thickness_mismatch_timoshenko/README.md` -- skeleton and boundaries
  for the planned thickness-mismatch / Timoshenko article.
- `README.md` -- repository structure, baseline analytic/FEM entry points, and
  high-level source-of-truth notes.
