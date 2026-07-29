# Refactoring Status

## Safety principle

Repository refactoring must improve ownership and navigation without silently
changing equations, determinants, unknown ordering, root multiplicity, branch
identity, FEM conventions, or generated scientific evidence. Each later stage
requires a separately reviewed scope.

## Verified pre-refactor snapshot

Stage 0 created a verified external local backup and the public
[pre-refactor inventory](pre_refactor_inventory/README.md). The public copy is
privacy-sanitized; its [sanitization note](pre_refactor_inventory/public_sanitization_note.md)
records what was removed while the complete evidence remains in the external
backup.

## Current stage

Stage 1 is the documentation and navigation pass. It creates the
[research index](../research_index.md), [generated results index](../results_index.md),
[archive policy](../archive_policy.md), and [script status map](../../scripts/STATUS.md).
It does not refactor scientific code or physically archive anything.

## Completed actions

- Stage 0: backup and inventory — completed.
- Stage 1: documentation and navigation — current.
- Public privacy review of the inventory — completed for this pass.
- Canonical status, result-directory, and preservation navigation — documented
  for manual review.

## Prohibited changes in this stage

No Python code, test, formula, determinant, root solver, branch tracker, FEM
model, generated result, cache, or private file may be changed. Files are not
moved, renamed, deleted, or physically archived. Scientific calculations are
outside this stage.

## Planned later stages

These are possible review stages, not commitments:

- Stage 2: soft archive classification — documentation only.
- Stage 3: common thickness-scaling helper — not started.
- Stage 4: effective rod properties interface — not started.
- Stage 5: anisotropic research scaffold — not started.

Stages 3--5 require new regression gates and, for anisotropy, an approved
physical model and isotropic-limit contract before implementation.

## Manual path-reference review

The inventory found 371 absent inline path references. Of these, 345 point
under generated `results/` and are handled by the generated-results policy;
absence from a fresh clone is expected. The remaining 26 references require
manual context review and are grouped below without changing or fabricating
any path.

| Source group | References | Count | Review needed |
| --- | --- | ---: | --- |
| Historical `CHANGELOG.md` entries | `paper_dorofeev_style/*` and one abbreviated `scripts/plot_flat_mu_bending_desc*` path | 14 | Preserve history; locate the old workspace/command before adding current navigation. |
| Project rules | local article figure/workspace paths | 2 | Determine whether they are optional local prerequisites; the rules file is intentionally unchanged in this pass. |
| Literature index | one named dissertation PDF | 1 | Locate the authoritative source or mark it unavailable in a separately approved literature update. |
| Thickness crossing status | one local article `main.tex` path | 1 | Resolve article-workspace ownership before editing the status note. |
| Veering figure index | one ambiguous relative data path and one local article source figure | 2 | Prefer the tracked curated copy; review the relative link separately. |
| Thickness script map | one future manifest shorthand and one explicitly proposed future shape CLI | 2 | Keep clearly historical/proposed; do not create files to silence the scanner. |
| Internal helper README | one shorthand `scripts/reports` reference | 1 | Clarify whether it denotes a concept or a real future directory. |
| Main scripts README | three local article-workspace references | 3 | Treat the workspace as an external/local prerequisite, not part of the tracked checkout. |

The exact scanner rows remain in
[`documentation_path_references.csv`](pre_refactor_inventory/manifests/documentation_path_references.csv).
No tracked Markdown link was broken in the pre-refactor snapshot.
