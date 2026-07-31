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

The residual refactoring audit is complete and recommends no further code
refactoring for the current scope. Stage 3, the common thickness-scaling
helper, is completed. The value-only
[`thickness_scaling.py`](../../src/my_project/analytic/common/thickness_scaling.py)
module now owns the duplicated mass-preserving denominator polynomial and
`tau1`/`tau2` divisions. The EB and Timoshenko entry points remain legacy
wrappers with their original validation, square-root backends, dataclass return
types, exception messages, and public imports.

## Completed actions

- Stage 0: backup and inventory — completed.
- Stage 1: documentation and navigation — completed.
- Stage 3: common thickness-scaling helper — completed.
- Residual duplicate, soft-archive, possible-orphan, path-reference, and module
  boundary audit — completed.
- Public privacy review of the inventory — completed for this pass.
- Canonical status, result-directory, and preservation navigation — documented
  for manual review.
- Pre/post regression gates passed for 98 factor outputs, 28 exception cases,
  public imports/signatures/dataclass fields/repr, and six fixed-`Lambda`
  matrices for each of EB and Timoshenko, all with zero exact mismatches.
- The local baseline remains identified by
  `pre-common-code-refactor-2026-07-29`.

## Residual refactoring audit

Status: **completed — no further code refactoring recommended**.

The [closing audit](residual_refactoring_audit.md) and its
[candidate](manifests/residual_refactor_candidates.csv),
[possible-orphan](manifests/orphan_script_audit.csv), and
[manual-reference](manifests/manual_path_reference_review.csv) manifests record
why the remaining similar functions retain different contracts, why soft
archive status is sufficient, and why two exact filename utilities are only
future candidates if their existing call sites are already being changed.

## Preserved Stage 3 boundaries

No determinant or coupling-matrix entry, formula sign or coefficient, unknown
ordering, root solver, branch tracker, FEM model, generated result, cache, or
private file changed. The refactor added only the shared value algebra, legacy
wrapper calls, targeted regression tests, and this status documentation. It
performed no eigenvalue root calculation or scientific result regeneration.

## Planned later stages

These are possible review stages, not commitments:

- Later soft-archive status harmonization — documentation only and only where
  the closing audit records a concrete navigation mismatch.
- Stage 3: common thickness-scaling helper — completed.
- Stage 4: effective rod properties interface — not started.
- Stage 5: anisotropic source assimilation is active as a separate diagnostic
  research line; no coupled-rod model or production API has been started.

Stage 4 and any future production anisotropic implementation require separately
approved regression gates, a physical model, and an isotropic-limit contract.

## Anisotropic direction

The [anisotropic direction](../anisotropic_rods/README.md) is now an
`active-diagnostic` source-assimilation line. Yartsev's 2024 monograph is
registered as its primary source, and the Chapter-2 single free-free rod gate
is complete. Coupled rods and a production anisotropic API remain not started.
This diagnostic work is not a justification for speculative abstractions in
the refactoring line; any future refactoring step must be evaluated
independently.

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
