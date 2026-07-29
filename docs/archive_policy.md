# Archive and Preservation Policy

> Archive means a change of maintenance status, not deletion of scientific
> history.

This policy applies to code, documentation, and generated scientific evidence.
It does not authorize moving or deleting any file. Candidate classifications
remain subject to scientific, dependency, and provenance review.

## Active/stable

Active or stable material is used by current workflows, has a documented
verification role, or forms part of the scientific source of truth. It may be
changed only through the normal review and testing process. Stable baselines
must keep their mathematical and numerical contracts explicit.

## Soft archive

A soft-archived file stays at its current path and receives a status such as
`historical`, `completed`, `closed-research`, or `superseded`. It links to its
successor or canonical report, is preserved for reproduction and provenance,
and should not gain unrelated new functionality. Soft archive is the default
for completed research workflows because it preserves commands, imports, and
historical evidence.

## Hard archive

A physical move is permitted only after all of these gates pass:

1. dependency and import audit;
2. documentation-reference audit;
3. a reviewed replacement already exists;
4. relevant regression tests pass;
5. the move uses `git mv`;
6. an old public CLI path remains available through a wrapper, or a migration
   note explicitly records the change;
7. numerical results and output contracts are unchanged.

Hard archive is not performed in the current documentation/navigation stage.
Zero importers alone is not evidence that a runnable scientific script is
unused.

## Generated/local results archive

The ignored `results/` tree combines caches, active diagnostics, completed
studies, negative results, and canonical closure evidence. Before any move, it
requires an external backup and SHA256 manifest. A canonical conclusion must
also be promoted into tracked documentation with its provenance and
reproduction entry point. Generated data and caches are never deleted merely
because they can probably be recomputed.

`results/_smoke/` is temporary and normally regenerable, but it may be cleaned
only by a separate, explicitly authorized task after checking that it contains
no unique report or active-process output. See the [generated results
index](results_index.md) for the current directory classification.

## Never delete automatically

The following are never automatic deletion candidates:

- verified theory and equations;
- recorded assumptions and literature source notes;
- determinant implementations and unknown-order/sign conventions;
- branch-tracking and spectrum-completeness machinery;
- the production FEM baseline and transform convention;
- canonical tests and regression fixtures;
- counterexamples and negative results;
- stage-closure reports;
- provenance and backup manifests.

The detailed pre-refactor candidate audit remains available in the
[inventory](refactoring/pre_refactor_inventory/archive_candidates.md). It is a
review aid, not an approved move plan.
