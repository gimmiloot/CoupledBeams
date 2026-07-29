# Pre-refactor snapshot and inventory

This directory is a read-only inventory of the repository state captured before
any structural refactoring. No production code, scientific formula,
determinant, solver, test, existing result, or existing documentation file was
changed while preparing it.

No file listed as an archive candidate has been declared unnecessary. Every
archive classification is a draft that requires a manual scientific and
dependency review. In particular, zero Python importers is not evidence that a
runnable research script is unused.

The SHA256 manifests under `manifests/` fix the tracked and key scientific
state at HEAD `df701f36569723444e7b131e7cfdad8c894889db`. The full pre-existing
worktree manifest, committed-ref bundle, ignored/untracked snapshot, and
per-file results manifest are stored in the local external backup
`CoupledBeams_pre_refactor_20260729_094757`; private data were not sent to an
external service.

Inventory contents:

- `pre_refactor_snapshot.md`: Git, environment, frozen hashes, and scientific
  anchors;
- `repository_inventory.md`: structure, composition, and ownership risks;
- `research_topics_inventory.md`: status of the visible research lines;
- `script_status_draft.md`: preliminary runnable-workflow status map;
- `archive_candidates.md`: retain/soft/hard/local/manual-review draft;
- `refactor_risk_register.md`: risks and gates for a later refactor;
- `documentation_consistency_audit.md`: path and status consistency findings;
- `anisotropic_preparation_notes.md`: preparation only, not an implementation;
- `pre_post_file_integrity_check.md`: final byte-integrity comparison;
- `public_sanitization_note.md`: public-copy privacy redactions, with full
  local evidence retained externally;
- `manifests/`: machine-readable inventories.

The next refactoring stage is permitted only after manual review of this
inventory and explicit decisions on canonical entry points, historical
workflows, generated results, and the proposed anisotropic data model.
