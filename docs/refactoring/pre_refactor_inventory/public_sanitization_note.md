# Public sanitization note

The pre-refactor inventory was reviewed before it became part of the public
documentation. The public copy now uses `<repository-root>` and
`<external-local-backup>` instead of absolute local paths. Individual paths
and example filenames inside `private/` were removed from the two aggregate
inventory CSV files; only counts, total size, category, and backup status
remain.

No credentials, tokens, local user names, email credentials, or private-file
contents were found in the inventory. Repository-relative scientific and
generated-result paths remain because they describe the public project.

The complete unsanitized evidence, including the original manifests and the
verified ignored-data snapshot, remains in the external local backup. That
backup was not modified by this sanitization. No scientific file, generated
result, source snapshot, or external snapshot manifest was changed; only the
derived public inventory copy was redacted.
