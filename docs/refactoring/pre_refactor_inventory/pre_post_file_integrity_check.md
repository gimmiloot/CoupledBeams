# Pre/post file integrity check

The baseline was captured before creating this inventory. It contains every
pre-existing tracked, untracked, and ignored worktree file reported by Git,
excluding `.git` internals.

| Check | Result |
| --- | ---: |
| Pre-existing files in baseline | 9,759 |
| Pre-existing files found after work | 9,759 |
| `pre_existing_files_deleted` | 0 |
| `pre_existing_files_moved` | 0 |
| `pre_existing_files_content_changed` | 0 |
| New files in final worktree | 19 |
| `new_files_outside_allowed_inventory_dir` | 0 |

Required machine-readable result:

```text
pre_existing_files_deleted = 0
pre_existing_files_moved = 0
pre_existing_files_content_changed = 0
new_files_outside_allowed_inventory_dir = 0
```

The pre-existing and post-existing manifests have the same SHA256 because
their sorted paths, sizes, mtimes, and content hashes are identical:

```text
d2f341b313f54304f6586479587afce8766bd125881ff8cdf130cd5bea2cf043
```

The final post manifest is stored in the external backup as
`post_existing_files_manifest.csv`. New files are exclusively the nine
Markdown audit documents, this integrity report, and nine CSV manifests under
`docs/refactoring/pre_refactor_inventory/`.

No pre-existing file was restored because no mismatch was detected. Git
metadata are outside the content-integrity scope; the verified Git bundle
separately protects committed refs.

