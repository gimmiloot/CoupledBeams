# Migration Log

## Canonical package after consolidation

The canonical `veering_bridge` package is now:

- `README.md`
- `foundations.md`
- `theorem_map.md`
- `lemma_statements.md`
- `proof_notes.md`
- `migration_log.md`

No internal `archive/` folder is kept as part of the canonical package.

## Consolidation policy

- old bridge notes were treated as working drafts, not as permanent canonical
  files;
- overlapping status files and registries were absorbed into
  `theorem_map.md`;
- clean theorem-language statements were absorbed into
  `lemma_statements.md`;
- proof attempts, cautions, failure modes, constructive comments, and
  ideal-model warnings were absorbed into `proof_notes.md`;
- definitions and conventions were absorbed into `foundations.md`;
- no new archive was created inside `docs/theory/veering_bridge/`.

## Old-file mapping

| old file(s) | absorbed into | disposition |
| --- | --- | --- |
| `README.md` | `README.md` | rewritten as the canonical overview |
| `definitions_and_language.md`, `assumptions_and_regular_domain.md`, `right_left_data_definitions.md` | `foundations.md` | removed after merge |
| `bridge_theorem_skeleton.md`, `bridge_lemma_hypotheses_update.md`, `conditional_bridge_lemma.md`, `claim_registry.md`, `conditional_bridge_registry.md` | `theorem_map.md`, `lemma_statements.md`, `proof_notes.md` | removed after merge |
| `packet_construction_lemma_simple_root.md`, `packet_construction_lemma_status.md`, `packet_construction_lemma_registry.md` | `lemma_statements.md`, `theorem_map.md`, `proof_notes.md` | removed after merge |
| `packet_construction_program.md`, `packet_construction_registry.md`, `packet_construction_next_step.md`, `packet_failure_modes.md`, `constructive_packet_workflow.md`, `constructive_packet_registry.md`, `constructive_bridge_next_step.md` | `proof_notes.md`, `theorem_map.md` | removed after merge |
| `complement_regular_lemma.md`, `complement_regular_status.md`, `complement_regular_registry.md` | `lemma_statements.md`, `theorem_map.md`, `proof_notes.md` | removed after merge |
| `schur_root_capture_lemma.md`, `schur_root_capture_status.md`, `schur_root_capture_registry.md` | `lemma_statements.md`, `theorem_map.md`, `proof_notes.md` | removed after merge |
| `lambda_normalization_lemma.md`, `lambda_normalization_status.md`, `lambda_normalization_registry.md` | `lemma_statements.md`, `theorem_map.md`, `proof_notes.md` | removed after merge |
| `first_order_reduced_model_lemma.md`, `first_order_reduced_model_status.md`, `first_order_reduced_model_registry.md` | `lemma_statements.md`, `theorem_map.md`, `proof_notes.md` | removed after merge |
| `frozen_model_root_comparison_lemma.md`, `frozen_model_root_comparison_status.md`, `frozen_model_root_comparison_registry.md` | `lemma_statements.md`, `theorem_map.md`, `proof_notes.md` | removed after merge |
| `symbolic_checks.md` | `proof_notes.md` | removed after merge |
| `archive/*.md` | `foundations.md`, `theorem_map.md`, `proof_notes.md`, `migration_log.md` | archive removed from canonical package |

## Files treated as duplicate or superseded

The following categories were judged non-canonical after consolidation:

- narrow `*_status.md` files;
- narrow `*_registry.md` files;
- overlapping theorem skeletons;
- overlapping packet-program notes;
- the local `archive/` subtree.

Their useful content now lives in the canonical files above.

## Canonical versus non-canonical outcome

- keep canonical:
  `README.md`, `foundations.md`, `theorem_map.md`, `lemma_statements.md`,
  `proof_notes.md`, `migration_log.md`;
- do not keep a markdown archive inside `docs/theory/veering_bridge/`;
- rely on git history for long-term archival rather than maintaining a second
  markdown forest inside the canonical package.
