# Documentation consistency audit

The path scanner inspected tracked Markdown, reStructuredText, LaTeX, and text
files under the requested roots. It recorded 904 distinct references:

- 60 Markdown links, all of which resolve in this checkout;
- 844 inline/code-style path references;
- 371 inline paths absent from the snapshot;
- 345 of the absent paths are under `results/` and are usually documented
  expected/generated outputs rather than broken navigation links;
- 26 absent non-results inline paths require context review; several are
  historical, proposed-future, shorthand, or optional article-workspace paths.

Therefore “371 absent path references” must not be simplified to “371 broken
links.” The exact classification is in
`manifests/documentation_path_references.csv`.

## Findings

| Document | Current statement | Newer conflicting evidence | Severity | Recommended later correction |
| --- | --- | --- | --- | --- |
| `docs/literature/source_index.md` | Lists `docs/literature/pdf/02.-Dissertatsiya_M.-K.-CHien-_21.03.25_.pdf` as a local PDF | The file is absent from the tracked and ignored snapshot | high | Locate the authoritative local source or mark it unavailable; do not substitute a different PDF silently |
| `README.md` | Presents `paper_thickness_mismatch_timoshenko/` as an existing skeleton | The directory is absent from this checkout and is not currently ignored | medium | After locating the authoritative workspace, state whether it is optional/external/ignored or remove the present-tense layout claim |
| `docs/project_rules.md` | Refers to `paper_thickness_mismatch_timoshenko/README.md` | That workspace and file are absent | medium | Resolve workspace ownership first; then use an explicit optional/local-workspace note |
| `docs/thickness_mismatch/frequency_crossing_verification_status.md` | Refers to `paper_thickness_mismatch_article/main.tex` | The article workspace is absent | medium | Restore/locate the workspace or record that the path is external to this checkout |
| `scripts/README.md` | `plot_article_fig3_with_fp_and_ff_refs.py` reuses `paper_dorofeev_style/generate_article_spectral_figures.py` | `paper_dorofeev_style/` is absent, so this documented command cannot be assumed runnable in this checkout | high | Add an explicit prerequisite and workspace-location check; do not copy article code into the diagnostic script without review |
| `docs/veering/figures/figure_index.md` | Gives `paper_dorofeev_style/figures/lambda_mu_beta15_r5mm.png` as source | External source path is absent, although a curated tracked copy exists in `docs/veering/figures/` | medium | Preserve the curated-copy provenance hash and mark the original workspace as unavailable locally |
| `docs/veering/figures/figure_index.md` | Says “see `data/dataset_index.md`” from a nested `figures/` file | The actual tracked path is `docs/veering/data/dataset_index.md`; inline code is ambiguous when interpreted relative to the current file | low | Later replace with an explicit relative Markdown link `../data/dataset_index.md` |
| `docs/theory/out_of_plane_eb_torsion.md` and thickness script map | Name `results/out_of_plane_fem_validation/` and other out-of-plane output directories | None of those generated directories is present in this snapshot | medium | State “expected output” versus “verified output present”; restore canonical generated reports before promoting validation claims |
| `README.md` | Tests section names only `tests/test_analytic_smoke.py` | 24 tracked test files now cover formulas, branch tracking, completeness, safe-prefix, FEM transforms, and out-of-plane code | low | Expand public test guidance after inventory approval; preserve warnings about heavy/root-solving tests |
| `scripts/analysis/thickness_mismatch/README.md` | Contains a future example `shapes/plot_beta_scan.py` | The proposed unified script does not exist | informational | Keep explicitly labeled as proposal; do not turn the example into an active command until implemented and regression-tested |
| `CHANGELOG.md` | Contains many paths under absent `paper_dorofeev_style/` | Changelog entries describe historical changes, not necessarily current checkout availability | informational | Do not rewrite historical changelog entries; add a current workspace-location note elsewhere if needed |
| safe-prefix plan and script map | Retain `not_ready_for_step3`, `ready_for_targeted_step3`, and `rule_B_safety_survives_cost_test_required` in chronological/Phase-I sections | Later Step-3A and negative cost results exist | informational | No factual correction needed: these sections are explicitly historical or Phase-I; keep the final closure prominent |

## Stale-status assessment

One active present-tense status conflict is counted: the root README describes
the absent `paper_thickness_mismatch_timoshenko/` skeleton as part of the
current layout. The incomplete root test overview is a coverage omission, not
a scientific status conflict.

The requested Step-3 and Rule-B checks found no unlabelled active contradiction:

- the general completeness `not_ready_for_step3` result and later gateway
  `ready_for_targeted_step3` result are different chronological workflows;
- the gateway truthfully says it wrote but did not execute the future manifest;
  the separate Step-3A row and closure document record the later execution;
- `rule_B_safety_survives_cost_test_required` is retained as an immutable
  Phase-I status, while the final sections clearly record
  `rule_S_cost_not_beneficial` and stage closure;
- the research plan explicitly labels older “remaining steps” as historical
  and superseded.

No documentation was repaired in this task.

