# Repository inventory

This is a baseline inventory of the 9,759 pre-existing worktree files. It
excludes `.git` internals and excludes the new files in this inventory itself.
Exact rows are in `manifests/top_level_inventory.csv`.

## Top-level structure

| Directory | Files | Bytes | Tracked / untracked / ignored | README / AGENTS anywhere below | Status | Deletion risk | Later action |
| --- | ---: | ---: | --- | --- | --- | --- | --- |
| `.agents/` | 0 | 0 | 0 / 0 / 0 | no / no | unknown local metadata | low | manual ownership review; do not infer purpose from an empty directory |
| `.codex/` | 1 | 31 | 1 / 0 / 0 | no / no | active configuration | medium | retain; keep separate from scientific code |
| `.pytest_cache/` | 6 | 35,884 | 0 / 0 / 6 | generated README / no | disposable cache | low | unchanged in this task; eligible for later cache cleanup only |
| `data/` | 2 | 0 | 2 / 0 / 0 | no / no | stable placeholders | medium | retain `.gitkeep`; define data ownership before adding datasets |
| `docs/` | 89 | 35,093,583 | 89 / 0 / 0 | nested README / nested AGENTS | active source and evidence | high | retain; separate canonical theory from research notes in navigation only after review |
| `notebooks/` | 1 | 0 | 1 / 0 / 0 | no / no | placeholder | medium | retain `.gitkeep`; no notebook workflow is present |
| `private/` | 52 | 5,549,997 | 0 / 0 / 52 | nested README / no | private local data | critical | retain locally; never publish or move without content-owner review |
| `results/` | 9,317 | 1,363,576,511 | 1 / 0 / 9,316 | no / no | generated and scientific evidence | critical | backup first; distinguish canonical closure evidence, active diagnostics, caches, and smoke output |
| `scripts/` | 205 | 8,752,061 | 157 / 0 / 48 | yes / no | mixed active code | high | keep imports stable; introduce statuses before any move |
| `skills/` | 4 | 9,378 | 4 / 0 / 0 | yes / no | active writing instructions | medium | retain; keep external-writing policy separate from code |
| `src/` | 24 | 179,398 | 16 / 0 / 8 | no / no | stable analytic/FEM source | critical | retain; mathematical audit required before refactor |
| `tests/` | 54 | 1,011,432 | 25 / 0 / 29 | no / no | stable tests plus bytecode | critical | retain tracked tests; caches are disposable only later |

The requested `paper_dorofeev_style/`,
`paper_thickness_mismatch_article/`, and
`paper_thickness_mismatch_timoshenko/` directories are absent from the current
checkout. The root README nevertheless advertises
`paper_thickness_mismatch_timoshenko/`; this is a documentation/manual-review
issue, not a reason to create the directory during this task.

## Composition and ownership observations

- `scripts/` mixes reusable helpers (`scripts/lib/`), presentation entry
  points, compatibility wrappers, research audits, heavy FEM workflows, and
  one-case plotting programs. It is the main mixed-ownership directory.
- `src/my_project/analytic/` contains both stable library modules and three
  historically runnable `Freq*.py` modules. Forty runnable-looking Python
  files are also imported by other code, so moving a CLI based only on its
  filename would break dependencies.
- `results/` contains both reproducible caches and nontrivial scientific
  closure evidence. The largest current subtree is
  `results/eb_vs_timoshenko_3d_validation/` at 814,249,750 bytes; it remains an
  active diagnostic, not an automatic archive target.
- `private/` is ignored but scientifically and operationally material. All 52
  files were included in the local external snapshot with hashes.
- The three requested article workspaces are not hidden by the current
  `.gitignore`; they are simply absent from this checkout. No conclusion about
  their status outside this checkout is possible.

## Empty directories and placeholders

The only physically empty pre-existing directory found outside `.git` was
`.agents/`.

Tracked `.gitkeep` placeholders must not be deleted:

```text
data/external/.gitkeep
data/input/.gitkeep
docs/literature/notes/.gitkeep
docs/theory/derivations/.gitkeep
notebooks/.gitkeep
results/.gitkeep
scripts/analysis/thickness_mismatch/audits/.gitkeep
scripts/analysis/thickness_mismatch/lambda_mu/.gitkeep
scripts/analysis/thickness_mismatch/maps/.gitkeep
scripts/analysis/thickness_mismatch/postprocess/.gitkeep
scripts/analysis/thickness_mismatch/shapes/.gitkeep
src/my_project/analytic/.gitkeep
src/my_project/fem/.gitkeep
src/my_project/io/.gitkeep
src/my_project/utils/.gitkeep
tests/.gitkeep
```

Some of these directories now also contain files; the placeholder remains
part of the committed state and is frozen for this task.

## README coverage

Only `scripts/` and `skills/` have immediate tracked top-level README files.
`docs/` and `private/` have nested README files but no immediate directory
README. `.pytest_cache/README.md` is generated cache documentation, not project
navigation. Immediate project README coverage is absent for `.codex/`,
`data/`, `docs/`, `notebooks/`, `private/`, `results/`, `src/`, and `tests/`.
This is an inventory finding only; no README was added outside this audit
directory.

## Ignored-path classification

| Category | Files | Bytes | Included in backup | Notes |
| --- | ---: | ---: | --- | --- |
| `generated-results` | 9,316 | 1,363,576,511 | yes | includes active, closed, historical, cache, plot, CSV, and report artifacts |
| `private-data` | 52 | 5,549,997 | yes | local-only; no external service used |
| `disposable-python-cache` | 91 | 4,086,091 | no | `.pytest_cache`, `__pycache__`, and `.pyc`; all paths listed externally |

No ignored files were found in the other requested categories
(`scientific-data`, `article-workspace`, `root-cache`, `temporary-backup`,
`editor-config`, `virtual-environment`, or `unknown-manual-review`). This says
only that no such ignored files are present in this checkout.

## Python/import inventory summary

| Measure | Count |
| --- | ---: |
| Tracked `.py` files under `src/`, `scripts/`, `tests/` | 181 |
| Runnable entry points | 130 |
| Compatibility wrappers | 7 |
| Reusable libraries | 20 |
| Tests | 24 |
| Import edges | 2,410 |
| Resolved repository edges | 800 |
| External edges | 1,595 |
| Unresolved best-effort edges | 14 |
| Dynamic edges | 1 |
| Runnable/compatibility files that are also imported | 40 |
| Runnable files with no importer but documented CLI/path reference | 96 |
| Possible orphan scripts after both checks | 1 |

The sole possible orphan candidate is
`scripts/analysis/thickness_mismatch/audits/audit_timoshenko_modes456_visualization.py`.
This is a manual-review flag only; its name and proximity to shape-construction
audits suggest scientific provenance that must be checked before any status or
move decision.

The most imported repository modules are
`formulas_thickness_mismatch.py` (50 importers), `formulas.py` (30),
`solvers.py` (24), `python_fem.py` (19), and the shared shape/tracking/
Timoshenko helpers (10–16). These are high-risk move targets.

## Results directory inventory

There are 23 immediate ignored result directories. The complete aggregate
table is `manifests/result_directory_summary.csv`; the per-file manifest is
stored only in the external backup.

- Closed-research evidence: `eb_epsilon_lower_envelope_step3a/`,
  `eb_rule_ab_exact_pareto/`, and `eb_rule_s_cost_break_even/`.
- Completed prerequisite diagnostics: epsilon pilot variants, corrected
  baseline thresholds, general completeness, and the branch-continuation
  gateway.
- Active diagnostics: EB/Timoshenko maps, mode-shape audits, fixed-epsilon and
  Stage-1 studies, and 3D-validation workspaces.
- Temporary smoke output: `_smoke/` (1,094 files, 205,691,561 bytes).

All 9,316 ignored result files were included in the external verified backup.

