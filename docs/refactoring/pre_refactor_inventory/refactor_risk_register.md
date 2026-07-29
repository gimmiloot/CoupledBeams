# Refactor risk register

Probability and severity are qualitative. Each high/critical row includes a
concrete gate for the next stage.

| Risk | Affected files | Consequence | Probability | Severity | Mitigation / required gate |
| --- | --- | --- | --- | --- | --- |
| Loss of ignored results | `results/` (9,316 files) | irrecoverable scientific evidence or expensive recomputation | high | critical | verified external per-file backup and restore spot-test before any move; preserve canonical reports with provenance |
| Loss of untracked scientific code | future untracked paths; current count 0 | silent loss outside Git bundle | medium | critical | re-run untracked/ignored manifest immediately before every structural stage; stop on new paths |
| Accidental determinant replacement | theory and analytic formula modules | wrong spectrum with superficially plausible output | medium | critical | freeze SHA256; symbolic entry-by-entry comparison; determinant spot checks and limits; explicit mathematical approval |
| Unknown ordering change | `equations.tex`, `formulas*.py`, shape reconstruction | eigenvector coefficients mapped to wrong fields | medium | critical | assert `Z=(A1,B1,A2,B2,P1,P2)` in tests/docs; matrix-column audit before refactor |
| Sorted frequency confused with descendant branch | branch tracker and plots | false crossings/veering/localization claims | high | critical | branch-ID regression fixtures; require `branch_id` and `current_sorted_index` as separate fields |
| Compatibility import breakage | 40 runnable-and-imported files; seven wrappers | documented commands/imports stop working | high | high | complete import-edge review; wrapper tests on old paths; no move until all 14 unresolved edges are classified |
| Negative results removed as “old” | veering assessments; safe-prefix closure/cost outputs | repeated failed research or incorrect positive claim | medium | critical | mark `closed-negative-result`; preserve report, input manifest, hashes, and decision wording |
| Diagnostic/article workflow confusion | scripts/analysis, absent paper workspaces | unreviewed diagnostic figure/claim enters article | medium | high | locate article workspaces; explicit promotion checklist from writing rules; separate output roots |
| Superseded thresholds reused | baseline/pilot/Rule A–D results | unsafe EB prefix or stale conclusion | high | critical | canonical-status registry; exact threshold source/hash in every consumer; reject legacy cache algorithm IDs |
| Root solver moved without multiplicity contract | `solvers.py`, completeness helpers, Timoshenko scans | missed close roots or collapsed cross-family multiplicity | high | critical | dedicated synthetic close-pair/multiplicity suite; K10/root-11 guards; SVD and cluster contract comparison |
| EB and Timoshenko helpers unified despite different assumptions | thickness/Timoshenko modules | hidden kinematic/constitutive model change | medium | critical | theory/code consistency audit; eta=0, thin limit, swap, cutoff, and energy tests before shared layer |
| Constitutive formulas duplicated | Timoshenko and solid-FEM diagnostics | divergent `G`, `kappa`, `A`, `I`, mass, or normalization | high | high | inventory all property constructors; define one reviewed data contract; compare dimensional units and outputs |
| Isotropic API frozen before anisotropic design | current scalar `E`, `EI`, `EA`, `GJ` code | future extension cannot express coupling/orientation cleanly | medium | high | decide orthotropy/general anisotropy and 2D/3D first; approve constitutive matrix/interface before implementation |
| Private-file leakage | `private/`, article workspaces | confidentiality breach | low | critical | local-only operations; exclude private paths from tracked manifests/navigation; manual content-owner gate |
| Backup placed inside repository | snapshot workflow | recursive bloat, accidental commit, self-copy | low | high | validate resolved backup path is outside Git root before writing; retain external-root assertion |
| Ignored article workspace archived without backup | absent but documented paper paths | loss of manuscript sources | medium | critical | discover actual workspace first; independent local manifest and owner approval before any move |
| Root caches overwritten by smoke test | result cache directories | provenance loss or contaminated “canonical” cache | medium | high | tests use OS temp only, `PYTHONDONTWRITEBYTECODE=1`, no default results output, pre/post hash check |
| FEM right-arm transform changed | `python_fem.py` and FEM diagnostics | incorrect orientation/stiffness/mass assembly | medium | critical | preserve `q_global=T@q_local`, `K=T K_local T.T`, `M=T M_local T.T`; run convention regression |
| Full 3D FEM implementation treated as analytic replica | 3D validation scripts/results | tuning to wrong physical joint model | medium | high | keep independent engineering benchmark role; mesh/joint/MAC audits and physical interpretation review |
| Generated-output path mistaken for broken source dependency | 371 absent inline refs, mostly results | unnecessary script/doc rewrite | high | medium | distinguish navigational Markdown links (0 broken) from absent expected generated outputs (345 results refs) |
| One possible orphan archived automatically | visualization audit | loss of one-off scientific provenance | medium | high | author review, search output signatures, execute only in a separately authorized later task |

The dominant refactor gates are therefore scientific, not stylistic: frozen
matrix/ordering checks, multiplicity/completeness regression, branch-identity
regression, and verified result/private backups.

